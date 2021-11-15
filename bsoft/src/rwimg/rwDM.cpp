/**
@file	rwDM.cpp
@brief	Functions for reading and writing Digital Micrograph files
@author Bernard Heymann
@date	Created: 20020619
@date 	Modified: 20170126
**/

#include "rwDM.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			readFixedDMHeader(ifstream* fimg, Bimage* p, int readdata);
int			readTagGroupData(ifstream* fimg, int dim_flag, Bimage* p, int readdata);
int			readTagGroupWithVersion(ifstream* fimg, Bimage* p, int readdata, int img_select);

int			show(0), keep(0), version(0), sb(0), endianness(1);
size_t		level(0);

/**
@brief	Reading a Digital Micrograph image file format.
@param	*p				the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int				error code (<0 means failure).
A 2D/3D image format used with CCD cameras in electron microscopy.
	File format extensions:  	.dm, .DM
	Two types: Fixed format (new) and the Macintosh format (old)
	Fixed format:
		Header size:				24 bytes (fixed).
		Byte order determination:	An endian flag: Must be 65535 or swap everything
		Data types: 				many.
	Macintosh format: 			Hermitian
		Header size:				8 bytes (fixed).
		Byte order determination:	Big-endian
		Data types: 				many.
**/
int 	readDM(Bimage* p, int readdata, int img_select)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	// Read a small block at the beginning to see if it is the fixed or tagged format
	char		buf[1024];
	fimg->read(buf, 4);
	if ( fimg->fail() ) return -2;
	
//	version = *((int *) buf);
	version = buf[0];

	if ( version < 1 || version > 4 ) version = buf[3];
	
	if ( show == 1 )
		cout << "Magic number = " << version << endl;
	
	fimg->seekg(0, ios::beg);
	
	switch ( version ) {
		case 0: readFixedDMHeader(fimg, p, readdata); break;
		case 3:
		case 4: readTagGroupWithVersion(fimg, p, readdata, img_select); break;
		default:
			cout << "Digital Micrograph format version " << version << " not supported!" << endl;
	}
	
	fimg->close();
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a Digital Micrograph map image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 2D/3D image format used in electron microscopy.
**/
int 	writeDM(Bimage* p)
{
 	DMhead* 	header = new DMhead;
	memset(header, 0, sizeof(DMhead));
	
	long	 	datatypesize = p->channels()*p->data_type_size();
	long		datasize = p->data_size();
	
	// Map the parameters
	header->endian = 65535;
	header->xSize = p->sizeX();
	header->ySize = p->sizeY();
	header->zSize = p->sizeZ();
	header->depth = datatypesize;
	switch ( p->data_type() ) {
        case UCharacter  : header->type = BINARY_DATA; break;
        case UShort : header->type = UNSIGNED_INT16_DATA; break;
        case Short  : header->type = SIGNED_INT16_DATA; break;
        case Integer    : header->type = SIGNED_INT32_DATA; break;
        case Float  : header->type = ( p->compound_type() == TComplex )? COMPLEX8_DATA: REAL4_DATA; break;
        case Double : header->type = ( p->compound_type() == TComplex )? COMPLEX16_DATA: REAL8_DATA; break;
        default : header->type = NULL_DATA;
    }
    
	p->data_offset(24);
	
    ofstream        fimg(p->file_name().c_str());
    if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, p->data_offset());
	fimg.write((char *)p->data_pointer(), datasize);
	
	fimg.close();
	
	delete header;
		
	return 0;
}


DataType	datatype_from_dm3_type(DMDataType dm3_type, Bimage* p)
{
	if ( show == 1 ) cout << "\tDMDataType=" << dm3_type << endl;
	
	p->compound_type(TSimple);
	p->channels(1);
	
	DataType		datatype = Unknown_Type;
	
	switch( dm3_type ) {
		case BINARY_DATA :         datatype = UCharacter; break;
		case UNSIGNED_INT8_DATA :  datatype = UCharacter; break;
		case SIGNED_INT8_DATA :    datatype = SCharacter; break;
		case RGB_DATA :
		case RGB_UINT8_DATA :
			datatype = UCharacter;
			p->channels(3);
			p->compound_type(TRGB);
			break;
		case OS_RGBA_UINT8_DATA :
			datatype = UCharacter;
			p->channels(4);
			p->compound_type(TRGBA);
			break;
		case UNSIGNED_INT16_DATA : datatype = UShort; break;
		case RGB_UINT16_DATA :
			datatype = UShort;
			p->channels(3);
			p->compound_type(TRGB);
			break;
		case RGBA_UINT16_DATA :
			datatype = UShort;
			p->channels(4);
			p->compound_type(TRGBA);
			break;
		case SIGNED_INT16_DATA :   datatype = Short; break;
		case UNSIGNED_INT32_DATA : datatype = UInteger; break;
		case SIGNED_INT32_DATA :   datatype = Integer; break;
		case REAL4_DATA :          datatype = Float; break;
		case RGBA_FLOAT32_DATA :
			datatype = Float;
			p->channels(4);
			p->compound_type(TRGBA);
			break;
		case REAL8_DATA :          datatype = Double; break;
		case RGB_FLOAT64_DATA :
			datatype = Double;
			p->channels(3);
			p->compound_type(TRGB);
			break;
		case RGBA_FLOAT64_DATA :
			datatype = Double;
			p->channels(4);
			p->compound_type(TRGBA);
			break;
		case COMPLEX8_DATA :       datatype = Float; p->compound_type(TComplex); break;
		case COMPLEX16_DATA :      datatype = Double; p->compound_type(TComplex); break;
		default : datatype = UCharacter;
	}
	
	p->data_type(datatype);
	if ( p->compound_type() == TComplex ) p->channels(2);
	
	return datatype;
}



int			readFixedDMHeader(ifstream* fimg, Bimage* p, int readdata)
{
	DMhead*		header = new DMhead;
	DMMachead*	macheader = (DMMachead *) header;
	
//	if ( fread( header, sizeof(DMhead), 1, fimg ) < 1 ) return -2;
	fimg->read((char *)header, sizeof(DMhead));
	if ( fimg->fail() ) return -2;
	
	if ( verbose & VERB_DEBUG || show )
		cout << "DEBUG readFixedDMHeader: macheader=" << macheader << 
				" width=" << macheader->width << " height=" << macheader->height << endl;
	
	p->data_offset(8);
	if ( ( macheader->width < 1 ) || ( macheader->height < 1 ) ) {
		macheader = (DMMachead *) ((char*)header + 6);
		p->data_offset(p->data_offset() + 6);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readFixedDMHeader: macheader=" << macheader << 
				" width=" << macheader->width << " height=" << macheader->height << endl;
	
	// Determine header type
    int     	i, fixed(1), sb(0);
/*	if ( macheader->width == 0 && macheader->height == -1 ) {
		fixed = 1;
		sb = 0;
	} else if ( macheader->width == -1 && macheader->height == 0 ) {
		fixed = 1;
		sb = 1;
	} else if ( systype(0) >= LittleIEEE ) sb = 1;
*/	
    // Swap bytes if necessary
    unsigned char*   	b = (unsigned char *) header;
    if ( sb ) {
		if ( fixed ) {
	    	if ( verbose & VERB_PROCESS )
		    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    		for ( i=0; i<24; i+=4 ) swapbytes(b+i, 4);
		} else {
			b = (unsigned char *) macheader;
	    	if ( verbose & VERB_PROCESS )
		    	cerr << "Warning: Swapping header byte order for 2-byte types" << endl;
    		for ( i=0; i<8; i+=2 ) swapbytes(b+i, 2);
		}
    }
    
	// Map the parameters
	p->images(1);
	p->channels(1);
	p->compound_type(TSimple);
	if ( fixed ) {
		p->size(header->xSize, header->ySize, header->zSize);
		p->data_offset(24);
		datatype_from_dm3_type(header->type, p);
	} else {
		p->size(macheader->width, macheader->height, 1);
    	switch( macheader->type ) {
        	case 6 :  p->data_type(UCharacter); break;	
        	case 14 : p->data_type(UCharacter); break;
        	case 9 :  p->data_type(SCharacter); break;
        	case 10 : p->data_type(UShort); break;
        	case 1 :  p->data_type(Short); break;
        	case 7 :  p->data_type(Integer); break;
        	case 2 :  p->data_type(Float); break;
        	case 3 :  p->data_type(Float); p->compound_type(TComplex); p->channels(2); break;
        	case 13 : p->data_type(Double); p->compound_type(TComplex); p->channels(2); break;
        	default : p->data_type(UCharacter);
		}
	}
	
	
	delete header;

	if ( readdata ) {
		fimg->read((char *)p->data_pointer(), p->alloc_size());
		if ( fimg->fail() ) return -3;
		if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	
	return sb;
}

int			dm3_type_length(int dm3_type)
{
	int			dtlen = 0;
	
	switch ( dm3_type ) {
		case 2: dtlen = 2; break;
		case 3: dtlen = 4; break;
		case 4: dtlen = 2; break;
		case 5: dtlen = 4; break;
		case 6: dtlen = 4; break;
		case 7: dtlen = 8; break;
		case 8: dtlen = 1; break;
		case 9: dtlen = 1; break;
		case 10: dtlen = 1; break;
		case 11: dtlen = 8; break;
		case 12: dtlen = 8; break;
		default:
			cerr << "Error: Data type " << dm3_type << " length not defined!" << endl;
	}
	
	return dtlen;
}

double		dm3_value(ifstream* fimg, int dm3_type)
{
	int			dtlen = dm3_type_length(dm3_type);
	if ( dtlen < 1 ) return 0;
	
	long		i, ivalue(0);
	double		dvalue(0);
	TypePointer	buf;
	buf.uc = new unsigned char[1024];

	fimg->read((char *)buf.uc, dtlen);

	if ( sb - endianness ) swapbytes(buf.uc, dtlen);
	
	switch ( dm3_type ) {
		case 2: ivalue = *buf.ss; break;
		case 3: ivalue = *buf.si; break;
		case 4: ivalue = *buf.us; break;
		case 5: ivalue = *buf.ui; break;
		case 6: dvalue = *buf.f; break;
		case 7: dvalue = *buf.d; break;
		case 8: ivalue = *buf.uc; break;
		case 9: ivalue = *buf.uc; break;
		case 10: ivalue = *buf.uc; break;
		case 11: ivalue = *buf.sl; break;
		case 12: ivalue = *buf.ul; break;
		default:
			cerr << "Error: Data type " << dm3_type << " not defined!" << endl;
	}
	
	if ( !dvalue ) dvalue = ivalue;
	
	if ( show == 1 ) cout << tab << dvalue;
	
	if ( show == 2 ) {
		for ( i=0; i<level; i++ ) cout << tab;
		if ( dm3_type == 6 || dm3_type == 7 )
			cout << "<real>" << dvalue << "</real>" << endl;
		else
			cout << "<integer>" << ivalue << "</integer>" << endl;
	}

	return dvalue;
}

unsigned long	dm_read_integer(ifstream* fimg, long len)
{
	unsigned short		sval(0);
	unsigned int		ival(0);
	unsigned long		lval(0);
	TypePointer	buf;
	buf.uc = new unsigned char[1024];
	
	if ( version < 4 && len > 4 ) len = 4;
	
	fimg->read((char *)buf.uc, len);
	
	if ( sb ) swapbytes(buf.uc, len);
	
	switch ( len ) {
		case 2: sval = *buf.us; lval = sval; break;
		case 4: ival = *buf.ui; lval = ival; break;
		case 8: lval = *buf.ul; break;
		default: lval = *buf.ul;
	}
	
	return lval;
}

int			tag_convert(unsigned char* tag)
{
	for ( ; *tag; tag++ ) if ( *tag > 127 ) *tag -= 64;
	
	return 0;
}

/*
	Image data is assigned to the data pointer based on the keep flag
	Associated properties are also kept
*/
int			readTag(ifstream* fimg, int dim_flag, Bimage* p, int readdata, int& notag)
{
	unsigned char	t, tag[128];
	unsigned short	len;
	size_t			i, k, data_size(0);
	size_t			nnum, size(0), narr, nel, dt, dtarr[128], dtlen;	// Change in v4
	double			val(0), arr[256];
	char			buf[1024];
	
	level++;
	
	// Common part of TagSubGroup and TagEntry records
	// Selector (1 byte): TagSubGroup=20, TagEntry=21
	// Length of label (2 bytes)
	// Label
	*fimg >> t;
	len = dm_read_integer(fimg, sizeof(unsigned short));
//	cout << "\tlen=" << len << endl;

	if ( len > 128 ) {
		cout << "\tlen=" << len << endl;
		cerr << "Error: tag length too long!" << endl;
		bexit(-1);
	}
	
	if ( len ) {
		fimg->read((char *)tag, len);
		tag[len] = 0;
		tag_convert(tag);
		notag = 0;
	} else {
		snprintf((char *)tag, 128, "%d", notag);
		notag++;
	}

	// Image list
	if ( strcmp((char *)tag, "ImageList") == 0 ) keep = 1;

	// Start of a new image specification
	if ( keep == 1 && p->data_offset() > 0 ) {
		if ( strcmp((char *)tag, "ImageTags") == 0 || strcmp((char *)tag, "ImageData") == 0 ) {
			p->next = p->copy_header();
			p = p->next;
			p->images(1);
			keep = 2;
		}
	}
	
	if ( show == 1 ) {
		cout << endl;
		for ( i=0; i<level; i++ ) cout << tab;
		cout << "Tag=" << (int)t << "\tlen=" << len << tab << tag;
	}
	if ( show == 2 ) {
		for ( i=0; i<level; i++ ) cout << tab;
		cout << "<key>" << tag << "</key>" << endl;
	}

	if ( version == 4 )	{	// Changed in v4
		size = dm_read_integer(fimg, sizeof(size_t));	// Size of data type & data
		if ( show == 1 )
			cout << "\tsize=" << size;
	}
	
	if ( t == 20 ) {			// TagSubGroup record
		if ( strcmp((char *)tag, "Dimensions") == 0 ) dim_flag = 1;
		readTagGroupData(fimg, dim_flag, p, readdata);
	} else if ( t == 21 ) {		// TagEntry record
		fimg->read((char *)buf, 4);		// String = "%%%%"
		buf[4] = 0;
		nnum = dm_read_integer(fimg, sizeof(size_t));	// Number of integers in array
		if ( show == 1 ) cout << tab << setw(4) << buf << "\tnnum=" << nnum;
		
		for ( i=0; i<nnum; i++ ) {
			dtarr[i] = dm_read_integer(fimg, sizeof(size_t));	// Data type array
			if ( show == 1 ) cout << "\tdt[" << i << "]=" << dtarr[i];
		}
		if ( show == 1 ) cout << endl;

		dt = dtarr[0];
		if ( dt <= 12 ) {
			val = dm3_value(fimg, dt);
		} else if ( dt == 15 ) {										// Struct
			narr = dtarr[2];
			nel = 1;
			if ( show == 1 )
				cout << "\tstruct=" << narr << "x" << nel;
			if ( show == 2 ) {
				for ( i=0; i<level; i++ ) cout << tab;
				cout << "<array>" << endl;
			}
			level++;
			for ( k=0; k<narr; k++ )
				arr[k] = dm3_value(fimg, dtarr[4+2*k]);
			level--;
			if ( show == 2 ) {
				for ( i=0; i<level; i++ ) cout << tab;
				cout << "</array>" << endl;
			}
		} else if ( dt == 18 ) {	// String
			fimg->read((char *)buf, dtarr[1]);
			if ( show == 1 ) cout << "\tstrlen=" << dtarr[1];
			if ( show == 2 ) {
				for ( i=0; i<level; i++ ) cout << tab;
				cout << "<string>" << buf << "</string>" << endl;
			}
		} else if ( dt == 20 ) {	// Data array
			if ( dtarr[1] == 15 ) narr = dtarr[3];
			else narr = 1;
			nel = dtarr[nnum-1];		// Number of elements
			if ( show == 1 )
				cout << "\tarr=" << narr << "x" << nel;
			if ( show == 2 ) {
				for ( i=0; i<level; i++ ) cout << tab;
				cout << "<array>" << endl;
			}
			level++;
			if ( dtarr[1] == 15 ) {	// Struct
				if ( show == 1 )
					cout << "\tstruct=" << narr << "x" << nel;
				if ( show == 2 ) {
					for ( i=0; i<level; i++ ) cout << tab;
					cout << "<array>" << endl;
				}
				level++;
				for ( i=0; i<nel; i++ )
					for ( k=0; k<narr; k++ )
						arr[k] = dm3_value(fimg, dtarr[5+2*k]);
				level--;
				if ( show == 2 ) {
					for ( i=0; i<level; i++ ) cout << tab;
					cout << "</array>" << endl;
				}
			} else {
				narr = 1;
				dtlen = dm3_type_length(dtarr[1]);
				data_size = nel*dtlen;
				// The big image data block is not printed for the plist format
				if ( strcmp((char *)tag, "Data") == 0 ) {
//					if ( nel > 1e6 ) {
//						keep = 1;
//						p->data_offset(fimg->tellg());
//					}
//					if ( show == 1 )
//						cout << "\tdata_size=" << data_size;
//					if ( readdata ) {
//						p->data_alloc(data_size);
//						fimg->read((char *)p->data_pointer(), data_size);
//					} else {
//						fimg->seekg(data_size, ios_base::cur);
//					}
//					fimg->seekg(data_size, ios_base::cur);
				} else {
					for ( i=0; i<nel; i++ ) {
						if ( i < 10 ) val = dm3_value(fimg, dtarr[1]);	// Element value
						else fimg->read((char *)buf, dtlen);			// Element value
					}
				}
			}
			level--;
			if ( show == 2 ) {
				for ( i=0; i<level; i++ ) cout << tab;
				cout << "</array>" << endl;
			}
		} else {
			cerr << "Error: Data type " << dt << " not defined!" << endl;
		}
		
//		if ( keep && strcmp((char *)tag, "ImageData") == 0 ) keep = 0;

		if ( keep ) {
//			show = 1;
//			if ( show == 1 ) cout << "\tdim_flag=" << dim_flag;
			if ( strcmp((char *)tag, "Data") == 0 ) {
				p->data_offset(fimg->tellg());
				if ( show == 1 )
					cout << "\tdata_size=" << data_size;
				fimg->seekg(data_size, ios_base::cur);
			}
			if ( strcmp((char *)tag, "DataType") == 0 ) {
				datatype_from_dm3_type((DMDataType)val, p);
				if ( show == 1 )
					cout << "\tdatatype=" << p->data_type() << "\tcompoundtype=" << 
						p->compound_type() << "\tc=" << p->channels() << "\tval=" << val;
			}
			if ( strcmp((char *)tag, "Minimum Value (counts)") == 0 ) {
				p->minimum(val);
			}
			if ( strcmp((char *)tag, "Maximum Value (counts)") == 0 ) {
				p->maximum(val);
			}
			if ( dim_flag == 1 ) {
				p->sizeX((long) val);
				dim_flag = 2;
			} else if ( dim_flag == 2 ) {
				p->sizeY((long) val);
				dim_flag = 3;
			} else if ( dim_flag == 3 ) {
				p->images((long) val);
				dim_flag = 0;
				if ( show == 1 ) cout << "\tx=" << p->sizeX() << "\ty=" << 
					p->sizeY() << "\tn=" << p->images();
			}
			if ( strcmp((char *)tag, "PixelDepth") == 0 ) {
				if ( show == 1 ) cout << tab << tag << "=" << val << endl;
			}
//			if ( strcmp((char *)tag, "Scale") == 0 ) {
//				p->sampling(val, val, 1);
//				if ( show == 1 ) cout << tab << tag << "=" << p->sampling(0) << endl;
//				cout << tab << tag << "=" << p->sampling(0) << endl;
//			}
			if ( strcmp((char *)tag, "Pixel Size (um)") == 0 ) {
				p->sampling(arr[0], arr[1], 1);
				if ( show == 1 ) cout << tab << tag << "=" << p->sampling(0) << endl;
//				cout << tab << tag << "=" << p->sampling(0) << endl;
			}
			if ( strcmp((char *)tag, "Actual Magnification") == 0 ) {
				if ( val ) {
					p->show_scale(val/1e4);
					if ( show == 1 ) cout << tab << tag << "=" << val << endl;
//					cout << tab << tag << "=" << val << endl;
				}
			}
			if ( strcmp((char *)tag, "Pixel Upsampling") == 0 ) {
				if ( show == 1 ) cout << tab << tag << "=" << arr[0] << "x" << arr[1] << endl;
			}
			if ( strcmp((char *)tag, "SourceSize_Pixels") == 0 ) {
				if ( show == 1 ) cout << tab << tag << "=" << arr[0] << "x" << arr[1] << endl;
			}
//			if ( strcmp(tag, "ImageIndex") == 0 ) {
//				if ( show == 1 ) cout << tab << tag << "=" << val << endl;
//			}
			if ( strstr((char *)tag, "Emission") ) {
				if ( show == 1 ) {
					cout << tab << tag << "=" << val << endl;
//					for ( i=0; i<len; i++ ) cout << (int)tag[i] << tab << tag[i] << endl;
				}
			}
//			show = 0;
		}
//	} else {
//		cerr << "Error: Undefined tag type! (" << t << ")" << endl;
	}
	
	level--;

	return 0;
}

	// TagGroupData contains 3 items
	// Sorted flag: 0,1; v3-v4 (1 byte)
	// Open flag: 0,1; v3-v4 (1 byte)
	// Number of records (TagEntry/TagSubGroup): v3 (4 byte), v4 (8 byte)
	// Content (TagGroupData)
int			readTagGroupData(ifstream* fimg, int dim_flag, Bimage* p, int readdata)
{
	unsigned char	sorted, open;
	size_t			i, ntag(0);
	int				notag(0);
	
	*fimg >> sorted;
	*fimg >> open;
	
	ntag = dm_read_integer(fimg, sizeof(size_t));	// Number of tags
	
	if ( show == 1 ) cout << "\tntag=" << ntag;
	if ( show == 2 ) {
		for ( i=0; i<level; i++ ) cout << tab;
		cout << "<dict>" << endl;
	}
	
	for ( i=0; i<ntag; i++ ) {
		readTag(fimg, dim_flag, p, readdata, notag);
		if ( dim_flag ) dim_flag++;
	}

	if ( show == 2 ) {
		for ( i=0; i<level; i++ ) cout << tab;
		cout << "</dict>" << endl;
	}
	
	return 0;
}

/*
	Main DM3/4 file database:
	
	Arrangement:
		The file header (TagGroupWithVersion) contains 4 items:
		Version: 1-4 (4 byte)
		Size of contained data (TagGroupData): v1-v3 (4 byte), v4 (8 byte)
		Endianness: 0,1 (4 byte)
		Content (TagGroupData)
	
	Important image information:
		Main block: "ImageList"
			May contain more than one image, usually a small display version
			followed by the real data
		Important elements:
			Pixel size (nm):		ImageData/Calibrations/Dimension/Scale
			Data element type:		ImageData/DataType
			Image size:				ImageData/Dimensions
			Data element size:		ImageData/PixelDepth
			Camera pixel size (um):	ImageTags/Acquisition/Device/CCD/Pixel Size (um)
			Dose rate:				ImageTags/Calibration/Dose rate/Calibration
			Magnification:			ImageTags/Microscope Info/Actual Magnification
			Acceleration voltage:	ImageTags/Microscope Info/Voltage
*/
int			readTagGroupWithVersion(ifstream* fimg, Bimage* p, int readdata, int img_select)
{
	if ( verbose & VERB_DEBUG_DM ) show = 2;
	
	keep = level = 0;
	
	size_t			file_length(0), val(0);
	char			buf[1024];
	int*			vp = (int *) buf;
	
	fimg->read(buf, 4);
	if ( fimg->fail() ) return -2;
	
//	version = *((int *) buf);
	version = *vp;

	if ( version > 100 ) {
		sb = 1;
		version = buf[3];
	}

	file_length = dm_read_integer(fimg, sizeof(size_t));	// File length
	file_length += 16;
	
	endianness = dm_read_integer(fimg, sizeof(unsigned int));	// Endianness

	if ( show == 1 ) {
		if ( sb ) cout << "Warning: Swapping header bytes" << endl;
		cout << "Version: " << version << endl;
		cout << "File length: " << file_length << endl;
		cout << "Endianness: " << endianness << endl;
	}
	
	if ( file_length <= 16 ) {
		cout << "Error: file length = " << file_length << endl;
		return error_show("File length specifier incorrect!\n", __FILE__, __LINE__);
	}
	
	if ( show == 2 ) {
		cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    	cout << "<!DOCTYPE plist PUBLIC \"-//Apple Computer//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">" << endl;
		cout << "<plist version=\"1.0\">" << endl;
	}

	// Set up default parameters
	p->size(1, 1, 1);
	p->channels(1);
	p->images(1);
	
	readTagGroupData(fimg, 0, p, readdata);
	
	// End of file
	val = dm_read_integer(fimg, sizeof(unsigned int));	// 500
	val = dm_read_integer(fimg, sizeof(unsigned int));	// Label size
	if ( val ) {
		fimg->read(buf, val);
		if ( show == 2 ) cout << "<" << buf << "/>" << endl;
	}
	
	if ( show == 1 ) cout << endl;
	if ( show == 2 ) cout << "</plist>" << endl;
	
	int				i;
	Bimage*			p1 = NULL;
	if ( verbose & VERB_FULL ) {
		for ( p1 = p, i = 0; p1; p1 = p1->next, ++i )
			cout << "Image " << i << " dimensions:             " << 
				p1->channels() << tab << p1->size() << tab << p1->images() << endl;
		cout << endl;
	}
	
	if ( p->next ) {
		if ( p->next->sizeX() > p->sizeX() ) *p = *(p->next);
		delete p->next;
		p->next = NULL;
	}
	
	if ( img_select >= 0 ) {
		p->data_offset(p->data_offset()+img_select*p->sizeX()*p->sizeY()*p->sizeZ()*p->channels()*p->data_type_size());
		p->images(1);
	}
	
	if ( p->show_scale() != 1 ) {
		p->sampling(p->sampling(0)[0]/p->show_scale(), p->sampling(0)[1]/p->show_scale(), 1);
		p->show_scale(1);
	}

	if ( readdata ) {
//		fimg->seekg(p->data_offset(), ios_base::beg);
//		p->data_alloc();
//		fimg->read((char *)p->data_pointer(), p->alloc_size());
//		if ( version == 4 ) sb = 0;
		sb -= endianness;
//		cout << "sb = " << sb << endl;
		p->read_data(fimg, 0, sb, 0, 0);
	}
	
	return version;
}


