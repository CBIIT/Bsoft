/**
@file	rwMRC.cpp
@brief	Functions for reading and writing MRC files
@author Bernard Heymann
@date	Created: 19990321
@date 	Modified: 20180712
**/

#include "rwMRC.h"
#include "rwCCP4.h"
#include "file_util.h"
#include "utilities.h"

#ifndef NOSYMOP
#include "rwsymop.h"
#endif

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

struct SERIhead {
	short		tilt;		// 100 x tilt angle in degrees
	short		stageX;		// stage coordinate x 25
	short		stageY;		// stage coordinate x 25
	short		mag;		// magnfication / 100
	short		intensity;	// intensity x 25000
	short		dose1;		// dose as two shorts
	short		dose2;
};

/**
@brief	Reading a MRC map image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@return	int			error code (<0 means failure).
A 3D image format used in electron microscopy.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.mrc
	The identifier is a 4-byte machine stamp (same as for CCP4 maps):
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = signed byte, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
**/
int 	readMRC(Bimage* p, int readdata, int img_select)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;
//	cerr << "MRC file = " << p->file_name() << endl;

	MRChead*		header = new MRChead;
	
	fimg->read((char *)header, MRCSIZE);
	if ( fimg->fail() ) return -2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readMRC: File opened" << endl;

	string			exttype(header->exttype, 0, 4);
	if ( exttype.length() < 1 || exttype[0] == ' ' )
		if ( strncmp(header->labels, "Seri", 4) == 0 ) exttype = "SERI";	// SerialEM
	
//	cout << header->labels << endl;	
//	cout << "length of exttype: " << exttype.length() << endl;

	if ( verbose & VERB_PROCESS ) {
		cout << "Header version:                 " << header->nversion << endl;
		cout << "Extended header type:           " << exttype << " (" << header->nsymbt << ")" << endl;
	}
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    long			i, sb(0);
    if ( abs(header->mode) > SWAPTRIG || abs(header->mapc) > SWAPTRIG ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		long		extent = MRCSIZE - 800; // exclude labels from swapping
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
    
	if ( header->nx < 1 || header->ny < 1 || header->nz < 1 || header->mode < 0 ) {
		cerr << "Error readMRC: The file " << p->file_name() << " is not a valid MRC image!" << endl;
		cerr << "\tnx = " << header->nx << " ny = " << header->ny << " nz = " << header->nz
			<< " mode = " << header->mode << endl;
		fimg->close();
		return -3;
	}
	
    // Convert VAX floating point types if necessary
/*    if ( header->amin > header->amax ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Converting VAX floating point" << endl;
 		for ( i=40; i<64; i+=4 )		// Unit cell parameters
    	    vax2ieee(b+i, sb);
		for ( i=76; i<88; i+=4 )		// Min, max, mean
    	    vax2ieee(b+i, sb);
		for ( i=96; i<220; i+=4 )		// Extra and origin
    	    vax2ieee(b+i, sb);
    }
*/
	string			ext(extension(p->file_name()));
	
	long			nimg(0);
    
	// Old versions: pre-2015
	if ( header->nversion < 20140 || header->nversion > 20209 ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readMRC: ---old version---" << endl;
		if ( strncmp(header->map, "STK ", 4) == 0 || ext == "st" || ext == "stk" ||
//				ext.contains("ali") || ext == "mrcs" || header->ispg < 1 ) {
				ext.find("ali") != string::npos || ext == "mrcs" || exttype == "SERI" ) {
//				ext.contains("ali") || ext == "mrcs" ) {
			nimg = header->nz;
			p->size(header->nx, header->ny, 1);
		} else if ( header->mz && header->mz < header->nz ) {
			nimg = header->nz/header->mz;
			p->size(header->nx, header->ny, header->mz);
		} else {
			nimg = 1;
			p->size(header->nx, header->ny, header->nz);
		}
	}
	// New versions: 2015+
	else {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readMRC: ---new version---" << endl;
		if ( header->nz%header->mz ) {
			cerr << "Error: The image " << p->file_name() << " incorrectly specifies the number of slices and volumes!" << endl;
			cerr << tab << "NZ = " << header->nz << "  MZ = " << header->mz << endl;
		}
		if ( exttype == "SERI" ) {
			p->size(header->nx, header->ny, 1);
			nimg = header->nz;
		} else {
			p->size(header->nx, header->ny, header->mz);
			nimg = header->nz/header->mz;
		}
	}

	p->channels(1);
	
	// Allocating the single sub-image and setting its origin
	if ( img_select > -1 ) {
		if ( img_select >= nimg ) img_select = nimg - 1;
		nimg = 1;
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readMRC: number of images = " << nimg << endl;
	p->images(nimg);
	
	switch ( header->mode ) {
		case 0:
			if ( header->imodFlags & 1 ) p->data_type(SCharacter);
			else p->data_type(UCharacter);
//			if ( ext == "mrc" && exttype != "SERI" ) p->data_type(SCharacter);
			break;
		case 1: p->data_type(Short); break;
		case 2: p->data_type(Float); break;
		case 3: p->data_type(Short); p->compound_type(TComplex); p->channels(2); break;
		case 4: p->data_type(Float); p->compound_type(TComplex); p->channels(2); break;
		case 6: p->data_type(UShort); break;
		case 7: p->data_type(Integer); break;
		case 16: p->data_type(UCharacter); p->compound_type(TRGB); p->channels(3); break;
		default: p->data_type(UCharacter); break;
	}

	p->minimum(header->amin);
	p->maximum(header->amax);
	p->average(header->amean);
	p->standard_deviation(header->arms);

	UnitCell	uc(header->a, header->b, header->c, header->alpha, header->beta, header->gamma);
	p->unit_cell(uc);
	
	if ( header->mx && header->my && header->mz )
		p->sampling(header->a/header->mx, header->b/header->my, header->c/header->mz);
	if ( header->ispg >= 400 ) p->space_group(header->ispg-400);
	else p->space_group(header->ispg);
	p->label(header->labels);

	p->origin(-header->xOrigin/p->sampling(0)[0], -header->yOrigin/p->sampling(0)[1], -header->zOrigin/p->sampling(0)[2]);	// New header

	if ( header->nsymbt > 0 ) {
		char*			ext_header = new char[header->nsymbt];
		fimg->read((char *)ext_header, header->nsymbt);
		if ( fimg->fail() ) return -1;
		if ( exttype == "JSON" ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readMRC: reading JSON header" << endl;
			string		file_name(p->file_name());	// JSON header may contain old file name
			JSparser	parser(ext_header, header->nsymbt);
			p->meta_data() = parser.parse();
			if ( img_select > -1 ) p->meta_data_retain_one_image(img_select);
			p->update_from_meta_data();
			p->file_name(file_name);
		} else if ( exttype == "SERI" ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readMRC: reading SerialEM header" << endl;
			SERIhead*	serihead = (SERIhead *) ext_header;
			JSvalue		arr(JSarray);
			JSvalue		img(JSobject);
			for ( i=0; i<header->nz; ++i ) if ( img_select < 0 || img_select == i ) {
				img["Tilt"] = serihead[i].tilt/100.0;
				img["StageX"] = serihead[i].stageX/25.0;
				img["StageY"] = serihead[i].stageY/25.0;
				img["Magnification"] = serihead[i].mag*100.0;
				img["Intensity"] = serihead[i].intensity/25000.0;
				img["Dose"] = serihead[i].dose1;
				arr.push_back(img);
			}
			p->meta_data()["image"] = arr;
//			cout << p->meta_data() << endl;
		}
		if ( ext_header ) delete[] ext_header;
	}
	
//	cout << "readMRC:" << endl << p->meta_data() << endl;
	
	p->data_offset(MRCSIZE + header->nsymbt);
	if ( header->mode%5 > 2 && header->mode%5 < 5 ) {
		fimg->seekg(0, ios::end);
		if ( (double) fimg->tellg() > p->data_offset() + 0.8*p->alloc_size() )
			p->sizeX(2*(p->sizeX() - 1));
		if ( header->mx%2 == 1 ) p->sizeX(p->sizeX() + 1);	// Quick fix for odd x-size maps
		fimg->clear();
		fimg->seekg(0, ios::beg);
	}
	
	delete header;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwMRC: Header setup done" << endl;

	long			readsize = p->alloc_size();
	unsigned char*	data = NULL;
	
	if ( img_select < 0 ) img_select = 0;
//	cerr << "MRC file = " << p->file_name() << endl;
	
	if ( readdata ) {
		p->data_alloc();
		if ( p->compound_type() == TSimple ) {
			fread_large(p->data_pointer(), readsize, p->data_offset() + img_select*readsize, fimg);
			if ( sb ) swapbytes(readsize, p->data_pointer(), p->data_type_size());
		} else {
			readsize = p->channels()*(p->sizeX()/2+1)*p->sizeY()*p->sizeZ()*p->images()*p->data_type_size();
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG rwMRC: readsize: " << readsize << endl;
			data = new unsigned char[readsize];
			fread_large(data, readsize, p->data_offset() + img_select*readsize, fimg);
			if ( sb ) swapbytes(readsize, data, p->data_type_size());
			p->unpack_transform(data, CentHerm);
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a MRC map image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D image format used in electron microscopy.
**/
int 	writeMRC(Bimage* p)
{
	p->color_to_simple();

	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
	string			ext(extension(p->file_name()));
    
    switch ( p->data_type() ) {
		case Bit: case UCharacter:
			p->change_type(UCharacter); break;
		case SCharacter:
			if ( ext == "mrc" ) p->change_type(SCharacter);
			if ( ext != "mrc" ) p->change_type(UCharacter);
			break;
    	case UShort: case Short:
			p->change_type(Short); break;
    	default:
			p->change_type(Float); break;
    }
	
	MRChead*	header = new MRChead;
	memset(header, 0, sizeof(MRChead));
	
	// Map the parameters
	UnitCell	uc = p->unit_cell();
//	strncpy(header->map, "MAP ", 4);
	memcpy(header->map, "MAP ", 4);
	header->nversion = 20140;
	set_CCP4_machine_stamp(header->machst);
	header->nx = p->sizeX();
	header->ny = p->sizeY();
	header->nz = p->images() * p->sizeZ();
	if ( p->compound_type() == TComplex ) header->nx = p->sizeX()/2 + 1;	// If a transform, physical storage is nx/2 + 1
	switch ( p->data_type() ) {
		case UCharacter: header->mode = 0; header->imodFlags = 0; break;
		case SCharacter: header->mode = 0; header->imodFlags |= 1; break;
		case UShort: header->mode = 6; break;
		case Short: header->mode = 1; if ( p->compound_type() == TComplex ) header->mode = 3; break;
		case Float: header->mode = 2; if ( p->compound_type() == TComplex ) header->mode = 4; break;
		default: header->mode = 0; break;
	}
	if ( p->image->origin()[0] >= 0 )
		header->nxStart = (int) -(p->image->origin()[0] + 0.5);
	else
		header->nxStart = (int) -(p->image->origin()[0] - 0.5);
	if ( p->image->origin()[1] >= 0 )
		header->nyStart = (int) -(p->image->origin()[1] + 0.5);
	else
		header->nyStart = (int) -(p->image->origin()[1] - 0.5);
	if ( p->image->origin()[2] >= 0 )
		header->nzStart = (int) -(p->image->origin()[2] + 0.5);
	else
		header->nzStart = (int) -(p->image->origin()[2] - 0.5);
	header->mx = (int) (uc.a()/p->sampling(0)[0] + 0.5);
	header->my = (int) (uc.b()/p->sampling(0)[1] + 0.5);
	header->mz = (int) (uc.c()/p->sampling(0)[2] + 0.5);
	header->mz = p->sizeZ();
	header->mapc = 1;
	header->mapr = 2;
	header->maps = 3;
	header->amin = p->minimum();
	header->amax = p->maximum();
	header->amean = p->average();
	header->arms = p->standard_deviation();
	header->a = uc.a();
	header->b = uc.b();
	header->c = uc.c();
//	header->xOrigin = p->image->origin()[0];
//	header->yOrigin = p->image->origin()[1];
//	header->zOrigin = p->image->origin()[2];
	header->xOrigin = -p->image->origin()[0]*p->sampling(0)[0];
	header->yOrigin = -p->image->origin()[1]*p->sampling(0)[1];
	header->zOrigin = -p->image->origin()[2]*p->sampling(0)[2];
	
	// This is a band-aid to overcome the limitations of the image format
	if ( fabs(uc.a() - p->sampling(0)[0]*header->mx) > 0.001 || fabs(uc.b() - p->sampling(0)[1]*header->my) > 0.001 ||
			fabs(uc.c() - p->sampling(0)[2]*header->mz) > 0.001 ) {
		header->a = p->sampling(0)[0]*header->mx;
		header->b = p->sampling(0)[1]*header->my;
		header->c = p->sampling(0)[2]*header->mz;
		if ( verbose & VERB_FULL )
			cerr << "Warning: Resetting the unit cell to: " << header->a << "," << header->b << "," << header->c << endl;
	}
	
	header->alpha = uc.alpha()*180/M_PI;
	header->beta = uc.beta()*180/M_PI;
	header->gamma = uc.gamma()*180/M_PI;
	if ( p->sizeZ() > 1 ) {
		header->ispg = p->space_group();
		if ( header->ispg < 1 ) header->ispg = 1;
		if ( p->images() > 1 ) header->ispg += 400;
	}
	
	int				nsym(0), jflag(1);
	Bstring			temp;
	char* 			symop = NULL;
	stringstream    ss;

	if ( jflag ) {
		ss << p->meta_data() << endl;
		header->nsymbt = ss.tellp();
//		cout << "Metadata size: " << header->nsymbt << endl;
//		strncpy(header->exttype, "JSON", 4);
		memcpy(header->exttype, "JSON", 4);
	} else {
#ifndef NOSYMOP
		if ( p->space_group() > 0 ) {
			symop = read_symop(temp, p->space_group(), nsym);
//			strncpy(header->exttype, "CCP4", 4);
			memcpy(header->exttype, "CCP4", 4);
		}
#endif
		header->nsymbt = nsym*80;
	}
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwMRC: nsymbt = " << header->nsymbt << endl;

	long 			i;
	
	header->nlabl = 10;
	strncpy(header->labels, p->label().c_str(), 799);
	header->labels[799] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG rwMRC: nlabl = " << header->nlabl << endl;
		for ( i=0; i<header->nlabl; i++ )
			if ( header->labels[i*80] ) cout << &header->labels[i*80] << endl;
	}
		
	p->data_offset(MRCSIZE + header->nsymbt);
	
	long			datatypesize = p->channels()*p->data_type_size();
	long   			datasize = (long)header->nx*header->ny*header->nz*datatypesize;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwMRC: offset=" << p->data_offset() << " typesize=" << datatypesize << " datasize=" << datasize << endl;

	unsigned char*	data = NULL;
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, MRCSIZE);
	
	if ( header->nsymbt ) {
		if ( jflag ) fimg << ss.rdbuf();
		else fimg.write((char *)symop, header->nsymbt);
	}
	
	if ( p->data_pointer() ) {
		if ( p->compound_type() < TComplex ) {
			fimg.write((char *)p->data_pointer(), datasize);
		} else {
			data = new unsigned char[datasize];
			p->pack_transform(data, CentHerm);
			fimg.write((char *)data, datasize);
			delete[] data;
		}
	}
	
	fimg.close();
	
	if ( symop ) delete[] symop;
	delete header;
		
	return 0;
}

