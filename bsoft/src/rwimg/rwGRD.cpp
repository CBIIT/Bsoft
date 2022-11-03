/**
@file	rwGRD.cpp
@brief	Functions for reading and writing GRD files
@author Bernard Heymann
@date	Created: 19990410
@date 	Modified: 20210304
**/

#include "rwGRD.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			grid_decompress(Bimage* p, long len, unsigned char* buf)
{
	if ( p->data_type() != UCharacter ) return -1;
	if ( !p->data_pointer() ) return -1;
	
	long			i, j, k, cnt;
	unsigned char	val;
	unsigned char*	data = p->data_pointer();
	
	for ( i=j=0; i<p->data_size() && j<len;  ) {
		cnt = buf[j++];
		val = buf[j++];
//		cout << cnt << tab << int(val) << endl;
		for ( k=0; k<cnt; ++k ) data[i++] = val;
	}
			
	return 0;
}

/*
	Run length encoding
	The data is encoded in byte pairs, where every pair
	consists of a count and a value.
	If the count reaches 255, it writes a pair and
	restarts the count.
*/
unsigned char*	grid_compress(Bimage* p, long& compress)
{
	compress = 0;
	if ( p->data_type() != UCharacter ) return NULL;
	if ( !p->data_pointer() ) return NULL;
	
	long			i, j, cnt(0);
	unsigned char	val;
	unsigned char*	data = p->data_pointer();
	unsigned char*	buf = new unsigned char[p->alloc_size()];
	
	for ( i=j=0, val = data[0]; i<p->data_size(); ++i ) {
		if ( val != data[i] || cnt > 254 ) {
//			cout << cnt << tab << int(val) << endl;
			buf[j++] = cnt;
			buf[j++] = val;
			val = data[i];
			cnt = 1;
		} else cnt++;
	}
//	cout << cnt << tab << int(val) << endl;
	buf[j++] = cnt;
	buf[j++] = val;
	
	compress = j;
			
	return buf;
}


/**
@brief	Reading a Basel GRD map image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@param 	img_select	image selection in multi-image file (-1 = all images).
@return	int			error code (<0 means failure).
A 3D image format used in electron microscopy.
	Header size:			512 bytes (fixed).
	File format extensions:  	.grd
	The identifier is a 4-byte magic number (not used).
	Machine identifier: 		1 = little endian, 2 = big endian.
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Old data types: 			1 = byte, 2 = float, 3 = complex float
						4 = 3-value vector (float), 5 = view (float)
						6 = RGB byte, 7 = RGB float
						8 = RGBA byte, 9 = RGBA float
						10 = multi
	The Bsoft data type and compound type are combined to give the mode:
		mode = 100*compound_type + data_type
**/
int 	readGRD(Bimage *p, int readdata, int img_select)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	GRDhead*	header = new GRDhead;
	
	fimg->read((char *)header, GRDSIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
	unsigned char*	b = (unsigned char *) header;
	int				i, sb = 0;
    if ( abs( header->mode ) > SWAPTRIG ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		int 	extent = GRDSIZE;	// Swap whole header
    	for ( i=0; i<extent; i+=4 ) swapbytes(b+i, 4);
    }

	// Map the parameters
	p->images(header->nn);
	if ( img_select > -1 ) {
		p->images(1);
		if ( img_select >= header->nn ) img_select = header->nn - 1;
	}

	p->size(header->nx, header->ny, header->nz);
	p->channels(header->channels);
	p->sampling(header->ux, header->uy, header->uz);
	
	if ( header->version < 99 ) {	// Old spec
		switch ( header->mode ) {
			case 1: p->data_type(UCharacter); break;
			case 2: p->data_type(Float); break;
			case 3: p->data_type(Float);
				p->compound_type(TComplex); p->channels(2);
				p->fourier_type(Standard); break;
			case 4: p->data_type(Float);
				p->compound_type(TVector3); p->channels(3); break;
			case 5: p->data_type(Float);
				p->compound_type(TView); p->channels(4); break;
			case 6: p->data_type(UCharacter);
				p->compound_type(TRGB); p->channels(3); break;
			case 7: p->data_type(Float);
				p->compound_type(TRGB); p->channels(3); break;
			case 8: p->data_type(UCharacter);
				p->compound_type(TRGBA); p->channels(4); break;
			case 9: p->data_type(Float);
				p->compound_type(TRGBA); p->channels(4); break;
			case 10: p->data_type(Float);
				p->compound_type(TMulti); break;
			default: p->data_type(UCharacter); break;
		}
	} else {
		p->compound_type(CompoundType(header->mode/100));
		if ( p->compound_type() == TComplex )
			p->fourier_type(Standard);
		p->data_type(DataType(header->mode%100));
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readGRD: mode=" << header->mode << " datatype=" << p->data_type() << endl;
	p->data_offset(GRDSIZE + header->doffset);
	p->minimum(header->amin);
	p->maximum(header->amax);
	p->average(header->amean);
	p->standard_deviation(header->arms);
	p->unit_cell(UnitCell(header->a, header->b, header->c, header->alpha, header->beta, header->gamma));
	
	// Allocating the single sub-image and setting its parameters
	p->background(0, header->bkg);
	p->origin(header->ox, header->oy, header->oz);
	p->view(header->vx, header->vy, header->vz, header->va);

	if ( header->doffset > 0 ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readGRD: reading JSON header:" << header->doffset << endl;
		char*			ext_header = new char[header->doffset];
		fimg->read((char *)ext_header, header->doffset);
		if ( fimg->fail() ) return -1;
		if ( verbose & VERB_DEBUG )
			cout << ext_header << endl;
		JSparser	parser(ext_header, header->doffset);
		p->meta_data() = parser.parse();
		if ( img_select > -1 ) p->meta_data_retain_one_image(img_select);
		p->update_from_meta_data();
		delete[] ext_header;
	}
	
	tm*			t = p->get_localtime();
	sscanf(header->datetime, "%4d%2d%2d%2d%2d%2d",
			&t->tm_year, &t->tm_mon, &t->tm_mday, 
			&t->tm_hour, &t->tm_min, &t->tm_sec);
    t->tm_mon -= 1;
    t->tm_year -= 1900;
	p->set_time(t);
	
	p->symmetry(header->symmetry);
	p->label(header->label);
	
	long			readsize = p->alloc_size();
	if ( img_select < 0 ) img_select = 0;
		
	if ( readdata ) {
		p->data_alloc();
		if ( header->rle ) {
			unsigned char*	buf = new unsigned char[header->rle];
			fimg->read((char *)buf, header->rle);
			grid_decompress(p, header->rle, buf);
			delete[] buf;
		} else {
			if ( fread_large(p->data_pointer(), readsize, p->data_offset() + img_select*readsize, fimg) < 0 ) return -1;
			if ( sb ) swapbytes(readsize, p->data_pointer(), p->data_type_size());
		}
	}
		
	fimg->close();
	delete fimg;

	delete header;

	return 0;
}

/**
@brief	Writing a Basel GRD map image file format.
@param	*p			the image structure.
@param	flags		flags to set output properties.
@return	int			error code (<0 means failure).
A 3D image format used in electron microscopy.
	flags:
		0	no compression
		1	rul-length compression (only byte data types)
**/
int 	writeGRD(Bimage *p, int flags)
{
	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			p->change_type(UCharacter); break;
    	default:
			p->change_type(Float); break;
    }

	long			compress(flags);
	unsigned char*	rle = NULL;
	
	if ( compress ) rle = grid_compress(p, compress);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: header size=" << sizeof(GRDhead) << endl;

	GRDhead*	header = new GRDhead;
	memset(header, 0, sizeof(GRDhead));
	
	header->magic = 42;
	header->version = 100;
	
	// Map the parameters
	UnitCell	uc = p->unit_cell();
	header->nx = p->size()[0];
	header->ny = p->size()[1];
	header->nz = p->size()[2];
	header->nn = p->images();
	header->channels = p->channels();
	header->ux = p->sampling(0)[0];
	header->uy = p->sampling(0)[1];
	header->uz = p->sampling(0)[2];
	
/*	switch ( systype(0) ) {
		case BigIEEE: header->origp = 2; break; 	// MIPS, Motorola
		case LittleIEEE: header->origp = 1; break;	// Alpha, Intel 
		case LittleVAX: header->origp = 1; break;	// VAX 
		default: break;
	}
	switch ( p->data_type() ) {
		case UCharacter:
			switch ( p->compound_type() ) {
				case TRGB: header->mode = 6; break;
				case TRGBA: header->mode = 8; break;
				default: header->mode = 1;
			}
			break;
		case Float:
			switch ( p->compound_type() ) {
				case TComplex: header->mode = 3; break;
				case TVector3: header->mode = 4; break;
				case TView: header->mode = 5; break;
				case TRGB: header->mode = 7; break;
				case TRGBA: header->mode = 9; break;
				case TMulti: header->mode = 10; break;
				default: header->mode = 2;
			}
			break;
		default: header->mode = 1; break;
	}
*/
	
	header->mode = 100*p->compound_type() + p->data_type();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: mode=" << header->mode << endl;

	stringstream    ss;
	ss << p->meta_data() << endl;
//	cout << p->meta_data() << endl;
//	ss.seekp(0, ios::end);
	header->doffset = ss.tellp();
//	string			jstr = ss.str();
//	cout << jstr << endl;
//	header->doffset = jstr.size();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: offset=" << header->doffset << endl;
	
	p->data_offset(GRDSIZE+header->doffset);
	header->amin = p->minimum();
	header->amax = p->maximum();
	header->amean = p->average();
	header->arms = p->standard_deviation();
	
	header->a = uc.a();
	header->b = uc.b();
	header->c = uc.c();
	header->alpha = uc.alpha()*180/M_PI;
	header->beta = uc.beta()*180/M_PI;
	header->gamma = uc.gamma()*180/M_PI;

	header->bkg = p->background(long(0));
	header->ox = p->image->origin()[0];
	header->oy = p->image->origin()[1];
	header->oz = p->image->origin()[2];
	header->vx = p->image->view()[0];
	header->vy = p->image->view()[1];
	header->vz = p->image->view()[2];
	header->va = p->image->view().angle();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: bkg=" << header->bkg << endl;
	
	tm*			t = p->get_localtime();
	snprintf(header->datetime, 16, "%4d%02d%02d%02d%02d%02d", t->tm_year+1900,
			t->tm_mon+1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: datetime=" << header->datetime << endl;
	
	strncpy(header->symmetry, p->symmetry().c_str(), 30);
	strncpy(header->label, p->label().c_str(), 200);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: label=" << header->label << endl;

	header->rle = compress;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGRD: rle=" << header->rle << endl;

	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, GRDSIZE);
	fimg << ss.rdbuf();
//	fimg << jstr << endl;

	if ( p->data_pointer() ) {
		if ( compress && rle )
			fimg.write((char *)rle, compress);
		else
			fimg.write((char *)p->data_pointer(), p->alloc_size());
	}
	
	fimg.close();
	
	delete header;
	if ( rle ) delete[] rle;
		
	return 0;
}

