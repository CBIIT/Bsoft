/**
@file	rwCCP4.cpp
@brief	Functions for reading and writing CCP4 files
@author Bernard Heymann
@date	Created: 19990410
@date 	Modified: 20120730
**/

#include "rwCCP4.h"
#include "UnitCell.h"
#include "file_util.h"
#include "utilities.h"

#ifndef NOSYMOP
#include "rwsymop.h"
#endif

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Setting a CCP4 style machine stamp.
@param	*machine_stamp	machine stamp string.
@return	int					error code (<0 means failure).
The 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
**/
int			set_CCP4_machine_stamp(char* machine_stamp)
{
	switch ( systype(0) ) {
		case BigIEEE:
			machine_stamp[0] = machine_stamp[1] = 17;
			break;
		case LittleIEEE:
			machine_stamp[0] = 68;
			machine_stamp[1] = 65;
			break;
		case LittleVAX:
			machine_stamp[0] = 34;
			machine_stamp[1] = 65;
			break;
		default:
			break;
	}
	
	return 0;
}

/**
@brief	Reading a CCP4 map image file format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 3D image format used in X-ray crystallography.
	Header size:				1024 bytes followed by the symmetry 
		operator table which is composed of 80 character lines, each 
		line for a symmetry operator.
	File format extensions:  	.map, .ccp, .ccp4
	The identifier is a 4-byte machine stamp:
				1	Big-endian IEEE 	17 17 00 00
                2	VAX 				34 65 00 00
				3	Cray				-
                4	Little-endian IEEE	68 65 00 00
                5	Convex				85 17 00 00
				6	Fijitsu VP			-
				(Note: not always implemented - so unreliable)
	Byte order determination:	Data type and third dimension values
								must be less than 256*256.
	Data types: 				0 = signed char, 1 = short, 2 = float,
								3 = complex short, 4 = complex float.
	Transform type: 			Centered hermitian
								The x-dimension contains the x-size
								of the full transform
**/
int 	readCCP4(Bimage* p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readCCP4: File opened" << endl;

	CCP4head*	header = new CCP4head;
	
	fimg->read((char *)header, CCP4SIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    long  			i, sb(0);
    if ( abs(header->mode) > SWAPTRIG || abs(header->mapc) > SWAPTRIG ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		long	 	extent = CCP4SIZE - 800; // exclude labels from swapping
    	for ( i=0; i<extent; i+=4 ) 
			if ( i < 208 || i > 212 ) swapbytes(b+i, 4);
    }
    
	if ( header->nx < 1 || header->ny < 1 || header->nz < 1 || header->mode < 0 ) {
		cerr << "Error readCCP4: The file " << p->file_name() << " is not a valid CCP4 image!" << endl;
		cerr << "Size:  " << header->nx << tab << header->ny <<
			tab << header->nz << endl;
		fimg->close();
		return -3;
	}
	
    // Convert VAX floating point types if necessary
    // Machine stamp for vax is: 34 65 00 00
    if ( ( header->machst[0] == 34 ) || ( header->amin > header->amax ) ) {
		if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Converting VAX floating point" << endl;
		for ( i=40; i<64; i+=4 )		// Unit cell parameters
    	    vax2ieee(b+i, sb);
		for ( i=76; i<88; i+=4 )		// Min, max, mean
    	    vax2ieee(b+i, sb);
		for ( i=100; i<204; i+=4 )		// Skew matrix and extra
    	    vax2ieee(b+i, sb);
		vax2ieee(b+216, sb); 			// RMS
    }
	
	// Map the parameters
	if ( strncmp(header->map, "STK ", 4) == 0 ) {
		p->images(header->nz);
		p->size(header->nx, header->ny, 1);
	} else {
		p->images(1);
		p->size(header->nx, header->ny, header->nz);
	}
	p->channels(1);
	switch ( header->mode ) {
		case 0: p->data_type(SCharacter); break;
		case 1: p->data_type(Short); break;
		case 2: p->data_type(Float); break;
		case 3: p->data_type(Short); p->compound_type(TComplex); p->channels(2); break;
		case 4: p->data_type(Float); p->compound_type(TComplex); p->channels(2); break;
		default: p->data_type(SCharacter); break;
	}
	p->data_offset(CCP4SIZE + header->nsymbt);
	if ( header->mode > 2 ) {
		fimg->seekg(0, ios::end);
		if ( (double) fimg->tellg() > p->data_offset() + 0.8*p->alloc_size() )
			p->sizeX(2*(p->sizeX() - 1));
		if ( header->mx%2 == 1 ) p->sizeX(p->sizeX() + 1);	// Quick fix for odd x-size maps
		fimg->clear();
		fimg->seekg(0, ios::beg);
	}
	p->minimum(header->amin);
	p->maximum(header->amax);
	p->average(header->amean);
	p->standard_deviation(header->arms);

	if ( header->mx && header->my && header->mz )
		p->sampling(header->a/header->mx, header->b/header->my, header->c/header->mz);
	p->space_group(header->ispg);
	UnitCell	uc(header->a, header->b, header->c, header->alpha, header->beta, header->gamma);
	p->unit_cell(uc);
	p->label(header->labels);
	
	// Allocating the single sub-image and setting its origin
	p->origin(-header->nxStart, -header->nyStart, -header->nzStart);
	
	delete header;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: Header setup done" << endl;

	long	readsize = p->alloc_size();
	unsigned char*	data = NULL;
	
	if ( readdata ) {
		p->data_alloc();
		if ( p->compound_type() == TSimple ) {
			fread_large(p->data_pointer(), readsize, p->data_offset(), fimg);
			if ( sb ) swapbytes(readsize, p->data_pointer(), p->data_type_size());
		} else {
			readsize = p->channels()*(p->sizeX()/2+1)*p->sizeY()*p->sizeZ()*p->data_type_size();
			data = new unsigned char[readsize];
			fread_large(data, readsize, p->data_offset(), fimg);
			if ( sb ) swapbytes(readsize, data, p->data_type_size());
//			img_unpack_transform(p, data, CentHerm);
			p->unpack_transform(data, CentHerm);
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: Data read" << endl;

	return 0;
}

/**
@brief	Writing a CCP4 map image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D image format used in X-ray crystallography.
**/
int 	writeCCP4(Bimage* p)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: Writing CCP4 file" << endl;
	
	p->color_to_simple();
	
//	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			if ( p->compound_type() == TComplex ) p->change_type(Short);
			else p->change_type(SCharacter);
			break;
    	case UShort: case Short:
			p->change_type(Short); break;
    	default:
			p->change_type(Float); break;
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: data type = " << p->data_type() << endl;

	CCP4head*	header = new CCP4head;
	memset(header, 0, sizeof(CCP4head));

	UnitCell		uc = p->unit_cell();
	
	// Map the parameters
	set_CCP4_machine_stamp(header->machst);
	header->nx = p->sizeX();
	header->ny = p->sizeY();
	if ( p->sizeZ() > 1 || p->images() == 1 ) {
		memcpy(header->map, "MAP ", 4);
		header->nz = p->sizeZ();
	} else {
		memcpy(header->map, "STK ", 4);
		header->nz = p->images();
	}
	if ( p->compound_type() == TComplex ) header->nx = p->sizeX()/2 + 1;	// If a transform, physical storage is nx/2 + 1
	switch ( p->data_type() ) {
		case SCharacter: header->mode = 0; break;
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
	header->mapc = 1;
	header->mapr = 2;
	header->maps = 3;
	header->amin = p->minimum();
	header->amax = p->maximum();
	header->amean = p->average();
	header->arms = p->standard_deviation();
//	cout << "Unitcell %g %g %g\n", uc.a(), uc.b(), uc.c());
	header->a = uc.a();
	header->b = uc.b();
	header->c = uc.c();
	
	// This is a band-aid to overcome the limitations of the image format
	if ( fabs(header->a - p->sampling(0)[0]*header->mx) > 0.001 || fabs(header->b - p->sampling(0)[1]*header->my) > 0.001 ||
			fabs(header->c - p->sampling(0)[2]*header->mz) > 0.001 ) {
		header->a = p->sampling(0)[0]*header->mx;
		header->b = p->sampling(0)[1]*header->my;
		header->c = p->sampling(0)[2]*header->mz;
		if ( verbose & VERB_FULL )
			cerr << "Warning: Resetting the unit cell to: " << header->a << "," << header->b << "," << header->c << endl;
	}
	
	long 			i,j,k;
	Matrix3			frac_mat = uc.skew_matrix();
	header->alpha = uc.alpha()*180/M_PI;
	header->beta = uc.beta()*180/M_PI;
	header->gamma = uc.gamma()*180/M_PI;
	for ( i=k=0; i<3; i++ ) for ( j=0; j<3; j++, k++ ) header->skwmat[k] = frac_mat[i][j];
//	for ( i=0; i<3; i++ ) header->skwtrn[i] = frac_mat[i+9];
	header->ispg = p->space_group();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: Writing labels" << endl;

	header->nlabl = 10;
	strncpy(header->labels, p->label().c_str(), 799);
	header->labels[799] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG rwCCP4: nlabl = " << header->nlabl << endl;
		for ( i=0; i<header->nlabl; i++ )
			if ( header->labels[i*80] ) cout << &header->labels[i*80] << endl;
	}
	
	int				nsym(0);
	Bstring			temp;
	char*			symop = NULL;
	
#ifndef NOSYMOP
	if ( p->space_group() > 0 ) symop = read_symop(temp, p->space_group(), nsym);
#endif

	header->nsymbt = nsym*80;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: nsymbt = " << header->nsymbt << endl;
	
	p->data_offset(CCP4SIZE + header->nsymbt);
	
	long			datatypesize = p->channels()*p->data_type_size();
	long			datasize = (long)header->nx*header->ny*header->nz*datatypesize;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: offset=" << p->data_offset() << " typesize=" << datatypesize << " datasize=" << datasize << endl;
	
	unsigned char*	data = NULL;
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, CCP4SIZE);
	
	if ( header->nsymbt ) fimg.write((char *)symop, header->nsymbt);
	
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
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwCCP4: Done writing" << endl;
	
	return 0;
}

