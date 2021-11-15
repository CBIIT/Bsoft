/**
@file	rwEM.cpp
@brief	Functions for reading and writing EM files
@author Bernard Heymann
@date	Created: 19990418
@date 	Modified: 20120321
**/

#include "rwEM.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a EM map image file format.
A 3D image format used in electron microscopy.
	Header size:				512 bytes (fixed).
	File format extensions:  	.em, .EM
	The identifier is a machine stamp in the first byte:
				0	OS-9
                1	VAX (VMS little endian)
                2	CONVEX (historical)
                3	SGI IRIX, Linux (big endian)
                4	SUN
                5	MAC
                6	DEC UNIX, Linux (little endian)
				(Note: not always implemented - so unreliable)
	Byte order determination:	The z-dimension value
								must be less than 256*256.
	Data types: 				1 = byte, 2 = short, 4 = int, 5 = float,
								8 = complex float.
	Transform type: 			Hermitian
								The x-dimension contains the x-size
								of the full transform
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readEM(Bimage* p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	EMhead*		header = new EMhead;
	
	fimg->read((char *)header, EMSIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if necessary
    unsigned char*	b = (unsigned char *) header;
    int				i, sb = 0;
    if ( abs(header->nz) > SWAPTRIG || header->nz < 1 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		int 	extent = EMSIZE;	// exclude labels from swapping
    	for ( i=4; i<16; i+=4 ) swapbytes(b+i, 4);		// Image size fields
		for ( i=96; i<extent; i+=4 ) swapbytes(b+i, 4); // Extras
    }
    
	// Map the parameters
	p->size(header->nx, header->ny, header->nz);
	p->images(1);
	p->channels(1);
	
    switch( header->stamp[3] ){
        case 1 : p->data_type(UCharacter); break;	    // Images
        case 2 : p->data_type(Short); break;
        case 4 : p->data_type(Integer); break;
        case 5 : p->data_type(Float); break;
        case 8 : p->data_type(Float);  		// Transform
			p->compound_type(TComplex); p->channels(2);
			break;
        default : p->data_type(UCharacter);
    }
    
	p->data_offset(EMSIZE);
	p->label(header->label);
	
	
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
			p->unpack_transform(data, Hermitian);
			delete[] data;
		}
	}
	
	fimg->close();
	delete fimg;
	
	delete header;
		
	return 0;
}

/**
@brief	Writing a EM map image file format.
A 3D image format used in electron microscopy.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeEM(Bimage* p)
{
	p->color_to_simple();
	
	if ( p->compound_type() == TComplex ) p->change_type(Float);
	
    switch ( p->data_type() ) {
		case Bit: case UCharacter: case SCharacter:
			p->change_type(UCharacter); break;
    	case UShort: case Short:
			p->change_type(Short); break;
    	case UInteger: case Integer:
			p->change_type(Integer); break;
    	default:
			p->change_type(Float); break;
    }
	
    switch ( p->data_type() ) {
    	case SCharacter: p->change_type(UCharacter); break;
    	case UShort: p->change_type(Short); break;
    	default: break;
    }
	
 	EMhead*		header = new EMhead;
	memset(header, 0, sizeof(EMhead));
	
	// Map the parameters
	header->nx = p->sizeX();
	header->ny = p->sizeY();
	header->nz = p->sizeZ();
	switch ( systype(0) ) {
		case BigIEEE: header->stamp[0] = 3; break;		// SGI machine stamp
		case LittleIEEE: header->stamp[0] = 6; break;	// Alpha, Intel
		case LittleVAX: header->stamp[0] = 1; break;	// VAX
		default: break;
	}
	switch ( p->data_type() ) {
        case UCharacter : header->stamp[3] = 1; break;		// Images
        case Short : header->stamp[3] = 2; break;
        case Integer   : header->stamp[3] = 4; break;
        case Float : header->stamp[3] = 5;
			if ( p->compound_type() == TComplex ) header->stamp[3] = 8;
			break;	// Transform
        default : header->stamp[3] = 1;
    }
    
	strncpy(header->label, p->label().c_str(), 80);
	
	p->data_offset(EMSIZE);
	
	// If a transform, then the physical storage in x is only half+1
	int				xstore = header->nx;
	if ( p->compound_type() == TComplex ) xstore = header->nx/2 + 1;
	
	int				datatypesize = p->channels()*p->data_type_size();
	long			datasize = (long) xstore*header->ny*header->nz*datatypesize;
	
	unsigned char*	data = NULL;
	
	ofstream		fimg(p->file_name().c_str());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, p->data_offset());

	if ( p->data_pointer() ) {
		if ( p->compound_type() < TComplex ) {
			fimg.write((char *)p->data_pointer(), datasize);
		} else {
			data = new unsigned char[datasize];
			p->pack_transform(data, Hermitian);
			fimg.write((char *)data, datasize);
			delete[] data;
		}
	}
	
	fimg.close();

	delete header;
		
	return 0;
}

