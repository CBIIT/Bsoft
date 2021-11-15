/**
@file	rwBIORAD.cpp
@brief	Functions for reading and writing BioRad confocal files
@author Bernard Heymann
@date	Created: 19990427@date Modified: 20120321
**/

#include "rwBIORAD.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			swap_header(BIORADhead*	header)
{
    unsigned char*	b = (unsigned char *) header;
    int				i;
	int				extent = BIORADSIZE - 10; // exclude magnification from swapping
	
	for ( i=0; i<extent; i+=2 ) swapbytes(b+i, 2);
	
	swapbytes(b+extent, 4);
	
	return 0;
}

/**
@brief	Reading a BioRad image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 3D image format intended for confocal light microscopy.
	Header size:				76 bytes (fixed).
	File format extension:  	.PIC.
	Identifier:					12345 (byte 54).
	Byte order determination:	Magic number.
	Data types: 				0 = short, 1 = byte.
**/
int			readBIORAD(Bimage* p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	
	
	BIORADhead*	header = new BIORADhead;

	fimg->read((char *)header, BIORADSIZE);
	if ( fimg->fail() ) return -2;
		
    // Determine byte order and swap bytes if from little-endian machine
    int     	sb = 0;
    if ( header->magic != 12345 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order" << endl;
    	sb = 1;
		swap_header(header);
    }
    if ( header->magic != 12345 ) {
		cerr << "Error: The file " << p->file_name() << " is not a BioRad file!" << endl;
		fimg->close();
		delete fimg;
		return -3;
	}
    
	// Map the parameters
	p->size(header->nx, header->ny, header->nz);
	p->images(1);
	p->channels(1);
	switch ( header->mode ) {
		case 1: p->data_type(UCharacter); break;
		case 0: p->data_type(UShort); break;
		default: p->data_type(UCharacter); break;
	}
	p->minimum(header->black);
	p->maximum(header->white);
	p->data_offset(BIORADSIZE);

	
	delete header;
		
	if ( readdata ) {
		p->data_alloc();
		fread_large(p->data_pointer(), p->alloc_size(), p->data_offset(), fimg);
		if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	delete fimg;

	return 0;
}

/**
@brief	Writing a BioRad image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D image format intended for confocal light microscopy.
**/
int			writeBIORAD(Bimage* p)
{
	if ( p->data_type() != UCharacter && p->data_type() != UShort ) p->change_type(UCharacter);
	p->color_to_simple();
	
	BIORADhead*	header = new BIORADhead;
	memset(header, 0, sizeof(BIORADhead));
	
	int			sb = 0;
	if ( systype(0) < LittleIEEE ) sb = 1;
	
	// Map the parameters
	header->nx = p->size()[0];
	header->ny = p->size()[1];
	header->nz = p->size()[2];
	switch ( p->data_type() ) {
		case UCharacter: header->mode = 1; break;
		case UShort: header->mode = 0; break;
		default: header->mode = 0; break;
	}
	header->black = (short) p->minimum();
	header->white = (short) p->maximum();
	header->magic = 12345;		// Magic number for the BioRad format
	p->data_offset(BIORADSIZE);
	
	if ( sb ) {
		swap_header(header);
		swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write((char *)header, p->data_offset());
	if ( p->data_pointer() ) fimg.write((char *)p->data_pointer(), p->alloc_size());
	
	fimg.close();
	
	delete header;
		
	return 0;
}

