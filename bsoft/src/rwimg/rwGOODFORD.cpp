/**
@file	rwGOODFORD.cpp
@brief	Functions for reading and writing Peter Goodford's GRID files
@author Bernard Heymann
@date	Created: 20000924
@date 	Modified: 20120321
**/

#include "rwGOODFORD.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading Peter Goodford's GRID map image file format.
A 3D image format used for electrostatic potential maps.
	Header size:				168 bytes (fixed).
	File format extensions:  	.pot
	The data is packed as slices, each slice with a header of size 20 bytes.
	Byte order determination:	File type and third dimension values
								must be less than 256*256.
	Data types: 				1 = float.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readGOODFORD(Bimage *p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	GOODFORDhead*	header = new GOODFORDhead;
	
	fimg->read((char *)header, GOODFORDSIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    int				i, sb = 0;
    if ( ( abs( header->grdflg ) > SWAPTRIG ) || ( abs(header->km) > SWAPTRIG ) ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
		int 	extent = GOODFORDSIZE;	// Swap whole header
		swapbytes(b, 4);
    	for ( i=76; i<extent; i+=4 ) swapbytes(b+i, 4);
    }
/*    
    // Convert VAX floating point types if necessary
    if ( header->h > 100 || header->h < 0.01 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Converting VAX floating point" << endl;
    	vax = 1;
    	vax2ieee(b+76, sb);
		for ( i=116; i<132; i+=4 )		// All floating point parameters
    	    vax2ieee(b+i, sb);
    }
*/    
	// Map the parameters
	p->size(header->im, header->jm, header->km);
	p->page_size(header->im, header->jm, 1);
	p->images(1);
	p->channels(1);
	switch ( header->grdflg ) {
		case 1: p->data_type(Float); break;
		default: p->data_type(Float); break;
	}
	p->data_offset(GOODFORDSIZE + GFSLICESIZE + sizeof(int));
	p->sampling(header->h, header->h, header->h);
	p->label(header->title);
	
	// Allocating the single sub-image and setting its origin
	p->origin(-header->ox, -header->oy, -header->oz);
	
	long 	pad = GFSLICESIZE + 2*sizeof(int);
	long	pagesize = p->sizeX()*p->sizeY()*p->data_type_size();
	long	offset = p->data_offset();
	unsigned char*	data = NULL;
	
	if ( readdata ) {
		p->data_alloc();
		for ( i=0, data=p->data_pointer(); i<p->sizeZ(); i++, data+=pagesize, offset+=pagesize+pad )
			fread_large(data, pagesize, offset, fimg);
		if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	delete fimg;
		
	delete header;
		
	return 0;
}

/**
@brief	Writing Peter Goodford's GRID map image file format.
A 3D image format used for electrostatic potential maps.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeGOODFORD(Bimage *p)
{
	p->color_to_simple();
	
	p->change_type(Float);
	
	GOODFORDhead*	header = new GOODFORDhead;
	memset(header, 0, sizeof(GOODFORDhead));
	
	GFslice_head*	slice_header = new GFslice_head;
	memset(slice_header ,0, sizeof(GFslice_head));
	
	// Map the parameters
	memset(header->title, ' ', 72);
	strncpy(header->title, p->label().c_str(), 72);
	header->im = slice_header->im = p->sizeX();
	header->jm = slice_header->jm = p->sizeY();
	header->km = p->sizeZ();
	switch ( p->data_type() ) {
		case Float: header->grdflg = 1; break;
		default: header->grdflg = 1; break;
	}
	p->data_offset(GOODFORDSIZE + GFSLICESIZE + sizeof(int));
	header->ox = -p->image->origin()[0];
	header->oy = -p->image->origin()[1];
	header->oz = -p->image->origin()[2];
	header->h = p->sampling(0)[0];
	header->pad1 = header->pad8 = 160;
	header->pad4 = header->pad6 = p->sizeZ();
	header->pad5 = 1;
	slice_header->pad1 = slice_header->pad2 = 12;
	
	long			i;
	int				pagesize = p->sizeX()*p->sizeY()*p->data_type_size();
	
	ofstream		fimg(p->file_name().c_str());
	if ( fimg.fail() )  return -1;

	unsigned char*	data = p->data_pointer();
	
	fimg.write((char *)header, GOODFORDSIZE);

	for ( i=0; i<p->sizeZ(); i++, data += pagesize ) {
		slice_header->k = i + 1;
		fimg.write((char *)slice_header, GFSLICESIZE);
		fimg.write((char *)&pagesize, sizeof(int));
		fimg.write((char *)data, pagesize);
		fimg.write((char *)&pagesize, sizeof(int));
	}
	
	fimg.close();

	delete header;
	delete slice_header;
		
	return 0;
}

