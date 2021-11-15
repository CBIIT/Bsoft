/**
@file	rwMFF.cpp
@brief	Functions for reading and writing What If MFF files
@author Bernard Heymann
@date	Created: 19990424
@date 	Modified: 20150130
**/

#include "rwMFF.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a MFF image format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 3D file format developed with the What If molecular dynamics package.
	Header size:				268 bytes (fixed). 
	File format extensions:  	.mff.
	Byte order determination:	The third dimension
								must be less than 256*256.
	Data type: 					float.
	Special features: The MFF files are written with double integers separating
		different records by What If.  
		This was solved by adding "padding" fields in the
		header and adding 8 bytes to every page (= section) offset.
		The data is stored as 2D pages separated by 8 bytes.
**/
int 	readMFF(Bimage* p, int readdata)
{
    ifstream*		fimg = new ifstream;
    fimg->open(p->file_name());
    if ( fimg->fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readMFF: File opened" << endl;

	MFFhead*	header = new MFFhead;
	
	fimg->read((char *)header, MFFSIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*	b = (unsigned char *) header;
    int				i, sb = 0, vax = 0;
    if ( abs(header->nz) > SWAPTRIG ) {
    	if ( verbose & VERB_PROCESS )
			cerr << "Warning: Swapping header byte order for 4-byte types" << endl;
    	sb = 1;
    	for ( i=84; i<MFFSIZE; i+=4 ) swapbytes(b+i, 4);
    }
    
    // Convert VAX floating point types if necessary
    if ( header->a > SWAPTRIG || header->b > SWAPTRIG || header->c > SWAPTRIG ) {
    	if ( verbose & VERB_PROCESS )
			cerr << "Warning: Converting VAX floating point" << endl;
    	vax = 1;
		for ( i=92; i<136; i+=4 )		// Unit cell parameters & step sizes
    	    vax2ieee(b+i, sb);
		for ( i=204; i<MFFSIZE; i+=4 )	// Origin & skew matrix
    	    vax2ieee(b+i, sb);
    }
    
	// Map the parameters
	p->size(header->nx, header->ny, header->nz);
	p->images(1);
	p->channels(1);
	p->data_type(Float);
	p->data_offset(MFFSIZE);
	
	UnitCell	uc(header->a, header->b, header->c, header->alpha, header->beta, header->gamma);
	p->unit_cell(uc);
	
	p->sampling(header->stepu, header->stepv, header->stepw);
	p->label(header->title);
	
	p->page_size(p->size()[0], p->size()[1], 1);
	
	// Allocating the single sub-image and setting its origin
	p->origin(-header->uvworg[0]/header->uvwmat[0][0],
		-header->uvworg[1]/header->uvwmat[1][1],
		-header->uvworg[2]/header->uvwmat[2][2]);

	delete header;
	
	if ( readdata )
		p->read_data( fimg, -1, sb, vax, 8 );

	fimg->close();
	delete fimg;

	return 0;
}

/**
@brief	Writeing a MFF image format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D file format developed with the What If molecular dynamics package.
**/
int 	writeMFF(Bimage* p)
{
	p->color_to_simple();
	p->change_type(Float);
	
	MFFhead*	header = new MFFhead;
	memset(header, 0, sizeof(MFFhead));
	
	// Map the parameters
	UnitCell	uc = p->unit_cell();
	header->nx = p->sizeX();
	header->ny = p->sizeY();
	header->nz = p->sizeZ();
    header->idiru = 1; 				// Column, row, section order        
    header->idirv = 2;     
    header->idirw = 3;
	header->uvwmat[0][0] = 1/uc.a();
	header->uvwmat[1][1] = 1/uc.b();
	header->uvwmat[2][2] = 1/uc.c();
	header->uvworg[0] = p->image[0].origin()[0]*header->uvwmat[0][0];
	header->uvworg[1] = p->image[0].origin()[1]*header->uvwmat[1][1];
	header->uvworg[2] = p->image[0].origin()[2]*header->uvwmat[2][2];
	header->stepu = p->sampling(0)[0];
	header->stepv = p->sampling(0)[1];
	header->stepw = p->sampling(0)[2];
	header->a = uc.a();
	header->b = uc.b();
	header->c = uc.c();
	header->alpha = uc.alpha()*180/M_PI;
	header->beta = uc.beta()*180/M_PI;
	header->gamma = uc.gamma()*180/M_PI;
	header->ndiva = (int) (uc.a()/p->sampling(0)[0]);     // Unit cell size in pixels
	header->ndivb = (int) (uc.b()/p->sampling(0)[1]);
	header->ndivc = (int) (uc.c()/p->sampling(0)[2]);
	strncpy(header->title, p->label().c_str(), 80);
	
	header->pad1    = 80;			// Padding fields
    header->pad2[0] = 80;
    header->pad2[1] = 24;
    header->pad3[0] = 24;
    header->pad3[1] = 12;
    header->pad4[0] = 12;
    header->pad4[1] = 12;
    header->pad5[1] = 12;
    header->pad5[0] = 12;
    header->pad5[1] = 12;
    header->pad6[0] = 12;
    header->pad6[1] = 12;
    header->pad7[0] = 12;
    header->pad7[1] = 12;
    header->pad8[0] = 12;
    header->pad8[1] = 36;
    header->pad9[0] = 36;
	
	p->data_offset(MFFSIZE);
	
	long 	datatypesize = p->data_type_size();
	long	pagesize = header->nx*header->ny*datatypesize;
    header->pad9[1] = pagesize;
	
	unsigned char*	data = p->data_pointer();
	
    ofstream        fimg(p->file_name());
    if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, p->data_offset());
	
	for ( int i=0; i<header->nz; i++, data += pagesize ) {
		fimg.write((char *)data, pagesize);
		if ( i < header->nz - 1 )
			fimg.write((char *)&pagesize, sizeof(unsigned long));
		else
			fimg.write((char *)&pagesize, sizeof(unsigned int));
	}
	
	fimg.close();
	
	delete header;
	
	return 0;
}

