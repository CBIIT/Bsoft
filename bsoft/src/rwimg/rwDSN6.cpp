/**
@file	rwDSN6.cpp
@brief	Functions for reading and writing DSN6 files
@author Bernard Heymann
@date	Created: 20011226
@date 	Modified: 20150130
**/

#include "rwDSN6.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a DSN6 map image file format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 3D image format used in X-ray crystallography.
	Header size:				512 bytes 
	File format extensions:  	.omap, .dsn6, .dn6
	Header byte order always big-endian, composed of short's.
	Byte order determination:	Scaling fields must be less than 256.
	Data type: 					unsigned char or byte.
**/
int 	readDSN6(Bimage* p, int readdata)
{
    ifstream*		fimg = new ifstream;
    fimg->open(p->file_name());
    if ( fimg->fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readDSN6: File opened" << endl;

	DSN6head*	header = new DSN6head;
	
	fimg->read((char *)header, DSN6SIZE);
	if ( fimg->fail() ) return -2;
	
    // Determine byte order and swap bytes if from little-endian machine
    unsigned char*   b = (unsigned char *) header;
	int				i;
//	int				sb = 0;
    if ( header->scale > 255 && header->cell_scale > 255 ) {
    	if ( verbose & VERB_PROCESS )
	    	cerr << "Warning: Swapping header byte order for 2-byte types" << endl;
//		sb = 2;
    	for ( i=0; i<40; i+=2 ) swapbytes(b+i, 2);
    }
    
	if ( header->x_extent < 1 || header->y_extent < 1 || header->z_extent < 1 ) {
		cerr << "Error readDSN6: The file " << p->file_name() << " is not a valid DSN6 image!" << endl;
		fimg->close();
		delete fimg;
		return -3;
	}
	
	// Map the parameters
	p->size(header->x_extent, header->y_extent, header->z_extent);
	p->images(1);
	p->channels(1);
	p->data_type(UCharacter);
	p->data_offset(DSN6SIZE);
	
	UnitCell	uc(header->a*1.0/header->cell_scale, header->b*1.0/header->cell_scale,
		header->c*1.0/header->cell_scale, header->alpha*1.0/header->cell_scale,
		header->beta*1.0/header->cell_scale, header->gamma*1.0/header->cell_scale);
	
	p->unit_cell(uc);
	
	if ( header->x_sampling && header->y_sampling && header->z_sampling )
		p->sampling(uc.a()/header->x_sampling, uc.b()/header->y_sampling, uc.c()/header->z_sampling);
	
	// Allocating the single sub-image and setting its origin
	p->origin(-header->x_start, -header->y_start, -header->z_start);

	float		min = (3 - header->plus)*header->scale*1.0/header->product;
	float		max = (253 - header->plus)*header->scale*1.0/header->product;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readDSN6: min=" << min << " max=" << max << endl;
	
	delete header;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwDSN6: Memory freed" << endl;

    p->page_size(8, 8, 8);	 // Page size = one brick
	
	if ( readdata)
		p->read_data( fimg, -1, 2, 0, 0 );
	
	p->minimum(0);
	p->maximum(255);

	fimg->close();
	delete fimg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwDSN6: Data read" << endl;

	return 0;
}

/**
@brief	Writing a DSN6 map image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D image format used in X-ray crystallography.
**/
int 	writeDSN6(Bimage* p)
{
	p->color_to_simple();
	
	float		min = p->minimum();
	float		max = p->maximum();
	
	p->change_type(UCharacter);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwDSN6: Writing DSN6 file" << endl;

	DSN6head*	header = new DSN6head;
	memset(header, 0, sizeof(DSN6head));
	
	// Map the parameters
	UnitCell	uc = p->unit_cell();
	header->x_extent = p->sizeX();
	header->y_extent = p->sizeY();
	header->z_extent = p->sizeZ();
	header->x_start = (short) (-p->image[0].origin()[0] - 0.5);
	header->y_start = (short) (-p->image[0].origin()[1] - 0.5);
	header->z_start = (short) (-p->image[0].origin()[2] - 0.5);
	header->x_sampling = (short) (uc.a()/p->sampling(0)[0] + 0.5);
	header->y_sampling = (short) (uc.b()/p->sampling(0)[1] + 0.5);
	header->z_sampling = (short) (uc.c()/p->sampling(0)[2] + 0.5);
	
	header->cell_scale = 100;
	header->a = (short) (uc.a()*header->cell_scale + 0.5);
	header->b = (short) (uc.b()*header->cell_scale + 0.5);
	header->c = (short) (uc.c()*header->cell_scale + 0.5);
	header->alpha = (short) (uc.alpha()*header->cell_scale*180/M_PI + 0.5);
	header->beta = (short) (uc.beta()*header->cell_scale*180/M_PI + 0.5);
	header->gamma = (short) (uc.gamma()*header->cell_scale*180/M_PI + 0.5);
	
	header->scale = 100;
	header->product = (short) (header->scale*(253 - 3)/(max - min) + 0.5);
	header->plus = (short) ((3*max - 253*min)/(max - min) + 0.5);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeDSN6: product=" << header->product << " plus=" << header->plus << endl;
	
	p->data_offset(DSN6SIZE);
    p->page_size(8, 8, 8);	 // Page size = one brick
	
    unsigned char*	b = (unsigned char *) header;
	long   			i, j, x, y, z, px, py, pz;
//	int				sb = 0;
	long			datatypesize = p->data_type_size();
	long			valuesize = datatypesize*p->channels();
	long			pagesize = p->page_size()[0]*p->page_size()[1]*p->page_size()[2]*valuesize;
	unsigned char*	data = p->data_pointer();
	unsigned char*	page = new unsigned char[pagesize];
	memset(page, 0, pagesize);
	
	if ( systype(0) == LittleIEEE ) {
//		sb = 1;
    	for ( i=0; i<40; i+=2 ) swapbytes(b+i, 2);
		if ( verbose & VERB_PROCESS )
			cout << "Swapping byte pairs for the DSN6 format" << endl;
	}
	
    ofstream        fimg(p->file_name().c_str());
    if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, DSN6SIZE);
	
	for ( pz=0; pz<p->sizeZ(); pz+=p->page_size()[2] )
		for ( py=0; py<p->sizeY(); py+=p->page_size()[1] )
			for ( px=0; px<p->sizeX(); px+=p->page_size()[0] ) {
				for ( x=0; x<p->page_size()[0] && x<p->sizeX()-px; x++ )
					for ( y=0; y<p->page_size()[1] && y<p->sizeY()-py; y++ )
						for ( z=0; z<p->page_size()[2] && z<p->sizeZ()-pz; z++ ) {
							i = (((z+pz)*p->sizeY() + y+py)*p->sizeX() + x+px)*valuesize;
							j = ((z*p->page_size()[1] + y)*p->page_size()[0] + x)*valuesize;
							memcpy(page+j, data+i, valuesize);
						}
//				if ( sb ) for ( i=0; i<pagesize; i+=2 ) swapbytes(page+i, 2);
				for ( i=0; i<pagesize; i+=2 ) swapbytes(page+i, 2);
				fimg.write((char *)page, pagesize);
			}
			
	fimg.close();
	
	delete[] page;
	
	delete header;
		
	return 0;
}

