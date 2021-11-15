/**
@file	rwBRIX.cpp
@brief	Functions for reading and writing BRIX files
@author Bernard Heymann
@date	Created: 19990424
@date 	Modified: 20150130
**/

#include "rwBRIX.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a BRIX image file format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int					error code (<0 means failure).
A 3D image format intended for use with the 'O' program.
	It has a text header of 512 bytes (fixed size).
	File format extension:  	.brx.
	Identifier:					':-)' (byte 0).
	Data type: 					byte, packed as 8x8x8 pages.
	The header specifies the following keywords:
		origin: 	origin
		extent: 	size in pixels
		grid:		unit cell size in pixels
		cell:		unit cell parameters in angstroms and degrees
		prod:		scale to convert original map to byte
		plus:		average 
		sigma:		standard deviation of original map
	The voxel size is given by the unit cell size divided by the grid size
	(Since the unit cell size may not fit exactly onto the grid, the voxel
	size as calculated from the unit cell size and grid may not be accurate).
**/
int 	readBRIX(Bimage* p, int readdata)
{
    ifstream*		fimg = new ifstream;
    fimg->open(p->file_name());
    if ( fimg->fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBRIX: File opened" << endl;

	char*	header = new char[BRIXSIZE];
	
	fimg->read((char *)header, BRIXSIZE);
	if ( fimg->fail() ) return -2;
	
	int				mx(0), my(0), mz(0);
	long			x(1), y(1), z(1);
    float			prod(1);
    int				plus(127);
    float			sigma(1);
	Vector3<double>	ori;
    
	// Allocating the single sub-image
	p->images(1);
	
    // extract the image attributes from the file header
    char    	*var;
	UnitCell	uc;
    for ( int i=0; i<BRIXSIZE; i++ ) header[i] = tolower(header[i]);
    if ( ( var = strstr(header,"origin") ) != NULL )  				// Origin
    	    sscanf(var,"origin %lf %lf %lf", &ori[0], &ori[1], &ori[2]);
    if ( ( var = strstr(header,"extent") ) != NULL ) 
    	    sscanf(var,"extent %ld %ld %ld", &x, &y, &z);		// Size
    if ( ( var = strstr(header,"grid") ) != NULL ) 
    	    sscanf(var,"grid %d %d %d", &mx, &my, &mz);	    		// Unit cell size
    if ( ( var = strstr(header,"cell") ) != NULL ) 
    	    sscanf(var,"cell %lf %lf %lf %lf %lf %lf",
    	    	&uc[0], &uc[1], &uc[2], &uc[3], &uc[4], &uc[5]);	// Unit cell parameters
    if ( ( var = strstr(header,"prod") ) != NULL ) 
    	    sscanf(var,"prod %f", &prod);	    	    			// Scale product
    if ( ( var = strstr(header,"plus") ) != NULL ) 
    	    sscanf(var,"plus %d", &plus);	    	    			// Scale shift
    if ( ( var = strstr(header,"sigma") ) != NULL ) 
    	    sscanf(var,"sigma %f", &sigma);	    	    			// RMS

	// Negating the origin to get the correct convention for Bsoft
	p->image->origin(-ori);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBRIX: image header parsed" << endl;

	p->size(x, y, z);
    p->page_size(8, 8, 8);	// Page size = one brick
	if ( p->sizeZ() < 8 ) p->page_size(8, 8, p->sizeZ());
    p->channels(1);
	p->data_type(UCharacter);
	p->data_offset(BRIXSIZE);
	p->minimum(0);
	p->maximum(255);
//	p->average(plus);
//	p->scale = prod;
	p->standard_deviation(sigma*prod);		// Scaled to correct value
	if ( mx && my && mz ) p->sampling(uc.a()/mx, uc.b()/my, uc.c()/mz);
	
	delete[] header;
		
	if ( readdata )
		p->read_data( fimg, -1, 0, 0, 0 );
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readBRIX: image read" << endl << endl;
	
	fimg->close();
	delete fimg;
	
	return 0;
}

/**
@brief	Writing a BRIX image file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 3D image format intended for use with the 'O' program.
**/
int 	writeBRIX(Bimage* p)
{
	p->change_type(UCharacter);
	p->color_to_simple();
	
	char*	header = new char[BRIXSIZE];
	memset(header, ' ', BRIXSIZE);
	
	// Ensure that the scaling is done correctly
//	if ( p->scale < 250/(p->max - p->min) || p->scale > 255.5/(p->max - p->min) ) 
//		img_rescale_to_min_max(p, 0, 255);
	
	// Map the parameters
	double			scale(1);
	Vector3<double>	u(p->sampling(0));
	UnitCell		uc = p->unit_cell();
	int 			mx = (int) (uc.a()/u[0]);
	int 			my = (int) (uc.b()/u[1]);
	int 			mz = (int) (uc.c()/u[2]);
	float			prod(scale);
	int 			plus = (int) p->minimum();
    float			sigma = p->standard_deviation()/scale;	// Sigma scaled to original value
	
	// This is a band-aid to overcome the limitations of the image format
	if ( fabs(uc.a() - u[0]*mx) > 0.001 || fabs(uc.b() - u[1]*my) || fabs(uc.c() - u[2]*mz) ) {
		uc.a(u[0]*mx);
		uc.b(u[1]*my);
		uc.c(u[2]*mz);
		if ( verbose & VERB_FULL )
			cerr << "Warning: Resetting the unit cell to: " << uc.a() << " " << uc.b() << " " << uc.c() << " A" << endl;
	}
	
    snprintf(header, BRIXSIZE, ":-) Origin%5.0f%5.0f%5.0f Extent%5ld%5ld%5ld Grid%5d%5d%5d "
    	"Cell %10.3f%10.3f%10.3f%10.3f%10.3f%10.3f "
    	"Prod%12.5f Plus%8d Sigma %12.5f",
    	-p->image->origin()[0], -p->image->origin()[1], -p->image->origin()[2], p->sizeX(), p->sizeY(), p->sizeZ(), mx, my, mz,
		uc.a(), uc.b(), uc.c(), uc.alpha()*180/M_PI, uc.beta()*180/M_PI, uc.gamma()*180/M_PI, prod, plus, sigma);

	p->data_offset(BRIXSIZE);
	p->page_size(8, 8, 8);
	if ( p->sizeZ() < 8 ) p->page_size(8, 8, p->sizeZ());

	long			i, j, x, y, z, px, py, pz;
	long   			datatypesize = p->data_type_size();
	long   			valuesize = datatypesize*p->channels();
	long   			pagesize = p->page_size()[0]*p->page_size()[1]*p->page_size()[2]*valuesize;
	unsigned char*	data = p->data_pointer();
	unsigned char*	page = new unsigned char[pagesize];
	memset(page, 0, pagesize);
	
    ofstream        fimg(p->file_name());
    if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, p->data_offset());
	
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
				fimg.write((char *)page, pagesize);
			}
			
	fimg.close();
	
	delete[] page;
	
	delete[] header;
		
	return 0;
}

