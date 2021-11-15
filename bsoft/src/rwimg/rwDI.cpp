/**
@file	rwDI.cpp
@brief	Functions for reading (only) Digital Instruments files
@author Bernard Heymann
@date	Created: 19990424
@date 	Modified: 20120321
**/

#include "rwDI.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a Digital Instruments image file format.
A two-image 2D image format used in Digital Instruments atomic 
		force microscopes.
	The header is text varying with the version of DI software, which may
		be 4,8,12,20K and is supposed to be read from the header itself.
	The header contains the following keywords (tags):
		image list: 			presumably one entry per image
		Data offset:			offset to start of data, may however be
								incorrect if the file was transferred
								with ftp as text
		Data length:			total length of data block
		Samps/line: 			samples per line (x-dimension)
		Number of lines:		lines (y-dimension)
		list end:				end of header
	File format extensions:  	.di, .DI
	Byte order:					Always little endian because the DI 
								software runs on Intel machines and
								there is no way to detect byte order from
								the text header.
	Data type: 					short.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int					error code (<0 means failure).
**/
int 	readDI(Bimage* p, int readdata, int img_select)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

    // extract the image attributes from the file header
    int				end_parse = 0;
	int				dataLength = 0;
	long			offset(0), x, y, n(0);
	char			aline[MAXLINELEN], *var;
	
//	while ( end_parse == 0 && fimg->getline(aline, MAXLINELEN) != NULL ) {
	while ( end_parse == 0 ) {
		fimg->getline(aline, MAXLINELEN);
		if ( !fimg->good() ) end_parse = 1;
    	if ( ( var = strstr(aline,"image list") ) != NULL )
    	    n++; 	    	    	    	    	    	    // Image
    	if ( ( var = strstr(aline,"Data offset") ) != NULL && offset < 1 ) 
    	    sscanf(var,"Data offset: %ld",&offset);	    	    // Offset
    	if ( ( var = strstr(aline,"Data length") ) != NULL )
    	    sscanf(var,"Data length: %d",&dataLength);	    	    // Length
    	if ( ( var = strstr(aline,"Samps/line") ) != NULL ) 
    	    sscanf(var,"Samps/line: %ld",&x);	    	    		// x
    	if ( ( var = strstr(aline,"Number of lines") ) != NULL ) 
    	    sscanf(var,"Number of lines: %ld",&y);	    	    // y
    	if ( ( var = strstr(aline,"list end") ) != NULL )
    	    end_parse = 1;	    	    	    	    	    	// End
    	if ( aline[0] != '\\' )
    	    end_parse = 1;	    	    	    	    	    	// End
	}
	
	p->size(x, y, 1);
	p->data_offset(offset);

    if ( y < 2 ) p->sizeY(dataLength/(x*n*sizeof(short)));
	
	long			i, imgstart(0), imgend(n-1);
	
	if ( img_select > -1 ) {
		if ( img_select >= n ) img_select = n - 1;
		p->images(1);
		imgstart = imgend = img_select;
	} else {
		p->images(n);
		img_select = -1;
	}
	
    p->channels(1);
	p->data_type(Short);
	p->minimum(-256*256);
	p->maximum(256*256-1);
	
	long			image_size = p->channels()*p->sizeX()*p->sizeY()*p->sizeZ()*p->data_type_size();
	unsigned char*	data = NULL;

	offset = p->data_offset() + imgstart*image_size;

	if ( readdata ) {
		p->data_alloc();
		for ( i=imgstart, data=p->data_pointer(); i<=imgend; i++, data+=image_size, offset+=image_size )
			fread_large(data, image_size, offset, fimg);
	}
	
	fimg->close();
	delete fimg;
	    
	return 0;
}

