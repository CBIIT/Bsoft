/**
@file	rwIP.cpp
@brief	Functions for reading and writing image plate reader files
@author Bernard Heymann
@date	Created: 20041110
@date 	Modified: 20111217
**/

#include "rwIP.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an image plate reader file format.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
An image format used in the Ditabis image plate reader.
	The header is text with tag-value pairs.
	Required tags:
		CREATED = time and date string
		HEADER = length of header, typically 2048
		XPIXEL = x size
		YPIXEL = y size
		BYTE PER PIXEL = 2,3,4
		XRESOLUTION = x pixel size in um
		YRESOLUTION = y pixel size in um
		THUMB-NAIL-ZOOM = zoom factor for preview, typically 10
	Additional tags:
		MAGNIFICATION = microscope magnification
		OFFSET = zero dose value
		OFFSET CORRECTION = flag to indicate correction
		GAIN = gain setting
		LASER = laser setting
		PARAMS = settings file
		FORMAT = format file
		CHANNEL = data channel
		MICRON-MARK = position, length and text of micron mark overlay
	Final tag:
		COMMENT = comment, always last, terminated by NULL.
	File format extensions:  	.IPL, .IPH, .IPR, .IPC
	Byte order:					Always little endian because the IP 
								software runs on Intel machines and
								there is no way to detect byte order from
								the text header.
	Data type: 					short.
**/
int 	readIP(Bimage* p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	
	
    // extract the image attributes from the file header
    int				bpp = 2, sb = 0;
	long			offset(IPSIZE), x, y;
	char			aline[MAXLINELEN] = "", *value;
	Vector3<double>	sampling(1,1,1);
	
	if ( systype(0) < LittleIEEE ) sb = 1;
	
	while ( !fimg->eof() && !strstr(aline, "COMMENT") ) {
    	fimg->getline(aline, MAXLINELEN);
		value = strchr(aline, '=');
		if ( value ) value += 1;
		if ( strstr(aline, "HEADER") )
			sscanf(value, "%ld", &offset);					// Header size
		if ( strstr(aline, "YPIXEL") )
    	    sscanf(value, "%ld", &x);	    	    	// x
		if ( strstr(aline, "XPIXEL") )
    	    sscanf(value, "%ld", &y);					// y
		if ( strstr(aline, "BYTE PER PIXEL") )
    	    sscanf(value, "%d", &bpp);						// bytes per pixel
		if ( strstr(aline, "XRESOLUTION") )
    	    sscanf(value, "%lf", &sampling[0]);	    	   	// x unit
		if ( strstr(aline, "YRESOLUTION") )
    	    sscanf(value, "%lf", &sampling[1]);	    	    // y unit
	}
	p->data_offset(offset);

	p->label(0);
	while ( !fimg->eof() && strchr(aline, '\n') ) {
    	fimg->getline(aline, MAXLINELEN);
		p->label(p->label() + aline);
	}
	
	p->size(x, y, 1);
	p->images(1);
	p->channels(1);
	p->data_type(UShort);
	if ( bpp > 2 ) p->data_type(Integer);
	p->sampling(sampling);
	
	
	if ( readdata ) {
		p->data_alloc();
		fimg->read((char *)p->data_pointer(), p->alloc_size());
		if ( fimg->fail() ) return -3;
		if ( sb ) swapbytes(p->alloc_size(), p->data_pointer(), p->data_type_size());
	}
	
	fimg->close();
	
	delete fimg;
		  
	return 0;
}

/**
@brief	Writing  an image plate reader file format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
A 2D image format used with the image plate reader.
**/
int			writeIP(Bimage* p)
{
	if ( p->data_type() < Integer )
		p->change_type(UShort);
	else
		p->change_type(Integer);
	
	char*		header = new char[IPSIZE];
	memset(header, 0, IPSIZE);
	p->data_offset(IPSIZE);

//	int			sb = 0;
//	if ( systype(0) < LittleIEEE ) sb = 1;
	
	snprintf(header, IPSIZE, "DITABIS micron Data File\r\n"
		"CREATED = %s"
		"HEADER = %ld\r\n"
		"YPIXEL = %ld\r\n"
    	"XPIXEL = %ld\r\n"
    	"BYTE PER PIXEL = %ld\r\n"
    	"XRESOLUTION = %g\r\n"
    	"YRESOLUTION = %g\r\n"
    	"THUMB-NAIL-ZOOM = 10\r\n"
    	"MAGNIFICATION = 1\r\n"
    	"OFFSET = 0\r\n"
    	"OFFSET CORRECTION = NO\r\n"
//    	"GAIN = 10000\r\n"
//    	"LASER = 50\r\n"
//    	"PARAMS = Standard.set\r\n"
//  	"FORMAT = Standard.fmt\r\n"
//    	"CHANNEL = PMT LOWSENS\r\n"
//    	"MICRON-MARK = \r\n"
    	"COMMENT\r\n%s",
		asctime(p->get_localtime()), p->data_offset(), p->sizeX(), p->sizeY(), 
			p->data_type_size(), p->sampling(0)[0], p->sampling(0)[1], p->label().c_str());

	ofstream		fimg(p->file_name());
	if ( fimg.fail() ) return -1;
	
	fimg.write((char *)header, p->data_offset());
	if ( p->data_pointer() ) fimg.write((char *)p->data_pointer(), p->alloc_size());
	
	fimg.close();
	
	delete[] header;
		
	return 0;
}

