/**
@file	rwPIC.cpp
@brief	Functions for reading and writing PIC BP and BQ files
@author Bernard Heymann
@date	Created: 20000412
@date 	Modified: 20120409
**/

#include "rwPIC.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading a PIC image file format(s).
A 2D image format used in electron microscopy.
	Header size:				variable - determined by title or label.
	File format extensions:  	.bp, .bq
	Byte order determination:	Title length must be less than 256.
	Data types: 				BP files = byte, BQ files = float,
								(Note: float is VMS G_FLOAT).
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@return	int					error code (<0 means failure).
**/
int 	readPIC(Bimage *p, int readdata)
{
	ifstream*		fimg = new ifstream(p->file_name());
	if ( fimg->fail() ) return -1;	

	short			titlen, nbr, nbc, sb = 0, pad;
	char			title[128];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIC: Reading PIC file = " << p->file_name() << endl;
	
    // Read and swap bytes if from different-endian machine
	fimg->read( (char *)&pad, sizeof(short) ); if ( fimg->fail() ) return -2;
	fimg->read( (char *)&titlen, sizeof(short) ); if ( fimg->fail() ) return -2;
	if ( titlen > 256 ) {
		sb = 1;
		swapbytes((unsigned char *)&titlen, sizeof(short));
	}
	fimg->read( (char *)&pad, sizeof(short) );
	if ( fimg->fail() ) return -2;
	
	if ( !pad ) {
		fimg->read( (char *)&pad, sizeof(short) );
		if ( fimg->fail() ) return -2;
	}
	
	if ( titlen > 59 ) titlen = 59;
	if ( titlen > 0 ) {
		fimg->read( title, 2*titlen );
		if ( fimg->fail() ) return -2;
		title[2*titlen+1] = 0;
	}
	
	fimg->read( (char *)&pad, sizeof(short) );
	if ( fimg->fail() ) return -2;
	
	if ( !pad ) {
		fimg->read( (char *)&pad, sizeof(short) );
		if ( fimg->fail() ) return -2;
	}
	
	fimg->read( (char *)&nbr, sizeof(short) );
	if ( fimg->fail() ) return -2;
	
	fimg->read( (char *)&pad, sizeof(short) );
	if ( fimg->fail() ) return -2;
	
	if ( !pad ) {
		fimg->read( (char *)&pad, sizeof(short) );
		if ( fimg->fail() ) return -2;
	}
	
	fimg->read( (char *)&nbc, sizeof(short) );
	if ( fimg->fail() ) return -2;
	
	if ( sb ) {
		swapbytes((unsigned char *)&nbr, sizeof(short));
		swapbytes((unsigned char *)&nbc, sizeof(short));
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readPIC: sb=" << sb << endl;
		cout << "DEBUG readPIC: title: " << title << " (length = " << 2*titlen << ")" << endl;
		cout << "DEBUG readPIC: size=" << nbc << "x" << nbr << endl;
		cout << "DEBUG readPIC: offset=" << 2*titlen + 8*sizeof(short) << endl;
	}
	
	// Map the parameters
	p->size(nbc, nbr, 1);
	p->images(1);
	p->channels(1);
	p->data_type(UCharacter);
//	if ( p->file_name().contains(".bq") || p->file_name().contains(".BQ") )
	if ( p->file_name().find(".bq") != string::npos || p->file_name().find(".BQ") != string::npos )
		p->data_type(Float);
	p->label(title);
	p->data_offset(2*titlen + 8*sizeof(short));
	
	
	long			datatypesize = p->data_type_size();
	long			datasize = p->alloc_size();
	long			linesize = p->sizeX()*datatypesize;
	long			i, y, y2, pagesize;
//	char			buf[128];
	char*			data = NULL;
	
	// At least some PIC files are written with a variable page size
	// with a maximum of 2042, padded with a short in between
	if ( readdata ) {
		data = (char *) p->data_alloc();
		for ( y=0; y<p->sizeY(); y++ ) {
			for ( i=1, y2=0; y2<p->sizeX(); y2+=2042, i++ ) {
				pagesize = 2042;
				if ( pagesize > linesize - y2 ) pagesize = linesize - y2;
				fimg->read( (char *)&pad, sizeof(short) );
				if ( fimg->fail() ) return -3;
				if ( !pad ) fimg->read( (char *)&pad, sizeof(short) );
				if ( fimg->fail() ) return -3;
				fimg->read( data+y*linesize+y2, pagesize );
				if ( fimg->fail() ) return -3;
			}
		}
	}

	fimg->close();
	delete fimg;

	if ( p->data_pointer() && p->data_type() == Float ) {
		if ( verbose > 0 )
			cout << "Warning: BQ files may be VMS G_float, not supported in Bsoft" << endl;
		for ( y=0; y<datasize; y+=4 )
    	    vax2ieee(p->data_pointer()+y, 1-sb);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIC: PIC file " << p->file_name() << " read" << endl;
	
	return 0;
}

/**
@brief	Writing a PIC image file format(s).
A 2D image format used in electron microscopy.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writePIC(Bimage *p)
{
	p->change_type(UCharacter);
	p->color_to_simple();
	
	short		titlen, nbr, nbc, pad = 3;
	char		title[120];
	for ( int i=0; i<120; ++i ) title[i] = ' ';
	
	// Map the parameters
	strncpy(title, p->label().c_str(), 120);
	titlen = strlen(title)/2 + 1;
	nbc = p->sizeX();
	nbr = p->sizeY();
		
	long	datatypesize = p->data_type_size();
	long	linesize = p->sizeX()*datatypesize;
	
	ofstream		fimg(p->file_name());
	if ( fimg.fail() )  return -1;
	
	fimg.write( (char *)&pad, sizeof(short) );
	fimg.write( (char *)&titlen, sizeof(short) );
	fimg.write( (char *)&pad, sizeof(short) );
	fimg.write( title, 2*titlen );
	fimg.write( (char *)&pad, sizeof(short) );
	fimg.write( (char *)&nbr, sizeof(short) );
	fimg.write( (char *)&pad, sizeof(short) );
	fimg.write( (char *)&nbc, sizeof(short) );
	
	// At least some PIC files are written with a variable page size
	// with a maximum of 2042, padded with a short in between
	long	y, y2, pagesize;
	char*			data = (char *) p->data_pointer();
	
	pad = 3;
	for ( y=0; y<p->sizeY(); y++ ) {
		for ( y2=0; y2<p->sizeX(); y2+=2042 ) {
			if ( p->sizeX() > 2042 ) pad = y2 + 1;
			pagesize = 2042;
			if ( pagesize > linesize - y2 ) pagesize = linesize - y2;
			fimg.write( (char *)&pad, sizeof(short) );
			fimg.write( data+y*linesize+y2, pagesize );
		}
	}
	
	fimg.close();
	
	return 0;
}

