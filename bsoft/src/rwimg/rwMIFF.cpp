/**
@file	rwMIFF.cpp
@brief	Functions for reading and writing ImageMagick MIFF files
@author Bernard Heymann
@date	Created: 19990321
@date 	Modified: 20120211
**/

#include "rwMIFF.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an Image Magick map image file format.
A 2D image format used for X-windows display.
	The header is text of variable size ending in a colon, ':'.
	A binary colour table may follow the header.
	Files can be directly concatenated for a multi-image file.
	File format extensions:  	.miff
	Byte order determination:	none.
	Data types: 				depth 8 = byte, depth 16 = short.
@param	*p			the image structure.
@param 	readdata		flag to activate reading of image data.
@param 	img_select		image selection in multi-image file (-1 = all images).
@return	int 				0.
**/
int 	readMIFF(Bimage* p, int readdata, int img_select)
{
    ifstream        fimg(p->file_name());
    if ( fimg.fail() ) return -1;
		
	int					depth(8), colors(0), rle(0);
	long 				i, j, k, c(3);
	long				n, nimg(0), imgsize(0), datatypesize(1);
	char				aline[MAXLINELEN], header[10000], *tag;
	unsigned char**		buf = new unsigned char*[10000];
	unsigned char*		lut = NULL;
	Vector3<long>		size(1,1,1), pagesize(1,1,1);
	
//	memset(aline, 0, MAXLINELEN);
//	memset(header, 0, 10000);
//	memset(buf, 0, 10000);
	
	p->compound_type(TRGB);
	
	// Read multiple images
	while ( fimg.getline(header, MAXLINELEN) ) {
		if ( strstr(header, "id=ImageMagick") ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readMIFF: image " << nimg << endl;
			
			while ( fimg.getline(aline, MAXLINELEN) && !strchr( aline, ':' ) )
				strcat( header, aline );
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readMIFF: " << header << endl;
			if ( ( tag = strstr( header, "columns" ) ) )
				sscanf( tag, "columns=%ld", &size[0] );
			if ( ( tag = strstr( header, "rows" ) ) )
				sscanf( tag, "rows=%ld", &size[1] );
			if ( ( tag = strstr( header, "page" ) ) )
				sscanf( tag, "page=%ldx%ld", &pagesize[0], &pagesize[1] );
			if ( ( tag = strstr( header, "depth" ) ) )
				sscanf( tag, "depth=%d", &depth );
			if ( ( tag = strstr( header, "colors" ) ) )
				sscanf( tag, "colors=%d", &colors );
			if ( ( tag = strstr( header, "compression=RunlengthEncoded" ) ) ) rle = 1;
			if ( rle )
				cerr << "Warning: RLE for MIFF files not supported!" << endl;
			switch ( depth ) {
				case 8: p->data_type(UCharacter); break;
				case 16: p->data_type(UShort); break;
				case 32: p->data_type(Integer); break;
				default: break;
			}
			if ( ( tag = strstr( header, "class" ) ) ) {	// Default is DirectClass
				if ( strstr(tag, "PseudoClass") ) {
					c = 1;
					lut = new unsigned char[3*colors];
//					fread( lut, 3*colors, 1, fimg );
					fimg.read((char *)lut, 3*colors);
				}
			}
			if ( nimg == 0 ) {
				datatypesize = p->data_type_size();
				imgsize = size[0]*size[1]*c;
				p->data_offset(fimg.tellg());
			}
			if ( readdata ) {
				buf[nimg] = new unsigned char[imgsize];
//				fread( buf[nimg], imgsize, datatypesize, fimg );
				fimg.read((char *)buf[nimg], imgsize*datatypesize);
			} else {
//				fseek(fimg, imgsize*datatypesize, SEEK_CUR);
				fimg.seekg(imgsize*datatypesize, ios_base::cur);
			}
			nimg++;
		}
	}
	
	fimg.close();
	
	p->size(size);
	p->page_size(pagesize);
	p->channels(3);
	
	if ( img_select > -1 ) {
		if ( img_select >= nimg ) img_select = nimg - 1;
		nimg = 1;
	} else img_select = -1;
	
	p->images(nimg);
	
	if ( !readdata ) {
		delete[] buf;
		if ( lut ) delete [] lut;
		return 0;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG rwMIFF: Packing " << p->images() << " 2D images of size " << p->size() << endl;

	imgsize = size[0]*size[1]*p->channels();
	
	p->data_alloc();
	
	for ( n=i=j=0; n<nimg; n++ ) {
		if ( img_select < 0 || img_select == n ) {
			if ( colors ) {
				for ( j=0; j<size[0]*size[1]; j++ )
					for ( c=0, k=3*buf[n][j]; c<3; c++, i++, k++ ) p->set(i, lut[k]);
			} else {
				for ( j=0; j<imgsize; j++, i++ ) p->set(i, buf[n][j]);
			}
		}
		delete [] buf[n];
	}
	delete [] buf;
	if ( lut ) delete [] lut;
	
	return 0;
}

/**
@brief	Writing an Image Magick map image file format.
A 2D image format used for X-windows display.
@param	*p			the image structure.
@return	int 				0.
**/
int 	writeMIFF(Bimage* p)
{
	p->simple_to_rgb();
	
    ofstream        fimg(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long   typesize = p->data_type_size();
	long   imgsize = p->sizeX()*p->sizeY()*p->channels()*typesize;
	long	z, n;
	unsigned char*	data = p->data_pointer();
	
	for ( n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++, data+=imgsize ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG writeMIFF: Writing image " << n << " slice " << z << endl;
			
			fimg << "id=ImageMagick" << endl;
			fimg << "class=DirectClass" << endl;
			fimg << "columns=" << p->size()[0] << "  rows=" << p->size()[1] << "  depth=" << typesize*8 << endl;
			fimg << "page=" << p->size()[0] << "x" << p->size()[1] << "+0+0" << endl;
			fimg << "\f" << endl << ":" << endl;
			
			fimg.write((char *)data, imgsize);
		}
	}
		
	fimg.close();
	
	return 0;
}
