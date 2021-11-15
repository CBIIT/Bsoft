/**
@file	rwTGA.h
@brief	Header file for reading and writing Truevision TGA files
@author Bernard Heymann
@date	Created: 20150811
@date	Modified: 20150815

	Format: Truevision image file format
**/

#include "rwimg.h"

#define TGASIZE     18		// Size of the TGA header (constant)

struct __attribute__((packed)) TGAhead {                // file header for TGA data
        char idlen;  	        //  0   ID length (0-255)
        char colmaptype;  	    //  1	Color map type: 0=none, 1=included
        char imgtype;   	  	//  2	Image type:
								//			0: No image data
								//			1: Uncompressed color-mapped
								//			2: Uncompressed true-color
								//			3: Uncompressed grayscale
								//			9: Run=length encoded color-mapped
								//			10: Run=length encoded true-color
								//			11: Run=length encoded grayscale
								//  3	Color map specification:
		short colmapindex;		//			index of first color map entry
		short colmaplength;		//			number of color map entries
        char colmapbits;	    //  		number of bits per entry (15,16,24,32)
        						//  8	Image specification:
		short xorigin;			//			left screen position
		short yorigin;			//			bottom screen position
		short width;			//			width of image
		short height;			//			height of image
		char pixeldepth;		//			number of bits per pixel (8,16,24,32)
		char imgdesc;			//			image descriptor:
								//				0-3: attribute bits per pixel (alpha/overlay)
								//				4: left-to-right ordering
								//				5: top-to-bottom ordering
								//				6,7: 0
} ;

struct __attribute__((packed)) TGAfoot {
	int extoff;				// Eextension area offset
	int devoff;				// Developer area offset
	char sig[16];			// TGA v2 signature
	char term[2];			// Terminator
} ;


// I/O prototypes
int 		readTGA(Bimage* p, int readdata);
int 		writeTGA(Bimage* p);
