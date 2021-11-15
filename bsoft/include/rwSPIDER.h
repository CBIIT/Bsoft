/**
@file	rwSPIDER.h
@brief	Header file for reading and writing SPIDER files
@author Bernard Heymann
@date	Created: 19990410
@date	Modified: 20010928

	Format: 3D image file format for the SPIDER package
**/

#include "rwimg.h"

#define SPIDERSIZE 1024	// Minimum size of the SPIDER header (variable)

struct SPIDERhead {         	// file header for SPIDER data
    	float nslice;   		//  0	    	slices in volume (image = 1)
    	float nrow;   	    	//  1	    	rows per slice
    	float irec;   	    	//  2	    	# records in file (unused)
    	float nhistrec;	    	//  3	    	(obsolete)
    	float iform;  	    	//  4	    	file type specifier
    	float imami;  	    	//  5	    	max/min flag (=1 if calculated)
    	float fmax; 	    	//  6	    	maximum
    	float fmin; 	    	//  7	    	minimum
    	float av;   	    	//  8	    	average
    	float sig;  	    	//  9	    	standard deviation (=-1 if not calculated)
    	float ihist;  	    	// 10	    	(obsolete)
    	float nsam;   	    	// 11	    	pixels per row
    	float labrec;	    	// 12	    	# records in header
    	float iangle; 	    	// 13	    	flag: tilt angles filled
    	float phi;  	    	// 14	    	tilt angles
    	float theta;	    	// 15
    	float gamma;	    	// 16	    	(=psi)
    	float xoff; 	    	// 17	    	translation
    	float yoff; 	    	// 18
    	float zoff; 	    	// 19
    	float scale;	    	// 20	    	scaling
    	float labbyt; 	    	// 21	    	# bytes in header
    	float lenbyt; 	    	// 22	    	record length in bytes (row length)
    	float istack; 	    	// 23	    	indicates stack of images
    	float inuse;  	    	// 24	    	indicates this image in stack is used (not used)
    	float maxim;  	    	// 25	    	max image in stack used
    	float imgnum;  	    	// 26	    	number of current image
    	float unused[2];		// 27-28    	(unused)
    	float kangle; 	    	// 29	    	flag: additional angles set
    	float phi1; 	    	// 30	    	additional angles
    	float theta1;	    	// 31
    	float psi1; 	    	// 32
    	float phi2; 	    	// 33
    	float theta2;	    	// 34
    	float psi2; 	    	// 35
    	float extra[175];		// 36-210   	extra (49-75 for Jose Maria's transforms)
    	char cdat[12];  		// 211-213  	creation date
    	char ctim[8];			// 214-215  	creation time
    	char ctit[160]; 		// 216-255  	title
} ;


// I/O prototypes
int 		readSPIDER(Bimage* p, int readdata, int img_select);
int 		writeSPIDER(Bimage* p);
