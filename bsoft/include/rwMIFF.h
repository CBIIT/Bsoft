/**
@file	rwimg.h
@brief	Header file for reading and writing ImageMagick MIFF files
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20111217

	Format: 2D image file format for the Image Magick package
**/

#include "rwimg.h"

// I/O prototypes
int 		readMIFF(Bimage* p, int readdata, int img_select);
int 		writeMIFF(Bimage* p);
