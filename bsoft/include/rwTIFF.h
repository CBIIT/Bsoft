/**
@file	rwTIFF.h
@brief	Header file for reading and writing TIFF files
@author Bernard Heymann
@date	Created: 19990509
@date	Modified: 20010412

	Format: 3D image file format (LIBTIFF)
**/

#include "rwimg.h"

// I/O prototypes
int 		readTIFF(Bimage* p, int readdata, int img_select);
int 		writeTIFF(Bimage* p, int istiled);
