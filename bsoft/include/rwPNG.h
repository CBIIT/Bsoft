/**
@file	rwPNG.h
@brief	Header file for reading and writing PNG files
@author Bernard Heymann
@date	Created: 20041223
@date	Modified: 20041223

	Format: PNG 2D image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readPNG(Bimage* p, int readdata);
int 		writePNG(Bimage* p);
