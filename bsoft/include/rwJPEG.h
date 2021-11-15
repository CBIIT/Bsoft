/**
@file	rwJPEG.h
@brief	Header file for reading and writing JPEG files
@author Bernard Heymann
@date	Created: 20030510
@date	Modified: 20210402

	Format: 2D image file format (JPEG)
**/

#include "rwimg.h"

// I/O prototypes
int 		readJPEG(Bimage* p, int readdata);
int 		writeJPEG(Bimage* p, int quality);
