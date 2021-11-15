/**
@file	rwDI.h
@brief	Header file for reading (only) Digital Instruments files
@author Bernard Heymann
@date	Created: 19990424
@date	Modified: 20010410

	Format: 2D image file format for AFM (Digital Instruments)
**/

#include "rwimg.h"

#define DISIZE    20480	// Maximum size of the DI header (4,8,12,20K)

// I/O prototypes
int 		readDI(Bimage* p, int readdata, int img_select);
