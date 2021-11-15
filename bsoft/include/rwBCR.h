/**
@file	rwBCR.h
@brief	Header file for reading and writing AFM BCR-STM files
@author Bernard Heymann
@date	Created: 20170214
@date	Modified: 20170214

	Format: BCR-STM image file format
**/

#include "rwimg.h"

#define BCRSIZE   2048	// Size of the BCR header (2x for unicode)

// I/O prototypes
int 		readBCR(Bimage* p, int readdata);
int 		writeBCR(Bimage* p);
