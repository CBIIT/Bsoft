/**
@file	rwASCII.h
@brief	Header file for reading and writing ASCII files
@author Bernard Heymann
@date	Created: 20000318
@date	Modified: 20111217

	Format: Generic ASCII image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readASCII(Bimage* p, int readdata);
int 		writeASCII(Bimage* p);
