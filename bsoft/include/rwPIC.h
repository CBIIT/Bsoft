/**
@file	rwPIC.h
@brief	Header file for reading and writing PIC's BP and BQ files
@author Bernard Heymann
@date	Created: 20000412
@date	Modified: 20111217

	Format: 2D image file format for the PIC package
**/

#include "rwimg.h"

// I/O prototypes
int 		readPIC(Bimage* p, int readdata);
int 		writePIC(Bimage* p);
