/**
@file	rwDX.h
@brief	Header file for reading and writing OpenDX files
@author Bernard Heymann
@date	Created: 20150319
@date	Modified: 20150319

	Format: OpenDX image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readDX(Bimage* p, int readdata);
int 		writeDX(Bimage* p);
