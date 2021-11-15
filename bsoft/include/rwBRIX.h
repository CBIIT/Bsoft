/**
@file	rwBRIX.h
@brief	Header file for reading and writing Brix files
@author Bernard Heymann
@date	Created: 19990424
@date	Modified: 20111217

	Format: 3D crystallographic image file format for the program 'O'
**/

#include "rwimg.h"

#define BRIXSIZE    512	// Size of the Brix header (constant)


// I/O prototypes
int 		readBRIX(Bimage* p, int readdata);
int 		writeBRIX(Bimage* p);
