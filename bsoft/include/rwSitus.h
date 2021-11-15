/**
@file	rwSitus.h
@brief	Header file for reading and writing Situs files
@author Bernard Heymann
@date	Created: 20150316
@date	Modified: 20150316

	Format: Situs image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readSitus(Bimage* p, int readdata);
int 		writeSitus(Bimage* p);
