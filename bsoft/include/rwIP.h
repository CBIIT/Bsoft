/**
@file	rwIP.h
@brief	Header file for reading and writing image plate reader files
@author Bernard Heymann
@date	Created: 20041110
@date	Modified: 20111217

	Format: 2D image file format for the image plate reader
**/

#include "rwimg.h"

#define IPSIZE	2048

// I/O prototypes
int 		readIP(Bimage* p, int readdata);
int			writeIP(Bimage* p);

