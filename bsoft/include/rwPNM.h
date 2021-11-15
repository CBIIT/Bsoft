/**
@file	rwPNM.h
@brief	Header file for reading and writing PNM files
@author Bernard Heymann
@date	Created: 20110317
@date	Modified: 20111217

	Format: 2D simple image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readPNM(Bimage* p, int readdata, int img_select);
int 		writePNM(Bimage* p);
