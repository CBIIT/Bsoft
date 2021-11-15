/**
@file	rwBrookhavenSTEM.h
@brief	Header file for reading and writing Brookhaven STEM files
@author Bernard Heymann
@date	Created: 20050729
@date	Modified: 20111217

	Format: 2D STEM image file format
**/

#include "rwimg.h"

#define BrookhavenSTEMSIZE    4096	// Text header size

// I/O prototypes
int 		readBrookhavenSTEM(Bimage* p, int readdata);
int 		writeBrookhavenSTEM(Bimage* p);
