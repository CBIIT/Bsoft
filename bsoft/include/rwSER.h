/**
@file	rwSER.h
@brief	Header file for reading and writing FEI Ser files
@author Bernard Heymann
@date	Created: 20150130
@date	Modified: 20150204

	Format: 1D and 2D spectral and image format for electron microscopy.
**/

#include "rwimg.h"


// I/O prototypes
int 		readSER(Bimage* p, int readdata);
int 		writeSER(Bimage* p);
