/**
@file	rwND2.h
@brief	Header file for reading Nikon ND2 files
@author Bernard Heymann
@date	Created: 20210628
@date	Modified: 20210628

	Format: Light microscopy image format for Nikon microscopes
**/

#include "rwimg.h"

// I/O prototypes
int 		readND2(Bimage* p, int readdata, int img_select);
