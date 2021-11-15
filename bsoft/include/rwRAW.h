/**
@file	rwRAW.h
@brief	Header file for reading and writing Raw files
@author Bernard Heymann
@date	Created: 19990724
@date	Modified: 20010410

	A Raw file is defined as consisting of only a block data with no additional info
	Format: Generic customizable 3D image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		readRAW(Bimage* p, int img_select);
int 		writeRAW(Bimage* p);
