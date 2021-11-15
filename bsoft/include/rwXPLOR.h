/**
@file	rwXPLOR.h
@brief	Header file for reading and writing XPLOR reflection files
@author Bernard Heymann
@date	Created: 19981229
@date	Modified: 20000715

	Format: 3D structure factor file format for the program XPLOR
**/

#include "rwimg.h"

// I/O prototypes
int 	readXPLOR(Bimage* p, int readdata);
int 	writeXPLOR(Bimage* p);
