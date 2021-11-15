/**
@file	rwHKL.h
@brief	Header file for reading and writing HKL reflection files as images
@author Bernard Heymann
@date	Created: 19981229
@date	Modified: 20000708

	Format: 3D generic structure factor file format
**/

#include "rwimg.h"

// I/O prototypes
int 	readHKL(Bimage* p, int readdata);
int 	writeHKL(Bimage* p);
