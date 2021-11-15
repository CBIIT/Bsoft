/**
@file	rwEER.h
@brief	Header file for reading EER files
@author Bernard Heymann
@date	Created: 20210624
@date	Modified: 20210624

	Format: Electron event representation file
**/

#include "rwimg.h"

// I/O prototypes
int 		readEER(Bimage* p, int readdata, int img_select, int supres);
