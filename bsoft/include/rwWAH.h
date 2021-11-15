/**
@file	rwWAH.h
@brief	Header file for reading and writing Wayne Hendrickson coordinate files
@author Bernard Heymann
@date	Created: 20050217
@date	Modified: 20050217

	Format: Atomic coordinate file format for PROLSQ
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readWAH(Bstring& filename, Bmolgroup* molgroup);
int 	writeWAH(Bstring& filename, Bmolgroup *molgroup);
