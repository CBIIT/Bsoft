/**
@file	rwEMBL.h
@brief	Header file for reading and writing EMBL sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20000729

	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readEMBL(Bstring& filename, Bmolgroup* molgroup);
int 	writeEMBL(Bstring& filename, Bmolgroup* molgroup);
