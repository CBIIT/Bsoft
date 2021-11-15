/**
@file	rwGenBank.h
@brief	Header file for reading and writing GenBank sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20030308

	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readGenBank(Bstring& filename, Bmolgroup* molgroup);
int 	writeGenBank(Bstring& filename, Bmolgroup* molgroup);
