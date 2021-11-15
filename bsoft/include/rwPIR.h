/**
@file	rwPIR.h
@brief	Header file for reading and writing PIR sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20000729
	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int			readPIR(Bstring& filename, Bmolgroup* molgroup);
int			writePIR(Bstring& filename, Bmolgroup* molgroup);
