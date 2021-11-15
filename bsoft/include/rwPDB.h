/**
@file	rwPDB.h
@brief	Header file for reading and writing PDB coordinate files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20000729
	Format: Protein data bank atomic coordinate file format
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readPDB(Bstring& filename, Bmolgroup* molgroup);
int 	writePDB(Bstring& filename, Bmolgroup* molgroup);
