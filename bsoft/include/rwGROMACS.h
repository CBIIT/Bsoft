/**
@file	rwGROMACS.h
@brief	Header file for reading and writing GROMACS coordinate files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20000729

	Format: Atomic coordinate file format for the GROMACS package
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readGROMACS(Bstring& filename, Bmolgroup* molgroup);
int 	writeGROMACS(Bstring& filename, Bmolgroup *molgroup);
