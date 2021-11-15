/**
@file	rwClustal.h
@brief	Header file for reading and writing Clustal sequence files
@author Bernard Heymann
@date	Created: 20030309
@date	Modified: 20030309

	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int 	readClustal(Bstring& filename, Bmolgroup* molgroup);
int 	writeClustal(Bstring& filename, Bmolgroup* molgroup);
