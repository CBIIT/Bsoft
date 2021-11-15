/**
@file	rwFASTA.h
@brief	Header file for reading and writing FASTA sequence files
@author Bernard Heymann
@date	Created: 20001112 
@date	Modified: 20001112

	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int			readFASTA(Bstring& filename, Bmolgroup* molgroup);
int			writeFASTA(Bstring& filename, Bmolgroup* molgroup);
