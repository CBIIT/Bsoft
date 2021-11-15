/**
@file	rwPhylip.h
@brief	Header file for reading and writing Phylip sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20030308

	Format: Protein sequence file format
**/

#include "rwmolecule.h"

// I/O prototypes
int			readPhylip(Bstring& filename, Bmolgroup* molgroup);
int			writePhylip(Bstring& filename, Bmolgroup* molgroup);
