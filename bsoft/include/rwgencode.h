/**
@file	rwgencode.h
@brief	Header file for reading a genetic code from a file
@author Bernard Heymann
@date	Created: 20030316
@date	Modified: 20060710
**/

#include "Bstring.h"

/* Function prototypes */
Bstring		get_genetic_code(Bstring& filename);
int 		write_genetic_code(Bstring& filename, Bstring& gc);
int			index_from_codon(const char codon[3]);


