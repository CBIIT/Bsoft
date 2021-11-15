/**
@file	rwsymop.h
@brief	Header file for reading and writing symmetry operators
@author Bernard Heymann
@date	Created: 19990509
@date	Modified: 20070621
**/

#include "symmetry.h"
#include "rwimg.h"
#include "View.h"

// Function prototypes in symmetry.c and rwsymop.c
float* 		read_symat(Bstring& filename, int spacegroup, int& nsym);
char* 		read_symop(Bstring& filename, int spacegroup, int& nsym);
int 		write_symat(Bstring& filename, int spacegroup);
int 		write_pointgroup(Bstring& filename, Bstring& symmetry_string, View ref_view);
int 		write_pointgroup(Bstring& filename, Bsymmetry& sym, View ref_view);

