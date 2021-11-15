/**
@file	rwmol_star.h
@brief	Read and write molecules in STAR format
@author Bernard Heymann
@date	Created: 19991113
@date	Modified: 20050908
**/

// Function prototypes
long 		read_mol_star(Bstring& filename, Bmolgroup *molgroup);
long 		write_mol_star(Bstring& filename, Bmolgroup *molgroup);

