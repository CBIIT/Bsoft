/**
@file	rwmol_text.h
@brief	Read and write molecules in plain text
@author Bernard Heymann
@date	Created: 20050419
@date	Modified: 20050419
**/

// Function prototypes
int 		read_mol_text(Bstring& filename, Bmolgroup *molgroup);
int 		write_mol_text(Bstring& filename, Bmolgroup *molgroup);

