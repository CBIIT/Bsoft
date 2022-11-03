/**
@file	rwmodel_pdb.h
@brief	Header file for reading and writing molecular model parameters
@author Bernard Heymann
@date	Created: 20211231
@date	Modified: 20220929
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_pdb(Bstring* file_list, Bstring& paramfile);
int			write_model_pdb(Bstring& filename, Bmodel* model, int split);


