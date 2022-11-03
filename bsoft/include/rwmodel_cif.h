/**
@file	rwmodel_cif.h
@brief	Read and write molecules in CIF format
@author Bernard Heymann
@date	Created: 19991113
@date	Modified: 20220103
**/

#include "rwmodel.h"

// Function prototypes
Bmodel*		read_model_cif(Bstring& filename, Bstring& paramfile);
long 		write_model_cif(Bstring& filename, Bmodel* model);

