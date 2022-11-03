/**
@file	rwmodel_xyz.h
@brief	Header file for reading and writing Kirkland's xyz atomic model parameters
@author Bernard Heymann
@date	Created: 20220426
@date	Modified: 20220426
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_xyz(Bstring* file_list, Bstring& paramfile);
int			write_model_xyz(Bstring& filename, Bmodel* model);


