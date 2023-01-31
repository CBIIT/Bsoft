/**
@file	rwmodel_cmm.h
@brief	Header file for reading and writing Chimera marker model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20221115
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_chimera(Bstring* file_list);
int			write_model_chimera(Bstring& filename, Bmodel* model, int split);


