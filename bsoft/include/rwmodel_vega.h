/**
@file	rwmodel_vega.h
@brief	Header file for reading and writing Vega model parameters
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20080408
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_vega(Bstring* file_list);
int			write_model_vega(Bstring& filename, Bmodel* model);


