/**
@file	rwmgSTAR.h
@brief	Header file for reading and writing micrograph parameters from and to the STAR format
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220113
**/

#include "mg_processing.h"
#include "mg_select.h"

// Function prototypes
int			read_project_star(Bstring& filename, Bproject* project, int flag);
int			write_project_star(Bstring& filename, Bproject* project, int mg_select, int rec_select);



