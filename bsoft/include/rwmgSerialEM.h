/**
@file	rwmgSerialEM.h
@brief	Converts between SerialEM MDOC files and a micrograph parameter file
@author	Bernard Heymann
@date	Created: 20190109
@date	Modified: 20190109
**/

#include "mg_processing.h"



// Function prototypes
int			read_project_serialem(Bstring& filename, Bproject* project, int flag);
int			write_project_serialem(Bstring& filename, Bproject* project);


