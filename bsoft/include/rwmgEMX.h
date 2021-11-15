/**
@file	rwmgEMX.h
@brief	Reads and writes micrograph exchange files
@author Bernard Heymann
@date	Created: 20130201@date	
@date	Modified: 20130430
**/

#include "mg_processing.h"

// Function prototypes
int			read_project_emx(Bstring& filename, Bproject* project, Bstring& xsdfile);
int			read_project_emx(Bstring* file_list, Bproject* project, Bstring& xsdfile);
int			write_project_emx(Bstring& filename, Bproject* project, int mg_select, int rec_select);

