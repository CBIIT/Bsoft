/**
@file	rwmgRELION.h
@brief	Header file for reading and writing micrograph parameters from and to the RELION STAR format
@author Bernard Heymann
@date	Created: 20061101
@date	Modified: 20180525
**/

#include "mg_processing.h"

// Function prototypes
int			read_project_relion(Bstring& filename, Bproject* project);
int			write_project_relion(Bstring& filename, Bproject* project, int mg_select, int rec_select);
int			project_split_particles(Bproject* project, Bstring partfile, 
				Bstring path, Bstring partext);



