/**
@file	rwmgSTAR.h
@brief	Header file for reading and writing micrograph parameters from and to the STAR format
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210330
**/

#include "mg_processing.h"
#include "mg_select.h"
#include "rwstar.h"

// Function prototypes
int			read_project_star(Bstring& filename, Bproject* project);
int			read_project_star2(Bstring& filename, Bproject* project);
int			write_project_star(const char* filename, Bproject* project, int mg_select, int rec_select);
int			write_project_star(Bstring& filename, Bproject* project, int mg_select, int rec_select);
int			write_project_star2(Bstring& filename, Bproject* project, int mg_select, int rec_select);
int			mg_star_update_tags(Bstar* star);



