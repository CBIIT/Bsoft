/**
@file	rwmg.h
@brief	Header file for reading and writing micrograph parameters
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210412
**/

#include "mg_processing.h"
#include "mg_select.h"

#ifndef _RWMG_
// Function prototypes
Bproject*	read_project(const char* filename, int flags=0);
Bproject*	read_project(Bstring& filename, int flags=0);
Bproject*	read_project(Bstring& filename, Bstring& xsdfile, int flags=0);
Bproject*	read_project(Bstring* file_list, int flags=0);
Bproject*	read_project(Bstring* file_list, Bstring& xsdfile, int flags=0);
long		append_project(Bproject* project, Bstring* file_list, Bstring& xsdfile, int flags);
Bparticle*	read_particle(Bstring& filename, Bparticle** partlist, FOMType fom_tag[NFOM]);
int			write_project(const char* filename, Bproject* project);
int			write_project(Bstring& filename, Bproject* project, 
				int mg_select, int rec_select);
int			write_project(Bstring& filename, Bproject* project, int flags=0);
long		write_particle_list(Bstring& filename, Bproject* project, int flags=0);
int			write_particle(Bstring& filename, Bparticle* part, int euler_flag, int omega_flag, FOMType fom_tag[NFOM]);
Bstring		ppx_filename(Bstring& id, int part_id);
int			ppx_exists(Bparticle* part, int flag);
int			ppx_check(Bparticle* part, FOMType fom_tag[NFOM]);
int			ppx_check(Bparticle* part);
Bparticle*	project_find_particle(Bproject* project, Bstring& fn);
long		project_update_from_ppx(Bproject* project);
long		project_list_ppx(Bproject* project, int flag);
int			project_split_write(Bstring& filename, Bproject* project);
int			project_split_field_write(Bproject* project);
#define _RWMG_
#endif



