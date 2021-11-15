/**
@file	rwmgXML.h
@brief	Reads and writes micrograph XML files
@author Bernard Heymann
@date	Created: 20050920
@date	Modified: 20130827
**/

#include "mg_processing.h"

// Function prototypes
int			read_project_xml(Bstring& filename, Bproject* project);
int			read_project_xml(Bstring* file_list, Bproject* project);
Bparticle*	read_particle_xml(Bstring& filename, Bparticle** partlist, FOMType fom_tag[NFOM]);
int			write_project_xml(Bstring& filename, Bproject* project, int mg_select, int rec_select);
int			write_particle_xml(Bstring& filename, Bparticle* part, int euler_flag, int omega_flag, FOMType fom_tag[NFOM]);

