/**
@file	ps_micrograph.h
@brief	Header file for postscript tools dealing with micrograph particle sets
@author Bernard Heymann
@date	Created: 20011127
@date	Modified: 20210706
**/

#include <fstream>
#include "mg_processing.h"
#include "mg_select.h"

//Function prototypes
int			ps_mg_origins(Bstring& filename, Bstring& title, Bproject* project);
int			ps_mg_particle_positions(Bstring& filename, Bstring& title, Bproject* project);
int 		ps_particle_views_origins(Bstring& filename, Bstring& title, 
					Bstring& symmetry_string, Bproject* project, int selection);
int 		ps_particle_phi_theta(Bstring& filename, Bstring& title, 
					Bproject* project, int selection);
int			ps_origins(Bstring& filename, Bstring& title, Bproject* project, int flags);
int			ps_origins(ofstream* fps, Bstring& title, Bproject* project, int flags);
int			ps_part_fom_histogram(Bstring& filename, Bproject* project);
int			ps_defocus_histogram(Bstring& filename, Bproject* project);
int			ps_astigmatism_plot(Bstring& filename, Bproject* project);
int			ps_class_average_fom(ofstream* fps, Bproject* project);

