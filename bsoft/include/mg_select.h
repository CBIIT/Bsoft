/**
@file	mg_select.h
@brief	Header file for micrograph processing
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210515
**/

#include "ctf.h"
#include "symmetry.h"
#include "marker.h"
#include "View.h"
#include "Euler.h"
#include "Bstring.h"

// Function prototypes
Bfield* 	field_find_id(Bfield* field, Bstring& field_id);
Bmicrograph*	field_find_micrograph(Bfield* field, int mg_select, int mg_index, double mg_ang);
Bmicrograph*	field_find_micrograph_n(Bfield* field, int n);
Bmicrograph*	field_find_micrograph_by_focus(Bfield* field, int focus_opt, int index);
Bmicrograph*	field_find_micrograph_by_rotang(Bfield* field, double rotang);
Bmicrograph*	field_find_micrograph_by_tiltang(Bfield* field, double tiltang);
Bmicrograph*	field_find_zero_tilt_mg(Bfield* field);
Bmicrograph*	field_find_low_tilt_mg_with_markers(Bfield* field);
long		project_count_fields(Bproject* project);
long		project_count_micrographs(Bproject* project);
long		project_count_mg_selected(Bproject* project);
long		project_count_reconstructions(Bproject* project);
long		project_count_rec_selected(Bproject* project);
long		project_count_mg_particles(Bproject* project);
long 		project_count_mg_part_selected(Bproject* project);
long 		project_count_mg_part_selected(Bproject* project, int num_select);
long		project_count_mg_groups(Bproject* project);
long		project_count_mg_groups_selected(Bproject* project);
long		project_count_rec_particles(Bproject* project);
long		project_count_rec_part_selected(Bproject* project);
long		project_count_rec_groups(Bproject* project);
long 		project_count_mg_filaments(Bproject* project);
long 		project_count_mg_filament_nodes(Bproject* project);
long 		project_count_rec_filaments(Bproject* project);
long 		project_count_rec_filament_nodes(Bproject* project);
long 		field_count(Bfield* field);
long		field_count_micrographs(Bfield* field);
long		field_count_mg_selected(Bfield* field);
long		field_count_particles(Bfield* field);
long 		micrograph_count(Bmicrograph* mg);
long		micrograph_count_particles(Bmicrograph* mg);
long 		particle_count(Bparticle* part);
long 		particle_count_selected(Bparticle* part);
long 		filament_count(Bfilament* fil);
long 		filament_node_count(Bfilament* fil);
long		project_maximum_selection(Bproject* project);
long		project_show_selection_numbers(Bproject* project);
long		project_show_selected(Bproject* project);
long		project_show_selected_parameters(Bproject* project, int show);
long		project_show_mg_parameter(Bproject* project, Bstring& tag);
long		project_show_part_parameter(Bproject* project, Bstring& tag);
long		project_show_fom_histogram(Bproject* project, long bins,
				double min, double max);
long		project_show_mag_histogram(Bproject* project, long bins, double increment);
int			project_select_field(Bproject* project, Bstring& field_id);
int			project_select_micrograph(Bproject* project, Bstring& mg_id);
long		project_select_with_particles(Bproject* project, long part_sel);
vector<pair<Bmicrograph*,double>>	project_mg_sort(Bproject* project, Bstring tag);
