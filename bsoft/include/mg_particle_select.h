/**
@file	mg_particle_select.h
@brief	Select particles
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20200517
**/

#include "mg_processing.h"
#include "symmetry.h"
#include "qsort_functions.h"

// Function prototypes
int_float* 	project_fom_order(Bproject* project, long& npart, int fom_index, int defocus_fit);
int_float* 	part_fom_order(Bparticle* partlist, long& npart, int fom_index);
long		part_reset_selection(Bproject* project, int flag=0);
long		part_unset_selection(Bproject* project);
long		part_invert_selection(Bproject* project);
long		part_select_micrograph(Bproject* project, Bstring& mgselect);
long		part_consolidate_selection(Bproject* project, int number);
long		part_set_selection(Bproject* project, int number);
long		part_set_sequential(Bproject* project);
long		part_set_FOM(Bproject* project, int fom_index, double fom);
long		part_deselect_from_list(Bproject* project, Bstring list);
long		part_deselect(Bproject* project, int fom_index, double fommin, double fommax=1e30);
long		part_deselect_redundant(Bparticle* partlist, double excl_dist, int part_select, int fom_index);
long		part_deselect_redundant(Bproject* project, double excl_dist, int part_select, int fom_index);
long		part_set_multi_maps(Bproject* project, int part_select, int nmaps);
long		part_set_filament_maps(Bproject* project);
long		part_reselect(Bproject* project, Bstring& tag, double reselect_min, double reselect_max);
double		part_series_comparison(Bproject* project, Bsymmetry& sym, double angle_cutoff);
long		part_select_FOM_groups(Bproject* project, int ngroups, int fom_index, int defocus_fit);
long		part_select_percentage(Bproject* project, double percentage, int fom_index, int defocus_fit);
long		part_select_best(Bproject* project, long number, int fom_index, int defocus_fit);
long		part_select_FOM_avg_std(Bproject* project, double factor, int fom_index);
long		part_select_random(Bproject* project, long number);
long		part_select_random_fraction(Bproject* project, double fraction);
long		part_select_random_group(Bproject* project, long number);
long		part_select_random_filaments(Bproject* project, int nmaps);
long		part_select_bootstrap(Bproject* project, int number);
long		part_select_random_within_view(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, int number);
long		part_select_maxsmooth(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, double threshfrac, double sigma, int fom_index);
long		part_select_best_within_view(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, int number, int fom_index);
long		part_select_to_group(Bproject* project);
long		part_select_group(Bproject* project, int group);
long		part_select_sets(Bproject* project, int size, int flag);
long		part_select_frames(Bproject* project, int frame_start, int frame_end);
long 		part_filament_direction(Bproject* project, double minpct);
long 		part_filament_direction(Bproject* project, Bproject* project2, double minpct);
long 		part_view_select(Bproject* project, View view, double angle);
long 		part_side_view_select(Bproject* project, double angle);
long 		part_euler_angle_select(Bproject* project, double* euler6);
long 		part_origin_select(Bproject* project, Vector3<double> origin, double distance);
long		part_set_first_view_in_series(Bproject* project);
long		part_set_best_view_in_series(Bproject* project, int fom_index);
long		part_select_series(Bproject* project, int size, int flag);
long		part_series_from_seed(Bproject* project, int flags);

double		part_fom_defocus_fit(Bproject* project, int fom_index, double& intercept, double& slope);
double		part_fom_defocus_fit_deselect(Bproject* project, int fom_index, double cutoff);
long		part_select_closest_to_focus(Bproject* project);
long		part_select_furthest_from_focus(Bproject* project);
long		part_delete_deselected(Bproject* project);
long		part_delete_deselected(Bparticle** partlist);
long		part_fix_defocus(Bproject* project, double max_dev);
vector<pair<Bparticle*,double>>	project_part_sort(Bproject* project, Bstring tag);
