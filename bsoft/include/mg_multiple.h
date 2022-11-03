/**
@file	mg_multiple.h
@brief	Selection of single particle parameters from multiple files for classification
@author Bernard Heymann
@date	Created: 20010319
@date	Modified: 20220831
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "symmetry.h"
#include "utilities.h"
#include "timer.h"

// Function prototypes
Bproject*	project_multi_merge(Bstring* file_list, int fom_index, int flags);
int			project_multi_add(Bproject* project, Bstring* file_list, int fom_index);
Bproject*	project_multi_add_particles(Bstring* file_list);
long 		project_multi_adjust_FOM(Bproject* project_list, int fom_index);
long 		project_multi_select_best_FOM(Bproject* project_list, double fom_cut, int fom_index, int fom_def_flag);
long 		project_multi_selection_stats(Bproject* project_list);
long 		project_multi_select_low_variance(Bproject* project_list, Bsymmetry& sym,
				double origin_dev, double view_dev, double angle_dev, double mag_dev);
long 		project_multi_select_low_difference(Bproject* project1, Bproject* project2, Bsymmetry& sym,
				double origin_err, double view_err, double angle_err, double mag_err);
long 		project_multi_select_low_rmsd(Bproject* project1, Bproject* project2,
				Bsymmetry& sym, double origin_rmsd, double view_rmsd,
				double angle_rmsd, double mag_rmsd, int flag);
long 		project_multi_selection_compare(Bproject* project1, Bproject* project2);

