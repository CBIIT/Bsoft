/**
@file	mg_particles.h
@brief	Manipulates single particle images.
@author  Bernard Heymann
@date	Created: 20080424
@date	Modified: 20220531
**/

#include "rwimg.h"

// Function prototypes 
int			project_align_particles(Bproject* project, int part_select,
				double nuavg, double nustd);
int			project_center_particles(Bproject* project, int part_select,
				double nuavg, double nustd);
int			project_set_particle_centers(Bproject* project, Bimage* pref,
				int part_select, double hires, double lores);
long		project_find_part_centers_in_mgs(Bproject* project,
				double hires, double lores, int filter_flag);
long		mg_find_part_centers(Bmicrograph* mg, 
				double hires, double lores, int filter_flag);
Bimage*		project_find_particle_centers(Bproject* project, int max_iter, 
				int part_select, double hires, double lores);
long		project_mask_particles(Bproject* project, Bimage* pmask, Bstring& partpath);
long		project_compare_particles(Bproject* project, Bproject* projcomp);
double		project_tilt_from_particle_defocus(Bproject* project);
long		project_set_particle_defocus_from_tilt(Bproject* project, double axis, double tilt);
Bimage*		project_correlation_sum(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel);
Bimage*		project_aberration_phase_difference(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel);
Bimage*		project_ewald_phase(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel);
Bimage*		project_ewald_correlation(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel);

