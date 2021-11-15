/**
@file	mg_extract.h
@brief	Header file for functions to extract particles from micrographs
@author Bernard Heymann
@date	Created: 20040406
@date	Modified: 20200616
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
int			particle_setup_filenames(Bparticle* part, Bstring filename, Bstring partpath);
int			filament_setup_filenames(Bfilament* fil, Bstring filename, Bstring filpath);
Vector3<double>*	vector3_spline_from_nodes(Bfilnode* fnode, long& nspline);
long		project_extract_particles(Bproject* project, double scale,
				int back_flag, int norm_flag, int fill_type, double fill, int mask_width, 
				int split, Bstring& partbase, Bstring& partpath, Bstring& partext);
long		micrograph_extract_particles(Bmicrograph* mg, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width);
long		reconstruction_extract_particles(Breconstruction* rec, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width);
Bparticle*	reconstruction_project_extract_particles(Breconstruction* rec, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width);
Bimage*		particle_extract(Bparticle* particles, Bbadarea* bad_areas, Bimage* p,
				Vector3<long> size, double scale, double bad_radius,
				int back_flag, int norm_flag, double fill, int mask_width);
Bimage*		micrograph_extract_gold(Bmicrograph* mg, Bimage* p, double radius);
Bimage*		reconstruction_extract_gold(Breconstruction* rec, Bimage* p, double radius);
Bimage*		marker_extract_gold(Bmarker* marker_list, Bimage* p, int img_num, double radius);
long		project_extract_filaments(Bproject* project, int filament_width, int axis, 
				Bstring& base, Bstring& path, Bstring& ext, int split);
long		micrograph_extract_filaments(Bmicrograph* mg, Bimage* p, double width, int axis);
long		reconstruction_extract_filaments(Breconstruction* rec, Bimage* p, double width, int axis);
Bimage*		filament_extract(Bfilament* filaments, Bimage* p, double width, int axis);
int			project_filaments_to_particles(Bproject* project, Vector3<long> box_size,
				double boxing_interval, double rise, double angle);
int			micrograph_filaments_to_particles(Bmicrograph* mg, Vector3<long> box_size, 
				double boxing_interval, double rise, double angle);
int			reconstruction_filaments_to_particles(Breconstruction* rec, Vector3<long> box_size, 
				double boxing_interval, double rise, double angle);
Bparticle*	filaments_to_particles(Bfilament* filaments,
				Vector3<double> pixel_size, Vector3<long> box_size,
				double boxing_interval, double rise, double angle);
int			project_mask_filament_particles(Bproject* project, int mask_width);
int			project_rotate_mask_filament_particles(Bproject* project, 
				int rotation_axis, int back_flag, int mask_width);

