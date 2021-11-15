/**
@file	mg_img_proc.h
@brief	Header file for image processing from micrograph structures
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210722
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
Bstring		image_type(Bstring type);
Bproject*	project_create_from_image(Bimage* p, Bstring type);
Bproject*	project_create_from_images(Bstring* file_list, Bstring type);
double		micrograph_intensity(Bmicrograph* mg, Bimage* p, int flag);
double		micrograph_intensity(Bmicrograph* mg, int flag=0);
int			project_mg_avg_intensities(Bproject* project);
int			project_catenate_micrographs(Bproject* project);
Bimage*		particle_read_img(Bparticle* part, int readflag);
int			particle_write_img(Bparticle* part, Bimage* p, int compression);
int			project_check_particles(Bproject* project);
Vector3<long>	micrograph_get_size(Bmicrograph* mg);
Vector3<double>	micrograph_get_nominal_origin(Bmicrograph* mg);
int			project_set_nominal_mg_origins(Bproject* project);
int			project_reset_origins(Bproject* project);
int			project_set_part_img_origins(Bproject* project);
long		project_delesect_edge_particles(Bproject* project);
int			project_get_part_box_size(Bproject* project);
Vector3<long>	particle_get_box_size(Bparticle* part);
int			project_write_particle_classes(Bproject* project);
int			project_trim_class_averages(Bproject* project, Bstring& list);
int			project_flip_particle_coordinates(Bproject* project, int flip);
int			project_set_views_from_images(Bproject* project);
int			project_set_views_in_images(Bproject* project);
int			project_bin_micrographs(Bproject* project, int bin, Bstring mgpath, Bstring partpath);
int			project_change_pixel_size(Bproject* project,
				Vector3<double> new_pixel_size,
				Bstring mgpath, Bstring partpath);



