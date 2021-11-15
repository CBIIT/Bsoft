/**
@file	mg_tomography.h
@brief	Header file for single particle image processing in Bsoft
@author Bernard Heymann
@date	Created: 20020416
@date	Modified: 20210706
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "Matrix3.h"

#ifndef _mg_tomo_
#define _mg_tomo_

// Function prototypes
int			mg_marker_init(Bmicrograph* mg, Bmarker* model);
int			mg_marker_update(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref, int update_location);
Vector3<double>	mg_location_from_3D_model(Vector3<double> loc, Matrix3 mat,
				Vector3<double> ori);
Vector3<double>	mg_location_from_3D_model(Vector3<double> loc, Matrix3 mat,
				Vector3<double> ori, Vector3<double> scale);
Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Vector3<double> ori3D,
				Matrix3 mat, Vector3<double> ori);
Vector3<double>	mg_location_from_3D_model(Vector3<double> loc3D, Vector3<double> ori3D,
				Matrix3 mat, Vector3<double> ori, Vector3<double> scale);
Bimage*		img_gold_particle(Vector3<long> size, double radius);
Bmarker*	img_find_gold_particles(Bimage* p, int img_select, double radius, long edge, double cutoff);
int			project_mass_normalize(Bproject* project, double avg, double std, int norm_type,
				DataType datatype, int setinputZslices, int setoutputZslices,
				double cutmin, double cutmax, double replace_threshold);
Bplot*		project_intensity_plot(Bproject* project);
double		project_fit_intensities(Bproject* project, Bplot* plot=NULL, int flag=0);
double		project_thickness_to_lambda_ratio(Bproject* project);
double		project_lambda(Bproject* project, double thickness);
double		project_thickness(Bproject* project, double lambda);
int			project_sort_markers_by_id(Bproject* project);
double		project_tomo_residuals(Bproject* project, int show);
double		mg_tomo_residuals(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref);
double		project_tomo_errors(Bproject* project);
int			project_calculate_model(Bproject* project);
int			project_generate_markers(Bproject* project);
long		img_erase_markers(Bimage* p, Bmarker* mark, double marker_radius);
Bimage*		mg_erase_markers(Bmicrograph* mg, double marker_radius);
int			project_erase_markers(Bproject* project, double marker_radius);
int			mg_reset_model(Bmicrograph* mg, Bmarker* model);
long		project_find_markers(Bproject* project, long edge, int add);
long		mg_find_markers(Bmicrograph* mg, long edge, int add);
Bimage*		mg_composite_particle(Bmicrograph* mg, Vector3<long> size);
long		project_count_markers(Bproject* project);
long		project_show_markers(Bproject* project);
long		project_show_errors(Bproject* project, double error_cutoff);
long		project_deselect_markers(Bproject* project, Bstring& deselect_list);
long		project_delete_markers(Bproject* project, Bstring& delete_list);
long		project_renumber_markers(Bproject* project);
long		project_mg_select(Bproject* project, Bstring& mg_select);
long		project_mg_exclude(Bproject* project, Bstring& mg_exclude);
int			project_set_tilt_axis(Bproject* project, double tilt_axis);
int			project_invert_tilt_axis(Bproject* project);
int			project_set_tilt_angles(Bproject* project, double tilt_start, double tilt_step);
int			project_mg_tilt_to_matrix(Bproject* project);
int			project_calculate_angles(Bproject* project);
int			project_invert_matrices(Bproject* project);
int			project_mg_marker_select(Bproject* project, double fom);
int			project_set_marker_radius(Bproject* project, double mark_radius);
int			project_check_markers(Bproject* project, int flags);
int			project_fix_markers(Bproject* project);
int			img_clear_extraneous_areas(Bimage* p, double tilt_axis, double tilt_angle,
				long thickness, double width);
int			micrograph_clear_extraneous_areas(Bmicrograph* mg, Bimage* p, long thickness, double width);
int			project_clear_extraneous_areas(Bproject* project, long thickness, double width);
double		project_marker_rotation_axis(Bproject* project);
int			project_merge_rec_markers(Bproject* project);

#endif

