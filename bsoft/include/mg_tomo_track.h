/**
@file	mg_tomo_track.h
@brief	Functions to track fiducial markers in a tomographic series of images
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20210112
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "rwimg.h"

// Function prototypes
double		project_find_tilt_axis(Bproject* project, double tilt,
				double axis_start, double axis_end, double axis_step,
				double hi_res, double lo_res, double shift_limit);
int			project_track_markers(Bproject* project, double hi_res, double lo_res, 
				double shift_limit, double thickness, int max_cycle,
				double target, int cc_type, int recenter, Bstring paramfile);
int			project_track_markers_dual(Bproject* project, double hi_res, double lo_res, 
				double shift_limit, double thickness, int max_cycle,
				double target, int cc_type, int recenter, Bstring paramfile);
int			project_refine_markers(Bproject* project, double hi_res, double lo_res);
int			project_refine_one_marker(Bproject* project, int id, double hi_res, double lo_res);
int			mg_refine_markers(Bmicrograph* mg, Bimage* pgold, double hi_res, double lo_res);
double		project_refine(Bproject* project, int iter, double tol, Bstring refop);
double		project_refine_z(Bproject* project);
double		project_refine(Bproject* project, int do_view, int do_origin, int do_scale);
double		project_tilt_axis_from_markers(Bproject* project);
double		mg_marker_shift(Bmicrograph* mg, Bimage* pgold, double hi_res, double lo_res, double shift_limit);
double		mg_marker_shift(Bmicrograph* mg, Bimage* pgold, double hi_res,
				double lo_res, double shift_limit, fft_plan planf, fft_plan planb);
int			project_transfer_seed(Bproject* project,
				double rot_start, double rot_end, double rot_step, 
				double hi_res, double lo_res, double shift_limit);
int			project_transform_dual(Bproject* project);
double		project_dual_zcompare(Bproject* project);

