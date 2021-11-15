/**
@file	mg_tomo_resol.h
@brief	Functions to assess the resoltion of a tomographic tilt series
@author	Bernard Heymann
@date	Created: 20031205
@date	Modified: 20200723
**/

#include "mg_processing.h"
#include "ps_plot.h"
#include "rwimg.h"

// Function prototypes
Bimage*		mg_tomo_resolution(Bproject* project, int micrograph_id, double hi_res, 
				double sampling_ratio, double scale, Vector3<long> size,
				double fast_angle, int action, double wiener, double cutoff, Bstring& psfile);
vector<Bplot*>	project_tomo_resolution(Bproject* project, double hi_res,
				double sampling_ratio, double scale, Vector3<long> size,
				double fast_angle, int action, double wiener, double cutoff);
Bplot*		project_tomo_particle_resolution(Bproject* project, double hi_res,
				double sampling_ratio, double fast_angle, double cutoff);
long		img_pack_2D_into_central_section(Bimage* p, Bimage* prec, Bimage* prec2,
				long ft_size, double scale, double hi_res,
				Matrix3 matr, Matrix3 mat, int inplane);
Bplot*		plot_tilt_resolution(Bproject* project);

