/**
@file	mg_orient.h
@brief	Library routines for single particle analysis
@author	Bernard Heymann and David M. Belnap
@date	Created: 20010403
@date	Modified: 20190826 (BH)
**/

#include "mg_processing.h"
#include "rwimg.h"

#define TEMP_PROJ_FILE	"temp_proj.spi"

// Function prototypes
Bimage*		project_prepare_2D_references(Bproject* project, long first,
				long number, int bin=1, int ctf_action=0, double wiener=0.2);
Bimage*		img_prepare_projections(Bstring& filename, Bstring& mask_file,
				int bin, Bsymmetry& sym,
				double theta_step, double phi_step, double side_ang);
int 		project_determine_orientations(Bproject* project, Bimage* proj, Bstring& mask_file,
				int bin, Bsymmetry& sym, int part_select, vector<double>& band,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags);
int 		project_determine_orientations2(Bproject* project, Bimage* proj, Bstring& mask_file,
int bin, Bsymmetry& sym, int part_select, vector<double>& band,
double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
double shift_limit, double angle_limit, double edge_radius, int flags);
int 		project_determine_origins(Bproject* project, Bimage* proj, int bin,
				Bsymmetry& sym, int part_select, double res_lo, double res_hi,
				double shift_limit, int flags);
double		img_cross_validate(Bimage* p, Bimage* pref, Bimage* pmask, fft_plan planf);
