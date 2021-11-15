/**
@file	mg_align.h
@brief	Header file for functions to align micrographs or coordinates from micrographs and apply the resultant transformation.
@author Bernard Heymann and Samuel Payne
@date	Created: 20000505
@date	Modified: 20211004
**/

#include "mg_processing.h"
#include "Bimage.h"
#include "Vector3.h"

// Function prototypes
int			mg_align_coordinates(Bproject* project, int refset);
int 		mg_align_micrographs(Bproject* project, int refset,
				Vector3<long> tile_size, double res_lo, double res_hi,
				double max_shift, int filter_flag, int refine_flag);
int         mg_align_feature_extraction(Bproject* project, int max_features,
				double res_low, double res_high, double thresh, int extract_method);
int			mg_apply_transform(Bmicrograph* mg_ref, Bmicrograph* mg_apply);
int			mg_merge_focal_series(Bproject* project, int use_old_origins);
double		project_write_aligned_averages(Bproject* project, Bimage* pgr,
				Bstring& imgfile, DataType datatype, Bstring& subset);
double		project_write_aligned_images(Bproject* project, Bimage* pgr,
				Bstring& imgfile, DataType datatype);
double		project_write_frame_sums(Bproject* project, Bimage* pgr,
				DataType datatype, Bstring& subset, double sampling_ratio, int flag);
double		project_align_frames(Bproject* project, int ref_img, long window, long step,
				Bimage* pgr, Bimage* pmask, Vector3<double> origin, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width,
				long bin, Bstring& subset, int flag);
double		project_align_series(Bproject* project, int ref_img, Bimage* pgr, 
				Bimage* pmask, Vector3<double> origin, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width,
				long bin, Bstring& subset, int flag);
Bimage*		mg_tomo_reconstruct2D(Bproject* project, long dimg,
				Vector3<long> size, double scale, double hi_res);
long		project_tomo_align(Bproject* project, long thickness, long iter, double dchange,
				long dimg, double resolution, double shift_limit, double edge_width, double gauss_width);
double		project_ssnr(Bproject* project);
int			project_frames_snr(Bproject* project, double res_hi,
				long window, Bstring& subset, double sampling_ratio, int flag);
int			project_frame_shift_analysis(Bproject* project, long window, double resolution);

