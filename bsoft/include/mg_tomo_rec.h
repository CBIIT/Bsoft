/**
@file	mg_tomo_rec.h
@brief	Functions to do a tomographic reconstruction
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20180629
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
Bimage*		project_tomo_reconstruct(Bproject* project, double hi_res,
				double scale, Vector3<long> size, int interp_type, int pad_factor,
				double edge_width, double marker_radius, 
				int fill_type, double fill, int action, double wiener);
Bimage*		project_fourier_reconstruction_slab(Bproject* project, double hi_res,
				double scale, Vector3<long> size, int slab_start, int slab_end,
				double marker_radius, int fill_type, double fill, int action, double wiener);
long		project_tomo_reconstruct_particles(Bproject* project, 
				double resolution, int interp_type, int pad_factor, 
				int ctf_action, double wiener, Bsymmetry& sym,
				Bstring& partbase, Bstring& partpath, Bstring& partext);
long		img_pack_2D_in_recip_space_slab(Bimage* p, Bimage* prec, 
				long zsize, long slab_start, float* weight, float* weight2,
				double hi_res, Matrix3 mat, double scale);
int			img_backtransform_slices(Bimage* p);
int			img_phase_shift_slab_to_origin(Bimage* p, int zsize, int slab_start);
int			mg_fft_write(Bproject* project, Vector3<int> size, double scale, int pad_factor, 
				DataType datatype, double marker_radius, int fill_type, double fill);
Bimage*		img_extract_ytile(Bstring* file_list, int ystart, int ysize);
Bimage*		img_backtransform_z_on_disk(Bstring* file_list, Bstring& recfile, 
				DataType datatype, double avg, double std, double cutmin, double cutmax);
int			img_backtransform_z_lines(Bimage* p);
Bimage*		project_missing_mask(Bproject* project, Vector3<long> size,
				Vector3<double> origin, double hi_res, double scale) ;


