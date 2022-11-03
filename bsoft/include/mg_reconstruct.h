/**
@file	mg_reconstruct.h
@brief	Functions for reconstruction
@author	Bernard Heymann
@date	Created: 20010403
@date	Modified: 20220722
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "fft.h"
#include "symmetry.h"

// Function prototypes
int			part_ft_size(int xsize, double scale, int pad_factor);
Bimage*		particle_reconstruct(Bparticle* partlist, Bsymmetry sym, int sym_mode,
				double hi_res, Vector3<double> scale, Vector3<double> sam, Vector3<long> size,
				int ft_size, fft_plan plan, int interp_type=0,
				int ctf_action=0, double wiener=0.2, int flags=0, int first=0);
Bimage*		img_reconstruction_sum_weigh(Bimage** pacc, int imap, int nmaps, int nthreads, double hi_res);
long		project_single_particle_reconstruction(Bproject* project, 
				Bstring& maskfile, Bsymmetry& sym,
				int num_select, double hi_res, Vector3<double> scale, Vector3<long> size, 
				int pad_factor, int interp_type, int ctf_action, double wiener, int flags);
Vector3<long>	project_set_reconstruction_size(Bproject* project, Vector3<double>& sam, double scale, int twoD_flag);
int			project_configure_for_reconstruction(Bproject* project, Bstring classes, long& nmaps, int& nthreads);
Bimage* 	project_reconstruct_2D(Bproject* project, Bstring file_name, int transform_output);
Bimage* 	project_reconstruct_2D_fast(Bproject* project, Bstring file_name);
Bimage*		project_reconstruct_3D(Bproject* project, long selnum, int calcfom, Vector3<long> size);
Bimage*		project_reconstruct_3D(Bproject* project, long selnum, Vector3<long> size, double resolution);
Bimage* 	project_back_projection(Bproject* project, long num_select,
				Vector3<long> map_size, Vector3<double> sam, double scale,
				double resolution, fft_plan planf, fft_plan planb);
Bimage*		img_backprojection_accumulate(Bimage** pacc, int imap,
				int nmaps, int nthreads);


