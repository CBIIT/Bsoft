/*
	img_radon.h
	Header file for Radon transform functions.
	Author: Bernard Heymann, Salvatore Lanzavecchia and Francesca Cantele
	Created: 20010519
	Modified: 20160728 (BH)
*/

#include "mg_processing.h"
#include "symmetry.h"
#include "rwimg.h"

/* Function prototypes */
Bimage *img_radon_transform(Bimage *p, int type, int nkernel, int kernel_power, int padd, int ntheta);
Bimage *img_radon_inverse_transform(Bimage *p, int type, int nkernel, int kernel_power, int padd);
int img_radon_pocs_filter(Bimage *p, int n_cyc_out, int n_cyc_in, double rad_3D, double rad_plane, int support, Bimage *pmask);
Bimage *img_radon_reconstruction(Bproject *project, Bsymmetry& sym, Bstring& file_mask, int rec_size, int ntheta, int table_size, double threshold, Vector3<double> origin, int nkernel, int kernel_power);
int		img_resize_to_next_power2(Bimage* p, int fill_type, double fill);
