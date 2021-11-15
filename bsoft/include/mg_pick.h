/**
@file	mg_pick.h
@brief	Header file for single particle picking functions.
@author Bernard Heymann
@date	Created: 20000505
@date	Modified: 20190912
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "marker.h"

// Function prototypes
Bparticle*	particles_from_peaks(Bimage* pcc, long bin, double excl_dist, 
				double part_ori, double& fommin, double fommax=1e30, 
				long maxnum=1000000, double pix_min=2, double pix_max=10);
Bparticle*	particles_pick_cc(Bimage* p, Bimage* ptemp, Bimage* pmask, 
				double hires, double lores, double fommin, double fommax, 
				double excl_dist, long bin=1);
Bparticle*	particles_pick_var(Bimage* p, long average_kernel, long var_kernel,
				double nsig, double part_ori, double excl_dist, long bin);
Bparticle*	particles_pick_var(Bimage* p, long average_kernel, long var_kernel,
				double cutmin, double cutmax, double part_ori, double excl_dist, long bin);
Bparticle*	particles_pick_cc(Bstring& filename, long img_num, Bimage* ptemp, 
				Bimage* pmask, double hires, double lores, double fommin, double fommax, 
				double excl_dist, long bin=1);
Bparticle*	particles_pick_var(Bstring& filename, long img_num, 
				long average_kernel, long var_kernel, double nsig, 
				double part_ori, double excl_dist, long bin);
double		project_pick_particles(Bproject* project, Bimage* ptemp, Bimage* pmask, 
				double hires, double lores, double fommin, double fommax, double excl_dist, long bin);
double		project_pick_particles(Bproject* project, long average_kernel, long var_kernel,
				double nsig, double part_ori, double excl_dist, long bin);
long		project_pick_particles(Bproject* project, double din, double dout,
				int avg_kernel, double ainc, int flags, int contrast);
long		project_pick_background(Bproject* project, long number, 
				long average_kernel, long var_kernel, double excl_dist);
long		project_pick_sym_axis(Bproject* project, Bsymmetry& sym, int sym_axis, double axis_dist);
double		project_extract_orient_particles(Bproject* project, Bstring& tempfile, 
				Bsymmetry& sym, double hires, double lores, long bin);

