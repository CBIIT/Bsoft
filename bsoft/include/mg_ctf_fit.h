/**
@file	mg_ctf_fit.h 
@brief	Header file for CTF (contrast transfer function) functions
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210813
**/

#include "rwimg.h"
#include "ctf.h"

// Function prototypes
Bimage*		img_ctf_radial_average(Bimage* p, long n, CTFparam& em_ctf);
double		img_ctf_fit(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires,
				double def_start=1e3, double def_end=2e5, double def_inc=1e3, int flag=0);
double		img_ctf_isotropy(Bimage* p, long n, double lores, double hires);
double		img_ctf_fit_baseline(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires);
double		img_ctf_fit_envelope(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires);
double		ctf_find_defocus(vector<double>& v, CTFparam& em_ctf,
				long rmin, long rmax, double step_size,
				double def_start, double def_end, double def_inc);
double		img_ctf_find_defocus(Bimage* p, long n, CTFparam& em_ctf,
				double lores, double hires, 
				double def_start=1e3, double def_end=2e5, double def_inc=1e3);
double		img_ctf_fit_astigmatism(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires);
double		img_water_ring_index(Bimage* p, long img_num, CTFparam& em_ctf);
double		img_water_ring_index(Bimage* prad);

