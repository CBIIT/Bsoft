/**
@file	mg_ctf_focal.h
@brief	Header file for CTF (contrast transfer function) functions
@author 	Bernard Heymann
@date	Created: 20000426
@date	Modified: 20230127
**/

#include "rwimg.h"
#include "ctf.h"

// Function prototypes
Bimage*		img_ctf_focal_series(CTFparam& cp, double def_min, double def_max, double def_inc,
				Vector3<long> size, Vector3<double> sam, double lores, double hires);
Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter);
Bimage*		img_fspace_extract_sphere(Bimage* p, double volt);

