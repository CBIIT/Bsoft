/**
@file	mg_ctf_sim.h
@brief	Header file for functions to simulate tilted micrographs
@author Bernard Heymann
@date	Created: 20150224
@date	Modified: 20200525
**/

#include "Bimage.h"
#include "ctf.h"

// Function prototypes
Bimage*		img_ttf_simulate(Bimage* pn, CTFparam& ctf, int action,
				double wiener, double tilt, double tilt_inc, double axis);

