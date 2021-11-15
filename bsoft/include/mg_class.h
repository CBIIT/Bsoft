/**
@file	mg_class.h
@brief	Header file for classification of raw single particle images with respect to multiple models
@author Bernard Heymann
@date	Created: 20010222
@date	Modified: 20151008
**/

#include "mg_processing.h"
#include "FSI_Kernel.h"

// Function prototypes
long	   mg_classify(Bproject* project, double resolution_hi, double resolution_lo,
					int fom_type, double fom_cut, FSI_Kernel* kernel, int ctf_apply, int img_out);

