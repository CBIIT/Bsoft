/**
@file	model_multifit.h
@brief	Header file for searching for a template in a map and returning multiple hits in a model.
@author Bernard Heymann
@date	Created: 20021027
@date	Modified: 20151127
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "View.h"


/* Function prototypes */
Bmodel*		model_from_densities(Bimage* p, Bimage* ptemp, View* view,
				double alpha, double alpha_step, double hires, double lores, Bimage* pmask, double threshold);
Bmodel*		model_from_densities_for_view(Bimage* p, Bimage* ptemp, View view,
				double hires, double lores, Bimage* pmask, double threshold);
Bmodel*		model_from_peaks(Bimage* p, double threshold, int wrap);

