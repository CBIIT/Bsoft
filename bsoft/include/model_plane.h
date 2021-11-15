/**
@file	model_plane.h
@brief	Library routines used for plane models.
@author Bernard Heymann
@date	Created: 20140925
@date	Modified: 20141008
**/

#include "rwmodel.h"

// Function prototypes
Vector3<double>	model_fit_plane(Bmodel* model);
Bmodel*		model_generate_from_plane_guide(Bmodel* guide, double separation, double sigma);

