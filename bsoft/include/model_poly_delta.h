/**
@file	model_poly_delta.h
@brief	Functions to generate deltagraph polyhedra.
@author Bernard Heymann
@date	Created: 20080103
@date	Modified: 20111229
**/

#include "rwmodel.h"

// Function prototypes 
//Bmodel*		model_delta_create_tube(int radius, int height);
Bmodel*		model_delta_create_tube(int h, int k, int height);
Bmodel*		model_delta_create_cylinder(int type, int radius, int height);

