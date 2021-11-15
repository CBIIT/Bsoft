/**
@file	spline.h 
@brief	Header file for spline interpolation
@author Bernard Heymann
@date	Created: 20020808
@date	Modified: 20151023
**/

#include "Vector3.h"

// Function prototypes
Vector3<double>*	vector3_catmull_rom_spline(long ncoord, Vector3<double>* coords, long& nspline);


