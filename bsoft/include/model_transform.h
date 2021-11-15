/**
@file	model_transform.h
@brief	Library routines used for model transformation
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20141027
**/

#include "rwmodel.h"
#include "Transform.h"
#include <vector>

// Function prototypes
long		model_center(Bmodel* model);
long		model_shift(Bmodel* model, Vector3<double> shift);
long		models_shift(Bmodel* model, Vector3<double> shift);
long		model_scale(Bmodel* model, Vector3<double> scale, Vector3<double> origin);
long		models_scale(Bmodel* model, Vector3<double> scale, Vector3<double> origin);
long		model_reflect(Bmodel* model, Vector3<double> normal, Vector3<double> origin);
long		models_reflect(Bmodel* model, Vector3<double> normal, Vector3<double> origin);
double		model_reflect_and_compare(Bmodel* model, Vector3<double> normal, Vector3<double> origin);
long		model_rotate(Bmodel* model, Matrix3 mat);
long		model_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin);
long		model_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin, Vector3<double> shift);
long		model_rotate(Bmodel* model, View2<float> view);
long		model_rotate(Bmodel* model, View2<float> view, Vector3<double> origin, Vector3<double> shift);
long		model_rotate(Bmodel* model, Transform t);
double		model_rotate_and_compare(Bmodel* model, Transform t);
long		model_adjust_for_binning(Bmodel* model, Vector3<long> bin);
long		model_align_to_guide(Bmodel* model, Bmodel* guide);
Transform	model_find_transform(Bmodel* model, Bmodel* refmod);

