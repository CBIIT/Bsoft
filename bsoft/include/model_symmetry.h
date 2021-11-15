/**
@file	model_symmetry.h
@brief	Library routines used for model symmetry operations
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20191101
**/

#include "rwmodel.h"
#include "View.h"
#include "Bstring.h"

// Function prototypes
long		model_find_asymmetric_unit(Bmodel* model, string& symmetry_string);
long 		model_apply_point_group(Bmodel* model, string& symmetry_string,
					Vector3<double> origin, View ref_view, int flags=0);
long		models_apply_point_group(Bmodel* model, string& symmetry_string,
					Vector3<double> origin, View ref_view, int flags=0);
long 		model_symmetrize(Bmodel* model, string& symmetry_string);
long		model_symmetry_related(Bmodel* model, string& symmetry_string);

