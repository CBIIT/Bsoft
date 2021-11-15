/**
@file	model_mask.h
@brief	Functions to generate masks from models
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20200329
**/

#include "rwimg.h"
#include "rwmodel.h"

// Function prototypes
Bimage*		model_create_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, Vector3<double> sam, double edge);
Bimage*		model_create_hull_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, Vector3<double> sam, int curv_flag, int fast);
Bimage*		model_create_shell_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, Vector3<double> sam, double shell_width, int curv_flag, int fast);
Bimage*		model_create_level_mask(Bmodel* model, Vector3<long> size,
				Vector3<double> origin, Vector3<double> sam);
Bimage*		img_extract_segments_using_model(Bimage* p, Bmodel* model, int multi_level);
int			img_add_model_to_mask(Bimage* p, Bmodel* model);
Bmodel*		model_from_multilevel_mask(Bimage* p);
Bimage*		model_create_projected_mask(Bmodel* model,
				Vector3<long> size, Vector3<double> ori,
				Vector3<double> sam, double dang, Bsymmetry& sym);

