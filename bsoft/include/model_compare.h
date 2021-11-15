/**
@file	model_compare.h
@brief	Functions to compare models and components
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20161011
**/

#include "rwmodel.h"
#include "Matrix.h"

// Function prototypes
long		model_component_number_difference(Bmodel* model1, Bmodel* model2);
long		model_maxnum_components(Bmodel* model);
double		model_compare(Bmodel* model1, Bmodel* model2);
Matrix		model_distance_matrix(Bmodel* model, int view_flag);
Matrix		model_distance_matrix(Bmodel* m1, Bmodel* m2);
Matrix		model_adjacency_matrix(Bmodel* model);
long		model_consolidate(Bmodel* model, double distance);
Bmodel*		models_consensus(Bmodel* model, double distance);

