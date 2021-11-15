/**
@file	model_color.h
@brief	Functions to color models.
@author Bernard Heymann
@date	Created: 20080206
@date	Modified: 20210319
**/

#include "rwmodel.h"

// Function prototypes 
int			model_color_uniformly(Bmodel* model, RGBA<float> color);
int			model_color_selected(Bmodel* model, RGBA<float> color);
int			model_color_by_order(Bmodel* model);
int			model_color_by_density(Bmodel* model);
int			model_color_by_fom(Bmodel* model);
int			model_color_by_selection(Bmodel* model);
int			model_color_selected_types(Bmodel* model, RGBA<float> rgba);
int			model_color_curvature(Bmodel* model);
int			model_color_chiral_vertices(Bmodel* model);

