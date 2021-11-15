/**
@file	model_views.h
@brief	Library routines used for analysing model component views
@author Bernard Heymann
@date	Created: 20081120
@date	Modified: 20141029
**/

#include "rwmodel.h"
#include "Bstring.h"

// Function prototypes
//View*		views_from_model(Bmodel* model);
list<View2<float>>	views_from_model(Bmodel* model);
//View*		views_from_models(Bmodel* model);
list<View2<float>>	views_from_models(Bmodel* model);
long		model_set_views(Bmodel* model, View2<float> view);
long		model_invert_views(Bmodel* model);
long		model_find_views(Bmodel* model, Bstring& reffile, Bstring& paramfile);
long		model_calculate_views(Bmodel* model, Bstring& mode);
long		model_calculate_local_views(Bmodel* model);
long		model_view_directions(Bmodel* model, int bin_width, int ref_flag);
int			component_hand(Bstring s);

