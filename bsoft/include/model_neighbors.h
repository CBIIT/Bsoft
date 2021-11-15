/**
@file	model_neighbors.h
@brief	Functions to manipulate model component neighbors
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20141006
**/

#include "rwmodel.h"

// Function prototypes 
int			model_set_neighbors(Bmodel* model, int number);
int			model_set_neighbors(Bmodel* model, int number, double distance);
int			model_neighbor_reciprocity(Bmodel* model);

