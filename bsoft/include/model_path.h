/**
@file	model_path.h
@brief	Library routines used for model processing
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20150208
**/

#include "rwmodel.h"
#include "Matrix.h"
#include "Bstring.h"

// Function prototypes
Matrix		model_shortest_path(Bmodel* model);
double		model_wiener_index(Bmodel* model);
Bmodel*		model_hamiltonian_cycle(Bmodel* model);

