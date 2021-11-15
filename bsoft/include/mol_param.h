/**
@file	mol_param.h
@brief	Header for functions to extract parameters from coordinate files.
@author Bernard Heymann
@date	Created: 20050304
@date	Modified: 20090226
**/

#include "rwmolecule.h"
#include "rwmd.h"

// Function prototypes
Batomtype*	molgroup_get_atom_types(Bmolgroup* molgroup, int elements);
Bmd*		md_calculate_parameters(Bmolgroup* molgroup, int elements, int show_bond, int show_angle);

