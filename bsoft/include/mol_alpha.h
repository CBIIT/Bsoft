/**
@file	mol_alpha.h
@brief	Header for functions to make and analyze alpha helices.
@author Bernard Heymann
@date	Created: 20050315
@date	Modified: 20060121
**/

#include "rwmolecule.h"

// Function prototypes
Bmolgroup*	molgroup_generate_alpha_helix(int length);
Bmolecule*	mol_generate_alpha_helix(int length);
int			molgroup_set_alpha_helix(Bmolgroup* molgroup, int helix_start, int helix_end);
int			molgroup_find_helical_axes(Bmolgroup* molgroup);
Vector3<double>	alpha_find_center(Bresidue* resfirst, Bresidue* reslast);
Vector3<double>	alpha_find_orientation(Bresidue* resfirst, Bresidue* reslast);
Vector3<double>	mol_find_alpha_orientation(Bmolecule* mol, int set_std);
Vector3<double>	point_on_helix_axis(Vector3<double> ca1, Vector3<double> ca2, Vector3<double> ca3);
Bmolgroup*	molgroup_consolidate_alpha(Bmolgroup* molgroup);

