/**
@file	mol_transform.h
@brief	Library routines used for atomic coordinate transformations
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20070614
**/

#include "Matrix3.h"
#include "Transform.h"
#include "View.h"
#include "rwmolecule.h"

// Function prototypes
int 		molgroup_coor_shift(Bmolgroup* molgroup, Vector3<double> shift);
int 		mol_coor_shift(Bmolecule* mol, Vector3<double> shift);
int 		molgroup_coor_rotate(Bmolgroup* molgroup, Transform t);
int 		mol_coor_rotate(Bmolecule* mol, Transform t);
Bmolgroup*	molgroup_rotate_to_view(Bmolgroup* molgroup, View view, Vector3<double> origin, Vector3<double> trans);
Bmolgroup*	molgroup_rotate_from_view(Bmolgroup* molgroup, View view, Vector3<double> origin, Vector3<double> trans);
Bmolecule*	mol_rotate_to_view(Bmolecule* mol, View view, Vector3<double> origin, Vector3<double> trans);
Bmolecule*	mol_rotate_from_view(Bmolecule* mol, View view, Vector3<double> origin, Vector3<double> trans);
int 		molgroup_coor_transform(Bmolgroup* molgroup, Transform t);
int 		mol_coor_transform(Bmolecule* mol, Transform t);
int 		molgroup_coor_invert(Bmolgroup* molgroup, Vector3<double> point);
int 		molgroup_resolve_pbc(Bmolgroup* molgroup);
int  		molgroup_pack_in_periodic_box(Bmolgroup* molgroup);
int 		molgroup_coor_shift_PBC(Bmolgroup* molgroup, Vector3<double> shift);
int 		molgroup_coor_shift_rotate_PBC(Bmolgroup* molgroup, Vector3<double> origin, 
				Matrix3 mat, Vector3<double> shift);
int 		molgroup_shift_to_center_of_mass(Bmolgroup* molgroup);
int 		mol_shift_to_center_of_mass(Bmolecule* mol);
int			molgroup_place_at_coordinates(Bmolgroup* molgroup, Vector3<double> location);

