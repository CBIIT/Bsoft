/**
@file	mol_compare.h
@brief	Library routines used to compare sets of atomic coordinates
@author Bernard Heymann
@date	Created: 20021020
@date	Modified: 20200917
**/

#include "Transform.h"
#include "Matrix.h"
#include "rwmolecule.h"

// Function prototypes
double		molgroup_rotate_and_compare(Bmolgroup* molgroup, Transform t);
Transform	molgroup_find_transformation(Bmolgroup* molgroup1, Bmolgroup* molgroup2);
Transform	mol_find_transformation(Bmolecule* mol1, Bmolecule* mol2, int offset);
double		molgroup_calculate_rmsd(Bmolgroup* molgroup, Bmolgroup* molgroup2);
double		mol_calculate_rmsd(Bmolecule* mol1, Bmolecule* mol2);
double		molgroup_calc_brute_rmsd(Bmolgroup* molgroup1, Bmolgroup* molgroup2);
Matrix		mol_distance_matrix(Bmolecule* m1, Bmolecule* m2);

