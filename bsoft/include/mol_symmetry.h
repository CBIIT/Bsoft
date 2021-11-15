/**
@file	mol_symmetry.h
@brief	Library routines used for symmetry operations on atomic coordinates
@author Bernard Heymann
@date	Created: 20021020
@date	Modified: 20150813
**/

#include "rwmolecule.h"
#include "rwresprop.h"
#include "symmetry.h"

// Function prototypes
int 		molgroup_apply_point_group(Bmolgroup* molgroup, Bsymmetry& sym, View ref_view);
int			molgroup_generate_helix(Bmolgroup* molgroup, View ref_view, 
				double helix_rise, double helix_angle, int gen_down, int gen_up);
int			molgroup_apply_symmetry_from_pdb(Bmolgroup* molgroup, Bstring& filename);
int			molgroup_apply_matrices_from_pdb(Bmolgroup* molgroup, Bstring& filename);
int 		molgroup_find_standard_view(Bmolgroup* molgroup, Bsymmetry& sym, View ref_view);
int 		molgroup_orient_to_standard_view(Bmolgroup* molgroup, Bsymmetry& sym, 
				View ref_view, Bresidue_matrix* simat);
double		molgroup_symmetry_RMSD(Bmolgroup* molgroup, Bsymmetry& sym);
double		molgroup_symmetry_B(Bmolgroup* molgroup, Bsymmetry& sym);
int 		molgroup_generate_crystal(Bmolgroup* molgroup, UnitCell unit_cell, Vector3<int> number);


