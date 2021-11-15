/**
@file	mol_bonds.h
@brief	Header for molecular dynamics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20180228
**/

#include "rwmd.h"
#include "rwmolecule.h"
#include "mol_util.h"
#include "utilities.h"

// Function prototypes
Bbond*		md_generate_bond_list(Bmolgroup* molgroup, Bmd* md);
Bbond*		md_generate_molecular_bond_list(Bmolgroup* molgroup, Bmd* md);
Bbond*		md_generate_bond_list_with_valence(Bmolgroup* molgroup, Bmd* md, int valence);
double		md_bond_list_set_parameters(Bbond* bondlist, Bbondtype* bondtype);
int			md_show_bonds(Bmolgroup* molgroup);
//int			md_show_angles(Bmolgroup* molgroup);
int			md_show_bond_stats(Bmolgroup* molgroup);
//Bangle*		md_generate_angle_list(Bmolgroup* molgroup, Bmd* md);
//int			md_angle_list_set_parameters(Bangle* anglelist, Bangletype* angletype);
double		md_find_bond_length(Batom* atom1, Batom* atom2, Bbondtype* bondtype);
Bbondtype*	md_find_bond_type(Batom* atom1, Batom* atom2, Bbondtype* bondtype);
int			md_show_bond_types(Bmolgroup* molgroup, Bbondtype* bondtype);
double		md_angle(Batom* atom1, Batom* atom2, Batom* atom3);
//Bangletype*	md_find_angle_type(Batom* atom1, Batom* atom2, Batom* atom3, Bangletype* angletype);
//double		md_find_angle(Batom* atom1, Batom* atom2, Batom* atom3, Bangletype* angletype);
double		md_calculate_deviations(Bmolgroup* molgroup, int wrap);
int			md_calculate_radial_deviation(Bmolgroup* molgroup);


