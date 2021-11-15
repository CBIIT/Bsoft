/**
@file	mol_md.h
@brief	Header for molecular dynamics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20060412
**/

#include "rwmd.h"
#include "rwmolecule.h"
#include "mol_util.h"
#include "utilities.h"

// Function prototypes
double		md_leapfrog(Bmolgroup* molgroup, Bmd* md, int max_iter, double velocitylimit);
int			md_zero_forces(Bmolgroup* molgroup);
double		md_bond_forces(Bmolgroup* molgroup, double Kbond, int wrap);
double		md_angular_forces(Bmolgroup* molgroup, double Kangle, int wrap);
double		md_nonbonded_forces(Bmolgroup* molgroup, Bmd* md);
int			atom_nonbonded_forces(Batom* atom1, Batom* atom2, Bmd* md, Vector3<double> box);
double		md_point_force(Bmolgroup* molgroup, Vector3<double> point, double Kpoint, double decay);



