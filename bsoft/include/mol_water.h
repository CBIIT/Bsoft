/**
@file	mol_water.h
@brief	Generating and managing water
@author Bernard Heymann
@date	Created: 20001014
@date	Modified: 20060122
**/

#include "rwmolecule.h"
#include "Vector3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bmolecule*	mol_generate_one_water(Bmolecule** mollist, char* watername, Vector3<double> Ocoord);
Bmolgroup*	molgroup_generate_regular_water(Vector3<double> size, int type);
Bmolgroup*	molgroup_generate_random_water(Vector3<double> size);
Bbond*		water_bond_list(Bmolgroup* molgroup);
Bangle*		water_angle_list(Bmolgroup* molgroup);
int			molgroup_calc_water_rdf(Bmolgroup* molgroup, double interval, double cutoff);

