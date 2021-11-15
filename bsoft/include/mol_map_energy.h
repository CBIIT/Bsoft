/**
@file	mol_map_energy.h
@brief	Headers of functions to calculate estimates of the fitting of molecules to maps.
@author Bernard Heymann
@date	Created: 20041230
@date	Modified: 20071224
**/

#include "rwmolecule.h"
#include "rwimg.h"

// Function prototypes 
double		molgroup_map_energy(Bmolgroup* molgroup, Bimage* map, double Kmap);
double		mol_map_energy(Bmolecule* mol, Bimage* map, double Kmap);
double		molgroup_map_correlation(Bmolgroup* molgroup, Bimage* map);
double		mol_map_correlation(Bmolecule* mol, Bimage* map);
double		molgroup_bond_fit_map_energy(Bmolgroup* molgroup, Bimage* map, double Kmap, int steps);

