/**
@file	mol_select.h
@brief	Library routines to select atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20211029
**/

#include "Vector3.h"
#include "rwmolecule.h"

// Function prototypes
long		molgroup_select(Bmolgroup* molgroup, Bstring selstr);
long		molgroup_atoms_selected(Bmolgroup* molgroup);
long		molgroup_select_all(Bmolgroup* molgroup);
long		molgroup_deselect_all(Bmolgroup* molgroup);
int			molgroup_select_chains(Bmolgroup* molgroup, Bstring chains);
int 		molgroup_coor_select_ring(Bmolgroup* molgroup, double rmin, double rmax);
int 		molgroup_coor_select(Bmolgroup* molgroup, Vector3<double> min, Vector3<double> max);
int 		molgroup_select(Bmolgroup* molgroup, long number);

