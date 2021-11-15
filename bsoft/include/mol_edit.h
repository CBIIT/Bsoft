/**
@file	mol_edit.h
@brief	Library routines used for atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20211029
**/

#include "rwmolecule.h"

// Function prototypes
long		molgroup_set_atom_types_to_elements(Bmolgroup* molgroup);
int			molgroup_remove_hydrogens(Bmolgroup* molgroup);
int			molgroup_add_disulfides(Bmolgroup* molgroup, double distance);
Bmolgroup**	molgroup_split_into_slices(Bmolgroup* molgroup, double slice_thickness, int& nslices);
int			molgroup_insert(Bmolgroup* molgroup, Bmolgroup* molinsert, double distance);
int			molgroup_randomize(Bmolgroup* molgroup, double random_max);
int			molgroup_randomize_B(Bmolgroup* molgroup, double B);
int			molgroup_random_displace_number(Bmolgroup* molgroup, long number, double stdev);
int			molgroup_remove_overlapping_atoms(Bmolgroup* molgroup, double mindist);
int			molgroup_bond_pseudo_atoms(Bmolgroup* molgroup, int atoms_per_bond, int wrap);
int			molgroup_prune_molecules(Bmolgroup* molgroup);
long		molgroup_delete_deselected_molecules(Bmolgroup* molgroup);
long		molgroup_prune_overlapping_atoms(Bmolgroup* molgroup, double mindist);
long		molgroup_delete_deselected_atoms(Bmolgroup* molgroup);
int			molgroup_untangle_molecules(Bmolgroup* molgroup, double sampling, double lambda);
int			molgroup_untangle_groups(Bmolgroup* molgroup, double sampling, double lambda);

