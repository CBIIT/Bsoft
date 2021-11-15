/**
@file	mol_util.h
@brief	Library routines used for atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20210512
**/

#include "Vector3.h"
#include "rwmolecule.h"
#include "json.h"

#ifndef _Latom_
#define _Latom_
/************************************************************************
@Object: struct Latom
@Description:
	A structure for use in a linked list of atoms.
@Features:
	This is used to group an arbitrary set of atoms.
*************************************************************************/
struct Latom {
	Latom*		next;		// Next atom 
	Batom*		atom;
} ;
#endif


// Function prototypes
int  		molgroup_Bfactors(Bmolgroup* molgroup);
int 		molgroup_print_sequence(Bmolgroup* molgroup);
int 		molgroup_set_radius(Bmolgroup* molgroup, double radius);
int			molgroup_rename(Bmolgroup* molgroup, char first_name);
int 		molgroup_residue_renumber(Bmolgroup* molgroup, int first);
int 		molgroup_atom_renumber(Bmolgroup* molgroup, int first);
int 		molgroup_coor_reset_occupancy(Bmolgroup* molgroup, int range_first, 
				int range_last, double occupancy);
double 		molgroup_weight_from_atoms(Bmolgroup* molgroup);
double 		molgroup_weight_from_sequence(Bmolgroup* molgroup);
Vector3<double> 	molgroup_center_of_mass(Bmolgroup* molgroup);
Vector3<double> 	molgroup_selected_center_of_mass(Bmolgroup* molgroup);
Vector3<double> 	molgroup_show_center_of_mass(Bmolgroup* molgroup);
Vector3<double> 	mol_center_of_mass(Bmolecule* mol);
Vector3<double> 	mol_show_center_of_mass(Bmolecule* mol);
Vector3<double> 	molgroup_principal_axes(Bmolgroup* molgroup, Vector3<double>* eigenvec);
Vector3<double> 	mol_principal_axes(Bmolecule* mol, Vector3<double>* eigenvec);
double		molgroup_sphericity(Bmolgroup* molgroup, double da);
double		molgroup_density(Bmolgroup* molgroup);
double		molgroup_volume(Bmolgroup* molgroup, Bstring& paramfile, int wrap);
int			molgroup_composition(Bmolgroup* molgroup, Bstring& paramfile);
JSvalue		molgroup_elements(Bmolgroup* molgroup);
JSvalue		molgroup_elements(Bmolgroup* molgroup, Bstring& paramfile);
JSvalue		protein_composition_default(double mass);
int			molgroup_radial_density(Bmolgroup* molgroup, double interval, double cutoff, int wrap);
Latom**		molgroup_atom_mesh_lists(Bmolgroup* molgroup, Vector3<int> size, Vector3<double> sampling);
vector<Batom*>	atom_get_array(Bmolgroup* molgroup);

