/**
@file	mol_monte.h
@brief	Headers of functions using a monte carlo metroplis algorithm to energy minimize molecular positions.
@author Bernard Heymann
@date	Created: 20041230
@date	Modified: 20071223
**/

#include "rwmolecule.h"
#include "rwimg.h"
#include "rwmodel.h"
#include "rwmd.h"


// Function prototypes 
Bmolgroup*	monte_carlo_metropolis(Bmolgroup* molgroup, Bmd* md, Bimage* map, 
				double beta, double max_angle, double max_shift, 
				long max_iter, int rigid, 
				double (Efunc)(Bmolgroup*, Bimage*, Bmd*),
				int (Tfunc)(Bmolgroup*, double, double));
Bmolgroup*	molgroup_generate_masked_grid_list(Bmolgroup* molgroup, Bimage* pmask,
				Vector3<double> grid_sampling, Bstring filename);
Bmodel*		molgroup_generate_masked_grid_list(Bmolgroup* molgroup, 
				Vector3<double> grid_sampling, Bimage* pmask);
Bmolgroup*	molgroup_generate_orientation_list(Bmolgroup* molgroup,
				double angle_step, Bstring filename, int whole);
Bmodel*		molgroup_generate_orientation_list(Bmolgroup* molgroup, double angle_step);
Bmolgroup*	mcm_molecule_list(Bmolgroup* molgroup, Bmd* md, Bimage* map, double beta, 
				double max_angle, double max_shift, long max_iter, int rigid);
int			mcm_molecule_groups(Bmolgroup* molgroup, Bmodel* model, Bmd* md, Bimage* map, double beta, 
				double max_angle, double max_shift, long max_iter, int rigid);
int			mcm_molecule_list(Bmolgroup* molgroup, Bmodel* model, Bmd* md, Bimage* map, double beta, 
				double max_angle, double max_shift, long max_iter, int rigid);
int			molgroup_set_box_to_map_boundaries(Bmolgroup* molgroup, Bimage* map);
long		molgroup_test_if_within_box(Bmolgroup* molgroup, Vector3<double> min, Vector3<double> max);
double		monte_rigid_body_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md);
double		monte_atom_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md);
double		monte_bond_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md);
double		molgroup_atom_overlap(Bmolgroup* molgroup, Bmd* md);
int			molgroup_rigid_body_transform(Bmolgroup* molgroup, double max_angle, double shift_std);
int			molgroup_move_atoms_down_energy(Bmolgroup* molgroup, double max_angle, double max_shift);

