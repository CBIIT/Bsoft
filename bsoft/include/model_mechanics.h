/**
@file	model_mechanics.h
@brief	Functions to do molecular mechanics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20210324
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "rwimg.h"

// Function prototypes 
double		model_mechanics(Bmodel* model, Bmodparam& md, int mm_type, int max_iter,
				double max_shift, double velocitylimit);
int			model_minimize(Bmodel* model, double max_shift);
double		model_verlet(Bmodel* model, double timestep, double Kfriction, double velocitylimit);
double		model_electrostatic_energy(Bmodel* model, Bmodparam& md);
double		model_distance_energy(Bmodel* model, Bmodparam& md);
double		model_grid_distance_energy(Bmodel* model, Bmodparam& md);
double		model_neighbor_distance_energy(Bmodel* model, Bmodparam& md);
double		model_soft_sphere_energy(Bmodel* model, double Kdistance, double dist_ref);
double		model_lennard_jones_energy(Bmodel* model, double Kdistance, double distance);
double		model_morse_energy(Bmodel* model, double Kdistance, double distance);
double		model_link_energy(Bmodel* model, Bmodparam& md);
double		model_link_energy(Bmodel* model, double Klink);
double		model_angle_energy(Bmodel* model, Bmodparam& md);
double		model_polygon_angle_energy(Bmodel* model, double Kangle);
double		model_polygon_energy(Bmodel* model, double Kpoly);
double		model_polygon_plane_energy(Bmodel* model, double Kpolyplane);
double		model_neighbor_plane_energy(Bmodel* model, double Kplane);
double		model_point_force(Bmodel* model, Vector3<double> point, double Kpoint, double decay);
double		model_radial_energy(Bmodel* model, Vector3<double> point, double d0, double Kd);
double		model_guide_energy(Bmodel* model, Bmodparam& md);
double		model_polyhedron_guide_energy(Bmodel* model, Bmodel* gmod, double Kguide);
double		model_map_energy(Bmodel* model, Bimage* map, double Kmap);
double		model_map_energy(Bmodel* model, Bimage* map, double Kmap, double sigma);
int			model_zero_forces(Bmodel* model);
int			model_calculate_deviations(Bmodel* model);
int			model_calculate_deviations(Bmodel* model, Bmodparam& md);
int			model_regularize(Bmodel* model, int reg_iter, double distance, 
				double Kdistance, double Klink, double Kangle, double Kpoly, 
				double Kpolyplane, double Kpoint, double decay);

