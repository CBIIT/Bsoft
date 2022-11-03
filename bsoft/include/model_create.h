/**
@file	model_create.h
@brief	Library routines used for model creation
@author Bernard Heymann
@date	Created: 20090714
@date	Modified: 20220210
**/

#include "rwmodel.h"
#include "symmetry.h"

// Function prototypes
Bcomponent*	model_add_component(Bmodel* model, Bstring &id, Bstring &type, Vector3<double> loc);
Bcomponent*	model_add_components(Bmodel* model, Bstring &id, Bstring &type, vector<Vector3<double>> loc);
Bmodel*		model_platonic(Bsymmetry& sym, double radius);
Bmodel*		model_helix(double radius, double helix_rise, double helix_angle, long ncomp);
Bmodel*		model_random(long ncomp, double comp_radius, double max_radius);
Bmodel*		model_random(long ncomp, double comp_radius, Vector3<double> min, Vector3<double> max);
Bmodel*		model_random_gaussian(long ncomp, double std);
Bmodel*		model_random_shell(long ncomp, double radius);
Bmodel*		model_random_shell(long ncomp, double radius, double separation);
Bmodel*		model_create_shell(long number, double radius, double distance);
Bmodel*		model_create_fibonacci_sphere(long number, double radius);
Bmodel*		model_create_dodecahedron(double radius, long divisions, double sphere_fraction);
Bmodel*		model_create_icosahedron(double radius, long divisions, double sphere_fraction);
Bmodel*		model_create_circle(double radius, double z, double distance);
Bmodel*		model_create_ellipse(Vector3<double> axes, double distance);
Bmodel*		model_create_ellipsoid(Vector3<double> axes, double distance);
Bmodel*		model_create_cylinder(Vector3<double> direction, double radius,
				double length, double distance);
Bmodel*		model_create_plane(double length, double width, double z, double distance);
Bmodel*		model_create_spindle(Vector3<double> direction, double radius,
				double separation, double packing);
Bmodel*		model_create_cubic_lattice(Vector3<long> lattice, double separation);
Bmodel*		model_create_hexagonal_lattice(Vector3<long> lattice, double separation);

