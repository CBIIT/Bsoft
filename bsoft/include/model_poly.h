/**
@file	model_poly.h
@brief	Functions to manipulate polyhedral coordinate files
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20100329
**/

#include "rwmodel.h"

// Function prototypes 
int			model_poly_faces(Bmodel* model);
int			model_poly_generate(Bmodel* model);
Bmodel*		model_poly_dual(Bmodel* model, int order);
int			model_vertex_types(Bmodel* model);
int			model_extended_vertex_types(Bmodel* model);
int			model_poly_analyze(Bmodel* model);
int			model_poly_links(Bmodel* model);
int			model_poly_angles(Bmodel* model);
double		model_poly_regularity(Bmodel* model);
double		model_poly_planarity(Bmodel* model);
double		model_poly_energy(Bmodel* model, double angle_ref);
int			model_poly_pentagon_adjacency(Bmodel* model);
Bstring		model_poly_find_symmetry(Bmodel* model, double find);
int			model_poly_hand(Bmodel* model);
int			model_poly_compare(Bmodel* model, Bmodel* refmodel);
vector<double>	model_poly_eigenvalues(Bmodel* model, int show);
vector<double>	model_poly_sphere_coor(Bmodel* model);


