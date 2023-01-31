/**
@file	model_util.h
@brief	Library routines used for model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20230120
**/

#include "rwmodel.h"
#include "rwmolecule.h"
#include "Bstring.h"

// Function prototypes
long		models_process(Bmodel* model, long (modfunc)(Bmodel*));
long		models_process(Bmodel* model, long i, long (modfunc)(Bmodel*, long));
long		models_process(Bmodel* model, double d, long (modfunc)(Bmodel*, double));
long		models_process(Bmodel* model, Bstring& str, long (modfunc)(Bmodel*, Bstring& str));
long		models_process(Bmodel* model, string& str, long (modfunc)(Bmodel*, string& str));
long		model_list(Bmodel* model);
long		model_list_comp(Bmodel* model);
long		model_replace_components(Bmodel* model, Bmodel* modref);
long		model_merge(Bmodel* model);
long		model_number_ids(Bmodel* model);
long		model_rename(Bmodel* model, char first_name);
long		model_rename_components(Bmodel* model);
double		model_mass(Bmodel* model);
long		model_mass_all(Bmodel* model);
Vector3<double>	model_center_of_mass(Bmodel* model);
Vector3<double>	models_center_of_coordinates(Bmodel* model);
Vector3<double>	model_geometric_median(Bmodel* model);
double		model_gyration_radius(Bmodel* model);
Vector3<double> 	model_principal_axes(Bmodel* model, Vector3<double>* eigenvec);
Vector3<double> 	model_principal_axes(Bmodel* model, Matrix& eigenvec);
long 		model_principal_axes(Bmodel* model);
long		model_radial_distribution(Bmodel* model, double interval);
Bmolgroup*	model_assemble(Bmodel* model, Bstring& paramfile, int separate);
Bmodel*		model_generate_com(Bmolgroup* molgroup);
long		model_update_centers_of_mass(Bmodel* model, Bmolgroup* molgroup);
long		model_average_components(Bmodel* model, int number);
Bcomponent**	component_get_array(Bmodel* model, long& ncomp);
Vector3<double>	component_plane(vector<Bcomponent*>& comparray, double& offset);
vector<Vector3<double>>	models_calculate_bounds(Bmodel* model);
vector<vector<Bcomponent*>>	model_component_grid(Bmodel* model, Vector3<long>& size,
				Vector3<double>& origin, Vector3<double>& sampling);
vector<Bcomponent*>	models_get_component_array(Bmodel* model);
vector<vector<Bcomponent*>>	model_split_into_slices(Bmodel* model, double bottom, double top, double thickness);

