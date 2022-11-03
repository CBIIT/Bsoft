/**
@file	scatter.h
@brief	Functions for calculating electron scattering profiles
@author Bernard Heymann
@date	Created: 20190521
@date	Modified: 20210515
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "ctf.h"
#include "string_util.h"
#include "utilities.h"

// Function prototypes 
vector<string>	all_elements(map<string,Bcomptype>& types);
vector<double>	calculate_scattering_curve(Bcomptype& ct,
				double ds, double scut);
map<string, vector<double>>		calculate_scattering_curves(JSvalue& el,
				map<string,Bcomptype>& types, double ds, double scut);
map<string, vector<double>>	calculate_scattering_curves(
				map<string,Bcomptype>& types, double ds, double scut);
double		scatter_curve_integral(Bcomptype& ct);
map<string, vector<double>>		calculate_potential_curves(JSvalue& el,
				map<string,Bcomptype>& types, double dr, double rmax);
double		elastic_cross_section(Bcomptype& ct, double volt);
double		elastic_cross_section_integrated(Bcomptype& ct, CTFparam& ctf);
double		elastic_cross_section(map<string,Bcomptype>& types, double volt);
double		elastic_cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf);
double		inelastic_cross_section_langmore(long Z, double volt);
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, double volt);
double		cross_section(map<string,Bcomptype>& types, double volt);
double		cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf);
double		cross_section_half_maximal_frequency(Bcomptype& ct);
double		cross_section_half_maximal_frequency(map<string,Bcomptype>& at);
Bmaterial	material_combine(vector<Bmaterial>& mlist, vector<double> fractions);
Bmaterial	material_ice(map<string,Bcomptype>& types);
//double		effective_mean_free_path(map<string,Bcomptype>& types, CTFparam& ctf);
double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf);
double		effective_mean_free_path(vector<Bmaterial>& material, vector<double> fractions, CTFparam& ctf);
double		signal_intensity(Bmaterial& material, double thickness, CTFparam& ctf);
vector<double>	signal_intensity(Bmaterial& material, double thick_step, double thick_max, CTFparam& ctf);
int			write_scattering_curves(Bstring& paramfile, Bstring& outfile, Bstring& selection, double resolution);
int			write_potential_curves(Bstring& paramfile, Bstring& outfile, Bstring& selection, double radius);
