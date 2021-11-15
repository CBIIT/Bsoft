/**
@file	scatter.cpp
@brief	Functions for calculating electron scattering profiles
@author Bernard Heymann
@date	Created: 20190521
@date	Modified: 20210311
**/

#include "scatter.h"
#include "rwmodel_param.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
	The intensity is composed of the direct and scattering components.
	The average for the scattering component (mu) is modified
	with aperture-based fraction k for the elastic part
	and the use of an energy filter for the inelastic part.
*/
double			signal_integrated(double emu, double imu,
		double ke, double ki, int ef)
{
	double			sig(1), f(1), mu(ke*emu), pmu(1);
	if ( !ef ) mu += ki*imu;
	
	for ( long n=1; n<=10; ++n ) {
		f *= n;			// Factorial
		pmu *= mu;		// Powers of average
		sig += pmu/f;
//		cout << n << tab << f << tab << pmu << tab << sig << endl;
	}

	sig *= exp(-emu-imu);
	
	return sig;
}

double			signal_integrated_new(double emu, double imu, double k, int ef)
{
	double			sig(0), f(1), pmu(1);
	
	// Elastic scattering
	for ( long n=1; n<=10; ++n ) {
		f *= n;			// Factorial
		pmu *= emu;		// Powers of average
		sig += pmu/f;
	}

	if ( ef )
		sig = (1 + k*sig)*exp(-emu-imu);
	else
		sig = 1 - (1 - k)*sig*exp(-emu-imu);
	
	return sig;
}

auto string_is_empty = []( const std::string &s )
{
	return s.size() < 1;
};

vector<string>	all_elements(map<string,Bcomptype>& types)
{
	long			mx(0);
	
	for ( auto ct: types )
		if ( mx < ct.second.index() ) mx = ct.second.index();
	
//	cout << "Maximum Z = " << mx << endl;

	vector<string>	el(mx);
	
	for ( auto ct: types )
		el[ct.second.index()-1] = ct.second.identifier();

	el.erase( remove_if( el.begin(), el.end(), string_is_empty ), el.end() );

	return el;
}

JSvalue			js_all_elements(map<string,Bcomptype>& types)
{
	JSvalue			el(JSobject);

	for ( auto ct: types )
		el[ct.first] = 1.0;
	
	return el;
}


JSvalue			js_composition(vector<string>& vs)
{
	JSvalue			el(JSobject);
	
	for ( auto es: vs )
		el[es] = 1.0;
		
	return el;
}

JSvalue			js_composition(string s)
{
	vector<string>	vs(split(s,','));
	return js_composition(vs);
}

/**
@brief 	Calculates an elastic electron scattering profile for an component type.
@param 	&ct				component type.
@param 	ds 				frequency space sampling increment.
@param 	scut 			maximum spacial frequency.
@return vector<double>		array with scattering profile.
**/
vector<double>		calculate_scattering_curve(Bcomptype& ct,
				double ds, double scut)
{
	if ( ds < 1e-10 ) ds = 0.01;
	
	long			i, j, k, imax(scut/ds + 0.5);
	double			s2, ds2(ds*ds/4);
    
    // Calculate the atomic scattering profiles
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG calculate_scattering_curve: z = " << ct.index() << endl;
		
	vector<double>	scat;
	vector<double>	coef(ct.coefficients());
	
	for ( i=0; i<imax; ++i ) {
		s2 = i*i*ds2;
		scat.push_back(coef[10]);
		for ( j=0, k=5; j<5; ++j, ++k )
			scat[i] += coef[j]*exp(-coef[k]*s2);
	}
		
	return scat;
}

/**
@brief 	Calculates elastic electron scattering profiles for a subset of component types.
@param	&el				elements to calculate curves for.
@param 	&types			component type.
@param 	ds 				frequency space sampling increment.
@param 	scut 			maximum spacial frequency.
@return map<string, vector<double>>	array with scattering profile.
**/
map<string, vector<double>>		calculate_scattering_curves(JSvalue& el,
				map<string,Bcomptype>& types, double ds, double scut)
{
	if ( el.size() < 1 ) el = js_all_elements(types);
	
	map<string, vector<double>>	scat;
	
	for ( auto ct: types ) {
		if ( el.exists(ct.first) ) {
			if ( verbose & VERB_FULL )
				cout << "Calculating scattering curve for " << ct.first << endl;
			vector<double>	scat1 = calculate_scattering_curve(ct.second, ds, scut);
			scat[ct.first] = scat1;
		}
	}
		
	return scat;
}

map<string, vector<double>>		calculate_scattering_curves(
				map<string,Bcomptype>& types, double ds, double scut)
{
	map<string, vector<double>>	scat;
	
	for ( auto ct: types ) {
		if ( verbose & VERB_FULL )
			cout << "Calculating scattering curve for " << ct.first << endl;
		vector<double>	scat1 = calculate_scattering_curve(ct.second, ds, scut);
		scat[ct.first] = scat1;
	}
		
	return scat;
}

/**
@brief 	Calculates the integral of the whole curve.
@param 	&ct				component type.
@return double			integral.
**/
double		scatter_curve_integral(Bcomptype& ct)
{
	long			i, j, k, l;
	double			cs(0);
	vector<double>	coef(ct.coefficients());
	
	for ( i=0, k=5; i<5; ++i, ++k )
		for ( j=0, l=5; j<5; ++j, ++l )
			cs += coef[i]*coef[j]/(coef[k] + coef[l]);
	
//	cout << ct->el << tab << cs << endl;
	
	return cs;
}

/**
@brief 	Calculates the full elastic cross section (no aperture).
@param 	&ct				component type.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		elastic_cross_section(Bcomptype& ct, double volt)
{
	double			b2 = beta2(volt);
	double			hf = 1e10*PLANCK/(EMASS*LIGHTSPEED);
	double			cs = scatter_curve_integral(ct);

	cs *= 4*M_PI*hf*hf/b2;
	
	return cs;
}

/**
@brief 	Integrates the elastic scattering curve up to the aperture cutoff.
@param 	&ct				component type.
@param	ctf				CTF and microscope parameters.
@return double			integral.
**/
double		elastic_cross_section_integrated(Bcomptype& ct, CTFparam& ctf)
{
	long			i(0);
	double			ds(0.01), f, cs(0);
	double			b2 = beta2(ctf.volt());
	double			hf = 1e10*PLANCK/(EMASS*LIGHTSPEED);
	double			scut(ctf.frequency_cutoff());

	vector<double>	scat = calculate_scattering_curve(ct, ds, scut);

	for ( auto it = scat.begin(); it != scat.end(); ++it, ++i ) {
		f = *it;
		cs += f*f*i;
	}
	
	cs *= 2*M_PI*hf*hf*ds*ds/b2;
	
	return cs;
}

/**
@brief 	Calculates the combined elastic cross section for a defined composition.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		elastic_cross_section(map<string,Bcomptype>& types, double volt)
{
	double			cs(0);

	for ( auto ct: types )
		cs += elastic_cross_section(ct.second, volt) * ct.second.component_count();

	return cs;
}


/**
@brief 	Calculates the combined elastic cross section for a defined composition.
@param 	&types			component types.
@param	ctf				CTF and microscope parameters.
@return double			integral.
**/
double		elastic_cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double			cs(0);

	for ( auto ct: types )
		cs += elastic_cross_section_integrated(ct.second, ctf) * ct.second.component_count();

	return cs;
}

/**
@brief 	Calculates the inelastic cross section for an component type.
@param 	Z				atomic number.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		inelastic_cross_section_langmore(long Z, double volt)
{
	double			b2 = beta2(volt);

	double			cs = 1.5e-4*
		(sqrt(Z)/b2)*log(2*b2*
		(ECHARGE*volt+EMASS*LIGHTSPEED*LIGHTSPEED)/(20*ECHARGE));
	
	if ( Z == 1 ) cs *= 0.154;	// Correction for hydrogen
	
	return cs;
}

/**
@brief 	Calculates the combined inelastic cross section for a defined composition.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, double volt)
{
	double			cs(0);

	for ( auto ct: types )
		cs += inelastic_cross_section_langmore(ct.second.index(), volt) * ct.second.component_count();
	
	return cs;
}

/**
@brief 	Calculates the combined total cross section for a defined composition.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			integral.

**/
double		cross_section(map<string,Bcomptype>& types, double volt)
{
	double			ecs = elastic_cross_section(types, volt);
	double			ics = inelastic_cross_section_langmore(types, volt);

	if ( verbose & VERB_PROCESS ) {
		cout << "Elastic cross section:          " << ecs << " A2" << endl;
		cout << "Inelastic cross section:        " << ics << " A2" << endl;
	}
	
	return ecs + ics;
}

/**
@brief 	Calculates the combined total cross section for a defined composition.
@param 	&types			component types.
@param	&ctf			CTF and microscope parameters.
@return double			integral.

	If the slit width is specified, the energy filter is assumed to be used
	and only the elastic cross section is returned.
**/
double		cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double			ecs = elastic_cross_section_integrated(types, ctf);
	double			ics = inelastic_cross_section_langmore(types, ctf.volt());

	if ( verbose ) {
		cout << "Elastic cross section:          " << ecs << " A2" << endl;
		cout << "Inelastic cross section:        " << ics << " A2" << endl;
	}
	
	if ( ctf.slit_width() > 1 ) return ecs;
	else return ecs + ics;
}

/**
@brief 	Calculates the half-maximal frequency for an component type.
@param 	&ct				component type.
@return double			half-maximal frequency.
**/
double		cross_section_half_maximal_frequency(Bcomptype& ct)
{
	long			i(0);
	double			ds(0.001), f, csi(0);
	double			cs = scatter_curve_integral(ct);
	vector<double>	scat = calculate_scattering_curve(ct, ds, 10);
	
	for ( auto it = scat.begin(); it != scat.end() && csi < 0.5*cs; ++it, ++i ) {
		f = *it * ds;
		csi += 0.5*f*f*i;
	}
	
	return ds*i;
}

/**
@brief 	Calculates the half-maximal frequency for an component type.
@param 	&types			component types.
@return double			half-maximal frequency.
**/
double		cross_section_half_maximal_frequency(map<string,Bcomptype>& types)
{
	double			c, n(0), k(0);

	for ( auto ct: types ) {
		c = ct.second.component_count();
		k += c * cross_section_half_maximal_frequency(ct.second);
//		cout << ct.second.identifier() << tab << ct.second.component_count() << tab << k << endl;
		n += c;
	}

	if ( n ) k /= n;
	else k = 1;
	
	return k;
}

/**
@brief 	Creates a newAssembles a component composition from a set of materials with fractional contributions.
@param 	&mlist				list of materials.
@param 	&fractions			fractional contributions.
@return Bmaterial				new combined material.
**/
Bmaterial	material_combine(vector<Bmaterial>& mlist, vector<double> fractions)
{
	long					i(0);
	double					d(0);
	Bmaterial				material;
	material.identifier("Combined");
	map<string,Bcomptype>&	nucomp = material.composition();
	map<string,Bcomptype>	types;
	
	for ( auto m: mlist ) {
		d += fractions[i] * m.dalton_per_angstrom3();
		map<string,Bcomptype>	comp = m.composition();
		for ( auto ct: comp ) {
			if ( nucomp.find(ct.first) == nucomp.end() ) {
				nucomp[ct.first] = ct.second;
				nucomp[ct.first].component_count(0);
			}
			nucomp[ct.first].component_count_add(fractions[i] * ct.second.component_count());
		}
		i++;
	}
	
	material.density(d, 1);
	
	return material;
}

/**
@brief 	Returns vitreous ice density.
@param	types			reference parameters.
@return double			density in molecules/A3.
**/
Bmaterial	material_ice(map<string,Bcomptype>& types)
{
	Bmaterial		ice;
	map<string,Bcomptype>&	comp = ice.composition();
	Bcomptype		ct;

	ice.identifier("Vitreous ice");
	ice.density(0.93,0);		// Density in g/cm3
	ice.unit(1);				// Density in Da/A3
	
	ct.identifier() = "H";
	ct.index(1);
	ct.component_count(2);
	ct.mass(1);
	comp["H"] = ct;
	
	ct.identifier() = "O";
	ct.index(8);
	ct.component_count(1);
	ct.mass(16);
	comp["O"] = ct;
	
	ice.update_parameters(types);
	
//	double			wdens(mdens/18);				// Density in waters/A3
/*
	if ( verbose ) {
		cout << "Ice:" << endl;
		cout << "Density:                        " << dens << " g/cm3" << endl;
		cout << "Density:                        " << mdens << " Da/A3" << endl;
		cout << "Density:                        " << wdens << " molecules/A3" << endl;
	}
*/
	return ice;
}

/**
@brief 	Returns gold density.
@return double			density in molecules/A3.
**/
double		gold_density()
{
	double			dens(19.32);					// Density in g/cm3
	double			mdens(dens*AVOGADRO/1.0e24);	// Density in Da/A3
	double			wdens(mdens/196.966);			// Density in atoms/A3

	if ( verbose ) {
		cout << "Gold:" << endl;
		cout << "Density:                        " << dens << " g/cm3" << endl;
		cout << "Density:                        " << mdens << " Da/A3" << endl;
		cout << "Density:                        " << wdens << " molecules/A3" << endl;
	}
	
	return wdens;
}

/**
@brief 	Calculates the total cross section for ice up to the aperture cutoff.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			cross section.
**/
double		ice_cross_section(map<string,Bcomptype>& types, double volt)
{
	Bmaterial		ice = material_ice(types);

	double			cs = cross_section(ice.composition(), volt);
		
	return cs;
}


/**
@brief 	Calculates the average effective mean free path.
@param 	&material		material with types.
@param	&ctf			microscope parameters.
@return double			EMFP average.
**/
double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf)
{
	double			thick_step(100), thick_max(5000), t(thick_step), emfp_avg(0);

	vector<double>	I = signal_intensity(material, thick_step, thick_max, ctf);

	for ( long i=0; i<I.size(); ++i, t+=thick_step )
		emfp_avg += -t/log(I[i]);
		
	emfp_avg /= I.size();

	return emfp_avg;
}


/*
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&material		material with types.
@param	&ctf			microscope parameters.
@return double			effective mean free path.
**/
/*double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf)
{
	return effective_mean_free_path(material.composition(), ctf);
	map<string,Bcomptype>& comp = material.composition();
	
	double			d(material.number_per_angstrom3());
	bool			ef(ctf.slit_width()>0);	// Energy filter
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
//	double			ki = sc/(sc+0.054);
	double			ki = sc/(sc+0.2);
//	double			ki = sc*sc/(sc*sc+0.45*0.45);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	double			ics = inelastic_cross_section_langmore(comp, ctf.volt());
//	cout << d << tab << ecs << tab << ics << tab << ke << tab << ki << endl;
//	double			emfp = (ef)? 1/(d*ke*ecs): 1/(d*(ke*ecs+ki*ics));
	double			emfp = 1/(d*ke*ecs);
	if ( !ef ) emfp += 1/(d*ki*ics);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Material:                       " << material.identifier() << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
	}
	
	return emfp;
}*/

/**
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&materials		list of materials with types.
@param 	&fractions		fractional contributions.
@param	&ctf			microscope parameters.
@return double			effective mean free path.
**/
double		effective_mean_free_path(vector<Bmaterial>& materials, vector<double> fractions, CTFparam& ctf)
{
	Bmaterial		m = material_combine(materials, fractions);
	
	return effective_mean_free_path(m, ctf);
}

/**
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	thickness		component types.
@param 	&types			component types.
@param	ctf				microscope parameters.
@return double			intensity.
**/
double		ice_intensity(map<string,Bcomptype>& types, double thickness, CTFparam& ctf)
{
	if ( thickness < 1 ) thickness = 1000;

	Bmaterial		ice = material_ice(types);
	
	return signal_intensity(ice, thickness, ctf);
}

/**
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&material		material with types.
@param 	thickness		specimen thickness.
@param	ctf				microscope parameters.
@return double			intensity.
**/
double		signal_intensity(Bmaterial& material, double thickness, CTFparam& ctf)
{
	if ( thickness < 1 ) thickness = 1000;
	
	material.show();

	map<string,Bcomptype>& comp = material.composition();
	
	double			d(material.number_per_angstrom3());
	
	int				ef(ctf.slit_width()>0);	// Energy filter
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
	double			ki = sc/(sc+0.054);
//	double			ki = sc/(sc+0.2);
//	double			ki = sc*sc/(sc*sc+0.45*0.45);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	double			ics = inelastic_cross_section_langmore(comp, ctf.volt());
	double			mfp = 1/(d*(ecs+ics));
	
	double			emu = ecs*d*thickness;
	double			imu = ics*d*thickness;
	double			i = signal_integrated(emu, imu, ke, ki, ef);
//	cout << i << endl;
	double			emfp = -thickness/log(i);
//	double			emfp = 1/(d*(ke*ecs+ki*ics));
	
	if ( verbose ) {
		cout << "Mean free path:                 " << mfp << " A" << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
		cout << "Half_maximal frequency:         " << sm << " 1/Å" << endl;
		cout << "Frequency cutoff:               " << sc << " 1/Å" << endl;
		cout << "Aperture effect:                " << ke << tab << ki << endl;
		cout << "Thickness:                      " << thickness << " A" << endl;
		cout << "Intensity:                      " << i << endl;
		cout << endl;
	}
	
	return i;
}

/**
@brief 	Calculates the expected intensity vs thickness of vitreous ice.
@param 	&material		material with types.
@param 	thick_step		specimen tickness step size (angstrom).
@param 	thick_max		maximum specimen thickness (angstrom).
@param	&ctf			microscope parameters.
@return vector<double>	array of intensities.
**/
vector<double>	signal_intensity(Bmaterial& material, double thick_step, double thick_max, CTFparam& ctf)
{
	if ( thick_step < 1 ) thick_step = 100;
	if ( thick_max < 1 ) thick_max = 5000;

	map<string,Bcomptype>& comp = material.composition();

	int				ef(ctf.slit_width()>0);	// Energy filter
	double			d(material.number_per_angstrom3());
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
	double			ki = sc/(sc+0.054);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	
	double			ics = inelastic_cross_section_langmore(comp, ctf.volt());
	
	double			t, emu, imu, i, emfp;
	vector<double>	I;

//	cout << "Aperture effect:                " << k << endl;

	if ( verbose & VERB_FULL )
		cout << "t(A)\tI\tEMFP" << endl;
	for ( t=thick_step; t<=thick_max; t+=thick_step ) {
		emu = ecs*d*t;
		imu = ics*d*t;
		i = signal_integrated(emu, imu, ke, ki, ef);
		emfp = -t/log(i);
		I.push_back(i);
		if ( verbose & VERB_FULL )
			cout << t << tab << i << tab << emfp << endl;
	}
	
	return I;
}

/**
@brief 	Calculates and writes atomic scattering profiles to a file.
@param 	&paramfile	file with scattering coefficients.
@param 	&outfile	file to write curves to.
@param 	&selection	element selection.
@return int			0.

	The scattering coefficients are obtained from an input parameter file.

**/
int			write_scattering_curves(Bstring& paramfile, Bstring& outfile, Bstring& selection)
{
	
	long			i, j, nscat(100);
	double			s2;
	double			scut(1.0), ds(scut/nscat);
    double			ds2(ds*ds/4);
    
 	map<string,Bcomptype> 	types = read_atom_properties(paramfile);

	vector<string>	vs;
    if ( selection.length() < 1 ) vs = all_elements(types);
    else vs = split(selection.str(), ',');

    JSvalue			el = js_composition(vs);
	
	map<string, vector<double>>	scat = calculate_scattering_curves(el, types, ds, scut);

	ofstream		fd(outfile.c_str());
	if ( fd.fail() ) return -1;

	if ( verbose ) {
		cout << "Writing scattering curves to:   " << outfile << endl;
		cout << "Elemental selection:            " << selection << endl << endl;
	}
	
	fd << "Scattering coefficients:" << endl;
	for ( auto ct: types )
		if ( el.exists(ct.first) )
			fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( j=0; j<5; j++ ) {
		fd << "a" << j+1;
		for ( auto ct: types )
			if ( el.exists(ct.first) )
				fd << tab << ct.second.coefficients()[j];
		fd << endl;
	}
	for ( j=0; j<5; j++ ) {
		fd << "b" << j+1;
		for ( auto ct: types )
			if ( el.exists(ct.first) )
				fd << tab << ct.second.coefficients()[j+5];
		fd << endl;
	}
	fd << "c";
	for ( auto ct: types )
		if ( el.exists(ct.first) )
			fd << tab << ct.second.coefficients()[10];
	fd << endl;
	
	fd << "Scattering curves:" << endl << "s";
	for ( auto ct: types )
		if ( el.exists(ct.first) )
			fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( i=0; i<nscat; i++ ) {
		s2 = i*i*ds2;
		fd << sqrt(s2);
		for ( auto ct: types )
			if ( el.exists(ct.first) )
				fd << tab << scat[ct.first][i];
		fd << endl;
	}
	fd << endl;
	
	fd.close();
	
	return 0;
}

