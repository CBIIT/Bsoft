/**
@file	bscatter.cpp
@brief	Program to calculate scattering cross sections
@author	Bernard Heymann
@date	Created: 20190521
@date	Modified: 20210311
**/

#include "molecule_to_map.h"
#include "mol_util.h"
#include "seq_util.h"
#include "rwmodel_param.h"
#include "rwresprop.h"
#include "ctf.h"
#include "scatter.h"
#include "json.h"
//#include "linked_list.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Constants
#define MAXSCAT 	1000			// Maximum number of atomic scattering data points

Bmaterial	protein_material_default();
Bmaterial	material_from_molgroup(Bmolgroup* molgroup);
//double		material_update_parameters(Bmaterial& material, map<string,Bcomptype>& atompar);
int			cross_sections(vector<string>& elements, map<string,Bcomptype>& atompar, CTFparam& ctf);
int			cross_sections(Bmaterial& material, CTFparam& ctf);
int			cross_section_half_maximal_frequencies(vector<string>& elements, map<string,Bcomptype>& atompar);
int			cross_section_half_maximal_frequencies(Bmaterial& material);
int			aperture_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> apser);
int			collection_angle_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> angser);
double		particle_snr(Bmaterial& material, double mass, double radius, double thickness, CTFparam& ctf);
Bplot*		particle_spectral_signal(Bmaterial& material, double mass,
				double radius, double thickness, CTFparam& ctf);

// Usage assistance
const char* use[] = {
" ",
"Usage: bscatter [options] input.pdb",
"-----------------------------------",
"Calculating scattering cross sections.",
"Elements can be specified with the -select option,",
"from a composition file with the -composition option,",
"or from an input sequence or coordinate file.",
"Otherwise a general elemental composition for proteins is used.",
" ",
"Actions:",
"-halfmaximum             Calculate the half-maximal frequency for each element",
"-series 30,50,70         Calculate EMFP using series of aperture sizes (um).",
"-angles 2.6,4.6,11.8     Calculate EMFP using series of collection angles (mrad).",
" ",
"Material parameters:",
"-material \"Vitreous ice\" Material name from properties file,",
"                         or coordinate file.",
"-elements H,C,O          Material elements.",
"-counts 10,3,2           Element counts.",
"-density 1.27            Material density (g/cm3).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-Volt 100                Acceleration voltage (default 300 kV).",
"-aperture 100            Aperture size (default 100 µm).",
"-focallength 2.5         Focal length (default 3.5 mm).",
"-slitwidth 25            Energy filter slit width (default 0 eV - no filter).",
"-thickness 1500          Specimen thickness (angstrom).",
"-mass 234m               Molecular weight (k=kDa, m=MDa).",
"-radius 284              Particle radius (angstrom).",
"-dose 25                 Electron dose (electrons/pixel).",
"-defocus 1.2,3.5,0.1     Defocus range an steps (µm).",
" ",
"Input:",
"-atoms atomprop.star     Input atom properties file.",
"-residues resprop.star   Input residue properties file.",
"-properties matter.star  Input material density and composition file.",
" ",
NULL
};


int		main(int argc, char** argv)
{
	// Initialize variables
	bool			halfmax(0);			// Flag to calculate half-maximal frequency
	double			volts(300000);	 	// In volts
	double			aperture(1e6);		// Aperture size in angstrom
	vector<double>	apser;				// Series of apertures
	vector<double>	angser;				// Series of collection angels
	Bstring			material_str;		// Material string (name or coordinate file)
	Bstring			elements;			// Elements
	Bstring			count_str;			// String with element counts
	double			density(0);			// Material density (g/cm3)
	double			focal(3.5e7);		// Focal length in angstrom
	double			slit(0);			// Energy filter slit width
	Bstring			atom_select("all");	// Selection
	double			thickness(0);		// Specimen thickness
	double			mass(0);			// Molecular weight
	double			radius(0);			// Particle radius
	double			dose(0);			// Fluence
	double			def_min(0), def_max(0), def_step(0.1); // Defocus range
	Bstring 		atompropfile;		// Atom properties file
	Bstring 		respropfile;		// Residue properties file
	Bstring 		propfile("material.star");	// Material density and composition file
	Bstring			paramfile;

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "halfmaximum" ) halfmax = 1;
		if ( curropt->tag == "series" )
			apser = curropt->value.split_into_doubles(",");
		if ( curropt->tag == "angles" )
			angser = curropt->value.split_into_doubles(",");
		if ( curropt->tag == "material" )
			material_str = curropt->value;
		if ( curropt->tag == "elements" )
			elements = curropt->value;
		if ( curropt->tag == "counts" )
			count_str = curropt->value;
		if ( curropt->tag == "density" )
			if ( ( density = curropt->value.real() ) < 0.01 )
				cerr << "-density: A density must be specified!" << endl;
		if ( curropt->tag == "Volt" ) {
			if ( ( volts = curropt->value.real() ) < 1 )
				cerr << "-Volt: A voltage must be specified!" << endl;
			else
				if ( volts < 1e3 ) volts *= 1e3;	// Assume kilovolts
		}
		if ( curropt->tag == "aperture" ) {
			if ( ( aperture = curropt->value.real() ) < 1 )
				cerr << "-aperture: An aperture diameter must be specified!" << endl;
			else {
				if ( aperture < 1e4 ) aperture *= 1e4;	// Assume µm
			}
		}
		if ( curropt->tag == "focallength" ) {
			if ( ( focal = curropt->value.real() ) < 1 )
				cerr << "-focallength: A focal length must be specified!" << endl;
			else {
				if ( focal < 1e3 ) focal *= 1e7;	// Assume mm
			}
		}
		if ( curropt->tag == "slitwidth" )
			if ( ( slit = curropt->value.real() ) < 1 )
				cerr << "-slitwidth: A slit width must be specified!" << endl;
		if ( curropt->tag == "thickness" )
			if ( ( thickness = curropt->value.real() ) < 1 )
				cerr << "-thickness: A specimen thickness must be specified!" << endl;
		if ( curropt->tag == "mass" )
			mass = get_option_mass(curropt->value);
		if ( curropt->tag == "radius" )
			if ( ( radius = curropt->value.real() ) < 1 )
				cerr << "-radius: A particle radius must be specified!" << endl;
		if ( curropt->tag == "dose" )
			if ( ( dose = curropt->value.real() ) < 1 )
				cerr << "-dose: A total electron dose per pixel must be specified!" << endl;
		if ( curropt->tag == "defocus" ) {
			if ( curropt->values(def_min, def_max, def_step) < 1 )
				cerr << "-defocus: At least one defocus value must be specified!" << endl;
			else {
				if ( def_min < 20 ) def_min *= 1e4;	// Assume µm
				if ( def_max < 20 ) def_max *= 1e4;	// Assume µm
				if ( def_max < def_min ) def_max = def_min;
				if ( def_step < 10 ) def_step *= 1e4;
			}
		}
		if ( curropt->tag == "atoms" )
			atompropfile = curropt->filename();
		if ( curropt->tag == "residues" )
			respropfile = curropt->filename();
		if ( curropt->tag == "properties" )
			propfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
/*
	vector<string>			el;
	if ( elements.length() ) {
		stringstream			ss(elements.c_str());
		string					es;
		while ( getline(ss, es, ',') )
			el.push_back(es);
		cout << elements << tab << el.size() << endl;
	}
*/

	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);
	
	Bmaterial			material;
	map<string,Bcomptype>&	comp = material.composition();
	Bmolgroup*			molgroup = NULL;
	
	material.identifier("Unknown");
	
//	cout << "-" << material_str << "-" << endl;
	
	if ( material_str.length() ) {
		if ( file_type(material_str) == Molecule ) {
			molgroup = read_molecule(material_str, atom_select, paramfile);
			material = material_from_molgroup(molgroup);
		} else {
			material_str = "\"" + material_str + "\"";
			cout << "Requested material: " << material_str << endl;
			map<string,Bmaterial>	mprop = read_material_properties(propfile);
//			for ( auto m: mprop )
//				cout << m.second.identifier() << tab << m.second.density(1) << endl;
			material = mprop[material_str.str()];
		}
	} else if ( elements.length() ) {
		vector<string>		el = split(elements.str(),',');
		vector<long>		cnt;
		if ( count_str.length() )
			cnt = count_str.split_into_integers(",");
		for ( long i=0; i<el.size(); ++ i ) {
			comp[el[i]].identifier(el[i]);
			if ( i < cnt.size() ) comp[el[i]].component_count(cnt[i]);
			else comp[el[i]].component_count(1);
		}
	} else {
		material = protein_material_default();
	}
	
	if ( density ) material.density(density, 0);
	
	if ( material.density(1) < 0.01 ) material.density(RHO, 1);
	
	material.unit(1);	// Convert the density
	
	if ( comp.size() < 1 ) comp["C"].component_count(1);
	
	material.update_parameters(atompar);

	if ( mass ) material.mass(mass);
	
	mass = material.mass();
	
	// Calculate radius from mass, assuming a globular/spherical protein
	if ( radius < 1 )
		radius = pow((3*mass)/(4*M_PI*material.density(1)), 1.0/3.0);
	
	if ( thickness < 1 ) thickness = 2*radius;
	
	if ( verbose )
		material.show();

	CTFparam		ctf(volts, 2.7, 0.07);
	ctf.objective_aperture(aperture);
	ctf.focal_length(focal);
	ctf.slit_width(slit);
	ctf.Cs(2.7);
	ctf.amp_shift(0.07);

	double			scut = ctf.frequency_cutoff();

	if ( elements.length() ) {
		if ( halfmax )
			cross_section_half_maximal_frequencies(material);
		else
			cross_sections(material, ctf);
		bexit(0);
	}
	
	if ( verbose ) {
		cout << "General parameters:" << endl;
		cout << "Acceleration voltage:           " << 1e-3*volts << " kV" << endl;
		cout << "Beta:                           " << beta(volts) << endl;
		cout << "Aperture:                       " << 1e-4*aperture << " µm" << endl;
		cout << "Focal length:                   " << 1e-7*focal << " mm" << endl;
		cout << "Slit width:                     " << slit << " eV" << endl;
		cout << "Frequency cutoff:               " << scut << " /A (" <<
			1/scut << " A)" << endl;
		cout << "Specimen density:               " << material.density(1) << " Da/A3" << endl;
		cout << "Specimen thickness:             " << thickness << " A" << endl;
		cout << "Dose/fluence:                   " << dose << " e/px" << endl;
		if ( def_min ) {
			cout << "Defocus range:                  " << def_min*1e-4 << " - " << def_max*1e-4 << " µm" << endl;
			cout << "Defocus step:                   " << def_step*1e-4 << " µm" << endl;
		}
		cout << endl;
	}

	Bmaterial	ice = material_ice(atompar);
	
	double		def_fac = defocus_factor(ctf, def_min, def_max, def_step);
	double		intensity = signal_intensity(ice, thickness, ctf);
	double		emfp = effective_mean_free_path(ice, ctf);
	double		snre = particle_snr(material, mass, radius, thickness, ctf);
	
	if ( dose ) {
		cout << "Total intensity:                " << dose << " e/px" << endl;
		cout << "Ice EMFP:                       " << emfp << endl;
		cout << "Direct beam attenuation:        " << intensity << endl;
		cout << "Defocus factor:                 " << def_fac << endl;
		cout << "Particle SNR:                   " << dose*def_fac*snre*intensity << endl << endl;
	}

//	ice_intensity(100, thickness, atompar, ctf);
	
	if ( apser.size() )
		aperture_series(thickness, ice, ctf, apser);

	if ( angser.size() )
		collection_angle_series(thickness, ice, ctf, angser);

//	Bplot*		plot = particle_spectral_signal(protcomp, mass,
//		radius, thickness, atompar, ctf);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

/**
@brief 	Default protein composition adjusted by mass.
@return Bmaterial			material.
**/
Bmaterial		protein_material_default()
{
	Bmaterial			material;
	map<string,Bcomptype>&	comp = material.composition();
	
	material.identifier("Protein");
	material.density(RHO, 1);
	
	// Atoms per thousand
	comp["H"].component_count(498);
	comp["C"].component_count(320);
	comp["N"].component_count(85);
	comp["O"].component_count(95);
	comp["S"].component_count(2);
	
	return material;
}

Bmaterial	material_from_molgroup(Bmolgroup* molgroup)
{
	Bmaterial			material;
	map<string,Bcomptype>&	comp = material.composition();
	material.density(RHO, 1);
	
	Bmolecule*			mol;
	Bresidue*			res;
	Batom*  			atom;

    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				if ( comp.find(atom->el) != comp.end() )
					comp[atom->el].component_count_increment();
				else
					comp[atom->el].component_count(1);
	
	return material;
}
/*
double		material_update_parameters(Bmaterial& material, map<string,Bcomptype>& atompar)
{
	double			mass(0);
	
	if ( verbose )
		cout << "Updating material parameters" << endl;
	
	map<string,Bcomptype>&	comp = material.composition();

	for ( auto &it: comp ) {
		Bcomptype&		ct = it.second;
		if ( atompar.find(it.first) != atompar.end() ) {
			Bcomptype&	at = atompar[it.first];
			ct.identifier(at.identifier());
			ct.index(at.index());
			ct.mass(at.mass());
			ct.charge(at.charge());
			ct.coefficients(at.coefficients());
			mass += ct.mass();
			cout << ct.identifier() << tab << ct.mass() << endl;
		}
	}

	if ( verbose )
		cout << "Mass:              " << mass << endl << endl;

//	material.show();
	
	return mass;
}
*/
double		elastic_cross_section_lenz(long Z, double volts)
{
	double			f = 1e10*PLANCK/(EMASS*LIGHTSPEED);
	double			b2 = beta2(volts);

	double			cs = (f*f/(M_PI*b2))*pow(Z, 4.0/3.0);
	
	return cs;
}

double		elastic_cross_section_langmore(long Z, double volts)
{
	double			b2 = beta2(volts), b = sqrt(b2);

	double			cs = 1.4e-4*(pow(Z, 1.5)/b2)*(1-0.26*Z/(137*b));
	
	return cs;
}



int			cross_sections(vector<string>& el, map<string,Bcomptype>& atompar, CTFparam& ctf)
{
	if ( el.size() < 1 || el[0] == "all" ) el = all_elements(atompar);
	
//	long			nscat(0), i, j, t;
//	double			ds(0.01), f, sci, cs, csi, csa, csin;
	double			sci, cs, csi, csa, csin;
	
	double			lambda = electron_wavelength(ctf.volt());
	double			scut(ctf.frequency_cutoff());
	double			b2 = beta2(ctf.volt());
//	double			hf = 1e10*PLANCK/(EMASS*LIGHTSPEED);
	
	if ( verbose ) {
		cout << "Cross sections:" << endl;
		cout << "Acceleration voltage:           " << ctf.volt()*1e-3 << " kV" << endl;
		cout << "Aperture diameter:              " << ctf.objective_aperture()*1e-4 << " µm" << endl;
		cout << "Focal length:                   " << ctf.focal_length()*1e-7 << " mm" << endl;
		cout << "Wavelength:                     " << lambda << endl;
		cout << "Beta:                           " << sqrt(b2) << endl;
		cout << "Cutoff frequency:               " << scut << " /Å" << endl;
//		cout << "Element\tZ\tcs\tcs(ap)\tcs(ap)\tcsin(lang)" << endl;
		cout << "Element\tZ\tsci\tcs\tcsb2\tcsin(lang)\tcsinb2" << endl;
	}
	
	for ( auto el1: el ) {
		if ( atompar.find(el1) != atompar.end() ) {
			Bcomptype&	at = atompar[el1];
			sci = scatter_curve_integral(at);
			cs = elastic_cross_section(at, ctf.volt());
			csi = elastic_cross_section_integrated(at, ctf);
			csa = cs*scut*scut/(0.44*0.44 + scut*scut);	// Only for carbon
			csin = inelastic_cross_section_langmore(at.index(), ctf.volt());
//			cout << at->el << tab << at->z << tab << cs << tab << csi << tab << csa << tab << csin << endl;
			cout << at.identifier() << tab << at.index() << tab << sci << tab
				<< cs << tab << cs*b2 << tab << csin << tab
				<< csin*b2 << endl;
		}
	}
	
	
	return 0;
}

int			cross_sections(Bmaterial& material, CTFparam& ctf)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	double			sci, cs, csi, csa, csin;
	double			lambda = electron_wavelength(ctf.volt());
	double			scut(ctf.frequency_cutoff());
	double			b2 = beta2(ctf.volt());
	
	if ( verbose ) {
		cout << "Cross sections:" << endl;
		cout << "Acceleration voltage:           " << ctf.volt()*1e-3 << " kV" << endl;
		cout << "Aperture diameter:              " << ctf.objective_aperture()*1e-4 << " µm" << endl;
		cout << "Focal length:                   " << ctf.focal_length()*1e-7 << " mm" << endl;
		cout << "Wavelength:                     " << lambda << endl;
		cout << "Beta:                           " << sqrt(b2) << endl;
		cout << "Cutoff frequency:               " << scut << " /Å" << endl;
//		cout << "Element\tZ\tcs\tcs(ap)\tcs(ap)\tcsin(lang)" << endl;
		cout << "Element\tZ\tsci\tcs\tcsb2\tcsin(lang)\tcsinb2" << endl;
	}
	
	for ( auto el1: comp ) {
		Bcomptype&	ct = el1.second;
		sci = scatter_curve_integral(ct);
		cs = elastic_cross_section(ct, ctf.volt());
		csi = elastic_cross_section_integrated(ct, ctf);
		csa = cs*scut*scut/(0.44*0.44 + scut*scut);	// Only for carbon
		csin = inelastic_cross_section_langmore(ct.index(), ctf.volt());
//		cout << at->el << tab << at->z << tab << cs << tab << csi << tab << csa << tab << csin << endl;
		cout << ct.identifier() << tab << ct.index() << tab << sci << tab
				<< cs << tab << cs*b2 << tab << csin << tab
				<< csin*b2 << endl;
	}
	
	
	return 0;
}



/*
	Calculates the half-maximal cross section to estimate the mid-cutoff
	frequency coefficient for each element to model the aperture effect.
 */
int			cross_section_half_maximal_frequencies(vector<string>& el, map<string,Bcomptype>& atompar)
{
	if ( el[0] == "all" ) el = all_elements(atompar);

	double			sm;

//	if ( verbose ) {
		cout << "Half-maximal frequencies:" << endl;
		cout << "Element\tZ\tsm" << endl;
//	}

	for ( auto el1: el ) {
		if ( atompar.find(el1) != atompar.end() ) {
			Bcomptype&	at = atompar[el1];
			sm = cross_section_half_maximal_frequency(at);
			cout << setprecision(3) << at.identifier() << tab << at.index() << tab << sm << endl;
		}
	}
	
	return 0;
}

int			cross_section_half_maximal_frequencies(Bmaterial& material)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	double			sm;

//	if ( verbose ) {
		cout << "Half-maximal frequencies:" << endl;
		cout << "Element\tZ\tsm" << endl;
//	}

	for ( auto el1: comp ) {
		Bcomptype&		ct = el1.second;
		sm = cross_section_half_maximal_frequency(ct);
		cout << setprecision(3) << ct.identifier() << tab << ct.index() << tab << sm << endl;
	}
	
	return 0;
}

int			aperture_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> apser)
{
//	map<string,Bcomptype>& 	comp = material.composition();
	
	double		oa, ang, Iunf, Ifil, lam_unf, lam_fil;
	
	cout << "Aper\ta(mr)\tLunf\tLfil\tIunf\tIfil\tln(Iu)\tln(If)\tIf/Iu\tln(If/Iu)" << endl;
	for ( auto it = apser.begin(); it != apser.end(); ++it ) {
		oa = *it;
		if ( oa < 1e4 ) oa *= 1e4;
		ctf.objective_aperture(oa);
		ang = 1e3*oa/(2*ctf.focal_length());
		ctf.slit_width(0);
		Iunf = signal_intensity(material, thickness, ctf);
		lam_unf = -thickness/log(Iunf);
		ctf.slit_width(20);
		Ifil = signal_intensity(material, thickness, ctf);
		lam_fil = -thickness/log(Ifil);
		cout << *it << tab << setprecision(4) << ang << tab
			<< lam_unf << tab << lam_fil << tab
			<< Iunf << tab << Ifil << tab
			<< -log(Iunf) << tab << -log(Ifil) << tab
			<< Ifil/Iunf << tab << -log(Ifil/Iunf) << endl;
	}

	return 0;
}

int			collection_angle_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> angser)
{
//	map<string,Bcomptype>& 	comp = material.composition();

	double		oa, ang, Iunf, Ifil, lam_unf, lam_fil;
	
	cout << "Aper\ta(mr)\tLunf\tLfil\tIunf\tIfil\tln(Iu)\tln(If)\tIf/Iu\tln(If/Iu)" << endl;
	for ( auto it = angser.begin(); it != angser.end(); ++it ) {
		ang = *it;
		if ( ang > 1 ) ang *= 1e-3;
		oa = 2*ang*ctf.focal_length();
		ctf.objective_aperture(oa);
		ctf.slit_width(0);
		Iunf = signal_intensity(material, thickness, ctf);
		lam_unf = -thickness/log(Iunf);
		ctf.slit_width(20);
		Ifil = signal_intensity(material, thickness, ctf);
		lam_fil = -thickness/log(Ifil);
		cout << oa*1e-4 << tab << setprecision(4) << ang*1e3 << tab
			<< lam_unf << tab << lam_fil << tab
			<< Iunf << tab << Ifil << tab
			<< -log(Iunf) << tab << -log(Ifil) << tab
			<< Ifil/Iunf << tab << -log(Ifil/Iunf) << endl;
	}

	return 0;
}


double		particle_snr(Bmaterial& material, double mass, double radius, double thickness, CTFparam& ctf)
{
	map<string,Bcomptype>& comp = material.composition();

	if ( mass < 1 ) mass = material.mass();

	// Calculate radius from mass, assuming a globular/spherical protein
	if ( radius < 1 ) radius = 50;
		radius = pow((3*mass)/(4*M_PI*material.density(1)), 1.0/3.0);

	if ( thickness < 1 ) thickness = 2*radius;

	double			cs = elastic_cross_section_integrated(comp, ctf);
	double			ics = inelastic_cross_section_langmore(comp, ctf.volt());
	double			emfp = effective_mean_free_path(material, ctf);
	double			snre = 0.5*cs/(M_PI*radius*radius);

	if ( verbose ) {
		cout << "Particle:" << endl;
		cout << "Molecular weight:               " << mass << " Da" << endl;
		cout << "Particle radius:                " << radius << " A" << endl;
		cout << "Density:                        " << material.density(0) << " g/cm3" << endl;
		cout << "Density:                        " << material.density(1) << " Da/A3" << endl;
		cout << "Density:                        " << material.density(2) << " #/A3" << endl;
		cout << "Occupancy:                      " << mass/(M_PI*radius*radius*thickness) << " Da/A3" << endl;
		cout << "Elastic cross section:          " << cs << " A2" << endl;
		cout << "Inelastic cross section:        " << ics << " A2" << endl;
		cout << "Particle SNR per electron:      " << snre << " /e" << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
		cout << endl;
	}
	
	return snre;
}

Bplot*		particle_spectral_signal(Bmaterial& material, double mass,
	double radius, double thickness, CTFparam& ctf)
{
	map<string,Bcomptype>& comp = material.composition();

	double			ds(0.5/radius);
	double			scut(ctf.frequency_cutoff());

	map<string, vector<double>>		sc = calculate_scattering_curves(comp, ds, scut);

	string			e;
	long			i, ns(sc.begin()->second.size());
	double			pf(ctf.lambda()/sqrt(1-beta2(ctf.volt())));
	double			w, area(M_PI * radius * radius);
	Bplot*			plot = new Bplot(1,ns,2);
	
	cout << "Prefactor: " << pf << endl;
	
	for ( i=0; i<ns; ++i ) {
		(*plot)[i] = ds*i;
		(*plot)[i+ns] = 0;
	}

	for ( auto it: comp ) {
		e = it.first;
		w = pf * it.second.component_count();
		cout << e << tab << w << endl;
		i = ns;
		for ( auto is = sc[e].begin(); is != sc[e].end(); ++is )
			(*plot)[i++] += *is * w;
	}
	
	vector<double>	dp = defocus_range_profile(ctf, ds, 5e3, 2e4, 1e3);
	
	for ( i=0; i<ns; ++i ) {
		(*plot)[i+ns] *= (*plot)[i+ns];
		(*plot)[i+ns] *= dp[i]/area;
	}

	cout << "s\tAmp\t" << endl;
	for ( i=0; i<ns; ++i )
		cout << (*plot)[i] << tab << (*plot)[i+ns] << tab << dp[i] << endl;

	return plot;
}
