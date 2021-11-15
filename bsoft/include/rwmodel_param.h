/**
@file	rwmodel_param.h
@brief	Header to read and write model dynamics parameters in STAR format
@author Bernard Heymann
@date	Created: 20100305
@date	Modified: 20210224
**/

#include "rwmodel.h"

#ifndef _Bmodparam_
#define _Bmodparam_
#define	DISTMAT_COMPTYPE1	"distance_matrix.component_type1"
#define	DISTMAT_COMPTYPE2	"distance_matrix.component_type2"
#define	DISTMAT_LLEN		"distance_matrix.link_length"
#define	DISTMAT_DIST		"distance_matrix.distance"
#define	DISTMAT_KL			"distance_matrix.link_constant"
#define	DISTMAT_KD			"distance_matrix.distance_constant"
#define	DISTMAT_POTTYPE		"distance_matrix.potential_type"

#define ANGLEMAT_COMPTYPE1	"angle_matrix.component_type1"
#define ANGLEMAT_COMPTYPE2	"angle_matrix.component_type2"
#define ANGLEMAT_COMPTYPE3	"angle_matrix.component_type3"
#define ANGLEMAT_ANGLE		"angle_matrix.angle"
#define ANGLEMAT_KA			"angle_matrix.angle_constant"

/************************************************************************
@Object: class Bmodparam
@Description:
	A structure used for model mechanics.
@Features:
	This defines variables and constants for doing model dynamics or
	other model mechanics.
*************************************************************************/
class Bmodparam {
public:
	string		comment;
	double		timestep;		// Integration time step for MD
	double		Kfriction;		// Friction coefficient
	double		Kdistance;		// Distance/Van der Waals force constant
	double		Kelec;			// Electrostatic force constant
	double		Klink;			// Link/bond force constant
	double		Kangle;			// Angular force constant
	double		Kpolyangle;		// Angular force constant for polygons
	double		Kpolygon;		// Polygon regularity force constant
	double		Kpolyplane;		// Polygon planarity force constant
	double		Kpoint;			// Point force constant
	double		Kradial;		// Radial force constant
	double		Kplane;			// Neighbor plane force constant
	double		Kguide;			// Polyhedron guide force constant
	double		Kmap;			// Density map force constant
	double		sepdist;		// Separation distance grid sampling
	double		cutoff;			// Distance cutoff for non-bonded calculations
	double		pointdecay;		// Decay constant for the point force
	double		radius;			// Radius for radial force
	Vector3<double>	point;		// Center of point or radial force
	Bmodel*		guide;			// Polyhedron guide model
	int			linksteps;		// Number of sampling intervals along a link
	int			wrap;			// Flag to turn periodic boundaries on
	double		sigma;			// Gaussian decay for density fitting
	double		Edistance;		// Non-linked distance energy
	double		Eelec;			// Electrostatic energy
	double		Elink;			// Link energy
	double		Eangle;			// Angle energy
	double		Epolyangle;		// Angle energy for polygons
	double		Epolygon;		// Polygon regularity energy
	double		Epolyplane;		// Polygon planarity energy
	double		Epoint;			// Point force energy
	double		Eradial;		// Radial force energy
	double		Eplane;			// Neighbor plane force energy
	double		Eguide;			// Polyhedron guide force energy
	double		Emap;			// Density map associated energy
	double		Ekin;			// Kinetic energy
	double		Epot;			// Potential energy
	map<string,Bcomptype>		comptype;
	vector<vector<Blinktype>>	linktype;
	vector<vector<vector<Bangletype>>>	angletype;
private:
	void	initialize() {
		comment = "?";
		timestep = 1;
		Kfriction = 1;
		Kdistance = 0;
		Kelec = 0;
		Klink = 0;
		Kangle = 0;
		Kpolyangle = 0;
		Kpolygon = 0;
		Kpolyplane = 0;
		Kpoint = 0;
		Kradial = 0;
		Kplane = 0;
		Kguide = 0;
		Kmap = 0;
		sepdist = 10;
		cutoff = 10;
		pointdecay = 0.01;
		radius = 0;
		linksteps = 3;
		wrap = 0;
		sigma = 0;
	}
public:
	Bmodparam()	{ initialize(); }
} ;
#endif

#ifndef _Bmaterial_
#define _Bmaterial_
#define	MATERIAL_NAME		"material.name"
#define	MATERIAL_DENSITY	"material.density"
#define	MATERIAL_DENS_UNIT	"material.density_unit"

class Bmaterial {
private:
	string				id;		// Material name
	double				den;	// Density in g/cm3
	int					unt;	// Density units: 0 = g/cm3, 1 = Da/Å3, 2 = number/Å3
	long				sel;	// Selection
	map<string,Bcomptype>	comp;	// Elemental composition
public:
	void		identifier(string s) { id = s; }
	string		identifier() { return id; }
	void		select(long i) { sel = i; }
	long		select() { return sel; }
	void		density(double d, int u) { den = d; unt = u; }
	double		density(int u) {
		if ( u == 0 ) return gram_per_cm3();
		else if ( u == 1 ) return dalton_per_angstrom3();
		else return number_per_angstrom3();
	}
	double		gram_per_cm3() {
		if ( unt == 1 ) return den*1.0e24/AVOGADRO;
		else if ( unt == 2 ) return den*mass()*1.0e24/AVOGADRO;
		return den;
	}
	double		dalton_per_angstrom3() {
		if ( unt == 0 ) return AVOGADRO*den/1.0e24;
		else if ( unt == 2 ) return den*mass();
		return den;
	}
	double		number_per_angstrom3() {
		if ( unt == 0 ) return AVOGADRO*den/(1.0e24*mass());
		else if ( unt == 1 ) return den/mass();
		return den;
	}
	double		unit(int u) {
		if ( u == unt ) return den;
		if ( u == 0 ) den = gram_per_cm3();
		else if ( u == 1 ) den = dalton_per_angstrom3();
		else den = number_per_angstrom3();
		unt = u;
		return den;
	}
	int			unit() { return unt; }
	map<string,Bcomptype>&	composition() { return comp; }
	Bcomptype&	operator[](string s) { return comp[s]; }
	long		number() {
		long		n(0);
		for ( auto i: comp ) n += i.second.component_count();
		return n;
	}
	double		mass() {
		double		m(0);
		for ( auto ct: comp ) m += ct.second.mass() * ct.second.component_count();
		return m;
	}
	void		mass(double m) {
		double		r(m/mass());
		for ( auto &ct: comp ) ct.second.component_count(ct.second.component_count()*r);
	}
	void		update_parameters(map<string,Bcomptype>& types) {
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			if ( types.find(it.first) != types.end() ) {
				Bcomptype&	at = types[it.first];
				ct.identifier(at.identifier());
				ct.index(at.index());
				ct.mass(at.mass());
				ct.charge(at.charge());
				ct.coefficients(at.coefficients());
			}
		}
	}
	void		show() {
		cout << "Material:                       " << id << endl;
		cout << "Density:                        " << den;
		if ( unt == 0 ) cout << " g/cm3" << endl;
		else if ( unt == 1 ) cout << " Da/A3" << endl;
		else cout << " #/A3" << endl;
		cout << "Composition:" << endl;
		cout << "Element\tCount\tMass" << endl;
		for ( auto ct: comp )
			cout << ct.second.identifier() << tab << ct.second.component_count()
				<< tab << ct.second.mass() << endl;
		cout << "Mass:                           " << mass() << " Da" << endl << endl;
	}
} ;

#endif

// Function prototypes
Bmodparam	model_param_generate(Bmodel* model);
int			model_param_generate(Bmodparam& md, Bmodel* model);
int			model_param_set_type_indices(Bmodel* model, Bmodparam& mp);
map<string,Bcomptype> 	read_atom_properties(Bstring& filename);
map<string,Bmaterial> 	read_material_properties(Bstring& filename);
int			write_material_properties(Bstring& filename, map<string,Bmaterial> material);
Bmodparam 	read_dynamics_parameters(Bstring& filename);
int			update_dynamics_parameters(Bmodparam& md, Bstring& filename);
int		 	write_dynamics_parameters(Bstring& filename, Bmodparam& md);

