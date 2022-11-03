/**
@file	model_create.cpp
@brief	A tool to expand models.
@author Bernard Heymann
@date	Created: 20090714
@date	Modified: 20220210
**/

#include "rwmodel.h"
#include "model_create.h"
#include "model_transform.h"
#include "model_symmetry.h"
#include "model_util.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Adds a component to a model.
@param 	*model			model list.
@param 	&id				model identifier.
@param 	&type			component type.
@param 	loc				component location.
@return Bcomponent*		new component.

	The component is added to the model indicated by the ID, the type and
	the location given.

**/
Bcomponent*	model_add_component(Bmodel* model, Bstring &id, Bstring &type, Vector3<double> loc)
{
	if ( !model ) return NULL;
	
	Bmodel*			mp = model;
	
	if ( id.length() ) {
		for ( mp = model; mp; mp = mp->next )
			if ( mp->identifier() == id.str() ) break;
		if ( !mp ) mp = model;
	}

	if ( type.length() < 1 ) {
		if ( model->type ) type = model->type->identifier();
		else type = "VER";
	}

	Bcomptype*		ct = model->add_type(type);
	
	long			n;
	Bcomponent*		comp = NULL;

	for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;

	string			cid = to_string(++n);
	
//	comp = component_add(&mp->comp, cid);
	comp = mp->add_component(cid);
	comp->location(loc);
	comp->type(ct);
	
	return comp;
}

/**
@brief 	Adds a set of components to a model.
@param 	*model			model list.
@param 	&id				model identifier.
@param 	&type			component type.
@param 	loc				component locations.
@return Bcomponent*		last new component.

	The component is added to the model indicated by the ID, the type and
	the location given.

**/
Bcomponent*	model_add_components(Bmodel* model, Bstring &id, Bstring &type, vector<Vector3<double>> loc)
{
	if ( !model ) return NULL;
	
	Bmodel*			mp = model;
	
	if ( id.length() ) {
		for ( mp = model; mp; mp = mp->next )
			if ( mp->identifier() == id.str() ) break;
		if ( !mp ) mp = model;
	}

	if ( type.length() < 1 ) {
		if ( model->type ) type = model->type->identifier();
		else type = "VER";
	}

	Bcomptype*		ct = model->add_type(type);
	
	long			n;
	Bcomponent*		comp = NULL;

	for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;

	string			cid = to_string(++n);
	
	for ( auto p: loc ) {
		comp = mp->add_component(cid);
		comp->location(p);
		comp->type(ct);
	}
	
	return comp;
}

/**
@brief 	Generates a platonic polyhedron.
@param 	*sym			symmetry.
@param 	radius			distance of each vertex from origin.
@return Bmodel*			new model.

	Polyhedral components are generated on symmetry axes.
	Symmetry designations supported:
		T-3		tetrahedron.
		O-2		truncated octahedron.
		O-3		cube.
		O-4		octahedron.
		I-2		truncated icosahedron.
		I-3		dodecahedron.
		I-5		icosahedron.

**/
Bmodel*		model_platonic(Bsymmetry& sym, double radius)
{
	string			id("1");
	Bmodel*			model = new Bmodel(id);
	model->symmetry(sym.label().str());
	
	string			comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
//	Bcomponent*		comp = component_add(&model->comp, id);
	Bcomponent*		comp = model->add_component(id);
	comp->type(ct);
	comp->radius(radius/10);

//	double			search_radius;
	
	if ( model->symmetry() == "T-3" ) {
		model->identifier() = "Tetrahedron";
//		search_radius = 2*radius;
		comp->location(Vector3<double>(radius/sqrt(3), radius/sqrt(3), radius/sqrt(3)));
	} else if ( model->symmetry() == "O-2" ) {
		model->identifier() = "Truncated octahedron";
//		search_radius = 1.2*radius;
		comp->location(Vector3<double>(0, -radius/sqrt(2), radius/sqrt(2)));
	} else if ( model->symmetry() == "O-3" ) {
		model->identifier() = "Cube";
//		search_radius = 1.5*radius;
		comp->location(Vector3<double>(radius/sqrt(3), radius/sqrt(3), radius/sqrt(3)));
	} else if ( model->symmetry() == "O-4" ) {
		model->identifier() = "Octahedron";
//		search_radius = 1.5*radius;
		comp->location(Vector3<double>(0, 0, radius));
	} else if ( model->symmetry() == "I-2" ) {
		model->identifier() = "Truncated icosahedron";
//		search_radius = 0.8*radius;
		comp->location(Vector3<double>(0, 0, radius));
	} else if ( model->symmetry() == "I-3" ) {
		model->identifier() = "Dodecahedron";
//		search_radius = radius;
		comp->location(Vector3<double>(radius/sqrt(3), radius/sqrt(3), radius/sqrt(3)));
	} else if ( model->symmetry() == "I-5" ) {
		model->identifier() = "Icosahedron";
//		search_radius = 1.2*radius;
		comp->location(Vector3<double>(0, radius/sqrt(2+GOLDEN), radius/sqrt(3-GOLDEN)));
	} else {
		cerr << "Error: Symmetry " << sym.label() << " not supported!" << endl;
		model_kill(model);
		return NULL;
	}
	
	if ( verbose ) {
		cout << "Generating a " << model->identifier() << ":" << endl;
		cout << "Symmetry:                       " << model->symmetry() << endl;
		cout << "Radius:                         " << radius << endl << endl;
	}

	comp->view(View2<float>(comp->location()));

	Vector3<double>	origin;
	View			ref_view;
	model_apply_point_group(model, model->symmetry(), origin, ref_view, 1);

	return model;
}

/**
@brief 	Generates a helix model.
@param 	radius			distance of each vertex from helical axis.
@param 	helix_rise		helical rise.
@param 	helix_angle		helical rotation angle.
@param 	ncomp			number of components.
@return Bmodel*			new model.
**/
Bmodel*		model_helix(double radius, double helix_rise, double helix_angle, long ncomp)
{
	long			i;
	double			start = -(ncomp/2.0 - 0.5)*helix_rise;
//	double			search_radius = 1.2*sqrt(helix_rise*helix_rise + 2*radius*radius*(1 - cos(helix_angle)));
	string			id("Helix");
	Bstring			symstr(Bstring(helix_rise, "H%.3g") + Bstring(helix_angle*180.0/M_PI, ",%.3g"));
	Bmodel*			model = new Bmodel(id);
	model->symmetry(symstr.str());
	string			comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
//	Bcomponent*		comp = component_add(&model->comp, id);
	Bcomponent*		comp = model->add_component(1);
	comp->type(ct);
	comp->location()[0] = radius;
	comp->location()[2] = start;
	comp->radius(radius/10);
	
	if ( verbose ) {
		cout << "Generating a " << model->identifier() << ":" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Helical rise and rotation:      " << helix_rise << " " << helix_angle*180.0/M_PI << endl;
		cout << "Number of components:           " << ncomp << endl << endl;
	}

	for ( i=1; i<ncomp; i++ ) {
//		id = to_string(i+1);
//		comp = component_add(&comp, id);
		comp = comp->add(i+1);
		comp->type(ct);
		comp->location()[0] = radius*cos(i*helix_angle);
		comp->location()[1] = radius*sin(i*helix_angle);
		comp->location()[2] = start + i*helix_rise;
		comp->radius(radius/10);
	}
	
	return model;
}

/**
@brief 	Generates components at random non-overlapping locations and with random views.
@param 	ncomp			number of components.
@param 	comp_radius		component radius.
@param 	max_radius		maximum radius of components.
@return Bmodel*			new model.

	If a new component overlaps within an existing component, as defined
	by the component radius, new random coordinates are generated for it.

**/
Bmodel*		model_random(long ncomp, double comp_radius, double max_radius)
{
	Bstring			id("Random");
	double			R = ncomp * pow(comp_radius/max_radius, 3);
		
	if ( verbose ) {
		cout << "Generating a random model:      " << id << endl;
		cout << "Number of components:           " << ncomp << endl;
		cout << "Component radius:               " << comp_radius << endl;
		cout << "Maximum model radius:           " << max_radius << endl;
		cout << "Packing density:                " << R << endl << endl;
	}
	
	if ( R > 0.35 ) {
		cerr << "Error: The packing density is too high!" << endl << endl;
		return NULL;
	}

	random_seed();
	
	long				i;
	double			inner_radius(max_radius - comp_radius);
	double			d;
	Bmodel*			model = new Bmodel(id);
	Bstring			comptype = "VER";
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;
	Bcomptype*		ct = model->add_type(comptype);
	
	for ( i=1; i<=ncomp; i++ ) {
//		id = Bstring(i, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(i);
		else comp = model->add_component(i);
		comp->type(ct);
		for ( d = 0; d < 2*comp_radius;  ) {
			comp->location(vector3_random(inner_radius));
			for ( d = 4*comp_radius, comp2 = model->comp; 
					comp2 != comp && d >= 2*comp_radius; comp2 = comp2->next )
				d = comp->location().distance(comp2->location());
		}
		comp->view().random_view();
		comp->radius(comp_radius);
	}
	
	return model;
}

/**
@brief 	Generates components at random non-overlapping locations and with random views.
@param 	ncomp			number of components.
@param 	comp_radius		component radius.
@param 	min				minimum bounds.
@param 	max				maximum bounds.
@return Bmodel*			new model.

	If a new component overlaps within an existing component, as defined
	by the component radius, new random coordinates are generated for it.

**/
Bmodel*		model_random(long ncomp, double comp_radius, Vector3<double> min, Vector3<double> max)
{
	Bstring			id("Random");
	Vector3<double>	cmin(min + comp_radius);
	Vector3<double>	cmax(max - comp_radius);
	double			R(ncomp * (4.0*M_PI/3.0)*comp_radius*comp_radius*comp_radius/(max - min).volume());
		
	if ( verbose ) {
		cout << "Generating a random model:      " << id << endl;
		cout << "Number of components:           " << ncomp << endl;
		cout << "Component radius:               " << comp_radius << endl;
		cout << "Minimum coordinates:            " << cmin << endl;
		cout << "Maximum coordinates:            " << cmax << endl;
		cout << "Packing density:                " << R << endl << endl;
	}
	
	if ( R > 0.35 ) {
		cerr << "Error: The packing density is too high!" << endl << endl;
		return NULL;
	}

	random_seed();
	
	long				i;
	double			d;
	Bmodel*			model = new Bmodel(id);
	Bstring			comptype = "VER";
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;
	Bcomptype*		ct = model->add_type(comptype);
	
	for ( i=1; i<=ncomp; i++ ) {
//		id = Bstring(i, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(i);
		else comp = model->add_component(i);
		comp->type(ct);
		for ( d = 0; d < 2*comp_radius;  ) {
			comp->location(vector3_random(cmin, cmax));
			for ( d = 4*comp_radius, comp2 = model->comp; 
					comp2 != comp && d >= 2*comp_radius; comp2 = comp2->next )
				d = comp->location().distance(comp2->location());
		}
		comp->view().random_view();
		comp->radius(comp_radius);
	}
	
	return model;
}

/**
@brief 	Generates components at random gaussian-distributed locations and with random views.
@param 	ncomp			number of components.
@param 	std				standard deviation of gaussian distribution.
@return Bmodel*			new model.

	Overlapping components are generated.

**/
Bmodel*		model_random_gaussian(long ncomp, double std)
{
	Bstring			id("RandomGaussian");
		
	if ( verbose ) {
		cout << "Generating a random gaussian model: " << id << endl;
		cout << "Number of components:           " << ncomp << endl;
	}
	
	random_seed();
	
	long				i;
	Bmodel*			model = new Bmodel(id);
	Bstring			comptype = "VER";
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = model->add_type(comptype);
	
	for ( i=1; i<=ncomp; i++ ) {
//		id = Bstring(i, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(i);
		else comp = model->add_component(i);
		comp->type(ct);
		comp->location(vector3_random_gaussian(0, std));
		comp->view().random_view();
	}

	double			R = model_gyration_radius(model);
	
	if ( verbose )
		cout << "Radius of gyration:             " << R << " A" << endl << endl;
	
	return model;
}

/**
@brief 	Generates random components on a shell.
@param 	ncomp			number of components.
@param 	radius			shell radius.
@return Bmodel*			new model.

	Overlapping components are generated.

**/
Bmodel*		model_random_shell(long ncomp, double radius)
{
	Bstring			id("RandomShell");
		
	if ( verbose ) {
		cout << "Generating a random shell model: " << id << endl;
		cout << "Number of components:           " << ncomp << endl;
		cout << "Radius:                         " << radius << endl;
	}
	
	random_seed();
	
	long				i;
	Bmodel*			model = new Bmodel(id);
	Bstring			comptype = "VER";
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = model->add_type(comptype);
	
	for ( i=1; i<=ncomp; i++ ) {
//		id = Bstring(i, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(i);
		else comp = model->add_component(i);
		comp->type(ct);
		comp->location(vector3_random_unit_sphere());
		comp->view(View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0));
		comp->scale(radius);
	}
	
	double			R = model_gyration_radius(model);
	
	if ( verbose )
		cout << "Radius of gyration:             " << R << " A" << endl << endl;
	
	return model;
}

/**
@brief 	Generates random components on a shell.
@param 	ncomp			number of components.
@param 	radius			shell radius.
@param 	separation		minimum vertex separation distance.
@return Bmodel*			new model.

	The vertices are created with a desired separation, which means that
	if too many vertices are requested with too large separation, the 
	function may never return. The following must therefore be true:
		n*(d/r)^2 < 8
	where
		n: number of vertices
		d: minimum separation distance
		r: polyhedron radius
	If this not true, the number of vertices is decreased.

**/
Bmodel*		model_random_shell(long ncomp, double radius, double separation)
{
	Bstring			id("RandomShell");
		
	if ( verbose ) {
		cout << "Generating a random shell model: " << id << endl;
		cout << "Number of components:           " << ncomp << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Separation:                     " << separation << endl;
	}
	
	random_seed();
	
	long				i;
	double			d, pd(ncomp*separation*separation/(radius*radius));

	if ( pd > 8 ) ncomp = 8*radius*radius/(separation*separation);
	
	Bmodel*			model = new Bmodel(id);
	Bstring			comptype = "VER";
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;
	Bcomptype*		ct = model->add_type(comptype);
	
	for ( i=1; i<=ncomp; i++ ) {
//		id = Bstring(i, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(i);
		else comp = model->add_component(i);
		comp->type(ct);
		comp->radius(separation/10);
		for ( d = 0; d < separation;  ) {
			comp->location(vector3_random_unit_sphere());
			comp->scale(radius);
			for ( d = 2*separation, comp2 = model->comp; 
					comp2 != comp && d >= separation; comp2 = comp2->next )
				d = comp->location().distance(comp2->location());
		}
		comp->view(View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0));
	}
	
	double			R = model_gyration_radius(model);
	model->FOM(R);
	
	return model;
}

/**
@brief 	Creates a spherical shell point model.
@param 	number			number of points.
@param 	radius			shell radius.
@param 	separation		distance between points.
@return Bmodel*			new shell model.

	Two of the three arguments need to be given (the other zero).
	If all three arguments are given, only the radius and distance is used.

**/
Bmodel*		model_create_shell(long number, double radius, double separation)
{
	// flag = 0: the number is set by generation
	// flag = 1: the radius must be derived
	// flag = 2: the separation must be derived
	int				flag(0);
	if ( radius <= 0 ) flag = 1;
	else if ( separation <= 0 ) {
		if ( number <= 0 ) separation = radius/10;
		else flag = 2;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_create_shell: flag=" << flag << " number=" << number << " radius=" << radius << " separation=" << separation << endl;
	
	string			id("Shell");
	Bmodel*			model = new Bmodel(id);
	
	Bcomponent*		comp = NULL;
	string			comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
	
	long			i(number), n, t;
	double			x, y, z;
	double			theta_step, theta, phi, phi_step, rad_sin_theta;
	double			lam(0.1), rn;

	random_seed();
	
	double			irm = lam/get_rand_max();
	
	if ( flag ) {
		if ( flag == 1 ) radius = separation*sqrt(number/12.5L);
		else separation = radius/sqrt(number/12.5L);
		for ( t=i=0; t<1000 && i != number; t++ ) {
			theta_step = sqrt(3.0/4.0)*separation/radius;
			for ( i=0, theta = 0; theta <= M_PI; theta += theta_step ) {
				rad_sin_theta = radius*sin(theta);
				n = (long) (TWOPI*rad_sin_theta/separation);
				phi_step = 1.5*TWOPI;
				if ( theta > 0 && theta < TWOPI )
					phi_step = TWOPI*1.0L/n;
				for ( phi = 0; phi < TWOPI - phi_step/2; phi += phi_step, i++ ) ;
			}
			if ( i != number ) {
				rn = irm*random();
				if ( flag == 1 ) radius *= 1 + (number - i)*rn/number;
				else separation *= 1 + (i - number)*rn/number;
			}
		}
	}
	
	theta_step = sqrt(3.0/4.0)*separation/radius;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_create_shell: number=" << number << " separation=" <<
			separation << " radius=" << radius << " theta_step=" << theta_step << endl;

	for ( i=0, theta = 0; theta <= M_PI; theta += theta_step ) {
		z = radius*cos(theta);
		rad_sin_theta = radius*sin(theta);
		n = (long) (TWOPI*rad_sin_theta/separation);
		phi_step = 1.5*TWOPI;
		if ( theta > 0 && theta < TWOPI )
			phi_step = TWOPI*1.0L/n;
		for ( phi = 0; phi < TWOPI - phi_step/2; phi += phi_step ) {
			x = cos(phi)*rad_sin_theta;
			y = sin(phi)*rad_sin_theta;
//			id = to_string(++i);
//			comp = component_add(&comp, id);
//			if ( !model->comp ) model->comp = comp;
			if ( comp ) comp = comp->add(++i);
			else comp = model->add_component(++i);
			comp->type(ct);
			comp->location(Vector3<double>(x, y, z));
			comp->view(View2<float>(x, y, z, 0));
			comp->view().normalize();
		}
	}
	number = i;
	
	double			cr(radius/sqrt(number));
	
	for ( comp = model->comp; comp; comp = comp->next )
		comp->radius(cr);
	
	if ( verbose ) {
		cout << "Generating a spherical shell:" << endl;
		cout << "Number of components:           " << number << endl;
		cout << "Separation distance:            " << separation << endl;
		cout << "Radius:                         " << radius << endl << endl;
	}

	return model;
}

/**
@brief 	Creates a Fibonacci sphere point model.
@param 	number			number of points.
@param 	radius			sphere radius.
@return Bmodel*			new sphere model.
**/
Bmodel*		model_create_fibonacci_sphere(long number, double radius)
{
	if ( radius <= 0 ) radius = 1;
	if ( number <= 0 ) number = 20;
	
	string			id("FibonacciSphere");
	Bmodel*			model = new Bmodel(id);
	
	Bcomponent*		comp = NULL;
	string			comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
	
	long			i;
	double			x, y, z, st, phi;
	double			ga(4*M_PI/(sqrt(5.0)+1.0));
	double			cr(radius/sqrt(number));
	RGBA<float>		color(0,0,1,1);

	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_create_fibonacci_sphere: number=" << number << 
			" radius=" << radius << endl;

	for ( i=0; i<number; ) {
		z = 1 - (2.0*i + 1.0)/number;
		st = sin(acos(z));
		phi = ga*i;
		x = cos(phi)*st;
		y = sin(phi)*st;
//		id = to_string(++i);
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(++i);
		else comp = model->add_component(++i);
		comp->type(ct);
		comp->location(Vector3<double>(x, y, z));
		comp->scale(radius);
		comp->view(View2<float>(x, y, z, 0));
		comp->view().normalize();
		comp->radius(cr);
		comp->color(color);
	}
	
	if ( verbose ) {
		cout << "Generating a Fibonacci sphere:" << endl;
		cout << "Number of components:           " << i << endl;
		cout << "Radius:                         " << radius << endl << endl;
	}

	return model;
}

long		model_subdivide(Bmodel* model, long division)
{
	long			i, j, n, ncomp;
	Bcomponent*		comp, *comp2, *comp3;
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) ;
	
	if ( division < 1 ) return ncomp;
	
	double			search_radius = 1/pow(2.0,division - 1.0);
	Bstring			id, comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
	
	if ( model->symmetry() == "I-3" ) search_radius *= 1.44;
	else if ( model->symmetry() == "I-5" ) search_radius *= 1.104;
	
	for ( i=0, n=ncomp, comp = model->comp; comp && i<ncomp; comp = comp->next, i++ ) {
		for ( j=i, comp2 = comp; comp2 && j<ncomp; comp2=comp2->next, j++ ) {
			if ( comp != comp2 ) {
				if ( comp->location().distance(comp2->location()) < search_radius ) {
//					id = Bstring(++n, "%ld");
//					comp3 = component_add(&comp2, id);
					comp3 = comp2->add(++n);
					comp3->type(ct);
					comp3->location((comp->location() + comp2->location())/2);
					comp3->view(View2<float>(comp3->location()));
				}
			}
		}
	}
	
	return n;
}

long		model_spherize(Bmodel* model, double sphere_fraction)
{
	if ( fabs(sphere_fraction) < 1e-3 ) return 0;
	
	double			v;
	Bcomponent*		comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		v = 1/(1-sphere_fraction*(1-comp->location().length()));
		comp->scale(v);
	}
	
	return 0;
}

/**
@brief 	Creates an dodecahedron point model.
@param 	radius			sphere radius.
@param 	divisions		number of divisions from a base dodecahedron.
@param 	sphere_fraction	spherical fraction: 0=dodecahedral, 1=spherical.
@return Bmodel*			new sphere model.

**/
Bmodel*		model_create_dodecahedron(double radius, long divisions, double sphere_fraction)
{
	Bsymmetry		sym("I-3");
	Bsymmetry		sym2("I-5");
	
	Bmodel*			model = model_platonic(sym, 1);
	if ( divisions > 1 ) model->next = model_platonic(sym2, 0.8);
	
	if ( divisions < 1 ) divisions = 1;
	if ( divisions > 4 ) {
		cout << "Warning: Divisions bigger than 4 cannot be tesselated!" << endl;
		cout << "	Divisions reset to 4" << endl << endl;
		divisions = 4;
	}
	
	if ( verbose ) {
		cout << "Generating a dodecahedral sphere:" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Divisions:                      " << divisions << endl << endl;
	}

	long			i, n;
	Vector3<double>	origin;
	Vector3<double>	scale(radius, radius, radius);

	
	for ( i=0; i<divisions; i++ ) {
		if ( i == 1 )
			n = model_merge(model);
		else
			n = model_subdivide(model, i);
		if ( verbose )
			cout << "Division " << i+1 << ": " << n << " components" << endl;
	}
	if ( verbose )
		cout << endl;
	
	model_spherize(model, sphere_fraction);
	
	model_scale(model, scale, origin);

	model_set_component_radius(model, radius/(4*divisions));

	return model;
}

/**
@brief 	Creates an icosahedron point model.
@param 	radius			sphere radius.
@param 	divisions		number of divisions from a base icosahedron.
@param 	sphere_fraction	spherical fraction: 0=dodecahedral, 1=spherical.
@return Bmodel*			new sphere model.

**/
Bmodel*		model_create_icosahedron(double radius, long divisions, double sphere_fraction)
{
	Bsymmetry		sym("I-5");
	Bmodel*			model = model_platonic(sym, 1);
	
	if ( divisions < 1 ) divisions = 1;
	if ( divisions > 4 ) {
		cout << "Warning: Divisions bigger than 4 cannot be tesselated!" << endl;
		cout << "	Divisions reset to 4" << endl << endl;
		divisions = 4;
	}
	
	if ( verbose ) {
		cout << "Generating an icosahedral sphere:" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Divisions:                      " << divisions << endl << endl;
	}

	long			i, n;
	Vector3<double>	origin;
	Vector3<double>	scale(radius, radius, radius);
	
	for ( i=0; i<divisions; i++ ) {
		n = model_subdivide(model, i);
		if ( verbose )
			cout << "Division " << i+1 << ": " << n << " components" << endl;
	}
	if ( verbose )
		cout << endl;
	
	model_spherize(model, sphere_fraction);
		
	model_scale(model, scale, origin);

	model_set_component_radius(model, radius/(4*divisions));

	return model;
}

/**
@author Jonathan Luo
@brief 	Calculates how far the slice separation is from ideal
@param 	axes			ellipsoid axis lengths
@param 	z1				previous slice z value
@param 	z2				next slice z value
@param 	distance		target separation distance
@return	double			the distance error (2.0*distance is target sum for the 2 distances)

	For the generalized ellipsoid, the distance between points of 2 slices varies
	if the x and y axes are not equal. The best slice will have a component some
	distance+error separated at the shorter axis and distance-error at the longer axis,
	distributing the error. Solve the ellipsoid equation for those two points and distances. 

**/
double		ellipsoid_slice_distance_error(Vector3<double> axes, double z1, double z2, double distance)
{
	double		temp1 = z1/axes[2], temp2 = z2/axes[2];
	double		a1 = sqrt(1 - temp1*temp1);
	double		a2 = sqrt(1 - temp2*temp2);
	
	return Vector3<double>(a1*axes[0],0,z1).distance(Vector3<double>(a2 * axes[0],0,z2))
			+ Vector3<double>(0,a1*axes[1],z1).distance(Vector3<double>(0,a2*axes[1],z2)) - 2.0 * distance;
}

/**
@author Jonathan Luo
@brief 	Calculates the next slice's z component to evenly space slices
@param 	axes			ellipsoid axis lengths
@param 	z1				previous slice's z
@param 	distance		separation distance
@return double			next slice's z-coordinate
**/
double 		ellipsoid_next_z(Vector3<double> axes, double z1, double distance)
{
	double 		z2(z1 - distance);
	if ( z2 < -axes[2] ) z2 = -axes[2];

	double 		epsilon(0.0001); // any smaller and it never exits while loop below!! //z precision
	double 		slice_error(0.1 * distance); //margin of error for the distance between points between slices
	
	if ( ellipsoid_slice_distance_error(axes, z1, z2, distance) < epsilon * distance )
		return z2; // z2 does not bracket the root with z1, it is the root.
	
	double		z3 = z1;
	double		midpoint(0);
	
	while (fabs(z3 - z2) > 2 * epsilon
			|| fabs(ellipsoid_slice_distance_error(axes,z1,midpoint,distance)) > slice_error) {
		midpoint = (z3 + z2) / 2;
		if ( ellipsoid_slice_distance_error(axes, z1, z2, distance)
				* ellipsoid_slice_distance_error(axes, z1, midpoint, distance) < 0 )
			z3 = midpoint;
		else if ( ellipsoid_slice_distance_error(axes, z1, z3, distance)
				 * ellipsoid_slice_distance_error(axes, z1, midpoint, distance) < 0 )
			z2 = midpoint;
		else break;
	}
	
	return midpoint;
}


/**
@author Jonathan Luo
@brief 	Calculates the next slice's x component to evenly space points
@param 	axes			the ellipse axis lengths and z-height
@param 	x1				???
@param 	distance		separation distance
@return double			next x-coordinate
	
	Solves the 4th order polynomial by bracketing.
**/
double		ellipse_next_x(Vector3<double> axes, double x1, double distance)
{
	double			distance2 = distance*distance;
	double			temp = axes[1] * x1 / axes[0];
	double			y1 = sqrt(axes[1]*axes[1] - temp*temp);
//	if (isnan(y1)) y1 = 0;	
	if ( !isfinite(y1) ) y1 = 0;	
	Vector3<double>	old(x1,y1,0), mid;
	
	if ( Vector3<double>(-axes[0],0,0).distance2(old) < distance2 ) return x1;

	double			right = x1;
	double			left = fmax(-axes[0], x1 - distance);
	temp = axes[1] * left / axes[0];
	double			f_left = sqrt(axes[1]*axes[1] - temp*temp);
	double			d2;
	
	// exit conditions (precision in x and distance)
	double			epsilon(0.0001);
	double			maxerror(0.1 * distance);
	
	if ( fabs(Vector3<double>(left,f_left,0).distance(old) - distance) < epsilon * distance )
		return left; //if left doesn't bracket the root, it is the root

	while ( fabs((right - left) / right) > 2 * epsilon 
			|| fabs(mid.distance(old) - distance) > maxerror ) {
		mid[0] = (right + left) / 2.0;
		temp = axes[1] * mid[0] / axes[0];
		mid[1] = sqrt(axes[1]*axes[1] - temp*temp);		
		d2 = mid.distance2(old);
		if ( d2 > distance2 ) {
			left = mid[0];
		} else if ( d2 < distance2 ) {
			right = mid[0];
		} else break;
	}
	
	return mid[0];
}

/**
@author Jonathan Luo
@brief 	Creates a circle with the given radius, height, and component separation.
@param 	radius			radius length.
@param 	z				height.
@param 	distance		separation distance.
@return Bmodel*			new circle model.

	Find angle associated with separation distance using Law of Cosines.

**/
Bmodel*		model_create_circle(double radius, double z, double distance)
{
	if ( radius <= 0 ) {
		cerr << "Error: The circle radius must given!" << endl;
		return NULL;
	}
	
	if ( distance <= 0 ) distance = radius/10;
	
	if ( verbose & VERB_FULL )
		cout << "Creating a circle of radius " << radius << " at z=" << z << 
			" with separation " << distance << endl << endl;
	
	string			id("Circle");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;

	string			comptype("CIR");
	Bcomptype*		ct = model->add_type(comptype);

	long			i, n(0);
	double			temp = distance/radius;
	double			angle = acos(1 - temp*temp/2);
	double			angle2;
	long			numpts = (long) (TWOPI / angle + 0.5);		
	angle = TWOPI / numpts;
	
	if ( radius * 2 < distance ) { //if radius is too small, just place one point at the center
		comp = model->add_component(1);
		comp->location(Vector3<double>(0,0,z));
	} else {
		for( i = 0; i < numpts; i++ ){
			angle2 = angle * i;
//			id = to_string(++n);
//			comp = component_add(&model->comp, id);
			if ( comp ) comp = comp->add(++n);
			else comp = model->add_component(++n);
			comp->location(Vector3<double>(radius * cos(angle2), radius * sin(angle2), z));
			comp->type(ct);
		}
	}
	
	return model;
}

/**
@author Jonathan Luo
@brief 	Creates an ellipse with given semiaxes and component separation
@param 	axes			ellipse axis x and y lengths and z height
@param 	distance		separation distance
@return Bmodel*			new ellipse model.
**/
Bmodel*		model_create_ellipse(Vector3<double> axes, double distance)
{
	if (fabs((axes[0] - axes[1]) / axes[1]) < 0.001)		// when eccentricity ~ 0, make circle
		return model_create_circle((axes[0] + axes[1])/2, axes[2], distance);

	//error margins, to be parameterized possibly
	double			derror(0.1 * distance);	//error for component separation distance
	double			lastpair(0.85 * distance); //accuracy for an added pair of coordinates' separation
	double			enderror(0.85 * distance); //accuracy for an endpoint distance to be added
	
	long			n = 1;
	double			temp;
	Vector3<double>	loc_old(axes[0], 0, axes[2]), loc_new(0, 0, axes[2]);
	
	string			id("Ellipse");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = model->add_component(1);
	string			comptype("ELL");
	Bcomptype*		ct = model->add_type(comptype);
	comp->type(ct);

	if ( axes[0] * 2 < distance && axes[1] * 2 < distance ) {	//if too small, just put 1 point in middle
		comp->location(loc_new);
	} else {
		comp->location(loc_old);
		while ( loc_old[0] > -axes[0] ) {
			loc_new[0] = ellipse_next_x(axes, loc_old[0], distance);
			temp = axes[1] * loc_new[0] / axes[0];
			loc_new[1] = sqrt(axes[1]*axes[1] - temp*temp);						
			if (fabs(loc_new.distance(loc_old) - distance) < derror //top point is distance from prev top point
					&& 2 * loc_new[1] > lastpair) {	//2 points upper and lower too close together
//				id = to_string(++n);
//				comp = component_add(&comp, id);
				comp = comp->add(++n);
				comp->type(ct);
				comp->location(loc_new);
				comp->location()[1] = -comp->location()[1];
//				id = Bstring(++n, "%ld");
//				comp = component_add(&comp, id);
				comp = comp->add(++n);
				comp->type(ct);
				comp->location(loc_new);
				loc_old = loc_new;
			} else {
				loc_new = Vector3<double>(-axes[0], 0, axes[2]);
				if ( loc_new.distance(loc_old) > enderror ) { // if not too bad, put a point at the end
//					id = to_string(++n);
//					comp = component_add(&comp, id);
					comp = comp->add(++n);
					comp->type(ct);
					comp->location(loc_new);
					break;
				} else break; //don't add endpoint, exit
			}
		}
	}
	
	return model;
}

/**
@author Jonathan Luo
@brief 	Creates an ellipsoid with given semiaxes and component separation
@param 	axes			ellipsoid axis lengths
@param 	distance		separation distance
@return Bmodel*			new ellipsoid model.
**/
Bmodel*		model_create_ellipsoid(Vector3<double> axes, double distance)
{
	long		swap_yx(0), swap_zx(0), swap_zy(0), swap_slice(0);
	double		distance2 = distance*distance;

	if ( verbose ) {
		cout << "Generating an ellipsoid:" << endl;
		cout << "Axes:                         " << axes << endl;
		cout << "Separation:                   " << distance << endl;
	}
	
	// arrange so semi_x > semi_y > semi_z
	if (axes[1] > axes[0]) {axes = Vector3<double>(axes[1],axes[0],axes[2]); swap_yx = 1;}
	if (axes[2] > axes[0]) {axes = Vector3<double>(axes[2],axes[1],axes[0]); swap_zx = 1;}
	if (axes[2] > axes[1]) {axes = Vector3<double>(axes[0],axes[2],axes[1]); swap_zy = 1;}
	
	// compare x-y and y-z eccentricities to set slice axis to z
	if (axes[1] * axes[1] < axes[0] * axes[2]) { //simplified ecc. comparison
		axes = Vector3<double>(axes[1],axes[2],axes[0]);
		swap_slice = 1;
	}

	
	long			nellipse(0);
	double			zend = -axes[2], dz = 2*axes[2], temp(0);
	Vector3<double>	ell_axis(0,0,axes[2]);

	string			id("Ellipsoid");
	Bmodel*			model = new Bmodel(id);
	Bmodel*			current = model;
//	Bcomponent*		comp = component_add(&model->comp, 1);
	Bcomponent*		comp = model->add_component(1);
	comp->location()[2] = -zend;
	
	while ( ell_axis[1]*ell_axis[1] + dz*dz > distance2 ) {		//y is shorter axis than x, use y to test
		ell_axis[2] = ellipsoid_next_z(axes, ell_axis[2], distance);
		temp = ell_axis[2]/axes[2];
		temp = sqrt(1 - temp*temp);
		ell_axis[0] = temp * axes[0];
		ell_axis[1] = temp * axes[1];
		current->next = model_create_ellipse(ell_axis, distance);
		current = current->next;
		dz = ell_axis[2] - zend;
		nellipse++;
	}

	if ( verbose ) cout << "Ellipses Generated:           " << nellipse << endl << endl;
	
	id = "Endpoint";
	double			derror(0.9*distance2); //allowed error in separation distance for an endpoint
	
	if ( ( ell_axis[0]*ell_axis[0] + dz*dz > derror ) && ( ell_axis[1]*ell_axis[1] + dz*dz > derror ) ) {
		current->next = new Bmodel(id);
//		comp = component_add(&current->next->comp, 1);
		comp = current->next->add_component(1);
		comp->location()[2] = zend;
	}
	
	model_merge(model);
	
	string			comptype("VER");
	Bcomptype*		ct = model->add_type(comptype);
	for ( comp = model->comp; comp; comp = comp->next) {
		comp->type(ct);
		comp->radius(distance/3);
	}
	
	Vector3<double>	x(1,0,0),y(0,1,0),z(0,0,1);
	Matrix3			mat(1);
	if(swap_slice) { mat = Matrix3(y,x) * mat; mat = Matrix3(x,z) * mat;}
	if(swap_zy) mat = Matrix3(z,y) * mat;
	if(swap_zx) mat = Matrix3(z,x) * mat;
	if(swap_yx) mat = Matrix3(y,x) * mat;
	
	if ( verbose )
		cout << mat << endl;

	model_rotate(model, mat);
	
	return model;
}

/**
@author Jonathan Luo
@brief 	Creates a cylinder with the given radius, length, and component separation.
@param 	direction		cylinder axis direction
@param 	radius			radius length.
@param 	length			cylinder length.
@param 	distance		separation distance.
@return Bmodel*			new cylinder model.

	Find angle associated with separation distance using Law of Cosines.

**/
Bmodel*		model_create_cylinder(Vector3<double> direction, double radius,
				double length, double distance)
{
	if ( radius <= 0 ) {
		cerr << "Error: The cylinder radius must given!" << endl;
		return NULL;
	}
	
	if ( length <= 0 ) {
		cerr << "Error: The cylinder length must given!" << endl;
		return NULL;
	}
	
	if ( distance <= 0 || distance > radius ) distance = radius/10;
	
	if ( direction.length() < 0.5 )
		direction = Vector3<double>(0,0,1);
		
	direction.normalize();
	
	if ( verbose & VERB_FULL )
		cout << "Creating a cylinder in direction " << direction << ", radius "
			<< radius << " and length " << length <<
			" with separation " << distance << endl << endl;
	
	string			id("Cylinder");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;

	string			comptype("CYL");
	Bcomptype*		ct = model->add_type(comptype);

	long			i, t(0), n(0);
	double			a, z, dz(distance/sqrt(2.0));
	double			temp(distance/radius);
	double			da(acos(1 - temp*temp/2));
	long			numpts = (long) (TWOPI / da + 0.5);
	da = TWOPI / numpts;
	
	for ( z = -length/2; z <= length/2; z += dz ) {
		a = ( t )? da/2: 0;
		for( i = 0; i < numpts; i++, a += da ){
//			id = to_string(++n);
//			comp = component_add(&model->comp, id);
			if ( comp ) comp = comp->add(++n);
			else comp = model->add_component(++n);
			comp->location(Vector3<double>(radius * cos(a), radius * sin(a), z));
			comp->type(ct);
		}
		t = 1 - t;
	}

	Vector3<double>	ref(0,0,1);
	Matrix3			mat = Matrix3(ref, direction);
	model_rotate(model, mat);
	
	return model;
}

/**
@brief 	Creates a plane with the given length, width, height, and component separation.
@param 	length			dimension in x.
@param 	width			dimension in y.
@param 	z				height.
@param 	distance		separation distance.
@return Bmodel*			new plane model.

	The components are placed in a hexagonal arrangement to simulate the
	most commonly expected close packingx.

**/
Bmodel*		model_create_plane(double length, double width, double z, double distance)
{
	if ( length <= 0 || width <= 0 ) {
		cerr << "Error: The plane length and width must given!" << endl;
		return NULL;
	}
	
	if ( distance <= 0 ) distance = length/10;
	
	if ( verbose & VERB_FULL ) {
		cout << "Creating a plane:" << endl;
		cout << "Size:                           " << length << " x " << width << endl;
		cout << "Height:                         " << z << endl;
		cout << "Separation:                     " << distance << endl;
	}
	
	string			id("Plane");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;

	string			comptype("PLN");
	Bcomptype*		ct = model->add_type(comptype);

	long			t(0), n(0);
	double			x, y, hl(length/2), hw(width/2), radius(distance/3);
	double			dx = distance;
	double			dy = distance/sqrt(2.0);
	
	for ( y = -hw; y <= hw; y += dy ) {
		t = 1 - t;
		for ( x = (t)? -hl: -hl+dx/2; x < hl; x += dx ) {
//			id = to_string(++n);
//			comp = component_add(&model->comp, id);
			if ( comp ) comp = comp->add(++n);
			else comp = model->add_component(++n);
			comp->location(Vector3<double>(x, y, z));
			comp->type(ct);
			comp->view(View2<float>(0,0,1,0));
			comp->radius(radius);
		}
	}
	
	if ( verbose & VERB_FULL )
		cout << "Number of components created:   " << n << endl << endl;
	
	return model;
}

/**
@brief 	Creates a spindle packed inside a sphere
@param 	direction		spindle axis direction
@param 	radius			sphere radius.
@param 	separation		separation distance between successive components.
@param 	packing			distance between spindle strands.
@return Bmodel*			new spindle model.

The spindle starts at the minimum Z-pole of the sphere, winding up along the
spherical wall to the other pole, and back down inside the first shell.
The process is repeated until the whole sphere is packed.

**/
Bmodel*		model_create_spindle(Vector3<double> direction, double radius,
				double separation, double packing)
{
	if ( radius <= 0 ) {
		cerr << "Error: The sphere radius must given!" << endl;
		return NULL;
	}
	
	if ( separation <= 0 ) {
		cerr << "Error: The component separation distance must given!" << endl;
		return NULL;
	}
	
	if ( packing <= 0 ) {
		cerr << "Error: The packing distance must given!" << endl;
		return NULL;
	}
	
	if ( direction.length() < 0.5 )
		direction = Vector3<double>(0,0,1);
		
	direction.normalize();
	
	if ( verbose & VERB_FULL )
		cout << "Creating a spindle in direction " << direction << ", radius "
			<< radius << " with separation " << separation <<
			" and packing " << packing << endl << endl;
	
	string			id("Spindle");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;

	string			comptype("SPN");
	Bcomptype*		ct = model->add_type(comptype);

	long			n(0);
	double			rp, da(0);
	double			a, b, amax;
	Vector3<double>	loc;
	
	for ( ; radius > packing; radius -= packing ) {
		amax = 3*M_PI*M_PI*radius/(2*packing);
		for ( a=0; a<amax; a+=da ) {
			b = M_PI*(a/amax - 0.5);
			loc = Vector3<double>(cos(a) * cos(b), sin(a) * cos(b), sin(b)) * radius;
//			id = to_string(++n);
//			comp = component_add(&model->comp, id);
			if ( comp ) comp = comp->add(++n);
			else comp = model->add_component(++n);
			comp->location(loc);
			comp->type(ct);
			cout << loc << endl;
			rp = sqrt(loc[0]*loc[0] + loc[1]*loc[1]);
			if ( rp > 0.1 ) da = separation/rp;
			else da = 1;
		}
	}

	Vector3<double>	ref(0,0,1);
	Matrix3			mat = Matrix3(ref, direction);
	model_rotate(model, mat);
	
	return model;
}


/**
@brief 	Creates a cubic lattice with the given length, width, height, and component separation.
@param 	lattice			lattice dimensions.
@param 	separation		distances between components.
@return Bmodel*			new lattice model.

	The components are placed in a cubic arrangement.

**/
Bmodel*		model_create_cubic_lattice(Vector3<long> lattice, double separation)
{
	if ( lattice.volume() <= 0 || separation <= 0 ) {
		cerr << "Error: The lattice size and step size must given!" << endl;
		return NULL;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Creating a cubic lattice:" << endl;
		cout << "Size:                           " << lattice << endl;
		cout << "Separation:                     " << separation << endl;
	}
	
	string			id("CubicLattice");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	string			comptype("CEL");
	Bcomptype*		ct = model->add_type(comptype);

	long			n(0), ix, iy, iz;
	double			x, y, z;
	double			radius(separation/5);

	for ( iz = 0, z = 0; iz < lattice[2]; ++iz, z += separation ) {
		for ( iy = 0, y = 0; iy < lattice[1]; ++iy, y += separation ) {
			for ( ix = 0, x = 0; ix < lattice[0]; ++ix, x += separation ) {
				if ( comp ) comp = comp->add(++n);
				else comp = model->add_component(++n);
				comp->location(Vector3<float>(x, y, z));
				comp->type(ct);
				comp->radius(radius);
			}
		}
	}

	Vector3<double>	c;
	
	if ( verbose & VERB_FULL )
		cout << "Number of components created:   " << n << endl << endl;
	
	return model;
}

/**
@brief 	Creates a hexagonal lattice with the given length, width, height, and component separation.
@param 	lattice			lattice dimensions.
@param 	separation		distances between components.
@return Bmodel*			new lattice model.

	The components are placed in a hexagonal arrangement to simulate the
	most commonly expected close packing.

**/
Bmodel*		model_create_hexagonal_lattice(Vector3<long> lattice, double separation)
{
	if ( lattice.volume() <= 0 || separation <= 0 ) {
		cerr << "Error: The lattice size and step size must given!" << endl;
		return NULL;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Creating a hexagonal lattice:" << endl;
		cout << "Size:                           " << lattice << endl;
		cout << "Separation:                     " << separation << endl;
	}
	
	string			id("HexagonalLattice");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	string			comptype("CEL");
	Bcomptype*		ct = model->add_type(comptype);

	long			n(0), ix, iy, iz;
	double			x, y, z;
	double			dx(separation);
	double			dy(separation*sqrt(3.0)/2);
	double			dz(separation*sqrt(2.0/3.0));
	double			radius(separation/5);
	
	for ( iz = 0, z = 0; iz < lattice[2]; ++iz, z += dz ) {
		for ( iy = 0, y = dy*(iz%3)/3.0; iy < lattice[1]; ++iy, y += dy ) {
			for ( ix = 0, x = 0.5*dx*((iy+iz+iz/3)%2); ix < lattice[0]; ++ix, x += dx ) {
				if ( comp ) comp = comp->add(++n);
				else comp = model->add_component(++n);
				comp->location(Vector3<float>(x, y, z));
				comp->type(ct);
				comp->radius(radius);
			}
		}
	}

	if ( verbose & VERB_FULL )
		cout << "Number of components created:   " << n << endl << endl;
	
	return model;
}
