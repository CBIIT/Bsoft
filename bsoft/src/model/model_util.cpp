/**
@file	model_util.cpp
@brief	Library routines used for model processing
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20220419
**/

#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_links.h"
#include "mol_transform.h"
#include "mol_compare.h"
#include "mol_util.h"
#include "symmetry.h"
#include "matrix_linear.h"
#include "matrix_util.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	modfunc			function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, long (modfunc)(Bmodel*))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		n += (modfunc)(mp);
	
	if ( verbose & VERB_PROCESS )
		cout << "Elements processed:             " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	i				an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, long i, long (modfunc)(Bmodel*, long))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		n += (modfunc)(mp, i);
	
	if ( verbose & VERB_PROCESS )
		cout << "Elements processed:             " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	d				an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, double d, long (modfunc)(Bmodel*, double))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		n += (modfunc)(mp, d);
	
	if ( verbose & VERB_PROCESS )
		cout << "Elements processed:             " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	&str			an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, Bstring& str, long (modfunc)(Bmodel*, Bstring& str))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		n += (modfunc)(mp, str);
	
	if ( verbose & VERB_PROCESS )
		cout << "Elements processed:             " << n << endl << endl;
	
	return n;
}

/**
@brief 	Process a list of models using the specified function.
@param 	*model			list of models.
@param 	&str			an argument.
@fn 	(modfunc)(Bmodel*)	function to be called.
@return long			aggregate number returned by function.

**/
long		models_process(Bmodel* model, string& str, long (modfunc)(Bmodel*, string& str))
{
	long 			n(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		n += (modfunc)(mp, str);
	
	if ( verbose & VERB_PROCESS )
		cout << "Elements processed:             " << n << endl << endl;
	
	return n;
}

/**
@brief 	Lists models in table form.
@param 	*model		model parameters.
@return long			number of models.
**/
long		model_list(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			nmod(0), ncomp, nct(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bstring			type_id, sym, fmap;

	cout << "Model\tType\tNcomp\tPG\tHand\tFOM\tSelect\tMap\tNumber" << endl;
	for ( mp = model; mp; mp = mp->next, nmod++ ) {
		for ( ncomp=0, comp = mp->comp; comp; comp = comp->next ) ncomp++;
		nct += ncomp;
		if ( mp->model_type().length() ) type_id = mp->model_type();
		else type_id = "?";
		if ( mp->symmetry().length() ) sym = mp->symmetry();
		else sym = "?";
		if ( mp->mapfile().length() ) fmap = mp->mapfile();
		else fmap = "?";
		cout << mp->identifier() << tab << type_id << tab << ncomp << tab << sym << tab << mp->handedness() << tab <<
			mp->FOM() << tab << mp->select() << tab << fmap << tab << mp->image_number() << endl;
	}
	cout << "Total\t\t" << nct << endl << endl;
	
	return nmod;
}

/**
@brief 	Lists models with component counts in table form.
@param 	*model			model parameters.
@return long			number of models.
**/
long		model_list_comp(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			i, nmod(0), nct(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		ct;
	Bstring			type_id, sym, fmap;
	Bstring*		comptypelist = NULL;
	Bstring*		comptype = NULL;
	Bstring*		type = NULL;
	
	for ( mp = model; mp; mp = mp->next, nmod++ ) if ( mp->select() ) {
		for ( ct = mp->type; ct; ct = ct->next ) {
			for ( type = comptypelist; type && *type != ct->identifier(); type = type->next ) ;
			if ( !type ) {
				comptype = string_add(&comptype, ct->identifier().c_str());
				if ( !comptypelist ) comptypelist = comptype;
				nct++;
			}
		}
	}
    
	int*			n = new int[nct];
	int*			ntot = new int[nct];
	for ( i=0; i<nct; i++ ) ntot[i] = 0;

	cout << "Model";
	for ( comptype = comptypelist; comptype; comptype = comptype->next )
		cout << tab << *comptype;
	cout << endl;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( i=0; i<nct; i++ ) n[i] = 0;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			for ( i=0, comptype = comptypelist; comptype && *comptype != comp->type()->identifier(); comptype = comptype->next, i++ ) ;
			if ( comptype ) n[i]++;
		}
		cout << mp->identifier();
		for ( i=0, comptype = comptypelist; comptype; comptype = comptype->next, i++ ) {
			ntot[i] += n[i];
			cout << tab << n[i];
		}
		cout << endl;
		nmod++;
	}

	cout << "Total";
	for ( i=0, comptype = comptypelist; comptype; comptype = comptype->next, i++ )
		cout << tab << ntot[i];
	cout << endl << endl;

	delete[] n;
	delete[] ntot;
	
	return nmod;
}

/**
@brief	Replaces all components in a model with reference components.
@param 	*model			model structure to be modified.
@param 	*modref			reference model.
@return long			number of components used in reference.

	Only the first model in the reference is used.
	All models will have identical sets of components.

**/
long		model_replace_components(Bmodel* model, Bmodel* modref)
{
	int				i;
	Bmodel*			mp;
	Bcomptype*		ct;
	Bcomptype*		ctref;
	Bcomponent*		comp;
	Bcomponent*		compref;
	Blink*			link;
	Blink*			linkref;
	Bpolygon*		poly;
	Bpolygon*		polyref;
	
	if ( verbose & VERB_FULL )
		cout << "Replacing components in a model from a reference model: " << modref->identifier() << endl;

	long	n = count_list((char *) modref->comp);
	
	for ( mp = model; mp; mp = mp->next ) {
		model_link_list_kill(mp);
		comp_type_list_kill(mp->type);
		component_list_kill(mp->comp);
		poly_list_kill(mp->poly);
		mp->type = ct = NULL;
		mp->comp = comp = NULL;
		mp->link = link = NULL;
		mp->poly = poly = NULL;
		for ( ctref = modref->type; ctref; ctref = ctref->next ) {
//			ct = (Bcomptype *) add_item((char **) &ct, sizeof(Bcomptype));
//			if ( !mp->type ) mp->type = ct;
//			component_type_copy(ctref, ct);
			if ( mp->type ) mp->type->add(ctref);
			else mp->type = new Bcomptype(ctref);
		}
		for ( compref = modref->comp; compref; compref = compref->next ) {
			comp = mp->add_component(compref);
			comp->type(mp->add_type(compref->type()->identifier()));
		}
		for ( linkref = modref->link; linkref; linkref = linkref->next ) {
			link = (Blink *) add_item((char **) &link, sizeof(Blink));
			if ( !modref->link ) modref->link = link;
			for ( compref = modref->comp, comp = mp->comp; compref; compref = compref->next, comp = comp->next ) {
				if ( linkref->comp[0] == compref ) link->comp[0] = comp;
				if ( linkref->comp[1] == compref ) link->comp[1] = comp;
			}
			link->length(linkref->length());
			link->radius(linkref->radius());
			link->color(linkref->color());
		}
		for ( polyref = modref->poly; polyref; polyref = polyref->next ) {
			poly = (Bpolygon *) add_item((char **) &poly, sizeof(Bpolygon));
			if ( !modref->poly ) modref->poly = poly;
			poly->normal(polyref->normal());
			poly->closed(polyref->closed());
			for ( compref = modref->comp, comp = mp->comp; compref; compref = compref->next, comp = comp->next ) {
				for ( i=0; i<polyref->comp.size() && polyref->comp[i]; i++ )
					if ( polyref->comp[i] == compref ) poly->comp.push_back(comp);
			}
		}
	}
	
	model_setup_links(model);
	
	return n;
}

/**
@brief 	Merges components from all models into one.
@param 	*model		model parameters.
@return long		number of components.
**/
long		model_merge(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			nct, ncomp, nlink;
	Bmodel*			mp;
	Bcomptype*		ct = model->type;
	Bcomponent*		comp = model->comp;
	Blink*			link = model->link;

	if ( comp ) for ( ; comp->next; comp = comp->next ) ;
	if ( link ) for ( ; link->next; link = link->next ) ;
	
	for ( mp = model->next; mp; mp = mp->next ) {
		if ( ct ) ct->next = mp->type;
		else model->type = ct = mp->type;
		mp->type = NULL;
		if ( ct ) for ( ; ct->next; ct = ct->next ) ;
		if ( comp ) comp->next = mp->comp;
		else model->comp = comp = mp->comp;
		mp->comp = NULL;
		if ( comp ) for ( ; comp->next; comp = comp->next ) ;
		if ( link ) link->next = mp->link;
		else model->link = link = mp->link;
		mp->link = NULL;
		if ( link ) for ( ; link->next; link = link->next ) ;
	}
	
	model_kill(model->next);
	model->next = NULL;
	
	for ( nct=0, ct = model->type; ct; ct = ct->next ) nct++;

	// Renumber all the components
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next )
		comp->identifier() = to_string(++ncomp);

	for ( nlink=0, link = model->link; link; link = link->next ) nlink++;

	if ( verbose & VERB_PROCESS ) {
		cout << "Merged component types:         " << nct << endl;
		cout << "Merged components:              " << ncomp << endl;
		cout << "Merged links:                   " << nlink << endl;
		cout << endl;
	}
		
	return ncomp;
}

/**
@brief 	Add a number to each model id.
@param 	*model		model parameters.
@return long		number of models.

	The intention is to give unique id's to models.

**/
long		model_number_ids(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			n(0);
	Bmodel*			mp;
	Bstring			id;
	
	for ( mp = model; mp; mp = mp->next ) {
		id = mp->identifier().c_str() + Bstring(++n, "_%04d");
		mp->identifier() = id.str();
	}
	
	return n;
}

/**
@brief 	Rename components.
@param 	*model		model parameters.
@return long		number of components.

	The number of links to a component determines its new name.
	Only the first model is processed.

**/
long		model_rename_components(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	long			i, n(0);
	Bcomponent*		comp;

	Bstring			ctstr[10];
	ctstr[0] = "NIL";
	ctstr[1] = "MON";
	ctstr[2] = "DI";
	ctstr[3] = "TRI";
	ctstr[4] = "TET";
	ctstr[5] = "PEN";
	ctstr[6] = "HEX";
	ctstr[7] = "HEP";
	ctstr[8] = "OCT";
	ctstr[9] = "NON";
	
	for ( comp = model->comp; comp; comp = comp->next, n++ ) {
		for ( i=0; comp->link[i]; i++ ) ;
		if ( i > 9 ) i = 9;
//		comp->type = model_add_type_by_id(model, ctstr[i]);
		comp->type(model->add_type(ctstr[i]));
	}
	
	return n;
}

/**
@brief 	Calculates the mass of a model from component masses.
@param 	*model		model parameters.
@return double		model mass.

	The component type masses must be provided.
	Only the first model in the list is processed.

**/
double		model_mass(Bmodel* model)
{
	double			mass = 0;

	if ( !model ) return mass;
	if ( !model->select() ) return mass;
	
	Bcomptype*		ct = NULL;
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next ) {
		ct = comp->type();
		if ( ct ) mass += ct->mass();
	}
	
	if ( verbose & VERB_FULL )
		cout << "Model: " << model->identifier() << "  Mass = " << mass << endl;
	
	return mass;
}

/**
@brief 	Calculates the masses of all the models in the list.
@param 	*model		linked list of model parameters.
@return long		number of selected models.

	The component type masses must be provided.

**/
long		model_mass_all(Bmodel* model)
{
	long			nmod(0);
	Bmodel*			mp;
	
	cout << "Model\tMass" << endl;
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		cout << mp->identifier() << tab << model_mass(mp) << endl;
		nmod++;
	}
	cout << endl;
	
	return nmod;
}


/**
@brief 	Calculates the center-of-mass of a model.
@param 	*model			model parameters.
@return Vector3<double>	center-of-mass.

	Only the first model in the list is processed.

**/
Vector3<double>	model_center_of_mass(Bmodel* model)
{
	Vector3<double>	com;

	if ( !model ) return com;
	if ( !model->select() ) return com;
	
	int				n = 0;
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		com += comp->location();
		n++;
	}
	
	com /= n;
	
	return com;
}

/**
@brief 	Calculates the center-of-mass of a list of models.
@param 	*model			model parameters.
@return Vector3<double>	center-of-mass.

**/
Vector3<double>	models_center_of_coordinates(Bmodel* model)
{
	Vector3<double>	com;

	if ( !model ) return com;
	if ( !model->select() ) return com;
	
	long		ncomp(0);
	Bmodel* 	mp;
	Bcomponent*	comp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			com += comp->location();
			ncomp++;
		}
	}
	
	com /= ncomp;
	
	return com;
}

/**
@brief 	Calculates the sum of distances of model components to a location.
@param 	*model			model parameters.
@param 	loc				reference location.
@return double			sum of distances.

	Only the first model in the list is processed.

**/
double			model_distance_sum(Bmodel* model, Vector3<double> loc)
{
	long			n(0);
	double			ds(0);
	Bcomponent*		comp;
	
	for ( comp = model->comp; comp; comp = comp->next, ++n )
		ds += loc.distance(comp->location());

	if ( n ) ds /= n;
	
	return ds;
}

/**
@brief 	Calculates the geometric median estimate of model components.
@param 	*model			model parameters.
@param 	pgm				previous geometric median estimate.
@return Vector3<double>	new geomtric median estimate.

	Only the first model in the list is processed.

**/
Vector3<double>	model_geometric_median_estimate(Bmodel* model, Vector3<double> pgm)
{
	double			d, ds(0);
	Vector3<double>	gm;
	Bcomponent*		comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		d = 1.0/pgm.distance(comp->location());
		if ( d > 0.001 ) {
			ds += d;
			gm += comp->location() * d;
		}
	}
	
	if ( ds ) gm /= ds;

	return gm;
}

/**
@brief 	Calculates the geometric median of model components.
@param 	*model			model parameters.
@return Vector3<double>	geometric median.

	Only the first model in the list is processed.
	Based on Weiszfeld’s method.

**/
Vector3<double>	model_geometric_median(Bmodel* model)
{
	long			i, iter(1000);
	double			tol(0.01), dd, ddd(1), ds;
	Vector3<double>	gm, pgm;
	
	pgm = model_center_of_mass(model);
	ddd = dd = pgm.length();
	ds = model_distance_sum(model, gm);
	if ( verbose & VERB_PROCESS ) {
		cout << "Iter\tGMx\tGMy\tGMz\t∆d\tDistSum" << endl;
		cout << 0 << tab << pgm << tab << ddd << tab << ds << endl;
	}
	
	for ( i=0; i<iter && ddd > tol; ++i ) {
		ddd = dd;
		gm = model_geometric_median_estimate(model, pgm);
		dd = gm.distance(pgm);
		ddd -= dd;
		pgm = gm;
		ds = model_distance_sum(model, gm);
		if ( verbose & VERB_PROCESS )
			cout << i+1 << tab << gm << tab << ddd << tab << ds << endl;
	}
	
	return gm;
}

/**
@brief 	Calculates the radius of gyration for a model.
@param 	*model		model parameters.
@return double		radius of gyration.

	Only the first model in the list is processed.

**/
double		model_gyration_radius(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	int				n;
	double			d, R;
	Bcomponent*		comp;
	
	Vector3<double>	com = model_center_of_mass(model);
	
	for ( n = 0, R = 0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		d = (comp->location() - com).length();
		R += d*d;
		n++;
	}
	
	if ( n ) R = sqrt(R/n);
	
	return R;
}

/**
@brief	Calculates the principal axes of a model.
@param 	*model			model structure.
@param 	*eigenvec		eigen vectors (can be NULL).
@return Vector3<double>	3-valued vector of principal axes.

	Only the first model in the list is processed.

**/
Vector3<double> 	model_principal_axes(Bmodel* model, Vector3<double>* eigenvec)
{
	Vector3<double>	eigenval;

	if ( !model  ) return eigenval;
	if ( !model->select() ) return eigenval;
	
	double			summass(0);
	Vector3<double>	loc, vec, vec2, vecx;
	Bcomponent*		comp;

//	Vector3<double>	com = model_center_of_mass(model);
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		loc = comp->location();
		vec += loc;					// Sums
		vec2 += loc*loc;		// Square sums
		vecx[0] += loc[0]*loc[1];	// Cross-term sums
		vecx[1] += loc[0]*loc[2];
		vecx[2] += loc[1]*loc[2];
		summass += 1;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Number of components:           " << summass << endl;
	}
	
	if ( summass < 1 ) {
		cout << "Error: No components found!" << endl << endl;
		return vec;
	}
	
	vec /= summass;
	vec2 /= summass;
	vecx /= summass;
	
	return principal_axes(vec, vec2, vecx, eigenvec);
}

Vector3<double> 	model_principal_axes(Bmodel* model, Matrix& eigenvec)
{
	Vector3<double>			pax;

	if ( !model  ) return pax;
	if ( !model->select() ) return pax;

	vector<Vector3<double>>	coor;
	Vector3<double>			loc;
	Bcomponent*				comp;

	for ( comp = model->comp; comp; comp = comp->next )
		if ( comp->select() ) {
			loc = comp->location();
			coor.push_back(loc);
		}
	
	if ( verbose & VERB_PROCESS )
		cout << "Model:                          " << model->identifier() << endl;
	
	if ( coor.size() < 1 ) {
		cout << "Error: No components found for model " <<
			model->identifier() << "!" << endl << endl;
		return pax;
	}
	
	return principal_axes(coor, eigenvec);
}

long		model_principal_axes(Bmodel* model)
{
	if ( !model  ) return -1;

	Vector3<double>		pax;
	Bmodel*				mp;
	Matrix				eigenvec(3,3);
	
	if ( verbose )
		cout << "Principal axes:" << endl << "Model\tMajor\tMiddle\tMinor" << endl;
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		pax = model_principal_axes(mp, eigenvec);
		if ( verbose )
			cout << mp->identifier() << tab << pax[0] << tab << pax[1] << tab << pax[2] << endl;
	}
	
	if ( verbose )
		cout << endl;
	
	return 0;
}


/**
@brief	Calculates the radial distribution function of a model.
@param 	*model		model structure.
@param 	interval	interval between bins.
@return long		0.

	Only the first model in the list is processed.

**/
long		model_radial_distribution(Bmodel* model, double interval)
{
	if ( interval <= 0 ) interval = 1;	
	if ( !model  ) return 0;
	if ( !model->select() ) return 0;
	
	Vector3<double>	min(1e10,1e10,1e10), max(-1e10,-1e10,-1e10);
	Bcomponent*		comp1, *comp2;

	for ( comp1 = model->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		min = min.min(comp1->location());
		max = max.max(comp1->location());
	}
	
	int				r, rmax = (int) (max.distance(min)/interval + 1);
	double			d;
	int*			rdf = new int[rmax];
	for ( r=0; r<rmax; r++ ) rdf[r] = 0;
	
	for ( comp1 = model->comp; comp1->next; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			r = (int) (d/interval + 0.5);
			if ( r < rmax ) rdf[r]++;
		}
	}
	
//	if ( verbose ) {
		cout << "Radius\tCount" << endl;
		for ( r=0; r<rmax; r++ )
			cout << r*interval << tab << rdf[r] << endl;
		cout << endl;
//	}
	
	delete[] rdf;
	
	return 0;
}

long		molgroup_write_into_grid(Bmolgroup* molgroup, Vector3<int> size, Vector3<double> min, double sampling, int* grid)
{
	long			i, x, y, z;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;

	long	vol = (long) size.volume();
	char*			tgrid = new char[vol];
	for ( i=0; i<vol; i++ ) tgrid[i] = 0;
	
	for ( mol=molgroup->mol; mol; mol=mol->next ) {
		for ( res=mol->res; res; res=res->next ) {
			for ( atom=res->atom; atom; atom=atom->next ) {
				x = (int) ((atom->coord[0] - min[0])/sampling);
				y = (int) ((atom->coord[1] - min[1])/sampling);
				z = (int) ((atom->coord[2] - min[2])/sampling);
				i = (z*size[1]+y)*size[0]+x;
				if ( i < 0 || i > vol ) {
					cerr << "Error in molgroup_write_into_grid: i=" << i << " (vol=" << vol << ")" << endl;
					return -1;
				}
				tgrid[i] = 1;
			}
		}
	}
	
	for ( i=0; i<vol; i++ ) grid[i] += tgrid[i];

	delete[] tgrid;
	
	return 0;
}

/**
@brief 	Concatenates selected molecules into one group.
@param 	*model		model parameters.
@param 	&paramfile	atomic parameter file.
@param 	separate	flag to generate separate molecule groups.
@return Bmolgroup*	list of molecule groups.

	Only the first model in the linked list is processed.

**/
Bmolgroup*	model_assemble(Bmodel* model, Bstring& paramfile, int separate)
{
	Bcomponent*		comp = NULL;
//	Bcomptype*		comptype = NULL;

	Bmolgroup*		mglist = NULL;
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		molgroup1 = NULL;
	Bmolecule*		mol = NULL;
    Bstring    		atom_select("all");
	Quaternion		q;
	Transform		t;

	long			i, nsel, nover;
	double			sampling = 5;
	Vector3<double>	min(model->comp->location()), max(model->comp->location());
	Vector3<int>	size;
	
	if ( verbose )
		cout << "Assembling components" << endl;
	
	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		min = min.min(comp->location());
		max = max.max(comp->location());
		nsel++;
	}
	
	Bstring 		fn(model->type->file_name());
	molgroup1 = read_molecule(fn, atom_select, paramfile);
	min -= molgroup1->box;
	max += molgroup1->box;
	molgroup_kill(molgroup1);
	for ( i=0; i<3; i++ ) size[i] = (int) ((max[i] - min[i])/sampling);
	
	if ( nsel < 1 ) {
		cerr << "Error: No components are selected!" << endl;
		return NULL;
	}
	
	long	vol = (long) size.volume();
	int*			grid = new int[vol];
	for ( i=0; i<vol; i++ ) grid[i] = 0;
	
	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		nsel++;
//		comptype = model_get_type(model, comp->type);
		fn = comp->type()->file_name();
		molgroup1 = read_molecule(fn, atom_select, paramfile);
		molgroup1->id = comp->identifier();
//		q = quaternion_from_view(comp->view);
		q = comp->view().quaternion();
//		t = transform_from_quaternion(q);
		t = Transform(q);
		t.origin = molgroup_center_of_mass(molgroup1);
		t.trans = comp->location() - t.origin;
		molgroup_coor_rotate(molgroup1, t);
		molgroup_stats(molgroup1);
		if ( molgroup_write_into_grid(molgroup1, size, min, sampling, grid) < 0 ) {
			error_show("Error in model_assemble", __FILE__, __LINE__);
			return NULL;
		}
		if ( separate ) {
			if ( !mglist ) mglist = molgroup = molgroup1;
			else {
				molgroup->next = molgroup1;
				molgroup = molgroup1;
			}
		} else {
			if ( molgroup ) {
				if ( mol ) {
					for ( ; mol->next; mol = mol->next ) ;
					mol->next = molgroup1->mol;
				} else {
					mol = molgroup->mol = molgroup1->mol;
				}
				molgroup1->mol = NULL;
				molgroup_kill(molgroup1);
			} else {
				mglist = molgroup = molgroup1;
				mol = molgroup->mol;
			}
		}
	}

	for ( i=nover=0; i<vol; i++ ) if ( grid[i] > 1 ) nover++;
	
	delete[] grid;
		
	if ( verbose ) {
		cout << "Components assembled:           " << nsel << endl;
		cout << "Molecule group overlap:         " << 
			sampling*sampling*sampling*nover << " A3 (" << nover*100.0/size.volume() << " %)" << endl << endl;
	}
	
	return mglist;
}

/**
@brief 	Calculates the centers-of-mass of molecule group components and generates a new model.
@param 	*molgroup	list of molecule groups.
@return Bmodel*		new model.

	Each molecule is assumed to be a component.

**/
Bmodel*		model_generate_com(Bmolgroup* molgroup)
{
	Bstring			id, path;
	Bstring			comptype("VER");
	Bmolgroup*		mg;
	Bmolecule*		mol;
	
	int				i, j, n=0;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Generating a centers-of-mass model" << endl << endl;
	
	for ( i=1, mg = molgroup; mg; mg = mg->next, i++ ) {
		mp = (Bmodel *) add_item((char **) &mp, sizeof(Bmodel));
		if ( !model ) model = mp;
		if ( mg->id.length() ) mp->identifier(mg->id);
		else mp->identifier() = to_string(i);
		comp = NULL;
		for ( j=1, mol = molgroup->mol; mol; mol = mol->next, j++, n++ ) {
			cout << "Adding molecule " << j << " as component" << endl;
//			comp = component_add(&comp, j);
//			if ( !mp->comp ) mp->comp = comp;
			if ( comp ) comp = comp->add(j);
			else mp->comp = comp = new Bcomponent(j);
			comp->location(mol_center_of_mass(mol));
			if ( mol->id.length() ) id = mol->id.no_space();
			else id = comptype;
//			comp->type = model_add_type_by_id_and_filename(mp, id, molgroup->filename, 0);
			comp->type(mp->add_type(id.str(), molgroup->filename.str(), 0));
		}
	}
	
	cout << "Models generated:               " << --i << endl;
	cout << "Components generated:           " << n << endl << endl;

	model_check(model, path);
	
	return model;
}

/**
@brief 	Updates the centers-of-mass of molecule group components.
@param 	*model		model parameters.
@param 	*molgroup	list of molecule groups.
@return long		number of selected components.

	The identifiers of the molecule groups must correspond to the component identifiers.

**/
long		model_update_centers_of_mass(Bmodel* model, Bmolgroup* molgroup)
{
	long			nsel(0);
	Bcomponent*		comp = NULL;
	Bmolgroup*		mg = NULL;
	
	if ( verbose )
		cout << "Updating component centers-of-mass" << endl << endl;

	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		nsel++;
		for ( mg = molgroup; mg; mg = mg->next ) if ( comp->identifier() == mg->id.str() ) break;
		if ( mg ) comp->location(molgroup_center_of_mass(mg));
	}
	
	return nsel;
}

/**
@brief     Averages sequential components.
@param 	*model		model structure to be modified.
@param 	number		number of components to average.
@return long		number of remaining components.

	Only the first component in each set with modified coordinates is kept.

**/
long		model_average_components(Bmodel* model, int number)
{
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		comp_avg;

	if ( verbose & VERB_PROCESS )
		cout << "Averaging every " << number << " of selected components" << endl;
		
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		comp_avg = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( !comp_avg ) {
				comp_avg = comp;
				comp_avg->select(1);
			} else {
				comp->select(0);
				comp_avg->shift(comp->location());
				comp_avg->select_increment();
				if ( comp_avg->select() == number ) {
					comp_avg->scale(1.0/comp_avg->select());
					comp_avg = NULL;
				}
			}
		} else comp->select(1);
		if ( comp_avg ) comp_avg->location(comp_avg->location() / comp_avg->select());
	}
	
	return model_delete_non_selected(&model);
}

/**
@brief	Generates an array of pointers to component structures.
@param 	*model			model structure.
@param 	&n				pointer to number of comps found.
@return Bcomponent**	array of pointers to components.
**/
Bcomponent**	component_get_array(Bmodel* model, long& n)
{
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) n++;

	Bcomponent**	comparray = new Bcomponent*[n];

	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->select() )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() )
			comparray[n++] = comp;

	return comparray;
}

/**
@brief 	Calculates a plane through an array of components.
@param 	**comparray		array of components.
@param 	&offset			offset from plane.
@return Vector3<double>	plane normal.

	A plane is fit through the polygon vertices and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	The polygon planarity is defined as the root-mean-square-deviation from 
	the fitted plane.

**/
Vector3<double>	component_plane(vector<Bcomponent*>& comparray, double& offset)
{
    long     			i, j, n;
	vector<double>		b(3);
	Matrix				a(3,3);
	Vector3<double>		loc, normal, center;

	for ( i=0; i<3; i++ ) b[i] = 0;
	
	center = 0;
	for ( n=0; n<comparray.size() && comparray[n]; n++ ) {
		loc = comparray[n]->location();
		center += loc;
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<=i; j++ )
				a[3][i] += loc[i]*loc[j];
			b[i] += loc[i];
		}
	}
	center /= n;
	
	for ( i=1; i<3; i++ )
		for ( j=0; j<i; j++ )
			a[3][j] = a[3][i];
	
	offset = 0;
	if ( a[0][0] == 0 ) {
		normal = Vector3<double>(1, 0, 0);
	} else if ( a[1][1] == 0 ) {
		normal = Vector3<double>(0, 1, 0);
	} else if ( a[2][2] == 0 ) {
		normal = Vector3<double>(0, 0, 1);
	} else {
		a.LU_decomposition(b);
		normal = Vector3<double>(b[0], b[1], b[2]);
		normal.normalize();
	}

	offset = center.scalar(normal);

	return normal;
}

/**
@brief	Calculates the bounds of a list of models.
@param 	*model	 	model list.
@return vector<Vector3<double>>	minimum and maximum bounds.

**/
vector<Vector3<double>>	models_calculate_bounds(Bmodel* model)
{
	vector<Vector3<double>> bounds;
	Vector3<double>			mn(1e30,1e30,1e30), mx(-1e30,-1e30,-1e30);
	Bmodel*					mp;
	
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp->calculate_bounds();
		mn = mn.min(mp->minimum());
		mx = mx.max(mp->maximum());
	}
	
	bounds.push_back(mn);
	bounds.push_back(mx);
	return bounds;
}

/**
@brief	Generates lists of atoms based on a grid.
@param 	*model	 	model list.
@param 	size		size of grid.
@param 	origin		origin of grid.
@param 	sampling	spacing in each dimension.
@return vector<vector<Bcomponent*>>	array of component arrays.

	The goal is to fit all the components within the grid boundaries.
	Components located outside the grid will be added to the edges.

**/
vector<vector<Bcomponent*>>	model_component_grid(Bmodel* model, Vector3<long>& size,
			Vector3<double>& origin, Vector3<double>& sampling)
{
	if ( sampling.volume() < 1 ) {
		cerr << "Error in model_component_grid: sampling must be specified!" << endl;
		bexit(-1);
	}

	long			i;
	Vector3<double>	mn(1e30,1e30,1e30), mx(-1e30,-1e30,-1e30), s1, loc;
	
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( size.volume() < 1 ) {
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			mp->calculate_bounds();
			mn = mn.min(mp->minimum());
			mx = mx.max(mp->maximum());
		}
		s1 = mx - mn;
		origin = -mn;
		size = s1/sampling;
	}

	long			gridvol = (long) size.volume();
	Vector3<long>	size1 = size - 1;
	vector<vector<Bcomponent*>>	grid(gridvol);

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			loc = (comp->location() + origin)/sampling;
			loc = loc.max(0);
			loc = loc.min(size1);
			i = (loc[2]*size[1] + loc[1])*size[0] + loc[0];
			grid[i].push_back(comp);
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_component_grid: Done!" << endl;
	
	return grid;
}

/**
@brief	Generates an array of pointers to model components.
@param 	*model 				model.
@return vector<Bcomponent*>	array of pointers to components.
**/
vector<Bcomponent*>	models_get_component_array(Bmodel* model)
{
	Bmodel*				mp;
	Bcomponent*			comp;
	vector<Bcomponent*>	carr;

    for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		for( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() )
			carr.push_back(comp);

	return carr;
}

/**
@brief	Splits the components of models into slices.
@param 	*model 				model.
@param 	bottom 				minimum coordinate in z.
@param 	top 					maximum coordinate in z.
@param 	thickness 			slice thickness.
@return vector<vector<Bcomponent*>>	sets of arrays of pointers to components.
**/
vector<vector<Bcomponent*>>	model_split_into_slices(Bmodel* model, double bottom, double top, double thickness)
{
	long			i, n((top - bottom)/thickness);
	vector<vector<Bcomponent*>>	comp_slice(n);
	
	if ( verbose ) {
		cout << "Splitting model into slices:" << endl;
		cout << "Bottom:                      " << bottom << " A" << endl;
		cout << "Top:                         " << top << " A" << endl;
		cout << "Slice thickness:             " << thickness << " A" << endl;
		cout << "Number of slices:            " << n << endl;
	}
	
	Bmodel*			mp;
	Bcomponent*		comp;

	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			i = (comp->location()[2] - bottom)/thickness;
			if ( i>=0 && i<n )
				comp_slice[i].push_back(comp);
		}
	}
	
	if ( verbose ) {
		long		ncomp(0);
		cout << "Slice\tComponents" << endl;
		for ( i=0; i<n; ++i ) {
			ncomp += comp_slice[i].size();
			cout << i+1 << tab << comp_slice[i].size() << endl;
		}
		cout << "Total number of components:  " << ncomp << endl << endl;
	}
	
	return comp_slice;
}
