/**
@file	model_views.cpp
@brief	Library routines used for analyzing the views of components in models
@author Bernard Heymann
@date	Created: 20081120
@date	Modified: 20210205
**/

#include "rwimg.h"
#include "model_views.h"
#include "model_util.h"
#include "model_neighbors.h"
#include "mol_compare.h"
#include "Matrix3.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Returns a list of component views.
@param 	*model			model parameters.
@return list<View2<float>>	list of views.

	Only the first model is processed.

**/
list<View2<float>>		views_from_model(Bmodel* model)
{
	Bcomponent*			comp;
	list<View2<float>>	v;

	for ( comp = model->comp; comp; comp = comp->next )
		if ( comp->select() > 0 )
			v.push_back(comp->view());
	
	return v;
}

/**
@brief 	Returns a list of component views.
@param 	*model		model parameters.
@return View*			list of views.

	All models are processed.

**/
/*View*		views_from_models(Bmodel* model)
{
	Bmodel*			mp;
	View*			view = NULL;
	View*			v = NULL;
	View*			v2 = NULL;

	for ( mp = model; mp; mp = mp->next ) {
		v = views_from_model(mp);
		if ( !view ) view = v2 = v;
		else {
			while ( v2->next ) v2 = v2->next;
			v2->next = v;
		}
	}

	return view;
}*/

list<View2<float>>	views_from_models(Bmodel* model)
{
	Bmodel*				mp;
	Bcomponent*			comp;
	list<View2<float>>	v;

	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next )
			if ( comp->select() > 0 )
				v.push_back(comp->view());

	return v;
}

/**
@brief 	Sets views.
@param 	*model		model parameters.
@param 	view		set view.
@return long			number of selected components.

	Each component view is set tot the given view.

**/
long		model_set_views(Bmodel* model, View2<float> view)
{
	long			nsel(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( verbose )
		cout << "Setting views" << endl << endl;

	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			nsel++;
			comp->view(view);
		}
	}
	
	return nsel;
}

/**
@brief 	Invert views.
@param 	*model		model parameters.
@return long			number of selected components.

	It calculates the inverse of each component view.
	Only the first model is processed.

**/
long		model_invert_views(Bmodel* model)
{
	long			nsel(0);
	Bmodel*			mp = model;
	Bcomponent*		comp;
	
	if ( verbose )
		cout << "Inverting views" << endl << endl;

	for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		nsel++;
		comp->view(comp->view().backward());
	}
	
	return nsel;
}

/**
@brief 	Finds the molecule views with respect to a reference.
@param 	*model		model parameters.
@param 	&reffile	reference molecule file name.
@param 	&paramfile	atomic parameter file.
@return long			number of molecules selected.

	The positioning of each molecule is based on the center of mass of the reference.

**/
long		model_find_views(Bmodel* model, Bstring& reffile, Bstring& paramfile)
{
	Bcomponent*		comp = NULL;
	Transform		t;
    Bstring    		atom_select("all");
	
	if ( !model->type ) {
		cerr << "Error: No component types found!" << endl;
		return -1;
	}
	
	long			nsel(0);
	Bstring			fn;
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		ref_molgroup = read_molecule(reffile, atom_select, paramfile);
	
	if ( verbose )
		cout << "Molecule\tDist\t\t\tAxis\t\t\tAngle" << endl;
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		nsel++;
		fn = comp->type()->file_name();
		molgroup = read_molecule(fn, atom_select, paramfile);
		t = molgroup_find_transformation(molgroup, ref_molgroup);
		comp->view(View2<float>(t.angle, t.axis));
		molgroup_kill(molgroup);
		if ( verbose )
			cout << comp->type()->file_name()  << tab << comp->location() << tab << comp->view() << endl;
	}
	
	molgroup_kill(ref_molgroup);
	
	return nsel;
}

/* The normal is defined as the average of all the rotation axes associated 
with the angles between the links */
Vector3<double>	component_normal(Bcomponent* comp)
{
	Vector3<double>		normal;

	if ( !comp->link[1] ) return normal;
	
	long					i, j;
	Vector3<double>		v1, v2;
	
	for ( i=1; i<comp->link.size(); i++ ) {
		for ( j=0; j<i; j++ ) {
			v1 = comp->link[i]->location() - comp->location();
			v2 = comp->link[j]->location() - comp->location();
			v1 = v1.cross(v2);
			if ( comp->location().angle(v1) > M_PI_2 ) v1 = -v1;
			v1.normalize();
			normal += v1;
		}
	}
	if ( normal.length2() > 1e-10 ) normal.normalize();
	
	return normal;
}

/* The view vector is the normal and the rotation angle is chosen such that
the first link falls in the xz plane when the component is rotated to the 
standard view */
View2<float>	component_view(Bcomponent* comp)
{
	Vector3<double>	w(comp->view().vector3());
	
	if ( comp->link[1] ) w = component_normal(comp);
	
	if ( w.length2() < 0.1 ) w = comp->location();
	
	long			i, im;
	Vector3<double>	v, u;
	double			a(0);

	for ( i=1, im=0; i<comp->link.size() && comp->link[i]; i++ ) {
		v = comp->location().normal(comp->link[im]->location(), comp->link[i]->location());
		a = w.angle(v);
		if ( a > M_PI_2 ) im = i;
	}

	w.normalize();
	u = comp->link[im]->location() - comp->location();
	u.normalize();
	v = w.cross(u);
	v.normalize();
	u = v.cross(w);
	a = angle_set_negPI_to_PI(atan2(w[1], w[0]) + atan2(v[2], -u[2]));
	
	View2<float>	view(w[0], w[1], w[2], a);

	return view;
}

/**
@brief 	Calculates the views associated with each component.
@param 	*model		model structure (views modified).
@param 	&mode		none=current, com, map, local.
@return long			number of selected components.

	The view for a vertex is calculated from the vertex vector and the first link.
	The vertex vector points away from the origin and is calculated as follows:
	Angles defined for the vertex:
		vv = normalized sum of the cross products for each pair of links
			constituting an angle
	Angles not defined:
		vv = component coordinates - origin
	The origin is defined as the current zero coordinates, the center-of-mass,
	or from the map.

**/
long		model_calculate_views(Bmodel* model, Bstring& mode)
{
	if ( mode.length() < 1 ) {
		cerr << "Error: The view calculation mode must be specified!" << endl << endl;
		bexit(-1);
	}
	
	if ( mode[0] == 'l' ) return model_calculate_local_views(model);
	
	Bmodel*			mp;
	Bcomponent*  	comp;
	Vector3<double>	normal;	
	Vector3<double>	origin;

	long			nsel(0);
	
	if ( verbose & VERB_PROCESS )
		cout << "Calculating views based on origin-component vectors" << endl << endl;
	
    for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		if ( mode[0] == 'c' )
			origin = model_center_of_mass(model);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG model_calculate_views: origin=" << origin << endl;
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
			normal = comp->location() - origin;
			if ( verbose & VERB_DEBUG )
				cout << "Old view: " << comp->view() << endl;
			comp->view(View2<float>(normal[0], normal[1], normal[2], 0));
			if ( comp->link.size() ) comp->view(component_view(comp));
			if ( verbose & VERB_DEBUG )
				cout << "New view: " << comp->view() << endl;
			nsel++;
		}
	}

	return nsel;
}

/**
@brief 	Calculates the views associated with each component based on neighbors.
@param 	*model		model structure (views modified).
@return long		number of selected components.

	The view for a vertex is calculated from the vertex vector and the first link.
	The vertex vector points away from the origin and is calculated as follows:
	Angles defined for the vertex:
		vv = normalized sum of the cross products for each pair of links
			constituting an angle
	Angles not defined:
		vv = component coordinates - origin
	The origin is defined as the current zero coordinates, the center-of-mass,
	or from the map.

**/
long		model_calculate_local_views(Bmodel* model)
{
	Bmodel*			mp;
	Bcomponent*  	comp;
	Vector3<double>	normal;	

	long			nsel(0);
	double			offset(0);
	
	if ( verbose & VERB_PROCESS )
		cout << "Calculating views based on neighbors" << endl << endl;
	
    for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		model_set_neighbors(mp, 3);
		for ( comp = mp->comp; comp; comp = comp->next ) {
			normal = component_plane(comp->neighbor, offset);
			if ( verbose & VERB_DEBUG )
				cout << "Old view: " << comp->view() << endl;
			comp->view(View2<float>(normal[0], normal[1], normal[2], 0));
			if ( comp->link[0] ) comp->view(component_view(comp));
			if ( verbose & VERB_DEBUG )
				cout << "New view: " << comp->view() << endl;
			if ( comp->select() > 0 ) nsel++;
		}
	}

	return nsel;
}

/**
@brief 	Analyzes view directions in a model.
@param 	*model		model parameters.
@param 	bin_width	bin width in degrees.
@param 	ref_flag	flag to select the reference vector (0=z-axis, 1=component location)
@return long		number of molecules selected.

	The angle between the component view and a reference vector is calculated.
	A histogram of the angles is constructed and reported.
	The reference vector can be the z-axis, or the component location.

**/
long		model_view_directions(Bmodel* model, int bin_width, int ref_flag)
{
	long				i, j, im, nsel(0);
	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomptype*		comptype = model->type;
	
	if ( !comptype ) {
		cerr << "Error: No component types found!" << endl;
		return -1;
	}

	for ( mp = model; mp; mp = mp->next )
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) nsel++;

//	Vector3<double>	com;
	Vector3<double>	vr(0,0,1), vc;
	Matrix3			mat;
	double			max;
	double			d, da;
	
	double			width = (M_PI/180.0) * (int) (1800.0/nsel);
	if ( bin_width ) width = bin_width * M_PI/180.0;
	
	int				bins = (int) (M_PI/width + 0.999);
	int				a3t[3] = {0,0,0};
	int*			a = new int[bins];
	int*			a3 = new int[3*bins];
	
	for ( i=0; i<bins; i++ ) a[i] = 0;
	for ( i=0; i<3*bins; i++ ) a3[i] = 0;
	
	if ( verbose )
		cout << "Model\tComp\tDist\tAngle\tvx\tvy\tvz" << endl;
	for ( mp = model; mp; mp = mp->next ) {
//		com = model_center_of_mass(mp);
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
//			vr = comp->location() - com;
			d = comp->location().length();
			if ( ref_flag ) vr = comp->location();
			vc = comp->view().vector3();
			da = vc.angle(vr);
			i = (int) (da/width);
			if ( i >= bins ) cout << "i too large: " << i << endl;
			a[i]++;
			for ( j=im=0, max = 0; j<3; j++ ) {
				if ( max < fabs(comp->location()[j]) ) {
					max = fabs(comp->location()[j]);
					im = j;
				}
			}
			a3[im*bins+i]++;
			mat = comp->view().matrix();
			mat = mat.transpose();
//			cout << mat << endl;
			comp->location(mat * vr);
			comp->location().normalize();
			comp->radius(0.05);
			if ( verbose )
				cout << mp->identifier() << tab << comp->identifier() << tab << d << tab << 
					da*180.0/M_PI << tab << comp->location()[0] << tab << 
					comp->location()[1] << tab << comp->location()[2] << endl;
		}
	}

	for ( i=0; i<bins; i++ )
		for ( j=0; j<3; j++ ) a3t[j] += a3[j*bins+i];
	
	if ( verbose ) {
		cout << "\nHistogram:\nAngle\tCount\tExpect\tRatio" << endl;
		for ( i=0; i<bins; i++ ) {
			da = nsel*(cos(width*i) - cos(width*(i+1)))/2.0;
			cout << (i+0.5)*width*180.0/M_PI << tab << a[i] << tab << da << tab << a[i]/da << endl;
		}
		cout << "\nHistogram for three directions of component locations:" << endl;
		cout << "Angle\txCount\tyCount\tzCount\txRatio\tyRatio\tzRatio" << endl;
		for ( i=0; i<bins; i++ ) {
			da = (cos(width*i) - cos(width*(i+1)))/2.0;
			cout << (i+0.5)*width*180.0/M_PI;
			for ( j=0; j<3; j++ ) cout << tab << a3[j*bins+i];
			for ( j=0; j<3; j++ ) cout << tab << a3[j*bins+i]/(a3t[j]*da);
			cout << endl;
		}
		cout << "Total";
		for ( j=0; j<3; j++ ) cout << tab << a3t[j];
		cout << endl << endl;
	}

	delete[] a;
	delete[] a3;
		
	return nsel;
}

/**
@brief 	Determines the hand of a component.
@param 	&s			string with encoded order.
@return int			the hand.

	Requirement: The string must be either a 3 or 6-digit code.
	The reverse of the string is generated and the canonical version 
	compared to the canonical version of the original string.
	The handedness is then returned as the sign of the comparison.
	The reverse of the 6-digit code is defined as reversing the
	first and last 3 digits separately.

**/
int			component_hand(Bstring s)
{
	int			h(0), len = s.length();
	Bstring		ss(s);
	Bstring		sh(s);
		
	if ( len < 3 ) return (0);
	else if ( len == 3 ) {
		ss = ss.canonical(1);
		sh = sh.swap(0,2);
		sh = sh.canonical(1);
	} else if ( len == 6 ) {
		ss = ss.canonical(2);
		sh = sh.swap(0,2);
		sh = sh.swap(3,5);
		sh = sh.canonical(2);
	}
	
	h = ss.compare(sh);
	
	if ( h > 0 ) h = 1;
	if ( h < 0 ) h = -1;
	
	return h;
}

