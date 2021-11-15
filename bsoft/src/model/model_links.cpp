/**
@file	model_links.cpp
@brief	Library routines used for model links processing
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20210830
**/

#include "model_links.h"
#include "model_transform.h"
#include "model_compare.h"
#include "model_util.h"
#include "matrix_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Set up the link list for each component.
@param 	*model		model parameters.
@return long			total number of links.

	Only the first model is processed.

**/
long		model_setup_links(Bmodel* model)
{
	if ( !model ) return 0;
	
	long			i, j, n(0);
	Blink*			link = NULL;
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_setup_links: id=" << model->identifier() << endl;
	
	for ( n=0, link = model->link; link; link = link->next, n++ ) {
		comp1 = link->comp[0];
		comp2 = link->comp[1];
		for ( i=0; i<comp1->link.size() && comp1->link[i] != comp2; i++ ) ;
		for ( j=0; j<comp2->link.size() && comp2->link[j] != comp1; j++ ) ;
		if ( i < comp1->link.size() && j < comp2->link.size() ) {
			if ( verbose & VERB_FULL )
				cerr << "Error: " << comp1->identifier() << " already linked to " << comp2->identifier() << 
					"! (" << i << " -" << j << ")" << endl;
		} else {
			comp1->link.push_back(comp2);
			comp2->link.push_back(comp1);
			comp1->flag.push_back(1);
			comp2->flag.push_back(1);	// flag=1 indicates link
		}
		if ( link->length() < 1e-6 ) link->length(comp1->location().distance(comp2->location()));
	}

	if ( model->link ) model->calculate_normals();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_setup_links: n=" << n << endl;
	
	return n;
}

/**
@brief 	Generates a link list for a model based on vertex separation.
@param 	*model		model parameters.
@param 	maxlength	maximum link length.
@return long			number of links.

	Only the first model is processed.

**/
long		model_link_list_generate(Bmodel* model, double maxlength)
{
	if ( !model ) return 0;
	if ( !model->comp ) return 0;
	
	long			n(0);
	double			linkrad(0.1*maxlength);
	double			d;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	Blink*			link = NULL;

	if ( verbose )
		cout << "Generating a model link list with maximum link length " << maxlength << endl;
	
	link = model->link;
	if ( link ) linkrad = link->radius();
	for ( comp = model->comp; comp && comp->next; comp = comp->next ) if ( comp->select() ) {
		for ( comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp->location().distance(comp2->location());
			if ( d <= maxlength ) {
				link = link_add(&link, comp, comp2, d, linkrad);
				if ( !model->link ) model->link = link;
				n++;
			}
		}
	}

	if ( verbose )
		cout << "Number of links generated       " << n << endl << endl;
	
	return  n;
}

/**
@brief 	Generates a link list for a model based on vertex separation between defined componet types.
@param 	*model		model parameters.
@param 	maxlength	maximum link length.
@param 	type1		first component type.
@param 	type2		second component type.
@param 	flag		1=closest components.
@return long			number of links.

	Only the first model is processed.
	A link is generated for the minimum distance if that is below the maximum length.
	If the flag is set, a link is generated only for the closest second component to the first.

**/
long		model_link_list_generate(Bmodel* model, double maxlength,
				Bstring& type1, Bstring& type2, int flag)
{
	if ( !model ) return 0;
	if ( !model->comp ) return 0;
	if ( type1.length() < 1 ) return model_link_list_generate(model, maxlength);
	if ( type2.length() < 1 ) type2 = type1;

	Bcomptype*		ct1 = model->find_type(type1.c_str());
	Bcomptype*		ct2 = model->find_type(type2.c_str());

	if ( !ct1 || !ct2 ) {
		if ( !ct1 ) cerr << "Error: No component type " << type1 << " found!" << endl;
		if ( !ct2 ) cerr << "Error: No component type " << type2 << " found!" << endl;
		bexit(-1);
	}

	int				closest(flag&1);	// Flag to find the closest second component to the first
	long			n(0);
	double			linkrad(0.1*maxlength);
	double			d, dmin(maxlength), d_avg(0), d_std(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	Bcomponent*		compsel = NULL;
	Blink*			link = model->link;

	if ( verbose ) {
		cout << "Generating a model link list between component types:" << endl;
		cout << "Maximum link length:            " << maxlength << endl;
		cout << "Component type 1:               " << ct1->identifier() << endl;
		cout << "Component type 2:               " << ct2->identifier() << endl;
	}
	
	link = model->link;
	if ( link ) linkrad = link->radius();
	for ( comp = model->comp; comp && comp->next; comp = comp->next ) {
		if ( comp->select() && comp->type() == ct1 ) {
			if ( closest ) {
				dmin = maxlength;
				compsel = NULL;
			}
			for ( comp2 = model->comp; comp2; comp2 = comp2->next ) {
				if ( comp != comp2 && comp2->select() && comp2->type() == ct2 ) {
					if ( model->link->find(comp, comp2) == NULL ) {
						d = comp->location().distance(comp2->location());
						if ( closest ) {
							if ( dmin > d ) {
								dmin = d;
								compsel = comp2;
							}
						} else if ( d <= maxlength ) {
							link = link_add(&link, comp, comp2, d, linkrad);
							if ( !model->link ) model->link = link;
							n++;
							d_avg += d;
							d_std += d*d;
						}
					}
				}
			}
			if ( closest && compsel ) {
				d = comp->location().distance(compsel->location());
				link = link_add(&link, comp, compsel, d, linkrad);
				if ( !model->link ) model->link = link;
				n++;
				d_avg += d;
				d_std += d*d;
			}
		}
	}
	
	if ( n ) {
		d_avg /= n;
		d_std = d_std/n - d_avg*d_avg;
		if ( d_std > 0 ) d_std = sqrt(d_std);
		else d_std = 0;
	}

	if ( verbose ) {
		cout << "Number of links generated:      " << n << endl;
		cout << "Average and standard deviation: " << d_avg << tab << d_std << endl << endl;
	}
	
	return  n;
}

/**
@brief 	Sets the reference link lengths.
@param 	*model		model parameters.
@param 	linklength	reference link length.
@return long			number of links.

	Only the first model is processed.

**/
long		model_set_link_length(Bmodel* model, double linklength)
{
	if ( !model ) return 0;
	
	long			nsel(0);
	Blink*			link;

	if ( verbose & VERB_FULL )
		cout << "Setting reference link lengths to " << linklength << endl << endl;
	
	for ( link = model->link; link; link = link->next ) if ( link->select() ) {
		link->length(linklength);
		nsel++;
	}
	
	return  nsel;
}

/**
@brief 	Set the display radius for all links to a specific value.
@param 	*model		model parameters.
@param 	linkrad		link display radius.
@return long			number of components selected.

	Only the first model is processed.

**/
long		model_set_link_radius(Bmodel* model, double linkrad)
{
	if ( !model ) return 0;
	
	long			nsel(0);
	Blink*			link = NULL;

	for ( link = model->link; link; link = link->next ) if ( link->select() ) {
		link->radius(linkrad);
		nsel++;
	}
	
	return  nsel;
}

Bcomponent*	model_linked_submodel(Bcomponent* comp, Bmodel* model, long& n, Vector3<float>& avg_loc)
{
	if ( !comp->select() ) return NULL;
	
	avg_loc += comp->location();
	n++;

	int				i;
	Bstring			id("1");
	Bcomponent*		comp_sub = NULL;
	Bcomponent*		comp_link = NULL;
	
	if ( model ) {
//		comp_sub = component_add(&model->comp, id);
//		component_copy(comp, comp_sub);
		comp_sub = model->add_component(comp);
	}
	
	comp->select(0);
	
	for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) if ( comp->link[i]->select() ) {
		comp_link = model_linked_submodel(comp->link[i], model, n, avg_loc);
		if ( comp_link ) link_add(&model->link, comp_sub, comp_link, 1, 1);
	}
	
	return comp_sub;
}

/**
@brief	Reduces a model to one component per linked set of components.
@param 	*model			model structure to be modified.
@param 	&submodname		sub-model file name.
@param 	flags			1=generate averaged types.
@return long				number of new components.

	A new reduced set of components is generated.
	Each old set of linked components is saved as a sub-model and referenced as a type.
	Optionally, the old sets of components are averaged and the 
	averaged model saved as a type.

**/
long		model_reduce_linked(Bmodel* model, Bstring& submodname, int flags)
{
	long			i, n, nn(0);
	Bstring			id;
	Bmodel*			mp;
	Bmodel*			model_sub = NULL;
	Bmodel*			mpt = NULL;
	Bmodel*			mps = NULL;
	Bcomponent*		comp;
	Bcomponent*		comp_list = NULL;
	Bcomponent*		comp_new = NULL;
	Bcomponent*		comp_sub = NULL;
	Bcomponent*		comp_temp = NULL;
	Bcomptype*		ct;
	Vector3<float>	avg_loc;
	Transform		t;
	
	if ( verbose & VERB_PROCESS )
		cout << "Reducing models and writing sub-models to file " << submodname << endl << endl;

		
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		comp_list = comp_new = NULL;
		for ( ct = mp->type; ct; ct = ct->next ) ct->component_count(0);
		for ( i=0, comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
//			id = Bstring(++i, "%d");
//			comp_new = component_add(&comp_new, id);
//			if ( !comp_list ) comp_list = comp_new;
//			component_copy(comp, comp_new);
			if ( comp_new ) comp_new = comp_new->add(comp);
			else comp_new = comp_list = new Bcomponent(comp);
			comp_new->select(0);
			comp_new->location(Vector3<float>(0,0,0));
			comp_new->type()->component_count_increment();
			mpt = new Bmodel(++i);
			mpt->model_type(comp_new->type()->identifier());
			n = 0;
			avg_loc = 0;
//			model_linked_submodel(comp, mpt, &comp_new->select(), &comp_new->location());
			model_linked_submodel(comp, mpt, n, avg_loc);
			avg_loc /= n;
			comp_new->location(avg_loc);
			comp_new->select(n);
			for ( mps = model_sub; mps && mps->model_type() != mpt->model_type(); mps = mps->next ) ;
			if ( !mps ) {
				if ( verbose )
					cout << "Adding sub-model " << mpt->model_type() << endl;
//				if ( !model_sub ) model_sub = mps = model_copy(mpt);
				if ( !model_sub ) model_sub = mps = mpt->copy();
				else {
					for ( mps = model_sub; mps->next; mps = mps->next ) ;
//					mps->next = model_copy(mpt);
					mps->next = mpt->copy();
					mps = mps->next;
				}
				comp_new->type()->file_name(submodname.str());
				comp_new->type()->image_number(id.integer());
				mps->add_type(comp_new->type()->identifier());
				for ( comp_sub = mps->comp; comp_sub; comp_sub = comp_sub->next ) {
					comp_sub->velocity(comp_sub->location());	// Keep sum in velocity vector
					comp_sub->select(1);
				}
			} else {	// Fit this model onto the type model
				if ( model_component_number_difference(mpt, mps) == 0 ) {
					t = model_find_transform(mpt, mps);
					model_rotate(mpt, t);
//					comp_new->view = view_from_angle_and_axis3(t.angle, t.axis);
					comp_new->view(View2<float>(t.angle, t.axis));
					if ( verbose )
						cout << "Fitting sub-model " << mpt->model_type() << ": RMSD = " << model_compare(mpt, mps) << endl;
					for ( comp_sub = mps->comp, comp_temp = mpt->comp; comp_sub && comp_temp; 
							comp_sub = comp_sub->next, comp_temp = comp_temp->next ) {
						comp_sub->velocity(comp_sub->velocity() + comp_temp->location());	// Sum in velocity vector
						comp_sub->select_increment();
						cout << comp_sub->identifier() << tab << comp_temp->identifier() << tab << comp_sub->location().distance(comp_temp->location()) << endl;
					}
				}
				model_kill(mpt);
			}
		}
		model_link_list_kill(mp);
		mp->link = NULL;
		component_list_kill(mp->comp);
		mp->comp = comp_list;
		nn += i;
	}

	for ( mps = model_sub; mps; mps = mps->next )
		for ( comp_sub = mps->comp; comp_sub; comp_sub = comp_sub->next )
			if ( comp_sub->select() )
				comp_sub->location(comp_sub->velocity()/comp_sub->select());

	write_model(submodname, model_sub);
	model_kill(model_sub);
	
	return nn;
}

/**
@brief 	Generates links such that each component has at least the given valency.
@param 	*model		model structure.
@param 	valency		the desired valency.
@return long			number of links.

	A distance matrix is calculated and a cutoff determined that ensures
	that each column has at least the given valency number below it.
	The links are then created based on the cutoff.
	Only the first model in the linked list is used.

**/
long		model_links_minimum_valency(Bmodel* model, long valency)
{
	if ( !model ) {
		cerr << "Error: No model selected!" << endl << endl;
		return(-1);
	}
	
	if ( verbose )
		cout << "Generating links to ensure minimum valency of " << valency << endl;

	Matrix			m = model_distance_matrix(model, 0);
	
	double			dcut = matrix_find_cutoff_for_number(m, valency);
	
	long			nl = model_link_list_generate(model, dcut);
	
	return nl;
}


