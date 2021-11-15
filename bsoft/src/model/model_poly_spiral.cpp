/**
@file	model_poly_spiral.cpp
@brief	Functions to generate polyhedra using the spiral algorithm.
@author Bernard Heymann
@date	Created: 20071127
@date	Modified: 20210124
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_poly_spiral.h"
#include "model_mechanics.h"
#include "model_links.h"
#include "model_transform.h"
#include "model_util.h"
#include "Matrix.h"
#include "math_util.h"
#include "linked_list.h"
#include "Vector3.h"
#include "utilities.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
@brief 	Adds a polygon to a model based on the spiral algorithm.

	If there are no components, a new polygon is generated.
	The next polygon is attached to the last component.

@param 	n				polygon order = number of vertices.
@param 	valence			vertex valence.
@param 	*model			model structure.
@return int				number of components saturated with respect to polygons.
**/
int			polygon_add(int n, int valence, Bmodel* model)
{
	int				i, j, mid(0), unsat, nsat;
	int				sat(valence - 1);
	Bstring			id;
	Bstring			comptype = "VER";
	Bcomponent*		comp = model->comp;
	Bcomponent*		comp2 = NULL;
	Bcomponent*		comp3 = NULL;
	Bcomponent*		complist[10];
	Blink*			link = model->link;

	if ( verbose & VERB_PROCESS )
		cout << "Adding polygon: " << n << ":" << endl;
	Bpolygon*		poly = (Bpolygon *) add_item((char **) &model->poly, sizeof(Bpolygon));
	
	poly->select(1);

	for ( comp2 = model->comp, comp = NULL; comp2; comp2 = comp2->next )
//		if (  !comp2->link[sat] ) comp = comp2;					// Find the last vertex with unsaturated valency
		if (  comp2->link.size() < valence ) comp = comp2;		// Find the last vertex with unsaturated valency
	
	if ( !model->comp ) {		// Create the first polygon
		for ( i=mid=0; i<n; i++ ) {
//			id = Bstring(++mid, "%d");
//			comp = component_add(&comp, id);
//			if ( !model->comp ) model->comp = comp;
			if ( comp ) comp = comp->add(++mid);
			else model->comp = comp = new Bcomponent(++mid);
			comp->type(model->add_type(comptype));
			comp->select(1);			// Member of one polygon
			if ( verbose & VERB_FULL )
				cout << "Adding vertex: " << id << endl;
			poly->comp.push_back(comp);
			if ( i ) {
				link = link_add(&link, comp2, comp, 1, 1);
				if ( !model->link ) model->link = link;
				if ( verbose & VERB_FULL )
					cout << "Adding link: " << comp2->identifier() << " - " << comp->identifier() << endl;
			}
			comp2 = comp;
		}
		link = link_add(&link, model->comp, comp, 1, 1);
		poly->closed(1);
		if ( verbose & VERB_FULL )
			cout << "Adding link: " << model->comp->identifier() << " - " << comp->identifier() << endl;
	} else if ( !comp ) {		// Add the last polygon
		for ( comp = model->comp; comp && comp->select() > sat; comp = comp->next ) ;	// Find the first vertex with less than 3 polygons
		if ( comp ) {
			poly->comp.push_back(comp);
			comp->select_increment();
			for ( i = 1; i<n && comp; i++ ) {
				for ( j=0; j<comp->link.size() && comp->link[j]->select() > sat; j++ ) ;	// Find the next vertex with less than 3 polygons
				if ( j < valence ) {
					comp = comp->link[j];
					poly->comp.push_back(comp);
					comp->select_increment();
				} else comp = NULL;		// Polygon not complete
			}
			if ( comp ) {
				poly->closed(1);
			}
		}
	} else {
//		for ( mid=0, comp2 = model->comp; comp2; comp2 = comp2->next ) mid++;
		for ( mid=0, comp2 = model->comp; comp2; comp2 = comp2->next ) {
			i = stol(comp2->identifier());
			if ( mid < i ) mid = i;
		}
		comp2 = comp;							// Polygon end component
		complist[0] = comp;
//		for ( j=0; comp2->link[j]; j++ ) ;		// Find the last link of the last vertex
		j = comp2->link.size();					// Find the last link of the last vertex
		j--;
		if ( j < 0 ) {
			cerr << "Error: " << comp2->identifier() << " does not have any links!" << endl;
			bexit(-1);
		}
		if ( verbose & VERB_FULL )
			cout << "Trace:\t\t" << comp2->identifier();
		for ( unsat = 0, i = 1; i<n && unsat == 0 && comp2; i++ ) {	// Trace to the start component
			if ( i > 1 )
				for ( j=0; j<comp2->link.size(); j++ )		// Find the next vertex with less than 3 polygons
					if ( comp2->link[j] != comp3 && comp2->link[j]->select() < valence ) break;
			if ( j > sat ) {
				if ( verbose & VERB_FULL ) {
					cerr << endl << "Error: No link with the required properties" << endl;
					cerr << "Component: " << comp2->identifier() << endl;
					for ( j=0; j<comp2->link.size(); j++ )
						cerr << "Link: " << comp2->link[j]->identifier() << "  Flag: " 
							<< comp2->flag[j] << "  Select: " << comp2->link[j]->select() << endl;
				}
				comp2 = NULL;
			} else {
				comp3 = comp2;
				comp2 = comp2->link[j];
				complist[i] = comp2;
//				if ( !comp2->link[sat] ) unsat = 1;	// Finish when a vertex with unsaturated valency is found
				if ( comp2->link.size() < valence ) unsat = 1;	// Finish when a vertex with unsaturated valency is found
				if ( verbose & VERB_FULL )
					cout << tab << comp2->identifier();
			}
		}
		if ( verbose & VERB_FULL )
			cout << endl;
		if ( comp2 ) {
			if ( verbose & VERB_FULL )
				cout << "Backtrace:";
			for ( j=i-1, i=0; j>=0; j--, i++ ) {
				poly->comp.push_back(complist[j]);
				poly->comp[i]->select_increment();
				if ( verbose & VERB_FULL )
					cout << tab << poly->comp[i]->identifier();
			}
			if ( verbose & VERB_FULL )
				cout << endl;
			comp3 = comp;
			for ( ; i<n; i++ ) {				// New vertices
//				id = Bstring(++mid, "%d");
//				comp3 = component_add(&comp3, id);
				comp3 = comp3->add(++mid);
				comp3->type(model->add_type(comptype));
				comp3->select(1);			// Member of one polygon
				if ( verbose & VERB_FULL )
					cout << "Adding vertex: " << id << endl;
				poly->comp.push_back(comp3);
				link = link_add(&link, comp, comp3, 1, 1);
				if ( verbose & VERB_FULL )
					cout << "Adding link: " << comp->identifier() << " - " << comp3->identifier() << endl;
				comp = comp3;
			}
			link = link_add(&link, comp, comp2, 1, 1);
			poly->closed(1);
			if ( verbose & VERB_FULL )
				cout << "Adding link: " << comp->identifier() << " - " << comp2->identifier() << endl;
		}
	}
	
	for ( nsat=0, comp = model->comp; comp; comp = comp->next ) {
		if ( comp->select() == valence ) nsat++;
		else if ( comp->select() > valence ) {
			if ( verbose & VERB_FULL )
				cerr << "Error: Valence too large for " << comp->identifier() << ": " << comp->select() << " > " << valence << endl;
			nsat += 100;
//			bexit(-1);
		}
	}
	
	if ( verbose & VERB_FULL ) {
		if ( poly->size() != n )
			cerr << "Error: Polygon incomplete: " << poly->size() << "/" << n << endl;
		if ( !poly->closed() )
			cerr << "Error: Polygon not closed" << endl;
	}
	
	return nsat;
}

int			seq_count_adj_pentagons(Bstring& seq)
{
	int				i, j, adj(0);
	
	for ( i=0, j=1; j<seq.length(); i++, j++ )
		if ( seq[i] == '5' && seq[j] == '5' ) adj++;

	return adj;
}

int			model_vertex_types_in_first_polygons(Bmodel* model, int npoly)
{
	int				i, j, n;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bpolygon*		poly;
	
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( n=0, poly = mp->poly; poly && n<npoly; poly = poly->next, n++ ) {
			for ( i=0; poly->comp[i]; i++ ) {
				if ( i+1 < poly->size() ) comp = poly->comp[i+1];
				else comp = poly->comp[0];
				for ( j=0; poly->comp[i]->link[j] != comp; j++ ) ;
				poly->comp[i]->flag[j] = poly->size();
			}
		}
	}
	
	return 0;
}

int			model_adjacent_pentagons(Bmodel* model)
{
	int				i, n, adj(0);
	Bcomponent*		comp;

	model_vertex_types(model);
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		for ( i=n=0; i<comp->link.size(); i++ ) if ( comp->flag[i] == 5 ) n++;
		if ( n > 1 ) adj++;
	}
	
	return adj;
}

int			model_adjacent_pentagons_in_first_polygons(Bmodel* model, int npoly)
{
	int				i, n, adj(0);
	Bcomponent*		comp;

	model_vertex_types_in_first_polygons(model, npoly);
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		for ( i=n=0; i<comp->link.size(); i++ ) if ( comp->flag[i] == 5 ) n++;
		if ( n > 1 ) adj++;
	}
	
	return adj;
}

/**
@brief 	Generates a polyhedron using the spiral algorithm.
@param 	&seq			polygon sequence.
@param 	valence			vertex valence.
@param 	requirements	polyhedron requirements.
@return Bmodel*			new model, NULL if generation failed.

	Polygons are added based on the given sequence.
	The success of the algorithm is checked using the indicated requirements:
		0		only a polyhedron consistency check is done
		1		the exact number of vertices must be obtained
		2		only a polyhedron with isolated pentagons is accepted
	The generation fails when an incorrect number of vertices are added
	or some of the vertices have incorrect valency.

**/
Bmodel*		model_poly_spiral(Bstring& seq, int valence, int requirements)
{
	int				npoly = seq.length();
	int				vertices = 20 + 2*(npoly - 12);
	
	if ( verbose ) {
		cout << "Generating a new model from a polyhedral sequence:" << endl;
		cout << "Polygon sequence:              " << seq << endl;
		cout << "Valence:                       " << valence << endl;
		cout << "Requirements:                  " << requirements << endl;
	}
	
	Bmodel*			model = new Bmodel(seq);
	
	int				i, nvert(1), nsat(0), adj(0), good(1);
	
	for ( i=1; nsat < nvert && i<=npoly; i++ ) {
		nvert = seq[i-1]-48;
		if ( verbose & VERB_PROCESS )
			cout << "Polygon " << i << ":\t" << nvert << endl;
		nsat = polygon_add(nvert, valence, model);
		nvert = model->component_count();
		if ( verbose & VERB_PROCESS )
			cout << "Saturation:\t" << nsat << "/" << nvert << endl;
	}

	if ( ( i = model_polyhedron_check(model, valence) ) ) good = 0;
	if ( requirements%2 == 1 ) if ( nvert != vertices ) good = 0;
	if ( requirements > 1 ) if ( model_adjacent_pentagons(model) ) good = 0;

	if ( !good ) {
		if ( verbose & VERB_FULL ) {
			cout << "Model rejected:";
			if ( nvert != vertices )
				cout << "\tvertices = " << nvert << endl;
			else if ( nsat > nvert )
				cout << "\tsaturation = " << nsat << " > " << nvert << endl;
			else if ( adj )
				cout << "\tadjacent pentagons = " << adj << endl;
			else
				cout << "\tcheck = " << i << endl;
			cout << seq << endl;
		}
		model_kill(model);
		model = NULL;
	}
	
	return model;
}

int			model_poly_spiral_regularize(Bmodel* model)
{
	model_set_link_length(model, 1);
	
	return model_regularize(model, 10000, 0, 0, 0.2, 0.001, 0, 0, 0.1, 0.01);
}

/**
@brief 	Checks a polyhedron for accuracy and completeness.
@param 	*model			model structure.
@param 	valence			vertex valence.
@return int				number of failed conditions.

	Every component must have the required number of links = valence.
	Every component must have the required number of polygons = valence.
	The polyhedron must adhere to Euler's formula:
		components + polygons - links = 2

**/
int			model_polyhedron_check(Bmodel* model, int valence)
{
	int				ncomp(0), nlink(0), npoly(0), err(0);
	Bmodel*			mp;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	Bpolygon*		poly = NULL;
	
	if ( verbose & VERB_PROCESS )
		cout << endl << "Checking the model:" << endl;
	
	for ( ncomp=nlink=npoly=0, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next, ncomp++ ) {
			if ( comp->select() != valence ) {
				if ( verbose & VERB_PROCESS )
					cerr << comp->identifier() << ": select = " << comp->select() << endl;
				err++;
			}
//			for ( j=0; comp->link[j]; j++ ) ;
			if ( comp->link.size() != valence ) {
				if ( verbose & VERB_PROCESS )
					cerr << comp->identifier() << ": valency = " << comp->link.size() << endl;
				err++;
			}
		}
		for ( link = mp->link; link; link = link->next, nlink++ ) ;
		for ( poly = mp->poly; poly; poly = poly->next, npoly++ )
			if ( !poly->closed() ) {
				if ( verbose & VERB_PROCESS )
					cerr << npoly << ": polygon not closed!" << endl;
				err++;
			}
	}
	
	if ( ncomp - nlink + npoly != 2 ) {
		if ( verbose & VERB_PROCESS )
			cerr << "Euler's formula incorrect: " << ncomp << " - " << nlink << " + " 
				<< npoly << " = " << ncomp - nlink + npoly << endl;
		err++;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of components:  " << ncomp << endl;
		cout << "Number of links:       " << nlink << endl;
		cout << "Number of polygons:    " << npoly << endl;
		cout << "Errors:                " << err << endl << endl;
	}
	
	return err;
}

/**
@brief 	Generates a polyhedron using a given sequence.
@param 	&seq			polygon sequence.
@param 	valence			vertex valence.
@param 	enantiomorph	flag to generate enantiomorphs.
@param 	requirements	polyhedron requirements.
@param 	nm				current number of models (before creating this one).
@return Bmodel*			new model, NULL if generation failed.

	A single model is generated based on the sequence.

**/
Bmodel*		model_poly_gen_sequence(Bstring& seq, int valence, int enantiomorph, int requirements, int nm)
{
	vector<double>	table;

	return model_poly_gen_sequence(seq, valence, enantiomorph, requirements, nm, table);
}

/**
@brief 	Generates a polyhedron using a given sequence.
@param 	&seq			polygon sequence.
@param 	valence			vertex valence.
@param 	enantiomorph	flag to generate enantiomorphs.
@param 	requirements	polyhedron requirements.
@param 	nm				current number of models (before creating this one).
@param 	&table			table for sets of eigenvalues.
@return Bmodel*			new model, NULL if generation failed.

	A single model is generated based on the sequence.
	A table is used to keep track of sets of eigenvalues of previous models
	to avoid generating redundant models.

**/
Bmodel*		model_poly_gen_sequence(Bstring& seq, int valence, int enantiomorph, int requirements, int nm, vector<double>& table)
{
	int				j(0), k;
	int				d(1), enant(1);
	Vector3<double>	origin, normal(1,0,0);
	Bstring			sym_label("C1");
	
	int				vertices = 20 + 2*(seq.length() - 12);

	Bstring			id = Bstring(vertices, "%d_");

	if ( verbose & VERB_PROCESS )
		cout << nm+1 << ": " << seq << endl;

	Bmodel*			model = model_poly_spiral(seq, valence, requirements);
	if ( !model ) return model;
	
	vector<double>	eigval = model_poly_sphere_coor(model);
	
	if ( table.size() ) {
		if ( nm )
			for ( j=0, d=1; j<nm*vertices && table[j] && d; j+=vertices )
				for ( k=0, d=0; k<vertices; k++ )
					if ( fabs(eigval[k] - table[j+k]) > 1e-6 ) d++;
		if ( d ) for ( k=0; k<vertices; j++, k++ )
			table[j] = eigval[k];
		else
			id += Bstring(j/vertices, "%04d");
	}
	
	if ( !d ) {
		model_kill(model);
		if ( verbose )
//			cout << "Model rejected because of redundancy" << endl;
			cout << id << tab << seq << endl;
		return NULL;
	}
	
	id = Bstring(count_list((char *)model->comp), "%d_");
	Bstring			id2(++nm, "%04d");
	id2 = id + id2;
	model->identifier(id2);
	model->model_type(model->identifier());
	model_poly_spiral_regularize(model);
	sym_label = model_poly_find_symmetry(model, 0.1);

	if ( verbose )
		cout << model->identifier() << tab << seq << endl;
	
	model_set_component_radius(model, 0.1);
	model_set_link_radius(model, 0.1);
	
	Bstring			filename = "temp/" + model->identifier() + ".cmm";
	write_model(filename, model);
	
	if ( enantiomorph ) {
		if ( sym_label.contains("s") || sym_label.contains("v") || 
				sym_label.contains("d") || sym_label.contains("h") ) enant = 0;
		if ( enant ) {
//			model->next = model_copy(model);
			model->next = model->copy();
			id2 = id + Bstring(++nm, "%04d");
			model->next->identifier(id2);
			model_reflect(model->next, normal, origin);
		}
	}
	
	return model;
}

/**
@brief 	All polyhedra are generated for a given number of vertices. 
@param 	vertices			number of vertices.
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_permutations(int vertices, int valence, int enantiomorph)
{
	long			nt = (long) pow(10, (4*(log(1.0L*vertices) - 3.1)));	// Maximum number of comparison table entries
	if ( verbose )
		cout << "Generating polyhedra with permutations:" << endl << 
			"Vertices = " << vertices << "  Table rows = " << nt << endl;
	if ( nt < 100 ) nt = 100;
	if ( nt > 1e8 ) {
		cerr << "The number of table entries is too large! (" << nt << ")" << endl;
		return NULL;
	}
	
	int				i, np(1), nm(0);
	vector<double>	table(nt*vertices, 0);
	Bstring			seq = "555555555555";			// Initial sequence
	for ( i=0; i<(vertices-20)/2; i++ ) seq += '6';
//	for ( i=0; i<nt*vertices; i++ ) table[i] = 0;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_poly_gen_permutations: " << i << ": " << seq << ": start" << endl;
	
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			new_mp = NULL;
	
	for ( i=1, nm=0; np; i++ ) {
		new_mp = model_poly_gen_sequence(seq, valence, enantiomorph, 1, nm, table);
		if ( new_mp ) {
			if ( !model ) model = new_mp;
			else mp->next = new_mp;
			mp = new_mp;
			nm++;
			if ( verbose )
				cout << mp->identifier() << " " << i << ": " << seq << ": " << mp->symmetry() << endl;
			if ( new_mp->next ) {
				mp = new_mp->next;
				nm++;
				if ( verbose )
					cout << mp->identifier() << " " << i << ": " << seq << ": " << mp->symmetry() << endl;
			}
		}
		np = next_permutation(seq);
	}

//	delete[] table;
	
	for ( nm=0, mp = model; mp; mp = mp->next ) nm++;
	if ( verbose )
		cout << "Models made:                    " << nm << "/" << --i << endl << endl;
	
	return model;
}

/**
@brief 	A cone tip with 5 pentagons and a base with 7 pentagons are use to generate polyhedra. 
@param 	tip					number of vertices in the tip.
@param 	body				number of vertices in the body.
@param 	base				number of vertices in the base.
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@param 	requirements		polyhedron requirements.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_cone(int tip, int body, int base, int valence, int enantiomorph, int requirements)
{
	if ( tip < 6 ) tip = 6;
	if ( body < 1 ) body = 1;
	if ( base < 15 ) base = 15;
	
	int				i;
	Bstring			stip;
	Bstring			sbody('6', body);
	Bstring			sbase;

	// Tip starting cases
	int				type = 4;
	switch ( type ) {
		case 1: stip = "5666665656565"; break;
		case 2: stip = "56666656566565"; break;
		case 3: stip = "56666665656565"; break;
		case 4: stip = "665656565665"; break;	// Equivalent of 1
		case 5: stip = "6666565656565"; break;
		default: stip = Bstring('5', 5) + Bstring('6', tip-5);
	}
	while ( stip.length() < tip ) stip += "6";

	for ( i=0; i<7; i++ ) sbase += "56";
	sbase += Bstring('6', base-14);

	return model_poly_gen_3part(stip, sbody, sbase, valence, enantiomorph, requirements);
}

/**
@brief 	Two icosahedral tips are set up and polyhedra generated by rotating the 2 tips. 
@param 	ttop				number of polygons between pentagons in the tip.
@param 	tbody				number of body rings.
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@param 	requirements		polyhedron requirements.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_lozenge(int ttop, int tbody, int valence, int enantiomorph, int requirements)
{
	int				i, j, h, k, vertices, nm, nt;
	int				top = 1 + (5*ttop*(ttop + 1))/2;
	int				dm = (5*ttop + 1)/2;
	int				body = tbody*dm;
	
	Bstring			stip, sbody, sbase;
	Bstring			stip_rpt, sbase_rpt;
	Bstring			id, seq;
	
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			new_mp = NULL;
	
	for ( nm=0, nt=1, i=0; i<ttop; i++ ) {
		stip_rpt = 0;
		for ( k=0; k<i; k++ ) stip_rpt += "6";
		stip_rpt += "5";
		for ( k++; k<ttop; k++ ) stip_rpt += "6";
		stip = "5";
		stip += Bstring('6', top - 1 - 5*ttop);
		for ( k=0; k<5; k++ ) stip += stip_rpt;
		for ( h=body; h<=body+2*dm; h++ ) {
			if ( h ) sbody = Bstring('6', body+h);
			for ( j=0; j<ttop; j++, nt++ ) {
				sbase_rpt = 0;
				for ( k=0; k<j; k++ ) sbase_rpt += "6";
				sbase_rpt += "5";
				for ( k++; k<ttop; k++ ) sbase_rpt += "6";
				sbase = 0;
				for ( k=0; k<5; k++ ) sbase += sbase_rpt;
				sbase += Bstring('6', top - 1 - 5*ttop);
				sbase += "5";
				seq = stip + sbody + sbase;
				new_mp = model_poly_gen_sequence(seq, valence, enantiomorph, requirements, nm);
				if ( new_mp ) {
					vertices = 20 + 2*(seq.length() - 12);
					id = Bstring(vertices, "%d_") + Bstring(++nm, "%04d");
					new_mp->identifier(id);
					cout << new_mp->identifier() << " " << nt << ": " << seq << ": " << new_mp->symmetry() << endl;
					if ( mp ) mp->next = new_mp;
					else model = new_mp;
					mp = new_mp;
				}
			}
		}
	}

	cout << "Models made:                   " << nm << "/" << --nt << endl << endl;

	return model;
}

/**
@brief 	An icosahedral tip and a 6-fold base is set up and polyhedra generated by 
	rotating the tip and base. 
@param 	ttop				number of polygons between pentagons in the tip.
@param 	tbody				number of body rings.
@param 	tbase				number of polygons between pentagons in the base
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@param 	requirements		polyhedron requirements.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_coffin(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements)
{
	int				i, j, h, k, vertices, nm, nt;
	int				top = 1 + (5*ttop*(ttop + 1))/2;
	int				dm = (5*ttop + 6*tbase + 1)/2;
	int				body = tbody*dm;
	int				base = 1 + 3*tbase*(tbase + 1);
	
	Bstring			stip, sbody, sbase;
	Bstring			stip_rpt, sbase_rpt;
	Bstring			id, seq;
	
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			new_mp = NULL;
	
	for ( nm=0, nt=1, i=0; i<ttop; i++ ) {
		stip_rpt = 0;
		for ( k=0; k<i; k++ ) stip_rpt += "6";
		stip_rpt += "5";
		for ( k++; k<ttop; k++ ) stip_rpt += "6";
		stip = "5";
		stip += Bstring('6', top - 1 - 5*ttop);
		for ( k=0; k<5; k++ ) stip += stip_rpt;
		for ( h=0; h<dm; h++ ) {
			sbody = Bstring('6', body+h);
			for ( j=0; j<tbase; j++, nt++ ) {
				sbase_rpt = 0;
				for ( k=0; k<j; k++ ) sbase_rpt += "6";
				sbase_rpt += "5";
				for ( k++; k<tbase; k++ ) sbase_rpt += "6";
				sbase = 0;
				for ( k=0; k<6; k++ ) sbase += sbase_rpt;
				sbase += Bstring('6', base - 6*tbase);
				seq = stip + sbody + sbase;
				new_mp = model_poly_gen_sequence(seq, valence, enantiomorph, requirements, nm);
				if ( new_mp ) {
					vertices = 20 + 2*(seq.length() - 12);
					id = Bstring(vertices, "%d_") + Bstring(++nm, "%04d");
					new_mp->identifier(id);
					cout << new_mp->identifier() << " " << nt << ": " << seq << ": " << new_mp->symmetry() << endl;
					if ( mp ) mp->next = new_mp;
					else model = new_mp;
					mp = new_mp;
				}
			}
		}
	}

	cout << "Models made:                   " << nm << "/" << --nt << endl << endl;

	return model;
}

/**
@brief 	An icosahedral tip and a 6-fold base is set up and polyhedra generated by 
	permuting the tip and base. 
@param 	ttop				number of polygons between pentagons in the tip.
@param 	tbody				number of body rings.
@param 	tbase				number of polygons between pentagons in the base
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@param 	requirements		polyhedron requirements.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_coffin_loose(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements)
{
	int				i;
	int				top = 1 + (5*ttop*(ttop + 1))/2;
	int				body = tbody*(5*ttop + 6*tbase + 1)/2;
	int				base = 1 + 3*tbase*(tbase + 1);
	
	Bstring			stip("5");
	Bstring			sbody('6', body);
	Bstring			sbase;
	Bstring			stip_rpt("5");
	Bstring			sbase_rpt("5");
	
	for ( i=1; i<ttop; i++ ) stip_rpt += "6";
	for ( i=1; i<tbase; i++ ) sbase_rpt += "6";

	stip += Bstring('6', top - 1 - 5*ttop);
	for ( i=0; i<5; i++ ) stip += stip_rpt;

	for ( i=0; i<6; i++ ) sbase += sbase_rpt;
	sbase += Bstring('6', base - 6*tbase);

	return model_poly_gen_3part(stip, sbody, sbase, valence, enantiomorph, requirements);
}

/**
@brief 	An icosahedral tip and a 6-fold base is set up and polyhedra generated by moving pentagons around. 
@param 	ttop				number of polygons between pentagons in the tip.
@param 	tbody				number of body rings.
@param 	tbase				number of polygons between pentagons in the base
@param 	valence				vertex valence.
@param 	enantiomorph		flag to generate enantiomorphs.
@param 	requirements		polyhedron requirements.
@return Bmodel*				new model, NULL if generation failed.
**/
Bmodel*		model_poly_gen_coffin_jiggle(int ttop, int tbody, int tbase, int valence, int enantiomorph, int requirements)
{
	int				i;
	int				top = 1 + (5*ttop*(ttop + 1))/2;
	int				body = tbody*(5*ttop + 6*tbase + 1)/2;
	int				base = 1 + 3*tbase*(tbase + 1);
	
	Bstring			stip("5");
	Bstring			sbody('6', body);
	Bstring			sbase;
	Bstring			stip_rpt("5");
	Bstring			sbase_rpt("5");
	
	for ( i=1; i<ttop; i++ ) stip_rpt += "6";
	for ( i=1; i<tbase; i++ ) sbase_rpt += "6";

	stip += Bstring('6', top - 1 - 5*ttop);
	for ( i=0; i<5; i++ ) stip += stip_rpt;

	for ( i=0; i<6; i++ ) sbase += sbase_rpt;
	sbase += Bstring('6', base - 6*tbase);
	
	Bstring			seq = stip + sbody + sbase;
	
	return model_poly_gen_move_pentagons(seq, valence, enantiomorph, requirements);
}

/**
@brief 	Generates several polyhedra by permuting the first and last parts of a three-part sequence.
@param 	stip			tip sequence (permuted).
@param 	sbody			body sequence (all hexagons).
@param 	sbase			base/end sequence (permuted).
@param 	valence			vertex valence.
@param 	enantiomorph	flag to generate enantiomorphs.
@param 	requirements	polyhedron requirements.
@return Bmodel*			new model, NULL if generation failed.

	A recursive algorithm is used to generate several models by shifting the
	positions of pentagons in a sequence.

**/
Bmodel*		model_poly_gen_3part(Bstring stip, Bstring sbody, Bstring sbase, int valence, int enantiomorph, int requirements)
{
	int				i, np(1), npb(1), nm(0), max(1000);
	int				nt(1000);			// Maximum number of comparison table entries
	Bstring			sbase_one;
	
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			new_mp = NULL;

	Bstring			seq = stip + sbody + sbase;
	int				vertices = 20 + 2*(seq.length() - 12);
//	double*			table = new double[nt*vertices];
	vector<double>	table(nt*vertices, 0);
//	for ( i=0; i<nt*vertices; i++ ) table[i] = 0;

	if ( verbose )
		cout << "Start: " << seq << endl;

	for ( i=0, nm=0; nm<max && np; np = next_permutation(stip) )
				if ( seq_count_adj_pentagons(stip) < 1 ) {
		sbase_one = sbase;
		for ( npb=1; nm<max && npb; npb = next_permutation(sbase_one) )
					if ( seq_count_adj_pentagons(sbase_one) < 1 ) {
			i++;
			seq = stip + sbody + sbase_one;
			new_mp = model_poly_gen_sequence(seq, valence, enantiomorph, requirements, nm, table);
			if ( new_mp ) {
				if ( !model ) model = new_mp;
				else mp->next = new_mp;
				mp = new_mp;
				nm++;
				cout << mp->identifier() << " " << i << ": " << seq << ": " << mp->symmetry() << endl;
				if ( new_mp->next ) {
					mp = new_mp->next;
					nm++;
					cout << mp->identifier() << " " << i << ": " << seq << ": " << mp->symmetry() << endl;
				}
			}
		}
	}

//	delete[] table;
	
	for ( nm=0, mp = model; mp; mp = mp->next ) nm++;

	cout << "Models made:                   " << nm << "/" << --nt << endl << endl;
	
	return model;
}

Bmodel*		move_pentagons(int len, int p[12], int n, int valence, int enantiomorph, int requirements, int* nm)
{
	int				i, j, pc[12];
	Bstring			seq;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			new_mp = NULL;
		
	if ( n < 12 ) {
		for ( i=0; i<12; i++ ) pc[i] = p[i];
		for ( i = -1; i < 2; i++ ) {
			pc[n] = p[n] + i;
			if ( pc[n] >= 0  && pc[n] < len ) {
				new_mp = move_pentagons(len, pc, n+1, valence, enantiomorph, requirements, nm);
				if ( new_mp ) {
					if ( !model ) model = mp = new_mp;
					else mp->next = new_mp;
					while ( mp->next ) {
						(*nm)++;
						mp = mp->next;
					}
				}
			}
		}
	} else {
		for ( i=j=0; i<len; i++ ) if ( i==p[j] ) {
			seq += "5";
			j++;
		} else {
			seq += "6";
		}
		model = model_poly_gen_sequence(seq, valence, enantiomorph, requirements, *nm);
		if ( model )
			cout << model->identifier() << " " << i << ": " << seq << ": " << model->symmetry() << endl;
	}
	
	return model;
}

/**
@brief 	Generates many polyhedrons using a given sequence and moving pentagons around.
@param 	&seq			polygon sequence.
@param 	valence			vertex valence.
@param 	enantiomorph	flag to generate enantiomorphs.
@param 	requirements	polyhedron requirements.
@return Bmodel*			new model, NULL if generation failed.

	A recursive algorithm is used to generate several models by shifting the
	positions of pentagons in a sequence.

**/
Bmodel*		model_poly_gen_move_pentagons(Bstring& seq, int valence, int enantiomorph, int requirements)
{
	int				i, j, p[12], nm(0);
	
	for ( i=j=0; i<seq.length() && j<12; i++ ) if ( seq[i] == '5' ) p[j++] = i;
	
	Bmodel*			model = move_pentagons(seq.length(), p, 0, valence, enantiomorph, requirements, &nm);
	
	return model;
}

