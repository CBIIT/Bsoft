/**
@file	model_path.cpp
@brief	Library routines used for model processing
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20150208
**/

#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_compare.h"
#include "mol_transform.h"
#include "mol_compare.h"
#include "mol_util.h"
#include "symmetry.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen



/**
@brief 	Calculates the shortest path between every pair of vertices.
@param 	*model			model structure.
@return Matrix					path matrix.

	The adjacency matrix is calculated as a starting point.
	Only the first model in the linked list is used.

**/
Matrix		model_shortest_path(Bmodel* model)
{
	Matrix			mat = model_adjacency_matrix(model);
	
	int				h, i, j, k, ncomp;
	double			s = 0;
	Bcomponent*		comp;
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next ) ncomp++;
	
	if ( verbose & VERB_FULL )
		cout << "Cycle: Coverage" << endl;
	for ( h=2; h<10 && s<ncomp*(ncomp-1); h++ ) {
		for ( i=0; i<ncomp; i++ ) {
			for ( j=0; j<ncomp; j++ ) {
				if ( i!=j && mat[i][j] && mat[i][j]<h ) {
					for ( k=0; k<ncomp; k++ ) {
						s = mat[i][j] + mat[k][j];
						if ( i!=k && j!=k && mat[k][j] && s==h && mat[k][i] == 0 )
							mat[k][i] = s;
					}
				}
			}
		}
		for ( i=0, s=0; i<ncomp; i++ )
			for ( j=0; j<ncomp; j++ ) if ( mat[i][j] ) s++;
		if ( verbose & VERB_FULL )
			cout << h << ":\t" << s << endl;
	}
	
//	if ( verbose & VERB_FULL )
//		print_matrix(ncomp, mat);
	
	return mat;
}

/**
@brief 	Calculates the Wiener index.
@param 	*model			model structure.
@return double					average Wiener index.

	The Wiener index is calculated as the sum of the elements of
	the shortest path matrix.
	Only the first model in the linked list is used.
Reference:  H. Wiener, J. Am. Chem. Soc., 1947, 69, 17.

**/
double		model_wiener_index(Bmodel* model)
{
	int				i, j, n, ncomp;
	double			wi, wia;
	Matrix			mat;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	if ( verbose & VERB_RESULT )
		cout << "Model\tWienerIndex" << endl;
	for ( n=0, wia=0, mp = model; mp; mp = mp->next, n++ ) {
		for ( ncomp = 0, comp = mp->comp; comp; comp = comp->next ) ncomp++;
		if ( ncomp ) {
			mat = model_shortest_path(mp);
			for ( i=0, wi=0; i<ncomp; i++ )
				for ( j=0; j<ncomp; j++ ) wi += (long) mat[i][j];
			wi /= 2;
			wia += wi;
		}
		if ( verbose & VERB_RESULT )
			cout << mp->identifier() << tab << wi << endl;
	}

	if ( n ) wia /= n;

	if ( verbose & VERB_RESULT )
		cout << "Avg:\t" << wia << endl << endl;
		
	return wia;
}

int			comp_link_check(Bcomponent* comp, Bcomponent* comp2)
{
	int			j;
	
	for ( j=0; j<comp->link.size() && comp->link[j]; j++ )
		if ( comp->link[j] == comp2 ) return 1;
	
	return 0;
}

int			comp_link_check(Bmodel* model, string& id1, string& id2)
{
	Bcomponent*		comp, *comp2;
	
	for ( comp = model->comp; comp && comp->identifier() != id1; comp = comp->next ) ;
	
	for ( comp2 = model->comp; comp2 && comp2->identifier() != id2; comp2 = comp2->next ) ;
	
	return comp_link_check(comp, comp2);
}

Blink*		model_link_add(Bmodel* model, string& id1, string& id2)
{
	Bcomponent*		comp, *comp2;
	
	for ( comp = model->comp; comp && comp->identifier() != id1; comp = comp->next ) ;
	
	for ( comp2 = model->comp; comp2 && comp2->identifier() != id2; comp2 = comp2->next ) ;
	
	return link_add(&model->link, comp, comp2);
}

int			model_check_path(Bmodel* model, Bcomponent** path)
{
	int				same(0), i, n;
	Bmodel*			modpath;
	Bcomponent*		comp;

	for ( n=0, comp = model->comp; comp; comp = comp->next ) n++;
	
	for ( modpath = model->next; modpath && same<n; modpath = modpath->next ) {
		for ( same=0, i=1; i<n; i++ )
			if ( comp_link_check(modpath, path[i-1]->identifier(), path[i]->identifier()) ) same++;
		if ( comp_link_check(modpath, path[0]->identifier(), path[n-1]->identifier()) ) same++;
	}
	
	if ( same < n ) return 1;
	
	return 0;
}

double		model_add_path(Bmodel* model, Bcomponent** path)
{
	int				i, n;
	double			len(0);
	Bcomponent*		comp;
	Blink*			link;
	
	for ( n=0, comp = model->comp; comp; comp = comp->next ) n++;
	
//	Bmodel*		modpath = model_copy(model);
	Bmodel*		modpath = model->copy();
	
	model_link_list_kill(modpath);
	modpath->link = NULL;
	
	for ( i=1; i<n; i++ ) {
		link = model_link_add(modpath, path[i-1]->identifier(), path[i]->identifier());
		len += link->length();
	}
	link = model_link_add(modpath, path[0]->identifier(), path[n-1]->identifier());
	len += link->length();
	
	for ( ; model->next; model = model->next ) ;
	model->next = modpath;
	
	return len;
}

int			comp_pick_next(int i, int n, Bcomponent* comp, Bcomponent** path, Bmodel* model)
{
	path[i] = comp;
	
	int				j, k, use;
	double			len;
	Bcomponent*		comp_link;

	for ( j=0; j<comp->link.size() && comp->link[j]; j++ ) {
		comp_link = comp->link[j];
		for ( use=1, k=0; k<i && use; k++ )		// Check if already in path
			if ( comp_link == path[k] ) use = 0;
		if ( use )							// Add to path and check next
			comp_pick_next(i+1, n, comp_link, path, model);
		else {
			if ( i+2 > n && comp_link_check(path[0], path[i]) ) {
				// Check if path already exists
				use = model_check_path(model, path);
				// If not, add as a new model
				if ( use ) {
					len = model_add_path(model, path);
					if ( verbose ) {
						cout << i+1 << ":";
						for ( k=0; k<=i; k++ ) cout << " " << path[k]->identifier();
						cout << " " << path[0]->identifier();
						cout << " " << len << endl;
					}
				}
			}
		}
	}
	
	return i;
}

/**
@brief 	Calculates the Hamiltonian cycles in a model.
@param 	*model			model structure.
@return Bmodel*					linked list of models with cycle links.

	A Hamiltonian cycle passes through each component once.
	The cycles are limited to the links in the input model.
	Redundant cycles are discarded.
	Only the first model in the linked list is used.
	Each of the resultant models has links forming a Hamiltonian cycle.
Reference:  H. Wiener, J. Am. Chem. Soc., 1947, 69, 17.

**/
Bmodel*		model_hamiltonian_cycle(Bmodel* model)
{
	int				i, n, m;
	Bcomponent*		comp;
	Blink*			link;

	for ( n=0, comp = model->comp; comp; comp = comp->next ) n++;
	for ( m=0, link = model->link; link; link = link->next ) m++;
	
	if ( verbose ) {
		cout << "Finding Hamiltonian cycles:" << endl;
		cout << "Components:                 " << n << endl;
		cout << "Starting links:             " << m << endl;
		cout << "Paths:" << endl;
	}
	
	Bcomponent**	path = new Bcomponent*[n];
	for ( i=0; i<n; i++ ) path[i] = NULL;

	comp_pick_next(0, n, model->comp, path, model);
	
	delete[] path;

	Bmodel*			modpath = model->next;
	model->next = NULL;
	
	if ( verbose )
		cout << endl;
	
	return modpath;
}

