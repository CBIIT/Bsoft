/**
@file	model_neighbors.cpp
@brief	Functions to manipulate model component neighbors
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20210324
**/

#include "model_neighbors.h"
#include "model_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Sets the desired number of components as neighbors for each component.
@param 	*model		model.
@param 	number		number of neighbors to set.
@return int			0.

	For each component, its neighbor list is set to the closest components.
	Only the first model in the linked list is used.

**/
int			model_set_neighbors(Bmodel* model, int number)
{
	if ( number > MAXLINK ) number = MAXLINK;
	
	int				i, imax(0);
	double			d1, dmax;
	vector<double>	d(number);
	Bcomponent*		comp;
	Bcomponent*		comp2;

	for ( comp = model->comp; comp; comp = comp->next ) {
		for ( i=0; i<number; ++i ) d[i] = 1e30;
		for ( comp2 = model->comp; comp2; comp2 = comp2->next ) if ( comp != comp2 ) {
			for ( i=0, dmax = 0; i<number; ++i ) if ( dmax < d[i] ) {
				dmax = d[i];
				imax = i;
			}
			d1 = comp->location().distance(comp2->location());
			if ( d1 < dmax ) {
				d[imax] = d1;
				if ( comp->neighbor.size() <= imax )
					comp->neighbor.push_back(comp2);
				else
					comp->neighbor[imax] = comp2;
			}
		}
	}

	model_neighbor_reciprocity(model);

	return 0;
}

int			model_set_neighbors(Bmodel* model, int number, double distance)
{
	if ( number > MAXLINK ) number = MAXLINK;
	
	int				i, j, imin, ncomp(0), nd(0);
	vector<int>		n(number+1);
	double			a, d1, dmin, dist(0);
	Vector3<double>	vz(0,0,1), v0[MAXLINK], v[MAXLINK], vc;
	vector<double>	d(number);
	Matrix3			mat;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	Vector3<double>	com = model_center_of_mass(model);

	if ( verbose ) {
		cout << "Setting up neighbors for model: " << model->identifier() << endl;
		cout << "Number:                         " << number << endl;
		cout << "Distance:                       " << distance << " A" << endl << endl;
	}
	
	// Set up the guide positions
	for ( i=0; i<number; i++ ) {
		a = TWOPI*i*1.0L/number;
		v0[i][0] = distance*cos(a);
		v0[i][1] = distance*sin(a);
	}

	for ( i=0; i<=number; i++ ) n[i] = 0;

	for ( comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		vc = comp->location() - com;
		mat = Matrix3(vz, vc);
		for ( i=0; i<number; i++ ) d[i] = 1e30;
		for ( i=0; i<number; i++ )	// Transform the position vectors for this component
			v[i] = mat * v0[i] + comp->location();
		for ( comp2 = model->comp; comp2; comp2 = comp2->next ) if ( comp != comp2 ) {
			for ( i=imin=0, dmin=1e30; i<number; i++ ) {	// Find the closest guide position
				d1 = comp2->location().distance(v[i]);
				if ( dmin > d1 ) {
					dmin = d1;
					imin = i;
				}
			}
			d1 = comp->location().distance(comp2->location());
			if ( d1 < d[imin] ) {
				d[imin] = d1;
				if ( comp->neighbor.size() <= imin )
					comp->neighbor.push_back(comp2);
				else
					comp->neighbor[imin] = comp2;
			}
		}
		for ( i=j=0; i<comp->neighbor.size(); i++ ) {	// Shift all neighbor references to the head of the array
			if ( comp->neighbor[i] ) {
				if ( i > j ) {
					comp->neighbor[j] = comp->neighbor[i];
					comp->neighbor[i] = NULL;
				}
				dist += comp->location().distance(comp->neighbor[j]->location());
				nd++;
				j++;
			}
		}
		n[j]++;
	}
	
	dist /= nd;
	
	if ( verbose ) {
		cout << "Neighbor numbers assigned:\nNumber\tCount" << endl;
		for ( i=0; i<=number; i++ )
			cout << i << tab << n[i] << endl;
		cout << "Average distance:               " << dist << " A" << endl << endl;
	}

	model_neighbor_reciprocity(model);

	return 0;
}

/**
@brief 	Counts how many neighbors of a component also has that component as a neighbor.
@param 	*model		model.
@return int			number of non-reciprocated neighbors.

	For each neighbor of component, the neighbor list of the neighbor is searched
	for the id of the first component.
	Only the first model in the linked list is used.

**/
int			model_neighbor_reciprocity(Bmodel* model)
{
	int				i, j, nr(0), nn(0);
	Bcomponent*		comp, *compnb;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_neighbor_reciprocity: neigbor list:" << endl;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		for ( i = 0; i<comp->neighbor.size(); i++ ) {
			compnb = comp->neighbor[i];
			for ( j=0; j<compnb->neighbor.size(); j++ )
				if ( compnb->neighbor[j]->identifier() == comp->identifier() ) break;
			if ( j<compnb->neighbor.size() ) nr++;
			else nn++;
			if ( verbose & VERB_DEBUG )
				cout << comp->identifier() << " - " << compnb->identifier() << endl;
		}
	}
	
	if ( verbose ) {
		cout << "Reciprocated neighbors:         " << nr << endl;
		cout << "Non-reciprocated neighbors:     " << nn << endl << endl;
	}
	
	return nn;
}



