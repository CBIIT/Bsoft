/**
@file	model_compare.cpp
@brief	Functions to compare models and components
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20161012
**/

#include "model_select.h"
#include "model_util.h"
#include "Matrix3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Compares the numbers of components in two models.
@param 	*model1			first model.
@param 	*model2			second model.
@return long			difference in the number of components.

	Only the first models in the lists are compared.

**/
long		model_component_number_difference(Bmodel* model1, Bmodel* model2)
{
	long			n1, n2;
	Bcomponent*		comp1, *comp2;
	
	for ( n1 = 0, comp1 = model1->comp; comp1; comp1 = comp1->next ) n1++;

	for ( n2 = 0, comp2 = model2->comp; comp2; comp2 = comp2->next ) n2++;

	return n1 - n2;
}

/**
@brief	Calculates the largest number of components in a model.
@param 	*model			model.
@return long			largest number of components in a model.
**/
long		model_maxnum_components(Bmodel* model)
{
	long	ncomp, n;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( ncomp = 0, mp = model; mp; mp = mp->next ) {
		for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;
		if ( ncomp < n ) ncomp = n;
	}

	return ncomp;
}

/**
@brief	Compares two models.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@return double			RMSD.

	Only the first models in the linked lists are compared.

**/
double		model_compare(Bmodel* model1, Bmodel* model2)
{
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	Bcomponent*		compsel;
	
	long	n(0), n1(0), n2(0);
	double			d, dmin, R(0);

//	Vector3<float>	com1 = model_center_of_mass(model1);
//	Vector3<float>	com2 = model_center_of_mass(model2);

//	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) comp1->location() -= com1;
//	for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) comp2->location() -= com2;

	if ( verbose ) {
		cout << "Comparing " << model1->identifier() << " with " << model2->identifier() << ":" << endl;
		cout << "Comp1\tComp2\tDmin" << endl;
	}
	
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		dmin = 1e30;
		compsel = NULL;
		for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			if ( dmin > d ) {
				dmin = d;
				compsel = comp2;
			}
		}
		if ( compsel ) {
			comp1->select(stol(compsel->identifier()));
			comp1->FOM(dmin);
		}
		n1++;
	}

	for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
		dmin = 1e30;
		compsel = NULL;
		for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
			d = comp1->location().distance(comp2->location());
			if ( dmin > d ) {
				dmin = d;
				compsel = comp1;
			}
		}
		if ( compsel ) {
			comp2->select(stol(compsel->identifier()));
			comp2->FOM(dmin);
		}
		n2++;
	}

	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = model2->comp; comp2 && comp1->select() != stoi(comp2->identifier()); comp2 = comp2->next ) ;
		if ( comp2 && comp2->select() == stoi(comp1->identifier()) ) {
			R += comp1->FOM()*comp1->FOM();
			n++;
			if ( verbose )
				cout << comp1->identifier() << tab << comp2->identifier() << tab << comp1->FOM() << endl;
		}
	}

	R = sqrt(R/n);
	
	if ( verbose ) {
		cout << "Components compared: " << n << endl;
		cout << "Model 1 components:  " << n1 << " (" << n*100.0/n1 << "%)" << endl;
		cout << "Model 2 components:  " << n2 << " (" << n*100.0/n2 << "%)" << endl;
		cout << "RMSD:                " << R << endl << endl;
	}
	
	return R;
}

/**
@brief 	Constructs the adjacency matrix for a model.
@param 	*model			model structure.
@return Matrix			adjacency matrix.

	The matrix contains ones for adjacent components and zero elsewhere. 
	Only the first model in the linked list is used.
	The component selections are reset.

**/
Matrix		model_adjacency_matrix(Bmodel* model)
{
	int				i, j, k, ncomp;
	Bcomponent*		comp;
	
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next, ncomp++ )
		comp->select(ncomp);
	
	Matrix			mat(ncomp,ncomp);

	for ( comp = model->comp; comp; comp = comp->next ) {
		i = comp->select();
		for ( k=0; k<comp->link.size() && comp->link[k]; k++ ) {
			j = comp->link[k]->select();
			mat[i][j] = 1;
		}
	}

	for ( comp = model->comp; comp; comp = comp->next )
		comp->select(1);
	
	return mat;
}

/**
@brief 	Calculates the distance matrix for a model.
@param 	*model			model structure.
@param 	view_flag		flag to include views.
@return Matrix			distance matrix.

	The matrix is calculated from the pairwise distances between selected 
	components, with the option to include the views.
	When the views are included, the euclidean distances are rescaled to
	the maximum so that their ranges are similar to the view differences. 
	Only the first model in the linked list is used.

**/
Matrix		model_distance_matrix(Bmodel* model, int view_flag)
{
	long			n, i, j;
	double			d, dmax;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( n=0, dmax = 0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		d = comp->location().length();
		if ( dmax < d ) dmax = d;
		n++;
	}
	
	if ( verbose )
		cout << "Calculating a " << n << " x " << n << " distance matrix:" << endl << endl;

	Matrix			dmat(n, n);
	
	for ( i=0, comp = model->comp; comp->next; comp = comp->next ) if ( comp->select() ) {
		for ( j=i+1, comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			if ( view_flag )
				dmat[i][j] = dmat[j][i] =
					comp->location().distance(comp2->location())/dmax + comp->view().residual(comp2->view());
			else
				dmat[i][j] = dmat[j][i] = comp->location().distance(comp2->location());
			j++;
		}
		i++;
	}

	if ( verbose & VERB_FULL )
		cout << dmat << endl;
	
	return dmat;
}

/**
@brief 	Calculates the distance matrix between two models.
@param 	*m1				first model structure.
@param 	*m2				second model structure.
@return Matrix			distance matrix.

	The matrix is calculated from the pairwise distances between selected 
	components. 
	Only the first model in each linked list is used.

**/
Matrix		model_distance_matrix(Bmodel* m1, Bmodel* m2)
{
	long			i, j, n1(0), n2(0);
	Bcomponent		*c1, *c2;
	
	n1 = m1->component_count();
	n2 = m2->component_count();
	
	Matrix			mat(n1, n2);

	for ( i=0, c1 = m1->comp; c1; c1 = c1->next, ++i )
		for ( j=0, c2 = m2->comp; c2; c2 = c2->next, ++j )
			mat[i][j] = c1->location().distance(c2->location());
	
	return mat;
}

/**
@brief 	Consolidates close components within a model.
@param 	*model			model structure.
@param 	distance		cutoff distance to consider components to be the same.
@return long			number of components retained.

	A matrix is calculated from the pairwise distances between selected components.
	Components closer to each other than the given distance are consolidated. 
	Only the first model in the linked list is used.

**/
long		model_consolidate(Bmodel* model, double distance)
{
	long			i, j, k, n(0);
	Bcomponent		*c1, *c2;
	Bcomponent*		c = NULL;
	Bcomponent*		c_list = NULL;
	string			cid;

	n = model->component_count();

	Matrix			mat = model_distance_matrix(model, 0);
	
	if ( verbose ) {
		cout << model->identifier() << tab << n << endl;
		mat.show_below_cutoff(distance);
	}
	
	for ( i=0; i<mat.rows(); ++i ) mat[i][i] = 2*distance;

	for ( i=n=0, c1 = model->comp; c1; c1 = c1->next, ++i ) if ( mat[i][i] ) {
		cid = to_string(++n);
//		c = component_add(&c_list, cid);
//		component_copy(c1, c);
		if ( c_list ) c = c_list->add(c1);
		else c = c_list = new Bcomponent(c1);
		c->identifier() = cid;
		for ( j=i+1, k=1, c2 = c1->next; c2; c2 = c2->next, ++j ) {
			if ( mat[i][j] <= distance ) {
				c->shift(c2->location());
				mat[j][j] = 0;
				k++;
			}
		}
		if ( k > 1 ) c->scale(1.0L/k);
		c->select(k);
	}
	
	component_list_kill(model->comp);
	model->comp = c_list;

	if ( verbose )
		cout << "Number of components:     " << n << endl << endl;
	
	return n;
}

/**
@brief 	Determines a consensus between models in a list.
@param 	*model			model list.
@param 	distance		cutoff distance to consider components to be the same.
@return long			number of components retained.

	A matrix is calculated from the pairwise distances between selected components
	from each apir of models.
	Components closer to each other than the given distance in different
	models are consolidated.
	The component selection field contains the number of contributing components.
	A new model containing the result is returned.

**/
Bmodel*		models_consensus(Bmodel* model, double distance)
{
	long			nmod(0), i, j, k, n;
	Bmodel			*m1, *m2;
	Bcomponent		*c1, *c2;
	
	for ( m1 = model; m1; m1 = m1->next ) nmod++;
	
	for ( m1 = model; m1; m1 = m1->next )
		model_consolidate(m1, distance);
	
	long			nmat(factorial(nmod)/(2*factorial(nmod-2)));
	
	Matrix*			mat = new Matrix[nmat];
				
	for ( i=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( m2 = m1->next; m2; m2 = m2->next, ++i ) {
			mat[i] = model_distance_matrix(m1, m2);
			if ( verbose & VERB_FULL ) {
				cout << m1->identifier() << " vs " << m2->identifier() << endl;
				mat[i].show_below_cutoff(distance);
			}
		}
	}
	
	models_process(model, model_reset_selection);

	for ( i=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( m2 = m1->next; m2; m2 = m2->next, ++i ) {
			for ( j=0, c1 = m1->comp; c1; c1 = c1->next, ++j ) if ( c1->select() ) {
//				c1->select(1);
				for ( k=0, c2 = m2->comp; c2; c2 = c2->next, ++k ) if ( c2->select() ) {
//					c2->select(1);
					if ( mat[i][j][k] <= distance ) {
						c1->shift(c2->location());
						c1->select_increment();
						c2->select(0);
					}
				}
//				c1->loc /= c1->select();
			}
		}
	}
	
	delete[] mat;
	
	Bmodel*			numod = new Bmodel(model->identifier());
	Bstring			cid;
	RGBA<float>		color;
	int*			ns = new int[nmod+1];
	for ( i=0; i<=nmod; i++ ) ns[i] = 0;
	
	numod->mapfile() = model->mapfile();
	numod->image_number(model->image_number());
	numod->comment(model->comment());
	Bcomptype*		ct = numod->add_type(model->type->identifier());

	for ( n=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( c1 = m1->comp; c1; c1 = c1->next ) if ( c1->select() ) {
			cid = Bstring(++n, "%d");
			c2 = numod->add_component(c1);
			c2->type(ct);
			c2->identifier(cid);
			c2->location(c2->location()/c2->select());
			c2->FOM(c2->select()*1.0L/nmod);
			ns[c2->select()]++;
			color.spectrum(c2->select(), nmod - 0.5, 0.5);
			c2->color(color);
			if ( verbose & VERB_FULL )
				cout << n << tab << c2->location() << tab << c2->select() << endl;
		}
	}
	
	if (verbose ) {
		cout << "Number\tCount\t%" << endl;
		for ( i=0; i<=nmod; i++ )
			cout << i << tab << ns[i] << tab << ns[i]*100.0/n << endl;
		cout << endl;
	}
	
	delete[] ns;
	
	return numod;
}



