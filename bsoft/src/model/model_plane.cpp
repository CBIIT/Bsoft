/**
@file	model_plane.cpp
@brief	Library routines used for plane models.
@author Bernard Heymann
@date	Created: 20140925
@date	Modified: 20141008
**/

#include "model_create.h"
#include "model_transform.h"
#include "model_views.h"
#include "model_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Fits a plane through a model.

	A plane is fit through the components and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	Only the first model is fit.

@param 	*model		model structure.
@return Vector3<double>		plane normal.
**/
Vector3<double>	model_fit_plane(Bmodel* model)
{
	int					i, j, n;
	double				offset, d, dstd(0), dmin(1e30);
	vector<double>		b(3);
	Matrix				a(3,3);
	Vector3<double>		cen, loc, normal;
	Bcomponent*			comp;

	for ( i=0; i<3; i++ ) b[i] = 0;
	
	for ( n=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		loc = comp->location();
		cen += loc;
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<=i; j++ )
				a[j][i] += loc[i]*loc[j];
			b[i] += loc[i];
		}
	}
	cen /= n;
	for ( i=1; i<3; i++ )
		for ( j=0; j<i; j++ )
			a[i][j] = a[j][i];
		
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

	offset = cen.scalar(normal);
	
	for ( n=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		d = comp->location().scalar(normal) - offset;
		dstd += d*d;
		if ( fabs(dmin) > fabs(d) ) dmin = d;
		comp->FOM(d);
	}

	dstd /= n;
	model->FOM(dstd);
	
	if ( dmin < 0 ) normal = -normal;
	for ( comp = model->comp; comp; comp = comp->next )
		comp->view(View2<float>(normal));
		
	return normal;
}

/**
@brief 	Calculates the distance between a location and the local guide model.

	Each guide component is used to calculate a distance of the location
	to the guide plane and a lateral distance to the guide component.
	The lateral distance is used to calculate a gaussian distributed weight
	which in turn is used to calculate a weighted average of the distance
	to the guide plane.

@param 	loc	location to be tested.
@param 	vec	direction vector associated with the location.
@param 	*guide		guide model.
@param 	sigma2		smoothness weighting parameter.
@return Bmodel*				new circle model.
**/
double		model_distance_to_guide(Vector3<double> loc, Vector3<double> vec, Bmodel* guide, double sigma2)
{
	double			dp, dl, dpa(0), w, wa(0);
	Vector3<double>	dvec;
	Bcomponent*		comp;

	for ( comp = guide->comp; comp; comp = comp->next ) {
		dvec = comp->location() - loc;
		dp = dvec.scalar(vec);			// distance to the general guide plane
		dl = dvec.length2() - dp*dp;	// lateral distance to the guide component
		w = exp(-dl/sigma2);
		dpa += dp*w;
		wa += w;
	}
	
//	cout << dpa << tab << wa << endl;
	
	if ( wa ) dpa /= wa;
	
	return dpa;
}

/**
@brief 	Adjusts a shell model to the faces of a polyhedral guide model.

	For each component, it is determined whether it is located inside or 
	outside the appropriate polyhedral face, and its coordinates adjust
	closer to the face by the indicated fraction.

@param 	*model			shell model.
@param 	*gmod			guide polyhedron model.
@param 	sigma			smoothness of model.
@return int						0.
**/
int			model_adjust_to_guide(Bmodel* model, Bmodel* gmod, double sigma)
{
	int				n;
	double			d, davg(0), dstd(0), sigma2(2*sigma*sigma);
	Vector3<double>	center;
	Bcomponent*		comp = NULL;

	if ( verbose ) {
		cout << "Adjusting the model " << model->identifier() << " to the guide " << gmod->identifier() << endl;
		cout << "Sigma:                          " << sigma << endl << endl;
	}
	
	if ( model->mapfile().length() < 1 && gmod->mapfile().length() ) model->mapfile() = gmod->mapfile();
	
	for ( n=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		comp->force(comp->view().vector3());
		d = model_distance_to_guide(comp->location(), comp->force(), gmod, sigma2);
		if ( !isfinite(d) ) cerr << "d not finite: " << d << endl;
		if ( verbose )
			cout << comp->identifier() << tab << comp->location() << tab << comp->force() << tab << d << endl;
		comp->shift(comp->force() * d);
		davg += d;
		dstd += d*d;
	}
	
	davg /= n;
	dstd /= n;
	dstd -= davg*davg;
	if ( dstd > 0 ) dstd = sqrt(dstd);
	else dstd = 0;
	
	if ( verbose )
		cout << "Average displacement:           " << davg << " (" << dstd << ")" << endl;

	return 0;
}


/**
@brief 	Generates a plane model and fits a plane through it.
@param 	*guide			model structure.
@param 	separation		separation between components for new model.
@param 	sigma			smoothness of model.
@return Bmodel* 				new model based on input model.

	The principal axes of a planar model are calculated and a new plane model
	generated with the given separation between components.
	The components of the new model are then adjusted to fit the guide model.

**/
Bmodel*		model_generate_from_plane_guide(Bmodel* guide, double separation, double sigma)
{
	if ( !guide ) {
		cerr << "Error: No guide model!" << endl;
		return NULL;
	}
	
	Vector3<double>* eigenvec = new Vector3<double>[3];
	Vector3<double>	pax = model_principal_axes(guide, eigenvec);

	Bmodel*			model = model_create_plane(2*pax[0], 2*pax[1], 0, separation);

	if ( !model ) {
		cerr << "Error: The model creation from " << guide->identifier() << " failed!" << endl;
		return NULL;
	}
	
	model->identifier() = guide->identifier();
	model->mapfile() = guide->mapfile();
	model->image_number(guide->image_number());
	model->comment(guide->comment());

	int				i;
	Vector3<double>	origin;
	Vector3<double>	shift = model_center_of_mass(guide);
	Matrix3			mat;
	
	for ( i=0; i<3; i++ ) mat[i] = eigenvec[i];
	
	model_rotate(model, mat, origin, shift);
	
	View2<float>	view(eigenvec[2][0], eigenvec[2][1], eigenvec[2][2], 0);
	model_set_views(model, view);
	
	delete[] eigenvec;

	model_adjust_to_guide(model, guide, sigma);
	
	return model;
}



