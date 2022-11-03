/**
@file	model_transform.cpp
@brief	Library routines used for model transformation
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20220215
**/

#include "model_transform.h"
#include "model_select.h"
#include "model_compare.h"
#include "model_util.h"
#include "Transform.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Shifts one models to its center of mass.
@param 	*model		model parameters.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_center(Bmodel* model)
{
	Vector3<double>	com = model_center_of_mass(model);

	return model_shift(model, -com);
}

/**
@brief 	Shifts a model.
@param 	*model	model parameters.
@param 	shift	translation vector.
@return long	number of components processed.

	Only the first model in the list is processed.

**/
long		model_shift(Bmodel* model, Vector3<double> shift)
{
	long			ncomp(0);
	Bcomponent*		comp;

	if ( verbose & VERB_FULL )
		cout << "Shifting " << model->identifier() << " by " << shift << endl;
	
	for ( comp = model->comp; comp; comp = comp->next, ncomp++ )
		comp->shift(shift);
	
	return ncomp;
}

/**
@brief 	Shifts all models.
@param 	*model			model parameters.
@param 	shift	translation vector.
@return long					number of components processed.

	All models in the list are processed.

**/
long		models_shift(Bmodel* model, Vector3<double> shift)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += mp->shift(shift);
	
	return ncomp;
}

/**
@brief 	Scales a model.
@param 	*model			model parameters.
@param 	scale			scale.
@param 	origin			model origin.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_scale(Bmodel* model, Vector3<double> scale, Vector3<double> origin)
{
	if ( scale == 1 ) return 0;
	
	long			ncomp(0);
	Bcomponent*		comp;
	Blink*			link;

	if ( verbose & VERB_FULL )
		cout << "Scaling " << model->identifier() << " by " << scale << endl;
	
	for ( link = model->link; link; link = link->next )
			link->length(link->length()/link->comp[0]->location().distance(link->comp[1]->location()));
	
	for ( comp = model->comp; comp; comp = comp->next, ncomp++ )
			comp->location((comp->location() - origin) * scale + origin);
	
	for ( link = model->link; link; link = link->next )
			link->length(link->length()*link->comp[0]->location().distance(link->comp[1]->location()));
	
	return ncomp;
}

/**
@brief 	Scales a model.
@param 	*model			model parameters.
@param 	scale	scale.
@param 	origin	model origin.
@return long					number of components processed.
**/
long		models_scale(Bmodel* model, Vector3<double> scale, Vector3<double> origin)
{
	if ( scale == 1 ) return 0;
	
	long			ncomp(0);
	Bmodel*			mp;

	if ( verbose & VERB_FULL )
		cout << "Scaling:                        " << scale << endl << endl;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_scale(mp, scale, origin);
	
	return ncomp;
}

/**
@brief	Reflects the model through a mirror plane.
@param 	*model			model structure.
@param 	normal	plane normal.
@param 	origin	model origin.
@return long					number of components processed.

	Only the first model in the list is processed.

**/
long		model_reflect(Bmodel* model, Vector3<double> normal, Vector3<double> origin)
{
	long			ncomp(0);
	double			d;
	Bcomponent*		comp;
	
	normal.normalize();
	
	if ( verbose & VERB_FULL )
		cout << "Reflecting " << model->identifier() << " through " << normal << endl;
	
	for ( comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		comp->location() -= origin;
		d = -2 * comp->location().scalar(normal);
		comp->shift(normal * d);
		comp->shift(origin);
	}

	return ncomp;
}

/**
@brief	Reflects the model through a mirror plane.
@param 	*model			model structure.
@param 	normal	plane normal.
@param 	origin	model origin.
@return long				number of components processed.
**/
long		models_reflect(Bmodel* model, Vector3<double> normal, Vector3<double> origin)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_scale(mp, normal, origin);

	return ncomp;
}

/**
@brief	Copies and reflects the model and compares it with the original.
@param 	*model			model structure.
@param 	normal	plane normal.
@param 	origin	model origin.
@return double					RMSD.

	Only the first model in the list is processed.

**/
double		model_reflect_and_compare(Bmodel* model, Vector3<double> normal, Vector3<double> origin)
{
	Bcomponent*		comp;
	Bcomponent*		complist = NULL;
	Bcomponent*		comp2 = NULL;
	
	long	n = 0;
	double			d, dmin, R = 0;

	normal.normalize();

	for ( comp = model->comp; comp; comp = comp->next ) {
		comp2 = (Bcomponent *) add_item((char **) &comp2, sizeof(Bcomponent));
		if ( !complist ) complist = comp2;
		comp2->location(comp->location() - origin);
		d = -2 * comp2->location().scalar(normal);
		comp2->shift(normal * d);
		comp2->shift(origin);
	}

	for ( comp = model->comp; comp; comp = comp->next ) {
		dmin = 1e30;
		for ( comp2 = complist; comp2; comp2 = comp2->next ) {
			d = comp->location().distance(comp2->location());
			if ( dmin > d ) dmin = d;
		}
		R += dmin*dmin;
		n++;
	}

	R = sqrt(R/n);
	
	kill_list((char *) complist, sizeof(Bcomponent));
	
	return R;
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@param 	origin			rotation origin.
@param 	shift			translation after rotation.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin, Vector3<double> shift)
{
	long			ncomp(0);
	Bcomponent*		comp;
	
	if ( verbose & VERB_FULL ) {
		cout << "Rotating " << model->identifier() << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Origin:                         " << origin << endl;
		cout << "Translation:                    " << shift << endl << endl;
	}
	
	Matrix3			vmat;
	
	for ( comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		comp->location(comp->location() - origin);
		comp->location(mat * comp->location());
		comp->location(comp->location() + shift);
		vmat = comp->view().matrix();
		vmat = vmat * mat;
		comp->view(View2<float>(vmat));
//		comp->view((View2<float>(vmat)).backward());
//		cout << comp->view << endl;
	}

	return 0;
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@param 	origin			rotation origin.
@param 	shift			translation after rotation.
@return long				number of components processed.
**/
long		models_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin, Vector3<double> shift)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_rotate(mp, mat, origin, shift);

	return ncomp;
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, Matrix3 mat)
{
	Vector3<double>		origin, shift;
	return model_rotate(model, mat, origin, shift);
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@return long				number of components processed.
**/
long		models_rotate(Bmodel* model, Matrix3 mat)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_rotate(mp, mat);

	return ncomp;
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@param 	origin			rotation origin.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin)
{
	Vector3<double>		shift;
	return model_rotate(model, mat, origin, shift);
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	mat				rotation matrix.
@param 	origin			rotation origin.
@return long				number of components processed.
**/
long		models_rotate(Bmodel* model, Matrix3 mat, Vector3<double> origin)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_rotate(mp, mat, origin);

	return ncomp;
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	view			view to rotate to.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, View2<float> view)
{
	Vector3<double>		origin, shift;
	Matrix3				mat = view.matrix();
	return model_rotate(model, mat, origin, shift);
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	view			view to rotate to.
@return long				number of components processed.
**/
long		models_rotate(Bmodel* model, View2<float> view)
{
	long			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_rotate(mp, view);

	return ncomp;
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	view			view to rotate to.
@param 	origin			rotation origin.
@param 	shift			translation after rotation.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, View2<float> view, Vector3<double> origin, Vector3<double> shift)
{
	Matrix3				mat = view.matrix();
	return model_rotate(model, mat, origin, shift);
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	view			view to rotate to.
@param 	origin			rotation origin.
@param 	shift			translation after rotation.
@return long				number of components processed.
**/
long		models_rotate(Bmodel* model, View2<float> view, Vector3<double> origin, Vector3<double> shift)
{
	Matrix3				mat = view.matrix();
	return models_rotate(model, mat, origin, shift);
}

/**
@brief	Rotates the model.
@param 	*model			model structure.
@param 	t				rotation operation.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		model_rotate(Bmodel* model, Transform t)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating:" << endl;
		cout << "Axis:                           " << t.axis << endl;
		cout << "Angle:                          " << t.angle*180.0/M_PI << endl;
		cout << "Origin:                         " << t.origin << endl;
		cout << "Translation:                    " << t.trans << endl << endl;
	}
	
	Matrix3			mat = Matrix3(t.axis, t.angle);
	Vector3<double>	shift = t.origin + t.trans;
	
//	cout << mat << endl;

	return model_rotate(model, mat, t.origin, shift);
}

/**
@brief	Rotates the models.
@param 	*model			model structure.
@param 	t				rotation operation.
@return long				number of components processed.

	Only the first model in the list is processed.

**/
long		models_rotate(Bmodel* model, Transform t)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating:" << endl;
		cout << "Axis:                           " << t.axis << endl;
		cout << "Angle:                          " << t.angle*180.0/M_PI << endl;
		cout << "Origin:                         " << t.origin << endl;
		cout << "Translation:                    " << t.trans << endl << endl;
	}
	
	Matrix3			mat = Matrix3(t.axis, t.angle);
	Vector3<double>	shift = t.origin + t.trans;
	
//	cout << mat << endl;

	return models_rotate(model, mat, t.origin, shift);
}

/**
@brief	Copies and rotates the model and compares it with the original.
@param 	*model			model structure.
@param 	t				rotation operation.
@return double					RMSD.

	Only the first model in the list is processed.

**/
double		model_rotate_and_compare(Bmodel* model, Transform t)
{
	Bcomponent*		comp;
	Bcomponent*		complist = NULL;
	Bcomponent*		comp2 = NULL;
	
	long	n = 0;
	double			d, dmin, R = 0;

	Matrix3			mat = Matrix3(t.axis, t.angle);
	Vector3<double>	shift = t.origin + t.trans;

	for ( comp = model->comp; comp; comp = comp->next ) {
		comp2 = (Bcomponent *) add_item((char **) &comp2, sizeof(Bcomponent));
		if ( !complist ) complist = comp2;
		comp2->location(comp->location() - t.origin);
		comp2->location(mat * comp2->location());
		comp2->location(comp2->location() + shift);
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		dmin = 1e30;
		for ( comp2 = complist; comp2; comp2 = comp2->next ) {
			d = comp->location().distance(comp2->location());
			if ( dmin > d ) dmin = d;
		}
		R += dmin*dmin;
		n++;
	}

	R = sqrt(R/n);
	
	kill_list((char *) complist, sizeof(Bcomponent));
	
	return R;
}

/**
@brief 	A binned model list is converted back to the unbinned version.
@param 	*model		the model list.
@param 	bin	3-value bin vector.
@return long				number of components processed.

	All component coordinates are converted back to the unbinned values.
	Only the first model is adjusted.

**/
long		model_adjust_for_binning(Bmodel* model, Vector3<long> bin)
{
	if ( bin[0]*bin[1]*bin[2] < 2 ) return 0;
	
	long			ncomp(0);
	Bcomponent*		comp = NULL;
	Vector3<double>	binf(bin[0], bin[1], bin[2]);
	
	if ( verbose )
		cout << "Unbinning:      " << binf << endl << endl;
	
	for ( comp = model->comp; comp; comp = comp->next, ncomp++ )
		comp->scale(binf);
	
	return ncomp;
}

/**
@brief 	A model is oriented to coincide with a guide model.
@param 	*model		model (rotated).
@param 	*guide		guide model.
@return long				number of components processed.

	The principal axes of the guide model is calculated and the model rotated.
	Only the first model is aligned.

**/
long		model_align_to_guide(Bmodel* model, Bmodel* guide)
{
	Vector3<double>	eigenvec[3];	
	
	model_principal_axes(guide, eigenvec);
	
	Matrix3		mat(eigenvec[0][0], eigenvec[0][1], eigenvec[0][2], 
					eigenvec[1][0], eigenvec[1][1], eigenvec[1][2], 
					eigenvec[2][0], eigenvec[2][1], eigenvec[2][2]);
	
	mat = mat.transpose();

	if ( verbose )
		cout << mat << endl;
	
	return model_rotate(model, mat);
}

/**
@brief 	A model is fitted to a reference model.
@param 	*model			model.
@param 	*refmod			template model.
@return Transform			transform.

	The components in the model and the reference must match exactly.
	Only the first model and template in the lists are processed.

**/
Transform	model_find_transform(Bmodel* model, Bmodel* refmod)
{
	Transform		t;
		
	if ( model_component_number_difference(model, refmod) ) return t;

	if ( verbose & VERB_PROCESS )
		cout << "Mapping model " << model->identifier() << " to reference " << refmod->identifier() << endl;
	
	long			i, j;
	Bcomponent*		comp, *compr;
	vector<double>	bx(4), by(4), bz(4), v(4,0);
	Matrix			a(4,4);
	
	v[3] = 1;
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	Vector3<double>	com = model_center_of_mass(model);
	Vector3<double>	comr = model_center_of_mass(refmod);
	Vector3<double>	loc, locr;
	
	for ( comp = model->comp, compr = refmod->comp; comp && compr; comp = comp->next, compr = compr->next ) {
		if ( verbose & VERB_FULL )
			cout << "Mapping component " << comp->identifier() << " to reference " << compr->identifier() << endl;
		loc = comp->location() - com;
		locr = compr->location() - comr;
		if ( verbose & VERB_FULL )
			cout << loc[2] << tab << locr[2] << tab << loc[2] - locr[2] << endl;
		v[0] = loc[0];
		v[1] = loc[1];
		v[2] = loc[2];
		for ( i=0; i<4; i++ ) {
			bx[i] += v[i]*locr[0];
			by[i] += v[i]*locr[1];
			bz[i] += v[i]*locr[2];
			for ( j=0; j<=i; j++ ) a[i][j] += v[i]*v[j];
		}
	}

	t = transform_matrix_solve(a, bx, by, bz, 0);
	t.trans = comr - com;
	t.origin = com;
		
	return t;
}

/**
@brief	Applies random displacements to a selected number of coordinates.
@param 	*model 		model to be modified.
@param 	number		number of coordinates to displace.
@param 	stdev		standard deviation of displacement.
@return int			0.
**/
int			model_random_displace_number(Bmodel* model, long number, double stdev)
{
	if ( stdev < 1e-30 ) return 0;
	
	model_select_random(model, number);

	Bmodel*				mp;
	Bcomponent*			comp;
	
	if ( verbose & VERB_PROCESS )
		cout << "Randomizing " << number << " coordinates to standard deviation " << stdev << " angstrom" << endl << endl;
	
    for ( mp = model; mp; mp = mp->next ) if ( mp->select() )
		for( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() )
			comp->location() += vector3_random_gaussian(0, stdev);
    
    model->select_all();
	
	return 0;
}
