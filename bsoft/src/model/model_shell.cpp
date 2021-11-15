/**
@file	model_shell.cpp
@brief	Library routines used for shell model processing
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20161013
**/

#include "model_create.h"
#include "model_shell.h"
#include "model_color.h"
#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Adds a shell to a model using selected components.
@param 	*model			model structure.
@param 	add_distance	distance to shift new components.
@param 	&nutype			new component type.
@return long			number of components in new shell.

	Each selected component is duplicated and shifted radially.
	Only the first model in the list is used.

**/
long		model_add_shell(Bmodel* model, double add_distance, Bstring& nutype)
{
	if ( !model ) return -1;
	if ( !model->comp ) return 0;
	
	Bcomptype*		ct = model->add_type(nutype);
	
	long			i, ncomp;
	Vector3<double>	vv;
	Bcomponent*		comp;
	Bcomponent*		comp_new = NULL;
	Bcomponent*		comp_list = NULL;

	for ( i=0, comp = model->comp; comp; comp = comp->next )
		if ( i < stol(comp->identifier()) ) i = stol(comp->identifier());
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		if ( comp_new ) comp_new = comp_new->add(comp);
		else comp_new = comp_list = new Bcomponent(comp);
		comp_new->identifier() = to_string(++i);
		comp_new->type(ct);
		vv = comp->view().vector3();
		comp_new->shift(vv * add_distance);
		ncomp++;
	}

	for ( comp = model->comp; comp->next; comp = comp->next ) ;
	if ( comp ) comp->next = comp_list;
	
	if ( verbose ) {
		cout << "Generated new shell" << endl;
		cout << "New type:                       " << nutype << endl;
		cout << "Distance added:                 " << add_distance << " A" << endl;
		cout << "Components:                     " << ncomp << endl << endl;
	}
	
	return ncomp;
}

/**
@brief 	Adjusts a shell model to the faces of a polyhedral guide model.
@param 	*model			shell model.
@param 	*gmod			guide polyhedron model.
@param 	fraction		fraction to adjust.
@param 	curv_flag		flag to indicate curved surface.
@return int				0.

	For each component, it is determined whether it is located inside or 
	outside the appropriate polyhedral face, and its coordinates adjust
	closer to the face by the indicated fraction.

**/
int			model_adjust_shell_to_guide(Bmodel* model, Bmodel* gmod, double fraction, int curv_flag)
{
	if ( fraction < 1e-3 || fraction > 1 ) fraction = 1;
	
	int				n;
	double			d;
	Vector3<double>	center;
	Bcomponent*		comp = NULL;
	Bcomponent*		gcomp = NULL;

	// Calculate the center of the vertices in the guide model
	for ( n=0, gcomp = gmod->comp; gcomp; gcomp = gcomp->next, n++ )
		center += gcomp->location();
	
	if ( n ) center /= n;
	
	if ( verbose ) {
		cout << "Adjusting the model " << model->identifier() << " to the guide polyhedron " << gmod->identifier() << endl;
		cout << "Fraction:                       " << fraction << endl;
		cout << "Guide model center:             " << center << endl << endl;
	}
	
	if ( model->mapfile().length() < 1 && gmod->mapfile().length() ) model->mapfile() = gmod->mapfile();
	
	// Calculate the component vectors in the guide model
	for ( gcomp = gmod->comp; gcomp; gcomp = gcomp->next )
		gcomp->force(gcomp->location());
//		gcomp->force(gcomp->location() - center);
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		d = fraction * model_inside_outside(comp->location(), gmod, curv_flag, 0);
		if ( !isfinite(d) ) cerr << "d not finite: " << d << endl;
		d += comp->location().normalize();
		comp->scale(d);
//		comp->shift(center;
	}

	return 0;
}

/**
@brief 	Generates a new set of models by vonverting each component to a shell model.
@param 	*model			model structure.
@param 	distance		distance between components.
@param 	&nutype			new component type.
@param	twod			flag to indicate a 2D model.
@return Bmodel*			new list of models.

	Each selected component is converted to a shell model.
	The radius is taken from the original component radius.
	The given distance defines the distance between the new
	components and their radii.

**/
Bmodel*		model_components_to_shells(Bmodel* model, double distance, Bstring& nutype, int twod)
{
	if ( !model ) return NULL;
	if ( !model->comp ) return NULL;
	
	if ( distance > 1.5*model->comp->radius() ) distance = 1.5*model->comp->radius();
	if ( nutype.length() < 1 ) nutype = "VER";
	
	if ( verbose ) {
		cout << "Generating new shell models from components:" << endl;
		cout << "Number of components:           " << model->component_count() << endl;
		cout << "New type:                       " << nutype << endl;
		cout << "Distance between components:    " << distance << " A" << endl;
		if ( twod ) cout << "Generating plane models" << endl;
	}
	
	long			ncomp(0);
	Bmodel*			mp = NULL;
	Bmodel*			numod = NULL;
	Bmodel*			numod_list = NULL;
	Bcomponent*		comp = NULL;
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( twod )
				numod = model_create_circle(comp->radius(), 0, distance);
			else
				numod = model_create_shell(0, comp->radius(), distance);
//			numod = model_add(&numod_list, numod);
			if ( numod_list ) numod_list->add(numod);
			else numod_list = numod;
			numod->identifier() = mp->identifier() + "_" + comp->identifier();
			model_set_type(numod, nutype);
			model_set_component_radius(numod, distance/2);
			model_shift(numod, comp->location());
			numod->mapfile() = mp->mapfile();
			numod->image_number(mp->image_number());
			model_color_uniformly(numod, comp->color());
			ncomp += numod->component_count();
		}
	}
	
	if ( verbose )
		cout << "New components:                 " << ncomp << endl << endl;
	
	return numod_list;
}

/**
@brief     Calculates how close the model is to a spherical shape.
@param 	*model				model structure.
@return double				sphericity.

	The deviation of vertices from the average radius is calculated.
	Only the first model in the list is processed.

**/
double		model_sphericity(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	int				ncomp;
	double			d, sph = 0, da = 0, ds = 0;	

	Vector3<double>	com = model_center_of_mass(model);
	
	Bcomponent*		comp;
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		d = (comp->location() - com).length();
		da += d;
		ds += d*d;
	}
	
	if ( ncomp ) {
		da /= ncomp;
		ds = ds/ncomp - da*da;
		if ( ds > 0 ) sph = 1 - sqrt(ds)/da;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Sphericity:                    " << sph << endl << endl;
	
	return sph;
}

/**
@brief     Calculates how close the model is to an ellipsoid shape.
@param 	*model				model structure.
@return double				ellipsoidicity.

	The principal axes are calculated first.
	Then the fit to an ellipsoid function is calculated.
	Only the first model in the list is processed.

**/
double		model_ellipsoidicity(Bmodel* model)
{
	if ( !model ) return 0;
	if ( !model->select() ) return 0;
	
	int				i, ncomp;
	double			d, ell = 0;
	Vector3<double>	c, c1, axis[3];	
	Vector3<double> 	alen = model_principal_axes(model, axis);
	
	for ( i=0; i<3; i++ ) axis[i].normalize();

	Vector3<double>	com = model_center_of_mass(model);
	
	Bcomponent*		comp;
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		c = comp->location() - com;
		for ( i=0; i<3; i++ ) c1[i] = c.scalar(axis[i]);
		c1 /= alen;
		d = 1 - c1.length2();
		ell += d*d;
	}
	
	ell = sqrt(1 - ell/ncomp);
	
	if ( verbose & VERB_PROCESS )
		cout << "Ellipsoidicity:                 " << ell << endl << endl;
	
	return ell;
}

/**
@brief 	Calculates the curvature associated with each link.
@param 	*model			model structure (views modified).
@return int				0.

	The difference angle in the normals or view vectors of each pair of
	vertices of a link represents the curvature for that link.
	Curvature at each component is calculated as the average of the 
	curvature of the attached links and set as the component FOM.
	The model curvature is calculated as the average component curvature
	and set as the model FOM.

**/
int			model_curvature(Bmodel* model)
{
	int				i, ncomp;
	double			amax;
	Bmodel*			mp;
	Bcomponent*		comp;
	Blink*			link;
	
	if ( verbose )
		cout << "Calculating curvature" << endl << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Model\tAverage\tMaximum" << endl;
	for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		mp->FOM(0);
		amax = 0;
		for ( comp = mp->comp; comp; comp = comp->next ) comp->FOM(0);
		for ( link = mp->link; link; link = link->next ) {
			link->FOM(link->comp[0]->view().angle(link->comp[1]->view()));
			link->comp[0]->FOM(link->comp[0]->FOM() + link->FOM());
			link->comp[1]->FOM(link->comp[1]->FOM() + link->FOM());
			if ( amax < link->FOM() ) amax = link->FOM();
		}
		for ( ncomp=0, comp = mp->comp; comp; comp = comp->next, ncomp++ ) {
			for ( i=0; i<comp->link.size() && comp->link[i]; i++ ) ;
			if ( i ) comp->FOM(comp->FOM() / i);
			mp->FOM(mp->FOM() + comp->FOM());
		}
		if ( ncomp ) mp->FOM(mp->FOM() / ncomp);
		if ( verbose & VERB_PROCESS )
			cout << mp->identifier() << tab << mp->FOM()*180.0/M_PI << tab << amax*180.0/M_PI << endl;
		mp->FOM(amax);
	}	
	
	return 0;
}

/**
@brief 	Determines if a vector is inside or outside model boundaries.
@param 	vec				vector to test for.
@param 	*model			model structure.
@param 	curv_flag		flag to indicate curved surface.
@return double			distance from interface (inside positive).

	Whether a point is inside or outside is based on the model structure, 
	where the intersection of a vector through the point with the plane through 
	the three closest vertices is the dividing point. The following set of 
	equations is solved:
		v/|v| = w0 * a0 + w1 * a1 + w2 * a2
	where v is the point to be judged, and w<n> are the weights for the 
	contributing vectors a<n>.
	The distance to the triangle plane intersecting point is:
		d = 1/(w0 + w1 + w2)
	The distance to the curved surface through the three vertices is:
		d = w0 * |a0| + w1 * |a1| + w2 * |a2|
	where the weights have been normalized:
		w0 + w1 + w2 = 1
	Requirement: The sampled component location vectors must be set up 
		in the component force vectors.
	Only the first model in the linked list is used.

**/
double		model_inside_outside_fast(Vector3<double> vec, Bmodel* model, int curv_flag)
{
	int				i, ib, j;
	double			testangle, d, angle[3] = {10,10,10};
	Vector3<double>	tvec[3], vec3(vec[0], vec[1], vec[2]);
	Matrix			mat(3,3);
	vector<double>	w(3);
	Bcomponent*		comp;

	double			vec_len = vec3.normalize();
	
	// Find the closest three vertices
	for ( ib = 0, comp = model->comp; comp && ib < 3; comp = comp->next ) {
		testangle = 0;
		for ( i=ib=0; i<3; i++ ) if ( testangle < angle[i] ) {	// Find the biggest angle
			testangle = angle[i];
			ib = i;
		}
		testangle = vec3.angle(comp->force());
		if ( testangle < 1e-30 ) {
			return comp->force().length() - vec_len;
		}
		if ( angle[ib] > testangle ) {	// Test to replace the biggest angle
			angle[ib] = testangle;
			tvec[ib] = comp->force();
		}
	}
	
	for ( i=0; i<3; i++ ) {
		for ( j=0; j<3; j++ )
			mat[i][j] = tvec[j][i];
		w[i] = vec3[i];
	}

//	double		det = mat.LU_decomposition(w);
	double		det = mat.singular_value_decomposition(w);
	if ( det < -1e20 ) {
		cerr << "Error in model_inside_outside_fast: vec=" << vec << " det=" << det << endl;
		cerr << mat << endl;
		return 0;
	}
	
	if ( w[0]<0 || w[1]<0 || w[2]<0 )
		if ( verbose & VERB_DEBUG )
			cerr << "Warning: Negative weights: " << 
				w[0] << " " << w[1] << " " << w[2] << " " << vec << endl;
	
//	for ( i=0, d=0; i<3; i++ ) if ( w[i] > 0 ) d += w[i]; else w[i] = 0;
	for ( i=0, d=0; i<3; i++ ) d += w[i];
	
	if ( curv_flag ) {
		if ( d ) for ( i=0; i<3; i++ ) w[i] /= d;
		for ( i=0, d=0; i<3; i++ ) d += w[i] * tvec[i].length();
	} else if ( d ) {
		d = 1/d;
	}
	
	return d - vec_len;
}

double		model_inside_outside(Vector3<double> vec, Bmodel* model, int curv_flag, int fast)
{
	if ( fast ) return model_inside_outside_fast(vec, model, curv_flag);
	
	int				i, j, iw, ib;
	double			a, d, dmin = 1e30;
//	double			w[3], wb[3] = {0,0,0};
	vector<double>	w(3), wb(3,0);
	Matrix			mat(3,3);
	Vector3<double>	tvec[3], bvec[3], vec3(vec[0], vec[1], vec[2]);
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	Bcomponent*		comp3;
	
	double			vec_len = vec3.normalize();
	
	if ( !isfinite(vec_len) )
		cerr << "Vector not finite: " << vec << " length=" << vec_len << endl;
	
	for ( ib = 0, a = M_PI/16; ib < 3 && a < M_PI + 0.1; a *= 2 ) {
		for ( comp1 = model->comp; comp1; comp1 = comp1->next ) {
			tvec[0] = comp1->force();
			if ( vec3.angle(tvec[0]) < a ) {
				for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) {
					tvec[1] = comp2->force();
					if ( vec3.angle(tvec[1]) < a ) {				
						for ( comp3 = comp2->next; comp3; comp3 = comp3->next ) {
							tvec[2] = comp3->force();
							if ( vec3.angle(tvec[2]) < a ) {
								for ( i=0; i<3; i++ ) {
//									cout << tvec[i] << endl;
									for ( j=0; j<3; j++ )
										mat[i][j] = tvec[j][i];
									w[i] = vec3[i];
								}
								if ( mat.check_for_singularity() < 1 ) {
									mat.LU_decomposition(w);
									for ( i=iw=0; i<3; i++ ) if ( w[i] >= 0 ) iw++;
									if ( iw == 3 ) {
										for ( i=0, d=0; i<3; i++ ) d += w[i];
										if ( d > 1e-20 && isfinite(d) ) {
											if ( dmin > d ) {
												dmin = d;
												for ( i=0; i<3; i++ ) {
													wb[i] = w[i];
													bvec[i] = tvec[i];
												}
												ib = iw;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	if ( !isfinite(wb[0]) || !isfinite(wb[1]) || !isfinite(wb[2]) )
		cerr << "Warning: Non-finite weights: " << wb[0] << " " << wb[1] << " " << wb[2] << endl;
	
	if ( wb[0]<0 || wb[1]<0 || wb[2]<0 )
		cerr << "Warning: Negative weights: " << wb[0] << " " << wb[1] << " " << wb[2] << endl;
	
	for ( i=0, d=0; i<3; i++ ) d += wb[i];
	
	if ( curv_flag ) {
		if ( d ) for ( i=0; i<3; i++ ) wb[i] /= d;
		for ( i=0, d=0; i<3; i++ ) d += wb[i] * bvec[i].length();
	} else if ( d ) {
		d = 1/d;
	} else cout << "Distance = " << d << endl;
	
	if ( !isfinite(d) )
		cerr << "Weights and distance: " << wb[0] << " " << wb[1] << " " << wb[2] << " " << d << endl;
	
	return d - vec_len;
}



