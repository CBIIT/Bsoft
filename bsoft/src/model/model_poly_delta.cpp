/**
@file	model_poly_delta.cpp
@brief	Functions to generate deltagraph polyhedra.
@author Bernard Heymann
@date	Created: 20080103
@date	Modified: 20180327
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_links.h"
#include "model_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates a new tubular deltagraph.
@param 	radius			radius in integer units.
@param 	height			height in integer units.
@return Bmodel*				new model.
**/
Bmodel*		model_delta_create_tube_old(int radius, int height)
{
	string			id("Tube");
	Bmodel*			model = new Bmodel(id);
	
	Bstring			hextype("HEX");
	int				i, nv;
	double			a, sa, ea;
	int				c = (int) (TWOPI*radius + 0.5);
	double			rad = c*1.0L/TWOPI;
	double			da = 1.0L/rad;
	Bcomponent*  	comp = NULL;
	Vector3<double>	loc;	
	
	if ( verbose & VERB_LABEL )
		cout << "Generating a tubular deltagraph" << endl;
		
	if ( verbose & VERB_PROCESS ) {
		cout << "Radius:                         " << rad << endl;
		cout << "Height:                         " << height << endl;
	}
	
	for ( nv=0, i=0; i<=height; i++ ) {
		loc[2] = (i-height/2.0)*sqrt(3.0)/2.0;
		sa = 0;
		if ( i%2 ) sa = da/2;
		ea = TWOPI - da/2 + sa;
		for ( a=sa; a<ea; a+=da ) {
			loc[0] = radius*cos(a);
			loc[1] = radius*sin(a);
//			id = Bstring(++nv, "%d");
//			comp = component_add(&comp, id);
//			if ( !model->comp ) model->comp = comp;
			if ( comp ) comp = comp->add(++nv);
			else model->comp = comp = new Bcomponent(++nv);
//			comp->type = model_add_type_by_id(model, hextype);
			comp->type(model->add_type(hextype));
			comp->location(loc);
			comp->radius(0.1);
		}
	}

	model_link_list_generate(model, 1.2);
	
	return model;
}

/**
@brief 	Generates a new tubular deltagraph.
@param 	h				first lattice parameter, must be 1+.
@param 	k				second lattice parameter, must be 0+.
@param 	height			height in integer units.
@return Bmodel*			new model.

	The lattice parameters determine the radius.

**/
Bmodel*		model_delta_create_tube(int h, int k, int height)
{
	if ( height < 1 ) height = 1;
	
	Bstring			id("Tube");
	Bmodel*			model = new Bmodel(id);
	
	Bstring			hextype("HEX");
	int				nv(0), close;
	double			c = sqrt(h*h + h*k + k*k);
	double			rad = c/TWOPI;
	double			bottom = -height/2.0, top = height/2.0;
	Bcomponent*  	comp = NULL;
	Bcomponent*  	comp2 = NULL;
	Vector3<double>	lat, start, loc;
	
	if ( verbose & VERB_LABEL )
		cout << "Generating a tubular deltagraph" << endl;
		
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Lattice parameters:             " << h << " " << k << endl;
		cout << "Radius:                         " << rad << endl;
		cout << "Height:                         " << height << endl;
	}
	
	double			a = atan2(k*sqrt(3.0)/2, h+k/2.0);
	Vector3<double>	u(cos(a), -sin(a), 0), v(cos(M_PI/3-a), sin(M_PI/3-a), 0);
	Vector3<double>	w = v - u;
	
	while ( w[0] > 0 ) w -= u;
	
	start[1] = bottom;
	while ( lat[1] < top ) {
		for ( lat = start; lat[0] < c; lat += u ) {
			if ( lat[0] >= 0 && lat[0] < c && lat[1] >= bottom && lat[1] <= top ) {
				a = lat[0]/rad;
				loc[0] = rad*cos(a);
				loc[1] = rad*sin(a);
				loc[2] = lat[1];
				for ( close=0, comp2 = model->comp; comp2 && !close; comp2 = comp2->next )
					if ( comp2->location().distance(loc) < 0.5 ) close = 1;
				if ( !close ) {
//					id = Bstring(++nv, "%d");
//					comp = component_add(&comp, id);
//					if ( !model->comp ) model->comp = comp;
					if ( comp ) comp = comp->add(++nv);
					else model->comp = comp = new Bcomponent(++nv);
//					comp->type = model_add_type_by_id(model, hextype);
					comp->type(model->add_type(hextype));
					comp->location(loc);
					comp->radius(0.1);
				}
			}
		}
		start += w;
	}
	
	if ( verbose )
		cout << "Number of components generated: " << nv << endl << endl;

	model_link_list_generate(model, 1.2);
	
	return model;
}

/**
@brief 	Generates a new deltagraph based on the type.
@param 	type			type of cylinder.
@param 	radius			radius in integer units.
@param 	height			height in integer units.
@return Bmodel*			new model.

	Types supported:
		1	5-fold ends.
		2	6-fold ends.
		4	starting cap.
		8	ending cap.

**/
Bmodel*		model_delta_create_cylinder(int type, int radius, int height)
{
	int				order = 6;
	if ( type & 2 ) order = 5;
	
	Bstring			id("Cylinder");
	Bmodel*			model = new Bmodel(id);
	model->symmetry("D6");
	if ( type & 2 ) model->symmetry("D5");
	
	Bstring			pentype("PEN");
	Bstring			hextype("HEX");
	int				i, nv(0);
	double			a, sa(0), ea, f;
	double			da = TWOPI/(1.0L*order);
	double			hra = da/(2*radius);
	Bcomponent*  	comp = NULL;
	Vector3<double>	start, vec;	
	
	if ( verbose & VERB_LABEL )
		cout << "Generating a cylindrical deltagraph" << endl;
	
	if ( type & 4 ) {
//		id = "1";
//		comp = component_add(&model->comp, id);
		comp = model->add_component(++nv);
//		if ( order == 5 ) comp->type = model_add_type_by_id(model, pentype);
//		else comp->type = model_add_type_by_id(model, hextype);
		if ( order == 5 ) comp->type(model->add_type(pentype));
		else comp->type(model->add_type(hextype));
		comp->location()[2] = -height*sqrt(3.0)/4.0 - 0.5*radius;
		ea = TWOPI - da/2.0;
		for ( i=1; i<=radius; i++ ) {
			start[2] = -height*sqrt(3.0)/4.0 - 0.5*(radius - i);
			for ( a=0; a<ea; a+=da ) {
				start[0] = i*cos(a);
				start[1] = i*sin(a);
				vec[0] = i*cos(a+da) - start[0];
				vec[1] = i*sin(a+da) - start[1];
				for ( f=0; f<0.95; f+=1.0L/i ) {
//					id = Bstring(++nv, "%d");
//					comp = component_add(&comp, id);
					comp = comp->add(++nv);
//					if ( i == radius && f == 0 ) comp->type = model_add_type_by_id(model, pentype);
//					else comp->type = model_add_type_by_id(model, hextype);
					if ( i == radius && f == 0 ) comp->type(model->add_type(pentype));
					else comp->type(model->add_type(hextype));
					comp->location(start + vec*f);
				}
			}
		}
	}
	
	for ( i=1; i<height; i++ ) {
		start[2] = (i-height/2.0)*sqrt(3.0)/2.0;
		sa = 0;
		if ( i%2 ) sa = hra;
		ea = TWOPI - hra + sa;
		for ( a=sa; a<ea; a+=da ) {
			start[0] = radius*cos(a);
			start[1] = radius*sin(a);
			if ( height%2 == 0 || start[2] < 0 ) {
				vec[0] = radius*cos(a+da-2*sa) - start[0];
				vec[1] = radius*sin(a+da-2*sa) - start[1];
				if ( sa && radius > 1 ) vec *= radius/(radius - 1.0);
			} else {
				if ( sa == 0 ) {
					start[0] = radius*cos(a+2*hra);
					start[1] = radius*sin(a+2*hra);
					vec[0] = radius*cos(a+da) - start[0];
					vec[1] = radius*sin(a+da) - start[1];
					if ( radius > 1 ) vec *= radius/(radius - 1.0);
				} else {
					vec[0] = radius*cos(a+da) - start[0];
					vec[1] = radius*sin(a+da) - start[1];
				}
			}
			for ( f=0; f<0.95; f+=1.0L/radius ) {
//				id = Bstring(++nv, "%d");
//				comp = component_add(&model->comp, id);
				comp = comp->add(++nv);
//				comp->type = model_add_type_by_id(model, hextype);
				comp->type(model->add_type(hextype));
				comp->location(start + vec*f);
			}
		}
	}
	
	if ( type & 8 ) {
		if ( sa > 0 ) sa = 0;
		else sa = hra;
		for ( i=radius; i>0; i-- ) {
			start[2] = height*sqrt(3.0)/4.0 + 0.5*(radius - i);
			ea = TWOPI - hra + sa;
			for ( a=sa; a<ea; a+=da ) {
				start[0] = i*cos(a);
				start[1] = i*sin(a);
				vec[0] = i*cos(a+da) - start[0];
				vec[1] = i*sin(a+da) - start[1];
				for ( f=0; f<0.95; f+=1.0L/i ) {
//					id = Bstring(++nv, "%d");
//					comp = component_add(&comp, id);
					comp = comp->add(++nv);
//					if ( i == radius && f == 0 ) comp->type = model_add_type_by_id(model, pentype);
//					else comp->type = model_add_type_by_id(model, hextype);
					if ( i == radius && f == 0 ) comp->type(model->add_type(pentype));
					else comp->type(model->add_type(hextype));
					comp->location(start + vec*f);
					if ( verbose & VERB_FULL ) cout << comp->identifier() << ":\t" << comp->location() << endl;
				}
			}
		}

//		id = Bstring(++nv, "%d");
//		comp = component_add(&comp, id);	
		comp = comp->add(++nv);
//		if ( order == 5 )  comp->type = model_add_type_by_id(model, pentype);
//		else comp->type = model_add_type_by_id(model, hextype);
		if ( order == 5 ) comp->type(model->add_type(pentype));
		else comp->type(model->add_type(hextype));
		comp->location()[2] = height*sqrt(3.0)/4.0 + 0.5*radius;
	}
	
	model_link_list_generate(model, 1.3);
	
	return model;
}

