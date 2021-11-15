/**
@file	model_color.cpp
@brief	Functions to color models.
@author Bernard Heymann
@date	Created: 20080206
@date	Modified: 20210319
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_shell.h"
#include "model_views.h"
#include "model_util.h"
#include "Color.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Colors model components uniformly.
@param 	*model		model to color.
@param 	color		4-value color vector ([0,1]).
@return int			0.

	Both selected and non-selected elements are colored.

**/
int			model_color_uniformly(Bmodel* model, RGBA<float> color)
{
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components uniformly: " << color[0] << " " << 
			color[1] << " " << color[2] << " " << color[3] << endl << endl;
	
	for ( m = model; m; m = m->next ) {
		for ( comp = m->comp; comp; comp = comp->next )
 			comp->color(color);
		for ( link = m->link; link; link = link->next )
 			link->color(color);
	}

	return 0;
}

/**
@brief 	Colors selected model components.
@param 	*model		model to color.
@param 	color		4-value color vector ([0,1]).
@return int			0.
**/
int			model_color_selected(Bmodel* model, RGBA<float> color)
{
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components by selection: " << color[0] << " " << 
			color[1] << " " << color[2] << " " << color[3] << endl << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() )
			comp->color(color);
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			if ( link->comp[0]->select() && link->comp[1]->select() )
				link->adopt_component_color();
		}
	}
	
	return 0;
}

/**
@brief 	Colors model components in order from red to blue.
@param 	*model		model to color.
@return int			0.
**/
int			model_color_by_order(Bmodel* model)
{
	int				i, ncomp;
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components by order" << endl << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( ncomp = 0, comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) ncomp++;
		for ( i=0, comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->color().spectrum(i, 0, ncomp);
			comp->color()[3] = 1;
			i++;
		}
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->adopt_component_color();
		}
	}
	
	return 0;
}

/**
@brief 	Colors model components by density.
@param 	*model		model to color.
@return int			0.

	Color assignments are:
		density = 0		blue
		density = 0.5	green
		density = 1		red

**/
int			model_color_by_density(Bmodel* model)
{
	double			dmin = 1e30, dmax = -1e30;
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( dmin > comp->density() ) dmin = comp->density();
			if ( dmax < comp->density() ) dmax = comp->density();
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components by density (" << dmin << " - " << dmax << ")" << endl << endl;

	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->color().spectrum(comp->density(), dmin, dmax);
			comp->color()[3] = 1;
		}
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->adopt_component_color();
		}
	}
	
	return 0;
}

/**
@brief 	Colors model components by fom.
@param 	*model		model to color.
@return int			0.

	Color assignments are:
		FOM = 0		blue
		FOM = 0.5	green
		FOM = 1		red

**/
int			model_color_by_fom(Bmodel* model)
{
	double			fmin = 1e30, fmax = -1e30;
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( fmin > comp->FOM() ) fmin = comp->FOM();
			if ( fmax < comp->FOM() ) fmax = comp->FOM();
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components by FOM (" << fmin << " - " << fmax << ")" << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->color().spectrum(comp->FOM(), fmin, fmax);
			comp->color()[3] = 1;
		}
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->adopt_component_color();
		}
	}
	
	return 0;
}

/**
@brief 	Colors model components by the selection number.
@param 	*model		model to color.
@return int			0.

	Color assignments are:
		FOM = 0		blue
		FOM = 0.5	green
		FOM = 1		red

**/
int			model_color_by_selection(Bmodel* model)
{
	long			smin(1000000), smax(0);
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( smin > comp->select() ) smin = comp->select();
			if ( smax < comp->select() ) smax = comp->select();
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring components by selection (" << smin << " - " << smax << ")" << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->color().spectrum(comp->select(), smin, smax);
			comp->color()[3] = 1;
		}
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->adopt_component_color();
		}
	}
	
	return 0;
}


/**
@brief 	Colors model components by fom.
@param 	*model		model to color.
@param 	rgba		RGBA color.
@return int					0.

	Color assignments are:
		FOM = 0		blue
		FOM = 0.5	green
		FOM = 1		red

**/
int			model_color_selected_types(Bmodel* model, RGBA<float> rgba)
{
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring selected component types" << endl << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			if ( comp->type()->select() )
				comp->color(rgba);
		}
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->adopt_component_color();
		}
	}
	
	return 0;
}

/**
@brief 	Colors model vertices and links as a function of curvature.
@param 	*model		model to color.
@return int			0.
**/
int			model_color_curvature(Bmodel* model)
{
	Bmodel*			m;
	Bcomponent*		comp;
	Blink*			link;

	model_curvature(model);
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring based on curvature" << endl << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( link = m->link; link; link = link->next ) if ( link->select() ) {
			link->color().spectrum(link->FOM(), 0, M_PI_4);
		}
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->color().spectrum(comp->FOM(), 0, M_PI_4);
		}
	}	
	
	return 0;
}

/**
@brief 	Colors model vertices as a function of chirality.
@param 	*model		model to color.
@return int			number of chiral vertices.

	Chiral vertices are colored blue for + and red for - handedness.

**/
int			model_color_chiral_vertices(Bmodel* model)
{
	int				h, n(0);
	RGBA<float>		red(1.0,0.0,0.0,1.0);
	RGBA<float>		blu(0.0,0.0,1.0,1.0);
	Bmodel*			m;
	Bcomponent*		comp;

	model_extended_vertex_types(model);
	
	if ( verbose & VERB_PROCESS )
		cout << "Coloring based on chirality." << endl;
	
	for ( m = model; m; m = m->next ) if ( m->select() ) {
		for ( comp = m->comp; comp; comp = comp->next ) if ( comp->select() ) {
			h = component_hand(comp->type()->identifier());
			if ( h > 0 ) comp->color(blu);
			else if ( h < 0 ) comp->color(red);
			if ( h ) n++;
		}
	}	
	
	if ( verbose )
		cout << "Number of chiral vertices:      " << n << endl << endl;
	
	return n;
}
