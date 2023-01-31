/**
@file	rwmodel_bild.cpp
@brief	Library routines to read and write Chimera BILD model parameters
@author 	Bernard Heymann
@date	Created: 20080521
@date	Modified: 20221115
**/


#include "rwmodel.h"
#include "model_neighbors.h"
#include "model_poly.h"
#include "Color.h"
#include "string_util.h"
#include "utilities.h"
#include "timer.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bcomponent*	component_add_if_not_exist(Vector3<double> loc, Bcomponent** comp, int &number)
{
	Bstring			id;
	Bcomponent*		comp1 = *comp;
	if ( comp1 ) comp1 = comp1->find_closest(loc);
		
//	for ( comp1 = *comp; comp1; comp1 = comp1->next )
//		if ( comp1->loc.distance(loc) < 0.001 ) break;
	
	if ( !comp1 ) {
//		id = Bstring(++number, "%d");
//		comp1 = component_add(comp, id);
		if ( *comp ) comp1 = (*comp)->add(++number);
		else *comp = comp1 = new Bcomponent(++number);
		comp1->location(loc);
	}
	
	return comp1;
}

/**
@brief 	Reads Chimera BILD model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
Bmodel*		read_model_bild(Bstring* file_list)
{
	int				i, j, nc, first_comment;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp1 = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	Bstring			id, line;
	Bstring*		filename;
	Bstring			*tokens, *token;
	ifstream		fmod;
	char			aline[1024];
	RGBA<float>		rgba(1,1,1,1);
	Vector3<double>	loc, vec;

	for ( i=1, filename = file_list; filename; filename = filename->next, i++ ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		first_comment = 1;
		fmod.open(filename->c_str());
		if ( fmod.fail() ) return NULL;
		id = Bstring(i, "%d");
//		mp = model_add(&mp, id);
//		if ( !model ) model = mp;
		if ( model ) mp = mp->add(i);
		else mp = new Bmodel(i);
		link = NULL;
		comp = NULL;
		nc = 0;
		while ( fmod.getline(aline, 1024) ) {
			line = aline;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG read_model_bild: " << line << endl;
			tokens = line.split();
			if ( tokens ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG read_model_bild: " << *tokens;
				if ( *tokens == ".comment" ) {
					if ( first_comment ) {
						for ( token = tokens->next; token; token = token->next )
							mp->comment() += token->str() + " ";
						first_comment = 0;
					}
				} else if ( *tokens == ".color" ) {
					if ( verbose & VERB_DEBUG )
						cout << " " << *(tokens->next);
					if ( isdigit((*tokens)[0]) )
						for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
							rgba[j] = token->real();
					else
						;
				} else if ( *tokens == ".transparency" ) {
					rgba[3] = tokens->next->real();
				} else if ( *tokens == ".dot" || *tokens == ".dotat" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->radius(1);
					comp->color(rgba);
				} else if ( *tokens == ".sphere" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					if ( token ) comp->radius(token->real());
					comp->color(rgba);
				} else if ( *tokens == ".m" || *tokens == ".move" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->color(rgba);
				} else if ( *tokens == ".mr" || *tokens == ".moverel" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						vec[j] = token->real();
					loc += vec;
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->color(rgba);
				} else if ( *tokens == ".d" || *tokens == ".draw" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->color(rgba);
					if ( comp && comp2 ) {
						link = link_add(&link, comp, comp2, 0, 1);
						if ( !mp->link ) mp->link = link;
						link->color(rgba);
					}
				} else if ( *tokens == ".dr" || *tokens == ".drawrel" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						vec[j] = token->real();
					loc += vec;
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->color(rgba);
					if ( comp && comp2 ) {
						link = link_add(&link, comp, comp2, 0, 1);
						if ( !mp->link ) mp->link = link;
						link->color(rgba);
					}
				} else if ( *tokens == ".v" || *tokens == ".vector" || *tokens == ".cylinder" ) {
					for ( j=0, token = tokens->next; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp = component_add_if_not_exist(loc, &mp->comp, nc);
					comp->color(rgba);
					for ( j=0; token && j<3; token = token->next, j++ )
						loc[j] = token->real();
					comp1 = component_add_if_not_exist(loc, &mp->comp, nc);
					comp1->color(rgba);
					if ( comp && comp1 ) {
						link = link_add(&link, comp, comp1, 0, 1);
						if ( !mp->link ) mp->link = link;
						link->color(rgba);
						if ( token ) link->radius(token->real());
						if ( link->radius() < 0.001 ) link->radius(1);
					}
				} else if ( *tokens == ".polygon" ) {
					comp1 = NULL;
					for ( j=0, token = tokens->next; token; token = token->next ) {
						loc[j++] = token->real();
						if ( j > 2 ) {
							j = 0;
							comp = component_add_if_not_exist(loc, &mp->comp, nc);
							comp->color(rgba);
							if ( comp && comp1 ) {
								link = link_add(&link, comp, comp1, 0, 1);
								if ( !mp->link ) mp->link = link;
								link->color(rgba);
							}
							comp1 = comp;
						}
					}
				}
				if ( verbose & VERB_DEBUG )
					cout << endl;
			}
			comp2 = comp;	// Previous component
		}
		fmod.close();
	}

	return model;
}

void		bild_color(ofstream& fbld, RGB<float> rgb)
{
	fbld << ".color " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
}

void		bild_color(ofstream& fbld, RGBA<float> rgba)
{
	fbld << ".color " << rgba[0] << " " << rgba[1] << " " << rgba[2] << endl;
}

void		bild_sphere(ofstream& fbld, Vector3<double> loc, double radius)
{
	fbld << ".sphere " << loc[0] << " " << loc[1] << " " << loc[2] << " " << radius << endl;
}

void		bild_cylinder(ofstream& fbld, Vector3<double> start, Vector3<double> end, double radius)
{
	fbld << ".cylinder " << start[0] << " " << start[1] << " " << start[2]
		<< " " << end[0] << " " << end[1] << " " << end[2]
		<< " " << radius << endl;
}

void		bild_polygon_start(ofstream& fbld, Vector3<double> vertex)
{
	fbld << ".polygon " << vertex[0] << " " << vertex[1] << " " << vertex[2];
}
void		bild_polygon_add(ofstream& fbld, Vector3<double> vertex)
{
	fbld << " " << vertex[0] << " " << vertex[1] << " " << vertex[2];
}
void		bild_polygon_end(ofstream& fbld, Vector3<double> vertex)
{
	fbld << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
}

void		bild_arrow(ofstream& fbld, Vector3<double> start, Vector3<double> end, double cylrad, double conerad, double rho)
{
	fbld << ".arrow " << start[0] << " " << start[1] << " " << start[2]
		<< " " << end[0] << " " << end[1] << " " << end[2]
		<< " " << cylrad << " " << conerad << " " << rho << endl;
}

/**
@brief 	Writes Chimera BILD model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_bild(Bstring& filename, Bmodel* model, int splt)
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	Bstring			onename;
	char			format[32];

	snprintf(format, 32, "_%%0%dd.", splt);

	ofstream		fmod;

	for ( n=1, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = filename.pre_rev('.') + Bstring(n, format) + filename.post_rev('.');
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( mp->comment().length() ) {
//			Bstring*	sa = mp->comment().split("\n");
//			for ( Bstring* s = sa; s; s = s->next )
//				if ( s->length() ) fmod << ".comment " << *s << endl;
			vector<string>	sa = split(mp->comment(), '\n');
			for ( auto s: sa )
				if ( s.size() ) fmod << ".comment " << s << endl;
		}
		for ( comp = mp->comp; comp; comp = comp->next ) {
			bild_color(fmod, comp->color());
			bild_sphere(fmod, comp->location(), comp->radius());
		}
		for ( link = mp->link; link; link = link->next ) {
			bild_color(fmod, link->color());
			bild_cylinder(fmod, link->comp[0]->location(), link->comp[1]->location(), link->radius());
		}
		fmod.close();
	}
	
	return n;
}

/**
@brief 	Converts a model into a representation with views.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	vec_type	0=views, 1=vectors, 2=velocity vectors.
@param 	color_type	flag to color components.
@return int			0, <0 on error.

	A sphere is drawn for every component with a cylinder indicating its view.

**/
int			model_to_bild_orientations(Bstring& filename, Bmodel* model, int vec_type, int color_type)
{
	int				n, selmax;
	double			vlen, vrad, dmin(1e10), dmax(-1e10), fommin(1e10), fommax(-1e10);
	RGB<float>		rgb(1,1,1);
	Vector3<double>	vv, av, v1, v2, v3;
	Matrix3			mat;
	Bcomponent*		comp;
	
	for ( n=selmax=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		if ( dmin > comp->density() ) dmin = comp->density();
		if ( dmax < comp->density() ) dmax = comp->density();
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
		if ( selmax < comp->select()) selmax = comp->select();
	}

	if ( verbose ) {
		cout << "Generating a bild views file:   " << filename << endl;
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select()) {
		vlen = comp->radius()*2;
		vrad = comp->radius()/5;
		switch ( vec_type ) {
			case 1: vv = comp->view().vector3(); break;
			case 2: vv = comp->force(); break;
			case 3: vv = comp->velocity(); break;
		}
		vv *= vlen;
		if ( vec_type < 2 ) {
			av = Vector3<double>(1,0,0);
			mat = comp->view().matrix();
//			mat = mat.transpose();
			av = mat * av;
			av *= 3*vrad;
		}
		v1 = comp->location()+ vv;
		v2 = comp->location()+ vv*0.75;
		v3 = v2 + av;
		switch ( color_type ) {
			case 1:
				rgb.spectrum(comp->select(), 0, selmax);
				break;
			case 2:
				rgb.spectrum(comp->density(), dmin, dmax);
				break;
			case 3:
				rgb.spectrum(comp->FOM(), fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
		bild_sphere(fbld, comp->location(), comp->radius());
		bild_arrow(fbld, comp->location(), v1, vrad, vrad, 0.75);
		if ( vec_type < 2 ) {
			bild_polygon_start(fbld, v1);
			bild_polygon_add(fbld, v2);
			bild_polygon_end(fbld, v3);
		}
		fbld << endl;		
	}
	
	fbld << endl;
	
	fbld.close();
		
	return 0;
}

/**
@brief 	Converts a model into a representation with views.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	color_type	flag to color components.
@return int			0, <0 on error.

	A sphere is drawn for every component with a cylinder indicating its view.

**/
int			model_to_bild_view_sphere(Bstring& filename, Bmodel* model, int color_type)
{
	int				n, selmax;
	double			vrad(TWOPI/360), dmin(1e10), dmax(-1e10), fommin(1e10), fommax(-1e10), fomrange(1);
	RGB<float>		rgb(1,1,1);
	Vector3<double>	vv, av, v1, v2, v3;
	Matrix3			mat;
	Bcomponent*		comp;
	
	for ( n=selmax=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		if ( dmin > comp->density() ) dmin = comp->density();
		if ( dmax < comp->density() ) dmax = comp->density();
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
		if ( selmax < comp->select()) selmax = comp->select();
	}
	
	fomrange = 5*(fommax - fommin);

	if ( verbose ) {
		cout << "Generating a bild views file:   " << filename << endl;
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	fbld << ".transparency 0.5" << endl;
	bild_sphere(fbld, Vector3<double>(0,0,0), 1);
	fbld << ".transparency 0.0" << endl;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select()) {
		vv = comp->view().vector3();
		av = Vector3<double>(1,0,0);
		mat = comp->view().matrix();
//		mat = mat.transpose();
		av = mat * av;
		av *= 3*vrad;
		if ( fomrange > 0 ) v1 = vv*(1 + vrad + (comp->FOM()- fommin)/fomrange);
		else v1 = vv*(1 + 10*vrad);
		v2 = vv*1.15;
		v3 = v2 + av;
		switch ( color_type ) {
			case 1:
				rgb.spectrum(comp->select(), 0, selmax);
				break;
			case 2:
				rgb.spectrum(comp->density(), dmin, dmax);
				break;
			case 3:
				rgb.spectrum(comp->FOM(), fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
//		bild_arrow(fbld, vv, v1, vrad, vrad, 0.75);
		bild_cylinder(fbld, vv, v1, vrad);
//		bild_polygon_start(fbld, v1);
//		bild_polygon_add(fbld, v2);
//		bild_polygon_end(fbld, v3);
		fbld << endl;		
	}
	
	fbld << endl;
	
	fbld.close();
		
	return 0;
}

/**
@brief 	Converts a model into a representation with views.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	color_type	flag to color components.
@return int			0, <0 on error.

	A sphere is drawn for every component with a cylinder indicating its view.

**/
int			model_to_bild_force_vectors(Bstring& filename, Bmodel* model, int color_type)
{
	int				n, selmax;
	double			vrad(TWOPI/360), dmin(1e10), dmax(-1e10);
	double			f, forcemax(0), fommin(1e10), fommax(-1e10), scale(1);
	RGB<float>		rgb(1,1,1);
	Vector3<float>	fv, v0;
	Bcomponent*		comp;
	
	for ( n=selmax=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		if ( dmin > comp->density() ) dmin = comp->density();
		if ( dmax < comp->density() ) dmax = comp->density();
		f = comp->force().length();
		if ( forcemax < f ) forcemax = f;
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
		if ( selmax < comp->select()) selmax = comp->select();
	}
	
	scale = 1/forcemax;
	vrad = scale/360;
//	fomrange = 5*(fommax - fommin);

	if ( verbose ) {
		cout << "Generating a bild force vector file:   " << filename << endl;
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	fbld << ".transparency 0.5" << endl;
	bild_sphere(fbld, Vector3<double>(0,0,0), 1);
	fbld << ".transparency 0.0" << endl;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select()) {
		fv = v0 = comp->force();
		fv *= scale;
		v0.normalize();
		fv += v0;
		switch ( color_type ) {
			case 1:
				rgb.spectrum(comp->select(), 0, selmax);
				break;
			case 2:
				rgb.spectrum(comp->density(), dmin, dmax);
				break;
			case 3:
				rgb.spectrum(comp->FOM(), fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
		bild_cylinder(fbld, v0, fv, vrad);
		fbld << endl;
	}
	
	fbld << endl;
	
	fbld.close();
		
	return 0;
}

/**
@brief 	Converts a model into a representation with view polygons.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	order		polygon order.
@param 	color_type	flag to color components.
@return int			0, <0 on error.

	A sphere is drawn for every component with a cylinder indicating its view.

**/
int			model_to_bild_view_polygons(Bstring& filename, Bmodel* model, int order, int color_type)
{
	int					i, n, selmax;
	double				dmin(1e10), dmax(-1e10), fommin(1e10), fommax(-1e10);
	RGB<float>			rgb(1,1,1);
	Vector3<double>		av;
	Matrix3				mat;
	Bcomponent*			comp;
	
	double				a(TWOPI*1.0L/order);
	Vector3<double>*	p = new Vector3<double>[order];
	
	for ( i=0; i<order; i++ )
		p[i] = Vector3<double>(cos(i*a), sin(i*a), 0);
	
	for ( n=selmax=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		if ( dmin > comp->density() ) dmin = comp->density();
		if ( dmax < comp->density() ) dmax = comp->density();
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
		if ( selmax < comp->select()) selmax = comp->select();
	}

	if ( verbose ) {
		cout << "Generating a bild view polygons file:   " << filename << endl;
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select()) {
		switch ( color_type ) {
			case 1:
				rgb.spectrum(comp->select(), 0, selmax);
				break;
			case 2:
				rgb.spectrum(comp->density(), dmin, dmax);
				break;
			case 3:
				rgb.spectrum(comp->FOM(), fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
		mat = comp->view().matrix();
		for ( i=0; i<order; i++ ) {
			av = comp->location()+ mat * p[i] * comp->radius();
			if ( i == 0 ) bild_polygon_start(fbld, av);
			else if ( i < order - 1 ) bild_polygon_add(fbld, av);
			else bild_polygon_end(fbld, av);
		}
		fbld << endl;
	}
	
	fbld << endl;
	
	fbld.close();
	
	delete[] p;
		
	return 0;
}

/**
@brief 	Converts a model into a bild polygon representation.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	color_type	flag to color polygons.
@return int			0, <0 on error.

	Color scheme types:
		1	spectrum/rainbow.
		2	fom.

**/
int			model_to_bild_polygons(Bstring& filename, Bmodel* model, int color_type)
{
	if ( !model->poly ) model_poly_generate(model);
	
	int				i, n, npoly;
	double			fom, fommin(1e10), fommax(-1e10);
	RGB<float>		rgb(1,1,1);
	Bcomponent*		comp;
	Bpolygon*		poly;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
	}

	for ( npoly=0, poly = model->poly; poly; poly = poly->next ) npoly++;
	
	if ( verbose ) {
		cout << "Generating a bild polygon file: " << filename << endl;
		switch ( color_type ) {
			case 1: cout << "Color spectrum range:           0 - " << npoly << endl; break;
			case 2: cout << "FOM range:                      " << fommin << " - " << fommax << endl; break;
			default: cout << "Homogeneous color" << endl;
		}
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	for ( n=0, poly = model->poly; poly; poly = poly->next, n++ ) {
		for ( i=0, fom=0; i<poly->comp.size() && poly->comp[i]; i++ ) fom += poly->comp[i]->FOM();
		if ( i ) fom /= i;
		switch ( color_type ) {
			case 1:
				rgb.spectrum(n, 0, npoly);
				break;
			case 2:
				rgb.spectrum(fom, fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
		fbld << ".polygon";
		for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ )
			bild_polygon_add(fbld, poly->comp[i]->location());
		fbld << endl;
	}
	fbld << endl;
	
	fbld.close();
		
	return 0;
}

/**
@brief 	Converts a model into a bild neighbor plane representation.
@param 	&filename	bild format file name.
@param 	*model		model structure.
@param 	color_type	flag to color planes.
@return int			0, <0 on error.

	Color scheme types:
		1	spectrum/rainbow.
		2	fom.

**/
int			model_to_bild_neighbor_planes(Bstring& filename, Bmodel* model, int color_type)
{
	model_set_neighbors(model, 8, 50);
	
	int				i, n, ncomp;
	double			fom, fommin(1e10), fommax(-1e10);
	RGB<float>		rgb;
	Bcomponent*		comp;
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		if ( fommin > comp->FOM()) fommin = comp->FOM();
		if ( fommax < comp->FOM()) fommax = comp->FOM();
	}

	if ( verbose ) {
		cout << "Generating a bild neighbor plane file: " << filename << endl;
		switch ( color_type ) {
			case 1: cout << "Color spectrum range:           0 - " << ncomp << endl; break;
			case 2: cout << "FOM range:                      " << fommin << " - " << fommax << endl; break;
			default: cout << "Homogeneous color" << endl;
		}
	}
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) {
		cerr << "Error: Unable to open " << filename << " for writing" << endl;
		return -1;
	}
	
	fbld << endl;
	fbld << ".color white" << endl;
	for ( n=0, comp = model->comp; comp; comp = comp->next, n++ ) {
		for ( i=0, fom=0; i<comp->neighbor.size() && comp->neighbor[i]; i++ )
			fom += comp->neighbor[i]->FOM();
		if ( i ) fom /= i;
		switch ( color_type ) {
			case 1:
				rgb.spectrum(n, 0, ncomp);
				break;
			case 2:
				rgb.spectrum(fom, fommin, fommax);
				break;
		}
		bild_color(fbld, rgb);
		fbld << ".polygon";
		for ( i=0; i<comp->neighbor.size() && comp->neighbor[i]; i++ )
			bild_polygon_add(fbld, comp->neighbor[i]->location());
		fbld << endl;
	}
	fbld << endl;
	
	fbld.close();
		
	return 0;
}

/*
BILD format

BILD is a simple text format that describes lines, polygons, and geometric primitives such as spheres, boxes, cylinders, and cones with commands (see example.bild and XYZ-axes.bild). Character strings can also be included. The commands are used to generate a VRML model in Chimera. The objects can be specified in absolute coordinates and/or transformed and scaled coordinates. See also: shape, define, geometric objects

This file type is indicated by the filename suffix .bild or .bld, or by using bild:filename in the command line.

Curved objects (spheres, cones, and cylinders) read from BILD format are described as perfectly smooth in exported files, but display within Chimera uses planar facets. The number of facets and the apparent smoothness can be increased by raising the subdivision level with the Effects tool or the set command.

BILD Commands

Square brackets denote optional parameters.

.arrow x1 y1 z1 x2 y2 z2 [r1 [r2 [rho]]] 
Draw an arrow from (x1, y1, z1) to (x2, y2, z2). An arrow consists of a cylinder stretching from the initial point to an intermediary junction, and a cone stretching from the junction to the final point. The radius of the cylinder is r1 (default 0.1), the radius of the base of the cone is r2 (default 4*r1), and the fraction of the total distance taken up by the cylinder is rho (default 0.75).
.box x1 y1 z1 x2 y2 z2 
Draw a box with opposite corners at (x1, y1, z1) and (x2, y2, z2).
.color name 
or
.color r g b 
Set the color of subsequently defined items. The name can be a built-in name, a name defined previously with colordef, or an integer that refers to the old BILD color wheel (0-65, inclusive). Alternatively, a color can be described by its red (r), green (g), and blue (b) components, each in the range 0-1, inclusive. Any transparency in the color is ignored, but transparency can be set separately.
.cmov x y z 
Define the starting location of the next character string. Lines in the BILD file that do not start with a period (.) are character strings to be displayed. See also .font.
.comment text 
User comment line (ignored during object creation).
.cone x1 y1 z1 x2 y2 z2 r [open] 
Draw a cone with a base of radius r centered at (x1, y1, z1) and a tip at (x2, y2, z2). If the keyword open is present, the base of the cone will be invisible.
.cylinder x1 y1 z1 x2 y2 z2 r [open] 
Draw a cylinder with radius r and bases centered at (x1, y1, z1) and (x2, y2, z2). If the keyword open is present, the bases of the cylinder will be invisible.
.dotat x y z 
or
.dot x y z 
Draw a sphere of unit radius centered at (x, y, z). The sphere center is treated as a vertex if there is a subsequent .draw, .drawrel, or .moverel command.
.draw x y z 
or
.d x y z 
Add another vertex to the current set of line segments. A line will be drawn from the previous vertex to this vertex at (x, y, z). There should be a prior .move, .moverel, .dotat, or .marker command (these initiate sets of line segments).
.drawrel dx dy dz 
or
.dr dx dy dz 
Add another vertex to the current set of line segments. A line will be drawn from the previous vertex at (x, y, z) to this vertex at (x + dx, y + dy, z + dz).
.font fontname pointsize [fontstyle] 
Set the font, font size, and font style of subsequent character strings. Lines in the BILD file that do not start with a period (.) are character strings to be displayed. Options for fontname include: Times, Helvetica, Courier, SERIF, SANS, TYPEWRITER. Options for fontstyle: plain, bold, italic, bold italic. See also .cmov.
.marker x y z 
Draw a box of unit cubic diagonal centered at (x, y, z), i.e., a box with opposite corners at (x – 0.5, y – 0.5, z – 0.5) and (x + 0.5, y + 0.5, z + 0.5). The box center is treated as a vertex if there is a subsequent .draw, .drawrel, or .moverel command.
.move x y z 
or
.m x y z 
Start a new set of line segments whose first vertex is at (x, y, z).
.moverel dx dy dz 
or
.mr dx dy dz 
Start a new set of line segments whose first vertex is at (x + dz, y + dy, z + dz), where (x, y, z) is the coordinate of the last vertex defined.
.polygon x1 y1 z1 x2 y2 z2 ... xN yN zN 
Draw a flat polygon with vertices at (x1, y1, z1), (x2, y2, z2), ..., (xN, yN, zN).
.pop 
Discard the most recent transformation (rotation, scaling, or translation) from the transformation stack.
.rotate angle axis 
or
.rot angle axis 
Rotate all subsequent coordinates by angle degrees about the given axis. The axis can be given as a single letter (x, y, or z) or as three numbers defining an arbitrary vector: xa ya za. This transformation is added to the top of the transformation stack.
.scale xscale [yscale [zscale]] 
Scale all subsequent coordinates by the given factor(s). The x coordinates will be scaled by xscale, y coordinates by yscale (equal to xscale by default), and z coordinates by zscale (equal to xscale by default). This transformation is added to the top of the transformation stack.
.sphere x y z r 
Draw a sphere centered at (x, y, z) with radius r.
.translate dx dy dz 
or
.tran dx dy dz 
Translate all subsequent coordinates by the specified amount. This transformation is added to the top of the transformation stack.
.transparency value 
Set the transparency of subsequently defined items. The value can range from 0.0 (not transparent) to 1.0 (completely transparent).
.vector x1 y1 z1 x2 y2 z2 
or
.v x1 y1 z1 x2 y2 z2 
Draw a line segment from (x1, y1, z1) to (x2, y2, z2). This command is a shorthand for
.m x1 y1 z1 
.d x2 y2 z2
UCSF Computer Graphics Laboratory / October 2013
*/
