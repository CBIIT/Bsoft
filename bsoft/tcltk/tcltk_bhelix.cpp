/**
@file	tcltk_bhelix.cpp
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20130726
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bhelix.h"
#include "mg_helix.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern Bimage* 	imglist;

// Internal function prototypes
Tcl_Obj*	layerline_count(Bmicrograph* mg);
Tcl_Obj*	layerline_ids(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_renumber(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_distance(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_fom(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_select(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_move(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_create(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_order(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_generate(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_mask(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_plot(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	layerline_bessel_function(double realsizeX, int objc, Tcl_Obj *CONST objv[]);


/*
@brief 	Functions for helical indexing.
	Layer line interface syntax:
		Bmg layerline [item] [action] [params]
	Actions:
		count						counts the number of layer lines
		ids							returns the identifiers of the layer lines
		location [id]				returns the location of a layer line
		fom [id]					returns the figure-of-merit of a layer line
		select [x] [y]				returns the layer line identifier at the given location
		move [id] [dx] [dy]			updates a layer line location by the given vector
		create [x] [y] [id]			creates a new layer line at the given location
		delete [all|mg|id]			deletes layer lines
		order [id] [dx] [dy]		sets or returns the order of a layer line
		generate [m]				generates all layer lines
		mask [w]					masks layer lines
		plot [id] [l]				retrieves one layer line's values as an array of length l
		bessel [o] [r] [l]			retrieves a bessel function of given order and radius as an array of length l
	Parameters:
		id			layer line identifier
		x			x location
		y			y location
		dx			change in x
		dy			change in y
		m			index limit
		w			mask width
		l			array length
		o			Bessel order
*/
Tcl_Obj*	do_layerline(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_layerline: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	if ( !item.contains("Field") && !item.contains("Micrograph") ) {
		cerr << "Error in do_layerline: Item " << item << " must be a micrograph!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2, base, delete_what;
	char				string[MAXLINELEN] = "";
	
	if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
		if ( field ) mg = field->mg;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next ) ;
			if ( mg ) break;
		}
	}

	if ( mg ) {
		filename = mg->fft;
		if ( filename.length() < 1 ) filename = mg->fps;
		if ( filename.length() < 1 ) filename = mg->fmg;
		if ( filename.length() < 1 ) filename = mg->fframe;
	}
	
	if ( filename.length() ) {
		base = filename.base();
		for ( p = imglist; p; p = p->next ) {
			filename2 = p->file_name();
			if ( filename2.base() == base ) break;
		}
	}
			
	if ( !p ) {
		returnObj = Tcl_NewObj();
		snprintf(string, MAXLINELEN, "Image %s not found!", filename.c_str());
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_layerline: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		returnObj = layerline_count(mg);
	} else if ( action == "ids" ) {
		returnObj = layerline_ids(mg, objc, objv);
	} else if ( action == "renumber" ) {
		returnObj = layerline_renumber(mg, objc, objv);
	} else if ( action == "distance" ) {
		returnObj = layerline_distance(mg, objc, objv);
	} else if ( action == "fom" ) {
		returnObj = layerline_fom(mg, objc, objv);
	} else if ( action == "select" ) {
		if ( p ) returnObj = layerline_select(mg, p, objc, objv);
	} else if ( action == "move" ) {
		returnObj = layerline_move(mg, objc, objv);
	} else if ( action == "create" ) {
		if ( p ) returnObj = layerline_create(mg, p, objc, objv);
	} else if ( action == "delete" ) {
		if ( mg ) returnObj = layerline_delete(mg, objc, objv);
		else returnObj = layerline_delete(project, objc, objv);
	} else if ( action == "order" ) {
		returnObj = layerline_order(mg, objc, objv);
	} else if ( action == "generate" ) {
		returnObj = layerline_generate(mg, objc, objv);
	} else if ( action == "mask" ) {
		if ( p ) returnObj = layerline_mask(mg, p, objc, objv);
	} else if ( action == "plot" ) {
		if ( p ) returnObj = layerline_plot(mg, p, objc, objv);
	} else if ( action == "bessel" ) {
		if ( p ) returnObj = layerline_bessel_function(p->sampling(0)[0]*p->sizeX(), objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	if ( !returnObj ) {
		returnObj = Tcl_NewObj();
		Tcl_SetIntObj(returnObj, 0);
	}

	return returnObj;
}

Tcl_Obj*	layerline_count(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nll(0);
	Blayerline*			line;

	if ( mg )
		for ( line = mg->layer; line; line = line->next ) nll++;
	
	Tcl_SetIntObj(returnObj, nll);
	
	return returnObj;
}

Tcl_Obj*	layerline_ids(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				fom_cut(-1);
	char				string[MAXLINELEN] = "";
	Blayerline*			line;

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	if ( mg ) {
		for ( line = mg->layer; line; line = line->next ) if ( line->fom >= fom_cut ) {
			snprintf(string, MAXLINELEN, " %d", line->number);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}
		
	return returnObj;
}

Tcl_Obj*	layerline_renumber(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id1, id2;
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id1);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id2);

	Blayerline*			line = NULL;
	
	if ( mg ) {
		for ( line = mg->layer; line && line->number != id2; line = line->next ) ;
		if ( line ) {
			id2 = id1;
		} else {
			for ( line = mg->layer; line && line->number != id1; line = line->next ) ;
			if ( line ) line->number = id2;
		}
	}

	Tcl_SetIntObj(returnObj, id2);
	
	return returnObj;
}

Tcl_Obj*	layerline_distance(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id;
	double				d(0);
	Blayerline*			line = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);

	if ( mg )
		for ( line = mg->layer; line && line->number != id; line = line->next ) ;

	if ( line ) d = line->distance;
	
	Tcl_SetDoubleObj(returnObj, d);
	
	return returnObj;
}

Tcl_Obj*	layerline_fom(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id;	
	double				fom(1);
	Blayerline*			line = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);

	if ( mg )
		for ( line = mg->layer; line && line->number != id; line = line->next ) ;

	if ( line ) fom = line->fom;
	
	Tcl_SetDoubleObj(returnObj, fom);
	
	return returnObj;
}

Tcl_Obj*	layerline_select(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					i(0), id(-1);
	unsigned long		n = mg->img_num;
	double				x(0), y(0), w(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &w);
	
	double				d, ds, dmin;
//	Vector3<float>		loc(x, y, 0);
	Blayerline*			line;

	if ( mg ) {
		ds = cos(mg->helix_axis)*(x - p->image[n].origin()[0]) + sin(mg->helix_axis)*(y - p->image[n].origin()[1]);
		dmin = 2*w;
		for ( line = mg->layer; line; line = line->next ) {
			d = fabs(ds - line->distance);
			if ( dmin > d ) {
				dmin = d;
				i = line->number;
			}
		}
		if ( dmin <= w ) id = i;
	}

	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	layerline_move(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(-1);
	double				d(0), dx(0), dy(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &dx);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &dy);

	Blayerline*			line;

	if ( mg ) {
		for ( line = mg->layer; line && line->number != id; line = line->next ) ;
		if ( line ) {
			d = cos(mg->helix_axis)*dx + sin(mg->helix_axis)*dy;
			line->distance += d;
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	layerline_create(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);
	unsigned long		n = mg->img_num;
	double				x, y;
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &x);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &y);

	Blayerline*			line = NULL;

	if ( mg )
		line = (Blayerline *) add_item((char **) &mg->layer, sizeof(Blayerline));

	if ( line ) {
		line->number = id;
		line->distance = fabs(cos(mg->helix_axis)*(x - p->image[n].origin()[0]) + sin(mg->helix_axis)*(y - p->image[n].origin()[1]));
		line->freq = fabs(line->distance/(p->sizeY()*mg->pixel_size[1]));
		line->amp = 1;
		line->fom = 1;
		line->sel = 1;
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	layerline_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), n(-1);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Blayerline*			line;

	if ( action == "all" || action == "mg" ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				kill_list((char *) mg->layer, sizeof(Blayerline));
				mg->layer = NULL;
			}
		}
	}

	if ( action == "fom" ) {
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);
		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &n);
			if ( n < 0 ) {
				if ( field ) for ( mg = field->mg; mg; mg = mg->next ) {
					for ( line = mg->layer; line; line = line->next )
						if ( line->fom < fom_cut ) line->sel = 0;
//					layerline_delete_non_selected(&mg->layer);
				}
			} else {
				if ( field ) for ( mg = field->mg; mg && mg->img_num != n; mg = mg->next ) ;
				if ( mg ) {
					for ( line = mg->layer; line; line = line->next )
						if ( line->fom < fom_cut ) line->sel = 0;
//					layerline_delete_non_selected(&mg->layer);
				}
			}
	} else {
		if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &id);
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &n);
		if ( id > 0 ) {
				if ( n < 0 ) {
					if ( field ) for ( mg = field->mg; mg; mg = mg->next ) {
						for ( line = mg->layer; line && line->number != id; line = line->next ) ;
						if ( line )
							remove_item((char **)&mg->layer, (char *)line, sizeof(Blayerline));
					}
				} else {
					if ( field ) for ( mg = field->mg; mg && mg->img_num != n; mg = mg->next ) ;
					if ( mg ) {
						for ( line = mg->layer; line && line->number != id; line = line->next ) ;
						if ( line )
							remove_item((char **)&mg->layer, (char *)line, sizeof(Blayerline));
					}
				}
		}
	}	
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	layerline_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Blayerline*			line;

	if ( action == "all" ) {
		kill_list((char *) mg->layer, sizeof(Blayerline));
		mg->layer = NULL;
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
		if ( id > 0 ) {
			if ( mg ) {
				for ( line = mg->layer; line && line->number != id; line = line->next ) ;
				if ( line )
					remove_item((char **)&mg->layer, (char *)line, sizeof(Blayerline));
			}
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	layerline_order(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(-1), order(0), set_order(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) {
		Tcl_GetIntFromObj(NULL, objv[5], &order);
		set_order = 1;
	}

	Blayerline*			line;

	if ( mg ) {
		for ( line = mg->layer; line && line->number != id; line = line->next ) ;
		if ( line ) {
			if ( set_order ) line->order = order;
			else order = line->order;
		}
	}
	
	Tcl_SetIntObj(returnObj, order);
	
	return returnObj;
}

Tcl_Obj*	layerline_generate(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					n(0), ind_lim(3);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &ind_lim);
	
	if ( mg )
		n = mg_generate_layer_lines(mg, ind_lim);

	Tcl_SetIntObj(returnObj, n);
	
	return returnObj;
}

Tcl_Obj*	layerline_mask(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					n(0);
	double				width(5);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &width);
	
	if ( mg )
		n = img_mask_layer_lines(p, mg->layer, mg->helix_axis, width);

	Tcl_SetIntObj(returnObj, n);
	
	return returnObj;
}

Tcl_Obj*	layerline_plot(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					l(0), len(p->sizeX()), i;
	Blayerline*			line;
	double*				plot;
	char				string[64];
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &l);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &len);
	
	if ( mg ) {
		for ( line = mg->layer; line && line->number != l; line = line->next ) ;
		if ( line ) {
			plot = img_extract_layer_line(p, line, mg->helix_axis, len);
			snprintf(string, 60, "%d", len);
			Tcl_AppendToObj(returnObj, string, strlen(string));
			for ( i=0; i<p->sizeX(); i++ ) {
				snprintf(string, 60, " %g", plot[i]);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
			delete[] plot;
		}
	}
	
	return returnObj;
}

Tcl_Obj*	layerline_bessel_function(double realsizeX, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					n(0), len(100);
	double				radius(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &n);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &radius);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &len);
	
	if ( radius < 1 ) return returnObj;
	
	int					x;
	char				string[64];
	double				j, scale = TWOPI*radius/realsizeX;
	
	snprintf(string, 60, "%d", len);
	Tcl_AppendToObj(returnObj, string, strlen(string));
	for ( x=0; x<len; x++ ) {
		j = jn(n, (x - len/2)*scale);
		snprintf(string, 60, " %g", j*j);
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}


