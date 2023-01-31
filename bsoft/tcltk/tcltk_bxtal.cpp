/**
@file	tcltk_bxtal.cpp
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20181030
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bxtal.h"
#include "mg_xtal.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern Bimage* 	imglist;

// Internal function prototypes
Tcl_Obj*	spot_count(Bmicrograph* mg);
Tcl_Obj*	spot_ids(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_reindex(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_location(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_fom(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_select(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_move(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_unitcell_vectors(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_generate(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	spot_mask(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);


/*
@brief 	Functions for 2D crystallographic indexing.
	Diffraction spot interface syntax:
		Bmg spot [item] [action] [params]
	Actions:
		count						counts the number of spots
		ids							returns the spot indices (3 per spot)
		reindex [h1] [k1] [l1] [h2] [k2] [l2]	re-indexes a spot
		location [h] [k] [l]		returns the location of a spot
		fom [h] [k] [l]				returns the FOM of a spot
		select [x] [y] [z]			returns the spot indices at the given location
		move [h] [k] [l] [dx] [dy] [dz]	updates a spot location by the given vector
		create [h] [k] [l] [x] [y] [z]	creates a new spot at the given location
		delete [all|mg|rec|h] [k] [l]	deletes spots
		unitcell					calculates unit cell vectors
		generate [m]				generates spots
		mask [r]					masks spots
	Parameters:
		h k l		Miller indices
		x			x location
		y			y location
		z			z location
		dx			change in x
		dy			change in y
		dz			change in z
		m			limit in reciprocal space for generating spots
		r			mask radius
*/

Tcl_Obj*	do_spot(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_spot: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	if ( !item.contains("Field") && !item.contains("Micrograph") ) {
		cerr << "Error in do_spot: Item " << item << " must be a micrograph!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2, base, delete_what;
	
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
	}
	
	if ( filename.length() ) {
		base = filename.base();
		for ( p = imglist; p; p = p->next ) {
			filename2 = p->file_name();
			if ( filename2.base() == base ) break;
		}
	}

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_spot: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		returnObj = spot_count(mg);
	} else if ( action == "ids" ) {
		returnObj = spot_ids(mg, objc, objv);
	} else if ( action == "reindex" ) {
		returnObj = spot_reindex(mg, objc, objv);
	} else if ( action == "location" ) {
		returnObj = spot_location(mg, objc, objv);
	} else if ( action == "fom" ) {
		returnObj = spot_fom(mg, objc, objv);
	} else if ( action == "select" ) {
		returnObj = spot_select(mg, objc, objv);
	} else if ( action == "move" ) {
		returnObj = spot_move(mg, objc, objv);
	} else if ( action == "create" ) {
		returnObj = spot_create(mg, objc, objv);
	} else if ( action == "delete" ) {
		returnObj = spot_delete(mg, objc, objv);
	} else if ( action == "unitcell" ) {
        returnObj = spot_unitcell_vectors(mg, objc, objv);
	} else if ( action == "generate" ) {
		returnObj = spot_generate(mg, p, objc, objv);
	} else if ( action == "mask" ) {
		if ( p ) returnObj = spot_mask(mg, p, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	if ( !returnObj ) {
		returnObj = Tcl_NewObj();
		Tcl_SetIntObj(returnObj, 0);
	}

	return returnObj;
}

Tcl_Obj*	spot_count(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nspot(0);
	Bstrucfac*			sf;

	if ( mg )
		for ( sf = mg->sf; sf; sf = sf->next ) nspot++;
	
	Tcl_SetIntObj(returnObj, nspot);
	
	return returnObj;
}

Tcl_Obj*	spot_ids(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				fom_cut(0);
	char				string[MAXLINELEN] = "";
	Bstrucfac*			sf;

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	if ( mg ) {
		for ( sf = mg->sf; sf; sf = sf->next ) if ( sf->fom >= fom_cut ) {
			snprintf(string, MAXLINELEN, " %d %d %d", sf->index[0], sf->index[1], sf->index[2]);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}
		
	return returnObj;
}

Tcl_Obj*	spot_reindex(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					h1(0), k1(0), l1(0), h2(0), k2(0), l2(0);	
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h1);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k1);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l1);
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &h2);
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &k2);
	if ( objc > 9 ) Tcl_GetIntFromObj(NULL, objv[9], &l2);

	Bstrucfac*			sf;
	Vector3<int>		id1(h1, k1, l1), id2(h2, k2, l2);

	if ( mg ) {
		for ( sf = mg->sf; sf && sf->index != id2; sf = sf->next ) ;
		if ( sf ) {
			id2 = id1;
		} else {
			for ( sf = mg->sf; sf && sf->index != id1; sf = sf->next ) ;
			if ( sf ) sf->index = id2;
		}
	}

	snprintf(string, MAXLINELEN, "%d %d %d", id2[0], id2[1], id2[2]);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}

Tcl_Obj*	spot_location(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					h(0), k(0), l(0);	
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l);

	Bstrucfac*			sf = NULL;
	Vector3<int>		id(h, k, l);

	if ( mg )
		for ( sf = mg->sf; sf && sf->index != id; sf = sf->next ) ;

	if ( sf ) {
		Vector3<double>		loc = sf->loc + mg->origin;
		snprintf(string, MAXLINELEN, "%f %f %f", loc[0], loc[1], loc[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
//		cout << id << tab << sf->loc << endl;
	}
	
	return returnObj;
}

Tcl_Obj*	spot_fom(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					h(0), k(0), l(0);	

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l);

	Bstrucfac*			sf = NULL;
	Vector3<int>		id(h, k, l);
	double				fom(1);

	if ( mg )
		for ( sf = mg->sf; sf && sf->index != id; sf = sf->next ) ;

	if ( sf ) fom = sf->fom;
	
	Tcl_SetDoubleObj(returnObj, fom);
	
	return returnObj;
}

Tcl_Obj*	spot_select(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	
	double				d, dmin, rad(0);
	char				string[64] = "";
	Vector3<int>		ind;
	Vector3<double>		loc(x, y, z);
	Bstrucfac*			sflist = NULL;
	Bstrucfac*			sf;

	if ( mg ) {
		sflist = mg->sf;
		rad = mg->sf_radius;
	}

	if ( sflist ) {
		loc -= mg->origin;
		dmin = 2*rad;
		for ( sf = sflist; sf; sf = sf->next ) {
			d = loc.distance(sf->loc);
			if ( dmin > d ) {
				dmin = d;
				ind = sf->index;
			}
		}
		if ( dmin <= rad ) {
			snprintf(string, 60, " %d %d %d", ind[0], ind[1], ind[2]);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	spot_move(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					h(0), k(0), l(0);
	double				dx(0), dy(0), dz(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &dx);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &dy);
	if ( objc > 9 ) Tcl_GetDoubleFromObj(NULL, objv[9], &dz);

	char				string[64] = "";
	Bstrucfac*			sflist = NULL;
	Bstrucfac*			sf;
	Vector3<int>		id(h, k, l);
	Vector3<double>		d(dx, dy, dz);

	if ( mg ) sflist = mg->sf;

	if ( sflist ) {
		for ( sf = sflist; sf && sf->index != id; sf = sf->next ) ;
		if ( sf ) {
			sf->loc += d;
			snprintf(string, 60, " %d %d %d", sf->index[0], sf->index[1], sf->index[2]);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	spot_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					h(0), k(0), l(0);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &x);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &y);
	if ( objc > 9 ) Tcl_GetDoubleFromObj(NULL, objv[9], &z);

	char				string[64] = "";
	Bstrucfac*			sf = NULL;
	Vector3<int>		id(h, k, l);

//	cout << "origin=" << mg->origin << endl;
	if ( mg )
		sf = (Bstrucfac *) add_item((char **) &mg->sf, sizeof(Bstrucfac));

	if ( sf ) {
		sf->index = id;
		sf->loc = Vector3<double>(x, y, z) - mg->origin;
		sf->amp = 1;
		sf->fom = 1;
		sf->sel = 1;
		snprintf(string, 60, " %d %d %d", sf->index[0], sf->index[1], sf->index[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
//		cout << id << tab << sf->loc << endl;
	}
	
	return returnObj;
}

Tcl_Obj*	spot_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !mg ) return returnObj;
	
	int					h(0), k(0), l(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	char				string[64] = "";
	Bstrucfac*			sf;

	if ( action == "all" ) {
		kill_list((char *) mg->sf, sizeof(Bstrucfac));
		mg->sf = NULL;
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &h);
		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &k);
		if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &l);

		Vector3<int>		id(h, k, l);
	
		if ( h!=0 || k!=0 || l!=0 ) {
			for ( sf = mg->sf; sf && sf->index != id; sf = sf->next ) ;
			if ( sf )
				remove_item((char **)&mg->sf, (char *)sf, sizeof(Bstrucfac));
		}
		snprintf(string, 60, "%d %d %d", h, k, l);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	spot_unitcell_vectors(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( mg ) mg_unitcell_vectors(mg);

	return returnObj;
}

Tcl_Obj*	spot_generate(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				resolution(20);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &resolution);

	if ( verbose ) {
		cout << "sampling=" << p->sampling(0)[0] << tab << p->sampling(0)[1] << endl;
		cout << "real_size=" << p->real_size() << endl;
	}
	
	if ( mg ) mg_generate_reflections(mg, p->real_size(), resolution);

	return returnObj;
}

Tcl_Obj*	spot_mask(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				radius(3);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &radius);
	
	if ( mg ) img_mask_reflections(p, mg->sf, radius);

	return returnObj;
}


