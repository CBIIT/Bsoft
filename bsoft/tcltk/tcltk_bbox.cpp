/**
@file	tcltk_bbox.cpp
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20200512
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bbox.h"
#include "tcltk_bimage.h"
#include "mg_img_proc.h"
#include "mg_select.h"
#include "mg_particle_select.h"
#include "mg_pick.h"
#include "mg_ctf.h"
#include "mg_extract.h"
#include "rwmg.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern Bimage* 	imglist;
extern Bimage*	imgtemp;

// Internal function prototypes
Tcl_Obj*	box_count(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_count(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_ids(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_location(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_fom(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_fom_min(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_fom_max(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_select(Bparticle* part, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_select(Bparticle* part, Bbadarea* bad, Vector3<float> rad, double badrad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_select_min(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_select_max(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_move(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_create(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_renumber(Bparticle* part, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_center(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_update_template(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_pickcc(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_pickvar(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	box_extract(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);

/*
@brief 	Functions for particle picking.
	Box interface syntax:
		Bmg box [item] [action] [params]
	Item: Must contain "Micrograph" or "Reconstruction"
	Actions:
		count [n] [|"bad"]				counts the number of particles or bad areas
		ids [fc|"bad"] [fi]				returns the identifiers of particles or bad areas
		location [id]					returns the location of a particle or bad area
		fom [id] [fi]					returns the fom of a particle with the indicated index
		fom_min							returns minimum FOM value
		fom_max							returns maximum FOM value
		select [id] [s]					sets the selection number of a particle
		select [x] [y] [z]				returns the identifier at the given location
		select_min						returns minimum selection value
		select_max						returns maximum selection value
		move [id] [dx] [dy] [dz]		updates a particle location by the given vector
		create [x] [y] [z] [n] [t] [s]	creates a new particle at the given location
		renumber [n]					renumbers particles
		delete [all|mg|rec|id] [n]		deletes particles or bad areas
		template						updates picked template
		pickcc [hr] [lr] [fmn] [fmx] [ed] [b]	picks particles by cross-correlation with a template
		pickvar [ak] [vk] [ed] [b] [ns]	picks particles from a variance map
		extract [f] [of] [bf] [nf] [ft] [fv] [mw] [sf] [p]	extracts particles from a micrograph
	Parameters:
		id			identifier
		t			box type: 1=particle, -1=bad area
		fc			FOM cutoff
		fmn			FOM minimum
		fmx			FOM maximum
		fi			FOM index (0-4)
		s			selection number
		x			x location
		y			y location
		z			z location
		n			image number
		dx			change in x
		dy			change in y
		dz			change in z
		r			particle radius
		e			gaussian edge for picking particles
		hr			high resolution limit
		lr			low resolution limit
		ak			averaging kernel edge size
		av			variance kernel edge size
		ed			exclusion distance for picking particles
		b			binning for picking particles
		ns			number of standard deviations for picking particles
		f			file name
		of			flag indicating odd box size
		bf			background correction flag
		nf			normalization flag
		ft			fill type
		fv			fill value
		mw			filament mask width
		sf			split flag
		p			particle path
*/
Tcl_Obj*	do_box(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_box: Item: " << item << " (" << item.length() << ")" << endl;

	if ( item != "all" && !item.contains("Field") && !item.contains("Micrograph") &&
			!item.contains("Reconstruction") ) {
		cerr << "Error in do_box: Item " << item << " must be a micrograph or reconstruction!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2, base, delete_what;
	double				badrad(0);
	Vector3<float>		radius;
	
	if ( item == "all" ) {
	} else if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
		if ( field ) mg = field->mg;
		if ( !mg )
			cerr << "Error in do_box: Field " << id << " not found!" << endl;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next ) ;
			if ( mg ) break;
		}
		if ( !mg )
			cerr << "Error in do_box: Micrograph " << id << " not found!" << endl;
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( !rec )
			cerr << "Error in do_box: Reconstruction " << id << " not found!" << endl;
	} else {
		cerr << "Error in do_box: Item " << item << " must be a field, micrograph or reconstruction!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}
	
	if ( mg ) {
		if ( mg->fmg.length() )
			filename = mg->fmg;
		else if ( mg->fframe.length() )
			filename = mg->fframe;
		part = mg->part;
		bad = mg->bad;
		radius = mg->box_size/2;
		badrad = mg->bad_radius;
	} else if ( rec ) {
		filename = rec->frec;
		part = rec->part;
		bad = rec->bad;
		radius = rec->box_size/2;
		badrad = rec->bad_radius;
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
		cout << "DEBUG do_box: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		if ( item == "all" )
			returnObj = box_count(project, objc, objv);
		else
			returnObj = box_count(part, bad, objc, objv);
	} else if ( action == "ids" ) {
		returnObj = box_ids(part, bad, objc, objv);
	} else if ( action == "location" ) {
		returnObj = box_location(part, bad, objc, objv);
	} else if ( action == "fom" ) {
		returnObj = box_fom(part, bad, objc, objv);
	} else if ( action == "fom_min" ) {
		returnObj = box_fom_min(project, objc, objv);
	} else if ( action == "fom_max" ) {
		returnObj = box_fom_max(project, objc, objv);
	} else if ( action == "select" ) {
		if ( objc < 7 )
			returnObj = box_select(part, objc, objv);
		else
			returnObj = box_select(part, bad, radius, badrad, objc, objv);
	} else if ( action == "select_min" ) {
		returnObj = box_select_min(project, objc, objv);
	} else if ( action == "select_max" ) {
		returnObj = box_select_max(project, objc, objv);
	} else if ( action == "move" ) {
		returnObj = box_move(part, bad, objc, objv);
	} else if ( action == "create" ) {
		if ( mg ) returnObj = box_create(mg, objc, objv);
		else if ( rec ) returnObj = box_create(rec, objc, objv);
	} else if ( action == "renumber" ) {
		returnObj = box_renumber(part, objc, objv);
	} else if ( action == "delete" ) {
		if ( objc > 4 ) delete_what = Tcl_GetStringFromObj(objv[4], NULL);
		if ( delete_what == "all" || delete_what == "mg" || delete_what == "rec" )
			returnObj = box_delete(project, objc, objv);
		else if ( mg )
			returnObj = box_delete(mg, objc, objv);
		else if ( rec )
			returnObj = box_delete(rec, objc, objv);
	} else if ( action == "center" ) {
		if ( p ) returnObj = box_center(mg, rec, p, objc, objv);
	} else if ( action == "template" ) {
		if ( p ) returnObj = box_update_template(mg, rec, p, objc, objv);
	} else if ( action == "pickcc" ) {
		if ( p ) returnObj = box_pickcc(mg, rec, p, objc, objv);
	} else if ( action == "pickvar" ) {
		if ( p ) returnObj = box_pickvar(mg, rec, p, objc, objv);
	} else if ( action == "extract" ) {
		if ( p ) returnObj = box_extract(mg, rec, p, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
		bexit(-1);
	}

	return returnObj;
}

Tcl_Obj*	box_count(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nbox = project_count_mg_particles(project) + project_count_rec_particles(project);
	
	Tcl_SetIntObj(returnObj, nbox);
	
	return returnObj;
}

Tcl_Obj*	box_count(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nbox(0);

	Bstring				select_bad;
	if ( objc > 4 ) select_bad = Tcl_GetStringFromObj(objv[4], NULL);

	if ( select_bad == "bad" ) nbox = count_list((char *) bad);
	else nbox = count_list((char *) part);
	
	Tcl_SetIntObj(returnObj, nbox);
	
	return returnObj;
}

Tcl_Obj*	box_ids(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fom_index(0);
	double				fom_cut(-1);
	Bstring				select_bad;
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &fom_index);
	if ( objc > 4 ) select_bad = Tcl_GetStringFromObj(objv[4], NULL);

	if ( select_bad == "bad" ) {
		for ( ; bad; bad = bad->next ) {
			snprintf(string, MAXLINELEN, " %d", bad->id);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	} else {
		for ( ; part; part = part->next ) if ( part->fom[fom_index] >= fom_cut ) {
			snprintf(string, MAXLINELEN, " %d", part->id);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}
			
	return returnObj;
}

Tcl_Obj*	box_location(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);	
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);

	if ( id > 0 ) for ( ; part && part->id != id; part = part->next ) ;
	else for ( ; bad && bad->id != id; bad = bad->next ) ;

	if ( id > 0 && part ) {
		snprintf(string, MAXLINELEN, "%f %f %f", part->loc[0], part->loc[1], part->loc[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	} else if ( id < 0 && bad ) {
		snprintf(string, MAXLINELEN, "%f %f %f", bad->loc[0], bad->loc[1], bad->loc[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	} else {
		snprintf(string, MAXLINELEN, "-1 -1 -1");
		Tcl_SetStringObj(returnObj, string, strlen(string));
	}
	
//	cout << "location for " << id << ": " << string << endl;
	
	return returnObj;
}

Tcl_Obj*	box_fom(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), fom_index(0);
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &fom_index);

	if ( id > 0 ) {
		for ( ; part && part->id != id; part = part->next ) ;
		if ( part ) {
			snprintf(string, MAXLINELEN, "%f", part->fom[fom_index]);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
//	cout << "fom for " << id << ": " << string << endl;
	
	return returnObj;
}

Tcl_Obj*	box_fom_min(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fom_index(0);
	double				fom_min(1e30);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( fom_min > part->fom[fom_index] ) fom_min = part->fom[fom_index];
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( fom_min > part->fom[fom_index] ) fom_min = part->fom[fom_index];
	}

	Tcl_SetDoubleObj(returnObj, fom_min);
	
	return returnObj;
}

Tcl_Obj*	box_fom_max(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fom_index(0);
	double				fom_max(-1e30);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( fom_max < part->fom[fom_index] ) fom_max = part->fom[fom_index];
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( fom_max < part->fom[fom_index] ) fom_max = part->fom[fom_index];
	}

	Tcl_SetDoubleObj(returnObj, fom_max);
	
	return returnObj;
}

Tcl_Obj*	box_select(Bparticle* part, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	int					id(0), sel(-1);

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &sel);
	
	if ( id > 0 ) {
		for ( ; part && part->id != id; part = part->next ) ;
		if ( part ) {
			if ( sel >= 0 ) part->sel = sel;
			else sel = part->sel;
		}
	}

	Tcl_SetIntObj(returnObj, sel);
	
	return returnObj;
}

Tcl_Obj*	box_select(Bparticle* part, Bbadarea* bad, Vector3<float> rad, double badrad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), i(0);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	
	double				d, dmin;
	Vector3<float>		loc(x, y, z), b, vs, ve;

	// Particle box
		dmin = 2*rad.length();
		for ( ; part; part = part->next ) {
			d = loc.distance(part->loc);
			if ( dmin > d ) {
				dmin = d;
				i = part->id;
				b = part->loc;
			}
		}
//		if ( dmin <= rad.length() ) id = i;
		if ( i ) {
			vs = b - rad;
			ve = b + rad;
			if ( loc.within(vs, ve) ) id = i;
		}
		
	// Bad area
		for ( ; bad; bad = bad->next ) {
			d = loc.distance(bad->loc);
			if ( dmin > d ) {
				dmin = d;
				i = bad->id;
			}
		}
		if ( id == 0 ) if ( dmin <= badrad ) id = i;
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_select_min(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				sel_min(1e30);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( sel_min > part->sel ) sel_min = part->sel;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( sel_min > part->sel ) sel_min = part->sel;
	}

	Tcl_SetDoubleObj(returnObj, sel_min);
	
	return returnObj;
}

Tcl_Obj*	box_select_max(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				sel_max(0);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( sel_max < part->sel ) sel_max = part->sel;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( sel_max < part->sel ) sel_max = part->sel;
	}

	Tcl_SetDoubleObj(returnObj, sel_max);
	
	return returnObj;
}

Tcl_Obj*	box_move(Bparticle* part, Bbadarea* bad, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);
	double				dx(0), dy(0), dz(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &dx);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &dy);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &dz);

	Vector3<float>		d = Vector3<float>(dx, dy, dz);

	if ( id > 0 ) {
		for ( ; part && part->id != id; part = part->next ) ;
		if ( part ) part->loc += d;
	} else if ( id < 0 ) {
		for ( ; bad && bad->id != id; bad = bad->next ) ;
		if ( bad ) bad->loc += d;
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), sel(1);
	double				x, y, z;
	Bstring				select_bad;
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	if ( objc > 7 ) select_bad = Tcl_GetStringFromObj(objv[7], NULL);
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &sel);

	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;

	if ( select_bad != "bad" ) {
		for ( part = mg->part; part; part = part->next )
			if ( id < part->id ) id = part->id;
		part = particle_add(&mg->part, ++id);
		part->loc = Vector3<float>(x, y, z);
		part->ori = mg->box_size/2;
		part->sel = sel;
	} else {
		for ( bad = mg->bad; bad; bad = bad->next )
			if ( id > bad->id ) id = bad->id;
		bad = (Bbadarea *) add_item((char **) &mg->bad, sizeof(Bbadarea));
		bad->id = --id;
		bad->loc = Vector3<float>(x, y, z);
	}
	
//	cout << "Creating a box " << id << " at " << x << "," << y << "," << z << endl;

	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_create(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), sel(1);
	double				x, y, z;
	Bstring				select_bad;
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	if ( objc > 7 ) select_bad = Tcl_GetStringFromObj(objv[7], NULL);
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &sel);

	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;

	if ( select_bad != "bad" ) {
		for ( part = rec->part; part; part = part->next )
			if ( id < part->id ) id = part->id;
		part = particle_add(&rec->part, ++id);
		part->loc = Vector3<float>(x, y, z);
		part->ori = rec->box_size/2;
		part->sel = sel;
	} else {
		for ( bad = rec->bad; bad; bad = bad->next )
			if ( id > bad->id ) id = bad->id;
		bad = (Bbadarea *) add_item((char **) &rec->bad, sizeof(Bbadarea));
		bad->id = --id;
		bad->loc = Vector3<float>(x, y, z);
	}
	
//	cout << "Creating a box " << id << " at " << x << "," << y << "," << z << endl;

	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_renumber(Bparticle* part, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					i;
	
	for ( i=0; part; part = part->next ) part->id = ++i;

	return returnObj;
}

Tcl_Obj*	box_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fom_index(0);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = project->rec;
	Bparticle*			part = NULL;
	
//	cout << "Box deletion action: " << action << endl;

	if ( action == "all" || action == "mg" ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				particle_kill(mg->part);
				mg->part = NULL;
				kill_list((char *) mg->bad, sizeof(Bbadarea));
				mg->bad = NULL;
			}
		}
	}

	if ( action == "all" || action == "rec" ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			particle_kill(rec->part);
			rec->part = NULL;
			kill_list((char *) rec->bad, sizeof(Bbadarea));
			rec->bad = NULL;
		}
	}
	
	if ( action == "fom" ) {
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
		if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &fom_index);
//		cout << "FOM cutoff: " << fom_cut << tab << project->select << endl;
		if ( project->select < 1 ) {
			for ( field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg; mg = mg->next ) {
					for ( part = mg->part; part; part = part->next )
						if ( part->fom[fom_index] < fom_cut )
							part->sel = 0;
					part_delete_deselected(&mg->part);
				}
			}
		} else {
			for ( rec = project->rec; rec; rec = rec->next ) {
				for ( part = rec->part; part; part = part->next )
					if ( part->fom[fom_index] < fom_cut )
						part->sel = 0;
				part_delete_deselected(&rec->part);
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	box_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), fom_index(0);
	double				fom_cut(0);
	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;
	
	Bstring				action;
	if ( objc > 4 ) {
		action = Tcl_GetStringFromObj(objv[4], NULL);
		Tcl_GetIntFromObj(NULL, objv[4], &id);
	}
	
//	cout << "Box deletion action: " << action << endl;

	if ( mg ) {
		if ( action == "fom" ) {
			if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
			if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &fom_index);
			for ( part = mg->part; part; part = part->next )
				if ( part->fom[fom_index] < fom_cut )
					part->sel = 0;
			part_delete_deselected(&mg->part);
		} else if ( action == "current" ) {
			particle_kill(mg->part);
			mg->part = NULL;
		} else if ( id > 0 ) {
			for ( part = mg->part; part && part->id != id; part = part->next ) ;
			if ( part )
				remove_item((char **)&mg->part, (char *)part, sizeof(Bparticle));
		} else if ( id < 0 ) {
			for ( bad = mg->bad; bad && bad->id != id; bad = bad->next ) ;
			if ( bad )
				remove_item((char **)&mg->bad, (char *)bad, sizeof(Bbadarea));
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), fom_index(0);
	double				fom_cut(0);
	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;
	
	Bstring				action;
	if ( objc > 4 ) {
		action = Tcl_GetStringFromObj(objv[4], NULL);
		id = action.integer();
	}
	
//	cout << "box_delete: " << action << tab << id << endl;

	if ( rec ) {
		if ( action == "fom" ) {
			if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
			if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &fom_index);
//			cout << "fom_cut = " << fom_cut << endl;
			for ( part = rec->part; part; part = part->next )
				if ( part->fom[fom_index] < fom_cut )
					part->sel = 0;
			part_delete_deselected(&rec->part);
		} else if ( action == "current" ) {
			particle_kill(rec->part);
			rec->part = NULL;
		} else if ( id > 0 ) {
			for ( part = rec->part; part && part->id != id; part = part->next ) ;
			if ( part )
				remove_item((char **)&rec->part, (char *)part, sizeof(Bparticle));
		} else if ( id < 0 ) {
			for ( bad = rec->bad; bad && bad->id != id; bad = bad->next ) ;
			if ( bad )
				remove_item((char **)&rec->bad, (char *)bad, sizeof(Bbadarea));
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	box_center(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	box_update_template(mg, rec, p, objc, objv);
	
	int					n(0);
	double				hires(p->sampling(0)[0]*4), lores(200), tlim(mg->box_size[0]/8);
	Bparticle*			part = NULL;
	Bimage*				p1;
	Vector3<long>		size;
	Vector3<double>		origin;
	
//	cout << "mg box size = " << mg->box_size << endl;
	
	if ( mg ) {
		size = mg->box_size;
		origin = mg->box_size/2;
		part = mg->part;
	} else if ( rec ) {
		size = rec->box_size;
		origin = rec->box_size/2;
		part = rec->part;
	}
	
	tlim = size[0]/8;
		
	for ( ; part; part = part->next ) {
		p1 = p->extract(0, part->loc - origin, size);
		p1->find_shift(imgtemp, NULL, hires, lores, tlim, 0, 1);
		part->loc += imgtemp->image->origin() - p1->image->origin();
		part->fom[0] = imgtemp->image->FOM();
		delete p1;
		n++;
	}
	
	Tcl_SetIntObj(returnObj, n);
	
	return returnObj;
}

Tcl_Obj*	box_update_template(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					n(0), img_num(0), nn;
	Bparticle*			part = NULL;
	Bimage*				p1;
	Vector3<long>		size;
	Vector3<double>		origin;
	
	delete imgtemp;
	imgtemp = NULL;

	if ( mg ) {
		img_num = mg->img_num;
		size = mg->box_size;
		origin = mg->box_size/2;
		origin[2] = 0;
		part = mg->part;
	} else if ( rec ) {
		size = rec->box_size;
		origin = rec->box_size/2;
		part = rec->part;
	}
	
	for ( ; part; part = part->next ) {
		p1 = p->extract(img_num, part->loc - origin, size, FILL_AVERAGE);
		p1->change_type(Float);
		if ( imgtemp ) {
			imgtemp->add(p1);
			delete p1;
		} else imgtemp = p1;
		n++;
	}
	
	if ( imgtemp ) {
		imgtemp->calculate_background();
		for ( nn=0; nn<imgtemp->images(); nn++ ) {
			p1 = imgtemp->extract(nn);
			p1->find_center(NULL, 1000, 20*p->sampling(0)[0], size[0]/10, 0, 1);
//			cout << "shift=" << origin - p1->image->origin() << endl;
			imgtemp->shift(nn, origin - p1->image->origin(), FILL_BACKGROUND);
			delete p1;
		}
		imgtemp->origin(origin);
		imgtemp->identifier("thetemp");
		imgtemp->rescale_to_avg_std(0, 1);
	}
	
	Tcl_SetIntObj(returnObj, n);
	
	return returnObj;
}

Tcl_Obj*	box_pickcc(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
	
	int					bin(1);
	double				hires(p->sampling(0)[0]*4), lores(200);
	double				fommin(0), fommax(1e30), excl_dist(0), threshold(0);
	Bimage*				pmask = NULL;
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &hires);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &lores);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &fommin);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &fommax);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &excl_dist);
	if ( objc > 9 ) Tcl_GetIntFromObj(NULL, objv[9], &bin);

	if ( imgtemp ) {
		if ( mg ) {
			mg->box_size = imgtemp->size();
			mg->part = particles_pick_cc(p, imgtemp, pmask, hires, lores, fommin, fommax, excl_dist, bin);
			threshold = mg->fom = imgtemp->image->FOM();
		} else if ( rec ) {
			rec->box_size = imgtemp->size();
			rec->part = particles_pick_cc(p, imgtemp, pmask, hires, lores, fommin, fommax, excl_dist, bin);
			threshold = rec->fom = imgtemp->image->FOM();
		}
	}
	
	Tcl_SetDoubleObj(returnObj, threshold);
	
	return returnObj;
}

Tcl_Obj*	box_pickvar(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
	
	int					bin(1), avg_kernel(21), var_kernel(301);
//	double				excl_dist(0), threshold(0), nsig(5);
	double				excl_dist(0), threshold(0), cutmin(0.1), cutmax(0.2);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &avg_kernel);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &var_kernel);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &excl_dist);
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &bin);
//	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &nsig);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &cutmin);
	if ( objc > 9 ) Tcl_GetDoubleFromObj(NULL, objv[9], &cutmax);

//	cout << "excl_dist=" << excl_dist << endl;
	
	if ( mg ) {
//		mg->part = particles_pick_var(p, avg_kernel, var_kernel, nsig,
//						mg->box_size[0]/2, excl_dist, bin);
		mg->part = particles_pick_var(p, avg_kernel, var_kernel, cutmin,
						cutmax, mg->box_size[0]/2, excl_dist, bin);
		threshold = mg->fom ;
	} else if ( rec ) {
//		rec->part = particles_pick_var(p, avg_kernel, var_kernel, nsig,
//						rec->box_size[0]/2, excl_dist, bin);
		rec->part = particles_pick_var(p, avg_kernel, var_kernel, cutmin,
						cutmax, rec->box_size[0]/2, excl_dist, bin);
		threshold = rec->fom;
	}
	
	Tcl_SetDoubleObj(returnObj, threshold);
	
	return returnObj;
}

Tcl_Obj*	box_extract(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return 0;
		
	Bstring			filename, partpath;
	int				back_flag(0), norm_flag(0), mask_width(0), split(0);
	int				filltype(2), npart(0);
	double			scale(1), fill(0);
	
	if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &scale);	/* Extraction scale */
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &back_flag);	/* Flag for background correction */
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &norm_flag);	/* Flag for normalization */
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &filltype);	/* Fill value */
	if ( objc > 9 ) Tcl_GetDoubleFromObj(NULL, objv[9], &fill);
	if ( objc > 10 ) Tcl_GetIntFromObj(NULL, objv[10], &mask_width);	/* Filament mask width */
	if ( objc > 11 ) Tcl_GetIntFromObj(NULL, objv[11], &split);		/* Flag to split particles into individual images */
	if ( objc > 12 ) partpath = Tcl_GetStringFromObj(objv[12], NULL);	/* Particle image path */

	if ( filltype == 1 ) fill = p->average();
	else if ( filltype == 2 ) {
		if ( fabs(p->background(long(0))) < 1e-6 ) p->calculate_background();
		fill = p->background(long(0));
	}
	
	if ( filename.length() ) {
		if ( mg ) {
			if ( split ) {
				particle_setup_filenames(mg->part, filename, partpath);
				mg->fpart = 0;
			} else {
				if ( partpath.length() ) {
					filename = partpath + "/" + filename;
					mkdir(partpath.c_str(), (mode_t)0755);
				}
				mg->fpart = filename;
			}
			npart = micrograph_extract_particles(mg, p, scale,
						back_flag, norm_flag, fill, mask_width);
		} else if ( rec ) {
			if ( split ) {
				particle_setup_filenames(rec->part, filename, partpath);
				rec->fpart = 0;
			} else {
				if ( partpath.length() ) {
					filename = partpath + "/" + filename;
					mkdir(partpath.c_str(), (mode_t)0755);
				}
				rec->fpart = filename;
			}
			npart = reconstruction_extract_particles(rec, p, scale,
						back_flag, norm_flag, fill, mask_width);
		}
	}
	
	filename = 0;
	partpath = 0;

	Tcl_SetIntObj(returnObj, npart);
	
	return returnObj;
}

