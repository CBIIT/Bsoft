/**
@file	tcltk_bmarker.cpp
@brief	A shared object to manage markers in micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20200513
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bmarker.h"
#include "mg_img_proc.h"
#include "mg_extract.h"
#include "mg_tomography.h"
#include "mg_tomo_track.h"
#include "rwmg.h"
#include "mg_tomography.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern Bimage* 	imglist;

// Internal function prototypes
Tcl_Obj*	marker_count(Bproject* project);
Tcl_Obj*	marker_count(Bfield* field);
Tcl_Obj*	marker_count(Bmicrograph* mg);
Tcl_Obj*	marker_count(Breconstruction* rec);
Tcl_Obj*	marker_ids(Bproject* project, double fom_cut);
Tcl_Obj*	marker_ids(Bfield* field, double fom_cut);
Tcl_Obj*	marker_ids(Bmicrograph* mg, double fom_cut);
Tcl_Obj*	marker_ids(Breconstruction* rec, double fom_cut);
Tcl_Obj*	marker_list(Bproject* project);
Tcl_Obj*	marker_list(Bfield* field);
Tcl_Obj*	marker_list(Bmicrograph* mg);
Tcl_Obj*	marker_list(Breconstruction* rec);
Tcl_Obj*	marker_location(Bmarker* mark);
Tcl_Obj*	marker_residual(Bproject* project);
Tcl_Obj*	marker_residual(Bmarker* mark);
Tcl_Obj*	marker_fom(Bmarker* mark);
Tcl_Obj*	marker_fom_maximum(Bmarker* mark);
Tcl_Obj*	marker_select(Bmarker* marklist, double rad, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_select_rectangle(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_get_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_set_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_clear_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_move(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_create(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_center(Bmarker* mark, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_generate_from_seed(Bproject* project);
Tcl_Obj*	marker_find(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_snap(Bmicrograph* mg, Bimage* p, long mid, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_extract(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_extract(Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_delete(Bfield* field, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_edge(Bmarker* mark, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	marker_accept(Bmarker* mark, int objc, Tcl_Obj *CONST objv[]);

Tcl_Obj*	do_tomo_thickness(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*		returnObj = Tcl_NewObj();
	
	int				adjust_tilt(0);
	double			lambda(2260);
	Bstring			psfile;
	
	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &lambda);		// Proportionality parameter
	if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &adjust_tilt);		// Tilt ddjustment flag
	if ( objc > 4 ) psfile = Tcl_GetStringFromObj(objv[4], NULL);		// Postscript file name

//	cout << "Lambda = " << lambda << tab << adjust_tilt << endl;
	
	Bplot*			plot = project_intensity_plot(project);
	if ( !plot ) return returnObj;
	
	double			tlr = project_fit_intensities(project, plot, adjust_tilt);
	double			thickness = lambda*tlr;
	
//	cout << "tlr = " << tlr << endl;

	Bstring			txt;
	txt = Bstring(thickness, "Thickness: %lg");
	plot->page(0).add_text(txt);

	if ( psfile.length() )
		ps_plot(psfile, plot);

	Tcl_SetDoubleObj(returnObj, thickness);
	
	return returnObj;
}

int			do_tomo_transfer_seed(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	double			rot_start(-180), rot_end(180);	// Rotation search limits
	double			rot_step(0);					// Rotation search step size for seed transfer
	double			hi_res(0), lo_res(1000);		// Default resolution range
	double			shift_limit(-1);				// Micrograph shift search limit
	int				refine(0);

	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &rot_start);
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &rot_end);
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &rot_step);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &hi_res);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &lo_res);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &shift_limit);
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &refine);
	
	rot_step *= M_PI/180.0;
	rot_start *= M_PI/180.0;
	rot_end *= M_PI/180.0;

	if ( rot_step > 0 ) {
		project_transfer_seed(project, rot_start, rot_end, rot_step, hi_res, lo_res, shift_limit);
		if ( refine )
			project_refine_markers(project, hi_res, lo_res);
	}
	
	return 0;
}

int			do_tomo_findaxis(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	double		axis(0), axis_start(-180), axis_end(180), axis_step(1);
	double		hi_res(0), lo_res(1000), shift_limit(-1);
	
	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &axis);			// Target tilt axis
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &axis_start);	// Starting angle to test
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &axis_end);		// Ending angle to test
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &axis_step);	// Angular step size
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &hi_res);		// High resolution limit
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &lo_res);		// Low resolution limit
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &shift_limit);	// Maximum shift limit

	axis *= M_PI/180.0;
	axis_start *= M_PI/180.0;
	axis_end *= M_PI/180.0;
	axis_step *= M_PI/180.0;

	project_find_tilt_axis(project, axis, axis_start, axis_end, axis_step, hi_res, lo_res, shift_limit);
	
	return 0;
}

int			do_tomo_track(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	int			iter(5), thickness(500), cc_type(0), recenter(0);
	double		hi_res(0), lo_res(500), shift_limit(-1);
	Bstring		paramfile("track.star");
	
	if ( objc > 2 ) Tcl_GetIntFromObj(NULL, objv[2], &iter);			// Maximum number of iterations
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &hi_res);		// High resolution limit
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &lo_res);		// Low resolution limit
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &shift_limit);	// Maximum shift limit
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &thickness);		// Maximum thickness
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &cc_type);			// Flag to cross correlate
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &recenter);		// Flag to recenter z coordinates

    ofstream 	out("track.log");
    streambuf 	*coutbuf = cout.rdbuf(); //save old buf
    cout.rdbuf(out.rdbuf()); //redirect std::cout to out!

	project_track_markers(project, hi_res, lo_res, shift_limit, thickness, cc_type, iter, 1, recenter, paramfile);

    cout.rdbuf(coutbuf); //reset to standard output again
	
	return 0;
}

int			do_tomo_refine(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
//	int			ref_view(0), ref_ori(0), ref_scale(0);
	int			ref_iter(1);
	double		hi_res(20), lo_res(1000), ref_tol(1e-4);
	Bstring		refop;

	project_check_markers(project, 14);
	project_tomo_residuals(project, 0);
	project_sort_markers_by_id(project);
	
	/* Type of operations */
	if ( objc > 2 ) refop = Tcl_GetStringFromObj(objv[2], NULL);

	if ( refop[0] == 'm' ) {
		if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &hi_res);	// High resolution limit
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &lo_res);	// Low resolution limit
		project_refine_markers(project, hi_res, lo_res);
//	} else if ( refop[0] == 'z' ) {
//		project_refine_z(project);
	} else {
//		if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &ref_view);	/* Do views */
//		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &ref_ori);		/* Do origins */
//		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &ref_scale);	/* Do scales */
//		project_refine(project, ref_view, ref_ori, ref_scale);
		if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &ref_iter);	// Iterations
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &ref_tol);	// Stopping condition
		project_refine(project, ref_iter, ref_tol, refop);
	}
	
	project_check_markers(project, 14);
	project_tomo_residuals(project, 0);
	project_calculate_angles(project);
	
	return 0;
}

int			do_tomo_refine_one(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	int			id(0);
	double		hi_res(20), lo_res(1000);
	
	/* Type of operations */
	if ( objc > 2 ) Tcl_GetIntFromObj(NULL, objv[2], &id);
	if ( id < 1 ) return 0;
	
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &hi_res);	/* High resolution limit */
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &lo_res);	/* Low resolution limit */
	project_refine_one_marker(project, id, hi_res, lo_res);
	
	project_check_markers(project, 14);
	project_tomo_residuals(project, 0);
	
	return 0;
}



/*
@brief 	Functions for tomographic marker manipulation.
	Marker interface syntax:
		Bmg marker [item] [action] [params]
	Actions:
		count							counts the number of markers
		ids [f]							returns the identifiers of markers selected above the given FOM cutoff
		list							returns a list of markers
		location [id]					returns the location and error vector of a marker
		z_range							returns the maximum difference in z coordinates
		residual [id]					returns marker residual or residual of all markers
		fom [id]						returns the FOM of a marker
		select [x] [y] [z]				returns the identifier at the given location
		selectrect [x1] [y1] [x2] [y2]	returns a list of identifiers within the rectangle
		getflag [id]					sets the selection flag of a marker
		setflag [id] [s]				sets the selection flag of a marker
		clearflag [id]					clears the selection flag of a marker
		move [id] [dx] [dy] [dz]		updates a marker location by the given vector
		create [x] [y] [z]				creates a new marker at the given location
		renumber						renumbers markers
		center							centers markers
		generate						generates missing markers from an existing seed
		find [r]						finds markers
		snap [id]						snaps a marker to a low density center
		extract [fn]					extracts markers
		delete [all|fom|id] [f]			deletes markers
		edge [e]						sets FOM based on edge parameter
		accept [id]						accept model location for marker
	Parameters:
		id			identifier
		x			x location
		y			y location
		z			z location
		f			FOM cutoff
		s			selection flag
		dx			change in x
		dy			change in y
		dz			change in z
		r			marker radius
		fn			new file name
		e			edge for marker selection
*/
Tcl_Obj*	do_marker(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_marker: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bmarker*			mark = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2, base, delete_what;
	int					mid(-1);
	double				radius(0), fom_cut(0);
	
	if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
		if ( field )
			filename = field->mg->fmg;
		else
			cerr << "Error in do_marker: Field " << id << " not found!" << endl;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next ) ;
			if ( mg ) break;
		}
		if ( mg ) {
			if ( mg->fmg.length() )
				filename = mg->fmg;
			else if ( mg->fframe.length() )
				filename = mg->fframe;
			mark = mg->mark;
			radius = mg->mark_radius;
		} else
			cerr << "Error in do_marker: Micrograph " << id << " not found!" << endl;
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( rec ) {
			filename = rec->frec;
			mark = rec->mark;
			radius = rec->mark_radius;
		} else
			cerr << "Error in do_marker: Reconstruction " << id << " not found!" << endl;
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
		cout << "DEBUG do_marker: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		if ( mg ) returnObj = marker_count(mg);
		else if ( field ) returnObj = marker_count(field);
		else if ( rec ) returnObj = marker_count(rec);
		else returnObj = marker_count(project);
	} else if ( action == "ids" ) {
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);
		if ( mg ) returnObj = marker_ids(mg, fom_cut);
		else if ( field ) returnObj = marker_ids(field, fom_cut);
		else if ( rec ) returnObj = marker_ids(rec, fom_cut);
		else returnObj = marker_ids(project, fom_cut);
	} else if ( action == "list" ) {
		if ( mg ) returnObj = marker_list(mg);
		else if ( field ) returnObj = marker_list(field);
		else if ( rec ) returnObj = marker_list(rec);
		else returnObj = marker_list(project);
	} else if ( action == "location" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &mid);
		for ( ; mark && mark->id != mid; mark = mark->next ) ;
		if ( mark ) returnObj = marker_location(mark);
	} else if ( action == "z_range" ) {
		rec = project->rec;
		if ( rec && rec->mark ) {
			Vector3<double>	mr = marker_range(rec->mark);
			mr *= rec->voxel_size;
			returnObj = Tcl_NewDoubleObj(mr[2]);
		}
	} else if ( action == "residual" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &mid);
		if ( mid < 0 ) {
			returnObj = marker_residual(project);
		} else {
			for ( ; mark && mark->id != mid; mark = mark->next ) ;
			if ( mark ) returnObj = marker_residual(mark);
		}
	} else if ( action == "fom" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &mid);
		for ( ; mark && mark->id != mid; mark = mark->next ) ;
		returnObj = marker_fom(mark);
	} else if ( action == "fom_maximum" ) {
		returnObj = marker_fom_maximum(mark);
	} else if ( action == "select" ) {
		if ( mark ) returnObj = marker_select(mark, radius, objc, objv);
	} else if ( action == "selectrect" ) {
		if ( mark ) returnObj = marker_select_rectangle(mark, objc, objv);
	} else if ( action == "getflag" ) {
		if ( mark ) returnObj = marker_get_selection_flag(mark, objc, objv);
	} else if ( action == "setflag" ) {
		if ( mark ) returnObj = marker_set_selection_flag(mark, objc, objv);
	} else if ( action == "clearflag" ) {
		if ( mark ) returnObj = marker_clear_selection_flag(mark, objc, objv);
	} else if ( action == "move" ) {
		if ( mark ) returnObj = marker_move(mark, objc, objv);
	} else if ( action == "create" ) {
		if ( mg ) returnObj = marker_create(mg, objc, objv);
		else if ( rec ) returnObj = marker_create(rec, objc, objv);
	} else if ( action == "renumber" ) {
		if ( mark ) markers_renumber(mark);
		else project_renumber_markers(project);
	} else if ( action == "center" ) {
		if ( mark ) returnObj = marker_center(mark, p, objc, objv);
	} else if ( action == "generate" ) {
		returnObj = marker_generate_from_seed(project);
	} else if ( action == "find" ) {
		if ( mg ) returnObj = marker_find(mg, p, objc, objv);
	} else if ( action == "snap" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &mid);
		if ( mg ) returnObj = marker_snap(mg, p, mid, objc, objv);
	} else if ( action == "extract" ) {
		if ( mg ) returnObj = marker_extract(mg, p, objc, objv);
		else if ( rec ) returnObj = marker_extract(rec, p, objc, objv);
	} else if ( action == "delete" ) {
		if ( mg ) returnObj = marker_delete(mg, objc, objv);
		else if ( field ) returnObj = marker_delete(field, objc, objv);
		else if ( rec ) returnObj = marker_delete(rec, objc, objv);
		else returnObj = marker_delete(project, objc, objv);
	} else if ( action == "edge" ) {
		if ( mark ) returnObj = marker_edge(mark, p, objc, objv);
	} else if ( action == "accept" ) {
		if ( mark ) returnObj = marker_accept(mark, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}
	
	if ( !returnObj ) {
		returnObj = Tcl_NewObj();
		Tcl_SetIntObj(returnObj, 0);
	}

	return returnObj;
}

Tcl_Obj*	marker_count(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nmark(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bmarker*			mark;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( mark = mg->mark; mark; mark = mark->next ) nmark++;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( mark = rec->mark; mark; mark = mark->next ) nmark++;
	}
	
	Tcl_SetIntObj(returnObj, nmark);
	
	return returnObj;
}

Tcl_Obj*	marker_count(Bfield* field)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nmark(0);
	Bmicrograph*		mg;
	Bmarker*			mark;

	for ( mg = field->mg; mg; mg = mg->next )
		for ( mark = mg->mark; mark; mark = mark->next ) nmark++;
	
	Tcl_SetIntObj(returnObj, nmark);
	
	return returnObj;
}

Tcl_Obj*	marker_count(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nmark(0);
	Bmarker*			mark;

	for ( mark = mg->mark; mark; mark = mark->next ) nmark++;
	
	Tcl_SetIntObj(returnObj, nmark);
	
	return returnObj;
}

Tcl_Obj*	marker_count(Breconstruction* rec)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nmark(0);
	Bmarker*			mark;

	for ( mark = rec->mark; mark; mark = mark->next ) nmark++;
	
	Tcl_SetIntObj(returnObj, nmark);
	
	return returnObj;
}

Tcl_Obj*	marker_ids(Bproject* project, double fom_cut)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = project->rec;
	Bmarker*			mark = NULL;

	if ( project->select < 1 ) {
		if ( field ) for ( mg = field->mg; mg && !mg->mark; mg = mg->next ) ;
		if ( mg ) {
			for ( mark = mg->mark; mark; mark = mark->next ) if ( mark->fom >= fom_cut ) {
				snprintf(string, MAXLINELEN, " %d", mark->id);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
		}
	} else {
		if ( rec ) {
			for ( mark = rec->mark; mark; mark = mark->next ) if ( mark->fom >= fom_cut ) {
				snprintf(string, MAXLINELEN, " %d", mark->id);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
		}
	}
		
	return returnObj;
}

Tcl_Obj*	marker_ids(Bfield* field, double fom_cut)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmicrograph*		mg = NULL;
	Bmarker*			mark = NULL;

	for ( mg = field->mg; mg && !mg->mark; mg = mg->next ) ;
	if ( mg ) {
		for ( mark = mg->mark; mark; mark = mark->next ) if ( mark->fom >= fom_cut ) {
			snprintf(string, MAXLINELEN, " %d", mark->id);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}
		
	return returnObj;
}

Tcl_Obj*	marker_ids(Bmicrograph* mg, double fom_cut)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmarker*			mark = NULL;

	for ( mark = mg->mark; mark; mark = mark->next ) if ( mark->fom >= fom_cut ) {
		snprintf(string, MAXLINELEN, " %d", mark->id);
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
		
	return returnObj;
}

Tcl_Obj*	marker_ids(Breconstruction* rec, double fom_cut)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmarker*			mark = NULL;

	for ( mark = rec->mark; mark; mark = mark->next ) if ( mark->fom >= fom_cut ) {
		snprintf(string, MAXLINELEN, " %d", mark->id);
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
		
	return returnObj;
}

Tcl_Obj*	marker_list(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bfield*				field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = project->rec;
	Bmarker*			mark = NULL;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( mark = mg->mark; mark; mark = mark->next ) {
					snprintf(string, MAXLINELEN, " %d %d %20.16f %20.16f %d", mg->img_num, mark->id, mark->res, mark->fom, mark->sel);
					Tcl_AppendToObj(returnObj, string, strlen(string));
				}
	} else {
		if ( rec ) {
			for ( mark = rec->mark; mark; mark = mark->next ) {
				snprintf(string, MAXLINELEN, " 0 %d %20.16f %20.16f %d", mark->id, mark->res, mark->fom, mark->sel);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
		}
	}
		
	return returnObj;
}

Tcl_Obj*	marker_list(Bfield* field)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmicrograph*		mg = NULL;
	Bmarker*			mark = NULL;

	for ( mg = field->mg; mg; mg = mg->next )
		for ( mark = mg->mark; mark; mark = mark->next ) {
			snprintf(string, MAXLINELEN, " %d %d %20.16f %20.16f %d", mg->img_num, mark->id, mark->res, mark->fom, mark->sel);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
		
	return returnObj;
}

Tcl_Obj*	marker_list(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmarker*			mark = NULL;

	for ( mark = mg->mark; mark; mark = mark->next ) {
		snprintf(string, MAXLINELEN, " %d %d %20.16f %20.16f %d", mg->img_num, mark->id, mark->res, mark->fom, mark->sel);
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
		
	return returnObj;
}

Tcl_Obj*	marker_list(Breconstruction* rec)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bmarker*			mark = NULL;

	for ( mark = rec->mark; mark; mark = mark->next ) {
		snprintf(string, MAXLINELEN, " 0 %d %20.16f %20.16f %d", mark->id, mark->res, mark->fom, mark->sel);
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
		
	return returnObj;
}

Tcl_Obj*	marker_location(Bmarker* mark)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";

	if ( mark ) {
		snprintf(string, MAXLINELEN, "%f %f %f %f %f %f", 
			mark->loc[0], mark->loc[1], mark->loc[2], 
			mark->err[0], mark->err[1], mark->err[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	marker_residual(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				res = project_tomo_residuals(project, 0);
	
	Tcl_SetDoubleObj(returnObj, res);
	
	return returnObj;
}

Tcl_Obj*	marker_residual(Bmarker* mark)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				res(0);
	
	if ( mark ) res = mark->res;
	
	Tcl_SetDoubleObj(returnObj, res);
	
	return returnObj;
}

Tcl_Obj*	marker_fom(Bmarker* mark)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				fom(0);

	if ( mark ) fom = mark->fom;
	
	Tcl_SetDoubleObj(returnObj, fom);
	
	return returnObj;
}

Tcl_Obj*	marker_fom_maximum(Bmarker* mark)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				fom_max(0);
	
	for ( ; mark; mark = mark->next )
		if ( fom_max < mark->fom ) fom_max = mark->fom;

	Tcl_SetDoubleObj(returnObj, fom_max);
	
	return returnObj;
}

Tcl_Obj*	marker_select(Bmarker* marklist, double radius, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), i(-1);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	
	double				d, dmin;
	Vector3<float>		loc(x, y, z);
	Bmarker*			mark;

	if ( marklist ) {
		dmin = 2*radius;
		for ( mark = marklist; mark; mark = mark->next ) {
			d = loc.distance(mark->loc);
			if ( dmin > d ) {
				dmin = d;
				i = mark->id;
			}
		}
		if ( dmin <= radius ) id = i;
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_select_rectangle(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				x1(0), y1(0), x2(0), y2(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x1);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y1);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &x2);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &y2);
	
	if ( x2 < x1 ) swap(x1, x2);
	if ( y2 < y1 ) swap(y1, y2);
	
	char				string[MAXLINELEN] = "";
	Bmarker*			mark;

	if ( marklist ) {
		for ( mark = marklist; mark; mark = mark->next ) {
			if ( mark->loc[0] >= x1 && mark->loc[0] <= x2 && 
					mark->loc[1] >= y1 && mark->loc[1] <= y2 ) {
				snprintf(string, MAXLINELEN, " %d", mark->id);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	marker_get_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), sel(-1);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	
	Bmarker*			mark;

	for ( mark = marklist; mark && mark->id != id; mark = mark->next ) ;
	
	if ( mark ) sel = mark->sel;

	Tcl_SetIntObj(returnObj, sel);
		
	return returnObj;
}

Tcl_Obj*	marker_set_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), sel(1), found(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &sel);
	
	Bmarker*			mark;

	for ( mark = marklist; mark && mark->id != id; mark = mark->next ) ;
	
	if ( mark ) mark->sel = found = sel;

	Tcl_SetIntObj(returnObj, found);
		
	return returnObj;
}

Tcl_Obj*	marker_clear_selection_flag(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), found(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	
	Bmarker*			mark;

	for ( mark = marklist; mark && mark->id != id; mark = mark->next ) ;
	
	if ( mark ) {
		mark->sel = 0;
		found = 1;
	}

	Tcl_SetIntObj(returnObj, found);
		
	return returnObj;
}

Tcl_Obj*	marker_move(Bmarker* marklist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);
	double				dx(0), dy(0), dz(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &dx);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &dy);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &dz);

	Bmarker*			mark;
	Vector3<float>		d = Vector3<float>(dx, dy, dz);

	if ( marklist ) {
		for ( mark = marklist; mark && mark->id != id; mark = mark->next ) ;
		if ( mark ) {
			mark->loc += d;
			mark->err -= d;
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

int			marker_create(Bmarker* mark, int objc, Tcl_Obj *CONST objv[])
{
//	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				x, y, z;
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);

	if ( mark ) {
		mark->loc = Vector3<float>(x, y, z);
		mark->fom = 1;
		mark->sel = 1;
	}
	
	return 0;
}

Tcl_Obj*	marker_create(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);

	Bmarker*			mark = NULL;

	if ( mg ) {
		for ( mark = mg->mark; mark; mark = mark->next )
			if ( id < mark->id ) id = mark->id;
		mark = (Bmarker *) add_item((char **) &mg->mark, sizeof(Bmarker));
	}

	if ( mark ) {
		mark->id = ++id;
		mark->img_num = mg->img_num;
		marker_create(mark, objc, objv);
	}

	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_create(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);

	Bmarker*			mark = NULL;

	if ( rec ) {
		for ( mark = rec->mark; mark; mark = mark->next )
			if ( id < mark->id ) id = mark->id;
		mark = (Bmarker *) add_item((char **) &rec->mark, sizeof(Bmarker));
	}

	if ( mark ) {
		mark->id = ++id;
		mark->img_num = 0;
		marker_create(mark, objc, objv);
	}
	
	Tcl_SetIntObj(returnObj, id);

	return returnObj;
}

Tcl_Obj*	marker_center(Bmarker* mark, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	Vector3<float>		center(p->sizeX()/2, p->sizeY()/2, p->sizeZ()/2);

	if ( mark ) markers_center(mark, center);

	return returnObj;
}

Tcl_Obj*	marker_generate_from_seed(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !project ) return returnObj;
	if ( !project->field ) return returnObj;
	if ( !project->field->mg ) return returnObj;

	if ( project->field->mg->origin[0] < 1 || project->field->mg->origin[1] < 1 )
		project_set_nominal_mg_origins(project);

	project_mg_tilt_to_matrix(project);
		
//	if ( !project->rec || !project->rec->mark )
		project_calculate_model(project);
	
	project_generate_markers(project);

	return returnObj;
}

Tcl_Obj*	marker_find(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	double				r(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &r);

	Bmarker*			mark = NULL;
	
	if ( mg ) {
		mark = img_find_gold_particles(p, mg->img_num, r, 0, 0);
		if ( mark ) {
			if ( !mg->mark ) mg->mark = mark;
			else markers_add(&mg->mark, mark, 2*mg->mark_radius, 2);
		}
	}
	
	int					nmark = count_list((char *) mg->mark);

	Tcl_SetIntObj(returnObj, nmark);
	
	return returnObj;
}

Tcl_Obj*	marker_snap(Bmicrograph* mg, Bimage* p, long mid, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	if ( mg ) {
		Bmarker* 		mark = NULL;
		for ( mark = mg->mark; mark && mark->id != mid; mark = mark->next ) ;
		if ( mark ) {
			long			i(0);
			double			v, avg(0), wcoor(0);
			Vector3<float>	acoor;	
			Vector3<long>	size = Vector3<long>(6*mg->mark_radius, 6*mg->mark_radius, 1);
			Vector3<long>	coor = mark->loc - size/2;
			Bimage*			pm = p->extract(mg->img_num, coor, size, 0, 0);
			for ( i=0; i<pm->image_size(); ++i ) avg += (*pm)[i];
			avg /= pm->image_size();
			for ( i=0; i<pm->image_size(); ++i ) {
				v = (*pm)[i];
				if ( (*pm)[i] < avg ) {
					acoor += pm->coordinates(i) * v;
					wcoor += v;
				}
			}
			acoor /= wcoor;
			mark->loc = acoor + coor;
			delete pm;
		}
	}
	
	Tcl_SetDoubleObj(returnObj, 0);
	
	return returnObj;
}

Tcl_Obj*	marker_extract(Bmicrograph* mg, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
		
	Bstring				filename;
	float				box_radius(0);
	Bimage*				part = NULL;
	
	/* Particle image file name */
	if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
	
	if ( filename.length() ) {
		if ( mg ) {
			part = micrograph_extract_gold(mg, p, box_radius);
			if ( part ) {
				mg->fpart = filename;
				write_img(filename, part, 0);
				delete part;
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	marker_extract(Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
		
	Bstring				filename;
	float				box_radius(0);
	Bimage*				part = NULL;
	
	/* Particle image file name */
	if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
	
	if ( filename.length() ) {
		if ( rec ) {
			part = reconstruction_extract_gold(rec, p, box_radius);
			if ( part ) {
				rec->fpart = filename;
				write_img(filename, part, 0);
				delete part;
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	marker_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bfield*				field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec;
	Bmarker*			mark;

	if ( action == "all" ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				kill_list((char *) mg->mark, sizeof(Bmarker));
				mg->mark = NULL;
			}
		}
		for ( rec = project->rec; rec; rec = rec->next ) {
			kill_list((char *) rec->mark, sizeof(Bmarker));
			rec->mark = NULL;
		}
	} else if ( action == "fom" ) {
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( mark = mg->mark; mark; mark = mark->next )
					if ( mark->fom < fom_cut ) mark->sel = 0;
					else if ( mark->sel < 1 ) mark->sel = 1;
				markers_delete_non_selected(&mg->mark);
			}
		}
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( mark = rec->mark; mark; mark = mark->next )
				if ( mark->fom < fom_cut ) mark->sel = 0;
				else if ( mark->sel < 1 ) mark->sel = 1;
			markers_delete_non_selected(&rec->mark);
		}
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
		if ( id > 0 ) {	// Delete all markers with this id for both micrographs and reconstructions
			for ( field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg; mg = mg->next ) {
					for ( mark = mg->mark; mark && mark->id != id; mark = mark->next ) ;
					if ( mark )
						remove_item((char **)&mg->mark, (char *)mark, sizeof(Bmarker));
				}
			}
			for ( rec = project->rec; rec; rec = rec->next ) {
				for ( mark = rec->mark; mark && mark->id != id; mark = mark->next ) ;
				if ( mark )
					remove_item((char **)&rec->mark, (char *)mark, sizeof(Bmarker));
			}
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_delete(Bfield* field, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	if ( !field ) return returnObj;

	int					id(0);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bmicrograph*		mg = NULL;
	Bmarker*			mark;

	if ( action == "all" ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			kill_list((char *) mg->mark, sizeof(Bmarker));
			mg->mark = NULL;
		}
	} else if ( action == "fom" ) {
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( mark = mg->mark; mark; mark = mark->next )
				if ( mark->fom < fom_cut ) mark->sel = 0;
				else if ( mark->sel < 1 ) mark->sel = 1;
			markers_delete_non_selected(&mg->mark);
		}
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
		if ( id > 0 ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( mark = mg->mark; mark && mark->id != id; mark = mark->next ) ;
				if ( mark )
					remove_item((char **)&mg->mark, (char *)mark, sizeof(Bmarker));
			}
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	if ( !mg ) return returnObj;

	int					id(0);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bmarker*			mark;

	if ( action == "all" ) {
		kill_list((char *) mg->mark, sizeof(Bmarker));
		mg->mark = NULL;
	} else if ( action == "fom" ) {
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
		for ( mark = mg->mark; mark; mark = mark->next )
			if ( mark->fom < fom_cut ) mark->sel = 0;
			else if ( mark->sel < 1 ) mark->sel = 1;
//		mark = mg->mark;
//		mg->mark = markers_copy_selected(mark);
//		kill_list((char *) mark, sizeof(Bmarker));
		markers_delete_non_selected(&mg->mark);
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
		if ( id > 0 ) {
			for ( mark = mg->mark; mark && mark->id != id; mark = mark->next ) ;
			if ( mark ) remove_item((char **)&mg->mark, (char *)mark, sizeof(Bmarker));
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	if ( !rec ) return returnObj;

	int					id(0);
	double				fom_cut(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bmarker*			mark;

	if ( action == "all" ) {
		kill_list((char *) rec->mark, sizeof(Bmarker));
		rec->mark = NULL;
	} else if ( action == "fom" ) {
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &fom_cut);
		for ( mark = rec->mark; mark; mark = mark->next )
			if ( mark->fom < fom_cut ) mark->sel = 0;
			else if ( mark->sel < 1 ) mark->sel = 1;
		markers_delete_non_selected(&rec->mark);
	} else {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
		if ( id > 0 ) {
			for ( mark = rec->mark; mark && mark->id != id; mark = mark->next ) ;
			if ( mark ) remove_item((char **)&rec->mark, (char *)mark, sizeof(Bmarker));
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	marker_edge(Bmarker* mark, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					edge(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &edge);

	float				xmin, xmax, ymin, ymax;

	if ( mark ) {
		xmin = edge;
		xmax = p->sizeX() - edge;
		ymin = edge;
		ymax = p->sizeY() - edge;
		for ( ; mark && mark->id; mark = mark->next ) {
			if ( mark->loc[0] < xmin || mark->loc[0] > xmax || mark->loc[1] < ymin || mark->loc[1] > ymax )
				mark->fom = -fabs(mark->fom);
			else
				mark->fom = fabs(mark->fom);
		}
	}
	
	return returnObj;
}

Tcl_Obj*	marker_accept(Bmarker* mark, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0);

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);

	for ( ; mark && mark->id != id; mark = mark->next ) ;

	if ( mark ) {
		mark->loc += mark->err;
		mark->err = 0;
	}
	
	return returnObj;
}



