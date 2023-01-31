/**
@file	tcltk_bfil.cpp
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20150617
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bfil.h"
#include "mg_extract.h"
#include "mg_helix.h"
#include "mg_select.h"
#include "rwmg.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern Bimage* 	imglist;

// Internal function prototypes
Tcl_Obj*	filament_count(Bproject* project);
Tcl_Obj*	filament_count(Bfield* field);
Tcl_Obj*	filament_count(Bmicrograph* mg);
Tcl_Obj*	filament_count(Breconstruction* rec);
Tcl_Obj*	filament_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_extract(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_profile(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_center(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	filament_to_particles(Bproject* project, int objc, Tcl_Obj *CONST objv[]);

Tcl_Obj*	node_count(Bproject* project);
Tcl_Obj*	node_count(Bfield* field);
Tcl_Obj*	node_count(Bmicrograph* mg);
Tcl_Obj*	node_count(Breconstruction* rec);
Tcl_Obj*	node_ids(Bfilament* fillist);
Tcl_Obj*	node_location(Bfilament* fillist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_select(Bfilament* fillist, double radius, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_move(Bfilament* fillist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_create(Bfilament** fillist, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	node_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[]);


/*
@brief 	Functions for filament picking.
	Filament interface syntax:
		Bmg filament [item] [action] [params]
	Actions:
		count							counts the number of filaments
		delete [all|mg|rec|fid]			deletes filaments
		extract [f] [w] [ax] [sf] [p]	extracts filaments from a micrograph
		profile [fid] [nid] [w]			plots a filament profile
		center [w]						centers nodes based on profiles
		toparticles [w] [b] [r] [a]		generates particles along filaments
	Parameters:
		fid			filament identifier
		nid			node identifier
		f			file name
		w			filament width
		ax			extracted filament axis
		sf			split flag
		p			filament path
		b			boxing interval
		r			helical rise
		a			helical rotation angle
*/
Tcl_Obj*	do_filament(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_filament: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	if ( !item.contains("Field") && !item.contains("Micrograph") && !item.contains("Reconstruction") ) {
		cerr << "Error in do_filament: Item " << item << " must be a micrograph or reconstruction!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2, base, delete_what;
	
	if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
		if ( field ) {
			mg = field->mg;
			filename = field->mg->fmg;
		} else
			cerr << "Error in do_filament: Field " << id << " not found!" << endl;
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
		} else
			cerr << "Error in do_filament: Micrograph " << id << " not found!" << endl;
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( rec ) {
			filename = rec->frec;
		} else
			cerr << "Error in do_filament: Reconstruction " << id << " not found!" << endl;
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
		cout << "DEBUG do_filament: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		if ( mg ) returnObj = filament_count(mg);
		else if ( field ) returnObj = filament_count(field);
		else if ( rec ) returnObj = filament_count(rec);
		else returnObj = filament_count(project);
	} else if ( action == "delete" ) {
		if ( objc > 4 ) delete_what = Tcl_GetStringFromObj(objv[4], NULL);
		if ( delete_what == "all" || delete_what == "mg" || delete_what == "rec" )
			returnObj = filament_delete(project, objc, objv);
		else if ( mg ) returnObj = filament_delete(mg, objc, objv);
		else if ( rec ) returnObj = filament_delete(rec, objc, objv);
		else returnObj = filament_delete(project, objc, objv);
	} else if ( action == "extract" ) {
		if ( p ) returnObj = filament_extract(mg, rec, p, objc, objv);
	} else if ( action == "profile" ) {
		if ( p ) returnObj = filament_profile(mg, rec, p, objc, objv);
	} else if ( action == "center" ) {
		if ( p ) returnObj = filament_center(mg, rec, p, objc, objv);
	} else if ( action == "toparticles" ) {
		returnObj = filament_to_particles(project, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	if ( !returnObj ) {
		returnObj = Tcl_NewObj();
		Tcl_SetIntObj(returnObj, 0);
	}

	return returnObj;
}

Tcl_Obj*	filament_count(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nfil(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				if ( mg->fil )
					nfil += count_list((char *) mg->fil);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			if ( rec->fil )
				nfil += count_list((char *) rec->fil);
	}
	
	Tcl_SetIntObj(returnObj, nfil);
	
	return returnObj;
}

Tcl_Obj*	filament_count(Bfield* field)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nfil(0);
	Bmicrograph*		mg;

	if ( field ) {
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mg->fil )
				nfil += count_list((char *) mg->fil);
	}
	
	Tcl_SetIntObj(returnObj, nfil);
	
	return returnObj;
}

Tcl_Obj*	filament_count(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nfil(0);

	if ( mg && mg->fil )
		nfil = count_list((char *) mg->fil);
	
	Tcl_SetIntObj(returnObj, nfil);
	
	return returnObj;
}

Tcl_Obj*	filament_count(Breconstruction* rec)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nfil(0);

	if ( rec && rec->fil )
		nfil = count_list((char *) rec->fil);
	
	Tcl_SetIntObj(returnObj, nfil);
	
	return returnObj;
}

Tcl_Obj*	filament_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), n(-1), all(0);

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	char				string[64] = "";
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec;
	Bfilament*			fil = NULL;
	
	if ( action ==  "all" || action ==  "mg" ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				filament_kill(mg->fil);
				mg->fil = NULL;
			}
		}
		all = 1;
	}

	if ( action ==  "all" || action ==  "rec" ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			filament_kill(rec->fil);
			rec->fil = NULL;
		}
		all = 1;
	}
	
	if ( !all ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
//		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &n);
	}
	
	if ( fid > 0 ) {
		if ( project->select < 1 ) {
			if ( n < 0 ) {
				if ( field ) for ( mg = field->mg; mg; mg = mg->next ) {
					for ( fil = mg->fil; fil && fil->id != fid; fil = fil->next ) ;
					if ( fil ) {
						kill_list((char *)fil->node, sizeof(Bfilnode));
						remove_item((char **)&mg->fil, (char *)fil, sizeof(Bfilament));
					}
				}
			} else {
				if ( field ) for ( mg = field->mg; mg && mg->img_num != n; mg = mg->next ) ;
				if ( mg ) {
					for ( fil = mg->fil; fil && fil->id != fid; fil = fil->next ) ;
					if ( fil ) {
						kill_list((char *)fil->node, sizeof(Bfilnode));
						remove_item((char **)&mg->fil, (char *)fil, sizeof(Bfilament));
					}
				}
			}
		} else {
			for ( rec = project->rec; rec; rec = rec->next ) {
				for ( fil = rec->fil; fil && fil->id != fid; fil = fil->next ) ;
				if ( fil ) {
					kill_list((char *)fil->node, sizeof(Bfilnode));
					remove_item((char **)&rec->fil, (char *)fil, sizeof(Bfilament));
				}
			}
		}
	}	
	
	snprintf(string, 60, "%d", fid);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}

Tcl_Obj*	filament_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0);
	Bfilament*			fil = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	
	if ( mg && fid > 0 ) {
		for ( fil = mg->fil; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil ) {
			kill_list((char *)fil->node, sizeof(Bfilnode));
			remove_item((char **)&mg->fil, (char *)fil, sizeof(Bfilament));
		}
	}
	
	Tcl_SetIntObj(returnObj, fid);
	
	return returnObj;
}

Tcl_Obj*	filament_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0);
	Bfilament*			fil = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	
	if ( rec && fid > 0 ) {
		for ( fil = rec->fil; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil ) {
			kill_list((char *)fil->node, sizeof(Bfilnode));
			remove_item((char **)&rec->fil, (char *)fil, sizeof(Bfilament));
		}
	}
	
	Tcl_SetIntObj(returnObj, fid);
	
	return returnObj;
}

int			filament_renumber(Bfilament* fil)
{
	int					nnode(0);
	Bfilnode*			fnode = NULL;

	for ( fnode = fil->node; fnode; fnode = fnode->next ) fnode->id = ++nnode;

	return nnode;
}

Tcl_Obj*	filament_extract(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
	
	Bstring				filename, filpath;
	int					nfil(0), split(0), axis(0);
	double				filament_width(10);
	
	if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &filament_width);
	if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &axis);	/* Axis selection: x=1, y=2, z=3 */
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &split);	/* Flag to split particles into individual images */
	if ( objc > 8 ) filpath = Tcl_GetStringFromObj(objv[8], NULL);
	
//	cout << filename << "\t" << filament_width << "\t" << axis << "\t" << split << endl;
	
	if ( filename.length() ) {
		if ( mg ) {
			if ( split ) {
				filament_setup_filenames(mg->fil, filename, filpath);
				mg->ffil = 0;
			} else {
				if ( filpath.length() ) {
					if ( filename.contains("/") ) filename = filename.post_rev('/');
					filename = filpath + "/" + filename;
				}
				mg->ffil = filename;
			}
			nfil = micrograph_extract_filaments(mg, p, filament_width, axis);
		} else if ( rec ) {
			if ( split ) {
				filament_setup_filenames(rec->fil, filename, filpath);
				rec->ffil = 0;
			} else {
				if ( filpath.length() ) {
					if ( filename.contains("/") ) filename = filename.post_rev('/');
					filename = filpath + "/" + filename;
				}
				rec->ffil = filename;
			}
			nfil = reconstruction_extract_filaments(rec, p, filament_width, axis);
		}
	}
	
	Tcl_SetIntObj(returnObj, nfil);

	return returnObj;
}

Tcl_Obj*	filament_profile(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
	
	int					fid(0), nid(0);
	double				filament_width(10);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &nid);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &filament_width);
	
	long				i, n(0), img_num(0);
	double*				profile = NULL;
	char				string[MAXLINELEN] = "";
	Bfilament* 			fil = NULL;
	
	if ( mg ) {
		fil = mg->fil;
		img_num = mg->img_num;
	} else if ( rec ) fil = rec->fil;
	
	if ( fil ) {
		for ( ; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil ) {
			profile = filament_profile(fil->node, p, img_num, nid, filament_width, n);
			snprintf(string, MAXLINELEN, "%ld", n);
			Tcl_SetStringObj(returnObj, string, strlen(string));
			if ( profile ) {
				for ( i=0; i<2*n; i++ ) {
					snprintf(string, MAXLINELEN, " %lf", profile[i]);
					Tcl_AppendToObj(returnObj, string, strlen(string));
				}
				delete[] profile;
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	filament_center(Bmicrograph* mg, Breconstruction* rec, Bimage* p, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !p ) return returnObj;
	
	int					filament_width(40);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &filament_width);
	
	long				n(0), img_num(0);
	double				cc(0);
	Bfilament* 			fil = NULL;
	
	if ( mg ) {
		fil = mg->fil;
		img_num = mg->img_num;
	} else if ( rec ) fil = rec->fil;
	
	if ( fil ) {
		cc += filaments_center(fil, p, img_num, filament_width);
		n++;
	}
	
	if ( n ) cc /= n;

	Tcl_SetDoubleObj(returnObj, cc);
	
	return returnObj;
}

Tcl_Obj*	filament_to_particles(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nfil(0);
	double				filament_width(10), boxing_interval(0), rise(1e30), angle(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &filament_width);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &boxing_interval);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &rise);		/* Helical rise */
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &angle);	/* Helical angle */
	
	Vector3<int>		box_size((int)filament_width, (int)filament_width, (int)filament_width);
	nfil = project_filaments_to_particles(project, box_size, boxing_interval, rise, angle);
	
	Tcl_SetIntObj(returnObj, nfil);

	return returnObj;
}

/*
	Filament node interface syntax:
		Bmg node [item] [action] [params]
	Actions:
		count						counts the number of filament nodes
		ids							returns the identifiers of filaments and nodes
		location [fid] [nid]		returns the location of a node
		select [x] [y] [z]			returns the filament and node identifiers at the given location
		move [fid] [nid] [dx] [dy] [dz]	updates a node location by the given vector
		create [x] [y] [z] [fid] [nid]	creates a new node at the given location and after the given node
		delete [all|mg|rec|fid] [nid]	deletes nodes
	Parameters:
		fid			filament identifier
		nid			node identifier
		x			x location
		y			y location
		z			z location
		dx			change in x
		dy			change in y
		dz			change in z
*/
Tcl_Obj*	do_node(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_node: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');

	if ( !item.contains("Field") && !item.contains("Micrograph") && !item.contains("Reconstruction") ) {
		cerr << "Error in do_node: Item " << item << " must be a micrograph or reconstruction!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	
	double				radius(10);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bfilament*			fil = NULL;
	Bstring				delete_what;
	
	if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next ) ;
			if ( mg ) break;
		}
		if ( mg ) {
			fil = mg->fil;
			radius = mg->fil_node_radius;
		} else cerr << "Error in do_node: Micrograph " << id << " not found!" << endl;
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( rec ) {
			fil = rec->fil;
			radius = rec->fil_node_radius;
		} else cerr << "Error in do_node: Reconstruction " << id << " not found!" << endl;
	}

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_node: Action: " << action << " (" << action.length() << ")" << endl;

	if ( action == "count" ) {
		if ( mg ) returnObj = node_count(mg);
		else if ( field ) returnObj = node_count(field);
		else if ( rec ) returnObj = node_count(rec);
		else returnObj = node_count(project);
	} else if ( action == "ids" ) {
		if ( fil ) returnObj = node_ids(fil);
	} else if ( action == "location" ) {
		if ( fil ) returnObj = node_location(fil, objc, objv);
	} else if ( action == "select" ) {
		if ( fil ) returnObj = node_select(fil, radius, objc, objv);
	} else if ( action == "move" ) {
		if ( fil ) returnObj = node_move(fil, objc, objv);
	} else if ( action == "create" ) {
		if ( mg ) returnObj = node_create(&mg->fil, objc, objv);
		else if ( rec ) returnObj = node_create(&rec->fil, objc, objv);
	} else if ( action == "delete" ) {
		if ( objc > 4 ) delete_what = Tcl_GetStringFromObj(objv[4], NULL);
		if ( delete_what == "all" || delete_what == "mg" || delete_what == "rec" )
			returnObj = node_delete(project, objc, objv);
		else if ( mg ) returnObj = node_delete(mg, objc, objv);
		else if ( rec ) returnObj = node_delete(rec, objc, objv);
		else returnObj = node_delete(project, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	if ( !returnObj ) {
		returnObj = Tcl_NewObj();	// Return object should be empty in some cases
//		Tcl_SetIntObj(returnObj, 0);
	}

	return returnObj;
}

Tcl_Obj*	node_count(Bproject* project)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nnode(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec;
	Bfilament*			fil = NULL;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( fil = mg->fil; fil; fil = fil->next )
					if ( fil->node ) nnode += count_list((char *) fil->node);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( fil = rec->fil; fil; fil = fil->next )
				if ( fil->node ) nnode += count_list((char *) fil->node);
	}
	
	Tcl_SetIntObj(returnObj, nnode);
	
	return returnObj;
}

Tcl_Obj*	node_count(Bfield* field)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nnode(0);
	Bmicrograph*		mg = NULL;
	Bfilament*			fil = NULL;

	if ( field )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( fil = mg->fil; fil; fil = fil->next )
				if ( fil->node ) nnode += count_list((char *) fil->node);
	
	Tcl_SetIntObj(returnObj, nnode);
	
	return returnObj;
}

Tcl_Obj*	node_count(Bmicrograph* mg)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nnode(0);
	Bfilament*			fil = NULL;

	if ( mg && mg->fil )
		for ( fil = mg->fil; fil; fil = fil->next )
			if ( fil->node ) nnode += count_list((char *) fil->node);
	
	Tcl_SetIntObj(returnObj, nnode);
	
	return returnObj;
}

Tcl_Obj*	node_count(Breconstruction* rec)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nnode(0);
	Bfilament*			fil = NULL;

	if ( rec && rec->fil )
		for ( fil = rec->fil; fil; fil = fil->next )
			if ( fil->node ) nnode += count_list((char *) fil->node);
	
	Tcl_SetIntObj(returnObj, nnode);
	
	return returnObj;
}

Tcl_Obj*	node_ids(Bfilament* fillist)
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

	if ( fillist ) {
		for ( fil = fillist; fil; fil = fil->next ) {
			for ( fnode = fil->node; fnode; fnode = fnode->next ) {
				snprintf(string, MAXLINELEN, " %d %d", fil->id, fnode->id);
				Tcl_AppendToObj(returnObj, string, strlen(string));
			}
		}
	}
	
	return returnObj;
}

Tcl_Obj*	node_location(Bfilament* fillist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), nid(0);	
	char				string[MAXLINELEN] = "";
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &nid);

	if ( fillist ) {
		for ( fil = fillist; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil )
			for ( fnode = fil->node; fnode && fnode->id != nid; fnode = fnode->next ) ;
		if ( fnode ) {
			snprintf(string, MAXLINELEN, "%f %f %f", fnode->loc[0], fnode->loc[1], fnode->loc[2]);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	node_select(Bfilament* fillist, double radius, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), nid(0);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	
	double				d, dmin;
	char				string[64] = "";
	Vector3<float>		loc(x, y, z);
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

//	cout << "Selection coordinates: " << loc << endl;
	if ( fillist ) {
		dmin = 2*radius;
		for ( fil = fillist; fil; fil = fil->next ) {
			for ( fnode = fil->node; fnode; fnode = fnode->next ) {
				d = loc.distance(fnode->loc);
				if ( dmin > d ) {
					dmin = d;
					fid = fil->id;
					nid = fnode->id;
				}
			}
		}
		if ( dmin <= radius ) {
			snprintf(string, 60, "%d %d", fid, nid);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	node_move(Bfilament* fillist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), nid(0);
	double				dx(0), dy(0), dz(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &nid);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &dx);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &dy);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &dz);

	Vector3<float>		d = Vector3<float>(dx, dy, dz);
	char				string[64] = "";
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

	if ( fillist ) {
		for ( fil = fillist; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil )
			for ( fnode = fil->node; fnode && fnode->id != nid; fnode = fnode->next ) ;
		if ( fnode ) {
			fnode->loc += d;
			snprintf(string, 60, "%d %d", fid, nid);
			Tcl_SetStringObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	node_create(Bfilament** fillist, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	int					fid(0), nid(0), fid_max(0), nid_max(0);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &fid);	// Filament
	if ( objc > 8 ) Tcl_GetIntFromObj(NULL, objv[8], &nid);	// Node

	char				string[64] = "";
	Vector3<float>		loc(x, y, z);
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;
	Bfilnode*			fnode2 = NULL;

	if ( !(*fillist) ) {
		fid = nid = 1;
		fil = filament_add(fillist, fid);
	}
	
	for ( fil = *fillist; fil && fil->id != fid; fil = fil->next )
		if ( fid_max < fil->id ) fid_max = fil->id;
	
	if ( !fil ) {
		fil = filament_add(fillist, ++fid_max);
		nid = 1;
	}
	
	for ( fnode = fil->node; fnode && fnode->id != nid; fnode = fnode->next )
		if ( nid_max < nid ) nid_max = nid;
	
	if ( !fnode ) {
		nid = ++nid_max;
		fnode = filament_node_add(&fil->node, nid);
	} else {
		fnode2 = fnode->next;
		fnode->next = NULL;
		fnode = filament_node_add(&fnode, ++nid);
		fnode->next = fnode2;
		for ( ; fnode2; fnode2 = fnode2->next ) fnode2->id = ++nid;
	}

	fnode->loc = loc;

	snprintf(string, 60, "%d %d", fil->id, fnode->id);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}

Tcl_Obj*	node_delete(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	Bstring				action = Tcl_GetStringFromObj(objv[4], NULL);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec;

	if ( action == "all" || action == "mg" ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				filament_kill(mg->fil);
				mg->fil = NULL;
			}
		}
	}

	if ( action == "all" || action == "rec" ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			filament_kill(rec->fil);
			rec->fil = NULL;
		}
	}
	
//	cout << "Deleting all" << endl;
	
	return returnObj;
}

Tcl_Obj*	node_delete(Bmicrograph* mg, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), nid(0);
	char				string[64] = "";
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &nid);
	
	if ( mg ) {
		for ( fil = mg->fil; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil ) {
			for ( fnode = fil->node; fnode && fnode->id != nid; fnode = fnode->next ) ;
			if ( fnode )
				remove_item((char **)&fil->node, (char *)fnode, sizeof(Bfilnode));
			filament_renumber(fil);
			if ( !fil->node )
				remove_item((char **)&mg->fil, (char *)fil, sizeof(Bfilament));
		}
	}
	
	snprintf(string, 60, "%d %d", fid, nid);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}


Tcl_Obj*	node_delete(Breconstruction* rec, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					fid(0), nid(0);
	char				string[64] = "";
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &fid);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &nid);
	
	if ( rec ) {
		for ( fil = rec->fil; fil && fil->id != fid; fil = fil->next ) ;
		if ( fil ) {
			for ( fnode = fil->node; fnode && fnode->id != nid; fnode = fnode->next ) ;
			if ( fnode )
				remove_item((char **)&fil->node, (char *)fnode, sizeof(Bfilnode));
			filament_renumber(fil);
			if ( !fil->node )
				remove_item((char **)&rec->fil, (char *)fil, sizeof(Bfilament));
		}
	}
	
	snprintf(string, 60, "%d %d", fid, nid);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}



