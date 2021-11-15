/**
@file	tcltk_bmg.cpp
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author Bernard Heymann
@date	Created: 20030813
@date	Modified: 20210311
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

//#include "tcltk_bmg.h"
#include "tcltk_bbox.h"
#include "tcltk_bfil.h"
#include "tcltk_bhelix.h"
#include "tcltk_bxtal.h"
#include "tcltk_bmarker.h"
#include "mg_processing.h"
#include "mg_multiple.h"
#include "mg_img_proc.h"
#include "mg_select.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "mg_tomography.h"
#include "mg_tags.h"
#include "rwmg.h"
#include "rwimg.h"
#include "rwmodel_param.h"
#include "scatter.h"
#include "linked_list.h"
#include "timer.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 			verbose;		// Level of output to the screen
extern Bproject*	project;
extern Bimage* 		imglist;

// Internal function prototypes
Tcl_Obj*	do_get(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
int			do_set(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_rps(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_emfp(int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_material_list();
Tcl_Obj*	do_ctf_fit(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_mg_sort(Bproject* project, int objc, Tcl_Obj *CONST objv[]);

/**
@brief 	Implements the "Bmg" command in Tcl/Tk to access micrograph parameter files through Bsoft.
@param 	*interp		a Tcl interpreter within Tcl.
@param 	objc		number of arguments passed (+1).
@param 	*objv[]		arguments passed as Tcl objects.
@return int			Tcl result.

	Bmg command syntax:
		Bmg <action> <arguments>.
		where:
			action			"create", "exists", "read", "add", "write", "kill", "get", "set", 
							"image_type", "rps", "ctf_fit", "update_matrices",
							"unitcell_vectors"	"refine" "box" "filament" "node" "layerline" "spot" "marker"
			arguments				action-specific arguments:
				"create"			<number_micrographs> <number_reconstructions>
				"read"				<filename>
				"add"				<filename>
				"write"				<filename>
				"get"				<property> [arguments]
				"set"				<property> <value>
				"rps"
				"ctf_fit"			<level> <lores> <hires>
				"sort"				<tag>
				"update_matrices"
				"unitcell_vectors"
				"findaxis"			<axis> <step> <start> <end>
				"track"				<iterations> <refine_markers>
				"refine"			<operation>
				"box"				[arguments]
				"filament"			[arguments]
				"node"				[arguments]
				"layerline"			[arguments]
				"spot"				[arguments]
				"marker"			[arguments]
				where:
					property	"active <flag>"
								"id <string>"
								"id_from_index <mg_num>"
								"img_num <mg_num>"
								"number_of_mg"
								"number_of_rec"
								"number_of_part"
								"field <string>"
								"filename <string> <imgtype>"
								"select <y/n>"
								"fom <mg_num>"
								"pixel_size <angstrom>"
								"size <x> <y> <z>"
								"origin <x> <y> <z> / <mg_num>"
								"part_origin <x> <y> <z> / <part_num>"
								"scale <x> <y> <z> / <mg_num>"
								"dose" <electrons/Å2>
								"intensity"
								"defocus <angstrom>"
								"defocus_deviation <angstrom>"
								"astigmatism_angle <radians>"
								"volt <volts>"
								"Cs <angstrom>"
								"amp_fac <fraction>"
								"focal_length <angstrom>"
								"aperture <angstrom>"
								"slit_width <volts>"
								"zero <angstrom>"
								"baseline <string>"
								"envelope <string>"
								"view <mg_num>"
								"axis [<radians>] <mg_num>"
								"tilt [<radians>] <mg_num>"
								"level [<radians>] <mg_num>"
								"unitcell <ux> <uy> <vx> <vy>"
								"helix_axis <radians>"
								"helix_rise <angstrom>"
								"helix_angle <radians>"
								"helix_radius <angstrom>"
								"box_size <x> <y> <z>"
								"bad_radius <pixels>"
								"filament_width <width>"
								"filament_node_radius"
								"marker_radius <pixels>"
	Return values:
		Each action may have a return value:
			"create"	(none)
			"exists"	0=no, 1=yes
			"read"		micrograph id
			"write"		(none)
            "kill"		(none)
			"get"		return value based on property
			"set"		modify micrograph property

**/
int 		project_processing(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
//	verbose = 967;

	Tcl_ResetResult(interp);
	
	if ( objc < 2 ) {
		Tcl_AppendResult(interp, "wrong # args", (char *)NULL);
		return TCL_ERROR;
	}

	Bstring				action = Tcl_GetStringFromObj(objv[1], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG BmgCmd: Action: " << action << " (" << action.length() << ")" << endl;
	
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	char				str[128] = " ";
	Tcl_SetStringObj(returnObj, str, 1);
	
//	int					nmg(1), nrec(0);
//	int					type(1);
	Bstring				imgtype("mg");
	Bstring				filename;
	Bstring*			files = NULL;
    Bfield*         	field = NULL;
	Bmicrograph* 		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;

	time_t				ti = time(NULL);
	
//	if ( project ) cout << "Flag active = %d\n", project->select);
	
	if ( action == "create" ) {
		if ( objc > 2 ) filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( objc > 3 ) imgtype = Tcl_GetStringFromObj(objv[3], NULL);
		if ( imglist ) {
			project = project_create_from_image(imglist, imgtype);
		} else {
//			string_add(&files, filename);
//			project = project_create_from_images(files, imgtype);
			project = project_create_from_images(&filename, imgtype);
//			string_kill(files);
		}
	} else if ( action == "exists" ) {
		if ( project ) Tcl_SetIntObj(returnObj, 1);
		else Tcl_SetIntObj(returnObj, 0);
	} else if ( action == "read" ) {
		if ( objc < 2 ) {
			Tcl_AppendResult(interp, "No file name given!", (char *)NULL);
			return TCL_ERROR;
		}
		filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG BmgCmd: File name: " << filename << " (" << filename.length() << ")" << endl;
		if ( project ) project_kill(project);
		project = read_project(filename);
		if ( project->field ) {
			mg = project->field->mg;	// Only use the first micrograph
			Tcl_SetStringObj(returnObj, mg->id.c_str(), mg->id.length());
		}
//		cout << "nrec=" << project_count_reconstructions(project) << endl;
	} else if ( action == "add" ) {
		if ( objc < 2 ) {
			Tcl_AppendResult(interp, "No file name given!", (char *)NULL);
			return TCL_ERROR;
		}
		filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG BmgCmd: File name: " << filename << " (" << filename.length() << ")" << endl;
		files = filename.split();
		project_multi_add(project, files, 0);
		string_kill(files);
		if ( project->field ) {
			mg = project->field->mg;	// Only use the first micrograph
			Tcl_SetStringObj(returnObj, (char *)mg->id.c_str(), mg->id.length());
		}
	} else if ( action == "write" ) {
		if ( objc < 2 ) {
			Tcl_AppendResult(interp, "No file name given!", (char *)NULL);
			return TCL_ERROR;
		}
//		cout << "nrec=" << project_count_reconstructions(project) << endl;
		filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG BmgCmd: File name: " << filename << " (" << filename.length() << ")" << endl;
		if ( !project ) {
			Tcl_AppendResult(interp, "No project in memory!", (char *)NULL);
			return TCL_ERROR;
		}
		project->comment += "\n# Written from Bshow\n# " + Bstring(asctime(localtime(&ti))) + "\n";
		write_project(filename, project);
	} else if ( action == "kill" ) {
		project_kill(project);
		project = NULL;
	} else if ( action == "ids" ) {
        for ( field = project->field; field; field = field->next ) {    
			snprintf(str, 128, "Field-of-view %s %ld ", field->id.c_str(), field->select);
			Tcl_AppendToObj(returnObj, str, strlen(str));
            for ( mg = field->mg; mg; mg = mg->next ) {
				snprintf(str, 128, "Micrograph %s %ld ", mg->id.c_str(), mg->select);
				Tcl_AppendToObj(returnObj, str, strlen(str));
				for ( part = mg->part; part; part = part->next ) {
					snprintf(str, 128, "Particle %d %ld ", part->id, part->sel);
					Tcl_AppendToObj(returnObj, str, strlen(str));
				}
			}
		}
		for ( rec = project->rec; rec; rec = rec->next ) {
			snprintf(str, 128, "Reconstruction %s %ld ", rec->id.c_str(), rec->select);
			Tcl_AppendToObj(returnObj, str, strlen(str));
			for ( part = rec->part; part; part = part->next ) {
				snprintf(str, 128, "Particle %d %ld ", part->id, part->sel);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		}
//	} else if ( action == "item" ) {
//		returnObj = do_item(project, objc, objv);
	} else if ( action == "get" ) {
		returnObj = do_get(project, objc, objv);
	} else if ( action == "set" ) {
		do_set(project, objc, objv);
	} else if ( action == "image_type" ) {
		imgtype = Tcl_GetStringFromObj(objv[2], NULL);
		imgtype = image_type(imgtype);
		Tcl_SetStringObj(returnObj, (char *)imgtype.c_str(), imgtype.length());
	} else if ( action == "rps" ) {
		returnObj = do_rps(project, objc, objv);
	} else if ( action == "ctf_fit" ) {
		returnObj = do_ctf_fit(project, objc, objv);
	} else if ( action == "sort" ) {
		returnObj = do_mg_sort(project, objc, objv);
	} else if ( action == "update_matrices" ) {
		project_set_nominal_mg_origins(project);
		project_mg_tilt_to_matrix(project);
	} else if ( action == "thickness" ) {
		returnObj = do_tomo_thickness(project, objc, objv);
	} else if ( action == "emfp" ) {
		returnObj = do_emfp(objc, objv);
	} else if ( action == "material" ) {
		returnObj = do_material_list();
	} else if ( action == "transfer" ) {
		do_tomo_transfer_seed(project, objc, objv);
	} else if ( action == "findaxis" ) {
		do_tomo_findaxis(project, objc, objv);
	} else if ( action == "track" ) {
		do_tomo_track(project, objc, objv);
	} else if ( action == "refine" ) {
		do_tomo_refine(project, objc, objv);
	} else if ( action == "refine_one" ) {
		do_tomo_refine_one(project, objc, objv);
	} else if ( action == "box" ) {
		returnObj = do_box(project, objc, objv);
	} else if ( action == "filament" ) {
		returnObj = do_filament(project, objc, objv);
	} else if ( action == "node" ) {
		returnObj = do_node(project, objc, objv);
	} else if ( action == "layerline" ) {
		returnObj = do_layerline(project, objc, objv);
	} else if ( action == "spot" ) {
		returnObj = do_spot(project, objc, objv);
	} else if ( action == "marker" ) {
		returnObj = do_marker(project, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}
	
	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}


Bstring		get_all_file_names(Bproject* project)
{
	Bstring				bstr;
	Bstring*			strlist = NULL;
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	Breconstruction*	rec = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Getting all file names" << endl;

	for ( field = project->field; field; field = field->next ) {
		mg = field->mg;
		if ( mg && mg->next && mg->fmg.length() && mg->fmg == mg->next->fmg ) {
			bstr = mg->fmg + ":mg:Field:" + field->id + " ";
			string_add(&strlist, bstr);
		} else {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->fmg.length() ) {
					bstr = mg->fmg + ":mg:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
				if ( mg->fframe.length() ) {
					bstr = mg->fframe + ":frame:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			}
		}
		if ( mg && mg->next && mg->fps.length() && mg->fps == mg->next->fps ) {
			bstr = mg->fps + ":ps:Field:" + field->id + " ";
			string_add(&strlist, bstr);
		} else {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->fps.length() ) {
					bstr = mg->fps + ":ps:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			}
		}
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->fpart.length() ) {
				bstr = mg->fpart + ":part:Micrograph:" + mg->id + " ";
				string_add(&strlist, bstr);
			} else for ( part = mg->part; part; part = part->next ) {
				if ( part->fpart.length() ) {
					bstr = part->fpart + ":part:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			}
			if ( mg->ffil.length() ) {
				bstr = mg->ffil + ":fil:Micrograph:" + mg->id + " ";
				string_add(&strlist, bstr);
			}
			if ( mg->fft.length() ) {
				bstr = mg->fft + ":ft:Micrograph:" + mg->id + " ";
				string_add(&strlist, bstr);
			}
		}
	}
	
	for ( rec = project->rec; rec ; rec = rec->next ) {
		if ( rec->frec.length() ) {
			bstr = rec->frec + ":rec:Reconstruction:" + rec->id + " ";
			string_add(&strlist, bstr);
		}
		if ( rec->fpart.length() ) {
			bstr = rec->fpart + ":part:Reconstruction:" + rec->id + " ";
			string_add(&strlist, bstr);
		} else for ( part = rec->part; part; part = part->next ) {
			if ( part->fpart.length() ) {
				bstr = part->fpart + ":part:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		}
		if ( rec->ffil.length() ) {
			bstr = rec->ffil + ":fil:Reconstruction:" + rec->id + " ";
			string_add(&strlist, bstr);
		}
		if ( rec->fps.length() ) {
			bstr = rec->fps + ":ps:Reconstruction:" + rec->id + " ";
			string_add(&strlist, bstr);
		}
		if ( rec->fft.length() ) {
			bstr = rec->fft + ":ft:Reconstruction:" + rec->id + " ";
			string_add(&strlist, bstr);
		}
	}

	for ( part = project->class_avg; part; part = part->next ) {
		if ( part->fpart.length() && part->fpart != bstr ) {
			bstr = part->fpart + ":part:Class:1 ";
			string_add(&strlist, bstr);
			bstr = part->fpart;
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Last file name: " << bstr << endl;
	
	bstr = string_catenate(strlist);
	
	string_kill(strlist);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Catenation done" << endl;
	
	return bstr;
}

Bstring		get_all_file_names(Bproject* project, Bstring imgtype)
{
	Bstring				bstr;
	Bstring*			strlist = NULL;
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	Breconstruction*	rec = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Getting all file names" << endl;

	for ( field = project->field; field; field = field->next ) {
		mg = field->mg;
		if ( mg && mg->next && mg->fmg.length() && mg->fmg == mg->next->fmg ) {
			if ( imgtype == "mg" ) {
				bstr = mg->fmg + ":mg:Field:" + field->id + " ";
				string_add(&strlist, bstr);
			}
		} else {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->fmg.length() ) {
					if ( imgtype == "mg" ) {
						bstr = mg->fmg + ":mg:Micrograph:" + mg->id + " ";
						string_add(&strlist, bstr);
					}
				}
				if ( mg->fframe.length() ) {
					if ( imgtype == "frame" ) {
						bstr = mg->fframe + ":frame:Micrograph:" + mg->id + " ";
						string_add(&strlist, bstr);
					}
				}
			}
		}
		if ( mg && mg->next && mg->fps.length() && mg->fps == mg->next->fps ) {
			if ( imgtype == "ps" ) {
				bstr = mg->fps + ":ps:Field:" + field->id + " ";
				string_add(&strlist, bstr);
			}
		} else {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->fps.length() ) {
					if ( imgtype == "ps" ) {
						bstr = mg->fps + ":ps:Micrograph:" + mg->id + " ";
						string_add(&strlist, bstr);
					}
				}
			}
		}
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->fpart.length() ) {
				if ( imgtype == "part" ) {
					bstr = mg->fpart + ":part:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			} else for ( part = mg->part; part; part = part->next ) {
				if ( part->fpart.length() ) {
					if ( imgtype == "part" ) {
						bstr = part->fpart + ":part:Micrograph:" + mg->id + " ";
						string_add(&strlist, bstr);
					}
				}
			}
			if ( mg->ffil.length() ) {
				if ( imgtype == "fil" ) {
					bstr = mg->ffil + ":fil:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			}
			if ( mg->fft.length() ) {
				if ( imgtype == "ft" ) {
					bstr = mg->fft + ":ft:Micrograph:" + mg->id + " ";
					string_add(&strlist, bstr);
				}
			}
		}
	}
	
	for ( rec = project->rec; rec ; rec = rec->next ) {
		if ( rec->frec.length() ) {
			if ( imgtype == "rec" ) {
				bstr = rec->frec + ":rec:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		}
		if ( rec->fpart.length() ) {
			if ( imgtype == "part" ) {
				bstr = rec->fpart + ":part:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		} else for ( part = rec->part; part; part = part->next ) {
			if ( part->fpart.length() ) {
				if ( imgtype == "part" ) {
					bstr = part->fpart + ":part:Reconstruction:" + rec->id + " ";
					string_add(&strlist, bstr);
				}
			}
		}
		if ( rec->ffil.length() ) {
			if ( imgtype == "fil" ) {
				bstr = rec->ffil + ":fil:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		}
		if ( rec->fps.length() ) {
			if ( imgtype == "ps" ) {
				bstr = rec->fps + ":ps:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		}
		if ( rec->fft.length() ) {
			if ( imgtype == "ft" ) {
				bstr = rec->fft + ":ft:Reconstruction:" + rec->id + " ";
				string_add(&strlist, bstr);
			}
		}
	}

	for ( part = project->class_avg; part; part = part->next ) {
		if ( part->fpart.length() && part->fpart != bstr ) {
			if ( imgtype == "part" ) {
				bstr = part->fpart + ":part:Class:1 ";
				string_add(&strlist, bstr);
			}
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Last file name: " << bstr << endl;
	
	bstr = string_catenate(strlist);
	
	string_kill(strlist);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_all_file_names: Catenation done" << endl;
	
	return bstr;
}

Tcl_Obj*	part_select(Bparticle* part, int objc, Tcl_Obj *CONST objv[])
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

Tcl_Obj*	part_fom(Bparticle* part, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					id(0), fom_index(0);
	char				str[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &fom_index);

	if ( id > 0 ) {
		for ( ; part && part->id != id; part = part->next ) ;
		if ( part ) {
			snprintf(str, MAXLINELEN, "%f", part->fom[fom_index]);
			Tcl_SetStringObj(returnObj, str, strlen(str));
		}
	}
	
//	cout << "fom for " << id << ": " << string << endl;
	
	return returnObj;
}


Tcl_Obj*	do_get(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	int					i, n(-1), mgi(0), nmg(0), nrec(0), npart(0), err(0);
	char				str[MAXLINELEN] = "";
	Bstring				bstr, imgtype, filename;
//	Vector3<double>		sam(1,1,1);
	Vector3<int>		size;
	View				v;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	Breconstruction*	rec = NULL;
	
	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);

	Bstring				property = Tcl_GetStringFromObj(objv[3], NULL);

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG do_get: Item: " << item << " (" << item.length() << ")" << endl;
		cout << "DEBUG do_get: ID: " << id << " (" << id.length() << ")" << endl;
		cout << "DEBUG do_get: Property: " << property << " (" << property.length() << ")" << endl;
	}
	
//	verbose = 512;
//	cout << "id = \"" << id << "\"" << endl;
	
	if ( item.contains("all" ) ) {
		field = project->field;
		if ( field ) mg = field->mg;
		rec = project->rec;
		nmg = project_count_micrographs(project);
		nrec = project_count_reconstructions(project);
	} else if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
//		if ( !field ) field = project->field;
		if ( !field ) err++;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next) mgi++;
			if ( mg ) break;
		}
		if ( !mg ) err++;
	} else if ( item.contains("Reconstruction") ) {
		project->select = 1;
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( !rec ) err++;
	} else if ( item.contains("Class") ) {
		part = project->class_avg;
	} else {
		cerr << "Error in do_get: item \"" << item << "\" not supported! (property: " << property << ")" << endl;
		return returnObj;
	}

	if ( err ) cerr << "Error in do_get: " << item << " not found! (property: " << property << ")" << endl;


	if ( !project || property.empty() ) {
		Tcl_SetIntObj(returnObj, 0);
	} else if ( property == "active" ) {
		Tcl_SetIntObj(returnObj, project->select);
	} else if ( property == "id" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &n);
		if ( project->select && rec ) {
			Tcl_SetStringObj(returnObj, (char *)rec->id.c_str(), rec->id.length());
		} else if ( mg ) {
			if ( !mg ) {
				cerr << "Error: No micrograph specified!" << endl;
				mg = field->mg;
			}
			if ( mg )
				Tcl_SetStringObj(returnObj, (char *)mg->id.c_str(), mg->id.length());
		} else if ( field ) {
			if ( n >= 0 ) {
				for ( i=0, mg = field->mg; mg && i != n; mg = mg->next, i++ ) ;
				if ( !mg ) {
					cerr << "Error: No micrograph specified!" << endl;
					mg = field->mg;
				}
				Tcl_SetStringObj(returnObj, (char *)mg->id.c_str(), mg->id.length());
			} else {
				Tcl_SetStringObj(returnObj, (char *)field->id.c_str(), field->id.length());
			}
		}
	} else if ( property == "id_from_index" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &n);
		if ( project->select ) {
			for ( i=1, rec = project->rec; rec; rec = rec->next, ++i )
				if ( i == n ) break;
			Tcl_SetStringObj(returnObj, (char *)rec->id.c_str(), rec->id.length());
		} else {
			for ( i=1, field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg; mg = mg->next, ++i )
					if ( i == n) break;
				if ( mg ) break;
			}
			Tcl_SetStringObj(returnObj, (char *)mg->id.c_str(), mg->id.length());
		}
	} else if ( property == "img_num" ) {
		if ( mg ) n = mg->img_num;
		else n = 0;
		Tcl_SetIntObj(returnObj, n);
	} else if ( property == "field" ) {
		if ( !field ) field = project->field;
		if ( field )
			Tcl_SetStringObj(returnObj, (char *)field->id.c_str(), field->id.length());
	} else if ( property == "number_of_item" ) {
		if ( project->select ) {
			if ( part ) n = project_count_rec_particles(project);
			else n = nrec;
			if ( !n ) n = project_count_reconstructions(project);
		} else if ( field ) {
			if ( part ) n = project_count_mg_particles(project);
			else n = nmg;
			if ( !n ) for ( n=0, mg = field->mg; mg; mg = mg->next) n++;
		}
		Tcl_SetIntObj(returnObj, n);
//		cout << "Number of items = " << n << endl;
	} else if ( property == "number_of_mg" ) {
		if ( !nmg && field ) for ( nmg=0, mg = field->mg; mg; mg = mg->next) nmg++;
		Tcl_SetIntObj(returnObj, nmg);
//		cout << "Number of micrographs = " << nmg << endl;
	} else if ( property == "number_of_rec" ) {
		for ( nrec=0, rec = project->rec; rec; rec = rec->next ) nrec++;
		Tcl_SetIntObj(returnObj, nrec);
	} else if ( property == "number_of_part" ) {
		if ( !part ) {
			if ( mg ) part = mg->part;
			if ( rec ) part = rec->part;
		}
		for ( npart=0; part; part = part->next ) npart++;
		Tcl_SetIntObj(returnObj, npart);
	} else if ( property == "index_of_item" ) {
		if ( rec ) {
			for ( i=0, rec = project->rec; rec && rec->id != id; rec = rec->next ) i++;
		} else {
			for ( i=0, field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg && mg->id != id; mg = mg->next) i++;
				if ( mg ) break;
			}
		}
		Tcl_SetIntObj(returnObj, ++i);
	} else if ( property == "image_filenames" ) {
		if ( objc > 4 ) {
			imgtype = Tcl_GetStringFromObj(objv[4], NULL);
			bstr = get_all_file_names(project, imgtype);
		} else {
			bstr = get_all_file_names(project);
		}
		Tcl_AppendToObj(returnObj, bstr.c_str(), bstr.length());
	} else if ( property == "filename" ) {
		imgtype = Tcl_GetStringFromObj(objv[4], NULL);
		if ( mg ) {
//			cout << "Get: id=" << mg->id << " file=" << mg->fps << " imgtype=" << imgtype << endl;
			if ( imgtype == "mg" )
				Tcl_SetStringObj(returnObj, (char *)mg->fmg.c_str(), mg->fmg.length());
			else if ( imgtype == "frame" )
				Tcl_SetStringObj(returnObj, (char *)mg->fframe.c_str(), mg->fframe.length());
			else if ( imgtype == "part" )
				Tcl_SetStringObj(returnObj, (char *)mg->fpart.c_str(), mg->fpart.length());
			else if ( imgtype == "fil" )
				Tcl_SetStringObj(returnObj, (char *)mg->ffil.c_str(), mg->ffil.length());
			else if ( imgtype == "ft" )
				Tcl_SetStringObj(returnObj, (char *)mg->fft.c_str(), mg->fft.length());
			else if ( imgtype == "ps"  )
				Tcl_SetStringObj(returnObj, (char *)mg->fps.c_str(), mg->fps.length());
		} else if ( rec ) {
			if ( imgtype == "rec" )
				Tcl_SetStringObj(returnObj, (char *)rec->frec.c_str(), rec->frec.length());
			else if ( imgtype == "part" )
				Tcl_SetStringObj(returnObj, (char *)rec->fpart.c_str(), rec->fpart.length());
			else if ( imgtype == "fil" )
				Tcl_SetStringObj(returnObj, (char *)rec->ffil.c_str(), rec->ffil.length());
			else if ( imgtype == "ft" )
				Tcl_SetStringObj(returnObj, (char *)rec->fft.c_str(), rec->fft.length());
			else if ( imgtype == "ps"  )
				Tcl_SetStringObj(returnObj, (char *)rec->fps.c_str(), rec->fps.length());
		} else if ( part ) {
			if ( imgtype == "part" )
				Tcl_SetStringObj(returnObj, (char *)part->fpart.c_str(), part->fpart.length());
		}
		imgtype = 0;
	} else if ( property == "select" ) {
		if ( project->select && rec ) {
			Tcl_SetIntObj(returnObj, rec->select);
		} else if ( mg ) {
			Tcl_SetIntObj(returnObj, mg->select);
		} else if ( field ) {
			Tcl_SetIntObj(returnObj, field->select);
		} else if ( part ) {
			returnObj = part_select(part, objc, objv);
		} else {
			Tcl_SetIntObj(returnObj, 0);
		}
	} else if ( property == "fom" ) {
		if ( project->select && rec ) {
			Tcl_SetDoubleObj(returnObj, rec->fom);
		} else if ( mg ) {
			Tcl_SetDoubleObj(returnObj, mg->fom);
		} else if ( part ) {
			returnObj = part_fom(part, objc, objv);
//			Tcl_SetDoubleObj(returnObj, part->fom[0]);
		} else {
			Tcl_SetDoubleObj(returnObj, 0);
		}
/*	} else if ( property == "pixel_size" ) {
		if ( objc > 4 ) imgtype = Tcl_GetStringFromObj(objv[4], NULL);
		if ( field && !mg ) mg = field->mg;
		if ( part ) sam = part->pixel_size;
		else if ( imgtype == "part" ) {
			if ( project->select && rec && rec->part )
				sam = rec->part->pixel_size;
			else if ( mg && mg->part )
				sam = mg->part->pixel_size;
		} else if ( mg ) {
			sam = mg->pixel_size;
		} else if ( rec ) {
			sam = rec->voxel_size;
		}
		snprintf(str, MAXLINELEN, "%g %g %g ", sam[0], sam[1], sam[2]);
		Tcl_AppendToObj(returnObj, str, strlen(str));*/
	} else if ( property == "pixel_size" ) {
		if ( project->select && rec ) {
			snprintf(str, 128, "%g %g %g ", rec->voxel_size[0],  rec->voxel_size[1],  rec->voxel_size[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( mg ) {
			snprintf(str, 128, "%g %g %g ", mg->pixel_size[0],  mg->pixel_size[1],  mg->pixel_size[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( field ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				snprintf(str, 128, "%g %g %g ", mg->pixel_size[0],  mg->pixel_size[1],  mg->pixel_size[2]);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		} else if ( part ) {
			for ( part = project->class_avg; part; part = part->next ) {
				snprintf(str, 128, "%g %g %g ", part->pixel_size[0],  part->pixel_size[1],  part->pixel_size[2]);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		}
	} else if ( property == "part_pixel_size" ) {
		if ( !part ) {
			if ( mg ) part = mg->part;
			else if ( rec ) part = rec->part;
		}
		if ( part ) {
			if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
			if ( filename.length() < 1 || ( mg && filename == mg->fpart ) ||
					( rec && filename == rec->fpart ) ) {
				for ( ; part; part = part->next ) {
					snprintf(str, 128, "%g %g %g ", part->pixel_size[0],  part->pixel_size[1],  part->pixel_size[2]);
					Tcl_AppendToObj(returnObj, str, strlen(str));
				}
			} else {
				for ( ; part; part = part->next ) if ( filename == part->fpart ) {
					snprintf(str, 128, "%g %g %g ", part->pixel_size[0],  part->pixel_size[1],  part->pixel_size[2]);
					Tcl_AppendToObj(returnObj, str, strlen(str));
					break;
				}
			}
		}
	} else if ( property == "size" ) {
		if ( !mg && field ) mg = field->mg;
		if ( mg ) {
			size = micrograph_get_size(mg);
			snprintf(str, 128, "%d %d %d ", size[0],  size[1],  size[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		}
	} else if ( property == "origin" ) {
		if ( project->select && rec ) {
			snprintf(str, 128, "%g %g %g ", rec->origin[0],  rec->origin[1],  rec->origin[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( mg ) {
			snprintf(str, 128, "%g %g %g ", mg->origin[0],  mg->origin[1],  mg->origin[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( field ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				snprintf(str, 128, "%g %g %g ", mg->origin[0],  mg->origin[1],  mg->origin[2]);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		} else if ( part ) {
			for ( part = project->class_avg; part; part = part->next ) {
				snprintf(str, 128, "%g %g %g ", part->ori[0],  part->ori[1],  part->ori[2]);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		}
	} else if ( property == "part_origin" ) {
		if ( !part ) {
			if ( mg ) part = mg->part;
			else if ( rec ) part = rec->part;
		}
		if ( part ) {
			if ( objc > 4 ) filename = Tcl_GetStringFromObj(objv[4], NULL);
			if ( filename.length() < 1 || ( mg && filename == mg->fpart ) ||
					( rec && filename == rec->fpart ) ) {
				for ( ; part; part = part->next ) {
					snprintf(str, 128, "%g %g %g ", part->ori[0],  part->ori[1],  part->ori[2]);
					Tcl_AppendToObj(returnObj, str, strlen(str));
				}
			} else {
				for ( ; part; part = part->next ) if ( filename == part->fpart ) {
					snprintf(str, 128, "%g %g %g ", part->ori[0],  part->ori[1],  part->ori[2]);
					Tcl_AppendToObj(returnObj, str, strlen(str));
					break;
				}
			}
		}
	} else if ( property == "scale" ) {
		if ( project->select && rec ) {
			snprintf(str, 128, "%g %g %g ", rec->scale[0],  rec->scale[1],  rec->scale[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( mg ) {
			snprintf(str, 128, "%g %g %g ", mg->scale[0],  mg->scale[1],  mg->scale[2]);
			Tcl_AppendToObj(returnObj, str, strlen(str));
		} else if ( field ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				snprintf(str, 128, "%g %g %g ", mg->scale[0],  mg->scale[1],  mg->scale[2]);
				Tcl_AppendToObj(returnObj, str, strlen(str));
			}
		} else {
			snprintf(str, 128, "1.0 1.0 1.0 ");
			Tcl_AppendToObj(returnObj, str, strlen(str));
		}
	} else if ( property == "dose" ) {
		if ( !mg ) mg = field->mg;
		if ( mg ) {
			Tcl_SetDoubleObj(returnObj, mg->dose);
		} else {
			if ( mg->intensity < 0.01 ) micrograph_intensity(mg, 1);
			Tcl_SetDoubleObj(returnObj, mg->intensity);
		}
	} else if ( property == "intensity" ) {
		if ( !mg ) mg = field->mg;
		if ( mg ) {
			if ( mg->intensity < 0.01 ) micrograph_intensity(mg, 1);
			Tcl_SetDoubleObj(returnObj, mg->intensity);
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "defocus_deviation" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->defocus_deviation());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->defocus_deviation());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "defocus" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->defocus_average());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->defocus_average());
//			cout << "Micrograph: " << mg->id << " defocus=" << mg->ctf->defocus_average() << endl;
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "astigmatism_angle" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->astigmatism_angle()*180/M_PI);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->astigmatism_angle()*180/M_PI);
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "volt" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->volt());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->volt());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "Cs" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->Cs());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->Cs());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "amp_fac" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, sin(rec->ctf->amp_shift()));
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, sin(mg->ctf->amp_shift()));
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "focal_length" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->focal_length());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->focal_length());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "aperture" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->objective_aperture());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->objective_aperture());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "slit_width" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, rec->ctf->slit_width());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetDoubleObj(returnObj, mg->ctf->slit_width());
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "baseline" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			bstr = rec->ctf->baseline_equation();
			Tcl_SetStringObj(returnObj, bstr.c_str(), bstr.length());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			bstr = mg->ctf->baseline_equation();
			Tcl_SetStringObj(returnObj, bstr.c_str(), bstr.length());
		} else {
			Tcl_SetStringObj(returnObj, "1", 1);
		}
	} else if ( property == "baseline_type" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetIntObj(returnObj, rec->ctf->baseline_type());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetIntObj(returnObj, mg->ctf->baseline_type());
		} else {
			Tcl_SetIntObj(returnObj, 1);
		}
	} else if ( property == "envelope" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			bstr = rec->ctf->envelope_equation();
			Tcl_SetStringObj(returnObj, bstr.c_str(), bstr.length());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			bstr = mg->ctf->envelope_equation();
			Tcl_SetStringObj(returnObj, bstr.c_str(), bstr.length());
		} else {
			bstr = "1.0*exp(0.0*$s2)";
			Tcl_SetStringObj(returnObj, bstr.c_str(), bstr.length());
		}
	} else if ( property == "envelope_type" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_SetIntObj(returnObj, rec->ctf->envelope_type());
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_SetIntObj(returnObj, mg->ctf->envelope_type());
		} else {
			Tcl_SetIntObj(returnObj, 1);
		}
	} else if ( property == "zero" ) {
		if ( project->select && rec ) {
			if ( rec->ctf )
				Tcl_SetDoubleObj(returnObj, rec->ctf->zero(1));
		} else if ( mg ) {
			if ( mg->ctf )
				Tcl_SetDoubleObj(returnObj, mg->ctf->zero(1));
		} else {
			Tcl_SetDoubleObj(returnObj, 0.0);
		}
	} else if ( property == "view" ) {
		if ( mg ) v = View(mg->matrix);
		else if ( rec ) v = rec->view;
		snprintf(str, 128, "%g %g %g %g ", v[0],  v[1],  v[2],  v.angle());
		Tcl_AppendToObj(returnObj, str, strlen(str));
	} else if ( property == "axis" ) {
		if ( !mg && field ) mg = field->mg;
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->tilt_axis*180.0/M_PI);
		else Tcl_SetDoubleObj(returnObj, 0);
	} else if ( property == "tilt" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->tilt_angle*180.0/M_PI);
		else Tcl_SetDoubleObj(returnObj, 0);
	} else if ( property == "level" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->level_angle*180.0/M_PI);
		else Tcl_SetDoubleObj(returnObj, 0);
	} else if ( property == "unitcell" ) {
		if ( mg )
			snprintf(str, 128, "%g %g %g %g ", mg->hvec[0],  mg->hvec[1],  mg->kvec[0],  mg->kvec[1]);
		Tcl_AppendToObj(returnObj, str, strlen(str));
	} else if ( property == "helix_axis" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->helix_axis*180.0/M_PI);
		else Tcl_SetDoubleObj(returnObj, 0.0);
	} else if ( property == "helix_rise" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->helix_rise);
		else Tcl_SetDoubleObj(returnObj, 0.0);
	} else if ( property == "helix_angle" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->helix_angle*180.0/M_PI);
		else Tcl_SetDoubleObj(returnObj, 0.0);
	} else if ( property == "helix_radius" ) {
		if ( mg ) Tcl_SetDoubleObj(returnObj, mg->helix_radius);
		else Tcl_SetDoubleObj(returnObj, 0.0);
	} else if ( property == "box_size" ) {
		if ( project->select && rec ) snprintf(str, 128, "%ld %ld %ld", rec->box_size[0], rec->box_size[1], rec->box_size[2]);
		else {
			if ( !mg ) mg = field->mg;
			if ( mg ) snprintf(str, 128, "%ld %ld %ld", mg->box_size[0], mg->box_size[1], mg->box_size[2]);
		}
		Tcl_SetStringObj(returnObj, str, strlen(str));
	} else if ( property == "bad_radius" ) {
		if ( project->select && rec ) Tcl_SetDoubleObj(returnObj, rec->bad_radius);
		else {
			if ( !mg ) mg = field->mg;
			if ( mg ) Tcl_SetDoubleObj(returnObj, mg->bad_radius);
		}
	} else if ( property == "filament_width" ) {
		if ( project->select && rec ) Tcl_SetDoubleObj(returnObj, rec->filament_width);
		else {
			if ( !mg ) mg = field->mg;
			if ( mg ) Tcl_SetDoubleObj(returnObj, mg->filament_width);
		}
	} else if ( property == "filament_node_radius" ) {
		if ( project->select && rec ) Tcl_SetDoubleObj(returnObj, rec->fil_node_radius);
		else {
			if ( !mg ) mg = field->mg;
			if ( mg ) Tcl_SetDoubleObj(returnObj, mg->fil_node_radius);
		}
	} else if ( property == "sf_radius" ) {
		if ( project->select && rec ) Tcl_SetDoubleObj(returnObj, rec->sf_radius);
		else {
			if ( !mg ) mg = field->mg;
			if ( mg ) Tcl_SetDoubleObj(returnObj, mg->sf_radius);
		}
	} else if ( property == "marker_radius" ) {
		if ( part ) Tcl_SetDoubleObj(returnObj, 0.0);
		else if ( project->select && rec ) Tcl_SetDoubleObj(returnObj, rec->mark_radius);
		else {
			if ( !mg && field ) mg = field->mg;
			if ( mg ) Tcl_SetDoubleObj(returnObj, mg->mark_radius);
		}
	} else {
		cerr << "Error in do_get: Property " << property << " not recognized!" << endl;
	}
	
	return returnObj;
}


int			do_set(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	if ( !project )
		return error_show("Error: No project defined!\n", __FILE__, __LINE__);
	
	int					update(0);
	
	int					i, n(-1), err(0);
	double				value;
	Vector3<int>		box;
	Vector3<double>		origin, sam(1,1,1), scale;
	
	Bstring				name, imgtype;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);

	Bstring				property = Tcl_GetStringFromObj(objv[3], NULL);

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG do_set: Item: " << item << " (" << item.length() << ")" << endl;
		cout << "DEBUG do_set: ID: " << id << " (" << id.length() << ")" << endl;
		cout << "DEBUG do_set: Property: " << property << " (" << property.length() << ")" << endl;
	}
	
	if ( item.contains("all" ) ) {
	} else if ( property == "id" ) {
		if ( item.contains("Field") ) {
			field = project->field;
		} else if ( item.contains("Micrograph") ) {
			field = project->field;
			if ( field ) mg = field->mg;
		} else if ( item.contains("Reconstruction") ) {
			project->select = 1;
			rec = project->rec;
		} else err++;
	} else if ( item.contains("Field") ) {
		for ( field = project->field; field && field->id != id; field = field->next ) ;
		if ( !field ) err++;
//		cout << "field id = " << project->field->id << endl;
	} else if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next) ;
			if ( mg ) break;
		}
		if ( !mg ) err++;
//		cout << "set mg id: " << mg->id << endl;
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( !rec ) rec = project->rec;
		if ( !rec ) err++;
	} else if ( item.contains("Class") ) {
		part = project->class_avg;
	} else {
		cerr << "Error in do_set: item \"" << item << "\" not supported! (property: " << property << ")" << endl;
		return update;
	}

	if ( err ) cerr << "Error in do_set: " << item << " not found! (property: " << property << ")" << endl;

	if ( !project || property.empty() || objc < 4 ) {
		return update;
	} else if ( property == "active" ) {
		Tcl_GetIntFromObj(NULL, objv[4], &project->select);
		if ( project->select && !project->rec ) rec = reconstruction_add(&project->rec, id);
	} else if ( property == "id" ) {
		id = Tcl_GetStringFromObj(objv[4], NULL);
//		cout << "Setting item " << item << " to id " << id << endl;
		if ( project->select ) {
			if ( rec ) {
				rec->id = id;
			} else {
				rec = project->rec;
				if ( !rec->next ) rec->id = id;
				else for ( n=0, rec = project->rec; rec; rec = rec->next, n++ ) {
					rec->id = id + Bstring(n, "_%03d");
				}
			}
		} else if ( !mg && field ) {
			field->id = id;
		} else if ( mg ) {
			mg = field->mg;
			if ( !mg->next ) mg->id = id;
			else for ( n=0, mg = field->mg; mg; mg = mg->next, n++ ) {
				mg->id = id + Bstring(n, "_%03d");
			}
		}
	} else if ( property == "filename" ) {
		name = Tcl_GetStringFromObj(objv[4], NULL);
		imgtype = Tcl_GetStringFromObj(objv[5], NULL);
//		cout << "filename=" << name << " imgtype=" << imgtype << endl;
		if ( mg ) {
//			cout << "Set: id=" << id << " file=" << name << " imgtype=" << imgtype << endl;
			if ( imgtype == "mg" ) mg->fmg = name;
			else if ( imgtype == "frame" ) mg->fframe = name;
			else if ( imgtype == "part" ) mg->fpart = name;
			else if ( imgtype == "fil" ) mg->ffil = name;
			else if ( imgtype == "ft" ) mg->fft = name;
			else if ( imgtype == "ps"  ) mg->fps = name;
		} else if ( field ) {
			if ( imgtype == "mg" ) {
				for ( i=0, mg = field->mg; mg; mg = mg->next, i++ ) {
					mg->fmg = name;
					mg->img_num = i;
				}
			} else if ( imgtype == "ps" ) {
				for ( mg = field->mg; mg; mg = mg->next )
					mg->fps = name;
			}
		} else if ( rec ) {
			if ( imgtype == "rec" ) rec->frec = name;
			else if ( imgtype == "part" ) rec->fpart = name;
			else if ( imgtype == "fil" ) rec->ffil = name;
			else if ( imgtype == "ft" ) rec->fft = name;
			else if ( imgtype == "ps"  ) rec->fps = name;
		} else if ( part ) {
			if ( imgtype == "part" ) part->fpart = name;
			if ( part == project->class_avg )
				for ( ; part; part = part->next ) part->fpart = name;
		}
		imgtype = 0;
		name = 0;
	} else if ( property == "select" ) {
		Tcl_GetIntFromObj(NULL, objv[4], &i);
		if ( project->select && rec ) rec->select = i;
		else if ( mg ) mg->select = i;
		else if ( field ) field->select = i;
		else if ( part ) part->sel = i;
	} else if ( property == "pixel_size" ) {
//		cout << objc << endl;
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &sam[0]);
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &sam[1]);
		if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &sam[2]);
//		cout << "pixel size = " << sam << endl;
		if ( project->select && rec )
			rec->voxel_size = sam;
		else if ( mg )
			mg->pixel_size = sam;
		else if ( field ) for ( mg = field->mg; mg; mg = mg->next )
			mg->pixel_size = sam;
	} else if ( property == "origin" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		origin[0] = value;
		Tcl_GetDoubleFromObj(NULL, objv[5], &value);
		origin[1] = value;
		Tcl_GetDoubleFromObj(NULL, objv[6], &value);
		origin[2] = value;
		if ( project->select && rec ) rec->origin = origin;
		else if ( mg ) mg->origin = origin;
		else if ( field ) for ( mg = field->mg; mg; mg = mg->next ) mg->origin = origin;
	} else if ( property == "scale" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		scale[0] = value;
		Tcl_GetDoubleFromObj(NULL, objv[5], &value);
		scale[1] = value;
		Tcl_GetDoubleFromObj(NULL, objv[6], &value);
		scale[2] = value;
		if ( project->select && rec ) rec->scale = scale;
		else if ( mg ) mg->scale = scale;
		else if ( field ) for ( mg = field->mg; mg; mg = mg->next ) mg->scale = scale;
	} else if ( property == "dose" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( item.contains("all" ) )
			project_set_dose(project, value);
		else if ( mg )
			mg->dose = value;
		else if ( field ) for ( mg = field->mg; mg; mg = mg->next )
			mg->dose = value;
	} else if ( property == "defocus_deviation" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->defocus_deviation(value);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->defocus_deviation(value);
		}
	} else if ( property == "defocus" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->defocus_average(value);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->defocus_average(value);
		}
	} else if ( property == "astigmatism_angle" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		value *= M_PI/180.0;
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->astigmatism_angle(value);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->astigmatism_angle(value);
		}
	} else if ( property == "volt" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_volts(project, value);
	} else if ( property == "Cs" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_Cs(project, value);
	} else if ( property == "amp_fac" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_amp_shift(project, asin(value));
	} else if ( property == "focal_length" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_focal_length(project, value);
	} else if ( property == "aperture" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_aperture(project, value);
	} else if ( property == "slit_width" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		project_set_slit_width(project, value);
	} else if ( property == "baseline_type" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_GetIntFromObj(NULL, objv[4], &i);
			rec->ctf->baseline_type(i);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_GetIntFromObj(NULL, objv[4], &i);
			mg->ctf->baseline_type(i);
//			cout << "Micrograph: " << mg->id << " baseline_type: " << mg->ctf->baseline_type() << endl;
		}
	} else if ( property == "baseline" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->parse_baseline_equation(Tcl_GetStringFromObj(objv[4], NULL));
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->parse_baseline_equation(Tcl_GetStringFromObj(objv[4], NULL));
		}
	} else if ( property == "envelope_type" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			Tcl_GetIntFromObj(NULL, objv[4], &i);
			rec->ctf->envelope_type(i);
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			Tcl_GetIntFromObj(NULL, objv[4], &i);
			mg->ctf->envelope_type(i);
//			cout << "Micrograph: " << mg->id << " envelope_type: " << mg->ctf->envelope_type() << endl;
		}
	} else if ( property == "envelope" ) {
		if ( project->select && rec ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->parse_envelope_equation(Tcl_GetStringFromObj(objv[4], NULL));
		} else if ( mg ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->parse_envelope_equation(Tcl_GetStringFromObj(objv[4], NULL));
		}
	} else if ( property == "axis" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		value *= M_PI/180.0;
		if ( field ) {
			if ( objc > 5 ) {
				Tcl_GetIntFromObj(NULL, objv[5], &n);
				for ( i=0, mg = field->mg; mg && i<n; mg = mg->next, i++ ) ;
				if ( mg ) mg->tilt_axis = value;
			} else {
				for ( mg = field->mg; mg; mg = mg->next ) mg->tilt_axis = value;
			}
		} else if ( mg ) {
			mg->tilt_axis = value;
		}
		project_mg_tilt_to_matrix(project);
	} else if ( property == "tilt" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		value *= M_PI/180.0;
//		cout << "Trying to set tilt angle for item " << item << " to " << value*180/M_PI << endl;
		if ( !mg && field ) {
			Tcl_GetIntFromObj(NULL, objv[5], &n);
			for ( i=0, mg = field->mg; mg && i<n; mg = mg->next, i++ ) ;
//			cout << "Trying to set tilt angle for mg " << n << " to " << value*180/M_PI << endl;
		}
//		if ( mg ) cout << "Set tilt for mg " << n << " to " << value*180/M_PI << endl;
		if ( mg ) mg->tilt_angle = value;
	} else if ( property == "level" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		value *= M_PI/180.0;
		if ( !mg && field ) {
			Tcl_GetIntFromObj(NULL, objv[5], &n);
			for ( i=0, mg = field->mg; mg && i<n; mg = mg->next, i++ ) ;
		}
		if ( mg ) mg->level_angle = value;
	} else if ( property == "unitcell" ) {
		if ( mg ) {
			Tcl_GetDoubleFromObj(NULL, objv[4], &value);
			mg->hvec[0] = value;
			Tcl_GetDoubleFromObj(NULL, objv[5], &value);
			mg->hvec[1] = value;
			Tcl_GetDoubleFromObj(NULL, objv[6], &value);
			mg->kvec[0] = value;
			Tcl_GetDoubleFromObj(NULL, objv[7], &value);
			mg->kvec[1] = value;
		}
	} else if ( property == "helix_axis" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( mg ) mg->helix_axis = value*M_PI/180.0;
	} else if ( property == "helix_rise" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( mg ) mg->helix_rise = value;
	} else if ( property == "helix_angle" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( mg ) mg->helix_angle = value*M_PI/180.0;
	} else if ( property == "helix_radius" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( mg ) mg->helix_radius = value;
	} else if ( property == "box_size" ) {
		if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &box[0]);
		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &box[1]);
		if ( objc > 6 ) Tcl_GetIntFromObj(NULL, objv[6], &box[2]);
		project_set_particle_box_size(project, box);
/*		if ( project->select && project->rec ) {
			for ( rec = project->rec; rec; rec = rec->next ) rec->box_size = box;
		} else if ( mg ) {
			box[2] = 1;
			mg->box_size = box;
		} else if ( field ) {
			box[2] = 1;
			for ( mg = field->mg; mg; mg = mg->next )
				mg->box_size = box;
		}*/
	} else if ( property == "bad_radius" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec )
			for ( rec=project->rec; rec; rec=rec->next ) rec->bad_radius = value;
		else if ( mg ) mg->bad_radius = value;
		else if ( field )
			for ( mg = field->mg; mg; mg = mg->next ) mg->bad_radius = value;
	} else if ( property == "filament_width" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec )
			for ( rec=project->rec; rec; rec=rec->next ) rec->filament_width = value;
		else if ( mg ) mg->filament_width = value;
		else if ( field )
			for ( mg = field->mg; mg; mg = mg->next ) mg->filament_width = value;
	} else if ( property == "filament_node_radius" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec )
			for ( rec=project->rec; rec; rec=rec->next ) rec->fil_node_radius = value;
		else if ( mg ) mg->fil_node_radius = value;
		else if ( field )
			for ( mg = field->mg; mg; mg = mg->next ) mg->fil_node_radius = value;
	} else if ( property == "sf_radius" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec )
			for ( rec=project->rec; rec; rec=rec->next ) rec->sf_radius = value;
		else if ( mg ) mg->sf_radius = value;
		else if ( field )
			for ( mg = field->mg; mg; mg = mg->next ) mg->sf_radius = value;
	} else if ( property == "marker_radius" ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &value);
		if ( project->select && rec )
			for ( rec=project->rec; rec; rec=rec->next ) rec->mark_radius = value;
		else if ( mg ) mg->mark_radius = value;
		else if ( field )
			for ( mg = field->mg; mg; mg = mg->next ) mg->mark_radius = value;
	} else {
		cerr << "Error in do_set: Property " << property << " not recognized!" << endl;
	}
	
	return update;
}

Tcl_Obj*	do_rps(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*		returnObj = Tcl_NewObj();

	int					err(0);
	long				img_num(0);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	CTFparam* 			em_ctf = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_rps: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG do_rps: Item: " << item << " (" << item.length() << ")" << endl;
		cout << "DEBUG do_rps: ID: " << id << " (" << id.length() << ")" << endl;
	}
		
	if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next) ;
			if ( mg ) break;
		}
		if ( !mg ) err++;
		else {
			img_num = mg->img_num;
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			em_ctf = mg->ctf;
			filename = mg->fps;
			if ( filename.contains("/") ) filename = filename.post_rev('/');
			for ( p = imglist; p; p = p->next ) {
				filename2 = p->file_name();
				if ( filename2.contains(filename) ) break;
			}
			if ( p ) p->sampling(mg->pixel_size);
			else err++;
		}
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( !rec ) err++;
		else {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			em_ctf = rec->ctf;
			filename = rec->fps;
			if ( filename.contains("/") ) filename = filename.post_rev('/');
			for ( p = imglist; p; p = p->next ) {
				filename2 = p->file_name();
				if ( filename2.contains(filename) ) break;
			}
			if ( p ) p->sampling(rec->voxel_size);
			else err++;
		}
	} else {
		cerr << "Error in do_rps: item \"" << item << "\" not supported!" << endl;
		return returnObj;
	}
	
	if ( !p ) {
		cerr << "Warning in do_rps: No image! " << endl;
		return returnObj;
	}
	
	if ( img_num >= p->images() ) img_num = 0;

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG do_rps: image id = " << p->identifier() << endl;
		cout << "DEBUG do_rps: image = " << p->file_name() << tab << img_num << endl;
	}

	if ( em_ctf->defocus_average() <= 1 ||
		em_ctf->defocus_average() >= 2e5 )
			em_ctf->defocus_average(2e4);
	
/*	if ( objc > 3 ) {
		Tcl_GetDoubleFromObj(NULL, objv[3], &defocus_avg);
		em_ctf->defocus_average(defocus_avg);
	}
		
	if ( objc > 4 ) {
		Tcl_GetDoubleFromObj(NULL, objv[4], &defocus_dev);
		em_ctf->defocus_deviation(defocus_dev);
	}
		
	if ( objc > 5 ) {
		Tcl_GetDoubleFromObj(NULL, objv[5], &astigmatism_angle);
		em_ctf->astigmatism_angle(astigmatism_angle);
	}
*/	
	char			str[128];
	Bimage*			prps = img_ctf_radial_average(p, img_num, *em_ctf);
	
	if ( !prps ) {
		snprintf(str, 128, "Warning: Radial average failed!");
		Tcl_AppendToObj(returnObj, str, strlen(str));
		return returnObj;
	}
	
	if ( prps->sizeX() < 1 ) {
		snprintf(str, 128, "Warning: Radial average has zero values!");
		Tcl_AppendToObj(returnObj, str, strlen(str));
		return returnObj;
	}
	
	long   			i;

	double			wri = img_water_ring_index(prps);
	
	for ( i=0; i<prps->sizeX(); i++ ) {
		if ( !isfinite((*prps)[i]) ) {
			cerr << "Warning in do_rps: Point " << i << " is a NaN!" << endl;
			Tcl_AppendToObj(returnObj, str, strlen(str));
			return returnObj;
		}
		snprintf(str, 128, " %g", (*prps)[i]);
		Tcl_AppendToObj(returnObj, str, strlen(str));
//		cout << i << tab << (*prps)[i] << endl;
	}
	
	delete prps;

	snprintf(str, 128, " %g", wri);
	Tcl_AppendToObj(returnObj, str, strlen(str));
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_rps: done" << endl;

	return returnObj;
}

Tcl_Obj*	do_emfp(int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*	returnObj = Tcl_NewObj();

	double				v(0), emfp(0);
	Bstring				material_str("Vitreous ice");
	CTFparam 			ctf;
	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &v);
	ctf.volt(v);
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &v);
	ctf.focal_length(v);
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &v);
	ctf.objective_aperture(v);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &v);
	ctf.slit_width(v);
	if ( objc > 6 ) material_str = Tcl_GetStringFromObj(objv[6], NULL);
		
	material_str = "\"" + material_str + "\"";
	
//	cout << "Material requested: " << material_str << endl;
	
	Bstring					propfile("material.star");
	map<string,Bmaterial>	mprop = read_material_properties(propfile);
	
	if ( mprop.size() < 1 ) {
		cerr << "Error: No material parameter file read!" << endl;
		return returnObj;
	}

	if ( mprop.find(material_str.str()) == mprop.end() ) {
		material_str = "Vitreous ice";
		material_str = "\"" + material_str + "\"";
	}
	
	Bmaterial		material = mprop[material_str.str()];

	cout << "Material set: " << material.identifier() << endl;
	
	material.show();
	
	ctf.show();
	
	emfp = effective_mean_free_path(material, ctf);
	
	cout << "EMFP = " << emfp << endl;

	char				str[MAXLINELEN] = "";
	snprintf(str, MAXLINELEN, "%g %s ", emfp, material.identifier().c_str());
	Tcl_AppendToObj(returnObj, str, strlen(str));

//	cout << "Return list: " << str << endl;

	return returnObj;
}

Tcl_Obj*	do_material_list()
{
	Tcl_Obj*	returnObj = Tcl_NewObj();

	Bstring					propfile("material.star");
	map<string,Bmaterial>	mprop = read_material_properties(propfile);
	
	char				str[MAXLINELEN] = "";
	
	for ( auto m: mprop ) {
		snprintf(str, MAXLINELEN, "%s ", m.first.c_str());
		Tcl_AppendToObj(returnObj, str, strlen(str));
	}
	
	return returnObj;
}

Tcl_Obj*	do_ctf_fit(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*	returnObj = Tcl_NewObj();

	int					err(0);
	int					level(0);		// 0=quick, 1=baseline, 2=defocus, 3=astigmatism
	long				img_num(0);
	double				lores(1e10), hires(0.1);
	double 				def_start(1e3), def_end(2e5), def_inc(1e3);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	CTFparam* 			em_ctf = NULL;
	Bimage*				p = NULL;
	Bstring				filename, filename2;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_ctf_fit: Item: " << item << " (" << item.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
//	cout << "id = " << id << endl;
	
	if ( item.contains("Micrograph") ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg && mg->id != id; mg = mg->next) ;
			if ( mg ) break;
		}
		if ( !mg ) err++;
		else {
			img_num = mg->img_num;
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			em_ctf = mg->ctf;
			filename = mg->fps;
			if ( filename.contains("/") ) filename = filename.post_rev('/');
			for ( p = imglist; p; p = p->next ) {
				filename2 = p->file_name();
				if ( filename2.contains(filename) ) break;
			}
			if ( p ) p->sampling(mg->pixel_size);
			else err++;
		}
	} else if ( item.contains("Reconstruction") ) {
		for ( rec = project->rec; rec && rec->id != id; rec = rec->next ) ;
		if ( !rec ) err++;
		else {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			em_ctf = rec->ctf;
			filename = rec->fps;
			if ( filename.contains("/") ) filename = filename.post_rev('/');
			for ( p = imglist; p; p = p->next ) {
				filename2 = p->file_name();
				if ( filename2.contains(filename) ) break;
			}
			if ( p ) p->sampling(rec->voxel_size);
			else err++;
		}
	} else {
		cerr << "Error in do_ctf_fit: item \"" << item << "\" not supported!" << endl;
		return returnObj;
	}

	if ( !p ) {
		cerr << "Warning in do_ctf_fit: No image! " << endl;
		return returnObj;
	}
	
	if ( img_num >= p->images() ) img_num = 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_rps: image id = " << p->identifier() << endl;

	if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &level);
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &lores);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &hires);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &def_start);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &def_end);
	if ( objc > 8 ) Tcl_GetDoubleFromObj(NULL, objv[8], &def_inc);

	if ( level == 0 || em_ctf->defocus_average() == 0 ) {
		em_ctf->defocus_deviation(0);
		em_ctf->astigmatism_angle(0);
		img_ctf_find_defocus(p, img_num, *em_ctf, lores, hires, def_start, def_end, def_inc);
		img_ctf_fit_baseline(p, img_num, *em_ctf, lores, hires);
		img_ctf_fit_envelope(p, img_num, *em_ctf, lores, hires);
	}
	
	if ( level == 1 )
		img_ctf_fit_baseline(p, img_num, *em_ctf, lores, hires);

	if ( level == 2 )
		img_ctf_fit_envelope(p, img_num, *em_ctf, lores, hires);

	if ( level == 3 )
		img_ctf_find_defocus(p, img_num, *em_ctf, lores, hires,
			em_ctf->defocus_average()/1.2, em_ctf->defocus_average()*1.2,
			em_ctf->defocus_average()*0.2);

	if ( level == 4 )
		img_ctf_fit_astigmatism(p, img_num, *em_ctf, lores, hires);

	double		wri = img_water_ring_index(p, img_num, *em_ctf);
	if ( mg ) mg->wri = wri;
	
	Bstring		base_eq = em_ctf->baseline_equation();
	Bstring		env_eq = em_ctf->envelope_equation();
	
//	cout << base_eq << endl;

	char		str[MAXLINELEN] = "";
	if ( base_eq.length() )
		snprintf(str, MAXLINELEN, "\"%s\" \"%s\" %g %g %g %g", base_eq.c_str(), env_eq.c_str(),
				em_ctf->defocus_average(), em_ctf->defocus_deviation(), em_ctf->astigmatism_angle()*180.0/M_PI, wri);

	Tcl_SetStringObj(returnObj, str, strlen(str));
	
	return returnObj;
}

Tcl_Obj*	do_mg_sort(Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();
	char				str[MAXLINELEN] = "";
	Bstring				mg_sort, tag;

	if ( objc > 3 ) mg_sort = Tcl_GetStringFromObj(objv[3], NULL);
	if ( mg_sort == "Defocus" ) tag = CTF_DEF_AVG;
	if ( mg_sort == "Water ring" ) tag = MICROGRAPH_WATER_RING;
	if ( mg_sort == "Particles" ) tag = PARTICLE;
	if ( mg_sort == "Intensity" ) {
		project_mg_avg_intensities(project);
		tag = MICROGRAPH_INTENSITY;
	}
	
	vector<pair<Bmicrograph*,double>>	mgv = project_mg_sort(project, tag);
	
	if ( mgv.size() < 1 ) return returnObj;
	
	for ( auto it = mgv.begin(); it != mgv.end(); ++it ) {
		snprintf(str, MAXLINELEN, "%s %20.16f ", it->first->id.c_str(), it->second);
		Tcl_AppendToObj(returnObj, str, strlen(str));
	}

	return returnObj;
}

