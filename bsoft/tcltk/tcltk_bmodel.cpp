/**
@file	tcltk_bmodel.cpp
@brief	A shared object to load Bsoft model files in TCL/Tk
@author Bernard Heymann
@date	Created: 20071002
@date	Modified: 20191010
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bmodel.h"
#include "rwmodel.h"
#include "rwimg.h"
#include "model_create.h"
#include "model_map.h"
#include "model_select.h"
#include "model_mask.h"
#include "model_color.h"
#include "model_views.h"
#include "mg_processing.h"
#include "Vector3.h"
#include "linked_list.h"
#include "Color.h"
#include "utilities.h"

// Declaration of global variables
extern int			verbose;		// Level of output to the screen
extern Bimage*		imglist;
extern Bmodel*		model;
extern Bproject*	project;

// Internal function prototypes
Tcl_Obj*	do_find(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_get(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
int			do_set(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
int			do_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
int			do_delete_non_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
int			do_components_to_particles(Bmodel* model, Bproject* project, int objc, Tcl_Obj *CONST objv[]);
int			do_extract_segments(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
int			do_create_shell(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_component_type(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_component(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_link(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);

Tcl_Obj*	comptype_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_list(Bmodel* model);
Tcl_Obj*	comptype_create_update(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_deselect(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	comptype_first_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);

Tcl_Obj*	component_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_location(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_distance(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_img_coords(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_fom(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_radius(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_move(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_create(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	component_color(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);

Tcl_Obj*	link_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_img_coords(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_create(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	link_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[]);

/**
@brief 	Implements the "Bmodel" command in Tcl/Tk to access model parameter files through Bsoft.
@param 	*interp		a Tcl interpreter within Tcl.
@param 	objc		number of arguments passed (+1).
@param 	*objv[]		arguments passed as Tcl objects.
@return int			Tcl result.

	Bmodel command syntax:
		Bmodel <action> <arguments>.
		where:
			action			"create", "exists", "read", "write", "kill", "get", "set", "delete"
							"delete_non_selected", "components_to_particles", "extract_segments",
							"comptype", "component", "link"
			arguments			action-specific arguments:
				"create"		<id>
				"read"			<filename>
				"write"			<filename>
				"get"			<property> [arguments]
				"set"			<property> <value>
				"delete"		<property>
				"delete_non_selected"	<property> <selection>
				"components_to_particles"
				"extract_segments"	<filename> [multilevel]
				"comptype"		[arguments]
				"component"		[arguments]
				"link"			[arguments]
				where:
					property	"id <string>"
								"map <string> <value>"
								"number <value>"
								"radius <value>"
								"linkradius <value>"
	Return values:
		Each action may have a return value:
			"create"	(none)
			"exists"	0=no, 1=yes
			"read"		model id
			"write"		(none)
            "kill"		(none)
			"get"		return value based on property
			"set"		modify model property

**/
int			model_processing(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
//	verbose = 967;

	Tcl_ResetResult(interp);
	
	if ( objc < 2 ) {
		Tcl_AppendResult(interp, "wrong # args", (char *)NULL);
		return TCL_ERROR;
	}

	Bstring				action = Tcl_GetStringFromObj(objv[1], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_processing: Action: " << action << " (" << action.length() << ")" << endl;
	
	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	char				cstr[128] = " ";
	Tcl_SetStringObj(returnObj, cstr, 1);
	
	Bstring				filename, id("1");

	time_t				ti = time(NULL);
	
	if ( action == "create" ) {
		model = model_from_images(imglist);
		id = "Model:" + model->identifier() ;
		Tcl_SetStringObj(returnObj, (char *)id.c_str(), id.length());
	} else if ( action == "exists" ) {
		if ( model ) Tcl_SetIntObj(returnObj, 1);
		else Tcl_SetIntObj(returnObj, 0);
	} else if ( action == "read" ) {
		if ( objc < 2 ) {
			Tcl_AppendResult(interp, "No file name given!", (char *)NULL);
			return TCL_ERROR;
		}
		filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG model_processing: File name: " << filename << " (" << filename.length() << ")" << endl;
		if ( model ) model_kill(model);
		model = read_model(filename);
		Tcl_SetStringObj(returnObj, (char *)model->identifier() .c_str(), model->identifier() .length());
	} else if ( action == "write" ) {
		if ( objc < 2 ) {
			Tcl_AppendResult(interp, "No file name given!", (char *)NULL);
			return TCL_ERROR;
		}
		filename = Tcl_GetStringFromObj(objv[2], NULL);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG BmodelCmd: File name: " << filename << " (" << filename.length() << ")" << endl;
		if ( !model ) {
			Tcl_AppendResult(interp, "No model in memory!", (char *)NULL);
			return TCL_ERROR;
		}
		string		com = "\n# Written from Bshow\n# ";
		model->comment() += com + asctime(localtime(&ti)) + "\n";
		write_model(filename, model);
	} else if ( action == "kill" ) {
		model_kill(model);
	} else if ( action == "find" ) {
		returnObj = do_find(model, objc, objv);
	} else if ( action == "get" ) {
		returnObj = do_get(model, objc, objv);
	} else if ( action == "set" ) {
		do_set(model, objc, objv);
	} else if ( action == "delete_non_selected" ) {
		do_delete_non_selected(model, objc, objv);
	} else if ( action == "delete" ) {
		do_delete(model, objc, objv);
	} else if ( action == "components_to_particles" ) {
		do_components_to_particles(model, project, objc, objv);
	} else if ( action == "extract_segments" ) {
		do_extract_segments(model, objc, objv);
	} else if ( action == "create_shell" ) {
		do_create_shell(model, objc, objv);
	} else if ( action == "comptype" ) {
		returnObj = do_component_type(model, objc, objv);
	} else if ( action == "component" ) {
		returnObj = do_component(model, objc, objv);
	} else if ( action == "link" ) {
		returnObj = do_link(model, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}
	
	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}

Tcl_Obj*	do_find(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					number(0);
	Bstring				name;
	Bmodel*				mp = model;

	if ( objc > 2 ) name = Tcl_GetStringFromObj(objv[2], NULL);
	if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &number);

	for ( mp = model; mp; mp = mp->next ) {
		if ( mp->mapfile().find(name.str()) != string::npos
				|| name.contains(mp->mapfile()) )
			if ( mp->image_number() == number ) break;
	}

	if ( !mp ) return returnObj;
	
	Bstring				bstr = mp->mapfile() + "::Model:" + mp->identifier() + " ";
	
	Tcl_AppendToObj(returnObj, bstr.c_str(), bstr.length());
	
	return returnObj;
}

Tcl_Obj*	do_get(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					err(0);
	Bmodel*				mp = model;
	Bstring				bstr, imgtype;
	
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
	
	if ( item.contains("all" ) ) {
		mp = model;
	} else if ( item.contains("Model") ) {
		for ( mp = model; mp && mp->identifier() != id.str(); mp = mp->next ) ;
		if ( !mp ) err++;
	} else {
		cerr << "Error in do_get: item \"" << item << "\" not supported!" << endl;
		return returnObj;
	}
	
	if ( !mp || property.empty() ) {
		Tcl_SetIntObj(returnObj, 0);
	} else if ( property == "id" ) {
		Tcl_SetStringObj(returnObj, (char *)mp->identifier() .c_str(), mp->identifier() .length());
	} else if ( property == "selection" ) {
		Tcl_SetIntObj(returnObj, mp->select());
	} else if ( property == "map" ) {
		Tcl_SetStringObj(returnObj, (char *)mp->mapfile().c_str(), mp->mapfile().length());
	} else if ( property == "mask" ) {
		Tcl_SetStringObj(returnObj, (char *)mp->maskfile().c_str(), mp->maskfile().length());
	} else if ( property == "number" ) {
		Tcl_SetIntObj(returnObj, mp->image_number());
	} else if ( property == "image_filenames" ) {
		for ( mp = model; mp; mp = mp->next ) {
			if ( mp->mapfile().length() && mp->mapfile() != "?" ) {
				bstr = mp->mapfile() + "::Model:" + mp->identifier() + " ";
				Tcl_AppendToObj(returnObj, bstr.c_str(), bstr.length());
//				cout << bstr << endl;
			}
		}
//	} else if ( property == "types" ) {
//		returnObj = component_type_list(model);
	} else {
		cerr << "Error in do_get: Property " << property << " not recognized!" << endl;
	}
	
	return returnObj;
}

int			set_component_radius(Bmodel* model, float radius)
{
	if ( radius < 0.001 ) return 0;
	
	Bcomponent*		comp = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG set_component_radius: radius = " << radius << endl;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() )
		comp->radius(radius);

	return 0;
}

int			set_link_radius(Bmodel* model, float radius)
{
	if ( radius < 0.001 ) return 0;
	
	Blink*			link = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG set_link_radius: radius = " << radius << endl;

	for ( link = model->link; link; link = link->next ) if ( link->select() )
		link->radius(radius);

	return 0;
}

int			do_set(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	if ( !model ) return 0;
	
	int					err(0);
	int					i, rtn(0);
	double				radius(1);
	Bstring				calc_views;
	Bmodel*				mp = model;
	
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
		mp = model;
	} else if ( item.contains("Model") ) {
		for ( mp = model; mp && mp->identifier() != id.str(); mp = mp->next ) ;
		if ( !mp ) err++;
	} else {
		cerr << "Error in do_get: item \"" << item << "\" not supported!" << endl;
		return rtn;
	}
	
	if ( property == "id" ) {
		if ( objc > 4 )
			mp->identifier() = Tcl_GetStringFromObj(objv[4], NULL);
	} else if ( property == "map" ) {
		if ( objc > 4 )
			mp->mapfile(Tcl_GetStringFromObj(objv[4], NULL));
		if ( objc > 5 ) {
			Tcl_GetIntFromObj(NULL, objv[5], &i);
			mp->image_number(i);
		}
	} else if ( property == "mask" ) {
		if ( objc > 4 )
			mp->maskfile(Tcl_GetStringFromObj(objv[4], NULL));
	} else if ( property == "number" ) {
		if ( objc > 4 ) {
			Tcl_GetIntFromObj(NULL, objv[4], &i);
			mp->image_number(i);
		}
	} else if ( property == "compradius" ) {
		if ( objc > 5 ) {
			component_radius(mp, objc, objv);
		} else if ( objc > 4 ) {
			Tcl_GetDoubleFromObj(NULL, objv[4], &radius);
			set_component_radius(mp, radius);
		}
	} else if ( property == "linkradius" ) {
		if ( objc > 4 ) {
			Tcl_GetDoubleFromObj(NULL, objv[4], &radius);
			set_link_radius(mp, radius);
		}
	} else if ( property == "views" ) {
		if ( objc > 4 )
			calc_views = Tcl_GetStringFromObj(objv[4], NULL);
		model_calculate_views(mp, calc_views);
	} else {
		cerr << "Error in do_set: Property " << property << " not recognized!" << endl;
	}
	
	return rtn;
}

int			do_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	if ( !model ) return 0;
	
	model_link_list_kill(model);
	component_list_kill(model->comp);
	
	model->comp = NULL;
	model->link = NULL;
	
	return 0;
}

int			do_delete_non_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	model_delete_non_selected(&model);

	return 0;
}

int			do_components_to_particles(Bmodel* model, Bproject* project, int objc, Tcl_Obj *CONST objv[])
{
	if ( !project ) return 0;
	if ( project->select < 1 ) return 0;
	
	Bstring				rec_id, filename2;
	Breconstruction*	rec = project->rec;
/*
	Bimage*				p = NULL;

	Bstring				b = base(model->mapfile());
	for ( p = imglist; p; p = p->next ) {
		filename2 = p->file_name();
		if ( filename2.base() == b ) break;
	}
*/
	Bimage*				p = imglist->find(model->mapfile());

	if ( !p ) return 0;
	
	if ( !project->rec ) {
		rec_id = "1";
		rec = reconstruction_add(&project->rec, rec_id);
		rec->frec = p->file_name();
	}
	
	int					id(0);
	Bcomponent*			comp = NULL;
	Bparticle*			part = NULL;
	
	particle_kill(rec->part);
	rec->part = NULL;
	
	if ( rec->box_size.volume() < 1 ) rec->box_size = Vector3<float>(10,10,10);
	if ( rec->bad_radius < 1 ) rec->bad_radius = 5;
	if ( rec->voxel_size < 0.001 ) rec->voxel_size = p->sampling(0)[0];
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		id = stol(comp->identifier());
		part = particle_add(&project->rec->part, id);
		part->group = stol(comp->type()->identifier());
		part->loc[0] = comp->location()[0]/p->sampling(0)[0] + p->image->origin()[0];
		part->loc[1] = comp->location()[1]/p->sampling(0)[1] + p->image->origin()[1];
		part->loc[2] = comp->location()[2]/p->sampling(0)[2] + p->image->origin()[2];
		part->ori = rec->box_size/2;
		part->view[2] = 1;
		part->fom[0] = comp->FOM();
		part->sel = comp->select();
	}
	
	return 0;
}

int			do_extract_segments(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	if ( !model ) return 0;
	
	Bstring				filename, filename2;
	int					multi_level(0);
	
	if ( objc > 2 ) filename = Tcl_GetStringFromObj(objv[2], NULL);
	if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &multi_level);
/*
	Bimage*				p = NULL;

	Bstring				b = base(model->mapfile());
	for ( p = imglist; p; p = p->next ) {
		filename2 = p->file_name();
		if ( filename2.base() == b ) break;
	}
*/
	Bimage*				p = imglist->find(model->mapfile());

	if ( !p ) {
		cerr << "Error: No map found!" << endl;
		return 0;
	}
	
	Bimage*				pseg = img_extract_segments_using_model(p, model, multi_level);
	
	write_img(filename, pseg, 0);
	
	delete pseg;
	
	return 0;
}

int			do_create_shell(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	if ( !model ) return 0;
	
	int				n(0), number(0), twod(0);
	double			radius(0), distance(0);
	
	if ( objc > 2 ) Tcl_GetIntFromObj(NULL, objv[2], &number);
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &radius);
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &distance);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &twod);

	if ( number > 0 ) n++;
	if ( radius > 0 ) n++;
	if ( distance > 0 ) n++;
	
	Bmodel*			numod = NULL;
	Bcomponent*		comp = NULL;
	
	if ( n > 1 ) {
		if ( twod )
			numod = model_create_circle(radius, 0, distance);
		else
			numod = model_create_shell(number, radius, distance);
		if ( model->comp ) {
			model->add(numod);
		} else {
			for ( comp = numod->comp; comp; comp = comp->next )
				model->add_component(comp);
			model_kill(numod);
		}
	}	
	
	return 0;
}

Tcl_Obj*	draw_components(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
/*
	Bimage*				p = NULL;

	Bstring				filename;
	Bstring				b = base(model->mapfile());
	for ( p = imglist; p; p = p->next ) {
		filename = p->file_name();
		if ( filename.base() == b ) break;
	}
*/
	Bimage*				p = imglist->find(model->mapfile());

	if ( !p ) return 0;

	int				slice(0), img_num(0);
	double			scale(1);
	
	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &scale);
	if ( objc > 3 ) Tcl_GetIntFromObj(NULL, objv[3], &slice);
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &img_num);

	if ( img_num >= p->images() ) img_num = 0;
	if ( slice >= p->sizeZ() ) slice = 0;
	if ( scale > 10 ) scale = 1;
	
	double				rad, dz;
	Vector3<long>		coor;
	char				string[MAXLINELEN] = "";
	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp; comp = comp->next ) {
		coor = p->image[img_num].image_coordinates(comp->location());
		rad = comp->radius()/p->image[img_num].sampling()[0];
		dz = fabs(coor[2] - slice);
		if ( dz <= rad ) {
			rad = sqrt(rad*rad - dz*dz);
			sprintf(string, "%ld %ld %g ", coor[0], coor[1], rad);
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}

	return returnObj;
}


/*
	Component type interface syntax:
		Bmodel comptype [item] [action] [params]
	Actions:
		count [n]						counts the number of component types
		countselected [n]				counts the number of selected component types
		ids [n]							returns the identifiers of component types
		list							returns a list of component types
		create [tid] [m] [fn] [in] [mn]	creates or updates a component type
		select [t]						sets the selection for a component type
		deselect [t]					unsets the selection for a component type
		firstselected					returns first selected component type
	Parameters:
		n			image number
		tid			type id
		mn			model number
		m			mass
		fn			file name
		in			image number
		t			type id list for selection/deselection
*/
Tcl_Obj*	do_component_type(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_component_type: action: " << action << " (" << action.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	Bmodel*				mp = model;

	if ( item.contains("Model") )
		for ( mp = model; mp && mp->identifier() != id.str(); mp = mp->next ) ;
	
	if ( model && !mp ) {
		cerr << "Error in do_component_type: item \"" << item << "\" not supported!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	if ( action == "count" ) {
		returnObj = comptype_count(mp, objc, objv);
	} else if ( action == "countselected" ) {
		returnObj = comptype_count_selected(mp, objc, objv);
	} else if ( action == "ids" ) {
		returnObj = comptype_ids(mp, objc, objv);
	} else if ( action == "list" ) {
		returnObj = comptype_list(mp);
	} else if ( action == "create" ) {
		returnObj = comptype_create_update(mp, objc, objv);
	} else if ( action == "select" ) {
		returnObj = comptype_select(mp, objc, objv);
	} else if ( action == "deselect" ) {
		returnObj = comptype_deselect(mp, objc, objv);
	} else if ( action == "firstselected" ) {
		returnObj = comptype_first_selected(mp, objc, objv);
	} else {
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}
	
	action = 0;

	return returnObj;
}

Tcl_Obj*	comptype_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nct = model->component_type_count();

	Tcl_SetIntObj(returnObj, nct);
	
	return returnObj;
}

Tcl_Obj*	comptype_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					nct(0);

	Bcomptype*			ct;

	for ( ct = model->type; ct; ct = ct->next ) if ( ct->select() ) nct++;
	
	Tcl_SetIntObj(returnObj, nct);
	
	return returnObj;
}

Tcl_Obj*	comptype_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN] = "";

	Bcomptype*			ct;

	for ( ct = model->type; ct; ct = ct->next ) if ( ct->select() ) {
		sprintf(string, " %ld", stol(ct->identifier()));
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	comptype_list(Bmodel* model)
{
	Tcl_Obj*		returnObj = Tcl_NewObj();

	char			string[MAXLINELEN];
	Bcomptype*		ct;
	
	model->update_type_counts();

	for ( ct = model->type; ct; ct = ct->next ) {
		sprintf(string, " %s %s %ld %g %g %g %ld", ct->identifier().c_str(),
			ct->file_name().c_str(), ct->image_number(),
			ct->component_count(), ct->mass(), ct->FOM(), ct->select());
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	comptype_create_update(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	int					num(-1);
	double				mass(0);
	Bstring				type_id("VER");
	Bcomptype*			ct;
	Bstring				filename;

	if ( objc > 4 ) type_id = Tcl_GetStringFromObj(objv[4], NULL);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &mass);
	if ( objc > 6 ) filename = Tcl_GetStringFromObj(objv[6], NULL);
	if ( objc > 7 ) Tcl_GetIntFromObj(NULL, objv[7], &num);

	ct = model->add_type(type_id);
	if ( mass > 0 ) ct->mass(mass);
	if ( filename.length() && filename != "?" ) ct->file_name(filename.str());
	if ( num >= 0 ) ct->image_number(num);
	
	filename = 0;
	
	return returnObj;
}

Tcl_Obj*	comptype_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	Bstring				selection_string;
	Bstring*			typelist, *type;
	Bcomptype*			ct = NULL;

	if ( objc > 4 ) selection_string = Tcl_GetStringFromObj(objv[4], NULL);

	if ( selection_string == "all" ) {
		for ( ct = model->type; ct; ct = ct->next ) ct->select(1);
	} else {
		typelist = selection_string.split(",");
		for ( ct = model->type; ct; ct = ct->next ) {
			for ( type = typelist; type; type = type->next )
				if ( ct->identifier() == type->str() ) ct->select(1);
		}
	}
	
	selection_string = 0;
	
	return returnObj;
}

Tcl_Obj*	comptype_deselect(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	Bstring				selection_string;
	Bstring*			typelist, *type;
	Bcomptype*			ct = NULL;

	if ( objc > 4 ) selection_string = Tcl_GetStringFromObj(objv[4], NULL);

	if ( selection_string == "all" ) {
		for ( ct = model->type; ct; ct = ct->next ) ct->select(1);
	} else {
		typelist = selection_string.split(",");
		for ( ct = model->type; ct; ct = ct->next ) {
			for ( type = typelist; type; type = type->next )
				if ( ct->identifier() == type->str() ) ct->select(0);
		}
	}
	
	selection_string = 0;
	
	return returnObj;
}

Tcl_Obj*	comptype_first_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	char				string[MAXLINELEN];
	Bcomptype*			ct = NULL;
	Bcomptype*			ctsel = NULL;

	for ( ct = model->type; ct && !ctsel; ct = ct->next )
		if ( ct->select() ) ctsel = ct;

	if ( ct ) {
		sprintf(string, " %s %s %ld %g %g %g %ld", ct->identifier().c_str(),
			ct->file_name().c_str(), ct->image_number(),
			ct->component_count(), ct->mass(), ct->FOM(), ct->select());
	} else {
		sprintf(string, " VER ? 0 0 1 1 1");
	}

	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}


/*
	Component interface syntax:
		Bmodel component [item] [action] [params]
	Actions:
		count [n]						counts the number of components
		countselected [n]				counts the number of selected components
		ids [n] [f]						returns the identifiers of components
		img_coords [n] [id]				returns one component's coordinates in the image
		img_coords [n] [f]				returns a list of component coordinates in the image
		location [id] [n]				returns the location of a component
		fom [id] [n]					returns the FOM of a component
		radius [id] [n]					returns the radius of a component
		select [x] [y] [z] [n]			returns the component identifier at the given location
		move [id] [n] [dx] [dy] [dz]	updates a component location by the given vector
		create [x] [y] [z] [n]			creates a new component at the given location
		delete [id] [n]					deletes a component
		color [id [n]					colors a component
	Parameters:
		id			component identifier
		x			x location
		y			y location
		z			z location
		n			image number
		f			FOM cutoff
		dx			change in x
		dy			change in y
		dz			change in z
*/
Tcl_Obj*	do_component(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_component: action: " << action << " (" << action.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	Bmodel*				mp = model;

	if ( item.contains("Model") )
		for ( mp = model; mp && mp->identifier() != id.str(); mp = mp->next ) ;
	
	if ( model && !mp ) {
		cerr << "Error in do_component: item \"" << item << "\" not supported!" << endl;
		cerr << "ID = " << id << endl;
//		for ( mp = model; mp; mp = mp->next ) cerr << mp->identifier() << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}
	
//	cout << mp->identifier() << tab << action << endl;

	if ( action == "count" ) {
		returnObj = component_count(mp, objc, objv);
	} else if ( action == "countselected" ) {
		returnObj = component_count_selected(mp, objc, objv);
	} else if ( action == "ids" ) {
		returnObj = component_ids(mp, objc, objv);
	} else if ( action == "location" ) {
		returnObj = component_location(mp, objc, objv);
	} else if ( action == "distance" ) {
		returnObj = component_distance(mp, objc, objv);
	} else if ( action == "img_coords" ) {
		returnObj = component_img_coords(mp, objc, objv);
	} else if ( action == "fom" ) {
		returnObj = component_fom(mp, objc, objv);
	} else if ( action == "radius" ) {
		returnObj = component_radius(mp, objc, objv);
	} else if ( action == "select" ) {
		returnObj = component_select(mp, objc, objv);
	} else if ( action == "move" ) {
		returnObj = component_move(mp, objc, objv);
	} else if ( action == "create" ) {
		returnObj = component_create(mp, objc, objv);
	} else if ( action == "delete" ) {
		returnObj = component_delete(mp, objc, objv);
	} else if ( action == "color" ) {
		returnObj = component_color(mp, objc, objv);
	} else {
		returnObj = Tcl_NewObj();
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	return returnObj;
}

Tcl_Obj*	component_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					ncomp = model->component_count();

	Tcl_SetIntObj(returnObj, ncomp);
	
	return returnObj;
}

Tcl_Obj*	component_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					ncomp(0);

	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) ncomp++;
	
	Tcl_SetIntObj(returnObj, ncomp);
	
	return returnObj;
}

Tcl_Obj*	component_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	double				fom_cut(0);
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->FOM()  >= fom_cut ) {
		sprintf(string, " %ld", stol(comp->identifier()));
		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	component_location(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id;
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);

	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
	
	if ( comp ) {
		sprintf(string, "%f %f %f", comp->location()[0], comp->location()[1], comp->location()[2]);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	component_distance(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id1, id2;
	double				dist(0);

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id1);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id2);

	Bcomponent*			comp = NULL;
	Bcomponent*			comp1 = NULL;
	Bcomponent*			comp2 = NULL;

	for ( comp = model->comp; comp && ( !comp1 || !comp2 ); comp = comp->next ) {
		if ( stoi(comp->identifier()) == id1 ) comp1 = comp;
		if ( stoi(comp->identifier()) == id2 ) comp2 = comp;
	}
	
	if ( comp1 && comp2 ) dist = comp1->location().distance(comp2->location());
	
	Tcl_SetDoubleObj(returnObj, dist);
	
	return returnObj;
}

Tcl_Obj*	component_img_coords(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;

	Bimage*				p = imglist->find(model->mapfile());

	if ( !p ) {
		cerr << "Error: Image " << model->mapfile() << " not found!" << endl;
		return returnObj;
	}
	
	int					id(0), n(model->image_number());
	double				fom_cut(0);
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	if ( fom_cut >= 1 ) id = (int) fom_cut;
	
	Bcomponent*			comp = NULL;
	Bstring				hexcol;

	if ( id > 0 ) {
		for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
		if ( comp ) {
			hexcol = comp->color().hex();
			sprintf(string, " %ld %g %g %g %g %s", stol(comp->identifier()),
				comp->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				comp->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				comp->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				comp->radius()/p->sampling(0)[0], hexcol.c_str());
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	} else {
		for ( comp = model->comp; comp; comp = comp->next ) if ( comp->FOM() >= fom_cut ) {
			hexcol = comp->color().hex();
			sprintf(string, " %ld %g %g %g %g %s", stol(comp->identifier()),
				comp->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				comp->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				comp->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				comp->radius()/p->sampling(0)[0], hexcol.c_str());
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
//	} else {
//		sprintf(string, " -1");
//		Tcl_AppendToObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	component_fom(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id, set(0);
	double				fom = 1;

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) {
		Tcl_GetDoubleFromObj(NULL, objv[5], &fom);
		set = 1;
	}

	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
	
	if ( comp ) {
		if ( set ) comp->FOM(fom);
		else fom = comp->FOM() ;
	}
	
	Tcl_SetDoubleObj(returnObj, fom);
	
	return returnObj;
}

Tcl_Obj*	component_radius(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id(-1), set(0);
	double				rad(1);

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) {
		Tcl_GetDoubleFromObj(NULL, objv[5], &rad);
		set = 1;
	}

//	Bcomponent*			comp = NULL;

//	for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
	
	string				s = to_string(id);
	Bcomponent*			comp = model->find_component(s);
	
	if ( comp ) {
		if ( set ) comp->radius(rad);
		else rad = comp->radius();
	}
	
	Tcl_SetDoubleObj(returnObj, rad);
	
	return returnObj;
}

Tcl_Obj*	component_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id(0);
	double				x(0), y(0), z(0);
	Bstring				selection_string, types, ids;
	double				d, dmin(1000000);
	Vector3<float>		loc;
	Bcomponent*			comp = NULL;
	Bcomponent*			compsel = NULL;

	if ( objc > 4 ) selection_string = Tcl_GetStringFromObj(objv[4], NULL);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG component_select: Model " << model->identifier() << ": selection_string=" << selection_string << endl;

	if ( selection_string.length() < 1 ) {
		return returnObj;
	} else if ( selection_string == "none" ) {
		model_unset_selection(model);
	} else if ( selection_string == "all" ) {
		model_reset_selection(model);
	} else if ( selection_string == "types" ) {
		if ( objc > 5 ) {
			types = "%.";
			model_select(model, types);
			types = "%";
			types += Tcl_GetStringFromObj(objv[5], NULL);
			types = types.replace(' ', ',');
//			cout << "types = " << types << endl;
			model_select(model, types);
		}
	} else if ( selection_string == "id" ) {
		if ( objc > 5 ) {
			model_unset_selection(model);
			ids = "@";
			ids += Tcl_GetStringFromObj(objv[5], NULL);
			ids = ids.replace(' ', ',');
//			cout << "ids = " << ids << endl;
			model_select(model, ids);
		}
	} else {
		if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
		if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
		if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
		loc = Vector3<float>(x, y, z);
//		cout << "location for selection = " << loc[0] << " " << loc[1] << " " << loc[2] << endl;
		for ( comp = model->comp; comp; comp = comp->next ) {
			d = loc.distance(comp->location());
			if ( dmin > d ) {
				dmin = d;
				compsel = comp;
			}
		}
		if ( compsel && ( dmin <= compsel->radius() ) ) {
			id = stoi(compsel->identifier());
			if ( compsel->select() ) compsel->select(0);	// Toggle the selection
			else compsel->select(1);
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	component_move(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id(0);
	double				dx(0), dy(0), dz(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &dx);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &dy);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &dz);
	
	Bcomponent*			comp = NULL;
	Vector3<float>		d(dx, dy, dz);

	for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
	
	if ( comp ) comp->location(comp->location() + d);
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	component_create(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					i(0);
	double				x(0), y(0), z(0);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	
//	Bstring				id("0");
	string				tid("VER");
	Bcomponent*			comp = NULL;
	Bcomptype*			ct = NULL;
	
	for ( i=0, comp = model->comp; comp; comp = comp->next )
		if ( i < stoi(comp->identifier()) ) i = stoi(comp->identifier());
	
//	id = Bstring(++i, "%d");
//	comp = component_add(&model->comp, id);
	if ( model->comp ) comp = model->comp->add(++i);
	else model->comp = comp = new Bcomponent(++i);
	for ( ct = model->type; ct && ct->select() < 1; ct = ct->next ) ;
	if ( ct ) comp->type(ct);
//	else comp->type = model_add_type_by_id(model, tid);
	else comp->type(model->add_type(tid));
	comp->location(Vector3<float>(x, y, z));
	comp->FOM(1);
	comp->select(1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG component_create: id=" << i << endl;
	
	Tcl_SetIntObj(returnObj, i);
	
	return returnObj;
}

Tcl_Obj*	component_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id(0);

	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id);
	
	Bcomponent*			comp = NULL;

	for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
	
	if ( comp ) {
		comp_associated_links_kill(comp, &model->link);
		comp->type(NULL);
		remove_item((char **)&model->comp, (char *)comp, sizeof(Bcomponent));
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}

Tcl_Obj*	component_color(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id(0);
	RGBA<float>			color;
	Bstring				selection_string, color_string;

	if ( objc > 4 ) selection_string = Tcl_GetStringFromObj(objv[4], NULL);
	if ( objc > 5 ) {
		color_string = Tcl_GetStringFromObj(objv[5], NULL);
		color = RGBA<float>(color_string);
	}
	
	Bcomponent*			comp = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG component_color: Select=" << selection_string << " color=" << color_string << endl;

	if ( selection_string.length() < 1 || selection_string == "none" ) {
		return returnObj;
	} else if ( selection_string == "all" ) {
		model_color_uniformly(model, color);
	} else if ( selection_string == "selected" ) {
		model_color_selected_types(model, color);
	} else if ( selection_string == "density" ) {
		model_color_by_density(model);
	} else if ( selection_string == "fom" ) {
		model_color_by_fom(model);
	} else {
		color_string = Tcl_GetStringFromObj(objv[4], NULL);
		color = RGBA<float>(color_string);
		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id);
		if ( id > 0 ) {
			for ( comp = model->comp; comp && stoi(comp->identifier()) != id; comp = comp->next ) ;
			if ( comp ) comp->color(color);
		} else {
			for ( comp = model->comp; comp; comp = comp->next )
				comp->color(color);
		}
	}
	
	Tcl_SetIntObj(returnObj, id);
	
	return returnObj;
}



/*
	Link interface syntax:
		Bmodel link [action] [params]
	Actions:
		count [n]						counts the number of links
		countselected [n]				counts the number of selected links
		ids [n] [f]						returns the link component identifiers
		img_coords [n] [id1] [id2]		returns one pair of link coordinates in the image
		img_coords [n] [f]				returns a list of link coordinates in the image
		select [x] [y] [z] [n] [r]		returns the link component identifiers at the given location
		create [id1] [id2] [n]			creates a new link between two components
		delete [id1] [id2] [n]			deletes a link between two components
	Parameters:
		id1			first link component identifier
		id2			second link component identifier
		x			x location
		y			y location
		z			z location
		n			image number
		f			FOM cutoff
		r			radius to select a link
*/
Tcl_Obj*	do_link(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = NULL;

	Bstring				item = Tcl_GetStringFromObj(objv[2], NULL);

	Bstring				action = Tcl_GetStringFromObj(objv[3], NULL);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG do_link: action: " << action << " (" << action.length() << ")" << endl;

	Bstring				id = item.post(':');
	id = id.remove(' ');
	item = item.pre(':');
	
	Bmodel*				mp = model;

	if ( item.contains("Model") )
		for ( mp = model; mp && mp->identifier() != id.str(); mp = mp->next ) ;
	
	if ( model && !mp ) {
		cerr << "Error in do_link: item \"" << item << "\" not supported!" << endl;
		returnObj = Tcl_NewObj();
		return returnObj;
	}

	if ( action == "count" ) {
		returnObj = link_count(mp, objc, objv);
	} else if ( action == "countselected" ) {
		returnObj = link_count_selected(mp, objc, objv);
	} else if ( action == "ids" ) {
		returnObj = link_ids(mp, objc, objv);
	} else if ( action == "img_coords" ) {
		returnObj = link_img_coords(mp, objc, objv);
	} else if ( action == "select" ) {
		returnObj = link_select(mp, objc, objv);
	} else if ( action == "create" ) {
		returnObj = link_create(mp, objc, objv);
	} else if ( action == "delete" ) {
		returnObj = link_delete(mp, objc, objv);
	} else {
		returnObj = Tcl_NewObj();
		cerr << "Error: Action " << action << " not recognized!" << endl;
	}

	return returnObj;
}

Tcl_Obj*	link_count(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					nlink(0);

	Blink*				link = NULL;

	for ( link = model->link; link; link = link->next ) nlink++;
	
	Tcl_SetIntObj(returnObj, nlink);
	
	return returnObj;
}

Tcl_Obj*	link_count_selected(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					nlink(0);

	Blink*				link = NULL;

	for ( link = model->link; link; link = link->next ) if ( link->select() ) nlink++;
	
	Tcl_SetIntObj(returnObj, nlink);
	
	return returnObj;
}

Tcl_Obj*	link_ids(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	double				fom_cut(0);
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	Blink*				link = NULL;

	for ( link = model->link; link; link = link->next )
		if ( link->comp[0]->FOM()  >= fom_cut && link->comp[1]->FOM()  >= fom_cut ) {
			sprintf(string, " %ld %ld", stol(link->comp[0]->identifier()), stol(link->comp[1]->identifier()));
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	
	return returnObj;
}

Tcl_Obj*	link_img_coords(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	Bimage*				p = imglist->find(model->mapfile());

	if ( !p ) {
		cerr << "Error: Image " << model->mapfile() << " not found!" << endl;
		return returnObj;
	}


	int					id1(0), id2(0), n(model->image_number());
	double				fom_cut(0);
	char				string[MAXLINELEN] = "";

	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &fom_cut);

	if ( fom_cut >= 1 ) {
		id1 = (int) fom_cut;
		if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id2);
	}
	
	Blink*				link = NULL;
	Bstring				hexcol;

	if ( id1 > 0 && id2 > 0 ) {
		for ( link = model->link; link; link = link->next )
			if ( ( stoi(link->comp[0]->identifier()) == id1 && stoi(link->comp[1]->identifier()) == id2 ) ||
				( stoi(link->comp[0]->identifier()) == id2 && stoi(link->comp[1]->identifier()) == id1 ) ) break;
		if ( link ) {
			hexcol = link->color().hex();
			sprintf(string, " %g %g %g %g %g %g %g %s",
				link->comp[0]->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				link->comp[0]->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				link->comp[0]->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				link->comp[1]->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				link->comp[1]->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				link->comp[1]->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				link->radius()/p->sampling(0)[0], hexcol.c_str());
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	} else {
		for ( link = model->link; link; link = link->next )
				if ( link->comp[0]->FOM()  >= fom_cut && link->comp[1]->FOM()  >= fom_cut ) {
			hexcol = link->color().hex();
			sprintf(string, " %g %g %g %g %g %g %g %s",
				link->comp[0]->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				link->comp[0]->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				link->comp[0]->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				link->comp[1]->location()[0]/p->sampling(0)[0] + p->image[n].origin()[0],
				link->comp[1]->location()[1]/p->sampling(0)[1] + p->image[n].origin()[1],
				link->comp[1]->location()[2]/p->sampling(0)[2] + p->image[n].origin()[2],
				link->radius()/p->sampling(0)[0], hexcol.c_str());
			Tcl_AppendToObj(returnObj, string, strlen(string));
		}
	}
	
	return returnObj;
}

Tcl_Obj*	link_select(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id1(0), id2(0);
	double				x(0), y(0), z(0), rad(20);
	
	if ( objc > 4 ) Tcl_GetDoubleFromObj(NULL, objv[4], &x);
	if ( objc > 5 ) Tcl_GetDoubleFromObj(NULL, objv[5], &y);
	if ( objc > 6 ) Tcl_GetDoubleFromObj(NULL, objv[6], &z);
	if ( objc > 7 ) Tcl_GetDoubleFromObj(NULL, objv[7], &rad);
	
	char				string[64] = "";
	double				d, d1, d2, dmin(1000000);
	Vector3<float>		loc(x, y, z);
	Blink*				link = NULL;
	Blink*				linksel = NULL;

	for ( link = model->link; link; link = link->next ) {
		d = link->comp[0]->location().distance(link->comp[1]->location());
		d1 = link->comp[0]->location().distance(loc);
		d2 = link->comp[1]->location().distance(loc);
		d = d1 + d2 - d;
		if ( dmin > d ) {
			dmin = d;
			linksel = link;
		}
	}
	
	if ( linksel && dmin <= rad ) {
		id1 = stoi(linksel->comp[0]->identifier());
		id2 = stoi(linksel->comp[1]->identifier());
		sprintf(string, "%d %d", id1, id2);
		Tcl_SetStringObj(returnObj, string, strlen(string));
	}
	
	return returnObj;
}

Tcl_Obj*	link_create(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id1(0), id2(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id1);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id2);
	
	char				string[64] = "";
	Bcomponent*			comp1 = NULL;
	Bcomponent*			comp2 = NULL;

	if ( id1 != id2 ) {
		for ( comp1 = model->comp; comp1 && stoi(comp1->identifier()) != id1; comp1 = comp1->next ) ;
		for ( comp2 = model->comp; comp2 && stoi(comp2->identifier()) != id2; comp2 = comp2->next ) ;
		if ( comp1 && comp2 )
			link_add(&model->link, comp1, comp2, 0, 1);
	}
	
	sprintf(string, "%d %d", id1, id2);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}

Tcl_Obj*	link_delete(Bmodel* model, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_Obj*			returnObj = Tcl_NewObj();

	if ( !model ) return returnObj;
	
	int					id1(0), id2(0);
	
	if ( objc > 4 ) Tcl_GetIntFromObj(NULL, objv[4], &id1);
	if ( objc > 5 ) Tcl_GetIntFromObj(NULL, objv[5], &id2);
	
	char				string[64] = "";
	Blink*				link = NULL;
	Bcomponent*			comp1 = NULL;
	Bcomponent*			comp2 = NULL;

	if ( id1 != id2 ) {
		for ( comp1 = model->comp; comp1 && stoi(comp1->identifier()) != id1; comp1 = comp1->next ) ;
		for ( comp2 = model->comp; comp2 && stoi(comp2->identifier()) != id2; comp2 = comp2->next ) ;
		if ( comp1 && comp2 )
			for ( link = model->link; link && link->comp[0] != comp1 && link->comp[1] != comp2; link = link->next ) ;
		if ( link )
			remove_item((char **)&model->link, (char *)link, sizeof(Blink));
	}
	
	sprintf(string, "%d %d", id1, id2);
	Tcl_SetStringObj(returnObj, string, strlen(string));
	
	return returnObj;
}

