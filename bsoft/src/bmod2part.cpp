/**
@file	bmod2part.cpp
@brief	Converting from a pseudo-atomic model to particle parameters
@author Bernard Heymann
@date	Created: 20070423
@date 	Modified: 20220223
**/

#include "rwmodel.h"
#include "model_select.h"
#include "model_util.h"
#include "rwmg.h"
#include "mg_processing.h"
#include "rwimg.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bproject*	project_part_from_model(Bmodel* model, Vector3<double> origin, 
				Vector3<double> sampling, Vector3<long> box_size);
Bmodel*		models_from_particles(Bproject* project, Bstring modtype);
Bmodel*		components_from_2D_particles(Bproject* project, Bstring comptype);
Bmodel*		components_from_3D_particles(Bproject* project, Bstring comptype);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmod2part [options] input.star",
"-------------------------------------",
"Converting between model components and reconstruction particle parameters.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-all                     Reset selection to all components before other selections.",
"-type comp               Type of conversion: part->models or part->components.",
"-settype VER             Set model or component types.",
"-associate TRS,trs.pdb   Associate a component type with a file name.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-sampling 2.3,2.3,1      Sampling (angstrom/voxel, one value sets all).",
"-box 100,83,120          Particle box size (default 50,50,50; one sets all).",
" ",
"Input:",
"-map image.spi           Input 3D map file.",
"-parameters param.star   Input atomic parameter file.",
" ",
"Output:",
"-output part.star        Output particle parameter or model file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int				use_rec(0);					// Flag to process reconstructions
	int 			reset(0);					// Keep selection as read from file
	Bstring			type;						// Type of conversion
	Bstring			set_type("PRT");			// Component type to set
	Bstring			associate_type;				// Component type
	Bstring			associate_file;				// Component type file name
	Vector3<double>	origin;						// Coordinate origin placement
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam;    					// Sampling in angstrom/voxel side
	Vector3<long>	box_size(50,50,50);			// Particle box size
	Bstring			mapfile;					// Input map file
	Bstring			paramfile;					// Input parameter file
	Bstring			outfile;					// Output parameter file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "type" ) type = curropt->value.lower();
		if ( curropt->tag == "settype" )
			set_type = curropt->value;
		if ( curropt->tag == "associate" ) {
			associate_type = curropt->value.pre(',');
			associate_file = curropt->value.post(',');
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "box" )
			box_size = curropt->size();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No model files specified!" << endl;
		bexit(-1);
	}

	Bmodel*			model = NULL;
	Bproject*		project = NULL;
	Bimage*			p = NULL;

	Vector3<double>		ori3, sam3;
	
	if ( file_type(*file_list) == Model ) {
		model = read_model(file_list, paramfile);
		if ( !model ) {
			cerr << "Error: Input file not read!" << endl;
			bexit(-1);
		}
		if ( reset ) models_process(model, model_reset_selection);
		if ( mapfile.length() ) model->mapfile(mapfile.str());
		if ( model->mapfile().length() ) {
			p = read_img(model->mapfile(), 0, -1);
			if ( p ) {
				if ( sam.volume() <= 0 ) sam = p->sampling(0);
				if ( origin.length() <= 0 ) origin = p->image->origin();
				else if ( set_origin == 2 ) origin = p->default_origin();
				delete p;
			}
		}
		ori3 = Vector3<double>(origin[0], origin[1], origin[2]);
		sam3 = Vector3<double>(sam[0], sam[1], sam[2]);
		project = project_part_from_model(model, ori3, sam3, box_size);
		if ( project && outfile.length() ) {
				write_project(outfile, project, 0, 0);
		}
	} else {
		project = read_project(file_list);
		if ( use_rec ) project->select = 1;
		Vector3<double>		ori3(origin[0], origin[1], origin[2]);
		if ( project->select && project->rec ) {
			if ( mapfile.length() ) project->rec->frec = mapfile;
			if ( project->rec->frec.length() ) {
				p = read_img(project->rec->frec, 0, -1);
				if ( p ) {
					if ( sam.volume() <= 0 ) sam = p->sampling(0);
					if ( origin.length() <= 0 ) origin = p->image->origin();
					else if ( set_origin == 2 ) origin = p->default_origin();
					delete p;
				}
			}
			ori3 = Vector3<double>(origin[0], origin[1], origin[2]);
			if ( origin.length() > 0 ) project->rec->origin = ori3;
			if ( sam[0] > 0 ) project_set_mg_pixel_size(project, sam);
			if ( type.length() && type[0] == 'm' ) model = models_from_particles(project, set_type);
			else model = components_from_3D_particles(project, set_type);
		} else if ( project->field && project->field->mg ) {
			ori3 = Vector3<double>(origin[0], origin[1], origin[2]);
			if ( ori3.volume() > 0 ) project_set_micrograph_origins(project, ori3);
			model = components_from_2D_particles(project, set_type);
			if ( mapfile.length() ) model->mapfile(mapfile.str());
		}
		if ( associate_file.length() )
			model_associate(model, associate_type, associate_file);
		if ( outfile.length() && model ) {
				write_model(outfile, model);
		}
	}
	
	string_kill(file_list);
	model_kill(model);
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

/**
@brief 	Converts component coordinates in a model to reconstruction particle coordinates.  
@param 	*model		model.
@param 	origin		map origin in voxel coordinates.
@param 	sampling	voxel size.
@param 	box_size	particle box size.
@return Bproject*	new project.
**/
Bproject*	project_part_from_model(Bmodel* model, Vector3<double> origin, 
				Vector3<double> sampling, Vector3<long> box_size)
{
	long				i;

	for ( i=0; i<3; i++ ) if ( sampling[i] <= 0 ) sampling[i] = 1;

	Bstring				b = base(model->mapfile());
	Bproject*			project = new Bproject;
	Breconstruction*	rec = reconstruction_add(&project->rec, b);
	Bparticle*			part = NULL;

	rec->frec = model->mapfile();
	rec->voxel_size = sampling;
	rec->origin = origin;
	rec->box_size = box_size;
	rec->bad_radius = box_size[0]/4;
	
	if ( verbose )
		cout << "Generating particles from model components" << endl << endl;
	
	Bmodel*			mp;
	Bcomponent*		comp;
	View2<float>	v;
	
	for ( i=0, mp = model; mp; mp = mp->next ) if ( mp->select() ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			part = particle_add(&part, ++i);
			if ( !rec->part ) rec->part = part;
			part->ori = box_size/2;
			part->loc = comp->location()/sampling + origin;
			v = comp->view().backward();
			part->view = View(v[0],v[1],v[2],v[3]);
			part->sel = comp->select();
			part->fom[0] = comp->FOM();
		}
	}
	
	return project;
}

/**
@brief 	Converts particle coordinates for a reconstruction into model component coordinates.  
@param 	*project	project with reconstruction.
@param	modtype		model tyoe string.
@return Bmodel*		new project.
**/
Bmodel*		models_from_particles(Bproject* project, Bstring modtype)
{
	long				nmod(0);
	Breconstruction*	rec = project->rec;
	Bparticle*			part = NULL;
	Bmodel*				model = NULL;
	Bmodel*				mp = NULL;
	Bstring				id;
	
	if ( verbose )
		cout << "Generating models from particles" << endl;
	
	for ( part = rec->part; part; part = part->next ) {
		id = Bstring(part->id, "%d");
//		mp = model_add(&mp, id);
//		if ( !model ) model = mp;
		if ( mp ) mp = mp->add(part->id);
		else model = mp = new Bmodel(part->id);
		mp->model_type(modtype.str());
		if ( part->fpart.length() ) {
			mp->mapfile(part->fpart.str());
			mp->image_number(0);
		} else if ( rec->fpart.length() ) {
			mp->mapfile(rec->fpart.str());
			mp->image_number(part->id-1);
		}
		mp->FOM(part->fom[0]);
		mp->select(part->sel);
//		comp->loc = (part->loc - rec->origin)*rec->voxel_size;
//		comp->view = part->view.backward();
		nmod++;
	}

	if ( verbose )
		cout << "Models generated:               " << nmod << endl << endl;
	
	return model;
}

/**
@brief 	Converts particle coordinates from micrographs into model component coordinates.  
@param 	*project	project with micrographs.
@param	comptype	component type string.
@return Bmodel*		new project.
**/
Bmodel*		components_from_2D_particles(Bproject* project, Bstring comptype)
{
	long				i;
	Vector3<double>		ori;
	Bfield*				field;
	Bmicrograph*		mg;
	Bparticle*			part = NULL;
	Bstring				id("Particles");
	Bmodel*				model = new Bmodel(id);
	Bcomponent*			comp = NULL;
	Bcomptype*			ct = model->add_type(comptype);

	if ( verbose )
		cout << "Generating model components from 2D particles" << endl;
	
	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			ori = mg->box_size/2.0;
			cout << "box_size = " << mg->box_size << endl;
			cout << "ori = " << ori << endl;
			for ( part = mg->part; part; part = part->next ) {
//				id = Bstring(++i, "%d");
//				comp = component_add(&comp, id);
//				if ( !model->comp ) model->comp = comp;
				if ( comp ) comp = comp->add(++i);
				else model->comp = comp = new Bcomponent(++i);
				comp->type(ct);
				if ( part->loc.length() )
					comp->location((part->loc - mg->origin)*mg->pixel_size);
				else
					comp->location((part->ori - ori)*part->pixel_size);
				View v(part->view.backward());
				comp->view(View2<float>(v[0],v[1],v[2],v[3]));
				comp->select(part->sel);
				comp->FOM(part->fom[0]);
			}
		}
	}

	if ( verbose )
		cout << "Components generated:           " << i << endl << endl;
	
	return model;
}

/**
@brief 	Converts particle coordinates for a reconstruction into model component coordinates.  
@param 	*project	project with reconstruction.
@param	comptype	component type string.
@return Bmodel*		new project.
**/
Bmodel*		components_from_3D_particles(Bproject* project, Bstring comptype)
{
	long				ncomp(0);
	Vector3<double>		ori;
	Breconstruction*	rec = project->rec;
	Bparticle*			part = NULL;
	Bmodel*				model = new Bmodel(rec->id);
	Bcomponent*			comp = NULL;
	Bstring				id;
	Bcomptype*			ct = model->add_type(comptype);

	model->mapfile(rec->frec.str());
	
	if ( verbose )
		cout << "Generating model components from 3D particles" << endl;
	
	ori = rec->box_size/2;
	
	for ( part = rec->part; part; part = part->next ) {
//		id = Bstring(part->id, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(part->id);
		else model->comp = comp = new Bcomponent(part->id);
		comp->type(ct);
		if ( part->loc.length() )
			comp->location((part->loc - rec->origin)*rec->voxel_size);
		else
			comp->location((part->ori - ori)*part->pixel_size);
		View v(part->view.backward());
		comp->view(View2<float>(v[0],v[1],v[2],v[3]));
		comp->select(part->sel);
		comp->FOM(part->fom[0]);
		ncomp++;
	}

	if ( verbose )
		cout << "Components generated:           " << ncomp << endl << endl;
	
	return model;
}

