/**
@file	bmodmap.cpp
@brief	A tool to genrate a map from a model.
@author Bernard Heymann
@date	Created: 20081112
@date 	Modified: 20110804
**/

#include "rwmodel.h"
#include "model_map.h"
#include "rwimg.h"
#include "model_select.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodmap [options] in.star",
"--------------------------------",
"Generates a map from a model.",
" ",
"Actions:",
"-all                     Reset selection to all components before other selections.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new image data type.",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-size 10,10,10           Map size, input map size otherwise (voxels).",
"-sampling 1              Sampling (default 1 angstrom/voxel).",
"-sigma 3.5               Density decay around a component (default 10).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output model.star       New model file.",
"-map density.map         New map file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int 			reset(0);						// Keep selection as read from file
	DataType 		nudatatype(Unknown_Type);		// Conversion to new type
	Vector3<double>	origin;				// Coordinate origin placement
	int				set_origin(0);					// Flag to set origin
	Vector3<long>	size;					// New map size
	Vector3<double>	sam;    				// Sampling in angstrom/voxel side
	double			sigma(0);						// Component density decay radius
	Bstring			paramfile;						// Input parameter file name
	Bstring			outfile;						// Output parameter file name
	Bstring			mapfile;						// Output map file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "size" )
        	size = curropt->size();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "sigma" )
			if ( ( sigma = curropt->value.real() ) < 1 )
				cerr << "-sigma: The component density decay radius must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
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

	Bmodel*		model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}

	if ( reset ) models_process(model, model_reset_selection);

	if ( set_origin == 2 ) origin = {double(size[0]/2), double(size[1]/2), double(size[2]/2)};

	Bimage*			map = NULL;
	
	if ( mapfile.length() && size.volume() ) {
		model_set_map_filenames(model, mapfile);
		map = img_from_model(model, origin, size, sam, sigma);
		if ( nudatatype == Unknown_Type ) nudatatype = map->data_type();
		map->change_type(nudatatype);
		write_img(mapfile, map, 0);
		delete map;
	}

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

