/**
@file	bdelta.cpp
@brief	Generating a deltagraph with minimal input parameters
@author Bernard Heymann
@date	Created: 20080103
@date 	Modified: 20180327
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_poly_delta.h"
#include "model_transform.h"
#include "model_links.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bdelta [options] in.star [in2.star ...]",
"----------------------------------------------",
"Generates and analyzes deltagraphs.",
"The creation options are mutually exclusive.",
" ",
"Actions:",
"-pentagonal              Create a pentagonal deltagraph.",
"-hexagonal               Create a hexagonal deltagraph.",
"-tube 8,3                Create a tubular deltagraph with the given lattice parameters.",
"-caps 2                  Add 1 or 2 caps.",
"-dual                    Calculate the dual of a polyhedron.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
"-linklength 150          Link length (default 1).",
"-radius 3                Radius in integer units (default 1).",
"-height 5                Height in integer units (default 1).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int				type(0);				// Type of deltagraph
	int				h(0), k(0);				// Tube lattice
	int				radius(0);				// Cylinder radius
	int				height(0);				// Cylinder height
	int				dual(0);				// Flag to calculate dual of polyhedron
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	double			linklength(1);			// Link length
    Bstring    		comp_select("all");
	Bstring			paramfile;				// Input parameter file
	Bstring			outfile;				// Output model parameter file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "pentagonal" ) type = 2;
		if ( curropt->tag == "hexagonal" ) type = 1;
		if ( curropt->tag == "caps" ) type |= 8 * curropt->value.integer() - 4;
		if ( curropt->tag == "tube" )
			if ( curropt->values(h, k) < 1 )
				cerr << "-tube: At least one integer must be specified!" << endl;
		if ( curropt->tag == "dual" ) dual = 1;
		if ( curropt->tag == "radius" )
			if ( ( radius = curropt->value.integer() ) < 1 )
				cerr << "-radius: A radius must be specified!" << endl;
		if ( curropt->tag == "height" )
			if ( ( height = curropt->value.integer() ) < 1 )
				cerr << "-height: A height must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 0.001 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 0.001 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 0.001 )
				cerr << "-linklength: The link length between vertices must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();

	Bmodel*			model = NULL;		
	Vector3<double>	scale(linklength, linklength, linklength), origin;
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		if ( radius > 0 && height > 0 && ( type & 3 ) ) {
			model = model_delta_create_cylinder(type, radius, height);
		} else if ( h ) {
			model = model_delta_create_tube(h, k, height);
		}
		if ( model ) {
			if ( linklength != 1 )
				model_scale(model, scale, origin);
			models_process(model, 0.1*linklength, model_set_link_radius);
		} else {
			cerr << "Error: No model files specified or generated!" << endl;
			bexit(-1);
		}
	} else {
		model = read_model(file_list, paramfile);		
		string_kill(file_list);
	}

	if ( !model ) {
		cerr << "Error: No model read or generated!" << endl;
		bexit(-1);
	}
	
	if ( compradius ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius ) models_process(model, linkradius, model_set_link_radius);

	model_poly_generate(model);

	if ( dual ) {
		Bmodel*		new_model = model_poly_dual(model, 0);
		if ( new_model ) {
			model_kill(model);
			model = new_model;
		}
	}
	
	if ( outfile.length() ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

