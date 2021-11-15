/**
@file	bpath.cpp
@brief	A program to analyze paths and cycles in a model.
@author Bernard Heymann
@date	Created: 20111109
@date	Modified: 20140325
**/

#include "rwmodel.h"
#include "model_path.h"
#include "model_select.h"
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
"Usage: bpath [options] in.star",
"------------------------------",
"Analyzes paths and cycles in a model.",
" ",
"Actions:",
"-all                     Reset selection to all components before other selections.",
"-cutoff 23.5             Assign links within this cutoff distance.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output path.star        Output model file with components linked along paths.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int 			reset(0);				// Keep selection as read from file
	double			cutoff(0);				// Cutoff distance to assign links
	int				minval(2);				// Minimum valency
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	Bstring			paramfile;				// Input parameter file name
	Bstring			outfile;				// Output parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "cutoff" )
			if ( ( cutoff = curropt->value.real() ) < 1 )
				cerr << "-cutoff: A cutoff distance must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 0.001 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 0.001 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
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

	Bmodel*		model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}

	if ( reset ) models_process(model, model_reset_selection);
	
	if ( cutoff ) models_process(model, cutoff, model_link_list_generate);
	else models_process(model, minval, model_links_minimum_valency);
	
	Bmodel*		modpath = NULL;
	modpath = model_hamiltonian_cycle(model);
	model_kill(model);
	model = modpath;
	
	if ( compradius ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius ) models_process(model, linkradius, model_set_link_radius);

	model_selection_stats(model);
	
	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


