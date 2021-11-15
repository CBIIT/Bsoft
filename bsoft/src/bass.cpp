/**
@file	bass.cpp
@brief	Program to assemble molecular components.
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20091102
**/

#include "rwmolecule.h"
#include "mol_edit.h"
#include "rwmodel.h"
#include "model_select.h"
#include "model_views.h"
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
"Usage: bass [options] in1.star [in2.star ...]",
"---------------------------------------------",
"Assemble molecular components specified in a model.",
" ",
"Actions for preparation:",
"-com                     Generate a new model from the centers-of-mass of molecules.",
"-all                     Reset selection to all components before other selections.",
" ",
"Selections:",
"-select #232@14          Select models and components.",
"-first 22                Select the first number of components.",
" ",
"Actions:",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-associate TRS,trs.pdb   Associate a component type with a file name.",
"-untangle 3.5,0.2        Eliminate overlaps by moving molecules apart (sampling and damping factor).",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-separate                Each model is defined as a separate molecule group.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components before other selections.",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-coordinates all.pdb     Output coordinate files.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			all(0);						// Keep selection as read from file
	int 			reset(0);					// Keep selection as ouput
	Bstring			mod_select;					// Model and component selection
	int				first(0);					// First number of components to select
	int				com(0);						// Flag to generate a centers-of-mass model
	Bstring			calc_views;					// Mode to calculate component views
	Bstring			associate_type;				// Component type
	Bstring			associate_file;				// Component file name
	int				separate(0);				// Flag to define separate molecule groups
	double			untangle(0);				// Untangling grid sampling
	double			lambda(0.1);				// Untangling damping factor
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
    Bstring    		atom_select("all");
	Bstring			paramfile;					// Input parameter file name
	Bstring			outfile;					// Output parameter file name
	Bstring			coorfile;					// Output coordinates file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "first" )
			if ( ( first = curropt->value.integer() ) < 1 )
				cerr << "-first: An integer must be specified!" << endl;
		if ( curropt->tag == "com" ) com = 1;
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "associate" ) {
			associate_file = curropt->value;
			associate_type = associate_file.pre(',');
			associate_file = associate_file.post(',');
		}
		if ( curropt->tag == "untangle" )
			if ( curropt->values(untangle, lambda) < 1 )
				cerr << "-untangle: Grid sampling in angstrom must be specified!" << endl;
		if ( curropt->tag == "separate" ) separate = 1;
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 1 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->filename();
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
	Bmolgroup*		molgroup = NULL;
	
	if ( com ) {
		molgroup = read_molecule(*file_list, atom_select, paramfile);
		model = model_generate_com(molgroup);
		molgroup_kill(molgroup);
		molgroup = NULL;
	} else {
		model = read_model(file_list, paramfile);		
		string_kill(file_list);
	}

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( all ) models_process(model, model_reset_selection);

	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( first ) model_select_first(model, first);
	
	if ( associate_file.length() )
		model_associate(model, associate_type, associate_file);
	
	if ( comprad > 0 ) models_process(model, comprad, model_set_component_radius);

	if ( linkrad > 0 ) models_process(model, linkrad, model_set_link_radius);

	if ( calc_views.length() ) model_calculate_views(model, calc_views);
	
	molgroup = model_assemble(model, paramfile, separate);
	
	if ( molgroup ) {
		if ( untangle ) {
			molgroup_untangle_groups(molgroup, untangle, lambda);
			model_update_centers_of_mass(model, molgroup);
		}
		if ( coorfile.length() ) {
			molecule_update_comment(molgroup, argc, argv);
			molgroup_list_write(coorfile, molgroup);
		}
		molgroup_list_kill(molgroup);
	}
		
	if ( reset ) models_process(model, model_reset_selection);

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

