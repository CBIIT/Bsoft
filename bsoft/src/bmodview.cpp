/**
@file	bmodview.cpp
@brief	Manipulates models.
@author Bernard Heymann
@date	Created: 20081120
@date 	Modified: 20141029
**/

#include "rwmodel.h"
#include "model_views.h"
#include "model_util.h"
#include "model_transform.h"
#include "model_symmetry.h"
#include "model_select.h"
#include "model_links.h"
#include "ps_model.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodview [options] in1.star [in2.star...]",
"------------------------------------------------",
"Manipulates models.",
" ",
"Actions:",
"-merge                   Merge models in different files rather than concatenate.",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-invert                  Invert views for selected components.",
"-setasu D6               Set components to within an asymmetric unit.",
"-directions 10,1         Histogram of component view directions with histogram bin width (degrees)",
"                         and a flag: 0=wrt z-axis, 1=wrt component locations.",
"-setfom 0.8              Set the FOM of selected components to a specific value.",
" ",
"Selections:",
"-all                     Reset selection to all models and components before other selections.",
"-select #232@14          Select models and components.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
"-combined                The views for all models are combined (-Postscript option).",
"-separate                Each model is defined as a separate molecule group (-coordinates option).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
//"-reference file.pdb      Reference coordinate file to set views.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-split 3                 Split models into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": model ID's are used as file names.",
"-coordinates all.pdb     Output coordinate files.",
"-Postscript plot.ps      Output postscript file with a plot of views.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			reset(0);					// Keep selection as read from file
	Bstring			mod_select;					// Model and component selection
	int				merge(0);					// Flag to merge models rather than concatenate
	Bstring			calc_views;					// Mode to calculate component views
	int				inv_views(0);				// Flag to invert views for selected components
	string			asu_sym;					// Point group string
	int				view_dir_bins(0);			// Flag and histogram bin width to assess view directions
	int				view_dir_flag(0);			// Flag to set the reference for view directions
	double			setfom(-1);					// Value to set FOM to
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
	Bstring			paramfile;					// Input parameter file name
//	Bstring			reffile;					// Reference file name 
	Bstring			outfile;					// Output parameter file name
	Bstring			coorfile;					// Output coordinates file name
	int				split(0);					// Output one big STAR file
	int				combined(0);				// Flag to combine models for Postscript output
	Bstring			ps_file;					// Postscript output
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "select" ) mod_select = curropt->value;
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "invert" ) inv_views = 1;
		if ( curropt->tag == "setasu" ) asu_sym = curropt->value.upper().str();
		if ( curropt->tag == "directions" )
			if ( curropt->values(view_dir_bins, view_dir_flag) < 1 )
				cerr << "-directions: A histogram bin width in degrees must be specified!" << endl;
		if ( curropt->tag == "setfom" )
			if ( ( setfom = curropt->value.real() ) < 0 )
				cerr << "-setfom: A FOM value must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) <= 0 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) <= 0 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "combined" ) combined = 1;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
		}
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
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

	Bmodel*			model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( reset ) models_process(model, model_reset_selection);

	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( merge ) model_merge(model);
	
	if ( comprad > 0 ) models_process(model, comprad, model_set_component_radius);

	if ( linkrad > 0 ) models_process(model, linkrad, model_set_link_radius);

	if ( setfom >= 0 ) model->set_component_fom(setfom);
	
	if ( calc_views.length() ) model_calculate_views(model, calc_views);
	
	if ( inv_views ) model_invert_views(model);
	
//	cout << "Symmetry " << asu_sym << endl;
	if ( asu_sym.size() ) models_process(model, asu_sym, model_find_asymmetric_unit);
	
	if ( view_dir_bins ) model_view_directions(model, view_dir_bins, view_dir_flag);
	
	if ( ps_file.length() )
		ps_model_symmetry_views(ps_file, model, asu_sym, combined);
	
	model_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( model && ( outfile.length() || split == 9 ) )
		write_model(outfile, model, split);

	model_kill(model);
		
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

