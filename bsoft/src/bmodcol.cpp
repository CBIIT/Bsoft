/**
@file	bmodcol.cpp
@brief	A tool to color models.
@author Bernard Heymann
@date	Created: 20080206
@date 	Modified: 20210319
**/

#include "rwmodel.h"
#include "model_color.h"
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
"Usage: bmodcol [options] in.star",
"--------------------------------",
"Colors models based on selection.",
" ",
"Actions for preparation:",
"-all                     Reset selection to all components before other selections.",
" ",
"Selections:",
"-select #232@14          Select models and components.",
"-valence 3               Select components with the indicated valence.",
"-polygons 6              Select polygonal components with the indicated order.",
" ",
"Coloring actions:",
"-color 0.1,0.5,0.2,1     Color selected components and links.",
"-spectrum                Color a model from red to blue by component order.",
"-density                 Color a model from red to blue by component density.",
"-fom                     Color a model from red to blue by component FOM.",
"-number                  Color a model from red to blue by component selection number.",
"-curvature               Color a model from red to blue by curvature.",
"-chiral                  Color chiral vertices red (-) and blue (+).",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components before other selections.",
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
"-output color.star       Color model file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int 			all(0);					// Keep selection as read from file
	int 			reset(0);				// Keep selection as ouput
	Bstring			mod_select;				// Model and component selection
	int				valence(0);				// Valency selector
	int				poly_order(0);			// Polygon selector
	RGBA<float>		color(-1,-1,-1,1);		// Component color specification
	int				spectrum(0);			// Flag to color models by component order
	int				density(0);				// Flag to color models by component density
	int				fom(0);					// Flag to color models by component FOM
	int				number(0);				// Flag to color models by component selection number
	int				curvature(0);			// Flag to color models by curvature
	int				chiral(0);				// Flag to color chiral vertices
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	Bstring			paramfile;				// Input parameter file name
	Bstring			outfile;				// Output parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "valence" )
			if ( ( valence = curropt->value.integer() ) < 1 )
				cerr << "-valence: A valency must be specified!" << endl;
		if ( curropt->tag == "polygons" )
			if ( ( poly_order = curropt->value.integer() ) < 1 )
				cerr << "-polygons: A polygon order must be specified!" << endl;
		if ( curropt->tag == "color" )
			if ( curropt->values(color[0], color[1], color[2], color[3]) < 3 )
				cerr << "-color: At least three color values must be specified!" << endl;
		if ( curropt->tag == "spectrum" ) spectrum = 1;
		if ( curropt->tag == "density" ) density = 1;
		if ( curropt->tag == "fom" ) fom = 1;
		if ( curropt->tag == "number" ) number = 1;
		if ( curropt->tag == "curvature" ) curvature = 1;
		if ( curropt->tag == "chiral" ) chiral = 1;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 1 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 1 )
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
	
	if ( all ) models_process(model, model_reset_selection);
	
	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( valence > 0 ) model_select_valence(model, valence);
	
	if ( poly_order > 0 ) model_select_polygons(model, poly_order);
	
	if ( spectrum ) model_color_by_order(model);
	else if ( curvature ) model_color_curvature(model);
	else if ( chiral ) model_color_chiral_vertices(model);
	else if ( density ) model_color_by_density(model);
	else if ( fom ) model_color_by_fom(model);
	else if ( number ) model_color_by_selection(model);
	else if ( color[0] >= 0 && color[1] >= 0 && color[2] >= 0 )
		model_color_selected(model, color);
	
	if ( reset ) models_process(model, model_reset_selection);
	
	if ( compradius ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius ) models_process(model, linkradius, model_set_link_radius);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

