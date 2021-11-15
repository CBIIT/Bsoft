/**
@file	bbild.cpp
@brief	Generates BILD output from a model.
@author Bernard Heymann
@date	Created: 20080521
@date	Modified: 20160616
**/

#include "rwmodel.h"
#include "rwmodel_bild.h"
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
"Usage: bbild [options] in1.star [in2.star ...]",
"----------------------------------------------",
"Generates BILD output from a model.",
" ",
"Actions:",
"-all                     Reset selection to all models before other selections.",
"-type poly               Type of BILD output: views, polygons or neighbors.",
"-linklength 56.3         Generate links with this maximum length.",
"-color spectrum          Color selection or polygons: spectrum=red->blue, density, fom.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
" ",
"Output:",
"-output file.bild        Output model file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			reset(0);				// Keep selection as read from file
	Bstring			bild_type("views");		// Default BILD output is views
	double			linklength(0);			// Link length for generating links
	int				color_type(0);			// Color scheme type: 1=spectrum, 2=fom
	double			comprad(0);				// Component display radius
	double			linkrad(0);				// Link display radius
	Bstring			outfile;				// Output parameter file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "type" ) bild_type = curropt->value.lower();
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 0.1 )
				cerr << "-linklength: A link length must be specified!" << endl;
		if ( curropt->tag == "color" ) {
			if ( curropt->value.contains( "spe") ) color_type = 1;
			if ( curropt->value.contains( "den") ) color_type = 2;
			if ( curropt->value.contains( "fom") ) color_type = 3;
		}
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 0.001 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 0.001 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the model parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bmodel*		model = read_model(file_list);
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( reset ) models_process(model, model_reset_selection);

	if ( linklength > 0 ) model_link_list_generate(model, linklength);

	if ( comprad > 0 ) models_process(model, comprad, model_set_component_radius);

	if ( linkrad > 0 ) models_process(model, linkrad, model_set_link_radius);

	model_selection_stats(model);
	
	// Write an output parameter format file if a name is given
	Bstring			ext;
	int				order(6);
	
	if ( outfile.length() ) {
		ext = outfile.extension();
		if ( ext.contains("bild") || ext.contains("bld") ) {
			if ( bild_type[0] == 'v' )	// Views
				model_to_bild_orientations(outfile, model, 1, color_type);
			else if ( bild_type[0] == 'f' )	// Force vectors
				model_to_bild_orientations(outfile, model, 2, color_type);
			else if ( bild_type[0] == 's' )	// Velocity vectors
				model_to_bild_orientations(outfile, model, 3, color_type);
			else if ( bild_type[0] == 'b' )	// View sphere
				model_to_bild_view_sphere(outfile, model, color_type);
			else if ( bild_type[0] == 'c' )	
				model_to_bild_view_polygons(outfile, model, order, color_type);
			else if ( bild_type[0] == 'p' )	// Polygons
				model_to_bild_polygons(outfile, model, color_type);
			else if ( bild_type[0] == 'n' )	// Neighbor planes
				model_to_bild_neighbor_planes(outfile, model, color_type);
			else
				write_model_bild(outfile, model);
		} else {
				write_model(outfile, model);
		}
	}
	
	model_kill(model);
		
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

