/**
@file	bplane.cpp
@brief	A tool to expand models.
@author Bernard Heymann
@date	Created: 20140925
@date 	Modified: 20141008
**/

#include "model_create.h"
#include "model_plane.h"
#include "model_transform.h"
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
"Usage: bplane [options] in.star",
"-------------------------------",
"Generates and manipulates plane models.",
" ",
"Actions for creation:",
"-plane 386,125           Generate a plane with the given length and width.",
"-setviews                Fit a plane to the model and set the component views.",
"-fit                     Use the input model as guide to generate a plane model.",
" ",
"Actions:",
"-view 0.5,-0.2,0.3,35    View to orient plane: vector=normal, angle=in-plane rotation.",
"-guide poly.star         Input guide model to adjust the plane to.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-separation 14.5         New component separation (default 10).",
"-sigma 150               Smoothness of the fit (default 100).",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
" ",
"Input:",
"-map density.map         Input map file.",
" ",
"Output:",
"-output newmod.star      Output model file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	double			length(0), width(0);		// Length and width of new plane model
	int				setviews(0);				// Flag to set the views based on a plane fit
	int				fit(0);						// Generate a fitted plane model
//	View			view;						// Orientation for plane
	View2<float>	view;						// Orientation for plane
	double			separation(10);				// New component separation
	double			sigma(100);					// Smoothness of fit to guide
//	int				normal(0);					// Move components only along view vector
	double			compradius(0);				// Component display radius
	double			linkradius(0);				// Link display radius
	Bstring			paramfile;					// Input parameter file name
	Bstring			guidefile;					// Input guide model file name
	Bstring			mapfile;					// Input map file name
	Bstring			outfile;					// Output parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "plane" )
			if ( curropt->values(length, width) < 2 )
				cerr << "-plane: Both length and width must be specified!" << endl;
		if ( curropt->tag == "setviews" ) setviews = 1;
		if ( curropt->tag == "fit" ) fit = 1;
		if ( curropt->tag == "view" ) {
			View 	tv = curropt->view();
			view = View2<float>(tv[0],tv[1],tv[2],tv[3]);
		}
		if ( curropt->tag == "separation" )
			if ( ( separation = curropt->value.real() ) < 0.01 )
				cerr << "-separation: The component separation must be specified!" << endl;
		if ( curropt->tag == "sigma" )
			if ( ( sigma = curropt->value.real() ) < 0.01 )
				cerr << "-sigma: A distance must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 1 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 1 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "guide" )
			guidefile = curropt->filename();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	Bmodel*			model = NULL;		
	Bmodel*			guide = NULL;

	if ( length > 0 && width > 0 ) {
		model = model_create_plane(length, width, 0, separation);
		if ( !model ) {
			cerr << "Error: The model creation failed!" << endl;
			bexit(-1);
		}
		model_rotate(model, view);
		model_set_views(model, view);
	} else {
		while ( optind < argc ) string_add(&file_list, argv[optind++]);
		if ( file_list ) {
			model = read_model(file_list, paramfile);
			string_kill(file_list);
		}
	}
	
	if ( guidefile.length() ) {
		guide = read_model(guidefile, paramfile);
		if ( !model ) model = guide;
	}

	if ( compradius ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius ) models_process(model, linkradius, model_set_link_radius);
	
	if ( mapfile.length() ) model->mapfile(mapfile.str());
	
	if ( setviews ) model_fit_plane(model);
	
	if ( fit ) {
		guide = model;
		model = model_generate_from_plane_guide(guide, separation, sigma);
		model_kill (guide);
		if ( !model ) {
			cerr << "Error: The model creation from a guide failed!" << endl;
			bexit(-1);
		}
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



