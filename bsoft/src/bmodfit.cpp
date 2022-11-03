/**
@file	bmodfit.cpp
@brief	A tool to fit models.
@author Bernard Heymann
@date	Created: 20220223
@date	Modified: 20220223
**/

//#include "model_create.h"
//#include "model_plane.h"
#include "model_transform.h"
//#include "model_views.h"
//#include "model_shell.h"
//#include "model_links.h"
//#include "model_util.h"
//#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			model_fit(Bmodel* model, Bmodel* refmod, string id);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodfit [options] in.star",
"--------------------------------",
"Calculate fits of model to a reference.",
" ",
"Actions:",
"-fit mod_id              Fit a model with the given id to a reference.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Input:",
"-reference ref.pdb       Input reference model file.",
" ",
"Output:",
"-output newmod.star      Output model file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	string			fit_id;						// Model ID to fit
	Bstring			paramfile;					// Input parameter file name
	Bstring			reffile;					// Input reference model file name
	Bstring			outfile;					// Output model file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "fit" ) fit_id = curropt->value.str();
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	Bmodel*			model = NULL;
	Bmodel*			refmod = NULL;

	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( file_list ) {
		model = read_model(file_list, paramfile);		
		string_kill(file_list);
	}
	
	cout << "Models: " << model->count() << endl;
	
	if ( reffile.length() ) {
		refmod = read_model(reffile, paramfile);
		if ( !refmod ) {
			error_show("Error: No reference model read!", __FILE__, __LINE__);
			bexit(-1);
		}
	}

	if ( fit_id.length() > 0 ) model_fit(model, refmod, fit_id);

	cout << "Models: " << model->count() << endl;
	
	// Write an output model file if a name is given
    if ( outfile.length() && model ) {
//		model_update_comment(model, argc, argv);
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Fits a model to a reference model.
@param 	model			model structure.
@param 	refmod			reference model.
@param 	id				model id to fit.
@return int				0.

**/
int			model_fit(Bmodel* model, Bmodel* refmod, string id)
{
	Bmodel*			mp = model;
	while ( mp && mp->identifier() != id ) mp = mp->next;
	
	if ( !mp ) {
		cerr << "Error: Model with id " << id << " not found!" << endl;
		bexit(-1);
	}
	
	Transform		t = model_find_transform(mp, refmod);

	models_rotate(model, t);
	
	return 0;
}

