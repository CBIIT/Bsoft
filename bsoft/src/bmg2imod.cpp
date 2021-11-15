/**
@file	bmg2imod.cpp
@brief	Converts between IMOD files and a micrograph parameter file
@author	Bernard Heymann
@date	Created: 20070501
@date	Modified: 20180322

**/

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "rwmgIMOD.h"
#include "rwimg.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			write_IMOD_parameters(Bstring& imodfile, Bproject* project);

// Usage assistance
const char* use[] = {
" ",
"Usage: bmg2imod [options] input.mrc input.xf input.tlt input.xyz",
"----------------------------------------------------------------",
"Converts between IMOD files and a micrograph parameter file.",
" ",
"Actions:",
"-topif                   Flag to convert the image file to multi-image PIF (default not).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5            Sampling (A/pixel).",
" ",
"Output:",
"-output file.star        Output micrograph parameter file.",
"-imod file.xf            Output IMOD parameter files.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int				flag(0);			// Flag to rewrite the image file
	Vector3<double>	sam;				// A/pixel
	Bstring			outfile;			// Output parameter file
	Bstring			imodfile;			// IMOD output parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "topif" ) flag = 1;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "imod" )
			imodfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list, flag);

	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No project generated!" << endl;
		bexit(-1);
	}

	if ( sam[0] > 0.1 )
		project_set_mg_pixel_size(project, sam);

	if ( outfile.length() )
		write_project(outfile, project, 0, 0);
	
	if ( imodfile.length() )
		write_project_imod(imodfile, project);
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

