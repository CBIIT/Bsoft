/**
@file	emgrand.cpp
@brief	Generate random values for micrograph parameters.
@author Eduardo Sanz-Garcia, David Belnap and Bernard Heymann
@date	Created: 20051011
@date	Modified: 20100128
**/

#include "mg_random.h"
#include "rwmg.h"
#include "symmetry.h"
#include "utilities.h"
#include "timer.h"
#include "options.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: emgrand [options] input.star [input.star]",
"------------------------------------------------",
"Randomizes particle parameters.",
" ",
"Actions:",
"-randomviews             Replace views with random orientations.",
"-helical                 Replace views with random orientations for a helix.",
"-symmetry C5             Replace views with random symmetry-related orientations.",
"-origins 1.4             Introduce random errors into origins with this standard deviation (pixels).",
"-views 3.1               Introduce random errors into views with this standard deviation (degrees).",
"-defocus 0.3             Introduce random errors into defocus with this standard deviation (um).",
"-magnification 0.1       Introduce random errors into magnifications with this standard deviation (fraction).",
"-setasu C5               Set the views to the asymmetric unit (must agree with -symmetry if used).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Output:",
"-output file.star        Output parameter file name.",
" ",
NULL
};


int		main(int argc, char **argv)
{
	// Initialize optional variables
	double			randomorigins(0);	// Standard deviation for origin errors (angstrom)
	double			randomviews(0);		// Standard deviation for view errors (degrees)
	double			randomdefocus(0);	// Standard deviation for defocus errors (angstrom)
	double			randommag(0);		// Standard deviation for magnification errors (fraction)
	Bstring			outfile;			// Output parameter file
	Bsymmetry		sym;				// No symmetry specified
	int			 	random(0);			// Flag to generate random views
	int				helical(0);			// Flag to generate random helical views
	int			 	randsym(0);			// Flag to generate random symmetry-related views
	int				asu(0);				// Flag to change views to ASU

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "randomviews" )
			random = 1;
		if ( curropt->tag == "helical" )
			helical = 1;
		if ( curropt->tag == "symmetry" ) {
			sym = curropt->symmetry();
			randsym = 1;
		}
		if ( curropt->tag == "origins" ) {
			if ( ( randomorigins = curropt->value.real() ) < 0.00001 )
				cerr << "-origins: A standard deviation must be specified!" << endl;
			else if ( randomorigins < 0 ) randomorigins = 0;
		}
		if ( curropt->tag == "views" ) {
			if ( ( randomviews = curropt->value.real() ) < 0.00001 )
				cerr << "-views: A standard deviation must be specified!" << endl;
			else {
				randomviews *= M_PI/180.0;		// Assume degrees
				if ( randomviews < 0 ) randomviews = 0;
			}
		}
		if ( curropt->tag == "defocus" ) {
			if ( ( randomdefocus = curropt->value.real() ) < 0.00001 )
				cerr << "-defocus: A standard deviation must be specified!" << endl;
			else {
				if ( randomdefocus < 100 ) randomdefocus *= 1e4;	// Assume um
				if ( randomdefocus < 0 ) randomdefocus = 0;
			}
		}
		if ( curropt->tag == "magnification" ) {
			if ( ( randommag = curropt->value.real() ) < 0.00001 )
				cerr << "-magnification: A standard deviation must be specified!" << endl;
			else
				if ( randommag < 0 ) randommag = 0;
		}
		if ( curropt->tag == "setasu" ) {
			sym = curropt->symmetry();
			asu = 1;
		}
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
	}
	option_kill(option);

	double       ti = timer_start();

	// Read parameter file(s)
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*			project = read_project(file_list);
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( random )
		project_random_views(project);

	if ( helical )
		project_random_helical_views(project);

	if ( sym.point() > 101 && randsym )
		project_random_symmetry_related_views(project, sym);

	if ( randomorigins )
		project_random_origins(project, randomorigins);

	if ( randomviews )
		project_random_views(project, randomviews);

	if ( randomdefocus )
		project_random_defocus(project, randomdefocus);

	if ( randommag )
		project_random_magnification(project, randommag);

	if ( sym.point() > 101 && asu )
		project_set_particle_asu_views(project, sym);

	// Output parameters
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}

	// Memory cleanup
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

