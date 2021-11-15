/**
@file	bpart.cpp
@brief	Calculates centers of single particle images.
@author Bernard Heymann
@date	Created: 20080424
@date	Modified: 20180102
**/

#include "mg_processing.h"
#include "mg_particles.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpart [options] input.star [input.star]",
"----------------------------------------------",
"Manipulates single particle images.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-all                     Reset selection to all particles before other selections.",
"-find 20                 Finds the centers of the particle images (iterations).",
"-shift                   Center the actual particle images (new images with \"_cen\" in names).",
"-align                   Output aligned particle images (new images with \"_aln\" in names).",
"-rescale -0.1,5.2        Rescale particle images to average and standard deviation.",
"-settilt -30,45          Set micrograph tilt, axis and particle defocus.",
"-calctilt                Calculate micrograph tilt parameters from particle defocus.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-select 14               Selection number of particles to process (default all selected).",
"-resolution 10,500       Resolution limits for cross-correlation (angstrom, default 20,1000).",
"-partpath dir/subdir     Set the particle file paths.",
" ",
"Input:",
"-mask mask.pif           Intput 3D mask file to mask particles.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-reference file.pif      Composite image reference file.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int				use_rec(0);				// Flag to process reconstructions
	int 			reset(0);				// Keep selection as read from file
	int 			find(0);				// Iterations to find the centers
	int				shift(0);				// Flag to shift the actual images
	int				align(0);				// Flag to align images
	double			nuavg(0), nustd(0); 	// Rescaling to average and stdev
	int				calc_tilt(0);			// Flag to calculate micrograph tilt parameters
	double			axis(0), tilt(0);		// Tilt axis and angle to set
	int				part_select(-1);		// Process all selected particles
	double			hires(20);				// High resolution limit
	double			lores(1000);			// Low resolution limit
	Bstring			partpath;				// Particle file path
	Bstring			maskfile;				// Particle mask file
	Bstring			outfile;				// Output parameter file
	Bstring			reffile;				// Composite image reference file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "find" )
			if ( ( find = curropt->value.integer() ) < 1 )
				cerr << "-find: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "shift" ) shift = 1;
		if ( curropt->tag == "align" ) align = 1;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "calctilt" ) calc_tilt = 1;
		if ( curropt->tag == "settilt" ) {
			if ( curropt->values(tilt, axis) < 1 )
				cerr << "-settilt: A tilt angle must be specified!" << endl;
			else {
				tilt *= M_PI/180.0;
				axis *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 0 )
				cerr << "-select: A selection number must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

	if ( hires > lores ) swap(hires, lores);
	
	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( use_rec ) project->select = 1;
	
	if ( reset ) part_reset_selection(project, 3);
	
	Bimage*			pref = NULL;
	Bimage*			pmask = NULL;
	
	if ( find ) 
		pref = project_find_particle_centers(project, find, part_select, hires, lores);
	
	if ( align )
		project_align_particles(project, part_select, nuavg, nustd);
	else if ( shift )
		project_center_particles(project, part_select, nuavg, nustd);

	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		if ( pmask ) {
			project_mask_particles(project, pmask, partpath);
			delete pmask;
		}
	}

	if ( calc_tilt )
		project_tilt_from_particle_defocus(project);

	if ( fabs(tilt) > 0.01 )
		project_set_particle_defocus_from_tilt(project, axis, tilt);
	
	if ( outfile.length() )
		write_project(outfile, project, 0, 0);
	
	if ( pref && reffile.length() )
		write_img(reffile, pref, 0);
	
	delete pref;
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

