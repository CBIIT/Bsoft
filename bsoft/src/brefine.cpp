/**
@file	brefine.cpp
@brief	Reciprocal space refinement of orientation parameters of particle images.
@author  Bernard Heymann
@date	Created: 20070115
@date	Modified: 20190207
**/

#include "mg_processing.h"
#include "mg_refine.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "rwFSC_XML.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: brefine [options] input.star [input.star]",
"------------------------------------------------",
"Refines the orientations of particle images with respect to a reference map.",
"Defocus refinement on non-CTF corrected images is turned on by the -defocus option.",
" ",
"Actions:",
"-all                     Reset selection to all particles before other selections (default not).",
"-symmetry C5             Search for the best symmetry-related orientations.",
"-magnification 0.01      Turns on magnification refinement: Maximum adjustment (default 0).",
"-defocus 0.002           Turns on defocus refinement: Standard deviation for Monte Carlo,",
"                         Step size for grid search (default 0 um).",
" ",
"Actions for the type of refinement:",
"-monte 80                Monte Carlo with a maximum number of iterations.",
"-grid 2.5,0.3            Grid search with angular step size and accuracy (both in degrees).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-select 14               Selection number of particles to process (default all selected).",
"-fomtype dpr             Resolution measurement type: frc (default), dpr.",
"-resolution 20,250       High and low resolution limits (angstrom).",
"-kernel 11,3             Interpolation kernel width and power (default 8,2)",
"-edge 125                Edge mask radius using previous particle origins (default not used).",
//"-threads 4               Set the number of threads (default maximum allowed).",
"-nothreads               Do not use threads (default parallel processing).",
" ",
"Parameters only for Monte Carlo refinement:",
"-shift 0.4               Random origin shift standard deviation (default 0.1 pixels).",
"-view 0.02               Random view shift standard deviation (default 0.01).",
"-angle 0.7               Random rotation angle maximum adjustment (default 0.5 degrees).",
" ",
"Parameters only for grid search:",
"-step 1.5,0.3            Grid search shift step size and accuracy (default 0.4,0.1 pixels).",
" ",
"Input:",
"-reference file.ext      Input 3D reference map file.",
"-mask mask.map           Real space 2D mask to be applied to particles.",
"-FSC file.xml            FSC curve used to weigh FOM contributions.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-ppx                     Write temporary particle parameter files to directory \"ppx\".",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int 			reset(0);					// Keep selection as read from file
	int				part_select(-1);			// Process all selected particles
	int				max_iter(0);				// Maximum Monte Carlo refining iterations
	double			alpha_step(0);				// Angular step size for grid search
	double			accuracy(0.5);				// Angular accuracy for grid search
	double			shift_step(0.4);			// Angular step size for grid search
	double			shift_accuracy(0.1);		// Angular accuracy for grid search
	double			edge_radius(0);				// Edge mask radius
	double			def_std(0);					// Random defocus standard deviation
	double			shift_std(0.1);				// Random origin shift standard deviation
	double			view_std(0.01);				// Random view shift standard deviation
	double			max_angle(M_PI/360);		// Random rotation angle maximum adjustment
	double			max_mag(0);					// Random magnification maximum adjustment
	int				fom_type(0);				// Default type is FRC
	double 			hi_res(0), lo_res(0); 		// Must be set > 0 to limit resolution
	int				kernel_width(8);			// Interpolation kernel width
	int				kernel_power(2);			// Interpolation kernel power
	int				nthreads(1);				// Number of threads to run
	int				nothreads(0);				// Flag to turn off threads
	Bstring			symmetry_string;			// No symmetry specified
	int				flags(0);					// Flags for processing options
	Bstring			reffile;					// Reference map file
	Bstring			maskfile;					// Mask to be applied to particles
	Bstring			fscfile;					// FSC curve
	Bstring			outfile;					// Output parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*			curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 0 )
				cerr << "-select: A selection number must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "monte" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-monte: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "grid" ) {
			if ( curropt->values(alpha_step, accuracy) < 1 )
				cerr << "-grid: An angular step size must be specified!" << endl;
			else {
				alpha_step *= M_PI/180.0;
				accuracy *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "fomtype" ) {
			if ( curropt->value[0] == 'f' ) fom_type = 0;
			if ( curropt->value[0] == 'd' ) fom_type = 1;
		}
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
			else
				if ( lo_res > 0 && lo_res < hi_res ) swap(hi_res, lo_res);
		}
		if ( curropt->tag == "kernel" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-kernel: A kernel width must be specified!" << endl;
		if ( curropt->tag == "edge" )
			if ( ( edge_radius = curropt->value.real() ) < 0 )
				cerr << "-edge: An edge radius must be specified!" << endl;
		if ( curropt->tag == "defocus" ) {
			if ( ( def_std = curropt->value.real() ) < 0.00001 )
				cerr << "-defocus: An defocus standard deviation must be specified!" << endl;
			else
				if ( def_std < 10 ) def_std *= 1e4;	// Assume um
		}
		if ( curropt->tag == "magnification" )
			if ( ( max_mag = curropt->value.real() ) < 0 )
				cerr << "-magnification: A maximum magnification adjustment must be specified!" << endl;
		if ( curropt->tag == "shift" )
			if ( ( shift_std = curropt->value.real() ) < 0.001 )
				cerr << "-shift: An origin shift standard deviation must be specified!" << endl;
		if ( curropt->tag == "view" )
			if ( ( view_std = curropt->value.real() ) < 0.001 )
				cerr << "-view: A view shift standard deviation must be specified!" << endl;
		if ( curropt->tag == "angle" ) {
			if ( ( max_angle = curropt->value.real() ) < 0.1 )
				cerr << "-angle: A maximum rotation angle adjustment must be specified!" << endl;
			else
				max_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "step" )
			if ( curropt->values(shift_step, shift_accuracy) < 1 )
				cerr << "-step: A shift step size must be specified!" << endl;
		if ( curropt->tag == "threads" ) {
			if ( curropt->value[0] == 'm' ) nthreads = 1000;
			else if ( ( nthreads = curropt->value.integer() ) < 1 )
				cerr << "-threads: A number of threads must be specified!" << endl;
		}
		if ( curropt->tag == "nothreads" ) nothreads = 1;
		if ( curropt->tag == "ppx" ) flags |= WRITE_PPX | CHECK_PPX;
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "FSC" )
			fscfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();

#ifdef HAVE_GCD
	if ( !nothreads ) {
		fftwf_init_threads();
		fftwf_plan_with_nthreads(system_processors());
		if ( verbose )
			cout << "Number of threads:              " << system_processors() << endl;
	}
#endif
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
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

	project->fom_tag[0] = FOM;
	project->fom_tag[1] = FOM_CV;

	if ( reset ) part_reset_selection(project, 3);
	
	vector<double>	weight;
	int				i, j;
	Bplot*			plot = NULL;
	if ( fscfile.length() ) {
		plot = xml_read_fsc(fscfile);
		if ( !plot ) {
			cerr << "Error: The FSC file " << fscfile << " was not read!" << endl;
			bexit(-1);
		}
		weight.resize(plot->rows());
		for ( i=0, j=plot->rows(); i<plot->rows(); i++, j++ )
			weight[i] = sqrt((*plot)[j]);
		delete plot;
	}
	
	if ( outfile.length() ) {
		if ( reffile.length() && ( max_iter || alpha_step > 0 ) ) {
			if ( nothreads ) {
				if ( mg_refine_orientations(project, reffile, maskfile, 
						symmetry_string, part_select, max_iter, 
						alpha_step, accuracy, shift_step, 
						shift_accuracy, fom_type, weight,
						hi_res, lo_res, kernel_width, kernel_power, 
						edge_radius, def_std, shift_std, 
						view_std, max_angle, max_mag, flags) < 0 )
					bexit(-1);
			} else {
				if ( project_refine_orientations(project, reffile, maskfile, 
						symmetry_string, part_select, max_iter, 
						alpha_step, accuracy, shift_step, 
						shift_accuracy, fom_type, weight,
						hi_res, lo_res, kernel_width, kernel_power, 
						edge_radius, def_std, shift_std, 
						view_std, max_angle, max_mag, flags) < 0 )
					bexit(-1);
			}
		}
		
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

#ifdef HAVE_GCD
	if ( !nothreads ) fftwf_cleanup_threads();
#endif

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

