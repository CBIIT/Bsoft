/**
@file	borient.cpp
@brief	Determines orientation angles and x,y origins of single particle images 
@author	Bernard Heymann and David M. Belnap
@date	Created: 20010403
@date	Modified: 20220906 (BH)
**/

#include "mg_orient.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_reconstruct.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern int	thread_limit;	// Thread limit

// Usage assistance
const char* use[] = {
" ",
"Usage: borient [options] input.star [input.star]",
"------------------------------------------------",
"Determines the orientations and origins of single particles with respect to a reference.",
"Options -resolution and -bands are mutually exclusive.",
" ",
"Actions:",
"-mode scc                Mode: 0,pps: polar power spectrum only (fast, default),",
"                         1,scc: polar power spectrum with selected cross correlation.",
"                         2,ccc: polar power spectrum with all cross correlation (slow).",
"                         3,ori: origin determination by cross correlation.",
"-CTF                     Apply CTF to projections (default not).",
"-contrast                Invert contrast for comparisons (default not).",
"-all                     Reset selection to all particles before other selections (default not).",
"-side 15                 Generate side view projections within the given angle from the equator.",
"-bin 2,1                 Bin particles and reference by the given kernel size (default 1,1).",
"-prepare 45,12           Prepare a set of particle images as 2D references (number and first image).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-select 14               Selection number of particles to process (default all selected).",
"-symmetry C5             Point group symmetry.",
"-angles 2.5,5.2          Step size for theta and phi, one value sets both (default 45 degrees).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given.)",
"-resolution 20,360       High and low resolution limits (angstrom).",
"-polarresolution 23.5    Resolution limit for in-plane angular alignment (angstrom).",
"-bands 250,-1,150,1,30,0 Bands for analysis in reciprocal space (pairs of numbers: angstrom,flag).",
"                         Band with 0: Not used - end flag.",
"                         Band with 1: Use for orientation determination.",
"                         Band with -1: Use for cross-validation.",
"-annuli 2,95             Real space annular limits (default 0,inf pixels).",
"-shiftlimit 3.5          Limit on origin shift relative to nominal center (default 10% of box edge size).",
"-rotationlimit 16.7      Limit on in-plane rotation change from input orientation (default not applied).",
"-edge 125                Edge mask radius using previous particle origins (default not used).",
"-fom 0.3                 Set FOM threshold for selection (default 0).",
#ifdef HAVE_GCD
"-nothreads               Do not use threads (default parallel processing).",
#endif
" ",
"Input:",
"-reference file.ext      Input 3D reference map or projection file.",
"-mask mask.mrc           Real space 2D mask to be applied to projections and particles.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-log                     Write particle log files to directory \"log\".",
"-ppx                     Write temporary particle parameter files to directory \"ppx\".",
"-TwoD file.ext           File name for 2D real space reconstruction.",
"-oriented                Output oriented images with file names ending with \"_proj.spi\".",
"                         (must be used with the -TwoD option).",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int 			reset(0);					// Keep selection as read from file
	long			part_select(-1);			// Process all selected particles
	Bstring			symmetry_string("C1");		// Point group string
	double			theta_step(M_PI_4);			// Angular step size for theta
	double			phi_step(M_PI_4);			// Angular step size for phi
	double			side_ang(-1);				// Side view variation angle
	int				bin(0);						// Data compression by binning - particles
	int				refbin(0);					// Data compression by binning - reference
	long			fref(0), nref(0);			// First and number of reference images to prepare
	Vector3<double>	sam;		    			// Units for the three axes (A/pixel)
	double 			res_lo(10000);				// Resolution limits for orientation finding
	double 			res_hi(0);					// Must be set > 0 to determine orientations
	double 			res_polar(0);				// Must be set > 0 to limit in-plane angular resolution 
	vector<double>	band;						// Reciprocal space bands
	long 			ann_min(0);					// Minimum annulus
	long 			ann_max(1000); 				// Maximum annulus, reset to maximum radius in image
	double			shift_limit(-1);			// Maximum shift from nominal box origin
	double			angle_limit(0);				// Maximum in-plane rotation from input orientation
	double			edge_radius(0);				// Edge mask radius
	double 			FOM_cut(0); 				// Threshold to accept orientations
	int				mode(0);					// Only polar power spectrum
	int				transform_output(0);		// Flag to output transformed images
#ifdef HAVE_GCD
	int				nothreads(0);				// Flag to turn off threads
#endif
	int				flags(FULL_ASU);			// Flags for processing options
	Bstring			outfile;					// Output STAR format file
	Bstring			ref_filename;				// Input 3D map or projection file name
	Bstring			mask_filename;				// Mask to be applied to projections and particles
	Bstring			twoD_recons_name;			// Output 2D reconstruction file

	thread_limit = 1;		// Default number of threads

	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" )
			reset = 1;
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 0 )
				cerr << "-select: A selection number must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified" << endl;
			else {
				theta_step *= M_PI/180.0;
				if ( i < 2 ) phi_step = theta_step;
				else phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "side" ) {
			if ( ( side_ang = curropt->value.real() ) < 0.001 )
				cerr << "-side: One angle must be specified!" << endl;
			else
				side_ang *= M_PI/180.0;
		}
		if ( curropt->tag == "bin" ) {
			if ( curropt->values(bin, refbin) < 1 )
				cerr << "-bin: An integer must be specified!" << endl;
			else {
				if ( bin < 2 ) bin = 0;
				if ( bin > 4 ) bin = 4;
				if ( refbin < 2 ) refbin = 0;
				if ( refbin > 4 ) refbin = 4;
			}
		}
		if ( curropt->tag == "prepare" )
			if ( curropt->values(nref, fref) < 1 )
				cerr << "-prepare: A number of particle images must be specified!" << endl;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(res_hi, res_lo) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
			else
				if ( res_lo < res_hi ) swap(res_hi, res_lo);
		}
		if ( curropt->tag == "polarresolution" )
			if ( ( res_polar = curropt->value.real() ) < 0.1 )
				cerr << "-polarresolution: A resolution limit must be specified!" << endl;
        if ( curropt->tag == "bands" ) {
			band = curropt->value.split_into_doubles(",");
			if ( band.size() < 2 )
				cerr << "-bands: At least two pairs of numbers must be specified!" << endl;
			band[band.size()] = 0;
		}
		if ( curropt->tag == "annuli" )
			if ( curropt->values(ann_min, ann_max) < 1 )
				cerr << "-annuli: At least a minimum annulus must be specified!" << endl;
		if ( curropt->tag == "shiftlimit" )
			if ( ( shift_limit = curropt->value.real() ) < 0.1 )
				cerr << "-shiftlimit: A maximum shift in pixels must be specified!" << endl;
		if ( curropt->tag == "rotationlimit" ) {
			if ( ( angle_limit = curropt->value.real() ) < 0.1 )
				cerr << "-rotationlimit: A maximum rotation angle must be specified!" << endl;
			else
				angle_limit *= M_PI/180.0;
		}
		if ( curropt->tag == "edge" )
			if ( ( edge_radius = curropt->value.real() ) < 0 )
				cerr << "-edge: An edge radius must be specified!" << endl;
		if ( curropt->tag == "CTF" ) flags |= APPLY_CTF;
		if ( curropt->tag == "multiple" ) flags |= MULTI_FILE;
		if ( curropt->tag == "fom" )
			if ( ( FOM_cut = curropt->value.real() ) < 0.0001 )
				cerr << "-fom: A FOM threshold must be specified!" << endl;
#ifdef HAVE_GCD
		if ( curropt->tag == "nothreads" ) nothreads = 1;
#endif
		if ( curropt->tag == "mode" ) {
			if ( isdigit(curropt->value[0]) ) {
				mode = curropt->value.integer();
			} else {
				if ( curropt->value.contains("pps") ) mode = 0;
				else if ( curropt->value.contains("scc") ) mode = 1;
				else if ( curropt->value.contains("ccc") ) mode = 2;
				else if ( curropt->value.contains("ori") ) mode = 3;
			}
			if ( mode < 0 || mode > 3 ) mode = 0;
			flags |= mode;
		}
		if ( curropt->tag == "oriented" ) transform_output = 1;
		if ( curropt->tag == "contrast" ) flags |= INVERT;
		if ( curropt->tag == "log" ) flags |= PART_LOG;
		if ( curropt->tag == "ppx" ) flags |= WRITE_PPX | CHECK_PPX;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "reference" )
			ref_filename = curropt->filename();
		if ( curropt->tag == "mask" )
			mask_filename = curropt->filename();
		if ( curropt->tag == "TwoD" )
			twoD_recons_name = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

#ifdef HAVE_GCD
	if ( !nothreads ) {
		fftwf_init_threads();
//		fftwf_plan_with_nthreads(system_processors());
		if ( verbose )
			cout << "Number of threads:              " << system_processors() << endl;
	}
#endif
	
	Bsymmetry		sym(symmetry_string);
	
	if ( res_hi <= 0 ) res_hi = 0.1;

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Warning: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);

	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	Bimage*			proj = NULL;
	if ( ref_filename.length() )
		proj = img_prepare_projections(ref_filename, mask_filename,
						refbin, sym, theta_step, phi_step, side_ang);
	else if ( nref )
		proj = project_prepare_2D_references(project, fref, nref, bin, flags & APPLY_CTF);

	ref_filename = 0;
	
	project->fom_tag[0] = FOM;
	project->fom_tag[1] = FOM_CV;

//	if ( project_count_mg_particles(project) < 1 )
//		project_set_views_from_images(project); // This must be done before operations on particles
		
	if ( reset ) part_reset_selection(project, 3);

	if ( sam[0] > 0 && sam[1] > 0 && sam[2] > 0 ) {
		if ( sam[0] < 0.1 ) sam = Vector3<double>(1,1,1);
		project_set_mg_pixel_size(project, sam);
	}

	if ( proj && outfile.length() ) {
		if ( mode == 3 ) {
			if ( project_determine_origins(project, proj, bin, sym,
					part_select, res_lo, res_hi, shift_limit, flags) < 0 ) {
				cerr << "Error: Alignment failed!" << endl;
				bexit(-1);
			}
		} else {
			if ( project_determine_orientations2(project, proj, mask_filename,
					bin, sym, part_select, band, res_lo, res_hi, res_polar,
					ann_min, ann_max, shift_limit, angle_limit, edge_radius, flags) < 0 ) {
				cerr << "Error: Alignment failed!" << endl;
				bexit(-1);
			}
		}
	}
	
	if ( proj ) delete proj;
	mask_filename = 0;
	
	if ( FOM_cut > 0 ) part_deselect(project, 0, FOM_cut);
	
	Bimage* 		prec = NULL;
	if ( twoD_recons_name.length() ) {
		if ( transform_output )
			prec = project_reconstruct_2D(project, twoD_recons_name, transform_output);
		else
			prec = project_reconstruct_2D_fast(project, twoD_recons_name);
		if ( prec ) {
			write_img(twoD_recons_name, prec, 0);
			delete prec;
		}
	} else if ( verbose > VERB_RESULT ) {
		project_show_class_averages(project);
	}

	time_t		t = time(NULL);
	project->comment += "# Finished at: ";
	project->comment += asctime(localtime(&t));
	project->comment += "\n";
	
	// Write an output STAR format file if a name is given
    if ( outfile.length() )
		write_project(outfile, project, 0, 0);

	project_kill(project);
	outfile = 0;
	symmetry_string = 0;

#ifdef HAVE_GCD
	if ( !nothreads ) fftwf_cleanup_threads();
#endif

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

