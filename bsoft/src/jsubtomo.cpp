/**
@file	jsubtomo.cpp
@brief	Finds particles (subtomograms) in a tomographic reconstruction and refines their origins and orientations
@author Juha Huiskonen (JTH)
@author	Bernard Heymann
@date	Created:  20071010
@date	Modified: 20080630 -rescale option.
@date	Modified: 20081124 option to write star files of extracted particles
@date	Modified: 20081128 rescaling of the map moved to main loop, output files will be float by default
@date	Modified: 20100503 changes to how extract mode handles particle origins
@date	Modified: 20100511 target map cropped for refinement, cropped map rescaled
@date	Modified: 20120123 (BH)
@date	Modified: 20120308
@date	Modified: 20120316 added an option to limit shifts along view vectors
@date	Modified: 20120329 added functional topp and topn options
@date	Modified: 20131129 (BH)
@date	Modified: 20150108 (BH) - incorporated into Bsoft
@date	Modified: 20150806 (BH)

	Based on the code by Bernard Heymann (BH) (bfind.c in bsoft)
**/

#include "rwimg.h"
#include "mg_subtomo.h"
#include "mg_processing.h"
#include "mg_particle_select.h"
#include "mg_extract.h"
#include "mg_img_proc.h"
#include "rwmg.h"
#include "ps_micrograph.h"
#include "Matrix.h"
#include "linked_list.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: jsubtomo [options] input.star/input.img [input.star/input.img]",
"---------------------------------------------------------------------",
"Finds particles (subtomograms) in a tomographic reconstruction and refines their orientations and origins.",
"Extracts particles from the tomograms and transforms them for averaging,",
"or applies the inverse transformations on the template to position it in the original tomogram.",
" ",
"Actions mode:",
"-mode refine                 Search, Refine or Extract.",
" ",
"                             Search:     Searches the particles in the input.image.",
"                             Refine:     Refines the orientations of the particles in the input.star.",
"                             Extract:    Writes selected particles to image files.",
" ",
"Actions:",
"-transform                   Transform the particle before extraction (use only in Extract mode).",
"-inverse                     Transform the template instead of the particle (use only in Extract mode).",
"-star                        Output a star file for each extracted particle (use only in Extract mode).",
"-wedges                      Output rotated missing wedge masks for each particle (use only in Extract mode).",
"-iterate 2                   Iterate angular refinement (default 1 iteration).",
"                             The original step sizes and limits are halved after every iteration (use only in Refine mode with angle limits).",
" ",
"Parameters:",
"-verbose 7                   Verbosity of output.",
"-datatype u                  Force writing of a new data type (default is float).",
"-sampling 1.5,1.5,1.5        Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 110,50,44            Origin for rotations of the template (default taken from template image).",
"-convert                     Convert the first particle in the star-file into a reconstruction.",
" ",
"Parameters for search and refinement:",
"-angles 8.8,2.5,5.2          Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-alphalimit 30               Angular limit for alpha (default: 180 degrees).",
"-thetaphilimit 5             Angular limit for theta & phi (default: 180 degrees, cannot be used with option symmetry).",
"-shiftlimit 5                Radial limit to search for shift in cross-correlation map, around the location given in input.star.",
"                             (pixels; default: whole map in Search mode, 1/2 template in Refine mode).",
"-resolution 15.6,200         Resolution limits for cross-correlation (default 0).",
"-symmetry I                  Symmetry of the template (default: no symmetry, thetaphilimit used instead).",
"-bin 3                       Binning of input, template and mask by given kernel size,",
"                             output image remains same size as input.",
"-refinepeaks                 Refine peaks from cc search to subpixel accuracy.",
" ",
"Parameters for search:",
"-threshold 0.5               Report only particles above the ccc threshold.",
"-maxhits 10                  Maximum number of particles to report for each orientation (default: all).",
" ",
"Parameters for refinement:",
"-zshiftlimit 10              Shift allowed in addition to 'shiftlimit' in the direction of the view vector (Use only in Refine mode, no xyshiftlimit).",
"-xyshiftlimit 10             Shift allowed in addition to 'shiftlimit' in the orthogonal to the view vector (Use only in Refine mode, no zshiftlimit).",
" ",
"Parameters for particle selection:",
"-all                         Reset selection to all particles before other selections.",
"-mindist 64                  Discard overlapping particles (default: size of the template in Search, selection number set to zero).",
"-minangledif 50              Discard particles with too similar orientations (selection number set to zero).",
"-fom 0.3                     Discard particles below the threshold (selection number set to zero).",
"-spheremindist 5             Discard particles too far from a sphere that best fits all the particles (selection number set to zero).",
"-topp 75.0                   Select only the best percentage of particles (default 100%).",
"-topn 100                    Select only a number of best particles from each reconstruction (default all).",
"-withindist 100,100,100,45   Select only particles within a set distance from a point.",
"-remove                      Remove all particles with selection number 0 after other selections.",
"-renumber                    Renumber particles.",
" ",
"Parameters for particle extraction:",
"-select 1                    Extract particles with this selection number (default: all > 0).",
"-size 50,50,50               New size for transformed particles (use only in Extract mode. one value sets all).",
"-partpath dir/subdir         Set the particle file paths (use only in Extract mode).",
"-partextension pif           Set the particle image file format (dafault: map. use only in Extract mode).",
" ",
"Input:",
"-Template image.map          Template to search for or to refine.",
"-mask mask_1.map             Mask in real space to apply on template. Use a soft edged mask.",
"-Mask mask_2.map             Mask in reciprocal space to apply during cross-correlation and for CCC normalization.",
" ",
"Output:",
"-output file.star            Output parameter file with orientations.",
"-ccmax ccmax.img             Output image with maximum cross-correlation coefficient value for each position (use only in Search mode).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Bstring			mode;		  				// Mode, default is none
	Bstring			rec_id;					// Reconstruction ID
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Vector3<double>	origin;					// Template origin
	Vector3<double>	origin2;				// Transformation origin
	Vector3<double>	origin3;				// Wedge mask origin
	int				set_origin(0);			// Flag to set origin
	int				set_sampling(0);		// Flag to set sampling
	double			alpha_step(M_PI_4);		// Angular step size for alpha
	double			theta_step(M_PI_4);		// Angular step size for theta
	double			phi_step(M_PI_4);		// Angular step size for phi
	double			alpha_limit(M_PI);		// Angular limit of alpha in search
	double			thetaphi_limit(M_PI);	// Angular limit of theta and phi in search
	double			hires(0), lores(0);		// Limiting resolution for cross-correlation
	double			shiftlimit(-1);			// Use default search radius
	double			shiftlimitz(0);			
	double			shiftlimitxy(0);			
	View			ref_view;				// Reference view
	Quaternion		quat;					
	int				applytransform(0);		// Rotate and shift in extraction (default not)
	int				extractparticle(1);		// 0 = template, 1 = particle
	int 			refinepeaks(0);			// No refinement of the cc peaks to subpixel accuracy
	Vector3<long>	bin = {1,1,1};			// No binning by default
	Bstring			template_file;			// Reference template
	Bstring			mask_file;				// Reciprocal space mask
	Bstring			mask_file2;				// Real space mask
	Bstring			ps_file;				// Postscript output
	Bstring			param_file;				// Parameter file
	Bstring			ccmax_file;				// CCC max value image file
	double			threshold(0.005);		// Threshold for cross correlation search
	int				maxhits(-1);			// Maximum nimber of hits in cross correlation search
	Vector3<long>	size;					// Boxsize for extracting particles
	Vector3<long>	size2;					// Boxsize for extracting particles
	Bstring			partpath;				// Particle file path
	Bstring			part_ext;				// Particle image file format
	double			randomlocations(0);		// Standard deviation for origin errors (angstrom)
	double			randomviews(0);			// Standard deviation for view errors (degrees)
	double			randommag(0);			// Standard deviation for magnification errors (fraction)
	int				iters(1);				// Number of iterations
	double			mindist(-1);			// Minimum distance between positions
	double			peakdist;
	double			minangledif(-1);		// Minimum difference in the angles
	double			spheremindist(-1);		// Minimum distance to a sphere fitted to the particles
	double			angledif;	
	double			fom_cutoff(-1);
	double			fom_percentage(-1);		// Percentage particles to accept
	int				num_to_select(-1);	 	// Number of particles to accept
	double			withindist_x(0);
	double			withindist_y(0);
	double			withindist_z(0);
	Vector3<double>	withindist_point;
	double			withindist_d(-1);
	int 			reset(0);				// Keep selection as read from file
	int				select(-1);
	int				renumber(0);
	int				convert(0);				// Promotes a subtomogram to a tomogram
	int 			remove(0);
	int				extractstar(0);
	int				extractwedges(0);
	int				ipart;
	int				i, optind;
	Bsymmetry   	sym;
	
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
	  	if ( curropt->tag == "mode" ) {
			mode = curropt->value.lower();
	  		if ( mode.length() < 1 ) {
	  			cerr << "-mode: A mode must be specified" << endl;
				bexit(-1);
			}
	  	}
		if ( curropt->tag == "inverse" )
			extractparticle = 0;
		if ( curropt->tag == "transform")
			applytransform = 1;
		if ( curropt->tag == "star")
			extractstar = 1;
		if ( curropt->tag == "wedges")
			extractwedges = 1;
		if ( curropt->tag == "sampling") {
			sam = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "origin") {
			origin = curropt->origin();
			set_origin = 1;
		}
		if ( curropt->tag == "angles") {
			if ( ( i = curropt->values(alpha_step, theta_step, phi_step) ) < 1 ) {
				cerr << "-angles: An angle step size must be specified" << endl;
				bexit(-1);
			}
			else {
				if ( i < 2 ) phi_step = theta_step = alpha_step;
				else if ( i < 3 ) phi_step = theta_step;

				if ( alpha_step < 0.01 ) { alpha_step=1; }
				if ( theta_step < 0.01 ) { theta_step=1; }
				if ( phi_step < 0.01 ) { phi_step=1; }

				alpha_step *= M_PI/180.0;
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "alphalimit" ) {
			if ( ( alpha_limit = curropt->value.real() ) < 0.001 ) {
				cerr << "-alphalimit: An angular limit for alpha must be specified!" << endl;
				bexit(-1);
			}
			else {
				alpha_limit *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "thetaphilimit" ) {
			if ( ( thetaphi_limit = curropt->value.real() ) < 0.001 ) {
				cerr << "-thetaphilimit: An angular limit for theta & phi must be specified!" << endl;
				bexit(-1);
			}
			else {
				thetaphi_limit *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(hires, lores) < 1 ) {
				cerr << "-resolution: A resolution must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "shiftlimit" ) {
			if ( ( shiftlimit = curropt->value.real() ) < 0.001 ) {
				cerr << "-shiftlimit: A shift search radius must be specified!" << endl;
				bexit(-1);
			}
			else {
				if ( shiftlimit < 0 ) { shiftlimit = 0; }
			}
		}
		if ( curropt->tag == "zshiftlimit") {
			if ( ( shiftlimitz = curropt->value.real() ) < 0.001 ) {
				cerr << "-zshiftlimit: A maximum shift in the direction of the view vector must be specified!" << endl;
				bexit(-1);
			}
			else {
				if ( shiftlimitz < 0 ) { shiftlimitz = 0; }
			}
		}
		if ( curropt->tag == "xyshiftlimit") {
			if ( ( shiftlimitxy = curropt->value.real() ) < 0.001 ) {
				cerr << "-xyshiftlimit: A maximum shift orthogonal to the view vector must be specified!" << endl;
				bexit(-1);
			}
			else {
				if ( shiftlimitxy < 0 ) { shiftlimitxy = 0; }
			}
		}
		if ( curropt->tag == "mindist") {
			if ( ( mindist = curropt->value.real() ) < 0.001 ) {
				cerr << "-mindist: A minimum distance between particles must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "spheremindist") {
			if ( ( spheremindist = curropt->value.real() ) < 0.001 ) {
				cerr << "-spheremindist: A minimum distance between particles and a fitted sphere!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "minangledif") {
			if ( ( minangledif = curropt->value.real() ) < 0.001 ) {
				cerr << "-minangledif: A minimum difference in the orientations between particles must be specified!" << endl;
				bexit(-1);
			}
			else {
				minangledif *= M_PI/180.0;
			}
			
		}
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "bin" ) {
			if ( ( bin[0] = curropt->value.integer() ) < 1 ) {
				cerr << "-bin: A bin size must be specified!" << endl;
				bexit(-1);
			} else {
				bin[2] = bin[1] = bin[0];
			}
		}
		if ( curropt->tag == "refinepeaks" ) {
			refinepeaks = 1;
		}
		if ( curropt->tag == "threshold" ) {
			if ( ( threshold = curropt->value.real() ) < 0.001 ) {
				cerr << "-threshold: A threshold must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "topp" ) {
			if ( ( fom_percentage = curropt->value.real() ) < 0.001 ) {
				cerr << "-topp: A percentage must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "topn" ) {
			if ( ( num_to_select = curropt->value.integer() ) < 1 ) {
				cerr << "-topn: A number must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "withindist" ) {
			if ( ( i = curropt->values(withindist_x, withindist_y, withindist_z, withindist_d) ) < 4 ) {
				cerr << "-withindist: Three coordinates and a maximum distance (in pixels) must be specified" << endl;
				bexit(-1);
			} else {
				withindist_point = Vector3<double>(withindist_x, withindist_y, withindist_z);
			}
		}
		if ( curropt->tag == "maxhits" ) {
			if ( ( maxhits = curropt->value.integer() ) < 1 ) {
				cerr << "-maxhits: The maximum number of hits must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "Template" ) {
			template_file = curropt->filename();
		}
		if ( curropt->tag == "mask" ) {
			mask_file2 = curropt->filename();
		}
		if ( curropt->tag == "Mask" ) {
			mask_file = curropt->filename();
		}
		if ( curropt->tag == "output" ) {
			param_file = curropt->filename();
		}
		if ( curropt->tag == "ccmax" ) {
			ccmax_file = curropt->filename();

		}
		if ( curropt->tag == "size" ) {
			size = curropt->size();
			if ( size.volume() < 1 ) {
				cerr << "-size: A box size must be specified." << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else {
				if ( partpath[-1] != '/' ) partpath += "/";
			}
		}
		if ( curropt->tag == "partextension" ) {
			part_ext = curropt->value;
			if ( part_ext.length() < 1 ) {
				cerr << "-partextension: The particle file extension must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "randomlocations" ) {
			if ( ( randomlocations = curropt->value.real() ) < 0.001 ) {
				cerr << "-randomlocations: A standard deviation must be specified!" << endl;
				bexit(-1);
			}
			else if ( randomlocations < 0 ) randomlocations = 0;
		}
		if ( curropt->tag == "randomviews" ) {
			if ( ( randomviews = curropt->value.real() ) < 0.001 )
				cerr << "-randomviews: A standard deviation must be specified!" << endl;
			else {
				randomviews *= M_PI/180.0;		// Assume degrees
				if ( randomviews < 0 ) randomviews = 0;
			}
		}
		if ( curropt->tag == "randommag" ) {
			if ( ( randommag = curropt->value.real() ) < 0.00001 ) {
				cerr << "-randommag: A standard deviation must be specified!" << endl;
				bexit(-1);
			}
			else {
				if ( randommag < 0 ) randommag = 0;
			}
		}
		if ( curropt->tag == "iterate" ) {
			if ( ( iters = curropt->value.integer() ) < 1 ) {
				cerr << "-iterate: The number of iterations must be specified!" << endl;
				bexit(-1);
			}
			else {
				if ( iters < 1 ) { iters = 1; }
			}	
		}

		if ( curropt->tag == "renumber" )
	        	renumber = 1;
		if ( curropt->tag == "fom" ) {
			if ( ( fom_cutoff = curropt->value.real() ) < 0.001 ) {
				cerr << "-fom: A threshold value must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "all" )	{
			reset = 1;
		}
		if ( curropt->tag == "remove" ) {
			remove = 1;
		}
		if ( curropt->tag == "select" ) {
			if ( ( select = curropt->value.integer() ) < 1 ) {
				cerr << "-select: A selection number must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "convert" ) {
	        	convert = 1;
		}
    	}
	option_kill(option);
	
	double		ti = timer_start();
	
	Vector3<double>		currshift, bestshift;
	
	Bimage*				p = NULL;
	Bimage*				p2 = NULL;
	Bimage*				p3 = NULL;
	Bimage*				ptemp = NULL;
	Bimage*				pmask = NULL;
	Bimage*				pmask2 = NULL;
	Bimage*				pmask3 = NULL;
	Bstring				filename;
	Bstring				starfilename;
	Bstring				wedgefilename;
	Vector3<double>		scale(1,1,1);
	Bproject*			project = NULL;
	Bproject*			partproject = NULL;
	Breconstruction*	rec = NULL;
	Breconstruction*	newrec = NULL;
	Bparticle*			part = NULL;
	Bparticle*			part2 = NULL;
	Bparticle*			newpart = NULL;
	Vector3<long>		psize, ptempsize;	
	long				npart = 0;
	Matrix3				mat;
	Vector3<double>		translate;
	Vector3<long>		obin(bin[0], bin[1], bin[2]);

	// Read template
	if ( template_file.length() ) {
		if ( verbose & VERB_PROCESS ) { cout << "Reading the template." << endl << endl; } 
		ptemp = read_img(template_file, 1, 0);
		ptemp->change_type(Float);
		ptempsize = ptemp->size();
	}

	// Template required for these modes
	if ( ( (mode == "search") || (mode == "refine") ) && ptemp == NULL ) {
		cerr << "Error: No template file specified!" << endl << endl;
		bexit(-1);
	}

	// Template required for search and extract mode with inverse option
	if ( (mode == "extract") && ( extractparticle == 0) && ptemp == NULL ) {
		cerr << "Error: No template file specified!" << endl << endl;
		bexit(-1);
	}

	// Read all the parameter files (or reconstructions)
	if ( verbose & VERB_PROCESS ) { cout << "Reading the input file(s)." << endl << endl; } 
	Bstring* file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl << endl;
		bexit(-1);
	}
	if ( file_type(*file_list) == Micrograph ) {

		project = read_project(file_list);
		
		// Promote a subtomogram to a tomogram
		if ( convert ) {
			for ( rec=project->rec; rec; rec=rec->next ) {
				rec->frec = rec->part->fpart;
				rec->part->fpart = NULL;
				rec->fpart = NULL;
				rec->symmetry = 0;
				// old origin becomes new location
				rec->part->loc = rec->part->ori;
				// new origin
				rec->part->ori = origin;
			}
		}
	}
	else {
		project = project_create_from_images(file_list, 0);
		for ( rec=project->rec; rec; rec=rec->next ) { // Remove the particle record which was created unnecessarily
			rec->fpart=NULL;
			if ( ptemp ) {
				rec->box_size = ptempsize;
			}
		}
	}

	string_kill(file_list);

	// Read reciprocal space mask (missing wedge mask)
	if ( mask_file.length() ) {
		if ( verbose & VERB_PROCESS ) { cout << "Reading the reciprocal space mask." << endl << endl; } 
		pmask = read_img(mask_file, 1, 0);
		if (!( pmask )) {
			cerr << "Error: No mask found!" << endl;
			bexit (-1);
		}
	}
	if ( (mode == "extract") && ( extractwedges == 1 ) && pmask == NULL ) {
		cerr << "Error: No reciprocal space mask file specified!" << endl << endl;
		bexit(-1);
	}

	// Read real space mask
	if ( mask_file2.length() ) {
		if ( verbose & VERB_PROCESS ) { cout << "Reading the real space mask." << endl << endl; } 
		pmask2 = read_img(mask_file2, 1, 0);
		if (!( pmask2 )) {
			cerr << "Error: No mask found!" << endl;
			bexit (-1);
		}
	}

	// Set the origin and bin the template
	if ( ptemp ) {
		if ( set_origin ) { ptemp->origin(origin); }
		if ( set_sampling ) { ptemp->sampling(sam); }
		if ( bin[0] > 1 ) {
			if ( verbose & VERB_PROCESS ) { cout << "Binning the template." << endl << endl; } 
			ptemp->bin(bin);
		}
		ptempsize = ptemp->size();
	}

	// Set the origin and bin the real space mask
	if ( pmask2 ) {
		if ( set_origin ) { pmask2->origin(origin); }
		if ( set_sampling ) { pmask2->sampling(sam); }
		if ( bin[0] > 1 ) {
			if ( verbose & VERB_PROCESS ) { cout << "Binning the real space mask." << endl << endl; } 
			pmask2->bin(bin);
		}
	}

	// Real space mask must be the same size as template
	if ( pmask2 && ptemp ) {
		if (ptemp->size() != pmask2->size()) {
			cerr << "Error: Real space mask and template must be of the same size!" << endl << endl;
			bexit(-1);
		}
	}

	if ( mode != "refine" ) {
		shiftlimitz = 0;
		shiftlimitxy = 0;
	}

	if ( shiftlimitz && shiftlimitxy ) {
		cerr << "Error: Both zshiftlimit and xyshiftlimit specified!" << endl << endl;
	}

	// Take binning into account in distance limits
	if ( bin[0] > 1 ) {
		mindist /= bin[0];
		shiftlimit /= bin[0];
		shiftlimitz /= bin[0];
		shiftlimitxy /= bin[0];
	}

	// Run in given mode for all the reconstructions in the project
	for ( rec=project->rec; rec; rec=rec->next ) {

		// Set origin of all particles
		if ( set_origin == 1 ) {
			if ( verbose & VERB_PROCESS )
				cout << endl << "Setting particle origins." << endl;
			for (part = rec->part; part; part=part->next ) { part->ori = origin;	}		
		}

		// Read the reconstruction
		if ( (mode == "search") || (mode == "refine") || ( (mode == "extract") && (extractparticle == 1)) ) {
			if ( verbose & VERB_PROCESS )
				cout << endl << "Reading the reconstruction." << endl << endl;
			p = read_img(rec->frec, 1, 0);
			if ( !p ) {
				cerr << "Error: No reconstruction found!" << endl;
				bexit (-1);
			}
			if ( set_sampling ) { p->sampling(sam); }
			if ( bin[0] > 1 ) { p->bin(bin); }
			p->change_type(Float);
		}

		// SEARCH MODE

		if ( mode == "search" ) {
			npart = reconstruction_search_subtomo(rec, p, ptemp, pmask, pmask2,
				alpha_step, theta_step, phi_step, alpha_limit, thetaphi_limit,
				hires, lores, shiftlimit, mindist,
				threshold, maxhits, bin, sym, refinepeaks, ccmax_file);
			
			if ( verbose ) {
				if ( npart > 0 )
					cout << endl << "Found " << npart << " particle(s), some of which maybe overlapping each other." << endl;
				else
					cout << endl << "No particles found." << endl;
			}
		}

		// REFINE MODE		

		if ( mode == "refine" ) {
			npart = reconstruction_refine_subtomo(rec, p, ptemp, pmask, pmask2,
						alpha_step, theta_step, phi_step, alpha_limit, thetaphi_limit,
						hires, lores, shiftlimit, shiftlimitz, shiftlimitxy,
						mindist, bin, sym, iters, refinepeaks, ccmax_file);

			if ( verbose & VERB_PROCESS )
				cout << "Refined " << npart << " particle(s)." << endl;
		}

		// Normalize particle views automatically

		if ( verbose & VERB_PROCESS )
			cout << endl << "Normalizing particle views." << endl;
		for ( part = rec->part; part; part=part->next ) part->view.normalize();
		
		// PARTICLE SELECTIONS

		// Select all particles
		if ( reset == 1 ) {
			if ( verbose & VERB_PROCESS )
				cout << endl << "Resetting particle selections." << endl;
			for ( part = rec->part; part; part=part->next ) part->sel = 1;
		}

		// Determine the ccc threshold from percentage
		if ( (fom_percentage >= 0) || (num_to_select >= 0) ) {
			int_float* 	rank = part_fom_order(rec->part, npart, 0);

			if (fom_percentage >= 0) { num_to_select = (int) (npart*fom_percentage/100.0); }
			if (num_to_select > npart) { num_to_select = npart; }

			fom_cutoff = rank[num_to_select].f;

			delete[] rank;
		}

		// Discard below the ccc threshold
		if ( fom_cutoff > 0 ) {
			if ( verbose & VERB_PROCESS ) {	cout << endl << "Discarding particles below threshold." << endl; }		
			for (part = rec->part; part; part=part->next ) {
				if ( part->fom[0] <= fom_cutoff ) { part->sel = 0; }
			}		
		}

		// Discard overlapping particles
		if ( mindist > 0 ) {

			if ( verbose & VERB_PROCESS ) {	cout << endl << "Selecting the best particle for each position." << endl; }
		
			// for each particle, check if all its followers are worse. If not, unselect
			part = rec->part;
			while ( part ) {
				if ( part->sel != 0 ) {
					part2 = part->next;
					while ( part2 ) {
						if ( part2->sel != 0 ) {
							peakdist = (part->loc-part2->loc).length();
							// mindist in binned units
							if ( peakdist < mindist * bin[0] ) {
								if ( part2->fom[0] > part->fom[0] ) {
									part->sel = 0;
								} else {
									part2->sel = 0;
								}
							}
						}
						part2 = part2->next;
					}
				}
				part = part->next;
			}
			if ( verbose & VERB_PROCESS ) { cout << endl << "Minimum distance between particles: " << mindist * bin[0] << endl; }
		}

		// Discard particles with too similar views
		Vector3<double>	part_view;
		Vector3<double>	part2_view;
		View		view;
		View		view2;

		if ( minangledif >  0 ) {

			if ( verbose & VERB_PROCESS ) {	cout << endl << "Selecting the best particle for each view." << endl; } 
		
			// for each particle, check if all its followers are worse. If not, unselect
			part = rec->part;
			while ( part ) {
				if ( part->sel != 0 ) {
					part2 = part->next;			
					while ( part2 ) {
						if ( part2->sel != 0 ) {
							view = part->view;
							view2 = part2->view;

							part_view = Vector3<double>(view[0],view[1],view[2]);
							part2_view = Vector3<double>(view2[0],view2[1],view2[2]);

							angledif = part_view.angle(part2_view);

							// cout << "part %d--part %d: " << part->id, part2->id, (angledif*180/M_PI));
							if ( angledif < minangledif ) {
								if ( part2->fom[0] > part->fom[0] ) {
									part->sel = 0;
									cout << "unselected part: " << part->id;
								} else {
									part2->sel = 0;
									cout << "unselected part: " << part2->id;
								}
							}
						}
						part2 = part2->next;
					}
				}
				part = part->next;
			}
			if ( verbose & VERB_PROCESS ) {
				cout << endl << "Minimum difference in the view angles between particles: " << (minangledif*180/M_PI) << endl;
			}
		}

		// Discard particles too far from a specified point
		if ( withindist_d > 0 ) {
			if ( verbose & VERB_PROCESS ) {	cout << endl << "Discarding particles further than " << withindist_d << " from the given point." << endl; }
			for (part = rec->part; part; part=part->next ) {
				cout << "ID " << part->id << ": distance: " << (part->loc-withindist_point).length() << endl;
				if ( ((part->loc-withindist_point).length()) > withindist_d ) { part->sel = 0; }
			}		
		}

		// Discard particles too far from a fitted sphere surface
		if ( spheremindist > 0 ) {
			Sphere sphere = {0,0,0,0};
			double	 distance;
			double	 sumdistances2(0);
			if ( verbose & VERB_PROCESS ) {	cout << endl << "Fitting a sphere to particles locations" << endl; }
			sphere = locations_fit_sphere(rec->part, 100, 0.1);
			if ( verbose & VERB_PROCESS ) {	cout << "Found a sphere with center " << sphere.A << " " << sphere.B << " " << sphere.C << " and radius " << sphere.R << endl; }
			if ( verbose & VERB_PROCESS ) {	cout << "Discarding particles further than " << spheremindist << " pix from the surface of the sphere." << endl; }
			for (part = rec->part; part; part=part->next ) {
				distance = (part->loc-Vector3<double>(sphere.A,sphere.B,sphere.C)).length()-sphere.R;
				sumdistances2 += distance*distance;
				if ( verbose & VERB_PROCESS ) { cout << "ID " << part->id << ": distance: " << distance << endl;	}
				if ( fabs(distance) > spheremindist ) { part->sel = 0; }
			}
			if ( verbose & VERB_PROCESS ) {	cout << "RMSD of all locations: " << sqrt(sumdistances2)/particle_count(rec->part) << " pix" << endl; }
		}

		// remove unselected particles if wanted (currently removed particle structures are not killed)
		if ( remove ) {
			part2=rec->part;
	
			while ( part2 ) {
				if ( part2->sel == 0 ) {
					if ( verbose & VERB_PROCESS ) {
						cout << "Removing particle: " << part2->id << endl;
					}
					if ( part ) {
						part->next = part2->next;
						part2=part->next;
					} 
					else {
						rec->part = part2->next;
						part2 = part2->next;
					}
				
				} else {
					part = part2;
					part2 = part2->next;
				}
			}
		}

		if ( renumber ) {
			for (part = rec->part, ipart=1; part; part=part->next, ipart++) {
				part->id = ipart;
			}
		}
		
		if ( verbose ) project_show_selected(project);

		// EXTRACT MODE

		if ( mode == "extract" ) {
						
			part = rec->part;
			npart = particle_count_selected(part);

			if ( part_ext.length() < 1 ) { part_ext = "map"; }

			if ( verbose & VERB_PROCESS ) {	cout << endl << "Extracting " << npart << " particle(s)." << endl; }

			if ( !applytransform ) { cerr << "Warning: -transform option not give given, particles will not be rotated, only shifted." << endl; }

			// give box size is relative to unbinned units
			size = size / obin;

			if ( extractwedges ) {
				if ( size.volume() == 0 ) {
					// size of the template by default
					if ( ptemp ) { size = ptemp->size(); }
					else { size = p->size(); }
				}

				// Scale the reciprocal space mask to match the size of the extracted particle
				if ( pmask->size() != size ) {
					if ( verbose & VERB_PROCESS ) {	cout << endl << "Scaling the reciprocal space mask to the size of the extracted particle." << endl; } 
//					img_transform(pmask, size, size/pmask->size(),
  //            		    Vector3<double>(0,0,0), Vector3<double>(0,0,0), ref_view.matrix(), FILL_BACKGROUND, 0);
					pmask->mask_fspace_resize(size);
				}
			}
	
			// extract the images
			filename = rec->fpart;
			if ( filename.length() < 3 )
				filename = rec->frec.pre_rev('.') + "_part." + part_ext;
			
			particle_setup_filenames(rec->part, filename, partpath);

			for ( part = rec->part; part; part = part-> next ) {

				if ( ((select < 0) && (part->sel > 0)) || (part->sel == select) ) {

					view = part->view;
					view.normalize();
				
					if ( extractparticle == 1 ) {		// Extract the particle (and the corresponding wedge mask)
//						p2 = p->copy();

						if ( ( size.volume() == 0 ) ) {
							// size of the template by default
							if ( ptemp ) { size = ptemp->size(); }
							else { size = p->size(); }
						}

//						origin2 = part -> loc / Vector3<double>( bin[0],bin[1],bin[2] );
						origin2 = part->loc / obin;
//						p2->origin(origin2);
						//pmask2->origin(origin2);

//						translate = ( part->ori - part->loc ) / Vector3<double>(bin[0],bin[1],bin[2]);
						translate = ( part->ori - part->loc ) / obin;

						// keep extracted particle in the center
//						translate = ( translate - ( part->ori - Vector3<double>(size[0],size[1],size[2]) / 2 ) ) / Vector3<double>(bin[0],bin[1],bin[2]);
//						translate = ( translate - ( part->ori - size) / 2 ) / bin;

						if ( applytransform ) { 				
							mat = view.matrix();
						} else {
							mat = ref_view.matrix();
						}
	
//   						p2 = img_transform_copy(p, size, scale, origin2, translate, mat, FILL_BACKGROUND, 0);
   						p2 = p->transform(size, scale, origin2, translate, mat, FILL_BACKGROUND, 0);

						if ( pmask2 ) {
							if ( verbose & VERB_PROCESS ) { cout << endl << "Masking the extracted particle." << endl; }
							p2->multiply(pmask2);
						}

						if ( extractwedges ) {
//							pmask3 = pmask->copy();
//							origin3 = pmask3->size()/2;
//							pmask3->origin(origin3);
//							img_transform(pmask3, size, scale, origin3, Vector3<double>(0,0,0), mat, FILL_BACKGROUND, 0);
							origin3 = pmask->size()/2;
							pmask3 = pmask->transform(size, scale, origin3, Vector3<double>(0,0,0), mat, FILL_BACKGROUND, 0);
							pmask3->origin(origin3);
							wedgefilename = part->fpart.pre_rev('.') + "_mask." + part_ext;
							pmask3->change_type(nudatatype);
							write_img(wedgefilename, pmask3, 0);
							delete pmask3;
						}
//						p2->origin(origin);
						p2->origin(part->ori/obin);
					}

					if ( extractparticle == 0 ) {		// Transform the template instead
						if ( applytransform ) { 
							mat = view.matrix();
							mat = mat.transpose();
						} else {
							mat = ref_view.matrix();
						}

//						origin = part -> ori / Vector3<double>( bin[0],bin[1],bin[2] );
						origin = part->ori / obin;

						p2 = ptemp->copy();

						p2->origin(origin);

						if ( pmask2 ) {
							if ( verbose & VERB_PROCESS ) { cout << "Masking the template before transformation." << endl; }
							p2->multiply(pmask2);
						}

						// If no size is specified no translation is applied and the origin will reflect the location of the particle
						// in the reconstruction
	
						// If size is specified, translation will be applied and origin is set to zero

						if ( size.volume() == 0 ) {				
//	   						img_transform(p2, ptempsize, scale, origin, Vector3<double>(0,0,0), mat, FILL_BACKGROUND, 0);
	   						p2->transform(scale, origin, Vector3<double>(0,0,0), mat, FILL_BACKGROUND, 0);
//							origin = -(part->loc - part->ori) / Vector3<double>(bin[0],bin[1],bin[2]);
							origin = -(part->loc - part->ori) / obin;
						} else {
//							translate = (part->loc - part->ori) / Vector3<double>(bin[0],bin[1],bin[2]);
							translate = (part->loc - part->ori) / obin;
//		   					img_transform(p2, size, scale, origin, translate, mat, FILL_BACKGROUND, 0);
		   					p3 = p2->transform(size, scale, origin, translate, mat, FILL_BACKGROUND, 0);
							delete p2;
							p2 = p3;
//							origin = Vector3<double>(0,0,0);
							origin = p2->size()/2;
						}

						p2->origin(origin);
					}

					p2->change_type(nudatatype);
					write_img(part->fpart, p2, 0);
					delete p2;
	
					if ( extractstar ) {
						partproject = project_create(0, 1);				
	
						newrec = partproject->rec;
						newrec->id=part->fpart.pre_rev('.');
						newrec->frec=part->fpart;
						newrec->voxel_size=rec->voxel_size;
						newrec->box_size=rec->box_size;
						newpart = particle_add(&newrec->part, 1);

						if ( set_origin == 1 ) {
							newpart->loc = origin;
							newpart->ori = origin;
						}
						else if ( size[0] * size[1] * size[2] > 0 ) {				
							newpart->loc = Vector3<double>(size[0]/2, size[1]/2, size[2]/2);
							newpart->ori = Vector3<double>(size[0]/2, size[1]/2, size[2]/2);
						}
						else {
							newpart->loc = part->loc;
							newpart->ori = part->ori;
						}

						if ( applytransform ) {
							newpart->view = ref_view;
						}
						else {
							newpart->view = view;
						}
	
						newpart->fom[0] = part->fom[0];
	
						starfilename = part->fpart.pre_rev('.') + ".star";

						if ( starfilename.length() ) { write_project(starfilename, partproject, 0, 0); }

						project_kill(partproject);
					}
				} 
			} 
		} 

		delete p;
		delete pmask2;
	} 

	delete pmask;
	delete ptemp;

	if ( param_file.length() )
		write_project(param_file, project, 0, 0);

	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
