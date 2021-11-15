/**
@file	blocfilt.cpp
@brief	Local resolution analysis
@author Giovanni Cardone and Bernard Heymann
@date	Created: 20070611
@date	Modified: 20170206 (BH)
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistence
const char* use[] = {
" ",
"Usage: blocfilt [options] input.img output.img",
"----------------------------------------------",
" ",
"Actions:",
"-box 20                  Kernel size for filtering according to local resolution (voxels).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-symmetry D5             Point group symmetry.",
"-edge 10                 Pixels at the edge of the whole map to exclude from analysis (default: half the size of the resolution box).",
" ",
"Input:",
"-Resolution locres.tif   Resolution map to use for local filtering.",
"-Mask mask.tif,2         Mask file to use for limiting the analysis to a defined region",
"                           and level value to use (optional, default: all but zero).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	int 			boxsize(0);					// Linear size of box designated by user
	Vector3<long> 	vedge(-1,-1,-1);			// Edge size to exclude
	Bstring			maskfile;					// Mask file name
	Bstring			locresfile;					// Local resolution file name
	int				mask_level(-1);				// Mask level (-1 => all but zero)
	int				setimg(-1);					// Select all images
	Bstring			symmetry_string;			// Point group

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "box" )
			if ( ( boxsize = curropt->value.integer() ) < 1 )
				cerr << "-box: A resolution box size must be specified!" << endl;
		if ( curropt->tag == "edge" )
			vedge = curropt->size();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "Mask" ) {
			maskfile = curropt->filename();
			if ( maskfile.contains(",") ) {
				mask_level = (maskfile.post(',')).integer();
				maskfile = maskfile.pre(',');
			}
		}
		if ( curropt->tag == "Resolution" )
			locresfile = curropt->filename();
	}
	option_kill(option);

	double		ti = timer_start();


	Bimage*		p1 = read_img(argv[optind++], 1, setimg);
	if ( p1 == NULL ) bexit(-1);


	if ( sam[0] > 0 && sam[1] > 0 && sam[2] > 0 ) {
		p1->sampling(sam);
	}

	if ( set_origin ) {
		if ( set_origin == 2 ) p1->origin(p1->default_origin());
		else p1->origin(origin);
	}

	if ( ! locresfile.length()) {
		cerr << "ERROR: Local resolution file needed! (option -Resolution)" << endl;
		bexit(1);
	}

	if ( ! boxsize ) {
		cerr << "ERROR: Box size value needed! (option -box)" << endl;
		bexit(1);
	}

	Bimage*			p = NULL;
	Bimage*			psymmask = NULL;
	Bimage*			pmask = NULL;
	Bimage*			resmap = NULL;
	Bsymmetry		sym(symmetry_string);

	if ( maskfile.length() ) pmask = read_img(maskfile, 1, -1);
	if ( locresfile.length() ) resmap = read_img(locresfile, 1, -1);

	if ( sym.point() > 101 ) {
		psymmask = p1->levelmask_asymmetric_units(sym, 1);
		maskfile = p1->file_name();
		maskfile = maskfile.pre_rev('.') + "_levelmask." + maskfile.post_rev('.');
		psymmask->file_name(maskfile.str());
		psymmask->mask_dilate(2);
		if ( pmask ) {
			pmask->multiply(psymmask);
			delete psymmask;
		} else {
			pmask = psymmask;
		}
		psymmask = NULL;
	}

//		write_img("locres_mask.map", pmask);

	p = p1->local_filter(pmask, mask_level, resmap, boxsize, vedge);
	if ( !p ) {
		cerr << "Error: No local-filtered map produced!" << endl;
		bexit(-1);
	}

	if ( sym.point() > 101 )
		p->replicate_asymmetric_unit(sym);

	// Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}

	if (pmask) delete pmask;
	if (resmap) delete resmap;
	delete p1;
	delete p;

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}


