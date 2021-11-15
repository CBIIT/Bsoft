/**
@file	blocres.cpp
@brief	Local resolution analysis
@author Giovanni Cardone and Bernard Heymann
@date	Created: 20070611
@date	Modified: 20171122 (BH)
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
"Usage: blocres [options] input1.img input2.img output.img",
"---------------------------------------------------------",
"The local (-box) and shell (-shell) resolution calculations are mutually exclusive.",
" ",
"Actions:",
"-box 20                  Kernel size for determining local resolution (pixels/voxels).",
"-shell 14                Shell width for determining radial resolution (pixels/voxels).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-step 5                  Interval between voxel samples or shells for resolution analysis (pixels, default: 1).",
"-radii 20,130            Radial limits for local or shell resolution (default 0,<x-dimension>).",
"-maxresolution 5         Maximum frequency available in the data (angstrom).",
"-cutoff 0.5,0.143        Resolution cutoffs for FSC (default 0.5).",
"-fill 200                Value to fill the background (non-masked regions; default 0).",
" ",
"Parameters for local resolution:",
"-pad 2                   Resolution box padding factor (0 = none, default: 1 (box) and 0 (shell)).",
"-edge 10                 Pixels at the edge of the whole map to exclude from analysis (default: half the size of the resolution box).",
"-taper edge              Apply a windowing aperture to box volumes before resolution analysis.",
"                           none = do not taper",
"                           edge = circular mask with smooth edges",
"                           hann = hanning window (default)",
"-nofill                  Do not fill values in between (if step > 1).",
"-symmetry D5             Point group symmetry.",
" ",
"Parameters for shell resolution:",
"-smooth                  Smooth the shell edge.",
" ",
"Input:",
"-Mask mask.tif,2         Mask file to use for limiting the analysis to a defined region",
"                           and level value to use (optional, default: all but zero).",
" ",
"Output:",
"-usedmask used.mrc       Output file name for mask after modifications.",
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
	int				shell_width(0);				// Shell thickness
	int				taper(2);					// Flag to taper boxes
	int				smooth(0);					// Flag to smooth the shell edge
	Vector3<long> 	vedge(-1,-1,-1);			// Edge size to exclude
	int				nofill(0);                	// Fill flag
	int				pad(-1);        	        // Padding factor
	double			maxres(0);					// Resolution limit
	double			cutoff[4] = {0.5,0,0,0};	// FSC cutoff values
	int				minrad(0);					// Minimum radius for shell resolution
	int				maxrad(0);					// Maximum radius for shell resolution
	int				setimg(-1);					// Select all images
	long			step(1);					// Sampling step
	Bstring			symmetry_string;			// Point group
	double			fill(0);					// Background fill value
	Bstring			maskfile;					// Input mask file name
	int				mask_level(-1);				// Mask level (-1 => all but zero)
	Bstring			usedmaskfile;				// Output mask file name after modifications to original mask
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "box" )
			if ( ( boxsize = curropt->value.integer() ) < 1 )
				cerr << "-box: A kernel box size must be specified!" << endl;
		if ( curropt->tag == "shell" )
			if ( ( shell_width = curropt->value.integer() ) < 1 )
				cerr << "-shell: A shell width must be specified!" << endl;
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
		if ( curropt->tag == "step" )
			if ( ( step = curropt->value.integer() ) < 1 )
				cerr << "-step: A sampling step size must be specified!" << endl;
		if ( curropt->tag == "pad" )
			if ( ( pad = curropt->value.integer() ) < 0 )
				cerr << "-pad: A padding factor must be specified!" << endl;
		if ( curropt->tag == "maxresolution" )
			if ( ( maxres = curropt->value.real() ) < 0.1 )
				cerr << "-maxresolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "taper" ) {
			if ( curropt->value[0] == 'n' ) taper = 0;
			if ( curropt->value[0] == 'h' ) taper = 2;
			if ( curropt->value[0] == 'e' ) taper = 1;
		}
		if ( curropt->tag == "fill" )
			if ( ( fill = curropt->value.real() ) < 1 )
				cerr << "-fill: A fill value must be specified!" << endl;
		if ( curropt->tag == "nofill" )
			nofill = 1;
		if ( curropt->tag == "edge" )
			vedge = curropt->size();
		if ( curropt->tag == "radii" )
			if ( curropt->values(minrad, maxrad) < 1 )
				cerr << "-radii: At least a minimum radius must be specified!" << endl;
		if ( curropt->tag == "smooth" ) smooth = 1;
		if ( curropt->tag == "cutoff" )
			if ( curropt->values(cutoff[0], cutoff[1], cutoff[2], cutoff[3]) < 1 )
				cerr << "-cutoff: At least one FSC threshold must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "Mask" ) {
			maskfile = curropt->filename();
			if ( maskfile.contains(",") ) {
				mask_level = (maskfile.post(',')).integer();
				maskfile = maskfile.pre(',');
			}
		}
		if ( curropt->tag == "usedmask" )
			usedmaskfile = curropt->filename();
	}
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*		p1 = read_img(argv[optind++], 1, setimg);
	if ( p1 == NULL ) bexit(-1);

	Bimage*		p2 = read_img(argv[optind++], 1, setimg);
	if ( p2 == NULL ) bexit(-1);

	if ( p1->size() != p2->size() ) {
		cerr << "Error: The sizes of the two images do not match!" << endl;
		bexit(1);
	}
	
	if ( sam.volume() > 0 ) {
		p1->sampling(sam);
		p2->sampling(sam);
	}
	
	if ( set_origin ) {
		if ( set_origin == 2 ) origin = p1->default_origin();
	} else {
		origin = p1->image->origin();
		if ( origin.length() < 1 ) origin = p1->default_origin();
	}

	p1->origin(origin);
	p2->origin(origin);
	
	if ( boxsize>0 && maskfile.length() && step>1 && !nofill) {
		cerr << "ERROR: Can not set a step size more than one when using a mask!" << endl;
		cerr << "Please select only one of the two options, otherwise set the nofill option" << endl;
		bexit(1);
	}
		
	Bimage*			p = NULL;
	Bimage*			psymmask = NULL;
	Bimage*			pmask = NULL;
	Bsymmetry		sym(symmetry_string);
	
	if ( boxsize > 0 ) {
		if ( maskfile.length() ) pmask = read_img(maskfile, 1, -1);
		else {
			pmask = p1->copy_header();
			pmask->data_type(UCharacter);
			pmask->data_alloc();
			pmask->fill(1);
		}
		
		if ( maxrad )
			pmask->mask_shell(origin, minrad, maxrad);
		
		if ( sym.point() > 101 ) {
			psymmask = p1->levelmask_asymmetric_units(sym, 1);
			psymmask->mask_dilate(2);
			pmask->multiply(psymmask);
			delete psymmask;
		}
		
		if ( usedmaskfile.length() )
			write_img(usedmaskfile, pmask, 0);

		if ( pad < 0 ) pad = 1;
		p = p1->fsc_local(p2, pmask, maxres, cutoff, mask_level, boxsize,
						pad, vedge, step, taper, fill);
		if ( !p ) {
			cerr << "Error: No local resolution map produced!" << endl;
			bexit(-1);
		}
		
		if ( step > 1 && !nofill ) p->interpolate_gaps(step);

		if ( sym.point() > 101 )
			p->replicate_asymmetric_unit(sym);
	} else if ( shell_width > 0 ) {
		pad = 0;
		p = p1->fsc_shell(p2, maxres, cutoff, shell_width, step, minrad, maxrad,
						pad, smooth, fill);
		if ( !p ) {
			cerr << "Error: No shell resolution map produced!" << endl;
			bexit(-1);
		}
	}
	
	// Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}

	if (pmask) delete pmask;
	delete p1;
	delete p2;
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

