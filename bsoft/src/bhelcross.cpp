/**
@file	bhelcross.cpp
@brief	Dealing with images with helical symmetry.
@author Bernard Heymann
@date	Created: 20121119
@date	Modified: 20160603
**/

#include "rwimg.h"
#include "symmetry.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bhelcross [options] input.img output.img",
"-----------------------------------------------",
"Calculates the 2D cross section for a helical filament.",
"The helix axis must coincide with the y axis.",
" ",
"Actions:",
"-crosssection            Calculate a cross section based on helical parameters.",
"-extrude 254             Extrude the cross section to the indicated length.",
"-symmetry C5             Apply point group symmetry.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for helix (default 0,0,0 or from image).",
"-helix 25.3,67           Helical rise and angle (default 10 angstrom, 45 degrees).",
" ",
"Parameters for calculating the cross section:",
"-resolution 15.6         High resolution limit for reconstruction (default 0).",
"-scale 1.5               Reconstruction scale (default 1).",
"-snradius 36             Radius to distinguish signal from noise (alternative to -kernel).",
"-kernel 9                Local variance kernel size for generating a foreground mask (default 11).",
" ",
"Parameters for alignment:",
"-alignresolution 20,360  High and low resolution limits (angstrom).",
"-annuli 12,60            Real space annular limits (default 0,inf pixels).",
"-shiftlimit 3.5          Limit on origin shift relative to nominal center (default 10% of box edge size).",
" ",
"Parameter for extrusion:",
"-fill 1.77               Fill value for background: average (default), background, or value.",
" ",
"Input:",
"-reference file.map      2D reference image file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			set_sampling(0);
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	int				cross(0);					// Calculate cross section
	int				extrude(0);					// Extrude cross section
	Bsymmetry		sym;					// Point group
	double			helix_rise(10);				// Rise per asymmetric unit
	double			helix_angle(M_PI/8);		// Rotation angle per asymmetric unit
	double			hires(0);					// Resolution limit
	double			scale(1);					// Reconstruction scale
	double			nuavg(0), nustd(0);			// For rescaling
	Vector3<double>	origin;			// Helix origin
	int				set_origin(0);				// Flag to set origin
	int 			calc_background(0);			// Flag to calculate background values
//	int 			norm_flag(0);				// Normalization flag
	double			snradius(0);				// Signal-noise radius
	int 			var_kernel(11);				// Variance filter kernel size
	double	 		res_lo(10000);				// Resolution limits for orientation finding
	double	 		res_hi(0);					// Must be set > 0 to determine orientations
	int 			ann_min(0);					// Minimum annulus
	int 			ann_max(1000); 				// Maximum annulus, reset to maximum radius in image
	double			shift_limit(-1);			// Maximum shift from nominal box origin
	double			fill(0);	 				// Fill value for resizing
	int 			fill_type(FILL_AVERAGE);	// Fill type for resizing
	Bstring			ref_file;					// Reference file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "crosssection" ) cross = 1;
		if ( curropt->tag == "extrude" )
			if ( ( extrude = curropt->value.integer() ) < 1 )
				cerr << "-extrude: An extrusion length must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "background" )
			calc_background = 1;
//		if ( curropt->tag == "normalize" )
//			norm_flag = 1;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" ) {
			sam = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "helix" ) {
			if ( curropt->values(helix_rise, helix_angle) < 2 )
				cerr << "-helix: Both rise and angle must be specified" << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "resolution" )
			if ( ( hires = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.0001 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "snradius" )
			if ( ( snradius = curropt->value.real() ) < 1 )
				cerr << "-snradius: A radius must be specified!" << endl;
		if ( curropt->tag == "kernel" )
			if ( ( var_kernel = curropt->value.integer() ) < 1 )
				cerr << "-kernel: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "alignresolution" ) {
			if ( curropt->values(res_hi, res_lo) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
			else
				if ( res_lo < res_hi ) swap(res_hi, res_lo);
		}
		if ( curropt->tag == "annuli" )
			if ( curropt->values(ann_min, ann_max) < 1 )
				cerr << "-annuli: At least a minimum annulus must be specified!" << endl;
		if ( curropt->tag == "shiftlimit" )
			if ( ( shift_limit = curropt->value.real() ) < 0.001 )
				cerr << "-shiftlimit: A maximum shift in pixels must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
 		if ( curropt->tag == "reference" )
			ref_file = curropt->filename();
   }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*			p = read_img(argv[optind++], 1, -1);
	Bimage*			pcs = NULL;
	Bimage*			pmask = NULL;
	
	if ( !p ) {
		cerr << "Error: No input file!" << endl;
		bexit(-1);
	}
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
	
	if ( set_sampling ) p->sampling(sam);

	double			snr(0), nsam(fabs(p->real_size()[1]*helix_angle/(TWOPI*helix_rise)));
//	cout << "nsam=" << nsam << endl;
	
	if ( calc_background ) p->calculate_background();
	
	if ( cross ) {
		pcs = p->helical_cross_section(helix_rise, helix_angle, scale, hires);
		delete p;
		p = pcs;
		if ( snradius ) {
			snr = p->snvariance(snradius);
		} else {
			pmask = p->variance_mask(var_kernel, 0.01, 1);
			snr = pmask->image->FOM();
			delete pmask;
			snr /= sqrt(nsam);
			cout << "Adjusted signal-to-noise ratio: " << snr << endl << endl;
		}
	}

	Bimage*		pref = NULL;
	if ( ref_file.length() ) {
		pref = read_img(ref_file, 1, 0);
		p->align2D(pref, ann_min, ann_max, res_lo, res_hi, shift_limit, 0);
		delete pref;
	}
	
	if ( sym.point() > 101 )
		p->symmetrize(sym, 1);
	
	if ( extrude > 1 )
		p->extrude_cross_section(extrude, helix_rise, helix_angle, fill_type, fill);

	if ( optind < argc ) {
		if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}




				
