/**
@file	bsymcomp.cpp
@brief	Finds the best orientation that fits a symmetrized template.
@author	Bernard Heymann
@date	Created: 20070516
@date	Modified: 20150806
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bsymcomp [options] input.img output.img",
"----------------------------------------------",
"Finds the best orientation that fits a symmetrized template.",
"The output is the rotated input image.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 110,50,44        Origin for rotations of the template (default taken from template image).",
"-angles 8.8,2.5,5.2      Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-shift 0.5,3,-2.3        Shift to apply before symmetrization.",
"-symmetry C5             Point group symmetry.",
"-fill 0.02               Fill value (default image background).",
" ",
"Input:",
"-Template image.map      Template with symmetry to search for.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// Template origin
	int				set_origin(0);				// Flag to set origin
	double			alpha_step(M_PI_4);			// Angular step size for alpha
	double			theta_step(M_PI_4);			// Angular step size for theta
	double			phi_step(M_PI_4);			// Angular step size for phi
    Vector3<double>  shift;						// Shift to apply before symmetrization
	int 			fill_type(FILL_BACKGROUND);	// Use background
	double			fill(0);					// Default fill is zero
	Bsymmetry		sym;					// Point group
	Bstring			template_file;				// Symmmetric template
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
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
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(alpha_step, theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified" << endl;
			else {
				if ( i < 2 ) phi_step = theta_step = alpha_step;
				else if ( i < 3 ) phi_step = theta_step;
				alpha_step *= M_PI/180.0;
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "shift" )
			shift = curropt->vector3();
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "Template" )
			template_file = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

	// Read image file
	Bimage*			p = read_img(argv[optind++], 1, 0);

	if ( !p ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	Bimage*			ptemp = NULL;

	if ( template_file.length() )
		ptemp = read_img(template_file, 1, 0);

	if ( !ptemp ) {
		cerr << "Error: No template file read!" << endl;
		bexit(-1);
	}

	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();		// Preserve the old type
	
	if ( sam.volume() > 0 ) {
		p->sampling(sam);
		ptemp->sampling(sam);
	}

	if ( set_origin ) {
		p->origin(origin);
	}
	if ( set_origin ) {
		if ( set_origin == 2 ) {
			p->origin(p->size()/2);
			ptemp->origin(ptemp->size()/2);
		} else {
			p->origin(origin);
			ptemp->origin(origin);
		}
	}
		
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) {
		if ( !p->background(long(0)) ) p->calculate_background();
		fill = p->background(long(0));
	}
	
	p->background(fill);
	
	Bimage*		prot = p->find_symmetric_view(ptemp, sym,
                    phi_step, theta_step, alpha_step, shift);

	// Write an output file if a name is given
	if ( optind < argc ) {
		prot->change_type(nudatatype);
		write_img(argv[optind], prot, 0);
	}
	
	delete p;
	delete ptemp;
	delete prot;

	if ( verbose & VERB_TIME )
		timer_report(ti);

	return 0;
}

