/**
@file	btube.cpp
@brief	Symmetrizes helical tubes with hexagonal lattices.
@author Bernard Heymann
@date	Created: 20131022
@date	Modified: 20131028
**/

#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: btube [options] input.img output.img",
"-------------------------------------------",
"Symmetrizes helical tubes with hexagonal lattices.",
" ",
"Actions:",
"-lattice 5,9             Lattice parameters (units along u and v vectors).",
"-normalize               Normalize symmetrized output, use for maps close to symmetry (default not).",
"-background              Calculate new background values.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for tube (default 0,0,0 or from image).",
"-constant 95             Lattice constant: distance between units in angstrom (default 1).",
"-zlimits 8,47            Range of slices along the helical axis to use (default 0,inf).",
"-radius 22               Radial limit for symmetrization in pixels (default from image).",
" ",
//"Output:",
//"-Postscript file.ps      Postscript output file name.",
//" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				h(0), k(0);						// Lattice parameters
	int 			calc_background(0);				// Flag to calculate background values
	int 			norm_flag(0);					// Normalization flag
	double			nuavg(0), nustd(0);			// For rescaling
	DataType 		nudatatype(Unknown_Type);		// Conversion to new type
	int 			set_sampling(0);
	Vector3<double> sampling;						// Units for the three axes (A/pixel)
	Vector3<double>	origin;							// Helix origin
	int				set_origin(0);					// Flag to set origin
	double			latconst(1);					// Lattice constant
	int				zmin(0), zmax(1000000);			// Limits along the helical axis
	double			radius(0);						// Radial limit for symmetrization
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "lattice" )
			if ( curropt->values(h, k) < 2 )
				cerr << "-lattice: Both lattice parameters must be specified" << endl;
		if ( curropt->tag == "normalize" ) norm_flag = 1;
		if ( curropt->tag == "background" ) calc_background = 1;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" ) {
			sampling = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "constant" )
			if ( ( latconst = curropt->value.real() ) < 0.001 )
				cerr << "-constant: A lattice constant must be specified" << endl;
		if ( curropt->tag == "zlimits" )
			if ( curropt->values(zmin, zmax) < 2 )
				cerr << "-zlimits: Both limits on the helical axis must be specified" << endl;
		if ( curropt->tag == "radius" )
			if ( ( radius = curropt->value.real() ) < 0.001 )
				cerr << "-radius: A radius must be specified" << endl;
   }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*			p = read_img(argv[optind++], 1, -1);
	
	if ( !p ) {
		cerr << "Error: No input file!" << endl;
		bexit(-1);
	}
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	if ( set_sampling ) p->sampling(sampling);
	
	if ( calc_background ) p->calculate_background();
	
	if ( h ) p->tube_symmetrize(h, k, latconst, zmin, zmax, radius, norm_flag);

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

