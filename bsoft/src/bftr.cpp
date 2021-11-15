/**
@file	bftr.cpp
@brief	Program for resizing an image using Fourier transformation
@author Bernard Heymann
@date	Created: 20040225
@date	Modified: 20210129
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bftr [options] file.spi file.mrc",
"---------------------------------------",
"Resizing an image using Fourier transformation.",
" ",
"Actions:",
"-scale 1.56              Isotropic scale (required to resize the data).",
"-translate 0.3,-4.1,3.3  Translation.",
"-filterextremes          Remove extreme regions (such as micrograph labels)",
" ",
"Parameters:",
"-verbose 7               Verbosity of output",
"-datatype u              Force writing of a new data type",
"-origin 0,-10,30         Set origin, used with -size option (default 0,0,0)",
"-size 22,33,14           New size (voxels, default from data)",
"-sampling 2,3.5,1        Sampling before rescaling (angstrom/voxel, a single value sets all three)",
"-resolution 25.3,200     Resolution limits (angstrom).",
" ",
"Input:",
"-reference ref.map       Reference for resizing.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
    Vector3<long> 	size;					// New size
	Vector3<double> shift;					// Translation
	double			scale(0);				// Scale, 0=no scaling
	Vector3<double>	origin;					// Origin
	int				set_origin(0);			// Flag for setting the origin
	Vector3<double>	sam;    				// Sampling before rescaling
	double			res_hi(0);				// High resolution limit
	double			res_lo(1e37); 			// Low resolution limit
	int				filter_extremes(0);		// Filter extreme regions
  	Bstring			reffile;				// Reference template file name
  
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
 		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.001 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "translate" )
			shift = curropt->vector3();
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(res_hi, res_lo) < 1 )
				cerr << "-resolution: The resolution limits must be specified!" << endl;
			else
				if ( res_hi > res_lo ) swap(res_hi, res_lo);
		}
 		if ( curropt->tag == "size" )
			size = curropt->size();
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
 		if ( curropt->tag == "filterextremes" )
			filter_extremes = 1;
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();
	
    Bimage*	 		p = NULL;
	Bimage*	 		pref = NULL;
	Bimage*	 		pt = NULL;
    
	if ( ( p = read_img(argv[optind++], 1, -1) ) == NULL ) {
		cerr << "File not read!" << endl;
		exit(-1);
	}
	
	if ( reffile.length() ) {
		if ( ( pref = read_img(reffile, 1, -1) ) == NULL ) {
			cerr << "File " << reffile << " not read!" << endl;
			exit(-1);
		}
	}
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	
	if ( sam.volume() > 0 ) p->sampling(sam);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	Vector3<long>		translate;
	if ( optind < argc ) {
		if ( filter_extremes )
			p->filter_extremes();
		if ( shift.length() )
			p->fspace_translate(shift);
		if ( scale ) {
			p->fspace_resize(scale, res_hi, res_lo);
		} else if ( reffile.length() ) {
			pt = p->fspace_resize(pref);
			delete p;
			p = pt;
		}
		if ( size.volume() > 0 )
			p->resize(size, translate, FILL_BACKGROUND, 0);
		p->change_type(nudatatype);
	    write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete pref;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

