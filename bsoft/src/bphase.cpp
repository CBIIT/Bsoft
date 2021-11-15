/**
@file	bphase.cpp
@brief	A program to examine phase differences between image pairs
@author Bernard Heymann
@date	Created: 20020217
@date	Modified: 20151006
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
"Usage: bphase [options] file1.img file2.img outfile.img",
"-------------------------------------------------------",
"Calculates the phase difference between two images.",
" ",
"Actions:",
"-cosine                  Calculate the cosine of the phase angle difference.",
"-amplitude               Weight the phase difference with the amplitude product.",
"-center                  Center phase difference map and output it.",
"-flip                    Flip phases based on phase difference map and modify second image.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype short          Force writing of a new data type (default float).",
"-origin 0.8,-10,15.7     Set the origin of both images (default from input images).",
"-resolution 4.5,130      Resolution (default 0 - 10000 angstrom).",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize all settings
	int				pd_type(0);					// Phase difference result type
	int 			center(0); 					// Flag to center phase difference map
	int 			flip(0);					// Flag to flip phases of second image based on phase difference map
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
    double			hires(0), lores(1e4);		// Limiting resolution range (hires must be > 0 to be set)
	    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "cosine" ) pd_type |= 1;
		if ( curropt->tag == "amplitude" ) pd_type |= 2;
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "flip" ) flip = 1;
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
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
    }
	option_kill(option);
    
	double		ti = timer_start();
	
	int 		dataflag = 0;
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p1 = read_img(argv[optind++], dataflag, -1);
	if ( p1 == NULL )  {
		cerr << "Error: No first input file read!" << endl;
		bexit(-1);
	}
	if ( !p1->data_pointer() ) {
		cerr << "Error: Data missing from image " << p1->file_name() << endl;
		bexit(-1);
	}
	
    Bimage* 	p2 = read_img(argv[optind++], dataflag, -1);
	if ( p2 == NULL )  {
		cerr << "Error: No second input file read!" << endl;
		bexit(-1);
	}
	if ( !p2->data_pointer() ) {
		cerr << "Error: Data missing from image " << p2->file_name() << endl;
		bexit(-1);
	}
	
    if ( !p1->check_if_same_size(p2) ) bexit(-1);

	if ( set_origin ) {
		if ( set_origin == 2 ) {
			p1->origin(p1->size()/2);
			p2->origin(p2->size()/2);
		} else {
			p1->origin(origin);
			p2->origin(origin);
		}
	}
	
	Bimage* 		pd = NULL;
	Bimage*			pout = NULL;
	
	pd = p1->phase_difference(p2, pd_type, hires, lores);
	pout = pd;
	
	if ( flip ) {
		p2->phase_flip(pd);
		pout = p2;
	} else {
		if ( center ) pout->center_wrap();
	}
	
    if ( optind < argc ) {
		pout->change_type(nudatatype);
    	write_img(argv[optind], pout, 0);
	}
	
	delete p1;
	delete p2;
	delete pd;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

