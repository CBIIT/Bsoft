/**
@file	bmedian.cpp
@brief	Program to filter images.
@author Bernard Heymann
@date	Created: 20010414
@date	Modified: 20110804
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
"Usage: bmedian [options] input.img output.img",
"---------------------------------------------",
"Filters images.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale data to average and standard deviation after filtering.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-kernel 5                Median filter kernel edge size (default 3).",
"-iterations 5            Number of iterations (default 1).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	double			nuavg(0), nustd(0); 		// Values for rescaling
	int 			median_kernel(3);			// Median filter kernel size
	int				i, iterations(1);			// Number of filtering iterations
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "kernel" )
			if ( ( median_kernel = curropt->value.integer() ) < 1 )
				cerr << "-kernel: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "iterations" )
			if ( ( iterations = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number of iterations must be specified!" << endl;
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();			// Preserve old type
	
	if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( sam.volume() ) p->sampling(sam);
	
	if (verbose )
		cout << "Iteration\tStDev" << endl;
	for ( i=0; i<iterations; i++ ) {
		p->filter_rank(median_kernel, 0.5);
		if ( verbose )
			cout << i+1 << tab << p->standard_deviation() << endl;
	}
	
	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
