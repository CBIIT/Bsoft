/**
@file	bbif.cpp
@brief	Program to denoise images by bilateral filtering.
@author Bernard Heymann and Giovanni Cardone
@date	Created: 20070301
@date	Modified: 20150802
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
"Usage: bbif [options] input.img output.img",
"------------------------------------------",
"Denoises by bilateral filtering.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale data to average and standard deviation after filtering.",
"-kernel gaussian         Kernel used for for the intensity space:",
"                           1=gaussian (default)",
"                           2=lorentz ",
"                           3=tukey",
"-size 5                  kernel size (odd number, default 6*spacesigma)",
"-spacesigma 1.5          Standard deviation of spatial filter (pixels, default 1).",
"-rangesigma 23.8         Standard deviation of range filter (default image standard deviation).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;					// Units for the three axes (A/pixel)
	double			nuavg(0), nustd(0); 	// Values for rescaling
	double			spacesigma(1);			// Standard deviation of space filter
	double			rangesigma(-1);			// Standard deviation of range filter
	int				kernel_type(1);			// kernel type for range filter
											// 1=gaussian; 2=lorentz; 3=tukey biweight
	int				kernel_size(0);			// kernel size
	
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
		if ( curropt->tag == "kernel" ) {
			if ( ( kernel_type = curropt->value.integer() ) < 1 ) {
				for ( int i=0; i<(int)curropt->value.length(); i++ )
					curropt->value[i] = tolower(curropt->value[i]);
				if ( curropt->value[0] == 'g' ) {
					kernel_type = 1;
				} else if ( curropt->value[0] == 'l' ) {
					kernel_type =2;
				} else if ( curropt->value[0] == 't' ) {
					kernel_type =3;
				}
			}
			if ( kernel_type < 1 || kernel_type > 3 ) kernel_type = 1;
		}
		if ( curropt->tag == "spacesigma" )
			if ( ( spacesigma = curropt->value.real() ) < 0.001 )
				cerr << "-spacesigma: a space sigma must be specified!" << endl;
		if ( curropt->tag == "rangesigma" )
			if ( ( rangesigma = curropt->value.real() ) < 0.001 )
				cerr << "-rangesigma: a range sigma must be specified!" << endl;
		if ( curropt->tag == "size" )
			if ( ( kernel_size = curropt->value.integer() ) < 1 )
				cerr << "-size: a kernel size must be specified!" << endl;
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
	
	if ( sam.volume() ) p->sampling(sam);
	
	p->filter_bilateral(spacesigma, rangesigma, kernel_type, kernel_size/2);
			
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


