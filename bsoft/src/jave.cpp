/**
@file	jave.cpp
@brief	Program to average files and weigh the average in Real/Fourier space
@author Juha Huiskonen
@author	Bernard Heymann
@date	Created: 20080813
@date	Modified: 20100607
@date	Modified: 20100607
@date	Modified: 20110513 Fourier space mask must be centered to the center of the volume (shifted automatically)
@date	Modified: 20120123 (BH)
@date	Modified: 20120308
@date	Modified: 20150108 (BH) - incorporated into Bsoft
@date	Modified: 20150725
**/

#include "rwimg.h"
#include "img_combine.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"
#include "mg_subtomo.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: jave [options] input.images",
"----------------------------------",
"Averages multiple images together.",
"The average can be weighed with a weighting function,",
"either in real or Fourier space.",
" ",
"Actions:",
"-sum                     Sum input files (default is to average).",
"-rescale -0.1,5.2        Rescale input images to average and standard deviation.",
"-divide                  Values / amplitudes are divided with the weighting function (default).",
"-multiply                Values / amplitudes are multiplied with the weighting function.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
//"-size 50,50,50           Size of the output image (one value sets all).",
"-translate               Translate images to a common origin.",
"-origin 0.8,-10,15.7     Set the origin.",
"-wiener 0.1              Wiener factor for dividing (default 0.2).",
" ",
"Input:",
"-weight input.img        Real space weighting function.",
"-Weight input.img        Fourier space weighting function.",
" ",
"Output:",
"-output output.img       Output image.",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> sampling;				// Units for the three axes (A/pixel)
	Vector3<double> origin;					// New image origin
	int				set_origin(0);			// Flag to set origin
	int				divide(1), multiply(0);
	Vector3<long>	size(0,0,0);			// Boxsize for extracting particles
	Bstring			outfile;		
	Bstring			fomfile;		
	double			nuavg(0), nustd(0); 	// Average and standard deviation for rescaling
	double 			wiener(0.2);
	int				setweight(0);
	int				flags(1);				// Flags to modify summation
	Vector3<double> shift;
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype")
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling")
			sampling = curropt->scale();
		if ( curropt->tag == "origin") {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "rescale")
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified" << endl; 
		if ( curropt->tag == "wiener" )
			if ( ( wiener = curropt->value.real() ) < 1e-10 )
				cerr << "-wiener: Wiener factor larger than 0 must be specified" << endl;
		if ( curropt->tag == "Weight") {
			fomfile = curropt->filename();
			setweight = 1;
		}
		if ( curropt->tag == "weight") {
			fomfile = curropt->filename();
			setweight = 2;
		}
		if ( curropt->tag == "size") {
			size = curropt->size();
			if ( size.volume() < 1 )
				cerr << "-size: A box size must be specified." << endl;
		}
		if ( curropt->tag == "sum") {
	        	flags -= 1;
		}
		if ( curropt->tag == "divide") {
	        	divide = 1;
	        	multiply = 0;
		}
		if ( curropt->tag == "multiply") {
	        	multiply = 1;
	        	divide = 0;
		}
		if ( curropt->tag == "translate") {
	        	flags |= 8;
		}
		if ( curropt->tag == "output")
        		outfile = curropt->filename();
	}

	option_kill(option);
		
	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();
	
	if ( optind >= argc ) {
		cerr << "No input files!" << endl;
		bexit(-1);
	}
	
	// Set up image file names
	int			nfiles(0);
	Bstring*	file_list = NULL;
	while ( optind < argc ) {
		string_add(&file_list, argv[optind++]);
		nfiles++;
	}

	// calculate average
	vector<double>	weight;
	Bimage* 		psum = img_add_weighed(file_list, weight, nuavg, nustd, flags);

	// Masking / amplitude weighting
	if ( setweight > 0 ) {
		Bimage*		pw = read_img(fomfile, 1, -1);

		// Weigh amplitudes in Fourier space with values read from a file
		if ( setweight == 1 ) {	
	
			psum->fft();

			if ( divide ) psum->divide(pw, wiener);
			if ( multiply ) psum->multiply(pw);

			psum->fft_back();
		}
		// Weigh values in real space with values read from a file
		if ( setweight == 2 ) {	

			if ( divide ) psum->divide(pw, wiener);
			if ( multiply ) psum->multiply(pw);
		}
		
		delete pw;
	}

	if ( sampling > 0 ) psum->sampling(sampling);
	
	if ( set_origin ) psum->origin(origin);
	if ( set_origin ) {
		if ( set_origin == 2 ) psum->origin(psum->size()/2);
		else psum->origin(origin);
	}

	if ( outfile.length() ) {
		psum->change_type(nudatatype);
		write_img(outfile, psum, 0);
	}	
	
	delete psum;
	string_kill(file_list);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
