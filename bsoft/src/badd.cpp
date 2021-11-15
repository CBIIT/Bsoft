/**
@file	badd.cpp
@brief	Program to catenate image files
@author Bernard Heymann
@date	Created: 20010505
@date	Modified: 20110804
**/

#include "rwimg.h"
#include "img_combine.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: badd [options] input.images",
"----------------------------------",
"Adds multiple images to generate a new composite image.",
"Any number of input images may be given and may include the wild card \"*\".",
"All images are converted to floating point.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale input images to average and standard deviation.",
"-average                 Output average image rather than the sum.",
"-std                     Output standard deviation image rather than variance image.",
"-center                  Center images before summation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-weights 0.1,0.3,0.9     Weights for input files (default 1 for every file).",
" ",
"Output:",
"-output output.img       Output image (default sum.map).",
"-fom fom.img             Output variance/standard deviation image.",
" ",
"Examples:",
"To sum three images with different weights:",
"badd -v 7 -rescale 0,1 -weights 0.5,1.2,0.87 -output new.pif in1.map in2.map in3.map",
" ",
"To output the standard deviation of the summed input images:",
"badd -v 7 -fom stdev.map -output sum.pif *.spi",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;			// New image origin
	int				set_origin(0);				// Flag to set origin
	Bstring			sumfile("sum.map");			// Default output image file name
	double			nuavg(0), nustd(0); 		// Average and standard deviation for rescaling
	Bstring			weight_string;
	int 			setweigh(0);				// Default weights are 1
	int				flags(0);					// Flags to modify summation
	Bstring			fomfile;					// Output FOM file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified" << endl;
		if ( curropt->tag == "average" ) flags |= 1;
		if ( curropt->tag == "std" ) flags |= 4;
		if ( curropt->tag == "center" ) flags |= 8;
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
		if ( curropt->tag == "weights" ) {
			weight_string = curropt->value;
			if ( weight_string.length() < 1 )
				cerr << "-weights: Weights for all files must be specified" << endl;
			else
				setweigh = 1;
		}
		if ( curropt->tag == "output" )
        	sumfile = curropt->filename();
		if ( curropt->tag == "fom" ) {
        	fomfile = curropt->filename();
			flags |= 2;
		}
    }
	option_kill(option);
	
	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();

	if ( optind >= argc ) {
		cerr << "No input files!" << endl;
		bexit(-1);
	}
	
	// Set up image file names and weights
	long			nfiles(0);
	Bstring*		file_list = NULL;
	while ( optind < argc ) {
		string_add(&file_list, argv[optind++]);
		nfiles++;
	}
	
	vector<double>	weight;
	if ( setweigh ) {
		weight = weight_string.split_into_doubles(",");
		if ( weight.size() < nfiles ) {
			cerr << "Error: The number of weights must equal the number of input files" << endl;
			bexit(-1);
		} 
	}

	Bimage*			psum = img_add_weighed(file_list, weight, nuavg, nustd, flags);
	
	if ( sam.volume() > 0 ) psum->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) psum->origin(psum->default_origin());
		else psum->origin(origin);
	}
	
	if ( sumfile.length() ) {
		psum->change_type(nudatatype);
		write_img(sumfile, psum, 0);
	}
	
	if ( fomfile.length() )
		write_img(fomfile, psum->next, 0);
	
	delete psum;
	string_kill(file_list);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
