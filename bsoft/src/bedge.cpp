/**
@file	bedge.cpp
@brief	A program to float an image in a uniform background
@author Bernard Heymann
@date	Created: 19980520
@date	Modified: 20150717
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
"Usage: bedge [options] infile.mrc outfile.mrc",
"---------------------------------------------",
"Smooths the edge of an image to a uniform background.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-truncate -0.5,1.2       Truncate data to minimum and maximum before edge operation.",
"-invert                  Invert density in the image before edge operation.",
"-size 4000,4000,10       Size of region inside edge (default 90% of image size).",
//"-origin 2500,1500,5      Edge origin (deprecated).",
"-start 2500,1500,5       Start of edge (default 5% from image origin).",
//"-gausswidth 80           Gaussian width (deprecated, use -width).",
"-width 8                 Gaussian width (default 0.1, sharp edge, <0 inverts selection).",
//"-edgeshape oval          Edge shape (deprecated, use -shape).",
"-shape oval              Edge shape: rectangular (default), oval, cylindrical.",
"-fill 150                Fill value for background: average (default), background, or value.",
" ",
NULL
};

int 		main(int argc, char* argv[])
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			cutmin(0), cutmax(0);		// Truncation
	bool 			setinvert(0);				// Flag to invert density
    int     		setedge(0);    	    		// Rectangular edge smoothing
	Bstring			shape("rectangular");
	Vector3<long>	size;				// Size of smoothing area
	Vector3<double>	start;			// Start of smoothing edge
	int				set_start(0);				// Flag to set start
    double	   		width(0.1);   	    		// Gaussian width/sigma
    double	   		fill(127);   	    		// Fill value default for byte
	int 			fill_type(FILL_AVERAGE);	// Fill type for edge smoothing

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "truncate" )
			if ( curropt->values(cutmin, cutmax) < 2 )
				cerr << "-truncate: Both min and max must be specified!" << endl;
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "size" )
 			if ( curropt->values(size[0], size[1], size[2]) < 3 )
				cerr << "-size: All three dimensions must be specified!" << endl;
		if ( curropt->tag == "start" || curropt->tag == "origin" ) {
 			if ( curropt->values(start[0], start[1], start[2]) < 3 )
				cerr << "-start: All three dimensions must be specified!" << endl;
			set_start = 1;
		}
		if ( curropt->tag == "width" || curropt->tag == "gausswidth" ) {
			width = curropt->value.real();
	        if ( fabs(width) < 0.001 )
				cerr << "-width: The Gaussian width must be specified." << endl;
		}
		if ( curropt->tag == "shape" || curropt->tag == "edgeshape" ) {
			shape = curropt->value;
			if ( shape.length() < 1 ) {
				cerr << "-shape: An edge shape must be specified." << endl;
			} else {
				if ( tolower(shape[0]) == 'o' ) setedge = 1;
				if ( tolower(shape[0]) == 'c' ) setedge = 2;
			}
		}
 		if ( curropt->tag == "fill" )
	        fill = curropt->fill(fill_type);
   }
	option_kill(option);
    
    /* Start timer */
	double			ti = timer_start();
	
    // Read the input file
    Bimage* 		p = read_img(argv[optind++], 1, -1);
	Vector3<double>	sz = {double(p->sizeX()), double(p->sizeY()), double(p->sizeZ())};
	
	if ( cutmin < cutmax )
		p->truncate_to_min_max(cutmin, cutmax);

	if ( setinvert ) p->invert();

	if ( !set_start )
		start = sz * 0.05;
	
	if ( size.volume() < 1 ) {
		sz = sz * 0.95 - start;
		size = {long(sz[0]), long(sz[1]), long(sz[2])};
	}
	size = size.max(1);
    
	p->edge(setedge, size, start, width, fill_type, fill);
	
    // Write an output file if a file name is given
    if ( argc > optind ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

