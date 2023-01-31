/**
@file	bdens.cpp
@brief	Program to calculate density in selected regions of a map
@author Bernard Heymann
@date	Created: 20140424
@date	Modified: 20230110
**/

#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistence */
const char* use[] = {
" ",
"Usage: bdens [options] input.img output.img",
"-------------------------------------------",
"Calculate density in selected regions of a map.",
" ",
"Actions:",
"-box 12,34,2,45,34,25    Calculate statistics within this box (start and size).",
"-oval 5,10,28,50,22,89   Calculate statistics within this oval (start and size).",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
"-select 1                Select an image (first image = 0; default: 0).",
" ",
"Input:",
"-Mask mask.map           Calculate statistics with a mask with any number of levels.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	int				stat_type(0);					// Type: 1=box, 2=oval
	Vector3<long>	start, size;					// Start and size of box or oval to use
	DataType 		nudatatype(Unknown_Type);	 	// Conversion to new type
	int 			setimg(0);						// Select first image
	Bstring			maskfile;						// Mask file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "box" ) {
			if ( curropt->box(start, size) < 6 )
				cerr << "-box: A starting position and size must be specified!" << endl;
			else
				stat_type = 1;
		}
		if ( curropt->tag == "oval" ) {
			if ( curropt->box(start, size) < 6 )
				cerr << "-oval: A starting position and size must be specified!" << endl;
			else
				stat_type = 2;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "select" )
			if ( ( setimg = curropt->value.integer() ) < 0 )
				cerr << "-select: An image number must be specified!" << endl;
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
    }
	option_kill(option);
        
	double		ti = timer_start();
	
    // Read the input file
    Bimage* 	p = read_img(argv[optind++], 1, -1);

	Bimage*		pmask = NULL;
	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, -1);
		p->stats_in_mask(setimg, pmask);
		delete pmask;
	} else if ( stat_type > 0 ) {
		vector<double>	stats = p->stats_in_shape(setimg, stat_type, start, start+size);
		cout << "Number of voxels:        " << stats[0] << endl;
		cout << "Minimum:                 " << stats[1] << endl;
		cout << "Maximum:                 " << stats[2] << endl;
		cout << "Average density:         " << stats[3] << endl;
		cout << "Standard deviation:      " << stats[4] << endl << endl;
	}
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return(0);
}

