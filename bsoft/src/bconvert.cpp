/**
@file	bconvert.cpp
@brief	Program to catenate image files
@author Bernard Heymann
@date	Created: 20070927
@date	Modified: 20170823
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
"Usage: bconvert [options] input.images",
"--------------------------------------",
"Converts a list of images into a new format, using the original names and replacing the extension.",
"Any number of input images may be given and may include the wild card \"*\".",
"(However, VMS does not support Unix style usage of wild cards).",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale input images to average and standard deviation.",
"-extension zyx           Change the file extension after writing.",
" ",
"Parameters:",
"-verbose 3               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-raw d=f#x=120,120,1     Formatting to reinterpret image files.",
" ",
"Output:",
"-output mrc              Output format (default pif).",
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
	Bstring			nuformat("pif");			// New image file format
	Bstring			nuext;						// New extension after writing
	Bstring			rawstring;
	int 			setrescale(0);
	double			nuavg(0), nustd(0);		// Rescale image
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both averge and standard deviation must be specified!" << endl;
			else
				setrescale = 1;
			if ( nustd <= 0 ) setrescale = 0;
		}
		if ( curropt->tag == "extension" ) nuext = curropt->value;
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
		if ( curropt->tag == "raw" ) {
			rawstring = curropt->value;
			if ( rawstring.length() < 1 )
				cerr << "-raw: A valid format reinterpretation string must be specified!" << endl;
		}
		if ( curropt->tag == "output" ) {
			nuformat = curropt->value;
			if ( nuformat.length() < 1 )
				cerr << "-output: An output format must be specified!" << endl;
		}
    }
	option_kill(option);

	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();
	
	// Read the images and convert them
	Bimage*			p = NULL;
	Bstring			filename, nuname;
	
	if ( nuformat[0] != '.' ) nuformat = "." + nuformat;
	if ( nuext.length() && nuext[0] != '.' ) nuext = "." + nuext;
	if ( rawstring.length() > 0 && rawstring[0] != '#' ) rawstring = "#" + rawstring;
	
	while ( optind < argc ) {
		filename = argv[optind++];
		filename += rawstring;
		p = read_img(filename, 1, -1);
		if ( p != NULL ) {
			if ( nudatatype == Unknown_Type )
				nudatatype = p->data_type();
			else if ( nudatatype > p->data_type() )
				p->change_type(nudatatype);
			if ( setrescale ) p->rescale_to_avg_std(nuavg, nustd);
			p->statistics();
			p->calculate_background();
			if ( sam.volume() > 0 ) p->sampling(sam);
			if ( set_origin ) {
				if ( set_origin == 2 ) p->origin(p->default_origin());
				else p->origin(origin);
			}
				filename = p->file_name();
			filename = filename.pre_rev('.') + nuformat;
			p->change_type(nudatatype);
			write_img(filename, p, 0);
			delete p;
			if ( nuext.length() ) {
				nuname = filename.pre_rev('.') + nuext;
				if ( rename(filename.c_str(), nuname.c_str()) )
					cerr << "File " << filename << " not renamed to " << nuname << endl;
			}
		}
	}
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

