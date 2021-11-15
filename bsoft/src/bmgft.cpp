/**
@file	bmgft.cpp
@brief	Disk-based 3D reconstruction for a tomography series
@author	Bernard Heymann
@date	Created: 20060110
@date	Modified: 20090401
**/

#include "mg_tomo_rec.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmgft [options] input.star [input2.star]",
"-----------------------------------------------",
"Fourier transforms micrographs from a tomographic tilt series.",
" ",
"Actions:",
"-removemarkers 14        Mask out markers with this radius (pixels) before transformation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-select 23,45,66         Micrograph selection (default all).",
"-size 1024,1024,200      Intended reconstruction size (default from micrograph).",
"-scale 0.5               Scale of intended reconstruction (default 1).",
"-pad 3                   Image padding factor (default 2).",
"-fill 127                Fill value for erasing/painting markers (default average).",
" ",
"Output:",
"-output file.star        Output STAR file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	double			marker_radius(0);				// Radius to mask out markers
	DataType 		nudatatype(Unknown_Type);		// Conversion to new type
	Bstring			mg_select;						// Micrograph selection, default all
	Vector3<long>	size;							// Default use micrograph dimension
	double			scale(1);						// Intended reconstruction scale
	int				pad_factor(2);					// Image padding factor
	int 			fill_type(FILL_AVERAGE);		// Type of fill value
	double			fill(0);						// Fill value for new areas
	Bstring			outfile;						// Output STAR file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "removemarkers" )
			if ( ( marker_radius = curropt->value.real() ) < 1 )
				cerr << "-removemarkers: A marker radius must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "select" )
			mg_select = curropt->value;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.0001 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: An integer factor must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);

	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( mg_select.length() ) project_mg_select(project, mg_select);
	
	int				n = mg_fft_write(project, size, scale, pad_factor, 
						nudatatype, marker_radius, fill_type, fill);
	
	if ( n && project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

