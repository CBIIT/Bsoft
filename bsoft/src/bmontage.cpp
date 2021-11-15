/**
@file	bmontage.cpp
@brief	A program to generate a montage from a 3D image or a set of 2D images
@author Bernard Heymann
@date	Created: 20001126
@date	Modified: 20170117
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
"Usage: bmontage [options] input.img output.img",
"----------------------------------------------",
"Generates a montage from a 3D image or a set of 2D images",
" ",
"Actions:",
"-slices 3,4              Montage a 3D image in columns and rows (0,0: automatically calculated).",
"-images 3,4              Montage 2D images in columns and rows (0,0: automatically calculated).",
"-flip y                  Flip panel order on the x and/or y axis.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-first 12                First slice in 3D or first 2D image (default: 0).",
"-skip 2                  Number of slices or images to skip (default: 0).",
"-pad 5                   Padding to add around slices or images (default 0).",
"-fill 127                Fill value for padding (default average).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setmontage(0);
	long 			first(0);					// First slice or image in the montage
	long			skip(0);					// Number slices or images to skip
	long 			columns(0), rows(0);		// Determine montage dimensions automatically
	long 			pad(0);						// Padding around images
	int				flip(0);					// Flag to flip panels on: 1=x axis, 2=y axis
	double			fill(0);	 				// Fill value for resizing
	int 			fill_type(FILL_BACKGROUND);	// Fill type for resizing
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "slices" ) {
			if ( curropt->values(columns, rows) < 2 )
				cerr << "-slices: Columns and rows must be specified!" << endl;
			else
				setmontage = 1;
		}
		if ( curropt->tag == "images" ) {
        	if ( curropt->values(columns, rows) < 2 )
				cerr << "-images: Columns and rows must be specified!" << endl;
			else
				setmontage = 2;
        }
		if ( curropt->tag == "flip" ) {
			if ( curropt->value.contains("x") ) flip |= 1;
			if ( curropt->value.contains("y") ) flip |= 2;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "first" )
			if ( ( first = curropt->value.integer() ) < 1 )
				cerr << "-first: A number must be specified!" << endl;
		if ( curropt->tag == "skip" )
			if ( ( skip = curropt->value.integer() ) < 1 )
				cerr << "-skip: A number must be specified!" << endl;
		if ( curropt->tag == "pad" )
			if ( ( pad = curropt->value.integer() ) < 1 )
				cerr << "-pad: A number of pixels must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
    }
	option_kill(option);
	
	double		ti = timer_start();

	// Read all the file names
	Bstring*		file_list = NULL;
	while ( optind < argc - 1 ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No image files specified!" << endl;
		bexit(-1);
	}

	Vector3<long>	nusize;
	Vector3<long>	translate = {pad, pad, 0};

	Bstring 		rawstring;
	Bimage*			p = img_catenate(file_list, rawstring, Unknown_Type,
							nusize, 0, fill_type, fill, 0, 0);
	
	// Read image file
//	int 		dataflag = 0;
//	if ( optind < argc - 1 ) dataflag = 1;
//	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL )  {
		cerr << "Error: No input file(s) read!" << endl;
		bexit(-1);
	}
	
	if ( !setmontage && first ) {
		if ( p->sizeZ() > 1 ) setmontage = 1;
		else setmontage = 2;
		columns = rows = 1;
	}
	
	if ( setmontage == 2 && skip == 0 ) skip = p->sizeZ() - 1;
	
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) {
		if ( fabs(p->background(long(0))) < 1e-6 )
			p->calculate_background();
		fill = p->background(long(0));
	}
	
	Bimage*		pm = NULL;
	
	p->background(fill);
	
	if ( pad ) {
		nusize = {p->sizeX()+2*pad, p->sizeY()+2*pad, p->sizeZ()};
		p->resize(nusize, translate, fill_type, fill);
	}

	if ( setmontage )
		pm = p->montage(first, columns, rows, skip, flip);
	
	if ( pm ) {
		delete p;
		p = pm;
	}
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

