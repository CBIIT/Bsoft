/**
@file	bhead.cpp
@brief	Program to modify image parameters in the header.
@author Bernard Heymann
@date	Created: 19990801
@date	Modified: 20210228
**/

#include "rwimg.h"
#include "rwsymop.h"
#include "utilities.h"
#include "options.h"
#include "versions.h"
#include "timer.h"

#include <iostream>
using namespace std;

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bhead [options] input.img output.img",
"-------------------------------------------",
"Reports on and modifies image parameters in the image header.",
" ",
"Image format extensions supported:",
"BioRad (pic), Brix (brx), CCP4 (map), Digital Instruments (di)",
"EM (em), Goodford (pot), GRD (grd), Imagic (img), MFF (mff)", 
"Image Magick (miff), MRC (mrc), PIC (bp), PIF (pif), Spider (spi)", 
"Suprim (spm), TIFF (tif), RAW (raw)",
" ",
"Tags to use with the -get option:",
"	size: image size",
"	channels: number of channels",
"	x: x-dimension",
"	y: y-dimension",
"	z: z-dimension",
"	images: number of images",
"	sampling or pixelsize: size if pixels in Å",
"	origin: image origin",
"	view: image view",
"	min: minimum",
"	max: maximum",
"	avg: average",
"	std: standard deviation",
"	stat: statistics",
"	text: header label",
"	json: header metadata",
" ",
"Background calculation flag:",
"	1=outside enclosing circle",
"	2=inside enclosing circle",
"	3=inside half-radius circle",
" ",
"Actions:",
"-information             Show sub-image information.",
"-moments 6               Show sub-image moments up to the indicated order.",
"-get size                Show selected information.",
"-recalculate             Force recalculation of statistics.",
"-images                  Interpret slices of a single 3D image as 2D images.",
"-slices                  Interpret multiple 2D images as z-slices of a single 3D image.",
"-splitchannels           Convert a multi-channel into a single channel image.",
"-combinechannels 4,View  Convert a single channel into a multi-channel image.",
"-background 1            Calculate new background values (image not modified).",
"-Background 2            Calculate new background values and mask the background.",
"-subtractbackground      Subtract the background.",
"-shiftbackground 1.6     Shift the background to this value.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-select 1                Select an image (first image = 0; default: all images).",
"-origin 0,-10.5,30       Set the origin.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-Sampling 1.5,1.5,1.5    Sets sampling and unit cell dimensions to fill image (A/pixel).",
"-setbackground 0.5       Set the background to this value.",
"-symmetry D6             Point group identifier.",
"-spacegroup 90           Set space group.",
"-unitcell 82,82.3,100,90,66.7,70 Set unit cell parameters (A and degrees).",
"-View 0.3,0.5,-0.2,35    Sets the view parameters for every image.",
"-label \"New title\"       Write a new title or label into the header.",
"-fixtype                 8 or 16 bit integer types are switched between signed and unsigned.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			show_sub_images(0);			// Flag to show sub-image info
	int 			show_moments(0);			// Flag to show sub-image moments
	Bstring			get;						// Tag for selected information
	long 			setimg(-1);					// Select all images
	Vector3<double>	origin;						// Origin
	int				set_origin(0);				// Flag to set origin
	int 			set_stats(0);
	int 			set_sampling(0);			// 1=set sampling, 2=set sampling and reset unitcell
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Bstring			symmetry_string;			// No symmetry specified
	int 			spacegroup(0);				// Space group
	UnitCell		uc;							// Unit cell parameters
	int 			calc_background(0);			// Flag to calculate background values
	int 			correct_background(0);		// Flag to correct backgrounds
	int 			set_background(0);			// Flag to set backgrounds
	int 			subtract_background(0);		// Flag to subtract backgrounds
	int 			shift_background(0);		// Flag to shift backgrounds
	double	 		newbackground(0);			// New background
	View			view;						// View parameters
	int				set_view(0);				// Flag to set the view
	Bstring			label;						// Header title or label
	int				znswitch(0);				// 0=not, 1=z2n, 2=n2z
	int				cnswitch(0);				// 0=not, 1=c2n, 2=n2c
	CompoundType	type(TSimple);
	int				fixtype(0);					// Flag to switch between signed and unsigned
	int 			dataflag(0);				// Flag to set if data must be read
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "information" )
			show_sub_images = 1;
		if ( curropt->tag == "moments" )
			if ( ( show_moments = curropt->value.integer() ) < 0 )
				cerr << "-moments: A number of orders must be specified!" << endl;
		if ( curropt->tag == "get" )
			get = curropt->value.lower();
		if ( curropt->tag == "recalculate" )
			set_stats = dataflag = 1;
		if ( curropt->tag == "select" )
			if ( ( setimg = curropt->value.integer() ) < 0 )
				cerr << "-select: An image number must be specified!" << endl;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" ) {
			sam = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "Sampling" ) {
			sam = curropt->scale();
			if ( sam[0] > 0 ) set_sampling = 2;
		}
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "spacegroup" )
			if ( ( spacegroup = curropt->value.integer() ) < 1 )
				cerr << "-spacegroup: A valid space group must be specified!" << endl;
		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
		if ( curropt->tag == "background" ) {
			dataflag = 1;
			if ( ( calc_background = curropt->value.integer() ) < 1 )
				cerr << "-background: A calculation type must be specified!" << endl;
		}
		if ( curropt->tag == "Background" ) {
			dataflag = 1;
			if ( ( correct_background = curropt->value.integer() ) < 1 )
				cerr << "-Background: A calculation type must be specified!" << endl;
		}
		if ( curropt->tag == "setbackground" ) {
        	newbackground = curropt->value.real();
			set_background = 1;
        }
		if ( curropt->tag == "subtractbackground" )
        	subtract_background = dataflag = 1;
		if ( curropt->tag == "shiftbackground" ) {
        	newbackground = curropt->value.real();
			shift_background = dataflag = 1;
        }
		if ( curropt->tag == "View" ) {
			view = curropt->view();
			set_view = 1;
		}
		if ( curropt->tag == "label" )
			label = curropt->value;
		if ( curropt->tag == "images" )
			znswitch = 1;
		if ( curropt->tag == "slices" )
			znswitch = 2;
		if ( curropt->tag == "splitchannels" )
			cnswitch = 1;
		if ( curropt->tag == "combinechannels" ) {
			if ( ( cnswitch = curropt->value.integer() ) < 2 )
				cerr << "-tochannels: The number of channels must be specified!" << endl;
			else
				if ( curropt->value.contains(",") )
					type = getcompoundtype(cnswitch, curropt->value.post(',').str());
		}
		if ( curropt->tag == "fixtype" )
			fixtype = 1;
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	if ( verbose & VERB_PROCESS ) systype(1);
	 
	if ( verbose & VERB_PROCESS ) cout << "Number of processors:           " << system_processors() << endl;
	
	if ( get.contains("lib") || get.contains("ver") )
		show_library_versions();
	if ( optind >= argc ) bexit(0);
	 
	// Read image file
	if ( set_stats || show_moments || calc_background || optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, setimg);
	if ( !p ) bexit(-1);
	
	if ( fixtype ) p->fix_type();
	
	if ( znswitch == 1 )
		if ( p->slices_to_images() < 0 ) bexit(-1);
	
	if ( znswitch == 2 )
		if ( p->images_to_slices() < 0 ) bexit(-1);
	
	if ( cnswitch == 1 ) {
		Bimage*		pnu = p->split_channels();
		delete p;
		p = pnu;
	}
	
	if ( cnswitch > 1 )
		p->combine_channels(cnswitch, type);
	
	if ( set_stats ) {
		p->statistics();
		if ( verbose ) show_sub_images = 1;
	}
	
	if ( calc_background ) {
		p->calculate_background(calc_background);
		if (verbose ) show_sub_images = 1;
	}
	
	if ( optind < argc ) {
		if ( set_sampling ) {
			p->sampling(sam);
			if ( set_sampling > 1 ) p->unit_cell(UnitCell());
		}
	
		if ( set_origin ) {
			if ( set_origin == 2 ) p->origin(p->default_origin());
			else p->origin(origin);
		}
		
		if ( correct_background ) {
			p->correct_background(correct_background);
			if (verbose ) show_sub_images = 1;
		}
		
		if ( set_background ) {
			p->background(newbackground);
			if (verbose ) show_sub_images = 1;
		}
		
		if ( shift_background ) {
			p->shift_background(newbackground);
			if (verbose ) show_sub_images = 1;
		}
		
		if ( subtract_background ) {
			p->subtract_background();
			if (verbose ) show_sub_images = 1;
		}
		
		if ( set_view ) p->view(view);
	
/*		if ( set_minmax ) {
			p->min = p->smin = min;
			p->max = p->smax = max;
		}
*/		
		if ( symmetry_string.length() ) p->symmetry(symmetry_string.str());
			
		if ( uc.check() ) p->unit_cell(uc);
	
		if ( spacegroup ) {
			Bstring			symfile;
			int				nsym(0);
			char*			symop = read_symop(symfile, spacegroup, nsym);
			if ( verbose & VERB_FULL ) {
				cout << "Symmetry operator file:         " << symfile << endl;
				cout << "Number of symmetry operators:   " << nsym << endl;
				cout << "Operators:" << endl;
				for ( int i=0; i<80*nsym; i+=80 )
					cout << &symop[i] << endl;
				cout << endl;
			}
			p->space_group(spacegroup);
		}
		
		if ( nudatatype != Unknown_Type ) p->change_type(nudatatype);
		
		if ( label.length() ) p->label(label.str());
	
		write_img(argv[optind], p, 0);
	}

	if ( show_sub_images ) p->subimage_information();
	
	if ( show_moments ) p->moments(show_moments);
	
	if ( get.length() ) p->get(get);
		
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
