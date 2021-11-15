/**
@file	bmissing.cpp
@brief	Generating and manipulating masks.
@author Bernard Heymann
@date	Created: 20030831
@date	Modified: 20180914
**/

#include "rwimg.h"
#include "rwmg.h"
#include "mg_tomo_rec.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmissing [options] input.img/mg.star output.img",
"------------------------------------------------------",
"Generates binary masks corresponding to reciprocal space regions.",
"The selected regions have value 1 and the omitted regions 0.",
"The default is a blank mask created of the same size as the input image.",
" ",
"Actions for mask creation:",
"-create 20,60,30         Create a new mask of this size in place of input image",
"                         or using the micrograph orientations in the parameter file.",
"-find wedge              Find the missing regions and create a new mask.",
" ",
"Actions for mask modification:",
"-wedge 45,-70,65         Missing wedge mask: tilt axis angle, first and last tilt angles.",
"-pyramid 40,-60,60,130,-55,60 Missing pyramid mask: 2 sets of tilt axis angle, first and last tilt angles.",
"-cone 65                 Missing cone mask: one tilt angle.",
"-invert                  Invert mask just before writing the image.",
"-dilate 2                Dilate mask a number of times.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for the mask (default 0,0,0).",
//"-fill 28.3               Fill value (default average).",
"-wrap                    Flag to wrap around image boundaries (default not).",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
" ",
"Parameters for micrograph orientations:",
"-scale 1.2               Scale or magnification of mask compared to original images (default 1).",
" ",
"Examples:",
"1.	Creating a new mask:",
"		bmissing -v 7 -create 100,100,50 -wedge 45,-60,65 -origin 50,50,25 mask.mrc",
" ",
"2.	Generating a mask from an image:",
"		bmissing -v 7 -pyramid 30,-30,70,120,-60,50 -wrap input.map mask.mrc",
" ",
"3.	Generating a mask from micrograph orientations:",
"		bmissing -v 7 -create 100,100,50 -wrap mg.star mask.pif",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	int				invert(0);					// Inversion flag
//	int				operation(1);				// Default operation is and
	Vector3<long>	nusize;						// Size of new mask
	Bstring			find;						// Type of missing regions to find in input image
	double			tilt_axis(0);				// Tilt axis angle
	double			tilt_neg(0);				// First or negative tilt angle
	double			tilt_pos(0);				// Last or positive tilt angle
	double			tilt_axis2(0);				// Tilt axis angle
	double			tilt_neg2(0);				// First or negative tilt angle
	double			tilt_pos2(0);				// Last or positive tilt angle
	double			tilt_cone(0);				// Cone tilt angle
	int				dilate(0);					// Number of times to dilate
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	ori;						// Mask origin
//	int				wrap(0);					// Flag to wrap
	double			scale(1);					// Scale or magnification of mask compared to original images
	double 			resolution(0); 				// Must be set > 0 to limit resolution
//	int 			fill_type = FILL_USER;		// Use 1
//	double			fill = 1;					// Default fill is 1
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" )
			nusize = curropt->size();
		if ( curropt->tag == "find" )
			find = curropt->value;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "wedge" ) {
			if ( ( i = curropt->values(tilt_axis, tilt_neg, tilt_pos) ) < 2 )
				cerr << "-wedge: At least the tilt axis and one angle must be specified" << endl;
			else {
				tilt_axis *= M_PI/180.0;
				tilt_neg *= M_PI/180.0;
				tilt_pos *= M_PI/180.0;
				if ( i < 3 )  tilt_pos = tilt_neg;
			}
		}
		if ( curropt->tag == "pyramid" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 4 )
				cerr << "-pyramid: At least the four angles must be specified" << endl;
			else {
				tilt_axis = d[0] * M_PI/180.0;
				tilt_neg = d[1] * M_PI/180.0;
				tilt_pos = d[2] * M_PI/180.0;
				tilt_axis2 = d[3] * M_PI/180.0;
				if ( d.size() < 5 ) tilt_neg2 = tilt_neg;
				else tilt_neg2 = d[4] * M_PI/180.0;
				if ( d.size() < 6 ) tilt_pos2 = tilt_pos;
				else tilt_pos2 = d[5] * M_PI/180.0;
			}
		}
		if ( curropt->tag == "cone" ) {
			if ( ( tilt_cone = curropt->value.real() ) < 1 )
				cerr << "-cone: A tilt angle must be specified" << endl;
			else
				tilt_cone *= M_PI/180.0;
		}
		if ( curropt->tag == "dilate" )
       	    if ( ( dilate = curropt->value.integer() ) < 1 )
				cerr << "-dilate: A number of times to dilate must be specified!" << endl;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				ori = curropt->origin();
				set_origin = 1;
			}
		}
//		if ( curropt->tag == "fill" )
//			fill = curropt->fill(fill_type);
//		if ( curropt->tag == "wrap" ) wrap = 1;
		if ( curropt->tag == "invert" ) invert = 1;
/*		if ( curropt->tag == "operation" ) {
			if ( curropt->value.contains("xor") ) operation = 3;
			else if ( curropt->value.contains("or") ) operation = 2;
			else if ( curropt->value.contains("and") ) operation = 1;
			else operation = 0;
		}*/
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.0001 )
				cerr << "-scale: A scale must be specified!" << endl;
	}
	option_kill(option);
	
	double		ti = timer_start();
	
    Bimage* 	p = NULL;
	Bproject*	project = NULL;
	Bstring		filename(argv[optind++]);
	
	if ( file_type(filename) == Image ) {
		p = read_img(filename, 1, 0);
		if ( p == NULL ) bexit(-1);
		filename = argv[optind];
//		nusize = p->size();
	} else if ( file_type(filename) == Micrograph ) {
		project = read_project(filename);
		project->select = 1;
		if ( nusize.volume() > 0 )
			project_set_particle_box_size(project, nusize);
		p = project_missing_mask(project, nusize, ori, resolution, scale);
		filename = argv[optind];
	} else if ( nusize.volume() > 0 ) {
		if ( nudatatype == Unknown_Type ) nudatatype = UCharacter;
		p = new Bimage(nudatatype, TSimple, nusize, 1);
	}

	if ( p == NULL ) bexit(-1);

	if ( !project && ( verbose & VERB_PROCESS ) )
		cout << "Creating a mask of size:            " << nusize << endl;

/*
	if ( nusize.volume() > 0 ) {	// Create the mask in place of an input image
		if ( access(filename.c_str(), F_OK) ) {
			if ( verbose & VERB_PROCESS )
				cout << "Creating a mask of size:            " << nusize << endl;
			if ( nudatatype == Unknown_Type ) nudatatype = UCharacter;
			p = new Bimage(nudatatype, TSimple, nusize, 1);
//			operation = 0;					// Default operation is replace
//			fill = 1;
		} else if ( file_type(filename) == Micrograph ) {
			if ( nusize.volume() < 1 ) {
				cerr << "Error: A mask size must be specified!" << endl;
				bexit(-1);
			}
			project = read_project(filename);
			p = project_missing_mask(project, nusize, ori, resolution, scale);
			project_kill(project);
			filename = argv[optind];
		} else {
			cerr << "Error: " << filename << " already exists!" << endl;
		}
	} else {
		p = read_img(filename, 1, -1);
		if ( p == NULL ) bexit(-1);
//		if ( fill_type == FILL_AVERAGE ) fill = p->avg;
//		if ( fill_type == FILL_BACKGROUND ) fill = p->image[0].background();
		filename = argv[optind];
	}
*/

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(ori);
		ori = p->image->origin();
	}
	
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	if ( !project ) {
		if ( find.length() )
			p->mask_missing_find(ori, resolution, find);
		else if ( tilt_pos2 )
			p->mask_missing_pyramid(ori, tilt_axis, tilt_axis2,
					tilt_neg, tilt_pos, tilt_neg2, tilt_pos2, resolution);
		else if ( tilt_pos )
			p->mask_missing_wedge(ori, tilt_axis, tilt_neg, tilt_pos, resolution);
		else if ( tilt_cone )
			p->mask_missing_cone(ori, tilt_cone, resolution);
	}
	
	if ( invert ) p->mask_invert();
	
	if ( dilate > 0 )
		p->mask_dilate(dilate);
	
	p->mask_stats();

	if ( filename.length() ) {
		p->change_type(nudatatype);
		write_img(filename, p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

