/**
@file	bmask.cpp
@brief	Generating and manipulating masks.
@author Bernard Heymann
@date	Created: 20030831
@date	Modified: 20190821
**/

#include "rwimg.h"
#include "symmetry.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmask [options] input.img output.img",
"-------------------------------------------",
"Generates and uses binary masks of the size given or of the input image.",
"The selected regions have value 1 and the omitted regions 0.",
"The default is a blank mask created of the same size as the input image.",
" ",
"Actions for mask creation:",
"-create 20,60,30         Create a new mask of this size in place of input image",
"-threshold -0.55,1.2     One or more thresholds to generate a mask.",
"-variance 5,1            Generate mask based on local variance in a kernel of a given size,",
"                         and flag to generate a background with value -1.",
" ",
"Actions for mask modification:",
"-invert                  Invert mask.",
//"-operation xor           Binary operation for mask modifications: replace/and/or/xor (default replace).",
"-dilate 2                Dilate mask a number of times.",
"-erode 3                 Erode mask a number of times.",
"-open 4                  Open (erode-dilate) mask a number of times.",
"-close 2                 Close (dilate-erode) mask a number of times.",
"-hole 3,56,23            Fill a hole indicated by the voxel.",
"-plane 0.3,0.5,-0.1,0,3,25 Mask on one side of a plane with the given normal and origin.",
"-rectangle 40,20,45      Rectangular mask with given length, width and rotation angle.",
"-radius 10,45            Generate a mask inside these radii.",
"-asu D3                  Generate an asymmetric unit mask.",
"-color size              Generate a RGBA image (random/size).",
"-size                    Convert level mask to region size map.",
" ",
"Actions for multi-level masks:",
"-collapse                Convert a multi-level mask to a binary mask.",
"-select 1,5-7            Select levels in a multi-level mask.",
"-voxel 23,241,193        Select the level including this voxel.",
"-level 134               Convert a mask to multiple levels, with the given voxels per level.",
"-bands 250,-1,150,1,30,0 Banded mask for reciprocal space (pairs of numbers: angstrom,value).",
"                         Band with 0: Not used.",
"                         Band with 1: Use for orientation determination.",
"                         Band with -1: Use for cross-validation.",
" ",
"Actions for mask application:",
"-apply                   Apply the mask to the input image (default output mask).",
"-rescaleapply 0.5,2.3    Rescale region inside mask and apply the mask to the input image.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for the mask (default from input image).",
"-fill 28.3               Fill value (default average).",
"-wrap                    Flag to wrap around image boundaries (default not).",
" ",
"Input:",
"-Mask file.img           Input mask to apply or get difference with main input image.",
" ",
"Examples:",
"1.	Creating a new mask:",
"		bmask -v 7 -create 100,100,50 -radius 10,20 -origin 50,50,25 mask.mrc",
" ",
"2.	Generating a mask from an image:",
"		bmask -v 7 -threshold 1.5 test.map mask.mrc",
" ",
"3.	Applying a mask to an image:",
"		bmask -v 7 -Mask mask.mrc test.map new.map",
" ",
"4.	Generating and applying a mask to an image:",
"		bmask -v 7 -apply -radius 50,80 -origin 100,100,100 test.map new.map",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	int				invert(0);					// Inversion flag
//	int				operation(0);				// Default operation is and
	Vector3<long>	nusize;						// Size of new mask
	vector<double>	threshold;					// Threshold array
	int 			var_kernel(0);				// Variance filter kernel size
	int				bgk_flag(0);				// Low variance threshold
	double			rect_length(0);				// Rectangular mask length
	double			rect_width(0);				// Rectangular mask width
	double			rect_angle(0);				// Rectangular mask rotation angle
	double			rad_min(0);					// Minimum radius
	double			rad_max(0);					// Maximum radius, not used if 0
	Bsymmetry		asu_sym;					// Asymmetric unit mask symmmetry
	Vector3<double>	origin;						// Mask origin
	Vector3<double>	plane_origin;				// Point on plane
	Vector3<double>	plane_normal;				// Plane normal
	int				set_origin(0);				// Flag to set origin
	int				wrap(0);					// Flag to wrap
	int 			apply_flag(0);				// Flag to apply mask
	double			nuavg(0), nustd(0);			// Rescaling to average and stdev
	int				dilate(0);					// Number of times to dilate
	int				erode(0);					// Number of times to erode
	int				open(0);					// Number of times to open
	int				close(0);					// Number of times to close
	Vector3<long>	fill_voxel;					// Voxel to fill from
	Bstring			color;						// Generate RGBA image: random or size-based colors
	int				convert_size(0);			// Flag to cenvert level mask to region size
	int				collapse(0);				// Collapse a multi-level mask to a binary mask
	Bstring			level_select;				// Selection of levels
	Vector3<long>	voxel(-1,0,0);				// Voxel to select a level
	int				level(0);					// Number of voxels per level
	vector<double>	band;						// Reciprocal space bands
	int 			fill_type(FILL_AVERAGE);	// Use average
	double			fill(0);					// Default fill is average
	Bstring			mask_file;					// Mask file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" )
			nusize = curropt->size();
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
 		if ( curropt->tag == "threshold" ) {
			threshold = curropt->value.split_into_doubles(",");
       	    if ( threshold.size() < 1 )
				cerr << "-threshold: At least one threshold must be specified!" << endl;
		}
		if ( curropt->tag == "variance" )
			if ( curropt->values(var_kernel, bgk_flag) < 1 )
				cerr << "-variance: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "plane" )
			if ( curropt->box(plane_normal, plane_origin) < 3 )
				cerr << "-plane: The plane normal must be specified" << endl;
		if ( curropt->tag == "rectangle" ) {
			if ( curropt->values(rect_length, rect_width, rect_angle) < 2 )
				cerr << "-rectangle: Both length and width must be specified" << endl;
			rect_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "radius" )
			if ( curropt->values(rad_min, rad_max) < 2 )
				cerr << "-radius: Both radii must be specified" << endl;
		if ( curropt->tag == "asu" )
			asu_sym = curropt->symmetry();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "wrap" )
			wrap = 1;
		if ( curropt->tag == "apply" )
			apply_flag = 1;
		if ( curropt->tag == "rescaleapply" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescaleapply: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "invert" )
			invert = 1;
/*		if ( curropt->tag == "operation" ) {
			if ( curropt->value.contains("xor") ) operation = 3;
			else if ( curropt->value.contains("or") ) operation = 2;
			else if ( curropt->value.contains("and") ) operation = 1;
			else operation = 0;
		}*/
		if ( curropt->tag == "dilate" )
       	    if ( ( dilate = curropt->value.integer() ) < 1 )
				cerr << "-dilate: A number of times to dilate must be specified!" << endl;
		if ( curropt->tag == "erode" )
       	    if ( ( erode = curropt->value.integer() ) < 1 )
				cerr << "-erode: A number of times to erode must be specified!" << endl;
		if ( curropt->tag == "open" )
       	    if ( ( open = curropt->value.integer() ) < 1 )
				cerr << "-open: A number of times to open must be specified!" << endl;
		if ( curropt->tag == "close" )
       	    if ( ( close = curropt->value.integer() ) < 1 )
				cerr << "-close: A number of times to close must be specified!" << endl;
		if ( curropt->tag == "hole" )
			if ( curropt->values(fill_voxel[0], fill_voxel[1], fill_voxel[2]) < 3 )
				cerr << "-hole: A voxel must be specified!" << endl;
		if ( curropt->tag == "color" )
			color = curropt->value;
		if ( curropt->tag == "size" )
			convert_size = 1;
		if ( curropt->tag == "collapse" )
			collapse = 1;
		if ( curropt->tag == "select" ) {
			level_select = curropt->value;
			if ( level_select.length() < 1 )
				cerr << "-select: Level numbers must be specified!" << endl;
		}
		if ( curropt->tag == "voxel" )
			if ( curropt->values(voxel[0], voxel[1], voxel[2]) < 3 )
				cerr << "-voxel: A voxel must be specified!" << endl;
		if ( curropt->tag == "level" )
       	    if ( ( level = curropt->value.integer() ) < 1 )
				cerr << "-level: The voxels per level must be specified!" << endl;
        if ( curropt->tag == "bands" ) {
			band = curropt->value.split_into_doubles(",");
			if ( band.size() < 2 )
				cerr << "-bands: At least two pairs of numbers must be specified!" << endl;
			band[band.size()] = 0;
		}
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
    Bimage* 	p = NULL;
	Bimage* 	pmask = NULL;
	if ( nusize.volume() > 0 ) {	// Create the mask in place of an input image
		if ( verbose & VERB_PROCESS )
			cout << "Creating a mask of size:            " << nusize << endl;
		if ( nudatatype == Unknown_Type ) nudatatype = UCharacter;
		pmask = new Bimage(nudatatype, TSimple, nusize, 1);
		if ( set_origin ) {
			if ( set_origin == 2 ) origin = pmask->default_origin();
			pmask->origin(origin);
		}
//		operation = 0;					// Default operation is replace
		pmask->fill(fill);
	} else {
		p = read_img(argv[optind++], 1, -1);
		if ( p == NULL ) bexit(-1);
		if ( set_origin ) {
			if ( set_origin == 2 ) origin = p->default_origin();
			p->origin(origin);
		}
		if ( sam.volume() > 0 ) p->sampling(sam);
		if ( fill_type == FILL_AVERAGE ) fill = p->average();
		if ( fill_type == FILL_BACKGROUND ) fill = p->background(long(0));
	}
	
	if ( !pmask ) {		// Generate the mask from other sources
		if ( mask_file.length() ) {
			pmask = read_img(mask_file, 1, -1);
			if ( pmask == NULL ) bexit(-1);
			apply_flag = 1;
		} else if ( threshold.size() > 0 ) {
			pmask = p->mask_by_thresholds(threshold);
		} else if ( var_kernel > 2 ) {
			pmask = p->variance_mask(var_kernel, 1e-6, bgk_flag);
		} else {
			pmask = p->copy();
		}
	}
	
	if ( level_select.length() ) pmask->levelmask_select(level_select, 1);

	if ( voxel[0] >= 0 ) pmask->levelmask_select(0, voxel);
	
	if ( collapse ) pmask->levelmask_collapse();

	if ( sam.volume() > 0 ) pmask->sampling(sam);

	if ( plane_normal.length() > 0 )
		pmask->mask_plane(plane_origin, plane_normal);
	
	if ( rect_length*rect_width > 0 )
		pmask->mask_rectangle(rect_length, rect_width, rect_angle, wrap);
	
	if ( rad_max )
		pmask->mask_shell(origin, rad_min, rad_max);
	
	if ( band.size() && band[0] )
		pmask->mask_fspace_banded(band);
	
	Bimage*			psymmask = NULL;
	if ( asu_sym.point() > 101 ) {
		psymmask = pmask->levelmask_asymmetric_units(asu_sym, 1);
		pmask->multiply(psymmask);
		delete psymmask;
	}
	
	if ( dilate > 0 )
		pmask->mask_dilate(dilate);
	
	if ( erode > 0 )
		pmask->mask_erode(erode);
	
	if ( open > 0 )
		pmask->mask_open(open);
	
	if ( close > 0 )
		pmask->mask_close(close);
	
	if ( fill_voxel.length() )
		pmask->mask_fill(fill_voxel);
	
	if ( invert ) pmask->mask_invert();
	
	if ( p ) {
		if ( nustd > 0 ) {
			p->rescale_to_avg_std(nuavg, nustd, pmask);
			p->multiply(pmask);
		} else if ( apply_flag ) {
			p->mask(pmask, fill);	// For output of image
		} else {
			delete p;
			p = pmask;		// For output of mask
		}
	} else {
		p = pmask;		// For output of mask
	}
	
	if ( level ) pmask->mask_split(level);

	if ( color.length() ) {
		if ( color[0] == 'r' ) pmask->simple_to_rgba();
		else if ( color[0] == 's' ) p = pmask->levelmask_color_by_size();
	} else if ( convert_size ) pmask->levelmask_region_size();
	
	pmask->mask_stats();

	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	if ( p != pmask ) delete pmask;
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
