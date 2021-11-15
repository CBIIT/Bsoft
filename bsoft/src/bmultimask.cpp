/**
@file	bmultimask.cpp
@brief	Generating and manipulating composite masks.
@author Bernard Heymann
@date	Created: 20091117
@date	Modified: 20200609
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "model_mask.h"
#include "model_links.h"
#include "model_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bmultimask [options] input.img output.img",
"------------------------------------------------",
"Manipulates composite/level masks.",
"The input image is considered to have multiple integer levels (even for floating point).",
" ",
"Mask creation action:",
"-threshold -0.55,1.2     One or more thresholds for a conditional mask.",
" ",
"Actions:",
"-invert                  Invert mask before other operations.",
"-dilate 2                Dilate mask a number of times.",
"-switch 3,6              Switch two levels.",
"-color size              Generate a RGBA image (random/level/size).",
"-size                    Convert level mask to region size map.",
"-level 134               Convert a mask to multiple levels, with the given voxels per level.",
"-bands 250,-1,150,1,30,0 Banded mask for reciprocal space (pairs of numbers: angstrom,value).",
"                         Band with 0: Not used.",
"                         Band with 1: Use for orientation determination.",
"                         Band with -1: Use for cross-validation.",
"-collapse                Convert a multi-level mask to a binary mask.",
"-interfaces 3            Calculates interfaces for the indicated region (-1 gives whole matrix).",
"-merge 25,4              Calculates a region interface matrix and deletes/merges",
"                         regions smaller than the given size and interface size.",
"-symmetry C5             Set all symmetry-related segments to the same index.",
" ",
"Selections:",
"-select 1,5-7            Select levels in a multi-level mask.",
"-combine 1,5-7           Combine levels to create a binary mask.",
"-voxel 23,241,193        Select the level including this voxel.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for the mask (default from input image).",
" ",
"Parameters for model generation:",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
" ",
"Input:",
"-Mask file.img           Template mask to select levels.",
" ",
"Output:",
"-Model model.star        Output model file.",
"-Postscript his.ps       Output plot file of level size histogram.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				invert(0);				// Inversion flag
	int				dilate(0);				// Flag to dilate the mask
	long			index1(-1000), index2(-1000);	// Indices to switch
	vector<double>	threshold;				// Threshold array
	Bstring			color;					// RGBA: random/level/size-based colors
	int				convert_size(0);		// Flag to cenvert level mask to region size
	Bstring			level_select;			// Selection of levels
	int				select_flag(1);			// 0=combine, 1=multi-level
	Vector3<long>	voxel(-1,0,0);			// Voxel to select a level
	int				level(0);				// Number of voxels per level
	vector<double>	band;					// Reciprocal space bands
	int				collapse(0);			// Collapse a multi-level mask to a binary mask
	int				reg_if(-1000000);		// Calculate region interfaces
	int				min_size(0);			// Minimum size and flag to calculate interface matrix
	int				min_if(1);				// Minimum interface size to accept
	Bsymmetry		sym;					// Point group
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> sampling;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;					// Mask origin
	int				set_origin(0);			// Flag to set origin
	double			comprad(0);				// Component display radius
	double			linkrad(0);				// Link display radius
	Bstring			mask_file;				// Secondary mask file
	Bstring			mod_file;				// Output model file
	Bstring			ps_file;				// Poscript level size histogram plot

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( curropt->tag == "threshold" ) {
			threshold = curropt->value.split_into_doubles(",");
       	    if ( threshold.size() < 1 )
				cerr << "-threshold: At least one threshold must be specified!" << endl;
		}
		if ( curropt->tag == "invert" )
			invert = 1;
		if ( curropt->tag == "dilate" )
       	    if ( ( dilate = curropt->value.integer() ) < 1 )
				cerr << "-dilate: A number of times to dilate must be specified!" << endl;
		if ( curropt->tag == "switch" )
			if ( curropt->values(index1, index2) < 2 )
				cerr << "-switch: Both indices must be specified!" << endl;
		if ( curropt->tag == "color" )
			color = curropt->value;
		if ( curropt->tag == "size" )
			convert_size = 1;
		if ( curropt->tag == "select" ) {
			level_select = curropt->value;
			if ( level_select.length() < 1 )
				cerr << "-select: Level numbers must be specified!" << endl;
			else
				select_flag = 1;
		}
		if ( curropt->tag == "combine" ) {
			level_select = curropt->value;
			if ( level_select.length() < 1 )
				cerr << "-combine: Level numbers must be specified!" << endl;
			else
				select_flag = 0;
		}
		if ( curropt->tag == "voxel" ) {
			voxel = curropt->size();
			if ( voxel.volume() < 0 )
				cerr << "-voxel: A voxel must be specified!" << endl;
		}
		if ( curropt->tag == "level" )
       	    if ( ( level = curropt->value.integer() ) < 1 )
				cerr << "-level: The voxels per level must be specified!" << endl;
        if ( curropt->tag == "bands" ) {
			band = curropt->value.split_into_doubles(",");
			if ( band.size() < 2 )
				cerr << "-bands: At least two pairs of numbers must be specified!" << endl;
			band[band.size()] = 0;
		}
		if ( curropt->tag == "collapse" )
			collapse = 1;
		if ( curropt->tag == "interfaces" )
        	if ( ( reg_if = curropt->value.integer() ) < -1000000 )
				cerr << "-interfaces: A region number must be specified!" << endl;
		if ( curropt->tag == "merge" )
        	if ( curropt->values(min_size, min_if) < 1 )
				cerr << "-merge: A minimum size must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "Model" )
			mod_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage* 	p = read_img(argv[optind++], 1, -1);

	if ( p == NULL ) bexit(-1);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	if ( sampling.volume() > 0 ) p->sampling(sampling);

	Bimage*		pmask = NULL;
	if ( threshold.size() > 0 ) {
		pmask = p->mask_by_conditional_thresholds(threshold);
		delete p;
		p = pmask;
	}
	
	if ( invert ) p->invert();

	if ( level_select.length() ) {
		if ( select_flag )
			p->levelmask_select(level_select, select_flag);
		else
			p->levelmask_combine(level_select);
	}
	
	if ( voxel[0] >= 0 ) p->levelmask_select(0, voxel);
	
	if ( index1 > 0 ) p->levelmask_switch(index1, index2);
	
	if ( band.size() && band[0] )
		p->mask_fspace_banded(band);
	
	if ( level ) p->mask_split(level);

	if ( reg_if > -1000000 ) p->mask_region_interfaces(reg_if);
	
	if ( min_size > 0 ) p->mask_merge_delete(min_size, min_if);

	if ( dilate ) p->levelmask_dilate(dilate);
	
	if ( sym.point() > 101 ) p->levelmask_symmetrize(sym);
	
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, -1);
		if ( sampling.volume() > 0 ) pmask->sampling(sampling);
		p->levelmask_select(pmask);
		delete pmask;
	}
	
	if ( collapse ) p->levelmask_collapse();

	p->levelmask_clean();
	
	if ( color.length() ) {
		if ( color[0] == 'r' ) p->levelmask_colorize();
		else if ( color[0] == 'l' ) p->simple_to_rgba();
		else if ( color[0] == 's' ) {
			pmask = p->levelmask_color_by_size();
			delete p;
			p = pmask;
		}
	} else if ( convert_size ) {
		p->levelmask_region_size();
	}

	Bplot*		plot = NULL;
	if ( ps_file.length() ) {
		plot = p->levelmask_size_histogram();
		ps_plot(ps_file, plot);
		delete plot;
	}
	
	p->mask_stats();

	if ( mod_file.length() ) {
		Bmodel*		model = model_from_multilevel_mask(p);
		if ( comprad ) models_process(model, comprad, model_set_component_radius);
		if ( linkrad ) models_process(model, linkrad, model_set_link_radius);
		write_model(mod_file, model);
		model_kill(model);
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

