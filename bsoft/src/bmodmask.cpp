/**
@file	bmodmask.cpp
@brief	Generating a mask from an atomic model
@author Bernard Heymann
@date	Created: 20060301
@date	Modified: 20200329
**/

#include "rwmodel.h"
#include "model_mask.h"
#include "model_links.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodmask [options] in1.star [in2.star...]",
"------------------------------------------------",
"Generates a mask from a model if an output image is specified.",
" ",
"Actions:",
"-type hull               Mask type: model, hull, shell, level.",
"-invert                  Invert the mask after creation.",
"-project 2.5,D3          Projected masks at the angular step size and symmetry.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
" ",
"Parameters for mask creation:",
"-datatype u              Force writing of a new data type (default float).",
"-size 100,120,80         Mask size (default from model).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin for the mask.",
"-edge 1.5                Edge width for model mask (default 1 angstrom).",
"-width 34.8              Shell mask width (default 1 angstrom).",
"-curved                  Flag to indicate curved hull or shell surface (default not).",
"-fast                    Use fast but less comprehensive algorithm for shell and hull masks.",
"-smooth 2                Edge smoothing in pixels (default 0).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output file.star        Output model file.",
"-mask mask.map           Output mask image.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	DataType 		nudatatype = Float;	// Mask data type
	int				mask_type(-1);		// 0=model, 1=hull, 2=shell, 3=level, 4=project
	double			dang(0);			// Projection angular step size
	Bsymmetry		sym;			// Point group for projections
	Vector3<long>	mask_size;			// Size for creating a mask
	Vector3<double>	sam;	   			// Map sampling
	Vector3<double>	origin;				// Mask origin
	double			edge(1);			// Edge width for model mask
	double			shell_width(1);		// Shell mask width
	int				curv_flag(0);		// Curved surface flag for hull and shell masks
	int				invert(0);			// Flag for mask inversion
	double			comprad(0);			// Component display radius
	double			linkrad(0);			// Link display radius
	int				smooth(0);			// Edge smoothing in pixels
	int				fast(0);			// Flag for fast algorithm for hull and shell masks
	Bstring			paramfile;			// Input parameter file
	Bstring			outfile;			// Output model file
	Bstring			maskfile;			// Output mask file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'm' ) mask_type = 0;
			if ( curropt->value[0] == 'h' ) mask_type = 1;
			if ( curropt->value[0] == 's' ) mask_type = 2;
			if ( curropt->value[0] == 'l' ) mask_type = 3;
		}
		if ( curropt->tag == "project" ) {
			if ( ( dang = curropt->value.real() ) < 0.5 )
				cerr << "-project: An angular step size must be specified!" << endl;
			else {
				dang *= M_PI/180.0;
				if ( curropt->value.contains(",") )
					sym = Bsymmetry(curropt->value.post(','));
				mask_type = 4;
			}
		}
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 1 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "invert" ) invert = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
 		if ( curropt->tag == "size" )
			mask_size = curropt->size();
		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "edge" )
			if ( ( edge = curropt->value.real() ) < 0.001 )
				cerr << "-edge: The edge width must be specified!" << endl;
		if ( curropt->tag == "width" )
			if ( ( shell_width = curropt->value.real() ) < 1 )
				cerr << "-width: The shell width must be specified!" << endl;
		if ( curropt->tag == "curved" ) curv_flag = 1;
		if ( curropt->tag == "smooth" )
			if ( ( smooth = curropt->value.integer() ) < 1 )
				cerr << "-smooth: The number of smoothing pixels must be specified!" << endl;
		if ( curropt->tag == "fast" ) fast = 1;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bmodel*		model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( comprad ) models_process(model, comprad, model_set_component_radius);

	if ( linkrad ) models_process(model, linkrad, model_set_link_radius);

	Bimage*			pmask = NULL;
	if ( maskfile.length() && mask_size.volume() && mask_type >= 0 ) {
		switch ( mask_type ) {
			case 0:
				pmask = model_create_mask(model, mask_size, origin, sam, edge);
				break;
			case 1:
				pmask = model_create_hull_mask(model, mask_size, origin, sam, curv_flag, fast);
				break;
			case 2:
				pmask = model_create_shell_mask(model, mask_size, origin, sam, shell_width, curv_flag, fast);
				break;
			case 3:
				pmask = model_create_level_mask(model, mask_size, origin, sam);
				break;
			case 4:
				pmask = model_create_projected_mask(model, mask_size, origin, sam, dang, sym);
				break;
		}
		if ( invert ) pmask->mask_invert();
		if ( smooth > 0 ) pmask->filter_average(2*smooth + 1);
		pmask->change_type(nudatatype);
		model->mapfile(maskfile.str());
    	write_img(maskfile, pmask, 0);
		delete pmask;
	}
	
	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

    model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

