/**
@file	bgraphseg.cpp
@brief	Segment images into density regions using a nondirected graph
@author Bernard Heymann
@date	Created: 20110318
@date	Modified: 20210303
**/

#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bgraphseg [options] input.img output.img",
"-----------------------------------------------",
"Segment images into density regions using a nondirected graph.",
" ",
"Actions:",
"-invert                  Invert density in the image before other operations.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation before other operations.",
"-gaussian 11,2.6         Gaussian smoothing filter: kernel size and sigma.",
"-type srm                Type of segmentation (simple, srm).",
"-minimum 2300            Join segments smaller than this size.",
"-voxel 23,241,193,1      Select the segment containing this voxel with an optional inversion flag.",
"-colorize                Generate random segment colors.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-connectivity            Flag to use full neighbor connectivity (default 2*dimension).",
"-complexity 1.6          Segmentation complexity: larger value gives more regions (default 1).",
" ",
"Input:",
"-Mask mask.tif           Mask to limit segmentation.",
" ",
"Output:",
"-mask mask.tif           Multi-level mask with segmented regions.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setinvert(0);				// Flag to invert density
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	long			gauss_kernel(0), sigma(0);	// Gaussian kernel size and sigma
	int				type(0);					// Flag and segmentation type
	int				connect_type(0);			// Flag for full connectivity
	long			min_size(0);				// Minimum size of a blob to keep
	Vector3<long>	voxel(-1,0,0);				// Voxel to select a segement
	int				colorize(0);				// Generate color mask
	int				mask_invert(0);				// Invert voxel-selected mask
	double			complexity(1);				// Complexity value
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Bstring			maskfile;					// Mask file
	Bstring			segmask;					// Multi-level mask output

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "gaussian" )
			if ( curropt->values(gauss_kernel, sigma) < 2 )
				cerr << "-gaussian: The kernel edge size and sigma must be specified!" << endl;
		if ( curropt->tag == "type" ) {
			if ( curropt->value.contains("sim") ) type = 1;
			if ( curropt->value.contains("srm") ) type = 2;
		}
		if ( curropt->tag == "connectivity" )
			connect_type = 1;
		if ( curropt->tag == "minimum" )
			if ( ( min_size = curropt->value.integer() ) < 1 )
				cerr << "-minimum: A minimu segment size must be specified!" << endl;
		if ( curropt->tag == "voxel" )
			if ( curropt->values(voxel[0], voxel[1], voxel[2], mask_invert) < 3 )
				cerr << "-voxel: A voxel must be specified!" << endl;
		if ( curropt->tag == "colorize" )
			colorize = 1;
		if ( curropt->tag == "complexity" )
			if ( ( complexity = curropt->value.real() ) < 0.001 )
				cerr << "-complexity: A complexity value must be specified!" << endl;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "mask" )
			segmask = curropt->filename();
    }
	option_kill(option);
        
	double			ti = timer_start();
	
    // Read the input file
    Bimage* 	p = read_img(argv[optind++], 1, 0);

	
	Bimage*		pmask = NULL;
	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		p->mask(pmask, p->average());
	}
	
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	
	if ( nudatatype > p->data_type() ) p->change_type(nudatatype);
		
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}

	if ( setinvert ) p->invert();

	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( sigma > 0 ) p->filter_gaussian(gauss_kernel, sigma);

	Bimage*			pseg = p;
	Bimage*			pmm = NULL;
	if ( type ) {
		GSgraph		g = p->graph_segment(type, connect_type, complexity, min_size);
		pseg = p->graph_segments_to_image(g);
		if ( segmask.length() ) pmm = p->graph_segments_to_mask(g);
	}
	
	if ( pmm && voxel[0] >= 0 ) {
		pmm->levelmask_select(0, voxel);
		if ( mask_invert ) pmm->mask_invert();
	}
	
	if ( pmm && colorize ) pmm->levelmask_colorize();
	
	// Write an output file if a file name is given
    if ( optind < argc )
		write_img(argv[optind], pseg, 0);
	
	// Write an output file if a file name is given
    if ( pmm && segmask.length() )
		write_img(segmask, pmm, 0);
	
	delete p;
	delete pmask;
	delete pmm;
	if ( p != pseg ) delete pseg;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

