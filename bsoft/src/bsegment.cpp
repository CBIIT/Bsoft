/**
@file	bsegment.cpp
@brief	Segment images into density regions
@author Bernard Heymann
@date	Created: 19981222
@date	Modified: 20160315
**/

#include "rwimg.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bsegment [options] input.img/input.star [output.img]",
"-----------------------------------------------------------",
"Analyzes images for segmentation.",
" ",
"Actions:",
"-invert                  Invert density in the image before other operations.",
"-rescale -0.1,5.2        Rescale extracted regions to average and standard deviation.",
"-blobs 0.45,14,201,1.4   Show blobs above threshold (0.45) and set outside minimum (14) ",
"                         and maximum (201) sizes to value (1.4).",
" ",
"Mask generating actions:",
"-mass 342k               Find threshold associated with this molecular weight (Da)",
"                         value can be given in Da, kDa (with k), MDa (with m) or Gda (with g).",
"                         Output binary mask (data type byte, values 0 and 1).",
"-threshold 1.54          Report mass calculated from the volume included at this threshold.",
"                         Output binary mask (data type byte, values 0 and 1).",
"-multiple 1.2            Generate multiple regions above this threshold.",
"-kmeans 12,25,0.8        K-means segmentation: K, iterations and density:distance ratio.",
"-innervolume 2.3,15      Calculate internal volume below this threshold, output mask frequency.",
" ",
"Selections:",
"-select 1,5-7            Select levels in a multi-level mask.",
"-combine 1,5-7           Combine levels to create a binary mask.",
"-voxel 15,50,33          Select a region around this voxel.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-sign -1                 Sign for thresholding (default 0=set by sign of threshold-average;",
"                         1=positive density; -1=negative density).",
"-rho 0.8                 Protein density (default 0.81 Da/A3, use with -mass).",
"-fill 0.02               Fill value for extraction: average (default), background, or value.",
" ",
"Input:",
"-Mask mask.tif           Input mask (must be integer and same size as input image).",
" ",
"Output:",
"-Region region.tif       File name for the excised region.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setinvert(0);				// Flag to invert density
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	double			threshold(-1e36);			// Threshold for defining regions
	long			nregions(0);				// Number of regions for K-means
	long			max_iter(10);				// Maximum kmeans iterations
	double			ratio(1);					// Density-distance ratio
	int				region_mask(0);				// Flag to generate a mask of regions
	Vector3<double> origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	int				sign(0);					// Sign for region mask threshold
	double			blob_threshold(-1e36);		// Threshold for defining blobs
	double			min_size(0);				// Minimum size of a blob to keep
	double			max_size(1000000000);		// Maximum size of a blob to keep
	double			setvalue(0);				// Value to set a deleted region to
	double			mol_weight(0);				// Molecular weight for threshold finding
	int				internal_volume(0);			// Flag to calculate internal volume
	int				mask_out_freq(0);			// Frequency to output inner volume mask
	double			rho(RHO);					// Protein density in Da/A3
	int 			nvox(0);					// Number of voxel coordinates
	Bstring			level_select;				// Selection of levels
	int				select_flag(1);				// 0=combine, 1=multi-level
	Vector3<long>	voxel;						// Voxel coordinates for region extraction
	int 			fill_type(FILL_AVERAGE);	// Fill type for smaller images
	double			fill(0);
	Bstring			mask_file = NULL;			// Mask file
	Bstring			region_file = NULL;			// Excised region file

	int				i, optind;
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
		if ( curropt->tag == "blobs" ) {
			if ( ( i = curropt->values(blob_threshold, min_size, max_size, setvalue )  ) < 1 )
				cerr << "-blobs: A threshold must be specified!" << endl;
			else
				if ( i < 4 ) setvalue = blob_threshold;
		}
		if ( curropt->tag == "mass" )
			mol_weight = curropt->real_units();
		if ( curropt->tag == "threshold" )
			if ( ( threshold = curropt->value.real() ) < -1e30 )
				cerr << "-threshold: A threshold must be specified!" << endl;
		if ( curropt->tag == "innervolume" ) {
			if ( curropt->values(threshold, mask_out_freq ) < 1 )
				cerr << "-innervolume: A threshold must be specified!" << endl;
			else
				internal_volume = 1;
		}
		if ( curropt->tag == "multiple" ) {
			if ( curropt->values(threshold, sign) < 1 )
				cerr << "-multiple: A threshold must be specified!" << endl;
			else
				region_mask = 1;
		}
		if ( curropt->tag == "kmeans" )
			if ( curropt->values(nregions, max_iter, ratio) < 1 )
				cerr << "-kmeans: A number of regions must be specified!" << endl;
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
			nvox = 1;
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sign" )
			if ( ( sign = curropt->value.integer() ) == 0 )
				cerr << "-sign: A sign value must be specified!" << endl;
		if ( curropt->tag == "rho" )
			if ( ( rho = curropt->value.real() ) < 0.001 )
				cerr << "-rho: A protein density value must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "Region" )
			region_file = curropt->filename();
    }
	option_kill(option);
        
	double			ti = timer_start();
	
    // Read the input file
	Bstring			filename(argv[optind++]);
    Bimage* 		p = read_img(filename, 1, -1);

	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
		
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}

	if ( setinvert ) p->invert();
	
	if ( blob_threshold > -1e29 )
		p->blobs(blob_threshold, min_size, max_size, setvalue, sign);

	if ( mol_weight > 0 )
		threshold = p->mass_threshold(0, mol_weight, rho);
	else if ( threshold )
		p->mass_at_threshold(0, threshold, rho);

	// Mask generation or reading
	Bimage*			pmask = NULL;
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, -1);
	} else if ( nregions > 1 ) {
		pmask = p->kmeans_segment(nregions, max_iter, ratio);
	} else if ( region_mask ) {
		pmask = p->regions(threshold, sign);
	} else if ( mol_weight > 0 ) {
		pmask = p->mask_by_threshold(threshold);
	} else if ( internal_volume ) {
		pmask = p->internal_volume(threshold, mask_out_freq);
	}
	
	// Mask level selection
	if ( pmask ) {
		if ( level_select.length() ) pmask->levelmask_select(level_select, select_flag);
		if ( nvox ) pmask->levelmask_select(0, voxel);
	}
	
	// Mask analysis
//	if ( pmask )
//		p->level_masked_stats(pmask);
	
	// Extract mask regions
	
	// Extract masked regions
	if ( pmask && region_file.length() ) {
		Bimage*		pex = p->level_mask_extract(pmask, fill_type, fill);
		pex->change_type(nudatatype);
		if ( nustd > 0 ) pex->rescale_to_avg_std(nuavg, nustd);
		write_img(region_file, pex, 0);
		delete pex;
	}

    // Write an output mask file if a file name is given
    if ( !pmask ) pmask = p;
    if ( optind < argc ) {
		filename = argv[optind];
		write_img(filename, pmask, 0);
	}
	
	delete p;
	if ( pmask && pmask != p ) delete pmask;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

