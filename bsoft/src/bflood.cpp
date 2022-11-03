/**
@file	bflood.cpp
@brief	Segment images into density regions
@author Bernard Heymann
@date	Created: 19981222
@date	Modified: 20170718
**/

#include "rwimg.h"
#include "rwmg.h"
#include "rwmodel.h"
#include "model_map.h"
#include "model_mask.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bflood [options] input.img/input.star [output.img]",
"---------------------------------------------------------",
"Segments images using the watershed or flooding algorithm.",
"Note: -series, -peaks and -threshold are mutually exclusive.",
" ",
"Actions:",
"-invert                  Invert density in the image before other operations.",
"-list                    List regions or regions (use with -threshold).",
"-color                   Generate a color image based on region size.",
"-series 1.3,4.7,0.05     Segment through a series of thresholds: start,finish,increment.",
" ",
"Actions for seeding:",
"-peaks 12,1.4,1,1        Find peaks within a kernel size, above a threshold,",
"                         with optional flooding and wrapping.",
"-threshold 3.4           Segment at a threshold.",
//"-addmodel                Add input model components as additional regions.",
" ",
"Actions for flooding:",
"-flood 0.6,-0.03         Segment by flooding to a low threshold,",
"                         with increment (default increment = range/100).",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
"-fill                    Fill in borders between regions or regions.",
"-multilevel              Flag for multi-level output for extracted segments.",
"-fill 0.02               Fill value for extraction: average (default), background, or value.",
" ",
"Input:",
"-Map map.pif             Replacement for the map specified in the model file.",
"-Mask mask.tif           Mask for flooding (use with -threshold, must be integer and same size as input image).",
"-Model mod.star          Seed model for flooding (use with -threshold).",
" ",
"Output:",
//"-output regions.star     Model of segments.",
"-Extract mask.pif        Extract regions marked by the model (see flag for multi-level).",
"-Region region.pif       Base name to generate numbered output file names.",
"-Postscript his.ps       Output plot file of level size histogram.",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setinvert(0);				// Flag to invert density
	int				seg_type(0);				// 0=none, 1=series, 2=flooding(peaks), 3=flooding(threshold)
	double	 		threshold_hi(0);			// Segmentation thresholds
	double	 		threshold_lo(0);
	double	 		threshold_inc(0);
	long			kernel_size(3);				// Kernel for peak definition
	int				flood(0);					// Flag to flood from peaks
	int				wrap(0);					// Flag to turn wrapping on for finding peaks
//	int				addmodel(0);
	int				list(0);
	int				color(0);
	int 			fill_borders(0);			// Borders between regions left
	Bstring			seed_file;					// Model file
	Bstring			paramfile;					// Input parameter file name
	Bstring			map_file;					// Map file
	Bstring			mask_file;					// Mask file
	Bstring			seg_file;					// Model of segments
	Bstring			ext_file;					// Output for extracted regions
	int				multilevel(0);				// Flag for multi-level extraction
	int 			fill_type(FILL_AVERAGE);	// Fill type for smaller images
	double			fill(0);
	Bstring			region_file;				// region file base name
	Bstring			ps_file;					// Poscript level size histogram plot

	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "series" ) {
			if ( ( i = curropt->values(threshold_hi, threshold_lo, threshold_inc ) ) < 1 )
				cerr << "-series: At least one threshold must be specified!" << endl;
			else {
				seg_type = 1;
				if ( i < 2 ) threshold_lo = threshold_hi;
				if ( i < 3 ) threshold_inc = threshold_lo - threshold_hi;
			}
		}
		if ( curropt->tag == "peaks" ) {
			if ( curropt->values(kernel_size, threshold_hi, flood, wrap ) < 1 )
				cerr << "-peaks: A kernel size must be specified!" << endl;
			else {
				seg_type = 2;
				if ( flood ) flood = 1;
			}
		}
		if ( curropt->tag == "threshold" ) {
			threshold_hi = curropt->value.real();
			seg_type = 3;
		}
		if ( curropt->tag == "flood" ) {
			if ( curropt->values(threshold_lo, threshold_inc ) < 1 )
				cerr << "-threshold: A lower threshold must be specified!" << endl;
			else
				flood = 2;
		}
		if ( curropt->tag == "list" ) list = 1;
		if ( curropt->tag == "color" ) color = 1;
//		if ( curropt->tag == "addmodel" ) addmodel = 1;
		if ( curropt->tag == "fill" ) fill_borders = 1;
		if ( curropt->tag == "multilevel" ) multilevel = 1;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "Model" )
			seed_file = curropt->filename();
		if ( curropt->tag == "Map" )
			map_file = curropt->filename();
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "output" )
			seg_file = curropt->filename();
		if ( curropt->tag == "Extract" )
			ext_file = curropt->filename();
		if ( curropt->tag == "Region" )
			region_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
    }
	option_kill(option);
        
	double			ti = timer_start();
	
    // Read the input file
	Bstring			filename(argv[optind++]);
	Bproject*		project = NULL;
	Bmodel*			model = NULL;
	Bmodel*			seed = NULL;
	
	if ( file_type(filename) == Micrograph ) {
		project = read_project(filename);
		if ( project->rec ) filename = project->rec->frec;
	} else if ( file_type(filename) == Model ) {
		model = read_model(filename, paramfile);
		filename = model->mapfile();
	}
	
	if ( map_file.length() ) filename = map_file;
	if ( model ) model->mapfile(filename.str());

	if ( filename.length() < 1 ) {
		cerr << "Error: No image file name specified!" << endl;
		bexit(-1);
	}

	if ( !model && ext_file.length() ) {
		cerr << "Error: A model file is required to extract segments!" << endl;
		bexit(-1);
	}
	
    Bimage* 	p = read_img(filename, 1, -1);
	if ( !p ) {
		cerr << "Error: No image file " << filename << " read!" << endl;
		bexit(-1);
	}
	
	if ( nudatatype == Unknown_Type ) {
		if ( color )
			nudatatype = UCharacter;
		else
			nudatatype = p->data_type();
	}
	
	if ( setinvert ) p->invert();
	
	Bimage*			pmask = p;
	Bimage*			pseed = NULL;
	Bimage*			pex = NULL;
	Bplot*			plot = NULL;
	Vector3<double>	origin(p->image->origin());
	Vector3<long>	size = {p->sizeX(), p->sizeY(), p->sizeZ()};
	Vector3<double>	sam(p->sampling(0));

	// Mask generation or reading
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, -1);
	} else if ( seed_file.length() ) {
		seed = read_model(seed_file, paramfile);
		pseed = img_from_model(seed, origin, size, sam, 1);
		pmask = pseed->regions(0.5, 1);
		delete pseed;
	} else if ( seg_type == 3 ) {
		pmask = p->regions(threshold_hi, 0);
	} else if ( seg_type == 2 ) {
		pmask = p->region_peaks(kernel_size, threshold_hi, flood, wrap);
	}
	
//	if ( pmask && model && addmodel )
//		img_add_model_to_mask(pmask, model);

	if ( seg_type == 1 ) {
		p->region_threshold_series(threshold_hi, threshold_lo, threshold_inc);
	} else if ( flood == 2 ) {
		p->region_flood(pmask, threshold_hi, threshold_lo, threshold_inc, fill_borders);
	}
	
	if ( list ) p->level_masked_stats(pmask);
	
	// Extract masked regions
	if ( pmask && model && ext_file.length() ) {
		pex = img_extract_segments_using_model(pmask, model, multilevel);
		write_img(ext_file, pex, 0);
		delete pex;
	}

	if ( pmask && region_file.length() ) {
		pex = p->level_mask_extract(pmask, fill_type, fill);
		pex->change_type(nudatatype);
		write_img(region_file, pex, 0);
		delete pex;
	}

	if ( ps_file.length() ) {
		plot = pmask->levelmask_size_histogram();
		ps_plot(ps_file, plot);
		delete plot;
	}
	if ( color ) {
		if ( p != pmask ) delete p;
		p = pmask->levelmask_color_by_size();
		delete pmask;
		pmask = p;
	}
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		pmask->change_type(nudatatype);
     	write_img(argv[optind], pmask, 0);
	}
	
	if ( p != pmask ) delete p;
	delete pmask;
	project_kill(project);
	model_kill(model);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

