/**
@file	bseries.cpp
@brief	Program to align and analyze series of images
@author	Bernard Heymann
@date	Created: 20040407
@date	Modified: 20210915
**/

#include "mg_processing.h"
#include "mg_align.h"
#include "mg_extract.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bseries [options] input.star [input.star]",
"------------------------------------------------",
"Aligns series of micrographs using consecutive pairwise cross-correlation.",
" ",
"Actions:",
"-align 15,3,2            Align to reference image (first=0), moving sum and interval (default 1,1).",
"-snr 5                   Estimate SSNR over a summation window.",
"-envelope 5,4.5          Calculate an envelope due to shift: frame window and resolution.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for tilt axis (default 0,0,0).",
"-ratio 2.5               Radial sampling to Cartesian sampling ratio (default 1).",
" ",
#include "use_dose.inc"
" ",
"Parameters for alignment:",
"-frames                  Flag to align micrograph frames.",
"-counts                  Flag to rescale images based on their counts histogram.",
"-resolution 900,300      High and low resolution limits for cross-correlation (default 0.1,1000 angstrom).",
"-bin 3                   Binning by the given kernel size to speed up alignment.",
"-shiftlimit 3.5          Limit on origin shift relative to nominal center (default 10% of box edge size).",
"-edge 23,12              Smooth the edge to a given width, with gaussian decay of a given width.",
"-subset 2-8,12           Subset of micrographs to align and average.",
" ",
"Input:",
"-Gainreference gr.tif    Gain reference to correct the input micrographs.",
"-Mask mask.tif           Mask file to use for cross-correlation.",
" ",
"Output:",
"-output file.star        Output parameter file, if -average, also output average parameter file.",
"-write sum               Output either aligned frames, sum or avg.",
"-image file.pif          Output multi-image file with aligned images.",
"-average file.pif        Output average of aligned images.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	DataType 		datatype(Unknown_Type);		// Conversion to new type
	int				ref_img(-1);				// Reference image for alignment, <0 means don't align
	int				window(1), step(1);			// Moving sum window for alignment
	int				frames(0);					// Flag to align micrograph frames
	int				flags(0);					// Flags: 1=rescale based on histogram; 2=weigh by dose; 4=write aligned frames; 8=write frame sum
	double			sampling_ratio(1);			// For SNR estimation
	double			snr_window(0);				// Summing window for SNR estimation
	int				snr_prog(0);				// Flag for progressive summing
	long			shift_window(1);			// Frame window for shift envelope calculation
	double			shift_res(10);				// Resolution for shift envelope calculation
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	double			hi_res(0), lo_res(1e10);	// Default resolution range
	double			shift_limit(-1);			// Maximum shift from nominal image origin
	Vector3<double>	origin;						// Tilt axis origin
	double			edge_width(0), gauss_width(0);	// Edge parameters
	int 			fill_type(FILL_BACKGROUND);
	double			fill(0);
	long		 	bin(1);						// Binning before alignment and analysis
	Bstring			subset;						// Subset of micrographs to average
	JSvalue			dose_frac(JSobject);		// Container for dose fractionation parameters
	Bstring			paramfile;					// Output parameter file
    Bstring			grfile;						// Gain reference
    Bstring			maskfile;					// Mask to use for cross-correlation
	Bstring			imgfile;					// Output multi-image file
	Bstring			avgfile;					// Output average image file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			datatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "align" )
			if ( curropt->values(ref_img, window, step) < 1 )
 				cerr << "-align: A reference image number must be specified." << endl;
		if ( curropt->tag == "snr" ) {
			if ( ( snr_window = curropt->value.integer() ) < 2 )
 				cerr << "-snr: A window larger than 1 must be specified." << endl;
			if ( snr_prog ) snr_prog = 2;
		}
		if ( curropt->tag == "envelope" )
			if ( curropt->values(shift_window, shift_res) < 1 )
 				cerr << "-envelope: A frame window size must be specified." << endl;
		if ( curropt->tag == "frames" )
			frames = 1;
		if ( curropt->tag == "counts" )
			flags |= 1;
#include "dose.inc"
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( hi_res > lo_res )
				swap(hi_res, lo_res);
        }
		if ( curropt->tag == "shiftlimit" )
			if ( ( shift_limit = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A maximum shift in pixels must be specified!" << endl;
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "ratio" )
			if ( ( sampling_ratio = curropt->value.real() ) < 0.1 )
				cerr << "-ratio: A ratio must be specified!" << endl;
		if ( curropt->tag == "edge" )
    	    if ( curropt->values(edge_width, gauss_width) < 1 )
				cerr << "-edge: An edge width must be specified." << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: An ineteger greater than zero must be specified!" << endl;
		if ( curropt->tag == "subset" )
			subset = curropt->value;
		if ( curropt->tag == "Gainreference" )
			grfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
		if ( curropt->tag == "write" ) {
			if ( curropt->value[0] == 'f' ) flags |= 4;
			if ( curropt->value[0] == 's' ) flags |= 8;
			if ( curropt->value[0] == 'a' ) flags |= 8;
		}
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
		if ( curropt->tag == "average" )
			avgfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);		
	string_kill(file_list);

	if ( !project ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

//	cout << dose_frac << endl;
	if ( dose_frac.size() )
		project_set_dose(project, dose_frac);

	Bimage*			pgr = NULL;
	if ( grfile.length() )
		pgr = read_img(grfile, 1, 0);
	
	Bimage*			pmask = NULL;
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, 0);
		
	if ( ref_img > -1 ) {
		if ( frames )
			project_align_frames(project, ref_img, window, step, pgr, pmask, origin, hi_res, lo_res,
				shift_limit, edge_width, gauss_width, bin, subset, flags);
		else
			project_align_series(project, ref_img, pgr, pmask, origin, hi_res, lo_res,
				shift_limit, edge_width, gauss_width, bin, subset, flags);
		delete pgr;
		pgr = NULL;
	}
	
	if ( snr_window > 1 )
		project_frames_snr(project, hi_res, snr_window, subset, sampling_ratio, (flags&1));
	
	if ( shift_window )
		project_frame_shift_analysis(project, shift_window, shift_res);

	if ( !frames && imgfile.length() )
		project_write_aligned_images(project, pgr, imgfile, datatype);
	
	if ( avgfile.length() )
		project_write_aligned_averages(project, pgr, avgfile, datatype, subset);

	if ( frames && flags & 8 )
		project_write_frame_sums(project, pgr, datatype, subset, sampling_ratio, flags);

    if ( paramfile.length() )
		write_project(paramfile, project, 0, 0);
	
	project_kill(project);
	delete pgr;
	delete pmask;

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

