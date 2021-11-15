/**
@file	bhisto.cpp
@brief	Image histograms
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210202
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bhisto [options] input.img output.img",
"--------------------------------------------",
"Calculates the histogram of an image and associated statistics.",
" ",
"Actions:",
"-bins 256                Calculate histogram over a number of bins.",
"-counts                  Calculate a counts histogram and adjust the image.",
"-fit gauss,2             Do a fit to the histogram (gaussian, poisson, counts).",
"                         Second argument indicates number of gaussians.",
"                         Counts option will convert the image to discrete counts.",
"-otsu                    Calculate the Otsu inter-set variance.",
"-multi 5                 Calculate multiple thresholds for the given number of clusters.",
"-percentiles             Calculate the percentiles for an image.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-select 1                Select an image (default all).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
" ",
"Output:",
"-Postscript histo.ps     Postscript output file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
    long     		bins(0);    	    	// Number of bins, also selects histogram
    long     		counts(0);    	    	// Flag to calculate counts histogram
	int 			setimg(-1); 			// Select all images
	Bstring			setfit;					// Fitting specification: gau,poi,cou
	int				fitprm(1);				// Fit parameter: gaussians or counting flags
	int				percentiles(0);			// Calculate the percentiles for an image
	int				multi(0);				// Number of histogram clusters
	Vector3<double>	sam;					// Units for the three axes (A/pixel)
	Bstring			psfile;					// Postscript output file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "bins" ) {
        	if ( ( bins = curropt->value.integer() ) < 1 )
				cerr << "-bins: The number of bins must be specified!" << endl;
			else if ( bins < 2 ) {
				bins = 10;
				cerr << "-bins: Too few bins, changed to 10!" << endl;
			} else if ( bins > 100000 ) {
				bins = 1000;
				cerr << "-bins: Too many bins, changed to 1000!" << endl;
			}
		}
		if ( curropt->tag == "counts" ) counts = 1;
		if ( curropt->tag == "fit" ) {
			setfit = curropt->value;
			if ( setfit.contains(",") ) {
				setfit = setfit.pre(',');
				fitprm = curropt->value.post(',').integer();
			}
		}
		if ( curropt->tag == "otsu" ) setfit = "otsu";
		if ( curropt->tag == "multi" )
        	if ( ( multi = curropt->value.integer() ) < 2 )
				cerr << "-multi: The number of clusters must be specified!" << endl;
		if ( curropt->tag == "percentiles" ) percentiles = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "select" )
        	if ( ( setimg = curropt->value.integer() ) < 0 )
				cerr << "-select: An image number must be specified!" << endl;
		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read image file
	int 		dataflag(0);
	if ( bins || counts || percentiles || optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, setimg);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( sam.volume() ) p->sampling(sam);
	
	Bplot*			plot = NULL;
	
    if ( bins ) {
		if ( setfit.length() ) {
			if ( setfit[0] == 'g' ) {
				plot = p->histogram_gauss_plot(bins, fitprm);
			} else if ( setfit[0] == 'p' ) {
				plot = p->histogram_poisson_fit(bins, 1);
	 		} else if ( setfit[0] == 'o' ) {
				plot = p->histogram_otsu_variance(bins);
	 		} else if ( setfit[0] == 'c' ) {
				plot = p->histogram_counts(3);
			} else {
				cerr << "Incorrect fitting specification: " << setfit << endl;
			}
		} else if ( multi > 1 ) {
			p->histogram_multi_thresholds(bins, multi);
		} else {
			plot = p->histogram(bins);
		}
	} else if ( counts ) {
		plot = p->histogram_counts(3);
	} else if ( percentiles ) {
		plot = p->percentiles();
	}

	if ( psfile.length() && plot )
		ps_plot(psfile, plot);
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	if ( plot ) delete plot;
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

