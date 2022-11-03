/**
@file	bresolve.cpp
@brief	Calculate resolution estimates and Fourier shell statistics
@author Bernard Heymann
@date	Created: 20000612
@date	Modified: 20220324
**/

#include "rwimg.h"
#include "rwFSC_XML.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bresolve [options] input1.img input2.img",
"-----------------------------------------------",
"Determines the resolution (FSC and DPR) from two images.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-variancemask 45         Mask from local variance of reference with given kernel size.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (angstrom/pixel, a single value sets all three).",
"-resolution 15           Resolution limit for output (angstrom).",
"-aligned 25              Resolution limit for alignment (angstrom).",
"-fsccut 0.8,0.5,0.143    Resolution cutoffs for FSC.",
"-dprcut 45,70            Resolution cutoffs for DPR.",
"-ratio 2.5               Radial sampling to Cartesian sampling ratio (default 1).",
" ",
"Input:",
"-map file.map            3D reference file name.",
"-mask mask.tif           Real space mask (same size as reference).",
" ",
"Output:",
"-output file.star        Output STAR format file.",
"-Postscript file.ps      Postscript output filename.",
"-XML file.xml            XML output filename.",
"-tsv file.tsv            Tab-separated values output filename.",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize all settings
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	long			var_kernel(0);				// Local variance kernel to generate mask
	Vector3<double>	sam;    					// Pixel size
	double			hi_res(0);					// Resolution limit for output
	double			aligned_res(0);				// Resolution limit for alignment
	vector<double>	fsccut{0.3,0,0,0};			// FSC cutoff values
	vector<double>	dprcut{0,0,0,0};			// DPR cutoff values
	double			sampling_ratio(1);			// Sampling ratio in frequency space
	Bstring			map_file;					// Reference map
	Bstring			ps_file;					// Output postscript file
	Bstring			xml_file;					// Output XML file
	Bstring			tsv_file;					// Output tsv file
	Bstring			real_mask_file;				// Real space mask file
	    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "variancemask" )
			if ( ( var_kernel = curropt->value.integer() ) < 3 )
				cerr << "-variancemask: A kernel edge size must be specified!" << endl;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( hi_res = curropt->value.real() ) < 1 )
				cerr << "-resolution: A resolution must be specified!" << endl;
		if ( curropt->tag == "aligned" )
			if ( ( aligned_res = curropt->value.real() ) < 1 )
				cerr << "-aligned: A resolution must be specified!" << endl;
		if ( curropt->tag == "fsccut" )
			if ( curropt->values(fsccut[0], fsccut[1], fsccut[2], fsccut[3]) < 1 )
				cerr << "-fsccut: At least one FSC threshold must be specified!" << endl;
		if ( curropt->tag == "dprcut" )
			if ( curropt->values(dprcut[0], dprcut[1], dprcut[2], dprcut[3]) < 1 )
				cerr << "-dprcut: A DPR threshold must be specified!" << endl;
		if ( curropt->tag == "ratio" )
			if ( ( sampling_ratio = curropt->value.real() ) < 0.1 )
				cerr << "-ratio: A ratio must be specified!" << endl;
		if ( curropt->tag == "map" )
			map_file = curropt->filename();
		if ( curropt->tag == "mask" )
			real_mask_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
		if ( curropt->tag == "XML" )
			xml_file = curropt->filename();
		if ( curropt->tag == "tsv" )
			tsv_file = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;

	
	Bimage*			p = read_img(argv[optind++], 1, -1);

	if ( map_file.length() < 1 ) map_file = argv[optind++];
	
	Bimage*			pr = read_img(map_file, 1, -1);
	
	if ( p->images() != pr->images() ) {
		cerr << "Error: The number of images should match!" << endl;
		bexit(-1);
	}
	
	Bplot*			plot = NULL;
		
	if ( p->compound_type() == TSimple ) {
	
		if ( nustd ) {
			p->rescale_to_avg_std(nuavg, nustd);
			pr->rescale_to_avg_std(nuavg, nustd);
		}
		
		Bimage*			pmask = NULL;
		if ( real_mask_file.length() ) {
			pmask = read_img(real_mask_file, 1, 0);
		} else if ( var_kernel > 1 ) {
			pmask = pr->variance_mask(var_kernel);
			pmask->filter_average(11);
			write_img("mask.pif", pmask, 0);
		}
		
		if ( pmask ) {
			p->multiply(pmask);
			pr->multiply(pmask);
			delete pmask;
		}

		Bimage*			pc = p->resolution_prepare(pr);
		delete p;
		delete pr;
		
		if ( pc->images() == 1 ) {
			plot = pc->fsc_dpr(hi_res, sampling_ratio, (dprcut[0] == 0));
			plot->resolution_display(fsccut, dprcut);
		} else {
			plot = pc->fsc(hi_res, sampling_ratio, fsccut);
		}
		
		delete pc;
		
	} else if ( p->compound_type() == TComplex ) {
	
		plot = p->fsc(pr, hi_res, sampling_ratio);
		plot->resolution_display(fsccut, dprcut);
		
		delete p;
		delete pr;
	}
	
	if ( ps_file.length() ) ps_plot(ps_file, plot);
	if ( xml_file.length() ) xml_write_fsc(xml_file, plot);
	if ( tsv_file.length() ) plot->write_tsv(tsv_file);
	delete plot;
	
#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

