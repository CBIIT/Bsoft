/**
@file	bfilwidth.cpp
@brief	Program to calculate moving average images.
@author Bernard Heymann
@date	Created: 20130204
@date	Modified: 20150618
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
"Usage: bfilwidth [options] img.pif out.pif",
"------------------------------------------",
"Analyzes the periodicity of a filament from width changes.",
"The filament must be aligned with the y axis and with positive contrast.",
" ",
"Actions:",
"-invert                  Invert density in the image before other operations.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation before other operations.",
"-average 34,1,1          Moving average window size (voxels).",
"-period 100,32,1         Periodic frame size (voxels).",
"-width 25                Window to fit a moving line for estimating the filament width.",
"-density 30,120          Density with filament width and number of histogram bins.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-resolution 15.6         Resolution limit for cross-correlation.",
"-limits 35,57            Limits on the filament width (default 0 - width of image).",
" ",
"Input:",
"-reference file.map      2D reference for cross-correlation.",
" ",
"Output:",
"-Postscript plot.ps      Plot of filament widths.",
" ",
NULL
};

int 		main(int argc, char** argv)
{
	// Initialize variables
//	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setinvert(0);				// Flag to invert density
	double			nuavg(0), nustd(0); 		// Values for rescaling
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;			// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<long>	average_window;				// Size of moving average window
	Vector3<long>	period;						// Size of periodic frame
	long			width(0);					// Moving window for filament width estimation
	long			filwidth(0), bins(100);		// Filament width and bins for density calculation
	double	 		res_lo(10000);				// High resolution limits
	double	 		res_hi(0);					// Low resolution limits
	int				lim_lo(0), lim_hi(0);		// Limits on the filament width
	Bstring			ref_file;					// Reference file name
	Bstring			psfile;						// Postscript output file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
//		if ( curropt->tag == "datatype" )
//			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "average" )
			average_window = curropt->size();
		if ( curropt->tag == "period" )
			period = curropt->size();
		if ( curropt->tag == "width" )
			if ( ( width = curropt->value.integer() ) < 1 )
				cerr << "-width: A window size must be specified!" << endl;
		if ( curropt->tag == "density" )
			if ( curropt->values(filwidth, bins) < 1 )
				cerr << "-density: A filament width must be specified!" << endl;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "resolution" )
			if ( ( res_hi = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "limits" )
			if ( curropt->values(lim_lo, lim_hi) < 2 )
				cerr << "-limits: Both low and high limits must be specified!" << endl;
 		if ( curropt->tag == "reference" )
			ref_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
	}
	option_kill(option);
	
	double		ti = timer_start();

	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);

	if ( sam.volume() > 0 ) p->sampling(sam);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}

	if ( setinvert ) p->invert();

	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( average_window.volume() ) p->filter_average(average_window);

	if ( period.volume() ) {
		Bimage*		pp = p->periodic_averaging(period);
		delete p;
		p = pp;
	}
	
	Bimage*		pref = NULL;
	Bimage*		pmask = NULL;
	if ( ref_file.length() ) {
    	pref = read_img(ref_file, 1, -1);
		if ( !pref ) bexit(-1);
		p->find_shift(pref, pmask, res_hi, res_lo, p->sizeY()/2, 0, 1);
		if ( verbose )
			cout << "New origin = " << p->image->origin() << endl;
		p->shift_wrap(pref->image->origin() - p->image->origin());
		p->origin(pref->image->origin());
		delete pref;
	}

	Bimage*		pd = NULL;
	Bplot*		plot = NULL;
	if ( width ) plot = p->filament_width(width, lim_lo, lim_hi);
	else if ( filwidth ) {
		pd = p->filament_density(filwidth);
		plot = pd->histogram_gauss_plot(bins, 1);
		if ( verbose )
			cout << "Gaussain fit: average = " << pd->average() << " standard deviation = " << pd->standard_deviation() << endl;
		delete pd;
	}

	if ( psfile.length() && plot )
		ps_plot(psfile, plot);

	if ( argc > optind ) {
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete plot;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

