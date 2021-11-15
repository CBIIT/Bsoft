/**
@file	bradsec.cpp
@brief	Generates radial sections from 3D images
@author Bernard Heymann
@date	Created: 20000620
@date	Modified: 20200509
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bplot*		img_plot_radial(Bimage* p, Bstring& title);

// Usage assistance
const char* use[] = {
" ",
"Usage: bradsec [options] input.img output.img",
"---------------------------------------------",
"Generates symmetry-adjusted radial sections from 3D images.",
"The symmetry for O must be either O-3 or O-4.",
"The symmetry for I must be either I-3 or I-5.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (default from image; a single value can be given).",
"-origin 0,-10,30         Set the origin for radial profile (default image origin).",
"-radii 14,123,2.68       Start, end and step size for radial shell calculation (default from input).",
"-symmetry C5             Point group symmetry.",
"-fraction 0.3            Spherical fraction (1=spherical shells, 0=symmetry-related shells).",
"-fill 127                Fill value for excluded region: average (default), background, or value.",
" ",
"Output:",
"-Postscript radsec.ps    Symmetry-adjusted profile.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    					// Sampling
	Vector3<double>	origin;						// Origin
	int				set_origin(0);				// Flag to set origin
	double			rad_start(0), rad_end(0);	// Radial limits
	double			rad_step(0);				// Radial step size
	Bsymmetry		sym;					// Point group
	double			spherical_fraction(0);		// Spherical fraction
	double			fill(0);	 				// Fill value for excluded region
	int 			fill_type(FILL_AVERAGE);	// Fill type for excluded region
	Bstring			outfile;					// Output image file
	Bstring			psfile;						// Output profile file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "radii" )
			if ( curropt->values(rad_start, rad_end, rad_step) < 2 )
				cerr << "-radii: A range and step size must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "fraction" )
			if ( ( spherical_fraction = curropt->value.real() ) <= 0 )
				cerr << "-fraction: A fraction must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	// Get the output image file name
	if ( optind < argc ) outfile = argv[optind];
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();		// Preserve the old type

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}

	Bimage*		prad = p->radial_sections(rad_start, rad_end, rad_step, spherical_fraction, sym, fill_type, fill);
	if ( !prad ) {
		cerr << "Error: No output image generated!" << endl;
		bexit(-1);
	}

	if ( psfile.length() ) {
		Bstring			title = "Symmetry-adjusted radial average: " + sym.label();
		Bimage*			prof = p->radial_symmetry_adjusted(rad_start, rad_end, rad_step, spherical_fraction, sym);
		Bplot*			plot = img_plot_radial(prof, title);
		if ( plot ) {
			ps_plot(psfile, plot);
			delete plot;
		}
		delete prof;
	}

	if ( outfile.length() ) {
		prad->change_type(nudatatype);
		write_img(outfile, prad, 0);
	}

	delete p;
	delete prad;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


Bplot*		img_plot_radial(Bimage* p, Bstring& title)
{
	long		i, k, ncol(p->images()*p->channels() + 1);
	
	Bplot*		plot = new Bplot(1, p->sizeX(), ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	plot->page(0).column(0).label("Distance(A)");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; ++i ) {
		plot->page(0).column(i).number(i);
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).label("Average");
		plot->page(0).column(i).axis(3);
	}
//	plot->page(0).axis(1).min(minrad);
//	plot->page(0).axis(1).max(maxrad);
	plot->page(0).axis(3).min(p->minimum());
	plot->page(0).axis(3).max(p->maximum());

	for ( i=0; i<plot->rows(); ++i )
		(*plot)[i] = p->sampling(0)[0]*i;

	for ( i=0, k=plot->rows(); i<p->data_size(); ++i, ++k )
		(*plot)[k] = (*p)[i];

	return plot;
}
