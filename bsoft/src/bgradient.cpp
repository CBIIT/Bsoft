/**
@file	bgradient.cpp
@brief	Calculating image gradients.
@author Bernard Heymann
@date	Created: 20201224
@date	Modified: 20210302
**/

#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bimage*		img_aniso_average(Bimage* p, long iter, long ksize, double w);

// Usage assistance
const char* use[] = {
" ",
"Usage: bgradient [options] input.img output.img",
"-----------------------------------------------",
"Calculates image gradients.",
" ",
"Actions:",
"-gaussian 2.4,5.1,20.4   Anisotropic gaussian filter (one value sets all).",
"-gradient 3x3            Gradient type: cd, 3x3, freq.",
"-magnitude               Convert gradient vectors to lengths.",
"-anisotropic 10,0.5      Anisotropic gradient smoothing: iterations and weight.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for the mask (default from input image).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype = Float;		// Conversion to new type
	Vector3<double> sampling;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;					// Mask origin
	int				set_origin(0); 			// Flag to set origin
	Vector3<double>	sigma;					// Gaussian sigma values
	int				gtype(0);				// Gradient type
	int				mag(0);					// Magnitude flag
	long			aniso_iter(0);			// Anisotropic gradient iterations
	double			aniso_k(1);				// Anisotropic kernel size
	double			aniso_w(1);				// Anisotropic gradient weight
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "gaussian" )
			sigma = curropt->scale();
		if ( curropt->tag == "gradient" ) {
			if ( curropt->value[0] == 'c' ) gtype = 1;
			if ( curropt->value[0] == '3' ) gtype = 2;
			if ( curropt->value[0] == 'f' ) gtype = 3;
		}
		if ( curropt->tag == "magnitude" ) mag = 1;
		if ( curropt->tag == "anisotropic" )
			if ( curropt->values(aniso_iter, aniso_w) < 1 )
				cerr << "-anisotropic: A number of iterations must be specified!" << endl;
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
    }
	option_kill(option);
	
	double		ti = timer_start();

    Bimage* 	p = read_img(argv[optind++], 1, -1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();

	if ( sampling.volume() > 0 ) p->sampling(sampling);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	Bimage*		pg = NULL;
	
	if ( gtype < 3 &&  sigma.length() ) {
		if ( p->fourier_type() == NoTransform ) p->fft();
		p->fspace_weigh_gaussian(0, sigma);
		p->fft_back();
	}
	
	if ( gtype == 1 )
		pg = p->gradient();
	else if ( gtype == 2 )
		pg = p->gradient3x3();
	else if ( gtype == 3 )
		pg = p->fspace_gradient(sigma);
	
	if ( pg && mag ) pg->vector_to_simple();

	if ( !pg && aniso_iter )
		pg = img_aniso_average(p, aniso_iter, aniso_k, aniso_w);

	if ( pg && optind < argc ) {
		pg->change_type(nudatatype);
		write_img(argv[optind], pg, 0);
	}
	
	delete pg;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


Bimage*		img_aniso_average(Bimage* p, long iter, long ksize, double w)
{
	Bimage*		pt = NULL;
	long		i;
	double		d(-1);

	if ( verbose ) {
		cout << "Anisotropic gradient smoothing:" << endl;
		cout << "Weight:                         " << w << endl << endl;
	}

	if ( verbose )
		cout << "Iter\tStDev" << endl;
	for ( i=0; i<iter && d<0; ++i ) {
		if ( verbose )
			cout << i << tab << p->standard_deviation() << endl;
		pt = p->aniso_average(ksize, w);
		d = pt->standard_deviation() - p->standard_deviation();
		delete p;
		p = pt;
	}
	if ( verbose )
		cout << i << tab << pt->standard_deviation() << endl << endl;

	return pt;
}
