/**
@file	bcyl.cpp
@brief	Program to manipulate cylindrical images.
@author Bernard Heymann
@date	Created: 20090309
@date	Modified: 20160326
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			img_plot_max_vs_angle(Bimage* p);

// Usage assistance
const char* use[] = {
" ",
"Usage: bcyl [options] input.img output.img",
"------------------------------------------",
"Unwraps an image to its cylindrical projection.",
" ",
"Actions:",
"-symmetrize              Symmetrizes a map cylindrically around the origin.",
"-unwrap 1.5              Unwraps an image to its cylindrical projection with the given angle.",
"-autocorrelate           Autocorrelates the unwrapped map and outputs the autocorrelation map.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation (after -truncate).",
"-shells                  Calculate cylindrical shells projecting along the y-axis.",
"-reslice yzx             Reslice (switch axes) after unwrapping.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 0,-10.5,30       Set the origin.",
"-resolution 4.5,130      Resolution range for correlation (default 0 - 1e6 angstrom).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	int				symmetrize(0);				// Flag to symmetrize cylindrically
	double			unwrap_angle(0);			// Angle for cylindrical unwrapping
	int				autocorrelate(0);			// Flag to autocorrelate the unwrapped map
	int				shells(0);					// Flag to calculate cylindrical shells
	int 			setrescale(0);
	double			nuavg(0), nustd(0);			// Rescaling to average and stdev
	int 			setreslice(0);
	Bstring			order("xyz");
	DataType 		nudatatype = Unknown_Type;	// Conversion to new type
	Vector3<double>	origin  {0,0,0};			// Origin
	int				set_origin(0);				// Flag to set origin
	double			hires(0), lores(1e6);		// Limiting resolution range (hires must be > 0 to be set)
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "symmetrize" ) symmetrize = 1;
		if ( curropt->tag == "unwrap" ) {
        	if ( ( unwrap_angle = curropt->value.real() ) < 0.1 )
				cerr << "-unwrap: An angle must be specified!" << endl;
			else
				unwrap_angle *= M_PI/180.0;
        }
		if ( curropt->tag == "autocorrelate" ) autocorrelate = 1;
		if ( curropt->tag == "shells" ) shells = 1;
		if ( curropt->tag == "rescale" ) {
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else
				setrescale = 1;
			if ( nustd <= 0 ) setrescale = 0;
		}
		if ( curropt->tag == "reslice" ) {
			if ( curropt->value.length() < 3 )
				cerr << "-reslice: The order of x, y and z must be specified!" << endl;
			else {
				order = curropt->value;
				setreslice = 1;
			}
		}
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
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read image file
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
	
	if ( symmetrize ) p->symmetrize_cylinder();
	
	Bimage*		pcyl = NULL;
	if ( unwrap_angle > 0 ) {
		pcyl = p->cartesian_to_cylindrical(p->sizeX()/2, TWOPI/unwrap_angle);
		delete p;
		p = pcyl;
		if ( autocorrelate ) {
			p->auto_correlate(hires, lores);
			img_plot_max_vs_angle(p);
		}
	}

	if ( shells ) p->cylindrical_shells();
	
	if ( optind < argc ) {
		if ( setrescale ) p->rescale_to_avg_std(nuavg, nustd);
		if ( setreslice ) p->reslice(order);
	
		p->change_type(nudatatype);
	
		// Copy the command line to the label field of the image structure
	
		write_img(argv[optind], p, 0);
	}
		
	delete p;
	order = 0;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

int			img_plot_max_vs_angle(Bimage* p)
{
	long			i, x, y(0), z, z1, z2, zm, n, dz(10);
	double			a, v, cc_max;
	
	for ( n=0; n<p->images(); n++ ) {
		z1 = 0; z2 = dz;
		cout << "Image " << n+1 << ":" << endl << "Angle\tz\tCC" << endl;
		for ( x=0; x<p->sizeX(); x++ ) {
			a = x*TWOPI/p->sizeX();
			for ( z=z1, zm=0, cc_max = p->minimum(); z<z2; z++ ) {
				i = p->index(0, x, y, z, n);
				v = (*p)[i];
				if ( cc_max < v ) {
					cc_max = v;
					zm = z;
				}
			}
			z1 = (zm > dz)? zm - dz: 0;
			z2 = zm + dz;
			cout << a*180.0/M_PI << tab << zm*p->sampling(0)[2] << tab << cc_max << endl;
		}
		cout << endl;
	}
	
	return 0;
}

