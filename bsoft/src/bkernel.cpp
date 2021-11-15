/**
@file	bkernel.cpp
@brief	Program to calculate filter kernels.
@author	Bernard Heymann
@date	Created: 20051102
@date	Modified: 20150910
**/

#include "rwimg.h"
#include "FSI_Kernel.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bkernel [options] output.krn",
"-----------------------------------",
"Calculates filter kernels.",
"The default kernel just calculates the mean.",
"The different types of kernel are mutually exclusive.",
" ",
"Actions:",
"-gaussian 1.7,8.4        Gaussian kernel with the indicated sigma and maximum.",
"-log 2.5,5               Laplacian of gaussian kernel with the indicated sigma and maximum.",
"-frequency 8,2           Frequency space kernel: width and power.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-size 5,7,9              Size of kernel (voxels).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	Vector3<long>	size;		// Kernel size
	double			sigma(0);			// Gaussian sigma
	double			log_sigma(0);		// LoG sigma
	double			max(1);				// Gaussian or LoG maximum
	int				kernel_width(0);	// Frequency space kernel width
	int				kernel_power(2);	// Frequency space kernel power
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "gaussian" )
        	if ( curropt->values(sigma, max) < 1 )
				cerr << "-gaussian: A sigma value must be specified!" << endl;
		if ( curropt->tag == "log" )
        	if ( curropt->values(log_sigma, max) < 1 )
				cerr << "-log: A sigma value must be specified!" << endl;
		if ( curropt->tag == "frequency" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-frequency: At least the kernel size must be specified!" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
    }
	option_kill(option);
	
	double		ti = timer_start();

	int			side;
	FSI_Kernel	ker;

	if ( size.volume() < 1 ) {
		if ( log_sigma > 0.001 ) {
			side = 2*( (int) (4*log_sigma + 0.9) ) + 1;
			size = {side, side, side};
		} else if ( sigma > 0.001 ) {
			side = 2*( (int) (4*sigma + 0.9) ) + 1;
			size = {side, side, side};
		} else if ( kernel_width ) {
			ker = FSI_Kernel(kernel_width, kernel_power);
			size = {ker.width(), ker.divisions() + 1, 1};
		} else {
			size = {3, 3, 3};
		}
	}
		
	Bimage*		p = new Bimage(Float, TSimple, size, 1);
	p->origin(p->default_origin());

	if ( log_sigma ) {
		p->kernel_laplacian_of_gaussian(log_sigma, max);
	} else if ( sigma > 0.001 ) {
		p->kernel_gaussian(sigma, max);
	} else if ( kernel_width ) {
		for ( long i=0; i<p->data_size(); ++i ) p->set(i, ker[i]);
	} else {
		p->fill(1);
	}
	
	if ( p && optind < argc ) {
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

