/**
@file	brandom.cpp
@brief	Program to generate random images
@author Bernard Heymann
@date	Created: 19990703
@date	Modified: 20220517
**/

#include "rwimg.h"
#include "random_numbers.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: brandom [options] [input.img] output.img",
"-----------------------------------------------",
"Generates random content images.",
"If an input image is not found a random output image is generated.",
"The default output data type is floating point.",
" ",
"Actions:",
"-type gauss              Distribution type (default: uniform, ",
"                         others: gaussian, poisson, logistical, spectral).",
"-distance gauss          Distance distribution type (default: uniform, ",
"                         others: gaussian)",
"-rescale -0.1,5.2        Rescale output data to average and standard deviation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Data type (default: byte).",
"-images 12               Number of images (default: 1).",
"-size 10,50,8            Size (default: 256,256,1).",
"-sampling 1.5,1.5,1.5    Sampling (default: 1,1,1; a single value can be given).",
"-origin 0,-10,30         Set the origin (default image origin).",
"-minmax -1.2,5.6         Set minimum and maximum of input and uniform random image (default: -1.732,1.732).",
"-avgstd 0.5,2.3          Set average and standard deviation of input and non-uniform random image (default: 0,1).",
"-alpha 2                 Spectral noise decay (default: white=0, pink=1, red/brown=2, blue=-1, violet=-2).",
"-snr 1.4,63              Signal-to-noise ratio and foreground radius (default: 1).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);		// Data type to be set
	Vector3<long>	size(256,256,1);
	int				set_origin(0);					// Flag to set origin
	Vector3<double> origin;							// Origin
	int 			set_sampling(0);
	Vector3<double> sampling(1,1,1);				// Units for the three axes (A/pixel)
	long 			nimg(1);						// Number of images
	double			min(-sqrt(3.0)), max(sqrt(3.0));
	double			avg(0), std(1);					// Rescaling input or random image to average and stdev
	int				set_avgstd(0);					// Flag to set average and standard deviation
	double			alpha(0);						// Spectral noise decay
	double			finalavg(0), finalstd(0);		// Rescaling output to average and stdev
	double			snr(1);							// Signal-to-noise ratio
	double			foreground_radius(0);			// Foreground radius
	int 			rand_type(0);					// Default uniform distribution
	int 			dist_type(0);					// Default uniform distribution
	
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
		if ( curropt->tag == "sampling" ) {
			sampling = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'u' ) rand_type = 0;
			else if ( curropt->value[0] == 'g' ) rand_type = 1;
			else if ( curropt->value[0] == 'p' ) rand_type = 2;
			else if ( curropt->value[0] == 'l' ) rand_type = 3;
			else if ( curropt->value[0] == 's' ) rand_type = 4;
			else
				cerr << "-type: Type " << curropt->value << " not supported! Default to uniform distribution" << endl;
		}
		if ( curropt->tag == "distance" ) {
			if ( curropt->value[0] == 'u' ) dist_type = 1;
			else if ( curropt->value[0] == 'g' ) dist_type = 2;
			else
				cerr << "-distance: Distance type " << curropt->value << " not supported! Default to uniform distribution" << endl;
		}
		if ( curropt->tag == "minmax" ) {
			if ( curropt->values(min, max) < 2 )
				cerr << "-minmax: Both min and max must be specified!" << endl;
			if ( min > max ) swap(min, max);
			if ( min == max ) {
				min -= 5;
				max += 5;
			}
		}
		if ( curropt->tag == "avgstd" ) {
			if ( curropt->values(avg, std) < 1 )
				cerr << "-avgstd: At least the average must be specified!" << endl;
			if ( std < 0 ) std = -std;
			if ( std == 0 ) std = 1;
			set_avgstd = 1;
		}
		if ( curropt->tag == "alpha" )
        	if ( ( alpha = curropt->value.real() ) < 0 )
				cerr << "-alpha: A positive spectral decay must be specified!" << endl;
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(finalavg, finalstd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "images" )
			if ( ( nimg = curropt->value.integer() ) < 1 )
				cerr << "-images: The number of images must be specified!" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "snr" ) {
			if ( curropt->values(snr, foreground_radius) < 1 )
				cerr << "-snr: A value must be specified!" << endl;
			else
				if ( snr <= 0 ) snr = 1;
		}
    }
	option_kill(option);
	
	if ( verbose & VERB_FULL )
		cout << "Random number range:            " << get_rand_max() << endl;
	
	double		ti = timer_start();
	
	double		scale, fg = 1e30;
	Bimage* 	p = NULL;
	Bimage* 	pran = NULL;
	
	if ( optind < argc - 1 ) {
		p = read_img(argv[optind++], 1, -1);
		p->change_type(Float);
		if ( rand_type > 0 ) {
			if ( set_avgstd ) p->rescale_to_avg_std(avg, std);
		} else
			p->rescale_to_min_max(min, max);
	}
		
	if ( nudatatype == Unknown_Type ) nudatatype = Float;
	
	if ( p ) {
		p->change_type(Float);
		nimg = p->images();
		size = p->size();
		if ( !set_avgstd ) {
			avg = p->average();
			std = p->standard_deviation();
		}
	} else {
		if ( size[0] < 1 ) size[0] = 256;
		if ( size[1] < 1 ) size[1] = 256;
		if ( size[2] < 1 ) size[2] = 1;
	}

	pran = new Bimage(Float, TSimple, size, nimg);

	if ( set_sampling ) pran->sampling(sampling);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) pran->origin(pran->size()/2);
		else pran->origin(origin);
	}
	
	switch ( rand_type ) {
		case 1:
			pran->noise_gaussian(avg, std);
			break;
		case 2:
			pran->noise_poisson(avg);
			break;
		case 3:
			pran->noise_logistical(avg, std);
			break;
		case 4:
			pran->noise_spectral(alpha);
			break;
		default:
			if ( dist_type < 1 ) pran->noise_uniform(min, max);
			break;
	}
	
	switch ( dist_type ) {
		case 1:
			pran->noise_uniform_distance(avg*pran->size().volume());
			break;
		case 2:
			pran->noise_gaussian_distance(avg*pran->size().volume(), std);
			break;
		default:
			break;
	}
	
	if ( p ) {
		if ( foreground_radius ) {
			if ( p->sizeZ() < 2 ) fg = M_PI*foreground_radius*foreground_radius;
			else fg = (M_PI*4.0/3.0)*foreground_radius*foreground_radius*foreground_radius;
		}
		if ( fg > size.volume() ) fg = size.volume();
		scale = sqrt(size.volume()/(snr*fg)) * p->standard_deviation()/pran->standard_deviation();
		p->add(pran, scale, -avg);
		delete pran;
		pran = p;
		pran->statistics();
		if ( finalstd ) pran->rescale_to_avg_std(finalavg, finalstd);
	}
	
	pran->change_type(nudatatype);
	
	if ( optind < argc )
		write_img(argv[optind], pran, 0);
	else
		write_img("temp.pif", pran, 0);
	
	delete pran;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

