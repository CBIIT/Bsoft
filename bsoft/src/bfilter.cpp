/**
@file	bfilter.cpp
@brief	Program to filter images.
@author Bernard Heymann
@date	Created: 20010414
@date	Modified: 20210127
**/

#include "rwimg.h"
#include "rwkernel.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bfilter [options] input.img output.img",
"---------------------------------------------",
"Filters images.",
" ",
"Actions:",
"-invert                  Invert density in the image.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation after filtering.",
"-minmax -0.5,1.2         Rescale data to minimum and maximum.",
"-average 7,5,3           Averaging/smoothing filter: kernel size.",
"-gaussian 11,2.6         Gaussian smoothing filter: kernel size and sigma.",
"-median 3                Median filter: kernel edge size.",
"-difference file.mrc     Filter by difference with the given image.",
"-peak 5                  Peak filter: kernel edge size.",
"-gradient                Gradient filter (3x3x3).",
"-laplacian               Laplacian filter (3x3x3).",
"-LoG 2.5                 Laplacian-of-Gaussian filter.",
"-DoG 1.5,3.1             Difference-of-Gaussians filter.",
"-denoise 2,0.4           Denoising filter: distance and density difference sigmas.",
"-rollingball 5,4.3       Rolling ball filter: radius and density scaling (default scaling=1).",
"-variance 11,1           Calculate a local variance image using the given size kernel,",
"                         with a flag to indicate standard deviation rather than variance.",
"-extremes his            Filter image extremes: his=histogram-based, mg=for micrograph.",
"-replacemaxima 1.5       Replace maxima above the given threshold with surrounding average.",
"-bandpass 25.3,200,0.02  Bandpass filter: resolution limits (angstrom) and band edge width (1/angstrom).",
"-frequency 25,30         Frequency filter: frequency (Å or 1/Å) and gaussian width (Å).",
"-gabor 25,35,45,30,50    Gabor filter: frequency location (Å or 1/Å) and gaussian widths (Å).",
"-anisotropic 3.4,6.7,8,2 Anisotropic Gaussian filter: three sigma values and gradient direction.",
"-amplitude 8.3           Amplitude filter: setting all amplitudes below the threshold to zero.",
"-phasesonly              Set all amplitudes to one.",
"-Bfactor 44,25.3         B-factor application: B-factor (A^2), high resolution limit (A,optional)",
"                         Multiplied in reciprocal space by exp(-B-factor/4 * s^2).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 10,-10,20        Origin for rotation (default 0,0,0).",
"-wrap                    Turn wrapping on (default off, use with -denoise option).",
" ",
"Input:",
"-kernel file.txt         Convolve image with a kernel specified in a text file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// Origin for rotation
	int				set_origin(0);				// Flag to set origin
	int 			setinvert(0);				// Flag to invert density
	double			nuavg(0), nustd(0); 		// Values for rescaling
	int				setminmax(0);				// Flag to rescale to min and max
	double			newmin(0), newmax(0); 		// Rescaling to new min and max
	Vector3<long>	average_kernel;				// Average filter kernel size
	long			gauss_kernel(0);			// Gaussian kernel size
	double			sigma(0);					// Gaussian sigma
	long 			median_kernel(0);			// Median filter kernel size
	Bstring			dif_file;					// File to calculate difference filter
	long 			peak_kernel(0);				// Peak filter kernel size
	int				gradient(0);				// Gradient filter
	int				laplacian(0);				// Laplacian filter
	double			log_sigma(0);				// Laplacian-of-Gaussian filter sigma
	int				dog1(0), dog2(0);			// Difference-of-Gaussians filter
	double			sigma1(0), sigma2(0);		// Denoising filter sigmas
	Vector3<double>	aniso;						// Anisotropic Gaussian filter
	int				anidir(0);					// Anisotropic gradient direction
	long			ball_radius(0);				// Rolling ball radius
	double			ball_scale(0);				// Denisty scale for rolling ball
//	int				wrap(0);					// Flag to wrap around
	long 			var_kernel(0);				// Variance filter kernel size
	int 			std_flag(0);				// Flag to claculate StDev rather than variance
	Bstring			kernel_file;				// Convolution kernel file
	int				dofft(0);					// Flag to do FFT
	double			res_hi(0);					// High resolution limit
	double			res_lo(0); 					// Low resolution limit
	double			Bfactor(0); 				// B-factor
	int 			bandpass(0);				// Fourier bandpass filter
	double			bandpass_width(0);			// Edge width for the bandpass filter
	double			frequency(0);				// Frequency filter
	Vector3<double>	freqloc;					// Gabor filter frequency location
	double			fsigma(0), psigma(0);		// Frequency sigmas
	double 			amp_threshold(0);			// Fourier amplitude filter
	int				phases_only(0);				// Flag to set all amplitudes to one
	int 			filter_extremes(0);			// Filter extremes, 1=histogram, 2=micrograph, 3=neighbor average
	double			exc_min(0), exc_max(0);		// Extreme limits
	double			replace_threshold(0);		// Threshold for replacing maxima
	int				replace_flag(0);			// Flag for replacing maxima
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
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
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "minmax" ) {
			if ( curropt->values(newmin, newmax) < 2 )
				cerr << "-minmax: Both min and max must be specified!" << endl;
			else
				setminmax = 1;
			if ( newmin == newmax ) setminmax = 0;
		}
		if ( curropt->tag == "average" ) {
			if ( ( i = curropt->values(average_kernel[0],
					average_kernel[1], average_kernel[2]) ) < 1 )
				cerr << "-average: The kernel edge size must be specified!" << endl;
			else
				if ( i < 2 ) average_kernel[2] = average_kernel[1] = average_kernel[0];
		}
		if ( curropt->tag == "gaussian" )
			if ( curropt->values(gauss_kernel, sigma) < 2 )
				cerr << "-gaussian: The kernel edge size and sigma must be specified!" << endl;
		if ( curropt->tag == "median" )
			if ( ( median_kernel = curropt->value.integer() ) < 1 )
				cerr << "-median: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "difference" )
			dif_file = curropt->filename();
		if ( curropt->tag == "peak" )
			if ( ( peak_kernel = curropt->value.integer() ) < 1 )
				cerr << "-peak: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "gradient" ) gradient = 1;
		if ( curropt->tag == "laplacian" ) laplacian = 1;
		if ( curropt->tag == "LoG" ) {
			if ( ( log_sigma = curropt->real() ) < 0.1 )
				cerr << "-LoG: A sigma value must be specified!" << endl;
			else
				dofft = 1;
		}
		if ( curropt->tag == "DoG" )
			if ( curropt->values(dog1, dog2) < 2 )
				cerr << "-DoG: Two sigma values must be specified!" << endl;
		if ( curropt->tag == "denoise" )
			if ( curropt->values(sigma1, sigma2) < 2 )
				cerr << "-denoise: Both denoising sigmas must be specified!" << endl;
		if ( curropt->tag == "rollingball" )
			if ( curropt->values(ball_radius, ball_scale) < 1 )
				cerr << "-rollingball: At least a ball radius must be specified!" << endl;
		if ( curropt->tag == "variance" )
			if ( curropt->values(var_kernel, std_flag) < 1 )
				cerr << "-variance: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "kernel" )
			kernel_file = curropt->filename();
		if ( curropt->tag == "extremes" ) {
			if ( curropt->value.contains("his") ) filter_extremes = 1;
			else if ( curropt->value.contains("mg") ) filter_extremes = 2;
			else if ( curropt->values(exc_min, exc_max) > 1 ) filter_extremes = 3;
			else cerr << "-extremes: The type of filtering must be specified!" << endl;
        }
		if ( curropt->tag == "replacemaxima" ) {
			if ( ( replace_threshold = curropt->value.real() ) < 0.0000001 )
				cerr << "-replacemaxima: A threshold must be specified!" << endl;
			else
				replace_flag = 1;
		}
		if ( curropt->tag == "bandpass" ) {
			if ( curropt->values(res_hi, res_lo, bandpass_width) < 1 )
				cerr << "-bandpass: A resolution limit must be specified!" << endl;
			else {
				if ( res_lo > 0 && res_hi > res_lo ) swap(res_hi, res_lo);
				if ( res_hi != res_lo ) bandpass = dofft = 1;
			}
		}
		if ( curropt->tag == "frequency" ) {
			if ( curropt->values(frequency, fsigma) < 1 )
				cerr << "-frequency: A frequency must be specified!" << endl;
			else
				dofft = 1;
		}
		if ( curropt->tag == "gabor" ) {
			if ( curropt->values(freqloc[0], freqloc[1], freqloc[2], fsigma, psigma) < 3 )
				cerr << "-gabor: A frequency location must be specified!" << endl;
			else
				dofft = 1;
		}
		if ( curropt->tag == "anisotropic" ) {
			if ( ( i = curropt->values(aniso[0], aniso[1], aniso[2], anidir) ) < 3 )
				cerr << "-anisotropic: Three values must be specified!" << endl;
			else {
				if ( i < 3 ) aniso[2] = 1;
				else if ( i < 2 ) aniso[1] = aniso[2] = aniso[0];
				dofft = 1;
			}
		}
		if ( curropt->tag == "amplitude" ) {
			if ( ( amp_threshold = curropt->value.real() ) < 0.0000001 )
				cerr << "-amplitude: A threshold must be specified!" << endl;
			else
				dofft = 1;
		}
		if ( curropt->tag == "phasesonly" )
			phases_only = dofft = 1;
		if ( curropt->tag == "Bfactor" ) {
			if ( curropt->values(Bfactor, res_hi) < 1 )
				cerr << "-Bfactor: A B-factor must be specified!" << endl;
			else
				dofft = 1;
		}
//		if ( curropt->tag == "wrap" ) wrap = 1;
    }
	option_kill(option);
	
	double		ti = timer_start();

#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	
	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();			// Preserve old type
	
	if ( sam.volume() ) p->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}

	if ( filter_extremes == 1 ) p->filter_extremes(1);
	
	if ( filter_extremes == 2 ) p->filter_extremes();
	
	if ( filter_extremes == 3 ) p->filter_extremes(exc_min, exc_max);

	if ( replace_flag ) p->replace_maxima(replace_threshold);
	
	if ( average_kernel[0] > 0 ) p->filter_average(average_kernel);
	
	if ( sigma > 0 ) p->filter_gaussian(gauss_kernel, sigma);
	
	if ( median_kernel > 0 ) p->filter_rank(median_kernel, 0.5);
	
	if ( dif_file.length() ) {
		Bimage*		p2 = read_img(dif_file, 1, -1);
		p->filter_by_difference(p2);
		delete p2;
	}
	
	if ( peak_kernel > 0 ) {
		Bimage*		pp = p->filter_peak(peak_kernel);
		delete p;
		p = pp;
	}

	Bimage*		pkernel = NULL;
	
	if ( var_kernel > 0 ) {
		if ( std_flag & 2 ) {
			Vector3<long> 	ksize(var_kernel, var_kernel, var_kernel);
			ksize = ksize.min(p->size());
			pkernel = new Bimage(Float, TSimple, ksize, 1);
			pkernel->origin(ksize/2);
			pkernel->sphere(ksize/2, var_kernel/2, 1, FILL_USER, 1);
//			write_img("pk.mrc", pkernel);
			p->variance(pkernel);
			delete pkernel;
		} else {
			p->variance(var_kernel);
		}
		if ( std_flag & 1 ) p->square_root();
	}

	if ( kernel_file.length() ) {
		pkernel = read_img(kernel_file, 1, 0);
		p->convolve(pkernel);
//		write_img("t.krn", pkernel);
		delete pkernel;
	}

	if ( gradient ) p->filter_ortho(0);
	if ( laplacian ) p->filter_ortho(1);
	
	if ( dog1 >0 && dog2 > 0 )
		p->filter_dog(dog1, dog2);
		
	if ( sigma1 > 0 && sigma2 > 0 )
		p->filter_bilateral(sigma1, sigma2, 1, 0);
	
	if ( ball_radius )
		p->filter_rolling_ball(ball_radius, ball_scale);
	
	if ( dofft ) {
		FourierType		old_fourier_type = p->fourier_type();
		if ( p->fourier_type() == NoTransform ) p->fft();
		
		if ( amp_threshold ) p->fspace_amp_threshold(amp_threshold);
	
		if ( phases_only ) p->fspace_amp_one();
		
		if ( Bfactor != 0 ) p->fspace_weigh_B_factor(Bfactor, res_hi);
	
		if ( log_sigma > 0 ) p->fspace_weigh_LoG(res_hi, log_sigma);
	
		if ( bandpass ) {
			if ( bandpass_width < 2 ) p->fspace_bandpass(res_hi, res_lo, bandpass_width);
			else p->fspace_butterworth_band(res_hi, res_lo, bandpass_width);
		}
		
		if ( frequency ) p->fspace_frequency_filter(frequency, fsigma);
	
		if ( freqloc.length() ) p->fspace_gabor_filter(freqloc, fsigma, psigma);
	
		if ( aniso[0] ) p->fspace_weigh_gaussian(0, aniso, anidir);
		
		if ( old_fourier_type == NoTransform ) p->fft_back();
	}
	
	if ( setinvert ) p->invert();
	
	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);

	if ( setminmax ) p->rescale_to_min_max(newmin, newmax);
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;

#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

