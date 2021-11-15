/**
@file	bfft.cpp
@brief	General FFT for n-dimensional data
@author Bernard Heymann
@date	Created: 19980805
@date 	Modified: 20210817
Implementing the FFTW library
**/

#include "rwimg.h"
#include "qsort_functions.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototype
int 		img_fft_times(int ndim, int minsize, int maxsize, int opt);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bfft [options] file.hkl file.mrc",
"---------------------------------------",
"Fast Fourier transforms to and from image and reflection formats.",
" ",
"Actions:",
"-inverse                 Backward transform (default forward).",
"-phaseshift              Shift phase by half of the size.",
"-powerspectrum           Calculate powerspectrum estimate (can be used with the -tile option).",
"-average                 Average multiple power spectra (can be used with the -tile option).",
"-logarithm               Calculate the logarithm of the power spectrum.",
"-color 1.5               Generate a phase-coloured power spectrum: amplitude scaling.",
//"-weigh                   Weigh reflections with FOM.",
//"-cutoff 0.2              Remove reflections with FOM < cutoff threshold.",
"-zeroorigin              Zero the transform origin.",
"-halfshift               Shift correlation map by half the size to put the origin in the middle.",
//"-cone 45                 Remove missing cone.",
"-test 2,75,278,1         Test FFT execution time: dimensions,size min and max, and optimization flag.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-resolution 4.5,130      Resolution (default 0 - 10000 angstrom).",
"-free 15                 % flagged for R-free calculations (default 10).",
"-origin 0,-10,30         Set origin, used with -size option (default 0,0,0).",
"-size 22,33,14           Size for transform or hkl input (voxels, default from data).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-unitcell 100,100,100,90,90,90 Unit cell parameters.",
"-symmetry 1              Space group.",
"-tile 1024,1024,1        Size of tiles for calculating a power spectrum.",
" ",
//"Output:",
//"-Postscript out.ps       Postscript output for radial logarithm of power spectrum (only with -logarithm option).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
    fft_direction	setdir(FFTW_FORWARD);		// Forward transform
	int				phase_shift(0);				// Flag to shift phase by half the size
	int				setpower(0);				// Flag to calculate a power spectrum
	int				power_flags(0);				// Power spectrum flags
	int				setfom(0);					// Flag to modify reflections and FOM
    
    Vector3<long> 	size;						// Size for extraction or HKL
	Vector3<double>	nuorigin;					// Origin for extraction or HKL
	int				set_origin(0);				// Flag for setting the origin
	Vector3<double>	sam;    					// Sampling
	int 			spacegroup(0);  	    	// Spacegroup
	UnitCell		uc(0,0,0,M_PI_2,M_PI_2,M_PI_2);		// Unit cell parameters
	Vector3<long>	tile_size;					// Size of power spectrum tiles
    
    double			hires(0), lores(1e4);		// Limiting resolution range (hires must be > 0 to be set)
	int 			zero_origin(0);				// Don't zero the transform origin
	int				halfshift(0);				// Flag to shift origin to middle
    double			cone(M_PI_2);    	    	// Missing cone (default none)
    int 			weigh(0);   	    	    // No weighting
    double			fomcut(0);   	    	    // Selection based on FOM cutoff
    double			pfree(10);   	    		// % Rfree flags
	
	double	 		colour_phase_scale(0);		// Scale for coloured phases
	
	int				test_dim(0);				// Number of dimensions to test execution times
	int				test_minsize(0);			// Minimum size to test execution times
	int				test_maxsize(0);			// Maximum size to test execution times
	int				test_opt(0);				// Optimization flag for testing execution times
	
	Bstring			ps_file;					// Postscript file for radial power spectrum

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "inverse" ) setdir = FFTW_BACKWARD;
		if ( curropt->tag == "phaseshift" ) phase_shift = 1;
		if ( curropt->tag == "powerspectrum" ) setpower = 1;
		if ( curropt->tag == "average" ) power_flags |= 2;
		if ( curropt->tag == "logarithm" ) { power_flags |= 8; }
		if ( curropt->tag == "halfshift" ) { power_flags |= 4; halfshift = 1; }
		if ( curropt->tag == "zeroorigin" ) { power_flags |= 1; zero_origin = 1; }
 		if ( curropt->tag == "color" )
			if ( ( colour_phase_scale = curropt->value.real() ) < 0.001 )
		    		cerr << "-color: A scale for the phase colours must be specified" << endl;
 		if ( curropt->tag == "test" )
        	if ( curropt->values(test_dim, test_minsize, test_maxsize, test_opt) < 2 )
				cerr << "-test: At least the number of dimensions and one size must be specified" << endl;
 		if ( curropt->tag == "cone" ) {
 			if ( ( cone = curropt->value.real() ) < 1 )
				cerr << "-cone: A cone angle must be specified" << endl;
			else
				setfom = 1;
			if ( cone < 0 ) cone = -cone;
			if ( cone > TWOPI ) cone *= M_PI/180;
		}
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
			else
				setfom = 1;
  		}
		if ( curropt->tag == "weigh" )
       		setfom = weigh = 1;
  		if ( curropt->tag == "cutoff" ) {
       		if ( ( fomcut = curropt->value.real() ) < 0.001 )
				cerr << "-cutoff: FOM cutoff must be specified" << endl;
			else
				setfom = 1;
			if ( fomcut < 0 ) fomcut = 0;
		}
  		if ( curropt->tag == "free" )
       		if ( ( pfree = curropt->value.real() ) < 0.001 )
				cerr << "-free: %% reflections flagged must be specified" << endl;
 		if ( curropt->tag == "origin" ) {
        	nuorigin = curropt->origin();
			set_origin = 1;
		}
 		if ( curropt->tag == "size" )
			size = curropt->size();
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
 		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
 		if ( curropt->tag == "symmetry" )
        	if ( ( spacegroup = curropt->value.integer() ) < 1 )
				cerr << "-symmetry: The space group number must be specified" << endl;
		if ( curropt->tag == "tile" )
			tile_size = curropt->size();
 		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	
	if ( test_dim && test_minsize )
		img_fft_times(test_dim, test_minsize, test_maxsize, test_opt);

	// Allocate memory for the structure factors and image parameters
    Bimage*	 		p = NULL;
	
    // Read the file
	if ( optind >= argc-1 ) {
		bexit(0);
	} else if ( ( p = read_img(argv[optind++], 1, -1) ) == NULL ) {
		cerr << "File not read!" << endl;
		bexit(-1);
	}
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	
	if ( p->data_type() < Float )
		p->change_type(Float);
	
	if ( spacegroup ) p->space_group(spacegroup);
	
	if ( sam.volume() > 0 ) p->sampling(sam);
		
	if ( set_origin ) p->origin(nuorigin);
		
	if ( uc.check() ) p->unit_cell(uc);
		
    // Apply symmetery and select specific reflections if it is a transform
	if ( p->compound_type() == TComplex ) {
		if ( setfom ) {
			if ( zero_origin ) p->zero_fourier_origin();
		}
	}
	
	if ( size.volume() > 0 )
		if ( p->fourier_type() > NoTransform )
			p->change_transform_size(size);
	
	Bimage*		pnu = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG bfft: power_flags=" << power_flags << endl;
	
	if ( setpower ) {
		if ( tile_size.volume() > 0 ) {
			power_flags |= 2;
			pnu = p->powerspectrum_tiled(0, tile_size, power_flags);
			if ( pnu ) {
				delete p;
				p = pnu;
			}
		} else {
			p->power_spectrum(power_flags);
		}
		p->information();
		p->powerspectrum_isotropy(0, lores, hires);
	} else {
		if ( phase_shift && p->fourier_type() == Standard )
			p->phase_shift_to_center();
		if ( p->fft(setdir) )
			return error_show("bfft", __FILE__, __LINE__);
		if ( phase_shift && p->fourier_type() == Standard )
			p->phase_shift_to_center();
	}
	
	if ( size.volume() > 0 )
		if ( p->fourier_type() > NoTransform )
			p->change_transform_size(size);
	
	if ( p->compound_type() == TComplex ) {
		if ( setfom ) {
			if ( zero_origin ) p->zero_fourier_origin();
		}
	
		if ( colour_phase_scale > 0 ) {
			pnu = p->intensities_phase_colored(colour_phase_scale);
			if ( pnu ) {
				delete p;
				p = pnu;
				nudatatype = p->data_type();
			}
			Bstring		filename("ps_color_wheel.png");
			pnu = new Bimage(UCharacter, TRGB, 512, 512, 1, 1);
			pnu->phase_colour_wheel();
			write_img(filename, pnu, 0);
			delete pnu;
		}
	}
	
	// Write an output file if a file name is given
	if ( optind < argc ) {
		p->change_type(nudatatype);
		if ( halfshift ) p->center_wrap();
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

/**
@brief 	Tests the execution times for a series of image sizes.
@param 	ndim		number of dimensions (1,2,3).
@param 	minsize		minimum image size.
@param 	maxsize		maximum image size.
@param 	opt			optimization (with FFTW_MEASURE).
@return int 		0.

	FFTW library (www.fftw.org).
	Blank complex floating point images are created and transformed.
	Only the call to the complex FFT function is timed.

**/
int 		img_fft_times(int ndim, int minsize, int maxsize, int opt)
{
	long			i, size, nprime;
	long*			prime = NULL;
    Bimage*			p = NULL;
	double			tvs, tvp, tvf;
	double			dtp, dte;
	fft_plan		plan;
	
	if ( ndim < 1 ) ndim = 1;
	if ( ndim > 3 ) ndim = 3;
	if ( minsize < 1 ) minsize = 1;
	if ( maxsize < minsize ) maxsize = minsize;
	
	long			n = maxsize - minsize + 1;
	int_float*		t = new int_float[n];

	cout << endl << "Timing the execution of fast Fourier transforms:" << endl;
	cout << "FFTW planner option:            " << opt << endl;
	cout << "Number of dimensions:           " << ndim << endl;
	cout << "Size range:                     " << minsize << " - " << maxsize << endl << endl;
	cout << "Size\tPtime\tEtime\tPrime factors" << endl;
	for ( size = minsize; size <= maxsize; size++ ) {
		switch ( ndim ) {
			case 1:
				p = new Bimage(Float, TComplex, size, 1, 1, 1);
				break;
			case 2:
				p = new Bimage(Float, TComplex, size, size, 1, 1);
				break;
			case 3:
				p = new Bimage(Float, TComplex, size, size, size, 1);
				break;
		}
		prime = prime_factors(size, nprime);
		tvs = getwalltime();
		plan = p->fft_setup(FFTW_FORWARD, opt);
		tvp = getwalltime();
		if ( p->fft(plan, 0) )
			return error_show("img_fft_times", __FILE__, __LINE__);
		tvf = getwalltime();
		fft_destroy_plan(plan);
		dtp = tvp - tvs;
		dte = tvf - tvp;
		cout << size << tab << fixed << setprecision(5) << dtp << tab << dte;
		for ( i=0; i<nprime; i++ ) cout << tab << prime[i];
		cout << endl;
		delete[] prime;
		delete p;
		t[size-minsize].i = size;
		t[size-minsize].f = dte;
	}
	cout << endl;
	
	qsort((void *) t, n, sizeof(int_float), QsortSmallToLargeIntFloat);
	
	cout << "Ordered times:" << endl;
	cout << "Size\tEtime" << endl;
	for ( i=0, size=minsize; i<n; i++, size++ ) cout << t[i].i << tab << t[i].f << endl;
	cout << endl;

	delete[] t;

	return 0;
}
