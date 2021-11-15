/**
@file	bnad.cpp
@brief	Image denoising by nonlinear anisotropic diffusion: Coherence and edge enhancing diffusion
@author Achilleas Frangakis
@author Bernard Heymann
@date	Created: 20020803
@date	Modified: 20161210 (BH)
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
"Usage: bnad [options] input.img output.img",
"------------------------------------------",
"Denoises by nonlinear anisotropic diffusion: Coherence and edge enhancing diffusion.",
"The default is edge enhancing diffusion.",
"Coherence enhancing diffusion is used if the coherence parameter is specified.",
"Diffusion is faster for higher lambda (1e-6 to 10 have been tested).",
" ",
"Actions:",
"-lambda 0.3              Lambda: edge enhancing diffusion (default 0.1).",
"-coherence 0.002         Coherence: coherence enhancing diffusion.",
"-sigma 2.5               Initial gaussian smoothing coefficient (default no smoothing).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-iterations 80           Iterations (default 1).",
"-slabsize 20             Slab size for piece-wise processing (default z size).",
"-output 5                Number of iterations between writing output maps (default 100).",
" ",
"Parameters for coherence enhancing diffusion:",
"-alpha 0.001             Alpha (default 0.001).",
"-rho 4                   Gaussian smoothing of structure tensor (default 6).",
" ",
NULL
};

int main (int argc, char **argv)
{
	DataType        nudatatype(Unknown_Type);     // Conversion to new type
	long   			i;             			/* loop variables */
	long   			imax(1);               	/* largest iteration number */
	long			slabsize(0);			// Slab size for threaded runs
	long   			iout(100);             	/* iterations between writing output maps */
	double	  		ht(0.1);               	/* time step size */
	double	  		sigma(0);              	/* noise scale - Presmoothing */
	double	  		lambda(0.1);           	/* lamda Parameter */
	double			alpha(0.001);          	/* linear diffusion fraction */
	double			rho(6);                	/* integration scale */
	double	  		C(0);                  	/* coherence parameter */

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "lambda" )
			if ( ( lambda = curropt->value.real() ) <= 0 )
				cerr << "-lambda: Lambda must be specified" << endl;
		if ( curropt->tag == "coherence" )
        	if ( ( C = curropt->value.real() ) <= 0 )
				cerr << "-coherence: Coherence must be specified" << endl;
		if ( curropt->tag == "sigma" )
        	if ( ( sigma = curropt->value.real() ) <= 0 )
				cerr << "-sigma: A value must be specified" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "iterations" )
        	if ( ( imax = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "slabsize" )
        	if ( ( slabsize = curropt->value.integer() ) < 3 )
				cerr << "-slabsize: A slab size must be specified!" << endl;
		if ( curropt->tag == "output" )
        	if ( ( iout = curropt->value.integer() ) < 1 )
				cerr << "-output: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "alpha" )
			if ( ( alpha = curropt->value.real() ) <= 0 )
				cerr << "-alpha: Alpha must be specified" << endl;
		if ( curropt->tag == "rho" )
			if ( ( rho = curropt->value.real() ) < 1 )
				cerr << "-rho: A value must be specified" << endl;
    }
	option_kill(option);
	
	double		ti = timer_start();

	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	
	p->change_type(Float);
	
	if ( slabsize < 1 ) slabsize = p->sizeZ();
	
	long		memreq(p->sizeX()*p->sizeY()*sizeof(float));
	if ( p->sizeZ() == 1 )
		memreq *= 3;
	else {
		memreq *= 6*slabsize;
		memreq *= system_processors();
	}
	memreq += 2*p->data_size()*sizeof(float);

	if ( verbose ) {
		cout << "Denoising by non-linear anisotropic diffusion: ";
		if ( C ) cout << "Coherence enhancing diffusion" << endl;
		else cout << "Edge enhancing diffusion" << endl;
		cout << "Image size:                     " << p->size() << endl;
		if ( p->sizeZ() == 1 )
			cout << "Images:                         " << p->images() << endl;
		else
			cout << "Slab thickness:                 " << slabsize << endl;
	}
	
	memory_check(memreq);
	
	if ( sigma > 0 ) {
		if ( verbose )
			cout << "Gaussian filter sigma:          " << sigma << endl;
		p->filter_gaussian((long) (6*sigma), sigma);
	}
	
	Bstring			filename;
	Bimage*			pd = NULL;

	time_t			t = time(NULL);
	
	double			pstd;
	
	if ( verbose ) {
		cout << "Iter\tAvg\tStd\t∆Std\tTime" << endl;
		cout << "0\t" << p->average() << tab << p->standard_deviation() << tab
			<< "0\t" << (unsigned long)(time(NULL) - t) << endl << flush;
	}
	
	for ( i=1; i<=imax; i++ ) {
		pstd = p->standard_deviation();
		if ( p->sizeZ() == 1 )
			pd = p->nad_2D(ht, lambda, C, alpha);
		else
			pd = p->nad(ht, slabsize, lambda, C, alpha);
		p->replace(pd);
		p->statistics();
		
		/* check minimum, maximum, mean, variance */
		if ( verbose )
			cout << i << tab << setprecision(4) << p->average() << tab <<
				p->standard_deviation() << tab <<
				pstd - p->standard_deviation() << tab <<
				(size_t)(time(NULL) - t) << endl << flush;
				
		if ( verbose & VERB_PROCESS )
			p->information();
		
		if ( i%iout == 0 && i != imax ) {
			filename = argv[optind];
			filename = filename.pre_rev('.') + Bstring(i, "_%04d.") + filename.post_rev('.');
			pd->change_type(nudatatype);
			write_img(filename, pd, 0);
		}
		
		delete pd;

		if ( verbose & VERB_TIME )
			timer_report(ti);
	
	}

	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

