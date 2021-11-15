/**
@file	fft.cpp
@brief	General FFT for n-dimensional data
@author Bernard Heymann
@date	Created: 19980805
@date	Modified: 20151002

		Implementing the FFTW library
**/

#include "fft.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Sets up a plan for fast Fourier transforms.
@param 	x			x dimension.
@param 	y			y dimension.
@param 	z			z dimension.
@param 	dir			direction of transformation (FFTW_FORWARD or FFTW_BACKWARD)
@param 	opt			optimization (0=FFTW_ESTIMATE, 1=FFTW_MEASURE, 2=FFTW_PATIENT, 3=FFTW_EXHAUSTIVE).
@return fft_plan 	FFTW plan.

	FFTW library (www.fftw.org).
	The size and direction determines the plan.
	Both FFTW versions 2 and 3 are supported.

**/
fft_plan	fft_setup_plan(long x, long y, long z, fft_direction dir, int opt)
{
	int					rank = 0;
	int					n[3] = {1, 1, 1};
	if ( z > 1 ) {
		rank = 3;
		n[0] = z;
		n[1] = y;
		n[2] = x;
	} else if ( y > 1 ) {
		rank = 2;
		n[0] = y;
		n[1] = x;
	} else {
		rank = 1;
		n[0] = x;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG fft_setup_plan: n=" << n[0] << "x" << n[1] << "x" << n[2] << " opt=" << opt << endl;

	int					flags = FFTW_ESTIMATE;
	
	switch ( opt ) {
		case 1: flags = FFTW_MEASURE; break;
		case 2: flags = FFTW_PATIENT; break;
		case 3: flags = FFTW_EXHAUSTIVE; break;
		default: flags = FFTW_ESTIMATE;
	}
	
	fft_complex*		in = NULL;
	fft_complex*		out = in;
	if ( opt )
		in = out = new fft_complex[x*y*z];
	
	fft_plan			plan = fftwf_plan_dft(rank, n, in, out, dir, flags);

	if ( opt )
		delete[] in;

	return plan;
}

fft_plan	fft_setup_plan(Vector3<long> size, fft_direction dir, int opt)
{
	return fft_setup_plan(size[0], size[1], size[2], dir, opt);
}

/**
@brief 	Deallocates a plan for fast Fourier transforms.
@param 	plan		FFTW plan.
@return int			0.

	FFTW library (www.fftw.org).

**/
int			fft_destroy_plan(fft_plan plan)
{
	fftwf_destroy_plan(plan);

	return 0;
}

/**
@brief 	Fast Fourier transforms a data array.
@param 	plan		Fourier transform plan.
@param 	*a			data array.
@return int 		error code.

	FFTW library (www.fftw.org).
	The transform plan encodes both the direction and size of the transform.
	The transformation is done in place and the resultant data are 
	returned within the original array.

**/
int			fftw(fft_plan plan, Complex<float>* a)
{
	fft_complex*	in = (fft_complex *) a;
	fft_complex*	out = in;
	
	fftwf_execute_dft(plan, in, out);
	
	return 0;
}


