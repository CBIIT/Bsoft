/**
@file	fft.h
@brief	General FFT for n-dimensional data
@author Bernard Heymann
@date	Created: 19980805
@date	Modified: 20151002
		Implementing the FFTW library
**/

#include "Complex.h"
#include "Vector3.h"

#include <fftw3.h>

// This must be updated together with the FFTW3 package
#define FFTW3_VERSION "3.3.6-pl2"

typedef fftwf_complex fft_complex;
typedef int fft_direction;
typedef fftwf_plan fft_plan;

#define MAX_RANK 	3	// Only up to 3D transforms

/* Function prototypes */
fft_plan	fft_setup_plan(long x, long y, long z, fft_direction dir, int opt);
fft_plan	fft_setup_plan(Vector3<long> size, fft_direction dir, int opt);
int			fft_destroy_plan(fft_plan plan);
int			fftw(fft_plan plan, Complex<float>* a);
