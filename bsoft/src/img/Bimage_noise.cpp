/**
@file	Bimage_noise.cpp
@brief	Functions for generating noise images
@author Bernard Heymann
@date	Created: 19990703
@date	Modified: 20150725
**/

#include "Bimage.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates an image with a uniform random distribution of densities.
@param 	rmin 	minimum density value.
@param 	rmax 	maximum density value.
@return int		0.

	An image with a given datatype, number of subimages, and size is
	generated with densities distributed uniformly in the range of the
	given minimum and maximum:
		density = random_value*(max - min) + min
	where random_value is between 0 and 1.
	The average and standard deviation are:
		average = (max + min)/2
		standard deviation = 0.5*sqrt(1/3)*(max - min).
	The output image is floating point.
	Statistics are calculated before returning.

**/
int 		Bimage::noise_uniform(double rmin, double rmax)
{
	if ( rmin < dtmin ) rmin = dtmin;
	if ( rmax > dtmax ) rmax = dtmax;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a random image with a uniform distribution:" << endl;
		cout << "Minimum and maximum:            " << rmin << " " << rmax << endl << endl;
	}
		
	long			i;
	double			range = (rmax - rmin)/get_rand_max();
	
	random_seed();
	
	for ( i=0; i<datasize; i++ ) add(i, random()*range + rmin);
	
	statistics();

	return 0;
}

/**
@brief 	Generates an image with a gaussian random distribution of densities.
@param 	ravg 		average.
@param 	rstd 		standard deviation.
@return int		0.

	An image with a given datatype, number of subimages, and size is
	generated with densities with a gaussian distribution with a given
	average and standard deviation:
		density = average + std_dev*sqrt(-2*log(random_value))*
						cos(2*PI*random_value);
	where random_value is between 0 and 1.
	The output image is floating point.
	Statistics are calculated before returning.

**/
int			Bimage::noise_gaussian(double ravg, double rstd)
{
	if ( rstd <= 0 ) {
		error_show("Error in Bimage::noise_gaussian: The standard deviation for a Gaussian distribution must be > 0!", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a random image with a gaussian distribution:" << endl;
		cout << "Average and standard deviation: " << ravg << " " << rstd << endl << endl;
	}
		
	long			i;

	random_seed();
	
	for ( i=0; i<datasize; i++ )
		add(i, random_gaussian(ravg, rstd));

	statistics();
	
	return 0;
}

/**
@brief 	Generates an image with a poisson random distribution of densities.
@param 	ravg 	average.
@return int		0.

	Algorithm taken from Numerical Recipes in C.
	The poisson distribution is given for j = 0,1,... by:
		        avg^j * exp(-avg)
		P(j) = -----------------
		               j!
	Note that only positive integer values are defined for j and sum(P(j)) = 1.
	An array of floating point numbers is generated with a poisson 
	distribution with a given average. The standard deviation is:
		std = sqrt(avg)
	If the average <= 0, the function exits.
	Statistics are calculated before returning.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
int			Bimage::noise_poisson(double ravg)
{
	if ( ravg <= 0 ) {
		error_show("Error in Bimage::noise_poisson: The average for a poisson distribution must be > 0!", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a random image with a poisson distribution:" << endl;
		cout << "Average:                        " << ravg << endl << endl;
	}
		
	long			i;

	random_seed();
	
	for ( i=0; i<datasize; i++ )
		add(i, random_poisson(ravg));

	statistics();
	
	return 0;
}

/**
@brief 	Generates an image with a gaussian random distribution of densities.
@param 	ravg 	average.
@param 	rstd 	standard deviation.
@return int		0.

	An image with a given datatype, number of subimages, and size is
	generated with densities with a logistical differential distribution 
	with a given average and standard deviation:
		density = average + (std_dev/golden)*ln(1/random_value - 1)
	where random_value is between 0 and 1 and:
		golden  = (sqrt(5) + 1)/2
	The output image is floating point.
	Statistics are calculated before returning.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
int			Bimage::noise_logistical(double ravg, double rstd)
{
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a random image with a logistical differential distribution:" << endl;
		cout << "Average and standard deviation: " << ravg << " " << rstd << endl << endl;
	}
		
	long			i;

	random_seed();
	
	for ( i=0; i<datasize; i++ )
		add(i, random_logistical(ravg, rstd));

	statistics();
	
	return 0;
}

/**
@brief 	Generates a noise map with a defined spectral decay.
@param 	alpha	spectral decay.
@return int		0.

	Uniform random phases are generated and the amplitudes are set to:
		amp = s^(-alpha/2).

**/
int			Bimage::noise_spectral(double alpha)
{
	long			i, j, nn, xx, yy, zz, cc, imgsize(x*y*z);
	double			x2, y2, z2, s2, w, sum;
	double			fac(-alpha/4);	// Factor of 2 because we work with s2, and another 2 for amplitude scaling
	double			range(2.0/get_rand_max());
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a spectral noise image:" << endl;
		cout << "Alpha:                          " << alpha << endl << endl;
	}

	for ( i=nn=0; nn<n; nn++ ) {
		sum = 0;
		for ( zz=0; zz<z; zz++ ) {
			z2 = zz*1.0L/z;
			if ( z2 >= 0.5 ) z2 -= 1;
			z2 *= z2;
			for ( yy=0; yy<y; yy++ ) {
				y2 = yy*1.0L/y;
				if ( y2 >= 0.5 ) y2 -= 1;
				y2 *= y2;
				for ( xx=0; xx<x; xx++, i++ ) {
					x2 = xx*1.0L/x;
					if ( x2 >= 0.5 ) x2 -= 1;
					x2 *= x2;
					s2 = x2 + y2 + z2;
					if ( s2 ) {
						w = pow(s2, fac);
						for ( cc=0, j=i*c; cc<c; cc++, j++ )
							set(j, w*(random()*range-1));
						sum += w*w;
					}
				}
			}
		}
		sum = sqrt(sum/(3*imgsize));
		for ( j=nn*imgsize; j<(nn+1)*imgsize; j++ ) set(j, complex(j)/sum);
	}
	
	fft_back();
	
	statistics();
	
	return 0;
}
