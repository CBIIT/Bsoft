/**
@file	random_numbers.cpp
@brief	Functions for generating random sequences
@author Bernard Heymann
@date	Created: 19990703
@date	Modified: 20151113
**/

#include "Vector3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
int			randseed(0);	// Flag to indicate if a random seed has been generated

/**
@brief 	Finds the maximum random number for a system.
@return	long 			the maximum random number.

	Loops through random numbers to determine if the maximum is 2 or 4 bytes.
**/
long 		get_rand_max()
{
	long 		i, rand_max(32767);
	
	for ( i=0; i<100; i++ )
		if ( random() > rand_max ) rand_max = 2147483647;
	
	return rand_max;
}

double		irm(1.0L/get_rand_max());

/**
@brief 	Gets a random seed and sets the flag.
@return	int 			flag.

	The random seed is obtained using the program pid.
	A flag is set to prevent the seed from being generated multiple times.
**/
long 		random_seed()
{
	if ( randseed ) return 1;
	
	srandom(getpid());
	
	randseed = 1;
	
	return 1;
}

int			random_array_uniform_chunk(float* r, long start, long end, double min, double range)
{
	for ( long i=start; i<end; i++ )
		r[i] = random()*range + min;
	
	return 0;
}

/**
@brief 	Generates a series with a uniform random distribution.
@param 	n				number of values.
@param 	min 			minimum value.
@param 	max 			maximum value.
@return float*				array with uniform random numbers.

	An array of floating point numbers is generated distributed uniformly 
	in the range of the given minimum and maximum:
		value = random_value*(max - min) + min
	where random_value is between 0 and 1.
	The average and standard deviation are:
		average = (max + min)/2
		standard deviation = 0.5*sqrt(1/3)*(max - min).

**/
float* 		random_array_uniform(long n, double min, double max)
{
	float* 			r = new float[n];
	 
	double			range = (max - min)/get_rand_max();
	
	if ( verbose & VERB_FULL ) {
		cout << "Generating a random series with a uniform distribution:" << endl;
		cout << "Size:                           " << n << endl;
		cout << "Minimum & maximum:              " << min << " " << max << endl;
	}
	
	random_seed();

	long			chunk_size = get_chunk_size(n);

#ifdef HAVE_GCD
	dispatch_apply((n - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		start = i*chunk_size;
		long		end = start + chunk_size;
		if ( end > n ) end = n;
		random_array_uniform_chunk(r, start, end, min, range);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<n; i+=chunk_size ) {
		long		end = i + chunk_size;
		random_array_uniform_chunk(r, i, end, min, range);
	}
#endif

	return r;
}

/**
@brief 	Generates a value with a gaussian random distribution.
@param 	avg 			average.
@param 	stdev 		standard deviation.
@return double				the random value.

	A floating point number is generated with a gaussian 
	distribution with a given average and standard deviation:
		value = average + std_dev*sqrt(-2*log(random_value))*
						cos(2*PI*random_value);
	where random_value is between 0 and 1.

**/
double 		random_gaussian(double avg, double stdev)
{
	if ( verbose & VERB_FULL ) {
		cout << "Generating a random value with a gaussian distribution:" << endl;
		cout << "Avg and std:                    " << avg << " " << stdev << endl;
	}
	
	double			r;
	
	r = -2*log(random()*irm);
	
	if ( r > 0 ) r = avg + stdev*sqrt(r)*cos(TWOPI*(random()*irm));
	else r = avg;
	
	if ( verbose & VERB_FULL )
		cout << "Random value:                   " << r << endl;

	return r;
}

int			random_array_gaussian_chunk(float* r, long start, long end, double avg, double stdev)
{
	double			v;

//	cout << start << tab << end << endl;
	
	for ( long i=start; i<end; i++ ) {
		v = -2*log(random()*irm);
		if ( v > 0 )
			r[i] = avg + stdev*sqrt(v)*cos(TWOPI*random()*irm);
	}
	
	return 0;
}

/**
@brief 	Generates a series with a gaussian random distribution of values.
@param 	n				number of values.
@param 	avg 			average.
@param 	stdev 		standard deviation.
@return float*				the array of values.

	An array of floating point numbers is generated with a gaussian 
	distribution with a given average and standard deviation:
		value = average + std_dev*sqrt(-2*log(random_value))*
						cos(2*PI*random_value);
	where random_value is between 0 and 1.

**/
float* 		random_array_gaussian(long n, double avg, double stdev)
{
	float* 			r = new float[n];
	 
	if ( verbose & VERB_FULL ) {
		cout << "Generating a random series with a gaussian distribution:" << endl;
		cout << "Size:                           " << n << endl;
		cout << "Avg and std:                    " << avg << " " << stdev << endl;
	}
		
	random_seed();
	
	long			chunk_size = get_chunk_size(n);

#ifdef HAVE_GCD
	dispatch_apply((n - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		start = i*chunk_size;
		long		end = start + chunk_size;
		if ( end > n ) end = n;
		random_array_gaussian_chunk(r, start, end, avg, stdev);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<n; i+=chunk_size ) {
		long		end = i + chunk_size;
		if ( end > n ) end = n;
		random_array_gaussian_chunk(r, i, end, avg, stdev);
	}
#endif
	
	return r;
}

/**

@brief 	Generates a value deviating from the average based on a poisson distribution.
@param 	avg 			average.
@return double				value.

	The poisson distribution is given for j = 0,1,... by:
		        avg^j * exp(-avg)
		P(j) = -----------------
		               j!
	Note that only positive integer values are defined for j and sum(P(j)) = 1.
	A value is generated with a poisson distribution with a given average.
	If the average <= 0, the return value is zero.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		random_poisson(double avg)
{
	double 			r = 0;
	
	if ( avg <= 0 ) return r;
	 
	double			sq, lnavg, g, t, y(0);

	if ( avg < 12.0 ) {
		g = exp(-avg);
		for ( r = -1, t = 1; t > g; r++ ) t *= random()*irm;
	} else {
		sq = sqrt(2.0*avg);
		lnavg = log(avg);
		g = avg*lnavg - lgamma(avg + 1.0L);
		do {
			for ( r = -1; r < 0; ) {
				y = tan(M_PI*irm*random());
				r = sq*y + avg;
			}
			r = floor(r);
			t = 0.9*(1 + y*y)*exp(r*lnavg - lgamma(r + 1.0L) - g);
		} while ( random()*irm > t );
	}
	
	return r;
}

int			random_array_poisson_chunk(float* r, long start, long end, double avg)
{
	long 			i;
	double			sq, lnavg, g, t, y;

	if ( avg < 12.0 ) {
		g = exp(-avg);
		for ( i=start; i<end; i++ )
			for ( r[i] = -1, t = 1; t > g; r[i]++ ) t *= random()*irm;
	} else {
		sq = sqrt(2.0*avg);
		lnavg = log(avg);
		g = avg*lnavg - lgamma(avg + 1.0L);
		for ( i=start; i<end; i++ ) {
			do {
				for ( r[i] = -1; r[i] < 0; ) {
					y = tan(M_PI*irm*random());
					r[i] = sq*y + avg;
				}
				r[i] = floor(r[i]);
				t = 0.9*(1 + y*y)*exp(r[i]*lnavg - lgamma(r[i] + 1.0L) - g);
			} while ( random()*irm > t );
		}
	}
	
	return 0;
}

/**

@brief 	Generates a series with a poisson random distribution of values.
@param 	n				number of values.
@param 	avg 			average.
@return float*				the array of values.

	The poisson distribution is given for j = 0,1,... by:
		        avg^j * exp(-avg)
		P(j) = -----------------
		               j!
	Note that only positive integer values are defined for j and sum(P(j)) = 1.
	An array of floating point numbers is generated with a poisson 
	distribution with a given average. The standard deviation is:
		std = sqrt(avg)
	If the average <= 0, the return array contains only zeroes.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
float*		random_array_poisson(int n, double avg)
{
	float* 			r = new float[n];
	
	if ( avg <= 0 ) return r;
	 
	random_seed();

	long			chunk_size = get_chunk_size(n);

#ifdef HAVE_GCD
	dispatch_apply((n - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		start = i*chunk_size;
		long		end = start + chunk_size;
		if ( end > n ) end = n;
		random_array_poisson_chunk(r, start, end, avg);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<n; i+=chunk_size ) {
		long		end = i + chunk_size;
		random_array_poisson_chunk(r, i, end, avg);
	}
#endif
	
	return r;
}

/**
@brief 	Generates a value with a logistical random distribution.
@param 	avg 			average.
@param 	stdev 			standard deviation.
@return double				the random value.

	A floating point number is generated with a logistical 
	differential distribution with a given average and standard deviation:
		value = average + (std_dev/golden)*ln(1/random_value - 1)
	where random_value is between 0 and 1 and:
		golden  = (sqrt(5) + 1)/2
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double 		random_logistical(double avg, double stdev)
{
	double			r, rand_max = get_rand_max();
	double			prefac = stdev/GOLDEN;
	
	if ( verbose & VERB_FULL ) {
		cout << "Generating a value with a logistical differential distribution:" << endl;
		cout << "Avg, std and range:             " << avg << " " << stdev << endl;
	}
	
	r = avg + prefac*log(rand_max/(random()+1.0L) - 1);
	
	return r;
}

int			random_array_logistical_chunk(float* r, long start, long end, double avg, double stdev)
{
	double			rand_max = get_rand_max();
	double			prefac = stdev/GOLDEN;

	for ( long i=start; i<end; i++ )
		r[i] = avg + prefac*log(rand_max/(random()+1.0L) - 1);
	
	return 0;
}

/**
@brief 	Generates an array with a logistical random distribution.
@param 	n				number of values.
@param 	avg 			average.
@param 	stdev 		standard deviation.
@return float*				array of values.

	An array of floating point numbers is generated with a logistical 
	differential distribution with a given average and standard deviation:
		value = average + (std_dev/golden)*ln(1/random_value - 1)
	where random_value is between 0 and 1 and:
		golden  = (sqrt(5) + 1)/2
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
float* 		random_array_logistical(long n, double avg, double stdev)
{
	float* 			r = new float[n];
	 
	if ( verbose & VERB_FULL ) {
		cout << "Generating a random series with a logistical differential distribution:" << endl;
		cout << "Size:                           " << n << endl;
		cout << "Avg, std and range:             " << avg << " " << stdev << endl;
	}
	
	random_seed();
	
	long			chunk_size = get_chunk_size(n);

#ifdef HAVE_GCD
	dispatch_apply((n - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		start = i*chunk_size;
		long		end = start + chunk_size;
		if ( end > n ) end = n;
		random_array_logistical_chunk(r, start, end, avg, stdev);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<n; i+=chunk_size ) {
		long		end = i + chunk_size;
		random_array_logistical_chunk(r, i, end, avg, stdev);
	}
#endif
	
	return r;
}

/**
@brief 	Generates a random vector on the unit sphere.
@return Vector3<double>	vector.

	A random vector is generated with a uniform distribution on the unit sphere.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
Vector3<double> 	vector3_random_unit_sphere()
{
	Vector3<double>	vec;
	
	vec[0] = random_gaussian(0, 1);
	vec[1] = random_gaussian(0, 1);
	vec[2] = random_gaussian(0, 1);
	vec.normalize();
	
	return vec;
}

/**
@brief 	Generates a random vector within a defined cube.
@param	min				minimum.
@param	max				maximum.
@return Vector3<double>	vector.

	Each vector element is set to a random value between the given minimum and maximum.

**/
Vector3<double> 	vector3_random(const double min, const double max)
{
	Vector3<double>	vec;
	double			range = (max - min)/get_rand_max();
	
	vec[0] = random()*range + min;
	vec[1] = random()*range + min;
	vec[2] = random()*range + min;
	
	return vec;
}

/**
@brief 	Generates a random vector within a defined cube.
@param	min				minimum vector.
@param	max				maximum vector.
@return Vector3<double>	vector.

	Each vector element is set to a random value between the given 
	minimum and maximum vectors.

**/
Vector3<double> 	vector3_random(Vector3<double> min, Vector3<double> max)
{
	Vector3<double>	vec;
	Vector3<double>	range = (max - min) / get_rand_max();
	
	vec[0] = random()*range[0] + min[0];
	vec[1] = random()*range[1] + min[1];
	vec[2] = random()*range[2] + min[2];
	
	return vec;
}

/**
@brief 	Generates a random vector within a defined sphere.
@param	length			maximum vector length.
@return Vector3<double>	vector.

	A random vector is generated, normalized and multiplied with
	a random value smaller than the given length.

**/
Vector3<double> 	vector3_random(const double length)
{
	Vector3<double>	vec = vector3_random_unit_sphere();
	
	vec *= random()*length/get_rand_max();
	
	return vec;
}

/**
@brief 	Generates a random vector within a random gaussian-distributed length.
@param 	avg				average of gaussian distribution.
@param 	stdev			standard deviation of gaussian distribution.
@return Vector3<double>	vector.

	A random vector is generated, normalized and multiplied with
	a random value derived from a gaussian distribution.

**/
Vector3<double> 	vector3_random_gaussian(double avg, double stdev)
{
	Vector3<double>	vec = vector3_random_unit_sphere();
	double			amp = random_gaussian(avg, stdev);
	
	vec *= amp;
	
	return vec;
}

/**
@brief 	Generates a random vector within a random gaussian-distributed length in the xy plane.
@param 	avg				average of gaussian distribution.
@param 	stdev			standard deviation of gaussian distribution.
@return Vector3<double>	vector.

	A random vector is generated, the z-component set to zero, normalized 
	and multiplied with a random value derived from a gaussian distribution.

**/
Vector3<double> 	vector3_xy_random_gaussian(double avg, double stdev)
{
	Vector3<double>	vec;
	double			amp = random_gaussian(avg, stdev);
	
	vec[0] = random_gaussian(0, 1);
	vec[1] = random_gaussian(0, 1);
	vec.normalize();
	vec *= amp;
	
	return vec;
}


double		halton_number(long i, long b)
{
	double		r(0), f(1);
	
	while ( i > 0 ) {
		f /= b;
		r += f * (i%b);
		i /= b;
	}
	
	return r;
}

/**
@brief 	Generates a Halton sequence.
@param 	n			length of sequence.
@param 	b			base.
@return double*		sequence.

	The Halton sequence is a pseudo-random array.

**/
double*		halton_sequence(long n, long b)
{
	long		i;
	double*		r = new double[n];
	
	for ( i=0; i<n; i++ )
		r[i] = halton_number(i, b);

	return r;
}

