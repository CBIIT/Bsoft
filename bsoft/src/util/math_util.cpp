/**
@file	math_util.cpp
@brief	Mathematics utility functions
@author Bernard Heymann 
@date	Created: 20030414
@date	Modified: 20151203
**/

#include "math_util.h" 
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Truncates a value to a specified number of decimal places.
@param 	value	value to be truncated.
@param 	places		number of decimal places.
@return int 			0.
**/
double		bfloor(double value, int places)
{
	int 		i;
	double		mult(1);
	
	if ( places < 0 ) places = 0;
	
	for ( i=0; i<places; i++ ) mult *= 10;
	
	return floor(mult*value)/mult;
}

/**
@brief 	Rounds a value to a specified number of decimal places.
@param 	value	value to be rounded.
@param 	places	number of decimal places.
@return double 	rounded value.
**/
double		bround(double value, int places)
{
	int 		i;
	double		mult(1);
	
	if ( places < 0 ) places = 0;
	
	for ( i=0; i<places; i++ ) mult *= 10;
	
	return floor(mult*value + 0.5)/mult;
}

/**

@brief 	Calculates the factorial of n.

	All values of n less than 1 returns 1.
	An exact calculation is done for 1 < n <= 50.
	The Lancos approximation is used for n > 50.
	Factorials of integers larger than 170 exceeds the capacity of a 
	double and causes program termination.
	The largest relative error is for 170: 1.22378e-13.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

@param 	n			integer.
@return double			factorial of n, <0 on error.
**/
double		factorial(int n)
{
	double		fn;
	
	if ( n > 170 ) {
		error_show("Error in factorial", __FILE__, __LINE__);
		cerr << "Factorial of %d exceeds maximum value of a double!" << n << endl;
		return -1;
	}
	
	if ( n <= 1 ) return 1;
	
	if ( n > 50 ) return exp(lgamma(n+1.0L));
	
	for ( fn = 1; n > 1; n-- ) fn *= n;
	
	return fn;
}

/**
@brief 	Calculates the number of combinations of size r in set of size n.
@param 	n			number in set.
@param 	r			number in subset.
@return double		number of combinations, <0 on error.

	All values of n less than 1 returns 1.
	An exact calculation is done for 1 < n <= 50.
	The Lancos approximation is used for n > 50.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		number_of_combinations(int n, int r)
{
	double		nc(0);
	
	if ( n <= 0 ) return nc;
	
	if ( n > 50 ) nc = floor(0.5 + exp(lgamma(n+1.0L) - lgamma(r+1.0L) - lgamma(n-r+1.0L)));
	
	else nc = factorial(n)/(factorial(r)*factorial(n-r));
	
	return nc;
}

/**
@brief 	Determines the index k'th value in the array.
@param 	*a			array.
@param 	n			number of array elements.
@param 	k			rank index to look for.
@return 0

	The array is partioned into 2 sides, with the left side lower or equal 
	to and the right side higher or equal to the k'th element.
	This is useful to determine the median without full sorting.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
int			partition(vector<double>& a, int n, int k)
{
	int			i(0), j(n - 1);
	
	while ( i < j ) {
		while ( i<k && a[i] <= a[k] ) i++;
		while ( j>k && a[j] >= a[k] ) j--;
		if ( i < j ) {
			swap(a[i], a[j]);
			if ( i == k || j==k ) {
				i = 0;
				j = n - 1;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Finds all the prime factor for the input number.
@param 	number		integer.
@param 	&n			number of prime factors.
@return long*		array of prime factors (can be NULL).

	Calculates the prime factors from the smallest to the largest.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
long*  	prime_factors(long number, long& n)
{
	long		divisor, pn[1024], i;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG prime_factors: for number " << number << ":   ";
	
	for ( n=0; number > 1; n++ ) {
		divisor = smallest_prime(number);
		if ( divisor < 1 ) break;
		pn[n] = divisor;
		number /= divisor;
	}
	
	if ( n < 1 ) return NULL;
	
	long*  	prime = new long[n];
	for ( i=0; i<n; i++ ) prime[i] = pn[i];
	
	if ( verbose & VERB_DEBUG ) {
		for ( i=0; i<n; i++ ) cout << " " << pn[i];
		cout << endl;
	}
	
	return prime;
}

/**
@brief 	Finds the smallest prime number factor of the input number.
@param 	number		integer.
@return long		smallest prime factor.

	Tries to divide the given positive integer number by 
	primes from 2 to the square root of the integer.
	Returns the first prime divisor found, which may be the 
	input number if it is prime.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
long		smallest_prime(long number)
{
	if ( number < 2 ) {
		cout << "Zero and one are neither composite nor prime numbers!" << endl;
		return 1;
	}
	
	if ( number%2 == 0 ) return 2;
	if ( number%3 == 0 ) return 3;

	/* Primes generated by alternatingly adding 2 or 4 */
	long			i, inc(2);
	double			max = sqrt((double)number);
	for( i=5; i<=max; i+=inc, inc=6-inc)
		if ( number%i == 0 ) return i;
 
	return number;
}

/**
@brief 	Produces the next permutation in lexical order.
@param 	&s				string of symbols to permute.
@return int				1 = success, 0 = no next permutation.

	To get all permutations, the first string needs
	to be ordered in ascending order.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
int				next_permutation(Bstring& s)
{
	int		i, j;
	long	len = s.length();
	
	for ( i=len-2; i>=0 && s[i] >= s[i+1]; i-- ) ;
	if ( i < 0 ) return 0;
	
	for ( j=len-1; j>=0 && s[i] >= s[j]; j-- ) ;
	s = s.swap(i, j);
	
	for ( i++, j=len-1; i<j; i++, j-- ) s = s.swap(i, j);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG next_permutation: " << s << endl;
	
	return 1;
}

/**
@brief 	Calculates Fisher's z-transform.
@param 	value			a value.
@return double			z value.

	Fisher's z-transform is given by:
		z = 0.5*log((1+v)/(1-v))
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		fishers_z_transform(double value)
{
	double		vmin(1e-30); // Will regularize the unusual case of complete correlation
	
	return 0.5*log((1 + value + vmin)/(1.0 - value + vmin));
}

/**

@brief 	Evaluates the continued fraction for the incomplete beta function.
@param 	a		first parameter.
@param 	b		second parameter.
@param 	x		domain variable [0,1].
@return double			

	Lentz's method.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double			betacf(double a, double b, double x)
{
	int		m,m2;
	double 	aa, c, d, del, h, qab, qam, qap;
	double 	fpmin(1.0e-30);	//Number near the smallest representative
							//floating-point number
	double 	eps(3.0e-7);	//Relative accuracy
	int 	maxit(100);		//Maximum allowed number of iterations
	
	qab=a+b;			// These q's will be used in factors that occur 
	qap=a+1.0;			// in the coefficients (6.4.6).
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < fpmin) d=fpmin;
	d=1.0/d;
	h=d;
	for (m=1; m<=maxit;m++){
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < fpmin) d=fpmin;
		c=1.0+aa/c;
		if (fabs(c) < fpmin) c=fpmin;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;			// Next step of the recurrence (the odd one).
		if (fabs(d) < fpmin) d=fpmin;
		c=1.0+aa/c;
		if (fabs(c) < fpmin) c=fpmin;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < eps) break;		// Are we done?
	}
	if (m > maxit) cerr << "Error: a or b too big, or maxit too small in betacf" << endl;
	return h;
}

/**
@brief 	Calculates the incomplete beta function Ix(a,b).
@param 	a		first parameter.
@param 	b		second parameter.
@param 	x		domain variable [0,1].
@return double			

	Limiting values: I0(a,b) = 0, I1(a,b) = 1.
	Symmetry relation: Ix(a,b) = 1 - I1-x(b,a).
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double			betai (double a, double b, double x)
{
	double 	bt;
	
	if(x < 0.0 || x > 1.0) cerr << "Error: Bad x in routine betai" << endl;
	if(x == 0.0 || x == 1.0) bt=0.0;
	else				//Factors in front of the continued fraction.
		bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
	if(x < (a+1.0)/(a+b+2.0))	//Use continued fraction directly.
		return bt*betacf(a,b,x)/a;
	else						//Use continued fraction after making the sym
		return 1.0-bt*betacf(b,a,1.0-x)/b;	//metry transformation
}

/**
@brief 	Finds a threshold that partitions an array into foreground and background clusters.
@param 	n		number of array elements.
@param 	v		array.
@return double	threshold.	

	Limiting values: I0(a,b) = 0, I1(a,b) = 1.
	Symmetry relation: Ix(a,b) = 1 - I1-x(b,a).
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		kmeans_threshold(long n, double* v)
{
	long			i, nfg, nbg;
	double			t, vfg(v[0]), vbg(v[1]), dvfg(1), dvbg(1), dt(0.001), sumfg, sumbg;
	
	for ( i=0; i<n && vfg == vbg; i++ );	// Ensure at least one value on each class
	
	while ( fabs(dvfg) > dt || fabs(dvbg) > dt ) {	// Exit on no more change
		nfg = nbg = 0;
		sumfg = sumbg = 0;
		t = (vbg + vfg)/2;
		for ( i=0; i<n; i++ ) {	// Sum into two classes
			if ( v[i] > t ) {
				sumfg += v[i];
				nfg++;
			} else {
				sumbg += v[i];
				nbg++;
			}
		}
		if ( nfg ) sumfg /= nfg;
		if ( nbg ) sumbg /= nbg;
		dvfg = sumfg - vfg;
		dvbg = sumbg - vbg;
		vfg = sumfg;
		vbg = sumbg;
	}
	
	return t;
}

