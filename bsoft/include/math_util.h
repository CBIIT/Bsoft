/**
@file	math_util.h 
@brief	Header file for general utilities 
@author Bernard Heymann 
@date	Created: 19990722
@date	Modified: 20220524
**/

#include <stdlib.h>
#include <cmath> 
#ifdef SUNOS
#include <ieeefp.h>
#endif

#include "Bstring.h"
#include "utilities.h" 

// Function prototypes 
double		bfloor(double value, int places);
double		bround(double value, int places);
double		sinc(double d);
double		factorial(int n);
double		number_of_combinations(int n, int r);
int			partition(vector<double>& a, int n, int k);
int			partition(double* a, int n, int k);
long*  		prime_factors(long number, long& n);
long		smallest_prime(long number);
int			next_permutation(Bstring& s);
double		fishers_z_transform(double value);
double		betacf(double a, double b, double x);
double		betai (double a, double b, double x);
double		kmeans_threshold(long n, double* v);


