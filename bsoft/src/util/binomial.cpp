/**
@file	binomial.cpp
@brief	Functions to generate and fit binomials
@author Daniel Nemecek and Bernard Heymann
@date	Created: 20091202
@date	Modified: 20190201
**/

#include "simplex.h"
#include "math_util.h" 
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function
double		Bnp_R(Bsimplex* simp);

/**
@author Daniel Nemecek and Bernard Heymann
@brief 	Calculates the probability of k occurrences in n cases.
@param 	n			number of cases.
@param 	p			probability of one occurrence [0-1].
@param 	k			occurrences (<=n).
@param 	w			weight [0-1].
@return double 		weighted probability.
**/
double		Bnpk(int n, double p, int k, double w)
{
	int			i, nk(n-k);
	double		p_1 = 1 - p;
	double		Pnpk(1);
	
	double		nc = number_of_combinations(n, k);
	
	for ( i=0; i<k; i++ ) Pnpk *= p;
	for ( i=0; i<nk; i++ ) Pnpk *= p_1;
	Pnpk *= nc * w;
	
	return Pnpk;
}
	
double		Bnp_R(Bsimplex& simp)
{
	int			i,j;
	double		R = 0;	// normalized standard deviation
	double		df;		// value of the minimization function
	double		fit, w;
	vector<double>&	f = simp.dependent_values();
	
	int			n = (int) simp.constant(0);
	int			nfit = (int) simp.constant(1);

	// Normalize weights
	for ( j=0, w=0; j<nfit; j++ ) w += simp.parameter(j);
	
	if ( w > 0 )
		for ( j=0; j<nfit; j++ ) simp.parameter(j, simp.parameter(j) / w);
	else
		for ( j=0; j<nfit; j++ ) simp.parameter(j, 1.0/nfit);

	// Calculate residual
	for ( i=0; i<simp.points(); i++ ) {
		fit = 0;
		for ( j=0; j<nfit; j++ )
			fit += Bnpk(n, simp.parameter(nfit+j), i, simp.parameter(j));
		df = f[i] - fit;
		R += df*df;
	}	
	
	R = sqrt(R/i);
	
	return R;
}


/**
@author Daniel Nemecek and Bernard Heymann
@brief 	Fits a distribution to one or more binomial distributions.
@param 	&distrib		distribution.
@param 	n				number of cases (number of values in distribution).
@param 	nfit			number of binomial curves to fit.
@param 	&coeff			fitted coefficients (weights and probabilities).
@return double 			residual.

	The binomial curve(s) is(are) fit using the downhill simplex method.

**/
double		binomial_fit(vector<double>& distrib, int n, int nfit, vector<double>& coeff)
{			
	int				i;
	vector<double>	x(n+1, 0);
	Bsimplex		simp(1, 2*nfit, 2, n+1, x, distrib);
	
	for ( i=0; i<2*nfit; i++ )
		simp.limits(i, 0, 1);			// probability/weight = [0,1]

	if ( nfit < 2 ) simp.limits(0, 1, 1);
	
	for ( i=0; i<nfit; i++ ) {
		simp.parameter(i, 1.0/nfit);		// initial weight
		simp.parameter(i+nfit, (i + 1)*1.0/(nfit + 1));	// initial probability
	}

	simp.constant(0, n);		// number of x-points
	simp.constant(1, nfit);		// number of fitted B(n,p)
	
	// call the simplex function that minimizes Bnp_R
	double			R = simp.run(10000, 0.001, Bnp_R);
	
	for ( i=0; i<2*nfit; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

