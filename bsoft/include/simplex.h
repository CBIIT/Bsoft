/**
@file	simplex.h
@brief	Nelder and Mead downhill simplex method for generalized parameter fitting
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210729

	Adapted from Numerical Recipes, 2nd edition, Press et al. 1992
	The function "funk" is user-defined and references the "Bsimplex" structure.
	It returns the "R" value and is called as:
		double	R = (funk)(simplex_struct);
**/

#include "utilities.h"

//using namespace std;

/**
@class	Bsimplex
@brief	Structure used in the downhill simplex method.

	Nelder and Mead downhill simplex method for generalized parameter fitting
	Adapted from Numerical Recipes, 2nd edition, Press et al. 1992
	The structure is set up to accommodate any number of variables, parameters,
	constants and points.
	The structure is very flexible in the sense that only some fields
	are absolutely required and with a fixed meaning for the simplex method.
	The required fields are:
		nparam, param, lo, hi
	The other fields may be recast and used as desired in the user function.
	Intended sizes:
		param		nparam.
		lo 			nparam.
		hi 			nparam.
		c			nconstant.
		x			npoint*nvar.
		fx			npoint.
	x or fx can be recast as a different pointer, as long as it is handled
	by the user before calling kill_simplex.
**/
#ifndef _Bsimplex_
class Bsimplex {
private:
	long			nvar;		// Number of variables
	long			nparam;		// Number of parameters
	long			nconstant;	// Number of constants
	long			npoint; 	// Number of function values
	double			yvar;		// Variance of dependent variable
	vector<double>	param;		// Parameter values
	vector<double>	lo;			// Lower limits on parameter values
	vector<double>	hi;			// Upper limits on parameter values
	vector<double>	c;			// Constant values
	vector<double>	x;			// Independent variables: npoint*nvar array
	vector<double>	fx;			// Function values
	double		calculate_dependent_variance() {
		double		yavg(0);
		yvar = 0;
		for ( auto v: fx ) {
			yavg += v;
			yvar += v*v;
		}
		yavg /= fx.size();
		yvar = yvar/fx.size() - yavg*yavg;
		return yvar;
	}
public:
	Bsimplex() { }
	Bsimplex(long nv, long np, long nc, long n, vector<double>& ax, vector<double>& ay);
	long		variables() { return nvar; }
	long		constants() { return nconstant; }
	double		constant(long i) { return c[i]; }
	void		constant(long i, double p) { c[i] = p; }
	void		constants(long n, vector<double>& p) {
		for ( long i=0; i<nconstant && i<n; i++ ) c[i] = p[i];
	}
	long		points() { return npoint; }
	long		parameters() { return nparam; }
	void		parameters(long n, vector<double>& p) {
		for ( int i=0; i<nparam && i<n; i++ ) param[i] = p[i];
	}
	void		parameter(long i, double p) { param[i] = p; }
	double		parameter(long i) { return param[i]; }
	void		limit_low(long i, double v) { lo[i] = v; }
	double		limit_low(long i) { return lo[i]; }
	void		limit_high(long i, double v) { hi[i] = v; }
	double		limit_high(long i) { return hi[i]; }
	void		limits(long i, double vlo, double vhi) { lo[i] = vlo; hi[i] = vhi; }
	void		limits_low(long n, vector<double>& p) {
		for ( long i=0; i<nparam && i<n; i++ ) lo[i] = p[i];
	}
	void		limits_high(long n, vector<double>& p) {
		for ( int i=0; i<nparam && i<n; i++ ) hi[i] = p[i];
	}
	vector<double>& independent_values() { return x; }
	vector<double>& dependent_values() { return fx; }
	double		dependent_variance() {
		if ( yvar <= 0 ) calculate_dependent_variance();
		return yvar;
	}
	double		run(long maxcycles, double tolerance, double (funk)(Bsimplex&));
	double		amotry(vector<double>& mp, vector<double>& R,
		long ihi, double fac, double (funk)(Bsimplex&));
	void		show() {
		cout << "Simplex status:" << endl;
		cout << "Variables: " << nvar << endl;
		cout << "Instances: " << npoint << endl;
		cout << "Parameters: " << nparam << endl;
		for ( long i=0; i<nparam; ++i )
			cout << i+1 << tab << param[i] << tab << lo[i] << tab << hi[i] << endl;
		cout << "Constants: " << nconstant << endl;
		for ( long i=0; i<nconstant; ++i )
			cout << c[i] << endl;
	}
} ;
#define _Bsimplex_
#endif

