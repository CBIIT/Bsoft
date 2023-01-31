/**
@file	simplex.cpp
@brief	Nelder and Mead downhill simplex method for generalized parameter fitting
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220301

	Adapted from Numerical Recipes, 2nd edition, Press et al. 1992
**/

#include "simplex.h"
#include "random_numbers.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Initializes a Bsimplex structure.
@param 	nv		number of independent variables.
@param 	np		number of parameters to determine.
@param 	nc		number of constants.
@param 	n	 	number of points or function values.
@param	&ax		array of independent variables (nv*n).
@param	&ay		array of dependent variables (n).

	The minimum setup of the simplex structure requires a number of
	parameters greater than 0
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
Bsimplex::Bsimplex(long nv, long np, long nc, long n, vector<double>& ax, vector<double>& ay)
{
	if ( np < 1 ) {
		error_show("Error in Bsimplex::Bsimplex: Zero parameters!", __FILE__, __LINE__);
		return;
	}

	if ( n > ax.size() || n > ay.size() ) {
		error_show("Error in Bsimplex::Bsimplex: The number of points is too large!", __FILE__, __LINE__);
		return;
	}

	if ( ax.size() != nv*ay.size() ) {
		error_show("Error in Bsimplex::Bsimplex: The vector sizes are different!", __FILE__, __LINE__);
		return;
	}

	nvar = nv;
	nparam = np;
	nconstant = nc;
	npoint = n;

	if ( nparam && nvar ) x = ax;
	
	if ( npoint ) fx = ay;
	
	param.resize(nparam, 0);
	lo.resize(nparam, 0);
	hi.resize(nparam, 0);
	
	if ( nconstant ) c.resize(nconstant, 0);
	
	calculate_dependent_variance();
}

/**
@brief 	Estimates parameters using a generalized search algorithm.
@param 	maxcycles 		maximum number of evaluation cycles to use.
@param 	tolerance 		absolute tolerance on the R function.
@fn		(funk)(Bsimplex&)	evaluation function returning an R value.
@param 	report 			interval to report R values, default 0 (no reporting).
@return double			final R value (such as a correlation index).

	Downhill simplex method of Nelder and Mead.
	The limits for each parameter must be set in "lo" and "hi".
	If the limits are not specified (i.e., the pointers are NULL) they are
	just ignored.
	The function "funk" is user-defined and references the "Bsimplex" structure.
	It returns the "R" value and is called as:
		double	R = (funk)(simplex_struct);
	The objective is to lower the R value calculated in the evaluation function.
	The tolerance or exit condition parameter is related to the evaluation function,
	and must therefore be supplied by the calling function.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		Bsimplex::run(long maxcycles, double tolerance, double (funk)(Bsimplex&), long report)
{
	long			i, j, k, npnt(nparam + 1), cycle(0);
	long			converged(0), ilo(0), inhi(1), ihi(2);
	long 			rand_max = get_rand_max();	// My own random maximum due to system differences
	double			irm(1.0L/rand_max);
	double			r, Rtry, Rsave, reltol(tolerance/10);
	
	if ( maxcycles < 1 ) maxcycles = 1;
	
	random_seed();
	
	if ( verbose & VERB_DEBUG ) {
		for ( i=0; i<10; i++ )
			cout << "DEBUG Bsimplex::run: Random test: " <<
					rand_max << " " << random() << " " << exp(random()*4.0/rand_max-2.0) << endl;
		cout << "DEBUG Bsimplex::run: Initial parameters:" << endl;
		for ( i=0; i<nparam; ++i )
			cout << i << tab << param[i] << tab << lo[i] << tab << hi[i] << endl;
	}
	
	// Initialize the R's
	vector<double>		R(npnt, 0);
//	for ( i=0; i<npnt; i++ ) R[i] = 0;
	R[ihi] = 1e30;

	for ( j=0; j<nparam; j++ ) {
		if ( param[j] < lo[j] ) param[j] = lo[j];
		if ( param[j] > hi[j] ) param[j] = hi[j];
	}
	
	// Initialize the simplex multiple point parameter array
	vector<double>			mp(npnt*nparam);
	for ( j=0; j<nparam; j++ ) mp[j] = param[j];	// Set up the first point
	for ( i=1; i<npnt; i++ ) {						// Set up the other points
		for ( j=0; j<nparam; j++ ) {
			k = i*nparam+j;
			r = random()*irm;
			if ( lo.size() && hi.size() ) { 		// Stay within limits
				mp[k] = lo[j] + (hi[j] - lo[j])*(0.25 + 0.5*r);
			} else {								// No limits are specified
				if ( param[j] == 0 )
					mp[k] = 2.0*r - 1.0;
				else
					mp[k] = param[j]*exp(4.0*r-2.0);
			}
		}
	}
	
	if ( verbose && report ) {
		cout << "Cycle";
		for ( j=0; j<nparam; j++ ) cout << "\tp" << j;
		cout << "\tR" << endl;
	}
	
	while ( !converged ) {
		for ( i=0; i<npnt; i++ ) {				// Set up simplex points
			for ( j=0; j<nparam; j++ ) param[j] = mp[i*nparam+j];
			R[i] = (funk)(*this);
		}
		ilo = ihi = inhi = 0;					// Rank the points
		for ( i=0; i<npnt; i++ ) {
			if ( R[i] <= R[ilo] ) ilo = inhi = i;
			if ( R[i] > R[ihi] ) ihi = i;
		}
		for ( i=0; i<npnt; i++ )
			if ( R[i] >= R[inhi] && i != ihi ) inhi = i;
			
		// Exit here if the process converged
		if ( cycle >= maxcycles || R[ihi] < tolerance ) {
			converged = 1;
		} else if ( (R[ihi]-R[ilo])/R[ilo] < reltol ) {
				// Set the first point to the lowest
			for ( j=0; j<nparam; j++ ) mp[j] = mp[ilo*nparam+j];
			for ( i=1; i<npnt; i++ ) {			// Randomize other points
				for ( j=0; j<nparam; j++ ) {		
					k = i*nparam+j;
					r = random()*irm;
					if ( lo[j] < hi[j] ) { // Stay within limits
						mp[k] = lo[j] + (hi[j] - lo[j])*(0.25 + 0.5*r);
					} else {					// If no limits are specified
						if ( mp[k] == 0 )
							mp[k] = 2.0*r - 1.0;
						else
							mp[k] = mp[k]*exp(4.0*r-2.0);
					}
				}
			}
		} else {
			// Reflect the simplex opposite the highest point
			Rtry = amotry(mp, R, ihi, -1.0, funk);
			// If it is better, try an additional extrapolation
			if ( Rtry <= R[ilo] )
				Rtry = amotry(mp, R, ihi, 2.0, funk);
			// If it is worse, contract around the lowest point
			else if ( Rtry >= R[inhi] ) {
				Rsave = R[ihi];
				Rtry = amotry(mp, R, ihi, 0.5, funk);
				if ( Rtry >= Rsave ) {
					for ( i=0; i<npnt; i++ ) {
						for ( j=0; j<nparam; j++ ) {
							k = i*nparam+j;
							param[j] = 0.5*(mp[k] + mp[ilo*nparam+j]);
							mp[k] = param[j];
						}
						R[i] = (funk)(*this);
					}
				}
			}
		}
		if ( verbose && report ) {
			if ( cycle%report == 0 ) {
				cout << cycle;
				for ( j=0; j<nparam; j++ ) cout << tab << setw(7) << mp[ilo*nparam+j];
				cout << tab << R[ilo] << endl;
			}
		}
		if ( ! isfinite(R[ilo]) ) bexit(-1);
		cycle++;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bsimplex::run: End conditions: " << cycle << " " << R[ihi] << " " << (R[ihi]-R[ilo])/(R[ihi]+R[ilo]) << endl;
		for ( i=0; i<npnt; i++ ) {			// Show simplex points
			cout << "DEBUG Bsimplex::run: Final: ";
			for ( j=0; j<nparam; j++ ) cout << " " << mp[i*nparam+j];
			cout << " %g" << R[i] << endl;
		}
	}
	
	for ( j=0; j<nparam; j++ ) param[j] = mp[ilo*nparam+j];
	
	Rsave = R[ilo];
	
	if ( verbose && report ) {
		cout << "Simplex cycles used:            " << cycle << endl;
		cout << "R range:                        " << setprecision(6) << R[ihi] - R[ilo] << endl << endl;
	}

	return Rsave;
}

/*
@brief 	Generates new points in parameter space for the simplex function.
@param 	*mp			multiple point matrix.
@param 	*R			R vector.
@param 	ihi 		index of highest R.
@param 	fac			direction and magnitude of simplex modification.
@fn		(funk)(Bsimplex *)	evaluation function returning an R value.
@return double 		an R value.

	Downhill simplex method of Nelder and Mead.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		Bsimplex::amotry(vector<double>& mp, vector<double>& R,
				long ihi, double fac, double (funk)(Bsimplex&))
{
	long			i, j, npnt(nparam + 1);
	double			psum, pnew, Rtry;
	double			fac1((1 - fac)/nparam);
	double			fac2(fac1 - fac);
	
	for ( j=0; j<nparam; j++ ) {
		psum = 0;
		for ( i=0; i<npnt; i++ ) psum += mp[i*nparam+j];
		param[j] = mp[ihi*nparam+j];
		pnew = psum*fac1 - mp[ihi*nparam+j]*fac2;
		if ( lo[j] < hi[j] ) {
			if ( pnew < lo[j] ) pnew = param[j];
			if ( pnew > hi[j] ) pnew = param[j];
		}
		param[j] = pnew;
	}
	
	Rtry = (funk)(*this);
	
	if ( Rtry <= R[ihi] ) {
		R[ihi] = Rtry;
		for ( j=0; j<nparam; j++ ) mp[ihi*nparam+j] = param[j];
	}
	
	return Rtry;
}

