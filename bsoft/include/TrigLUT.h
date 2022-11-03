/**
@file	TrigLUT.h
@brief	Trigonometric look-up tables
@author Bernard Heymann
@date	Created: 20220709
@date	Modified: 20220713
**/

#ifndef _TrigLUT_
#define _TrigLUT_

#include "utilities.h"

/************************************************************************
@Object: class TrigLUT
@Description:
	Look-up table class for trigonomtric functions.
@Features:
	One-parameter functions:
		cos, sin
	Two-parameter function:
		atan2
*************************************************************************/
class TrigLUT {
private:
	double			inc;		// Increment
	double			inv_inc;	// Increment
	double			max;		// Maximum value
	vector<double>	v;			// Lookup-table values
public:
	TrigLUT(long n, double (*func)(double)) {
		long		i;
		double		a;
		if ( fabs(func(M_PI/4.0) - sin(M_PI/4.0)) < 1e-10 )
			max = TWOPI;
		else
			max = 1;
		inc = max/n;
		inv_inc = 1/inc;
//		cout << "TrigLUT: " << n << "\t" << inc << endl;
		v.resize(n+1);
//		v.reserve(n+1);
		for ( i=0, a=0; i<=n; ++i, a+=inc ) v[i] = func(a);
//		cout << TWOPI << "\t" << v.size() << endl;
	}
	double	operator[](size_t i) const {
  		try {
    		return v.at(i);
  		}
		catch (const out_of_range& oor) {
			cerr << "Out of range error: " << oor.what() << '\n';
			return 0;
  		}
	}
	double	operator()(double a) {
		a = angle_set_0_to_TWOPI(a);
		double		f = a*inv_inc;
		long		j = f;
		f -= j;
		return (1-f)*v[j] + f*v[j+1];
	}
	double	operator()(double y, double x) {
		double		ax(fabs(x)), ay(fabs(y));
		if ( ax == 0 && ay == 0 ) return 0;
		long		j;
		double		r, f, a;
		if ( ay < ax ) {
			r = ay/ax;
			f = r*inv_inc;
			j = f;
			f -= j;
			a = (1-f)*v[j] + f*v[j+1];
		} else {
			r = ax/ay;
			f = r*inv_inc;
			j = f;
			f -= j;
			a = M_PI_2 - ((1-f)*v[j] + f*v[j+1]);
		}
		if ( x < 0 ) a = M_PI - a;
		if ( y < 0 ) a = -a;
		return a;
	}
	double	accuracy(double (*func)(double)) {
		double		a, d, ds(0);
		for ( a=inc/2; a<max; a+=inc ) {
			d = func(a) - (*this)(a);
			ds += d*d;
		}
		ds = sqrt(ds/v.size());
		return ds;
	}
	void	list(double (*func)(double)) {
		for ( double a=inc/2; a<max; a+=inc )
			cout << a << "\t" << func(a) << "\t" << (*this)(a) << endl;
	}
};
#endif



