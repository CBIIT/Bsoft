/**
@file	zernike.h
@brief	A class for handling Zernike polynomials
@author Bernard Heymann
@date	Created: 20220105
@date	Modified: 20220107
**/

//#include <vector>
//#include <iostream>
#include "math_util.h"
#include "utilities.h"

//using namespace std;

#ifndef _ZERNIKE_
#define _ZERNIKE_

class Zernike_Coefficient {
private:
	long			n,m,k;
	double			v;
public:
	Zernike_Coefficient(long nn, long mm, long kk) : n(nn), m(mm), k(kk) {
		v = (1 - 2*(k%2)) * factorial(n-k)
			/ (factorial(k) *
			  factorial((n+m)/2 - k) *
			  factorial((n-m)/2 - k));
	}
	void		indices(long& nn, long& mm, long& kk) {
		nn = n; mm = m; kk = k;
	}
	double		value() { return v; }
	bool		is(long nn, long mm, long kk) {
		return (nn==n && abs(mm)==m && kk==k);
	}
} ;

class Zernike {
private:
	vector<Zernike_Coefficient>	coef;
	void		prepare(long nmax) {
		for ( long n=0; n<nmax; ++n ) {
			for ( long m=0; m<=n; ++m ) {
				if ((n - m) % 2 == 1) continue;
				for ( long k=0; k<=(n-m)/2; ++k )
					coef.push_back(Zernike_Coefficient(n,m,k));
			}
		}
	}
public:
	Zernike(long nmax) { prepare(nmax); }
	long		maximum_order() { return long(sqrt(1.5*coef.size())); }
	long		number_of_coefficients() { return coef.size(); }
	double		coefficient(long i) {
		if ( i < coef.size() ) return coef[i].value();
		else return 0;
	}
	double		coefficient(long n, long m, long k) {
		return coefficient(index(n,m,k));
	}
	long		index(long n, long m, long k) {
		if ((n - m) % 2 == 1) {
			cerr << "Error: This combination of m (" << m << ") and n (" << n << ") is not allowed!" << endl;
			return -1;
		}
		for ( long i=0; i<number_of_coefficients(); ++i )
			if ( coef[i].is(n,m,k) ) return i;
		cerr << "Error: This combination of m (" << m << ") and n (" << n << ") is not present!" << endl;
		return -1;
	}
	long		indices(long i, long& n, long& m) {
		long		k(0);
		coef[i].indices(n, m, k);
		return k;
	}
	double		polynomial(long n, long m, double r) {
		if ( n > maximum_order() )
			cerr << "Error: the order (" << n << ") is greater than the maximum order (" << maximum_order() << ")!" << endl;
		if ((n - m) % 2 == 1) return 0.0;
		double		p(0);
		for ( long k=0; k<=(n-m)/2; ++k )
			p += coefficient(n,m,k) * pow(r, n-2*k);
		return p;
	}
	vector<double>	decompose_weight(long n, long m, double w) {
		if ( n > maximum_order() )
			cerr << "Error: the order (" << n << ") is greater than the maximum order (" << maximum_order() << ")!" << endl;
		m = abs(m);
		vector<double>	wv(n+1,0);
		if ((n - m) % 2 == 1) return wv;
		for ( long k=0; k<=(n-m)/2; ++k )
			wv[n-2*k] = coefficient(n,m,k) * w;
		return wv;
	}
	double		Z(long n, long m, double r, double p) {
		if ( m < 0 )
			return polynomial(n, -m, r) * sin(-m*p);
		else
			return polynomial(n, m, r) * cos(m*p);
	}
	double		Zxy(long n, long m, double x, double y) {
		return Z(n, m, sqrt(x*x + y*y), atan2(y,x));
	}
	void		show() {
		cout << "Zernike coefficients:" << endl;
		cout << "#\tn\tm\tk\tc" << endl;
		long		i, n, m, k;
		for ( i=0; i<coef.size(); ++i ) {
			k = indices(i, n, m);
			cout << i << tab << n << tab << m << tab << k << tab << coefficient(i) << endl;
		}
	}
} ;

#endif
