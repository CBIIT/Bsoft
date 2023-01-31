/**
@file	moving_average.cpp
@brief	Functions for moving average calculations
@author 	Bernard Heymann
@date	Created: 20000430
@date	Modified: 20221129
**/

#include "moving_average.h"
#include "utilities.h" 
 
// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates a moving average over an array of data.
@param 	number		number of values in the array.
@param 	*x			the array.
@param 	window		sliding window length.
@return double*		the moving average array.

	All data points within a sliding window are averaged.
	The window moves over the ends of the array and averages only the 
	available points.
	A new array is allocated and the moving averages returned.

**/
vector<double>	moving_average(long number, double* x, long window)
{
	vector<double>	v(number);
	
	for ( long i=0; i<number; ++i ) v[i] = x[i];
	
	return moving_average(v, window);
}

/**
@brief 	Calculates a moving average over an array of data.
@param 	&x				the array.
@param 	window			sliding window length.
@return vector<double>	the moving average array.

	All data points within a sliding window are averaged.
	The window moves over the ends of the array and averages only the
	available points.
	A new array is allocated and the moving averages returned.

**/
vector<double>	moving_average(vector<double>& x, long window)
{
	long			n(x.size());
	long 			i, j, hw(window/2);
	vector<double>	mov_avg(n,0.0);
	
	if ( verbose & VERB_PROCESS )
		cout << "Calculating moving average with window of " << window << " points" << endl << endl;
	
	// Calculate the moving sum
	for ( i=j=0; i<n; ++i ) {
		if ( i>0 ) mov_avg[i] = mov_avg[i-1];
		if ( j-window >= 0 ) {
			mov_avg[i] -= x[j-window];
			if ( j>=n ) j++;
		}
		while ( j<i+hw && j<n) {
			mov_avg[i] += x[j];
			j++;
		}
	}
	
	// Normalize the moving sum
	for ( i=0; i<n; ++i ) {
		j = window;
		if ( i-hw < 0 ) j += i - hw;
		if ( i+hw >= n ) j -= i + hw - n;
		mov_avg[i] /= j;
	}
		
	return mov_avg;
}

/**
@brief 	Calculates a moving average over an array of complex data.
@param 	&x				the array.
@param 	window			sliding window length.
@return Complex<float>*	the moving average complex array.

	All data points within a sliding window are averaged.
	The window moves over the ends of the array and averages only the
	available points.
	A new array is allocated and the moving averages returned in it.

**/
vector<Complex<float>>	moving_average_complex(vector<Complex<float>>& x, long window)
{
	long 					i, j(0), number(x.size()), halfwindow(window/2);
	vector<Complex<float>>	mov_avg(number);
	for ( i=0; i<number; i++ ) mov_avg[i] = 0;
	
	if ( verbose & VERB_PROCESS )
		cout << "Calculating complex moving average with window of " << window << " points" << endl << endl;
	
	// Calculate the moving sum
	for ( i=0; i<number; i++ ) {
		if ( i>0 ) mov_avg[i] = mov_avg[i-1];
		if ( j-window >= 0 ) {
			mov_avg[i] -= x[j-window];
			if ( j>=number ) j++;
		}
		while ( j<i+halfwindow && j<number) {
			mov_avg[i] += x[j];
			j++;
		}
	}
	
	// Normalize the moving sum
	for ( i=0; i<number; i++ ) {
		j = window;
		if ( i-halfwindow < 0 ) j += i - halfwindow;
		if ( i+halfwindow >= number ) j -= i + halfwindow - number;
		mov_avg[i] /= j;
	}
		
	return mov_avg;
}

/**
@brief 	Calculates a moving polynomial fit over an array of data.
@param 	order			polynomial order.
@param 	number			number of values in the array.
@param 	*x				the array.
@param 	window			sliding window length.
@return vector<double>	the moving polynomial fit array.

	All data points within a sliding window are fit to a polynomial.
	The window moves over the ends of the array and fits only the 
	available points.
	A new array is allocated and the moving polynomial fit values returned.

**/
vector<double>	moving_polynomial(long order, long number, double* x, long window)
{
	vector<double>	v(number);
	
	for ( long i=0; i<number; ++i ) v[i] = x[i];
	
	return moving_polynomial(order, v, window);
}
/*
vector<double>	moving_polynomial(long order, long number, double* x, long window)
{
	vector<double>	m;

	if ( order < 1 ) {
		cerr << "Error in moving_polynomial: the order must be > 0!" << endl;
		return m;
	}
	
	long			hw = window/2, i, j, is, i1, i2, w;
	double			coeff[order+1], v[window], p;
	m.resize(number, 0);

	for ( i=0; i<window; i++ ) v[i] = i;

	for ( i=0; i<number; i++ ) {
		is = 0;
		i1 = i - hw;
		if ( i1 < 0 ) { is = -i1; i1 = 0; }
		i2 = i + hw;
		if ( i2 >= number ) i2 = number - 1;
		w = i2 - i1 + 1;
		fit_polynomial(w, &v[is], &x[i1], order, coeff);
		for ( j=0, m[i]=0, p=1; j<=order; j++, p*=v[hw] ) m[i] += coeff[j]*p;
	}

	return m;
}
*/

/**
@brief 	Calculates a moving polynomial fit over an array of data.
@param 	order			polynomial order.
@param 	&x				the array.
@param 	window			sliding window length.
@return vector<double>	the moving polynomial fit array.

	All data points within a sliding window are fit to a polynomial.
	The window moves over the ends of the array and fits only the
	available points.
	A new array is allocated and the moving polynomial fit values returned.

**/
vector<double>	moving_polynomial(long order, vector<double>& x, long window)
{
	vector<double>	m;

	if ( order < 1 ) {
		cerr << "Error in moving_polynomial: the order must be > 0!" << endl;
		return m;
	}
	
	long			number(x.size());
	long			hw = window/2, i, j, is, i1, i2, w;
	double			coeff[order+1], v[window], p;
	m.resize(number, 0);

	for ( i=0; i<window; i++ ) v[i] = i;

	for ( i=0; i<number; i++ ) {
		is = 0;
		i1 = i - hw;
		if ( i1 < 0 ) { is = -i1; i1 = 0; }
		i2 = i + hw;
		if ( i2 >= number ) i2 = number - 1;
		w = i2 - i1 + 1;
		fit_polynomial(w, &v[is], &x[i1], order, coeff);
		for ( j=0, m[i]=0, p=1; j<=order; j++, p*=v[hw] ) m[i] += coeff[j]*p;
	}

	return m;
}

/**
@brief 	Calculates a moving fit of the local gradient in an array of data.
@param 	&x				the array.
@param 	window			sliding window length.
@return vector<double>	the moving gradient fit array.

	All data points within a sliding window are fit to a line.
	The window moves over the ends of the array and fits only the
	available points.

**/
vector<double>	moving_gradient(vector<double>& x, long window)
{
	long			number(x.size());
	long			hw(window/2), i, is, i1, i2, w;
	double			coeff[2], v[window];
	vector<double>	m(number, 0);

	for ( i=0; i<window; i++ ) v[i] = i;

	for ( i=0; i<number; i++ ) {
		is = 0;
		i1 = i - hw;
		if ( i1 < 0 ) { is = -i1; i1 = 0; }
		i2 = i + hw;
		if ( i2 >= number ) i2 = number - 1;
		w = i2 - i1 + 1;
		fit_polynomial(w, &v[is], &x[i1], 1, coeff);
		m[i] = coeff[1];
	}

	return m;
}

/**
@brief 	Calculates a moving fit of the local curvature in an array of data.
@param 	&x				the array.
@param 	window			sliding window length.
@return vector<double>	the moving curvature fit array.

	All data points within a sliding window are fit to a third-order polynomial.
	The window moves over the ends of the array and fits only the
	available points.

**/
vector<double>	moving_curvature(vector<double>& x, long window)
{
	long			number(x.size());
	long			hw(window/2), i, is, i1, i2, w;
	double			coeff[3], v[window];
	vector<double>	m(number, 0);

	for ( i=0; i<window; i++ ) v[i] = i;

	for ( i=0; i<number; i++ ) {
		is = 0;
		i1 = i - hw;
		if ( i1 < 0 ) { is = -i1; i1 = 0; }
		i2 = i + hw;
		if ( i2 >= number ) i2 = number - 1;
		w = i2 - i1 + 1;
		fit_polynomial(w, &v[is], &x[i1], 3, coeff);
		m[i] = coeff[2];
	}

	return m;
}

