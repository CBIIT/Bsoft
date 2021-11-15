/**
@file	matrix_linear.cpp
@brief	Solving sets of linear equations through matrix algebra
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20190201
**/
 
#include "Matrix.h" 
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Does a linear least squares fit between two vectors.
@param 	n1			the starting index in each vector (usually 0).
@param 	n2			the final index in each vector.
@param 	*x			x vector (at least n2+1 elements).
@param 	*y			y vector (at least n2+1 elements).
@param 	&a			the intercept.
@param 	&b			the slope.
@return double		the correlation index.

	The two input vectors must have elements between indices n1 and n2.

**/
double		linear_least_squares(int n1, int n2, double *x, double *y, double& a, double& b)
{
    int     	i, n;
    double  	sx, sx2, sy, sxy, sd, dy, d, denom;
	
	// Initial values in case the function returns prematurely
	a = 0;
	b = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG linear_least_squares: Starting fit" << endl;
    
	if ( n2 < n1 ) {
		n = n1;
		n1 = n2;
		n2 = n;
	}
	if ( n1 < 0 ) n1 = 0;
	n = n2 - n1 + 1;
	
    sx = 0;
    sx2 = 0;
    sy = 0;
    sxy = 0;
    for ( i=n1; i<=n2; i++ ) {
    	sx += x[i];
    	sx2 += x[i]*x[i];
    	sy += y[i];
    	sxy += x[i]*y[i];
    }
	
	denom = n*sx2 - sx*sx;
	if ( fabs(denom) < 1e-30 ) return 0;
	
    a = (sx2*sy-sx*sxy)/denom;
    b = (n*sxy-sx*sy)/denom;
	
    sy = sy/n;
    sd = 0;
    dy = 0;
    for ( i=n1; i<=n2; i++ ) {
		d = y[i] - sy;
    	sd += d*d;
		d = a + b*x[i] - y[i];
    	dy += d*d;
    } 
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG linear_least_squares: a = " << a << " b = " << b << endl;
    
    return 1-dy/sd;
}

double		linear_least_squares(int n1, int n2, vector<double>& x, vector<double>& y, double& a, double& b)
{
    int     	i, n;
    double  	sx, sx2, sy, sxy, sd, dy, d, denom;
	
	// Initial values in case the function returns prematurely
	a = 0;
	b = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG linear_least_squares: Starting fit" << endl;
    
	if ( n2 < n1 ) {
		n = n1;
		n1 = n2;
		n2 = n;
	}
	if ( n1 < 0 ) n1 = 0;
	n = n2 - n1 + 1;
	
    sx = 0;
    sx2 = 0;
    sy = 0;
    sxy = 0;
    for ( i=n1; i<=n2; i++ ) {
    	sx += x[i];
    	sx2 += x[i]*x[i];
    	sy += y[i];
    	sxy += x[i]*y[i];
    }

//	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG linear_least_squares: sx=" << sx << tab << " sx2=" << sx2 << endl;
		cout << "DEBUG linear_least_squares: sy=" << sy << tab << " sxy=" << sxy << endl;
//	}
	
	denom = n*sx2 - sx*sx;
	if ( fabs(denom) < 1e-30 ) return 0;
	
    a = (sx2*sy-sx*sxy)/denom;
    b = (n*sxy-sx*sy)/denom;
	
    sy = sy/n;
    sd = 0;
    dy = 0;
    for ( i=n1; i<=n2; i++ ) {
		d = y[i] - sy;
    	sd += d*d;
		d = a + b*x[i] - y[i];
    	dy += d*d;
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG linear_least_squares: a = " << a << " b = " << b << endl;
    
    return 1-dy/sd;
}

/**
@brief 	Fits a data set to a polynomial function.
@param 	n			number of data points.
@param 	&x			x array (at least order+1 values).
@param 	&y			y array (at least order+1 values).
@param 	order		polynomial order.
@param 	&coeff		array in which coefficients are returned (order+1 values)
 	               		(if NULL, no coefficients returned).
@return double		the deviation.

	A polynomial of any order is fitted to the data using a least squares.
	The polynomial is defined as:
		f(x) = a0 + a1*x + a2*x^2 + ...
	The number of coefficients returned is the order plus one.
	The deviation is defined as:
		R = sqrt(sum(y - f(x))^2/n)

**/
double		fit_polynomial(int n, vector<double>& x, vector<double>& y, int order, vector<double>& coeff)
{
	int				i, j, k;
	int				nt(order + 1);
	Matrix			a(nt,nt);
	vector<double>	b(nt);
	vector<double>	v(nt);
	double			f, df, R(0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG fit_polynomial: n=" << n << " order=" << order << endl;
	
	for ( i=0; i<nt; i++ ) b[i] = v[i] = 0;
	
	v[0] = 1;
	for ( i=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0; j<nt; j++ ) b[j] += v[j]*y[i];
		for ( j=0; j<nt; j++ )
			for ( k=0; k<=j; k++ ) a[j][k] += v[j]*v[k];
	}
	for ( j=0; j<nt-1; j++ )
		for ( k=j+1; k<nt; k++ ) a[j][k] = a[k][j];

	a.LU_decomposition(b);
//	a.singular_value_decomposition(b);
	
	for ( i=0, R=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0, f=0; j<nt; j++ ) f += v[j]*b[j];
		df = f - y[i];
		R += df*df;
		if ( verbose & VERB_DEBUG )
			cout << x[i] << " " << y[i] << " " << f << " " << df << endl;
	}
	
	R = sqrt(R/n);
	
	if ( coeff.size() >= nt ) for ( i=0; i<nt; i++ ) coeff[i] = b[i];
	
	if ( verbose & VERB_FULL ) {
		cout << "Polynomial: f(x) = " << b[0];
		for ( i=1; i<nt; i++ ) cout << " + " << b[i] << " x^" << i;
		cout << endl << "R = " << R << endl << endl;
	}
	
	return R;
}

/**
@brief 	Fits a data set to a polynomial function.
@param 	n			number of data points.
@param 	*x			x array (at least order+1 values).
@param 	*y			y array (at least order+1 values).
@param 	order		polynomial order.
@param 	*coeff		array in which coefficients are returned (order+1 values)
 	               		(if NULL, no coefficients returned).
@return double		the deviation.

	A polynomial of any order is fitted to the data using a least squares.
	The polynomial is defined as:
		f(x) = a0 + a1*x + a2*x^2 + ...
	The number of coefficients returned is the order plus one.
	The deviation is defined as:
		R = sqrt(sum(y - f(x))^2/n)

**/
double		fit_polynomial(int n, double* x, double* y, int order, double* coeff)
{
	int				i, j, k;
	int				nt = order + 1;
	Matrix			a(nt,nt);
	vector<double>	b(nt);
	vector<double>	v(nt);
	double			f, df, R=0;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG fit_polynomial: n=" << n << " order=" << order << endl;
	
	for ( i=0; i<nt; i++ ) b[i] = v[i] = 0;
	
	v[0] = 1;
	for ( i=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0; j<nt; j++ ) b[j] += v[j]*y[i];
		for ( j=0; j<nt; j++ )
			for ( k=0; k<=j; k++ ) a[j][k] += v[j]*v[k];
	}
	for ( j=0; j<nt-1; j++ )
		for ( k=j+1; k<nt; k++ ) a[j][k] = a[k][j];

	a.LU_decomposition(b);
//	a.singular_value_decomposition(b);
	
	for ( i=0, R=0; i<n; i++ ) {
		for ( j=1; j<nt; j++ ) v[j] = v[j-1]*x[i];
		for ( j=0, f=0; j<nt; j++ ) f += v[j]*b[j];
		df = f - y[i];
		R += df*df;
		if ( verbose & VERB_DEBUG )
			cout << x[i] << " " << y[i] << " " << f << " " << df << endl;
	}
	
	R = sqrt(R/n);
	
	if ( coeff ) for ( i=0; i<nt; i++ ) coeff[i] = b[i];
	
	if ( verbose & VERB_FULL ) {
		cout << "Polynomial: f(x) = " << b[0];
		for ( i=1; i<nt; i++ ) cout << " + " << b[i] << " x^" << i;
		cout << endl << "R = " << R << endl << endl;
	}

	return R;
}


/**
@brief 	Solves for fitting a plane through a model.
@param 	a				3x3 matrix with cross-terms.
@param 	b				3 vector with averages.
@return Vector3<double>	plane normal.

	A plane is fit through the components and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.

**/
Vector3<double>	fit_plane(Matrix a, vector<double> b)
{
	Vector3<double>		normal;
		
	if ( a[0][0] == 0 ) {
		normal = Vector3<double>(1, 0, 0);
	} else if ( a[1][1] == 0 ) {
		normal = Vector3<double>(0, 1, 0);
	} else if ( a[2][2] == 0 ) {
		normal = Vector3<double>(0, 0, 1);
	} else {
		a.LU_decomposition(b);
		normal = Vector3<double>(b[0], b[1], b[2]);
		normal.normalize();
	}

	return normal;
}
