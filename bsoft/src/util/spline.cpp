/**
@file	spline.cpp
@brief	Functions to calculate spline curves
@author Bernard Heymann
@date	Created: 20020808
@date	Modified: 20151023
**/

#include "Matrix.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Vector3<double>	vector3_catmull_rom_interpolate(Vector3<double>* pnt, double t)
{
    double			s = 1.0 - t;
    double			t2 = t * t;
    double			t3 = t2 * t;
	double			c0 = -t*s*s;
	double			c1 = 2 - 5*t2 + 3*t3;
	double			c2 = t * ( 1 + 4*t - 3*t2 );
	double			c3 = t2*s;
	Vector3<double> intpnt;

	intpnt[0] = 0.5 * ( c0*pnt[0][0] + c1*pnt[1][0] + c2*pnt[2][0] - c3*pnt[3][0] );
	intpnt[1] = 0.5 * ( c0*pnt[0][1] + c1*pnt[1][1] + c2*pnt[2][1] - c3*pnt[3][1] );
	intpnt[2] = 0.5 * ( c0*pnt[0][2] + c1*pnt[1][2] + c2*pnt[2][2] - c3*pnt[3][2] );

    return intpnt;
}

/**
@brief 	Calculates a 2D/3D spline curve using the Catmull-Rom algorithm.
@param 	ncoord			number of coordinates.
@param 	*coords 			node or point coordinates.
@param 	&nspline		pointer to number of values in spline curve.
@return Vector3<double>*	spline curve.

	A Catmull-Rom spline curve is defined for 4 points {p0,p1,p2,p3} by:
		spline = 0.5*(-t*(1-t)^2)*p0 + (2-5*t^2+3*t^3)*p1 +
		         t*(1+4*t-3*t^2)*p2 - t^2*(1-t)*p3
	where 0 <= t <= 1 is the fractional distance between points p1 and p2.
	Each dimension is interpolated separately.
	The end segments of the spline are defined as straight lines.

**/
Vector3<double>*	vector3_catmull_rom_spline(long ncoord, Vector3<double>* coords, long& nspline)
{
	if ( ncoord < 1 ) return NULL;
	
	long				i, ns;
	double				d, s, t;
	Vector3<double>		dx;
	
	// Calculate the length of the spline curve in voxels
	d = 0;
	for ( i=1; i<ncoord; i++ )
		d += (coords[i] - coords[i-1]).length();
	nspline = (long) (d+1);
	
	Vector3<double>*	spline = new Vector3<double>[nspline];
	
	// Calculate the spline curve by interpolation
	ns = 0;
	t = 0;
	for ( i=1; i<ncoord; i++ ) {
		dx = coords[i] - coords[i-1];
		d = dx.normalize();
		if ( i==1 || i==ncoord-1 ) {
			for ( s=t; s<d; s+=1 ) {
				spline[ns] = (dx * s) + coords[i-1];
				ns++;
			}
			t = s - d;
		} else {
			for ( s=t; s<d; s+=1 ) {
				spline[ns] = vector3_catmull_rom_interpolate(&coords[i-2], s/d);
				ns++;
			}
			t = s - d;
		}
	}
	
	return spline;
}

/*
double		catmull_rom_interpolate(double* pnt, double t)
{
    double		s = 1.0 - t;
    double		t2 = t * t;
    double		t3 = t2 * t;
	double	 	intpnt;

	intpnt = 0.5 * ( -t*s*s*pnt[0] + ( 2 - 5*t2 + 3*t3 ) * pnt[1] +
                  t * ( 1 + 4*t - 3*t2 ) * pnt[2] - t2*s*pnt[3] );

    return intpnt;
}
*/

double 		tps_base_func(double r)
{
	if ( r <= 0.0 )
		return 0.0;
	else
		return r*r * log(r);
}

vector<double>	 thin_plate_splines(vector< Vector3<double> > points, double lambda)
{
	long		np(points.size());
	vector<double>	w(np+3,0);

	if ( points.size() < 3 )
		return w;

	long		i, j;
	double		a(0);
	Matrix 		L(np+3, np+3);

	for ( i=0; i<np; ++i ) {
 		for ( j=i+1; j<np; ++j ) {
			Vector3<double> pt_i = points[i];
			Vector3<double> pt_j = points[j];
			pt_i[2] = pt_j[2] = 0;
			double elen = (pt_i - pt_j).length();
			L[i][j] = L[j][i] = tps_base_func(elen);
			a += elen * 2;
		}
	}
	a /= (double)(np*np);

	for ( i=0; i<np; ++i ) {
		L[i][i] = lambda * a*a;

		L[i][np] = L[np][i] = 1.0;
    	L[i][np+1] = L[np+1][i] = points[i][0];
    	L[i][np+2] = L[np+2][i] = points[i][1];
	}

	for ( i=0; i<np; ++i )
    	w[i] = points[i][2];

	L.LU_decomposition(w);
	
	return w;
}

double		tps_interpolate(Vector3<double>& loc, vector< Vector3<double> > points, vector<double> w)
{
	long			np(points.size());
	double			val = w[np] + w[np+1]*loc[0] + w[np+2]*loc[1];
	Vector3<double>	d;
	
	for ( long i=0; i<np; ++i ) {
		d = loc - points[i];
		d[2] = 0;
		val += w[i] * tps_base_func(d.length());
	}
	
	loc[2] = val;
	
	return val;
}

