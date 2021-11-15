/**
@file	UnitCell.h
@brief	Header file for unit cell functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20150103
**/

#include "Matrix3.h"
#include "utilities.h"

#ifndef _UnitCell_
/************************************************************************
@Object: struct UnitCell
@Description:
	Unit cell parameters.
@Features:
	The convention is as follows:
		a lies on the x-axis, 
		b lies in the x,y plane. 
*************************************************************************/
class UnitCell {
private:
	double   data[6];
public:
	UnitCell() { for ( int i=0; i<3; i++ ) data[i] = 0; for ( int i=3; i<6; i++ ) data[i] = M_PI_2; }
	UnitCell(double* v) { for ( int i=0; i<6; i++ ) data[i] = v[i]; }
	UnitCell(double a, double b, double c, double alf, double bet, double gam) {
		data[0] = a; data[1] = b; data[2] = c; 
		data[3] = alf; data[4] = bet; data[5] = gam;
	}
//	UnitCell operator=(const UnitCell& uc) { for ( int i=0; i<6; i++ ) data[i] = uc.data[i]; return *this; }
	double&	operator[](int i) { if ( i < 0 ) i = 0; if ( i > 5 ) i = 5; return data[i]; }
	bool	check() { double e=1; for ( int i=0; i<6; i++ ) e *= data[i]; return e > 0; }
	void	size(double a, double b, double c) { data[0] = a; data[1] = b; data[2] = c; }
	void	a(double d) { data[0] = d; }
	void	b(double d) { data[1] = d; }
	void	c(double d) { data[2] = d; }
	void	alpha(double d) { data[3] = d; }
	void	beta(double d) { data[4] = d; }
	void	gamma(double d) { data[5] = d; }
	double	a() { return data[0]; }
	double	b() { return data[1]; }
	double	c() { return data[2]; }
	double	alpha() { return data[3]; }
	double	beta() { return data[4]; }
	double	gamma() { return data[5]; }
	void	set_angle_range() {
		for ( int i=3; i<6; i++ ) {
			while ( data[i] < -M_PI ) data[i] += TWOPI;
			while ( data[i] > M_PI ) data[i] -= TWOPI;
		}
	}
	void	degrees_to_radians() { for ( int i=3; i<6; i++ ) data[i] *= M_PI/180.0; }
	void	radians_to_degrees() { for ( int i=3; i<6; i++ ) data[i] *= 180.0/M_PI; }
	double	volume() { return data[0]*data[1]*data[2]*(1 - cos(data[3])*cos(data[3])
						- cos(data[4])*cos(data[4]) - cos(data[5])*cos(data[5])
						+ 2*cos(data[3])*cos(data[4])*cos(data[5])); }
	Matrix3	skew_matrix();
	Matrix3	skew_rotation(int invert);
} ;
#define _UnitCell_
#endif


