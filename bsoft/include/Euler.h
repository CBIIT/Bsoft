/**
@file	Euler.h
@brief	Euler object
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20150106
**/

#ifndef _Euler_
#define _Euler_

#include "Vector3.h"
#include "Matrix3.h"
#include "Quaternion.h"
#include "View.h"
#include "utilities.h"

/************************************************************************
@Object: class Euler
@Description:
	Euler angle parameter structure.
@Features:
	Three Euler angles:
		psi:	rotation around z
		theta:	rotation around y
		phi:	rotation around z
*************************************************************************/
class Euler {
private:
	double	data[3];
	void	check() {
		for ( int i = 0; i<3; i++ ) data[i] = angle_set_negPI_to_PI(data[i]);
	}
public:
	Euler() { for ( int i = 0; i<3; i++ ) data[i] = 0; }
	Euler(const Euler& e) { for ( int i = 0; i<3; i++ ) data[i] = e.data[i]; }
	Euler(const double psi, const double theta, const double phi) {
		data[0] = psi; data[1] = theta; data[2] = phi;
	}
	Euler(View& v) {
		data[1] = acos(v[2]);
		if ( fabs(v[0]) > 1e-6 || fabs(v[1]) > 1e-6 )
		data[2] = atan2(v[1], v[0]);
		data[0] = angle_set_negPI_to_PI(v.angle() - data[2]);
		if ( fabs(data[0]) < TRIGPRECISION )  data[0] = 0.0L;
	}
	Euler(Matrix3& m) {
		View		v = View(m);
		*this = Euler(v);
	}
	Euler operator=(const Euler& e) {
		for ( int i = 0; i<3; i++ ) data[i] = e.data[i];
		return *this;
	}
	double&	operator[](int i) { if ( i < 0 ) i = 0; if ( i > 2 ) i = 2; return data[i]; }
	double	psi()	{ return data[0]; }
	double	theta() { return data[1]; }
	double	phi()	{ return data[2]; }
	void	psi(const double d)		{ data[0] = d; }
	void	theta(const double d)	{ data[1] = d; }
	void	phi(const double d)		{ data[2] = d; }
	View	view() {
		View		v(cos(data[2])*sin(data[1]), sin(data[2])*sin(data[1]), 
					cos(data[1]), angle_set_negPI_to_PI(data[0] + data[2]));
		return v;
	}
	Matrix3	matrix() {
		double		phi_2 = data[2]/2.0L;
		double		theta_2 = data[1]/2.0L;
		double		psi_2 = data[0]/2.0L;
		Quaternion	q1(cos(phi_2), 0, 0, sin(phi_2));
		Quaternion	q2(cos(theta_2), 0, sin(theta_2), 0);
		Quaternion	q3(cos(psi_2), 0, 0, sin(psi_2));
		Quaternion	q = q1 * q2 * q3;
		return Matrix3(q);
	}
} ;
/*
ostream& operator<<(ostream& output, Euler& euler) {
	output.setf(ios::fixed, ios::floatfield);
	output.precision(4);
	output << "{" << euler.psi() << "," << euler.theta() << "," << euler.phi() << "}" << endl;
	return output;
}
*/
#endif



