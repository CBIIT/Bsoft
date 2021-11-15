/**
@file	Quaternion.h
@brief	Class for quaternions
@author Bernard Heymann
@date	Created: 20051017
@date	Modified: 20150115
**/

#ifndef _Quaternion_
#define _Quaternion_

#include <cmath>

#include "Vector3.h"

template <typename Type> class Vector3;

/************************************************************************
@Object: class Quaternion
@Description:
	Quaternion class for orientations and rotations.
@Features:
	The internal variables are an array of 4 double precision floating point numbers.
*************************************************************************/
class Quaternion {
	double			s;
	Vector3<double>	v;
	void	check() {
		if ( fabs(s) < TRIGPRECISION ) s = 0;
		for ( long i=0; i<3; i++ )
			if ( fabs(v[i]) < TRIGPRECISION ) v[i] = 0;
	}
	template <typename Type>
	void	from_axis_angle(Vector3<Type>& axis, double angle) {
		double			a = angle/2.0L;
		Vector3<double>	unitaxis = axis;
		unitaxis.normalize();
		s = cos(a);
		v = unitaxis*sin(a);
		check();
	}
public:
	Quaternion() : s(0), v(0,0,0) { }
	Quaternion(const double scalar) : s(scalar), v(0,0,0) {
		check();
	}
	Quaternion(const double scalar, const double x, const double y, const double z) :
		s(scalar), v(x,y,z) {
		check();
	}
	Quaternion(const Quaternion& q) : s(q.s), v(q.v) {
		check();
	}
	template <typename Type>
	Quaternion(const Vector3<Type>& vector) :
		s(0), v(vector) {
		check();
	}
	// Note the order - compare Quaternion(axis, angle)
	template <typename Type>
	Quaternion(const double scalar, const Vector3<Type>& vector) :
		s(scalar), v(vector) {
		check();
	}
	template <typename Type>
	// Note the order - compare Quaternion(s,v)
	Quaternion(Vector3<Type>& axis, double angle) {
		from_axis_angle(axis, angle);
	}
	template <typename Type1, typename Type2>
	Quaternion(Vector3<Type1>& from_vec, Vector3<Type2>& to_vec) {
		Vector3<double>	v1(from_vec);
		Vector3<double>	v2(to_vec);
		v1.normalize();
		v2.normalize();
		// Calculate the vector product to get the axis of rotation
		Vector3<double>	axis = v1.cross(v2);
		double			angle = v1.angle(v2);
		if ( axis.length2() < 1e-12 ) {
			v2 = Vector3<double>(1,0,0);
			axis = v1.cross(v2);
			if ( axis.length2() < 1e-12 ) {
				v2 = Vector3<double>(0,1,0);
				axis = v1.cross(v2);
			}
		}
		from_axis_angle(axis, angle);
	}
	Quaternion	operator=(const Quaternion q) {
		s = q.s; v = q.v;
		return *this;
	}
	Quaternion	operator-() {
		s = -s; v = -v;
		return *this;
	}
	Quaternion	operator+=(const Quaternion& q) {
		s += q.s; v += q.v;
		return *this;
	}
	Quaternion	operator+(const Quaternion& q) {
		Quaternion	qn(*this);
		return qn += q;
	}
	Quaternion	operator-=(const Quaternion& q) {
		s -= q.s; v -= q.v;
		return *this;
	}
	Quaternion	operator-(const Quaternion& q) {
		Quaternion	qn(*this);
		return qn -= q;
	}
	Quaternion	operator*=(const Quaternion& q) {
		*this = *this * q;
		return *this;
	}
	Quaternion	operator*(const Quaternion& q) {
		Quaternion	qn;
		qn.s = s*q.s - v.scalar(q.v);
		qn.v[0] = s*q.v[0] + v[0]*q.s + v[1]*q.v[2] - v[2]*q.v[1];
		qn.v[1] = s*q.v[1] + v[1]*q.s + v[2]*q.v[0] - v[0]*q.v[2];
		qn.v[2] = s*q.v[2] + v[2]*q.s + v[0]*q.v[1] - v[1]*q.v[0];
		return qn;
	}
	Quaternion	operator/=(const double d) {
		double	div = 1/d;
		if ( !isfinite(div) ) cerr << "division by zero!" << endl;
		s *= div;
		v *= div;
		return *this;
	}
	Quaternion	operator/(const double d) {
		Quaternion	q(*this);
		q /= d;
		return q;
	}
	double&	operator[](long i) {
		if ( i<1 ) return s;
		if ( i>3 ) i=3;
		return v[i-1];
	}
	double	scalar() { return s; }
	double	angle() { return 2*acos(s); }
	Vector3<double> axis() {
		Vector3<double>		ax(v);
		ax.normalize();
		return ax;
	}
	double	norm() { return sqrt(norm2()); }
	double	norm2() { return s*s + v.length2(); }
	double	normalize() {
		double	size = norm();
		if ( !size ) {
			cerr << "quaternion size = " << size << endl;
			exit(-1);
		}
		*this /= size;
		return size;
	}
	Quaternion conj() { return Quaternion(s, -v); }
	double	invert() {
		double	size2 = norm2();
		double	div = 1/size2;
		s *= div;
		v *= -div;
		return sqrt(size2);
	}
	Quaternion inverse() {
		Quaternion	q(*this);
		q.invert();
		return q;
	}
	Quaternion	rotate(const Quaternion& point) {
		Quaternion	q(*this);
		q.invert();
		return *this * point * q;
	}
} ;

ostream& operator<<(ostream& output, Quaternion& q);

template <typename Type>
inline Vector3<Type>	vector3_rotate(Vector3<Type>& vec, Vector3<Type>& axis, double angle) 
{
	Quaternion		p(vec);
	Quaternion		q(axis, angle);
	
	Quaternion		nq = q.rotate(p);
	
	Vector3<Type>	nuvec = nq.axis();
	
	nuvec *= vec.length();
	
	return nuvec;
}

#endif



