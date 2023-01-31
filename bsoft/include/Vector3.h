/**
@file	Vector3.h
@brief	Class for 3-value vectors
@author Bernard Heymann
@date	Created: 20000501
@date	Modified: 20190708
**/


#include <iostream>
#include <cmath>
#include <vector>

#include "math_util.h"
#include "utilities.h"

using namespace std;

#ifndef _Vector3_
#define _Vector3_
/**
@class Vector3
@brief	Vector class for 3-value vectors used in 3D space.
	The internal variables are an array of 3 numbers.
**/
template <typename Type>
class Vector3 {
private:
	Type	data[3];
public:
	Vector3() { data[0] = 0; data[1] = 0; data[2] = 0; }
	Vector3(const Vector3& v) { data[0] = v.data[0]; data[1] = v.data[1]; data[2] = v.data[2]; }
	Vector3(const Type x, const Type y, const Type z) { data[0] = x; data[1] = y; data[2] = z; }
	template <typename T2>
	Vector3(const vector<T2>& v) { for ( size_t i=0; i<3 && i<v.size(); i++ ) data[i] = v[i]; }
	Vector3	operator=(const Type d) { data[0] = d; data[1] = 0; data[2] = 0; return *this; }
	Vector3	operator=(const Vector3& v) {
		data[0] = v.data[0]; data[1] = v.data[1]; data[2] = v.data[2];
		return *this;
	}
	template <typename T2>
	Vector3	operator=(const vector<T2>& v) {
		for ( size_t i=0; i<3 && i<v.size(); i++ ) data[i] = v[i];
		return *this;
	}
	Vector3	operator-() {
		Vector3<Type>	vn(-data[0], -data[1], -data[2]);
		return vn;
	}
	Vector3	operator+=(const double d) {
		for ( int i=0; i<3; i++ ) data[i] += Type(d);
		return *this;
	}
	Vector3	operator+(const double d) {
		Vector3<Type>	vn(Type(data[0]+d), Type(data[1]+d), Type(data[2]+d));
		return vn;
	}
	Vector3	operator+=(const Vector3& v) {
		for ( int i=0; i<3; i++ ) data[i] += v.data[i];
		return *this;
	}
	Vector3	operator+(const Vector3& v) {
		Vector3<Type>	vn(data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2]);
		return vn;
	}
	Vector3	operator-=(const double d) {
		for ( int i=0; i<3; i++ ) data[i] -= Type(d);
		return *this;
	}
	Vector3	operator-(const double d) {
		Vector3<Type>	vn(Type(data[0]-d), Type(data[1]-d), Type(data[2]-d));
		return vn;
	}
	Vector3	operator-=(const Vector3& v) {
		for ( int i=0; i<3; i++ ) data[i] -= v.data[i];
		return *this;
	}
	Vector3	operator-(const Vector3& v) {
		Vector3<Type>	vn(data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2]);
		return vn;
	}
	Vector3	operator*=(const double d) {
		for ( int i=0; i<3; i++ ) data[i] = Type(data[i] * d);
		return *this;
	}
	Vector3	operator*(const double d) {
		Vector3<Type>	vn(Type(data[0] * d), Type(data[1] * d), Type(data[2] * d));
		return vn;
	}
	Vector3	operator*=(const Vector3& v) {	// Note that the order is important for vectors of different types
		for ( int i=0; i<3; i++ ) data[i] *= v.data[i];
		return *this;
	}
	Vector3	operator*(const Vector3& v) {
		Vector3<Type>	vn(data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2]);
		return vn;
	}
	Vector3 operator/=(const double d) {
		if ( fabs(d) < 1e-30 ) return *this;
		double	div = 1/d;
		for ( int i=0; i<3; i++ ) data[i] *= div;
		return *this;
	}
	Vector3	operator/(const double d) {
		Vector3<Type>	vn(Type(data[0] / d), Type(data[1] / d), Type(data[2] / d));
		return vn;
	}
	Vector3	operator/=(const Vector3& v) {
		for ( int i=0; i<3; i++ ) data[i] /= v.data[i];
		return *this;
	}
	Vector3	operator/(const Vector3& v) {
		Vector3<Type>	vn(data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2]);
		return vn;
	}
	bool	operator==(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] == v.data[i] );
		return ( e == 3 );
	}
	bool	operator==(const double d) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] == d );
		return ( e == 3 );
	}
	bool	operator!=(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] == v.data[i] );
		return ( e != 3 );
	}
	bool	operator>(const double d) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] > d );
		return ( e == 3 );
	}
	bool	operator<(const double d) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] < d );
		return ( e == 3 );
	}
	bool	operator>(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] > v.data[i] );
		return ( e == 3 );
	}
	bool	operator<(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] < v.data[i] );
		return ( e == 3 );
	}
	bool	operator>=(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] >= v.data[i] );
		return ( e == 3 );
	}
	bool	operator<=(const Vector3& v) {
		int		i, e(0);
		for ( i=0; i<3; i++ ) e += ( data[i] <= v.data[i] );
		return ( e == 3 );
	}
	Type	operator[](size_t i) const { while ( i>2 ) i-=3; return data[i]; }
	Type&	operator[](size_t i) { while ( i>2 ) i-=3; return data[i]; }
	template <typename T2> operator Vector3<T2>() const {
		return Vector3<T2>((T2)data[0], (T2)data[1], (T2)data[2]);
	}
	vector<double>	array() {
		return {data[0], data[1], data[2]};
	}
	Type	min() {
		Type d(data[0]);
		if ( d > data[1] ) d = data[1];
		if ( d > data[2] ) d = data[2];
		return d;
	}
	Type	max() {
		Type d(data[0]);
		if ( d < data[1] ) d = data[1];
		if ( d < data[2] ) d = data[2];
		return d;
	}
	Vector3	min(const double d) {
		Vector3		vn(*this);
		for ( int i=0; i<3; i++ ) if ( vn.data[i] > d ) vn.data[i] = (Type)d;
		return vn;
	}
	Vector3	max(const double d) {
		Vector3		vn(*this);
		for ( int i=0; i<3; i++ ) if ( vn.data[i] < d ) vn.data[i] = (Type)d;
		return vn;
	}
	Vector3	min(const Vector3& v) {
		Vector3		vn(*this);
		for ( int i=0; i<3; i++ ) if ( vn.data[i] > v.data[i] ) vn.data[i] = v.data[i];
		return vn;
	}
	Vector3	max(const Vector3& v) {
		Vector3		vn(*this);
		for ( int i=0; i<3; i++ ) if ( vn.data[i] < v.data[i] ) vn.data[i] = v.data[i];
		return vn;
	}
	Vector3	abs() {
		Vector3		vn(*this);
		for ( int i=0; i<3; i++ ) if ( data[i] < 0 ) vn.data[i] = -data[i];
		return vn;
	}
	Vector3	floor(int places) {
		Vector3		vn;
		for ( int i=0; i<3; i++ ) vn.data[i] = bfloor(data[i], places);
		return vn;
	}
	Vector3	round(int places) {
		Vector3		vn;
		for ( int i=0; i<3; i++ ) vn.data[i] = bround(data[i], places);
		return vn;
	}
	Vector3	remainder(int divisor) {
		Vector3		vn;
		for ( int i=0; i<3; i++ ) vn.data[i] = fmod(data[i], divisor);
		return vn;
	}
	double	length2() { return (double)data[0]*data[0] + (double)data[1]*data[1] + (double)data[2]*data[2]; }
	double	length() { return sqrt(length2()); }
	void	length(double d) {
		double		len(length());
		if ( len ) *this *= d/len;
	}
	double	distance(const Vector3& v) { return (*this - v).length(); }
	double	distance2(const Vector3& v) { return (*this - v).length2(); }
	double	distance_along_vector(Vector3 v) {
		return fabs(scalar(v)/v.length());
	}
	double	distance_from_line2(Vector3 v1, Vector3 v2) {
		Vector3		v1d(*this - v1);
		Vector3		v2d(*this - v2);
		if ( v1 == v2 ) return v1d.length2();
		Vector3		vd = v2 - v1;
		return v1d.cross(v2d).length2()/vd.length2();
	}
	double	distance_from_line(Vector3 v1, Vector3 v2) {
		return sqrt(distance_from_line2(v1, v2));
	}
	Vector3	closest_point_on_line(Vector3 v1, Vector3 v2) {
		Vector3		v1d(*this - v1);
		if ( v1 == v2 ) return v1d.length2();
		Vector3		vd = v2 - v1;
		return v1 + v1d * (vd * vd)/vd.length2();
	}
	double	position_relative_to_line(Vector3 v1, Vector3 v2) {
		if ( v1 == v2 ) return 0;
		Vector3		v1d(*this - v1);
		Vector3		vd = v2 - v1;
		return v1d.scalar(vd)/vd.length2();
	}
	Vector3	square_root() {
		Vector3		v(*this);
		for ( int i=0; i<3; i++ ) {
			if ( v.data[i] > 0 ) v.data[i] = sqrt(data[i]);
			else v.data[i] = 0;
		}
		return v;
	}
	double	sum() { return (double)data[0] + (double)data[1] + (double)data[2]; }
	double	volume() { return (double)data[0]*data[1]*data[2]; }
	double	normalize() {
		double	len = length();
		if ( len < 1e-30 ) {
//			len = 1;
			data[2] = 1;
		} else {
			double	div = 1.0/len;
			for ( int i=0; i<3; i++ ) data[i] *= div;
		}
		return len;
	}
	double	scalar(const Vector3& v) const {
		return (double)data[0]*v.data[0] + (double)data[1]*v.data[1] + (double)data[2]*v.data[2];
	}
	Vector3	cross(const Vector3& v) {
		Vector3		vn((double)data[1]*v.data[2] - (double)data[2]*v.data[1], 
						(double)data[2]*v.data[0] - (double)data[0]*v.data[2], 
						(double)data[0]*v.data[1] - (double)data[1]*v.data[0]);
		return vn;
	}
	template <typename T2>
	double	angle(Vector3<T2>& v) {
		Vector3		v1(*this);
		Vector3		v2(v);
		v1.normalize();
		v2.normalize();
		double		prod = v1.scalar(v2);
		if ( prod > 1 ) prod = 1;
		if ( prod < -1 ) prod = -1;
		return(acos(prod));
	}
	Vector3	normal(Vector3& v1, Vector3& v2) {
		Vector3 	edge1 = *this - v1;
		Vector3 	edge2 = *this - v2;
		Vector3 	normal = edge1.cross(edge2);
		normal.normalize();
		return(normal);
	}
	bool finite() {
		int			e(0);
		for ( int i=0; i<3; i++ ) e += isfinite(data[i]);
		return ( e == 3 );
	}
	bool notfinite() {
		int			e(0);
		for ( int i=0; i<3; i++ ) e += isfinite(data[i]);
		return ( e != 3 );
	}
	template <typename T1, typename T2>
	bool within(Vector3<T1>& v1, Vector3<T2>& v2) {
		int			i, e(0);
		for ( i=0; i<3; i++ ) e += (data[i] >= v1[i]);
		for ( i=0; i<3; i++ ) e += (data[i] <= v2[i]);
		return ( e == 6 );
	}
//	friend Vector3<float>	operator=(Vector3<double>& vd);
//	friend Vector3<double>	operator=(Vector3<float>& vf);
} ;

template <typename Type>
ostream& operator<<(ostream& output, Vector3<Type> vec) {
	output.setf(ios::fixed, ios::floatfield);
	output << vec[0] << tab << vec[1] << tab << vec[2];
	return output;
}

template <typename Type>
ostream& operator<<(ostream& output, vector<Type> vec) {
	output.setf(ios::fixed, ios::floatfield);
	output << vec[0];
	for ( auto it = vec.begin()+1; it != vec.end(); ++it )
		cout << tab << *it;
	return output;
}
/*
template <typename T1, typename T2>
//inline vector<T1>		operator=(vector<T1> v, Vector3<T2> v3)
inline void		operator=(vector<T1>& v, Vector3<T2> v3)
//inline vector<T1>::operator=(Vector3<T2> v3)
{
	for ( long i = 0; i<3; ++i )
		v[i] = v3[i];
	
	return v;
}

template <typename T> operator vector<T>() const {
		return Vector3<T2>((T2)data[0], (T2)data[1], (T2)data[2]);
	}

template <typename T1, typename T2> operator vector<T1>(Vector3<T2> v3) 
{
	vector<T1>		v(3);
	
	for ( long i = 0; i<3; ++i )
		v[i] = v3[i];
	
	return v;
}
*/
/*
template <typename T1, typename T2>
inline operator vector<T1>(vector<T2> v1) {
	vector<T1>		v(v1.size());
	
	for ( long i = 0; i<v1.size(); ++i )
		v[i] = T1(v1[i]);
	
	return v;
}
*/
/*
template <typename T1, typename T2>
inline void				operator=(vector<T1>& v1, vector<T2> v2)
{
	v1.reserve(v2.size());
	
	for ( long i = 0; i<v1.size() && i<v2.size(); ++i )
		v1[i] = T2(v2[i]);
}
*/
template <typename T>
inline vector<T>		operator-(vector<T> v1)
{
	vector<T>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = -v1[i];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		operator+(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v[i] = v1[i] + v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		operator+(vector<T1> v1, T2 v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = v1[i] + v2;
	
	return v;
}

template <typename T1, typename T2>
inline void				operator+=(vector<T1>& v1, vector<T2> v2)
{
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v1[i] += v2[i];
}

template <typename T1, typename T2>
inline void				operator+=(vector<T1>& v1, T2 v2)
{
	for ( size_t i = 0; i<v1.size(); ++i )
		v1[i] += v2;
}


template <typename T1, typename T2>
inline vector<T1>		operator-(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v[i] = v1[i] - v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		operator-(vector<T1> v1, Vector3<T2>& v2)
{
	for ( size_t i = 0; i<v1.size() && i<3; ++i )
		v1[i] -= v2[i];
	
	return v1;
}

template <typename T1, typename T2>
inline vector<T1>		operator-(vector<T1> v1, T2 v2)
{
	vector<T1>		v(v1.size());
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = v1[i] - v2;
	
	return v;
}

template <typename T1, typename T2>
inline void				operator-=(vector<T1>& v1, vector<T2> v2)
{
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v1[i] -= v2[i];
}

template <typename T1, typename T2>
inline void				operator-=(vector<T1>& v1, T2 v2)
{
	for ( size_t i = 0; i<v1.size(); ++i )
		v1[i] -= v2;
}

template <typename T1, typename T2>
inline vector<T1>		operator*(vector<T1> v1, T2 v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = v1[i] * v2;
	
	return v;
}

template <typename T1, typename T2>
inline void				operator*=(vector<T1>& v1, T2 v2)
{
	for ( size_t i = 0; i<v1.size(); ++i )
		v1[i] *= v2;
}

template <typename T1, typename T2>
inline vector<T1>		operator*(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v[i] = v1[i] * v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline void				operator*=(vector<T1>& v1, vector<T2> v2)
{
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v1[i] *= v2[i];
}

template <typename T1, typename T2>
inline vector<T1>		operator/(vector<T1> v1, T2 v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = v1[i] / v2;
	
	return v;
}

template <typename T1, typename T2>
inline void				operator/=(vector<T1>& v1, T2 v2)
{
	for ( size_t i = 0; i<v1.size(); ++i )
		v1[i] /= v2;
}

template <typename T1, typename T2>
inline vector<T1>		operator/(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v[i] = v1[i] / v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline void				operator/=(vector<T1> v1, vector<T2> v2)
{
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		v1[i] /= v2[i];
}

template <typename T>
inline double			length2(vector<T> v)
{
	double			len(0);
	
	for ( auto it = v.begin(); it != v.end(); ++it )
		len += *it * *it;
	
	return len;
}

template <typename T>
inline double			length(vector<T> v)
{
	double			len(length2(v));
	
	return sqrt(len);
}

template <typename T1, typename T2>
inline double			distance(vector<T1>& v1, vector<T2>& v2)
{
	vector<T1>		v = v1 - v2;
	
	return length(v);
}


template <typename T>
inline T				volume(vector<T> v)
{
	T				vol(1);
	
	for ( auto it = v.begin(); it != v.end(); ++it )
		vol *= *it;
	
	return vol;
}

template <typename T>
inline vector<T>		absolute(vector<T> v1)
{
	vector<T>		v(v1.size());
	
	for ( size_t i = 0; i<v1.size(); ++i )
		v[i] = fabs(v1[i]);

	return v;
}

template <typename T1, typename T2>
inline vector<T1>		vecmin(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		if ( v1[i] > v2[i] ) v[i] = v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		vecmin(vector<T1> v1, T2 v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		if ( v1[i] > v2 ) v[i] = v2;
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		vecmax(vector<T1> v1, vector<T2> v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size() && i<v2.size(); ++i )
		if ( v1[i] < v2[i] ) v[i] = v2[i];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		vecmax(vector<T1> v1, T2 v2)
{
	vector<T1>		v = v1;
	
	for ( size_t i = 0; i<v1.size(); ++i )
		if ( v1[i] < v2 ) v[i] = v2;
	
	return v;
}

template <typename T1, typename T2>
inline double			scalar(vector<T1>& v1, vector<T2>& v2)
{
	double			d(0);
	
	for ( size_t i = 0; i<v1.size(); ++i )
		d += v1[i] * v2[i];
	
	return d;
}

template <typename T1, typename T2>
inline vector<T1>		cross(vector<T1>& v1, vector<T2>& v2)
{
	vector<T1>		v(v1.size());
	
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2],
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
	
	return v;
}

template <typename T1, typename T2>
inline vector<T1>		angle(vector<T1>& v1, vector<T2>& v2)
{
	double		prod = scalar(v1, v2);
	if ( prod > 1 ) prod = 1;
	if ( prod < -1 ) prod = -1;
	return(acos(prod));
}

/**
@brief 	Divides a value by each of the elements of a 3-value vector.
@param 	d				value to be divided.
@param 	vec				3-value vector.
@return Vector3<Type>	new 3-value vector.

	All elements of the vector must be non-zero.

**/
template <typename Type>
inline Vector3<Type>	operator/(double d, Vector3<Type> vec)
{
	if ( vec[0] == 0 || vec[1] == 0 || vec[2] == 0 )
		cerr << "Error: Zero elements in vector " << vec << endl;
	
	return Vector3<Type>(d/vec[0], d/vec[1], d/vec[2]);
}

template <typename Type>
inline Vector3<Type>	vector3_scalar_range(Vector3<Type> vec, double min, double max)
{
	if ( max < min ) return(vec);
	
	if ( vec[0] < min ) vec[0] = min;
	if ( vec[0] > max ) vec[0] = max;
	if ( vec[1] < min ) vec[1] = min;
	if ( vec[1] > max ) vec[1] = max;
	if ( vec[2] < min ) vec[2] = min;
	if ( vec[2] > max ) vec[2] = max;
	
	return vec;
}

template <typename T1, typename T2>
Vector3<T1>	vector3_origin_to_shift(Vector3<T1> origin, Vector3<T2> size)
{
	Vector3<double>	h(size/2);
	
	while ( origin[0] < -h[0] ) origin[0] += size[0];
	while ( origin[1] < -h[1] ) origin[1] += size[1];
	while ( origin[2] < -h[2] ) origin[2] += size[2];
	while ( origin[0] >  h[0] ) origin[0] -= size[0];
	while ( origin[1] >  h[1] ) origin[1] -= size[1];
	while ( origin[2] >  h[2] ) origin[2] -= size[2];
	
	return origin;
}

template <typename T1, typename T2>
inline Vector3<T1>	vector3_set_PBC(Vector3<T1> coord, Vector3<T2> box)
{
	while ( coord[0] < 0 ) coord[0] += box[0];
	while ( coord[1] < 0 ) coord[1] += box[1];
	while ( coord[2] < 0 ) coord[2] += box[2];
	while ( coord[0] >= box[0] ) coord[0] -= box[0];
	while ( coord[1] >= box[1] ) coord[1] -= box[1];
	while ( coord[2] >= box[2] ) coord[2] -= box[2];
	
	return coord;
}

template <typename T1, typename T2, typename T3>
inline Vector3<T1>	vector3_difference_PBC(Vector3<T1> v1, Vector3<T2> v2, Vector3<T3> box)
{
	Vector3<T1>		d = v1 - v2;
	
	if ( box[0] - d[0] < d[0] ) d[0] -= box[0];
	if ( box[0] + d[0] < -d[0] ) d[0] += box[0];
	if ( box[1] - d[1] < d[1] ) d[1] -= box[1];
	if ( box[1] + d[1] < -d[1] ) d[1] += box[1];
	if ( box[2] - d[2] < d[2] ) d[2] -= box[2];
	if ( box[2] + d[2] < -d[2] ) d[2] += box[2];
	
	return d;
}


#endif


