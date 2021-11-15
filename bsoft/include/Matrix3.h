/**
@file	Matrix3.h
@brief	3x3 Matrix manipulation functions
@author	Bernard Heymann
@date	Created: 20000501
@date	Modified: 20210515
*/

#ifndef _Matrix3_
#define _Matrix3_

#include "Vector3.h"
#include "Quaternion.h"
#include "Matrix.h" 


class Matrix3
{
private:
	Vector3<double>	d[3];
	void	check() {
		for ( long i=0; i<3; i++ ) {
			for ( long j=0; j<3; j++ ) {
				if ( fabs(d[i][j]) < 1e-12 ) d[i][j] = 0;
				if ( d[i][j] > 1 ) d[i][j] = 1;
				if ( d[i][j] < -1 ) d[i][j] = -1;
			}
		} 
	}
	void		from_quaternion(Quaternion q) {
		d[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    	d[0][1] = 2*q[1]*q[2] - 2*q[0]*q[3];
    	d[0][2] = 2*q[1]*q[3] + 2*q[0]*q[2];
    	d[1][0] = 2*q[1]*q[2] + 2*q[0]*q[3];
    	d[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    	d[1][2] = 2*q[2]*q[3] - 2*q[0]*q[1];
    	d[2][0] = 2*q[1]*q[3] - 2*q[0]*q[2];
    	d[2][1] = 2*q[2]*q[3] + 2*q[0]*q[1];
    	d[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
		check();
//		for ( long i=0; i<3; i++ ) cout << d[i] << endl;
	}
public:
	Matrix3() { }
	// Diagonal matrix
	Matrix3(const double v) { for ( long i=0; i<3; i++ ) d[i][i] = v; }
	Matrix3(const double* v) {
		long		i, j, k;
		for ( i=k=0; i<3; i++ )
			for ( j=0; j<3; j++, k++ )
				d[i][j] = v[k];
	}
	Matrix3(const vector<double>& v) {
		long		i, j, k;
		for ( i=k=0; i<3; i++ )
			for ( j=0; j<3; j++, k++ )
				d[i][j] = v[k];
	}
	Matrix3(const double d0, const double d1, const double d2,
			const double d3, const double d4, const double d5,
			const double d6, const double d7, const double d8) {
		d[0][0] = d0; d[0][1] = d1; d[0][2] = d2;
		d[1][0] = d3; d[1][1] = d4; d[1][2] = d5;
		d[2][0] = d6; d[2][1] = d7; d[2][2] = d8;
	}
	Matrix3(const Matrix3& m) { for ( long i=0; i<3; i++ ) d[i] = m.d[i]; }
	// Matrix from a quaternion
	Matrix3(Quaternion q) {
		from_quaternion(q);
	}
	// Matrix from a rotation axis and angle
	template <typename T>
	Matrix3(Vector3<T> axis, double angle) {
		Quaternion		q(axis, angle);
		from_quaternion(q);
	}
	// Matrix from a tilt angle and axis
	Matrix3(double tilt, double axis) {
		double			tilt_2 = tilt/2.0L;
		double			sin_tilt = sin(tilt_2);
	 	Quaternion		q(cos(tilt_2), cos(axis)*sin_tilt, sin(axis)*sin_tilt, 0);
		from_quaternion(q);
	}
	template <typename Type1, typename Type2>
	Matrix3(Vector3<Type1>& from_vec, Vector3<Type2>& to_vec) {
		Quaternion	q(from_vec, to_vec);
		from_quaternion(q);
	}
/*	Quaternion	quaternion() const {
		Quaternion	q;
		double		t = 1 + data[0] + data[4] + data[8];
		if ( t > 1 ) {
			t = 0.5*sqrt(t);
			q.s(t);
			t = 0.25/t;
			q.x((data[7] - data[5])*t);
			q.y((data[2] - data[6])*t);
			q.z((data[3] - data[1])*t);
//			puts("** 1 **");
		} else if ( data[0] > data[4] && data[0] > data[8] ) {
			t = 0.5*sqrt(1 + data[0] - data[4] - data[8]);
			q.x(t);
			t = 0.25/t;
			q.s(-(data[5] - data[7])*t);
			q.y((data[1] + data[3])*t);
			q.z((data[2] + data[6])*t);
//			puts("** 2 **");
		} else if ( data[4] > data[8] ) {
			t = 0.5*sqrt(1 - data[0] + data[4] - data[8]);
			q.y(t);
			t = 0.25/t;
//			q.s(-(data[2] - data[6])*t);
			q.s((data[2] - data[6])*t);
			q.x((data[1] + data[3])*t);
			q.z((data[5] + data[7])*t);
//			puts("** 3 **");
		} else {
			t = 0.5*sqrt(1 - data[0] - data[4] + data[8]);
			q.z(t);
			t = 0.25/t;
			q.s(-(data[1] - data[3])*t);
			q.x((data[2] + data[6])*t);
			q.y((data[5] + data[7])*t);
//			puts("** 4 **");
		}
		q.normalize();
		return q;
	}*/
	Quaternion	quaternion() const {
		Quaternion	q;
		double		t = 1 + d[0][0] + d[1][1] + d[2][2];
		if ( t > 1 ) {
			t = 0.5*sqrt(t);
			q[0] = t;
			t = 0.25/t;
			q[1] = (d[2][1] - d[1][2])*t;
			q[2] = (d[0][2] - d[2][0])*t;
			q[3] = (d[1][0] - d[0][1])*t;
//			puts("** 1 **");
		} else if ( d[0][0] > d[1][1] && d[0][0] > d[2][2] ) {
			t = 0.5*sqrt(1 + d[0][0] - d[1][1] - d[2][2]);
			q[1] = t;
			t = 0.25/t;
			q[0] = -(d[1][2] - d[2][1])*t;
			q[2] = (d[0][1] + d[1][0])*t;
			q[3] = (d[0][2] + d[2][0])*t;
//			puts("** 2 **");
		} else if ( d[1][1] > d[2][2] ) {
			t = 0.5*sqrt(1 - d[0][0] + d[1][1] - d[2][2]);
			q[2] = t;
			t = 0.25/t;
//			q[0] = -(d[0][2] - d[2][0])*t;
			q[0] = (d[0][2] - d[2][0])*t;
			q[1] = (d[0][1] + d[1][0])*t;
			q[3] = (d[1][2] + d[2][1])*t;
//			puts("** 3 **");
		} else {
			t = 0.5*sqrt(1 - d[0][0] - d[1][1] + d[2][2]);
			q[3] = t;
			t = 0.25/t;
			q[0] = -(d[0][1] - d[1][0])*t;
			q[1] = (d[0][2] + d[2][0])*t;
			q[2] = (d[1][2] + d[2][1])*t;
//			puts("** 4 **");
		}
		q.normalize();
		return q;
	}
	Vector3<double>&	operator[](long i) { while ( i>2 ) i-=3; return d[i]; }
//	This asignment operator does not work!
//	Matrix3 operator=(const Matrix3& m) { Matrix3 mn(m); return mn; }
	Matrix3 operator=(const Matrix3& m) {
		for ( long i=0; i<3; i++ ) d[i] = m.d[i];
		return *this;
	}
	Matrix3 operator+=(const Matrix3& m) {
		for ( long i=0; i<3; i++ ) d[i] += m.d[i];
		return *this;
	}
	Matrix3 operator+(const Matrix3& m) {
		Matrix3		mn(*this);
		return mn += m;
	}
	Matrix3 operator-=(const Matrix3& m) {
		for ( int i=0; i<3; i++ ) d[i] -= m.d[i];
		return *this;
	}
	Matrix3 operator-(const Matrix3& m) {
		Matrix3		mn(*this);
		return mn -= m;
	}
	Matrix3 operator-() {
		Matrix3		mn;
		for ( int i=0; i<3; i++ ) mn.d[i] = -d[i];
		return mn;
	}
	Matrix3 operator*=(const Matrix3& m) {
		*this = *this * m;
		return *this;
	}
	Matrix3 operator*(const Matrix3& m) {
		long			i, j, k;
		Matrix3			mn;
		for ( i=0; i<3; i++ )
			for ( j=0; j<3; j++ )
				for ( k=0; k<3; k++ )
					mn.d[i][j] += d[i][k] * m.d[k][j];
		return mn;
	}
/*	Matrix3 operator*(const Matrix3& m) const {
		long			i, j;
		Matrix3			mt = m.transpose();
		Matrix3			mn;
		for ( i=0; i<3; i++ )
			for ( j=0; j<3; j++ )
				mn.d[i][j] = d[i].scalar(m.d[j]);
		return mn;
	}*/
	Matrix3 operator*=(const double v) {
		for ( long i=0; i<3; i++ ) d[i] *= v;
		return *this;
	}
	Matrix3 operator*(const double v) {
		Matrix3			mn;
		for ( long i=0; i<3; i++ ) mn.d[i] = d[i] * v;
		return mn;
	}
	Matrix3 operator/=(const double v) {
		double			div = 1/v;
		return *this *= div;
	}
	Matrix3 operator/(const double v) {
		double			div = 1/v;
		return *this * div;
	}
	bool	operator==(const Matrix3& m) {
		int		e = 0;
		for ( long i=0; i<3; i++ )
			e += ( d[i] == m.d[i] );
		return ( e == 3 );
	}
	template <typename T>
	Vector3<T>	operator*(Vector3<T>& v) {
		Vector3<T>		vp;
		for ( long i=0; i<3; i++ ) vp[i] = v.scalar(d[i]);
		return vp;
	}
	template <typename T>
	vector<T>	operator*(vector<T>& v) {
		vector<T>		vp(3);
		Vector3<T>		vo(v[0], v[1], v[2]);
		for ( long i=0; i<3; i++ ) vp[i] = vo.scalar(d[i]);
		return vp;
	}
	template <typename T>
	Matrix3	operator/=(Vector3<T>& v) {
		for ( long i=0; i<3; i++ ) d[i] /= v[i];
		return *this;
	}
	template <typename T>
	Matrix3	operator/(Vector3<T>& v) {
		Matrix3			mn(*this);
		return mn /= v;
	}
	double	determinant() {
		return d[0][0]*(d[1][1]*d[2][2]-d[1][2]*d[2][1])
				- d[0][1]*(d[1][0]*d[2][2]-d[1][2]*d[2][0])
				+ d[0][2]*(d[1][0]*d[2][1]-d[1][1]*d[2][0]);
	}
	double	trace() {
		return d[0][0] + d[1][1] + d[2][2];
	}
	Matrix3			transpose() const {
		Matrix3			m;
		for ( long i=0; i<3; i++ )
			for ( long j=0; j<3; j++ )
				m[i][j] = d[j][i];
		return m;
	}
/*	void	normalize() {
		long		i;
		double		s;
		for ( i=0, s=0; i<9; i+=3 ) s += d[i]*d[i];	// First column
		s = sqrt(s);
		for ( i=0; i<9; i+=3 ) d[i] /= s;
		for ( i=0, s=0; i<9; i+=3 ) s += d[i]*d[i+1];	// Second column
		for ( i=0; i<9; i+=3 ) d[i+1] -= s*d[i];
		for ( i=1, s=0; i<9; i+=3 ) s += d[i]*d[i];
		s = sqrt(s);
		for ( i=1; i<9; i+=3 ) d[i] /= s;
		for ( i=0, s=0; i<9; i+=3 ) s += d[i]*d[i+2];	// Third column
		for ( i=0; i<9; i+=3 ) d[i+2] -= s*d[i];
		for ( i=1, s=0; i<9; i+=3 ) s += d[i]*d[i+1];
		for ( i=1; i<9; i+=3 ) d[i+1] -= s*d[i];
		for ( i=2, s=0; i<9; i+=3 ) s += d[i]*d[i];
		s = sqrt(s);
		for ( i=2; i<9; i+=3 ) d[i] /= s;
	}*/
	void	normalize() {
		long		i;
		double		s;
		for ( i=0, s=0; i<3; i++ ) s += d[i][0]*d[i][0];	// First column
		s = sqrt(s);
		for ( i=0; i<3; i++ ) d[i][0] /= s;
		for ( i=0, s=0; i<3; i++ ) s += d[i][0]*d[i][1];	// Second column
		for ( i=0; i<3; i++ ) d[i][1] -= s*d[i][0];
		for ( i=0, s=0; i<3; i++ ) s += d[i][1]*d[i][1];
		s = sqrt(s);
		for ( i=0; i<3; i++ ) d[i][1] /= s;
		for ( i=0, s=0; i<3; i++ ) s += d[i][0]*d[i][2];	// Third column
		for ( i=0; i<3; i++ ) d[i][2] -= s*d[i][0];
		for ( i=0, s=0; i<3; i++ ) s += d[i][1]*d[i][2];
		for ( i=0; i<3; i++ ) d[i][2] -= s*d[i][1];
		for ( i=0, s=0; i<3; i++ ) s += d[i][2]*d[i][2];
		s = sqrt(s);
		for ( i=0; i<3; i++ ) d[i][2] /= s;
	}
	// Tilt axis, tilt angle, level angle
	vector<double> tilt_angles() {
		Quaternion		q = quaternion();
		vector<double>	v;
		v.push_back(atan2(q[2], q[1]));
		v.push_back(2*acos(q[0]));
		v.push_back(2*asin(q[3]));
		return v;
	}
	Matrix3		rotation(const Matrix3& mat) {
		Quaternion	q1 = quaternion();
		Quaternion	q2 = mat.quaternion();
		Quaternion	q = q1 * q2.inverse();
		return Matrix3(q);
	}
	double		angle(const Matrix3& mat) {
		Quaternion	q1 = quaternion();
		Quaternion	q2 = mat.quaternion();
		Quaternion	q = q1 * q2.inverse();
		return q.angle();
	}
	Vector3<double>	plane_normal(vector<double>& b) {
		Vector3<double>		normal;
		if ( d[0][0] == 0 ) {
			normal = Vector3<double>(1, 0, 0);
		} else if ( d[1][1] == 0 ) {
			normal = Vector3<double>(0, 1, 0);
		} else if ( d[2][2] == 0 ) {
			normal = Vector3<double>(0, 0, 1);
		} else {
			LU_decomposition(b);
			normal = Vector3<double>(b[0], b[1], b[2]);
			normal.normalize();
		}
		return normal;
	}
	Vector3<double>	plane_normal(Vector3<double>& b) {
		vector<double>	bb = {b[0],b[1],b[2]};
		Vector3<double>	n = plane_normal(bb);
		b = Vector3<double>(bb[0],bb[1],bb[2]);
		return n;
	}
	Matrix3		LU_decomposition(vector<double>& b) {
		long		i, j;
		Matrix		mat(3,3);
		Matrix3		mat3;
		for ( i=0; i<3; i++ ) for ( j=0; j<3; j++ ) mat[i][j] = d[i][j];
		mat.LU_decomposition(b);
		for ( i=0; i<3; i++ ) for ( j=0; j<3; j++ ) mat3[i][j] = mat[i][j];
		return mat3;
	}
	Matrix3		singular_value_decomposition() {
		long		i, j;
		Matrix		mat(3,3);
		Matrix3		mat3;
		for ( i=0; i<3; i++ ) for ( j=0; j<3; j++ ) mat[i][j] = d[i][j];
		mat.singular_value_decomposition();
		for ( i=0; i<3; i++ ) for ( j=0; j<3; j++ ) mat3[i][j] = mat[i][j];
		return mat3;
	}
};

template <typename Type>
inline Matrix3		operator*(Vector3<Type>& v, Matrix3& m)
{
	Matrix3			mn;
	for ( long i=0; i<3; i++ )
		mn[i] = m[i] * v[i];
	
	return mn;
}

/*template <typename Type>
inline Matrix3		operator/=(Matrix3& m, Vector3<Type>& v)
{
//	Type*			d = (Type *) &v;
	for ( long i=0; i<3; i++ )
			m[i] /= v[i];

	return m;
}
*/
/*template <typename Type>
inline Matrix3		operator/(const Matrix3& m, Vector3<Type>& v)
{
	Matrix3			mn(m);
	return mn /= v;
}
*/
#endif

ostream& operator<<(ostream& output, Matrix3& m);

// Function prototypes
int 		matrix3_show_hp(Matrix3& mat);


