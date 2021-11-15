/**
@file	Transform.h 
@brief	Class for a generalized transform object 
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20150526
**/

#ifndef _Transform_
#define _Transform_

#include "Vector3.h"
#include "Quaternion.h"
#include "Matrix3.h"
#include "Matrix.h"

/************************************************************************
@Object: class Transform
@Description:
	Transformation parameter structure.
@Features:
	This defines a full affine transformation.
*************************************************************************/
class Transform {
public:
	Transform*		next;		// Next in linked list
	Vector3<double>	origin;		// Rotation origin
	Vector3<double>	axis;		// Rotation axis
	double			angle;		// Rotation angle
	Vector3<double>	trans;		// Translation
	Vector3<double>	scale;		// Scale
	double			fom;		// A figure-of-merit associated with the transform
private:
	void		from_quaternion(Quaternion q) {
		next = NULL; origin=0; scale = Vector3<double>(1, 1, 1); trans = 0;
		axis = q.axis();
		axis.normalize();
		angle = 2*acos(q[0]);
		fom = 0;
	}
public:
	Transform() { next=NULL; origin=0; scale = Vector3<double>(1, 1, 1);
		trans = 0; axis = Vector3<double>(0, 0, 1); angle = 0; fom = 0; }
	Transform(Quaternion& q) {
		from_quaternion(q);
	}
	Transform(Matrix3& m) {
		Quaternion		q = m.quaternion();
		from_quaternion(q);
	}
	Vector3<double>	operator*(Vector3<double> v) {
		Matrix3		m(Matrix3(axis, angle));
		v -= origin;
		v = m * v;
		v += origin + trans;
		return v;
	}
} ;
#endif

// Function prototypes 
Transform	transform_matrix_solve(Matrix a, vector<double>& bx, 
				vector<double>& by, vector<double>& bz, int flag);
Transform	transform_find(int n, Vector3<double>* vector1, Vector3<double>* vector2, int shift);

