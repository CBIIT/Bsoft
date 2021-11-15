/**
@file	Transform.cpp
@brief	Transform manipulation functions
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20150526
**/
 
#include "Transform.h" 
#include "Matrix.h"
#include "matrix_linear.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Solves a 3D geometry fit matrix.
@param 	a				nxn matrix.
@param 	bx				x-dimension n-vector.
@param 	by				y-dimension n-vector.
@param 	bz				z-dimension n-vector.
@param	flag			type of decomposition: 0=LUD, 1=SVD.
@return Transform		transform structure.

	The matrix is inverted and solutions for the three vectors generated.

**/
Transform	transform_matrix_solve(Matrix a, vector<double>& bx, 
				vector<double>& by, vector<double>& bz, int flag)
{
	long			n(a.rows());
	int				i, j, c;
	int				ft = (n>3)? 1: 0;
	Matrix3			mat;
	
	for ( i=0; i<n-1; i++ )
		for ( j=i+1; j<n; j++ ) a[i][j] = a[j][i];
	
	c = a.check_for_singularity();
	
	if ( verbose & VERB_FULL )
		cout << "Singularity check: " << c << endl;
	
	if ( verbose & VERB_FULL )
		cout << a;

	if ( c ) {
		n--;
		i = 0;
		j = c;
		while ( ( j = j >> 1 ) ) i++;	// I'm not sure this works correctly!
		if ( verbose & VERB_FULL )
			cout << "Deleting row " << i << endl;
		a = a.delete_row_column(i);
	}
	
	if ( verbose & VERB_FULL )
		cout << a;

	if ( flag ) a.singular_value_decomposition();
	else a.LU_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	a.multiply_in_place(bz);

	if ( verbose & VERB_FULL )
		cout << a;

	switch ( c ) {
		case 0:
			mat = Matrix3(bx[0],bx[1],bx[2],by[0],by[1],by[2],bz[0],bz[1],bz[2]);
			break;
		case 1:
			mat = Matrix3(0,0,1,0,by[1],by[2],0,bz[1],bz[2]);
			break;
		case 2:
			mat = Matrix3(bx[0],0,bx[2],0,1,0,bz[0],0,bz[2]);
			break;
		case 4:
			mat = Matrix3(bx[0],bx[1],0,by[0],by[1],0,0,0,1);
			break;
		default: break;
	}

//	cout << mat << endl;
	mat.normalize();
//	cout << mat << endl;
	
	Transform		t(mat);
	
	if ( ft ) t.trans = Vector3<double>(bx[3], by[3], bz[3]);
	
	return t;
}

/**
@brief	Determines the transformation between two sets of vectors.
@param 	n					number of vectors in each set.
@param 	*vector1			first vector set.
@param 	*vector2			second vector set.
@param 	shift				flag to indicate determining shift.
@return Transform			transform structure.

	The two sets must have the same number of vectors.

**/
Transform	transform_find(int n, Vector3<double>* vector1, Vector3<double>* vector2, int shift)
{
	long			m = 3 + shift;
	long			h, i, j;
//	double			bx[4], by[4], bz[4], v[4] = {0,0,0,1};
	vector<double>	bx(4), by(4), bz(4), v(4,0);
	Matrix			a(m,m);
	
	v[3] = 1;
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	for ( h=0; h<n; h++ ) {
		v[0] = vector1[h][0];
		v[1] = vector1[h][1];
		v[2] = vector1[h][2];
		for ( i=0; i<m; i++ ) {
			bx[i] += v[i]*vector2[h][0];
			by[i] += v[i]*vector2[h][1];
			bz[i] += v[i]*vector2[h][2];
			for ( j=0; j<=i; j++ ) a[i][j] += v[i]*v[j];
		}
	}
	
	return transform_matrix_solve(a, bx, by, bz, 0);
}

