/**
@file	Matrix3.cpp
@brief	Matrix manipulation functions
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20150111
**/
 
#include "Matrix3.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, Matrix3& m) {
	output.setf(ios::fixed, ios::floatfield);
//	output.precision(4);
	for ( long i=0; i<3; i++ )
		output << m[i] << endl;
	return output;
}

/**
@brief 	Print out a matrix.
@param 	mat			3x3 matrix.
@return int 			0.

	Matrix elements less than 1e-6 are set to zero for display.

**/
int 		matrix3_show_hp(Matrix3& mat)
{
	cout << "Matrix:" << endl;
	cout << mat[0][0] << "," << mat[0][1] << "," << mat[0][2] << " (" << 
		sqrt(mat[0][0]*mat[0][0]+mat[0][1]*mat[0][1]+mat[0][2]*mat[0][2]) << ")" << endl;
	cout << mat[1][0] << "," << mat[1][1] << "," << mat[1][2] << " (" << 
		sqrt(mat[1][0]*mat[1][0]+mat[1][1]*mat[1][1]+mat[1][2]*mat[1][2]) << ")" << endl;
	cout << mat[2][0] << "," << mat[2][1] << "," << mat[2][2] << " (" << 
		sqrt(mat[2][0]*mat[2][0]+mat[2][1]*mat[2][1]+mat[2][2]*mat[2][2]) << ")" << endl;
	cout << " (" << sqrt(mat[0][0]*mat[0][0]+mat[1][0]*mat[1][0]+mat[2][0]*mat[2][0]) << "," << 
		sqrt(mat[0][1]*mat[0][1]+mat[1][1]*mat[1][1]+mat[2][1]*mat[2][1]) << "," << 
		sqrt(mat[0][2]*mat[0][2]+mat[1][2]*mat[1][2]+mat[2][2]*mat[2][2]) << ")" << endl;
	
	cout << "Determinant = " << mat.determinant() << endl << endl;
	
	return 0;
}

