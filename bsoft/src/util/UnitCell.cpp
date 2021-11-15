/**
@file	UnitCell.cpp
@brief	Unit cell functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20150111
**/

#include "UnitCell.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the skew matrix from the unit cell parameters.
@return Matrix3 				the 3x3 skew matrix.

	Derived from the X-plor source rotate.s.
	New coordinates are obtained by r'(i)=sum_j matrix(i,j)*r(j)
	The convention to setup the matrices is as follows:
		a lies on the x-axis, 
		b lies in the x,y plane. 

**/
Matrix3		UnitCell::skew_matrix()
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG UnitCell::skew_matrix: Calculating the fractionalization matrix" << endl;
	
    if ( alpha() < 0 ) alpha(-alpha());
    if ( beta() < 0 ) beta(-beta());
    if ( gamma() < 0 ) gamma(-gamma());
	if ( alpha() > M_PI || beta() > M_PI || gamma() > M_PI )
		degrees_to_radians();
	
	Matrix3		skew;
	
    double   	cterm = (cos(beta())*cos(gamma())-
							cos(alpha()))/
    	    	    	(sin(beta())*sin(gamma()));
    double   	sterm = sqrt(1.0 - cterm*cterm);
    
    skew[0][0] = 1.0/a();
    skew[0][1] = -cos(gamma())/(sin(gamma())*a());
    skew[0][2] = -(cos(gamma())*sin(beta())*cterm+
				cos(beta())*sin(gamma()))/
     			(sin(beta())*sterm*sin(gamma())*a());
    skew[1][1] = 1.0/(sin(gamma())*b());
    skew[1][2] = cterm/(sterm*sin(gamma())*b());
    skew[2][2] = 1.0/(sin(beta())*sterm*c());
	
	for ( size_t i=0; i<3; i++ )
		for ( size_t j=0; j<3; j++ )
			if ( skew[i][j] < 1e-6 ) skew[i][j] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG UnitCell::skew_matrix: Fractionalization matrix: " << endl;
		cout << skew << endl;
	}
	
    return skew;
}

Matrix3		UnitCell::skew_rotation(int invert)
{
	Matrix3		skew = skew_matrix();
	Matrix3		mat(1);
	
	mat[0][0] = skew[0][0]*a();
	mat[0][1] = skew[0][1]*a();
	mat[0][2] = skew[0][2]*a();
	mat[1][1] = skew[1][1]*b();
	mat[1][2] = skew[1][2]*b();
	mat[2][2] = skew[2][2]*c();
	
	if ( invert ) { 	// Invert the skew matrix
		mat[0][0] = 1/mat[0][0];
		mat[1][1] = 1/mat[1][1];
		mat[2][2] = 1/mat[2][2];
		mat[0][1] = -mat[0][1]*mat[0][0]*mat[1][1];
		mat[1][2] = -mat[1][2]*mat[1][1]*mat[2][2];
		mat[0][2] = mat[0][1]*mat[1][2]/mat[1][1] - mat[0][0]*mat[2][2]*mat[0][2];
	}
	
	return mat;
}
