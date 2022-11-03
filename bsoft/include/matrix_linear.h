/**
@file	matrix_linear.h 
@brief	Header file for matrix functions 
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20220216
**/

#include "Matrix.h"

// Function prototypes
double		linear_least_squares(int n1, int n2, double *x, double *y, double& a, double& b);
double		linear_least_squares(int n1, int n2, vector<double>& x, vector<double>& y, double& a, double& b);
double		fit_polynomial(int n, double* x, double* y, int order, double* coeff);
double		fit_polynomial(int n, vector<double>& x, vector<double>& y, int order, vector<double>& coeff);
Vector3<double>	fit_plane(Matrix a, vector<double> b);
vector<double>	fit_conic(vector<Vector3<double>> v);
vector<double>	fit_ellipse(vector<Vector3<double>> v);

