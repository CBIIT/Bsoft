/**
@file	moving_average.h
@brief	Header file moving average calculations
@author 	Bernard Heymann
@date	Created: 20000430
@date	Modified: 20221129
**/

#include "Complex.h"
#include "matrix_linear.h"
#include "utilities.h"

// Function prototypes 
vector<double>	moving_average(long number, double* x, long window);
vector<double>	moving_average(vector<double>& x, long window);
vector<Complex<float>>	moving_average_complex(vector<Complex<float>>& x, long window);
vector<double>	moving_polynomial(long order, long number, double* x, long window);
vector<double>	moving_polynomial(long order, vector<double>& x, long window);
vector<double>	moving_gradient(vector<double>& x, long window);
vector<double>	moving_curvature(vector<double>& x, long window);

