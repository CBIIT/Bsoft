/**
@file	matrix_util.h
@brief	Utility functions for matrices.
@author Bernard Heymann
@date	Created: 20010723
@date	Modified: 20190911
**/

#include "Matrix.h"
#include "rwimg.h"

// Function prototypes
int			matrix_normalize(Matrix& m);
double		matrix_find_cutoff_for_number(Matrix m, int n);
int 		matrix_log_1_R(Matrix matrix);
int*		matrix_find_linear_sequence(Matrix, int window);
int			matrix_permute(int i, int j, Matrix, int* order, 
				double* best_R, int* best_order);
double		matrix_calc_dist_fit(Matrix, int* order);
//Bimage*		img_from_matrix(Matrix& matrix, long scale);

