/**
@file	dynamic_programming.h
@brief	Library functions for dynamic programming.
@author Bernard Heymann
@date	Created: 20050622
@date	Modified: 20110411
**/

#include "Matrix.h"

// Function prototypes
long		dp_matrix_scoring(Matrix mat, double gapopen, double gapextend);
int*		dp_matrix_backtrace(Matrix mat, double gapopen, double gapextend, long& length);

