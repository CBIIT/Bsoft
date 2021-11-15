/**
@file	dynamic_programming.cpp
@brief	Library functions for dynamic programming.
@author Bernard Heymann
@date	Created: 20050622
@date	Modified: 20110810
**/

#include "dynamic_programming.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the scoring matrix in dynamic programming.
@param 	mat			matrix.
@param 	gapopen		gap opening penalty.
@param 	gapextend	gap extension penalty.
@return long 			index of the maximum in the matrix.

	Implementation of the Needleman-Wunsch algorithm.
	The input matrix can be any form of similarity or coincidence matrix, and
	is modified by dynamic programming in preparation for backtracing.

**/
long		dp_matrix_scoring(Matrix mat, double gapopen, double gapextend)
{
	int			m = mat.columns(), n = mat.rows();
	int			t, pt(0);
	long		i, j, k, kmax(0);
	double		value, max, maxall(0);
	
	for ( k=j=0; j<n; j++ ) {
		for ( i=0; i<m; i++, k++ ) {
			if ( i==0 || j==0 ) max = mat[j][i];
			else max = mat[j][i] + mat[j-1][i-1];
			value = 0;
			t = 0;
			if ( i > 0 ) {
				if ( pt == 1 ) value = mat[j][i-1] - gapextend;
				else value = mat[j][i-1] - gapopen;
			}
			if ( max < value ) {
				max = value;
				t = 1;
			}
			if ( j > 0 ) {
				if ( pt == 2 ) value = mat[j-1][i] - gapextend;
				else value = mat[j-1][i] - gapopen;
			}
			if ( max < value ) {
				max = value;
				t = 2;
			}
			mat[j][i] = max;
			if ( maxall < max ) {
				maxall = max;
				kmax = k;
			}
			pt = t;
		}
	}
	
	return kmax;
}

/**
@brief 	Backtraces the scoring matrix in dynamic programming.
@param 	mat			matrix.
@param 	gapopen		gap opening penalty.
@param 	gapextend	gap extension penalty.
@param 	*length		length of alignment.
@return int* 			alignment array with indices.

	Implementation of the Needleman-Wunsch algorithm.
	The scoring matrix for dynamic programming is backtraced into an integer 
	array holding the indices from the corresponding sequences, or -1 for gaps. 
	The array is double the length of the alignment and holds two sets of 
	indices, the first in the first half and the second in last half.
	The length of the alignment is returned in the pointer of the last argument.

**/
int*		dp_matrix_backtrace(Matrix mat, double gapopen, double gapextend, long& length)
{
	int			m = mat.columns(), n = mat.rows();
	int			t, pt(0);
	long		i, j, k, h, hbeg(0), hend, imax(0), jmax(0);
	long		tlen = m + n;
	double		max(0), tmax(0), value;
	int*		aln = new int[2*tlen];
	
	for ( h=0; h<tlen*2; h++ ) aln[h] = -1;
	
	for ( k=j=0; j<n; j++ ) {
		for ( i=0; i<m; i++, k++ ) {
			if ( tmax < mat[j][i] ) {
				tmax = mat[j][i];
				imax = i;
				jmax = j;
			}
		}
	}
	
	// Fill in the trailing pieces
	hend = m - imax;
	if ( n - jmax > hend ) hend = n - jmax;
	hend = tlen - hend;
	for ( h=hend, i=imax; h<tlen && i<m; h++, i++ ) aln[h] = i;
	for ( h=hend + tlen, j=jmax; h<2*tlen && j<n; h++, j++ ) aln[h] = j;
	
	// Trace back to get the central aligned pieces
	for ( h=hend, i=imax, j=jmax; i>0 && j>0; h-- ) {
		k = m*j + i;
		mat[j][i] = 2*tmax;
		t = 0;
		max = mat[j-1][i-1];
		if ( pt == 1 ) value = mat[j][i-1] - gapextend;
		else value = mat[j][i-1] - gapopen;
		if ( max < value ) {
			max = value;
			t = 1;
		}
		if ( pt == 2 ) value = mat[j-1][i] - gapextend;
		else value = mat[j-1][i] - gapopen;
		if ( max < value ) {
			max = value;
			t = 2;
		}
		if ( t < 2 ) {
			aln[h] = i;
			i--;
		}
		if ( t%2 == 0 ) {
			aln[h+tlen] = j;
			j--;
		}
		pt = t;
	}
	k = m*j + i;
	mat[j][i] = 2*tmax;
	
	// Fill in the leading pieces
	for ( ; h >= 0; i--, j--, h-- ) {
		if ( i >= 0 ) aln[h] = i;
		if ( j >= 0 ) aln[h+tlen] = j;
		if ( i >= 0 || j >= 0 ) hbeg = h;
	}
	
	// Generate the return array
	length = tlen - hbeg;
	int*		align = new int[2*length];
	
	for ( i=0, h=hbeg; i<length; i++, h++ ) {
		align[i] = aln[h];
		align[i+length] = aln[h+tlen];
	}

	delete[] aln;
	
	return align;
}
