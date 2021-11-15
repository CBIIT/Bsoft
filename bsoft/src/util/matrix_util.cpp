/**
@file	matrix_util.cpp
@brief	Utility functions for matrices.
@author Bernard Heymann
@date	Created: 20010723
@date	Modified: 20150124
**/

#include "matrix_util.h"
#include "Matrix.h"
#include "matrix_linear.h"
#include "math_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Normalizes a matrix.
@param 	&m		matrix.
@return int				0.

	The rows and columns are alternatively iteratively normalized until
	the error is small enough.

**/
int			matrix_normalize(Matrix& m)
{
	long			i, r, c;
	double			err(1);
	double*			rw = new double[m.rows()];
	double*			cw = new double[m.columns()];
	
	if ( verbose & VERB_FULL )
		cout << "Cycle\tError" << endl;
	for ( i=0; i<100 && err > 1e-20; i++ ) {
		for ( r=0; r<m.rows(); r++ ) rw[r] = 0;	// Row scaling
		for ( r=0; r<m.rows(); r++ )
			for ( c=0; c<m.columns(); c++ )
				rw[r] += m[r][c];
		for ( r=0; r<m.rows(); r++ ) rw[r] = 1/rw[r];
		for ( r=0; r<m.rows(); r++ )
			for ( c=0; c<m.columns(); c++ )
				m[r][c] *= rw[r];

		for ( c=0; c<m.columns(); c++ ) cw[c] = 0;	// Column scaling
		for ( r=0; r<m.rows(); r++ )
			for ( c=0; c<m.columns(); c++ )
				cw[c] += m[r][c];
		for ( c=0; c<m.columns(); c++ ) cw[c] = 1/cw[c];
		for ( r=0; r<m.rows(); r++ )
			for ( c=0; c<m.columns(); c++ )
				m[r][c] *= cw[c];

		for ( err=0, r=0; r<m.rows(); r++ ) err += rw[r];
		for ( c=0; c<m.columns(); c++ ) err += cw[c];
		err = fabs(err/(m.rows() + m.columns()) - 1);
		if ( verbose & VERB_FULL )
			cout << i+1 << tab << err << endl;
	}

	delete[] rw;
	delete[] cw;
	
	return 0;
}

/**
@brief 	Calculates a cutoff to ensure that the given number of elements fall below it.
@param 	m			pairwise square matrix.
@param 	n				number of elements below cutoff.
@return double 				cutoff.

	Each column is scanned and the cutoff adjusted so that at least the
	given number of elements will fall below it.

**/
double		matrix_find_cutoff_for_number(Matrix m, int n)
{
	n=2;
	long			j, k, c, r;
	double			dcut(0);
	double*			d = new double[n];
	
	for ( c=0; c<m.columns(); c++ ) {
		for ( j=0; j<n; j++ ) d[j] = 1e30;
		for ( r=0; r<m.rows(); r++ ) if ( c != r ) {
			for ( j=0; j<n && d[j] < m[r][c]; j++ ) ;
			if ( j < n ) {
				for ( k=n-1; k>j; k-- ) d[k] = d[k-1];
				d[j] = m[r][c];
			}
		}
		if ( dcut < d[n-1] ) dcut = d[n-1];
	}
	
	delete[] d;
	
	return dcut;
}

/**
@brief 	Calculates the logarithmic derivative of an R-factor in a matrix.
@param 	matrix		pairwise square matrix.
@return int 				0.
**/
int 		matrix_log_1_R(Matrix matrix)
{
	long 		i, j;
	
	for ( i=0; i<matrix.rows(); i++ ) {
		for ( j=0; j<matrix.columns(); j++ ) {
			if ( matrix[i][j] < 1 )
				matrix[i][j] = -log(1 - matrix[i][j]);
			else
				matrix[i][j] = 1e30;
		}
	}
	
	return 0;
}

/**
@brief 	Finds a linear sequence of elements based on a pairwise matrix.
@param 	matrix		pairwise square matrix.
@param 	window			range of elements to permute.
@return int*				the best order, NULL on error.
**/
int*		matrix_find_linear_sequence(Matrix matrix, int window)
{
	if ( matrix.columns() != matrix.rows() ) {
		cerr << "Error: the matrix is not square! (" << matrix.columns() << "x" << matrix.rows() << ")" << endl;
		return NULL;
	}
	
	int 		i, j, n = matrix.columns();
	double		best_R = 1e30;
	int*		order = new int[n];
	int*		best_order = new int[n];
	
	if ( window < 2 ) window = 2;
	if ( window > n ) window = n;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Finding the best sequence from " << n << " entities:" << endl;
		cout << "Window:                         " << window << endl;
		cout << "Number of permutations:         " << (n-window+1)*factorial(window) << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Finding the best sequence from " << n << " entities with a window of " << window << endl;
	
	for ( j=0; j<n; j++ ) best_order[j] = j; 				// Set up the initial order
	for ( i=0; i<n-window+1; i++ ) {
		for ( j=0; j<n; j++ ) order[j] = best_order[j]; 	// Keep the best order
		matrix_permute(i+1, i+window, matrix, order, &best_R, best_order);
	}

	if ( verbose & VERB_RESULT ) {
		cout << "Best R = " << best_R << endl << "Best order:";
		for ( j=0; j<n; j++ ) cout << " " << best_order[j]+1;
		cout << endl;
	}
		
	delete[] order;
	
	return best_order;
}

/**
@brief 	Permutes the order of a pairwise matrix and calculate a distance criterion.
@param 	i				start of window within the order array.
@param 	j				end of window within the order array.
@param 	matrix		pairwise square matrix.
@param 	*order			n-value array holding the current order.
@param 	*best_R		R-factor for the best order.
@param 	*best_order 	n-value array holding the best order.
@return int 				0.

	A recursive function to pass through all possible permutations
	within a given window.

**/
int			matrix_permute(int i, int j, Matrix matrix, int* order, 
				double* best_R, int* best_order)
{
	long 		k, m, n(matrix.columns());
	int* 		this_order;
	double		R;
  
	if( i == j ) {
		R = matrix_calc_dist_fit(matrix, order);
		if ( *best_R > R ) {
			*best_R = R;
			for ( m=0; m<n; m++ ) best_order[m] = order[m];
		}
	} else {
		this_order = new int[n];
		for ( m=0; m<n; m++ ) this_order[m] = order[m];
		for( k=i; k<=j; k++) {
			matrix_permute(i+1, j, matrix, order, best_R, best_order);
			if( k != j ) {
				swap(this_order[i-1], this_order[k]);
				for ( m=0; m<n; m++ ) order[m] = this_order[m];
			}
		}
		delete[] this_order;
	}
		
	return 0;
}


/**
@brief 	Calculates an R-factor for a specific order of elements in a pairwise matrix.
@param 	matrix		pairwise square matrix.
@param 	*order			n-value array holding the order.
@return double 				the R-factor.

	A new matrix is composed with the lower triangle taken from the input
	matrix in the given order, and the upper triangle calculated as the 
	sum of the diagonal elements between the indices of the element. 
	I.e, if i and j are the indices of an upper triangle element, then 
	it is the sum of all the diagonal elements from i to j:
		up(i,j) = sum(mat(k,k+1)) for i<=k<j
	The difference between the lower and upper triangles is then calculated
	as:
		          1      up(i,j) - upmin   lo(j,i) - lomin
		R = sum(----- * (--------------- - ---------------) )
		        i - j     upmax - upmin     lomax - lomin
	where upmin, upmax, lomin, and lomax are the extremes in each triangle.

**/
double		matrix_calc_dist_fit(Matrix matrix, int* order)
{
	long 		i, j, k, n = matrix.columns();
	double		min(1), max(0), new_min(1), new_max(0);
	Matrix		new_mat(n,n);
	
	// Set up the lower triangle with the original distances
	for ( i=1; i<n; i++ ) {
		for ( j=0; j<i; j++ ) {
			new_mat[i][j] = matrix[order[i]][order[j]];
			if ( min > new_mat[i][j] ) min = new_mat[i][j];
			if ( max < new_mat[i][j] ) max = new_mat[i][j];
		}
	}
	
	// Set up the second upper diagonal
	for ( j=1; j<n; j++ )
		new_mat[j-1][j] = new_mat[j][j-1];
	
	// Calculate the other upper diagonals
	for ( i=0; i<n-1; i++ ) {
		for ( j=i+2; j<n; j++ ) {
			for ( k=i+1; k<=j; k++ )
				new_mat[i][j] += new_mat[k][k-1];
			if ( new_min > new_mat[i][j] ) new_min = new_mat[i][j];
			if ( new_max < new_mat[i][j] ) new_max = new_mat[i][j];
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG matrix_calc_dist_fit: Order:" << endl;
		for ( i=0; i<n; i++ ) cout << " " << order[i];
		cout << endl << "DEBUG matrix_calc_dist_fit: New matrix:" << endl;
		cout << new_mat << endl;
	}
	
	// Calculate a figure-of-merit weighted by diagonal
	double		diff, weight, wsum(0), R(0);
	double		scale = 1/(max - min);
	double		new_scale = 1/(new_max - new_min);
	double		offset = min/(max - min) - new_min/(new_max - new_min);
	for ( i=2; i<n; i++ ) {
		for ( j=0; j<i-1; j++ ) {
			diff = scale*new_mat[i][j] - new_scale*new_mat[j][i] - offset;
			weight = 1.0/(i - j);
			R += weight*diff*diff;
			wsum += weight;
		}
	}
	R = sqrt(R/wsum);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "R = " << R << endl;
		for ( j=0; j<n; j++ ) cout << " " << order[j]+1;
		cout << endl;
	}
	
	return R;
}

/**
@brief 	Converts a matrix to an image.
@param 	matrix		matrix.
@param 	scale		integer scale for enlarging the image.
@return Bimage* 	new image.
**/
/*Bimage*		img_from_matrix(Matrix& matrix, long scale)
{
	if ( scale < 1 ) scale = 1;
	if ( scale > 10 ) scale = 10;
	
	Bimage*		pmat = new Bimage(Float, TSimple, scale*matrix.columns(), scale*matrix.rows(), 1, 1);
	
	long		i, j, k, x, y;
	
	for ( k=y=0; y<pmat->sizeY(); ++y ) {
		i = y/scale;
		for ( x=0; x<pmat->sizeX(); ++x, ++k ) {
			j = x/scale;
			pmat->set(k, matrix[i][j]);
		}
	}
	
	return pmat;
}
*/

