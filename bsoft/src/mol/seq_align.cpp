/**
@file	seq_align.cpp
@brief	Library routines to generate and analyze dot plots
@author Bernard Heymann
@date	Created: 20001029
@date	Modified: 20180622
**/

#include "rwmolecule.h"
#include "rwresprop.h"
#include "dynamic_programming.h"
#include "seq_align.h"
#include "seq_util.h"
#include "linked_list.h"
#include "Matrix3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Aligns two sequences.
@param 	*mol1			first molecule structure.
@param 	*mol2			second molecule structure.
@param 	gapopen			gap opening penalty.
@param 	gapextend			gap extension penalty.
@param 	*simat	residue similarity matrix.
@return long					length of the aligned sequences.

	Two sequences are aligned and returned in the sequence strings
	of the molecule structures.

**/
long		seq_pair_align(Bmolecule* mol1, Bmolecule* mol2, double gapopen, 
				double gapextend, Bresidue_matrix* simat)
{	
	if ( verbose & VERB_PROCESS ) {
		cout << "Aligning two sequences:" << endl;
		cout << mol1->id << " with " << mol1->nres << " residues" << endl;
		cout << mol2->id << " with " << mol2->nres << " residues" << endl;
		cout << "Gap opening penalty:            " << gapopen << endl; 
		cout << "Gap extension penalty:          " << gapextend << endl; 
		cout << endl;
	}
	
	Matrix			dp = seq_dot_plot(mol1, mol2, simat);
	
	long			imax = dp_matrix_scoring(dp, gapopen, gapextend);
	
	long			i(imax%dp.rows()), j(imax/dp.rows());
	long			alnlen(0);
	double			score = dp[j][i];
	
	int*			indaln = dp_matrix_backtrace(dp, gapopen, gapextend, alnlen);

//	Bstring			alnseq1(alnlen, '-');
//	Bstring			alnseq2(alnlen, '-');
	Bstring			alnseq1('-', alnlen);
	Bstring			alnseq2('-', alnlen);

	for ( i=0, j=alnlen; i<alnlen; i++, j++ ) {
		if ( indaln[i] > -1 ) alnseq1[i] = mol1->seq[indaln[i]];
		if ( indaln[j] > -1 ) alnseq2[i] = mol2->seq[indaln[j]];
	}
	
	delete[] indaln;
	
	mol1->seq = alnseq1;
	mol2->seq = alnseq2;

	if ( verbose & VERB_PROCESS ) {
		cout << "Alignment length:               " << alnlen << endl;
		cout << "Alignment score:                " << score << endl << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_pair_align:"  << endl << alnseq1 << endl << alnseq2 << endl;
	
	return alnlen;
}

/**
@brief 	Finds the best offset for aligning two sequences without gaps.
@param 	*mol1		first molecule structure.
@param 	*mol2		second molecule structure.
@param 	*nres		number of similar residues in best diagonal.
@param 	*simat		residue similarity matrix.
@return int			the offset of sequence 2 with respect to sequence 1.

	The best offset is determined as the diagonal of a dot plot with
	the highest value.

**/
int			seq_find_best_offset(Bmolecule* mol1, Bmolecule* mol2, 
				long& nres, Bresidue_matrix* simat)
{
	nres = 0;
	if ( mol1->nres < 10 || mol2->nres < 10 ) return 0;
	
	int				offset(0);
	long   			k, l;
	Matrix			dot = seq_dot_plot(mol1, mol2, simat);
	vector<double>	diag(mol1->nres+mol2->nres, 0);
	
	for ( k=0; k<mol1->nres; k++ ) {
		for ( l=0; l<mol2->nres; l++ ) {
			if ( dot[l][k] > 0 ) diag[mol2->nres-l-1+k] += 1;
		}
	}
	
	for ( k=0; k<mol1->nres+mol2->nres; k++ ) {
		if ( nres < diag[k] ) {
			nres = (long) diag[k];
			offset = k - mol2->nres + 1;
		}
	}
	
	if ( verbose & VERB_FULL ) {
		cout << mol1->id << tab << mol2->id << ":\t" << offset << tab << nres << endl;
		if ( offset > 0 ) cout << &mol1->seq[offset] << endl;
		else cout << mol1->seq << endl;
		if ( offset < 0 ) cout << &mol2->seq[-offset] << endl;
		else cout << mol2->seq << endl;
		cout << endl;
	}
	
	return offset;
}

/**
@brief 	Calculates a dot plot for two sequences.
@param 	*mol1			first molecule structure.
@param 	*mol2			second molecule structure.
@param 	*simat	residue similarity matrix.
@return Matrix					matrix with dot plot.
**/
Matrix		seq_dot_plot(Bmolecule* mol1, Bmolecule* mol2, Bresidue_matrix* simat)
{
	if ( verbose & VERB_FULL )
	    cout << "Dot plot between sequences of sizes " << mol1->nres << " and " << mol2->nres << endl;
	
	long			nres = simat->n;
	Bstring			code1 = simat->c;
	float*			sim_mat = simat->m;
    
    int				k, l, ii, jj;
	
	Matrix			dot_plot(mol1->nres, mol2->nres);

	for ( k=0; k<mol1->nres; k++ ) {
		ii = code1.index(mol1->seq[k]);
		for ( l=0; l<mol2->nres; l++ ) {
			jj = code1.index(mol2->seq[l]);
			if ( ii > -1 && jj > -1 )
	        	dot_plot[l][k] = sim_mat[ii*nres+jj];
        }
	}

	return dot_plot;
}

/**
@brief 	Calculates a moving average along the diagonals of a dot plot.
@param 	dot_plot	dot plot matrix.
@param 	window			moving average window size.
@return Matrix				matrix with moving average dot plot.
**/
Matrix		seq_dot_plot_mov_avg(Matrix dot_plot, int window)
{
	int				isize = dot_plot.columns(), jsize = dot_plot.rows();
	
	if ( window < 1 ) window = 1;
	if ( window > 1000 ) window = 1000;
	int				i, j, w, start, end;
	int				shortest = isize;
	if ( shortest > jsize ) shortest = jsize;
	int				low_half = window/2;
	int				up_half = window - low_half - 1;
	
	Matrix			mov_avg(isize, jsize);
	
	for ( j=0; j<jsize; j++ ) { 		// Diagonals in upper triangle
		for ( i=0; i<jsize-j && i<shortest; i++ ) {
			start = i - low_half;
			if ( start < 0 ) start = 0;
			end = i + up_half;
			if ( end >= shortest ) end = shortest - 1;
			if ( end >= jsize-j ) end = jsize - j - 1;
			for ( w=start; w<=end; w++ )
				mov_avg[i+j][i] += dot_plot[w+j][w];
		}
	}
	
	for ( i=1; i<isize; i++ ) { 		// Diagonals in lower triangle
		for ( j=0; j<isize-i && j<shortest; j++ ) {
			start = j - low_half;
			if ( start < 0 ) start = 0;
			end = j + up_half;
			if ( end >= shortest ) end = shortest - 1;
			if ( end >= isize-i ) end = isize - i - 1;
			for ( w=start; w<=end; w++ )
				mov_avg[j][i+j] += dot_plot[w][i+w];
		}
	}
	
	return mov_avg;
}

/**
@brief 	Finds the row and column maxima in a dot plot.
@param 	dot_plot	dot plot matrix.
@return int 			number of maxima in both rows and columns.
**/
int			seq_dot_plot_interpret(Matrix dot_plot)
{
	int			isize = dot_plot.columns(), jsize = dot_plot.rows();
	
	int 			i, j, k;
	double			max;
	int				nmask[4] = {0,0,0,0};
	vector<char>	mask(isize*jsize, 0);
	
	for ( j=0; j<jsize; j++ ) {
		max = 0;
		for ( i=0; i<isize; i++ )
			if ( max < dot_plot[j][i] )
				max = dot_plot[j][i];
		for ( i=0; i<isize; i++ )
			if ( max == dot_plot[j][i] ) mask[j*isize+i] += 1;
	}
	
	for ( i=0; i<isize; i++ ) {
		max = 0;
		for ( j=0; j<jsize; j++ )
			if ( max < dot_plot[j][i] )
				max = dot_plot[j][i];
		for ( j=0; j<jsize; j++ )
			if ( max == dot_plot[j][i] ) mask[j*isize+i] += 1;
	}
	
	for ( j=1; j<jsize-1; j++ ) {
		for ( i=1; i<isize-1; i++ ) {
			if ( mask[j*isize+i] > 1 ) {
				if ( mask[(j-1)*isize+i-1] ) mask[j*isize+i] = 3;
				if ( mask[(j+1)*isize+i+1] ) mask[j*isize+i] = 3;
			}
			k = mask[j*isize+i];
			nmask[k] += 1;
		}
	}
	
	cout << "Dot plot mask sums:" << endl << "n\tCount" << endl;
	for ( i=0; i<4; i++ )
		cout << i << tab << nmask[i] << endl;
	cout << endl;
	
	return nmask[3];
}

/**
@brief 	Finding the best segments in a dot plot.
@param 	dot_plot	dot plot matrix.
@return double 			threshold to identify unique segments.
**/
double		seq_dot_plot_best_segments(Matrix dot_plot)
{
	int				isize = dot_plot.columns(), jsize = dot_plot.rows();
	
	int 		i, j, n, nr=1, nc=1, nr1, nc1;
	double		threshold, tmin, tmax, max = 0;
	
	for ( i=0; i<isize; i++ )
		for ( j=0; j<jsize; j++ )
			if ( max < dot_plot[i][j] ) max = dot_plot[i][j];
	
	threshold = max/2;
	tmin = 0;
	tmax = max;
	
	cout << "Thresh\tUnique\tMulti" << endl;
	while ( nc+nr>0 || ( threshold - tmin > 1e-3 && tmax - threshold > 1e-3 ) ) {
		for ( j=nr=nr1=0; j<jsize; j++ ) {
			for ( i=n=0; i<isize; i++ )
				if ( dot_plot[j][i] > threshold ) n++;
			if ( n > 1 ) nr++;
			else if ( n == 1 ) nr1++;
		}
	
		for ( i=nc=nc1=0; i<isize; i++ ) {
			for ( j=n=0; j<jsize; j++ )
				if ( dot_plot[j][i] > threshold ) n++;
			if ( n > 1 ) nc++;
			else if ( n == 1 ) nc1++;
		}
		cout << threshold << tab << nc1+nr1 << tab << nc+nr << endl;
		if ( nc+nr > 0 ) {
			tmin = threshold;
			threshold = (threshold + tmax)/2;
		} else {
			tmax = threshold;
			threshold = (threshold + tmin)/2;
		}
	}
	
	cout << "Threshold:                      " << threshold << endl;
	cout << "Unique assignments:             " << nc1 << " " << nr1 << endl;
	cout << "Multiple assignments:           " << nc << " " << nr << endl << endl;
	
	return threshold;
}

/**
@brief 	Finds the best offset for aligning two sequences without gaps.
@param 	*mol1		first molecule structure.
@param 	*mol2		second molecule structure.
@param 	threshold		threshold to identify segments.
@param 	dot_plot	dot plot.
@return int 				the offset of sequence 2 with respect to sequence 1.

	The best offset is determined as the diagonal of a dot plot with
	the highest value.

**/
int			seq_show_segments(Bmolecule* mol1, Bmolecule* mol2, double threshold, Matrix dot_plot)
{
	int				isize = dot_plot.columns(), jsize = dot_plot.rows();
	
	int				i, j, notfound;
	
	for ( i=0; i<isize; i++ ) {
		cout << mol1->seq[i];
		notfound = 1;
		for ( j=0; j<jsize && notfound; j++ ) {
			if ( dot_plot[j][i] > threshold ) {
				cout << tab << mol2->seq[j] << tab << j;
				notfound = 0;
			}
		}
		cout << endl;
	}
	
	return 0;
}
