/**
@file	Matrix.h 
@brief	Generalized matrix class
@author Bernard Heymann
@date	Created: 20000501
@date	Modified: 20180723
**/

#include "random_numbers.h"
#include <iostream>
#include <fstream>

using namespace std;

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#ifndef _Matrix_
#define _Matrix_

class Matrix
{
private:
	long		m, n, len;
	double*		d;
	double**	p;
	void		init() {
		len = m*n;
		d = new double[len];
		p = new double*[m];
		for ( long i=0; i<m; i++ ) p[i] = d + i*n;
		for ( long i=0; i<len; i++ ) d[i] = 0;
//		cout << d << "\t" << &d[8] << "\t" << p[2] << endl;
	}
	void		clear() {
		delete[] d; delete[] p;
		d = NULL; p = NULL;
	}
public:
	Matrix() : m(0), n(0) { d = NULL; p = NULL; }
	Matrix(const Matrix& mat) : m(mat.m), n(mat.n) {
		init();
		for ( long i=0; i<len; i++ ) d[i] = mat.d[i];
	}
	Matrix(long rows, long cols) : m(rows), n(cols) {
		init();
	}
	Matrix(Bstring& filename) {
		m = n = 0;
		ifstream		fmat(filename.c_str());
		if ( fmat.fail() ) return;
		char			buf[10000];
		long		 	i(0), j;
		Bstring			sline, *tokens, *token;
		while ( fmat.getline(buf, 1024) && strncmp(buf, "matrix_", 7) ) ;
		if ( verbose & VERB_DEBUG )
			cout << buf << endl;
		while ( fmat.getline(buf, 10000) ) {
			sline = buf;
			tokens = sline.split();
					for ( j=0, token = tokens; token; token = token->next ) ++j;
//			cout << j << tab << sline << endl;
			if ( tokens ) {
				if ( *tokens == "_rows" ) {
					m = tokens->next->integer();
					if ( m && n ) { init(); }
//					cout << "Rows: " << m << endl;
				} else if ( *tokens == "_columns" ) {
					n = tokens->next->integer();
					if ( m && n ) { init(); }
//					cout << "Columns: " << n << endl;
				} else {
					for ( j=0, token = tokens; token; token = token->next ) ++j;
					if ( verbose & VERB_DEBUG )
						cout << "Row " << i+1 << ": " << j << endl;
					if ( i < m && j >= n ) {
						for ( j=0, token = tokens; j<n; j++, token = token->next )
							(*this)[i][j] = token->real();
						i++;
					}
				}
				string_kill(tokens);
			}
		}
		fmat.close();
		if ( m < 1 || n < 1 ) {
			cerr << "Error: The matrix size is zero!" << endl;
			bexit(-1);
		}
		if ( verbose & VERB_DEBUG )
			cout << filename << " read" << endl;
	}
	~Matrix() { clear(); }
	
	void	write(Bstring& filename) {
		ofstream        fmat(filename.c_str());
		if ( fmat.fail() ) return;
		long	 		i, j;
		fmat << "# Matrix written by Bsoft\n\ndata_\n\nmatrix_" << endl;
		fmat << "_rows" << tab << m << endl;
		fmat << "_columns" << tab << n << endl;
		for ( i=0; i<m; i++ ) {
			fmat << (*this)[i][0];
			for ( j=1; j<n; j++ ) fmat << " " << (*this)[i][j];
			fmat << endl;
		}
		fmat << endl;
		fmat.close();
	}
	
	double*	data() { return d; }
//	double*	row(long i) { while ( i>=m ) i -= m; return p; }
	
	class Row {
		private:
			long	n;
			double*	p;
		public:
			Row(long len, double* row) : n(len), p(row) { }
			double&	operator[](long i) {
				while ( i>=n ) i -= n;
				return p[i];
			}
	};
	
	Row		operator[](long i) {
		while ( i>=m ) i -= m;
		return Row(n, p[i]);
	}

	Matrix	operator=(const Matrix mat) {
		clear();
		m = mat.m; n = mat.n;
		init();
		for ( long i=0; i<len; i++ ) d[i] = mat.d[i];
		return *this;
	}
	Matrix	operator-() {
		Matrix		mat(*this);
		for ( long i=0; i<len; i++ ) mat.d[i] = -mat.d[i];
		return mat;
	}
	Matrix	operator+=(Matrix mat) {
		return *this + mat;
	}
	Matrix	operator+(Matrix mat) {
		long		i, j;
		Matrix		numat(m, n);
		if ( m != mat.m || n != mat.n ) {
			cerr << "Matrices not the same size!" << endl;
			return numat;
		}
		for ( i=0; i<n; i++ )
			for ( j=0; j<m; j++ )
				numat[i][j] = (*this)[i][j] + mat[i][j];
		return numat;
	}
	Matrix	operator-=(Matrix mat) {
		return *this - mat;
	}
	Matrix	operator-(Matrix mat) {
		long		i, j;
		Matrix		numat(m, n);
		if ( m != mat.m || n != mat.n ) {
			cerr << "Matrices not the same size!" << endl;
			return numat;
		}
		for ( i=0; i<n; i++ )
			for ( j=0; j<m; j++ )
				numat[i][j] = (*this)[i][j] - mat[i][j];
		return numat;
	}
	Matrix	operator*=(Matrix mat) {
		return *this * mat;
	}
	Matrix	operator*(Matrix& mat) {
		long		i, j, k;
		Matrix		numat(mat.m, n);
		if ( m != mat.n ) {
			cerr << "Matrix rows not equal to second vector columns!" << endl;
			return numat;
		}
		for ( i=0; i<n; i++ )
			for ( j=0; j<mat.m; j++ )
				for ( k=0; k<m; k++ ) numat[i][j] += (*this)[i][k]*mat[k][j];
		return numat;
	}
	vector<double>	operator*(vector<double>& vec) {
		long			i, j;
		vector<double>	nuvec(m,0);
		if ( n != vec.size() ) {
			cerr << "Matrix columns not equal to vector size!" << endl;
			return nuvec;
		}
		for ( i=0; i<m; i++ )
			for ( j=0; j<n; j++ )
				nuvec[i] += (*this)[i][j]*vec[j];
		return nuvec;
	}

	long	rows() { return m; }
	long	columns() { return n; }
	long	size() { return len; }

	void	show_below_cutoff(double d) {
		for ( long i=0; i<m; ++i ) {
			cout << i;
			for ( long j=0; j<n; ++j )
				if ( (*this)[i][j] <= d )
					cout << tab << j << "," << (*this)[i][j];
			cout << endl;
		}
	}

	void	swap_rows_columns(long rc1, long rc2) {
		long			i;
		for ( i=0; i<m; i++ ) swap((*this)[i][rc1], (*this)[i][rc2]);
		for ( i=0; i<n; i++ ) swap((*this)[rc1][i], (*this)[rc2][i]);
	}
	
	Matrix	delete_row_column(long rc) {
		long			i, j, k, l;
		Matrix			mat(m-1,n-1);
		for ( i=k=0; i<m; i++ ) {
			for ( j=l=0; j<n; j++ ) {
				if ( i != rc && j != rc ) mat[k][l] = p[i][j];
				if ( j != rc ) l++;
			}
			if ( i != rc ) k++;
		}
		return mat;
	}

	Matrix	transpose() const {
		Matrix			mat(n,m);
		for ( long i=0; i<n; i++ )
			for ( long j=0; j<m; j++ )
				mat[i][j] = p[j][i];
		return mat;
	}

	void	fill(double v) { for ( long i=0; i<len; i++ ) d[i] = v; }

	int		check_for_singularity() {
		int				c(0);
		long			i, j, k;
		double			max;
		for ( i=k=0; i<m; i++ ) {
			for ( j=0, max=0; j<n; j++, k++ ) if ( max < fabs(d[k]) ) max = fabs(d[k]);
			if ( max < 1e-37 ) c |= 1 << i;
		}
		return c;
	}
	/**
		The rows and columns are alternatively iteratively normalized until
		the error is small enough.
	**/
	void			normalize()
	{
		long			i, j, r, c;
		double			err(1);
		vector<double>	rw(m);
		vector<double>	cw(n);
	
		if ( verbose & VERB_FULL )
			cout << "Cycle\tError" << endl;
		for ( i=0; i<100 && err > 1e-20; i++ ) {
			for ( r=0; r<m; r++ ) rw[r] = 0;	// Row scaling
			for ( j=r=0; r<m; r++ )
				for ( c=0; c<n; c++, j++ )
					rw[r] += d[j];
			for ( r=0; r<m; r++ ) rw[r] = 1/rw[r];
			for ( j=r=0; r<m; r++ )
				for ( c=0; c<n; c++, j++ )
					d[j] *= rw[r];

			for ( c=0; c<n; c++ ) cw[c] = 0;	// Column scaling
			for ( j=r=0; r<m; r++ )
				for ( c=0; c<n; c++, j++ )
					cw[c] += d[j];
			for ( c=0; c<n; c++ ) cw[c] = 1/cw[c];
			for ( j=r=0; r<m; r++ )
				for ( c=0; c<n; c++, j++ )
					d[j] *= cw[c];

			for ( err=0, r=0; r<m; r++ ) err += rw[r];
			for ( c=0; c<n; c++ ) err += cw[c];
			err = fabs(err/(m + n) - 1);
			if ( verbose & VERB_FULL )
				cout << i+1 << "\t" << err << endl;
		}
	}
	void		randomize() {
		random_seed();
		long		i, j;
		double		rm = INT_MAX/4, irm = 10.0L/INT_MAX;
		for ( i=0; i<m; i++ )
			for ( j=0; j<n; j++ )
				p[i][j] = irm*(random() - rm);
	}
//	int			multiply_in_place(double* vec);
	int			multiply_in_place(vector<double>& vec);
	double		determinant() { return LU_decomposition(); }
	double 		LU_decomposition();
/*	double 		LU_decomposition(double* b) {
		double	det = LU_decomposition();
		if ( b ) multiply_in_place(b);
		return det;
	}*/
	double 		LU_decomposition(vector<double>& b) {
		double	det = LU_decomposition();
		multiply_in_place(b);
		return det;
	}
	double 		singular_value_decomposition();
/*	double 		singular_value_decomposition(double* b) {
		singular_value_decomposition();
		if ( b ) multiply_in_place(b);
		return 0;
	}*/
	double 		singular_value_decomposition(vector<double>& b) {
		singular_value_decomposition();
		multiply_in_place(b);
		return 0;
	}
	int			jrotate(double s, double tau, long i, long j, long k, long l);
//	double*		jacobi_rotation_old();
	vector<double>	jacobi_rotation();
//	long		jacobi_rotation_row(long ip, long n, long i, double* val, double* z, Matrix* vec);
//	double*		jacobi_rotation_parallel();
//	void 		eigen_sort(double* val);
	void 		eigen_sort(vector<double>& val);
};

ostream& operator<<(ostream& output, Matrix& mat);

#endif

// Function prototypes 
Vector3<double> principal_axes(Vector3<double> avg, Vector3<double> avg2, Vector3<double> avgx, Vector3<double>* eigenvec);
Vector3<double> principal_axes(vector<Vector3<double>>& coor, Matrix& eigenvec);
vector<double>	dsyevq3(Matrix& A);
void 			dsytrd3(Matrix& A, double d[3], double e[2]);
