/**
@file	cluster.cpp
@brief	Clustering functions
@author Bernard Heymann 
@date	Created: 20070417
@date	Modified: 20210508
**/

#include "cluster.h"
#include "utilities.h"

//#ifndef DBL_MAX
//#define DBL_MAX 1.7976e308
//#endif

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generate clusters using a K-means algorithm.
@param 	n			number of data elements.
@param 	*data		floating point array.
@param 	k			number of classes.
@return	vector<long>	vector of cluster memberships.


**/
vector<long>	k_means(long n, float* data, long k)
{
	long			i, j, ik, jk, imax(100), nun[k], change(1);
	double			min, max, d, dmin, avg[k], sum[k];
	vector<long>	sel(n,0);
	
	min = max = data[0];
	for ( j=0; j<n; j++ ) {
		sel[j] = k;
		if ( min > data[j] ) min = data[j];
		if ( max < data[j] ) max = data[j];
	}

	for ( ik=0; ik<k; ik++ ) avg[ik] = (ik + 0.5)*(max - min)/k + min;

//	cout << avg[0] << tab << avg[1] << endl;
	
	for ( i=0; i<imax && change; i++ ) {
		change = 0;
		for ( ik=0; ik<k; ik++ ) {
			sum[ik] = 0;
			nun[ik] = 0;
		}
		for ( j=0; j<n; j++ ) {
			dmin = 1e37;
			for ( ik=jk=0; ik<k; ik++ ) {
				d = fabs(data[j] - avg[ik]);
				if ( dmin > d ) {
					dmin = d;
					jk = ik;
				}
			}
			sum[jk] += data[j];
			nun[jk]++;
			change += (sel[j] != jk);
			sel[j] = jk;
		}
		for ( ik=0; ik<k; ik++ ) {
			if ( nun[ik] ) avg[ik] = sum[ik]/nun[ik];
			else avg[ik] = (ik + 0.5)*(max - min)/k + min;
		}
//		cout << i << tab << avg[0] << tab << nun[0] << tab << avg[1] << tab << nun[1] << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "K-means clustering:" << endl;
		for ( ik=0; ik<k; ik++ )
			cout << "Class " << ik << ":          " << avg[ik] << " (" << nun[ik] << ")" << endl;
		cout << endl;
	}
	
	return sel;
}


/**
@brief 	Generate clusters from a similarity matrix using affinity propagation.  
@param 	s			n x n similarity matrix.
@param 	maxit		maximum iterations.
@param 	convit		convergence iterations.
@param 	lambda		damping factor.
@param 	&ncluster	number of clusters.
@return	vector<long>	vector of cluster memberships.

	Frey, B. J. and D. Dueck (2007). "Clustering by passing messages between 
	data points." Science 315(5814): 972-6.
	The clustering algorithm proceeds iteratively for a given maximum iterations.
	If the cluster solution does not change for a given number of iterations
	(convit), convergence is assumed and the function finishes.
	The damping factor should be ~0.5, or at least between 0.1 and 0.9.
	For larger damping factors, more convergence iterations are required.
	The diagonal of the matrix should be set to preferences for exemplars
	(typically the average or median of the similarities).

**/
/*long*		affin_prop_clustering(Matrix s, long maxit, long convit, double lambda, long& ncluster)
{
	if ( s.rows() < 2 ) {
		cerr << "Error: No input similarity matrix!" << endl;
		return NULL;
	}
	
	if ( lambda < 0.1 ) lambda = 0.1;
	if ( lambda > 0.9 ) lambda = 0.9;
	if ( convit < 2 ) convit = 2;
	if ( maxit < 2*convit ) maxit = 2*convit;
	
	int				n = s.rows();
	
	char			t;
	long			it, cit, i, j, k, nc;
	double			v;

	double*			a = new double[n*n];
	double*			r = new double[n*n];
	double*			mx1 = new double[n];
	double*			mx2 = new double[n];
	double*			srp = new double[n];
	char*			c = new char[n];
	long*			idx = new long[n];
	
//	for ( i=0; i<n*n; i++ ) s[i] = s[i] + (1e-16*s[i] + 1e-306)*(rand()/((double)RAND_MAX+1));

	for ( i=0; i<n*n; i++ ) a[i] = r[i] = 0;
	for ( j=0; j<n; j++ ) c[j] = 0;
	
	// The loop exits on two conditions:
	// 1. Exceeding the set maximum number of iterations
	// 2. No change in the cluster solution for the set number of convergence iterations
	//		(i.e., when cit >= convit)
	for ( it=cit=0; it<maxit && cit<convit; it++, cit++ ) {
		// Calculate responsibilities
		for ( i=j=0; j<n; j++ ) {
			mx1[j] = mx2[j] = -DBL_MAX;
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( mx1[j] < v ) {
					mx2[j] = mx1[j];
					mx1[j] = v;
				} else if ( mx2[j] < v ) {
					mx2[j] = v;
				}
			}
		}
		for ( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( v == mx1[j] ) v = mx2[j];
				else v = mx1[j];
				r[i] = lambda*r[i] + (1-lambda)*(s[j][k] - v);
			}
		}
		// Calculate availabilities
		for( k=0; k<n; k++ ) srp[k] = 0;
		for ( i=j=0; j<n; j++ )
			for ( k=0; k<n; k++, i++ )
				if( r[i] > 0 ) srp[k] += r[i];
		for( i=k=0; k<n; k++, i+=n+1 ) if ( r[i] < 0 ) srp[k] += r[i];
		for ( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) if ( j != k ) {
				v = srp[k];
				if( r[i] > 0 ) v -= r[i];
				if ( v < 0 ) a[i] = lambda*a[i]+(1-lambda)*v; 
				else a[i] = lambda*a[i];
			}
		}
		for( i=k=0; k<n; k++, i+=n+1 )
			a[i] = lambda*a[i] + (1-lambda)*(srp[k] - r[i]);
		// Update exemplars and test exit conditions
		for( i=j=0, nc=0; j<n; j++, i+=n+1 ) {
			if ( a[i] + r[i] > 0 ) t=1; else t=0;
			if ( c[j] != t ) cit = 0;	// Reset the convergence counter
			c[j] = t;
			nc += c[j];
		}
		if ( nc < 1 ) cit = 0;	// Reset the convergence counter
		if ( verbose & VERB_FULL )
			cout << it << tab << nc << endl;
	}

    double				netsim, dpsim(0), expref(0);
	for( j=0, nc=0; j<n; j++ ) nc += c[j];
	
	if ( nc > 0 ) {
		for( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				if ( c[k] ) a[i] = 0;
				else a[i] = -DBL_MAX;
			}
		}
		for ( i=j=0; j<n; j++ ) {
			mx1[j] = -DBL_MAX;
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( mx1[j] < v ) {
					mx1[j] = v;
					idx[j] = k;
				}
			}
			if( c[j] ) idx[j] = j;
		}
		dpsim = expref = 0;
		for( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				if ( idx[j] == k ) {
					if ( j == k ) expref += s[j][k];
					else dpsim += s[j][k];
				}
			}
		}
	}
	
    netsim = dpsim + expref;

	long*	m = new long[n];

	for ( j=0; j<n; j++ ) m[j] = 0;

	for ( j=0; j<n; j++ ) m[idx[j]]++;

	for ( j=ncluster=0; j<n; j++ ) if ( m[j] ) ncluster++;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Number of points:               " << n << endl;
		cout << "Number of iterations:           " << it << endl;
		cout << "Number of clusters:             " << nc << endl;
		cout << "Fitness (net similarity):       " << netsim << endl;
		cout << "  Similarities to exemplars:    " << dpsim << endl;
		cout << "  Preferences of exemplars:     " << expref << endl << endl;
		cout << "Cluster\tNumber" << endl;
		for ( j=0; j<n; j++ ) if ( j == idx[j] )
			cout << j+1 << tab << m[j] << endl;
		cout << endl;
		cout << "Point\tCluster\tSimilarity" << endl;
		for ( j=0; j<n; j++ )
			cout << j+1 << tab << idx[j]+1 << tab << s[idx[j]][j] << endl;
		cout << endl;
	}
	
	delete[] a;
	delete[] r;
	delete[] mx1;
	delete[] mx2;
	delete[] srp;
	delete[] c;
	delete[] m;

	return idx;
}
*/
vector<long>	affin_prop_clustering(Matrix s, long maxit, long convit, double lambda, long& ncluster)
{
	if ( s.rows() < 2 ) {
		cerr << "Error: No input similarity matrix!" << endl;
		return vector<long>(0);
	}
	
	if ( lambda < 0.1 ) lambda = 0.1;
	if ( lambda > 0.9 ) lambda = 0.9;
	if ( convit < 2 ) convit = 2;
	if ( maxit < 2*convit ) maxit = 2*convit;
	
	int				n = s.rows();
	
	char			t;
	long			it, cit, i, j, k, nc;
	double			v;

	vector<double>	a(n*n,0);
	vector<double>	r(n*n,0);
	vector<double>	mx1(n,0);
	vector<double>	mx2(n,0);
	vector<double>	srp(n,0);
	vector<char>	c(n,0);
	vector<long>	idx(n,0);
	
	// The loop exits on two conditions:
	// 1. Exceeding the set maximum number of iterations
	// 2. No change in the cluster solution for the set number of convergence iterations
	//		(i.e., when cit >= convit)
	for ( it=cit=0; it<maxit && cit<convit; it++, cit++ ) {
		// Calculate responsibilities
		for ( i=j=0; j<n; j++ ) {
			mx1[j] = mx2[j] = -DBL_MAX;
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( mx1[j] < v ) {
					mx2[j] = mx1[j];
					mx1[j] = v;
				} else if ( mx2[j] < v ) {
					mx2[j] = v;
				}
			}
		}
		for ( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( v == mx1[j] ) v = mx2[j];
				else v = mx1[j];
				r[i] = lambda*r[i] + (1-lambda)*(s[j][k] - v);
			}
		}
		// Calculate availabilities
		for( k=0; k<n; k++ ) srp[k] = 0;
		for ( i=j=0; j<n; j++ )
			for ( k=0; k<n; k++, i++ )
				if( r[i] > 0 ) srp[k] += r[i];
		for( i=k=0; k<n; k++, i+=n+1 ) if ( r[i] < 0 ) srp[k] += r[i];
		for ( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) if ( j != k ) {
				v = srp[k];
				if( r[i] > 0 ) v -= r[i];
				if ( v < 0 ) a[i] = lambda*a[i]+(1-lambda)*v;
				else a[i] = lambda*a[i];
			}
		}
		for( i=k=0; k<n; k++, i+=n+1 )
			a[i] = lambda*a[i] + (1-lambda)*(srp[k] - r[i]);
		// Update exemplars and test exit conditions
		for( i=j=0, nc=0; j<n; j++, i+=n+1 ) {
			if ( a[i] + r[i] > 0 ) t=1; else t=0;
			if ( c[j] != t ) cit = 0;	// Reset the convergence counter
			c[j] = t;
			nc += c[j];
		}
		if ( nc < 1 ) cit = 0;	// Reset the convergence counter
		if ( verbose & VERB_FULL )
			cout << it << tab << nc << endl;
	}

    double				netsim, dpsim(0), expref(0);
	for( j=0, nc=0; j<n; j++ ) nc += c[j];
	
	if ( nc > 0 ) {
		for( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				if ( c[k] ) a[i] = 0;
				else a[i] = -DBL_MAX;
			}
		}
		for ( i=j=0; j<n; j++ ) {
			mx1[j] = -DBL_MAX;
			for ( k=0; k<n; k++, i++ ) {
				v = a[i] + s[j][k];
				if ( mx1[j] < v ) {
					mx1[j] = v;
					idx[j] = k;
				}
			}
			if( c[j] ) idx[j] = j;
		}
		dpsim = expref = 0;
		for( i=j=0; j<n; j++ ) {
			for ( k=0; k<n; k++, i++ ) {
				if ( idx[j] == k ) {
					if ( j == k ) expref += s[j][k];
					else dpsim += s[j][k];
				}
			}
		}
	}
	
    netsim = dpsim + expref;

	vector<long>	m(n,0);

	for ( j=0; j<n; j++ ) m[idx[j]]++;

	for ( j=ncluster=0; j<n; j++ ) if ( m[j] ) ncluster++;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Number of points:               " << n << endl;
		cout << "Number of iterations:           " << it << endl;
		cout << "Number of clusters:             " << nc << endl;
		cout << "Fitness (net similarity):       " << netsim << endl;
		cout << "  Similarities to exemplars:    " << dpsim << endl;
		cout << "  Preferences of exemplars:     " << expref << endl << endl;
		cout << "Cluster\tNumber" << endl;
		for ( j=0; j<n; j++ ) if ( j == idx[j] )
			cout << j+1 << tab << m[j] << endl;
		cout << endl;
		cout << "Point\tCluster\tSimilarity" << endl;
		for ( j=0; j<n; j++ )
			cout << j+1 << tab << idx[j]+1 << tab << s[idx[j]][j] << endl;
		cout << endl;
	}
	
	return idx;
}
