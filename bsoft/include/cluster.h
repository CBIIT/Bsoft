/**
@file	cluster.h 
@brief	Header file for clustering functions 
@author Bernard Heymann 
@date	Created: 20070417
@date	Modified: 20210508
**/

#include "Matrix.h"

// Function prototypes 
vector<long>	k_means(long n, float* data, long k);
vector<long>	affin_prop_clustering(Matrix s, long maxit, long convit, double lambda, long& ncluster);

