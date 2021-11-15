/**
@file	binomial.h
@brief	Functions to generate and fit binomials 
@author Daniel Nemecek and Bernard Heymann
@date	Created: 20091202
@date	Modified: 20190201
**/


// Function prototypes 
double		Bnpk(int n, double p, int k, double w);
double		binomial_fit(vector<double>& distrib, int n, int nfit, vector<double>& coeff);

