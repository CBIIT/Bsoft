/**
@file	histogram.h
@brief	Header file for calculating and fitting histograms
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20221004
**/

#include "simplex.h"
#include "Bplot.h"

// Function prototypes 
vector<long>	histogram(vector<double> data, long bins, double& scale, double& offset);
Bplot*		plot_convert_to_histogram(Bplot* plot, long bins, long hiscol);
vector<double>	histogram_thresholds(vector<long> h, long number);
double		histogram_gaussian_R(Bsimplex& simp);
//vector<double>	histogram_gauss_fit(vector<long> histo, double scale, double offset, long ngauss);
//Bplot* 		histogram_gauss_plot(vector<long> histo, double scale, double offset, long ngauss);
int 		plot_histogram_fit(Bplot* plot, long ngauss);
