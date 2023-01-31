/**
@file	histogram.cpp
@brief	Library routines to calculate histograms for images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20221004
**/

#include "histogram.h"
#include "ps_plot.h"
#include "cluster.h"
#include "moving_average.h"
#include "matrix_linear.h"
#include "qsort_functions.h"
#include "string_util.h"
#include "utilities.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the histogram froma data array.
@param	data			vector of data.
@param 	bins			number of bins in the histogram.
@param 	&scale			scale.
@param 	&offset			offset.
@return vector<long>		histogram.

	A histogram of an image is calculated with a given number of bins.
	Multiple channels are output as successive one-dimensional arrays.
	The image data is not affected.
	The statistics for the input image must be correctly calculated.

**/
vector<long>	histogram(vector<double> data, long bins, double& scale, double& offset)
{
	scale = 1;
	offset = 0;

	vector<long>	histo(bins, 0);
	
	if ( data.size() < 2 ) return histo;
	
    long			i;
    double			max(data[0]), min(data[0]);
    
	for ( auto v: data ) {
//		cout << v << endl;
		if ( min > v ) min = v;
		if ( max < v ) max = v;
	}
	
	scale = (bins - 1)*1.0/(max - min);
	offset = 0.5 - scale*min;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG histogram: range=" << min << "-"
			<< max << " bins=" << bins << " scale=" << scale
			<< " offset=" << offset << endl;

	for ( auto v: data ) {
		i = (long) (scale*v + offset);
		if ( i >= 0 && i < bins ) histo[i]++;
	}
	
//	for ( i=0; i<bins*c; i++ ) cout << i << tab << histo[i] << endl;
	
	return histo;
}

/*
vector<double>	read_text_column(Bstring filename, long column)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_text_column: file=" << filename << " col=" << column << endl;
		
	vector<double>		data;

	ifstream			ftxt(filename.c_str());
	if ( ftxt.fail() ) return data;
	
	long				i;
	string				s;
	
	while ( !ftxt.eof() ) {
		getline(ftxt, s);
		vector<string>	vs = split(s);
		for ( i=0; i<vs.size() && i<column; ++i ) ;
		if ( i<vs.size() && i == column ) data.push_back(to_real(vs[i]));
	}
	
	ftxt.close();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_text_column: rows read=" << data.size() << endl;

	return data;
}
*/

Bplot*		plot_convert_to_histogram(Bplot* plot, long bins, long hiscol)
{
	long			i, j, ncol(2);
	double			scale(1), offset(0);
	vector<double>	data(plot->rows());
	
	for ( i=0, j=bins*hiscol; i<bins; ++i, ++j )
		data[i] = (*plot)[j];
	
	vector<long>	histo = histogram(data, bins, scale, offset);

	double			min(-offset/scale);
	double			max((bins-1)/scale + min);
	double			bin_half = (max - min)/(4*bins);

	Bstring			title("Histogram"), txt;
	Bplot*			hisplot = new Bplot(1, bins, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Bin");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).element_size(0.5);
	plot->page(0).column(1).label("Count");
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).color(0.5,0.5,0.5);
	plot->page(0).axis(1).min(min - bin_half);
	plot->page(0).axis(1).max(max + bin_half);
	
	for ( i=0, j=bins; i<bins; ++i, ++j ) {
		(*plot)[i] = min + i/scale;
		(*plot)[j] = histo[i];
	}
	
	return hisplot;
}


double		find_maximum_sigma(Matrix& H, long i, long level, long depth, double sum, double& maxsum, vector<double>& t)
{
	if ( level >= depth ) {
		sum += H[i+1][H.rows()-1];
//		cout << tab << i << tab << sum << endl;
		return sum;
	}
	
	long		j, k(0);
	double		thissum(0), thismax(0);

	level++;

	for ( j=i+1; j<H.rows(); ++j ) {
//		cout << tab << j << tab << H[i+1][j];
		thissum = find_maximum_sigma(H, j, level, depth, sum + H[i+1][j], maxsum, t);
		if ( thismax <= thissum ) {
			thismax = thissum;
			k = j;
//			cout << tab << level << tab << j << tab << thismax << endl;
		}
	}
	
	if ( maxsum <= thismax ) {
		maxsum = thismax;
		t[level-1] = k;
	}

	return thismax;
}

/**
@brief 	Calculates multiple thresholds from a histogram.
@param 	h				histogram.
@param 	number			number of clusters (one more than thresholds).
@return vector<double> 		thresholds.

	Reference: PS.Liao, TS.Chen, and PC. Chung,
           Journal of Information Science and Engineering, vol 17, 713-727 (2001)
**/
vector<double>	histogram_thresholds(vector<long> h, long number)
{
	long			i, j, bins(h.size());

	// Set up matrices
	Matrix			P(bins,bins);
	Matrix			S(bins,bins);
	Matrix			H(bins,bins);
	
	// Diagonals
	for ( i=0; i<bins; ++i ) {
		P[i][i] = h[i];
		S[i][i] = h[i] * i;
	}
	
	// Second rows
	for ( i=1, j=2; i<bins-1; ++i, ++j ) {
		P[1][j] = P[1][i] + h[j];
		S[1][j] = S[1][i] + j*h[j];
	}
	
	// Propagate sums
	for ( i=2; i<bins; ++i ) {
		for ( j=i+1; j<bins; ++j ) {
			P[i][j] = P[1][j] - P[1][i-1];
			S[i][j] = S[1][j] - S[1][i-1];
		}
	}
	
	// H
	for ( i=1; i<bins; ++i )
		for ( j=i+1; j<bins; ++j )
			if ( P[i][j] ) H[i][j] = S[i][j]*S[i][j]/P[i][j];
	
	long				depth(number-1);
	double				maxsum(0);
	vector<double>		t(depth,0);
	
	find_maximum_sigma(H, 0, 0, depth, 0, maxsum, t);
	
	if ( verbose )
		cout << "Best FOM:                        " << maxsum << endl << endl;
	
	return t;
}


double		histogram_gaussian_R(Bsimplex& simp)
{
	long			i, j, ng(simp.constant(0));
	double			v, df, R(0);
	vector<double>	amp(ng), pos(ng), invsig(ng);
	for ( i=j=0; i<ng; ++i ) {
		amp[i] = simp.parameter(j++);
		pos[i] = simp.parameter(j++);
		invsig[i] = 1/simp.parameter(j++);
	}
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		for ( j=0, df=0; j<ng; j++ ) {
			v = (x[i] - pos[j])*invsig[j];
			df += amp[j]*exp(-0.5*v*v);
		}
		df -= f[i];
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}

/**
@brief 	Fits a gaussian function to a histogram of an image.
@param 	plot			plot with histogram.
@param 	ngauss			number of gaussians.
@return vector<double> 	array with 3 values for each gaussian, final value R.

	The return has (in order) the amplitude, offset,
	and sigma value for each gaussian.
	The last value is the RMSD, scaled as:
		RMSD*bins/sum

**/
vector<double>	plot_histogram_fit_gaussian(Bplot* plot, long ngauss)
{
	long			i, j, k, l;
	long			bins(plot->rows()), hsum(0);
	double			hmax(0);
	vector<double>	mh(bins);
	vector<double>	xx(bins), fx(bins);
	for ( i=0, j=bins; i<bins; ++i, ++j ) {
		xx[i] = (*plot)[i];
		fx[i] = (*plot)[j];
		hsum += fx[i];
		if ( hmax < fx[i] ) hmax = fx[i];
		mh[i] = fx[i];
	}

	// Attempting to find peaks based on the gradient change
	long			w(bins/20);
	if ( w < 5 ) w = 5;
	mh = moving_gradient(mh, w);
	
	vector<double>	gmax(ngauss);
	double			ex;
	for ( i=j=k=l=0, ex=0; i<bins; i++ ) {
		if ( ex < mh[i] ) {
			ex = mh[i];
			j = i;
		} else if ( mh[i] < 0 && mh[j] > 10 ) {
			gmax[k++] = i;
			ex = 0;
			j = i;
		}
	}
	
	// 3 parameters: amplitude, position, sigma
	Bsimplex	simp(1, 3*ngauss, 1, bins, xx, fx);
	simp.constant(0, ngauss);

	for ( i=j=0; j<ngauss; j++, i++ ) {
		simp.parameter(i, hmax);		// Amplitude
		simp.limits(i, hsum/bins, 5*hmax);
//		cout << tab << simp.parameter(i);
		i++;
		simp.parameter(i, gmax[j]);		// Position
		if ( ngauss == 1 ) {
			simp.limits(i, 0, bins-1);
		} else if ( j == 0 ) {
			simp.limits(i, 0, (gmax[0]+gmax[1])/2);
		} else if ( j < ngauss-1 ) {
			simp.limits(i, (gmax[j-1]+gmax[j])/2, (gmax[j]+gmax[j+1])/2);
		} else {
			simp.limits(i, (gmax[0]+gmax[1])/2, bins-1);
		}

		i++;
		simp.parameter(i, bins*0.1/ngauss);	// Sigma
		simp.limits(i, 3, bins*0.5);
//		cout << tab << simp.parameter(i) << endl;
	}

	double			R = simp.run(5000, 1e-6, histogram_gaussian_R);
	R = bins*R/hsum;
	
	vector<double>	gauss;
	for ( j=0; j<3*ngauss; j++ ) gauss.push_back(simp.parameter(j));
	gauss.push_back(R);

	if ( verbose & VERB_PROCESS ) {
		cout << "Gaussian fit over " << bins << " bins:" << endl;
		cout << "Amp\tOffset\tSigma\tSIP" << endl;
		for ( j=0; j<ngauss; j++ )
			cout << gauss[3*j] << tab << gauss[3*j+1] << tab << gauss[3*j+2] << tab << bins*gauss[3*j]/hsum << endl;
		cout << "Sampling-independent RMSD:      " << R << endl << endl;
	}

	return gauss;
}

/**
@brief 	Fits a gaussian function to a histogram and plots it.
@param 	plot			histogram plot.
@param 	ngauss			number of gaussians.
@return int 				0.

	One or more gaussians are fit to a histogram plot.
	The gaussian curves are added to columns.

**/
int 		plot_histogram_fit(Bplot* plot, long ngauss)
{
	long			bins(plot->rows());
	long			i, j, k, m, ncol(ngauss+plot->columns());
	double			v;
	Bstring			txt;

	vector<double>	gauss = plot_histogram_fit_gaussian(plot, ngauss);
	
	plot->add_columns(ngauss);
	plot->page(0).add_columns(ngauss);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Bin");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(1).label("Count");
	for ( i=2; i<ncol; i++ )
		plot->page(0).column(i).label("Gaussian");
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).color(0.5,0.5,0.5);
	if ( ngauss > 0 ) plot->page(0).column(2).color(1,0,0);
	if ( ngauss > 1 ) plot->page(0).column(3).color(0,1,0);
	if ( ngauss > 2 ) plot->page(0).column(4).color(0,0,1);
	
	if ( ngauss ) {
		for ( j=0; j<ngauss; j++ ){
			txt = Bstring(j+1, "Gaussian %d:") + Bstring(gauss[3*j], " a=%g") +
				Bstring(gauss[3*j+1], " u=%g") + Bstring(gauss[3*j+2], " s=%g");
			plot->page(0).add_text(txt);
		}
		txt = Bstring(gauss[3*ngauss], "R: %g");
		plot->page(0).add_text(txt);
	}
	
	for ( i=0, j=bins; i<bins; i++, j++ ) {
		for ( k=0, m=i+2*bins; k<3*ngauss; k+=3, m+=bins ) {
			v = ((*plot)[i] - gauss[k+1])/gauss[k+2];
			(*plot)[m] = gauss[k]*exp(-0.5*v*v);
		}
	}

	return 0;
}
