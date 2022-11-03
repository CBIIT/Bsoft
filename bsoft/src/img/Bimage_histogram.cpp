/**
@file	Bimage_histogram.cpp
@brief	Library routines to calculate histograms for images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210304
**/

#include "Bimage.h"
#include "histogram.h"
#include "cluster.h"
#include "moving_average.h"
#include "matrix_linear.h"
#include "qsort_functions.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the histogram of an image.
@param 	bins			number of bins in the histogram.
@param 	&scale			scale.
@param 	&offset			offset.
@return vector<long>	histogram.

	A histogram of an image is calculated with a given number of bins.
	Multiple channels are output as successive one-dimensional arrays.
	The image data is not affected.
	The statistics for the input image must be correctly calculated.

**/
vector<long>	Bimage::histogram(long bins, double& scale, double& offset)
{
	scale = 1;
	offset = 0;

	vector<long>	histo(bins*c, 0);
	
	if ( !data_pointer() ) return histo;
	
    long			i, j, cc;
//	for ( i=0; i<bins*c; i++ ) histo[i] = 0;
	
	if ( datatype < Float ) {
		scale = 1/(ceil((max - min + 1)/bins));
		offset = -scale*min;
	} else {
		scale = (bins - 1)*1.0/(max - min);
		offset = 0.5 - scale*min;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::histogram: range=" << min << "-"
			<< max << " bins=" << bins << " scale=" << scale
			<< " offset=" << offset << endl;

	for ( i=0; i<n*image_size();  ) {
		for ( cc=0; cc<c; cc++, i++ ) {
			j = (long) (scale*(*this)[i] + offset);
			if ( j >= 0 && j < bins ) histo[cc*bins+j]++;
		}
	}
	
//	for ( i=0; i<bins*c; i++ ) cout << i << tab << histo[i] << endl;
	
	return histo;
}

/**
@brief 	Calculates the histogram of an image.
@param 	bins		number of bins in the histogram.
@return Bplot*		0.

	A histogram of an image is calculated with a given number of bins.
	Multiple channels are output as successive one-dimensional arrays.
	The image data is not affected.
	The statistics for the input image must be correctly calculated.
	If the postscript file name is given, a postscript plot is produced.

**/
Bplot* 		Bimage::histogram(long bins)
{
    long		  	i, j, k, cc, npage, ncol;
	double			Hi, H;
	double			bin_half = (max - min)/(4*bins);
	double			unit_mass = RHO * sampling(0).volume();
	vector<long>	hsum(bins*c, 0);
    
	double			scale, offset;
	vector<long>	histo = histogram(bins, scale, offset);
	
	long			hmax(0);
	for ( i=0; i<bins; i++ ) if ( hmax < histo[i] ) hmax = histo[i];

	for ( i=0; i<bins*c; i++ ) hsum[i] = histo[i];
    for ( i=bins-1; i>0; i-- )
		for ( cc=0; cc<c; cc++ )
			hsum[cc*bins+i-1] += hsum[cc*bins+i];
	
	if ( c > 1 ) {
		npage = 1;
		ncol = c + 1;
	} else {
		npage = 2;
		ncol = 3;
	}
	
	if ( verbose & VERB_LABEL ) {
    	cout << "Histogram" << endl;
    	if ( c == 1 )
			cout << "Bin\tLevel\tNumber\tSum\tVolume\tInformation" << endl;
    	else {
			cout << "Bin\tLevel";
			for ( cc=0; cc<c; cc++ ) cout << "\tChannel" << cc+1;
			cout << endl;
		}
	}

	Bstring			title("Histogram");
	Bplot*			plot = new Bplot(npage, bins, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(c+1);
	for ( i=0; i<=c; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Bin");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<=c; i++ ) {
		plot->page(0).column(i).type(1);
		plot->page(0).column(i).label("Count");
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	if ( c > 1 ) {
		plot->page(0).column(1).color(1,0,0);
		plot->page(0).column(2).color(0,1,0);
		plot->page(0).column(3).color(0,0,1);
	}
	plot->page(0).axis(1).min(min - bin_half);
	plot->page(0).axis(1).max(max + bin_half);
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(hmax);
	
	if ( c < 2 ) {
		title = "Cumulative mass";
		plot->page(1).title(title);
		plot->page(1).columns(2);
		plot->page(1).column(0).number(0);
		plot->page(1).column(1).number(2);
		plot->page(1).column(0).label("Bin");
		plot->page(1).column(0).axis(1);
		plot->page(1).column(1).type(1);
		plot->page(1).column(1).label("Mass(MDa)");
		plot->page(1).column(1).axis(3);
		plot->page(1).column(1).element_size(0.5);
		plot->page(1).axis(1).min(min - bin_half);
		plot->page(1).axis(1).max(max + bin_half);
		plot->page(1).axis(3).min(0);
		plot->page(1).axis(3).max(hsum[0] * unit_mass);
	}
	
    for ( i=0; i<bins; i++ ) {
		(*plot)[i] = min + i/scale;
		for ( cc=0, j=i, k=bins+i; cc<c; cc++, j+=bins, k+=bins ) (*plot)[k] = histo[j];
	}

	if ( c == 1 ) {
		H = 0;
    	for ( i=0, k=2*bins; i<bins; i++, k++ ) {
			(*plot)[k] = hsum[i] * unit_mass;
			Hi = 0;
			if ( histo[i] )
				Hi = (-histo[i]*1.0/hsum[0]) * log(histo[i]*1.0/hsum[0]) / log(2.0);
			H += Hi;
	    	cout << i+1 << tab << setw(12) << (*plot)[i] << tab << histo[i] << tab << hsum[i] << tab <<
				sampling(0).volume()*hsum[i] << tab << Hi << endl;
		}
	
		if ( verbose & VERB_LABEL ) {
	    	cout << "Voxels & number in histogram:   " << datasize << " " 
				<< hsum[0] << " (" << hsum[0]*100.0/datasize << " %)" << endl;
	    	cout << "Voxel volume:                   " << sampling(0).volume() << " A/voxel" << endl;
    		cout << "Informational entropy and rate: " << H << " " << log(1.0*bins)/log(2.0) - H << endl << endl;
		}
	} else {
    	for ( i=0; i<bins; i++ ) {
	    	cout << i+1 << tab << (*plot)[i];
	    	for ( cc=0, j=0; cc<c; cc++, j+=bins ) cout << tab << histo[j];
	    	cout << endl;
		}
	}
	
    return plot;
}

/**
@brief 	Finds the peaks in a quantized image from the histogram.
@param	flags		1=plot, 2=convert image.
@return Bplot*		plot of fits.

	A histogram of an image with quantized data is calculated.
	The peaks are determined and the image converted.

**/
Bplot*		Bimage::histogram_counts(int flags)
{
	long			i, j, k, bins(1000), imax(0);
 	double			phi, scale, offset, hmax(0);
 	Complex<double>	cv;
	
	if ( datatype < Float ) bins = max - min + 1;
	
	vector<long>	h = histogram(bins, scale, offset);
	
	if ( verbose )
		cout << "Estimating count scale over " << bins << " bins" << endl;
	
	vector<Complex<double>>	ft(bins,Complex<double>(0,0));
	
	// Fourier transform and find the peak
	for ( i=bins/100; i<bins/2; ++i ) {
		for ( j=0; j<100; ++j ) {
			phi = -TWOPI*i*j*1.0L/bins;
			cv = Complex<double>(cos(phi), sin(phi));
			ft[i] += cv*h[j]*1.0L/max;
		}
		ft[i] /= bins;
//		cout << i << tab << h[i] << tab << ft[i].power() << endl;
		if ( hmax < ft[i].power() ) {
			imax = i;
			hmax = ft[i].power();
		}
	}
	
	if ( verbose ) {
		cout << "Peak frequency:                 " << bins/imax << endl;
		cout << "Peak power:                     " << ft[imax].power() << endl;
	}
	
	imax = bins/imax;

	long			np(bins/imax+1), ninc(0), thr(datasize/1e4);
	double			stot(0), inc(0), d;
	double			amp[np], s[np], sv[np], sv2[np];
	
	for ( j=0; j<np; j++ ) s[j] = sv[j] = sv2[j] = amp[j] = 0;
	
	for ( i=0, j=0; i<bins && j<np; i++ ) {
		j = (i+imax/2)/imax;
		if ( j < np ) {
			s[j] += h[i];
			sv[j] += h[i]*i;
			sv2[j] += h[i]*i*i;
			if ( amp[j] < h[i] ) amp[j] = h[i];
		}
	}

	for ( j=0; j<np; j++ ) stot += s[j];
	
	if ( verbose )
		cout << "Peak\t%\tAmp\tAvg\tStDev" << endl;
	for ( j=k=0; j<np; j++ ) if ( s[j] ) {
		sv[j] /= s[j];
		sv2[j] = sv2[j]/s[j] - sv[j]*sv[j];
		if ( j && amp[j] > thr ) {
			inc += sv[j]/j;
			ninc++;
		}
		if ( verbose )
			cout << setprecision(4) << j << tab << 100*s[j]/stot << tab << 
				amp[j] << tab << sv[j] << tab << sqrt(sv2[j]) << endl;
	}
	
	if ( ninc ) {
		inc /= ninc;
		if ( verbose )
			cout << "Counts scale:                   " << inc << " (" << ninc << ")" << endl;
		if ( flags & 2 ) {
			if ( datatype < Float ) scale = inc;
			else scale /= inc;
			for ( i=0; i<datasize; i++ ) set(i, floor((*this)[i]/scale + 0.5));
			statistics();
			if ( verbose )
				cout << "New variance/average:           " << poisson_statistics_check() << endl;
		}
		if ( verbose )
			cout << endl;
	}

	long			ncol(1+np);
	Bstring			title("Peaks");
	Bplot*			plot = NULL;

	if ( flags & 1 ) {
		plot = new Bplot(1, bins, ncol);
		plot->title(title);
	
		plot->page(0).title(title);
		plot->page(0).columns(ncol);
		for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Bin");
		plot->page(0).column(0).axis(1);
		plot->page(0).column(1).type(1);
		plot->page(0).column(1).label("Count");
		plot->page(0).column(1).axis(3);
		plot->page(0).column(1).element_size(0.5);
		plot->page(0).column(1).color(0.5,0.5,0.5);
		for ( j=1, i=2; j<np; j++, i++ ) {
			plot->page(0).column(i).type(2);
			plot->page(0).column(i).label("Count");
			plot->page(0).column(i).axis(3);
			plot->page(0).column(i).element_size(0.5);
			plot->page(0).column(i).color(1,0,0);
		}

		for ( i=0; i<bins; i++ ) {
			(*plot)[i] = i;
			(*plot)[i+bins] = h[i];
			for ( j=1, k=2*bins; j<np; j++, k+=bins )
				if ( sv2[j] ) {
					d = sv[j] - i;
					(*plot)[i+k] = amp[j] * exp(-0.5*d*d/sv2[j]);
				} else
					(*plot)[i+k] = 0;
		}
	}

    return plot;
}
/*
Bplot*		Bimage::histogram_counts(int flags)
{
	long			i, j, k, bins(1000), imax(0), imax2(0);
 	double			scale, offset, hmax(0), hmax2(0);
	
	if ( datatype < Float ) bins = max - min + 1;
	
	vector<long>	h = histogram(bins, scale, offset);

	// Find highest peak
	for ( i=0, imax=0; i<bins; ++i ) if ( hmax < h[i] ) {
		hmax = h[i];
		imax = i;
	}
	hmax2 = hmax/2;
	
	if ( imax < 1 ) {
		cerr << "Warning: First bin in histogram is already the highest!" << endl;
		return NULL;
	}
	
	// Track to start of peak
	for ( i=imax; i && h[i] >= hmax2; --i ) ;
	
	// Track to end of previous peak
	for ( ; i && h[i] < hmax2; --i ) ;
	
	// Find previous peak
	for ( ; i && h[i] > hmax/2; --i ) if ( hmax2 < h[i] ) {
		hmax2 = h[i];
		imax2 = i;
	}
	
	imax -= imax2;
	
	if ( imax < 1 ) {
		cerr << "Warning: First bin in histogram is already the highest!" << endl;
		cerr << "scale = " << scale << " offset = " << offset << endl;
		cerr << "imax = " << imax << " imax2 = " << imax2 << endl;
		return NULL;
	}
	
	long			np(bins/imax+1), ninc(0), thr(datasize/1e4);
	double			stot(0), inc(0), d;
	double			amp[np], s[np], sv[np], sv2[np];
	
	for ( j=0; j<np; j++ ) s[j] = sv[j] = sv2[j] = amp[j] = 0;
	
	for ( i=0, j=0; i<bins && j<np; i++ ) {
		j = (i+imax/2)/imax;
		if ( j < np ) {
			s[j] += h[i];
			sv[j] += h[i]*i;
			sv2[j] += h[i]*i*i;
			if ( amp[j] < h[i] ) amp[j] = h[i];
		}
	}

	for ( j=0; j<np; j++ ) stot += s[j];
	
	if ( verbose )
		cout << "Peak\t%\tAmp\tAvg\tStDev" << endl;
	for ( j=k=0; j<np; j++ ) if ( s[j] ) {
		sv[j] /= s[j];
		sv2[j] = sv2[j]/s[j] - sv[j]*sv[j];
		if ( j && amp[j] > thr ) {
			inc += sv[j]/j;
			ninc++;
		}
		if ( verbose )
			cout << setprecision(4) << j << tab << 100*s[j]/stot << tab <<
				amp[j] << tab << sv[j] << tab << sqrt(sv2[j]) << endl;
	}
	
	if ( ninc ) {
		inc /= ninc;
		cout << "Counts scale:                   " << inc << " (" << ninc << ")" << endl;
		if ( flags & 2 ) {
			if ( datatype < Float ) scale = inc;
			else scale /= inc;
			for ( i=0; i<datasize; i++ ) set(i, floor((*this)[i]/scale + 0.5));
			statistics();
		}
	}

	long			ncol(1+np);
	Bstring			title("Peaks");
	Bplot*			plot = NULL;

	if ( flags & 1 ) {
		plot = new Bplot(1, bins, ncol);
		plot->title(title);
	
		plot->page(0).title(title);
		plot->page(0).columns(ncol);
		for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Bin");
		plot->page(0).column(0).axis(1);
		plot->page(0).column(1).type(1);
		plot->page(0).column(1).label("Count");
		plot->page(0).column(1).axis(3);
		plot->page(0).column(1).element_size(0.5);
		plot->page(0).column(1).color(0.5,0.5,0.5);
		for ( j=1, i=2; j<np; j++, i++ ) {
			plot->page(0).column(i).type(2);
			plot->page(0).column(i).label("Count");
			plot->page(0).column(i).axis(3);
			plot->page(0).column(i).element_size(0.5);
			plot->page(0).column(i).color(1,0,0);
		}

		for ( i=0; i<bins; i++ ) {
			(*plot)[i] = i;
			(*plot)[i+bins] = h[i];
			for ( j=1, k=2*bins; j<np; j++, k+=bins )
				if ( sv2[j] ) {
					d = sv[j] - i;
					(*plot)[i+k] = amp[j] * exp(-0.5*d*d/sv2[j]);
				} else
					(*plot)[i+k] = 0;
		}
	}

    return plot;
}
*/

/**
@brief 	Calculates the percentiles from the histogram of an image.
@return Bplot* 		plot of the percentiles.

	A histogram of an image is calculated with 10000 bins.
	The percentiles are calculated from the running sum of the histogram
	and returned in a plot.

**/
Bplot*		Bimage::percentiles()
{
	long			bins(10000), pbins(101);
    long		  	i, j, hsum;
	double			invsize(100.0/datasize);
 	double			scale, offset, fraction, d;
	
	vector<long>	histo = histogram(bins, scale, offset);

	vector<double>	perc(pbins);

	perc[0] = min;
	perc[100] = max;

	long			ncol(3);
	Bstring			title("Percentiles");
	Bplot*			plot = new Bplot(1, pbins, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("%");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(1);
	plot->page(0).column(1).label("Bin");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).element_size(0.5);
	plot->page(0).column(2).type(2);
	plot->page(0).column(2).label("Value");
	plot->page(0).column(2).axis(3);
	plot->page(0).column(2).element_size(0.5);
	plot->page(0).column(2).color(1,0,0);
//	plot->page(0).axis(1).min(min - bin_half);
//	plot->page(0).axis(1).max(max + bin_half);
//	plot->page(0).axis(3).min(prad->minimum());
//	plot->page(0).axis(3).max(prad->maximum());

	if ( verbose ) {
		cout << "%\tBin\tValue" << endl;
		cout << 0 << tab << 0 << tab << perc[0] << endl;
	}
	(*plot)[2*pbins] = perc[0];
	for ( i=j=1, hsum = histo[0]; i<bins; hsum += histo[i], i++ ) {
		fraction = hsum*invsize;
		for ( d = (fraction>j)? 1.0/(fraction - j): 1; fraction > j; j++ ) {
			perc[j] = min + (i - d*(fraction - j))/scale;
			if ( verbose )
				cout << j << tab << i << tab << perc[j] << endl;
			(*plot)[j] = j;
			(*plot)[pbins+j] = i;
			(*plot)[2*pbins+j] = perc[j];
		}
	}
	if ( verbose )
		cout << 100 << tab << i << tab << perc[100] << endl;
	(*plot)[100] = 100;
	(*plot)[pbins+100] = i;
	(*plot)[2*pbins+100] = perc[100];
	
    return plot;
}

/**
@brief 	Calculates minimum and maximum thresholds for truncation.
@param 	&tmin		minumum threshold.
@param 	&tmax		maxumum threshold.
@return int 		0.

	A histogram of an image is calculated.  The first minimum in the first
	quarter of the histogram and the last minimum in the last quarter
	are taken to define the small and large outliers.

**/
int			Bimage::histogram_minmax(double& tmin, double& tmax)
{
	if ( max <= min ) statistics();
	
	// Calculate the histogram and remove large numbers of white and black pixels
	long			i, j, bins(256);
	if ( (max - min)/std < 25 ) bins = 10*(max - min)/std;
	if ( datatype <= SCharacter ) bins = (long) (max - min + 1.1);
	
	long			ihi, imin, nmin, imax(bins/2), nmax(0), isig(0);
	double			g, sigma(8.0/(bins*bins));
    double			scale, offset;
    
	vector<long>	hist = histogram(bins, scale, offset);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Estimating extreme thresholds from the histogram:" << endl;
		cout << "Bins:                           " << bins << endl;
	}
	
	imax = (avg - min)*scale;
	nmax = hist[imax];
	
	for ( i=imax/2, j=(bins-imax)/2; i<j; i++ ) {
		if ( nmax < hist[i] ) { 	// Find the maximum in the middle part of the histogram	
			nmax = hist[i];
			imax = i;
		}
	}
    
	while ( imax+isig < bins - 1 && imax-isig > 0 && 
			hist[imax] < 2*hist[imax+isig] && hist[imax] < 2*hist[imax-isig] )
		isig++;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::histogram_minmax: Gaussian parameters:" << endl;
		cout << "DEBUG Bimage::histogram_minmax: max=" << nmax << " imax=" << imax << " sigma=" << isig << endl;
	}
	
	if ( isig > 0 ) sigma = -0.5/(isig*isig);
	ihi = 0;
	for ( i=0; i<imax/2; i++ ) {	// Find the highest peak in the first part
		g = nmax*exp(sigma*(imax - i)*(imax - i));
		if ( hist[i] > 20*g ) ihi = i;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::histogram_minmax: highest peak in first quarter = " << ihi << endl;
	
	imin = ihi;
	nmin = x*y*z;
	for ( i=ihi; i<imax/2; i++ ) {	// Find the minimum after the highest peak in the first part
		if ( nmin > hist[i] ) {
			nmin = hist[i];
			imin = i;
		}
	}
    	
	ihi = bins - 1;
	for ( i=bins-1, j=(bins-imax)/2; i>j; i-- ) {	// Find the highest peak in the last part
		g = nmax*exp(sigma*(imax - i)*(imax - i));
		if ( hist[i] > 20*g ) ihi = i;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::histogram_minmax: highest peak in last quarter = " << ihi << endl;
	
	imax = ihi;
	nmin = x*y*z;
	for ( i=(bins-imax)/2; i<ihi; i++ ) {	// Find the minumum before the highest peak in the last part
		if ( nmin >= hist[i] ) {
			nmin = hist[i];
			imax = i;
		}
	}
	
//	delete[] hist;
	
	if ( imin >= imax ) {
		cerr << "Error in Bimage::histogram_minmax: extremes not found!" << endl;
		return -1;
	}
	
	tmin = imin/scale + min;
	tmax = imax/scale + min;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::histogram_minmax: imin=" << imin << " imax=" << imax << " tmin=" << tmin << " tmax=" << tmax << endl;
    	
	if ( verbose & VERB_PROCESS )
		cout << "Thresholds:                      " << tmin << tab << tmax << endl << endl;
    	
	return 0;
}

/**
@brief 	Calculates the inter-set variance of the bisection of a historgram using the method of Otsu.
@param 	bins		number bins in histogram.
@return double 		threshold.

	Reference: NOBUYUKI OTSU, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS, VOL. SMC-9, NO. 1, JANUARY 1979
**/
Bplot*		Bimage::histogram_otsu_variance(long bins)
{
	double			scale, offset, num(x*y*z), num2(num*num);
	vector<long>	h = histogram(bins, scale, offset);
	vector<double>	v = otsu_variance(h);
	
	if ( verbose )
		cout << "Calculating the Otsu inter-set variance with " << bins << " bins" << endl << endl;
	
	long			ncol(3), i, j, k;
	Bstring			title("Otsu inter-set variance");
	Bplot*			plot = NULL;

//	if ( flags & 1 ) {
		plot = new Bplot(1, bins, ncol);
		plot->title(title);
	
		plot->page(0).title(title);
		plot->page(0).columns(ncol);
		for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Bin");
		plot->page(0).column(0).axis(1);
		plot->page(0).column(1).type(1);
		plot->page(0).column(1).label("Probability");
		plot->page(0).column(1).axis(3);
		plot->page(0).column(1).element_size(0.5);
		plot->page(0).column(1).color(0.5,0.5,0.5);
		plot->page(0).column(2).type(2);
		plot->page(0).column(2).label("Variance");
		plot->page(0).column(2).axis(4);
		plot->page(0).column(2).element_size(0.5);
		plot->page(0).column(2).color(1.0,0,0);

		for ( i=0, j=bins, k=2*bins; i<bins; ++i, ++j, ++k ) {
			(*plot)[i] = (i - offset)/scale;
			(*plot)[j] = h[i]/num;
			(*plot)[k] = v[i]/num2;
		}
//	}

    return plot;
}

/**
@brief 	Calculates the threshold from a histogram according to Otsu.
@param 	bins		number bins in histogram.
@return double 		threshold.

	Reference: NOBUYUKI OTSU, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS, VOL. SMC-9, NO. 1, JANUARY 1979
**/
double			Bimage::otsu_threshold(long bins)
{
	double			scale, offset;
	vector<long>	h = histogram(bins, scale, offset);

    long			i, total(x*y*z);
	double			sum(0), sumB(0), wB(0), wF(0), mB, mF, mx(0), bt(0), t1(0), t2(0);
	
    for ( i = 1; i < bins; ++i )
        sum += i * h[i];
	
    for ( i = 0; i < bins; ++i ) {
        wB += h[i];
        if (wB == 0)
            continue;
        wF = total - wB;
        if (wF == 0)
            break;
        sumB += i * h[i];
        mB = sumB / wB;
        mF = (sum - sumB) / wF;
        bt = wB * wF * (mB - mF) * (mB - mF);
        if ( bt >= mx ) {
            t1 = i;
            if ( bt > mx ) {
                t2 = i;
            }
            mx = bt;
        }
    }
	
    return ( t1 + t2 ) / (2 * scale) + min;
}

/**
@brief 	Calculates the inter-set variance of the bisection of a historgram using the method of Otsu.
@param 	h				histogram.
@return vector<double> 		variance.

	Reference: NOBUYUKI OTSU, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS, VOL. SMC-9, NO. 1, JANUARY 1979
**/
vector<double>	Bimage::otsu_variance(vector<long> h)
{
    long			i, bins(h.size()), total(x*y*z);
	double			sum(0), sumB(0), wB(0), wF(0), mB, mF;
	vector<double>	v(bins, 0);
	
    for ( i = 1; i < bins; ++i )
        sum += i * h[i];
	
    for ( i = 0; i < bins; ++i ) {
        wB += h[i];
        wF = total - wB;
		sumB += i * h[i];
        if ( wB && wF ) {
			mB = sumB / wB;
			mF = (sum - sumB) / wF;
			v[i] = wB * wF * (mB - mF) * (mB - mF);
		}
    }
	
    return v;
}

/**
@brief 	Calculates multiple thresholds from a histogram.
@param 	bins			number bins in histogram.
@param 	number			number of clusters (one more than thresholds).
@return vector<double> 		thresholds.

	Reference: PS.Liao, TS.Chen, and PC. Chung,
           Journal of Information Science and Engineering, vol 17, 713-727 (2001)
**/
vector<double>	Bimage::histogram_multi_thresholds(long bins, long number)
{
	double			scale, offset;
	vector<long>	h = histogram(bins, scale, offset);

	vector<double>	t = histogram_thresholds(h, number);
	
	if ( verbose )
		cout << "Thresholds:" << endl;
	for ( auto it=t.begin(); it!=t.end(); ++it ) {
		*it = (*it - offset)/scale;
		if ( verbose )
			cout << *it << endl;
	}

	return t;
}

/**
@brief 	Fits a gaussian function to a histogram of an image.
@param 	bins			number of bins in the histogram.
@param 	ngauss			number of gaussians.
@return vector<double> 		array with 3 values for each gaussian.

	A histogram of an image is calculated with a given number of bins.
	The return has (in order) the gaussian amplitude, the offset, and the sigma value.

**/
vector<double>	Bimage::histogram_gauss_fit(long bins, long ngauss)
{
	change_type(Float);
	
	long			i, j, nh(bins);
	long			thresh((100*datasize)/bins);
	double			hmax(0), havg[ngauss], hstd[ngauss], w[ngauss];

	double			scale, offset;
	vector<long>	histo = histogram(bins, scale, offset);
	
	// Attempting to find peaks by K-means clustering
	vector<long>	sel = k_means(datasize, (float *)data_pointer(), ngauss);
	
	for ( j=0; j<ngauss; j++ ) havg[j] = hstd[j] = w[j] = 0;
	for ( i=0; i<datasize; i++ ) {
		havg[sel[i]] += (*this)[i];
		hstd[sel[i]] += (*this)[i]*(*this)[i];
		w[sel[i]] += 1;
	}
	for ( j=0; j<ngauss; j++ ) {
		havg[j] /= w[j];
		hstd[j] /= w[j];
		hstd[j] -= havg[j]*havg[j];
		if ( hstd[j] > 0 ) hstd[j] = sqrt(hstd[j]);
		else hstd[j] = 0;
	}
	
	sel.resize(bins);
	
	// Selection of bins to fit
	for ( i=nh=0; i<bins; i++ ) {
		if ( histo[i] < thresh ) {
			sel[i] = 1;
			nh++;
		} else {
			sel[i] = 0;
		}
	}
	
	if ( nh < 20 && nh < bins/2 )
		cerr << "Warning: Too few bins selected for fitting! (" << nh << "/" << bins << ")" << endl;

	// 3 parameters: amplitude, position, sigma
	vector<double>	xx(bins), fx(nh);
	for ( i=j=0; i<bins; i++ ) if ( sel[i] ) {
		xx[j] = min + i/scale;
		fx[j] = histo[i];
		if ( hmax < histo[i] ) hmax = histo[i];
		j++;
	}

	Bsimplex	simp(1, 3*ngauss, 1, nh, xx, fx);
	simp.constant(0, ngauss);
	
	for ( i=j=0; j<ngauss; j++, i++ ) {
		simp.parameter(i, hmax);		// Amplitude
		simp.limits(i, datasize/bins, 5*hmax);
//		cout << tab << simp.parameter(i);
		i++;
		simp.parameter(i, havg[j]);		// Position
/*		if ( hstd[j] ) {
			simp.limits(i, havg[j] - 5*hstd[j], havg[j] + 5*hstd[j]);
			if ( simp.limit_low(i) < min ) simp.limit_low(i, min);
			if ( simp.limit_high(i) > max ) simp.limit_high(i, max);
		} else {
			simp.limits(i, min, max);
		}*/
		if ( ngauss == 1 ) {
			simp.limits(i, min, max);
		} else if ( j == 0 ) {
			simp.limits(i, min, (havg[j] + havg[j+1])/2);
		} else if ( j < ngauss-1 ) {
			simp.limits(i, (havg[j-1] + havg[j])/2, (havg[j] + havg[j+1])/2);
		} else {
			simp.limits(i, (havg[j-1] + havg[j])/2, max);
		}
//		cout << tab << simp.parameter(i);
		i++;
		if ( hstd[j] ) {
			simp.parameter(i, hstd[j]);	// Sigma
			simp.limits(i, hstd[j]/3, 3*hstd[j]);
		} else {
			simp.parameter(i, std);				// Sigma
			simp.limits(i, 1/scale, 0.5*(max - min)/n);
		}
//		cout << tab << simp.parameter(i) << endl;
	}

	double			R = simp.run(1000, 0.0001, histogram_gaussian_R);

	vector<double>	gauss;
	for ( j=0; j<3*ngauss; j++ ) gauss.push_back(simp.parameter(j));

	if ( verbose & VERB_PROCESS ) {
		cout << "Gaussian fit over " << nh << " bins:" << endl;
		cout << "Amp\tOffset\tSigma\tSIP" << endl;
		for ( j=0; j<ngauss; j++ )
			cout << gauss[3*j] << tab << gauss[3*j+1] << tab << gauss[3*j+2] << tab << bins*gauss[3*j]/datasize << endl;
		cout << "Sampling-independent RMSD:      " << bins*R/datasize << endl << endl;
	}

	image->FOM(R);

	return gauss;
}

vector<double>	Bimage::histogram_gauss_fit2(long bins, long ngauss)
{
	change_type(Float);
	
	long			i, j, k, l;
	double			hmax(0);
	double			scale, offset;
	vector<long>	histo = histogram(bins, scale, offset);
	vector<double>	mh(bins);
	for ( i=0; i<bins; i++ ) mh[i] = histo[i];

	long			w(bins/20);
	if ( w < 5 ) w = 5;
	double			damp(3/scale);
	
	// Attempting to find peaks based on the gradient change
	mh = moving_gradient(mh, w);
	
	vector<double>	mx(ngauss);
	double			ex;
	for ( i=j=k=l=0, ex=0; i<bins; i++ ) {
		if ( ex < mh[i] ) {
			ex = mh[i];
			j = i;
		} else if ( mh[i] < 0 && mh[j] > 10 ) {
			mx[k++] = i;
			ex = 0;
			j = i;
		}
	}
	
	for ( auto it=mx.begin(); it!=mx.end(); ++it )
		*it = (*it - offset)/scale;

//	for ( auto it=mx.begin(); it!=mx.end(); ++it )
//		cout << *it << endl;

	// 3 parameters: amplitude, position, sigma
	vector<double>	xx(bins), fx(bins);
	for ( i=j=0; i<bins; i++ ) {
		xx[j] = min + i/scale;
		fx[j] = histo[i];
		if ( hmax < histo[i] ) hmax = histo[i];
		j++;
	}

	Bsimplex	simp(1, 3*ngauss, 1, bins, xx, fx);
	simp.constant(0, ngauss);

	for ( i=j=0; j<ngauss; j++, i++ ) {
		simp.parameter(i, hmax);		// Amplitude
		simp.limits(i, datasize/bins, 5*hmax);
//		cout << tab << simp.parameter(i);
		i++;
		simp.parameter(i, mx[j]);		// Position
		if ( ngauss == 1 ) {
			simp.limits(i, min+damp, max-damp);
		} else if ( j == 0 ) {
			simp.limits(i, min+damp, (mx[0]+mx[1])/2);
		} else if ( j < ngauss-1 ) {
			simp.limits(i, (mx[j-1]+mx[j])/2, (mx[j]+mx[j+1])/2);
		} else {
			simp.limits(i, (mx[0]+mx[1])/2, max-damp);
		}

		i++;
		simp.parameter(i, std/ngauss);	// Sigma
		simp.limits(i, damp, std);
//		cout << tab << simp.parameter(i) << endl;
	}

	double			R = simp.run(5000, 1e-6, histogram_gaussian_R);
	R = bins*R/datasize;
	
	vector<double>	gauss;
	for ( j=0; j<3*ngauss; j++ ) gauss.push_back(simp.parameter(j));
	gauss.push_back(R);

	if ( verbose & VERB_PROCESS ) {
		cout << "Gaussian fit over " << bins << " bins:" << endl;
		cout << "Amp\tOffset\tSigma\tSIP" << endl;
		for ( j=0; j<ngauss; j++ )
			cout << gauss[3*j] << tab << gauss[3*j+1] << tab << gauss[3*j+2] << tab << bins*gauss[3*j]/datasize << endl;
		cout << "Sampling-independent RMSD:      " << R << endl << endl;
	}

	image->FOM(R);

	return gauss;
}

/**
@brief 	Fits a gaussian function to a histogram of an image.
@param 	bins		number of bins in the histogram.
@param 	ngauss		number of gaussians.
@return Bplot* 		plot of the histogram and the gaussian fit.

	A histogram of an image is calculated with a given number of bins.

**/
Bplot* 		Bimage::histogram_gauss_plot(long bins, long ngauss)
{
	long			i, j, k, m, ncol(ngauss+2);
	double			v, bin_half = (max - min)/(4*bins);

	vector<double>	gauss = histogram_gauss_fit2(bins, ngauss);

	double			scale, offset;
	vector<long>	histo = histogram(bins, scale, offset);
	
	// Set the average and standard deviation to the first gaussian parameters
	average(gauss[1]);
	standard_deviation(gauss[2]);

	Bstring			title("Gaussian fit"), txt;
	Bplot*			plot = new Bplot(1, bins, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
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
	plot->page(0).column(2).color(1,0,0);
	if ( n > 1 ) plot->page(0).column(3).color(0,1,0);
	if ( n > 2 ) plot->page(0).column(4).color(0,0,1);
	plot->page(0).axis(1).min(min - bin_half);
	plot->page(0).axis(1).max(max + bin_half);
	
	for ( j=0; j<ngauss; j++ ){
		txt = Bstring(j+1, "Gaussian %d:") + Bstring(gauss[3*j], " a=%g") +
			Bstring(gauss[3*j+1], " u=%g") + Bstring(gauss[3*j+2], " s=%g");
		plot->page(0).add_text(txt);
	}
	txt = Bstring(image->FOM(), "R: %g");
	plot->page(0).add_text(txt);
	
	for ( i=0, j=bins; i<bins; i++, j++ ) {
		(*plot)[i] = min + i/scale;
		(*plot)[j] = histo[i];
		for ( k=0, m=i+2*bins; k<3*ngauss; k+=3, m+=bins ) {
			v = ((*plot)[i] - gauss[k+1])/gauss[k+2];
			(*plot)[m] = gauss[k]*exp(-0.5*v*v);
		}
	}

//	delete[] histo;
	
	return plot;
}
/*
double		simplex_poisson_R(Bsimplex& simp)
{
	long		i, j, n;
	double		v, R(0);
//	double*		df = new double[simp.points()];
	vector<double>	df(simp.points());
	
	simp.parameter(1, floor(simp.parameter(1)+0.5));	// Enforce integer
	
	for ( i=0, n=0; i<simp.points(); i++ ) {
		if ( simp.fx[i] >= 0 ) {
			v = simp.x[i] + simp.parameter(1);
			df[i] = fabs(simp.fx[i] - simp.parameter(0) - v*log(simp.parameter(2)) +
				simp.parameter(2) + lgamma(v+1.0L));
			if ( !isfinite(df[i]) ) {
				error_show("Error in img_fit_poisson_to_his_R", __FILE__, __LINE__);
				cerr << "i=" << i << "\tfx[i]=" << simp.fx[i] << "\tamp=" << 
					simp.parameter(0) << "\tpos=" << simp.parameter(1)
					<< "\tlambda=" << simp.parameter(2) << endl;
				return 1e37;
			}
			n++;
		} else df[i] = -1;
	}
	
//	qsort((void *) df, simp.points(), sizeof(double), 
//		(int (*)(const void *, const void *)) QsortSmallToLargeDouble);
	sort(df.begin(), df.end());
	
//	if ( simp->npoint < 100 ) j = simp->npoint - 3;
//	else j = (simp->npoint*97)/100;
	j = simp.points() - 3;
	
	for ( i=0, R=0; i<j; i++ ) if ( df[i] > -1 ) {
		R += df[i]*df[i];
		n++;
	}
	
	if ( n ) R = sqrt(R/n);
	else R = 1e37;
			
//	delete[] df;
	
	return R;
}
*/
double		simplex_poisson_R(Bsimplex& simp)
{
	long			i;
	double			fact(1), pl(1), df, R(0);
	vector<double>&	f = simp.dependent_values();
	double			ep(simp.parameter(0)*exp(-simp.parameter(1)));
	
//	simp.parameter(1, floor(simp.parameter(1)+0.5));	// Enforce integer
	
	for ( i=0; i<simp.points(); i++ ) {
		if ( i ) {
			fact *= i;
			pl *= simp.parameter(1);
		}
		df = f[i] - pl*ep/fact;
		R += df*df;
	}
	
	if ( i ) R = sqrt(R/i);
	
	return R;
}

/**
@brief 	Fits a poisson function to a histogram of an image.
@param 	bins		number of bins in the histogram.
@param	flag		1=plot.
@return Bplot* 		plot of the histogram and the poisson fit.

	A histogram of an image is calculated with a given number of bins.
	The histogram is fit to the Poisson function:
		           lambda * exp(-lambda)
		f(x) = Amp ---------------------
		                    x!
	by first converting it to a linear form:
		ln(f(x)) = ln(Amp) + x*ln(lambda) - lambda - lgamma(x+1)
	and then using the downhill simplex method to find the solution.

**/
/*Bplot* 		Bimage::histogram_poisson_fit(long bins, int flag)
{
	long			i, j, m, nn, ni;
	double			hmax(0), sum(0);
	double			bin_half = (max - min)/(4*bins);

	double			scale, offset;
	long*			histo = histogram(bins, scale, offset);
	
	// Count the number of positive bins and readjust the histogram to have only positive bins
	for ( i=0, j=0, nn=0, ni=0; i<bins; i++ ) {
		if ( histo[i] > 0 ) {
			ni += j;
			nn++;
			j = 0;
		}
		j++;
	}
//	if ( verbose & VERB_FULL )
		cout << "Positive bins = " << nn << "   Interval size = " << ni*1.0/nn << endl;
	delete[] histo;
	
	if ( bins < nn ) bins = nn;
	else bins = (bins*nn)/ni;

	histo = histogram(bins, scale, offset);
	
//	if ( verbose & VERB_FULL )
		cout << "Bins reset to " << bins << ", scale = " << scale << ", offset = " << offset << endl << endl;

	// 3 parameters: amplitude, position, sigma
//	double*		x = new double[bins];
//	double*		fx = new double[bins];
	vector<double>	xx(bins), fx(bins);
	for ( i=0; i<bins; i++ ) xx[i] = min + i/scale;
	for ( i=1; i<bins-1; i++ ) {
		if ( histo[i] < 1 ) fx[i] = log(0.5);
		else fx[i] = log(1.0*histo[i]);
		sum += fx[i];
		if ( hmax < fx[i] ) hmax = fx[i];
		cout << i << tab << xx[i] << tab << histo[i] << tab << fx[i] << endl;
	}
	fx[0] = fx[1];
	fx[bins-1] = log(0.5);

	Bsimplex	simp(1, 3, 0, bins, xx, fx);
	
	simp.parameter(0, log(datasize*1.0));				// ln(Amplitude)
	simp.limits(0, log(datasize/2.0), log(datasize+1.0));

	simp.parameter(0, 0);							// Position in bins
	simp.limits(1, 0, 0);
	
	simp.parameter(2, std*std);						// Lambda in bins
	if ( simp.parameter(2) < 0.001 ) simp.parameter(2, 0.001);
	simp.limits(2, 0.001, 10*std);
	if ( simp.hi[2] < 1 ) simp.hi[2] = 10*avg;
	
	double			R = simp.run(1000, 0.001, simplex_poisson_R);
	R = exp(R);

//	delete[] x;
//	delete[] fx;
	
	double*			poisson = new double[3];
	poisson[0] = exp(simp.parameter(0));
	poisson[1] = simp.parameter(1);
	poisson[2] = simp.parameter(2);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Poisson function fit over " << bins << " bins:" << endl;
		cout << "Lambda:                         " << poisson[2] << endl;
//		cout << "Offset:                         " << poisson[1] << endl;
		cout << "Included voxels:                " << poisson[0] << " (" << 100*poisson[0]/datasize << "%)" << endl;
		cout << "Sampling-independent RMSD:      " << bins*R/datasize << endl << endl;
	}

	// Set the average and standard deviation to the poisson parameters
	average(poisson[2]);
	standard_deviation(sqrt(poisson[2]));
	image->FOM(R);

	long			ncol(3);
	Bstring			title("Poisson fit"), txt;
	Bplot*			plot = NULL;

	if ( flag ) {
		plot = new Bplot(1, bins, ncol);
		plot->title(title);
	
		plot->page(0).title(title);
		plot->page(0).columns(ncol);
		for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Bin");
		plot->page(0).column(0).axis(1);
		plot->page(0).column(1).type(1);
		plot->page(0).column(1).label("Count");
		plot->page(0).column(1).axis(3);
		plot->page(0).column(1).element_size(0.5);
		plot->page(0).column(1).color(0.5,0.5,0.5);
		plot->page(0).column(2).type(2);
		plot->page(0).column(2).label("Count");
		plot->page(0).column(2).axis(3);
		plot->page(0).column(2).element_size(0.5);
		plot->page(0).column(2).color(1,0,0);
		plot->page(0).axis(1).min(min - bin_half);
		plot->page(0).axis(1).max(max + bin_half);
//		plot->page(0).axis(3).min(prad->minimum());
//		plot->page(0).axis(3).max(prad->maximum());
	
		txt = Bstring(poisson[0], "Poisson: a=%g") + Bstring(poisson[2], " u=%g");
		plot->page(0).add_text(txt);
		txt = Bstring(R, "R: %g");
		plot->page(0).add_text(txt);
	
		for ( i=0, j=bins, m=2*bins; i<bins; i++, j++, m++ ) {
			(*plot)[i] = i;
			(*plot)[j] = histo[i];
			(*plot)[m] = poisson[0]*pow(poisson[2],(double)i)*exp(-poisson[2])/factorial(i);
			if ( !isfinite((*plot)[m]) ) (*plot)[m] = 0;
		}
	}

	delete[] histo;
	delete[] poisson;
		
	return plot;
}
*/
Bplot* 		Bimage::histogram_poisson_fit(long bins, int flag)
{
	if ( std <= 0 ) statistics();
	
	if ( bins > avg + 5*std ) bins = avg + 5*std;
	if ( bins > 150 ) bins = 150;
	if ( bins < 5 ) bins = 5;
	
	long			i, j, m, imax(0);
	double			hmax(0), scale(bins*1.0/(avg+5*std));
		
	vector<double>	b(bins, 0);
	vector<double>	h(bins, 0);
	
	for ( i=0; i<datasize; ++i ) {
		j = (long) (scale*(*this)[i] + 0.5);
		if ( j >= 0 && j < bins ) h[j]+=1;
	}
	
	for ( i=0; i<bins; ++i ) {
		b[i] = i;
		if ( hmax < h[i] ) {
			hmax = h[i];
			imax = i;
		}
//		cout << i << tab << b[i] << tab << h[i] << endl;
	}
	if ( imax < 1 ) imax = 1;
	
//	bexit(0);

	cout << "scale=" << scale << " imax=" << imax << endl;
	
	Bsimplex	simp(1, 2, 0, bins, b, h);
	
	simp.parameter(0, hmax);			// Amplitude
	simp.limits(0, hmax, 20*hmax);

	simp.parameter(1, imax);			// Lambda in bins
	if ( avg < imax ) simp.parameter(1, avg);
//	simp.limits(1, 0.5*imax, 2*imax);
	simp.limits(1, avg/2, 2*imax);
	
	double			R = simp.run(1000, 0.001, simplex_poisson_R);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Poisson function fit over " << bins << " bins:" << endl;
		cout << "Amplitude:                      " << simp.parameter(0) << endl;
		cout << "Lambda:                         " << simp.parameter(1) << endl;
		cout << "Sampling-independent RMSD:      " << bins*R/datasize << endl << endl;
	}

	// Set the average and standard deviation to the poisson parameters
	average(simp.parameter(1));
	standard_deviation(sqrt(simp.parameter(1)));
	image->FOM(R);

	long			ncol(3);
	Bstring			title("Poisson fit"), txt;
	Bplot*			plot = NULL;

	if ( flag ) {
		plot = new Bplot(1, bins, ncol);
		plot->title(title);
	
		plot->page(0).title(title);
		plot->page(0).columns(ncol);
		for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Bin");
		plot->page(0).column(0).axis(1);
		plot->page(0).column(1).type(1);
		plot->page(0).column(1).label("Count");
		plot->page(0).column(1).axis(3);
		plot->page(0).column(1).element_size(0.5);
		plot->page(0).column(1).color(0.5,0.5,0.5);
		plot->page(0).column(2).type(2);
		plot->page(0).column(2).label("Count");
		plot->page(0).column(2).axis(3);
		plot->page(0).column(2).element_size(0.5);
		plot->page(0).column(2).color(1,0,0);
		plot->page(0).axis(1).min(0);
		plot->page(0).axis(1).max(bins);
//		plot->page(0).axis(3).min(prad->minimum());
//		plot->page(0).axis(3).max(prad->maximum());
	
		txt = Bstring(simp.parameter(0), "Poisson: a=%g") + Bstring(simp.parameter(1), " u=%g");
		plot->page(0).add_text(txt);
		txt = Bstring(R, "R: %g");
		plot->page(0).add_text(txt);
	
		for ( i=0, j=bins, m=2*bins; i<bins; i++, j++, m++ ) {
			(*plot)[i] = i;
			(*plot)[j] = h[i];
			(*plot)[m] = simp.parameter(0)*pow(simp.parameter(1),(double)i)*exp(-simp.parameter(1))/factorial(i);
			if ( !isfinite((*plot)[m]) ) (*plot)[m] = 0;
		}
	}

	return plot;
}

