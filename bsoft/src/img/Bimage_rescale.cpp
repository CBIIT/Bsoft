/**
@file	Bimage_rescale.cpp
@brief	Library routines to rescale images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210110
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Rescales the image data with a given multiplier and offset.
@param 	scale		multiplier.
@param 	shift		addition or offset.
@return int			error code.

	The new data is calculated as:
		new_datum = datum*scale + shift
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int			Bimage::rescale(double scale, double shift)
{
	int				err(0);
    long		   	i;
    double			v1;
   
	if ( verbose & VERB_FULL )
    	cout << "Scale and shift:                " << scale << " " << shift << endl << endl;
	
	// Adjust the background values
	for ( i=0; i<n; i++ )
		background(i, background(i)*scale + shift);
//		image[i].background(image[i].background()*scale + shift);

	for ( i=0; i<datasize; i++ ) {
		v1 = (*this)[i]*scale + shift;
		set(i, v1);
	}
	
	if ( ( err = statistics() ) )
		cerr << tab << "in Bimage::rescale" << endl;
	
	return err;
}

/**
@brief 	Rescales the image data with a given multiplier and offset.
@param	nn			sub-image.
@param 	scale		multiplier.
@param 	shift		addition or offset.
@return int			error code.

	Requirement: The images must have the same size.

**/
int			Bimage::rescale(long nn, double scale, double shift)
{
	int				err(0);
	long			i, j, imgsize(c*image_size());
	double			v1;
	
	if ( verbose & VERB_FULL )
    	cout << "Scale and shift:                " << scale << " " << shift << endl << endl;
	
	// Adjust the background value
//	image[nn].background(image[nn].background()*scale + shift);
	image[nn].background(background(nn)*scale + shift);

	for ( i=0, j=nn*imgsize; i<imgsize; i++, j++ ) {
		v1 = (*this)[j] * scale + shift;
		set(j, v1);
	}
	
//	if ( ( err = statistics() ) )
//		cerr << tab << "in Bimage::rescale" << endl;
	
	return err;
}

/**
@brief 	Rescales the image data to a given minimum and maximum.
@param 	numin		new minimum.
@param 	numax		new maximum.
@return int			error code.

	The new data is calculated as:
		                            new_max - new_min
		new_datum = (datum - min) * ----------------- + new_min
		                                max - min
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::rescale_to_min_max(double numin, double numax)
{
	if ( numin >= numax ) {
		cerr << "Error in rescaling: minimum (" << numin << ") is greater than maximum ("
			<< numax << ") for " << file_name() << endl;
		return -1;
	}
	
    double			scale = (numax - numin)/(max - min);
	double			shift = numin - min*scale;
    
	if ( verbose & VERB_FULL )
	    cout << "Rescaling to:                   " << numin << " " << numax << endl;
	
	return rescale(scale, shift);
}

/**
@brief 	Rescales the image data to a given minimum and maximum.
@param	nn			sub-image.
@param 	numin		new minimum.
@param 	numax		new maximum.
@return int			error code.

	The new data is calculated as:
		                            new_max - new_min
		new_datum = (datum - min) * ----------------- + new_min
		                                max - min
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::rescale_to_min_max(long nn, double numin, double numax)
{
	if ( numin >= numax ) {
		cerr << "Error in rescaling: minimum (" << numin << ") is greater than maximum ("
			<< numax << ") for " << file_name() << endl;
		return -1;
	}
	
    double			scale = (numax - numin)/(image[nn].maximum() - image[nn].minimum());
	double			shift = numin - image[nn].minimum()*scale;
    
	if ( verbose & VERB_FULL )
	    cout << "Rescaling to:                   " << numin << " " << numax << endl;
	
	return rescale(nn, scale, shift);
}

/**
@brief 	Rescales the image data to a given average and standard deviation.
@param 	nuavg		new average.
@param 	nustd		new standard deviation.
@return int			error code.

	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::rescale_to_avg_std(double nuavg, double nustd)
{
	if ( nustd < 0 ) {
		cerr << "Warning: Cannot use a negative standard deviation to scale to! (" << nustd << ")" << endl;
		return -1;
	}
	
	if ( std < 1e-30 ) if ( statistics() ) {
		cerr << tab << "in Bimage::rescale_to_avg_std" << endl;
		return -1;
	}
	
    double			scale = (std>1e-30)? nustd/std: 1;
	double			shift = nuavg - avg*scale;
    
	if ( verbose & VERB_FULL )
	    cout << "Rescaling to average and stdev: " << nuavg << " " << nustd << endl;
	
	return rescale(scale, shift);
}

/**
@brief 	Rescales the image data to a given average and standard deviation.
@param	nn			sub-image.
@param 	nuavg		new average.
@param 	nustd		new standard deviation.
@return int			error code.

	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::rescale_to_avg_std(long nn, double nuavg, double nustd)
{
	if ( nustd < 0 ) {
		cerr << "Warning: Cannot use a negative standard deviation to scale to! (" << nustd << ")" << endl;
		return -1;
	}
	
	if ( std < 1e-30 ) if ( statistics() ) {
		cerr << tab << "in Bimage::rescale_to_avg_std" << endl;
		return -1;
	}
	
    double			scale = (image[nn].standard_deviation()>1e-30)? nustd/image[nn].standard_deviation(): 1;
	double			shift = nuavg - image[nn].average()*scale;
    
	if ( verbose & VERB_FULL )
	    cout << "Rescaling to average and stdev: " << nuavg << " " << nustd << endl;
	
	return rescale(nn, scale, shift);
}

/**
@brief 	Rescales the image data to a given average and standard deviation.
@param 	nuavg		new average.
@param 	nustd		new standard deviation.
@param	*pmask		statistical calculations limited to the masked region.
@return int			error code.

	The new data is calculated as:
		                            new_std_dev
		new_datum = (datum - avg) * ----------- + new_avg
		                              std_dev
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::rescale_to_avg_std(double nuavg, double nustd, Bimage* pmask)
{
	if ( nustd < 0 ) {
		cerr << "Warning: Cannot use a negative standard deviation to scale to! (" << nustd << ")" << endl;
		return -1;
	}
	
	double			regavg(0), regstd(0);
	
	statistics(pmask, regavg, regstd);
	
	if ( regstd == 0 ) {
		cerr << "Error in Bimage::rescale_to_avg_std: Region standard deviation is zero for "
			<< file_name() << "!" << endl;
		return -1;
	}
	
    double			scale = nustd/regstd;
	double			shift = nuavg - regavg*scale;
    
	if ( verbose & VERB_FULL )
	    cout << "Rescaling to average and stdev: " << nuavg << " " << nustd << endl;
	
	return rescale(scale, shift);
}

/**
@brief 	Truncates image data to a given minimum and maximum.
@param 	minim		minimum.
@param 	maxim		maximum.
@param 	setmin		value to set voxels smaller than minimum.
@param 	setmax		value to set voxels larger than maximum.
@return int 		0.

	All values smaller than the new minimum are set to the new minimum
	and all values larger than the new maximum are set to the new maximum.
	In cases where the given minimum and maximum are outside the data
	type ranges, they are reset to the data type range limits.
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::truncate(double minim, double maxim, double setmin, double setmax)
{
	if ( !data_pointer() ) return -1;
	
    long   			i;
    double  	    v1;
    
	// Make sure the extremes don't exceed the data
	if ( setmin < min ) setmin = min;
	if ( setmax > max ) setmax = max;
	if ( minim < dtmin ) minim = dtmin;
	if ( maxim > dtmax ) maxim = dtmax;
	if ( setmin < dtmin ) setmin = dtmin;
	if ( setmax > dtmax ) setmax = dtmax;
	
	// Adjust the background values
	for ( i=0; i<n; i++ ) {
//		if ( image[i].background() < minim ) image[i].background(setmin);
//		if ( image[i].background() > maxim ) image[i].background(setmax);
		if ( background(i) < minim ) background(i, setmin);
		if ( background(i) > maxim ) background(i, setmax);
	}
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Truncating to:                  " << minim << " " << maxim << endl;
	    cout << "Min and max replacement values: " << setmin << " " << setmax << endl << endl;
	}
	
    for ( i=0; i<datasize; i++ ) {
 		v1 = (*this)[i];
		if ( v1 < minim ) v1 = setmin;
		if ( v1 > maxim ) v1 = setmax;
		set(i, v1);
	}

	statistics();

	return 0;
}

/**
@brief 	Truncates image data to a given minimum and maximum.
@param 	minim		minimum.
@param 	maxim		maximum.
@return int 		0.

	All values smaller than the new minimum are set to the new minimum
	and all values larger than the new maximum are set to the new maximum.
	In cases where the given minimum and maximum are outside the data
	type ranges, they are reset to the data type range limits.
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::truncate_to_min_max(double minim, double maxim)
{
	return truncate(minim, maxim, minim, maxim);
}

/**
@brief 	Sets voxels in image data exceeding a given minimum and maximum to
	the average.
@param 	minim		minimum.
@param 	maxim		maximum.
@return int 			0.

	All values smaller than the new minimum or larger than the new 
	maximum are set to the image average.
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::truncate_to_avg(double minim, double maxim)
{
	return truncate(minim, maxim, avg, avg);
}

/**
@brief 	Sets voxels in image data exceeding a given minimum and maximum to
	the image background.
@param 	minim		minimum.
@param 	maxim		maximum.
@return int 		0.

	All values smaller than the new minimum or larger than the new 
	maximum are set to the background value.
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::truncate_to_background(double minim, double maxim)
{
	if ( !data_pointer() ) return -1;
	
    long   			i, nn, n1, imgsize = x*y*z*n;
    double  	    v1, bkg;
    
	// Make sure the extremes don't exceed the data
	if ( minim < dtmin ) minim = dtmin;
	if ( maxim > dtmax ) maxim = dtmax;
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Truncating to background:" << endl;
	    cout << "Range:                          " << minim << " " << maxim << endl;
	    cout << "Image\tBackground" << endl;
	}
	
	// Adjust the background values
	for ( nn=0; nn<n; nn++ ) {
//		if ( image[nn].background() < minim ) image[nn].background(minim);
//		if ( image[nn].background() > maxim ) image[nn].background(maxim);
//		if ( verbose & VERB_PROCESS )
//			cout << nn+1 << tab << image[nn].background() << endl;
		bkg = background(nn);
		if ( bkg < minim ) image[nn].background(minim);
		if ( bkg > maxim ) image[nn].background(maxim);
		if ( verbose & VERB_PROCESS )
			cout << nn+1 << tab << background(nn) << endl;
	}
	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	for ( nn=0; nn<n; nn++ ) {
//		bkg = image[nn].background();
		bkg = background(nn);
		for ( i=nn*imgsize, n1=nn+imgsize; i<n1; i++ ) {
			v1 = (*this)[i];
			if ( v1 < minim ) v1 = bkg;
			if ( v1 > maxim ) v1 = bkg;
			set(i, v1);
		}
	}
    
	statistics();

	return 0;
}

/**
@brief 	Converts a full gray scale image to a limited level image.
@param 	nlevels		number of levels.
@return int 		0.

	The dynamic range of the image is decreased to a given number of gray
	scale levels.
	The new data replaces the old data.
	Image statistics are recalculated.

**/
int 		Bimage::limit_levels(int nlevels)
{
	if ( !data_pointer() ) return -1;
	
	if ( nlevels > 255 ) {
		cerr << "Error: There should be less than 255 gray levels (" << nlevels << ")" << endl;
		return -1;
	}

	if ( verbose & VERB_PROCESS )
	    cout << "Limiting levels to:             " << nlevels << endl << endl;
	
    long   i;
    double			unit = 255.0/(nlevels - 1);
    double			scale = 0.999*nlevels/(max - min);
	
    unsigned char* 	nudata = new unsigned char[datasize];

	if ( verbose & VERB_LABEL )
	    cout << "Restricting to " << nlevels << " levels." << endl << endl;
	
	for ( i=0; i<datasize; i++ )
		nudata[i] = (unsigned char) (unit*floor(((*this)[i] - min)*scale));
    
	data_type(UCharacter);

	data_assign((unsigned char *) nudata);
	
	statistics();
	
	return 0;
}

/**
@brief 	Normalizes a sub-image to a desired average and standard deviation.
@param	imgnum		sub-image number.
@param 	average		desired average.
@param 	stdev		desired standard deviation (if 0, use defaults).
@param 	norm_type	type of determining the effective average and standard deviation:
 						0=simple, 1=Gaussian, 2=Poisson.
@param	bins		number of histogram bins required to fit distributions.
@return int			0.

	The effective average and standard deviation for each image is obtained
	in one of three ways:
		0.		The simple avergae and standard devaition for the image.
		1.		Gaussian fit of the histogram.
		2.		Poisson fit of the histogram.
	A histogram of an image is calculated with a given number of bins.
	The histogram is fit to a Gaussian or Poisson function with exclusion of a
	small number of bins in the histogram (defined as outliers).
	The effective average and standard deviation are used to 
	rescale the data for each image.

**/
int			Bimage::normalize(long imgnum, double average, double stdev, int norm_type, long bins)
{
	Bplot*          plot = NULL;
	
	Bimage*			ptemp = extract(imgnum);

	ptemp->statistics();
		
	ptemp->truncate_to_min_max(ptemp->average() - 10*ptemp->standard_deviation(), ptemp->average() + 10*ptemp->standard_deviation());
	
	switch ( norm_type ) {
		case 1:
			plot = ptemp->histogram_gauss_plot(bins, 1);
			break;
		case 2:
			plot = ptemp->histogram_poisson_fit(bins, 1);
			break;
		default:
			ptemp->statistics();
	}
	
	ptemp->rescale_to_avg_std(average, stdev);
	
//	if ( img_rescale_to_avg_std(ptemp, average, stdev) )
//		cerr << "Error: The extrema for image " << imgnum+1 << " are the same: " <<
//			ptemp->minimum() << " " << ptemp->maximum() << endl;
	
	replace(imgnum, ptemp);
	
	delete ptemp;
	if ( plot ) delete plot;

	return 0;
}

/**
@brief 	Normalizes a set of images to a desired average and standard deviation.
@param 	average		desired average.
@param 	stdev		desired standard deviation (if 0, use defaults).
@param 	norm_type	type of determining the effective average and standard deviation:
 						0=simple, 1=Gaussian, 2=Poisson.
@return int			0.

	The effective average and standard deviation for each image is obtained
	in one of three ways:
		0.		The simple avergae and standard devaition for the image.
		1.		Gaussian fit of the histogram.
		2.		Poisson fit of the histogram.
	A histogram of an image is calculated with a given number of bins.
	The histogram is fit to a Gaussian or Poisson function with exclusion of a
	small number of bins in the histogram (defined as outliers).
	The effective average and standard deviation are used to 
	rescale the data for each image.

**/
int			Bimage::normalize(double average, double stdev, int norm_type)
{
	int				bins = (int) max;
	if ( bins < 256 ) bins = 256;
	if ( bins > 1024 ) bins = 1024;
	
	if ( stdev <= 0 ) {			// Defaults
		switch ( datatype ) {
			case UCharacter: average = 127; stdev = 40; break;
			case SCharacter: average = 0; stdev = 40; break;
			case UShort: average = USHRT_MAX/2; stdev = 256; break;
			case Short: average = 0; stdev = 256; break;
			case UInteger: average = 100000; stdev = 1000; break;
			case Integer: average = 0; stdev = 1000; break;
			case ULong: average = 100000; stdev = 1000; break;
			case Long: average = 0; stdev = 1000; break;
			default: average = 0; stdev = 1; break;
		}
	}

	change_type(Float);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Mass normalizing ";
		if ( file_name().length() ) cout << file_name() << endl;
		else cout << endl;
		cout << "Normalizing type:               ";
		switch ( norm_type ) {
			case 1: cout << "Gaussian"; break;
			case 2: cout << "Poisson"; break;
			default: cout << "Simple";
		}
		cout << endl << "Average and standard deviation: " << average << " " << stdev << endl;
	}

#ifdef HAVE_GCD
//	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		normalize(nn, average, stdev, norm_type, bins);
//		dispatch_sync(myq, ^{
//			if ( verbose & VERB_FULL )
//				cout << "Image " << nn << " normalized" << endl;
//		});
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; nn++ ) {
		normalize(nn, average, stdev, norm_type, bins);
//	#pragma omp critical
//		{
//			if ( verbose & VERB_FULL )
//				cout << "Image " << nn << " normalized" << endl;
//		}
	}
#endif

	calculate_background();
	statistics();

	return 0;
}

/**
@brief 	Normalizes by subtracting local average and dividing by local standard deviation.
@param 	kernel_size	size of kernel edge.
@return int 		0.

	The local average and standard deviation within a kernel is calculated 
	and used to normalize the image.
	The convolution is threaded if compiled with GCD or OpenMP.

**/
int			Bimage::normalize_local(long kernel_size)
{
	Vector3<long>		kernel(kernel_size, kernel_size, kernel_size);
	
	return normalize_local(kernel);
}

/**
@brief 	Normalizes by subtracting local average and dividing by local standard deviation.
@param 	kernel		size of kernel edge.
@return int 		0.

	The local average and standard deviation within a kernel is calculated 
	and used to normalize the image.
	The convolution is threaded if compiled with GCD or OpenMP.

**/
int			Bimage::normalize_local(Vector3<long> kernel)
{
	change_type(Float);
	
	kernel = kernel.max(1);
	kernel = kernel.min(size());

	Bimage*				pc = copy();
	
	if ( verbose )
		cout << "Normalizing with a kernel size of " << kernel.volume() << endl << endl;

	long			i, nn;
	pc->next = new Bimage(Float, TSimple, size(), n);
	
	// Set up the square block
	for ( i=0; i<datasize; i++ ) pc->next->set(i, (*pc)[i] * (*pc)[i]);

	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn << endl;
		for ( i=0; i<3; i++ )
			pc->line_sums(nn, i, kernel[i]);
	}

	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	double		v, s;
	for ( i=0; i<datasize; i++ ) {
		v = (*pc)[i];
		s = (*pc->next)[i] - v*v;
		if ( s > 0 ) {
			s = sqrt(s);
			v = ((*this)[i] - v)/s;
		} else v = 0;
		set(i, v);
	}
	
	delete pc;

	statistics();
	
	return 0;
}

/**
@brief Calculates the square of an image.

	Values less than zero are set to zero.
	Values greater than the data type maximum are set to the maximum.

**/
void		Bimage::square()
{
	double			v;
	
	for ( long j=0; j<datasize; j++ ) {
		v = (*this)[j];
		v *= v;
		set(j, v);
	}
	
	statistics();
}

/**
@brief Calculates the square root of an image.

	Values less than zero are set to zero.
	Values greater than the data type maximum are set to the maximum.

**/
void		Bimage::square_root()
{
	double			v;
	
	for ( long j=0; j<datasize; j++ ) {
		v = (*this)[j];
		if ( v > 0 ) v = sqrt(v);
		else v = 0;
		set(j, v);
	}
	
	statistics();
}

/**
@brief 	Calculates the logarithm of the image data.

	The image is first converted to floating point, or intensities
	for complex data. The logarithm is calculated to place the minmum at
	zero and scale it to the standard deviation:
					data - min
		new_data = log ------------- + 1
		               std
	The new data replaces the old data.
	Image statistics are recalculated.
**/
void		Bimage::logarithm()
{
	if ( compoundtype == TComplex ) complex_to_intensities();

	if ( datatype < Float ) change_type(Float);
	
	double			norm = 1.0/std;
	
	for ( long i=0; i<datasize; i++ )
		set(i, log(((*this)[i] - min)*norm + 1));
		
	statistics();
}

/*void		Bimage::logarithm()
{
	if ( compoundtype == TComplex ) complex_to_intensities();

	if ( datatype < Float ) change_type(Float);
	
	double			fudge = (avg - min)/1000;
	double			norm = 1.0/(avg - min + fudge);
	
	for ( long i=0; i<datasize; i++ )
		set(i, log(((*this)[i] - min + fudge)*norm));
		
	statistics();
}*/

/**
@brief 	Calculates the exponential of the image data.

	The image is first converted to floating point.
	The new data replaces the old data.
	Image statistics are recalculated.
**/
void		Bimage::exponential()
{
	if ( compoundtype == TComplex ) complex_to_intensities();

	if ( datatype < Float ) change_type(Float);
	
	for ( long i=0; i<datasize; i++ )
		set(i, exp((*this)[i]));
		
	statistics();
}

/**
@brief 	Calculates and corrects for a linear gradient across an image.
@return int			0.

	The linear gradient across an image is calculated as:
		density(x,y,z) = b0 + b1*x + b2*y + b3*z
	The image is converted to floating point and corrected for the gradient.

**/
int 		Bimage::gradient_correction()
{
	change_type(Float);
	
    long			i, j, nn, xx, yy, zz;
	double			v;
	vector<double>	b(4);
	Matrix			a(4,4), a2(4,4);
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Correcting the image linear gradient:" << endl;
    	cout << "Image\tb0\tb1\tb2\tb3" << endl;
	}
	
	// Precalculate the a matrix which depends only on the size
	// Keep it like this in case we want to use a mask later
	for ( zz=0; zz<z; zz++ ) {
		for ( yy=0; yy<y; yy++ ) {
			for ( xx=0; xx<x; xx++ ) {
				a[0][0] += 1;
				a[0][1] += xx;
				a[0][2] += yy;
				a[0][3] += zz;
				a[1][1] += xx*xx;
				a[1][2] += xx*yy;
				a[1][3] += xx*zz;
				a[2][2] += yy*yy;
				a[2][3] += yy*zz;
				a[3][3] += zz*zz;
			}
		}
	}
	for ( i=1; i<4; i++ ) {
		for ( j=0; j<i; j++ )
			a[i][j] = a[j][i];
		if ( a[i][i] < 1 ) a[i][i] = 1; // The diagonals cannot be zero
	}
	
	for ( nn=0; nn<n; nn++ ) {
		for ( i=0; i<4; i++ ) b[i] = 0;
		a2 = a;
		for ( i = nn*x*y*z, zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				for ( xx=0; xx<x; xx++, i++ ) {
					v = (*this)[i];
					b[0] += v;
					b[1] += v*xx;
					b[2] += v*yy;
					b[3] += v*zz;
				}
			}
		}
		a2.LU_decomposition(b);
		if ( verbose & VERB_PROCESS )
    		cout << n << tab << b[0] << tab << b[1] << tab << b[2] << tab << b[3] << endl;
		for ( i = nn*x*y*z, zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				for ( xx=0; xx<x; xx++, i++ ) {
					v = (*this)[i] - (b[0] + b[1]*xx + b[2]*yy + b[3]*zz);
					set(i, v);
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Corrects for a quadric surface.
@param 	param		7-value array of parameters.
@return int 		error code.

**/
int			Bimage::quadric_correct(vector<double> param)
{
	long			i, xx, yy, zz, nn;
	double			dx, dy, dz, v;
	
	param[0] -= average();
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
//			dz = -image[nn].origin()[2] + zz;
			dz = -image[nn].origin()[2] + zz;
			for ( yy=0; yy<y; yy++ ) {
//				dy = -image[nn].origin()[1] + yy;
				dy = -image[nn].origin()[1] + yy;
				for ( xx=0; xx<x; xx++, i++ ) {
//					dx = -image[nn].origin()[0] + xx;
					dx = -image[nn].origin()[0] + xx;
					v = param[0] + dx*param[1] + dy*param[2] + dz*param[3]
						+ dx*dx*param[4] + dy*dy*param[5] + dz*dz*param[6];
					add(i, -v);
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Fits the whole image to a quadric surface.
@return vector<double>	7-value array of parameters.

	A quadric surface is defined as:
	v = a0 + a1*dx + a2*dy + a3*dz + a4*dx^2 + a5*dy^2 + a6*dz^2
	
**/
vector<double>	Bimage::quadric_fit()
{
	long			i, j, nn, nterm(7);
	long			xx, yy, zz;
	double			dx, dy, dz;
	Matrix			a(nterm,nterm);
	vector<double>	b(nterm);
	vector<double>	vec(nterm);

	for ( nn=0; nn<n; nn++ ) {
		for ( i=0; i<nterm; i++ ) {	// Clear the matrix and vector
			for ( j=0; j<nterm; j++ ) a[i][j] = 0;
			b[i] = 0;
		}
		vec[0] = 1;
		for ( zz=0; zz<z; zz++ ) {
			dz = -image[nn].origin()[2] + zz;
			vec[3] = dz;
			vec[6] = dz*dz;
			for ( yy=0; yy<y; yy++ ) {
				dy = -image[nn].origin()[1] + yy;
				vec[2] = dy;
				vec[5] = dy*dy;
				for ( xx=0; xx<x; xx++ ) {
					dx = -image[nn].origin()[0] + xx;
					vec[1] = dx;
					vec[4] = dx*dx;
					i = index(0,xx,yy,zz,nn);
					for ( j=0; j<nterm; j++ ) b[j] += vec[j]*(*this)[i];
					for ( i=0; i<nterm; i++ )
						for ( j=0; j<=i; j++ ) a[i][j] += vec[i]*vec[j];
				}
			}
		}
		for ( i=0; i<nterm-1; i++ )
			for ( j=i+1; j<nterm; j++ ) a[i][j] = a[j][i];
		for ( i=0; i<nterm; i++ )
			if ( fabs(a[i][i]) < 1e-37 ) a[i][i] = 1;
		a.LU_decomposition(b);
		if ( verbose & VERB_PROCESS )
			cout << "Quadric coefficients: " <<
				b[0] << tab << b[1] << tab << b[2] << tab << b[3] << tab << b[4] << tab << b[5] << tab << b[6] << endl << endl;
	}

	return b;
}

/**
@brief 	Calculates a thickness based on intensities with respect to a reference.
@param	reference		reference intensity.
@param	emfp			proportionaility coefficient.
@return vector<double>	7-value array of parameters.

	The thickness for each pixel is:
		t = emfp * ln(Iref/Ipix)
	If Ipix <= 0, t = 0
	
**/
Bimage*		Bimage::thickness(double reference, double emfp)
{
	if ( reference < 1e-30 ) reference = maximum();
	
	long		i, ds(x*y*z*n);
	double		v;
	Bimage*		pthick = new Bimage(Float, TSimple, size(), n);
	pthick->sampling(sampling(0));
	
	if ( verbose ) {
		cout << "Calculating a thickness image:" << endl;
		cout << "Proportionality coefficient:    " << emfp << endl;
		cout << "Reference intensity:            " << reference << endl << endl;
	}
	
	for ( i=0; i<ds; ++i ) {
		v = (*this)[i];
		if ( v > 0 && v < reference )
			pthick->set(i, emfp*log(reference/v));
		else
			pthick->set(i, 0);
	}
	
	return pthick;
}

