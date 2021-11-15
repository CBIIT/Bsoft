/**
@file	Bimage_fspace.cpp
@brief	Library routines used for modifying reciprocal space amplitudes
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210727
**/

#include "Bimage.h"
#include "ctf.h"
#include "matrix_linear.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the maximum radius in frequency space from a given resolution.
@param	resolution		high resolution limit.
@param	sampling_ratio	frequency space sampling (default 1 pixel/sample).
@return long			maximum radius.
	
	The maximum radius is either the physical image size divided by the resolution,
	or Nyquist frequency.
**/
long		Bimage::fspace_maximum_radius(double resolution, double sampling_ratio)
{
	double			rad_scale(real_size().max()/sampling_ratio);

	long			maxrad(size().max());
	
	if ( resolution > 0.1 )
		maxrad = (long) (rad_scale/resolution);
	else
		maxrad = (long) (maxrad/2.0);

	return maxrad;
}


/**
@brief 	Calculates the background for a Fourier transform.
@return int			0.
	The background is taken to be the area outside Nyquest.
**/
int			Bimage::fspace_background()
{
	long			nn, xx, yy, zz, i, nb;
	double			bkg;
	Vector3<long>	bmin(0.4*x, 0.4*y, 0.4*z), bmax(0.6*x, 0.6*y, 0.6*z);
	bmax = bmax.min(1);
	
	if ( verbose )
		cout << "Calculating the Fourier transform background" << endl;
	
	for ( nn=0; nn<n; ++nn ) {
		nb = 0;
		bkg = 0;
		for ( zz=bmin[2]; zz<bmax[2]; ++zz ) {
			for ( yy=bmin[1]; yy<bmax[1]; ++yy ) {
				for ( xx=bmin[0]; xx<bmax[0]; ++xx ) {
					i = index(xx, yy, zz, nn);
					bkg += complex(i).power();
					nb++;
				}
			}
		}
		if ( nb ) bkg /= nb;
		image[nn].background(bkg);
		if ( verbose )
			cout << nn << tab << scientific << bkg << endl;
	}
	
	return 0;
}

/**
@brief 	Calculates the complex value at an image location by kernel-based interpolation.
@param 	img_num				sub-image number.
@param 	m					location in image.
@param	kernel				interpolation kernel.
@return Complex<double>		interpolated value.
	The kernel lookup table must be precalculated.
**/
Complex<double>	Bimage::fspace_interpolate(long img_num, Vector3<double> m, FSI_Kernel* kernel)
{
	Complex<double>	value(0,0);

	if ( compoundtype != TComplex ) return value;
	
	long			xx, yy, zz, kx, ky, kz, mx, my, mz, jx, jy, jz;
	long			nx(kernel->width()), ny(kernel->width()), nz(kernel->width());
	long			ix = (long) floor(m[0]);
	long			iy = (long) floor(m[1]);
	long			iz = (long) floor(m[2]);
	long			ikx(kernel->half_width());
	long			iky(kernel->half_width());
	long			ikz(kernel->half_width());
	double			wyz;

	if ( x > 1 ) {
		ikx = ((long) ((m[0] - ix)*kernel->divisions() + 0.5))*kernel->width();
		ix -= kernel->half_width();
	} else nx = 1;
	
	if ( y > 1 ) {
		iky = ((long) ((m[1] - iy)*kernel->divisions() + 0.5))*kernel->width();
		iy -= kernel->half_width();
	} else ny = 1;
	
	if ( z > 1 ) {
		ikz = ((long) ((m[2] - iz)*kernel->divisions() + 0.5))*kernel->width();
		iz -= kernel->half_width();
	} else nz = 1;
	
	for ( zz=0, kz=ikz, mz=iz; zz<nz; ++zz, kz++, mz++ ) {
		jz = mz;
		while ( jz < 0 ) jz += z;
		while ( jz >= z ) jz -= z;
		jz = (img_num*z + jz)*y;
		for ( yy=0, ky=iky, my=iy; yy<ny; ++yy, ky++, my++ ) {
			jy = my;
			while ( jy < 0 ) jy += y;
			while ( jy >= y ) jy -= y;
			jy = (jz + jy)*x;
			wyz = (*kernel)[ky]*(*kernel)[kz];
			for ( xx=0, kx=ikx, mx=ix; xx<nx; ++xx, kx++, mx++ ) {
				jx = mx;
				while ( jx < 0 ) jx += x;
				while ( jx >= x ) jx -= x;
				value += complex(jy + jx) * ((*kernel)[kx]*wyz);
			}
		}
	}
	
	if ( !value.is_finite() )
		cerr << "Warning: Interpolated complex value not finite! (" << value << ")" << endl;

    return value;
}

/**
@brief 	Translates an image in frequency space to avoid interpolation.
@param 	shift		3-value real space shift vector.
@return int			error code.
**/
int			Bimage::fspace_translate(Vector3<double> shift)
{
	if ( verbose & VERB_FULL ) {
		cout << "Translating in frequency space:" << endl;
		cout << "Shift:                          " << shift << endl;
	}
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();

	for ( long nn=0; nn<n; nn++ )
		fspace_translate(nn, shift);
	
	if ( old_transform == NoTransform )
		fft_back();
	
	return 0;
}


/**
@brief 	Translates an image in frequency space to avoid interpolation.
@param	nn			sub-image to transform.
@param 	shift		3-value real space shift vector.
@return int			0.
**/
int			Bimage::fspace_translate(long nn, Vector3<double> shift)
{
	if ( verbose & VERB_FULL ) {
		cout << "Translating in frequency space:" << endl;
		cout << "Shift:                          " << shift << endl;
	}
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();

	phase_shift(nn, shift);
	
	if ( old_transform == NoTransform )
		fft_back();
	
	return 0;
}

/**
@brief 	Resizes an image in frequency space to avoid interpolation.
@param 	scale			isotropic scaling.
@param 	res_hi			high resolution limit.
@param	res_lo			low resolution limit.
@return int				0.
**/
int			Bimage::fspace_resize(double scale, double res_hi, double res_lo)
{
	Vector3<long>		nusize = {x,y,z};
	if ( x > 1 ) nusize[0] = (long) (x*scale);
	if ( y > 1 ) nusize[1] = (long) (y*scale);
	if ( z > 1 ) nusize[2] = (long) (z*scale);
	
	Vector3<long>		translate = {0,0,0};
	Vector3<long>		olsize = {x,y,z};
	if ( x > 1 ) olsize[0] = (long) (nusize[0]/scale);
	if ( y > 1 ) olsize[1] = (long) (nusize[1]/scale);
	if ( z > 1 ) olsize[2] = (long) (nusize[2]/scale);
	
	Vector3<double>		error;
	if ( x > 1 ) error[0] = olsize[0]*scale - nusize[0];
	if ( y > 1 ) error[1] = olsize[1]*scale - nusize[1];
	if ( z > 1 ) error[2] = olsize[2]*scale - nusize[2];
	
	if ( verbose ) {
//	if ( verbose & VERB_PROCESS ) {
		cout << "Rescaling in frequency space:" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "New size:                       " << nusize << " voxels" << endl;
		cout << "Error:                          " << error << " voxels" << endl << endl;
	}
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform ) {
		calculate_background();
		resize(olsize, translate, FILL_BACKGROUND, 0);
		fft();
	}

	change_transform_size(nusize);
	
	if ( res_hi > 0 ) fspace_bandpass(res_hi, res_lo, 0);
	
	if ( old_transform == NoTransform )
		fft_back();
	
	return 0;
}

/**
@brief 	Resizes an image in frequency space to avoid interpolation.
@param 	pref			reference.
@return Bimage*			resized image.
	The image is first resized to approximate the real size of
	the reference map.
**/
Bimage*		Bimage::fspace_resize(Bimage* pref)
{
	if ( image->origin().length() < 1 ) origin(size()/2);
	
	double			r(pref->real_size()[0]/real_size()[0]);
	long			nx(x * r + 0.5);
//	r = nx/(x * r);

	Vector3<long>	nusize(nx, nx, nx);
	Vector3<double> scale(1,1,1), origin(image->origin());
//	Vector3<double> translate(nusize/2 - image->origin());
//	Vector3<double> translate(image->origin() * (nx * r/x - 1));
//	Vector3<double> translate((nusize - size())/2);
	Vector3<double> translate(nusize/2 - size()/2);
	Matrix3 		mat(1);
	
	if ( verbose & VERB_RESULT ) {
		cout << "Fitting " << file_name() << " to " << pref->file_name() << endl;
		cout << "Intermediate size:              " << nusize << endl;
		cout << "Intermediate scale:             " << scale << endl;
		cout << "Intermediate translation:       " << translate << endl << endl;
	}
		
	if ( verbose & VERB_PROCESS ) {
		cout << "Comparing map " << file_name() << endl;
		cout << "Map size:                       " << size() << endl;
		cout << "Map origin:                     " << image->origin() << endl;
		cout << "Map sampling:                   " << sampling(0) << endl;
		cout << "Map real size:                  " << real_size() << endl << endl;
		cout << "Reference:                      " << pref->file_name() << endl;
		cout << "Reference size:                 " << pref->size() << endl;
		cout << "Reference origin:               " << pref->image->origin() << endl;
		cout << "Reference sampling:             " << pref->sampling(0) << endl;
		cout << "Reference real size:            " << pref->real_size() << endl << endl;
	}
	
	Bimage*			pt = transform(nusize, scale, origin, translate, mat, FILL_BACKGROUND);

//	write_img("pt.pif", pt);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Translate:                      " << translate << endl;
		cout << "New map size:                   " << pt->size() << endl;
		cout << "New map origin:                 " << pt->image->origin() << endl;
		cout << "New map sampling:               " << pt->sampling(0) << endl;
		cout << "New map real size:              " << pt->real_size() << endl << endl;
	}
	
	pt->fft();

	pt->fspace_bandpass(2*pt->image->sampling()[0]);	// Cut out the corners
	
	Vector3<long>	sz = {pref->sizeX(), pref->sizeY(), pref->sizeZ()};
	pt->change_transform_size(sz);
	
	pt->fft_back();

	return pt;
}

/**
@brief 	Sets all amplitudes to one.
@return int 		0.

	If the image is a not a Fourier transform, it is transformed, modified,
	and backtransformed. If the image is a Fourier transform, it is just
	modified. The resultant image is floating point for real space or
	complex for reciprocal space.

**/
int			Bimage::fspace_amp_one()
{
	if ( !d.uc ) {
		cerr << "Error: No image data to process!" << endl << endl;
		return -1;
	}
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();
	
	if ( verbose & VERB_PROCESS )
		cout << "Setting all amplitudes to one" << endl << endl;

    long			i, ds(x*y*z*n);
    double			amp;
	Complex<double>	cv;

    for ( i=0; i<ds; i++ ) {
		cv = complex(i);
		if ( ( amp = cv.amp() ) )
			set(i, cv/amp);
	}
		
	if ( old_transform == NoTransform )
		fft_back();
	
	statistics();
	
	return 0;
}

/**
@brief 	Filters the amplitudes of the Fourier transform of an image.
@param 	threshold	Miminum amplitude to accept.
@return int 			0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered. The filtering sets all amplitudes below the given threshold
	to zero. The resultant image is floating point for real space or
	complex for reciprocal space.

**/
int			Bimage::fspace_amp_threshold(double threshold)
{
	if ( !d.uc ) {
		cerr << "Error: No image data to process!" << endl << endl;
		return -1;
	}

	if ( threshold <= 0 ) return 0;
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();
	
	if ( verbose & VERB_PROCESS )
		cout << "Applying an amplitude filter with a threshold of " << threshold << endl;

    long			i, ds(x*y*z*n);
    long			nk(0), nz(0);
	double			maxamp(0), zeroamp(1e-5*threshold);
	double			amp, ampkept(0), ampzeroed(0);
	Complex<double>	cv;

    for ( i=0; i<ds; i++ ) {
		amp = complex(i).amp();
		if ( maxamp < amp ) maxamp = amp;
		if ( amp < threshold ) {
			ampzeroed += amp;
			if ( amp > zeroamp ) nz++;
			set(i, cv);
		} else {
			ampkept += amp;
			nk++;
		}
	}
		
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of amplitudes kept:      " << nk << " (" << nk*100.0/ds << " %)" << endl;
		cout << "Number of amplitudes zeroed:    " << nz << " (" << nz*100.0/ds << " %)" << endl;
		cout << "Maximum amplitude:              " << maxamp << endl;
		cout << "Amplitude kept sum:             " << ampkept << endl;
		cout << "Amplitude zeroed sum:           " << ampzeroed << endl << endl;
	}
		
	if ( old_transform == NoTransform )
		fft_back();
	
	statistics();
	
	return 0;
}

/**
@brief 	Change the amlitudes to their square roots.
@return int			0.
**/
int			Bimage::fspace_sqrt_amp()
{
	if ( fouriertype != NoTransform ) {
		cerr << "Error: File " << file_name() << " must be a real space map!" << endl;
		return -1;
	}
	
	fft();

	if ( verbose & VERB_PROCESS )
		cout << "Changing the amplitudes to their square roots" << endl << endl;
	
	long			i, ds(x*y*z*n);
	double			d;
	Complex<float>	cv;

	for ( i=0; i<ds; i++ ) {
		cv = complex(i);
		d = sqrt(cv.amp());
		set(i, cv/d);
	}
    
	fft_back();
	
	return 0;
}

/**
@brief 	Change the amlitudes to their squares.
@return int			0.
**/
int			Bimage::fspace_square_amp()
{
	if ( fouriertype != NoTransform ) {
		cerr << "Error: File " << file_name() << " must be a real space map!" << endl;
		return -1;
	}
	
	fft();

	if ( verbose & VERB_PROCESS )
		cout << "Changing the amplitudes to their squares" << endl << endl;
	
	long			i, ds(x*y*z*n);
	double			d;
	Complex<float>	cv;

	for ( i=0; i<ds; i++ ) {
		cv = complex(i);
		d = cv.amp();
		set(i, cv*d);
	}
    
	fft_back();
	
	return 0;
}

/**
@brief 	Applies a bandpass filter to an image.
@param 	res_hi		high resolution limit.
@param 	res_lo		low resolution limit.
@param 	width		gaussian width of edge.
@return int 		0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered. The filtering sets all values with frequencies above the
	given high resolution limit and below the given low resolution limit
	to zero.

**/
int			Bimage::fspace_bandpass(double res_hi, double res_lo, double width)
{
	fft_plan		planf = fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = fft_setup(FFTW_BACKWARD, 0);

	fspace_bandpass(res_hi, res_lo, width, planf, planb);

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	return 0;
}

/**
@brief 	Applies a bandpass filter to an image.
@param 	res_hi		high resolution limit.
@param 	res_lo		low resolution limit.
@param 	width		gaussian width of edge.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return int 		0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered. The filtering sets all values with frequencies above the
	given high resolution limit and below the given low resolution limit
	to zero.

**/
int			Bimage::fspace_bandpass(double res_hi, double res_lo, double width, fft_plan planf, fft_plan planb)
{
	if ( !d.uc ) {
		cerr << "Error: No image data to process!" << endl << endl;
		return -1;
	}
	
	if ( width < 1e-6 ) width = 1e-6;
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft(planf, 1);
	
	if ( res_lo > 0 && res_hi > res_lo ) swap(res_hi, res_lo);
	if ( res_hi == res_lo ) res_lo = 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a bandpass filter:" << endl;
		cout << "High and low resolution limits: " << res_hi << " - ";
		if ( res_lo > 0 ) cout << res_lo << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Gaussian edge:                  " << width << endl << endl;
	}

	long			i, nn, xx, yy, zz;
    double   		edge_lo, edge_hi, f, test_width(10*width), a(GOLDEN/fabs(width));
	double			s, sx2, sy2, sz2;
	double			shi = 1/res_hi;
	double			slo = (res_lo > 0)? 1/res_lo : 0;
	Vector3<long>	h((size()-1)/2);
	Vector3<double> iscale(1/real_size());
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			sz2 = zz;
			if ( zz > h[2] ) sz2 -= z;
			sz2 *= iscale[2];
			sz2 *= sz2;
			for ( yy=0; yy<y; ++yy ) {
				sy2 = yy;
				if ( yy > h[1] ) sy2 -= y;
				sy2 *= iscale[1];
				sy2 *= sy2;
				for ( xx=0; xx<x; ++xx, i++ ) {
					sx2 = xx;
					if ( xx > h[0] ) sx2 -= x;
					sx2 *= iscale[0];
					sx2 *= sx2;
					s = sqrt(sx2 + sy2 + sz2);
					if ( slo - s > test_width ) edge_lo = 1e300;
					else edge_lo = exp(a*(slo - s));
					if ( s - shi > test_width ) edge_hi = 1e300;
					else edge_hi = exp(a*(s - shi));
					f = (edge_hi > edge_lo)? 1/(1 + edge_hi): 1/(1 + edge_lo);
					set(i, complex(i) * f);
				}
			}
		}
	}
	
	if ( old_transform == NoTransform ) {
		fft(planb, 1);
		complex_to_real();
		fourier_type(NoTransform);
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Applies a frequency filter to an image.
@param 	freq		frequency or inverse frequency.
@param 	sigma		gaussian envelope width around frequency.
@return int 		0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered.
	The filter imposes a guassian envelope at the given frequency.
	If the frequency value is greater than one, it is assumed to be given
	as the inverse.

**/
int			Bimage::fspace_frequency_filter(double freq, double sigma)
{
	fft_plan		planf = fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = fft_setup(FFTW_BACKWARD, 0);

	fspace_frequency_filter(freq, sigma, planf, planb);

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	return 0;
}

/**
@brief 	Applies a frequency filter to an image.
@param 	freq		frequency or inverse frequency.
@param 	sigma		gaussian envelope width around frequency.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return int 			0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered.
	The filter imposes a guassian envelope at the given frequency.
	If the frequency value is greater than one, it is assumed to be given
	as the inverse.

**/
int			Bimage::fspace_frequency_filter(double freq, double sigma, fft_plan planf, fft_plan planb)
{
	if ( !d.uc ) {
		cerr << "Error: No image data to process!" << endl << endl;
		return -1;
	}
	
	if ( freq < 1e-6 ) return -1;
	if ( freq > 1 ) freq = 1/freq;	// Assume real space dimensions
	if ( sigma < 2*image->sampling()[0] ) sigma = 0.25/freq;
	if ( sigma < 2*image->sampling()[0] ) sigma = 2*image->sampling()[0];
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a frequency filter:" << endl;
		cout << "Frequency:                      " << freq << " 1/Å (" << 1/freq << " Å)" << endl;
		cout << "Gaussian sigma:                 " << sigma << " Å" << endl << endl;
	}

	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft(planf, 1);
	
	long			i, nn, xx, yy, zz;
	double			d;
	Vector3<double>	s;
	Vector3<double> fscale(1.0L/x, 1.0L/y, 1.0L/z);

	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			s[2] = fscale[2]*zz;
			if ( s[2] >= 0.5 ) s[2] -= 1;
			s[2] /= image->sampling()[2];
			for ( yy=0; yy<y; ++yy ) {
				s[1] = fscale[1]*yy;
				if ( s[1] >= 0.5 ) s[1] -= 1;
				s[1] /= image->sampling()[1];
				for ( xx=0; xx<x; ++xx, ++i ) {
					s[0] = fscale[0]*xx;
					if ( s[0] >= 0.5 ) s[0] -= 1;
					s[0] /= image->sampling()[0];
					d = s.length() - freq;
					d *= M_PI*sigma;
					d = exp(-2*d*d);
					set(i, complex(i) * d);
				}
			}
		}
	}

	if ( old_transform == NoTransform ) {
		fft(planb, 1);
		complex_to_real();
		fourier_type(NoTransform);
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Applies a Gabor filter to an image.
@param 	freq		frequency or inverse frequency location.
@param 	fsigma		gaussian envelope width in the direction of the frequency vector.
@param 	psigma		gaussian envelope width perpendicular to the frequency vector.
@return int 			0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered.
	The filter imposes a guassian envelope at the given frequency location.
	If the frequency vector size is greater than one, it is assumed to be given
	as the inverse.

**/
int			Bimage::fspace_gabor_filter(Vector3<double> freq, double fsigma, double psigma)
{
	fft_plan		planf = fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = fft_setup(FFTW_BACKWARD, 0);

	fspace_gabor_filter(freq, fsigma, psigma, planf, planb);

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	return 0;
}

/**
@brief 	Applies a Gabor filter to an image.
@param 	freq		frequency or inverse frequency direction location.
@param 	fsigma		gaussian envelope width in the direction of the frequency vector.
@param 	psigma		gaussian envelope width perpendicular to the frequency vector.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return int 			0.

	If the image is a not a Fourier transform, it is transformed, filtered,
	and backtransformed. If the image is a Fourier transform, it is just
	filtered.
	The filter imposes a guassian envelope at the given frequency location.
	If the frequency vector size is greater than one, it is assumed to be given
	as the inverse.

**/
int			Bimage::fspace_gabor_filter(Vector3<double> freq, double fsigma, double psigma, fft_plan planf, fft_plan planb)
{
	if ( !d.uc ) {
		cerr << "Error: No image data to process!" << endl << endl;
		return -1;
	}
	
	if ( freq.length() < 1e-6 ) return -1;
	if ( freq.length() > 1 ) {		// Assume real space dimensions
		if ( freq[0] ) freq[0] = 1/freq[0];
		if ( freq[1] ) freq[1] = 1/freq[1];
		if ( freq[2] ) freq[2] = 1/freq[2];
	}
	if ( fabs(freq[0]) < 1e-6 ) freq[0] = 1e-6;
	if ( fabs(freq[1]) < 1e-6 ) freq[1] = 1e-6;
	if ( fabs(freq[2]) < 1e-6 ) freq[2] = 1e-6;
	
	if ( fsigma < 2*image->sampling()[0] ) fsigma = 0.25/freq.length();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a frequency filter:" << endl;
		cout << "Frequency:                      " << freq << " 1/Å (" << 1/freq << " Å)" << endl;
		cout << "Sigma in frequency direction:   " << fsigma << " Å" << endl;
		cout << "Sigma perpendiular:             " << psigma << " Å" << endl << endl;
	}

	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft(planf, 1);
	
	long			i, nn, xx, yy, zz;
	double			d, w;
	Vector3<double>	s, v, v0;
	Vector3<double> fscale(1.0L/x, 1.0L/y, 1.0L/z);

	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			s[2] = fscale[2]*zz;
			if ( s[2] >= 0.5 ) s[2] -= 1;
			s[2] /= image->sampling()[2];
			for ( yy=0; yy<y; ++yy ) {
				s[1] = fscale[1]*yy;
				if ( s[1] >= 0.5 ) s[1] -= 1;
				s[1] /= image->sampling()[1];
				for ( xx=0; xx<x; ++xx, ++i ) {
					s[0] = fscale[0]*xx;
					if ( s[0] >= 0.5 ) s[0] -= 1;
					s[0] /= image->sampling()[0];
					v = s - freq;
					v *= M_PI*fsigma;
					if ( psigma ) {
						w = M_PI*psigma * s.distance_from_line(freq, v0);
						d = exp(-2*(v.scalar(v) + w*w));
					} else {
						d = exp(-2*v.scalar(v));
					}
					set(i, complex(i) * d);
				}
			}
		}
	}

	if ( old_transform == NoTransform ) {
		fft(planb, 1);
		complex_to_real();
		fourier_type(NoTransform);
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Calculates the radial power spectrum from a Fourier transform.
@param 	resolution		high resolution limit.
@param	sampling_ratio	frequency space sampling (default 1 pixel/sample).
@return Bimage*			radial power spectrum in the form of a 1D image.

	A radial average of a 2D or 3D Fourier transform is calculated.  
	An interpolative method is used where the value of 
	a voxel is distributed between the two nearest radial annuli.
	The final sum in an annulus is normalized by the number of voxels 
	contributing to the annulus sum.

**/
Bimage* 	Bimage::fspace_radial_power(double resolution, double sampling_ratio)
{
	int				cmplx(0);
	
//	Bstring			ct = compound_type_string();
//	cout << ct << endl;

	if ( compoundtype == TComplex && fouriertype == Standard ) {
		cmplx = 1;
	} else if ( compoundtype != TSimple ) {
		cerr << "Error: File " << file_name() << " must be a Fourier transform or power spectrum!" << endl;
		return NULL;
	}
	
	check_resolution(resolution);
	
	long			i, j, nn, xx, yy, zz, iradius;
	double			radius, f, f1, rx, ry, rz, v;
	double			rad_scale(real_size()[0]/sampling_ratio);
	Vector3<long>	h((size()-1)/2);
	Vector3<double>	freq_scale(1/real_size());
	long			maxrad = fspace_maximum_radius(resolution, sampling_ratio);
	
	Bimage*			prad = new Bimage(Float, TSimple, maxrad, 1, 1, n);
	prad->sampling(1.0/rad_scale, 1, 1);
	
	vector<float>	radstd(n*maxrad);
	vector<float>	num(maxrad);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a radial power spectrum:" << endl;
		cout << "Radial scale:                   " << rad_scale << " A" << endl;
		cout << "Radial step size:               " << 1.0/rad_scale << " 1/A" << endl;
		cout << "Resolution limit:               " << resolution << " A" << endl;
		cout << "Number of values:               " << maxrad << endl << endl;
	}
	
	for ( nn=0; nn<n; nn++ ) {
		for ( i=0; i<maxrad; i++ ) num[i] = 0;
		for ( zz=0; zz<z; ++zz ) {
			rz = zz;
			if ( zz > h[2] ) rz -= z;
			rz *= freq_scale[2];
			rz *= rz;
			for ( yy=0; yy<y; ++yy ) {
				ry = yy;
				if ( yy > h[1] ) ry -= y;
				ry *= freq_scale[1];
				ry *= ry;
				for ( xx=0; xx<x; ++xx ) {
					rx = xx;
					if ( xx > h[0] ) rx -= x;
					rx *= freq_scale[0];
					rx *= rx;
					radius = rad_scale*sqrt(rx + ry + rz);
					iradius = (long) radius;
					if ( iradius < maxrad ) {
						f = radius - iradius;
						f1 = 1 - f;
						i = index(xx, yy, zz, nn);
						num[iradius] += f1;
						j = nn*maxrad + iradius;
						if ( cmplx ) v = complex(i).power();
						else v = (*this)[i];
						prad->add(j, f1*v);
						radstd[j] += f1*v*v;
						iradius++;
						if ( iradius < maxrad ) {
							j++;
							num[iradius] += f;
							prad->add(j, f*v);
							radstd[j] += f*v*v;
						}
					}
				}
			}
		}
		for ( i=0; i<maxrad; i++ ) {
			if ( num[i] ) {
				j = nn*maxrad + i;
				prad->set(j, (*prad)[j]/num[i]);
				radstd[j] = radstd[j]/num[i] - (*prad)[j]*(*prad)[j];
				if ( radstd[j] > 0 )
					radstd[j] = sqrt(radstd[j]);
				else
					radstd[j] = 0;
			}
		}
		prad->image[nn].origin(0.0,0.0,0.0);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fspace_radial_power: Calculating statistics" << endl;
	
	prad->statistics();

	if ( verbose & VERB_FULL ) {
		cout << endl << "Radial power spectrum:" << endl;
		cout << "Radius";
		for ( nn=0; nn<prad->images(); nn++ ) cout << " " << nn+1 << " (std)";
		cout << endl;
		for ( xx=0; xx<prad->sizeX(); ++xx ) {
			cout << xx*prad->sampling(0)[0];
			for ( nn=0, i=xx; nn<prad->images(); nn++, i+=prad->sizeX() )
				cout << " " << (*prad)[i] << " " << radstd[i];
			cout << endl;
		}
	}
	
	return prad;
}

/**
@brief  Calculates the radial power spectrum from a Fourier transform.
@param  nn    		sub-image number.
@param  maxrad	maximum radius (i.e., high resolution limit).
@param  flag      	flag to calculate power (0) or amplitude (1), and not normalize (2).
@return double*     	radial power spectrum in the form of a 1D image.
 
A radial average of a 2D or 3D Fourier transform is calculated.
 An interpolative method is used where the value of
 a voxel is distributed between the two nearest radial annuli.
 The final sum in an annulus is normalized by the number of voxels
 contributing to the annulus sum.
 
**/
vector<double>	Bimage::fspace_radial(long nn, long maxrad, int flag)
{
	bool			amp(flag & 1), norm(1 - (flag & 2));
	long			i, xx, yy, zz, iradius, iradius2;
	double			radius, fraction, fraction2, v;
	double			rx, ry, rz;
	Vector3<long>	h((size()-1)/2);
	double			rad_scale(real_size()[0]);
	Vector3<double>	freq_scale(1.0/real_size());
	
	vector<double>	num(maxrad,0);
	vector<double>	rps(maxrad,0);
	
	for ( zz=0; zz<z; ++zz ) {
		rz = zz;
		if ( zz > h[2] ) rz -= z;
		rz *= freq_scale[2];
		rz *= rz;
		for ( yy=0; yy<y; ++yy ) {
			ry = yy;
			if ( yy > h[1] ) ry -= y;
			ry *= freq_scale[1];
			ry *= ry;
			for ( xx=0; xx<x/2; ++xx ) {
				rx = xx;
				if ( xx > h[0] ) rx -= x;
				rx *= freq_scale[0];
				rx *= rx;
				i = index(xx, yy, zz, nn);
				radius = rad_scale*sqrt(rx + ry + rz);
				iradius = (long) radius;
				iradius2 = iradius + 1;
				if ( iradius2 < maxrad ) {
					fraction = radius - iradius;
					fraction2 = 1.0 - fraction;
					if ( compoundtype == TComplex ) {
						if ( amp )
							v = complex(i).amp();
						else
							v = complex(i).power();
					} else {
						v = (*this)[i];
					}
					num[iradius] += fraction2;
					num[iradius2] += fraction;
					rps[iradius] += fraction2*v;
					rps[iradius2] += fraction*v;
				}
			}
		}
	}

	if ( norm ) {
		rps[0] = 0;
		for ( i=0; i<maxrad; i++ )
			if ( num[i] ) rps[i] /= num[i];
	}
	
	return rps;
}

int 		Bimage::fspace_weigh(vector<double>& curve, int flag)
{
	long			maxrad(curve.size());
	long			i, nn;
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();
	
	for ( nn=0; nn<n; nn++ ) {
		vector<double>	rps = fspace_radial(nn, maxrad, flag);
		for ( i=0; i<maxrad; i++ )
			if ( rps[i] ) rps[i] = curve[i]/rps[i];
		fspace_scale(nn, rps);
	}
	
	if ( old_transform == NoTransform )
		fft_back();
	
	return 0;
}

int 		Bimage::fspace_scale(vector<double>& scale, Bimage* pmask)
{
	long			nn;
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();
	
	for ( nn=0; nn<n; nn++ )
		fspace_scale(nn, scale, pmask);
	
	if ( old_transform == NoTransform )
		fft_back();
	
	return 0;
}

int 		Bimage::fspace_scale(long nn, vector<double>& scale, Bimage* pmask)
{
	long			maxrad(scale.size());
	int				use_this;
	long			i, j, xx, yy, zz, iradius, iradius2;
	double			rx, ry, rz;
	Vector3<long>	h((size()-1)/2);
	double			radius, fraction, fraction2, w;
	Complex<double>	cv;
	
	double			rad_scale = real_size()[0];
	Vector3<double>	freq_scale = 1.0/real_size();
	
	for ( j=zz=0; zz<z; ++zz ) {
		rz = zz;
		if ( zz > h[2] ) rz -= z;
		rz *= freq_scale[2];
		rz *= rz;
		for ( yy=0; yy<y; ++yy ) {
			ry = yy;
			if ( yy > h[1] ) ry -= y;
			ry *= freq_scale[1];
			ry *= ry;
			for ( xx=0; xx<x; ++xx, j++ ) {
				rx = xx;
				if ( xx > h[0] ) rx -= x;
				rx *= freq_scale[0];
				rx *= rx;
				i = index(xx, yy, zz, nn);
				radius = rad_scale*sqrt(rx + ry + rz);
				iradius = (long) radius;
				iradius2 = iradius + 1;
				use_this = 1;
				if ( pmask && (*pmask)[j] < 1 ) use_this = 0;
				if ( iradius2 >= maxrad ) use_this = 0;
				if ( use_this ) {
					fraction = radius - iradius;
					fraction2 = 1.0 - fraction;
					w = scale[iradius]*fraction2 + scale[iradius2]*fraction;
					set(i, complex(i) * w);
				} else {
					set(i, cv);
				} 
			} 
		}
	}
	
	return 0;
} 

/**
@brief 	Determines the overall B-factor of an image.
@param 	res_hi		high resolution limit.
@return double 		B-factor.

	The input image must a real space image. A radial power spectrum
	is calculated and fitted to the linearized version of the function:
		f^2 = scale * F^2 * exp(-B_factor/2 * s^2)
	where:
		f:	scattering profile for carbon
		F:	radial average amplitude
		s:	reciprocal space distance
		scale:	arbitrary scale
	The linear form of the function is:
		4*(log(f) - log(F)) = 2*log(scale) - B * s^2

**/
double		Bimage::fspace_fit_B_factor(double res_hi)
{
	if ( !d.uc ) {
		cerr << "Error: Cannot estimate the B-factor - the data block is empty!" << endl << endl;
		return -1;
	}
	
	check_resolution(res_hi);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "B-factor estimation:" << endl;
		cout << "High resolution limit:          " << res_hi << endl;
	}

	if ( compoundtype != TComplex ) fft();
	
	Bimage* 	prps = fspace_radial_power(res_hi);
	
	long 		i;
	double		realsize(real_size()[0]);
	double		frac_res2(1.0/(realsize));
	frac_res2 *= frac_res2;
	double		res_lo(200);
	long 		minrad((long) (realsize/res_lo));
	long 		maxrad((long) (realsize/res_hi));
	double		scat, logScale, Bfactor(-1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fspace_fit_B_factor: minrad=" << minrad << " maxrad=" << maxrad << endl;
	
	double*		s2 = new double[maxrad];
	double*		logF = new double[maxrad];
	
	for ( i=0; i<maxrad; i++ ) {
		s2[i] = i*i*frac_res2;
    	scat = 1.494*exp(-23.22*s2[i]) + 0.937*exp(-3.79*s2[i]);
		logF[i] = 4*(log(scat) - log((*prps)[i]));
		if ( verbose & VERB_FULL )
			cout << i << tab << sqrt(s2[i]) << tab << 1/sqrt(s2[i]) << tab << (*prps)[i] << tab << scat << endl;
    }
	
	double		CI = linear_least_squares(minrad, maxrad, s2, logF, logScale, Bfactor);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "B-factor:                       " << Bfactor << endl;
		cout << "Correlation index:              " << CI << endl << endl;
	} else if ( verbose & VERB_LABEL ) {
		cout << "B-factor:                       " << Bfactor << endl << endl;
	}
	
	delete[] s2;
	delete[] logF;
	delete prps;
	
	return Bfactor;
}

/**
@brief 	Weighs a transform with a ramp.
@param 	resolution	high resolution limit.
@param 	planf		2D forward Fourier transform plan.
@param 	planb		2D backward Fourier transform plan.
@return int 		0.

	Requirements:
		The data must be complex float and the FOM block must be allocated.

**/
int 		Bimage::fspace_weigh_ramp(double resolution, fft_plan planf, fft_plan planb)
{
	check_resolution(resolution);
	
	fft(planf, 1);
	
	long			i, xx, yy, zz, nzero(0), nuse(0);
	Vector3<long>	h(size()/2);
	Vector3<long>	k;
	Vector3<double>	s;
	double			s2max(1.0/(resolution*resolution)), weight(0), s2;
	Vector3<double>	realsize(real_size());
	
	if ( verbose & VERB_PROCESS )
		cout << "Filtering with a ramp to a resolution of " << resolution << " A" << endl;

	for ( i=zz=0; zz<z; ++zz ) {
		k[2] = (zz>h[2])? z - zz: zz;
		s[2] = k[2]/realsize[2];
		for ( yy=0; yy<y; ++yy ) {
			k[1] = (yy>h[1])? y - yy: yy;
			s[1] = k[1]/realsize[1];
			for ( xx=0; xx<x; ++xx, i++ ) {
				k[0] = (xx>h[1])? x - xx: xx;
				s[0] = k[0]/realsize[0];
				s2 = s.length2();
				if ( s2 > s2max ) {
					weight = 0;
					nzero++;
				} else {
					weight = sqrt(s2);
					nuse++;
				}
				set(i, complex(i) * weight);
			}
		}
	}

	if ( verbose & VERB_FULL )
		cout << "Weighed fractions: zeroed=" << nzero*1.0/i << " used=" << nuse*1.0/i << endl << endl;

	fft(planb, 1);
	
	complex_to_real();
	
	statistics();
	
	return 0;
}

/**
@brief 	Weighs a transform with a ramp.
@param 	resolution	high resolution limit.
@param 	axis		tilt axis for single tilt series (radians).
@param 	planf		2D forward Fourier transform plan.
@param 	planb		2D backward Fourier transform plan.
@return int 		0.

	Requirements:
		The data must be complex float and the FOM block must be allocated.

**/
int 		Bimage::fspace_weigh_ramp(double resolution, double axis, fft_plan planf, fft_plan planb)
{
	check_resolution(resolution);
	
	fft(planf, 1);
	
	long			i, xx, yy, zz, nzero(0), nuse(0);
	Vector3<long>	h(size()/2);
	Vector3<long>	k;
	Vector3<double>	s, a(cos(axis), sin(axis), 0);
	double			s2max(1.0/(resolution*resolution)), weight(0), s2;
	Vector3<double>	realsize(real_size());
	
	if ( verbose & VERB_PROCESS )
		cout << "Filtering with a ramp to a resolution of " << resolution << " A" << endl;

	for ( i=zz=0; zz<z; ++zz ) {
		k[2] = (zz>h[2])? z - zz: zz;
		for ( yy=0; yy<y; ++yy ) {
			k[1] = (yy>h[1])? y - yy: yy;
			s[1] = k[1]/realsize[1];
			for ( xx=0; xx<x; ++xx, i++ ) {
				k[0] = (xx>h[1])? x - xx: xx;
				s[0] = k[0]/realsize[0];
				s2 = s.length2();
				if ( s2 > s2max ) {
					weight = 0;
					nzero++;
				} else {
					weight = sqrt(s2)*sin(axis-atan2(s[1],s[0]));
					nuse++;
				}
				set(i, complex(i) * weight);
			}
		}
	}

	if ( verbose & VERB_FULL )
		cout << "Weighed fractions: zeroed=" << nzero*1.0/i << " used=" << nuse*1.0/i << endl << endl;

	fft(planb, 1);
	
	complex_to_real();
	
	statistics();
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with B-factor (gaussian) curve.
@param 	B			B factor.
@param 	resolution	high resolution limit.
@return int			0.

	The image is Fourier transformed and weighed with a gaussian curve:
		Fnew = F*exp(-B*s2/4)

**/
int 		Bimage::fspace_weigh_B_factor(double B, double resolution)
{
	check_resolution(resolution);
	
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	double			rad_scale = real_size()[0];
	
	//	maxrad is the main limit for the length of all the arrays, 
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + rad_scale/resolution);
	
	if ( verbose ) {
		cout << "Weighing amplitudes with a B factor:" << endl;
		cout << "B factor:                       " << B << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	long			i;
	vector<double>	gcurve(maxrad);
	double			s2, fac(-B/4);
	
	for ( i=0; i<maxrad; i++ ) {
		s2 = i/rad_scale;
		s2 *= s2;
		gcurve[i] = exp(fac*s2);
	}
	
	fspace_scale(gcurve);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with the carbon scattering curve.
@param 	resolution	high resolution limit.
@return int			0.

	The image is Fourier transformed and the radial power spectrum calculated.
	The ratio between the C curve and the average amplitudes in each shell
	is calculated and used to rescale the amplitudes of the image.

**/
int 		Bimage::fspace_weigh_C_curve(double resolution)
{
	check_resolution(resolution);
	
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	double			rad_scale = real_size()[0];
	
	//	maxrad is the main limit for the length of all the arrays, 
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + rad_scale/resolution);
	
	if ( verbose ) {
		cout << "Weighing amplitudes with a C cross-section:" << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	vector<double>	ccurve = C_curve(maxrad, 1/rad_scale);
	 
	fspace_weigh(ccurve, 1);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with a Laplacian-of-Gaussian function.
@param 	resolution	high resolution limit.
@param 	sigma		gaussian sigma.
@return int			0.

	The image is Fourier transformed and the radial power spectrum calculated.
	The ratio between the C curve and the average amplitudes in each shell
	is calculated and used to rescale the amplitudes of the image.

**/
int 		Bimage::fspace_weigh_LoG(double resolution, double sigma)
{
	check_resolution(resolution);
	
	//	maxrad is the main limit for the length of all the arrays,
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + real_size()[0]/resolution);
	
	if ( verbose ) {
		cout << "Weighing amplitudes with a Laplacian-of-Gaussian function:" << endl;
		cout << "Sigma:                          " << sigma << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	long			i;
	double			s2, freq_step(1/real_size()[0]), c1(-TWOPI*TWOPI), c2(-2*M_PI*M_PI*sigma*sigma);
	vector<double>	curve(maxrad);
	
	for ( i=0; i<maxrad; i++ ) {
		s2 = i*freq_step;
		s2 *= s2;
		curve[i] = c1*s2*exp(c2*s2);
	}
	
	fspace_weigh(curve, 1);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with a given RPS curve.
@param 	*plot		RPS plot.
@param 	resolution	high resolution limit.
@return int			0.

	The image is Fourier transformed and the radial power spectrum calculated.
	The ratio between the square root of the RPS curve and the average 
	amplitudes in each shell is calculated and used to rescale the amplitudes 
	of the image.

**/
int 		Bimage::fspace_weigh_RPS_curve(Bplot* plot, double resolution)
{
	check_resolution(resolution);

	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the x dimension.
	double			rad_scale = real_size()[0];
	double			scale = 1/((*plot)[1] * rad_scale);
	
	//	maxrad is the main limit for the length of all the arrays, 
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + rad_scale/resolution);
	
	if ( verbose ) {
		cout << "Weighing amplitudes with an RPS curve:" << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "rad_scale=" << rad_scale << endl;
		cout << "(*plot)[1]=" << (*plot)[1] << endl;
		cout << "scale=" << scale << endl;
	}
	
	long			i, j;
	double			s;
	vector<double>	rps(maxrad);

	for ( i=plot->rows(); i<2*plot->rows(); ++i ) (*plot)[i] = sqrt((*plot)[i]);
	
	for ( i=0; i<maxrad; ++i ) {
		s = scale*i;
		j = long(s);
		s -= j;
		j += plot->rows();
		rps[i] = (1-s)*(*plot)[j] + s*(*plot)[j+1];
		if ( verbose & VERB_DEBUG )
			cout << i << tab << j-plot->rows() << tab << i/rad_scale << tab << (j-plot->rows())*(*plot)[1] << tab << rps[i] << endl;
	}
	
	fspace_weigh(rps, 1);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with a given FSC curve.
@param 	*plot		FSC plot.
@param 	resolution	high resolution limit.
@return int			0.

	The image is Fourier transformed and the radial power spectrum calculated.
	The ratio between the square root of the FSC curve and the average 
	amplitudes in each shell is calculated and used to rescale the amplitudes 
	of the image.

**/
int 		Bimage::fspace_weigh_FSC_curve(Bplot* plot, double resolution)
{
	check_resolution(resolution);
	
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	double			rad_scale = real_size()[0];
	
	//	maxrad is the main limit for the length of all the arrays, 
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + rad_scale/resolution);
	if ( maxrad > plot->rows() ) maxrad = plot->rows();
	
	if ( verbose ) {
		cout << "Weighing amplitudes with an FSC curve:" << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	long			i, j;
	vector<double>	fsc(maxrad);
	
	for ( i=0, j=plot->rows(); i<maxrad; i++, j++ )
		fsc[i] = sqrt((*plot)[j]);
	 
	fspace_weigh(fsc, 1);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with an anisotropic Gaussian function.
@param 	nn			image index.
@param 	sigma		Gaussian sigma values.
@param	dir			derivative direction: 0=none, 1=x, 2=y, 3=z.
@return int			0.

	The image must be a complex Fourier transform.
	An anisotropic weight function is calculated at each voxel and applied to the transform.

**/
int 		Bimage::fspace_weigh_gaussian(long nn, Vector3<double> sigma, int dir)
{
	if ( dir < 0 || dir > 3 ) dir = 0;
	if ( z == 1 && dir > 2 ) dir = 0;
	
	long			i, xx, yy, zz, d1(dir-1);
	double			rx, ry, rz, w;
	Vector3<long>	h((size()-1)/2);
	Vector3<double>	s;
	Vector3<double>	fac = sigma * M_PI;
	Vector3<double>	freq_scale = 1.0/real_size();
	Complex<double>	cv;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Anisotropic Gaussian filter: " << endl;
		cout << "Widths:                         " << sigma << endl;
		if ( dir )
			cout << "Gradient direction:             " << dir << endl;
		cout << endl;
	}
	
	for ( i=nn*x*y*z, zz=0; zz<z; ++zz ) {
		s[2] = zz;
		if ( zz > h[2] ) s[2] -= z;
		s[2] *= freq_scale[2];
		rz = s[2]*fac[2];
		rz *= -2*rz;
		for ( yy=0; yy<y; ++yy ) {
			s[1] = yy;
			if ( yy > h[1] ) s[1] -= y;
			s[1] *= freq_scale[1];
			ry = s[1]*fac[1];
			ry *= -2*ry;
			for ( xx=0; xx<x; ++xx, ++i ) {
				s[0] = xx;
				if ( xx > h[0] ) s[0] -= x;
				s[0] *= freq_scale[0];
				rx = s[0]*fac[0];
				rx *= -2*rx;
				cv = complex(i) * exp(rx + ry + rz);
				if ( dir ) {
					w = TWOPI*s[d1];
					cv = Complex<double>(-cv.imag()*w, cv.real()*w);
				}
				set(i, cv);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates a image with orthogonal gradients encded in 3-value vectors.
@param 	sigma		Gaussian sigma values.
@return int			0.

	The image is Fourier transformed if needed.
	An anisotropic weight function is calculated at each voxel and applied to the transform.

**/
Bimage*		Bimage::fspace_gradient(Vector3<double> sigma)
{
	if ( fouriertype == NoTransform ) fft();

	long			i, j, k;
	long			nd = ( z > 1 )? 3: 2;
	Bimage*			pt;

	if ( verbose & VERB_PROCESS )
		cout << "Calculating a frequency space gradient" << endl << endl;

	Bimage*			pg = new Bimage(Float, TVector3, size(), 1);
	pg->sampling(sampling(0));
	pg->origin(image->origin());

	for ( i=0; i<nd; ++ i ) {
		pt = copy();
		pt->fspace_weigh_gaussian(0, sigma, i+1);
		pt->fft_back();
		for ( j=0, k=i; j<image_size(); ++j, k+=3 )
			pg->set(k, (*pt)[j]);
		delete pt;
	}

	return pg;
}

vector<double>	critical_exposure_curve(long size, double real_size, int flag)
{
	long			i;
	double			s, ds(1/real_size);
//	double			a(0.245), b(-1.665), c(2.81);	// Grant
	double			a(0.1), b(-2), c(3);			// Simplified
	double			deb(6), des(150), ec(30);		// Exponential decay

	if ( verbose & VERB_PROCESS ) {
		if ( flag )
			cout << "Critical exposure:              "
				<< deb << " + " << des << "exp(-" << ec << "s)"
				<< " e/A2" << endl << endl;
		else
			cout << "Critical exposure:              "
				<< a << "s^" << b << " + " << c
				<< " e/A2" << endl << endl;
	}
	
	vector<double>	critexp(size, 0);

	for ( i=0, s=0; i<size; ++i, s+=ds ) {
		if ( flag )
			critexp[i] = 2*(deb + des*exp(-ec*s));	// Exponential decay
		else
			critexp[i] = 2*(a*pow(s, b) + c);	// Grant
		if ( !isfinite(critexp[i]) ) critexp[i] = 1e10;
	}

	return critexp;
}

vector<double>	radiation_damage_curve(vector<double> dose, long size, double real_size)
{
	vector<double>	rdcurve(size, 0);

	if ( dose.size() < 1 ) return rdcurve;
//	if ( dose.size() < 2 ) dose.push_back(3);
	if ( dose.size() < 2 ) {
		dose.push_back(6);
		dose.push_back(150);
		dose.push_back(30);
	} else if ( dose.size() < 3 ) dose.push_back(0.1);

	double			s, is2, v, de(3), D(real_size);
	
	rdcurve[0] = 1;
	for ( long i=0; i<rdcurve.size(); ++i ) {
		s = i/D;
		if ( dose.size() == 3 ) {
			if ( i ) is2 = 1/s;					// Inverse frequency
			else is2 = 1;
			is2 *= is2;
			de = 2*(dose[1] + dose[2]*is2);		// Critical dose at frequency s (Grant)
		} else if ( dose.size() == 4 ) {
			de = 2*(dose[1] + dose[2]*exp(-dose[3]*s)); // Critical dose at frequency s
		}
		v = de*(1 - exp(-dose[0]/de));
		rdcurve[i] = v*v/dose[0];
	}

	return rdcurve;
}

int 		Bimage::fspace_weigh_dose(long nn, double dose_per_frame, vector<double> critdose)
{
	long 			maxrad(critdose.size());
	long			i;
	double			dose(0);
	
	// Accumulated dose in e/Å2
	if ( dose_per_frame > 1e-3 )
		for ( i=0; i<=nn; ++i ) dose += dose_per_frame;
	else
		for ( i=0; i<=nn; ++i ) dose += image[i].average()/image[i].sampling().volume();
	
	if ( dose < 1e-3 ) {
		cerr << "Error: The dose is too low! (" << dose << ")" << endl;
		return -1;
	}
		
	vector<double>	scale(maxrad);

	for ( i=0; i<maxrad; ++i )
		scale[i] = exp(-dose/critdose[i]);

	fspace_scale(nn, scale);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with the accumulated dose.
@param 	dose_per_frame	electron dose per frame.
@param	flag			0=Grant, 1=exponential decay.
@return int				0.

	The image is Fourier transformed and the amplitudes weighed using 
	the formula of Grant and Grigorieff (2015) or an exponential decay curve.

**/
int 		Bimage::fspace_weigh_dose(double dose_per_frame, int flag)
{
	long			maxrad = (long) (size().length()/2+2);

//	long			i;
//	double			s, ds(1.0/(real_size()[0]));
//	double			a(0.245), b(-1.665), c(2.81);	// Grant
//	double			a(0.1), b(-2), c(3);			// Simplified
//	double			deb(6), des(150), ec(30);		// Exponential decay

	if ( verbose ) {
		cout << "Weighing amplitudes with the accumulated dose:" << endl;
		cout << "Dose per frame:                 " << dose_per_frame << " e/Å2" << endl;
//		cout << "Critical dose:                  " << a << "*s^" << b << " + " << c << " e/A2" << endl << endl;
	}
/*
	vector<double>	critdose(maxrad);

	for ( i=0, s=0; i<maxrad; ++i, s+=ds ) {
		if ( flag )
			critdose[i] = 2*(deb + des*exp(-ec*s));	// Exponential decay
		else
			critdose[i] = 2*(a*pow(s, b) + c);	// Grant
		if ( !isfinite(critdose[i]) ) critdose[i] = 1e10;
	}
*/
	vector<double>	critexp = critical_exposure_curve(maxrad, real_size()[0], flag);

//	for ( i=0, s=0; i<maxrad; ++i, s+=ds )
//		cout << s << tab << critdose[i] << endl;
	
	FourierType		old_transform = fouriertype;
	
	if ( fouriertype == NoTransform )
		fft();
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		fspace_weigh_dose(nn, dose_per_frame, critexp);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; ++nn )
		fspace_weigh_dose(nn, dose_per_frame, critexp);
#endif
	
	if ( old_transform == NoTransform )
		fft_back();
		
	statistics();
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with the accumulated dose.
@param 	dose			array containing accumulated dose and parameters.
@return int				0.

	The amplitudes are weighed using:
	1. with 2 parameters: the formula of Grant and Grigorieff (2015).
	2. with 3 parameters: exponential decay.

**/
int			Bimage::fspace_weigh_accumulated_dose(vector<double> dose)
{
	if ( dose.size() < 1 ) return 0;

	vector<double>	rdcurve = radiation_damage_curve(dose, size().max(), real_size()[0]);

	fspace_weigh(rdcurve);
	
	return 0;
}

/**
@brief 	Weighs an image's amplitudes with the radial power spectrum of another.
@param 	*pref		reference image.
@param 	*pmask		reciprocal space mask (0 & 1, indicating inclusion of structure factors).
@param 	resolution	high resolution limit.
@return int			0.

	The two images must be of the same dimensions.
	The radial power spectra of the two images are calculated
	and used to calculate the ratio of the second to the first in each shell.
	This ratio is then used to rescale the amplitudes of the first image.

**/
int 		Bimage::fspace_weigh(Bimage* pref, Bimage* pmask, double resolution)
{
	if ( pmask ) pmask->change_type(UCharacter);
	
	check_resolution(resolution);
	
	//	maxrad is the main limit for the length of all the arrays,
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + real_size()[0]/resolution);
	
	if ( verbose ) {
		cout << "Weighing amplitudes:" << endl;
		cout << "Reference map:                  " << pref->file_name() << endl;
		cout << "Resolution:                     " << resolution << " A (" << maxrad << " pixels)" << endl << endl;
	}
	
	long			i, nn;

	if ( pref->compound_type() == TSimple && pref->fourier_type() == NoTransform )
		pref->fft();
	
	vector<double>	ref = pref->fspace_radial(0, maxrad, 1);
	
	if ( compoundtype == TSimple && fouriertype == NoTransform )
		fft();
	
	for ( nn=0; nn<n; nn++ ) {
		vector<double>		rps = fspace_radial(nn, maxrad, 1);
		for ( i=0; i<maxrad; i++ )
			if ( rps[i] ) rps[i] = ref[i]/rps[i];
		fspace_scale(nn, rps, pmask);
	}
	
	fft_back();

	return 0; 
} 


/**
@brief 	Normalizes an image's amplitudes.
@return int			0.

**/
int 		Bimage::fspace_normalize()
{ 
	if ( fouriertype != Standard ) {
		cerr << "Error: File " << file_name() << " must be a Fourier transform!" << endl;
		return -1; 
	}
	
	complex_power();
	
	long			nn, i, j, ds(x*y*z);
	double			scale;
	
	for ( nn=i=0; nn<n; ++nn ) {
		scale = 1/sqrt(image[nn].FOM());
		for ( j=0; j<ds; i++, j++ ) set(i, complex(i)*scale);
	}
		
	return 0;
} 

/**
@brief 	Normalizes an image's amplitudes.
@param 	*pmask		reciprocal space mask (0 & 1, indicating inclusion of structure factors).
@param 	resolution	high resolution limit.
@param  flag        flag to calculate power (0) or amplitude (1).
@return int			0.

**/
int 		Bimage::fspace_normalize_radial(Bimage* pmask, double resolution, int flag)
{ 
	if ( fouriertype != NoTransform ) {
		cerr << "Error: File " << file_name() << " must be a real space map!" << endl;
		return -1; 
	}
	
	if ( pmask ) pmask->change_type(UCharacter);
	
	check_resolution(resolution);
	
	//	maxrad is the main limit for the length of all the arrays,
	//	and is derived from the resolution set in the transform
	long			maxrad = (long) (2 + real_size()[0]/resolution);

	fft();
	
	long			i, nn;

	for ( nn=0; nn<n; nn++ ) {
		vector<double>	rps = fspace_radial(nn, maxrad, flag);
		for ( i=0; i<maxrad; i++ )
			if ( rps[i] ) rps[i] = 1/sqrt(rps[i]);
		fspace_scale(nn, rps, pmask);
	}
	
	fft_back();
	
	return 0;
} 

/**
@brief 	Sets the image to positive definite.
@return int			0.
Scales the amplitudes to give a zero (DC) value of one.
**/
int 		Bimage::fspace_positive()
{ 
	int				tflag(0);
	long			i(0), ds(x*y*z);
	double			fmax(0), scale;

	fmax = ds * (avg-min);
	
	if ( compoundtype != TComplex ) {
		tflag = 1;
		fft();
	}

	if ( (*this)[i] > SMALLFLOAT ) fmax = (*this)[i];
	if ( verbose ) {
		cout << "Converting to a positive definite image:" << endl;
		cout << "F0:                  " << fmax << endl << endl;
	}
	
	scale = 1/fmax;
	
	set(0, fmax);
	set(1, 0);
	
	for ( i=0; i<ds; i++ ) set(i, complex(i) * scale);
	
	if ( tflag ) fft_back();
	
	return 0;
} 

/**
@brief 	Checks Friedel symmetry.
@return double 		overall RMSD or residual.

	The differences between the complex and polar forms of Friedel-related 
	voxels are calculated and accumulated as squared sums weighted by their
	average intensities. The residuals are then calculated as 
	root-mean-square-deviations.

**/
double		Bimage::friedel_check()
{
	simple_to_complex();
	
	long			i, j, nn, xx, yy, zz, xf, yf, zf;
	double			d, w, s(0), R(0), Ri(0), Rr(0), Ra(0), Rp(0);
	Complex<double>	v, vf;

	for ( nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			zf = (zz)? z - zz: 0;
			for ( yy=0; yy<y; ++yy ) {
				yf = (yy)? y - yy: 0;
				for ( xx=0; xx<(x+1)/2; ++xx ) {
					xf = (xx)? x - xx: 0;
					i = index(xx, yy, zz, nn);
					j = index(xf, yf, zf, nn);
					v = complex(i);
					vf = complex(j);
					w = (v.power() + vf.power())/2;
					d = v.real() - vf.real();
					Rr += d*d*w;
					d = v.imag() + vf.imag();
					Ri += d*d*w;
					d = v.amp() - vf.amp();
					Ra += d*d*w;
					d = angle_set_negPI_to_PI(v.phi() + vf.phi());
					Rp += d*d*w;
					s += w;
				}
			}
		}
	}
	
	Rr = sqrt(Rr/s);
	Ri = sqrt(Ri/s);
	Ra = sqrt(Ra/s);
	Rp = sqrt(Rp/s);
	R = sqrt(Rr*Rr + Ri*Ri);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Friedel symmetry residuals:" << endl;
		cout << "Real:                           " << Rr << endl;
		cout << "Imaginary:                      " << Ri << endl;
		cout << "Amplitude:                      " << Ra << endl;
		cout << "Phase:                          " << Rp*180.0/M_PI << " degrees" << endl;
		cout << "Overall:                        " << R << endl << endl;
	}
	
	return R;
}

/**
@brief 	Applies Friedel symmetry.
@return int 		0.

	The Friedel-related voxels are converted to polar form and the 
	in difference in amplitude and phase calculated (taking into account
	the phases should differ by sign). The voxels are then set to the
	average and its conjugate.

**/
int 		Bimage::friedel_apply()
{
	long			i, j, nn, xx, yy, zz, xf, yf, zf, nR;
	double			amp, amps, phi, phis, damp, dphi;
	double			ampavg, phiavg, fomavg;
	double			Ramp, Rphi;
	double			tol(0.01);
	Complex<double>	v, vf;
	float*			fom = NULL;
	if ( next ) fom = (float *) next->data_pointer();

	if ( verbose & VERB_PROCESS ) {
		cout << "Applying Friedel symmetry:" << endl;
		cout << "Image\tRefl\tR(amp)\tR(phi)" << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Applying Friedel symmetry" << endl << endl;
	
	for ( nn=0; nn<n; nn++ ) {
		nR = 0;
		Ramp =  Rphi = 0;
		for ( zz=0; zz<z; ++zz ) {
			zf = (zz)? z - zz: 0;
			for ( yy=0; yy<y; ++yy ) {
				yf = (yy)? y - yy: 0;
				for ( xx=0; xx<(x+1)/2; ++xx ) {
					xf = (xx)? x - xx: 0;
					i = index(xx, yy, zz, nn);
					j = index(xf, yf, zf, nn);
					v = complex(i);
					vf = complex(j);
					amp = v.amp();
					amps = vf.amp();
					phi = v.phi();
					phis = vf.phi();
					damp = amp - amps;
					dphi = phi + phis;
					dphi = angle_set_negPI_to_PI(dphi);
					if ( verbose & VERB_DEBUG )
						if ( fabs(dphi) > tol ) 
							cout << "amp dphi: " << xx << " " << yy << " "
									<< zz << " " << j << " " << damp << " " << dphi << endl;
					ampavg = (amp + amps)/2;
					phiavg = phi - dphi/2;
					v = Complex<double>(ampavg*cos(phiavg), ampavg*sin(phiavg));
					set(i, v);
					set(j, v.conj());
			    	if ( fom ) {
						fomavg = (fom[i] + fom[j])/2;
						fom[i] = fomavg;
						fom[j] = fomavg;
					}
					Ramp += damp*damp;
					Rphi += dphi*dphi;
					nR++;
				}
			}
		}
		if ( nR ) {
			Ramp = sqrt(Ramp/nR);
			Rphi = sqrt(Rphi/nR);
		}
		if ( verbose & VERB_PROCESS )
			cout << nn+1 << tab << nR << tab << Ramp << tab << Rphi*180/M_PI << endl;
	}
	if ( verbose )
		cout << endl;
	
	return 0;
}

/**
@brief 	Calculates the cosine of the phase difference between two images.
@param 	*p			real space reference image.
@param	type		0=phase angle, 1=cos(phase angle), 2=scale by amplitude product
@param 	res_hi 		upper resolution limit.
@param 	res_lo 		lower resolution limit.
@return Bimage* 		phase difference image.

	Both images are Fourier transformed and the cosine of the phase
	difference calculated.

**/
Bimage* 	Bimage::phase_difference(Bimage* p, int type, double res_hi, double res_lo)
{
	if ( res_hi <= 0 ) res_hi = 0.1;
	if ( res_lo <= 0 ) res_lo = 1e10;
	if ( res_lo < res_hi + 1 ) res_lo = res_hi + 1;
	
	Bimage* 		ppd = pack_two_in_complex(p);
	
	ppd->fft();
	
//	if ( verbose & VERB_FULL )
	if ( verbose ) {
		if ( !(type&1) ) cout << "Phase angle";
		else if ( type&1 ) cout << "Cosine phase angle";
		cout << " difference" << endl;
		if ( type & 2 ) cout << "With amplitude scaling" << endl;
		cout << "Resolution range:               "
			<< res_hi << " - " << res_lo << endl << endl;
	}
	
	long	   		i, j, xx, yy, zz, nn;
	long			ix, iy, iz;
	Vector3<double>	shift;
	double			smax2 = 1.0/(res_hi*res_hi);
	double			smin2 = 1.0/(res_lo*res_lo);
	double			dy, dz, sx2, sy2, sz2, s2, w;
	Vector3<double>	scale(1.0/ppd->real_size());
	Complex<float>	temp1, temp2;
	
	float*			dphi = new float[datasize];
	
	if ( verbose )
		cout << "Image\tShift" << endl;
	for ( nn=0; nn<ppd->images(); nn++ ) {
		ppd->image[nn].origin(image[nn].origin() - p->image[nn].origin());
		shift = ppd->image[nn].origin()/ppd->size();
		if ( verbose )
			cout << nn+1 << tab << shift << endl;
		for ( zz=0; zz<ppd->z; ++zz ) {
			iz = -zz;
			if ( iz < 0 ) iz += ppd->z;
			dz = zz;
			if ( zz > (ppd->z - 1)/2 ) dz -= ppd->z;
			sz2 = dz*scale[2];
			sz2 *= sz2;
			for ( yy=0; yy<ppd->y; ++yy ) {
				iy = -yy;
				if ( iy < 0 ) iy += ppd->y;
				dy = yy;
				if ( yy > (ppd->y - 1)/2 ) dy -= ppd->y;
				sy2 = dy*scale[1];
				sy2 *= sy2;
				for ( xx=0; xx<ppd->x/2+1; ++xx ) {
					ix = -xx;
					if ( ix < 0 ) ix += ppd->x;
					sx2 = xx*scale[0];
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					if ( s2 >= smin2 && s2 <= smax2 ) {
						i = ppd->index(xx,yy,zz,nn);
						j = ppd->index(ix,iy,iz,nn);
						temp1 = ppd->complex(i).unpack_first(ppd->complex(j));
						temp2 = ppd->complex(i).unpack_second(ppd->complex(j));
						dphi[i] = temp1.phi() - temp2.phi() - TWOPI*(xx*shift[0] + yy*shift[1] + zz*shift[2]);
						if ( !(type & 1) ) dphi[i] = angle_set_negPI_to_PI(dphi[i]);
						if ( type & 1 ) dphi[i] = cos(dphi[i]);
						if ( type & 2 ) {
							w = temp1.power()*temp2.power();
							dphi[i] *= w/(w+1);
						}
						dphi[j] = dphi[i];
					}
				}
			}
		}
	}
	
	ppd->data_type(Float);
	ppd->compound_type(TSimple);
	ppd->channels(1);
	ppd->fourier_type(NoTransform);
	ppd->data_assign((unsigned char *) dphi);
	
	return ppd;
}

/**
@brief 	Calculates the average of the absolute phase difference between two images within given resolution shells.
@param 	*p				real space reference image.
@param 	res_hi 			upper resolution limit.
@param 	res_lo 			lower resolution limit.
@param 	weighting		weighting type: 0, none; 1, p1-amp; 2, p2-amp; 3, both amps
@return double	 		average of the cosine of the phase difference.

	Both images are Fourier transformed and the absolute phase
	difference calculated and averaged.

**/
double	 	Bimage::average_phase_difference(Bimage* p, double res_hi, double res_lo, int weighting)
{
	if ( res_hi <= 0 ) res_hi = 0.1;
	if ( res_lo < res_hi + 1 ) res_lo = res_hi + 1;
	
	Bimage* 		ppd = pack_two_in_complex(p);
	
	ppd->fft();
	
	long   		i, j, xx, yy, zz, nn;
	long			ix, iy, iz;
	Vector3<double>	shift;
	double			smax2 = 1.0/(res_hi*res_hi);
	double			smin2 = 1.0/(res_lo*res_lo);
	double			dy, dz, sx2, sy2, sz2, s2;
	Vector3<double>	scale(1.0/ppd->real_size());
	Complex<float>	temp1, temp2;
	
	double			amp, dphi, sum_amp(0), avg_pd(0);
	
	for ( nn=0; nn<ppd->images(); nn++ ) {
//		ppd->image[nn].origin(image[nn].origin() - p->image[nn].origin());
		ppd->image[nn].origin(image[nn].origin() - p->image[nn].origin());
		shift = ppd->image[nn].origin()/ppd->size();
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << shift << endl;
		for ( zz=0; zz<ppd->z; ++zz ) {
			iz = -zz;
			if ( iz < 0 ) iz += ppd->z;
			dz = zz;
			if ( zz > (ppd->z - 1)/2 ) dz -= ppd->z;
			sz2 = dz*scale[2];
			sz2 *= sz2;
			for ( yy=0; yy<ppd->y; ++yy ) {
				iy = -yy;
				if ( iy < 0 ) iy += ppd->y;
				dy = yy;
				if ( yy > (ppd->y - 1)/2 ) dy -= ppd->y;
				sy2 = dy*scale[1];
				sy2 *= sy2;
				for ( xx=0; xx<ppd->x/2+1; ++xx ) {
					ix = -xx;
					if ( ix < 0 ) ix += ppd->x;
					sx2 = xx*scale[0];
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					if ( s2 >= smin2 && s2 <= smax2 ) {
						i = ppd->index(xx,yy,zz,nn);
						j = ppd->index(ix,iy,iz,nn);
						temp1 = ppd->complex(i).unpack_first(ppd->complex(j));
						temp2 = ppd->complex(i).unpack_second(ppd->complex(j));
						dphi = temp1.phi() - temp2.phi() - TWOPI*(xx*shift[0] + yy*shift[1] + zz*shift[2]);
						dphi = angle_set_negPI_to_PI(dphi);
						switch ( weighting ) {
							case 1:
								amp = temp1.amp();
								break;
							case 2:
								amp = temp2.amp();
								break;
							case 3:
								amp = temp1.amp()*temp2.amp();
								break;
							default:						
								amp = 1e30;
						}
						amp /= amp + 1;
						sum_amp += amp;
						avg_pd += amp*fabs(dphi);
					}
				}
			}
		}
	}
	
	if ( sum_amp ) avg_pd /= sum_amp;
	else avg_pd = M_PI;
	
	delete ppd;
	
	return avg_pd;
}

/**
@brief 	Flips the phases of an image based on a phase difference map. 
@param 	*pd			reciprocal space phase difference map.
@return int			0.
**/
int			Bimage::phase_flip(Bimage* pd)
{
	fft();
	
	long		i, ds(x*y*z*n);
	
	for ( i=0; i<ds; i++ )
		if ( (*pd)[i] < 0 ) set(i, -complex(i));

	fft_back();
	
	return 0;
}

