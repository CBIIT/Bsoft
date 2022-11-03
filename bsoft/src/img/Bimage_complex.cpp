/**
@file	Bimage_complex.cpp
@brief	Routines to convert complex data sets
@author Bernard Heymann
@date	Created: 19990424
@date	Modified: 20220715
**/
	
#include "Bimage.h"
#include "timer.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief A simple image is converted to a complex image.

	The input image is written into the real part of the complex image.

**/
void		Bimage::simple_to_complex()
{
	if ( compoundtype == TComplex ) return;
	
	if ( compoundtype == TRGB || compoundtype == TRGBA || compoundtype == TCMYK ) color_to_simple();
	
	if ( compoundtype != TSimple ) {
		cerr << "Error: Conversion from compound types to complex type not supported!" << endl;
		exit(-1);
	}
	
	long			ds(x*y*z*n);
	Complex<float>*	cdata = new Complex<float>[ds];
	
	for ( long	 j=0; j<ds; j++ ) cdata[j] = (*this)[j];
	
	data_type(Float);
	compoundtype = TComplex;
	c = 2;
	
	data_assign((unsigned char *) cdata);
	
	statistics();
}

/**
@brief Two sub-images converted to a complex image.

**/
void		Bimage::two_to_complex()
{
	if ( compoundtype == TComplex ) return;
	
	if ( compoundtype != TSimple ) {
		cerr << "Error: Conversion from compound types to complex type not supported!" << endl;
		exit(-1);
	}
	
	if ( verbose )
		cout << "Converting from two images to a complex image" << endl << endl;
	
	long			i, j, is(x*y*z);
	Complex<float>*	cdata = new Complex<float>[is];
	
	for ( i=0, j=is; i<is; ++i, ++j ) cdata[i] = Complex<float>((*this)[i], (*this)[j]);
	
	data_type(Float);
	compoundtype = TComplex;
	n = 1;
	c = 2;
	
	data_assign((unsigned char *) cdata);
	
	statistics();
}

/**
@brief A multi-channel image is converted to a set of complex images.

	The input image channels are written into the real part of the complex image.

**/
void		Bimage::multi_channel_to_complex()
{
	if ( compoundtype == TComplex ) return;
	
	change_type(Float);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::multi_channel_to_complex: from " << compoundtype << endl;
	
	long			i, j, k, nn, cc, isz(x*y*z);
	Complex<float>*	cdata = new Complex<float>[c*isz*n];
	
	for ( nn=i=0; nn<n; ++nn )
		for ( j=0; j<isz; ++j )
			for ( cc=0, k=nn*isz+j; cc<c; ++cc, ++i, k+=isz )
				cdata[k] = (*this)[i];

	Bsub_image*		nusub = new Bsub_image[n*c];
	
	if ( image ) {
		for ( nn=j=0; nn<n; ++nn )
			for ( cc=0; cc<c; ++cc, ++j )
				nusub[j] = image[nn];
		delete[] image;
	}

	image = nusub;
	
	compoundtype = TComplex;
	n *= c;
	c = 2;
	
	data_assign((unsigned char *) cdata);
	
	statistics();
}

/**
@brief 	A complex image is split into two simple images.
@return	Bimage*		new image with the imaginary part.

	The imaginary part is written into the new image.
	The complex image is replaced by its real part.

**/
Bimage*		Bimage::complex_split()
{
	if ( compoundtype != TComplex ) {
		cerr << "Error in complex_split: The image must be complex!" << endl;
		bexit(-1);
	}
	
	Bimage*			pim = copy_header();
	pim->compound_type(TSimple);
	pim->channels(1);
	pim->data_alloc();
	
	long			ds(x*y*z*n);
	
	for ( long j=0; j<ds; ++j )
		pim->set(j, complex(j).imag());
	
	complex_to_real();

	pim->statistics();

	return pim;
}

/**
@brief A phase image is converted to a complex image.

**/
void		Bimage::phase_to_complex()
{
	simple_to_complex();
	
	long			ds(x*y*z*n);
	
	for ( long j=0; j<ds; j++ )
		set(j, Complex<double>(cos((*this)[j]), sin((*this)[j])));
}

/**
@brief The real part of a complex image is written to a simple image.
**/
void		Bimage::complex_to_real()
{
	if ( compoundtype != TComplex ) return;
	
	long				j, k, ds(x*y*z*n);
	float*				fdata = new float[ds];
	
	for ( j=k=0; j<ds; j++, k+=2 ) fdata[j] = (*this)[k];
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief The imaginary part of a complex image is written to a simple image.
**/
void		Bimage::complex_to_imaginary()
{
	if ( compoundtype != TComplex ) return;
	
	long				j, k, ds(x*y*z*n);
	float*			fdata = new float[ds];
	
	for ( j=0, k=1; j<ds; j++, k+=2 ) fdata[j] = (*this)[k];
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief The intensities from a complex image is written to a simple image.
**/
void		Bimage::complex_to_intensities()
{
	if ( compoundtype != TComplex ) return;
	
	long			j, ds(x*y*z*n);
	float*			fdata = new float[ds];
	
	for ( j=0; j<ds; j++ ) fdata[j] = (complex(j)).power();
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief The amplitudes from a complex image is written to a simple image.

**/
void		Bimage::complex_to_amplitudes()
{
	if ( compoundtype != TComplex ) return;
	
	long			j, ds(x*y*z*n);
	float*			fdata = new float[ds];
	
	for ( j=0; j<ds; j++ ) fdata[j] = (complex(j)).amp();
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief The signed amplitudes from a complex image is written to a simple image.

**/
void		Bimage::complex_to_signed_amplitudes()
{
	if ( compoundtype != TComplex ) return;
	
	long				j, ds(x*y*z*n);
	float*			fdata = new float[ds];
	
	for ( j=0; j<ds; j++ ) {
		fdata[j] = (complex(j)).amp();
		if ( (complex(j)).phi() < 0 ) fdata[j] = -fdata[j];
	}
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief The phases from a complex image is written to a simple image.

**/
void		Bimage::complex_to_phases()
{
	if ( compoundtype != TComplex ) return;
	
	long			j, ds(x*y*z*n);
	float*			fdata = new float[ds];
	
	for ( j=0; j<ds; j++ ) fdata[j] = (complex(j)).phi();
//	for ( j=0; j<ds; j++ ) fdata[j] = (complex(j)).phi_lut();
	
	data_type(Float);
	compoundtype = TSimple;
	c = 1;
	
	data_assign((unsigned char *) fdata);
	
	statistics();
}

/**
@brief Calculates the complex conjugate image.

**/
void 		Bimage::complex_conjugate()
{
	if ( compoundtype != TComplex ) return;
	
	long			j;
	
	for ( j=1; j<datasize; j+=2 ) set(j, -(*this)[j]);
}

/**
@brief Calculates the power in a complex image.
@return double			power.

	The first element is set to zero.
	The image is scaled by the remaining total power, i.e., the sum of the intensities.

**/
double 		Bimage::complex_power()
{
	if ( compoundtype != TComplex ) return 0;
	
	long			nn, i, j, ds(x*y*z);
	double			pwr(0), pwrn;
	
	for ( nn=i=0; nn<n; nn++ ) {
		pwrn = 0;
		for ( j=0; j<ds; i++, j++ ) pwrn += (complex(i)).power();
		pwr += pwrn;
		image[nn].FOM(pwrn);
	}
	
	return pwr;
}

/**
@brief Normalizes the power in a complex image.
@return double		power.

	The first element is set to zero.
	The image is scaled by the remaining total power, i.e., the sum of the intensities.

**/
double 		Bimage::complex_normalize()
{
	if ( compoundtype != TComplex ) return 0;
	
	set(0, 0); set(1, 0);
	
	long			nn, i, j, ds(x*y*z);
	double			pwr = complex_power(), pwrn, f;
	
	for ( nn=0; nn<n; nn++ ) {
		pwrn = image[nn].FOM();
		if ( pwrn ) {
			f = 1/sqrt(pwrn);
			for ( j=0, i=nn*ds; j<ds; i++, j++ ) set(i, complex(i) * f);
		}
	}
	
	statistics();
	
	return pwr;
}

/**
@brief 	Inverts a complex image.
@return int			0.

	The image is scaled by the remaining total power, i.e., the sum of the intensities.

**/
int 		Bimage::complex_invert()
{
	if ( compoundtype != TComplex ) return 0;
	
	for ( long i=0; i<data_size(); ++i )
		set(i, -complex(i));
	
	statistics();
	
	return 0;
}

/**
@brief 	Converts a complex image to a simple image.
@param 	conv		a flag for converting the complex image.
@return int			0.

**/
int 		Bimage::complex_convert(ComplexConversion conv)
{
	if ( compoundtype != TComplex ) return 0;
	
	switch ( conv ) {
		case NoConversion: break;
		case Real: complex_to_real(); break;
		case Imaginary: complex_to_imaginary(); break;
		case Amplitude: complex_to_amplitudes(); break;
		case Intensity: complex_to_intensities(); break;
	}
	
	if ( conv ) fourier_type(NoTransform);
	
	return 0;
}

/**
@brief 	Phase shifts a complex image.
@param 	shift		three-value real space shift vector.
@return int			error code.

	A real space translation with wrapping is equivalent to phase shifting
	in reciprocal space.

**/
int 		Bimage::phase_shift(Vector3<double> shift)
{
	if ( compoundtype != TComplex ) return -1;
	
	if ( shift.length() == 0 ) return 0;
	
	long				nn;
	
	for ( nn=0; nn<n; nn++ )
		phase_shift(nn, shift);
	
	return 0;
}

/**
@brief 	Phase shifts a complex sub-image.
@param	nn			sub-image to transform.
@param 	shift		three-value real space shift vector.
@return int			error code.

	A real space translation with wrapping is equivalent to phase shifting
	in reciprocal space.

**/
int 		Bimage::phase_shift(long nn, Vector3<double> shift)
{
	if ( compoundtype != TComplex ) return -1;
	
	if ( shift.length() == 0 ) return 0;
	
	if ( verbose & VERB_FULL )
		cout << "Translate within unit cell by:  " << shift << endl;
	
	long				j(nn*x*y*z), xx, yy, zz;
	long				h, k, l;
	double				phi(0);
	Vector3<double>		half((x - 1)/2, (y - 1)/2, (z - 1)/2);
	Vector3<double>		t(shift/size() * MIN2PI);
	Complex<double>		cv;
	
	for ( zz=0; zz<z; zz++ ) {
		l = zz;
		if ( l > half[2] ) l -= (long)z;
		for ( yy=0; yy<y; yy++ ) {
			k = yy;
			if ( k > half[1] ) k -= (long)y;
			for ( xx=0; xx<x; xx++, j++ ) {
				h = xx;
				if ( h > half[0] ) h -= (long)x;
				phi = h*t[0] + k*t[1] + l*t[2];
				cv = complex(j);
				cv.shift_phi(phi);
				set(j, cv);
			}
		}
	}
	
	image[nn].origin(image[nn].origin() + shift);
	
	return 0;
}

/**
@brief Phase shifts a set of reflections to the image origin.
@return int 			0.

	A real space translation with wrapping is equivalent to phase shifting
	in reciprocal space. The phases are shifted based on the embedded
	sub-image origins.

**/
int 		Bimage::phase_shift_to_origin()
{
	if ( compoundtype != TComplex ) return -1;
	
	if ( verbose & VERB_FULL )
		cout << "Translate within unit cell to phase origin:" << endl;
	
	long				j, nn, xx, yy, zz;
	long				h, k, l;
	double				skl, sl, phi;
	Vector3<double>		shift, half((x - 1)/2, (y - 1)/2, (z - 1)/2);
	Complex<double>		cv;
	
	if ( verbose & VERB_FULL )
		cout << "Image\tox\toy\toz" << endl;
	for ( j=nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin()[0] << tab << image[nn].origin()[1] << tab << image[nn].origin()[2] << endl;
		shift = (image[nn].origin()/size()) * TWOPI;
		for ( zz=0; zz<z; zz++ ) {
			l = zz;
			if ( l > half[2] ) l -= (long)z;
			sl = l*shift[2];
			for ( yy=0; yy<y; yy++ ) {
				k = yy;
				if ( k > half[1] ) k -= (long)y;
				skl = sl + k*shift[1];
				for ( xx=0; xx<x; xx++, j++ ) {
					h = xx;
					if ( h > half[0] ) h -= (long)x;
					phi = h*shift[0] + skl;
					cv = complex(j);
					cv.shift_phi(phi);
					set(j, cv);
				}
			}
		}
//		image[nn].origin(image[nn].origin() + shift*size());
		image[nn].origin(0.0,0.0,0.0);
	}
	if ( verbose & VERB_FULL )
		cout << endl;
	
	return 0;
}

/**
@brief 	Phase shifts a set of reflections to the nominal center of the image origin.
@return int 		0.

	A real space translation with wrapping is equivalent to phase shifting
	in reciprocal space. The phases are shifted based on the embedded
	sub-image origins.

**/
int 		Bimage::phase_shift_to_center()
{
	Vector3<double>	shift(x/2, y/2, z/2);

	phase_shift(shift);
	
	return 0;
}

/**
@brief 	Calculates the product of a complex image with a simple image.
@param 	*p			simple image.
@return int			error code.

	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_multiply(Bimage* p)
{
	if ( !d.uc ) return -1;
	
	if ( compoundtype != TComplex || p->compound_type() != TSimple ) return -1;
	
	change_type(Float);
	p->change_type(Float);
		
	if ( verbose & VERB_FULL )
		cout << "Multiplying the complex image" << endl << endl;
	
	long				i, ds(x*y*z*n);
	
	for ( i=0; i<ds; i++ ) set(i, complex(i) * (*p)[i]);
	
	return 0;
}

/**
@brief 	Calculates the complex product of two complex images.
@param 	*p			complex image.
@return int			error code.

	Complex product:
		(a + ib)*(c + id) = (a*c - b*d) + i(b*c + a*d).
	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_product(Bimage* p)
{
	if ( !d.uc ) return -1;
	
	if ( compoundtype != TComplex || p->compound_type() != TComplex ) return -1;
	
	change_type(Float);
	p->change_type(Float);
		
	if ( verbose & VERB_FULL )
		cout << "Calculating the complex product" << endl << endl;
	
	long				i, ds(x*y*z*n);
	
	for ( i=0; i<ds; i++ ) set(i, complex(i) * p->complex(i));
	
	return 0;
}

/**
@brief 	Calculates the product of a complex image with the conjugate of a second.
@param 	*p			complex image.
@param 	norm		normalize.
@return int			error code.

	Complex conjugate product:
		(a + ib)*(c - id) = (a*c + b*d) + i(b*c - a*d).
	Normalization is determined by the norm flag:
		1	square root of the power sum product of the two images
		2	the power sum of the first image
	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_conjugate_product(Bimage* p, int norm)
{
	if ( !d.uc ) return -1;
	
	if ( compound_type() != TComplex || p->compound_type() != TComplex ) return -1;
	
	change_type(Float);
	p->change_type(Float);
		
	if ( verbose & VERB_FULL )
		cout << "Calculating the complex conjugate product" << endl << endl;
	
	long				i, j, nn, ds(x*y*z);
	double				sum1(0), sum2(0), scale(0);
	Complex<double>		cv;
	
    for ( i=nn=0; nn<n; ++nn ) {
		sum1 = sum2 = 0;
		for ( j=0; j<ds; ++j, ++i ) {
			cv = p->complex(i).conj();
			if ( norm ) {
				sum1 += complex(i).power();
				sum2 += cv.power();
			}
			set(i, complex(i) * cv);
		}
		if ( norm ) {
			if ( norm == 1 ) {
				scale = sum1*sum2;
				if ( scale > 0 ) scale = 1.0/sqrt(scale);
			} else if ( norm == 2 ) {
				scale = sum1;
				if ( scale > 0 ) scale = 1.0/scale;
			}
			if ( scale > 0 ) {
				multiply(nn, scale);
			} else {
				cerr << "Error in Bimage::complex_conjugate_product: 	scaling failed!" << endl;
				return -1;
			}
		}
	}

	return 0;
}

/**
@brief 	Calculates the product of a complex image with the conjugate of a second.
@param 	*p			complex image.
@return int			error code.

	Complex conjugate product:
		(a + ib)*(c - id) = (a*c + b*d) + i(b*c - a*d).
	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
Bimage*		Bimage::complex_conjugate_product_one2many(Bimage* p)
{
	if ( !d.uc ) return NULL;
	
	if ( compound_type() != TComplex || p->compound_type() != TComplex ) return NULL;

	if ( !check_if_same_image_size(p) ) return NULL;
	
	Bimage*			pc = p->copy();
		
	if ( verbose & VERB_FULL )
		cout << "Calculating the complex conjugate product" << endl << endl;
	
	long			i, j, nn, ds(x*y*z);
	
	for ( nn=j=0; nn<p->n; nn++ )
		for ( i=0; i<ds; i++, j++ ) pc->set(j, complex(i) * (p->complex(j)).conj());
	
	return pc;
}

/**
@brief 	Applies a mask to a complex image.
@param 	*pmask		mask (converted to floating point).
@return int			error code.

	Wherever the mask is <= 0, the data is set to zero.
	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_apply_mask(Bimage* pmask)
{
	if ( !d.uc ) return -1;
	if ( !pmask ) return -1;
	
	pmask->change_type(Float);
	
	if ( verbose & VERB_FULL )
		cout << "Masking a complex image" << endl;
	
	long				i, ds(x*y*z*n);
	Complex<double>	cv;
	
	for ( i=0; i<ds; i++ )
		if ( (*pmask)[i] <= 0 ) set(i, cv);
	
	return 0;
}

/**
@brief 	Applies a negative mask to a complex image.
@param 	*pmask		mask (converted to floating point).
@return int			error code.

	Wherever the mask is >= 0, the data is set to zero.
	Requirement: The two images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_apply_negative_mask(Bimage* pmask)
{
	if ( !d.uc ) return -1;
	if ( !pmask ) return -1;
	
	if ( verbose & VERB_FULL )
		cout << "Masking a complex image" << endl;
	
	long				i, ds(x*y*z*n);
	Complex<double>	cv;
	
	for ( i=0; i<ds; i++ )
		if ( (*pmask)[i] >= 0 ) set(i, cv);
	
	return 0;
}

/**
@brief 	Applies a dual mask to a complex image.
@param 	*pmask		mask (converted to floating point).
@return int			error code.

	The mask is applied to an image, generating two images, the first 
	with values retained where the mask is positive, and the second with
	values where the mask is negative.
	The input image is replaced and the second is created and linked to
	the first.
	Requirement: The two input images must be the same size.
	No statistics are calculated.

**/
int 		Bimage::complex_apply_dual_mask(Bimage* pmask)
{
	if ( !d.uc ) return -1;
	if ( !pmask ) return -1;
	
	pmask->change_type(Float);
	
	if ( verbose & VERB_FULL )
		cout << "Masking a complex image with a dual mask" << endl;
	
	Bimage*			p = copy();
	next = p;
	
	long				i, ds(x*y*z*n);
	Complex<double>	cv;
	
	for ( i=0; i<ds; i++ ) {
		if ( (*pmask)[i] <= 0 ) set(i, cv);
		if ( (*pmask)[i] >= 0 ) p->set(i, cv);
	}
	
	return 0;
}

/**
@brief 	Bandpasses a Fourier transform.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return int			error code.

	Requirement: a complex Fourier transform.
	All data outside the band limits are set to zero.

**/
int 		Bimage::complex_bandpass(double hires, double lores)
{
	if ( !d.uc ) return -1;

	if ( compoundtype != TComplex ) return -1;
	
	change_type(Float);
	
	long	 		i, nn, xx, yy, zz;
	long 			ix, iy, iz;
	double			sx2, sy2, sz2, s2;
	
	check_resolution(hires);
	if ( lores <= 0 ) lores = 1e30;
	if ( lores < hires ) swap(lores, hires);
	double			s2hi = 1/(hires*hires);
	double			s2lo = 1/(lores*lores);
	Vector3<double>	iscale(1.0/real_size());
	Vector3<long>	h((x - 1)/2, (y - 1)/2, (z - 1)/2);
	Complex<double>	cv;
	
	if ( verbose & VERB_FULL )
		cout << "Bandpass limit to:              " << 1/sqrt(s2hi) << " - " << 1/sqrt(s2lo) << " Å" << endl;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::complex_bandpass: hires=" << hires << " lores" << lores << endl;
	
    for ( nn=i=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			iz = zz;
			if ( zz > h[2] ) iz -= (long)z;
			sz2 = iz*iscale[2];
			sz2 *= sz2;
			for ( yy=0; yy<y; yy++ ) {
				iy = yy;
				if ( yy > h[1] ) iy -= (long)y;
				sy2 = iy*iscale[1];
				sy2 *= sy2;
				for ( xx=0; xx<x; xx++, i++ ) {
					ix = xx;
					if ( xx > h[0] ) ix -= (long)x;
					sx2 = ix*iscale[0];
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					if ( s2 > s2hi || s2 < s2lo ) set(i, cv);
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Packs two real images into one complex image.
@param 	*p			second image.
@return Bimage*		the new complex image, NULL on error.

	Two real space images are packed into the real and imaginary parts
	of a new data block with conversion from the original non-complex type.
	The header values of the first image are adopted for the new image.
	No statistics are calculated.

**/
Bimage*		Bimage::pack_two_in_complex(Bimage* p)
{
	if ( !d.uc ) {
		error_show("Error in Bimage::pack_two_in_complex", __FILE__, __LINE__);
		cerr << "No data for " << file_name() << endl;
		return NULL;
	}
	
	if ( !p->data_pointer() ) {
		error_show("Error in Bimage::pack_two_in_complex", __FILE__, __LINE__);
		cerr << "No data for " << p->file_name() << endl;
		return NULL;
	}
	
    if ( compound_type() > TSimple ) {
		error_show("Error in Bimage::pack_two_in_complex", __FILE__, __LINE__);
		cerr << "Data for " << file_name() << " may not be compound!" << endl;
		return NULL;
	}
	
    if ( p->compound_type() > TSimple ) {
		error_show("Error in Bimage::pack_two_in_complex", __FILE__, __LINE__);
		cerr << "Data for " << p->file_name() << " may not be compound!" << endl;
		return NULL;
	}
	
	check_if_same_size(p);
		
    long			   	i, j;
	
	Bimage*				pc = new Bimage(Float, TComplex, size(), n);

	if ( !pc ) {
		error_show("Bimage::pack_two_in_complex", __FILE__, __LINE__);
		return NULL;
	}

	pc->fourier_type(NoTransform);
	for ( i=0; i<n; i++ ) pc->image[i] = image[i];

	if ( verbose & VERB_FULL )
    	cout << "Packing two real images into a complex image" << endl;

	for ( i=j=0; i<datasize; i++ ) {
		pc->set(j++, (*this)[i]);
		pc->set(j++, (*p)[i]);
	}
	
	pc->check();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::pack_two_in_complex: std=" << pc->standard_deviation() << endl;
	
	return pc;
}

/**
@brief 	Unpacks a complex transform obtained from two real images.
@return Bimage*		transform of second image.

	The complex image must be a Fourier transform obtained from two real 
	space images which were packed into the real and imaginary parts
	of the complex image before Fourier transformation.
	The input image is used to hold the transform of the first image.
	A new image is created to hold the transform of the second image.
	Both these images are complex floating point.
	Note: Images with even dimensions cannot be unpacked exactly because
		the inverse is not present when x=nx/2 or y=ny/2 or z=nz/2.

**/
Bimage* 	Bimage::unpack_combined_transform()
{
	if ( compoundtype != TComplex ) return NULL;
	
	change_type(Float);
	
	int				use;
	long	   		i, j, xx, yy, zz, nn, xh(x/2 + 1);
	long				ix, iy, iz;
	Vector3<long>	h(size()/2);
	Complex<float>	temp1, temp2, cv;
	
	Bimage* 		p = copy();

	if ( verbose & VERB_FULL )
		cout << "Unpacking a combined Fourier transform" << endl << endl;
	
	for ( nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) { 
			iz = (zz>0)? z - zz: 0;
			for ( yy=0; yy<y; yy++ ) { 
				iy = (yy>0)? y - yy: 0;
				for ( xx=0; xx<xh; xx++ ) {
					ix = (xx>0)? x - xx: 0;
					use = 1;
					if ( xx == ix ) {
						if ( yy > h[1] ) use = 0;
						if ( yy == iy && zz > h[2] ) use = 0;
					}
					if ( use ) {
						i = index(xx, yy, zz, nn);
						j = index(ix, iy, iz, nn);
						temp1 = complex(i).unpack_first(complex(j));
						temp2 = complex(i).unpack_second(complex(j));
						set(i, temp1);
						set(j, temp1.conj());
						p->set(i, temp2);
						p->set(j, temp2.conj());
					}
				}
			}
		}
	}
	
	return p;
}

/**
@brief 	Calculates the complex conjugate product of a complex image resulting from combining and Fourier transforming two real space images.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms, yielding 
	the correlation coefficient when the product is backtransformed into 
	the cross-correlation map.

**/
int 		Bimage::combined_complex_product()
{
	if ( !d.uc ) return -1;
	
	if ( compound_type() != TComplex ) return -1;
	
	change_type(Float);

	int 			use;
	long	 		i, j, nn, xx, yy, zz, xh(x/2 + 1);
	long	 		ix, iy, iz;
	Vector3<long>	h(size()/2);
		
	double			scale(1), sum1, sum2;
	Complex<double>	temp1, temp2;
	
	if ( verbose & VERB_FULL )
		cout << "Calculating the complex conjugate product" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::combined_complex_product: F0=" << complex(0).real() << " " << complex(0).imag() << endl;

	double			tvs = getwalltime(), tvf;
	
    for ( nn=0; nn<n; nn++ ) {
		sum1 = sum2 = 0;
		for ( zz=0; zz<z; zz++ ) {
			iz = (zz>0)? z - zz: 0;
			for ( yy=0; yy<y; yy++ ) {
				iy = (yy>0)? y - yy: 0;
				for ( xx=0; xx<xh; xx++ ) {
					ix = (xx>0)? x - xx: 0;
					use = 1;
					if ( xx == ix ) {
						if ( yy > h[1] ) use = 0;
						if ( yy == iy && zz > h[2] ) use = 0;
					}
					if ( use ) {
						i = index(xx, yy, zz, nn);
						j = index(ix, iy, iz, nn);
						temp1 = complex(i).unpack_first(complex(j));
						temp2 = complex(i).unpack_second(complex(j));
						set(i, temp1 * temp2.conj());
						set(j, complex(i).conj());
						sum1 += temp1.power();
						sum2 += temp2.power();
					}
				}
			}
		}
		if ( sum1 && sum2 ) {
			scale = 0.5/sqrt(sum1*sum2);
			multiply(nn, scale);
		}
	}

	if ( !isfinite(scale) ) {
		error_show("Error in Bimage::combined_complex_product", __FILE__, __LINE__);
		cerr << "Scale not finite = " << scale << endl;
		return -9;
	}
	
	if ( verbose & VERB_TIME ) {
		tvf = getwalltime();
		cout << "CC time: " << tvf - tvs << endl;
	}
	
	return 0;
}

/**
@brief 	Calculates the complex conjugate product of a complex image resulting from combining and Fourier transforming two real space images.
@param 	*pmask		binary mask (only 0 and 1), NULL if not desired.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms, yielding 
	the correlation coefficient when the product is backtransformed into 
	the cross-correlation map.

**/
int 		Bimage::combined_complex_product(Bimage* pmask)
{
	return combined_complex_product(0, 0, pmask);
}


/**
@brief 	Calculates the complex conjugate product of a complex image resulting from combining and Fourier transforming two real space images.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms within
	the resolution limits, yielding the correlation coefficient when
	the product is backtransformed into the cross-correlation map.

**/
/*int 		Bimage::combined_complex_product(double hires, double lores)
{
	return combined_complex_product(hires, lores, NULL);
}
*/
/**
@brief 	Calculates the complex conjugate product of a complex image resulting from combining and Fourier transforming two real space images.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		binary mask (only 0 and 1), NULL if not desired.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	The result is scaled by the total power of the two transforms within
	the resolution limits, yielding the correlation coefficient when
	the product is backtransformed into the cross-correlation map.

**/
int 		Bimage::combined_complex_product(double hires, double lores, Bimage* pmask)
{
	if ( !d.uc ) return -1;
	
	if ( compound_type() != TComplex ) return -1;
	
	change_type(Float);
	
	int 			use;
	long	 		i, j, nn, xx, yy, zz, xh(x/2 + 1);
	long	 		ix, iy, iz;
	double			sx2(0), sy2(0), sz2(0), s2;
	
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	double			s2hi = (hires > 0)? 1/(hires*hires): 0;
	double			s2lo = (lores > 0)? 1/(lores*lores): 0;

	Vector3<double>	iscale(1.0/real_size());
	Vector3<long>	h(size()/2);
	
	if ( verbose & VERB_FULL ) {
		cout << "Calculating the complex conjugate product:" << endl;
		cout << "Resolution range:               ";
		if ( hires > 0 ) cout << 1/sqrt(s2hi) << " - ";
		else cout << "0 - ";
		if ( lores > 0 ) cout << 1/sqrt(s2lo) << " A" << endl;
		else cout << "inf A" << endl;
		if ( pmask ) cout << "Mask:                   " << pmask->file_name() << endl;
	}
		
	double			scale(1), sum1, sum2;
	Complex<double>	temp1, temp2, cv(0,0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::combined_complex_product: F0 = " << complex(0).real() << " " << complex(0).imag() << endl;
	
	int				err(0);
	double			tvs = getwalltime(), tvf;
		
    for ( nn=0; nn<n; nn++ ) {
		sum1 = sum2 = 0;
		for ( zz=0; zz<z; zz++ ) {
			iz = (zz>0)? z - zz: 0;
			if ( s2hi > 0 ) {
				sz2 = zz;
				if ( zz > h[2] ) sz2 -= (double)z;
				sz2 *= iscale[2];
				sz2 *= sz2;
			}
			for ( yy=0; yy<y; yy++ ) {
				iy = (yy>0)? y - yy: 0;
				if ( s2hi > 0 ) {
					sy2 = yy;
					if ( yy > h[1] ) sy2 -= (double)y;
					sy2 *= iscale[1];
					sy2 *= sy2;
				}
				for ( xx=0; xx<xh; xx++ ) {
					ix = (xx>0)? x - xx: 0;
					use = 1;
					if ( xx == ix ) {
						if ( yy > h[1] ) use = 0;
						if ( yy == iy && zz > h[2] ) use = 0;
					}
					if ( use ) {
						i = index(xx, yy, zz, nn);
						j = index(ix, iy, iz, nn);
						if ( s2hi > 0 ) {
							sx2 = xx;
							if ( xx > h[0] ) sx2 -= (double)x;
							sx2 *= iscale[0];
							sx2 *= sx2;
							s2 = sx2 + sy2 + sz2;
							if ( s2 < s2lo || s2 > s2hi ) use = 0;
						}
						if ( pmask && (*pmask)[i] <= 0 ) use = 0;
						if ( use ) {
							temp1 = complex(i).unpack_first(complex(j));
							temp2 = complex(i).unpack_second(complex(j));
							set(i, temp1 * temp2.conj());
							set(j, complex(i).conj());
// Alternative calculation
//							temp1 = complex(i);
//							temp2 = complex(j);
//							temp.real((temp1.re()*temp2.im() + temp1.im()*temp2.re()));
//							temp.imag(temp1.power() - temp2.power());
//							set(i, temp);
//							set(j, temp.conj());
							sum1 += temp1.power();
							sum2 += temp2.power();
						} else {
							set(i, cv);
							set(j, cv);
						}
					}
				}
			}
		}
		scale = sum1*sum2;
		if ( scale > 0 ) {
			scale = 0.5/sqrt(scale);
			multiply(nn, scale);
/*			statistics();
			if ( image->maximum() > 1 ) {
				cerr << "Error in Bimage::combined_complex_product: maximum too large" << image->maximum() << endl;
				bexit(-1);
			}*/
		} else {
			cerr << "Error in Bimage::combined_complex_product: scaling failed!" << endl;
			err--;
		}
	}

	if ( verbose & VERB_TIME ) {
		tvf = getwalltime();
		cout << "CC(BP) time: " << tvf - tvs << endl;
	}
	
	if ( err ) {
		error_show("Error in Bimage::combined_complex_product", __FILE__, __LINE__);
		cerr << "Number of scaling errors: " << err << " (" << n << ")" << endl;
	}
	
	return err;
}

/**
@brief 	Calculates the phase product of a complex image resulting from combining and Fourier transforming two real space images.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		binary mask (only 0 and 1), NULL if not desired.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the
	phase product .
	The result is scaled by the size, yielding the phase the correlation coefficient when
	the product is backtransformed into the phasecorrelation map.

**/
int			Bimage::combined_phase_product(double hires, double lores, Bimage* pmask)
{
	if ( !d.uc ) return -1;
	
	if ( compound_type() != TComplex ) return -1;
	
	change_type(Float);
	
	int 			use;
	long	 		i, j, nn, xx, yy, zz, xh(x/2 + 1);
	long	 		ix, iy, iz;
	double			sx2(0), sy2(0), sz2(0), s2;
	
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	double			s2hi = (hires > 0)? 1/(hires*hires): 0;
	double			s2lo = (lores > 0)? 1/(lores*lores): 0;

	Vector3<double>	iscale(1.0/real_size());
	Vector3<long>	h(size()/2);
	
	if ( verbose & VERB_FULL ) {
		cout << "Calculating the complex phase product:" << endl;
		cout << "Resolution range:               ";
		if ( hires > 0 ) cout << 1/sqrt(s2hi) << " - ";
		else cout << "0 - ";
		if ( lores > 0 ) cout << 1/sqrt(s2lo) << " A" << endl;
		else cout << "inf A" << endl;
		if ( pmask ) cout << "Mask:                   " << pmask->file_name() << endl;
	}
	
	double			sd(1e-30);
	Complex<double>	temp1, temp2, temp, cv(0,0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::combined_complex_product: F0 = " << complex(0).real() << " " << complex(0).imag() << endl;
	
	int				err(0);
	double			tvs = getwalltime(), tvf;
	
    for ( nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			iz = (zz>0)? z - zz: 0;
			if ( s2hi > 0 ) {
				sz2 = zz;
				if ( zz > h[2] ) sz2 -= (double)z;
				sz2 *= iscale[2];
				sz2 *= sz2;
			}
			for ( yy=0; yy<y; yy++ ) {
				iy = (yy>0)? y - yy: 0;
				if ( s2hi > 0 ) {
					sy2 = yy;
					if ( yy > h[1] ) sy2 -= (double)y;
					sy2 *= iscale[1];
					sy2 *= sy2;
				}
				for ( xx=0; xx<xh; xx++ ) {
					ix = (xx>0)? x - xx: 0;
					use = 1;
					if ( xx == ix ) {
						if ( yy > h[1] ) use = 0;
						if ( yy == iy && zz > h[2] ) use = 0;
					}
					if ( use ) {
						i = index(xx, yy, zz, nn);
						j = index(ix, iy, iz, nn);
						if ( s2hi > 0 ) {
							sx2 = xx;
							if ( xx > h[0] ) sx2 -= (double)x;
							sx2 *= iscale[0];
							sx2 *= sx2;
							s2 = sx2 + sy2 + sz2;
							if ( s2 < s2lo || s2 > s2hi ) use = 0;
						}
						if ( pmask && (*pmask)[i] <= 0 ) use = 0;
						if ( use ) {
							temp1 = complex(i).unpack_first(complex(j));
							temp2 = complex(i).unpack_second(complex(j));
							temp = temp1 * temp2.conj();
							set(i, temp/(temp.amp()+sd));
							set(j, complex(i).conj());
						} else {
							set(i, cv);
							set(j, cv);
						}
					}
				}
			}
		}
		multiply(nn, 0.5L/size().volume());
	}

	if ( verbose & VERB_TIME ) {
		tvf = getwalltime();
		cout << "CC(BP) time: " << tvf - tvs << endl;
	}
	
	if ( err ) {
		error_show("Error in Bimage::combined_phase_product", __FILE__, __LINE__);
		cerr << "Number of scaling errors: " << err << " (" << n << ")" << endl;
	}
	
	return err;
}

/**
@brief 	Calculates the complex conjugate product of a complex image resulting from combining and Fourier transforming two real space images.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return int			error code.

	Requirement: Fourier transform of two images packed into one complex
		data block with the function Bimage::pack_two_in_complex and then
		transformed with the function Bimage::fft.
	The Friedel relationships in transforms from real space images are
	exploited to transform two images simultaneously and then extract
	the individual transforms from the complex data set.
	This function extracts the individual transforms and calculates the 
	complex conjugate product used in cross-correlation.
	An implicit mask is assumed: Only non-zero values in both transforms are considered.
	The result is scaled by the total power of the two transforms within
	the resolution limits, yielding the correlation coefficient when
	the product is backtransformed into the cross-correlation map.

**/
int 		Bimage::combined_complex_product_implicit_mask(double hires, double lores)
{
	if ( !d.uc ) return -1;
	
	if ( compound_type() != TComplex ) return -1;
	
	change_type(Float);
	
	int 			use;
	long	 		i, j, nn, xx, yy, zz, xh(x/2 + 1);
	long	 		ix, iy, iz;
	double			sx2(0), sy2(0), sz2(0), s2, tp1, tp2;
	
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	double			s2hi = (hires > 0)? 1/(hires*hires): 0;
	double			s2lo = (lores > 0)? 1/(lores*lores): 0;

	Vector3<double>	iscale(1.0/real_size());
	Vector3<long>	h(size()/2);
	
	if ( verbose & VERB_FULL ) {
		cout << "Calculating the complex conjugate product:" << endl;
		cout << "Resolution range:               ";
		if ( hires > 0 ) cout << 1/sqrt(s2hi) << " - ";
		else cout << "0 - ";
		if ( lores > 0 ) cout << 1/sqrt(s2lo) << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Implicit mask" << endl;
	}
		
	double			scale(1), sum1, sum2;
	Complex<double>	temp1, temp2, cv;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::combined_complex_product: F0 = " << complex(0).real() << " " << complex(0).imag() << endl;
	
	int				err(0);
	double			tvs = getwalltime(), tvf;
		
    for ( nn=0; nn<n; nn++ ) {
		sum1 = sum2 = 0;
		for ( zz=0; zz<z; zz++ ) {
			iz = (zz>0)? z - zz: 0;
			if ( s2hi > 0 ) {
				sz2 = zz;
				if ( zz > h[2] ) sz2 -= (double)z;
				sz2 *= iscale[2];
				sz2 *= sz2;
			}
			for ( yy=0; yy<y; yy++ ) {
				iy = (yy>0)? y - yy: 0;
				if ( s2hi > 0 ) {
					sy2 = yy;
					if ( yy > h[1] ) sy2 -= (double)y;
					sy2 *= iscale[1];
					sy2 *= sy2;
				}
				for ( xx=0; xx<xh; xx++ ) {
					ix = (xx>0)? x - xx: 0;
					use = 1;
					if ( xx == ix ) {
						if ( yy > h[1] ) use = 0;
						if ( yy == iy && zz > h[2] ) use = 0;
					}
					i = index(xx, yy, zz, nn);
					j = index(ix, iy, iz, nn);
					if ( use ) {
						if ( s2hi > 0 ) {
							sx2 = xx;
							if ( xx > h[0] ) sx2 -= (double)x;
							sx2 *= iscale[0];
							sx2 *= sx2;
							s2 = sx2 + sy2 + sz2;
							if ( s2 < s2lo || s2 > s2hi ) use = 0;
						}
						if ( use ) {
							temp1 = complex(i).unpack_first(complex(j));
							tp1 = temp1.power();
							if ( tp1 < SMALLFLOAT ) {
								use = 0;
							} else {
								temp2 = complex(i).unpack_second(complex(j));
								tp2 = temp2.power();
								if ( tp2 < SMALLFLOAT ) use = 0;
							}
						}
						if ( use ) {
							set(i, temp1 * temp2.conj());
							set(j, complex(i).conj());
							sum1 += temp1.power();
							sum2 += temp2.power();
						} else {
							set(i, cv);
							set(j, cv);
						}
					}
				}
			}
		}
		scale = sum1*sum2;
		if ( scale > 0 ) {
			scale = 0.5/sqrt(scale);
			multiply(nn, scale);
		} else {
			cerr << "Error in Bimage::combined_complex_product: scaling failed!" << endl;
			err--;
		}
	}

	if ( verbose & VERB_TIME ) {
		tvf = getwalltime();
		cout << "CC(BP) time: " << tvf - tvs << endl;
	}
	
	if ( err ) {
		error_show("Error in Bimage::combined_complex_product", __FILE__, __LINE__);
		cerr << "Number of scaling errors: " << err << " (" << n << ")" << endl;
	}
	
	return err;
}

/**
@brief 	Merges the amplitudes from one map with the phases of another.
@param 	*pamp		amplitude image (simple or complex).
@return double		RMSD of amplitudes.

	The amplitude image can be a floating point image or a complex image.
	The phase image must be complex and its amplitudes are replaced by
	the values from the amplitude image.
	No statistics are calculated.

**/
double		Bimage::merge_amplitudes_and_phases(Bimage* pamp)
{
	if ( !d.uc ) return -1;
	if ( !pamp ) return -1;

	if ( compound_type() != TComplex ) return -1;
	
	change_type(Float);
	
	pamp->change_type(Float);
	
	if ( verbose & VERB_FULL )
		cout << "Merging amplitudes and phases" << endl << endl;
	
	long				i, ds(x*y*z*n);
	double			a, ar, d, R(0);
	Complex<double>	cv;
	
	for ( i=0; i<ds; i++ ) {
		cv = complex(i);
		a = cv.amp();
		if ( pamp->compound_type() == TSimple ) ar = (*pamp)[i];
		else ar = (pamp->complex(i)).amp();
		d = a - ar;
		R += d*d;
		cv.amp(ar);
		set(i, cv);
	}
	
	R = sqrt(R/ds);
		
	return R;
}

/**
@brief 	Keeps selected phases and replaces amplitudes and other phases from a reference transform.
@param 	*pref		reference Fourier transform.
@param 	res_hi		high resolution limit.
@param 	res_lo		low resolution limit.
@return double		RMSD of amplitudes.

	The input transform is considered to be modified in some way (such as solvent flattening).
	Only the phases in the specified resolution shell are kept, while all
	the other phases and amplitudes within the high resolution limit is
	retrieved from the reference transform.
	No statistics are calculated.

**/
double		Bimage::merge_amplitudes_and_phases(Bimage* pref, double res_hi, double res_lo)
{
	if ( !d.uc ) return -1;
	if ( !pref ) return -1;
	
	if ( pref->compound_type() != TComplex ) {
		cerr << "Error: The reference map must be a Fourier transform!" << endl << endl;
		return -1;
	}
	
	change_type(Float);
	
	pref->change_type(Float);
	
	check_resolution(res_hi);
	if ( res_lo < res_hi ) res_lo = res_hi;
	
	if ( verbose & VERB_FULL )
		cout << "Merging amplitudes and phases in the range " << res_hi << " - " << res_lo << " Å" << endl << endl;
	
	double			minrad = real_size()[0]/res_lo;
	double			maxrad = real_size()[0]/res_hi;
	if ( maxrad < minrad + 1 ) maxrad = minrad + 1;
	double			minrad2 = minrad*minrad, maxrad2 = maxrad*maxrad;
	
	long	 			i, nn, xx, yy, zz, na(0);
	double 			x2, y2, z2, r2;
	Vector3<long>	h((x - 1)/2, (y - 1)/2, (z - 1)/2);
		
	double			a, ar, d, R(0);
	Complex<double> cv;
	
    for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			z2 = zz;
			if ( zz > h[2] ) z2 -= (long)z;
			z2 *= z2;
			for ( yy=0; yy<y; yy++ ) {
				y2 = yy;
				if ( yy > h[1] ) y2 -= (long)y;
				y2 *= y2;
				for ( xx=0; xx<x; xx++, i++ ) {
					x2 = xx;
					if ( xx > h[0] ) x2 -= (long)x;
					x2 *= x2;
					r2 = x2 + y2 + z2;
					cv = (*this)[i];
					if ( r2 > maxrad2 ) {			// Low-pass filter
						cv = 0;
					} else if ( r2 < minrad2 ) {	// Keep original low-frequency data
						cv = pref->complex(i);
					} else {						// Reset to the original amplitudes
						a = cv.amp();	
						ar = (pref->complex(i)).amp();
						d = a - ar;
						R += d*d;
						na++;
						cv.amp(ar);
					}
					set(i, cv);
				}
			}
		}
	}

	if ( na ) R = sqrt(R/na);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::merge_amplitudes_and_phases: minrad=" << minrad 
			<< " maxrad=" << maxrad << " n=" << n << " R" << R << endl;
		
	return R;
}

/**
@brief 	Generates a power spectrum with phases colored according to a color wheel.
@param 	scale 		amplitude scaling.
@return Bimage*		color power spectrum.

	A polar data type image (such as a Fourier transform) is converted
	to indicate the phases as colors. The primary colors are located
	0 degrees (red), 120 degrees (green) and -120 degrees (blue). All
	three color values are down-weighted based on the amplitude. The
	weighting is calculated as:
		weight = amplitude/(average+2*standard_deviation)
	The origin specified in the image is used to shift the phases.
	The scale is multiplied with the amplitude and cut off at one
	to give the user the ability to enhance the image.
	Default:
		scale = 1/(average_amplitude + standard_deviation_of_amplitude) 

**/
Bimage*		Bimage::intensities_phase_colored(double scale)
{
	if ( compound_type() != TComplex ) {
		cerr << "Error: The image should be complex!" << endl;
		return NULL;
	}
	
	change_type(Float);
	
	statistics();
	
    long				i, ds(x*y*z*n);
	double				amp_ratio;
	
	if ( scale <= 0 ) scale = 1/(avg + std);
	else scale /= avg + std;
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Generating a phase-colored power spectrum" << endl;
	    cout << "Scale:                          " << scale << endl << endl;
	} else if ( verbose & VERB_LABEL )
	    cout << "Generating a phase-colored power spectrum" << endl;
	
	scale *= 255;
	
	// Shift the phases
	phase_shift_to_origin();

	Bimage*     		ps = new Bimage(UCharacter, TRGB, size(), n);
	
	RGB<unsigned char>	rgb;

	for ( i=0; i<ds; i++ ) {
		amp_ratio = (complex(i)).amp()*scale;			// Weigh with amplitudes
		if ( amp_ratio > 0 ) {
			if ( amp_ratio > 255 ) amp_ratio = 255;
			rgb.phase((complex(i)).phi(), amp_ratio);
//			rgb.phase((complex(i)).phi_lut(), amp_ratio);
			ps->set(i, rgb);
		}
	}

	ps->center_wrap();
	
	return ps;
}


