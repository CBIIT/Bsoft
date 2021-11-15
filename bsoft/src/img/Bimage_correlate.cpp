/**
@file	Bimage_correlate.cpp
@brief	Cross correlation functions
@author Bernard Heymann
@date	Created: 19980805
@date	Modified: 20200218

		Implemented using the FFTW library
**/

#include "Bimage.h"
#include "rwimg.h"
#include "matrix_linear.h"
#include "math_util.h"
#include "qsort_functions.h"
#include "timer.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates correlation coefficient between two images.
@param 	*p			image to correlate with.
@return double 		correlation coefficient, -1 if not run.

	The correlation between two images is calculated and normalized as:
		           sum((image1 - avg1)*(image2 - avg2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	.
	Both images are converted to floating point.
	Only the first image is used.

**/
double		Bimage::correlate(Bimage* p)
{
	if ( c > 1 || p->channels() > 1 ) {
		cerr << "Error: Only single channel images are supported for correlation!" << endl;
		return -1.0;
	}
	
	if ( n > 1 || p->images() > 1 ) {
		cerr << "Error: Only single images are supported for correlation!" << endl;
		return -1.0;
	}
	
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::correlate", __FILE__, __LINE__);
		return -1;
	}
	
    long   		i;
	double			v1, v2, cc(0), sx(0), sy(0), sx2(0), sy2(0), sxy(0);
	
	for ( i=0; i<datasize; i++ ) {
		v1 = (*this)[i];
		v2 = (*p)[i];
		sx += v1;
		sy += v2;
		sx2 += v1*v1;
		sy2 += v2*v2;
		sxy += v1*v2;
	}
	
	double			varx = sx2 - sx*sx/datasize;
	double			vary = sy2 - sy*sy/datasize;
	
	if ( varx > 0 && vary > 0 )
		cc = (sxy - sx*sy/datasize)/sqrt(varx*vary);
	
	return cc;
}

/**
@brief 	Calculates a correlation coefficient between two images.
@param 	*p			second image.
@param 	rmin		minimum radius (pixel units).
@param 	rmax		maximum radius (pixel units).
@param	*pmask		mask (can be NULL).
@param	flag		flag to modify second image.
@return double 		correlation coefficient, -1 if error.

	The correlation between two images is calculated and normalized as:
		         sum((image1 - avg1)*(image2 - avg2))
		CC = -------------------------------------------------
		     sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)

**/
double		Bimage::correlate(Bimage* p, double rmin, double rmax, Bimage* pmask, int flag)
{
	if ( n > 1 || p->images() > 1 ) {
		cerr << "Error: Only single images are supported for correlation!" << endl;
		return -1.0;
	}
	
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::correlate", __FILE__, __LINE__);
		return -1;
	}
	
	if ( std <= 0 || p->std <= 0 ) return -1;
	
	if ( rmin > x/2 ) rmin = 0;
	if ( rmax <= rmin ) rmax = size().max();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a correlation map:" << endl;
		cout << "Origin:                         " << image->origin() << endl;
		cout << "Radii:                          " << rmin << " - " << rmax << endl;
		if ( pmask )
		cout << "Mask:                           " << pmask->file_name() << endl;
	}
	
	int				use_this;
    long			i, xx, yy, zz, ns(0);
	double			v1, v2, dx, dy, dz, d2, rmin2(rmin*rmin), rmax2(rmax*rmax);
    double  		avg1(0), avg2(0), std1(0), std2(0), c12(0), denom, cc(0);
	
	for ( i=zz=0; zz<z; zz++ ) {
		dz = (double)zz - image->origin()[2];
		dz *= dz;
		for ( yy=0; yy<y; yy++ ) {
			dy = (double)yy - image->origin()[1];
			dy *= dy;
			for ( xx=0; xx<x; xx++, i++ ) {
				dx = (double)xx - image->origin()[0];
				dx *= dx;
				d2 = dx + dy + dz;
				use_this = 1;
				if ( pmask && (*pmask)[i] < 1 ) use_this = 0;
				if ( d2 < rmin2 || d2 > rmax2 ) use_this = 0;
				if ( use_this ) {
					v1 = (*this)[i];
					v2 = (*p)[i];
					avg1 += v1;
					avg2 += v2;
					std1 += v1*v1;
					std2 += v2*v2;
					c12 += v1*v2;
					ns++;
				} else if ( flag ) {
					p->set(i, 0);
				}
			}
		}
	}
	
	if ( ns ) {
		avg1 /= ns;
		avg2 /= ns;
		std1 = std1/ns - avg1*avg1;
		std2 = std2/ns - avg2*avg2;
		if ( std1 > 0 ) std1 = sqrt(std1);
		else std1 = 0;
		if ( std2 > 0 ) std2 = sqrt(std2);
		else std2 = 0;
		c12 /= ns;
	}
	
	denom = std1*std2;
	if ( denom <= 0 ) {
//		cerr << "Error in Bimage::correlate: One or both standard deviations are zero!" << endl;
		return 0;
	}
	
	denom = 1/denom;

	cc = (c12 - avg1*avg2) * denom;
	
	if (flag ) {
		for ( i=0; i<datasize; i++ ) if ( (*p)[i] )
			p->set(i, ((*this)[i] - avg1)*((*p)[i] - avg2)*denom);
		p->statistics();
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Averages:                       " << avg1 << " " << avg2 << endl;
		cout << "Standard deviations:            " << std1 << " " << std2 << endl;
		cout << "Correlation:                    " << cc << " (" << ns << ")" << endl << endl;
	}
	
	return cc;
}

/**
@brief 	Rotates a copy of an image and correlates it with the original.
@param 	axis		rotation axis.
@param 	angle		rotation angle.
@return double 		correlation coefficient, -1 if error.

	This mainly used to determine the symmetry of a map.

**/
double		Bimage::rotate_correlate(Vector3<double> axis, double angle)
{
	Bimage*			p = rotate(size(), axis, angle);
	
	double			CC = correlate(p);
	
	delete p;
	
	return CC;
}

/**
@brief 	Calculates an R factor between two images.
@param 	*p			second image.
@return double 		R factor, -1 if not run.

	The difference between two images is calculated and normalized as:
		                      sum(image1 - image2)^2
		R = sqrt(-------------------------------------------------)
		         sqrt(sum(image1 - avg1)^2 * sum(image2 - avg2)^2)
	Both images are converted to floating point.

**/
double		Bimage::R_factor(Bimage* p)
{
	if ( n > 1 || p->images() > 1 ) {
		cerr << "Error: Only single images are supported for R factor calculation!" << endl;
		return -1.0;
	}
	
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::R_factor", __FILE__, __LINE__);
		return -1;
	}
	
    long			i;
	double			R(0);
	double			sx(0), sy(0), sx2(0), sy2(0), sxy(0);
	
	for ( i=0; i<datasize; i++ ) {
		sx += (*this)[i];
		sy += (*p)[i];
		sx2 += (*this)[i]*(*this)[i];
		sy2 += (*p)[i]*(*p)[i];
		sxy += (*this)[i]*(*p)[i];
	}
	
	R = sqrt((sx2 - 2*sxy + sy2)/sqrt((sx2 - sx*sx/datasize)*(sy2 - sy*sy/datasize)));
	
	return R;
}


/**
@brief 	Calculates an autocorrelation map by Fast Fourier transformation.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return int 		0.

	FFTW library (www.fftw.org).
	A multi-image 1D, 2D and 3D data set is transformed forward, the 
	transform multiplied by its complex conjugate, followed by backward 
	transformation and rescaling by 1/(N*N). Data beyond the resolution 
	set in the image structure are zeroed. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	The resultant image is Float.

**/
int 		Bimage::auto_correlate(double hires, double lores)
{
	if ( hires > lores ) swap(hires, lores);
	
	long	   		nn, i, xx, yy, zz;
	long			ix, iy, iz;
	double			s2, sx2, sy2, sz2, smin2(0), smax2(0), sum, scale;
	Vector3<double>	u(image->sampling());
	Vector3<double>	realsize(real_size());
	
	if ( lores > 0 ) smin2 = 1/(lores*lores);
	
	if ( hires > 0 ) {
		smax2 = 1/(hires*hires);
	} else {
		if ( x > 1 ) smax2 += 1.0/(4.0*u[0]*u[0]);
		if ( y > 1 ) smax2 += 1.0/(4.0*u[1]*u[1]);
		if ( z > 1 ) smax2 += 1.0/(4.0*u[2]*u[2]);
		hires = 1/sqrt(smax2);
	}

	if ( verbose & VERB_FULL ) {
		cout << "Auto-correlation:" << endl;
		if ( lores > 0 )
			cout << "Resolution range:               " << hires << " - " << lores << " A" << endl;
		else
			cout << "No resolution limits" << endl;
		cout << endl;
	}
		
	fft();
	
	Complex<double>		cv;
//	set(0, cv);
	
    for ( i=nn=0; nn<n; nn++ ) {
		sum = 0;
		for ( zz=0; zz<z; zz++ ) {
			iz = zz;
			if ( zz > (z - 1)/2 ) iz -= z;
			sz2 = iz/realsize[2];
			sz2 *= sz2;
			for ( yy=0; yy<y; yy++ ) {
				iy = yy;
				if ( yy > (y - 1)/2 ) iy -= y;
				sy2 = iy/realsize[1];
				sy2 *= sy2;
				for ( xx=0; xx<x; xx++, i++ ) {
					ix = xx;
					if ( xx > (x - 1)/2 ) ix -= x;
					sx2 = ix/realsize[0];
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					if ( s2 >= smin2 && s2 <= smax2 ) {
						cv = (complex(i)).power();
						sum += cv.real();
					} else {
						cv = 0;
					}
					set(i, cv);
				}
			}
		}
		if ( sum ) {
			scale = sqrt(x*y*z)/sum;
//			cout << "scale=" << scale << endl;
			multiply(nn, scale);
		}
		image[nn].origin(0.0,0.0,0.0);
	}
	
//	cout << "backtransforming" << endl;

	fft_back();
	
	return 0;
}

/**
@brief 	Calculates a cross-correlation map by Fast Fourier transformation.
@param 	*p			second image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param	planf		forward transform plan.
@param	planb		backward transform plan.
@return Bimage* 	cross-correlation image.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is unpacked before the first transform is multiplied with the complex 
	conjugate of the second transform. This is then back-transformed to 
	obtain the cross-correlation map in real space.
	Zero-valued data in the transforms are implicitly excluded.
	The low resolution limit can be 0, in which case no limits are applied.
	The resultant cross-correlation image data type is floating point.

**/
Bimage* 	Bimage::cross_correlate(Bimage* p, double hires, double lores,
				fft_plan planf, fft_plan planb)
{
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Cross-correlation:" << endl;
		if ( lores > 0 || hires > 0 ) {
			cout << "Resolution range:               " << hires << " - ";
			if ( lores > 0 ) cout << lores << " A" << endl;
			else cout << "inf A" << endl;
		} else
			cout << "No resolution limits" << endl;
		cout << endl;
	}
	
	pc->fft(planf, 0);

//	pc->set(0, Complex<double>(0,0));

	pc->combined_complex_product_implicit_mask(hires, lores);
	
	pc->fft(planb, 0);
	
	pc->complex_to_real();	

	long		nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(p->image[nn].origin() - image[nn].origin());
	
	return pc;
}

/**
@brief 	Calculates a cross-correlation map by Fast Fourier transformation.
@param 	*p			second image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		binary mask (only 0 and 1), NULL if not desired.
@return Bimage* 	cross-correlation image.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is unpacked before the first transform is multiplied with the complex 
	conjugate of the second transform. This is then back-transformed to 
	obtain the cross-correlation map in real space.
	The low resolution limit can be 0, in which case no limits are applied.
	The resultant cross-correlation image data type is floating point.

**/
Bimage* 	Bimage::cross_correlate(Bimage* p, double hires, double lores, Bimage* pmask)
{
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Cross-correlation:" << endl;
		if ( lores > 0 || hires > 0 ) {
			cout << "Resolution range:               " << hires << " - ";
			if ( lores > 0 ) cout << lores << " A" << endl;
			else cout << "inf A" << endl;
		} else
			cout << "No resolution limits" << endl;
		if ( pmask )
			cout << "With a mask:                    " << pmask->file_name() << endl;
		cout << endl;
	}
	
	pc->fft(FFTW_FORWARD, 0);
	
//	pc->set(0, Complex<double>(0,0));
	
	pc->combined_complex_product(hires, lores, pmask);
	
	pc->fft(FFTW_BACKWARD, 0);
	
	pc->complex_to_real();	

	long		nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(p->image[nn].origin() - image[nn].origin());
	
//	write_img("pcc.mrc", pc, 0);
	
	return pc;
}

/**
@brief 	Calculates a coefficient from a Fourier correlation transform given a shift.
@param 	shift		real space shift.
@return double	 	correlation coefficient.

	The input image is a cross-correlation transform.
	A brute force backtransform integration is done for the given shift to 
	calculate the corresponding correlation coefficient.

**/
double			Bimage::correlation_coefficient(Vector3<double> shift)
{
	long			i;
	Vector3<long>	k;
	double			phi, cc(0);
	Complex<double>	cv;
	
	shift = shift*TWOPI;
	shift /= size();
	
	for ( i=k[2]=0; k[2]<z; ++k[2] ) {
		for ( k[1]=0; k[1]<y; ++k[1] ) {
			for ( k[0]=0; k[0]<x; ++k[0], ++i ) {
				phi = (shift * k).sum();
				cv = complex(i);
				cc += cv.real()*cos(phi) - cv.imag()*sin(phi);
//				cout << k << tab << phi << endl;
			}
		}
	}
	
	return cc;			
}

/**
@brief 	Finds the shift by brute force backtransformation for selected shifts.
@param 	shift_limit	maximum	real space shift.
@return Vector3<double>	best shift.

	The input image is a cross-correlation transform.
	A brute force integration is done for the given shift to calculate
	the corresponding correlation coefficient.
	The correlation coefficient is return in the FOM field of the sub-image.

**/
Vector3<double>	Bimage::find_shift_in_transform(double shift_limit)
{
	long			i, nn(0);
	double			ds(1), cc(0), ccmax(0);
	Vector3<double>	sh, bestshift;
	Vector3<long>	smn, smx(shift_limit,shift_limit,shift_limit);
	smx = smx.min(size()/2);
	smn = -smx;
	Vector3<long>	size = smx - smn + 1;
	if ( verbose & VERB_FULL )
		cout << size << tab << smn << tab << smx << endl;
	
	Bimage*			pc = new Bimage(Double, TSimple, size, 1);
	
	for ( i=0, sh[2]=smn[2]; sh[2]<=smx[2]; sh[2]+=ds ) {
		for ( sh[1]=smn[1]; sh[1]<=smx[1]; sh[1]+=ds ) {
			for ( sh[0]=smn[0]; sh[0]<=smx[0]; sh[0]+=ds, ++i ) {
				cc = correlation_coefficient(sh);
				if ( ccmax < cc ) {
					ccmax = cc;
					bestshift = sh;
					if ( verbose & VERB_FULL )
						cout << setprecision(4) << ds << tab << sh << tab << cc << endl;
				}
				pc->set(i, cc);
			}
		}
	}
		
	pc->origin(bestshift+size/2);
	bestshift += pc->fit_peak();

	image[nn].origin(bestshift);
	image[nn].FOM(ccmax);
	
//	write_img("cck.mrc", pc, 0);
	
	delete pc;
	
	return bestshift;
}

/**
@brief 	Calculates a cross-correlation map by Fast Fourier transformation.
@param 	*p			second image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	shift_limit	maximum	real space shift.
@return Bimage* 	cross-correlation image.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is unpacked before the first transform is multiplied with the complex 
	conjugate of the second transform. This is then back-transformed to 
	obtain the cross-correlation map in real space.
	The low resolution limit can be 0, in which case no limits are applied.
	The resultant cross-correlation image data type is floating point.

**/
Bimage* 	Bimage::cross_correlate_fspace(Bimage* p, double hires, double lores, double shift_limit)
{
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Cross-correlation:" << endl;
		if ( lores > 0 || hires > 0 ) {
			cout << "Resolution range:               " << hires << " - ";
			if ( lores > 0 ) cout << lores << " A" << endl;
			else cout << "inf A" << endl;
		} else
			cout << "No resolution limits" << endl;
//		if ( pmask )
//			cout << "With a mask:                    " << pmask->file_name() << endl;
		cout << endl;
	}
	
	pc->fft(FFTW_FORWARD, 0);
	
//	pc->set(0, Complex<double>(0,0));
	
	pc->combined_complex_product(hires, lores, NULL);
	
	Vector3<double>	shift = pc->find_shift_in_transform(shift_limit);

	long		nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(shift);
	
//	write_img("pcc.mrc", pc, 0);
	
	return pc;
}


Bimage* 	Bimage::cross_correlate(Bimage* p, double hires, double lores,
				Bimage* pmask, fft_plan planf, fft_plan planb)
{
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Cross-correlation:" << endl;
		if ( lores > 0 || hires > 0 ) {
			cout << "Resolution range:               " << hires << " - ";
			if ( lores > 0 ) cout << lores << " A" << endl;
			else cout << "inf A" << endl;
		} else
			cout << "No resolution limits" << endl;
		if ( pmask )
			cout << "With a mask:                    " << pmask->file_name() << endl;
		cout << endl;
	}
	
	pc->fft(planf, 0);

//	pc->set(0, Complex<double>(0,0));

	pc->combined_complex_product(hires, lores, pmask);
	
	pc->fft(planb, 0);
	
	pc->complex_to_real();	
	
	if ( pc->image->maximum() > 1 ) {
		cerr << "Error in cross_correlate: maximum > 1: " << image->maximum() << endl;
		bexit(-1);
	}

	long		nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(p->image[nn].origin() - image[nn].origin());
	
///	write_img("pccp.mrc", pc, 0);
	
	return pc;
}

/**
@brief 	Calculates a cross-correlation map by Fast Fourier transformation.
@param 	*p			second image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@return Bimage* 	cross-correlation image.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is bandpass filtered to the given resolution limits and unpacked.
	The two transforms are normalized, and multiplied in two ways:
		prod1 = data1 * data2.conjugate
		prod2 = data1 * data2
	The second product represents a cross-correlation with a 180 degree 
	rotation of the second image.
	Both are back-transformed to obtain the two cross-correlation map in real space.
	The one with the highest maximum is selected to return.
	The resultant cross-correlation image data type is floating point.
	The angle in the cross-correlation image is set to zero if the first
	map is selected, otherwise it is set to PI.

**/
Bimage* 	Bimage::cross_correlate_two_way(Bimage* p, double hires,
				double lores, fft_plan planf, fft_plan planb)
{
	if ( hires > lores ) swap(hires, lores);
	
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	check_resolution(hires);
	
	if ( verbose & VERB_FULL ) {
		cout << "Two-way cross-correlation:" << endl;
		cout << "Resolution range:               " << hires << " - " << lores << " A" << endl;
		cout << endl;
	}
	
	pc->fft(planf, 0);

	Complex<double>		cv;
	pc->set(0, cv);

	pc->complex_bandpass(hires, lores);
	
	Bimage* 	pc2 = pc->unpack_combined_transform();
	
	pc->complex_normalize();
	pc2->complex_normalize();

	Bimage* 	pc1 = pc->copy();

	pc->complex_conjugate_product(pc2);
	pc1->complex_product(pc2);
	
//	write_img("pc.map", pc, 0);
//	bexit(0);

	delete pc2;
	
	pc->fft(planb, 0);
	pc1->fft(planb, 0);
	
	pc->complex_to_real();	
	pc1->complex_to_real();
	
//	cout << pc->maximum() << tab << pc1->maximum() << endl;

	if ( pc1->maximum() > pc->maximum() ) {
		delete pc;
		pc = pc1;
		pc->image->view_angle(M_PI);
	} else {
		delete pc1;
		pc->image->view_angle(0);
	}
	
	pc->fourier_type(NoTransform);

	long		nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(p->image[nn].origin() - image[nn].origin());
	
	return pc;
}

/**
@brief 	Calculates a masked cross-correlation map by Fast Fourier transformation.
@param 	*p			second image.
@param 	*pmask		binary mask (only 0 and 1), NULL if not desired.
@return Bimage* 	cross-correlation image.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D real space data sets are 
	packed into a complex data set and transformed forward. The transform
	is unpacked and masked with mask image before the first transform is
	multiplied with the complex conjugate of the second transform. This is 
	then back-transformed to obtain the cross-correlation map in real space.
	The mask must be composed of 0 and 1, and is converted to floating point.
	The mask can be omitted (NULL).
	The resultant cross-correlation image data type is floating point.

**/
Bimage* 	Bimage::cross_correlate_validate(Bimage* p, Bimage* pmask)
{
	Bimage* 	pc = pack_two_in_complex(p);
	if ( !pc ) return NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Cross-correlation with validation:" << endl;
		if ( pmask )
			cout << "With a mask:                    " << pmask->file_name() << endl;
		cout << endl;
	}
	
	pc->fft(FFTW_FORWARD, 0);

	Complex<double>		cv;
	pc->set(0, cv);

	long	nn;
	for ( nn=0; nn<pc->images(); nn++ )
		pc->image[nn].origin(p->image[nn].origin() - image[nn].origin());
		
	if ( pmask ) pc->complex_apply_dual_mask(pmask);
	
	pc->combined_complex_product();
	pc->next->combined_complex_product();
	
	pc->fft(FFTW_BACKWARD, 0);
	pc->next->fft(FFTW_BACKWARD, 0);
	
	pc->complex_to_real();	
	pc->next->complex_to_real();	
	
	return pc;
}

/*
@brief 	Searches a 2D/3D density map for a template using a specific view.
@param 	*p				the image.
@param 	*pref			the template to be searched for.
@param 	view			view.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	search_radius	radius for shift search.
@param 	*pmask			mask for cross-correlation (ignored if NULL).
@param 	&cc				the best correlation coefficient
@param 	planf			FFT forward plan.
@param 	planb			FFT backward plan.
@return Vector3<double>	best shift to impose on the image.

	The template is rotated to the view and cross-correlated to find
	a set of high-scoring fits.
	The views must be calculated externally to allow for custom sets.
**/
Vector3<double>	Bimage::rotate_cross_correlate(Bimage* pref, View view,
				double hires, double lores, double search_radius, Bimage* pmask, 
				double& cc, fft_plan planf, fft_plan planb)
{
	Vector3<double>	shift(-1e37,-1e37,-1e37);

	if ( !check_if_same_size(pref) ) {
		error_show("Bimage::rotate_cross_correlate", __FILE__, __LINE__);
		return shift;
	}
	
	Matrix3			mat = view.matrix();
	mat = mat.transpose();
	
	Bimage*			prot = pref->rotate(pref->size(), mat);
	
	Bimage*			pcc = cross_correlate(prot, hires, lores, pmask, planf, planb);

	pcc->find_peak(search_radius, 0);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::rotate_cross_correlate: found origin = " << pcc->image->origin() << endl;
		cout << "DEBUG Bimage::rotate_cross_correlate: cc = " << pcc->image->FOM() << endl;
	}

	pcc->refine_peak();
	
//	if ( verbose & VERB_DEBUG )
//		cout << "DEBUG Bimage::rotate_cross_correlate: refined origin = " << pcc->image->origin() << endl;
	
	shift = -pcc->image->origin();
	cc = pcc->image->FOM();

	delete prot;
	delete pcc;
	
	return shift;
}

/*
@brief 	Rotates a 2D reference image and cross-correlates with an image.
@param 	*pref		reference 2D image.
@param 	angle		rotation angle.
@param 	res_hi		high resolution limit.
@param 	res_lo		low resolution limit.
@param 	shift_limit	maximum shift from nominal origin of box.
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@return double 		correlation coefficient.

	The reference image is first rotated by the given angle.
	Two cross-correlations are done to account for a 180° ambiguity in
	the angle when it is derived from a polar power spectrum.
	The shift is then determined from the best cross-correlation map.
	The image origin is written into the image structure.
	The correlation coefficient return is the cross-correlation peak.
	The input images must be equal-sized square 2D images.
**/
double		Bimage::rotate_cross_correlate_two_way(Bimage* pref, 
				double angle, double res_hi, double res_lo, double shift_limit,
				fft_plan planf, fft_plan planb)
{
//	cout << pref->image->origin() << endl;
	
	Bimage*			prot = pref->rotate(pref->size(), fabs(angle));
	
	Bimage*			pcc = cross_correlate_two_way(prot, res_hi, res_lo, planf, planb);
	
	delete prot;
	
	pcc->origin(0.0, 0.0, 0.0);

	pcc->find_peak(shift_limit, 0);
	
	pcc->refine_peak();
	
	if ( angle > 0 ) image->view_angle(angle + pcc->image->view_angle());
	else image->view_angle(angle - pcc->image->view_angle());
	
//	cout << pcc->image->origin() << tab << pcc->image->FOM() << endl;
	
	if ( fabs(pcc->image->view_angle()) > M_PI_2 ) {
//		image->origin(pref->image->origin() - pcc->image->origin());
		origin(pref->image->origin() - pcc->image->origin());
	} else {
//		image->origin(pref->image->origin() + pcc->image->origin());
		origin(pref->image->origin() + pcc->image->origin());
	}
	
	image->FOM(pcc->image->FOM());
	
	delete pcc;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::rotate_cross_correlate_two_way: fom=" << image->FOM() << endl;
		cout << "DEBUG Bimage::rotate_cross_correlate_two_way: ori=" << image->origin() << endl;
	}
	
	return image->FOM();
}


/**
@brief 	Calculates a cross-correlation map to find the shift for the pair of images.
@param 	*pref		reference image.
@param 	*pmask		binary mask (only 0 and 1).
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	radius		search radius (if < 1, default 1e30).
@param 	sigma		attenuation around radius.
@param 	refine_flag set to refine shift to subpixel resolution.
@return int 		error code.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D data sets are transformed 
	forward, the first transform multiplied by the complex conjugate of
	the second transform, followed by backward transformation and 
	rescaling by 1/(N*N). Data beyond the resolution set in the first 
	image structure are not used. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	A shift vector for each pair of images is calculated to
	determine the cross-correlation peak to sub-pixel resolution.
	Note: The first image is the reference and the shift returned is to
		transform the second to fit the first.

**/
int			Bimage::find_shift(Bimage* pref, Bimage* pmask, double hires,
				double lores, double radius, double sigma, int refine_flag)
{
	if ( verbose & VERB_FULL ) {
		if ( refine_flag )
			cout << "Finding shift by cross-correlation and polynomial fitting" << endl << endl;
		else
			cout << "Finding shift by cross-correlation" << endl << endl;
	}
	
	Bimage* 		pc = cross_correlate(pref, hires, lores, pmask);
	if ( !pc ) return -1;
	
	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc, 0);
	
	long			nn;
	pc->find_peak(radius, sigma);

	if ( refine_flag ) pc->refine_peak();
	
	for ( nn=0; nn<pc->images(); nn++ ) {
		image[nn].origin(pref->image[nn].origin() + pc->image[nn].origin());
		image[nn].FOM(pc->image[nn].FOM());
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Image\t   ox\t   oy\t   oz\t  CC" << endl;
		for ( nn=0; nn<n; nn++ )
			cout << nn+1 << tab << image[nn].origin()
				<< tab << image[nn].FOM() << endl;
		cout << endl;
	}
	
	delete pc;
	
	return 0;
}

Vector3<double>	Bimage::find_shift(Bimage* pref, Bimage* pmask, double hires,
				double lores, double radius, double sigma, int refine_flag, double& cc)
{
	if ( verbose & VERB_FULL ) {
		if ( refine_flag )
			cout << "Finding shift by cross-correlation and polynomial fitting" << endl << endl;
		else
			cout << "Finding shift by cross-correlation" << endl << endl;
	}

	Vector3<double>	shift;
	
	Bimage* 		pc = cross_correlate(pref, hires, lores, pmask);
	if ( !pc ) return shift;
	
	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc, 0);
	
	pc->find_peak(radius, sigma);

	if ( refine_flag ) pc->refine_peak();
	
	shift = pc->image->origin();
	
	cc = pc->image->FOM();
	
	delete pc;
	
	return shift;
}

/**
@brief 	Calculates a cross-correlation map to find the shift for the pair of images.
@param 	*pref		reference image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	radius		search radius (if < 1, default 1e30).
@param 	sigma		attenuation around radius.
@param 	refine_flag set to refine shift to subpixel resolution.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return Vector3<double>	shift.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D data sets are transformed 
	forward, the first transform multiplied by the complex conjugate of
	the second transform, followed by backward transformation and 
	rescaling by 1/(N*N). Data beyond the resolution set in the first 
	image structure are not used. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	A shift vector for each pair of images is calculated to
	determine the cross-correlation peak to sub-pixel resolution.
	Note: The first image is the reference and the shift returned is to
		transform the second to fit the first.

	Only the first sub-image shift is calculated.

**/
Vector3<double>	Bimage::find_shift(Bimage* pref, double hires,
				double lores, double radius, double sigma, int refine_flag,
				fft_plan planf, fft_plan planb)
{
	if ( verbose & VERB_FULL ) {
		if ( refine_flag )
			cout << "Finding shift by cross-correlation and polynomial fitting" << endl << endl;
		else
			cout << "Finding shift by cross-correlation" << endl << endl;
	}
	
	Vector3<double>	shift;
	
	Bimage* 		pc = cross_correlate(pref, hires, lores, planf, planb);
	if ( !pc ) return shift;
	
	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc, 0);
	
	pc->find_peak(radius, sigma);

	if ( refine_flag ) pc->refine_peak();
	
	shift = pc->image->origin();
	
	image->FOM(pc->image->FOM());
	
	delete pc;
	
	return shift;
}

/**
@brief 	Calculates a cross-correlation map to find the shift for the pair of images.
@param 	*pref		reference image.
@param 	*pmask		binary mask (only 0 and 1).
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	radius		search radius (if < 1, default 1e30).
@param 	sigma		attenuation around radius.
@param 	refine_flag set to refine shift to subpixel resolution.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@param 	&cc			correlation coefficient
@return Vector3<double>	shift.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D data sets are transformed 
	forward, the first transform multiplied by the complex conjugate of
	the second transform, followed by backward transformation and 
	rescaling by 1/(N*N). Data beyond the resolution set in the first 
	image structure are not used. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	A shift vector for each pair of images is calculated to
	determine the cross-correlation peak to sub-pixel resolution.
	Note: The first image is the reference and the shift returned is to
		transform the second to fit the first.

	Only the first sub-image shift is calculated.

**/
Vector3<double>	Bimage::find_shift(Bimage* pref, Bimage* pmask, double hires,
				double lores, double radius, double sigma, int refine_flag,
				fft_plan planf, fft_plan planb, double& cc)
{
	if ( verbose & VERB_FULL ) {
		if ( refine_flag )
			cout << "Finding shift by cross-correlation and polynomial fitting" << endl << endl;
		else
			cout << "Finding shift by cross-correlation" << endl << endl;
	}
	
	Vector3<double>	shift;
	
	Bimage* 		pc = cross_correlate(pref, hires, lores, pmask, planf, planb);
	if ( !pc ) return shift;
	
	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc, 0);
	
	pc->find_peak(radius, sigma);

	if ( refine_flag ) pc->refine_peak();
	
	shift = pc->image->origin();
	
	cc = pc->image->FOM();
	
	if ( cc > 1 ) {
		cerr << "Error in find_shift: cc > 1: " << cc << endl;
		bexit(-1);
	}
	
	delete pc;
	
	return shift;
}

/**
@brief 	Calculates a cross-correlation map to find the shift for the pair of images.
@param 	nn			sub-image to align.
@param 	*pref		reference image.
@param 	*pmask		binary mask for cross-correlation (only 0 and 1).
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit.
@param 	shift_limit	search radius (if < 1, default 1e30).
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return double		correlation coefficient.

	The shift is returned as the origin of each sub-image relative to
	the reference origin.
	
**/
Vector3<double>	Bimage::find_shift(long nn, Bimage* pref, Bimage* pmask,
				double hi_res, double lo_res, double shift_limit, 
				fft_plan planf, fft_plan planb)
{
	double			cc(0);
	Bimage*			p1 = extract(nn);
	
	Vector3<double>	shift = p1->find_shift(pref, pmask, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);

	delete p1;
	
	image[nn].origin(pref->image->origin()+shift);
	
	image[nn].FOM(cc);
	
	return shift;
}

/**
@brief 	Finds one or more matches to a template.
@param 	*ptemp		template image.
@param 	*pmask		reciprocal space mask.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	bin			level of image binning.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.
@return Bimage*		cross-correlation map.

	A template is cross-correlated with the input image including
	bandpass filtering to target the size of the template particle.

**/
Bimage*		Bimage::find_template(Bimage* ptemp, Bimage* pmask,
				double hires, double lores, int bin, fft_plan planf, fft_plan planb)
{
	Bimage*			pc = bin_copy(bin);

	Bimage*			pt = ptemp->bin_copy(bin);
	pt->pad(pc->size(), FILL_USER, 0);

	Bimage*			pm = NULL;
	if ( pmask ) pm = pmask->pad_copy(pc->size(), FILL_USER, 0);
	
	Vector3<double>	shift(pt->image->origin());
	
	if ( verbose & VERB_PROCESS ) {
//		cout << "Number of templates:             " << pt->images() << endl;
		cout << "Template origin:                 " << pt->image->origin() << endl;
		cout << "Template avg & std:              " << pt->average() << 
			tab << pt->standard_deviation() << endl;
	}
	
	Bimage*			pcc = pc->cross_correlate(pt, hires, lores, pm, planf, planb);

	delete pc;
	delete pt;
	delete pm;

	pcc->shift_wrap(shift);
	
	pcc->statistics();
		
	return pcc;
}

/**
@brief 	Finds the center of mass of an image by cross-correlation with its inverse.
@param 	*pmask		binary mask (only 0 and 1).
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	radius		search radius (if < 1, default 1e30).
@param 	sigma		attenuation around radius.
@param 	refine_flag set to refine shift to subpixel resolution.
@return int			0.

	FFTW library (www.fftw.org).
	Two equally sized multi-image 1D, 2D and 3D data sets are transformed 
	forward, the first transform multiplied by the complex conjugate of
	the second transform, followed by backward transformation and 
	rescaling by 1/(N*N). Data beyond the resolution set in the first 
	image structure are not used. Therefore the correct setting 
	of units and resolution in the image are required. Defaults for the 
	units are usually 1 Angstrom/voxel and a zero resolution would
	include the whole image (i.e., no resolution limitation).
	A shift vector for each pair of images is calculated to
	determine the cross-correlation peak to sub-pixel resolution.
	Note: The first image is the reference and the shift returned is to
		transform the second to fit the first.

**/
int			Bimage::find_center(Bimage* pmask, double hires,
				double lores, double radius, double sigma, int refine_flag)
{
	if ( verbose & VERB_FULL ) {
		if ( refine_flag )
			cout << "Finding center by cross-correlation and polynomial fitting" << endl << endl;
		else
			cout << "Finding center by cross-correlation" << endl << endl;
	}
	
	Bimage*				pc = copy();
	
	pc->fft(FFTW_FORWARD, 0);
	
	Complex<double>		cv;
	pc->set(0, cv);
	if ( pmask ) pc->complex_apply_mask(pmask);
	
	if ( hires > sampling(0)[0] ) pc->complex_bandpass(hires, lores);

	pc->complex_normalize();

	pc->complex_product(pc);
	
	pc->fft(FFTW_BACKWARD, 0);
	
	pc->complex_to_real();

	if ( verbose & VERB_DEBUG )
		write_img("cc.map", pc, 0);
	
	long				nn;
	
	pc->origin(0.0,0.0,0.0);
	pc->find_peak(radius, sigma);

	if ( refine_flag ) pc->refine_peak();
	
	for ( nn=0; nn<n; nn++ ) {
		image[nn].origin(pc->image[nn].origin()/2 + size()/2);
		image[nn].FOM(pc->image[nn].FOM());
	}
	
	delete pc;
	
	if ( verbose & VERB_FULL ) {
		cout << "Image\t   ox\t   oy\t   oz\t  CC" << endl;
		for ( nn=0; nn<n; nn++ )
			cout << nn+1 << tab << image[nn].origin()
				<< tab << image[nn].FOM() << endl;
		cout << endl;
	}
	
	return 0;
}

/**
@brief 	Rotates and find shift by cross-correlation.
@param 	mat			rotation matrix.
@param 	hires		high resolution limit in angstroms.
@param 	lores		low resolution limit in angstroms.
@param 	radius		search radius.
@param 	sigma		soften search radius cutoff.
@param 	refine_flag	refine shift.
@param 	planf		forward Fourier transform plan.
@param 	planb		backward Fourier transform plan.k
@param 	&cc			correlation coefficient
@return Vector3<double>	origin.

	The point group symmetry operations are applied to an image with an
	orientation defined by the reference symmetry axis (default {0,0,1}). 

**/
Vector3<double>	Bimage::rotate_find_shift(Matrix3 mat,
				double hires, double lores, double radius, double sigma,
				int refine_flag, fft_plan planf, fft_plan planb, double& cc)
{
	if ( radius <= 0 ) radius = x/4;

	Bimage*			p = rotate(size(), mat);

	Vector3<double>	shift = find_shift(p, NULL, hires, lores, radius, sigma, 1, planf, planb, cc);

	Vector3<double>	origin;
	mat = Matrix3(1) - mat;
	mat = mat.singular_value_decomposition();
	origin = mat * shift + image->origin();
	
	delete p;
	
	return origin;
}

/**
@brief 	Finds the peak in an image to the nearest voxel.
@param 	radius		search radius (if < 1, default 1e30).
@param 	sigma		attenuation around radius.
@return int 		error code.

	An image is searched for the global maximum (typically used to find
	the shift vector in a cross-correlation map).
	The peak vectors returned are in actual pixel coordinates (no wrapping).
	The location around which to search is taken from the image origins.

**/
int			Bimage::find_peak(double radius, double sigma)
{
	if ( !d.uc ) return -1;
	
//	information();
//	write_img("cc.map", p, 0);
//	bexit(0);
	
	if ( radius < 0 ) radius = 1e30;
	else if ( radius < 1 ) radius = 0.9;
//	if ( sigma < 1e-30 ) sigma = 2e-30;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::find_peak: radius=" << radius << endl;
	
	long				i, nn, xx, yy, zz;
	long				xh(x/2), yh(y/2), zh(z/2);
	double				x2, y2, z2, r;
	double				value, maxval, fac(GOLDEN/sigma);
	Vector3<double>		ori, peak;
	
	if ( verbose & VERB_FULL ) {
	    cout << "Finding the brightest voxel:" << endl;
	    cout << "Radius & sigma:                 " << radius << " " << sigma << endl;
	    cout << "Origin:                         " << image->origin() << endl;
		cout << "n\tx\ty\tz\tv" << endl;
	}
	
	for ( i=nn=0; nn<n; nn++ ) {
		ori = image[nn].origin();
		if ( x <= 1 ) ori[0] = 0;
		if ( y <= 1 ) ori[1] = 0;
		if ( z <= 1 ) ori[2] = 0;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::find_peak: ori=" << ori << endl;
		maxval = -1e37;
		peak = ori;
		for ( zz=0; zz<z; zz++ ) {
			z2 = ori[2] - zz;
			if ( z2 < -zh ) z2 += z;
			if ( z2 >  zh ) z2 -= z;
			z2 *= z2;
			for ( yy=0; yy<y; yy++ ) {
				y2 = ori[1] - yy;
				if ( y2 < -yh ) y2 += y;
				if ( y2 >  yh ) y2 -= y;
				y2 *= y2;
				for ( xx=0; xx<x; xx++, i++ ) {
					x2 = ori[0] - xx;
					if ( x2 < -xh ) x2 += x;
					if ( x2 >  xh ) x2 -= x;
					x2 *= x2;
					r = sqrt(x2 + y2 + z2) - radius;
					value = (*this)[i];
					if ( sigma > 1e-30 ) value *= 1/(1 + exp(fac*r));
					if ( r <= 0 && maxval < value ) {
						maxval = value;
						peak = Vector3<double>(double(xx), double(yy), double(zz));
					}
//					if ( r < 0 )
//						cout << xx << tab << yy << tab << zz << tab << r << tab << value << tab << maxval << endl;
				}
			}
		}
		image[nn].FOM(maxval);
		image[nn].origin(peak);
		if ( verbose & VERB_DEBUG ) {
			cout << "DEBUG Bimage::find_peak: peak=" << peak << endl;
			cout << "DEBUG Bimage::find_peak: maxval=" << maxval << endl;
		}
		if ( maxval > 1 ) {
			cerr << "Error in find_peak: maxval > 1: " << maxval << endl;
			bexit(-1);
		}
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin() << tab << image[nn].FOM() << endl;
    }
	if ( verbose & VERB_FULL ) cout << endl;
	
	return 0;
}

/**
@brief 	Fits an elliptic parabole to locate the position of the peak to sub-voxel resolution.
@return Vector3<double> shift from input origin.

	The peak is expected to be near the middle of the image, close to the input origin.
	The function fits a 10-parameter parabola to the image.
	The shift is retrieved from the parameters by matrix inversion.
	The refined peak returned is the offset from the input origin.
**/
Vector3<double>	Bimage::fit_peak()
{
//	cout << "origin = " << image->origin() << endl;

//	write_img("cck.mrc", this, 0);
	
	long			i, j, k, nterm(9);
	long			xx, yy, zz;
	double			vm(get(0, image->origin()));
	Matrix			a(nterm,nterm);

	vector<double>	b(nterm);
	vector<double>	vec(nterm);
	
	for ( i=0; i<image_size(); ++i ) add(i, -vm);
	
	for ( j=0; j<nterm; ++j ) {	// Clear the matrix and vector
		for ( k=0; k<nterm; ++k ) a[j][k] = 0;
		b[j] = 0;
	}
	
	for ( i=zz=0; zz<z; ++zz ) {
		vec[2] = zz - image->origin()[2];
		vec[8] = vec[2]*vec[2];
		for ( yy=0; yy<y; ++yy ) {
			vec[1] = yy - image->origin()[1];
			vec[5] = vec[1]*vec[2];
			vec[7] = vec[1]*vec[1];
			for ( xx=0; xx<x; ++xx, ++i ) {
				vec[0] = xx - image->origin()[0];
				vec[3] = vec[0]*vec[1];
				vec[4] = vec[0]*vec[2];
				vec[6] = vec[0]*vec[0];
				for ( j=0; j<nterm; ++j ) b[j] += vec[j]*(*this)[i];
				for ( j=0; j<nterm; ++j )
					for ( k=0; k<=j; ++k ) a[j][k] += vec[j]*vec[k];
			}
		}
	}
		
	for ( j=0; j<nterm; ++j )
		for ( k=j+1; k<nterm; ++k ) a[j][k] = a[k][j];
	for ( j=0; j<nterm; ++j )
		if ( fabs(a[j][j]) < 1e-37 ) a[j][j] = 1;

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::fit_peak: Input matrix:" << endl;
		cout << a;
		cout << "DEBUG Bimage::fit_peak: Input vector:";
		for ( j=0; j<nterm; ++j ) cout << tab << b[j];
		cout << endl;
	}
	
//	a.LU_decomposition(b);
	a.singular_value_decomposition(b);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::fit_peak: Ouput matrix:" << endl;
		cout << a;
		cout << "DEBUG Bimage::fit_peak: Coefficients: ";
		for ( j=0; j<nterm; ++j ) cout << tab << b[j];
		cout << endl;
	}
	
	Matrix			c(3,3);
	vector<double>	d(3);
	for ( j=0; j<3; ++j ) d[j] = -b[j];
	c[0][0] = 2*b[6];
	c[0][1] = c[1][0] = b[3];
	c[0][2] = c[2][0] = b[4];
	c[1][1] = 2*b[7];
	c[1][2] = c[2][1] = b[5];
	c[2][2] = 2*b[8];
	
	c.LU_decomposition(d);
/*
	vector<double>		ev = c.jacobi_rotation();
	c.eigen_sort(ev);
	cout << c << endl;
	for ( auto it = ev.begin(); it != ev.end(); ++it ) cout << tab << *it;
	cout << endl;
*/	
	Vector3<double>	peak;
	
	double			cc = vm + b[0]*d[0] + b[1]*d[1] + b[2]*d[2]
		+ b[3]*d[0]*d[1] + b[4]*d[0]*d[2] + b[5]*d[1]*d[2]
		+ b[6]*d[0]*d[0] + b[7]*d[1]*d[1] + b[8]*d[2]*d[2];	

	if ( cc > 1 )
		cc = vm;
	else 
		peak = Vector3<double>(d[0], d[1], d[2]);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::fit_peak: Shift: ";
		for ( j=0; j<3; ++j ) cout << tab << d[j];
		cout << tab << cc << endl;
	}

	image->FOM(cc);
	if ( cc > 1 ) {
		cerr << "Error in fit_peak: cc = " << cc << endl;
		cerr << "peak = " << peak << endl;
		write_img("pf.mrc", this, 0);
		bexit(-1);
	}
	
/*
	double			R(0);
	Bimage*			pt = copy();
	for ( i=zz=0; zz<z; ++zz ) {
		d[3] = zz - image->origin()[2];
		for ( yy=0; yy<y; ++yy ) {
			d[2] = yy - image->origin()[1];
			for ( xx=0; xx<x; ++xx, ++i ) {
				d[1] = xx - image->origin()[0];
					cc = vm + b[0]*d[0] + b[1]*d[1] + b[2]*d[2]
						+ b[3]*d[0]*d[1] + b[4]*d[0]*d[2] + b[5]*d[1]*d[2]
						+ b[6]*d[0]*d[0] + b[7]*d[1]*d[1] + b[8]*d[2]*d[2];	
				pt->set(i, cc);
				cc -= (*this)[i];
				R += cc*cc;
			}
		}
	}
	write_img("cckc.mrc", pt, 0);
	delete pt;
	R /= size().volume();
	cout << "R2 = " << R << endl;
*/	
	
	return peak;
}

/**
@brief 	Refines the position of a peak to sub-voxel resolution.
@return int 		error code.

	The sub-voxel resolution peak in the vicinity of a voxel is defined 
	by fitting a 2D/3D second order function around the voxel.
	(typically used to find the shift vector in a cross-correlation map).

**/
int			Bimage::refine_peak_new()
{
	if ( !d.uc ) return -1;
	
	change_type(Float);
	
	long			nn, kernel_size(2);
	Vector3<long>	ksize(2*kernel_size+1,2*kernel_size+1,2*kernel_size+1);
	Vector3<double>	shift;
	Matrix3 		mat(1);

	if ( verbose & VERB_FULL ) {
	    cout << "Finding the shift in each image" << endl;
		cout << "n\tx\ty\tz" << endl;
	}
	
	for ( nn=0; nn<n; nn++ ) {
//		cout << "refine_peak: origin = " << image[nn].origin() << endl;
		Bimage*		pc = extract_wrap(nn, image[nn].origin(), ksize, ksize/2, mat);
		shift = image[nn].origin() + pc->fit_peak();
		shift = vector3_origin_to_shift(shift, size());
		image[nn].origin(shift);
		image[nn].FOM(pc->image->FOM());
		delete pc;
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << shift << endl;
	}
	if ( verbose & VERB_FULL ) cout << endl;

	return 0;
}

/**
@brief 	Refines the position of a peak to sub-voxel resolution.
@return int 		error code.

	The sub-voxel resolution peak in the vicinity of a voxel is defined 
	by fitting a 2D/3D second order function around the voxel.
	(typically used to find the shift vector in a cross-correlation map).

**/
int			Bimage::refine_peak()
{
	if ( !d.uc ) return -1;
	
	change_type(Float);
	
	long			i, j, nn, kernel_size(5), nterm(7);
	long			xx, yy, zz, xlo, xhi, ylo, yhi, zlo, zhi, ix, iy, iz;
	Vector3<long>	int_shift;
	Vector3<double>	shift;
	Matrix			a(nterm,nterm);
	vector<double>	b(nterm);
	vector<double>	vec(nterm);

	if ( verbose & VERB_FULL ) {
	    cout << "Finding the shift in each image" << endl;
		cout << "n\tx\ty\tz\tv" << endl;
	}
	
	for ( nn=0; nn<n; nn++ ) {
		shift = image[nn].origin();
		int_shift[0] = (long) shift[0];
		int_shift[1] = (long) shift[1];
		int_shift[2] = (long) shift[2];
		xlo = (long) (shift[0] - kernel_size);
		xhi = (long) (shift[0] + kernel_size);
		if ( xhi - xlo >= x ) {
			xlo = 0;
			xhi = x - 1;
		}
		ylo = (long) (shift[1] - kernel_size);
		yhi = (long) (shift[1] + kernel_size);
		if ( yhi - ylo >= y ) {
			ylo = 0;
			yhi = y - 1;
		}
		zlo = (long) (shift[2] - kernel_size);
		zhi = (long) (shift[2] + kernel_size);
		if ( zhi - zlo >= z ) {
			zlo = 0;
			zhi = z - 1;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::refine_peak: kernel limits = " << 
					xlo << tab << xhi << tab << ylo << tab << yhi << tab << zlo << tab << zhi << endl;
		for ( i=0; i<nterm; i++ ) {	// Clear the matrix and vector
			for ( j=0; j<nterm; j++ ) a[i][j] = 0;
			b[i] = 0;
		}
		vec[0] = 1;
		for ( zz=zlo; zz<=zhi; zz++ ) {
			iz = zz;
			if ( iz < 0 ) iz += z;
			if ( iz >= (long)z ) iz -= z;
			vec[3] = zz;
			vec[6] = zz*zz;
			for ( yy=ylo; yy<=yhi; yy++ ) {
				iy = yy;
				if ( iy < 0 ) iy += y;
				if ( iy >= (long)y ) iy -= y;
				vec[2] = yy;
				vec[5] = yy*yy;
				for ( xx=xlo; xx<=xhi; xx++ ) {
					ix = xx;
					if ( ix < 0 ) ix += x;
					if ( ix >= (long)x ) ix -= x;
					i = index(0,ix,iy,iz,nn);
					vec[1] = xx;
					vec[4] = xx*xx;
					if ( !isfinite((*this)[i]) ) {
						error_show("Error in Bimage::refine_peak", __FILE__, __LINE__);
						cerr << ix << " " << iy << " " << iz << " " << (*this)[i] << endl;
						return -9;
					}
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
		if ( verbose & VERB_DEBUG ) {
			cout << "DEBUG Bimage::refine_peak: Input matrix:" << endl;
			cout << a;
			cout << "DEBUG Bimage::refine_peak: Input vector:";
			for ( i=0; i<nterm; i++ ) cout << tab << b[i];
			cout << endl;
		}
		a.LU_decomposition(b);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::refine_peak: Coefficients: " <<
				b[0] << tab << b[1] << tab << b[2] << tab << b[3] << tab << b[4] << tab << b[5] << tab << b[6] << endl;
		if ( b[4] ) shift[0] = -0.5*b[1]/b[4];
		if ( b[5] ) shift[1] = -0.5*b[2]/b[5];
		if ( b[6] ) shift[2] = -0.5*b[3]/b[6];
		if ( fabs(shift[0] - int_shift[0]) > 2 ) shift[0] = image[nn].origin()[0];
		if ( fabs(shift[1] - int_shift[1]) > 2 ) shift[1] = image[nn].origin()[1];
		if ( fabs(shift[2] - int_shift[2]) > 2 ) shift[2] = image[nn].origin()[2];
		shift = vector3_origin_to_shift(shift, size());
		image[nn].origin(shift);
		if ( verbose & VERB_FULL )
			cout << n+1 << tab << shift << endl;
	}
	if ( verbose & VERB_FULL ) cout << endl;

	return 0;
}

/**
@brief 	Finds peaks in a map with periodic boundaries.
@param 	kernelsize		size of kernel side (must be ≥3).
@return Bimage*			image with peaks.

	A peak is defined as the local maximum within a specified kernel, where the
	maximum is also larger than a given threshold. This assumes that the region
	around a peak decreases monotonically on the scale length of the size of
	the kernel. With kernel edge of 3 (+-1 in all directions), the values around
	the peak shows a strict monotonicity with distance from the peak maximum.
	The peak positions and values are returned in the new image.

**/
Bimage*		Bimage::find_peaks(long kernelsize)
{
	if ( kernelsize < 3 ) kernelsize = 3;
	if ( kernelsize%2 == 0 ) kernelsize++;

    long			i, ii, ip, xx, yy, zz, nn, np(0);
	long			xlo, xhi, ylo, yhi, zlo, zhi, kx, ky, kz, ix, iy, iz;
	double			max;
	Vector3<long>	hk(size()/2);
	
	hk = hk.min(kernelsize/2);

	if ( verbose & VERB_PROCESS )
		cout << "Finding peaks:" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::find_peaks: hk=" << hk << endl;
	
	Bimage*			pfom = new Bimage(Float, TSimple, size(), n);
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			zlo = zz - hk[2];
			zhi = zz + hk[2];
			for ( yy=0; yy<y; yy++ ) {
				ylo = yy - hk[1];
				yhi = yy + hk[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					xlo = xx - hk[0];
					xhi = xx + hk[0];
					max = (*this)[i];
					if ( max > 0 ) {
						// Highest voxel search in the kernel
						ip = i;
						for ( kz=zlo; kz<=zhi; kz++ ) {
							iz = kz;
							if ( iz < 0 ) iz += z;	// Periodic boundaries
							if ( iz >= z ) iz -= z;
							for ( ky=ylo; ky<=yhi; ky++ ) {
								iy = ky;
								if ( iy < 0 ) iy += y;	// Periodic boundaries
								if ( iy >= y ) iy -= y;
								for ( kx=xlo; kx<=xhi; kx++ ) {
									ix = kx;
									if ( ix < 0 ) ix += x;	// Periodic boundaries
									if ( ix >= x ) ix -= x;
									ii = index(ix, iy, iz, nn);
									if ( max < (*this)[ii] ) {
										max = (*this)[ii];
										ip = ii;
									}
								}
							}
    				    }
						if ( ip == i ) {
							pfom->set(i, max);
							np++;
						}
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Peaks found:                    " << np << endl << endl;
	
	pfom->statistics();
	
	return pfom;
}

/**
@brief 	Finds the peaks in a cross-correlation map to find template matches.
@param 	excl_dist			exclusion distance between peaks.
@param 	&ncoor				number of coordinates.
@param 	&threshold_min		minimum threshold to choose peaks.
@param 	&threshold_max		maximum threshold to choose peaks.
@param 	pix_min				minimum peak width.
@param 	pix_max				maximum peak width.
@return Vector3<double>* 	list of coordinates.

	The map is searched in increments of the particle radius to identify
	peaks above the threshold and within a box the size of the
	particle radius. The identified peaks are further examined to eliminate 
	ones that are too close to a higher scoring peak. The acceptable distance
	between peaks is set to 1.8 times the particle radius.

**/
Vector3<double>*	Bimage::find_peaks(double excl_dist, long& ncoor, 
		double& threshold_min, double& threshold_max, double pix_min, double pix_max)
{
	if ( threshold_min < min ) threshold_min = min;
	if ( threshold_max > max ) threshold_max = max;
	
	long			rad = (long) excl_dist/2;

	Vector3<long>	tile(size());
	tile = tile.min(rad);
	
	Vector3<long>	start;
	start = start.max(tile);
	start = start.min(size()-1);
	
	Vector3<long>	end(size()-tile);
	end = end.max(1);
	
	Vector3<double>	ori(rad, rad, rad);
	ori = ori.min(size());
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Finding peaks in the cross-correlation map:" << endl;
		cout << "Exclusion distance:             " << excl_dist << endl;
		cout << "Thresholds:                     " << threshold_min << " - " << threshold_max << endl;
		cout << "Tile size:                      " << tile << endl;
		cout << "Starting coordinates:           " << start << endl;
		cout << "Ending coordinates:             " << end << endl;
	}
	
	long			i, k, l, nc, xx, yy, zz, tx, ty, tz, nn;
	long			nt = ((x-1)/tile[0]+1)*((y-1)/tile[1]+1)*((z-1)/tile[2]+1);
	double			max(-1e30), v;
	Vector3<long>*	ct = new Vector3<long>[nt];
	double*			cc = new double[nt];
	
	if ( verbose & VERB_FULL ) 
		cout << "Allocated number of coordinates: " << nt << endl;

	for ( nn=0, nc=0; nn<n; nn++ ) {
		for ( zz=start[2]; zz<end[2]; zz+=rad ) {
			for ( yy=start[1]; yy<end[1]; yy+=rad ) {
				for ( xx=start[0]; xx<end[0]; xx+=rad ) {
					cc[nc] = -1e30;
					for ( tz=zz; tz<zz+tile[2] && tz<end[2]; tz++ ) {
						for ( ty=yy; ty<yy+tile[1] && ty<end[1]; ty++ ) {
							for ( tx=xx; tx<xx+tile[0] && tx<end[0]; tx++ ) {
								i = index(0,tx,ty,tz,nn);
								v = (*this)[i];
								if ( v > threshold_min && v <= threshold_max && cc[nc] < v ) {
									cc[nc] = v;
									ct[nc][0] = tx;
									ct[nc][1] = ty;
									ct[nc][2] = tz;
								}
							}
						}
					}
					if ( cc[nc] > threshold_min && cc[nc] <= threshold_max ) {
//						sigma = peak_sigma(nn, ct[nc], 15);
//						cout << nc << tab << cc[nc] << tab << ct[nc] << tab << sigma << endl;
//						if ( sigma > pix_min && sigma < pix_max ) {
//							if ( max < cc[nc] ) max = cc[nc];
							nc++;
//						}
					} else cc[nc] = -1e30;
				}
			}
		}
	}

	if ( threshold_min < 1e-6 ) threshold_min = kmeans_threshold(nc, cc);
	image->FOM(threshold_min);

	if ( verbose & VERB_PROCESS ) {
		cout << "Number of peaks found:          " << nc << endl;
		cout << "Maximum CC:                     " << max << endl;
		cout << "Threshold:                      " << threshold_min << endl;
	}
	
	for ( k=0, ncoor=0; k<nc; k++ ) if ( cc[k] > threshold_min ) {
		for ( l=k+1; l<nc; l++ ) if ( cc[l] > threshold_min ) {
			if ( (ct[k] - ct[l]).length() < excl_dist ) {
				if ( cc[k] > cc[l] ) cc[l] = threshold_min - 1;
				else cc[k] = threshold_min - 1;
			}
		}
		if ( cc[k] > threshold_min ) ++ncoor;
	}
	
	Vector3<double>*	coor = new Vector3<double>[ncoor];
	
	for ( i=k=0; k<nc; k++ ) if ( cc[k] >= threshold_min )
		coor[i++] = ct[k];

	delete[] ct;
	delete[] cc;

	if ( verbose & VERB_PROCESS )
		cout << "Number of peaks found:          " << ncoor << endl << endl;
	
	return coor;
}

/**
@author Bernard Heymann and Xiange Zheng
@brief 	Calculates a confidence level to associate with a cross-correlation peak.
@param	nn			sub-image.
@return double		P value.

	The highest two peak maxima are found and used to calculate a
	confidence level for the hypothesis that the second highest peak is
	different from the first. The two maxima are found by a kernel search
	which implicitly assumes that the values associated with a peak
	monotonically decays with distance from the peak.
	Fisher's z-transform is used to calculate the confidence level:
		P = erfc((zmax - zmax2)/(2*z_sigma))
	where
		zmax, zmax2: z-transforms of the highest two peaks.
		z_sigma: standard deviation of the z-transform of the cc map.
	The P value is returned
**/
double		Bimage::ccmap_confidence(long nn)
{
	long			i, j, imax(0), imgsize(x*y*z);
	double			v, vmax(min), vmax2(min), zavg(0), zstd(0), pval(0);
	Vector3<long>	coor, cmax, cmax2;
	
	for ( i=nn*imgsize, j=0; j<imgsize; i++, j++ ) {
		v = (*this)[i];
		if ( vmax < v ) {
			vmax2 = vmax;
//			imax2 = imax;
			vmax = v;
			imax = i;
		} else if ( vmax2 < v ) {
			coor = coordinates(i);
			cmax = coordinates(imax);
			if ( coor.distance(cmax) > 1 ) {
				vmax2 = v;
//				imax2 = i;
			}
		}
		v = fishers_z_transform(v);
		zavg += v;
		zstd += v*v;
	}
	
	zavg /= imgsize;
	zstd = zstd/imgsize - zavg*zavg;
	if ( zstd <= 0 ) {
		zstd = 0;
	} else {
		zstd = sqrt(zstd);
		pval = erfc(0.5 * (fishers_z_transform(vmax) - fishers_z_transform(vmax2))/zstd);
	}
	if ( !isfinite(pval) ) pval = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::ccmap_confidence: z std:\t" << zstd << tab << 1/(sqrt(imgsize - 3.0)) << endl;
		cout << "DEBUG Bimage::ccmap_confidence: maxima:\t" << vmax << tab << vmax2 << endl;
	}
	
    return pval;
}

/**
@brief 	Calculates a sigma value for a cross-correlation peak.
@param	nn			sub-image.
@param	coor		coordinates in the image.
@param	kernel_size	kernel size to fit the gaussian.
@return double		sigma value.

	The peak is assumed to be close to a gaussian.
	The sigma is calculated for each kernel location and averaged:
		sigma = sqrt(rad^2/2(ln(vmax)-ln(v)))
		
**/
double		Bimage::peak_sigma(long nn, Vector3<long> coor, long kernel_size)
{
	long			xx, yy, zz, nv(0);
	double			sigma(0), dx2, dy2, dz2, r2, v, vmax(get(nn,coor));
	Vector3<long>	k(kernel_size, kernel_size, kernel_size);
	k = k.min(size());
	Vector3<long>	h(k/2), cmin(coor-h), cmax(coor+h+1);
	cmin = cmin.max(0);
	cmax = cmax.min(size());
	
//	cout << vmax << tab << cmin << tab << cmax << endl;
	
	for ( zz=cmin[2]; zz<cmax[2]; zz++ ) {
		dz2 = coor[2] - zz;
		dz2 *= dz2;
		for ( yy=cmin[1]; yy<cmax[1]; yy++ ) {
			dy2 = coor[1] - yy;
			dy2 *= dy2;
			for ( xx=cmin[0]; xx<cmax[0]; xx++ ) {
				dx2 = coor[0] - xx;
				dx2 *= dx2;
				r2 = dx2 + dy2 + dz2;
				if ( r2 ) {
					v = get(nn, xx, yy, zz);
//					cout << r2 << tab << vmax << tab << v << tab << log(vmax/v) << endl;
					if ( v > 0 && v < vmax ) {
						sigma += r2/log(vmax/v);
						nv++;
					}
				}
			}
		}
	}

	if ( nv ) sigma = sqrt(sigma/(2*nv));
//cout << nv << tab << sigma << endl;	
//bexit(-1);	
	return sigma;
}
