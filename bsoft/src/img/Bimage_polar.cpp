/**
@file	Bimage_polar.cpp
@brief	Library routines for polar and spherical transformations and calculations
@author Bernard Heymann
@date	Created: 19990904
@date	Modified: 20200509
**/

#include "Bimage.h"
#include "simplex.h"
#include "utilities.h"
#include "Matrix.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the radial average of an image.
@param 	minrad		minimum radius in voxels.
@param 	maxrad		maximum radius in voxels.
@param 	rad_step	step size in voxels.
@param	wrap		flag to wrap the data.
@return Bimage*		radial average in the form of a 1D image.

	A radial average of a 2D or 3D image is calculated between a minimum
	and maximum radius.  An interpolative method is used where the value of 
	a voxel is distributed between the two nearest radial annuli.
	The final sum in an annulus is normalized by the number of voxels 
	contributing to the annulus sum.
	The origins within the sub-image structures are used.

**/
Bimage* 	Bimage::radial(long minrad, long maxrad, double rad_step, int wrap)
{
	return radial(minrad, maxrad, rad_step, 1, 0, NULL, wrap);
}

/**
@brief 	Calculates the radial average of an image.
@param 	minrad		minimum radius in voxels.
@param 	maxrad		maximum radius in voxels.
@param 	rad_step	step size in voxels.
@param	*pmask		mask to limit calculation.
@param	wrap		flag to wrap the data.
@return Bimage*		radial average in the form of a 1D image.

	A radial average of a 2D or 3D image is calculated between a minimum
	and maximum radius.  An interpolative method is used where the value of 
	a voxel is distributed between the two nearest radial annuli.
	The final sum in an annulus is normalized by the number of voxels 
	contributing to the annulus sum.
	The origins within the sub-image structures are used.

**/
Bimage* 	Bimage::radial(long minrad, long maxrad, double rad_step, Bimage* pmask, int wrap)
{
	return radial(minrad, maxrad, rad_step, 1, 0, pmask, wrap);
}

/**
@brief 	Calculates the radial average of an image.
@param 	minrad		minimum radius in voxels.
@param 	maxrad		maximum radius in voxels.
@param 	rad_step	step size in voxels.
@param 	ellipticity	ratio of major and minor axes.
@param 	angle 		angle of major axis.
@param	*pmask		mask to limit calculation.
@param	wrap		flag to wrap the data.
@return Bimage*		radial average in the form of a 1D image.

	A radial average of a 2D or 3D image is calculated between a minimum
	and maximum radius.  An interpolative method is used where the value of 
	a voxel is distributed between the two nearest radial annuli.
	The final sum in an annulus is normalized by the number of voxels 
	contributing to the annulus sum.
	The origins within the sub-image structures are used.

**/
Bimage* 	Bimage::radial(long minrad, long maxrad, double rad_step,
				double ellipticity, double angle, Bimage* pmask, int wrap)
{
	long 			i, j, nn, xx, yy, zz, ir, nrad(0);
	double			dx, dy, dz, x2, y2, z2, a, r, f, f1;
	
	if ( rad_step < 0.1 ) rad_step = 1;
	
	Vector3<long>	h(size()/2);
	Vector3<double>	origin;
	Vector3<double>	scale(1/rad_step, 1/rad_step, 1/rad_step);
	if ( z < 2 ) scale[2] = 1;
	
	if ( maxrad < minrad ) swap(minrad, maxrad);
	if ( maxrad < 1 ) maxrad = size().max()/(2*rad_step);
	else if ( maxrad > size().max()/(2*rad_step) ) maxrad = size().max()/(2*rad_step);
	
	Bimage*			prad = new Bimage(Float, TSimple, maxrad, 1, 1, n);
	prad->sampling(rad_step*image->sampling()[0], 1, 1);

	vector<float>	num(maxrad);
	
	if ( pmask )
		if ( pmask->data_type() != UCharacter || pmask->minimum() != 0 || pmask->maximum() != 1 )
			pmask->to_mask();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a radial profile:" << endl;
		cout << "Minimum and maximum pixels:     " << minrad << " " << maxrad << endl;
		cout << "Radial step size:               " << rad_step << endl;
		cout << "Minimum and maximum radius:     " << minrad*image->sampling()[0] << " " << maxrad*image->sampling()[0] << " A" << endl;
		cout << "Radial step size:               " << rad_step*image->sampling()[0] << " A" << endl;
		if ( pmask )
			cout << "Mask file:                      " << pmask->file_name() << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Image\toriX\toriY\toriZ" << endl;
	for ( i=nn=0; nn<n; nn++ ) {
		for ( j=0; j<maxrad; j++ ) num[j] = 0;
		if ( x > 1 ) origin[0] = image[nn].origin()[0];
		if ( y > 1 ) origin[1] = image[nn].origin()[1];
		if ( z > 1 ) origin[2] = image[nn].origin()[2];
//		if ( (vector3_origin_to_shift(origin, size())).length() < x/10 ) wrap = 1;
//		else wrap = 0;
		if ( verbose & VERB_PROCESS )
			cout << nn+1 << tab << origin << endl;
		dz = -origin[2];
		if ( wrap && dz < -h[2] ) dz += z;
		for ( zz=0; zz<z; zz++, dz++ ) {
			if ( wrap && dz > h[2] ) dz -= z;
			z2 = dz * scale[2];
			z2 *= z2;
			dy = -origin[1];
			if ( wrap && dy < -h[1] ) dy += y;
			for ( yy=0; yy<y; yy++, dy++ ) {
				if ( wrap && dy > h[1] ) dy -= y;
				y2 = dy * scale[1];
				y2 *= y2;
				dx = -origin[0];
				if ( wrap && dx < -h[0] ) dx += x;
				for ( xx=0; xx<x; xx++, dx++, i++ ) if ( ( !pmask ) || ( pmask && (*pmask)[i] ) ) {
					if ( wrap && dx > h[0] ) dx -= x;
					x2 = dx * scale[0];
					x2 *= x2;
					r = sqrt(x2 + y2 + z2);
					if ( fabs(ellipticity - 1) > 0.001 ) {
						a = atan2(dy, dx) - angle;
						r *= sqrt((ellipticity*ellipticity+1)/
							(2*(ellipticity*ellipticity*cos(a)*cos(a)+
							sin(a)*sin(a))));
					}
					ir = (long) r;
					if ( ( ir >= minrad ) && ( ir < maxrad ) ) {
						f = r - ir;
						f1 = 1 - f;
						j = ir + nn*maxrad;	// Add for multiple images
						num[ir] += f1;
						prad->add(j, f1*(*this)[i]);
						ir++;
						if ( ir < maxrad ) {
							j++;
							num[ir] += f;
							prad->add(j, f*(*this)[i]);
						}
						nrad++;
					}
				}
			}
		}
		for ( j=0, ir = nn*maxrad; j<maxrad; j++, ir++ ) {
			if ( num[j] ) {
				prad->set(ir, (*prad)[ir]/num[j]);
			}
		}
		prad->image[nn].origin(0.0,0.0,0.0);
	}

	if ( verbose & VERB_PROCESS )
		cout << "Pixels in radial averages:      " << nrad << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage:radial: prad size = " << prad->sizeX() << endl;
	
	prad->statistics();
	
	return prad;
}

/**
@brief 	Calculates the symmetry-adjusted radial average of an image.
@param 	rad_start		minimum radius in voxels.
@param 	rad_end			maximum radius in voxels.
@param 	rad_step		step size in voxels.
@param 	spherical_fraction	ratio of major and minor axes.
@param	*sym			point group symmetry.
@return Bimage*			radial average in the form of a 1D image.

	A radial average of a 2D or 3D image is calculated between a minimum
	and maximum radius.  An interpolative method is used where the value of
	a voxel is distributed between the two nearest radial annuli.
	The final sum in an annulus is normalized by the number of voxels
	contributing to the annulus sum.
	The symmetry is used to adjust the radius to follow the contours of a shell.
	The origins within the sub-image structures are used.

**/
Bimage*			Bimage::radial_symmetry_adjusted(double rad_start,
					double rad_end, double rad_step,
					double spherical_fraction, Bsymmetry& sym)
{
	if ( rad_start < 0 ) rad_start = 0;
	if ( rad_start >= z/2 ) rad_start = z/2;
	if ( rad_end <= 0 ) rad_end = z/2;
	if ( rad_end >= z/2 ) rad_end = z/2 - 1;
	if ( rad_end < rad_start ) rad_end = rad_start;
	if ( rad_step <= 0 ) rad_step = 1;

	if ( set_radial_distances(spherical_fraction, sym) < 0 )
		return NULL;
	
	long			i, j, nn, ir;
	long			nr = (long) ((rad_end - rad_start)/rad_step + 1);
	double			r, f, f1;

	if ( verbose ) {
		cout << "Calculating radial profiles:" << endl;
		cout << "Radial range:                   " << rad_start << " - " << rad_end << endl;
		cout << "Radial step size:               " << rad_step << endl;
		cout << "Profile size:                   " << nr << endl << endl;
	}

	Bimage*			prad = new Bimage(Float, TSimple, nr, 1, 1, n);
	prad->sampling(sampling(0)[0] * rad_step, 1, 1);

	vector<double>	num(nr*n);
	
	for ( i=nn=0; nn<n; ++nn ) {
		for ( j=0; j<image_size(); ++j, ++i ) {
			r = ((*next)[i] - rad_start)/rad_step;
			ir = (long)r;
			f = r - ir;
			f1 = 1 - f;
			ir += nn*nr;
			if ( ir < nr ) {
				prad->add(ir, f1*(*this)[i]);
				num[ir] += f1;
				ir++;
				if ( ir < nr ) {
					prad->add(ir, f*(*this)[i]);
					num[ir] += f;
				}
			}
		}
	}

	for ( i=0; i<prad->data_size(); ++i )
		if ( num[i] ) prad->set(i, (*prad)[i]/num[i]);

	return prad;
}

/*
@brief Least squares function for the radial_fit function.
@param	&simp		a simplex structure.
@return double		the least squares R value.

	A radial profile is fitted to a reference profile using an iterative
	simplex down-hill method.  The equations solved are:
		xref = x*m
		yref(xref) = y(x)*a + b
	where	m is the dimension scaling (or magnification).
			a is the amplitude scaling.
			b is the amplitude shift.
	The least squares function is:
		         sum(weight*(yref(xref) - y(x)*a - b)^2)
		R = sqrt(---------------------------------------)
		                       sum(weight)
**/
double		radial_R(Bsimplex& simp)
{
	vector<double>&	y = simp.dependent_values();
	vector<double>&	yref = simp.independent_values();
	
	int			i, iref, n(0);
	double		fraction, xref, df, R(0), sum(0), weight;
	
	for ( i=0; i<simp.points(); i++ ) {
		xref = simp.parameter(0)*i;
		iref = (int) xref;
		if ( iref >= 0 && iref < simp.points() - 1 ) {
			fraction = xref - iref;
			df = (1.0 - fraction)*yref[iref] + fraction*yref[iref + 1]
					 - simp.parameter(1)*y[i] - simp.parameter(2);
//			weight = i*(fabs(yref[iref]) + fabs(y[i]));
			weight = i;
//			weight = fabs(y[i]);
			sum += weight;
			R += df*df*weight;
			n++;
		}
	}
	if ( sum && n > 0.2*simp.points() ) R = sqrt(R/sum);
	else R = 1e30;				// Out of range
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG radial_R: " <<  simp.parameter(0) << " " << 
				simp.parameter(1) << " " << simp.parameter(2) << " R=" << R << endl;
			
	return R;
}

/**
@brief 	Fits a radial profile to a reference radial profile.
@param 	*pref	reference radial profile as a 1D image.
@return double*	3-value result vector: dimension scale, amplitude scale and amplitude shift.

	A radial profile is fitted to a reference profile using an iterative
	simplex down-hill method.  The equations solved are:
		x(new) = x(old)*m
		y(new) = y(old)*a + b
	where	m is the dimension scaling (or magnification).
			a is the amplitude scaling.
			b is the amplitude shift.

**/
double*		Bimage::radial_fit(Bimage* pref)
{
	long 		i, npoint(x);
	if ( x < pref->sizeX() ) npoint = pref->sizeX();
	
	// 3 parameters: Dimension scaling, amplitude scaling, amplitude shift
	vector<double>	ax(npoint), fx(npoint);
	for ( i=0; i<npoint; i++ ) ax[i] = (*pref)[i];
	for ( i=0; i<npoint; i++ ) fx[i] = (*this)[i];

	Bsimplex	simp(1, 3, 0, npoint, ax, fx);
	
	// Dimension scaling estimated from the input units
	simp.parameter(0, image->sampling()[0]/pref->sampling(0)[0]);
	
	// Amplitude scaling estimated from the standard deviations
	simp.parameter(1, 1);
	if ( std != 0 && pref->standard_deviation() != 0 )
		simp.parameter(1, pref->standard_deviation()/std);
	else if ( max - min > 0 && pref->maximum() - pref->minimum() > 0 )
		simp.parameter(1, (pref->maximum() - pref->minimum())/(max - min));
	
	// Amplitude shift estimated from the averages
	simp.parameter(2, pref->average() - avg*simp.parameter(1));
//	if ( simp.parameter(2] < 1e-6 && max - min > 0 && pref->maximum() - pref->minimum() > 0 )
//		simp.parameter(2, (pref->maximum() - pref->minimum() + max - min)/2.0;
	
	// Parameter limits
	simp.limits(0, 0.8*simp.parameter(0), 1.2*simp.parameter(0));
	simp.limits(1, 0.1*simp.parameter(1), 10*simp.parameter(1));
	simp.limits(2, simp.parameter(2) - 10*pref->standard_deviation(), simp.parameter(2) + 10*pref->standard_deviation());

	if ( verbose & VERB_PROCESS ) {
		cout << "Radial fit initial parameters:" << endl;
		cout << "Dimension scale:                " << simp.parameter(0) << endl;
		cout << "Amplitude scale and shift:      " << simp.parameter(1) << " " << simp.parameter(2) << endl;
	}	

	// Do the downhill simplex iterative fit
	double		R = simp.run(10000, 1e-12*pref->standard_deviation(), radial_R) / pref->standard_deviation();
	
	int 		idat;
	double		fraction, xdat, f;

	if ( verbose & VERB_PROCESS ) {
		cout << "Radial fit final parameters:" << endl;
		cout << "Dimension scale:                " << simp.parameter(0) << endl;
		cout << "Relative dimension scale:       " << simp.parameter(0)*pref->sampling(0)[0]/image->sampling()[0] << endl;
		cout << "Amplitude scale and shift:      " << simp.parameter(1) << " " << simp.parameter(2) << endl;
		cout << "R:                              " << R << endl << endl;
		cout << "Radial averages:" << endl;
		cout << "Radius(A) " << file_name() << " " << pref->file_name() << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		for ( i=0; i<pref->sizeX(); i++ ) {
			xdat = i/simp.parameter(0);
			idat = (int) xdat;
			if ( idat >= 0 && idat < x - 1 ) {
				fraction = xdat - idat;
				f = simp.parameter(1)*((1.0 - fraction)*(*this)[idat]
							+ fraction*(*this)[idat+1]) + simp.parameter(2);
			} else {
				f = 0;
			}
			cout << pref->sampling(0)[0]*i << " " << f << " " << (*pref)[i] << endl;
		}
	}
	
	double*		result = new double[simp.parameters()];
	for ( i=0; i<simp.parameters(); i++ ) result[i] = simp.parameter(i);
	
	return result;
}

/**
@brief 	Generates a full 2D or 3D image from a radial profile in a 1D image.
@param 	nusize		size of image to expand to.
@param 	origin		origin for radial profile.
@return Bimage*		new 2D or 3D image.

	It assumes the resultant image is square or cubic.

**/
Bimage* 	Bimage::radial_to_full(Vector3<long> nusize, Vector3<double> origin)
{
	if ( nusize.volume() < 1 ) return NULL;
	
	Bimage*			p = copy_header();
	p->data_type(Float);
	p->size(nusize);
	if ( p->sizeX() < 2 ) origin[0] = 0;
	if ( p->sizeY() < 2 ) origin[1] = 0;
	if ( p->sizeZ() < 2 ) origin[2] = 0;
	p->data_alloc();
	
	long			i, nn, xx, yy, zz, ir, nexp(0);
	double			r, f, v;
	Vector3<double>	d;
	
	if ( verbose ) {
		cout << "Calculating a full image from a radial profile:" << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << "Origin:                         " << origin << endl;
	}
	
	for ( nn=i=0; nn<n; nn++ ) {
		for ( zz=0; zz<p->sizeZ(); zz++ ) {
			d[2] = zz - origin[2];
			for ( yy=0; yy<p->sizeY(); yy++ ) {
				d[1] = yy - origin[1];
				for ( xx=0; xx<p->sizeX(); xx++, i++ ) {
					d[0] = xx - origin[0];
					r = d.length(); 
					ir = (long) r;
					v = (*p)[i];
					if ( ir < x ) {
						f = 1 - r + ir;
						v += f*(*this)[ir];
						ir++;
						f = 1 - f;
						if ( ir < x ) v += f*(*this)[ir];
						else v += f*background(nn);
						nexp++;
					} else {
						v = background(nn);
					}
					p->set(i, v);
				}
			}
		}
		p->image[nn].origin(origin);
	}

	if ( verbose )
		cout << "Pixels generated:               " << nexp << endl << endl;
	
	p->statistics();
	
	return p;
}


/**
@brief 	Calculates an image with spherical coordinates.
@param 	nannuli 	number of annuli.
@param 	nphi		number of phi angles.
@param 	ntheta		number of theta angles.
@return Bimage*		spherical image.

	The image is converted to polar form with a specified number of annuli
	and angles in one (a 2D image) or two (a 3D map) directions. The 
	angular convention is consistent with the Euler angles used for 3D
	image processing. 
	A point vector p = {x,y,z} with respect to the image origin (typically
	the center of the image) is converted to polar coordinates as follows:
		The size of the vector gives the annulus:
			|p| = sqrt(x*x+y*y+z*z)
		Phi is the rotation angle around the z-axis, starting at the x-axis:
			phi = atan(y/x)
		Theta is the angle between the positive z-axis and the point vector:
			theta = acos(z/|p|)
	The new dimensions are mapped as follows with their maximum ranges:
		|p|   ===> x_dimension (0 - max(x_size,y_size,z_size)
		phi   ===> y_dimension (0 - 2*PI)
		theta ===> z_dimension (0 - PI)
	The sampling within these ranges are given by the calling function.
	The sampling must be isotropic.
	The origins within the sub-image structures are used.
	The interpolation routine actually calculates the old cartesian 
	coordinates for each set of spherical coordinates.

**/
Bimage*		Bimage::cartesian_to_spherical(long nannuli, long nphi, long ntheta)
{
	long   			j, nn, ann, phi, theta;
	double			maxrad, scale;
	
	maxrad = size().max()/2;
	
	if ( nannuli < 1 ) nannuli = (long) maxrad;
	if ( nphi < 1 ) nphi = 90;
	if ( y < 2 ) nphi = 1;
	if ( ntheta < 1 ) ntheta = 50;
	if ( z < 2 ) ntheta = 1;
	
	double			dannulus(maxrad*1.0/nannuli);	// Step size in each dimension
	double			dphi(TWOPI/nphi);
	double			dtheta(M_PI/ntheta);
	double			cos_phi, sin_phi, cos_theta, sin_theta;
	Vector3<double>	coor;
	Vector3<double>	u(image->sampling());

	Bimage*			psph = new Bimage(Float, TSimple, nannuli, nphi, ntheta, n);
	psph->sampling(dannulus*u[0], dphi*180.0/M_PI, dtheta*180.0/M_PI);
	u = psph->image->sampling();

	if ( verbose & VERB_PROCESS ) {
		cout << "Spherical image size:           " << psph->size() << endl;
		cout << "Step sizes:                     " << u[0] << " A\t"
			<< u[1] << " degree\t" << u[2] << " degree" << endl;
	}
	
	// Interpolate cartesian points to get spherical voxel values
	if ( verbose & VERB_FULL )
		cout << "Image\tOrigin" << endl;
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin() << endl;
		for ( theta=0; theta<ntheta; theta++ ) {
			cos_theta = dannulus*cos(dtheta*theta);
			sin_theta = dannulus*sin(dtheta*theta);
			for ( phi=0; phi<nphi; phi++ ) {
				cos_phi = cos(dphi*phi);
				sin_phi = -sin(dphi*phi);
				for ( ann=0; ann<nannuli; ann++ ) {
					j = psph->index(ann, phi, theta, nn);
					scale = ann*sin_theta;
					coor = Vector3<double>(scale*cos_phi, scale*sin_phi, ann*cos_theta);
					coor += image[nn].origin();
					psph->set(j, interpolate(coor, nn, 0));
				}
			}
		}
		psph->image[nn].origin(0.0, nphi/2, ntheta/2);
	}
	
	if ( verbose & VERB_PROCESS ) cout << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::cartesian_to_spherical: Calculating statistics" << endl;
	
	statistics();

	return psph;
}

/**
@brief 	Calculates an image with cylindrical coordinates.
@param 	nannuli 	number of annuli.
@param 	nphi		number of phi angles.
@param	flag		switch: annulus-angle-z (0), angle-annulus-z (1), z-annulus-angle (2)
@return Bimage*		cylindrical image.

	The image is converted to cylindrical form with a specified number of 
	annuli and angles. 
	A point vector p = {x,y,z} with respect to the image origin (typically
	the center of the image) is converted to cylindrical coordinates as follows:
		The distance from the z-origin gives the annulus:
			|p| = sqrt(x*x+y*y)
		Phi is the rotation angle around the z-axis, starting at the x-axis:
			phi = atan(y/x)
		The z-coordinate remains unchanged.
	The new dimensions are mapped as follows with their maximum ranges:
		|p|   ===> x_dimension (0 - max(x_size,y_size,z_size)
		phi   ===> y_dimension (0 - 2*PI)
		z     ===> z_dimension
	The sampling within these ranges are given by the calling function.
	The sampling must be isotropic.
	The origins within the sub-image structures are used.
	The interpolation routine actually calculates the old cartesian 
	coordinates for each set of cylindrical coordinates.

**/
Bimage*		Bimage::cartesian_to_cylindrical(long nannuli, long nphi, int flag)
{
	long   			j, nn, zz, ann, phi;
	double			maxrad;
	
	maxrad = x/2;			// Maximum radius taken from largest x and y dimension
	if ( maxrad < y/2 ) maxrad = y/2;
	
	if ( nannuli < 1 ) nannuli = (long) maxrad;
	if ( nphi < 1 ) nphi = 90;
	if ( y < 2 ) nphi = 1;
	
	double			dannulus(maxrad*1.0/nannuli);	// Step size in each dimension
	double			dphi(TWOPI/nphi);
	double			cos_phi, sin_phi;
	Vector3<double>	coor;
	Vector3<double>	u(image->sampling());
	
	Bimage*			pcyl;
	
	switch ( flag ) {
		case 0:
			pcyl = new Bimage(Float, TSimple, nannuli, nphi, z, n);
			pcyl->sampling(dannulus*u[0], dphi*180.0/M_PI, u[2]);
			break;
		case 1:
			pcyl = new Bimage(Float, TSimple, nphi, nannuli, z, n);
			pcyl->sampling(dphi*180.0/M_PI, dannulus*u[0], u[2]);
			break;
		case 2:
			pcyl = new Bimage(Float, TSimple, z, nannuli, nphi, n);
			pcyl->sampling(u[2], dannulus*u[0], dphi*180.0/M_PI);
			break;
		default:
			cerr << "Bimage::cartesian_to_cylindrical: flag = " << flag << " not supported!" << endl;
			return NULL;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Cylindrical image size:         " << pcyl->size() << endl;
		cout << "Step sizes:                     " << dannulus*image->sampling()[0] << " A\t"
			<< dphi*180.0/M_PI << " degree" << endl;
	}
	
	// Interpolate cartesian points to get spherical voxel values
	if ( verbose & VERB_FULL )
		cout << "Image\tOrigin" << endl;
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin() << endl;
		for ( zz=0; zz<z; zz++ ) {
			coor[2] = zz;
			for ( phi=0; phi<nphi; phi++ ) {
				cos_phi = cos(dphi*phi);
				sin_phi = sin(dphi*phi);
				for ( ann=0; ann<nannuli; ann++ ) {
					switch ( flag ) {
						case 0: j = pcyl->index(ann, phi, zz, nn); break;
						case 1: j = pcyl->index(phi, ann, zz, nn); break;
						case 2: j = pcyl->index(zz, ann, phi, nn); break;
						default: j=0;
					}
					coor[0] = ann*cos_phi + image[nn].origin()[0];
					coor[1] = ann*sin_phi + image[nn].origin()[1];
					pcyl->set(j, interpolate(coor, nn, 0));
				}
			}
		}
		if ( flag < 1 )
//			pcyl->image[nn].origin(0, nphi/2, image[nn].origin()[2]);
			pcyl->image[nn].origin(0.0, nphi/2, image[nn].origin()[2]);
		else
//			pcyl->image[nn].origin(nphi/2, 0, image[nn].origin()[2]);
			pcyl->image[nn].origin(nphi/2, 0.0, image[nn].origin()[2]);
	}
	
	if ( verbose & VERB_PROCESS ) cout << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::cartesian_to_cylindrical: Calculating statistics" << endl;
	
	statistics();

	return pcyl;
}

/**
@brief 	Calculates an image with cylindrical coordinates by integration.
@param 	nangles 	number of angles in each annulus.
@param 	ann_min 	minimum annulus (pixels).
@param 	ann_max 	maximum annulus (pixels).
@param 	dann	 	width of annulus (pixels).
@param 	zmin		minimum z (pixels).
@param 	zmax		maximum z (pixels).
@param 	zinc		increment in z (pixels).
@return Bimage*		cylindrical image.

	The image is converted to cylindrical form by integrating blocks 
		with a defined annular width and thickness in z at each angle. 
	The resultant image contains lines corresponding to integrated blocks
		covering 360° of angle.
	The sampling must be isotropic.
	The origins within the sub-image structures are used.
	The interpolation routine actually calculates the old cartesian 
	coordinates for each set of cylindrical coordinates.

**/
Bimage* 	Bimage::polar_transform(long nangles, long ann_min, long ann_max,
				long dann, long zmin, long zmax, long zinc)
{
	change_type(Float);
	
	if ( dann < 1 ) dann = 1;
	if ( ann_min >= x/2 ) ann_min = x/2;
	if ( ann_max < 1 || ann_max >= x/2 ) ann_max = x/2;
	if ( ann_max >= y/2 ) ann_max = y/2;
	if ( ann_min > ann_max ) swap(ann_min, ann_max);
	if ( ann_min == ann_max ) ann_max = ann_min + dann;
	
	if ( zinc < 1 ) zinc = 1;
	if ( zmin > zmax ) swap(zmin, zmax);
	if ( zmin == zmax ) zmax = zmin + zinc;
	if ( zmax >= z ) zmax = z - 1;

	long			i, ann, ang, r, zz, iz;
	long			nannuli((ann_max - ann_min)/dann + 1);
	long			nz((zmax - zmin)/zinc + 1);
	double			a, cosa, sina;
	double			dang(M_PI*2.0/nangles);
	Vector3<double>	coor;
	
	Bimage* 		p = new Bimage(Float, TSimple, nangles, nannuli, nz, 1);

	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a polar image from " << file_name() << endl;
		cout << "Image size:                     " << nangles << " x " << nannuli << " x " << nz << endl;
		cout << "Angular step size:              " << dang*180/M_PI << " degrees" << endl;
		cout << "Annular range:                  " << ann_min << " - " << ann_max << " pixels" << endl;
		cout << "Annular step size:              " << dann << " pixels" << endl;
		cout << "Z range:                        " << zmin << " - " << zmax << " pixels" << endl;
		cout << "Z step size:                    " << zinc << " pixels" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Calculating a polar image" << endl << endl;

	for ( i=0, zz=zmin; zz<=zmax; zz+=zinc ) {
		for ( ann=ann_min; ann<=ann_max; ann+=dann ) {
			for ( ang=0; ang<nangles; ang++, i++ ) {
				a = dang*ang;
				cosa = cos(a);
				sina = sin(a);
				for ( iz=zz; iz<zz+zinc; iz++ ) {
					for ( r=ann; r<ann+dann; r++ ) {
						coor[0] = r*cosa + image->origin()[0];
						coor[1] = r*sina + image->origin()[1];
						coor[2] = iz;
						p->add(i, interpolate(coor));
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::polar_transform: Calculating statistics" << endl;
	
	p->statistics();
	
	return p;
}

/**
@brief 	Calculates the polar power spectrum of a 2D transform amplitude or intensity image.
@param 	resolution	high resolution limit.
@param 	num_angle	number of angles to use for interpolation.
@return Bimage* 	polar power spectrum.

	The length of each annulus is number of angles times the radial offset.
	The samples on each annulus are calculated by bilinear interpolation.
	The maximum annulus calculated depends on the resolution limit specified
	in the input image.

**/
Bimage*		Bimage::polar_power_spectrum(double resolution, long num_angle)
{
	check_resolution(resolution);
	
	long				max_rad = (long) (x*image->sampling()[0]/resolution);
	if ( max_rad > x/2 ) max_rad = x/2;
	long				max_ang = num_angle*max_rad;
	Bimage*				pps = new Bimage(Double, TSimple, max_ang, max_rad, 1, n);
	
    long				j, nn, r, a;
	double				ang, xx, yy, fac;

	if ( verbose & VERB_FULL )
		cout << "Calculating a radial power spectrum: " << max_rad << " annuli x " << max_ang << " angles" << endl;

	for ( nn=0; nn<n; nn++ ) {
		for ( r=1; r<max_rad; r++ ) {
			max_ang = num_angle*r;
			fac = TWOPI/max_ang;
			for ( a=0; a<max_ang; a++ ) {
				j = (nn*pps->sizeY() + r)*pps->sizeX() + a;
				ang = a*fac;
				xx = r*cos(ang);
				yy = r*sin(ang);
				pps->set(j, interpolate_wrap(xx, yy, 0, nn));
			}
		}
	}

	return pps;
}

/**
@brief 	Calculates the 1D power spectrum of each line in an image.
@param 	plan		fft plan.
@return int			0.

	Each line is transformed to a power spectrum in place.

**/
int			Bimage::line_powerspectra(fft_plan plan)
{
	long				i, is, isx, xx, yy, zz;
	double				sum;
	Complex<float>*		data = new Complex<float>[x];
	
//	cout << "Doing 1D transforms" << endl;
	
	for ( is=zz=0; zz<z; zz++ ) {
		for ( yy=0; yy<y; yy++, is++ ) {
			isx = is*x;
			for ( xx=0, i=isx; xx<x; xx++, i++ )
				data[xx] = Complex<float>((*this)[i], 0);
			fftw(plan, data);
			set(isx, 0);	// Zero order
			for ( xx=1, i=isx+1, sum=0; xx<x; xx++, i++ ) {
				set(i, data[xx].power());
				sum += (*this)[i];
			}
			if ( sum ) for ( xx=1, i=isx+1; xx<x; xx++, i++ )
				set(i, (*this)[i]/sum);
		}
	}
	
	delete[] data;

	return 0;
}

/**
@brief 	Calculates an image with slices representing radial shell projections.
@return int		0.

	The 3D image is converted so that the sections contain 2D projections
	of radial shells. The projection in slice z is defined for:
		z >= sqrt((x-xo)^2 + (y-yo)^2) + zo
	The first projection is placed at the z origin and radiates out into
	both positive and negative directions.
	Sampling in x and y is not changed.
	The sampling must be isotropic.

**/
int		 	Bimage::radial_shells()
{
	change_type(Float);
	
	long   		i, j, nn, xx, yy, zz;
	long			rz;
	double			dx2, dy2, dz, dz2, d2, oldz, fraction;
		
	float*			shell = new float[datasize];
	
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << "Converting to radial shells" << endl << endl;
	
	// Interpolate cartesian points to get spherical voxel values
	if ( verbose & VERB_FULL )
		cout << "Image\tOrigin" << endl;
	for ( j=nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin() << endl;
		for ( zz=0; zz<z; zz++ ) {
			dz = zz - image[nn].origin()[2];
			dz2 = dz*dz;
			for ( yy=0; yy<y; yy++ ) {
				dy2 = yy - image[nn].origin()[1];
				dy2 *= dy2;
				for ( xx=0; xx<x; xx++, j++ ) {
					dx2 = xx - image[nn].origin()[0];
					dx2 *= dx2;
					d2 = dz2 - dy2 - dx2;
					if ( d2 >= 0 ) {
						if ( dz < 0 ) oldz = image[nn].origin()[2] - sqrt(d2);
						else oldz = image[nn].origin()[2] + sqrt(d2);
						rz = (long) oldz;
						if ( rz >= 0 && rz < z ) {
							fraction = 1 + rz - oldz;
							i = index(0,xx,yy,rz,nn);
							shell[j] = fraction*(*this)[i];
							rz++;
							if ( rz < z ) {
								fraction = 1 - fraction;
								i += y*x;
								shell[j] += fraction*(*this)[i];
							}
						}
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_FULL ) cout << "" << endl;

	data_assign((unsigned char *) shell);
	
	statistics();

	return 0;
}

/**
@brief 	Calculates an image with slices representing cylindrical shell projections.
@return int		0.

	The 3D image is converted so that the sections contain 2D projections
	of cylindrical shells. The projection in slice i is defined for:
		i >= sqrt(x^2 + y^2)
	The first projection is placed at the z origin and radiates out into
	both positive and negative directions.
	Sampling in x and y is not changed.
	The sampling must be isotropic.

**/
int		 	Bimage::cylindrical_shells()
{
	change_type(Float);
	
	long   			i, j, nn, xx, yy, zz;
	long			ry;
	double			dx, dy, dx2, dy2, d2, oldy, fraction;
		
	float*			shell = new float[datasize];
	
	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << "Converting to cylindrical shells" << endl << endl;
	
	// Interpolate cartesian points to get cylindrical voxel values
	if ( verbose & VERB_FULL )
		cout << "Image\tOrigin" << endl;
	for ( j=nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << nn+1 << tab << image[nn].origin() << endl;
		for ( zz=0; zz<z; zz++ ) {
//			dz = zz - image[nn].origin()[2];
//			dz2 = dz*dz;
			for ( yy=0; yy<y; yy++ ) {
				dy = yy - image[nn].origin()[1];
				dy2 = dy*dy;
				for ( xx=0; xx<x; xx++, j++ ) {
					dx = xx - image[nn].origin()[0];
					dx2 = dx*dx;
					d2 = dy2 - dx2;
					if ( dy2 >= dx2 ) {
						if ( dy < 0 ) oldy = image[nn].origin()[1] - sqrt(d2);
						else oldy = image[nn].origin()[1] + sqrt(d2);
						ry = (long) oldy;
						if ( ry >= 0 && ry < y ) {
							fraction = 1 + ry - oldy;
							i = index(0,xx,ry,zz,nn);
							shell[j] = fraction*(*this)[i];
							ry++;
							if ( ry < y ) {
								fraction = 1 - fraction;
								i += x;
								shell[j] += fraction*(*this)[i];
							}
						}
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_FULL ) cout << endl;

	data_assign((unsigned char *) shell);
	
	statistics();

	return 0;
}

/**
@brief 	Calculates the radius relative to the closest symmetry axis.
@param 	coor			coordinates in an image (voxels).
@param 	vs				symmetry axes.
@return double			symmetry-adjusted radial distance.

	The scalar production of the coordinate vector and each symmetry axis is calculated.
	The maximal value is returned as the symmetry-adjusted radial distance.

**/
double		radius_wrt_sym_axes(Vector3<double> coor, vector<Vector3<double>> vs)
{
	double				r, rmax(0);
	Vector3<double>		v(coor);
	v.normalize();
	
	for ( auto it = vs.begin(); it != vs.end(); ++it ) {
		r = v.scalar(*it);
		if ( rmax < r ) rmax = r;
	}

	return rmax;
}

/**
@brief 	Calculates the symmetry-adjusted radial distances.
@param 	spherical_fraction	fraction to impose spericity.
@param 	sym					symmetry.
@return int					0.

	The radius relative to the given symmetry is calculated for each voxel
	and adjusted with the provided spherical fraction..

**/
int			Bimage::set_radial_distances(double spherical_fraction, Bsymmetry& sym)
{
	change_type(Float);
	
	long				i, xx, yy, zz, nn;
	double				r, rf;
	Vector3<double>		coor;

	vector<Vector3<double>>	vs = symmetry_get_axes(sym);
	if ( vs.size() < 1 ) return -1;
	
	if ( sym.point() < 102 || vs.size() < 2 ) spherical_fraction = 1;
	
	if ( verbose ) {
		cout << "Calculating radial distances:" << endl;
		if ( sym.point() > 101 )
			cout << "Symmetry:                       " << sym.label() << " (" << vs.size() << ")" << endl;
		cout << "Spherical fraction:             " << spherical_fraction << endl;
	}

	next = new Bimage(Float, TSimple, size(), n);
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			coor[2] = zz - image[nn].origin()[2];
			for ( yy=0; yy<y; yy++ ) {
				coor[1] = yy - image[nn].origin()[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					coor[0] = xx - image[nn].origin()[0];
					r = coor.length();
					rf = 1;
					if ( vs.size() > 1 ) {
						rf = radius_wrt_sym_axes(coor, vs);
						rf = rf + spherical_fraction*(1 - rf);
						r = r*rf;
					}
					next->set(i, r);
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Calculates an image with slices giving the radial section projections.
@param 	prad				new image with radial sections.
@param 	nn					sub-image.
@param 	zr					radial section radius.
@param 	rad_start			starting radius (voxels).
@param 	rad_step			radial step size (voxels).
@param 	fill				value to fill in excluded regions.
@return int					0.

	The next linked image contains the symmetry-adjusted radius for each voxel.
	The section is then calculated for this radius with interpolation between the
	original z-slices.
**/
int			Bimage::radial_section(Bimage* prad, long nn,
				long zr, double rad_start, double rad_step, double fill)
{
	long			j, j1, xx, yy, zz, zj;
	long			slice_size = x*y;
	long			image_size = slice_size*z;
	long			i = (long) (nn*prad->size().volume() + zr*slice_size);
	double			r = rad_step*zr + rad_start;
	double			dr, drp;
	
	for ( yy=0; yy<y; yy++ ) {
		for ( xx=0; xx<prad->sizeX(); xx++, i++ ) {
			j = nn*image_size + yy*x + xx;
			for ( zz=0, zj=z+1, drp=1; zz<z && zj>=z; zz++, j+=slice_size ) {
				dr = (*next)[j] - r;
				if ( dr > 0 && drp <= 0 ) zj = zz;
				drp = dr;
			}
			if ( zj < z ) {
				j1 = nn*image_size + zj*slice_size + yy*x + xx;
				j = j1 - slice_size;
				dr = ((*next)[j1] - r)/((*next)[j1] - (*next)[j]);
				prad->set(i, dr*(*this)[j]);
				dr = 1 - dr;
				prad->add(i, dr*(*this)[j1]);
			} else {
				prad->set(i, fill);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Calculates an image with slices giving the radial section projections.
@param 	rad_start			starting radius (voxels).
@param 	rad_end				ending radius (voxels).
@param 	rad_step			radial step size (voxels).
@param 	spherical_fraction	fraction of spherical section (requires symmetry).
@param 	*sym				symmetry for non-spherical sections.
@param 	fill_type			FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill				value to fill in excluded regions.
@return Bimage*				image with radial sections.

	The 3D image is converted so that the sections contain 2D projections
	of radial shells. The projection in slice i is defined for:
		i >= sqrt(x^2 + y^2)
	Non-spherical sections can be generated by specifying the symmetry and
	the fraction of spherical nature (1=spherical, 0=based on symmetry).
	Sampling in x and y is not changed.
	The sampling in the input map must be isotropic.

**/
Bimage*		Bimage::radial_sections(double rad_start, double rad_end,
				double rad_step, double spherical_fraction,
				Bsymmetry& sym, int fill_type, double fill)
{
	if ( rad_start < 0 ) rad_start = 0;
	if ( rad_start >= z/2 ) rad_start = z/2;
	if ( rad_end <= 0 ) rad_end = z/2;
	if ( rad_end >= z/2 ) rad_end = z/2 - 1;
	if ( rad_end < rad_start ) rad_end = rad_start;
	if ( rad_step <= 0 ) rad_step = 1;

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) {
		if ( fabs(background(long(0))) < 1e-6 )
			calculate_background();
		fill = background(long(0));
	}

	if ( set_radial_distances(spherical_fraction, sym) < 0 )
		return NULL;
	
	long			nn, nz = (long) ((rad_end - rad_start)/rad_step + 1);

	if ( verbose ) {
		cout << "Calculating radial shells:" << endl;
		cout << "Radial range:                   " << rad_start << " - " << rad_end << endl;
		cout << "Radial step size:               " << rad_step << endl;
		cout << "New image size:                 " << x << " x " << y << " x " << nz << " x " << n << endl << endl;
	}

	Bimage*			prad = new Bimage(Float, TSimple, x, y, nz, n);
	prad->sampling(sampling(0)[0], sampling(0)[1], sampling(0)[2] * rad_step);
	
	for ( nn=0; nn<n; nn++ ) {
		prad->image[nn] = image[nn];
		prad->image[nn].origin(prad->image[nn].origin()[0], prad->image[nn].origin()[1], -rad_start);
#ifdef HAVE_GCD
		dispatch_apply(nz, dispatch_get_global_queue(0, 0), ^(size_t zz){
			radial_section(prad, nn, zz, rad_start, rad_step, fill);
		});
#else
#pragma omp parallel for
		for ( long zz=0; zz<nz; zz++ )
			radial_section(prad, nn, zz, rad_start, rad_step, fill);
#endif
	}
	
	prad->statistics();
	
	return prad;
}

/**
@brief 	Calculates the coverage in each radial shell.
@param 	threshold		density threshold to distinguish for- and background.
@param 	rad_step		radial step size (voxels).
@return int				0.

**/
Bimage*		Bimage::radial_coverage(double threshold, double rad_step)
{
	if ( rad_step <= 0 ) rad_step = 1;
	
	to_mask(threshold);
	
	Bimage* 	prad = radial(0, x/2, rad_step);

	if ( verbose ) {
		cout << "Calculating radial shell coverage:" << endl;
		cout << "Threshold:                      " << threshold << endl;
		cout << "Radial step size:               " << rad_step << endl;
		cout << "Radius\tCoverage" << endl;
	}

	long		i, xx, nn;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( xx=0; xx<prad->x; ++xx, ++i )
			cout << xx << tab << (*prad)[i] << endl;
	}
	
	delete prad;
	
	return 0;
}

