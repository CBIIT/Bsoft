/**
@file	Bimage_filter.cpp
@brief	Library routines used for filtering images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20160510
**/

#include "Bimage.h"
#include "matrix_linear.h"
#include "math_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

long		kmin(long i, long hk)
{
	return (i < hk)? hk - i: 0;
}

long		kmax(long i, long km, long m)
{
	long		hk = km/2;
	return (i + km < m + hk)? km: m + hk - i;
}

/**
@brief 	Generates a gaussian kernel image.
@param 	sigma		gaussian width.
@param 	max			kernel maximum.
@return int			0.
**/
int			Bimage::kernel_gaussian(double sigma, double max)
{
	if ( sigma < 0.001 ) sigma = 1;
	if ( fabs(max) < 1e-10 ) max = 1;

	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a gaussian kernel:" << endl;
		cout << "Size:                           " << size() << endl;
		cout << "Sigma:                          " << sigma << endl;
		cout << "Maximum:                        " << max << endl << endl;
	}
	
	Vector3<double>	h(size()/2);
	origin(h);

	long			i, xx, yy, zz;
	double			x2, y2, z2, s2;
	double			invsigma2(-0.5/(sigma*sigma));
	
	for ( i=zz=0; zz<z; zz++ ) {
		z2 = zz - h[2];
		z2 *= z2;
		for ( yy=0; yy<y; yy++ ) {
			y2 = yy - h[1];
			y2 *= y2;
			for ( xx=0; xx<x; xx++, i++ ) {
				x2 = xx - h[0];
				x2 *= x2;
				s2 = x2 + y2 + z2;
				set(i, max*exp(s2*invsigma2));
			}
		}
	}

	return 0;
}

/**
@brief 	Generates a laplacian-of-gaussian kernel.
@param 	sigma		gaussian width.
@param 	max			kernel maximum.
@return Bimage*		kernel image.
**/
int			Bimage::kernel_laplacian_of_gaussian(double sigma, double max)
{
	if ( sigma < 0.001 ) sigma = 1;
	if ( max < 1e-10 ) max = 1;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a laplacian-of-gaussian kernel:" << endl;
		cout << "Size:                           " << size() << endl;
		cout << "Sigma:                          " << sigma << endl;
		cout << "Maximum:                        " << max << endl << endl;
	}
	
	Vector3<double>	h(size()/2);
	origin(h);

	long			i, xx, yy, zz, ndim(0);
	if ( x > 1 ) ndim++;
	if ( y > 1 ) ndim++;
	if ( z > 1 ) ndim++;
	double			x2, y2, z2, s2;
	double			invsigma2 = 1/(sigma*sigma);
	double			invsigma4 = invsigma2*invsigma2;
	double			scale = max/(ndim*invsigma2);
	
	for ( i=zz=0; zz<z; zz++ ) {
		z2 = zz - h[2];
		z2 *= z2;
		for ( yy=0; yy<y; yy++ ) {
			y2 = yy - h[1];
			y2 *= y2;
			for ( xx=0; xx<x; xx++, i++ ) {
				x2 = xx - h[0];
				x2 *= x2;
				s2 = x2 + y2 + z2;
				set(i, scale*(s2*invsigma4 - ndim*invsigma2)*exp(-0.5*s2*invsigma2));
			}
		}
	}

	return 0;
}


int			Bimage::average_line_sum(long nn, long i, long ik, long nk)
{
	if ( size()[ik] < 2 ) return 0;
	
	long			ii(i*x), j, j1, j2, jc, cc, inc(1);
	long			nl(size()[ik]), ll(nl*c), hk(nk/2), nv(0);
	vector<double>	s(c,0);
	
	if ( ik == 1 ) {
		ii = (i%x) + (i/x) * x * y;
		inc = x;
	} else if ( ik == 2 ) {
		ii = i;
		inc = x * y;
	}
	
	ii += nn * x * y * z;
	
	double			step = inc*hk;
	vector<double>	v(ll);
	
//	cout << ii << tab << ll << tab << hk << tab << inc << tab << s[0] << endl;
	// Initial sum
	for ( j=ii, j1=0; j1<hk; j+=inc, j1++, nv++ )
		for ( cc=0, jc=j*c; cc<c; cc++, jc++ )
			s[cc] += (*this)[jc];
//	cout << "---" << endl;
	// Main sums
	for ( j=ii, j1=0; j1<nl; j+=inc, j1++ ) {
		if ( j1 < nl-hk ) {	// Add the leading values
			j2 = j + (long) step;
			for ( cc=0, jc=j2*c; cc<c; cc++, jc++ )
				s[cc] += (*this)[jc];
			nv++;
		}
		for ( cc=0, jc=j1*c; cc<c; cc++, jc++ )
			v[jc] = s[cc]/nv;
		if ( j1 >= hk ) {	// Subtract the trailing values
			j2 = j - (long) step;
			for ( cc=0, jc=j2*c; cc<c; cc++, jc++ )
				s[cc] -= (*this)[jc];
			nv--;
		}
	}

	for ( j=ii, j1=0; j1<ll; j+=inc )
		for ( cc=0, jc=j*c; cc<c; cc++, jc++, j1++ )
			set(jc, v[j1]);
	
	return 0;
}

int			Bimage::average_line_sums(long nn, long ik, long nk)
{
	if ( size()[ik] < 2 ) return 0;
	
	long		i, nls(1);
	
	// Count the number of lines
	for ( i=0; i<3; i++ ) if ( i != ik ) nls *= size()[i];
	
	if ( verbose & VERB_PROCESS )
		cout << tab << "Axis " << ik << " of size " << nls << endl;
	
	// Do each line
#ifdef HAVE_GCD
	dispatch_apply(nls, dispatch_get_global_queue(0, 0), ^(size_t i){
		average_line_sum(nn, i, ik, nk);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nls; i++ )
		average_line_sum(nn, i, ik, nk);
#endif

	return 0;
}

/**
@brief 	Applies an averaging filter to an image.
@param 	kernel_size		length of kernel edge (typically 3).
@return int 			0, <0 if error.

	A kernel of a given size is passed over the image and the average
	value within the kernel assigned to the central voxel.

**/
int			Bimage::filter_average(long kernel_size)
{
	Vector3<long>	ksize(kernel_size, kernel_size, kernel_size);
	return filter_average(ksize);
}

int			Bimage::filter_average(Vector3<long> k)
{
	change_type(Float);
	
	if ( k[0]%2 == 0 ) k[0] += 1;
	if ( k[1]%2 == 0 ) k[1] += 1;
	if ( k[2]%2 == 0 ) k[2] += 1;
	
	k = k.min(size());
	
	if ( verbose & VERB_LABEL ) {
		cout << "Applying an averaging filter:" << endl;
		cout << "Kernel size:                    " << k << " (" << k.volume() << ")" << endl << endl;
	}

	long			i, nn;
	
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn << endl;
		for ( i=0; i<3; i++ )
			average_line_sums(nn, i, k[i]);
	}
	
	if ( verbose & VERB_PROCESS ) cout << endl;
	
	statistics();
	
	return 0;
}

/**
@brief 	Calculates the central difference gradient image.
@return Bimage*		an image with 3-value gradient vectors.

	The central differences are calculated in the three orthogonal
	directions and written into 3-value vectors.

**/
Bimage*		Bimage::gradient()
{
	long			i, nn, xx, yy, zz, cc, xy(x*y);
	Vector3<float>	v1, v2;

	if ( verbose & VERB_PROCESS )
		cout << "Calculating a central difference gradient" << endl << endl;
	
	Bimage*			pg = new Bimage(Float, TVector3, size(), 1);
	pg->sampling(sampling(0));
	pg->origin(image->origin());

	for ( i=nn=0; nn<n; ++nn ) {
		for ( zz=0; zz<z; ++zz ) {
			for ( yy=0; yy<y; ++yy ) {
				for ( xx=0; xx<x; ++xx, ++i ) {
					for ( cc=0; cc<3; ++cc )
						v1[cc] = v2[cc] = (*this)[i];
					if ( xx ) v1[0] = (*this)[i-1];
					if ( xx < x-1 ) v2[0] = (*this)[i+1];
					if ( yy ) v1[1] = (*this)[i-x];
					if ( yy < y-1 ) v2[1] = (*this)[i+x];
					if ( zz ) v1[2] = (*this)[i-xy];
					if ( zz < z-1 ) v2[2] = (*this)[i+xy];
					pg->set(i, (v2 - v1)/2);
				}
			}
		}
	}
	
	pg->statistics();
	
	return pg;
}

Vector3<double>	Bimage::gradient_voxel(long i)
{
	long			xx, yy, zz, iz, iy, ix, jx, jy, jz;
	long			xy(x*y);
	double			len, v(0);
	Vector3<double>	g, w;
	Vector3<long>	coor = coordinates(i);
	Vector3<long>	coor2;
	
	zz = (coor[2])? -1: 0;
	for ( iz=coor[2]+zz; iz<z && zz<2; ++zz, ++iz ) {
		coor2[2] = zz;
		jz = iz*xy;
		yy = (coor[1])? -1: 0;
		for ( iy=coor[1]+yy; iy<y && yy<2; ++yy, ++iy ) {
			coor2[1] = yy;
			jy = jz + iy*x;
			xx = (coor[0])? -1: 0;
			for ( ix=coor[0]+xx; ix<x && xx<2; ++xx, ++ix ) {
				coor2[0] = xx;
				len = coor2.length();
				jx = jy + ix;
				if ( len ) {
					v = (*this)[jx]/len;
					g[0] += v*xx;
					g[1] += v*yy;
					g[2] += v*zz;
					w[0] += fabs(xx);
					w[1] += fabs(yy);
					w[2] += fabs(zz);
				}
			}
		}
	}
	if ( w[0] ) g[0] /= w[0];
	if ( w[1] ) g[1] /= w[1];
	if ( w[2] ) g[2] /= w[2];
	
	return g;
}

/**
@brief 	Calculates the difference gradient image in a 3x3 kernel.
@return Bimage*		an image with 3-value gradient vectors.

	The central differences are calculated in the three orthogonal
	directions and written into 3-value vectors.

**/
Bimage*		Bimage::gradient3x3()
{
	Bimage*			pg = new Bimage(Float, TVector3, size(), 1);
	pg->sampling(sampling(0));
	pg->origin(image->origin());

	if ( verbose & VERB_PROCESS )
		cout << "Calculating a 3x3 kernel gradient" << endl << endl;

#ifdef HAVE_GCD
	dispatch_apply(image_size(), dispatch_get_global_queue(0, 0), ^(size_t i){
		pg->set(i, gradient_voxel(i));
	});
#else
#pragma omp parallel for
	for ( long i=0; i<image_size(); ++i ) {
		pg->set(i, gradient_voxel(i));
	}
#endif
	
	pg->statistics();
	
	return pg;
}

int			Bimage::gaussian_line_sum(long nn, long i, long ik, Bimage* pk)
{
	if ( size()[ik] < 2 ) return 0;
	
	long		ii(i*x), j, j1, j2, jm, jc, ji, cc, inc(1), nl(size()[ik]), ll(nl*c), hk(pk->x/2);
	double		w;
	
	if ( ik == 1 ) {
		ii = (i%x) + (i/x) * x * y;
		inc = x;
	} else if ( ik == 2 ) {
		ii = i;
		inc = x * y;
	}
	
	ii += nn * x * y * z;
	jm = ii + nl*inc;
	
	vector<double>	v(ll);

	// Main sums
	for ( j1=0; j1<nl; j1++ ) {
		j2 = ( j1 < hk )? hk - j1: 0;
		j = ii;
		if ( j2 < 1 ) j += (j1 - hk)*inc;
		for ( cc=0, jc=j1*c; cc<c; cc++, jc++ )
			v[jc] = 0;
		for ( w=0; j2<pk->x && j<jm; j+=inc, j2++ ) {
			for ( cc=0, jc=j1*c, ji=j*c; cc<c; cc++, jc++, ji++ )
				v[jc] += (*pk)[j2] * (*this)[ji];
			w += (*pk)[j2];
		}
		for ( cc=0, jc=j1*c; cc<c; cc++, jc++ )
			v[jc] /= w;
	}
	
	// Transfer back
	for ( j=ii, j1=0; j1<ll; j+=inc )
		for ( cc=0, jc=j*c; cc<c; cc++, jc++, j1++ )
			set(jc, v[j1]);
	
	return 0;
}

int			Bimage::gaussian_line_sums(long nn, long ik, Bimage* pk)
{
	if ( size()[ik] < 2 ) return 0;
	
	long		i, nls(1);
	
	// Count the number of lines
	for ( i=0; i<3; i++ ) if ( i != ik ) nls *= size()[i];
	
	if ( verbose & VERB_PROCESS )
		cout << tab << "Axis " << ik << " of size " << nls << endl;
	
	// Do each line
#ifdef HAVE_GCD
	dispatch_apply(nls, dispatch_get_global_queue(0, 0), ^(size_t i){
		gaussian_line_sum(nn, i, ik, pk);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nls; i++ )
		gaussian_line_sum(nn, i, ik, pk);
#endif

	return 0;
}

/**
@brief 	Applies a gaussian filter to an image.
@param 	kernel_size		length of kernel edge.
@param	sigma			gaussian decay.
@return int 				0, <0 if error.

	The image is comvolved with a Gaussian kernel.
	If the sigma value is zero, it is set to a sixth of the kernel size.

**/
int			Bimage::filter_gaussian(long kernel_size, double sigma)
{
	if ( sigma < 1e-30 ) sigma = kernel_size/6.0;
	
	change_type(Float);

	if ( kernel_size%2 == 0 ) kernel_size += 1;
	
	Bimage*			pk = new Bimage(Float, TSimple, kernel_size, 1, 1, 1);
	pk->kernel_gaussian(sigma, 1);

	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a gaussian filter:" << endl;
		cout << "Kernel size:                    " << pk->x << endl;
		cout << "Sigma:                          " << sigma << endl << endl;
	}
	
	long			i, nn;
	
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn << endl;
		for ( i=0; i<3; i++ )
			gaussian_line_sums(nn, i, pk);
	}
	
	if ( verbose & VERB_PROCESS ) cout << endl;
	
	delete pk;
	
	statistics();
	
	return 0;
}

/**
@brief 	Weighs an image with a sinc function.
@return int			0.

	This filter compensates for trilinear interpolation during reconstruction.
	The origin must be properly defined.

**/
int			Bimage::filter_sinc()
{
	long			i, xx, yy, zz;
	double			r, w;
	Vector3<double>	d;
	
	if ( verbose )
		cout << "Weighing with a sinc function" << endl << endl;
	
	for ( i=zz=0; zz<z; zz++ ) {
		d[2] = (-image->origin()[2] + zz)/z;
		for ( yy=0; yy<y; yy++ ) {
			d[1] = (-image->origin()[1] + yy)/y;
			for ( xx=0; xx<x; xx++, i++ ) {
				d[0] = (-image->origin()[0] + xx)/x;
				r = M_PI*d.length();
				if ( r ) {
					w = r/sin(r);
					set(i, w * w * (*this)[i]);
				}
			}
		}
	}
	
	return 0;
}

int			Bimage::convolve_chunk(Bimage* pkernel, float* nudata, long i, long len)
{
	long			ic, j, jy, jz, jn, m, nn, xx, yy, zz, cc, kx, ky, kz, ix, iy, iz;
	Vector3<long>	hk(pkernel->image->origin());
	long			xmin, ymin, zmin, xmax, ymax, zmax;
	double			v, w;
	
	coordinates(i, cc, xx, yy, zz, nn);
	
	if ( cc != 0 )
		cerr << "Warning: chunk not on channel boundary! (" << i << ")" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::convolve_chunk: " << i << tab << cc << tab << xx
			<< tab << yy << tab << zz << tab << nn << endl;

	xmin = kmin(xx, hk[0]);
	xmax = kmax(xx, pkernel->x, x);
	ymin = kmin(yy, hk[1]);
	ymax = kmax(yy, pkernel->y, y);
	zmin = kmin(zz, hk[2]);
	zmax = kmax(zz, pkernel->z, z);
	
	m = i+len;
	if ( m > datasize ) m = datasize;
	
	for ( ; i<m; xx++, i+=c ) {
		if ( xx >= x ) {
			xx = 0;
			yy++;
			if ( yy >= y ) {
				yy = 0;
				zz++;
				if ( zz >= z ) {
					zz = 0;
					nn++;
				}
				zmin = kmin(zz, hk[2]);
				zmax = kmax(zz, pkernel->z, z);
			}
			ymin = kmin(yy, hk[1]);
			ymax = kmax(yy, pkernel->y, y);
		}
		xmin = kmin(xx, hk[0]);
		xmax = kmax(xx, pkernel->x, x);
		for ( cc=0, ic=i; cc<c; cc++, ic++ ) nudata[ic] = 0;
		w = 0;
		jn = nn*z;
		for ( kz=zmin; kz<zmax; kz++ ) {
			jz = (jn + zz + kz - hk[2])*y;
			iz = kz*pkernel->y;
			for ( ky=ymin; ky<ymax; ky++ ) {
				jy = (jz + yy + ky - hk[1])*x;
				iy = (iz + ky)*pkernel->x;
				for ( kx=xmin; kx<xmax; kx++ ) {
					j = (jy + xx + kx - hk[0])*c;
					ix = iy + kx;
					v = (*pkernel)[ix];
					for ( cc=0, ic=i; cc<c; cc++, ic++, j++ )
						nudata[ic] += (*this)[j]*v;
					w += v;
				}
			}
		}
		if ( w ) for ( cc=0, ic=i; cc<c; cc++, ic++ )
			nudata[ic] /= w;
	}
	
	return 0;
}

/**
@brief 	Convolves an image with an arbitrary size convolution filter.
@param 	*pkernel	kernel encoded as an image.
@return int 		0, <0 on error.

	The kernel is multiplied with each area surrounding the current voxel
	with wrapping to avoid image edge effects.
	The convolution is threaded if compiled with OpenMP.

**/
int			Bimage::convolve(Bimage* pkernel)
{
	if ( !pkernel || !pkernel->data_pointer() ) {
		cerr << "Error in Bimage::convolve: No kernel specified!" << endl << endl;
		return -1;
	}
	
	if ( pkernel->sizeX() > x || pkernel->sizeY() > y || pkernel->sizeZ() > z ) {
		error_show("Error in Bimage::convolve", __FILE__, __LINE__);
		cerr << "Error: Kernel size does not fit the image size:" << endl;
		cerr << "   Kernel size = " << pkernel->size() << endl;
		cerr << "   Image size  = " << size() << endl;
		return -1;
	}
	
	change_type(Float);
	pkernel->change_type(Float);

	if ( verbose & VERB_PROCESS ) {
		cout << "Convolving with a kernel:       " << pkernel->file_name() << endl;
		cout << "Kernel size:                    " << pkernel->size() << endl << endl;
	}

	long			i, xx, yy, chunk_size(get_chunk_size(datasize, c));
	
	float*			nudata = new float[datasize];

//	if ( verbose & VERB_FULL ) {
		cout << "Kernel:" << endl;
		for ( i=yy=0; yy<pkernel->sizeY()*pkernel->sizeZ(); yy++ ) {
			for ( xx=0; xx<pkernel->sizeX(); xx++, i++ )
				cout << tab << (*pkernel)[i];
			cout << endl;
		}
		cout << endl;
//	}
	
#ifdef HAVE_GCD
	dispatch_apply((datasize - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		convolve_chunk(pkernel, nudata, i*chunk_size, chunk_size);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<datasize; i+=chunk_size )
		convolve_chunk(pkernel, nudata, i, chunk_size);
#endif

	data_assign((unsigned char *) nudata);
	
	statistics();

	return 0;
}

/**
@brief 	Convolves the image with an orthogonal kernel with wrapping.
@param	type		type of kernel: 0=gradient magnitude, 1=laplacian.
@return int 		0, <0 on error.

	The gradient kernel is:
		0  0  0    0 -1  0    0  0  0
		0 -1  0   -1  0  1    0  1  0
		0  0  0    0  1  0    0  0  0

	The Laplacian filter kernel for a 3D volume is:
		0  0  0    0  1  0    0  0  0
		0  1  0    1 -6  1    0  1  0
		0  0  0    0  1  0    0  0  0
	For 1D and 2D the central value is -2 and -4 respectively.

**/
int 		Bimage::filter_ortho(int type)
{
    change_type(Float);
	
	long	   		i, xx, yy, zz, nn;
	double			value;	
	long   			slicesize(x*y*c);
	float*			nudata = new float[datasize];
	long			ksize(1);
	if ( x > 1 ) ksize += 2;
	if ( y > 1 ) ksize += 2;
	if ( z > 1 ) ksize += 2;
	vector<double>	kernel(ksize);
	
	if ( type == 0 ) {
		kernel[0] = 0;
		for ( i=1; i<7; i+=2 )
			if ( ksize > i ) {
				kernel[i] = -1;
				kernel[i+1] = 1;
			}
	} else if ( type == 1 ) {
		kernel[0] = 0;
		for ( i=1; i<7; i+=2 )
			if ( ksize > i ) {
				kernel[0] -= 2;
				kernel[i] = 1;
				kernel[i+1] = 1;
			}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Convolving with an orthogonal kernel:" << endl;
		cout << "Kernel type:                    ";
		switch ( type ) {
			case 0: cout << "gradient"; break;
			case 1: cout << "laplacian"; break;
			default: break;
		}
		cout << endl << endl;
	}

    for ( i=nn=0; nn<n; nn++ ) {
	    for ( zz=0; zz<z; zz++ ) {
    		for ( yy=0; yy<y; yy++ ) {
    	    	for ( xx=0; xx<x; xx++, i++ ) {
					value = kernel[0]*(*this)[i];
					if ( x > 1 ) {
						if ( xx > 0 ) value += kernel[1]*(*this)[i-1];
						else value += kernel[1]*(*this)[i+x-1];
						if ( xx < x - 1 ) value += kernel[2]*(*this)[i+1];
						else value += kernel[2]*(*this)[i+1-x];
						if ( type == 0 ) {
							nudata[i] = value*value;
							value = 0;
						}
					}
					if ( y > 1 ) {
						if ( yy > 0 ) value += kernel[3]*(*this)[i-x];
						else value += kernel[3]*(*this)[i+slicesize-x];
						if ( yy < y - 1 ) value += kernel[4]*(*this)[i+x];
						else value += kernel[4]*(*this)[i+x-slicesize];
						if ( type == 0 ) {
							nudata[i] += value*value;
							value = 0;
						}
					}
					if ( z > 1 ) {
						if ( zz > 0 ) value += kernel[5]*(*this)[i-slicesize];
						else value += kernel[5]*(*this)[i+datasize-slicesize];
						if ( zz < z - 1 ) value += kernel[6]*(*this)[i+slicesize];
						else value += kernel[6]*(*this)[i+slicesize-datasize];
						if ( type == 0 ) {
							nudata[i] += value*value;
							value = 0;
						}
					}
					if ( type == 1 ) nudata[i] = value;
					else nudata[i] = sqrt(nudata[i]);
				}
			}
		}
	}
	
	data_assign((unsigned char *) nudata);
	
	statistics();
		
	return 0;
}

/**
@brief 	Convolves the image with a difference of gausians kernel.
@param	sigma1		sigma for the inner gaussian.
@param	sigma2		sigma for the outer gaussian.
@return int 		0, <0 on error.
**/
int			Bimage::filter_dog(double sigma1, double sigma2)
{
	if ( sigma1 > sigma2 ) swap(sigma1, sigma2);

	long			k(long(6*sigma2+1));
	Vector3<long>	ksize(k, k, k);
	ksize = ksize.min(size());
	
	double			m1(1/(sqrt(TWOPI)*sigma1));
	double			m2(1/(sqrt(TWOPI)*sigma2));
	
	Bimage*			pk = new Bimage(Float, TSimple, ksize, 1);
	Bimage*			pk2 = new Bimage(Float, TSimple, ksize, 1);
	
	pk->kernel_gaussian(sigma1, m1);
	pk2->kernel_gaussian(sigma2, m2);
	pk->subtract(pk2);
	
	delete pk2;

	convolve(pk);
	
	delete pk;

	statistics();
	
	return 0;
}

int			Bimage::filter_bilateral_chunk(Bimage* pkernel, double sigma2,
				int kernel_type, float* nudata, long i, long len, int first)
{
	long			nn, xx, yy, zz, cc;
	long			ic, j, jy, jz, jn, m, kx, ky, kz, ix, iy, iz, cnt(0);
	Vector3<long>	hk(pkernel->image->origin());
	long			xmin, ymin, zmin, xmax, ymax, zmax;

	double			diff2, f, ref, acc, acc_num, acc_denom;
	double			gauss_neginv2sig2(0), lorentz_invsig2(0), tukey_invsig2(0);
	if ( kernel_type == 1) gauss_neginv2sig2 = -0.5/(sigma2*sigma2);
	else if ( kernel_type == 2) lorentz_invsig2 = 1./(sigma2*sigma2);
	else if ( kernel_type == 3) tukey_invsig2 = 0.2/(sigma2*sigma2);
	
	coordinates(i, cc, xx, yy, zz, nn);
	
	if ( cc != 0 )
		cerr << "Warning: chunk not on channel boundary! (" << i << ")" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::filter_bilateral_chunk: " << i << tab << cc << tab << xx
			<< tab << yy << tab << zz << tab << nn << endl;

	xmin = kmin(xx, hk[0]);
	xmax = kmax(xx, pkernel->x, x);
	ymin = kmin(yy, hk[1]);
	ymax = kmax(yy, pkernel->y, y);
	zmin = kmin(zz, hk[2]);
	zmax = kmax(zz, pkernel->z, z);
	
	m = i+len;
	if ( m > datasize ) m = datasize;
	
	for ( ; i<m; xx++, i+=c ) {
		if ( xx >= x ) {
			xx = 0;
			yy++;
			if ( yy >= y ) {
				yy = 0;
				zz++;
				if ( zz >= z ) {
					zz = 0;
					nn++;
				}
				zmin = kmin(zz, hk[2]);
				zmax = kmax(zz, pkernel->z, z);
			}
			ymin = kmin(yy, hk[1]);
			ymax = kmax(yy, pkernel->y, y);
		}
		xmin = kmin(xx, hk[0]);
		xmax = kmax(xx, pkernel->x, x);
		for ( cc=0, ic=i; cc<c; cc++, ic++ ) nudata[ic] = 0;
		ref = (*this)[i];
		acc_num = acc_denom = 0.0;
		jn = nn*z;
		for ( kz=zmin; kz<zmax; kz++ ) {
			jz = (jn + zz + kz - hk[2])*y;
			iz = kz*pkernel->y;
			for ( ky=ymin; ky<ymax; ky++ ) {
				jy = (jz + yy + ky - hk[1])*x;
				iy = (iz + ky)*pkernel->x;
				for ( kx=xmin; kx<xmax; kx++ ) {
					j = (jy + xx + kx - hk[0])*c;
					ix = iy + kx;
					diff2 = ref - (*this)[j];
					diff2 *= diff2;
					acc = 0;
					switch (kernel_type) {
						case 1:			// gaussian
							acc = exp(gauss_neginv2sig2*diff2);
							break;
						case 2:			// lorentz
							acc = 1/(1+lorentz_invsig2*diff2);
							break;
						case 3:			// tukey
							if (diff2 <= tukey_invsig2) {
								f = 1 - tukey_invsig2*diff2;
								acc = 0.5*f*f;
							};
							break;
					}
					acc *= (*pkernel)[ix];
					acc_num += acc*(*this)[j];
					acc_denom += acc;
				}
			}
		}
		for ( cc=0, ic=i; cc<c; cc++, ic++ )
			if (acc_denom > 0) nudata[ic] = acc_num/acc_denom;
		if ( first && ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) ) )
			cerr << "Complete:                       " << setprecision(3)
							<< (++cnt)*100.0/len << " %    \r" << flush;
	}
	
	return 0;
}

/**
@brief 	Denoise an image with combined gaussian distance and density difference kernel.
@param 	sigma1			sigma for distance weighting function.
@param 	sigma2			sigma for density difference weighting function.
@param 	kernel_type		kernel type for range filter.
@param 	kernel_radius	kernel radius.
@return int 			0.

	The kernel is multiplied with each area surrounding the current voxel.
	Kernel types:
		1.	Gaussian
		2.	Lorentzian
		3.	Tukey

**/
int 		Bimage::filter_bilateral(double sigma1, double sigma2,
				int kernel_type, long kernel_radius)
{
	// Initialize parameters
	if ( sigma1 <= 0 ) sigma1 = 1;
	if ( kernel_radius <= 0 ) kernel_radius = (long) (sigma1*3);
	if ( sigma2 <= 0 ) sigma2 = std;
		
	Vector3<long>	ksize(kernel_radius*2+1, kernel_radius*2+1, kernel_radius*2+1);
	ksize = ksize.min(size());

	Bimage*			pkernel = new Bimage(Double, TSimple, ksize, 1);
	pkernel->kernel_gaussian(sigma1, 1);
	
    change_type(Float);
	
	if ( verbose & VERB_LABEL ) {
		cout << "Denoising with a bilateral kernel:" << endl;
		cout << "Sigma values:                   " << sigma1 << " " << sigma2 << endl;
		cout << "Kernel type:                    " << kernel_type << endl;
		cout << "Kernel size:                    " << ksize << endl << endl;
	}

	long			chunk_size(get_chunk_size(datasize, c));
	
	float*			nudata = new float[datasize];
	
#ifdef HAVE_GCD
	dispatch_apply((datasize - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		filter_bilateral_chunk(pkernel, sigma2, kernel_type, nudata, i*chunk_size, chunk_size, i==0);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<datasize; i+=chunk_size )
		filter_bilateral_chunk(pkernel, sigma2, kernel_type, nudata, i, chunk_size, i==0);
#endif

	delete pkernel;

	data_assign((unsigned char *) nudata);
	
	statistics();
		
	return 0;
}

/**
@brief 	Apply a rolling ball filter.
@param 	radius		radius of rolling ball.
@param 	scale		density scale.
@return int 		0.
**/
int 		Bimage::filter_rolling_ball(long radius, double scale)
{
	if ( radius < 1 ) radius = 1;
	if ( fabs(scale) < 1e-30 ) scale = 1;
	
	long   			i, xx, yy, zz, nn, kx, ky, kz;
	long			ix, iy, iz, jx, jy, jz;
	double			r2(radius*radius);
	double			dx2, dy2, dz2, diff, bcenter;
	
	Vector3<long>	ksize(radius*2+1, radius*2+1, radius*2+1);
	ksize = ksize.min(size());
	Vector3<long>	hk(ksize/2);
	Vector3<long>	lo, hi;
	
	long			ball_size = (long) ksize.volume();
	vector<double>	ball(ball_size);
	
	for ( i=zz=0; zz<ksize[2]; zz++ ) {
		dz2 = (double)zz - hk[2];
		dz2 *= dz2;
		for ( yy=0; yy<ksize[1]; yy++ ) {
			dy2 = (double)yy - hk[1];
			dy2 *= dy2;
			for ( xx=0; xx<ksize[0]; xx++, i++ ) {
				dx2 = (double)xx - hk[0];
				dx2 *= dx2;
				diff = r2 - dx2 - dy2 - dz2;
				if ( diff > 0 ) ball[i] = sqrt(diff)/scale;
				else ball[i] = 0;
			}
		}
	}
	
    change_type(Float);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Rolling ball filter:" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Scaling:                        " << scale << endl << endl;
	}

	vector<float>	back(datasize);
	for ( i=0; i<datasize; i++ ) back[i] = (*this)[i];
	
    for ( i=nn=0; nn<n; nn++ ) {
	    for ( zz=0; zz<z; zz++ ) {
			lo[2] = kmin(zz, hk[2]);
			hi[2] = kmax(zz, ksize[2], z);
			for ( yy=0; yy<y; yy++ ) {
				lo[1] = kmin(yy, hk[1]);
				hi[1] = kmax(yy, ksize[1], y);
				for ( xx=0; xx<x; xx++, i++ ) {
					lo[0] = kmin(xx, hk[0]);
					hi[0] = kmax(xx, ksize[0], x);
					if ( scale > 0 ) bcenter = max;
					else bcenter = min;
					for ( kz=lo[2]; kz<hi[2]; kz++ ) {
						iz = (nn*z + zz + kz - hk[2])*y;
						jz = kz*ksize[1];
						for ( ky=lo[1]; ky<hi[1]; ky++ ) {
							iy = (iz + yy + ky - hk[1])*x;
							jy = (jz  + ky)*ksize[0];
							for ( kx=lo[0]; kx<hi[0]; kx++ ) {
								ix = iy + xx + kx - hk[0];
								jx = jy + kx;
								diff = (*this)[ix] - (bcenter + ball[jx]);
								if ( ball[jx] > 0 ) {
									if ( diff < 0 ) bcenter += diff;	// Moving the ball down
								} else if ( ball[jx] < 0 ) {
									if ( diff > 0 ) bcenter += diff;	// Moving the ball up
								}
							}
						}
					}
					for ( kz=lo[2]; kz<hi[2]; kz++ ) {
						iz = (nn*z + zz + kz - hk[2])*y;
						jz = kz*ksize[1];
						for ( ky=lo[1]; ky<hi[1]; ky++ ) {
							iy = (iz + yy + ky - hk[1])*x;
							jy = (jz  + ky)*ksize[0];
							for ( kx=lo[0]; kx<hi[0]; kx++ ) {
								ix = iy + xx + kx - hk[0];
								jx = jy + kx;
								diff = bcenter + ball[jx];
								if ( ball[jx] > 0 ) {
									if ( back[ix] > diff ) back[ix] = diff;	// Assigning the background
								} else if ( ball[jx] < 0 ) {
									if ( back[ix] < diff ) back[ix] = diff;	// Assigning the background
								}
							}
						}
					}
				}
			}
		}
	}

	for ( i=0; i<datasize; i++ ) set(i, (*this)[i] - back[i]);
	
	statistics();
		
	return 0;
}

int			Bimage::filter_rank_chunk(long kernel_size, double rank, float* nudata, long i, long len)
{
	Vector3<long>	blocksize(kernel_size, kernel_size, kernel_size);
	blocksize = blocksize.min(size());
	Vector3<long>	hk(blocksize/2);
	long			blocktotalsize = (long) blocksize.volume();
	vector<double>	block(blocktotalsize);

	long			nn, xx, yy, zz, cc;
	long			j, jy, jz, jn, m, kx, ky, kz;
	long			xmin, ymin, zmin, xmax, ymax, zmax;
	long			nblock, imed;
	
	coordinates(i, cc, xx, yy, zz, nn);
	
	if ( cc != 0 )
		cerr << "Warning: chunk not on channel boundary! (" << i << ")" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::filter_rank_chunk: " << i << tab << cc << tab << xx
			<< tab << yy << tab << zz << tab << nn << endl;

	xmin = kmin(xx, hk[0]);
	xmax = kmax(xx, blocksize[0], x);
	ymin = kmin(yy, hk[1]);
	ymax = kmax(yy, blocksize[1], y);
	zmin = kmin(zz, hk[2]);
	zmax = kmax(zz, blocksize[2], z);
	
	m = i+len;
	if ( m > datasize ) m = datasize;
	
	for ( ; i<m; xx++, i+=c ) {
		if ( xx >= x ) {
			xx = 0;
			yy++;
			if ( yy >= y ) {
				yy = 0;
				zz++;
				if ( zz >= z ) {
					zz = 0;
					nn++;
				}
				zmin = kmin(zz, hk[2]);
				zmax = kmax(zz, blocksize[2], z);
			}
			ymin = kmin(yy, hk[1]);
			ymax = kmax(yy, blocksize[1], y);
		}
		xmin = kmin(xx, hk[0]);
		xmax = kmax(xx, blocksize[0], x);
		nblock = 0;
		jn = nn*z;
		for ( kz=zmin; kz<zmax; kz++ ) {
			jz = (jn + zz + kz - hk[2])*y;
			for ( ky=ymin; ky<ymax; ky++ ) {
				jy = (jz + yy + ky - hk[1])*x;
				for ( kx=xmin; kx<xmax; kx++ ) {
					j = (jy + xx + kx - hk[0])*c;
					block[nblock++] = (*this)[j];
				}
			}
		}
		imed = nblock*rank;
		partition(block, nblock, imed);
		nudata[i] = block[imed];
	}

	return 0;
}

/**
@brief 	Applies a median filter to an image.
@param 	kernel_size	length of kernel edge (typically 3).
@param	rank		which value in the kernel to retain.
@return int 		0, <0 if error.

	A kernel of a given size is passed over the image and the median
	value within the kernel assigned to the central voxel.

**/
int			Bimage::filter_rank(long kernel_size, double rank)
{
	if ( compoundtype != TSimple ) {
		cerr << "Error: The rank filter only operates on single value images!" << endl << endl;
		return -1;
	}
	
	if ( kernel_size%2 == 0 ) kernel_size += 1;
	
//	long			i, j, jz, jy, nn, xx, yy, zz, kx, ky, kz;
	Vector3<long>	blocksize(kernel_size, kernel_size, kernel_size);
	blocksize = blocksize.min(size());
	long			blocktotalsize = (long) blocksize.volume();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a rank filter:" << endl;
		cout << "Rank:                           " << (long)(rank*blocktotalsize) << endl;
		cout << "Kernel size:                    " << blocksize << " (" << blocksize.volume() << ")" << endl << endl;
	}

	long			i, chunk_size(get_chunk_size(datasize, c));
	float*			nudata = new float[datasize];

#ifdef HAVE_GCD
	dispatch_apply((datasize - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		filter_rank_chunk(kernel_size, rank, nudata, i*chunk_size, chunk_size);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<datasize; i+=chunk_size )
		filter_rank_chunk(kernel_size, rank, nudata, i, chunk_size);
#endif
    
	for ( i=0; i<datasize; i++ ) set(i, nudata[i]);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::filter_rank: Done" << endl;
	
	delete[] nudata;
	
	statistics();
	
	return 0;
}

/**
@brief 	Finds the peaks in an image.
@param 	kernel_size	length of kernel edge (typically 3).
@return int 		0, <0 if error.

	A kernel of a given size is passed over the image and if the central
	voxel is the maximum, it is kept, otherwise it is set to background.

**/
Bimage*		Bimage::filter_peak(long kernel_size)
{
	if ( compoundtype != TSimple ) {
		cerr << "Error: The average filter only operates on single value images!" << endl << endl;
		return NULL;
	}
	
	if ( kernel_size%2 == 0 ) kernel_size += 1;
	if ( kernel_size < 3 ) kernel_size = 3;
	
	long				i, j, nn, xx, yy, zz, kx, ky, kz;
	long				zlo, zhi, ylo, yhi, xlo, xhi;
	Vector3<long>		blocksize(kernel_size, kernel_size, kernel_size);
	blocksize = blocksize.min(size());
	Vector3<long>		h(blocksize/2);
	long				blocktotalsize = (long) blocksize.volume();
	double				bkg, v, v2;
	
	Bimage*				ppeak = copy();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying a peak filter:" << endl;
		cout << "Block size:                     " << blocksize << " (" << blocktotalsize << ")" << endl << endl;
	}

	// Do the peak filtering
	for ( i=nn=0; nn<n; nn++ ) {
		bkg = background(nn);
		for ( zz=0; zz<z; zz++ ) {
			zlo = (zz > h[2])? zz - h[2]: 0;
			zhi = zz + h[2];
			if ( zhi >= z ) zhi = z - 1;
			for ( yy=0; yy<y; yy++ ) {
				ylo = (yy > h[1])? yy - h[1]: 0;
				yhi = yy + h[1];
				if ( yhi >= y ) yhi = y - 1;
				for ( xx=0; xx<x; xx++, i++ ) {
					xlo = (xx > h[0])? xx - h[0]: 0;
					xhi = xx + h[0];
					if ( xhi >= x ) xhi = x - 1;
					v = (*this)[i];
					for ( kz=zlo; kz<=zhi; kz++ ) {
						for ( ky=ylo; ky<=yhi; ky++ ) {
							for ( kx=xlo; kx<=xhi; kx++ ) {
								j = index(kx, ky, kz, nn);
								v2 = (*this)[j];
								if ( v < v2 ) v = bkg;
							}
						}
					}
					ppeak->set(i, v);
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::filter_peak: Done" << endl;
	
	ppeak->statistics();
	
	return ppeak;
}

/**
@brief 	Calculates an average within a periodic frame.
@param 	period		size of periodic frame.
@return Bimage*		new image.
**/
Bimage*		Bimage::periodic_averaging(Vector3<double> period)
{
	change_type(Float);
	
	period = period.min(size());

	Bimage*			pnu = copy_header();
	pnu->size(period);
	pnu->page_size(period);
	pnu->data_alloc_and_clear();
	vector<float>	w(pnu->image_size());
	
	long			i, j, k, xo, yo, zo, nn, xx, yy, zz, cc;
	
	if ( verbose ) {
		cout << "Calculating a periodic average:" << endl;
		cout << "Period:                         " << period << endl << endl;
	}
	
	for ( i=nn=0; nn<n; nn++ ) {	
		for ( k=0; k<pnu->image_size(); k++ ) w[k] = 0;
		for ( zo=zz=0; zo<z; zo++, zz++ ) {
			if ( zz >= pnu->sizeZ() ) zz = 0;
			for ( yo=yy=0; yo<y; yo++, yy++ ) {
				if ( yy >= pnu->sizeY() ) yy = 0;
				for ( xo=xx=0; xo<x; xo++, xx++ ) {
					if ( xx >= pnu->sizeX() ) xx = 0;
					k = pnu->index(xx, yy, zz);
					w[k] += 1;
					j = pnu->index(0, xx, yy, zz, nn);
					for ( cc=0; cc<c; cc++, j++, i++ )
						pnu->add(j, (*this)[i]);
				}
			}
		}
		for ( k=0, j=nn*c*pnu->image_size(); k<pnu->image_size(); k++ )
			if ( w[k] ) for ( cc=0; cc<c; cc++, j++ )
				pnu->set(j, (*pnu)[j]/w[k]);
	}
	
	pnu->origin(Vector3<long>(period/2));
	
	return pnu;
}

double		Bimage::aniso_voxel(Bimage* pg, long i, long ksize, double w)
{
	long			xx, yy, zz, iz, iy, ix, jx, jy, jz;
	long			xy(x*y);
	double			dx2, dy2, dz2, d, dk2(ksize*ksize);
	double			v(0), wone, wsum(0), gm;
	Vector3<double>	g = pg->vector3(i);
	Vector3<long>	coor = coordinates(i);

	zz = (coor[2]>ksize)? -ksize: -coor[2];
	for ( iz=coor[2]+zz; iz<z && zz<=ksize; ++zz, ++iz ) {
		dz2 = zz*zz;
		jz = iz*xy;
		yy = (coor[1]>ksize)? -ksize: -coor[1];
		for ( iy=coor[1]+yy; iy<y && yy<=ksize; ++yy, ++iy ) {
			dy2 = yy*yy;
			jy = jz + iy*x;
			xx = (coor[0]>ksize)? -ksize: -coor[0];
			for ( ix=coor[0]+xx; ix<x && xx<=ksize; ++xx, ++ix ) {
				dx2 = xx*xx;
				d = dx2 + dy2 + dz2;
				if ( d <= dk2 ) {
					jx = jy + ix;
					gm = sqrt(fabs(g[0]*xx + g[1]*yy + g[2]*zz));
//					wone = 1/(1 + w * fabs(g[0]*xx + g[1]*yy + g[2]*zz));	// Scalar magnitude
//					wone = 1/(1 + w * gm);	// Scalar magnitude
					wone = 1/(1 + exp(w*gm*(2*sqrt(d)-1)));	// Scalar magnitude
					v += (*this)[jx] * wone;
					wsum += wone;
				}
			}
		}
	}
	v /= wsum;
	
	return v;
}

/**
@brief 	Calculates an anisotropic average within a kernel based on the local gradient.
@param 	ksize		kernel radius.
@param 	w			gradient weight (0 reverts to isotropic).
@return Bimage*		new image.
**/
Bimage*		Bimage::aniso_average(long ksize, double w)
{
	Bimage*			pg = gradient();
	
	Bimage*			pd = copy();
	pd->change_type(Float);
	
#ifdef HAVE_GCD
	dispatch_apply(image_size(), dispatch_get_global_queue(0, 0), ^(size_t i){
		pd->set(i, aniso_voxel(pg, i, ksize, w));
	});
#else
#pragma omp parallel for
	for ( long i=0; i<image_size(); ++i ) {
		pd->set(i, aniso_voxel(pg, i, ksize, w));
	}
#endif

	delete pg;
	
	pd->statistics();

	return pd;
}

/**
@brief 	Calculates an anisotropic average within a kernel based on the local gradient.
@param 	*p			image to compare with.
@return int			0.
**/
int			Bimage::filter_by_difference(Bimage* p)
{
	long			i, exponent(2);
	double			v, d, f(1), b(avg);
	
	for ( i=0; i<image_size(); ++i ) {
		v = (*this)[i];
		d = fabs(v - (*p)[i]);
		if ( exponent == 2 ) d *= d;
		else if ( exponent > 2 ) d = pow(d, exponent);
		f = 1/(1 + d);
		set(i, f*v + (1-f)*b);
	}
	
	statistics();
	
	return 0;
}

