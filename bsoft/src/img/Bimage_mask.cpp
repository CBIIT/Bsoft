/**
@file	Bimage_mask.cpp
@brief	Methods for binary mask creation and manipulation.
@author Samuel Payne and Bernard Heymann
@date	Created: 20010710
@date	Modified: 20200821 (BH)
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@author	Samuel Payne
@brief	Masks an image.
@param 	*pmask		binary mask.
@param 	fill		value to use where mask is zero.
@return long		0.

	If the mask value for that pixel is 0, the image pixel is changed
	to the fill value.  Otherwise it is left alone.

**/
long 		Bimage::mask(Bimage* pmask, double fill)
{
    long		i, j, cc, ds(x*y*z*n);
	
	for ( i=j=0; i<ds; i++ )
		for ( cc=0; cc<c; cc++, j++ )
			if ( (*pmask)[i] < 0.5 ) set(j, fill);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask: datasize=" << ds << " type=" << datatype << endl;
    
	statistics();

	return 0;
}

/**
@brief 	Change the image to a mask.
@param	threshold	set mask above this value.
@return long			number of mask voxels.

	The input image is effectively thresholded and a mask with 0 and 1 generated.
	The new data type is unsigned char/byte.

**/
long		Bimage::to_mask(double threshold)
{
	if ( compoundtype != TSimple ) {
		cerr << "Error: Conversion from compound types to mask not supported!" << endl;
		return 0;
	}
	
	if ( datatype == UCharacter && min == 0 && max == 1 ) return 0;
	
	long				j, nm(0);
	unsigned char*		mask = new unsigned char[datasize];
	
	for ( j=0; j<datasize; j++ )
		if ( (*this)[j] > threshold ) { mask[j] = 1; nm++; }
		else mask[j] = 0;
	
	data_type(UCharacter);
	compoundtype = TSimple;
	c = 1;
	min = 0;
	max = 1;
	
	data_assign(mask);
	
	statistics();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::to_mask: max=" << max << endl;
	
	return nm;
}

/**
@author	Samuel Payne & Bernard Heymann
@brief 	Generates a mask based on an image at a given threshold.
@param 	threshold	gray scale level.
@return Bimage*		unsigned char mask.

	A binary mask is generated using the given threshold.
	Image statistics are recalculated.

**/
Bimage*		Bimage::mask_by_threshold(double threshold)
{
	if ( datatype == Bit ) return copy();
	
    long				i;
	
	Bimage* 			pmask = copy_header();
	pmask->compound_type(TSimple);
	pmask->data_type(UCharacter);
	pmask->data_alloc();
	
	if ( verbose & VERB_LABEL )
	    cout << "Thresholding to:                " << threshold << endl << endl;
    
	for ( i=0; i<datasize; i++ ) {
		if ( (*this)[i] >= threshold ) pmask->set(i, 1);
		else pmask->set(i, 0);
	}

	pmask->statistics();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_by_threshold: pmask->max=" << pmask->maximum() << endl;
	
	return pmask;
}

/**
@brief 	Generates a mask based on an image at given thresholds.
@param 	threshold	array of gray scale levels.
@return Bimage*		unsigned char mask.

	A multi-level mask is generated using the given threshold values.
	Image statistics are recalculated.

**/
Bimage*		Bimage::mask_by_thresholds(vector<double> threshold)
{
	if ( datatype == Bit ) return copy();
	
	double				v;
	
	Bimage* 			pmask = copy_header();
	pmask->compound_type(TSimple);
	pmask->data_type(UCharacter);
	pmask->data_alloc_and_clear();
	
	if ( verbose & VERB_LABEL )
	    cout << "Thresholding to:                " << threshold << endl << endl;
    
	for ( long i=0; i<datasize; i++ ) {
		v = (*this)[i];
		for ( size_t j=0; j<threshold.size(); j++ )
			if ( v >= threshold[j] ) pmask->set(i, j+1);
	}

	pmask->statistics();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_by_threshold: pmask->max=" << pmask->maximum() << endl;
	
	return pmask;
}

/**
@brief 	Generates a mask based on a conditional hierarchy of thresholds.
@param 	threshold	array of gray scale levels.
@return Bimage*		unsigned char mask.

	A multi-level mask is generated using the given threshold values.
	Image statistics are recalculated.

**/
Bimage*		Bimage::mask_by_conditional_thresholds(vector<double> threshold)
{
	if ( datatype == Bit ) return copy();
	
	long			t, i, j, ts(threshold.size());

	Bimage*			pmask = mask_by_threshold(threshold[0]);
	
	if ( ts < 2 ) return pmask;

	pmask->multiply(ts);

    for ( t=1; t<ts; ++t ) {
		if ( verbose & VERB_LABEL )
			cout << "Conditional thresholding to:    " << threshold[t] << endl << endl;
		for ( i=0; i<datasize; ++i ) if ( !(*pmask)[i] ) {
			if ( (*this)[i] >= threshold[t] ) {
				j = kernel_max(i, 1);
				if ( (*pmask)[j] ) pmask->set(i, ts-t);
			}
		}
	}

	pmask->statistics();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_by_conditional_thresholds: pmask->max=" << pmask->maximum() << endl;
	
	return pmask;
}

/**
@brief 	Calculates statistics for a mask.
@return long			voxels in positive levels.

	The mask can be any type, but the regions in the mask are rounded to
	the nearest integer for counting statistics.

**/
long		Bimage::mask_stats()
{
//	if ( max < 1 ) statistics();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_stats: max=" << max << endl;
	
	long			i, j, neg, sum, psum;
	long			nlev = (long) (max + 1.5);
	if ( nlev < 2 ) nlev = 2;
	
	double			value, volsum;
	double			voxel_volume(sampling(0).volume());
	vector<long>	h(nlev, 0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_stats: nlev=" << nlev << endl;
    
	if ( verbose )
		cout << "Calculating mask statistics:" << endl;
	
	for ( i=neg=0; i<datasize; i++ ) {
		value = (*this)[i];
		if ( value < 0 ) neg++;
		else {
			j = (long) (value + 0.5);
			h[j]++;
		}
	}

	for ( i=1, psum=0; i<nlev; i++ ) psum += h[i];
	sum = psum + neg + h[0];
	volsum = psum*voxel_volume;
	
	cout << "Region\tCount\t%\tVolume(A3)" << endl;
	if ( neg ) cout << "-\t" << neg << tab << neg*100.0/sum << tab << neg*voxel_volume << endl;
	for ( i=0; i<nlev; i++ ) if ( h[i] )
		cout << i << tab << h[i] << tab << h[i]*100.0/sum << tab << h[i]*voxel_volume << endl;
	cout << endl;
	
	cout << "Count sum:                      " << psum << endl;
	cout << "Volume sum:                     " << volsum << " A3" << endl << endl;
	
	return psum;
}

/**
@brief 	Inverts a mask.
@return long		number of mask voxels.

	The mask is assumed to be in range [0,1].
		new_mask = 1 - old_mask.
	Image statistics of the mask are recalculated.

**/
long		Bimage::mask_invert()
{
	to_mask();
	
	long				i, nm(0);
	
	if ( verbose & VERB_PROCESS )
		cout << "Inverting mask" << endl << endl;
	
	for ( i=0; i<datasize; i++ ) {
		if ( (*this)[i] > 0.5 ) set(i, 0);
		else { set(i, 1); nm++; }
	}
	
	statistics();
	
	return nm;
}

/**
@brief 	Combines two masks with different operations.
@param 	*p			second mask.
@param 	operation		combining operation.
@return Bimage*				unsigned char mask.

	The input images are assumed to be a masks of 0's and 1's and are
	modified according to the operation requested:
		0   val1 = val2
		1   val1 = val1 and val2
		2   val1 = val1 or  val2
		3   val1 = val1 xor val2
	where val1 and val2 are the data in the two images. 
	Image statistics of the first image are recalculated.

**/
long		Bimage::mask_combine(Bimage* p, int operation)
{
	to_mask();
	p->to_mask();
	
	long			i;
	
//	if ( max < 1 || p->maximum() < 1 ) operation = 2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_combine: op=" << operation << endl << endl;
	
	switch ( operation ) {
		case 1:
			for ( i=0; i<datasize; i++ ) d.uc[i] = (d.uc[i] & p->d.uc[i]);
			break;
		case 2:
			for ( i=0; i<datasize; i++ ) d.uc[i] = (d.uc[i] | p->d.uc[i]);
			break;
		case 3:
			for ( i=0; i<datasize; i++ ) d.uc[i] = (d.uc[i] ^ p->d.uc[i]);
			break;
		default:
			for ( i=0; i<datasize; i++ ) d.uc[i] = p->d.uc[i];
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Zeroes everything outside the mask and return the excised feature.
@param 	*pmask			mask.
@return Bimage*			excised feature.

	The mask must be the same size as the image.
**/
Bimage*		Bimage::mask_extract(Bimage* pmask)
{
	long		i, j, cc, ds(x*y*z*n);
	Bimage*		pnu = copy();
	
	for ( i=j=0; i<ds; i++ )
		for ( cc=0; cc<c; cc++, j++ )
			if ( (*pmask)[i] < 0.5 )
				pnu->set(j, 0);
	
	pnu->statistics();
	
	return pnu;
}

int			Bimage::max_in_kernel(long ksize)
{
	long				i, j, nn, xx, yy, zz, kx, ky, kz;
	double				v, m;
	
	Vector3<long>		k(ksize, ksize, ksize);
	k = k.min(size());

	Vector3<long>		kh(k/2);
	Vector3<long> 	lo, hi;
	
	Bimage*				p = new Bimage(datatype, compoundtype, size(), n);
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			lo[2] = (zz > kh[2])? zz - kh[2]: 0;
			hi[2] = (zz+kh[2] < z)? zz + kh[2]: z - 1;
			for ( yy=0; yy<y; yy++ ) {
				lo[1] = (yy > kh[1])? yy - kh[1]: 0;
				hi[1] = (yy+kh[1] < y)? yy + kh[1]: y - 1;
				for ( xx=0; xx<x; xx++, i++ ) {
					lo[0] = (xx > kh[0])? xx - kh[0]: 0;
					hi[0] = (xx+kh[0] < x)? xx + kh[0]: x - 1;
					m = (*this)[i];
					for ( kz=lo[2]; kz<=hi[2]; kz++ ) {
						for ( ky=lo[1]; ky<=hi[1]; ky++ ) {
							for ( kx=lo[0]; kx<=hi[0]; kx++ ) {
								j = index(kx, ky, kz, nn);
								v = (*this)[j];
								if ( m < v ) m = v;
							}
						}
					}
					p->set(i, m);
				}
			}
		}
	}
	
	for ( i=0; i<datasize; i++ ) set(i, (*p)[i]);
	
	delete p;
	
	statistics();
	
	return 0;
}

/**
@author	Samuel Payne & Bernard Heymann
@brief 	Dilates or erodes a binary mask.
@param 	dir			0=erode, 1=dilate.
@return long 		masked voxels.

	Traditional 3^dim kernel dilation.  Any pixel with a value of 1 turns 
	all of its neighbors to a value of 1.

**/
long		Bimage::mask_dilate_erode(unsigned char dir)
{
	to_mask();
	
	long				i, j, nn, xx, yy, zz, kx, ky, kz;
	
	Vector3<long>		k(3, 3, 3);
	k = k.min(size());

	Vector3<long>		kh(k/2);
	Vector3<long> 	lo, hi;
	
	unsigned char*		mask = new unsigned char[datasize];
	for ( i=0; i<datasize; i++ ) mask[i] = d.uc[i];
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			lo[2] = (zz > kh[2])? zz - kh[2]: 0;
			hi[2] = (zz+kh[2] < z)? zz + kh[2]: z - 1;
			for ( yy=0; yy<y; yy++ ) {
				lo[1] = (yy > kh[1])? yy - kh[1]: 0;
				hi[1] = (yy+kh[1] < y)? yy + kh[1]: y - 1;
				for ( xx=0; xx<x; xx++, i++ ) {
					lo[0] = (xx > kh[0])? xx - kh[0]: 0;
					hi[0] = (xx+kh[0] < x)? xx + kh[0]: x - 1;
					if ( mask[i] == dir ) {
						for ( kz=lo[2]; kz<=hi[2]; kz++ ) {
							for ( ky=lo[1]; ky<=hi[1]; ky++ ) {
								for ( kx=lo[0]; kx<=hi[0]; kx++ ) {
									j = index(kx, ky, kz, nn);
									set(j, dir);
								}
							}
						}
					}
				}
			}
		}
	}
	
	delete[] mask;
	
	statistics();
	
	return (long) (avg*datasize);
}

/**
@author	Samuel Payne & Bernard Heymann
@brief 	Dilates a binary mask.
@param 	times		the number of times to dilate the mask.
@return long 		masked voxels.

	Traditional 3^dim kernel dilation.  Any pixel with a value of 1 turns 
	all of its neighbors to a value of 1.

**/
long		Bimage::mask_dilate(long times)
{
	long		t, nv(0);
	
	for ( t=0; t<times; t++ )
		nv = mask_dilate_erode(1);
	
	return nv;
}

/**
@author	Samuel Payne & Bernard Heymann
@brief 	Erodes a binary mask.
@param 	times		the number of times to erode the mask
@return long 		masked voxels.

	Traditional 3^dim kernel erosion.  If all the neighboring pixels have
	value of 1, then that pixel is left at 1. Otherwise it is changed to 0.

**/
long		Bimage::mask_erode(long times)
{
	long		t, nv(0);
	
	for ( t=0; t<times; t++ )
		nv = mask_dilate_erode(0);
	
	return nv;
}

/**
@brief 	Opens a binary mask.
@param 	times		the number of times to erode and dilate the mask.
@return long 		masked voxels.

	Opening a mask is an erosion followed by a dilation.

**/
long 		Bimage::mask_open(int times)
{
	long		t, nv(0);
	
	for ( t=0; t<times; t++ ) {
		mask_erode(1);
		nv = mask_dilate(1);
	}
	
	return nv;
}

/**
@brief 	Closes a binary mask.
@param 	times		the number of times to dilate and erode the mask.
@return long 		masked voxels.

	Closing a mask is a dilation followed by an erosion.

**/
long 		Bimage::mask_close(int times)
{
	long		t, nv(0);
	
	for ( t=0; t<times; t++ ) {
		mask_dilate(times);
		nv = mask_erode(times);
	}
	
	return nv;
}

/**
@brief 	Fills an empty part of a mask indicated by the given voxel.
@param 	voxel		point from which to fill the mask.
@return long 		masked voxels.

**/
long 		Bimage::mask_fill(Vector3<long> voxel)
{
	return region_assign(this, index(voxel, 0), 1, 0.5, -1);
}

/**
@brief 	Generates a mask on one side of a plane.
@param 	origin		any point on plane.
@param 	normal		plane normal.
@return long		0.

	The given image is used to generate a binary mask of the same size.

**/
long		Bimage::mask_plane(Vector3<double> origin, Vector3<double> normal)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_plane: normal=" << normal << endl;
		
	normal.normalize();
	
	long   		i, nn, xx, yy, zz, total(0);
	Vector3<double>	point;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a planar mask:" << endl;
		cout << "Normal:                         " << normal << endl;
		cout << "Point on plane:                 " << origin << endl;
	}
	
    for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			point[2] = zz - origin[2];
			for ( yy=0; yy<y; yy++ ) {
				point[1] = yy - origin[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					point[0] = xx - origin[0];
					if ( normal.scalar(point) > 0 ) {
						set(i, 1);
						total++;
					} else {
						set(i, 0);
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << total << endl << endl;
	
	return total;
}

/**
@brief 	Generates a mask along an axis within a 2D image.
@param 	length		length along axis (pixels).
@param 	width		width perpendicular to axis (pixels).
@param 	rect_angle	angle from x-axis (radians).
@param 	wrap        flag to wrap around boundaries.
@return long		0.

	The input image is assumed to be a mask of 0's and 1's and is modified
	according to the operation requested.
	Only an axis in the xy plane is used, generating a rectangular mask
	of given length and width around the axis, with its center at the
	origins given in the image structure.

**/
long		Bimage::mask_rectangle(double length, double width,
				double rect_angle, int wrap)
{
//	change_type(UCharacter) ;
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating a mask along an axis" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Length and width:               " << length << " x " << width << " pixels" << endl;
		cout << "Axis angle:                     " << rect_angle*180.0/M_PI << " degrees" << endl;
	}
	
	long   		i, nn, xx, yy, zz, total(0);
	double			xo, yo, xr, yr;
	double			cos_ang(cos(rect_angle)), sin_ang(sin(rect_angle));
	double			half_length(length/2.0), half_width(width/2.0);
	
    for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				yo = yy - image[nn].origin()[1];
				if ( wrap ) if ( yo > (y - 1)/2 ) yo -= y;
				for ( xx=0; xx<x; xx++, i++ ) {
					xo = xx - image[nn].origin()[0];
					if ( wrap ) if ( xo > (x - 1)/2 ) xo -= x;
					xr = xo*cos_ang + yo*sin_ang;
					yr = -xo*sin_ang + yo*cos_ang;
					if ( xr >= -half_length && xr <= half_length && 
							yr >= -half_width && yr <= half_width ) {
						set(i, 1);
						total++;
					} else {
						set(i, 0);
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << total << endl << endl;
	
	return total;
}

/**
@brief 	Symetrizes a mask.
@param 	sym			point group symmetry.
@return long			maximum level index.

**/
long		Bimage::mask_symmetrize(Bsymmetry& sym)
{
	long			i, j, nn, xx, yy, zz, cc;
	
	View			ref = view_symmetry_reference(sym);
	Vector3<double>	vr(ref.vector3());
	Vector3<double>	v, vt, va;
	double			da, minda;
	
	long			nsym(0);
	vector<Matrix3>	m = sym.matrices();
	nsym = m.size();

	if ( verbose & VERB_PROCESS ) {
		cout << "Symmetrizing a mask:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Origin:                         " << image->origin() << endl << endl;
	}
	
	int*				mask = new int[data_size()];
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			v[2] = zz - image[nn].origin()[2];
			for ( yy=0; yy<y; yy++ ) {
				v[1] = yy - image[nn].origin()[1];
				for ( xx=0; xx<x; xx++ ) {
					v[0] = xx - image[nn].origin()[0];
					minda = 1e30;
					for ( j=0; j<nsym; j++ ) {
						vt = m[j]*v;
						da = vt.angle(vr);
						if ( minda > da ) {
							minda = da;
							va = vt;
						}
					}
					va += image[nn].origin();
					for ( cc=0; cc<c; cc++, i++ ) {
						mask[i] = get(nn, va, cc);
					}
				}
			}
		}
	}

	data_type(Integer);
	
	data_assign((unsigned char *) mask);

	statistics();
		
	return max;
}

/**
@brief 	Resizes a reciprocal space mask.
@param 	nusize			new mask size.
@return Vector3<double> scale.

	The mask is assumed to be centered at (0,0,0).
	The image is resized by inserting or removing
	rows or columns in the middle of the data set.

**/
Vector3<double>	Bimage::mask_fspace_resize(Vector3<long> nusize)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_fspace_resize: size=" << nusize << endl;

	Vector3<double>	scale(nusize[0]*1.0L/x, nusize[1]*1.0L/y, nusize[2]*1.0L/z);

	if ( nusize.volume() < 1 ) return scale;
    if ( nusize == size() ) return scale;
	
	if ( compoundtype != TSimple ) {
		cerr << "Error: Bimage::mask_fspace_resize: Only simple numbers supported!" << endl;
		return scale;
    }
	
    long   		i, j, nn, xx, yy, zz, iy, iz;
    long     	    xold, yold, zold, zerox, zeroy, zeroz;
	Vector3<long>	hold, h, d;
    long   		ds(n*nusize.volume());
	unsigned char*	nudata = new unsigned char[ds];
	
	hold = (size() - 1)/2;
	h = (nusize - 1)/2;
	d = nusize - size();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_fspace_resize: datasize=" << ds << endl;
	
	if ( verbose & VERB_FULL )
	    cout << "Changing mask to size:          " << nusize << " voxels" << endl << endl;

    for ( j=nn=0; nn<n; nn++ ) {
	    for ( zz=0; zz<nusize[2]; zz++ ) {
    		zold = zz;
			if ( zz > h[2] ) zold -= d[2];
			zeroz = 0;
			if ( zz > hold[2] && zz <= d[2] + hold[2] ) zeroz += 1;
			iz = (nn*z + zold)*y;
    		for ( yy=0; yy<nusize[1]; yy++ ) {
			 	yold = yy;
				if ( yy > h[1] ) yold -= d[1];
				zeroy = zeroz;
				if ( yy > hold[1] && yy <= d[1] + hold[1] ) zeroy += 1;
				iy = (iz + yold)*x;
				for ( xx=0; xx<nusize[0]; xx++, j++ ) {
					xold = xx;
					if ( xx > h[0] ) xold -= d[0];
					zerox = zeroy;
					if ( xx > hold[0] && xx <= d[0] + hold[0] ) zerox += 1;
					if ( zerox < 1 ) {
						i = iy + xold;
						nudata[j] = ((*this)[i] > 0.5);
					}
				}
	    	}
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::mask_fspace_resize: origin[" << n << "]=" << image[nn].origin() << endl;
	
    }
    
	data_type(UCharacter);
	sampling(sampling(0)/scale);
	size(nusize);
	page_size(nusize);

    data_assign((unsigned char *) nudata);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::mask_fspace_resize: sampling=" << sampling(0) << endl;
	
	return scale;
}

/**
@brief 	Sets up default frequency space bands to generate a mask.
@param 	res_lo			low resolution limit (angstrom).
@param 	res_hi			high resolution limit (angstrom).
@return vector<double>	array of band specifications.

	The band argument is a list of pairs of values, each pair indicating
	a resolution shell (in angstrom) and a flag indicating whether the 
	following shells should be:
		0: excluded - also the high resolution limit
		1: included in the FOM
		-1: included in the cross-validation FOM
	The high resolution limit for the mask is set in the image structure.

**/
vector<double>	Bimage::fspace_default_bands(double res_lo, double res_hi)
{
	vector<double>	band(6,0);
	
	band[0] = res_lo;
	band[1] = 1;
	band[2] = res_hi;
	band[3] = -1;
	band[4] = 0.9*res_hi;
	band[5] = 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG fspace_default_bands:";
		for ( auto it=band.begin(); it!=band.end(); ++it )
			cout << tab << *it;
		cout << endl;
	}
	
	return band;
}

/**
@brief 	Generate reciprocal space mask based on a specification of bands.
@param 	&band			array of number pairs, each a shell and value.
@return int				0.

	The band argument is a list of pairs of values, each pair indicating
	a resolution shell (in angstrom) and a flag indicating whether the 
	following shells should be:
		0: excluded - also the high resolution limit
		1: included in the FOM
		-1: included in the cross-validation FOM
	The high resolution limit for the mask is set in the image structure.

**/
int			Bimage::mask_fspace_banded(vector<double>& band)
{
	if ( band.size() < 2 || band[0] < 1 ) return -1;
	
	change_type(SCharacter);
	fill(0);
	
	long			i, j, xx, yy, zz, h[3] = {0,0,0};
	long			ix, iy, iz;
	double			sx2, sy2, sz2, s2;
	double			sband[100] = {1e-30, 0};
	Vector3<double>	scale(1.0/real_size());
	
//	cout << "scale = " << scale << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a reciprocal space banded mask:" << endl;
		cout << "Pixel size:                     " << sampling(0) << endl;
		cout << "Band\tValue" << endl;
	}
	
	// The reciprocal space specification starts with an excluded region, followed
	// by whatever the band specification is set at.
	for ( j=0; j<2*band.size() && band[j] > 0; j += 2 ) {
		if ( verbose & VERB_PROCESS )
			cout << band[j] << tab << band[j+1] << endl;
		sband[j+2] = 1.0/band[j];
		sband[j+2] *= sband[j+2];
		sband[j+3] = band[j+1];
	}
	sband[j+2] = 0;
	
	for ( i=zz=0; zz<z; zz++ ) {
		iz = zz;
		if ( zz > (z - 1)/2 ) iz -= z;
		sz2 = iz*scale[2];
		sz2 *= sz2;
		for ( yy=0; yy<y; yy++ ) {
			iy = yy;
			if ( yy > (y - 1)/2 ) iy -= y;
			sy2 = iy*scale[1];
			sy2 *= sy2;
			for ( xx=0; xx<x; xx++, i++ ) {
				ix = xx;
				if ( xx > (x -1 )/2 ) ix -= x;
				sx2 = ix*scale[0];
				sx2 *= sx2;
				s2 = sx2 + sy2 + sz2;
				for ( j=0; sband[j] > 0 && s2 > sband[j]; j += 2 ) ;
				if ( j > 0 ) j -= 2;
				set(i, sband[j+1]);
				if ( sband[j+1] < 0 ) h[0]++;
				else if ( sband[j+1] > 0 ) h[2]++;
				else h[1]++;
			}
		}
	}
	
	j = h[0] + h[1] + h[2];
	if ( verbose & VERB_PROCESS ) {
		cout << "Mask partitions:" << endl;
		for ( i=0; i<3; i++ )
			cout << i-1 << tab << h[i] << tab << h[i]*100.0/j << "%" << endl;
		cout << endl;
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Generates a mask with a missing wedge.
@param 	ori			image origin.
@param 	tilt_axis	tilt axis angle.
@param 	tilt_neg	negative tilt angle.
@param 	tilt_pos	positive tilt angle.
@param 	resolution	high resolution limit.
@return long		number of mask voxels.

	The given image is used to generate a binary mask of the same size.
	If the input origin is {0,0,0}, the origin in the image is used as 
	the center for the missing part.

**/
long		Bimage::mask_missing_wedge(Vector3<double> ori, double tilt_axis,
				double tilt_neg, double tilt_pos, double resolution)
{
	change_type(UCharacter);
	clear();
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating a mask with a missing wedge" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Tilt axis angle:                " << tilt_axis*180.0/M_PI << " degrees" << endl;
		cout << "Tilt angle range:               " << tilt_neg*180.0/M_PI << " - " << tilt_pos*180.0/M_PI << " degrees" << endl;
		cout << "Origin:                         " << ori << endl;
	}
	
	long   			i, nn, xx, yy, zz, cc, nv(0);
	double			d, dx, dy, dz, dx2, dy2, dz2, dn, dp;
	double			ax(cos(tilt_axis));
	double			ay(sin(tilt_axis));
	double			an(tan(fabs(tilt_neg) - M_PI_2));
	double			ap(tan(M_PI_2 - fabs(tilt_pos)));

	double			max_rad(1);
	if ( resolution > image->sampling()[0] )
		max_rad = image->sampling()[0]/resolution;
	double			max_rad_sq(max_rad*max_rad);	// Maximum radius squared on map scale

    for ( i=nn=0; nn<n; nn++ ) {
//		image[nn].origin(origin);
		image[nn].origin(ori);
		for ( zz=0; zz<z; zz++ ) {
			dz = ((double)zz - ori[2])/z;
			if ( dz < -0.5 ) dz += 1;
			if ( dz >= 0.5 ) dz -= 1;
			dz2 = dz*dz;
			dn = dz*an;
			dp = dz*ap;
			if ( dz < 0 ) swap(dn, dp);
			for ( yy=0; yy<y; yy++ ) {
				dy = ((double)yy - ori[1])/y;
				if ( dy < -0.5 ) dy += 1;
				if ( dy >= 0.5 ) dy -= 1;
				dy2 = dy*dy;
				for ( xx=0; xx<x; xx++, i+=c ) {
					dx = ((double)xx - ori[0])/x;
					if ( dx < -0.5 ) dx += 1;
					if ( dx >= 0.5 ) dx -= 1;
					dx2 = dx*dx;
					if ( dx2 + dy2 + dz2 < max_rad_sq ) {
						d = dx*ay - dy*ax;
						if ( d < dn || d > dp ) {
							for ( cc=0; cc<c; cc++ ) set(i+cc, 1);
							nv += c;
						}
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << nv << endl << endl;
	
	return nv;
}

/**
@brief 	Generates a mask with a missing pyramid.
@param 	ori			image origin.
@param 	tilt_axis1	tilt axis angle 1.
@param 	tilt_axis2	tilt axis angle 2.
@param 	tilt_neg1	negative tilt angle for axis 1.
@param 	tilt_pos1	positive tilt angle for axis 1.
@param 	tilt_neg2	negative tilt angle for axis 2.
@param 	tilt_pos2	positive tilt angle for axis 2.
@param 	resolution	high resolution limit.
@return long		number of mask voxels.

	The given image is used to generate a binary mask of the same size.
	If the input origin is {0,0,0}, the origin in the image is used as 
	the center for the missing part.

**/
long		Bimage::mask_missing_pyramid(Vector3<double> ori, double tilt_axis1,
				double tilt_axis2, double tilt_neg1, double tilt_pos1,
				double tilt_neg2, double tilt_pos2, double resolution)
{
	change_type(UCharacter);
	clear();
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating a mask with a missing pyramid" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Tilt axis angle 1:              " << tilt_axis1*180.0/M_PI << " degrees" << endl;
		cout << "Tilt angle range 1:             " << tilt_neg1*180.0/M_PI << " - " << tilt_pos1*180.0/M_PI << " degrees" << endl;
		cout << "Tilt axis angle 2:              " << tilt_axis2*180.0/M_PI << " degrees" << endl;
		cout << "Tilt angle range 2:             " << tilt_neg2*180.0/M_PI << " - " << tilt_pos2*180.0/M_PI << " degrees" << endl;
		cout << "Origin:                         " << ori << endl;
	}
	
	long   			i, nn, xx, yy, zz, cc, nv(0);
	double			d1, d2, dx, dy, dz, dx2, dy2, dz2, dn1, dp1, dn2, dp2;
	double			ax1(cos(tilt_axis1));
	double			ay1(sin(tilt_axis1));
	double			ax2(cos(tilt_axis2));
	double			ay2(sin(tilt_axis2));
	double			an1(tan(fabs(tilt_neg1) - M_PI_2));
	double			ap1(tan(M_PI_2 - fabs(tilt_pos1)));
	double			an2(tan(fabs(tilt_neg2) - M_PI_2));
	double			ap2(tan(M_PI_2 - fabs(tilt_pos2)));

	double			max_rad(1);
	if ( resolution > image->sampling()[0] ) 
		max_rad = image->sampling()[0]/resolution;
	double			max_rad_sq(max_rad*max_rad);	// Maximum radius squared on map scale
	
    for ( i=nn=0; nn<n; nn++ ) {
		image[nn].origin(ori);
		for ( zz=0; zz<z; zz++ ) {
			dz = ((double)zz - ori[2])/z;
			if ( dz < -0.5 ) dz += 1;
			if ( dz >= 0.5 ) dz -= 1;
			dz2 = dz*dz;
			dn1 = dz*an1;
			dp1 = dz*ap1;
			if ( dn1 > dp1 ) swap(dn1, dp1);
			dn2 = dz*an2;
			dp2 = dz*ap2;
			if ( dn2 > dp2 ) swap(dn2, dp2);
			for ( yy=0; yy<y; yy++ ) {
				dy = ((double)yy - ori[1])/y;
				if ( dy < -0.5 ) dy += 1;
				if ( dy >= 0.5 ) dy -= 1;
				dy2 = dy*dy;
				for ( xx=0; xx<x; xx++, i+=c ) {
					dx = ((double)xx - ori[0])/x;
					if ( dx < -0.5 ) dx += 1;
					if ( dx >= 0.5 ) dx -= 1;
					dx2 = dx*dx;
					if ( dx2 + dy2 + dz2 < max_rad_sq ) {
						d1 = dx*ay1 - dy*ax1;
						d2 = dx*ay2 - dy*ax2;
						if ( d1 < dn1 || d1 > dp1 || d2 < dn2 || d2 > dp2 ) {
							for ( cc=0; cc<c; cc++ ) set(i+cc, 1);
							nv += c;
						}
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << nv << endl << endl;
	
	return nv;
}

/**
@brief 	Generates a mask with a missing cone.
@param 	ori			image origin.
@param 	mis_ang		cone angle from xy plane.
@param 	resolution	high resolution limit.
@return long		number of mask voxels.

	The given image is used to generate a binary mask of the same size.
	If the input origin is {0,0,0}, the origin in the image is used as 
	the center for the missing part.

**/
long		Bimage::mask_missing_cone(Vector3<double> ori,
						double mis_ang, double resolution)
{
	change_type(UCharacter);
	clear();
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating a mask with a missing cone" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Cone angle:                     " << mis_ang*180.0/M_PI << " degrees" << endl;
		cout << "Origin:                         " << ori << endl;
	}
	
	long   			i, nn, xx, yy, zz, cc, nv(0);
	double			dx, dy, dz, dx2, dy2, dz2, d2;
	
	double			max_rad(1);
	if ( resolution > image->sampling()[0] ) 
		max_rad = image->sampling()[0]/resolution;
	double			max_rad_sq(max_rad*max_rad);	// Maximum radius squared on map scale

    for ( i=nn=0; nn<n; nn++ ) {
		image[nn].origin(ori);
		for ( zz=0; zz<z; zz++ ) {
			dz = ((double)zz - ori[2])/z;
			if ( dz < -0.5 ) dz += 1;
			if ( dz >= 0.5 ) dz -= 1;
			dz2 = dz*dz;
			dz = fabs(dz);
			d2 = dz*tan(M_PI_2 - mis_ang);
			d2 *= d2;
			for ( yy=0; yy<y; yy++ ) {
				dy = ((double)yy - ori[1])/y;
				if ( dy < -0.5 ) dy += 1;
				if ( dy >= 0.5 ) dy -= 1;
				dy2 = dy*dy;
				for ( xx=0; xx<x; xx++, i+=c ) {
					dx = ((double)xx - ori[0])/x;
					if ( dx < -0.5 ) dx += 1;
					if ( dx >= 0.5 ) dx -= 1;
					dx2 = dx*dx;
					if ( dx2 + dy2 + dz2 < max_rad_sq ) {
						if ( dx2 + dy2 > d2 ) {
							for ( cc=0; cc<c; cc++ ) set(i+cc, 1);
							nv += c;
						}
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << nv << endl << endl;
	
	return nv;
}

//#include "rwimg.h"

/**
@brief 	Generates a missing mask from an example image.
@param 	ori			image origin.
@param 	resolution	high resolution limit.
@param 	mis_type	missing region type: wedge.
@return long		number of mask voxels.

	The given image is used to generate a binary mask of the same size.
	If the input origin is {0,0,0}, the origin in the image is used as 
	the center for the missing part.

**/
long		Bimage::mask_missing_find(Vector3<double> ori, double resolution, Bstring& mis_type)
{
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Finding missing regions and creating a missing mask" << endl;
	
	power_spectrum(0);
	
	long			i, xx, yy, nz(14), hx(x/2), hy(y/2), nr(x), na(180), ia, ir, nv(0);
	double			v, a, sa(M_PI*1.0/na), vmin, vmax, ax, amin(0), amax(0), r, f;
	Vector3<long>	st(0,0,nz), sz(x,y,1);
	Vector3<double>	d;
	
	Bimage*			pex = extract(0, st, sz);

	Bimage*			ps = new Bimage(Float, TSimple, nr, na, 1, 1);
	
	for ( i=yy=0; yy<y; ++yy ) {
		d[1] = (yy < hy)? yy: yy - y;
		for ( xx=0; xx<x; ++xx, ++i ) {
			d[0] = (xx < hx)? xx: xx - x;
			v = (*pex)[i];
			for ( ia=0; ia<na; ++ia ) {
				a = sa*ia;
				r = cos(a)*d[1] - sin(a)*d[0];
//				r = cos(a)*d[0] + sin(a)*d[1];
				if ( r < 0 ) r += nr;
				ir = (long)r;
				f = r - ir;
				if ( ir < 0 ) ir += nr;
				if ( ir >= 0 && ir < nr )
					ps->add(ia*nr + ir, (1-f)*v);
				ir++;
				if ( ir >= nr ) ir -= nr;
				if ( ir >= 0 && ir < nr )
					ps->add(ia*nr + ir, f*v);
			}
		}
	}
	
//	write_img("ps.pif", ps, 0);

	ax = 0;
	vmin = 1e30;
	for ( ia=0; ia<na; ++ia ) {
		v = 0;
		for ( ir=0, i=ia*nr; ir<nz; ++ir, ++i )
			v += (*ps)[i];
		if ( vmin > v ) {
			vmin = v;
			ax = ia;
		}
//		cout << ia << tab << v << endl;
	}
	ia = ax;
	ax = M_PI*ax/na;
	i = ia*nr;
	vmax = 0;
	vmin = 0;
	v = (*ps)[i];
	for ( ir=0; ir<nr; ++ir, ++i ) {
		f = (*ps)[i] - v;
//		cout << ir << tab << f << endl;
		if ( vmin > f ) {
			vmin = f;
			amin = ir;
		}
		if ( vmax < f ) {
			vmax = f;
			amax = ir;
		}
		v = (*ps)[i];
	}
	amin = -atan2(-nz, nr - amin);
	amax = -atan2(nz, amax);
	
	if ( verbose & VERB_FULL ) {
		cout << "Axis angle:                     " << ax*180.0/M_PI << endl;
		cout << "Minimum tilt angle:             " << amax*180.0/M_PI << endl;
		cout << "Maximum tilt angle:             " << amin*180.0/M_PI << endl;
	}
	
	delete pex;
	delete ps;

//	change_type(UCharacter);
//	clear();
	mask_missing_wedge(ori, ax, amax, amin, resolution);
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of pixels in mask:       " << nv << endl << endl;
	
	return nv;
}

/**
@brief 	Packs a 2D mask into a 3D reciprocal space volume.  
@param 	mat				affine orientation matrix.
@param 	hi_res			high resolution limit.
@param 	scale			scale.
@return	long			0.

	The rotation matrix is used to determine the plane in reciprocal space
	to set as one.
	Both the high resolution limit and the scale must correspond to the 
	associated reconstruction.

**/
long			Bimage::mask_pack_plane(Matrix3 mat, double hi_res, double scale)
{
	check_resolution(hi_res);
	
	long 			j;
	double 			w, d2;
	Vector3<long>	coor;
	Vector3<double>	m, dist, iv;
	
	Vector3<double>	invsize(1.0/x, 1.0/y, 1.0/z);
	Vector3<double> vscale(1/scale, 1/scale, z*1.0/(scale*x));
//	mat = mat.transpose();
	mat = vscale * mat;	// Matrix includes scaling going from image to map
	
	if ( verbose & VERB_FULL ) {
		cout << "Packing a mask into reciprocal space up to " << hi_res << " A resolution" << endl;
		cout << "Transformation matrix:" << endl;
		cout << mat << endl;
		cout << endl;
	}
	
	double			max_rad(image->sampling()[0]/hi_res);
	if ( max_rad > 0.5 ) max_rad = 0.5;
	double			max_rad_sq = max_rad*max_rad;	// Maximum radius squared on map scale
//	cout << "max_rad_sq = " << max_rad_sq << endl;
	double			ymin = floor(-max_rad*x);
	double			ymax = -ymin;
	double			xmin = floor(-max_rad*x);
	double			xmax = -xmin;

	for ( iv[1]=ymin; iv[1]<=ymax; iv[1]+=1 ) {
		for ( iv[0]=xmin; iv[0]<=xmax; iv[0]+=1 ) {
			m = mat * iv;
			d2 = (m * invsize).length2();
			if ( d2 < max_rad_sq ) {
				coor = Vector3<long>((long) floor(m[0] + 0.5),
					(long) floor(m[1] + 0.5),
					(long) floor(m[2] + 0.5));
				dist = m - coor;
				w = 1 - dist.length();
				if ( w > 1e-6 ) {
					j = index_wrap(coor);
					set(j, 1);
				}
			}
		}
	}

	return 0;
}


/**
@brief 	Calculates a threshold for a local variance image.
@param 	lowvar			low variance threshold.
@return double 			threshold.

	A threshold is calculated to define the foreground and background.
	Excluded regions are defined by a low variance parameter.

**/
double	 	Bimage::variance_threshold(double lowvar)
{
	double			varmax(max);
	if ( varmax <= lowvar ) varmax = 1;

	long			i;
	long			nfg(0), nbg(0);
	double			snr(0), t(0), v, vfg(varmax), vbg(0);
	double			dvfg(varmax), dvbg(varmax), dt(lowvar), sumfg, sumbg;
	
	while ( fabs(dvfg) > dt || fabs(dvbg) > dt ) {
		nfg = nbg = 0;
		sumfg = sumbg = 0;
		t = (vbg + vfg)/2;
		for ( i=0; i<datasize; i++ ) {
			v = (*this)[i];
			if ( v > dt ) {
				if ( v > t ) {
					sumfg += v;
					nfg++;
				} else {
					sumbg += v;
					nbg++;
				}
			}
		}
		if ( nfg ) sumfg /= nfg;
		if ( nbg ) sumbg /= nbg;
		dvfg = sumfg - vfg;
		dvbg = sumbg - vbg;
		vfg = sumfg;
		vbg = sumbg;
	}
	
	snr = vfg/vbg - 1;
	image->FOM(snr);
	
	if ( verbose ) {
		cout << "Background variance:            " << vbg << " (" << nbg << " voxels, " << nbg*100.0/datasize << " %)" << endl;
		cout << "Foreground variance:            " << vfg << " (" << nfg << " voxels, " << nfg*100.0/datasize << " %)" << endl;
		cout << "Signal-to-noise ratio:          " << snr << endl << endl;
	}
	
	return t;
}

/**
@brief 	Calculates a mask based on local variance.
@param 	kernel_size		size of kernel edge.
@param 	lowvar			low variance threshold.
@param 	bkg_flag		flag to generate a background with value -1.
@return Bimage* 		binary mask.

	The local variance within a kernel is calculated and used to generate the mask.
	A threshold is calculated to define the foreground and background.
	Excluded regions are defined by a low variance parameter.
	Foreground is set to 1, background is set to -1, excluded areas are set to 0.

**/
Bimage* 	Bimage::variance_mask(long kernel_size, double lowvar, int bkg_flag)
{
	Bimage*			pvar = copy();
	
	pvar->variance(kernel_size);

	long			i;
	double			t(pvar->variance_threshold(lowvar)), v;
	
	Bimage* 		pmask = copy_header();
	pmask->data_type(SCharacter);
	pmask->data_alloc_and_clear();
	
	if ( verbose ) {
	    cout << "Generating a variance mask:" << endl;
	    cout << "Thresholding local variance to: " << t << endl;
	    cout << "Low variance threshold:         " << lowvar << endl << endl;
	}
    
	if ( bkg_flag == 0 ) {			// Set only foreground
		for ( i=0; i<datasize; i++ )
			if ( (*pvar)[i] > t ) pmask->set(i, 1);			// Foreground
	} else for ( i=0; i<datasize; i++ ) {
		v = (*pvar)[i];
		if ( v > t ) pmask->set(i, 1);				// Foreground
		else if ( v >= lowvar ) pmask->set(i, -1);	// Background
	}
	
	image->FOM(pvar->image->FOM());
	pmask->image->FOM(pvar->image->FOM());

	delete pvar;
	
	pmask->statistics();
	
//	write_img("t.pif", pmask);
	
	return pmask;
}

/**
@brief 	Converts a mask to a tiled multilevel mask.
@param 	step		isotropic tile edge size.
@return Bimage*		new mask.

	The old data is unmodified.

**/
Bimage*		Bimage::tile_mask(long step)
{
	Bimage*		pmask = copy_header(1);
	pmask->data_type(Integer);
	pmask->compound_type(TSimple);
	pmask->channels(1);
	pmask->data_alloc();

	long		i, xx, yy, zz, tx, ty, tz, ti;
	long		sx = (x - 1)/step + 1;
	long		sy = (y - 1)/step + 1;
	
	for ( i=zz=0; zz<z; zz++ ) {
		tz = zz/step;
		for ( yy=0; yy<y; yy++ ) {
			ty = yy/step;
			for ( xx=0; xx<x; xx++, i++ ) {
				tx = xx/step;
				ti = (tz*sy + ty)*sx + tx;
				pmask->set(i, ti);
			}
		}
	}
	
	pmask->statistics();

	return pmask;
}

/**
@brief 	Converts a mask to a multilevel mask with a given number of voxels per level.
@param 	voxels_per_level	number of voxels per level.
@return long				number of voxels retained.

	The new data replaces the old data.
	Image statistics are recalculated.

**/
long			Bimage::mask_split(long voxels_per_level)
{
	if ( c > 1 ) {
		cerr << "Error: Multi-channel masks cannot be converted to multi-level masks!" << endl;
		return 0;
	}
	
    long			i, k, nv;
    double			value, threshold(0.5);
	int*			mask = new int[datasize];
	
	if ( threshold > max/2 ) threshold = max/2;
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Splitting a mask into multiple levels:" << endl;
	    cout << "Voxels per level:               " << voxels_per_level << endl;
	    cout << "Threshold:                      " << threshold << endl;
	}
	
	for ( i=nv=0, k=1; i<datasize; i++ ) {
		value = (*this)[i];
		if ( value > threshold ) {
			mask[i] = k;
			nv++;
		} else mask[i] = 0;
		if ( nv >= voxels_per_level ) {
			nv = 0;
			k++;
		}
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Levels generated:               " << k << endl << endl;
    
	data_type(Integer);

	data_assign((unsigned char *) mask);

	statistics();
	
	return nv;
}

/**
@brief 	Collapse a multi-level mask to a binary mask.
@return long			number of voxels retained.

	All non-zero segments are converted to ones.

**/
long		Bimage::levelmask_collapse()
{
	if ( verbose & VERB_LABEL )
		cout << "Collapsing a multi-level mask with " << max << " levels" << endl << endl;
	
    long			i, nv;
    int				value;
	unsigned char*	mask = new unsigned char[datasize];
	
	for ( i=nv=0; i<datasize; i++ ) {
		value = (int) (*this)[i];
		if ( value > 0 ) {
			mask[i] = 1;
			nv++;
		} else mask[i] = 0;
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Number of voxels retained:      " << nv << " (" << nv*100.0/datasize << " %)" << endl;
    
	data_type(UCharacter);

	data_assign(mask);

	statistics();
	
	return nv;
}

/**
@brief 	Adds a bilevel mask to a level mask.
@param 	pmask		bilevel mask to add.
@param	add_level	the level index to add.
@return long		number of selected regions.

	Where the bilevel mask is one, the level mask is set to given level value.

**/
long		Bimage::levelmask_add(Bimage* pmask, int add_level)
{
	long		i;

	if ( verbose & VERB_PROCESS )
		cout << "Add mask " << pmask->file_name() << " at level " << add_level << endl << endl;
	
	for ( i=0; i<data_size(); ++i )
		if ( (*pmask)[i] > 0.5 ) set(i, add_level);
	
	statistics();
	
	return 0;
}

/**
@brief 	Dilates a level mask.
@param	times		number of dilation operations.
@return long 		0.

	Traditional 3^dim kernel dilation.  Any pixel with a value of 1 turns 
	all of its neighbors to a value of 1.

**/
long		Bimage::levelmask_dilate(int times)
{
	int			t;
	
	for ( t=0; t<times; t++ )
		levelmask_dilate();
	
	return 0;
}

/**
@brief 	Dilates a level mask.
@return long 		0.

	Traditional 3^dim kernel dilation.  Any pixel with a value of 1 turns 
	all of its neighbors to a value of 1.

**/
long		Bimage::levelmask_dilate()
{
	long			i, j, im, nn, xx, yy, zz, kx, ky, kz;
	long			vmax(max+1), v;
	
	Vector3<long>	k(3, 3, 3);
	k = k.min(size());

	Vector3<long>	kh(k/2);
	Vector3<long> 	lo, hi;
	
	if ( verbose & VERB_PROCESS )
		cout << "Dilating a level mask" << endl << endl;
	
	int*			lev = new int[vmax];
	
	unsigned char*	mask = new unsigned char[data_size()];
	for ( i=0; i<data_size(); i++ ) mask[i] = (*this)[i];
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			lo[2] = (zz > kh[2])? zz - kh[2]: 0;
			hi[2] = (zz+kh[2] < z)? zz + kh[2]: z - 1;
			for ( yy=0; yy<y; yy++ ) {
				lo[1] = (yy > kh[1])? yy - kh[1]: 0;
				hi[1] = (yy+kh[1] < y)? yy + kh[1]: y - 1;
				for ( xx=0; xx<x; xx++, i++ ) {
					lo[0] = (xx > kh[0])? xx - kh[0]: 0;
					hi[0] = (xx+kh[0] < x)? xx + kh[0]: x - 1;
					if ( mask[i] < 1 ) {
						for ( j=0; j<vmax; j++ ) lev[j] = 0;
						for ( kz=lo[2]; kz<=hi[2]; kz++ ) {
							for ( ky=lo[1]; ky<=hi[1]; ky++ ) {
								for ( kx=lo[0]; kx<=hi[0]; kx++ ) {
									j = index(kx, ky, kz, nn);
									lev[mask[j]]++;
								}
							}
						}
						for ( j=1, v=im=0; j<vmax; j++ ) if ( v < lev[j] ) {
							im = j;
							v = lev[j];
						}
						if ( v ) set(i, im);
					}
				}
			}
		}
	}
	
	delete[] mask;
	delete[] lev;
	
	statistics();
	
	return 0;
}

/**
@brief 	Retains the selected levels in a multi-level mask.
@param 	&select_list	comma_separated list of levels to select.
@param	flag			0=binary mask, 1=multi-level mask.
@return long				number of voxels retained.

	The new data replaces the old data.
	The result is either a binary mask (flag=0) or a multi-level mask where
	the level indices changed to reflect the new range of selected levels (flag=1).
	Image statistics are recalculated.

**/
long		Bimage::levelmask_select(Bstring& select_list, int flag)
{
	if ( verbose & VERB_LABEL )
		cout << "Selecting levels " << select_list << endl;
	
    long			i, nv;
	long			nlev = (long) max + 1;
    int				value;

	vector<int>		lev = select_numbers(select_list, nlev);
	
	if ( flag )
		for ( i=1, nv=0; i<nlev; ++i )
			if ( lev[i] ) lev[i] = ++nv;

	for ( i=nv=0; i<datasize; ++i ) {
		value = long((*this)[i] + 0.5);
		if ( value >= 0 && value < nlev && lev[value] ) {
//			cout << tab << value;
			set(i, lev[value]);
			nv++;
		} else set(i, 0);
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Number of voxels retained:      " << nv << " (" << nv*100.0/datasize << " %)" << endl << endl;
    
	statistics();
	
	return nv;
}

/**
@brief 	Combines the selected levels in a multi-level mask and renumber.
@param 	&select_list	comma_separated list of levels to select.
@return long				number of voxels retained.

	The new data replaces the old data.
	The result is a multi-level mask where the selected levels are combined into one.
	Image statistics are recalculated.

**/
long		Bimage::levelmask_combine(Bstring& select_list)
{
	if ( verbose & VERB_LABEL )
		cout << "Combining levels " << select_list << endl;
	
    long			i, nv, first(0);
	long			nlev = (long) max + 1;
    int				value;

	vector<int>		lev = select_numbers(select_list, nlev);
	
	for ( i=1, nv=first=0; i<nlev; ++i ) {
		if ( first ) {
			if ( lev[i] ) lev[i] = first;
			else lev[i] = ++nv;
		} else {
			if ( lev[i] ) first = nv+1;
			lev[i] = ++nv;
		}
	}

	for ( i=nv=0; i<datasize; ++i ) {
		value = long((*this)[i] + 0.5);
		if ( value >= 0 && value < nlev && lev[value] ) {
//			cout << tab << value;
			set(i, lev[value]);
			nv++;
		} else set(i, 0);
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Number of voxels retained:      " << nv << " (" << nv*100.0/datasize << " %)" << endl << endl;
    
	statistics();
	
	return nv;
}

/**
@brief 	Retains the selected levels in a multi-level mask.
@param 	nn			image number.
@param 	voxel		voxel with level to select.
@return long			number of voxels retained.

	The new data replaces the old data.
	Image statistics are recalculated.

**/
long		Bimage::levelmask_select(long nn, Vector3<long> voxel)
{
	if ( verbose & VERB_LABEL )
		cout << "Selecting level from voxel:     " << voxel << endl;
	
    long			i, nv;
    long			level, value;
	
	i = index(voxel, nn);
	level = (long) (*this)[i];
	
	for ( i=nv=0; i<datasize; i++ ) {
		value = (long) (*this)[i];
		if ( fabs(value - level) < 0.5  ) {
			set(i, 1);
			nv++;
		} else set(i, 0);
	}
	
	if ( verbose & VERB_PROCESS )
	    cout << "Number of voxels retained:      " << nv << " (" << nv*100.0/datasize << " %)" << endl << endl;
    
	statistics();
	
	return nv;
}

/**
@brief 	Selects regions overlapping a mask.
@param 	pmask		overlap template mask.
@return	long			number of selected regions.

	The input mask can be of any form. All non-zero parts of this mask is used.

**/
long		Bimage::levelmask_select(Bimage* pmask)
{
	long			i;
	long			v;
	long			nsel(0), nlevel((long) (max + 1));
	vector<int>		levels(nlevel, 0);
	
	if ( verbose ) {
		cout << "Selecting levels overlapping with a mask:" << endl;
		cout << "Mask file:                      " << pmask->file_name() << endl;
		cout << "Input levels:                   " << nlevel-1 << endl;
	}
	
	for ( i=0; i<datasize; i++ ) {
		v = (long) (*this)[i];
		if ( (*pmask)[i] > 0.5 && v > 0 )
			levels[v]++;
	}

	for ( i=0; i<nlevel; i++ ) if ( levels[i] > 0 ) nsel++;
	
	if ( verbose ) {
		cout << "Selected levels:                " << nsel << endl << endl;
		cout << "Level\tCount" << endl;
		for ( i=0; i<nlevel; i++ ) if ( levels[i] )
			cout << i << tab << levels[i] << endl;
		cout << endl;
	}
	
	for ( i=0; i<datasize; i++ ) {
		v = long( (*this)[i] + 0.5 );
		if ( v < 1 || levels[v] < 1 )
			set(i, 0);
	}
	
	return nsel;
}

/**
@brief 	Switches two segments in a multi-level mask.
@param 	index1		first index.
@param 	index2		second index.
@return	long			number of levels.

	If the two indices are not found, the mask is not modified.

**/
long		Bimage::levelmask_switch(long index1, long index2)
{
	long			i, v;
	
	if ( verbose )
		cout << "Switching levels " << index1 << " and " << index2 << endl << endl;
	
	for ( i=0; i<datasize; ++i ) {
		v = long((*this)[i] + 0.5);
		if ( v == index1 ) set(i, index2);
		else if ( v == index2 ) set (i, index1);
	}

	statistics();
	
	return max;
}

/**
@brief 	Caclulates statistics for all the regions defined by a multi-level mask.
@param 	*pmask		multi-level mask map.
@return int 			0.

	The input data is replaced by the average for each region.
	A linked image is created to hold the variance for each region.
**/
int		 	Bimage::level_masked_stats(Bimage* pmask)
{
	long			i, j, nr(pmask->maximum()+1);
	double			v;
	vector<double>	avg(nr, 0);
	vector<double>	var(nr, 0);
	vector<double>	num(nr, 0);

	for ( i=0; i<image_size(); ++i ) {
		j = long((*pmask)[i] + 0.5);
		if ( j >= 0 && j < nr ) {
			v = (*this)[i];
			avg[j] += v;
			var[j] += v*v;
			num[j] += 1;
		}
	}
	
	cout << "Region\tNum\tAvg\tVar" << endl;
	for ( j=0; j<nr; ++j ) {
		if ( num[j] ) {
			avg[j] /= num[j];
			var[j] /= num[j];
			var[j] -= avg[j] * avg[j];
		}
		cout << j << tab << num[j] << tab << avg[j] << tab << var[j] << endl;
	}
	
	cout << endl;
	
//	Bimage*		pavg = next = copy();
//	Bimage*		pvar = next->next = copy();
	Bimage*		pavg = this;
	Bimage*		pvar = next = copy();

	for ( i=0; i<image_size(); ++i ) {
		j = (*pmask)[i];
		if ( j >= 0 && j < nr ) {
			pavg->set(i, avg[j]);
			pvar->set(i, var[j]);
		}
	}

	return 0;
}

/**
@brief 	Symetrizes a multi-level mask.
@param 	sym			point group symmetry.
@return long			maximum level index.

**/
long		Bimage::levelmask_symmetrize(Bsymmetry& sym)
{
	long			i, j, k, sr, si, xx, yy, zz;
	
	Vector3<double>	v, vt;
	double			d, dm;
	
	long			nsym(0);
	vector<Matrix3>	m = sym.matrices();
	nsym = m.size();

	if ( verbose & VERB_PROCESS ) {
		cout << "Symmetrizing a multi-level mask:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Origin:                         " << image->origin() << endl << endl;
	}
	
	vector<int>				symmap(max+1, 0);
	vector<Vector3<double>>	segcom(max+1);
	vector<int>				segnum(max+1, 0);

	for ( i=zz=0; zz<z; zz++ ) {
		v[2] = zz - image[0].origin()[2];
		for ( yy=0; yy<y; yy++ ) {
			v[1] = yy - image[0].origin()[1];
			for ( xx=0; xx<x; xx++, ++i ) {
				v[0] = xx - image[0].origin()[0];
				k = int(get(0, xx, yy, zz) + 0.5);
				if ( k > 0 ) {
					segcom[k] += v;
					segnum[k]++;
				}
			}
		}
	}
	
	for ( k=1; k<=max; ++k )
		segcom[k] /= segnum[k];

	if ( verbose & VERB_FULL ) {
		cout << "Seg\tx\ty\tz\td" << endl;
		for ( k=1; k<=max; ++k )
			cout << setprecision(3) << k << tab << segcom[k] << tab << segcom[k].length() << endl;
	}
	
	if ( verbose & VERB_FULL )
		cout << "i\titest\tisym\tinew\tx\ty\tz\td" << setprecision(3) << endl;
	for ( i=1, si=0; i<=max; ++i ) if ( symmap[i] < 1 ) {
		si++;
		if ( segcom[i].length() < 2 ) {		// Segments centered at the origin
			symmap[i] = si;
			if ( verbose & VERB_FULL )
				cout << i << "\t0\t" << i << tab << si << tab << segcom[i] << tab << segcom[i].length() << endl;
		} else {
			for ( j=0; j<nsym; j++ ) {
				v = m[j]*segcom[i];
				dm = 1e30;
				for ( k=1, sr=i; k<=max; ++k ) {
					d = segcom[k].distance(v);
					if ( dm > d ) {
						dm = d;
						sr = k;
					}
				}
				if ( dm < 2 ) symmap[sr] = si;
				if ( verbose & VERB_FULL )
					cout << i << tab << j << tab << sr << tab << si << tab << segcom[sr] << tab << dm << endl;
			}
		}
	}
	
	for ( i=0; i<data_size(); ++i ) {
		j = int((*this)[i] + 0.5);
		if ( j > 0 ) set(i, symmap[j]);
	}

	data_type(Integer);
	
	statistics();
		
	return max;
}

/**
@brief 	Extracts all the regions associated with a multi-level mask.
@param 	*pmask		multi-level mask map.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value of edge voxels.
@return Bimage* 		image with excised regions.

	The feature index is found in a feature map and the corresponding voxels
	in the density maps are extracted.

**/
Bimage* 	Bimage::level_mask_extract(Bimage* pmask, int fill_type, double fill)
{
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = image->background();

	long			i, j, k, cc;
	long			imgsize(x*y*z);
	
    Bimage*			pex = new Bimage(datatype, compoundtype, size(), (long) pmask->maximum());
	pex->sampling(sampling(0));
	for ( i=0; i<pex->images(); i++ ) pex->image[i] = image[0];
	if ( fill ) pex->fill(fill);
	
	if ( verbose & VERB_PROCESS )
		cout << "Excising " << pex->images() << " regions" << endl << endl;
	
	for ( i=0; i<imgsize; i++ ) if ( (*pmask)[i] > 0 ) {
		j = (long((*pmask)[i] - 1)*imgsize + i)*c;
		for ( cc=0, k=i*c; cc<c; cc++, j++, k++ ) pex->set(j, (*this)[k]);
	}
	
	return pex;
}

/**
@brief 	Calculates the average size of regions in a level mask.
@return double			average feature size.
**/
double	 	Bimage::levelmask_average_region_size()
{
	long			i;
	long			v, nr(0);
	long			nl((long) (max + 1.9));
	vector<int>		f(nl, 0);
	double			a(0);
	
	for ( i=0; i<nl; i++ ) f[i] = 0;
	
	for ( i=0; i<datasize; i++ ) {
		v = (long) (*this)[i];
		if ( v > 0 ) f[v]++;
	}

	for ( i=0; i<nl; i++ ) if ( f[i] ) {
		a += f[i];
		nr++;
//		cout << nr << tab << f[i] << endl;
	}
	
	if ( nr ) a /= nr;
	
	if ( verbose & VERB_FULL ) {
		cout << "Number of regions:              " << nr << endl;
		cout << "Average region size:            " << a << endl << endl;
	}
	

	return a;
}

/**
@brief 	Removes empty levels from a multi-level mask.
@return long		0.

	The mask is converted to an integer mask.

**/
long		Bimage::levelmask_clean()
{
	change_type(Integer);
	
	long			i, j, m;
	
	for ( i=m=0; i<datasize; i++ )
		if ( m < (*this)[i] )
			m = (long) (*this)[i];
	
	m++;
	
	vector<int>		map(m, 0);
	
	for ( i=0; i<datasize; i++ ) if ( (*this)[i] > 0 ) map[(int)(*this)[i]] += 1;
	
	for ( i=j=1; i<m; i++ ) if ( map[i] ) map[i] = j++;

	for ( i=0; i<datasize; i++ ) if ( (*this)[i] > 0 ) set(i, map[(int)(*this)[i]]);
	
	statistics();

	return 0;
}

/**
@brief 	Colorizes a multi-level mask with random color assignments.
@return	long		number of levels.

	A lookup table (LUT) is calculated for the range of gray-scale values
	and random colors assigned to each value.
	The image is then converted to RGB using the LUT.

**/
long		Bimage::levelmask_colorize()
{
	random_seed();
	
	long				i, cc, nc((long) (max - min + 1));
	vector<RGB<unsigned char>>	lut(nc);

	if ( verbose & VERB_LABEL ) {
	    cout << "Coloring a multi-level mask:" << endl;
	    cout << "Levels:                         " << nc << endl << endl;
	}
	
	for ( auto it = lut.begin(); it != lut.end(); ++it )
		it->random_color();
	
    RGB<unsigned char>* nudata = new RGB<unsigned char>[datasize];
	
	for ( i=0; i<datasize; i++ ) {
		cc = (long) (*this)[i] - min;
		nudata[i] = lut[cc];
	}
	
    data_type(UCharacter);
    compound_type(TRGB);
    c = 3;
	
    data_assign((unsigned char *) nudata);

	return 0;
}

/**
@brief 	Color a multi-level mask based on the volumes of regions.
@return Bimage*		new colored image.

**/
Bimage*		Bimage::levelmask_color_by_size()
{	
	Bimage*		pcol = new Bimage(UCharacter, TRGB, size(), n);
	pcol->sampling(sampling(0));
	pcol->origin(image->origin());
	RGB<unsigned char>*	color = (RGB<unsigned char> *) pcol->data_pointer();

	if ( verbose )
		cout << "Coloring mask levels by size" << endl;
		
	long				i, j, nseg((long)(max+ 1.9));
	long				vmin(datasize), vmax(0);
	double				vavg(0), vstd(0);
	RGB<unsigned char>	black(0,0,0);
	vector<int>			vol(nseg, 0);

	for ( i=0; i<datasize; i++ )
		if ( (*this)[i] > -1 ) vol[(long)(*this)[i]]++;
	
	for ( i=1; i<nseg; i++ ) {
		if ( vmin > vol[i] ) vmin = vol[i];
		if ( vmax < vol[i] ) vmax = vol[i];
		vavg += vol[i];
		vstd += vol[i] * vol[i];
	}
	vavg /= nseg - 1;
	vstd /= nseg - 1;
	vstd -= vavg * vavg;
	vstd = sqrt(vstd);
	if ( vmin < vavg - 3*vstd ) vmin = vavg - 3*vstd;
	if ( vmax > vavg + 3*vstd ) vmax = vavg + 3*vstd;
	
	if ( verbose )
		cout << "Coloring range:                 " << vmin << " - " << vmax << endl;
	
	for ( i=0; i<datasize; i++ ) {
		j = (*this)[i];
		if ( j > 0 )
			color[i].spectrum(vol[j], vmin, vmax);
		else
			color[i] = black;
	}
	
	pcol->statistics();
	
	return pcol;
}

/**
@brief 	Convert a mask to reflect region sizes.
@return int			0.

**/
int			Bimage::levelmask_region_size()
{
	change_type(Integer);
	
	if ( verbose )
		cout << "Generating a region size level mask" << endl;
		
	long				i, j, nseg((long)(max+1.9));
	vector<int>			vol(nseg, 0);

	for ( i=0; i<datasize; i++ )
		if ( (*this)[i] > -1 ) vol[(long)(*this)[i]]++;

	for ( i=0; i<datasize; i++ ) {
		j = (*this)[i];
		if ( j > 0 )
			set(i, vol[j]);
		else
			set(i, 0);
	}
	
	statistics();
		
	return 0;
}

/**
@brief 	Color a multi-level mask based on the volumes of regions.
@return Bplot*		region size histogram.

**/
Bplot*		Bimage::levelmask_size_histogram()
{
	long				i, j, nseg((long)(max+ 1.9));
	long				vmin(datasize), vmax(0);
	double				vavg(0), vstd(0);
	vector<int>			vol(nseg, 0);

	if ( verbose ) {
		cout << "Generating a mask level size histogram" << endl;
		cout << "Number of regions:              " << nseg << endl << endl;
	}

	for ( i=0; i<datasize; i++ )
		if ( (*this)[i] > -1 ) vol[(long)(*this)[i]]++;
	
	for ( i=1; i<nseg; i++ ) {
//		cout << i << tab << vol[i] << endl;
		if ( vmin > vol[i] ) vmin = vol[i];
		if ( vmax < vol[i] ) vmax = vol[i];
		vavg += vol[i];
		vstd += vol[i] * vol[i];
	}
	vavg /= nseg - 1;
	vstd /= nseg - 1;
	vstd -= vavg * vavg;
	vstd = sqrt(vstd);
//	cout << "Stats:" << tab << vmin << tab << vmax << tab << vavg << tab << vstd << endl;
	if ( vmin < vavg - 3*vstd ) vmin = vavg - 3*vstd;
	if ( vmax > vavg + 3*vstd ) vmax = vavg + 3*vstd;
//	cout << "Stats:" << tab << vmin << tab << vmax << tab << vavg << tab << vstd << endl;
	
	long			bins(256);
	double			scale = bins*1.0/(vmax - vmin);
	vector<int>		h(bins, 0);

	for ( i=1; i<nseg; i++ ) {
		j = (long) (vol[i] - vmin)*scale;
		if ( j < 0 ) j = 0;
		if ( j >= bins ) j = bins - 1;
		h[j]++;
	}
		
	long			ncol(2);
	Bstring			title("Histogram of level sizes");
	Bplot*			plot = new Bplot(1, bins, ncol);
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
	plot->page(0).column(1).color(0.7,0.7,0.7);
	plot->page(0).axis(1).min(vmin);
	plot->page(0).axis(1).max(vmax);
//	plot->page(0).axis(3).min(0);
//	plot->page(0).axis(3).max(hmax);
		
	for ( i=0, j=bins; i<bins; i++, j++ ) {
		(*plot)[i] = vmin + i/scale;
		(*plot)[j] = h[i];
	}
	
	return plot;
}

/**
@brief 	Calculates the interfaces between regions.
@param 	img_num			sub-image number.
@return Matrix			interface matrix.

	The mask must be an integer image.
	An interface matrix is calculated to count the number of connected
	voxels between every pair of regions.
	Each value i in a row j gives the number of voxels of the region j that
	are adjacent to one or more voxels in region i.
	Row 0 and column 0 are excluded from the calculation.

**/
Matrix		Bimage::mask_interface_matrix(int img_num)
{
	long			ns(img_num), ne(img_num);
	if ( img_num < 0 ) {
		ns = 0;
		ne = n - 1;
	}
	
	long 			i, j, xx, yy, zz, nn, kx, ky, kz, m(max+1.9);
	long			vi, vj;
	long			xlo, xhi, ylo, yhi, zlo, zhi;
	long			imgsize(x*y*z);
	long			ds(imgsize*n);

	for ( i=m=0; i<ds; i++ ) if ( m < (*this)[i] ) m = (long) (*this)[i];
	m++;
	
	if ( verbose & VERB_PROCESS )
		cout << "Calculating an interface matrix of " << m << " x " << m << endl;
	
	vector<int>		tmp(m,0);
	Matrix			ifmat = Matrix(m,m);
	
	for ( nn=ns, i=nn*imgsize; nn<=ne; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			zlo = (zz>0)? zz-1: zz;
			zhi = (zz<z-1)? zz+1: zz;
			for ( yy=0; yy<y; yy++ ) {
				ylo = (yy>0)? yy-1: yy;
				yhi = (yy<y-1)? yy+1: yy;
				for ( xx=0; xx<x; xx++, i++ ) {
					xlo = (xx>0)? xx-1: xx;
					xhi = (xx<x-1)? xx+1: xx;
					vi = (int) (*this)[i];
					if ( vi >= 0 ) {
						for ( j=0; j<m; j++ ) tmp[j] = 0;
						for ( kz=zlo; kz<=zhi; kz++ ) {
							for ( ky=ylo; ky<=yhi; ky++ ) {
								for ( kx=xlo; kx<=xhi; kx++ ) {
									j = index(0,kx,ky,kz,nn);
									vj = (int) (*this)[j];
									if ( vj >= 0 && vi != vj )
										tmp[vj] = 1;
								}
							}
						}
						for ( j=0; j<m; j++ )
							if ( tmp[j] ) ifmat[vi][j] += 1;
					}
				}
			}
		}
	}

	return ifmat;
}

/**
@brief 	Reports the whole interface matrix or interfaces for one region only.
@param 	reg_num			region to report for (<0 whole interface matrix).
@return long			0.

	The mask is converted to an integer mask.
	An interface matrix is calculated to count the number of connected
	voxels between every pair of regions.

**/
long		Bimage::mask_region_interfaces(int reg_num)
{
	long			i, j, m, nadj(0);
	
	Matrix			ifmat = mask_interface_matrix(0);
	m = ifmat.rows();
	
	if ( reg_num >= 0 ) {
		cout << "Interfaces for region " << reg_num << ":" << endl;
		cout << "Region\tVin\tVadj" << endl;
		for ( i=0; i<m; i++ ) {
			if ( i && ifmat[reg_num][i] ) nadj++;
			if ( ifmat[reg_num][i] && ifmat[i][reg_num] )
				cout << i << tab << ifmat[reg_num][i] << tab << ifmat[i][reg_num] << endl;
		}
		cout << endl;
		cout << "Interfaced regions:             " << nadj << endl;
	} else {
		cout << "Interface matrix:" << endl;
//		if ( verbose & VERB_FULL ) {
			cout << "Region";
			for ( i=0; i<m; i++ ) cout << tab << i;
			cout << endl;
			for ( j=0; j<m; j++ ) {
				cout << j;
				for ( i=0; i<m; i++ ) cout << tab << ifmat[j][i];
				cout << endl;
			}
//		}
	}
	cout << endl;
	
	return 0;
}

/**
@brief 	Calculates the interfaces between regions and deletes/merges small ones.
@param 	min_size	minimum size to accept regions, small ones are deleted/merged.
@param 	min_if		minimum interface size to consider a small region connected.
@return long		0.

	The mask is converted to an integer mask.
	An interface matrix is calculated to count the number of connected
	voxels between every pair of regions.
	Regions that are smaller than the minimum size are either deleted or
	merged with other regions.
	If the maximum interface size for a small region is less than the minimum
	specified, the region is deleted, otherwise, it is added to that
	neighboring region with which it shares the biggest interface.

**/
long		Bimage::mask_merge_delete(long min_size, long min_if)
{
	if ( min_if < 1 ) min_if = 1;
	
//    change_type(Integer);

	long 			i, j, k, l, nm, nd, end, nn;
	long			imgsize(x*y*z);
	long			ds(imgsize*n);
	int				m, imax, max;

	for ( i=m=0; i<ds; i++ ) if ( m < (*this)[i] ) m = (long)(*this)[i];
	m++;
	
	vector<int>		cnt(m,0);
	vector<int>		map(m,0);
	Matrix			ifmat;
	
	if ( verbose ) {
		cout << "Eliminating small regions by deletion and merging:" << endl;
		cout << "Minimum size:                   " << min_size << " (" << sampling(0).volume()*min_size << " A3)" << endl;
		cout << "Minimum interface area:         " << min_if << " (" << sampling(0)[0]*sampling(0)[1]*min_if << " A2)" << endl;
		cout << "Regions:                        " << m << endl;
	}
	
	for ( nn=0; nn<n; nn++ ) {
		ifmat = mask_interface_matrix(nn);
		for ( j=0; j<m; j++ ) cnt[j] = 0;
		for ( i=nn*imgsize, j=i+imgsize; i<j; i++ ) cnt[(int)(*this)[i]] += 1;
		if ( verbose & VERB_FULL ) {
			cout << "Region";
			for ( k=1; k<m; k++ ) cout << tab << k;
			cout << endl;
			for ( j=1; j<m; j++ ) {
				cout << j;
				for ( k=1; k<m; k++ ) cout << tab << ifmat[j][k];
				cout << endl;
			}
		}
		for ( j=1; j<m; j++ ) {		// Find the largest interface for small regions
			map[j] = j;
			if ( cnt[j] < min_size ) {
				max = min_if;
				if ( max > cnt[j]/2 ) max = cnt[j]/2;	// Change the interface cutoff for really small regions
				if ( max < 1 ) max = 1;
				for ( k=1, imax=0; k<m; k++ ) {
					if ( max <= ifmat[j][k] ) {
						max = ifmat[j][k];
						imax = k;
					}
				}
				map[j] = imax;
			}
		}
		if ( verbose & VERB_FULL ) {
			cout << "Region\tPointer" << endl;
			for ( j=1; j<m; j++ ) cout << j << tab << map[j] << endl;
		}
		for ( j=1; j<m; j++ ) if ( map[j] != j ) {	// Resolve linked references
			for ( k=1, cnt[0] = j, end=0; k<m && end==0; k++ ) {
				cnt[k] = map[cnt[k-1]];
				for ( l=0; l<k && cnt[k] != cnt[l]; l++ ) ;	// Detect loops and self-reference
				if ( l < k ) {		// Set new references
					end = cnt[k];
					for ( l=0; l<k; l++ ) map[cnt[l]] = end;
				}
			}
		}
		for ( j=1, nm=nd=0; j<m; j++ ) {		// Re-assign small regions
			if ( map[j] != j ) {
				if ( map[j] < 1 ) nd++;
				for ( k=0, l=n*imgsize; k<imgsize; k++, l++ )
					if ( (int)(*this)[l] == j ) set(l, map[j]);
				if ( verbose & VERB_FULL )
					cout << j << " -> " << map[j] << endl;
			} else nm++;
		}
	}
	
	if ( verbose ) {
		cout << "New regions:                    " << nm << endl;
		cout << "Regions deleted:                " << nd << endl << endl;
	}
	
	return 0;
}

/**
@brief 	Applies a soft mask to a sub-image.
@param 	nn			sub-image.
@param	pmask		soft mask: [0,1].
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value of edge voxels.
@return int			0.

	The image is multiplied with the mask, filling in the remainder 
	with the fill value.
	Requirement: The mask must range from 0 to 1.

**/
long		Bimage::apply_soft_mask(long nn, Bimage* pmask, int fill_type, double fill)
{
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

	long		i, j;
	
	for ( i=0, j=nn*image_size(); i<image_size(); ++i, ++j )
		set(j, (*pmask)[i] * (*this)[j] + (1 - (*pmask)[i]) * fill);
	
	return 0;
}


