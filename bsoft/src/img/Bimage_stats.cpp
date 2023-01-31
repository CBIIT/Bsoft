/**
@file	Bimage_stats.cpp
@brief	Functions to calculate statistics on image regions
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210806
**/

#include "Bimage.h"
#include "matrix_linear.h"
#include "Complex.h"
#include "Vector3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the statistics for an image.
@return long			number of errors.
**/
long			Bimage::statistics()
{
	if ( !d.uc ) return 0;
	
    long   				nn, notfin(0);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::statistics: Data size:  " << c << " x " << x << 
			" x " << y << " x " << z << " x " << n << " = " <<  datasize << endl;
		cout << "DEBUG Bimage::statistics: datatype=" << datatype << endl;
		cout << "DEBUG Bimage::statistics: min=" << min << " max=" << max << endl;
		cout << "DEBUG Bimage::statistics: ave=" << avg << " std=" << std << endl;
	}

#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t i){
   		statistics(i);
	});
#else
#pragma omp parallel for
    for ( long i=0; i<n; i++ )
   		statistics(i);
#endif
	
	if ( notfin ) {
		cerr << "Error in Bimage:statistics: " << notfin << " values not finite!" << endl;
		cerr << tab << "Image: " << file_name() << endl;
	}
	
	min = image->minimum();
	max = image->maximum();
	avg = std = 0;
	for ( nn=0; nn<n; nn++ ) {
		avg += image[nn].average();
		std += image[nn].standard_deviation()
			*image[nn].standard_deviation();
		if ( min > image[nn].minimum() ) min = image[nn].minimum();
		if ( max < image[nn].maximum() ) max = image[nn].maximum();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::statistics: avg[" << nn << "] = " << image[nn].average() << " std[" << nn << "] = " << image[nn].standard_deviation() << endl;
	}
	
	// Dubious fix!!!
	if ( datatype <= Short && compoundtype == TRGB ) {
		min = dtmin;
		max = dtmax;
	}

    avg = avg/n;
    std = std/n;
    if ( std > 0 ) std = sqrt(std);
	else std = 0;
	show_minimum(minimum());
	show_maximum(maximum());

	if ( verbose & VERB_STATS ) {
		cout << "Data size:                      " << c << " x " << x << 
			" x " << y << " x " << z << " x " << n << " = " <<  datasize << endl;
	    cout << "Min, max, avg, std:             " << min << " " << max << " " << avg << " " << std << endl; 
	    if ( next ) cout << "FOM maximum:                    " << next->maximum() << endl;
		cout << endl << "Image\tMin\tMax\tAvg\tStd" << endl;
		for ( nn=0; nn<n; nn++ )
			cout << nn+1 << tab << image[nn].minimum() << tab <<
				image[nn].maximum() << tab <<
				image[nn].average() << tab <<
				image[nn].standard_deviation() << endl;
		cout << endl;
	}

	if ( next ) next->statistics();

	return notfin;
}

/**
@brief 	Calculates the statistics for a sub-image.
@param	img_num			sub-image number.
@return long			number of errors.
**/
long			Bimage::statistics(long img_num)
{
	if ( !d.uc ) return 0;

	long				j, k, notfin(0);
	double				imin(DBL_MAX), imax(-DBL_MAX), iavg(0), istd(0), v, pwr;
	Complex<double>		cv;
	long				imagesize = x*y*z;
	if ( compoundtype > TComplex ) imagesize *= c;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::statistics: image = " << img_num << endl;

	for ( j=0, k=img_num*imagesize; j<imagesize; k++, j++ ) {
		if ( compoundtype != TComplex ) {
			v = (*this)[k];
			if ( isfinite(v) ) {
				if ( v < imin ) imin = v;
				if ( v > imax ) imax = v;
				iavg += v;
				istd += v*v;
			} else {
//				cerr << "Error in Bimage::statistics: value not finite: k=" << k << endl;
				notfin++;
				set(k, 0);
			}
		} else if ( compoundtype == TComplex ) {
			cv = complex(k);
			pwr = cv.power();
			if ( isfinite(pwr) ) {
				if ( pwr < imin ) imin = pwr;
				if ( pwr > imax ) imax = pwr;
				iavg += sqrt(pwr);
				istd += pwr;
			} else {
//				cerr << "Error in Bimage::statistics: complex value not finite: k=" << k << endl;
				notfin++;
				cv = 0;
				set(k, cv);
			}
		}
	}
	if ( imin < dtmin ) imin = dtmin;
	if ( imax > dtmax ) imax = dtmax;
	iavg /= imagesize;
	istd = istd/imagesize - iavg*iavg;
	if ( istd > 0 ) istd = sqrt(istd);
	else istd = 0;
	image[img_num].average(iavg);
	image[img_num].standard_deviation(istd);
	image[img_num].minimum(imin);
	image[img_num].maximum(imax);
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::statistics: min[" << img_num << "] = "
			<< imin << " max[" << img_num << "] = " << imax << endl;
	
	return notfin;
}

/**
@brief 	Calculates the statistics for a region in an image.
@param	*pmask			region.
@param	&regavg			region average to be calculated.
@param	&regstd			region standard deviation to be calculated.
@return long			number of voxels.
**/
long			Bimage::statistics(Bimage* pmask, double& regavg, double& regstd)
{
	if ( !d.uc ) return 0;
	
	regavg = regstd = 0;
	
	long			i, nv(0);
	
	for ( i=0; i<datasize; i++ ) if ( (*pmask)[i] > 0.5 ) {
		regavg += (*this)[i];
		regstd += (*this)[i] * (*this)[i];
		nv++;
	}
	
	regavg /= nv;
	regstd /= nv;
	if ( regstd - regavg*regavg > 0 )
		regstd = sqrt(regstd - regavg*regavg);
	else
		regstd = 0;
	
	return nv;
}

/**
@brief 	Checks whether the statistics conform to a Poisson distribution.
@return double		variance-to-average scale.

	A warning is issued when the variance/average differs more than 5% from one.
**/
double			Bimage::poisson_statistics_check()
{
	if ( std <= 0 ) statistics();
	
	double		var(std*std/avg);
	
	if ( fabs(var-1) > 0.05 )
		cerr << "Warning: Poisson statistics error! (Scale = " << var << ")" << endl;
	
	return var;
}

/**
@brief 	Calculates the statistics for an image within given radii from a location.
@param 	nn				sub-image number.
@param	loc				center of shell.
@param 	rad_min    		minimum radius (pixel units).
@param 	rad_max    		maximum radius (pixel units).
@return vector<double>	statistical values (none means failure).

	If a voxel lies within the specified radii, it is included in
	the statistical calculations.

	Return vector:	num, min, max, avg, std

**/
vector<double>	Bimage::stats_within_radii(long nn, Vector3<double> loc,
				double rad_min, double rad_max)
{
	vector<double>		stats;

	if ( !d.uc ) {
		cerr << "Error: No data for image " << file_name() << " in memory!" << endl;
		return stats;
	}
	
    long				i, xx, yy, zz, cc, npx(0);
	double				dx, dy, dz, d2, rmax2(rad_max*rad_max), rmin2(rad_min*rad_min);
	double				value, sum(0), ssum(0), w(0);
	double				vmin(1e37), vmax(-1e37), vavg(0), vstd(0);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::stats_within_radii: Data size: " << size() << "x" << c << "x" << n << "=" << datasize << endl;
		cout << "DEBUG Bimage::stats_within_radii: datatype=" << datatype << endl;
		cout << "DEBUG Bimage::stats_within_radii: min=" << min << " max=" << max << endl;
		cout << "DEBUG Bimage::stats_within_radii: ave=" << avg << " std=" << std << endl;
	}

	for ( zz=0; zz<z; zz++ ) {
		dz = (double)zz - loc[2];
		dz *= dz;
		if ( dz <= rmax2 ) for ( yy=0; yy<y; yy++ ) {
			dy = (double)yy - loc[1];
			dy *= dy;
			if ( dy <= rmax2 ) for ( xx=0; xx<x; xx++ ) {
				dx = (double)xx - loc[0];
				dx *= dx;
				d2 = dx + dy + dz;
				if ( d2 <= rmax2 && d2 >= rmin2 ) {
					value = 0;
					i = index(xx,yy,zz,nn);
					for ( cc=0; cc<c; cc++, i++ ) {
						value = (*this)[i];
						sum += value;
						ssum += value*value;
						npx += 1;
						if ( vmin > value ) vmin = value;
						if ( vmax < value ) vmax = value;
					}
				}
			}
		}
	}

	if ( npx ) {
		vavg = sum/npx;
		vstd = ssum/npx - vavg*vavg;
		if ( vstd > 0 ) vstd = sqrt(vstd);
		else vstd = 0;
		stats.push_back(npx);
		stats.push_back(vmin);
		stats.push_back(vmax);
		stats.push_back(vavg);
		stats.push_back(vstd);
	}
	
	if ( verbose & VERB_STATS ) {
	    cout << "Data size: " << size() << "x" << c << "x" << n << "=" << w << endl;
	    cout << "Min, max, avg, std:             " << 
				vmin << " " << vmax << " " << vavg << " " << vstd << endl << endl;
	}
	
	return stats;
}

/**
@brief 	Calculates the statistics for an image within the given box.
@param 	nn				sub-image number.
@param 	type			type of selection: 1=rectangle, 2=ellipse.
@param 	start			starting coordinates.
@param 	end				ending coordinates.
@return vector<double>	statistical values (none means failure).

	If a voxel lies within the specified box, it is included in
	the statistical calculations.
	
	Return vector:	num, min, max, avg, std

**/
vector<double>	Bimage::stats_in_shape(long nn, int type,
					Vector3<long> start, Vector3<long> end)
{
	int					inc;
	long				i, xx, yy, zz, npx(0);
	double				a(1), b(1), c(1), fx, fy, fz;
	double				value, vmin(1e30), vmax(-1e30), vavg(0), vstd(0), sum(0), sum2(0);
	vector<double>		stats;

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::stats_in_shape: start=" << start << endl;
		cout << "DEBUG Bimage::stats_in_shape: end=" << end << endl;
	}

	if ( ( start[0] >= x ) && ( end[0] >= x ) ) return stats;
	if ( ( start[1] >= y ) && ( end[1] >= y ) ) return stats;
//	if ( ( start[2] >= z ) && ( end[2] >= z ) ) return stats;
	start = start.max(0);
	end = end.max(0);
	if ( start[0] >= x ) start[0] = x - 1;
	if ( start[1] >= y ) start[1] = y - 1;
	if ( start[2] >= z ) start[2] = z - 1;
	if ( end[0] >= x ) end[0] = x - 1;
	if ( end[1] >= y ) end[1] = y - 1;
	if ( end[2] >= z ) end[2] = z - 1;	
	if ( start[0] > end[0] ) swap(start[0], end[0]);
	if ( start[1] > end[1] ) swap(start[1], end[1]);
	if ( start[2] > end[2] ) swap(start[2], end[2]);
	
	Vector3<long>		box = end - start;
	Vector3<float>		center;
	if ( type == 2 ) {
		if ( box[0] ) end[0] += 1;
		if ( box[1] ) end[1] += 1;
		if ( box[2] ) end[2] += 1;
	}
	box = end - start;
	center = start + box/2.0;
	if ( box[0] ) a = 2.0/box[0];
	if ( box[1] ) b = 2.0/box[1];
	if ( box[2] ) c = 2.0/box[2];

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::stats_in_shape: start=" << start << endl;
		cout << "DEBUG Bimage::stats_in_shape: end=" << end << endl;
		cout << "DEBUG Bimage::stats_in_shape: center=" << center << endl;
	}
	
	for ( npx=0, zz=start[2]; zz<=end[2]; zz++ ) {
		fz = (zz - center[2]) * c;
		fz *= fz;
		for ( yy=start[1]; yy<=end[1]; yy++ ) {
			fy = (yy - center[1]) * b;
			fy *= fy;
			for ( xx=start[0]; xx<=end[0]; xx++ ) {
				fx = (xx - center[0]) * a;
				fx *= fx;
				inc = 1;
				if ( type == 2 && fx + fy + fz > 1 ) inc = 0;
				if ( inc ) {
					i = index(0,xx,yy,zz,nn);
					value = (*this)[i];
					sum += value;
					sum2 += value*value;
					if ( vmin > value ) vmin = value;
					if ( vmax < value ) vmax = value;
					npx++;
				}
			}
		}
	}
	
	if ( npx ) {
		vavg = sum/npx;
		vstd = (sum2 - sum*sum/npx)/npx;
		if ( vstd > 0 ) vstd = sqrt(vstd);
		else vstd = 0;
		stats.push_back(npx);
		stats.push_back(vmin);
		stats.push_back(vmax);
		stats.push_back(vavg);
		stats.push_back(vstd);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::stats_in_shape: npx=" << npx << " avg=" << vavg << " vstd=" << std << endl;
	
	return stats;
}



/*
@brief 	Determines if a vector is inside or outside a polyhedral structure.
@param 	vec			vector to test for.
@param 	ndim		number of dimensions.
@param 	nvert		number of vertices.
@param 	*poly		set of vertices.
@return double		distance from closest edge (inside positive).

	Whether a point is inside or outside is based on the polygon edges.
	First the normal to an edge is calculated and tested whether it points
	inside or outside the polygon by finding the distance to the center:
		dist_cent = normal . edge_vertex
	which must be positive. Then the distance to the given point is:
		dist = normal . (edge_vertex - point)
	The minimum distance gives the distance of the point to the closest edge.
	If this distance is positive, the point is inside the polygon.
	Requirement: The polygon vertices must be given in order of connectivity.
**/
double		vector3_inside_outside(Vector3<double> vec, int ndim, int nvert, Vector3<double>* poly)
{
	int				i, j;
	double			d, dmin;
	Vector3<double>	normal;
	
	// Find the closest edge
	dmin = 1e30;
	for ( i=0, j=1; i<nvert; i++, j++ ) {
		if ( j >= nvert ) j = 0;
		normal = poly[i] - poly[j];
		normal = Vector3<double>(-normal[1], normal[0], 0);
		normal.normalize();
		if ( normal.scalar(poly[i]) < 0 ) normal = -normal;
		d = normal.scalar(poly[i] - vec);
		if ( dmin > d ) dmin = d;
	}

	return dmin;
}

/**
@brief 	Calculates the statistics for an image within the given polyhedron.
@param 	nn				sub-image number.
@param 	nvert			number of polygon vertices.
@param 	*poly			array of polygon vertices.
@return vector<double>	statistical values (none means failure).

	If a voxel lies within the specified polyhedron, it is included in
	the statistical calculations.

	Return vector:	num, min, max, avg, std

**/
vector<double>	Bimage::stats_in_poly(long nn, int nvert, Vector3<double>* poly)
{
	int					ndim, inc;
	long				i, xx, yy, zz, cc, npx(0);
	long				imgsize(x*y*z*c);
	double				value, d, dmin(1e30), dmax(0);
	double				vmin(0), vmax(0), sum(0), sum2(0), vavg(0), vstd(0);
	Vector3<double>		center, var, vec;
	vector<double>		stats;
	
	for ( i=0; i<nvert; i++ ) center += poly[i];
	
	center /= nvert;
	
	for ( i=0; i<nvert; i++ ) {
		poly[i] -= center;
		var += poly[i]*poly[i];
		d = poly[i].length();
		if ( dmin > d ) dmin = d;
		if ( dmax < d ) dmax = d;
	}
	
	for ( i=0, ndim=0; i<3; i++ ) if ( var[i] > 0 ) ndim++;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::stats_in_poly: center=" << center << endl;
		cout << "DEBUG Bimage::stats_in_poly: ndim=" << ndim << " dmin=" << dmin << " dmax=" << dmax << endl;
	}
	
	for ( i = nn*imgsize, zz=0; zz<z; zz++ ) {
		vec[2] = zz - center[2];
		for ( yy=0; yy<y; yy++ ) {
			vec[1] = yy - center[1];
			for ( xx=0; xx<x; x++ ) {
				vec[0] = xx - center[0];
				d = vec.length();
				inc = 0;
				if ( d <= dmax ) {
					if ( dmax == 0 ) inc = 1;
					else if ( vector3_inside_outside(vec, ndim, nvert, poly) > 0 )
						inc = 1;
				}
				for ( cc=0; cc<c; cc++, i++ ) if ( inc ) {
					value = (*this)[i];
					sum += value;
					sum2 += value*value;
					if ( vmin > value ) vmin = value;
					if ( vmax < value ) vmax = value;
					npx++;
				}
			}
		}
	}
	
	vavg = vstd = 0;
	if ( npx ) {
		vavg = sum/npx;
		vstd = (sum2 - sum*sum/npx)/npx;
		if ( vstd > 0 ) vstd = sqrt(vstd);
		else vstd = 0;
		stats.push_back(npx);
		stats.push_back(vmin);
		stats.push_back(vmax);
		stats.push_back(vavg);
		stats.push_back(vstd);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::stats_in_poly: npx=" << npx << " vavg=" << vavg << " vstd=" << vstd << endl;
	
	return stats;
}

/**
@brief 	Calculates the statistics for an image for each level in a mask.
@param 	nn			sub-image number.
@param 	*pmask		mask with two or more levels
@return long			number of levels (0 means failure).
**/
long		Bimage::stats_in_mask(long nn, Bimage* pmask)
{
	if ( !pmask ) return 0;
	
	long			nlevel = (long) (pmask->maximum() + 1.5);
	long			lev, i, j, imgsize(image_size());
	double			vavg, vstd;
	long*			cnt = new long[nlevel];
	double*			sum = new double[nlevel];
	double*			sum2 = new double[nlevel];

	for ( lev=0; lev<nlevel; lev++ ) sum[lev] = sum2[lev] = cnt[lev] = 0;

	for ( i=0, j=nn*imgsize; i<imgsize; i++, j++ ) {
		lev = (long) ((*pmask)[i] + 0.5);
		cnt[lev]++;
		sum[lev] += (*this)[j];
		sum2[lev] += (*this)[j] * (*this)[j];
	}
	
	cout << "Level\tCount\tSum\tAverage\tStDev" << endl;
	for ( lev=0; lev<nlevel; lev++ ) {
		vavg = vstd = 0;
		if ( cnt[lev] ) {
			vavg = sum[lev]/cnt[lev];
			vstd = sum2[lev]/cnt[lev] - vavg*vavg;
			if ( vstd > 0 ) vstd = sqrt(vstd);
			else vstd = 0;
		}
		cout << lev << tab << cnt[lev] << tab << sum[lev] << tab << vavg << tab << vstd << endl;
	}
	
    return nlevel;
}

int			Bimage::kernel_sums(long nn, long i, long ik, long nk)
{
	if ( size()[ik] < 2 ) return 0;

	if ( !next )
		next = new Bimage(Float, TSimple, size(), n);
	
	long		ii(i*x), j, j1, j2, inc(1);
	long		nl(size()[ik]), hk(nk/2), nv(0);
	double		s(0), s2(0);
	
	if ( ik == 1 ) {
		ii = (i%x) + (i/x) * x * y;
		inc = x;
	} else if ( ik == 2 ) {
		ii = i;
		inc = x * y;
	}
	
	ii += nn * x * y * z;
	
	double			step = inc*hk;
	vector<double>	v(nl);
	vector<double>	v2(nl);
	
	// Initial sums
	for ( j=ii, j1=0; j1<hk; j+=inc, j1++, nv++ ) {
		s += (*this)[j];
		s2 += (*next)[j];
	}
	
	// Main sums
	for ( j=ii, j1=0; j1<nl; j+=inc, j1++ ) {
		if ( j1 < nl-hk ) {	// Add the leading values
			j2 = j + (long) step;
			s += (*this)[j2];
			s2 += (*next)[j2];
			nv++;
		}
		v[j1] = s/nv;
		v2[j1] = s2/nv;
//		cout << j << " " << s << " " << s2 << endl;
		if ( j1 >= hk ) {	// Subtract the trailing values
			j2 = j - (long) step;
			s -= (*this)[j2];
			s2 -= (*next)[j2];
			nv--;
		}
	}

	for ( j=ii, j1=0; j1<nl; j+=inc, j1++ ) {
		set(j, v[j1]);
		next->set(j, v2[j1]);
	}
	
//	cout << i << ": " << ii << " " << inc << " " << n << " " << s << " " << s2 << endl;
	
	return 0;
}

int			Bimage::line_sums(long nn, long ik, long nk)
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
		kernel_sums(nn, i, ik, nk);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nls; i++ )
		kernel_sums(nn, i, ik, nk);
#endif

	return 0;
}

/**
@brief 	Calculates the local variance within the given kernel.
@param 	kernel_size		size of kernel edge.
@param	flag			selects type of output.
@return int 				0.

	The local variance within a kernel is calculated.
	The calculation is threaded if compiled with GCD or OpenMP.
	Flag values:
	0	variance
	1	standard deviation
	2	poisson excess variance

**/
int			Bimage::variance(long kernel_size, int flag)
{
	if ( kernel_size < 3 ) kernel_size = 3;
	if ( kernel_size > 999 ) kernel_size = 999;
	
	Vector3<long>		hk(kernel_size/2, kernel_size/2, kernel_size/2);
	if ( kernel_size > x ) hk[0] = x/2;
	if ( kernel_size > y ) hk[1] = y/2;
	if ( kernel_size > z ) hk[2] = z/2;
	Vector3<long>		k = hk*2 + 1;
	
	return variance(k, flag);
}

/**
@brief 	Calculates the local variance within the given kernel.
@param 	kernel_size		3-value size of kernel edge.
@param	flag			selects type of output.
@return int 				0.

	The local variance within a kernel is calculated.
	The calculation is threaded if compiled with GCD or OpenMP.
	Flag values:
	0	variance
	1	standard deviation
	2	poisson excess variance

**/
int			Bimage::variance(Vector3<long> kernel_size, int flag)
{
	change_type(Float);

	kernel_size = kernel_size.max(1);
	kernel_size = kernel_size.min(size());

	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating local variance:" << endl;
		cout << "Kernel size:                    " << kernel_size << endl << endl;
	}

	long			i, nn;
	next = new Bimage(Float, TSimple, size(), n);
	
	// Set up the square block
	for ( i=0; i<datasize; i++ ) next->set(i, (*this)[i] * (*this)[i]);

	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn << endl;
		for ( i=0; i<3; i++ )
			line_sums(nn, i, kernel_size[i]);
	}

	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	double		v;
	for ( i=0; i<datasize; i++ ) {
		v = (*this)[i];
		if ( flag < 2 ) {
			v = (*next)[i] - v*v;
			if ( v < 0 ) v = 0;
			if ( flag == 1 ) v = sqrt(v);
		} else {
			v = (*next)[i] - v*v - v;
		}
		set(i, v);
	}
	
	delete next;
	next = NULL;

	statistics();
	
	return 0;
}

int			Bimage::kernel_sums(long nn, long i, Bimage* pweight)
{
	i += nn*image_size();
	
	long			j, k, xx, yy, zz, kx, ky, kz;
	double			v, s(0), s2(0), w, ws(0);
	Vector3<long>	hk(pweight->size()/2);
	Vector3<long>	lo = kernel_low(i, hk);
	Vector3<long>	hi = kernel_high(i, hk);
	Vector3<long>	start = hk + lo - coordinates(i);
	start = start.max(0);
	
//	cout << coordinates(i) << tab << lo << tab << hi << tab << start << endl;

	for ( zz=lo[2], kz=start[2]; zz<=hi[2]; zz++, kz++ ) {
		for ( yy=lo[1], ky=start[1]; yy<=hi[1]; yy++, ky++ ) {
			kx = start[0];
			k = pweight->index(kx,ky,kz);
			xx = lo[0];
			j = index(0,xx,yy,zz,nn);
			for ( ; xx<=hi[0]; xx++, k++, j++ ) {
				w = (*pweight)[k];
				if ( w > 1e-30 ) {
					v = (*next)[j];
					ws += w;
					s += w*v;
					s2 += w*v*v;
				}
			}
		}
	}
	
	if ( ws ) {
		s /= ws;
		s2 /= ws;
		v = s2 - s*s;
		if ( v < 0 ) v = 0;
	} else v = 0;
	
	set(i, v);

	return 0;
}

/**
@brief 	Calculates the local variance weighed with the given image.
@param 	pweight		weight image.
@return int 		0.

	The local variance is calculated using the given image as kernel weights.
	The calculation is threaded if compiled with GCD or OpenMP.

**/
int			Bimage::variance(Bimage* pweight)
{
	if ( !pweight ) return -1;
	
	change_type(Float);

	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating local variance:" << endl;
		cout << "Weight image size:              " << pweight->size() << endl << endl;
	}

	long			nn, imgsize(x*y*z);
	next = copy();
	
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn << endl;
		// Do each voxel
#ifdef HAVE_GCD
		dispatch_apply(imgsize, dispatch_get_global_queue(0, 0), ^(size_t i){
			kernel_sums(nn, i, pweight);
		});
#else
#pragma omp parallel for
		for ( long i=0; i<imgsize; i++ )
			kernel_sums(nn, i, pweight);
#endif
	}

	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	delete next;
	next = NULL;

	statistics();
	
	return 0;
}


