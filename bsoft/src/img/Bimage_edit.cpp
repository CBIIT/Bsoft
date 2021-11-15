/**
@file	Bimage_edit.cpp
@brief	Library routines used for editing image contents
@author Bernard Heymann
@date	Created: 19980520
@date	Modified: 20210105
**/

#include "Bimage.h"
#include "matrix_linear.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


double		rectangle_edge(Vector3<double> d, Vector3<double> hr)
{
	d = d.abs() - hr;
	
	if ( hr[2] > 1 && d[1] < d[2] ) d[1] = d[2];
	if ( hr[1] > 1 && d[0] < d[1] ) d[0] = d[1];
	
	return d[0];
}

double		oval_edge(Vector3<double> d, Vector3<double> hr)
{
	double			f(0);
	
	Vector3<double>	de(d/hr);
	double			ld(d.length());
	
	if ( ld < 1 ) return -100;
	
	double			lde(de.length());
	
	if ( lde > 1e-20 )
		f = ld*(1.0 - 1.0/lde);
	
	return f;
}

double		cylinder_edge(Vector3<double> d, Vector3<double> hr)
{
	double			f(0), ze(fabs(d[2])-hr[2]);
	
	d[2] = 0;
	
	double			ld(d.length());
	
	if ( ld < 1 ) return ze;
	
	Vector3<double>	de(d/hr);
	double			lde(de.length());

	if ( lde > 1e-20 )
		f = ld*(1.0 - 1.0/lde);

	if ( f < ze ) f = ze;
	
	return f;
}

double		shell_edge(Vector3<double> d, double mind, double maxd)
{
	double			f(0);
	
	double			ld(d.length());
	
	if ( mind < 0.001 && ld < 0.001 ) return 1;
	
	mind = ld - mind ;
	maxd -= ld;
	
	if ( mind < maxd ) f = mind;
	else f = maxd;
	
	return f;
}


/**
@brief 	Creates a shape in an image and fills it with a constant value.
@param	type		type of shape: 0=rectangle, 1=oval, 2=cylinder
@param 	rect		three-value size of the area to be filled.
@param 	start 		three-value start of area.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	The edge of the area is smoothed with a function:
					    v_old(x,y,z) + fill*exp(1.618*dist/width)
		v_new(x,y,z) = ------------------------------------------
		                       1 + exp(1.618*dist/width)
	where	fill is the constant fill value.
			dist is the distance to the rectangular boundary defined by 
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.

**/
int			Bimage::shape(int type, Vector3<long> rect, Vector3<double> start,
				double width, int fill_type, double fill)
{
	long 			nn;
	
	for ( nn=0; nn<n; nn++ )
		shape(nn, type, rect, start, width, fill_type, fill);
	
	return 0;
}

/**
@brief 	Creates a shape in an image and fills it with a constant value.
@param	type		type of edge: 0=rectangle, 1=oval, 2=cylinder
@param 	nn			sub-image.
@param 	rect		three-value size of the area to be filled.
@param 	start 		three-value start of area.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	The edge of the area is smoothed with a function:
					    v_old(x,y,z) + fill*exp(1.618*dist/width)
		v_new(x,y,z) = ------------------------------------------
		                       1 + exp(1.618*dist/width)
	where	fill is the constant fill value.
			dist is the distance to the rectangular boundary defined by 
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.

**/
int			Bimage::shape(long nn, int type, Vector3<long> rect, Vector3<double> start,
				double width, int fill_type, double fill)
{
	if ( fabs(width) < 0.001 ) {
		if ( width < 0 ) width = -0.001;
		else width = 0.001;
	}
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

	if ( rect[0] < 1 ) rect[0] = x;
	if ( rect[1] < 1 ) rect[1] = y;
	if ( rect[2] < 1 ) rect[2] = z;
	if ( z < 2 ) {
		start[2] = 0;
		rect[2] = z;
	}

    long     			i, xx, yy, zz, cc;
    double   			f(0), edge(3*width), a(-GOLDEN/width);
	Vector3<double>   	d, hr(rect/2);
    Vector3<double>   	cx(start[0] + (rect[0] - 1)/2.0, start[1] + (rect[1] - 1)/2.0, start[2] + (rect[2] - 1)/2.0);
	Vector3<long> 		lo(start - edge);
	Vector3<long> 		hi(start + rect + edge);

	hr = hr.max(1);
	lo = lo.max(0);
	if ( hi[0] >= x ) hi[0] = x - 1;
	if ( hi[1] >= y ) hi[1] = y - 1;
	if ( hi[2] >= z ) hi[2] = z - 1;
	
//	cout << hr << endl;
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Filling a shape:" << endl;
		cout << "Shape:                          " << type << endl;
    	cout << "Start:                          " << setprecision(3) << start << endl;
    	cout << "Size:                           " << rect << endl;
    	cout << "Width and fill value:           " << width << " " << fill << endl << endl;
	}
	
	
	for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
		d[2] = zz - cx[2];
		for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
			d[1] = yy - cx[1];
			for ( xx=lo[0]; xx<=hi[0]; ++xx ) {
				d[0] = xx - cx[0];
				if ( type == 0 )
					f = rectangle_edge(d, hr);
				else if ( type == 1 )
					f = oval_edge(d, hr);
				else if ( type == 2 )
					f = cylinder_edge(d, hr);
				f *= a;
				if ( f > 50 ) edge = 1e30;
				else edge = exp(f);
				if ( !isfinite(edge) ) {
					cerr << xx << " " << yy << " " << zz << ": Value too large or not finite!: " << edge << endl;
					edge = 1e30;
				}
				i = index(xx, yy, zz, nn);
				for ( cc=0; cc<c; cc++, ++i )
					set(i, ((*this)[i] + fill*edge)/(1+edge));
			}
		}
	}
	
	statistics();
    
    return 0;
}

/**
@brief 	Creates a line in an image and fills it with a constant value.
@param 	start 		three-value start of line.
@param 	end			three-value end of line.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	The edge of the area is smoothed with a function:
					    v_old(x,y,z) + fill*exp(1.618*dist/width)
		v_new(x,y,z) = ------------------------------------------
		                       1 + exp(1.618*dist/width)
	where	fill is the constant fill value.
			dist is the distance to the rectangular boundary defined by
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.

**/
int			Bimage::line(Vector3<double> start, Vector3<double> end, double width, int fill_type, double fill)
{
	if ( fabs(width) < 0.001 ) {
		if ( width < 0 ) width = -0.001;
		else width = 0.001;
	}
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(long(0));

    long     			i, xx, yy, zz, cc;
    double   			f(0), t(0), edge(3*width), a(-GOLDEN/width);
	Vector3<double>   	v;

	if ( verbose & VERB_PROCESS ) {
	    cout << "Drawing a line:" << endl;
    	cout << "Start:                          " << setprecision(3) << start << endl;
    	cout << "End:                            " << end << endl;
    	cout << "Width and fill value:           " << width << " " << fill << endl << endl;
	}
	
	fill *= 2;
	
	for ( i=zz=0; zz<z; ++zz ) {
		v[2] = zz;
		for ( yy=0; yy<y; ++yy ) {
			v[1] = yy;
			for ( xx=0; xx<x; ++xx ) {
				v[0] = xx;
				t = v.position_relative_to_line(start, end);
				if ( t < 0 ) f = (v - start).length();
				else if ( t > 1 ) f = (v - end).length();
				else f = v.distance_from_line(start, end);
				f *= a;
				if ( f > 50 ) edge = 1e30;
				else edge = exp(f);
				if ( !isfinite(edge) ) {
					cerr << xx << " " << yy << " " << zz << ": Value too large or not finite!: " << edge << endl;
					edge = 1e30;
				}
				for ( cc=0; cc<c; cc++, ++i )
					set(i, ((*this)[i] + fill*edge)/(1+edge));
			}
		}
	}

	statistics();
    
    return 0;
}

/**
@brief 	Creates a mask with the edge approaching zero.
@param	type		type of edge: 0=rectangle, 1=oval, 2=cylinder
@param 	rect		three-value size of the area to be masked.
@param 	start 		three-value start for mask.
@param 	width		gaussian width of smoothing function.
@return Bimage*		new soft mask.

	The edge of the image is smoothed with a function:
					                1
		v_new(x,y,z) = 	-------------------------
						1 + exp(1.618*dist/width)
	where	dist is the distance to the boundary defined by 
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.
	With negative width values, the area filled is outside the shape.

**/
Bimage*	Bimage::edge_mask(int type, Vector3<long> rect, 
				Vector3<double> start, double width)
{
	if ( fabs(width) < 0.001 ) {
		if ( width < 0 ) width = -0.001;
		else width = 0.001;
	}
	
	if ( rect[0] < 1 ) rect[0] = x;
	if ( rect[1] < 1 ) rect[1] = y;
	if ( rect[2] < 1 ) rect[2] = z;
	if ( z < 2 ) {
		start[2] = 0;
		rect[2] = z;
	}

    long     			i, xx, yy, zz, cc;
    Vector3<double>   	d, hr(rect/2);
    Vector3<double>   	cx(start[0] + (rect[0] - 1)/2.0, start[1] + (rect[1] - 1)/2.0, start[2] + (rect[2] - 1)/2.0);
    double   			f(0), edge, a(GOLDEN/width);

	hr = hr.max(1);
	
//	cout << hr << endl;
	
	if ( verbose & VERB_PROCESS ) {
	    cout << "Creating an edge mask:" << endl;
		cout << "Shape:                          " << type << endl;
    	cout << "Start:                          " << setprecision(3) << start << endl;
    	cout << "Size:                           " << rect << endl;
    	cout << "Width:                          " << width << endl << endl;
	}
	
	Bimage*			pmask = new Bimage(Float, TSimple, size(), 1);
	
    for ( i=0, zz=0; zz<z; ++zz ) {
		d[2] = (double)zz - cx[2];
		for ( yy=0; yy<y; ++yy ) {
			d[1] = (double)yy - cx[1];
			for ( xx=0; xx<x; ++xx ) {
				d[0] = (double)xx - cx[0];
				if ( type == 0 )
					f = rectangle_edge(d, hr);
				else if ( type == 1 )
					f = oval_edge(d, hr);
				else if ( type == 2 )
					f = cylinder_edge(d, hr);
				f *= a;
				if ( f > 50 ) edge = 1e30;
				else edge = exp(f);
				if ( !isfinite(edge) ) {
					cerr << xx << " " << yy << " " << zz << ": Value too large or not finite!: " << edge << endl;
					edge = 1e30;
				}
				for ( cc=0; cc<c; cc++, ++i )
					pmask->set(i, 1/(1+edge));
			}
		}
	}
	
	pmask->statistics();
    
    return pmask;
}

/**
@brief 	Smooths the image edge with a soft rectangular function.
@param	type		type of edge: 0=rectangle, 1=oval, 2=cylinder
@param 	rect		three-value size of the area to be smoothed.
@param 	start 		three-value start for smoothing.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value of edge voxels.
@return int			0.

	The edge of the image is smoothed with a function:
					    v_old(x,y,z) + fill*exp(1.618*dist/width)
		v_new(x,y,z) = ------------------------------------------
		                       1 + exp(1.618*dist/width)
	where	fill is the desired edge value.
			dist is the distance to the rectangular boundary defined by 
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.
	With negative width values, the area filled is outside the shape.

**/
int			Bimage::edge(int type, Vector3<long> rect, Vector3<double> start,
				double width, int fill_type, double fill)
{
	Bimage*		pmask = edge_mask(type, rect, start, width);
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		apply_soft_mask(nn, pmask, fill_type, fill);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; nn++ )
		apply_soft_mask(nn, pmask, fill_type, fill);
#endif

	delete pmask;

	statistics();
	
	return 0;
}

/**
@brief 	Smooths the image edge with a soft rectangular function.
@param 	nn			sub-image.
@param	type		type of edge: 0=rectangle, 1=oval, 2=cylinder
@param 	rect		three-value size of the area to be smoothed.
@param 	start 		three-value start for smoothing.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value of edge voxels.
@return int			0.

	The edge of the image is smoothed with a function:
					    v_old(x,y,z) + fill*exp(1.618*dist/width)
		v_new(x,y,z) = ------------------------------------------
		                       1 + exp(1.618*dist/width)
	where	fill is the desired edge value.
			dist is the distance to the rectangular boundary defined by 
				the input size and start
			width is the gaussian width (softness)
	With very small values of the gaussian width, the edge approaches a
	step function.
	With negative width values, the area filled is outside the shape.

**/
int			Bimage::edge(long nn, int type, Vector3<long> rect, Vector3<double> start,
				double width, int fill_type, double fill)
{
	Bimage*		pmask = edge_mask(type, rect, start, width);

	apply_soft_mask(nn, pmask, fill_type, fill);

    return 0;
}

/**
@brief 	Apply Hanning taper window to the image.
@param 	fill		value of edge voxels.
@return int			0.

	Along each dimension, the image is multiplied with a function:
		v_new(i) = fill + ( v_old(i) - fill )* 0.5 * ( 1 - cos( 2*PI*i/(n-1) )) i=0,...,n-1
	where	fill is the desired edge value.
			n is the size of the image

**/
int			Bimage::hanning_taper(double fill)
{
	
    long     			i, nn, xx, yy, zz, cc;
	double				vol, v;
    Vector3<double>   	cx, han;

	cx[0] = M_PI * 2.0L / ( x - 1);
	cx[1] = M_PI * 2.0L / ( y - 1);
	cx[2] = M_PI * 2.0L / ( z - 1);
	
	if ( verbose & VERB_FULL ) {
	    cout << "Hann tapering window:" << endl;
    	cout << "Size:                    " << size() << endl;
    	cout << "Fill value:              " << fill << endl << endl;
	}
	
    for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			han[2] = 0.5 * (1. - cos(cx[2]*zz));
	    	for ( yy=0; yy<y; ++yy ) {
 				han[1] = 0.5 * (1. - cos(cx[1]*yy));
				for ( xx=0; xx<x; ++xx ) {
					han[0] = 0.5 * (1. - cos(cx[0]*xx));
					vol = han.volume();
    		    	for ( cc=0; cc<c; cc++, ++i ) {
						v = ((*this)[i] - fill)*vol + fill;
						set(i, v);
					}
				}
    		}
		}
    }
	
    return 0;
}

/**
@brief 	Fills a sphere within an image with a uniform value.
@param 	center		three vector center of sphere.
@param 	radius		sphere radius.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels within a sphere at a given location and with a given radius 
	are increased by a given fill value.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::sphere(Vector3<double> center, double radius,
				double width, int fill_type, double fill)
{
	Vector3<long>		rect((long) (2*radius+0.5), (long) (2*radius+0.5), (long) (2*radius+0.5));
	rect = rect.min(size());
	
	center -= rect/2;
	
	return shape(1, rect, center, width, fill_type, fill);
}

/**
@brief 	Fills a cylinder within an image with a uniform value.
@param 	center		three vector center of cylinder.
@param 	radius		cylinder radius.
@param 	height		cylinder heigth.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels within a cylinder at a given location are set to a given fill value.
	The height is along the z-direction.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::cylinder(Vector3<double> center, double radius,
				double height, double width, int fill_type, double fill)
{
	Vector3<long>		cyl(long(2*radius+0.5), long(2*radius+0.5), long(height+0.5));
	cyl = cyl.min(size());

	center[0] -= radius;
	center[1] -= radius;
	center[2] -= height/2;
	
	return shape(2, cyl, center, width, fill_type, fill);
}

/**
@brief 	Fills a gaussian sphere within an image with a uniform value.
@param 	nn			sub-image.
@param 	center		center of sphere.
@param 	sigma		Gaussian sigma value.
@param 	amp			amplitude.
@return int			0.

	All voxels within a sphere at a given location and with a given radius 
	are increased by a given fill value.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::gaussian_sphere(long nn, Vector3<double> center, double sigma, double amp)
{
	if ( sigma < 1 ) sigma = 1;
	if ( amp == 0 ) amp = 1;
	
	if ( verbose ) {
		cout << "Calculating a Gaussian sphere:" << endl;
		cout << "Center:                         " << center << endl;
		cout << "Sigma:                          " << sigma << endl;
		cout << "Amplitude:                      " << amp << endl << endl;
	}
	
    long     			i, xx, yy, zz, cc;
	double				x2, y2, z2, d2, v;
	double				fac(-0.5/(sigma*sigma));
	
    for ( i=zz=0; zz<z; ++zz ) {
		z2 = zz - center[2];
		z2 *= z2;
		for ( yy=0; yy<y; ++yy ) {
			y2 = yy - center[1];
			y2 *= y2;
			for ( xx=0; xx<x; ++xx ) {
				x2 = xx - center[0];
				x2 *= x2;
				d2 = x2 + y2 + z2;
				v = amp*exp(fac*d2);
				for ( cc=0; cc<c; cc++, ++i )
					add(i, v);
			}
		}
    }
	    
	return 0;
}

/**
@brief 	Fills a shell within an image with a uniform value.
@param 	center		center of shell.
@param 	minrad		minimum radius of shell.
@param 	maxrad		maximum radius of shell.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels within a shell at a given location and within given radii
	are set to a given fill value.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::shell(Vector3<double> center, double minrad,
				double maxrad, double width, int fill_type, double fill)
{
	long 			nn;
	
	for ( nn=0; nn<n; nn++ )
		shell(nn, center, minrad, maxrad, width, fill_type, fill);
	
	return 0;
}

/**
@brief 	Fills a shell within an image with a uniform value.
@param 	nn			sub-image.
@param 	center		center of shell.
@param 	minrad		minimum radius of shell.
@param 	maxrad		maximum radius of shell.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels within a shell at a given location and within given radii
	are set to a given fill value.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::shell(long nn, Vector3<double> center, double minrad,
				double maxrad, double width, int fill_type, double fill)
{
	if ( minrad > maxrad ) swap(minrad, maxrad);
	if ( maxrad - minrad < 1 ) return 0;
	
	if ( fabs(width) < 0.1 ) {
		if ( width < 0 ) width = -0.1;
		else width = 0.1;
	}
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

    long     			i, xx, yy, zz, cc;
	double				f, edge, pad(5*width), mrp(maxrad + pad), a(GOLDEN/fabs(width));
	Vector3<long> 		lo(center - mrp);
	Vector3<long> 		hi(center + mrp);
	Vector3<double>		d;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Filling shell:" << endl;
		cout << "Center:                         " << center << endl;
		cout << "Radii:                          " << minrad << " - " << maxrad << endl;
    	cout << "Gaussian width and fill value:  " << width << " " << fill << endl << endl;
	}
	
//	lo = center - (maxrad + pad);
	lo = lo.max(0);
//	hi = center + (maxrad + pad);
	if ( hi[0] >= x ) hi[0] = x - 1;
	if ( hi[1] >= y ) hi[1] = y - 1;
	if ( hi[2] >= z ) hi[2] = z - 1;
/*	
	cout << lo << tab << hi << endl;
	
	for ( xx=-2*maxrad; xx<2*maxrad; ++xx ) {
		d[0] = xx;
		f = shell_edge(d, minrad, maxrad);
		f *= a;
		if ( f > 40 ) edge = 1e30;
		else edge = exp(f);		
		cout << xx << tab << f << tab << edge << endl;
	}
*/	
    for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
		d[2] = zz - center[2];
		for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
			d[1] = yy - center[1];
			for ( xx=lo[0]; xx<=hi[0]; ++xx ) {
				d[0] = xx - center[0];
				f = shell_edge(d, minrad, maxrad);
				f *= a;
				if ( f > 40 ) edge = 1e30;
				else edge = exp(f);
				i = index(0, xx, yy, zz, nn);
				for ( cc=0; cc<c; cc++, ++i )
					set(i, ((*this)[i] + fill*edge)/(1+edge));
			}
		}
    }
	
	return 0;
}

/**
@brief 	Fills a shell within an image with a uniform value with wrapping.
@param 	nn			sub-image.
@param 	center		center of shell.
@param 	minrad		minimum radius of shell.
@param 	maxrad		maximum radius of shell.
@param 	width		gaussian width of smoothing function.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels within a shell at a given location and within given radii
	are set to a given fill value.
	The new data replaces the old data.
	The default center is {0,0,0}.

**/
int			Bimage::shell_wrap(long nn, Vector3<double> center, double minrad,
				double maxrad, double width, int fill_type, double fill)
{
	if ( minrad > maxrad ) swap(minrad, maxrad);
	if ( maxrad - minrad < 1 ) return 0;
	
	if ( fabs(width) < 0.1 ) {
		if ( width < 0 ) width = -0.1;
		else width = 0.1;
	}
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

    long     			i, xx, yy, zz, cc;
	long				ix, iy, iz;
	double				f, edge, mrw(maxrad + 3*width), a(GOLDEN/fabs(width));
	Vector3<long> 		sh(size()/2);
	Vector3<long> 		lo(center - mrw);
	Vector3<long> 		hi(center + mrw);
	Vector3<double>		d;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Filling shell with wrapping:" << endl;
		cout << "Center:                         " << center << endl;
		cout << "Radii:                          " << minrad << " - " << maxrad << endl;
    	cout << "Gaussian width and fill value:  " << width << " " << fill << endl << endl;
	}
	
//	lo = center - (maxrad + 3*width);
//	hi = center + (maxrad + 3*width);
	
    for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
		d[2] = zz - center[2];
		if ( d[2] > sh[2] ) d[2] -= z;
		if ( d[2] < -sh[2] ) d[2] += z;
		iz = (zz<0)? zz + z: zz;
		for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
			d[1] = yy - center[1];
			if ( d[1] > sh[1] ) d[1] -= y;
			if ( d[1] < -sh[1] ) d[1] += y;
			iy = (yy<0)? yy + y: yy;
			for ( xx=lo[0]; xx<=hi[0]; ++xx ) {
				d[0] = xx - center[0];
				if ( d[0] > sh[0] ) d[0] -= x;
				if ( d[0] < -sh[0] ) d[0] += x;
				ix = (xx<0)? xx + x: xx;
				f = shell_edge(d, minrad, maxrad);
				f *= a;
				if ( f > 50 ) edge = 1e30;
				else edge = exp(f);
				i = index(ix, iy, iz, nn);
				for ( cc=0; cc<c; cc++, ++i )
					set(i, ((*this)[i] + fill*edge)/(1+edge));
			}
		}
    }
	
	return 0;
}

/**
@brief 	Generates a bar between start and end points.
@param 	start		start of bar.
@param 	end			end of bar.
@param 	width		width of bar.
@param 	edge_width	bar edge width.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	A bar is generated with the given fill value from the start point
	to the end point and with a defined width and gaussian edge.

**/
int			Bimage::bar(Vector3<double> start, Vector3<double> end,
				double width, double edge_width, int fill_type, double fill)
{
	long 			nn;
	
	for ( nn=0; nn<n; nn++ )
		bar(nn, start, end, width, edge_width, fill_type, fill);
	
	return 0;
}

/**
@brief 	Generates a bar between start and end points.
@param 	nn			sub-image.
@param 	start		start of bar.
@param 	end			end of bar.
@param 	width		width of bar.
@param 	edge_width	bar edge width.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	A bar is generated with the given fill value from the start point
	to the end point and with a defined width and gaussian edge.

**/
int			Bimage::bar(long nn, Vector3<double> start, Vector3<double> end,
				double width, double edge_width, int fill_type, double fill)
{
	if ( edge_width < 0.001 ) edge_width = 0.001;
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

	if ( verbose & VERB_PROCESS ) {
		cout << "Generating a bar" << endl;
    	cout << "Start:                          " << start << endl;
    	cout << "End:                            " << end << endl;
		cout << "Fill value:                     " << fill << endl;
    	cout << "Width and edge width:           " << width << " " << edge_width << endl << endl;
	}
	
	long				i, xx, yy, zz, cc;
	double				t, d, half_width(width/2.0), edge;
    double   			a(GOLDEN/fabs(edge_width));
	
	Vector3<double>		vec(end - start);
	Vector3<double>		u, w;
	double				vlen(vec.length());
	double				vpow(vlen*vlen);
	Vector3<double> 	lo, hi;
	
	lo = start.min(end);
	lo -= 3*edge_width;
	lo = lo.max(0);
	hi = start.max(end);
	hi += 3*edge_width;
	if ( hi[0] >= x ) hi[0] = x - 1;
	if ( hi[1] >= y ) hi[1] = y - 1;
	if ( hi[2] >= z ) hi[2] = z - 1;
	
    for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
		u[2] = zz;
		for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
			u[1] = yy;
			for ( xx=lo[0]; xx<=hi[0]; ++xx, ++i ) {
				u[0] = xx;
				w = start - u;
				t = -w.scalar(vec)/vpow;
				if ( t >= 0 && t <= 1 ) {
					d = (vec.cross(w)).length()/vlen;
					edge = exp(a*(half_width - d));
					if ( !isfinite(edge) || edge > 1e30 ) edge = 1e30;
					i = index(xx, yy, zz, nn);
					for ( cc=0; cc<c; cc++, ++i )
						set(i, ((*this)[i] + fill*edge)/(1+edge));
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Generates a quadric surface over the whole image.
@param 	param		7-value array of parameters.
@return int			0.

**/
int			Bimage::quadric(double* param)
{
	long			i, xx, yy, zz, nn;
	double			dx, dy, dz, v;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			dz = -image[nn].origin()[2] + zz;
			for ( yy=0; yy<y; ++yy ) {
				dy = -image[nn].origin()[1] + yy;
				for ( xx=0; xx<x; ++xx, ++i ) {
					dx = -image[nn].origin()[0] + xx;
					v = param[0] + dx*param[1] + dy*param[2] + dz*param[3]
						+ dx*dx*param[4] + dy*dy*param[5] + dz*dz*param[6];
					set(i, v);
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates a chirp image.
@param 	freq_scale		frequency scale.
@param 	freq_shift		frequency shift (radians).
@return int				0.

**/
int			Bimage::chirp(double freq_scale, double freq_shift)
{
	double			rf(freq_scale/real_size()[0]);
	long			i, xx, yy, zz, nn;
	double			v;
	Vector3<double>	d;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			d[2] = -image[nn].origin()[2] + zz;
			for ( yy=0; yy<y; ++yy ) {
				d[1] = -image[nn].origin()[1] + yy;
				for ( xx=0; xx<x; ++xx, ++i ) {
					d[0] = -image[nn].origin()[0] + xx;
					v = cos(rf*d.length2()+freq_shift);
					set(i, v);
				}
			}
		}
	}
	
	return 0;
}


/**
@brief 	Fill the voxels that are not calculated.
@param 	step		step increment for voxels with proper values.
@return int 		0, <0 if error.

	The step size used to calculate sparse voxels is here used to determine
	a kernel size for the filling operation. For each sparse voxel, all the
	neighbors within the kernel are filled with the center value.

**/
int			Bimage::fill_gaps(long step)
{
	if ( compoundtype != TSimple ) {
		cerr << "Error: The gap filter only operates on single value images!" << endl << endl;
		return -1;
	}
	
	long				kernel_size(step);
	double				fill(200), v;
	if ( kernel_size%2 == 0 ) kernel_size += 1;
	if ( kernel_size < 3 ) kernel_size = 3;
	
	// Calculate the new size and set up
	long				i, j, nn, xx, yy, zz, kx, ky, kz;
	Vector3<long>		blocksize(kernel_size,kernel_size,kernel_size);
	blocksize = blocksize.min(size());
	long				blocktotalsize(blocksize.volume());
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Filling the gaps:" << endl;
		cout << "Block size:                     " << blocksize[0] << " x "
			<< blocksize[1] << " x " << blocksize[2] << " = " << blocktotalsize << endl;
	}

    change_type(Float);
	
	// Do the fill filtering
	long 				zlo, zhi, ylo, yhi, xlo, xhi;
	for ( nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz+=step ) {
			zlo = zz - blocksize[2]/2;
			if ( zlo < 0 ) zlo = 0;
			zhi = zz + blocksize[2]/2;
			if ( zhi >= z ) zhi = z - 1;
			for ( yy=0; yy<y; yy+=step ) {
				ylo = yy - blocksize[1]/2;
				if ( ylo < 0 ) ylo = 0;
				yhi = yy + blocksize[1]/2;
				if ( yhi >= y ) yhi = y - 1;
				for ( xx=0; xx<x; xx+=step ) {
					xlo = xx - blocksize[0]/2;
					if ( xlo < 0 ) xlo = 0;
					xhi = xx + blocksize[0]/2;
					if ( xhi >= x ) xhi = x - 1;
					i = index(0, xx, yy, zz, nn);
					v = (*this)[i];
					if ( v < fill ) {
						for ( kz=zlo; kz<=zhi; kz++ ) {
							for ( ky=ylo; ky<=yhi; ky++ ) {
								for ( kx=xlo; kx<=xhi; kx++ ) {
									j = index(0, kx, ky, kz, nn);
									set(j, v);
								}
							}
						}
					}
				}
			}
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fill_gaps: Done" << endl;
	
	statistics();
	
	return 0;
}

/**
@brief 	Interpolate the voxels that are not calculated.
@param 	step		step increment for voxels with proper values.
@return int 		0, <0 if error.

	Sparse voxels calculated on a regular grid with the given step size are 
	used to fill intermediate voxels with linearly interpolated values.

**/
int			Bimage::interpolate_gaps(long step)
{
	if ( compoundtype != TSimple ) {
		cerr << "Error: The gap filter only operates on single value images!" << endl << endl;
		return -1;
	}
	
	// Calculate the new size and set up
	long			i, nn, xx, yy, zz;
	long 			x1(0), x2(0), y1(0), y2(0), z1(0), z2(0);
	double			fx, fx1, fy, fy1, fz, fz1, v[8], vn;
	for ( i=0; i<8; ++i ) v[i] = 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Interpolating the gaps:" << endl;
		cout << "Step size:                      " << step << endl;
	}

    change_type(Float);
	
	// Do the fill filtering	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			if ( zz%step == 0 ) {
				z1 = zz;
				z2 = z1 + step;
			}
			fz = (zz - z1)*1.0/step;
			fz1 = 1 - fz;
			for ( yy=0; yy<y; ++yy ) {
				if ( yy%step == 0 ) {
					y1 = yy;
					y2 = y1 + step;
				}
				fy = (yy - y1)*1.0/step;
				fy1 = 1 - fy;
				for ( xx=0; xx<x; ++xx, ++i ) {
					if ( xx%step == 0 ) {
						x1 = xx;
						x2 = x1 + step;
						v[0] = get(nn, x1, y1, z1);
						v[1] = get(nn, x2, y1, z1);
						v[2] = get(nn, x1, y2, z1);
						v[3] = get(nn, x2, y2, z1);
						v[4] = get(nn, x1, y1, z2);
						v[5] = get(nn, x2, y1, z2);
						v[6] = get(nn, x1, y2, z2);
						v[7] = get(nn, x2, y2, z2);
					}
					fx = (xx - x1)*1.0/step;
					fx1 = 1 - fx;
					vn = fx1 * fy1 * fz1 * v[0];
					vn += fx * fy1 * fz1 * v[1];
					vn += fx1 * fy * fz1 * v[2];
					vn += fx * fy * fz1 * v[3];
					vn += fx1 * fy1 * fz * v[4];
					vn += fx * fy1 * fz * v[5];
					vn += fx1 * fy * fz * v[6];
					vn += fx * fy * fz * v[7];
					set(i, vn);
				}
			}
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::interpolate_gaps: Done" << endl;
	
	statistics();
	
	return 0;
}
