/**
@file	Bimage_montage.cpp
@brief	Functions to montage images for display
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20170117
**/

#include "Bimage.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			Bimage::montage_one(Bimage* pm, long zz, long mi, long cols, long rows, int flip)
{
	long	i, j, xx, yy, ix, iy, cc;
	long	ny = mi/cols;
	long	nx = mi - ny*cols;
	
	if ( flip & 1 ) nx = cols - nx - 1;
	if ( flip & 2 ) ny = rows - ny - 1;
	
	nx *= x;
	ny *= y;
	
	for ( j=zz*y*x*c, yy=0, iy=ny*pm->x; yy<y; yy++, iy+=pm->x ) {
		for ( xx=0, ix=nx; xx<x; xx++, ix++ ) {
			i = (iy + ix)*c;
			for ( cc=0; cc<c; cc++, i++, j++ ) pm->set(i, (*this)[j]);
		}
	}
	
	return 0;
}

/**
@brief 	Rearranges an image into a montage of 2D slices for display.
@param 	cols		columns in montage.
@param 	rows		rows in montage.
@param 	first		first slice in montage.
@param 	skip		number of slices to skip.
@param 	flip		flip the order of panels on: 1=x axis, 2=y axis.
@return Bimage*		montaged image.

	The slices of a 3D image are packed into a 2D montage.
	The background value for the image is used for empty regions.

**/
Bimage*		Bimage::montage(int first, int cols, int rows, int skip, int flip)
{
	if ( first < 0 ) first = 0;
	if ( first > n*z ) first = n*z - 1;
	
	long			incr(1);
	if ( skip > 0 ) incr = skip + 1;
	
	if ( cols < 1 ) {
		cols = (int) (sqrt((1.0*n*z - first)/incr) + 0.99);
		if ( cols*x > 1500 ) cols = 1500/x;
		if ( cols < 1 ) cols = 1;
		if ( cols > (n*z - first)/incr ) cols = (n*z - first)/incr;
	}
	if ( rows < 1 )
		rows = ((n*z - first)/incr - 1)/cols + 1;
	
	if ( verbose & VERB_LABEL )
	    cout << "Creating a montage of " << cols << " x " << rows << " slices" << endl << endl;
	
	Bimage*			pm = new Bimage(datatype, compoundtype, cols*x, rows*y, 1, 1);
	pm->sampling(sampling(0));
	pm->fill(background(long(0)));
	
#ifdef HAVE_GCD
	dispatch_apply(cols*rows, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		zz(first+i*incr);
		if ( zz < n*z ) montage_one(pm, zz, i, cols, rows, flip);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<cols*rows; i++ ) {
		long		zz(first+i*incr);
		if ( zz < n*z ) montage_one(pm, zz, i, cols, rows, flip);
	}
#endif
	
	pm->statistics();
	
	return pm;
}

