/**
@file	Bimage_background.cpp
@brief	Library routines for manipulating image backgrounds
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20191002
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the background for one sub-image.
@param 	nn			sub-image.
@param	flag		flag to specify where to calculate the background: [0,3]
@return int 		0, <0 on error.

	The background is taken as the average of the values for:
	flag=0,1	outside the circle or sphere enclosed by the image
	flag=2		inside the circle or sphere enclosed by the image
	flag=3		inside the half-radius circle in the center of the image.

**/
int			Bimage::calculate_background(long nn, int flag)
{
	if ( !data_pointer() ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
    long				i, j, xx, yy, zz, cc, nm(0);
	double				x2(0), y2(0), z2(0), r;
	double				bkg(0);
	Vector3<long>		h(x/2, y/2, z/2);
	Vector3<double>		ori(image[nn].origin());
	
	i = (long) (nn*x*y*z);
	for ( zz=0; zz<z; zz++ ) {
		if ( z > 1 ) {
			z2 = (zz - ori[2])/h[2];
			z2 *= z2;
		}
		for ( yy=0; yy<y; yy++ ) {
			if ( y > 1 ) {
				y2 = (yy - ori[1])/h[1];
				y2 *= y2;
			}
			for ( xx=0; xx<x; xx++, i++ ) {
				if ( x > 1 ) {
					x2 = (xx - ori[0])/h[0];
					x2 *= x2;
				}
				r = x2 + y2 + z2 - 1;
				if ( flag == 2 ) r = -r;
				else if ( flag == 3 ) r = -r - 0.5;
				if ( r > 0 ) {
					if ( compoundtype == TComplex ) {
						bkg += complex(i).power();
						nm++;
					} else {
						for ( cc=0, j=i*c; cc<c; cc++, j++ ) {
							bkg += (*this)[j];
							nm++;
						}
					}
				}
			}
		}
	}

	if ( nm ) image[nn].background(bkg/nm);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::calculate_background: img=" << nn << " nm=" << nm <<
			" bkg=" << image[nn].background() << endl;
	
	return 0;
}

/**
@brief 	Calculates the background for each sub-image.
@param	flag		flag to specify where to calculate the background: [0,3]
@return int 		0, <0 on error.

	The background is taken as the average of the values outside the
	circle or sphere enclosed by the image.

**/
int			Bimage::calculate_background(int flag)
{
	if ( !data_pointer() ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
	if ( verbose & VERB_FULL )
		cout << "Calculating the background values of all images" << endl << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		calculate_background(nn, flag);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; nn++ )
		calculate_background(nn, flag);
#endif

	return 0;
}

/**
@brief 	Calculates the background for one sub-image.
@param 	*pmask		foreground mask.
@param 	nn			sub-image.
@param	flag		flag to specify where to calculate the background: [0,1]
@return int 		0, <0 on error.

	If the mask is specified, it is used to define the part of the 
	image to be used for background calculation (where the mask is zero).

**/
int			Bimage::calculate_background(Bimage* pmask, long nn, int flag)
{
	if ( !data_pointer() ) return -1;

	if ( !pmask ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
	long			i(nn*x*y*z), j, k, m(0), cc, nm(0);
	double			bkg(0);
	
	if ( pmask->n == n ) m = i;
	
	if ( flag ) {
		for ( k=0; k<image_size(); ++k, ++i, ++m ) if ( (*pmask)[m] > 0.5 ) {
			for ( cc=0, j=i*c; cc<c; cc++, j++ ) bkg += (*this)[j];
			nm++;
		}
	} else {
		for ( k=0; k<image_size(); ++k, ++i, ++m ) if ( (*pmask)[m] < 0.5 ) {
			for ( cc=0, j=i*c; cc<c; cc++, j++ ) bkg += (*this)[j];
			nm++;
		}
	}
	
	if ( nm ) image[nn].background(bkg/nm);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_calculate_background: img=" << nn << " nm=" << nm << " bkg=" << image[nn].background() << endl;

	return 0;
}


/**
@brief 	Calculates the background for each sub-image.
@param 	*pmask		foreground mask.
@param	flag		flag to specify where to calculate the background: [0,1]
@return int 		0, <0 on error.

	If the mask is specified, it is used to define the part of the 
	image to be used for background calculation (where the mask is zero).

**/
int			Bimage::calculate_background(Bimage* pmask, int flag)
{
	if ( !data_pointer() ) return -1;

	if ( !pmask ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
	if ( verbose & VERB_FULL )
		cout << "Calculating the background values of all images" << endl << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		calculate_background(pmask, nn, flag);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; nn++ )
		calculate_background(pmask, nn, flag);
#endif

	return 0;
}

/**
@brief 	Corrects the background for one sub-image.
@param 	nn			sub-image.
@param	flag		flag to specify where to calculate the background: [0,3]
@return int 		0, <0 on error.

	The background is taken as the average of the values outside the
	circle or sphere enclosed by the image.

**/
int			Bimage::correct_background(long nn, int flag)
{
	if ( !data_pointer() ) return -1;
	
	if ( data_type() == Bit ) return 0;

	calculate_background(nn, flag);
	
    long				i, j, xx, yy, zz, cc, nm(0);
	double				x2(0), y2(0), z2(0);
	double				bkg(image[nn].background());
	Vector3<long>		h(x/2, y/2, z/2);
	Vector3<double>		ori(image[nn].origin());
	
//	cout << "origin=" << ori << endl;

	i = (long) (nn*x*y*z);
	for ( zz=0; zz<z; zz++ ) {
		if ( z > 1 ) {
			z2 = (zz - ori[2])/h[2];
			z2 *= z2;
		}
		for ( yy=0; yy<y; yy++ ) {
			if ( y > 1 ) {
				y2 = (yy - ori[1])/h[1];
				y2 *= y2;
			}
			for ( xx=0; xx<x; xx++, i++ ) {
				if ( x > 1 ) {
					x2 = (xx - ori[0])/h[0];
					x2 *= x2;
				}
				if ( x2 + y2 + z2 > 1 ) {
					for ( cc=0, j=i*c; cc<c; cc++, j++ ) {
						set(j, bkg);
						nm++;
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::correct_background: img=" << nn << " bkg=" << bkg << endl;

	return 0;
}


/**
@brief 	Corrects the background for each sub-image.
@param	flag		flag to specify where to calculate the background: [0,1]
@return int 		0, <0 on error.

	The background is taken as the average of the values outside the
	circle or sphere enclosed by the image.

**/
int			Bimage::correct_background(int flag)
{
	if ( !data_pointer() ) return -1;

	if ( data_type() == Bit ) return 0;
		
	if ( verbose & VERB_FULL )
		cout << "Correcting the background values of all images" << endl << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		correct_background(nn, flag);
	});
#else
#pragma omp parallel for
	for ( long nn=0; nn<n; nn++ )
		correct_background(nn, flag);
#endif

	if ( statistics() )
		cerr << tab << "in Bimage::correct_background" << endl;
	
	return 0;
}

/**
@brief 	Corrects the background for each sub-image.
@param 	*pmask		foreground mask.
@param	flag		flag to specify where to calculate the background: [0,3]
@return int 		0.

	The background is taken as the average of the values outside the
	circle or sphere enclosed by the image.
	If the mask is specified, it is used to define the foreground of the 
	image. The correction then includes the area outside the circle or 
	sphere as well as the area where the mask is zero.

**/
int			Bimage::correct_background(Bimage* pmask, int flag)
{
	if ( !data_pointer() ) return -1;

	if ( !pmask ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
	if ( flag > 1 ) calculate_background(flag);
	else calculate_background(pmask, flag);

	long			i, j, k, m(0), nn, cc, nm;
	double			bkg;
		
	if ( verbose & VERB_FULL )
		cout << "Correcting the background values of all images" << endl << endl;
	
	for ( i=nn=0; nn<n; nn++ ) {
		if ( pmask->n == n ) m = i;
		else m = 0;
		bkg = image[nn].background();
		nm = 0;
		for ( k=0; k<image_size(); ++k, ++i, ++m ) if ( (*pmask)[m] < 0.5 ) {
			for ( cc=0, j=i*c; cc<c; cc++, j++ ) set(j, bkg);
			nm++;
		}
	}

	if ( statistics() )
		cerr << tab << "in Bimage::correct_background" << endl;
	
	return 0;
}

/**
@brief 	Subtracts the background for each sub-image.
@return int 		0, <0 on error.

	The background is taken from the sub-image structures.

**/
int			Bimage::subtract_background()
{
	return Bimage::shift_background(0);
}

/**
@brief 	Sets the background for each sub-image to a given value.
@param 	bkg			new background value.
@return int 		0, <0 on error.

	The background is taken as the average of the values outside the
	circle or sphere enclosed by the image.

**/
int			Bimage::shift_background(double bkg)
{
	if ( !data_pointer() ) return -1;
	
	if ( data_type() == Bit ) return 0;
	
    long    			i, j, nn, imagesize = x*y*z;
	if ( compound_type() != TComplex ) imagesize *= c;
	
	double				shift, v;
	Complex<double>		cv;

	if ( verbose & VERB_FULL )
		cout << "Shifting the background values of all images to " << bkg << endl;
	
    for ( nn=i=0; nn<n; nn++ ) {
//		shift = bkg - image[nn].background();
//		image[nn].background(bkg);
		shift = bkg - image[nn].background();
		image[nn].background(bkg);
	    for ( j=0; j<imagesize; i++, j++ ) {
			if ( compound_type() != TComplex ) {
				v = (*this)[i];
				v += shift;
				set(i, v);
			} else {
				cv = complex(i);
				cv += shift;
				set(i, cv);
			}
		}
	}
	
	if ( statistics() )
		cerr << tab << "in Bimage::shift_background" << endl;
	
	return 0;
}

