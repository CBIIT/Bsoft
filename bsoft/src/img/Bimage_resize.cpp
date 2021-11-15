/**
@file	Bimage_resize.cpp
@brief	Library routines to resize images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20150618
**/

#include "Bimage.h"
#include "Matrix.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Resizes without interpolation or rescaling.
@param 	nusize		new image size three-value vector.
@param 	translate	three-value translation vector.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill		value to fill in new regions.
@return int			0.

	An image is resized with translation and filling of new regions with 
	a given value.

**/
int			Bimage::resize(Vector3<long> nusize, Vector3<long> translate,
					int fill_type, double fill)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::resize: " << nusize << endl;

	if ( nusize == size() ) return 0;
	
	if ( nusize[0] < 1 ) nusize[0] = x;
	if ( nusize[1] < 1 ) nusize[1] = y;
	if ( nusize[2] < 1 ) nusize[2] = z;
	
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_MIN ) fill = min;
	if ( fill_type == FILL_MAX ) fill = max;

	long			i, j, xx, yy, zz, nn;
	long			oldx, oldy, oldz;
	long			elementsize = c*data_type_size();
	long			nualloc = (long) nusize.volume()*n*elementsize;
	long			fomsize = (long) nusize.volume()*n*sizeof(float);
	unsigned char*	nudata = new unsigned char[nualloc];
	float*			fom = NULL;
	float*			nufom = NULL;
	if ( next ) {
		fom = (float *) next->data_pointer();
		nufom = new float[fomsize];
		for ( i=0; i<fomsize; i++ ) nufom[i] = 0;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Resizing:" << endl;
		cout << "Shift:                          " << translate << endl;
		cout << "New size:                       " << nusize << endl;
		if ( fill_type != FILL_BACKGROUND )
			cout << "Fill value:                     " << fill << endl;
		cout << endl;
	}
	
	TypePointer		fp;

	for ( i=nn=0; nn<n; nn++ ) {
		image[nn].origin(image[nn].origin() + translate);
//		if ( verbose & VERB_PROCESS )
//			cout << "Image " << nn+1 << " new origin:             " << image[nn].origin() << endl;
		if ( fill_type == FILL_BACKGROUND ) {
			fill = background(nn);
//			if ( verbose & VERB_PROCESS )
//				cout << "Image " << nn+1 << " fill value:             " << fill << endl;
		}
		fp = fill_value(fill);
		for ( zz=0; zz<nusize[2]; zz++ ) {
			oldz = (long)zz - translate[2];
			for ( yy=0; yy<nusize[1]; yy++ ) {
				oldy = (long)yy - translate[1];
				for ( xx=0; xx<nusize[0]; xx++, i++ ) {
					oldx = (long)xx - translate[0];
					if ( oldx < 0 || oldx >= x || oldy < 0 ||
							oldy >= y || oldz < 0 || oldz >= z ) {
						memcpy(nudata+i*elementsize, fp.uc, elementsize);
					} else {
						j = index(oldx, oldy, oldz, nn);
						memcpy(nudata+i*elementsize, d.uc+j*elementsize, elementsize);
						if ( fom ) nufom[i] = fom[j];
					}
				}
			}
		}
		delete[] fp.uc;
	}
	
//	if ( verbose & VERB_PROCESS )
//		cout << endl;
	
	size(nusize);
	data_assign(nudata);

	if ( fom ) next->data_assign((unsigned char*) nufom);
	
	statistics();
	
	return 0;
}

/**
@brief 	Resizes without interpolation or rescaling.
@param 	nusize		new image size three-value vector.
@param 	translate	three-value translation vector.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill		value to fill in new regions.
@return Bimage*		resized image.

	An image is resized with translation and filling of new regions with 
	a given value.

**/
Bimage*		Bimage::resize_copy(Vector3<long> nusize, Vector3<long> translate,
					int fill_type, double fill)
{
	if ( !data_pointer() ) return NULL;
	
	if ( nusize[0] < 1 ) nusize[0] = x;
	if ( nusize[1] < 1 ) nusize[1] = y;
	if ( nusize[2] < 1 ) nusize[2] = z;
	
	long			i, j, cc, xx, yy, zz, nn;
	long			oldx, oldy, oldz;

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_MIN ) fill = min;
	if ( fill_type == FILL_MAX ) fill = max;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Resizing:" << endl;
		cout << "Shift:                          " << translate << endl;
		cout << "New size:                       " << nusize << endl;
		if ( fill_type != FILL_BACKGROUND )
			cout << "Fill value:                     " << fill << endl;
		cout << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Resizing" << endl << endl;

	Bimage*			pnu = copy_header();
	pnu->size(nusize);
	pnu->page_size(nusize);
	pnu->data_alloc();
	
	for ( i=nn=0; nn<n; nn++ ) {
		pnu->image[nn].origin(image[nn].origin() + translate);
		pnu->image[nn].origin(image[nn].origin() + translate);
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn+1 << " new origin:             " << pnu->image[nn].origin() << endl;
		if ( fill_type == FILL_BACKGROUND ) {
			fill = background(nn);
			if ( verbose & VERB_PROCESS )
				cout << "Image " << nn+1 << " fill value:             " << fill << endl;
		}
		for ( zz=0; zz<pnu->sizeZ(); zz++ ) {
			oldz = (long)zz - translate[2];
			for ( yy=0; yy<pnu->sizeY(); yy++ ) {
				oldy = (long)yy - translate[1];
				for ( xx=0; xx<pnu->sizeX(); xx++ ) {
					oldx = (long)xx - translate[0];
					if ( oldx < 0 || oldx >= x || oldy < 0 ||
							oldy >= y || oldz < 0 || oldz >= z ) {
						for ( cc=0; cc<c; cc++, i++ ) pnu->set(i, fill);
					} else {
						j = index(0, oldx, oldy, oldz, nn);
						for ( cc=0; cc<c; cc++, i++, j++ ) pnu->set(i, (*this)[j]);
					}
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	if ( pnu->statistics() )
		cerr << tab << "in Bimage::resize" << endl;

	return pnu;
}

/**
@brief 	Resizes without interpolation or rescaling.
@param 	nusize		new image size three-value vector.
@param 	translate	three-value translation vector.
@return int			0.

	An image is resized with translation and filling of new regions with 
	a given value.

**/
int			Bimage::resize_wrap(Vector3<long> nusize, Vector3<long> translate)
{
	if ( nusize[0] < 1 ) nusize[0] = x;
	if ( nusize[1] < 1 ) nusize[1] = y;
	if ( nusize[2] < 1 ) nusize[2] = z;
	
	long			i, j, xx, yy, zz, nn;
	long			oldx, oldy, oldz;
	long			elementsize = c*data_type_size();
	long			nualloc = (long) nusize.volume()*n*elementsize;
	long			fomsize = (long) nusize.volume()*n*sizeof(float);
	unsigned char*	nudata = new unsigned char[nualloc];
	float*			fom = NULL;
	float*			nufom = NULL;
	if ( next ) {
		fom = (float *) next->data_pointer();
		nufom = new float[fomsize];
		for ( i=0; i<fomsize; i++ ) nufom[i] = 0;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Resizing with wrapping:" << endl;
		cout << "Shift:                          " << translate << endl;
		cout << "New size:                       " << nusize << endl;
		cout << endl;
	}
	
	for ( i=nn=0; nn<n; nn++ ) {
//		image[nn].origin(image[nn].origin() + translate);
		image[nn].origin(image[nn].origin() + translate);
//		if ( verbose & VERB_PROCESS )
//			cout << "Image " << nn+1 << " new origin:             " << image[nn].origin() << endl;
		for ( zz=0; zz<nusize[2]; zz++ ) {
			oldz = (long)zz - translate[2];
			while ( oldz < 0 ) oldz += z;
			while ( oldz >= z ) oldz -= z;
			for ( yy=0; yy<nusize[1]; yy++ ) {
				oldy = (long)yy - translate[1];
				while ( oldy < 0 ) oldy += y;
				while ( oldy >= y ) oldy -= y;
				for ( xx=0; xx<nusize[0]; xx++, i++ ) {
					oldx = (long)xx - translate[0];
					while ( oldx < 0 ) oldx += x;
					while ( oldx >= x ) oldx -= x;
					j = index(oldx, oldy, oldz, nn);
					memcpy(nudata+i*elementsize, d.uc+j*elementsize, elementsize);
					if ( fom ) nufom[i] = fom[j];
				}
			}
		}
	}
	
//	if ( verbose & VERB_PROCESS )
//		cout << endl;
	
	size(nusize);
	data_assign(nudata);

	if ( fom ) {
		next->size(nusize);
		next->data_assign((unsigned char*) nufom);
	}
	
	statistics();
	
	return 0;
}

/**
@brief 	Resizes with wrapping without interpolation or rescaling.
@param 	nusize		new image size three-value vector.
@param 	translate	three-value translation vector.
@return Bimage*		resized image.

	An image is resized with translation and filling of new regions with 
	a given value.

**/
Bimage*		Bimage::resize_wrap_copy(Vector3<long> nusize, Vector3<long> translate)
{
	if ( !data_pointer() ) return NULL;
	
	if ( nusize[0] < 1 ) nusize[0] = x;
	if ( nusize[1] < 1 ) nusize[1] = y;
	if ( nusize[2] < 1 ) nusize[2] = z;
	
	long			i, j, cc, xx, yy, zz, nn;
	long			oldx, oldy, oldz;

	if ( verbose & VERB_PROCESS ) {
		cout << "Resizing with wrapping:" << endl;
		cout << "Shift:                          " << translate << endl;
		cout << "New size:                       " << nusize << endl;
		cout << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Resizing" << endl << endl;

	Bimage*			pnu = copy_header();
	pnu->size(nusize);
	pnu->page_size(nusize);
	pnu->data_alloc();
	
	for ( i=nn=0; nn<n; nn++ ) {
		pnu->image[nn].origin(image[nn].origin() + translate);
		if ( verbose & VERB_PROCESS )
			cout << "Image " << nn+1 << " new origin:             " << pnu->image[nn].origin() << endl;
		for ( zz=0; zz<pnu->sizeZ(); zz++ ) {
			oldz = (long)zz - translate[2];
			while ( oldz < 0 ) oldz += z;
			while ( oldz >= z ) oldz -= z;
			for ( yy=0; yy<pnu->sizeY(); yy++ ) {
				oldy = (long)yy - translate[1];
				while ( oldy < 0 ) oldy += y;
				while ( oldy >= y ) oldy -= y;
				for ( xx=0; xx<pnu->sizeX(); xx++ ) {
					oldx = (long)xx - translate[0];
					while ( oldx < 0 ) oldx += x;
					while ( oldx >= x ) oldx -= x;
					j = index(0, oldx, oldy, oldz, nn);
					for ( cc=0; cc<c; cc++, i++, j++ ) pnu->set(i, (*this)[j]);
				}
			}
		}
	}

	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	if ( pnu->statistics() )
		cerr << tab << "in Bimage::resize" << endl;

	return pnu;
}

/**
@brief 	Pads an image to a new size with a given fill value.
@param 	sz			new size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value to use when padding the image.
@return int			0.

	The image is enlarged with padding only on one side in each
	dimension with an input size greater than one.
	The data type is preserved.
**/
int			Bimage::pad(long sz, int fill_type, double fill)
{
	Vector3<long>	nusize(sz,sz,sz);
	if ( x < 2 ) nusize[0] = 1;
	if ( y < 2 ) nusize[1] = 1;
	if ( z < 2 ) nusize[2] = 1;
	
	Vector3<long>	translate;
	
	return resize(nusize, translate, fill_type, fill);
}

/**
@brief  Pads an image to a new size with a given fill value.
@param 	sz			new size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value to use when padding the image.
@return int			0.

 The image is enlarged with padding only on one side in each
 dimension with an input size greater than one.
 The data type is preserved.
**/
int			Bimage::pad(Vector3<long> sz, int fill_type, double fill)
{
	Vector3<long>	translate;
	
	return resize(sz, translate, fill_type, fill);
}

/**
@brief 	Pads an image to a new size with a given fill value.
@param 	sz			new size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value to use when padding the image.
@return Bimage*		resized image.

	The image is enlarged with padding only on one side in each
	dimension with an input size greater than one.
	The data type is preserved.
**/
Bimage*		Bimage::pad_copy(long sz, int fill_type, double fill)
{
	Vector3<long>	nusize(sz,sz,sz);
	if ( x < 2 ) nusize[0] = 1;
	if ( y < 2 ) nusize[1] = 1;
	if ( z < 2 ) nusize[2] = 1;
	
	Vector3<long>	translate;
	
	return resize_copy(nusize, translate, fill_type, fill);
}

/**
@brief  Pads an image to a new size with a given fill value.
@param 	sz			new size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value to use when padding the image.
@return Bimage*		resized image.

 The image is enlarged with padding only on one side in each
 dimension with an input size greater than one.
 The data type is preserved.
**/
Bimage*		Bimage::pad_copy(Vector3<long> sz, int fill_type, double fill)
{
	Vector3<long>	translate;
	
	return resize_copy(sz, translate, fill_type, fill);
}

/**
@brief 	Shrinks an image to a new size with wrapping of the excluded edges.
@param 	nusize		new image size three-value vector.
@param 	translate	three-value translation vector.
@return int 		0.

	An image is resized to a smaller size with translation and wrapping
	of the excluded edges to complete periodic boundaries. The image is
	first converted to floating point to prevent overflows of smaller 
	data types. The new image is finally converted back to the original
	data type.
	The new data replaces the old data.

**/
int 		Bimage::shrink_wrap(Vector3<long> nusize, Vector3<long> translate)
{
	if ( !data_pointer() ) return -1;
	
	if ( translate[0] == 0 && translate[1] == 0 && translate[2] == 0 &&
			nusize[0] == x && nusize[1] == y && nusize[2] == z )
				return 0;
	
	long			i, j, cc, xx, yy, zz, nn;
	long			nux, nuy, nuz;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Shrinking and wrapping:" << endl;
		cout << "New size:                       " << nusize << endl;
		cout << "New origin:                     " << translate << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Shrinking and wrapping" << endl << endl;

	DataType 		olddatatype = data_type();
	
	long			nudatasize = nusize[0]*nusize[1]*nusize[2]*n*c;
	float*			nudata = new float[nudatasize];
	for ( i=0; i<nudatasize; i++ ) nudata[i] = 0;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			nuz = zz + translate[2];
			if ( nuz < 0 ) nuz += nusize[2];
			if ( nuz >= nusize[2] ) nuz -= nusize[2];
			for ( yy=0; yy<y; yy++ ) {
				nuy = yy + translate[1];
				if ( nuy < 0 ) nuy += nusize[1];
				if ( nuy >= nusize[1] ) nuy -= nusize[1];
				for ( xx=0; xx<x; xx++ ) {
					nux = xx + translate[0];
					if ( nux < 0 ) nux += nusize[0];
					if ( nux >= nusize[0] ) nux -= nusize[0];
					j = (((nn*nusize[2] + nuz)*nusize[1] + nuy)*nusize[0] + nux)*c;
					for ( cc=0; cc<c; cc++, i++, j++ )
						nudata[j] += (*this)[i];
				}
			}
		}
	}
	
	size(nusize);
	page_size(nusize);
	data_type(Float);
	data_assign((unsigned char *) nudata);

	if ( statistics() )
		cerr << tab << "in Bimage::shrink_wrap" << endl;

	change_type(olddatatype);
	
	return 0;
}

/**
@brief 	Enlarges an image by an inetger scale.
@param 	scale		three-value scale.
@return int 		0.

	An image is enlarged by integer amounts.
	The new data replaces the old data.

**/
int 		Bimage::enlarge(Vector3<long> scale)
{
	if ( !data_pointer() ) return -1;
	
	Vector3<long>	nusize = size()*scale;
	
	long			i, j, xx, yy, zz, nn;
	long			oldx, oldy, oldz;
	long			elementsize = c*data_type_size();
	
	if ( verbose ) {
		cout << "Enlarging:" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "New size:                       " << nusize << endl << endl;
	}
	
	long			nudatasize = nusize[0]*nusize[1]*nusize[2]*n*elementsize;
	char*			nudata = new char[nudatasize];

	for ( i=nn=0; nn<n; nn++ ) {
		image[nn].origin(image[nn].origin()*scale);
		for ( zz=0; zz<nusize[2]; zz++ ) {
			oldz = zz/scale[2];
			for ( yy=0; yy<nusize[1]; yy++ ) {
				oldy = yy/scale[1];
				for ( xx=0; xx<nusize[0]; xx++, i++ ) {
					oldx = xx/scale[0];
					j = index(oldx, oldy, oldz, nn);
					memcpy(nudata+i*elementsize, d.uc+j*elementsize, elementsize);
				}
			}
		}
	}
	
	size(nusize);
	page_size(nusize);
	data_assign((unsigned char *) nudata);

	if ( statistics() )
		cerr << tab << "in Bimage::enlarge" << endl;

	return 0;
}

