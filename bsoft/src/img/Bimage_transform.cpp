/**
@file	Bimage_transform.cpp
@brief	Library routines for transforming images
@author Bernard Heymann
@date	Created: 19990904
@date	Modified: 20190214
**/

#include "Bimage.h"
#include "math_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Transforms an image by translation, rotation, scaling and skewing, in place.
@param 	scale			3-value scale factor vector to apply.
@param 	origin			3-value origin for rotation and skewing.
@param 	translate		3-value vector for translation after transformation.
@param 	mat				3x3 rotation or skewing matrix.
@param 	fill_type		fill type for filling empty regions.
@param 	fill			value to fill in empty regions.
@return 0 				error code.

	A number of transformation types are combined in this function for
	efficiency. The basic operation is the conversion of the image data 
	from the original data block to a new data block, applying a 
	transformation matrix and translation vector with a specified origin 
	followed by an optional translation:
		y = R*x + t
	where	R is the rotation/skewing/scaling matrix.
			t is the translation vector.
	Note: The image is converted to floating point to improve interpolation.

**/
int			Bimage::transform(Vector3<double> scale,
				Vector3<double> origin, Vector3<double> translate,
				Matrix3 mat, int fill_type, double fill)
{
	long			nn;
	Bimage*			p;
	
	for ( nn=0; nn<n; nn++ ) {
		p = transform(nn, size(), scale, origin, translate, mat, fill_type, fill);
		replace(nn, p);
		delete p;
	}
	
	return 0;
}

/**
@brief 	Transforms an image by translation, rotation, scaling and skewing, returning a new image.
@param 	nusize			3-value new image size.
@param 	scale			3-value scale factor vector to apply.
@param 	origin			3-value origin for rotation and skewing.
@param 	translate		3-value vector for translation after transformation.
@param 	mat				3x3 rotation or skewing matrix.
@param 	fill_type		fill type for filling empty regions.
@param 	fill			value to fill in empty regions.
@return Bimage* 		new image.

	A number of transformation types are combined in this function for
	efficiency. The basic operation is the conversion of the image data 
	from the original data block to a new data block, applying a 
	transformation matrix and translation vector with a specified origin 
	followed by an optional translation:
		y = R*x + t
	where	R is the rotation/skewing/scaling matrix.
			t is the translation vector.
	Note: The image is converted to floating point to improve interpolation.

**/
Bimage*		Bimage::transform(Vector3<long> nusize, Vector3<double> scale,
				Vector3<double> origin, Vector3<double> translate,
				Matrix3 mat, int fill_type, double fill)
{
	long			nn;
	Bimage*			p;

	Bimage*			pt = NULL;

	if ( n == 1 ) {
		pt = transform(0, nusize, scale, origin, translate, mat, fill_type, fill);
	} else {
		pt = new Bimage(datatype, compoundtype, nusize, n);
		for ( nn=0; nn<n; nn++ ) {
			p = transform(nn, nusize, scale, origin, translate, mat, fill_type, fill);
			pt->replace(nn, p);
			pt->unit_cell(p->unit_cell());
			delete p;
		}
	}

	for ( nn=0; nn<n; nn++ )
		pt->image[nn].sampling(sampling(nn)/scale.abs());

	pt->statistics();
	
	return pt;
}

void		Bimage::transform_voxel(long i, Bimage* pt, long nn, Vector3<double> oldorigin, 
				Vector3<double> nuorigin, Matrix3 affmat, double fill)
{
	Vector3<double>	nu = pt->coordinates(i) - nuorigin;
	
	Vector3<double> old = affmat * nu + oldorigin;
	
	for ( long cc=0; cc<c; cc++, i++ )
		pt->add(i, interpolate(cc, old, nn, fill));
}

/**
@brief 	Transforms a sub-image by translation, rotation, scaling and skewing, returning a single new image.
@param	nn				sub-image to process.
@param 	nusize			3-value new image size.
@param 	scale			3-value scale factor vector to apply.
@param 	origin			3-value origin for rotation and skewing.
@param 	translate		3-value vector for translation after transformation.
@param 	mat				3x3 rotation or skewing matrix.
@param 	fill_type		fill type for filling empty regions.
@param 	fill			value to fill in empty regions.
@return Bimage* 		new image.

	A number of transformation types are combined in this function for
	efficiency. The basic operation is the conversion of the image data 
	from the original data block to a new data block, applying a 
	transformation matrix and translation vector with a specified origin 
	followed by an optional translation:
		y = R*x + t
	where	R is the rotation/skewing/scaling matrix.
			t is the translation vector.
	Note: The image is converted to floating point to improve interpolation.

**/
Bimage*		Bimage::transform(long nn, Vector3<long> nusize, Vector3<double> scale,
				Vector3<double> origin, Vector3<double> translate,
				Matrix3 mat, int fill_type, double fill)
{
	if ( !d.uc ) return NULL;
	
	if ( compoundtype > TSimple ) {
		error_show("Bimage::transform", __FILE__, __LINE__);
		cerr << "Error: Interpolation not done on compound data types" << endl;
		return NULL;
	}
	
	if ( fabs(scale[0]) < 0.001 || x < 2 ) scale[0] = 1;
	if ( fabs(scale[1]) < 0.001 || y < 2 ) scale[1] = 1;
	if ( fabs(scale[2]) < 0.001 || z < 2 ) scale[2] = 1;
	
	if ( nusize[0] < 1 ) nusize[0] = (long) (fabs(scale[0])*x);
	if ( nusize[1] < 1 ) nusize[1] = (long) (fabs(scale[1])*y);
	if ( nusize[2] < 1 ) nusize[2] = (long) (fabs(scale[2])*z);
	nusize = nusize.max(1);
	
	long			i, cc, xx, yy, zz;
	Vector3<double>	old, nu, oldorigin, nuorigin;
	View			view;
	Matrix3			viewmat;
	Matrix3			affmat = mat/scale;
	Vector3<double>	nusampling(sampling(0)/scale.abs()), nurealsize(nusampling*nusize);
	
//	change_type(Float);
	
	oldorigin = origin;
	nuorigin = origin + translate;
	
	if ( verbose & VERB_FULL ) {
		cout << endl << "Geometric transformation:" << endl;
		cout << "Transformation origin:          " << oldorigin << endl;
		cout << "New origin:                     " << nuorigin << endl;
		cout << "Scales:                         " << scale << endl;
		cout << "New size:                       " << nusize << " voxels" << endl;
		cout << "                                " << nurealsize << " A" << endl;
		cout << "Matrix:" << endl << mat;
		cout << "Fill value:                     " << fill << endl << endl;
	}
	
	affmat = affmat.transpose();
	
//	cout << affmat << endl;
//	cout << "New size:                       " << nusize << " voxels" << endl;

	Bimage*			pt = new Bimage(datatype, compoundtype, nusize, 1);
	pt->sampling(nusampling);

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = image[nn].background();

	// Note: the matrix is used in a back-calculation of old coordinates corresponding to new
	for ( i=zz=0; zz<pt->z; zz++ ) {
		nu[2] = (double)zz - nuorigin[2];
		for ( yy=0; yy<pt->y; yy++ ) {
			nu[1] = (double)yy - nuorigin[1];
			for ( xx=0; xx<pt->x; xx++ ) {
				nu[0] = (double)xx - nuorigin[0];
				old = affmat * nu + oldorigin;
//				if ( old[1] < 0 )
//					cout << old << endl;
				for ( cc=0; cc<c; cc++, i++ )
					pt->add(i, interpolate(cc, old, nn, fill));
			}
		}
	}
/*
#ifdef HAVE_GCD
	dispatch_apply(pt->data_size()/c, dispatch_get_global_queue(0, 0), ^(size_t i){
		transform_voxel(i*c, pt, nn, oldorigin, nuorigin, affmat, fill);
	});
#else
#pragma omp parallel for
	for ( i=0; i<pt->data_size(); i+=c )
		transform_voxel(i, pt, nn, oldorigin, nuorigin, affmat, fill);
#endif
*/	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::transform: view=" << image[nn].view() << endl;
	
	pt->image->origin(nuorigin);
	view = image[nn].view();
	viewmat = view.matrix();
//	cout << viewmat << endl;
	viewmat = viewmat * mat.transpose();
//	cout << viewmat << endl;
	view = View(viewmat);
	pt->image->view(view);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::transform: view=" << view << endl;
	
	UnitCell	uc(ucell);
	uc.size(nurealsize[0], nurealsize[1], nurealsize[2]);
	pt->unit_cell(uc);
	
	pt->statistics();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::transform: interpolation done" << endl;
	
	return pt;
}

/**
@brief 	Rotates an image using parameters defined in the image in place.
@return 0 				error code.

	The image is rotated according to the view and origin encoded
	for the sub-image.

**/
int			Bimage::rotate()
{
	Vector3<double> 	scale(1,1,1);
	Vector3<double> 	translate;

	Matrix3				mat = image->view().matrix();
	mat = mat.transpose();
	
	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image around the z-axis by the given angle in place.
@param 	angle		rotation angle in radians.
@return 0 				error code.

	The image is rotated according to the axis and angle given.

**/
int			Bimage::rotate(double angle)
{
	Vector3<double> axis(0,0,1);
	Matrix3			mat = Matrix3(axis, angle);
	
	Vector3<double> scale(1,1,1);
	Vector3<double> translate;

	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using an axis and angle in place.
@param 	axis		3-value rotation axis.
@param 	angle		rotation angle in radians.
@return 0 				error code.

	The image is rotated according to the axis and angle given.

**/
int			Bimage::rotate(Vector3<double> axis, double angle)
{
	Matrix3			mat = Matrix3(axis, angle);
	
	Vector3<double> scale(1,1,1);
	Vector3<double> translate;

	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image to a specified view in place.
@param 	*view		view vector and angle.
@return 0 				error code.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.

**/
int			Bimage::rotate(View view)
{
	Vector3<double>	translate;
	
	return rotate(translate, view);
}

/**
@brief 	Rotates an image to a specified view in place.
@param 	translate	translation after rotation.
@param 	*view		view vector and angle.
@return 0 			error code.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.

**/
int			Bimage::rotate(Vector3<double> translate, View view)
{
	Vector3<double>	scale(1,1,1);
	Matrix3			mat = view.matrix();

	change_type(Float);
	
	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using the specified matrix.
@param 	mat			3x3 rotation matrix.
@return 0 			error code.

	An image is rotated according to the input matrix.
	The input image origin is the rotation center.

**/
int			Bimage::rotate(Matrix3 mat)
{
	Vector3<double>	scale(1,1,1), translate;

	change_type(Float);
	
	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using the specified shift and matrix.
@param 	translate	translation after rotation.
@param 	mat			3x3 rotation matrix.
@return 0 			error code.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.

**/
int			Bimage::rotate(Vector3<double> translate, Matrix3 mat)
{
	Vector3<double>	scale(1,1,1);

	change_type(Float);
	
	return transform(scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using parameters defined in the image.
@param 	nusize		3-value new image size.
@return Bimage*		rotated image.

	The image is rotated according to the view and origin encoded
	for the sub-image.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize)
{
	Vector3<double> 	scale(1,1,1);
	Vector3<double> 	translate;

	Matrix3				mat = image->view().matrix();
	mat = mat.transpose();
	
	return transform(nusize, scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image around the z-axis by the given angle.
@param 	nusize		3-value new image size.
@param 	angle		rotation angle in radians.
@return Bimage*		rotated image.

	The image is rotated according to the axis and angle given.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize, double angle)
{
	Vector3<double> axis(0,0,1);
	Matrix3			mat = Matrix3(axis, angle);
	
	Vector3<double> scale(1,1,1);
	Vector3<double> translate;

	return transform(nusize, scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using an axis and angle.
@param 	nusize		3-value new image size.
@param 	axis		3-value rotation axis.
@param 	angle		rotation angle in radians.
@return Bimage*		rotated image.

	The image is rotated according to the axis and angle given.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize, Vector3<double> axis, double angle)
{
	Matrix3			mat = Matrix3(axis, angle);
	
	Vector3<double> scale(1,1,1);
	Vector3<double> translate;

	return transform(nusize, scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image to a specified view.
@param 	nusize		3-value new image size.
@param 	*view		view vector and angle.
@return Bimage*		rotated image.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize, View view)
{
	Vector3<double>	translate;
	
	return rotate(nusize, translate, view);
}

/**
@brief 	Rotates an image to a specified view.
@param 	nusize		3-value new image size.
@param 	translate	translation after rotation.
@param 	*view		view vector and angle.
@return Bimage*		rotated image.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize, Vector3<double> translate, View view)
{
	Vector3<double>	scale(1,1,1);
	Matrix3			mat = view.matrix();

	change_type(Float);
	
	return transform(nusize, scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image using the specified matrix.
@param	nusize		new image size.
@param 	mat			3x3 rotation matrix.
@return Bimage*		rotated image.

	An image is rotated according to the input matrix.
	The input image origin is the rotation center.

**/
Bimage*		Bimage::rotate(Vector3<long> nusize, Matrix3 mat)
{
	Vector3<double>	scale(1,1,1), translate;

	change_type(Float);
	
	return transform(nusize, scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);
}

/**
@brief 	Rotates an image to a specified view and adds it to another image.
@param 	*p			the image to rotate and add (modified).
@param 	origin		origin in input image.
@param 	*view		view vector and angle.
@return int			0.

	An image is rotated according to the input view vector and angle. 
	The input image origin is the rotation center.
	The rotated image is translated to the input origin.
	The result is added to the current image.

**/
int			Bimage::rotate_and_add(Bimage* p, Vector3<double> origin, View view)
{
	p->rotate(image->origin() - origin, view);
	
	add(p);

	return 0;
}

/**
@brief 	Calculates multiple copies oriented according to views.
@param	views		orientations for new images.
@return Bimage*		new set of images.

	Each image in a Bimage structure is rotated to the corresponding view.

**/
Bimage*		Bimage::orient(View* views)
{
	View*		v;
	long		nviews(0);
	
	for ( v = views; v; v = v->next ) nviews++;
	
	if ( nviews < 1 ) return NULL;
	
	change_type(Float);
	
	Bimage*		pori = copy(nviews);
	
	long		nn(0);
	Bimage*		p;	
	
	for ( v = views; v; v = v->next, ++nn ) {
		p = rotate(size(), *v);
		pori->replace(nn, p);
		delete p;
		pori->image[nn].view(*v);
		pori->image[nn].origin(image->origin());
	}
	
	return pori;
}

/**
@brief 	Inverts/mirrors each image through its origin.
@return int			number of images.

	Each image in a Bimage structure is inverted.

**/
int			Bimage::mirror()
{
	Vector3<double>	scale(1,1,1), shift;
	Matrix3			mat(-1,0,0,0,-1,0,0,0,-1);
	Bimage*			p1;
	
	change_type(Float);
	
	long			nn;
	
	for ( nn=0; nn<n; nn++ ) {
		p1 = transform(nn, size(), scale, image[nn].origin(), shift, mat, FILL_BACKGROUND, 0);
		replace(nn, p1);
		delete p1;
	}
	
	return nn;
}

/**
@brief Shifts an image.
@param 	vec			3-value real space shift vector.
@param 	fill_type	fill type for filling empty regions.
@param 	fill		value to fill in empty regions.
@return int			error code.
**/
int 		Bimage::shift(Vector3<double> vec, int fill_type, double fill)
{
	for ( long nn=0; nn<n; nn++ )
		shift(nn, vec, fill_type, fill);
	
	return 0;
}

/**
@brief Shifts one sub-image.
@param 	nn			sub-image.
@param 	vec			3-value real space shift vector.
@param 	fill_type	fill type for filling empty regions.
@param 	fill		value to fill in empty regions.
@return int			error code.
**/
int 		Bimage::shift(long nn, Vector3<double> vec, int fill_type, double fill)
{
	if ( vec.length() == 0 ) return 0;
	
	if ( verbose & VERB_FULL )
		cout << "Translate by:                   " << vec << endl;

	long			i, j, xx, yy, zz, cc;
	long			imgsize = c*image_size();
	vector<float>	temp(imgsize);

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) fill = background(nn);

	for ( j=zz=0; zz<z; zz++ ) {
		for ( yy=0; yy<y; yy++ ) {
			for ( xx=0; xx<x; xx++ ) {
				for ( cc=0; cc<c; cc++, j++ ) {
					temp[j] = interpolate(cc, xx-vec[0], yy-vec[1], zz-vec[2], nn, fill);
				}
			}
		}
	}
	
	for ( i=nn*imgsize, j=0; j<imgsize; i++, j++ ) set(i, temp[j]);
	
	image[nn].origin(image[nn].origin() + vec);
	
	return 0;
}

/**
@brief Shifts an image.
@param 	vec		3-value real space shift vector.
@return int		error code.
**/
int 		Bimage::shift_wrap(Vector3<double> vec)
{
	for ( long nn=0; nn<n; nn++ )
		shift_wrap(nn, vec);
	
	return 0;
}

/**
@brief Shifts one sub-image with wrapping.
@param 	nn		sub-image.
@param 	vec		3-value real space shift vector.
@return int		error code.
**/
int 		Bimage::shift_wrap(long nn, Vector3<double> vec)
{
	if ( vec.length() == 0 ) return 0;
	
	if ( verbose & VERB_FULL )
		cout << "Translate by:                   " << vec << endl;

	long			i, j, xx, yy, zz, cc;
	long			imgsize(c*image_size());
	Vector3<double>	loc;
	vector<float>	temp(imgsize);

//	cout << "start interpolating " << nn << endl;
	for ( j=zz=0; zz<z; zz++ ) {
		loc[2] = -vec[2] + zz;
		for ( yy=0; yy<y; yy++ ) {
			loc[1] = -vec[1] + yy;
			for ( xx=0; xx<x; xx++ ) {
				loc[0] = -vec[0] + xx;
				for ( cc=0; cc<c; cc++, j++ )
					temp[j] = interpolate_wrap(cc, loc, nn);
			}
		}
	}
//	cout << "done interpolating " << nn << endl;
	
	for ( i=nn*imgsize, j=0; j<imgsize; i++, j++ ) set(i, temp[j]);
	
	image[nn].origin(image[nn].origin() + vec);
	
	return 0;
}

/**
@brief 	Centers an image.
@param 	fill_type	fill type for filling empty regions.
@param 	fill		value to fill in empty regions.
@return int			error code.
**/
int 		Bimage::center(int fill_type, double fill)
{
	Vector3<double>	origin(size()/2);
	
	for ( long nn=0; nn<n; nn++ )
		shift(nn, origin - image[nn].origin(), fill_type, fill);
	
	return 0;
}

/**
@brief 	Centers an image with wrapping.
@return int			error code.
**/
int 		Bimage::center_wrap()
{
	Vector3<double>	origin(size()/2);
	
	for ( long nn=0; nn<n; nn++ )
		shift_wrap(nn, origin - image[nn].origin());
	
	return 0;
}

/**
@brief	Shifts the origins to zero with wrapping.
@return int			error code.
**/
int 		Bimage::zero_origin()
{
	for ( long nn=0; nn<n; nn++ )
		shift_wrap(nn, -image[nn].origin());
	
	return 0;
}

/**
@brief	Scales and image to the given reference.
@param	pref		reference/template image.
@return	Bimage*		scaled image.
**/
Bimage*		Bimage::scale_to_same_size(Bimage* pref)
{
	if ( !pref ) return NULL;
	
	if ( image->origin().length() < 1 ) origin(size()/2);
	
//	double			scale(sampling()[0]/pref->sampling()[0]);
//	Vector3<double> vscale(scale, scale, scale);
//	Vector3<double> origin(image->origin()), translate((pref->size() - size())/2);
//	Matrix3 		mat(1);

//	return transform(pref->size(), vscale, origin, translate, mat, FILL_BACKGROUND);
	return fspace_resize(pref);
}

/**
@brief	Scales and image to the given reference, searching for the correct scale.
@param	pref		reference/template image.
@param	pmask		mask to limit correlation calculation (can be NULL).
@param	scalemin	minimum scaling to start search.
@param	scalemax	maximum scaling to end search.
@param	step		scaling search step size.
@return	Bimage*		scaled image.
**/
Bimage*		Bimage::scale_to_reference(Bimage* pref, Bimage* pmask, double scalemin, 
				double scalemax, double step)
{
	if ( image->origin().length() < 1 ) origin(size()/2);
	
	double			scale, scalebest(1), cc, ccbest(0);
	Vector3<double> vscale, origin(image->origin()), translate((pref->size() - size())/2);
	Matrix3 		mat(1);
	Bimage*			pt;
	
	if ( verbose ) {
		cout << "Finding the best scaling to the reference map:" << endl;
		cout << "Scale range:                    " << scalemin << " - " << scalemax << endl;
		cout << "Initial scale step size:        " << step << endl << endl;
		cout << "Scale\tSampling\tCC" << endl;
	}
	
	while ( step > 0.0001 ) {
		for ( scale=scalemin; scale<=scalemax+0.001; scale+=step ) {
			vscale = Vector3<double>(scale, scale, scale);
			pt = transform(pref->size(), vscale, origin, translate, mat, FILL_BACKGROUND);
			cc = pt->correlate(pref, 0, pref->sizeX()/2, pmask);
			if ( ccbest < cc ) {
				ccbest = cc;
				scalebest = scale;
			}
			cout << scale << tab << pt->sampling(0)[0] << tab << cc << endl;
			delete pt;
		}
		scalemin = scalebest - 0.75*step;
		scalemax = scalebest + 0.75*step;
		step /= 2;
	}
	
	if ( verbose ) {
		cout << "Best scale:                     " << scalebest << " (" 
			<< ccbest << ")" << endl;
		cout << "Sampling correction:            " 
			<< pref->image->sampling()[0]*scalebest/image->sampling()[0] << endl << endl;
	}
	
//	vscale = Vector3<double>(scalebest, scalebest, scalebest);

//	pt = transform(pref->size(), vscale, origin, translate, mat, FILL_BACKGROUND);
//	pt->sampling(pref->sampling());
	
	sampling(pref->sampling(0) * scalebest);
	pt = fspace_resize(pref);
	
	return pt;
}

/**
@brief	Scales and image to the given reference, searching for the correct scale.
@param	pref		reference/template image.
@param	pmask		mask to limit correlation calculation (can be NULL).
@return	Bimage*		scaled image.
**/
Bimage*		Bimage::scale_to_reference(Bimage* pref, Bimage* pmask)
{
	double			step(pref->image->sampling()[0]/pref->x);
	double			scalemin(image->sampling()[0]/pref->image->sampling()[0] - 10*step);
	double			scalemax(image->sampling()[0]/pref->image->sampling()[0] + 10*step);
	
	return scale_to_reference(pref, pmask, scalemin, scalemax, step);
}

