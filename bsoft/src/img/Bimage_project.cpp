/**
@file	Bimage_project.cpp
@brief	Functions for projections
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20190207
**/

#include "Bimage.h"
#include "linked_list.h"
#include "symmetry.h"
#include "Complex.h"
#include "Matrix3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Projects a 3D image to a 2D image down one of the three major axes.
@param 	axis		axis of projection.
@param 	flags		1=scale projection, 2=minimum, 4=maximum.
@return Bimage*		projection image (floating point).

	The sums of the z-planes are accumulated into a new floating point data 
	block.  This block is then rescaled and converted back to the original 
	data type.
	The new data replaces the old data.

**/
Bimage*		Bimage::project(char axis, int flags)
{
	bool				scale_flag(flags & 1);
	bool				mnproj(flags & 2);
	bool				mxproj(flags & 4);
    long				i, j, xx, yy, zz, cc, nn;
	double				v, scale;
	Vector3<double>		u(image->sampling());

	Bimage*				proj = NULL;
	if ( axis == 'x' ) {
		proj = new Bimage(Float, compoundtype, z, y, 1, n);
		scale = 1.0L/x;
		proj->sampling(u[2], u[1], 1);
	} else if ( axis == 'y' ) {
		proj = new Bimage(Float, compoundtype, x, z, 1, n);
		scale = 1.0L/y;
		proj->sampling(u[0], u[2], 1);
	} else {
		proj = new Bimage(Float, compoundtype, x, y, 1, n);
		scale = 1.0L/z;
		proj->sampling(u[0], u[1], 1);
	}
	
    long			   ds(n*x*y*c);
    
    if ( mxproj ) proj->fill(min);
    else if ( mnproj ) proj->fill(max);
	
	if ( verbose & VERB_LABEL )
    	cout << "Projecting along " << axis <<  " (" << flags << ")" << endl << endl;
 	
	// First calculate projection
	for ( i=nn=0; nn<n; nn++ ) {
		if ( axis == 'x' ) {
			proj->image[nn].origin(image[nn].origin()[2], image[nn].origin()[1], 0.0);
		} else if ( axis == 'y' ) {
			proj->image[nn].origin(image[nn].origin()[0], image[nn].origin()[2], 0.0);
		} else {
			proj->image[nn].origin(image[nn].origin()[0], image[nn].origin()[1], 0.0);
		}
		for ( zz=0; zz<z; zz++ )
			for ( yy=0; yy<y; yy++ )
				for ( xx=0; xx<x; xx++ )
					for ( cc=0; cc<c; cc++, i++ ) {
						if ( axis == 'x' )
							j = ((nn*y + yy)*z + zz)*c + cc;
						else if ( axis == 'y' )
							j = ((nn*z + zz)*x + xx)*c + cc;
						else
							j = ((nn*y + yy)*x + xx)*c + cc;
						v = (*this)[i];
						if ( mxproj ) {
							if ( (*proj)[j] < v ) proj->set(j, v);
						} else if ( mnproj ) {
							if ( (*proj)[j] > v ) proj->set(j, v);
						} else {
							proj->add(j, v);
						}
					}
	}
	
	if ( scale_flag ) for ( i=0; i<ds; i++ ) proj->set(i, (*proj)[i] * scale);
	
	proj->statistics();
	
	return proj;
}

/**
@brief 	Rotates a 3D map and projects it along the z-axis.
@param 	mat				3x3 rotation or skewing matrix.
@param 	translate		3-value vector for translation after transformation.
@param 	radial_cutoff	spherical cutoff to apply.
@param	norm_flag		flag to normalize projection.
@return Bimage*			new 2D projection image.

	The 3D map is rotated around its center as origin and the data 
	integrated along the z-direction after subtraction of the
	background (which must be calculated before this function).
	The resultant 2D image is translated if the translation vector is 
	non-zero. A radial cutoff can be applied to decrease the computation 
	time. A value of zero or less sets the default to the length of the 
	x-axis (i.e., no effective cutoff).
	The rotation origin is obtained from the map origin.

**/
Bimage*		Bimage::rotate_project(Matrix3 mat, Vector3<double> translate,
					double radial_cutoff, int norm_flag)
{
	long			i, j, xx, yy, zz;
	long			xo, yo, xp, yp, yw;
	double			oldz_sq, oldyz_sq;
	double			xf, yf, maxlen(0);
	Vector3<double>	old, newvec;
	
	if ( radial_cutoff < 1 ) radial_cutoff = x;
	double			rad_sq = radial_cutoff*radial_cutoff;
	
	Vector3<double>	oldorigin(image->origin());

	Vector3<double>	neworigin = oldorigin + translate;
	neworigin[2] = 0;
		
	Matrix3			matinv = mat.transpose();
	
	if ( verbose & VERB_FULL ) {
		cout << "Rotating around the map center and projecting:" << endl;
		cout << "Rotation origin:                " << oldorigin << endl;
		cout << matinv << endl;
		if ( translate[0] != 0 || translate[1] != 0 || translate[2] != 0 )
			cout << "Translation:                    " << translate << endl;
		cout << "Radial cutoff:                  " << radial_cutoff << " pixels" << endl;
		cout << endl;
	}

	Bimage*			proj = new Bimage(Float, TSimple, x, y, 1, 1);
	proj->sampling(image->sampling()[0], image->sampling()[1], 1);
	proj->origin(neworigin[0], neworigin[1], 0);
	
	double			bkg(background(long(0)));
	double			value;
	long			slice_size(x*y);
	
	for ( zz=0; zz<z; zz++ ) {
		old[2] = zz - oldorigin[2];
		oldz_sq = old[2]*old[2];
		if ( oldz_sq <= rad_sq ) for ( yy=0; yy<y; yy++ ) {
			old[1] = yy - oldorigin[1];
			oldyz_sq = oldz_sq + old[1]*old[1];
			if ( oldyz_sq <= rad_sq ) for ( xx=0; xx<x; xx++ ) {
				old[0] = xx - oldorigin[0];
				if ( old[0]*old[0] + oldyz_sq <= rad_sq ) {
					newvec = (matinv * old) + neworigin;
					yo = (long) newvec[1];
					if ( ( yo >= 0 ) && ( yo < y - 1 ) ) {
						xo = (long) newvec[0];
						if ( ( xo >= 0 ) && ( xo < x - 1 ) ) {
							if ( maxlen < fabs(newvec[2]) ) maxlen = fabs(newvec[2]);
							yf = newvec[1] - yo;
							xf = newvec[0] - xo;
							i = (zz*y + yy)*x + xx;
							value = (*this)[i] - bkg;
							for ( yp=yo; yp<yo+2; yp++ ) {
								yw = yp*x;
								yf = 1 - yf;
								for ( xp=xo; xp<xo+2; xp++ ) {
									j = yw + xp;
									xf = 1 - xf;
									proj->add(j, value*xf*yf);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if ( norm_flag )
		for ( i=0; i<slice_size; i++ ) proj->set(i, (*proj)[i] / maxlen);
	
	proj->statistics();
	
	return proj;
}

/**
@brief 	Calculates a set of projections from a 3D density map.
@param 	*view			linked list of views.
@param	norm_flag		flag to normalize projection.
@return Bimage* 		projections as sub-images.

	A set of projections is calculated according to a list of views.

**/
Bimage* 	Bimage::project(View* view, int norm_flag)
{
	if ( !view ) {
		error_show("Error in Bimage::project: No views defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	calculate_background();
	
	long			i;
	Vector3<double>	translate;

	long			nviews = count_list((char *)view);
	Vector3<double>	u(image->sampling());
	u[2] = 1;
	
	View*			vlist = new View[nviews];
	View*			v;
	
	for ( v=view, i=0; v; v = v->next, i++ ) vlist[i] = *v;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating projections:" << endl;
		cout << "Number of projections:          " << nviews << endl;
		if ( norm_flag ) cout << "with normalization" << endl;
	}

	Bimage*			proj = new Bimage(Float, TSimple, x, y, 1, nviews);
	proj->file_name(file_name());

#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t i){
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	one_proj = rotate_project(mat, translate, x/2.0, norm_flag);
		one_proj->shift_background(0);
		proj->replace(i, one_proj);
		proj->origin(i, one_proj->image->origin());
		proj->sampling(i, u);
		proj->image[i].view(vlist[i]);
		delete one_proj;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nviews; i++ ) {
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	one_proj = rotate_project(mat, translate, x/2.0, norm_flag);
		one_proj->shift_background(0);
		proj->replace(i, one_proj);
		proj->origin(i, one_proj->image->origin());
		proj->sampling(i, u);
		proj->image[i].view(vlist[i]);
		delete one_proj;
	}
#endif

	delete[] vlist;
	
	proj->statistics();
	
//	cout << proj->meta_data() < endl;
	
	return proj;
}

/**
@brief 	Calculates a central section of a 3D fourier transform.
@param 	mat				3x3 rotation or skewing matrix.
@param 	resolution		high resolution limit.
@param 	*kernel			frequency space interpolation kernel.
@return Bimage*			new 2D central section.

	The orientation of the central section is defined by a rotation matrix
	and the interpolation is done with a reciprocal space kernel.
	The rotation origin is obtained from the map origin.

**/
Bimage*		Bimage::central_section(Matrix3 mat, double resolution, FSI_Kernel* kernel)
{
	if ( image->sampling()[0] < 1e-10 ) sampling(1,1,1);
	if ( resolution < 1e-10 ) resolution = image->sampling()[0];
	
	long		 	i, xx, yy;
	long			hx((x - 1)/2), hy((y - 1)/2);
	Vector3<double>	m, iv;
	double			maxrad(real_size()[0]/resolution), maxrad2(maxrad*maxrad), rx2, ry2;
	double			zscale(z*1.0L/x);

	if ( verbose & VERB_FULL ) {
		cout << "Calculating a central section:" << endl;
		cout << "Resolution limit:               " << resolution << " A (" << maxrad << ")" << endl;
		cout << mat << endl << endl;
	}

	Bimage*			proj = new Bimage(Float, TComplex, x, y, 1, 1);
	proj->fourier_type(Standard);
	proj->image[0] = image[0];
	proj->sampling(image->sampling()[0], image->sampling()[1], 1);
	proj->origin(image->origin()[0], image->origin()[1], 0);
	
	for ( i=yy=0; yy<y; yy++ ) {
		iv[1] = yy;
		if ( iv[1] > hy ) iv[1] -= y;
		ry2 = iv[1]*iv[1];
		for ( xx=0; xx<x; xx++, i++ ) {
			iv[0] = xx;
			if ( iv[0] > hx ) iv[0] -= x;
			rx2 = iv[0]*iv[0];
			if ( rx2 + ry2 <= maxrad2 ) {
				m = mat * iv;
				m[2] *= zscale;
				proj->set(i, fspace_interpolate(0, m, kernel));
			}
		}
	}
		
	proj->statistics();
	
	return proj;
}

/**
@brief 	Calculates a set of projections as central sections from a 3D fourier transform.
@param 	*view			linked list of views.
@param 	resolution		high resolution limit.
@param 	*kernel			frequency space interpolation kernel.
@return Bimage* 		projections as sub-images.

	The map is Fourier transformed and shifted to its phase origin.
	For each view, a central section is calculated using reciprocal space interpolation.
	All the projections are phase shifted to a central origin and back-transformed.

**/
Bimage*     Bimage::project(View* view, double resolution, FSI_Kernel* kernel)
{
	if ( !view ) {
		error_show("Error in Bimage::project: No views defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::project: transforming and shifting origin" << endl;
	
	fft();

	phase_shift_to_origin();

//	friedel_check();
	
	long			i;
	long			nviews = count_list((char *)view);
	View*			vlist = new View[nviews];
	View*			v;
	
	for ( v=view, i=0; v; v = v->next, i++ ) vlist[i] = *v;
		
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating projections:" << endl;
		cout << "Number of projections:          " << nviews << endl;
	}

	Bimage*			proj = new Bimage(Float, TComplex, x, y, 1, nviews);
	proj->file_name(file_name());
	proj->sampling(image->sampling()[0], image->sampling()[1], 1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::project: origin after phase shift: " << image->origin() << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t i){
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	psec = central_section(mat, resolution, kernel);
		psec->phase_shift_to_center();
		proj->replace(i, psec);
		proj->image[i].view(vlist[i]);
		delete psec;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nviews; i++ ) {
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	psec = central_section(mat, resolution, kernel);
		psec->phase_shift_to_center();
		proj->replace(i, psec);
		proj->image[i].view(vlist[i]);
		delete psec;
	}
#endif

	delete[] vlist;
	
//	for ( i=0; i<proj->images(); i++ ) cout << proj->image[i].origin() << endl;

	proj->friedel_check();
	
    proj->fft_back();
	
	proj->statistics();
	
	return proj;
}

/*
@author Bernard Heymann and Peter Leong
@brief 	Back-projects a 2D image into a 3D volume.
@param 	*p			a 2D image.
@param 	resolution	high resolution limit (angstrom).
@param 	axis		tilt axis for single tilt series (radians).
@param 	planf		2D forward Fourier transform plan.
@param 	planb		2D backward Fourier transform plan.
@return int 		0.

	Requirements:
		The image must be 2D with the orientation parameters in
			the sub-image structure.
		The orientation parameters are defined as those giving the projection,
			and is reversed for backprojection.
	The 2D image is converted to floating point and rescaled to an average
	of zero and a standard deviation of one. The 3D map is traversed and
	the density added from the 2D image projected according to the orientation
	parameters. The origin of the new map is in the center of the volume.

**/
int 		Bimage::back_project(Bimage* p, double resolution, double axis, 
				fft_plan planf, fft_plan planb)
{
	if ( !p->data_pointer() ) return 0;
	
	p->change_type(Float);
	if ( axis ) p->fspace_weigh_ramp(resolution, axis, planf, planb);
	else p->fspace_weigh_ramp(resolution, planf, planb);
	
	p->statistics();
	p->rescale_to_avg_std(0, 1);
	p->calculate_background();
	
	long			i, xx, yy, zz;
	double			fill(p->background(long(0)));
	Vector3<double>	oldorigin(p->image->origin()), neworigin(image->origin());
	Vector3<double>	oldcoor, newcoor;
	Vector3<double> scale(p->image->magnification(), p->image->magnification(), 1);
			
	// The rotation matrix is transposed because the orientation parameters
	// give the orientation as rotated from the reference view
	Matrix3			mat = p->image->view().matrix();
	mat = mat.transpose();
	mat = scale*mat;
//	cout << mat << endl;
	
	if ( verbose & VERB_LABEL )
		cout << "Backprojecting" << endl << endl;
	
	for ( i=zz=0; zz<z; zz++ ) {
		newcoor[2] = zz - neworigin[2];
		for ( yy=0; yy<y; yy++ ) {
			newcoor[1] = yy - neworigin[1];
			for ( xx=0; xx<x; xx++, i++ ) {
				newcoor[0] = xx - neworigin[0];
				oldcoor = mat * newcoor + oldorigin;
				oldcoor[2] = 0;
				add(i, p->interpolate(oldcoor, 0, fill));
			}
		}
	}
	
	return 0;
}

