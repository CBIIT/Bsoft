/**
@file	Bimage_symmetry.cpp
@brief	Symmetry function library for crystallography
@author Bernard Heymann
@date	Created: 19990509
@date	Modified: 20191104
**/

#include "Bimage.h"
#include "rwsymop.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Applies point group symmetry to an image.
@param 	sym			point group.
@param 	ref_view	reference view vector and rotation angle.
@param 	flag		flag to normalize after symmetrization.
@return double		symmetry FOM.

	The point group symmetry operations are applied to an image with an
	orientation defined by the reference symmetry axis (default {0,0,1}). 
	The fill value is taken from image's background value.
	The symmetry FOM is taken as the ratio of the power after to before
	symmetrization. 

**/
double 		Bimage::symmetrize(Bsymmetry sym, View ref_view, int flag)
{
	if ( sym.point() < 102 ) return 0;
	
	long			theorder(sym.order());
	long 			i, m;
	double			pwr(std*std);
	
	if ( z == 1 ) {
		if ( sym.label()[0] != 'C' ) {
			error_show("Error in Bimage::symmetrize: Only cyclic point group symmetries can be applied to 2D images!", __FILE__, __LINE__);
			return -1;
		}
		ref_view = View(0, 0, 1, ref_view.angle());
	}

	calculate_background();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying symmetry " << sym.label() << ":" << endl;
		cout << "Transformation origin:          " << image->origin() << endl;
		cout << "Reference view:                 " << ref_view << endl;
		cout << "Fill value:                     " << background(long(0)) << endl << endl;
	}
	
	change_type(Float);

	Matrix3			mat = ref_view.matrix();
	
	Vector3<double>	new_axis(ref_view.vector3());

	Bimage**		plist;
	for ( i=0; i<sym.operations(); ++i ) {
		m = sym[i].order();
		if ( verbose & VERB_PROCESS )
			cout << "Order[" << i+1 << "] = " << m << endl;
		new_axis = mat * sym[i].axis();
		plist = new Bimage*[m];
		plist[0] = copy();
#ifdef HAVE_GCD
		dispatch_apply(m-1, dispatch_get_global_queue(0, 0), ^(size_t k){
			long 	j = k + 1;
			plist[j] = rotate(size(), new_axis, j*M_PI*2.0/m);
		});
#else
#pragma omp parallel for
		for ( long j=1; j<m; j++ )
			plist[j] = rotate(size(), new_axis, j*M_PI*2.0/m);
#endif
		sum(m, plist);
		for ( long j=0; j<m; j++ ) delete plist[j];
		delete[] plist;
	}
	if ( verbose )
		cout << endl;
	
	if ( flag ) multiply(1.0L/theorder);
	
	statistics();
	calculate_background();
	
	double		f(std*std/pwr);
	if ( !flag ) f /= theorder;
	
	if ( verbose & VERB_PROCESS )
		cout << "Symmetry FOM:                   " << f <<
			" (" << n << ")" << endl << endl;
			
	return f;
}

/**
@brief 	Checks for the requested symmetries.
@param 	&check_string	string with requested point groups.
@return double			symmetry correlation coefficient.

	Requirement: The input map must be in standard orientation for the symmetry.
	For each different symmetry operation, the map is rotated and compared by
	real space correlation to the original. 

**/
double 		Bimage::check_point_group(Bstring& check_string)
{
	check_string = check_string.upper();
	
	Bstring*		sym_list = check_string.split(",");
	Bstring*		sym;
	Bstring			sym_best("none");
	
	int				order;
	double			angle, CC1, CC2, CC3, CCavg(0), CCbest(0);
	Vector3<double>	axis;
	
	cout << "Group\tAxis1\tCC1\tAxis2\tCC2\tAxis3\tCC3\tCCavg" << endl;
	for ( sym = sym_list; sym; sym = sym->next ) {
		cout << *sym << tab;
		if ( (*sym)[0] == 'C' ) {
//			sscanf(sym->c_str(), "C%d", &order);
			order = sym->substr(1,10).integer();
			axis = Vector3<double>(0,0,1);
			angle = M_PI*2.0L/order;
			CCavg = CC1 = rotate_correlate(axis, angle);
			cout << order << tab << CC1 << "\t-\t-\t-\t-\t" << CCavg << endl;
		} else if ( (*sym)[0] == 'D' ) {
//			sscanf(sym->c_str(), "D%d", &order);
			order = sym->substr(1,10).integer();
			axis = Vector3<double>(0,0,1);
			angle = M_PI*2.0L/order;
			CC1 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1,0,0);
			angle = M_PI;
			CC2 = rotate_correlate(axis, angle);
			CCavg = (CC1 + CC2)/2;
			cout << order << tab << CC1 << "\t2\t" << CC2 << "\t-\t-\t" << CCavg << endl;
		} else if ( (*sym)[0] == 'T' ) {
			axis = Vector3<double>(0,0,1);
			angle = M_PI;
			CC1 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1,0,0);
			angle = M_PI;
			CC2 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1/sqrt(3.0),1/sqrt(3.0),1/sqrt(3.0));
			angle = M_PI*2.0L/3.0L;
			CC3 = rotate_correlate(axis, angle);
			CCavg = (CC1 + CC2 + CC3)/3;
			cout << "2\t" << CC1 << "\t2\t" << CC2 << "\t3\t" << CC3 << tab << CCavg << endl;
		} else if ( (*sym)[0] == 'O' ) {
			axis = Vector3<double>(0,0,1);
			angle = M_PI/2.0L;
			CC1 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1,0,0);
			CC2 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1/sqrt(3.0),1/sqrt(3.0),1/sqrt(3.0));
			angle = M_PI*2.0L/3.0L;
			CC3 = rotate_correlate(axis, angle);
			CCavg = (CC1 + CC2 + CC3)/3;
			cout << "4\t" << CC1 << "\t4\t" << CC2 << "\t3\t" << CC3 << tab << CCavg << endl;
		} else if ( (*sym)[0] == 'I' ) {
			axis = Vector3<double>(0,0,1);
			angle = M_PI;
			CC1 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1/sqrt(3.0),1/sqrt(3.0),1/sqrt(3.0));
			angle = M_PI*2.0L/3.0L;
			CC2 = rotate_correlate(axis, angle);
			axis = Vector3<double>(1.0L/sqrt(GOLDEN*GOLDEN + 1),GOLDEN/sqrt(GOLDEN*GOLDEN + 1.0L),0);
			angle = M_PI*2.0L/5.0L;
			CC3 = rotate_correlate(axis, angle);
			CCavg = (CC1 + CC2 + CC3)/3;
			cout << "2\t" << CC1 << "\t3\t" << CC2 << "\t5\t" << CC3 << tab << CCavg << endl;
		}
		if ( CCbest < CCavg ) {
			CCbest = CCavg;
			sym_best = *sym;
		}
	}
	
	cout << endl << "Best symmetry:                  " << sym_best << endl;
	cout << "Correlation coefficient:        " << CCbest << endl << endl;
	
	return CCbest;
}

Vector3<double>	vector3_closest_to_intersection(Vector3<double> origin1,
					Vector3<double> origin2, Vector3<double> axis1, Vector3<double> axis2)
{
	Vector3<double>	a1 = axis1;
	Vector3<double>	a2 = axis2;
	double			a = a1.scalar(a1);
	double			b = a1.scalar(a2);
	double			c = a2.scalar(a2);
	double			d = a1.scalar(origin1 - origin2);
	double			e = a2.scalar(origin1 - origin2);
	
	double			s, t, denom = a*c - b*b;

	if ( denom < 1e-30 ) {		// Parallel lines
		s = 0;
		t = d/b;
	} else {
		s = (b*e - c*d)/denom;
		t = (a*e - b*d)/denom;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Intersection distance:          " << (origin1 + a1 * s).distance(origin2 + a2 * t) << endl;
	
	Vector3<double>	point = (origin1 + a1 * s + origin2 + a2 * t)/2;
	
	return point;
}

/**
@brief 	Finds the orientation for an image with a cyclic point group symmetry.
@param 	sym			point group.
@param 	binfac		binning for faster searching (limited to 1,2,3).
@param 	hires		high resolution limit in angstroms.
@param 	lores		low resolution limit in angstroms.
@return double		symmetry correlation coefficient.

	The cyclic point group symmetry operation is applied to an image using
	the reference symmetry axis, aslo the 2D rotation axis (default {0,0,1}). 

**/
double 		Bimage::find_cyclic_point_group(Bsymmetry& sym,
				int binfac, double hires, double lores)
{
	if ( sym.point() < 102 || sym.point() > 200 ) return 0;

	if ( binfac < 1 ) binfac = 1;
	if ( binfac > 4 ) binfac = 4;
	
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	if ( hires < 2*binfac*image->sampling()[0] ) hires = 2*binfac*image->sampling()[0];
	if ( lores < hires ) lores = real_size()[0];
	
	Bimage*			pb = bin_copy(binfac);

	pb->calculate_background();

	fft_plan		planf = pb->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = pb->fft_setup(FFTW_BACKWARD, 1);

	double			cc;
	Vector3<double>	axis(0,0,1), scale(1,1,1), translate;	
	Matrix3			mat(axis, sym[0].angle());

	if ( verbose ) {
		cout << endl << "Finding orientation based on symmetry " << sym.label() << ":" << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Bin factor:                     " << binfac << endl;
	}

	Vector3<double>	origin = pb->rotate_find_shift(mat,
				hires, lores, pb->sizeX()/4, 0, 1, planf, planb, cc);

	origin *= binfac;

	delete pb;
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	translate = image->origin() - origin;
	transform(scale, origin, translate, mat, FILL_BACKGROUND, 0);

	if ( verbose ) {
		cout << "Best origin:                    " << origin << " " << cc << endl;
		cout << "Best shift:                     " << translate << endl << endl;
	}
	
	statistics();

	return cc;
}

/**
@brief 	Finds the orientation for an image with a specific point group symmetry.
@param 	sym			point group.
@param 	angle_step	angular step size for search.
@param 	binfac		binning for faster searching (limited to 1,2,3).
@param 	hires		high resolution limit in angstroms.
@param 	lores		low resolution limit in angstroms.
@param 	flags		flag to search only for minor axes.
@return double		symmetry correlation coefficient.

	The point group symmetry operations are applied to an image with an
	orientation defined by the reference symmetry axis (default {0,0,1}). 

**/
double 		Bimage::find_point_group(Bsymmetry& sym, double angle_step,
				int binfac, double hires, double lores, int flags)
{
	if ( sym.point() < 102 ) return 0;

	if ( ( flags & 1 ) && sym.point() < 200 )
		return find_cyclic_point_group(sym, binfac, hires, lores);
	
	if ( binfac < 1 ) binfac = 1;
	if ( binfac > 4 ) binfac = 4;
	
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	if ( hires < 2*binfac*image->sampling()[0] ) hires = 2*binfac*image->sampling()[0];
	if ( lores < hires ) lores = real_size()[0];
	
	long			i, op0(0), op1(0), grid_points;
	double			theta(M_PI_2), phi, phi_step, angle(M_PI), bestangle(0);
	double			bestcc(-1);
	Vector3<double>	ref_axis(0,0,1), axis1, axis2, axis3, bestaxis(ref_axis), bestaxis2(1,0,0);
	Vector3<double>	bestorigin, bestorigin2, scale(1,1,1), translate(0,0,0);
	Matrix3			mat;

	Bimage*			pb = bin_copy(binfac);

	pb->calculate_background();
	
	bestorigin = bestorigin2 = pb->image->origin();

	// Select the operation with the highest order
	for ( i=0, angle = TWOPI; i<sym.operations(); i++ )
		if ( sym[i].angle() < angle ) {
			angle = sym[i].angle();
			op0 = i;
		}

	// Select the operation with the second highest order
	for ( i=0, angle = TWOPI; i<sym.operations(); i++ )
		if ( sym[i].angle() < angle && i != op0 ) {
			angle = sym[i].angle();
			op1 = i;
		}

	// Count the number of points in the grid search
	if ( z > 1 ) {
		for ( grid_points=0, theta=0; theta<=M_PI_2+0.001; theta += angle_step ) {
			phi_step = TWOPI;
			if ( theta > 0 ) phi_step = angle_step/sin(theta);
			for ( phi=0; phi<TWOPI; phi += phi_step ) grid_points++;
		}
	} else {	// 2D
		phi_step = angle_step;
		for ( grid_points=0, phi=0; phi<TWOPI; phi += phi_step ) grid_points++;
	}
	
	// Set up the first grid of axes
	int					nax = (grid_points > 9)? grid_points: 9;
	Vector3<double>*	axis = new Vector3<double>[nax];
	Vector3<double>*	origin = new Vector3<double>[nax];
	double*				cc = new double[nax];
	
	if ( z > 1 ) {
		for ( i=0, theta=0; theta<=M_PI_2+0.001; theta += angle_step ) {
			phi_step = TWOPI;
			if ( theta > 0 ) phi_step = angle_step/sin(theta);
			for ( phi=0; phi<TWOPI; phi += phi_step, i++ ) {
				axis[i][0] = cos(phi)*sin(theta);
				axis[i][1] = sin(phi)*sin(theta);
				axis[i][2] = cos(theta);
			}
		}
	} else {	// 2D
		phi_step = angle_step;
		for ( i=0, phi=0; phi<TWOPI; phi += phi_step, i++ ) {
			axis[i][0] = axis[i][1] = 0;
			axis[i][2] = 1;
		}
	}

	fft_plan		planf = pb->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = pb->fft_setup(FFTW_BACKWARD, 1);
	
	if ( verbose ) {
		cout << endl << "Finding orientation based on symmetry " << sym.label() << ":" << endl;
		cout << "Angular step size:              " << angle_step*180.0/M_PI << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Bin factor:                     " << binfac << endl;
		cout << "Grid points:                    " << grid_points << endl;
		cout << "Major axis (" << sym[op0].order() << "-fold):" << endl;
	}

	if ( verbose & VERB_LABEL )
		cout << "#\tax\tay\taz\tCC" << endl;

	long				nmax(grid_points);
	double				da(angle_step), amin, amax;
	Quaternion			q1, q, qv;
	if ( flags & 1 ) da = 0;
	while ( da > 0.00175 ) {		// 0.1 degree limit
		if ( verbose )
			cout << "Current angular step: " << da*180.0/M_PI << "\r" << flush;

#ifdef HAVE_GCD
		dispatch_apply(nmax, dispatch_get_global_queue(0, 0), ^(size_t i){
			Matrix3			mat = Matrix3(axis[i], sym[op0].angle());
			origin[i] =	pb->rotate_find_shift(mat,
				hires, lores, pb->sizeX()/4, 0, 1, planf, planb, cc[i]);
		});
#else
#pragma omp parallel for
		for ( i=0; i<nmax; i++ ) {
			Matrix3			mat = Matrix3(axis[i], sym[op0].angle());
			origin[i] =	pb->rotate_find_shift(mat,
				hires, lores, pb->sizeX()/4, 0, 1, planf, planb, cc[i]);
		}
#endif

		for ( i=0; i<nmax; i++ ) {
			if ( bestcc < cc[i] ) {
				bestcc = cc[i];
				bestaxis = axis[i];
				bestorigin = origin[i];
			}
			if ( verbose & VERB_LABEL )
				cout << i+1 << tab << axis[i] << tab << cc[i] << endl;
			if ( verbose & VERB_TIME )
				cout << "Grid points done:  " << i << " (" << i*100.0/nmax << "%)  Best: " << bestcc << "\r" << flush;
		}
		axis1 = bestaxis.cross(ref_axis);
		if ( axis1.length() < 0.0001 ) axis1[0] = 1;
		axis1.normalize();
		axis2 = bestaxis.cross(axis1);
		if ( axis2.length() < 0.0001 ) axis2[1] = 1;
		axis2.normalize();
		qv = Quaternion(0, bestaxis);
		da /= 2;				// Contract and calculate small grid
		if ( z > 1 ) {
			for ( i=0, theta = -da; theta <= da; theta += da ) {
				q1 = Quaternion(axis1, theta);
				for ( phi = -da; phi <= da; phi += da, i++ ) {
					q = Quaternion(axis2, phi);
					q *= q1;
					q = q.rotate(qv);
					axis[i] = q.axis();
				}
			}
			nmax = i;
		} else break;
	}

	delete[] axis;
	delete[] origin;
	delete[] cc;
	
	if ( verbose ) {
		cout << "Best major axis:                " << bestaxis << " " << bestcc << endl;
		cout << "Best origin:                    " << bestorigin*binfac << endl;
	}

	if ( sym.operations() > 1 ) {
		if ( verbose )
			cout << "Minor axis (" << sym[op1].order() << "-fold):" << endl;
		if ( verbose & VERB_LABEL )
			cout << "#\tax\tay\taz\tCC" << endl;
		bestcc = 0;
		mat = Matrix3(bestaxis, sym[op0].axis());
		axis2 = mat * sym[op1].axis();
		da = angle_step;
		amin = 0;
		amax = M_PI;
		
		// Count the number of points in the grid search
		for ( grid_points=0, angle=amin; angle<amax; angle+=da ) grid_points++;
	
		// Set up the first grid of axes
		nax = (grid_points > 9)? grid_points: 9;
		axis = new Vector3<double>[nax];
		origin = new Vector3<double>[nax];
		cc = new double[nax];
	
		while ( da > 0.00175 ) {		// 0.1 degree limit
			for ( angle=amin, i=0; angle<amax; angle+=da, i++ ) {
				mat = Matrix3(bestaxis, angle);
				axis[i] = mat * axis2;
			}
			nax = i;
		
#ifdef HAVE_GCD
			dispatch_apply(nax, dispatch_get_global_queue(0, 0), ^(size_t i){
				Matrix3			mat = Matrix3(axis[i], sym[op1].angle());
				origin[i] =	pb->rotate_find_shift(mat,
					hires, lores, pb->sizeX()/4, 0, 1, planf, planb, cc[i]);
			});
#else
#pragma omp parallel for
			for ( i=0; i<nax; i++ ) {
				Matrix3			mat = Matrix3(axis[i], sym[op1].angle());
				origin[i] =	pb->rotate_find_shift(mat,
					hires, lores, pb->sizeX()/4, 0, 1, planf, planb, cc[i]);
			}
#endif

			for ( angle=amin, i=0; angle<amax; angle+=da, i++ ) {
				if ( bestcc < cc[i] ) {
					bestcc = cc[i];
					bestaxis2 = axis[i];
					bestangle = angle;
					bestorigin2 = origin[i];
				}

				if ( verbose & VERB_LABEL )
					cout << i+1 << tab << axis[i] << tab << cc[i] << endl;
				if ( verbose & VERB_TIME )
					cout << "Angles done:  " << i+1 << " (" << angle*100.0/M_PI << "%)  Best: " << bestcc << "\r" << flush;
			}
			da /= 2;
			amin = bestangle - da;
			amax = bestangle + da;
		}
		bestorigin = vector3_closest_to_intersection(bestorigin, bestorigin2, bestaxis, bestaxis2);
		if ( verbose ) {
			cout << "Best minor axis:                " << bestaxis2 << " " << bestcc << endl;
			cout << "Best origin:                    " << bestorigin2*binfac << endl;
			cout << "Best intersection origin:       " << bestorigin*binfac << endl;
		}
		
		delete[] axis;
		delete[] origin;
		delete[] cc;
	}

	delete pb;
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	bestorigin *= binfac;
//	bestorigin = bestorigin*binfac + 0.5;
	
	mat = Matrix3(bestaxis, sym[op0].axis());
	translate = image->origin() - bestorigin;
	Matrix3		mat2 = Matrix3(sym[op0].axis(), -bestangle);
	if ( sym.operations() > 1 ) mat = mat2 * mat;
	transform(scale, bestorigin, translate, mat, FILL_BACKGROUND, 0);
	view(View(mat));

	if ( verbose ) {
		cout << "Best view:                      " << image->view() << endl;
		cout << "Best shift:                     " << translate << endl << endl;
	}
	
	statistics();
	
	return bestcc;
}

/**
@brief 	Rotates to a symmetry axis.
@param 	*sym		symmetry structure.
@param 	axis		desired symmetry axis order.
@param 	axis_flag	view modifier.
@return long 		0.
**/
long		Bimage::rotate_to_axis(Bsymmetry& sym, long axis, long axis_flag)
{
	Matrix3			mat = symmetry_rotate_to_axis(sym, axis, axis_flag);
		
	rotate(mat);
			
	return 0;
}

/**
@brief 	Changes the symmetry order for cyclic and dihedral maps.
@param 	*symold		original symmetry.
@param 	*symnu		new symmetry.
@param 	radius		radius in voxels of repeating units.
@param 	z_slope		slope along z to adjust radius.
@return int 		0.

	If the radius is zero, the change in symmetry is done with proportional
	shifts that distorts the units in the map.
	If the radius is positive, it is used to determine the radial shift for
	each unit, giving much less distortion than the proportional shifts.
	Variation of the radius of the specimen along z can be handled with
	a non-zero z slope parameter.

**/
int			Bimage::change_symmetry(Bsymmetry& symold, Bsymmetry& symnu,
				double radius, double z_slope)
{
	long			i, xx, yy, zz, nn;
	double			dx, dy, dz, dxo, dyo, dzo, da, r, ro(1), ca, sa, s;
	double			rr(symold[0].order()*1.0L/symnu[0].order());
	double			af(1/rr - 1);
	
	float*			nudata = new float[datasize];
	
	if ( verbose ) {
		cout << "Changing the symmetry from " << symold.label() << " to " << symnu.label() << ":" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Z-slope:                        " << z_slope << endl << endl;
	}
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			dz = zz - image[nn].origin()[2];
			s = -(radius + dz*z_slope)*af;
			dzo = zz;
			if ( dzo >= z - 1 ) dzo = z - 1.001;
			for ( yy=0; yy<y; yy++ ) {
				dy = yy - image[nn].origin()[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					nudata[i] = 0;
					dx = xx - image[nn].origin()[0];
					da = af*atan2(dy, dx);
					ca = cos(da);
					sa = sin(da);
					if ( radius > 0 ) {
						r = sqrt(dx*dx + dy*dy);
						if ( r > 0 ) ro = (1 + s/r);
						else ro = 1;
						dxo = ro*(dx*ca - dy*sa) + image[nn].origin()[0];
						dyo = ro*(dx*sa + dy*ca) + image[nn].origin()[1];
					} else {
						dxo = rr*(dx*ca - dy*sa) + image[nn].origin()[0];
						dyo = rr*(dx*sa + dy*ca) + image[nn].origin()[1];
					}
					if ( ro>0 && dxo>=0 && dxo<=x-1 && dyo>=0 && dyo<=y-1 )
						nudata[i] = interpolate(dxo, dyo, dzo, nn, background(nn));
				}
			}
		}
	}
	
	data_type(Float);
	
	data_assign((unsigned char *) nudata);

	statistics();
	
	return 0;
}

/**
@brief 	Finds the symmetry equivalent orientation of a particle with respect to a template.
@param 	*pref			template to test against (same size as input image).
@param 	*pmask			mask to limit correlation calculation.
@param 	*sym			symmetry structure.
@return Matrix3			transformed image.

	For point groups up to tetrahedral symmetry, there is an ambiguity in the
	orientation of a particle even though it is oriented according to the 
	standard orientation for the symmetry.
	A template is used to select the desired orientation and the image is transformed.
	Both the input image and template must be in the standard orientation
	for the point group.

**/
Matrix3		Bimage::symmetry_equivalent(Bimage* pref, Bimage* pmask, Bsymmetry& sym)
{
	Matrix3			mat(1);
	if ( sym.point() < 102 ) return mat;
	if ( sym.point() > 600 ) return mat;
	if ( !pref ) return mat;

	if ( image->origin()[0] < 1 || image->origin()[1] < 1 || image->origin()[2] < 1 ) origin(size()/2);	
	
	if ( sym.point() < 200 )
		return symmetry_equivalent_cyclic(pref, pmask, sym);
	
	long			nstep(1), step;
	double			dang(0), angle(0);
	
	if ( sym.point() < 200 ) {							// Cyclic
		dang = M_PI/180.0;								// 1 degree step size
		nstep = 360/sym[0].order();
	} else if ( sym.point() == 202 ) {					// D2
		nstep = 6;
	} else if ( sym.point() < 300 ) {					// Dihedral
		dang = M_PI*1.0L/sym[0].order();
		if ( sym[0].order()%2 == 0 ) nstep = 2;
		else nstep = 4;
	} else if ( sym.point() < 400 ) {					// Tetrahedral
		dang = M_PI_2;
		nstep = 2;
	} else if ( sym.point() < 500 ) {					// Octahedral
		dang = M_PI_2;
		nstep = 1;
	} else {											// Icosahedral
		dang = M_PI_2;
		nstep = 4;
	}
	
	double			scale(sampling(0)[0]/pref->sampling(0)[0]);
	Vector3<double>	axis(0,0,1), best_axis(0,0,1);
	Vector3<double> vscale(scale, scale, scale);
	Vector3<double> ori(image->origin()), translate((pref->size() - size())/2);
	Bimage*			pt;

	double			ccbest(-1);
	double			cc, best_angle(0);
	Quaternion		q(axis, M_PI_2), q2;
	
	if ( verbose ) {
		cout << "Finding the symmetry equivalent for " << sym.label() << endl;
		cout << "Axis\t\t\tAngle\tCC" << setprecision(4) << endl;
	}

	for ( step = 0; step < nstep; ++step ) {
		if ( sym.point() == 202 ) {						// D2
			switch ( step ) {
				case 1: angle = M_PI_2; axis = Vector3<double>(0,0,1); break;
				case 2: angle = M_PI_2; axis = Vector3<double>(0,1,0); break;
				case 3: angle = M_PI_2; axis = Vector3<double>(1,0,0); break;
				case 4: angle = M_PI/1.5L; axis = Vector3<double>(1,1,1); break;
				case 5: angle = -M_PI/1.5L; axis = Vector3<double>(1,1,1); break;
				default: angle = 0;
			}
			axis.normalize();
		} else if ( sym.point() < 300 ) {				// Dihedral
			switch ( step ) {
				case 1: angle = dang; break;
				case 2: angle = M_PI_2; break;
				case 3: angle = M_PI_2 - dang; break;
				default: angle = 0;
			}
			axis.normalize();
		} else if ( sym.point() == 532 ) {				// Icosahedral
			switch ( step ) {
				case 1: angle = dang; break;
				case 2: angle = atan(1/GOLDEN); axis = Vector3<double>(1,0,0); break;
				case 3: q2 = Quaternion(axis, angle); q = q2*q;
					axis = q.axis(); angle = q.angle(); break;
				default: angle = 0;
			}
			axis.normalize();
		} else {										// Cyclic, tetra-, octahedral
			angle = step*dang;
		}
		mat = Matrix3(axis, angle);
		pt = transform(pref->size(), vscale, ori, translate, mat, FILL_BACKGROUND);
		cc = pt->correlate(pref, 0, pref->sizeX()/2, pmask);
		if ( ccbest < cc ) {
			ccbest = cc;
			best_angle = angle;
			best_axis = axis;
//			if ( verbose & VERB_PROCESS )
//				cout << "Better CC:                      " << cc << " (" << angle*180.0/M_PI << ")" << endl;
		}
		delete pt;
		if ( verbose )
			cout << axis << tab << angle*180.0/M_PI << tab << cc << endl;
	}
	
	if ( verbose ) {
		cout << "Best axis and angle:            " << best_axis << tab << best_angle*180.0/M_PI << endl;
		cout << "Best CC:                        " << ccbest << endl;
		cout << "Fixing symmetry " << sym.label() << " orientation" << endl << endl;
	}
	
	mat = Matrix3(best_axis, best_angle);

	return mat;
}

/**
@brief 	Finds the cyclic symmetry equivalent orientation of a particle with respect to a template.
@param 	*pref			template to test against (same size as input image).
@param 	*pmask			mask to limit correlation calculation.
@param 	*sym			symmetry structure.
@return Matrix3			transformed image.

	For point groups up to tetrahedral symmetry, there is an ambiguity in the
	orientation of a particle even though it is oriented according to the 
	standard orientation for the symmetry.
	A template is used to select the desired orientation and the image is transformed.
	Both the input image and template must be in the standard orientation
	for the point group.

**/
Matrix3		Bimage::symmetry_equivalent_cyclic(Bimage* pref, Bimage* pmask, Bsymmetry& sym)
{
	Matrix3			mat(1);
	if ( sym.point() < 102 ) return mat;
	if ( sym.point() > 200 ) return mat;	// Not suitable for groups other than cyclic
	if ( !pref ) return mat;
	
	if ( image->origin()[0] < 1 || image->origin()[1] < 1 || image->origin()[2] < 1 ) origin(size()/2);	

	double			hires(2*pref->sampling(0)[0]);
	if ( hires < 20 ) hires = 20;		// Starting resolution

	double			accuracy(M_PI*0.1/180.0);		// Angular accuracy ~ 0.1 degrees
	if ( accuracy < 3.0/pref->sizeX() ) accuracy = 3.0/pref->sizeX();	// or 2/3 Nyquest
	
	double			dang(hires/pref->real_size()[0]);		// Initial angular step size
	double			angle, angmin(0), angmax(TWOPI/sym[0].order()+0.01);
	double			scale(sampling(0)[0]/pref->sampling(0)[0]);
	Vector3<double>	axis(0,0,1);
	Vector3<double> vscale(scale, scale, scale);
	Vector3<double> ori(image->origin()), translate((pref->size() - size())/2);

	double			ccbest(-1), cc, best_angle(0);
	Vector3<double>	shift, best_shift;
	
	if ( verbose )
		cout << "Finding the symmetry equivalent for " << sym.label() << setprecision(4) << endl;

	Bimage*			pt = transform(pref->size(), vscale, ori, translate, mat, FILL_BACKGROUND);

//	long 			nannuli(pref->x/2);
	Bimage*			ptcyl = pt->symmetrize_cylinder(1);
	Bimage*			prcyl = pref->symmetrize_cylinder(1);
	
	delete pt;
	
	best_shift = ptcyl->find_shift(prcyl, NULL, 4*pref->image->sampling()[0], 0, pref->x/4, 0, 1, ccbest);

	if ( verbose ) {
		cout << "Shift\t\t\tAngle\tCC" << endl;
		cout << best_shift << tab << 0 << tab << ccbest << endl;
	}
	
	ptcyl->reslice("-xyz");

	shift = ptcyl->find_shift(prcyl, NULL, 4*pref->image->sampling()[0], 0, pref->x/4, 0, 1, cc);

	if ( verbose )
		cout << shift << tab << 180 << tab << cc << endl;

	if ( ccbest < cc ) {
		ccbest = cc;
		best_shift = shift;
		if ( verbose )
			cout << "Map is upside down - rotating 180 degrees around x axis" << endl << endl;
		axis = Vector3<double>(1,0,0);
		rotate(axis, M_PI);
		axis = Vector3<double>(0,0,1);
//		write_img("t.pif", this, 0);
	}
	
	delete ptcyl;
	delete prcyl;
	
	while ( dang >= accuracy ) {
		ccbest = -1;
		cout << "∆ang = " << dang*180.0/M_PI << tab << "res = " << hires << endl;
		cout << "Shift\t\t\tAngle\tCC" << endl;
		for ( angle = angmin; angle <= angmax; angle += dang ) {
			mat = Matrix3(axis, angle);
			pt = transform(pref->size(), vscale, ori, translate, mat, FILL_BACKGROUND);
//			cc = pt->correlate(pref, 0, pref->sizeX()/2, pmask);
			if ( pmask ) pt->multiply(pmask);
			shift = pt->find_shift(pref, NULL, hires, 0, 10, 0, 1, cc);
			if ( ccbest < cc ) {
				ccbest = cc;
				best_angle = angle;
				best_shift = shift;
//				if ( verbose & VERB_PROCESS )
//					cout << "Better CC:                      " << cc << " (" << angle*180.0/M_PI << ")" << endl;
			}
			delete pt;
			if ( verbose )
				cout << shift << tab << angle*180.0/M_PI << tab << cc << endl;
		}
		dang *= 0.75;
		angmin = best_angle - 2*dang;
		angmax = best_angle + 2.1*dang;
		hires = dang*pref->real_size()[0];
	}
	
	best_shift *= scale;
	
	if ( verbose ) {
		cout << "Best shift and angle:           " << best_shift << tab << best_angle*180.0/M_PI << endl;
		cout << "Best CC:                        " << ccbest << endl;
		cout << "Fixing symmetry " << sym.label() << " orientation" << endl << endl;
	}
	
	mat = Matrix3(axis, best_angle);
	origin(size()/2 + best_shift);

	return mat;
}

/**
@brief 	Calculates a multi-level mask to indicate asymmetric units.
@param 	*sym		symmetry structure.
@param 	index		asu index (<1 means all).
@return Bimage*		multi-level mask.

	A reference view for each point group is generated such that it is
	located close to the center of the canonical asymmetric unit.
	A set of symmetry-related reference views are then generated from 
	the original reference.
	The distance of each voxel is then calculated to each reference
	view, and assigned to asymmetric unit with the closest reference view.
	If the given index is larger than one, a mask is generated only for the
	asymmetric unit with that index. 

**/
Bimage*		Bimage::levelmask_asymmetric_units(Bsymmetry& sym, int index)
{
	View			ref = view_symmetry_reference(sym);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating an asymmetric unit mask:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Asymmetric unit:                ";
		if ( index > 0 ) cout << index << endl;
		else cout << "all" << endl;
		cout << "Reference vector:               " << ref << endl << endl;
	}
	
	Bimage* 		pmask = copy_header();
	pmask->data_type(UCharacter);
	pmask->data_alloc_and_clear();
	
	if ( sym.point() < 102 ) {
		pmask->fill(1);
		pmask->statistics();
		return pmask;
	}
	
	long			i, j, m, nn, xx, yy, zz;
	double			dist, mindist;
	Vector3<double>	testvec, symvec;
	View*			v;
	View*			symview = symmetry_get_all_views(sym, ref);

	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			testvec[2] = (double)zz - pmask->image[nn].origin()[2];
			for ( yy=0; yy<y; yy++ ) {
				testvec[1] = (double)yy - pmask->image[nn].origin()[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					testvec[0] = (double)xx - pmask->image[nn].origin()[0];
					mindist = 1e30;
					for ( j=m=1, v=symview; v; v=v->next, j++ ) {
						symvec = v->vector3();
						dist = testvec.distance(symvec);
						if ( mindist > dist ) {
							mindist = dist;
							m = j;
						}
					}
					if ( index < 1 ) pmask->set(i, m);
					else if ( index == m ) pmask->set(i, 1);
				}
			}
		}
	}
	
	kill_list((char *) symview, sizeof(View));
	
	pmask->statistics();
	
	return pmask;
}

/**
@brief 	Calculates a full map from one asymmetric unit.
@param 	*sym		symmetry structure.
@return int			0.

	A reference view for each point group is generated such that it is
	located close to the center of the canonical asymmetric unit.
	For each voxel in the target map, a set of symmetry-related views 
	are generated and the one closest to the reference view used
	to determine the corresponding voxels within the asymmetric unit.
	The new voxel value is calculated by trilinear interpolation of the
	voxels in the asymmetric unit. 

**/
int			Bimage::replicate_asymmetric_unit(Bsymmetry& sym)
{
	change_type(Float);
	
	View			ref = view_symmetry_reference(sym);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Replicating an asymmetric unit:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Reference vector:               " << ref << endl << endl;
	}
	
	if ( sym.point() < 102 ) return 0;
	
	long			i, j, cc, nn, xx, yy, zz, nsym(0);
	double			da, minda;
	Vector3<double>	vr(ref.vector3());
	Vector3<double>	v, vt, va;
//	Matrix3*		m = symmetry_get_all_matrices(sym, nsym);
	vector<Matrix3>	m = sym.matrices();
	nsym = m.size();

	float*			nudata = new float[datasize];

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
					for ( cc=0; cc<c; cc++, i++ )
						nudata[i] = interpolate(cc, va, nn, background(nn));
				}
			}
		}
	}
	
//	delete[] m;
	
	data_assign((unsigned char *) nudata);
	
	statistics();
	
	return 0;
}

/**
@brief 	Finds the view that on symmetrizing fits best to a symmetric template.  
@param 	*ptemp			symmetric template.
@param 	&sym			symmetry.
@param 	phi_step		phi angle step size (radians).
@param 	theta_step		theta angle step size (radians).
@param 	alpha_step		rotation angle step size (radians).
@param 	shift    		shift to impose before symmetrization.
@return Bimage* 		image rotated to the best view.

	The orientation parameters, view vector, angle of rotation and origin,
	of each image is packed into 3D reciprocal space.
	An image is used in the reconstruction if its selection flag has been set.
	The fill value is taken from image's background value.  

**/
Bimage*		Bimage::find_symmetric_view(Bimage* ptemp, Bsymmetry& sym,
				double phi_step, double theta_step, double alpha_step, Vector3<double> shift)
{
	if ( verbose ) {
		cout << "Finding the best symmetric view:" << endl;
		cout << "Theta step size:                " << theta_step*180.0/M_PI << endl;
		cout << "Phi step size:                  " << phi_step*180.0/M_PI << endl;
		cout << "Alpha step size:                " << alpha_step*180.0/M_PI << endl;
		cout << "Shift:                          " << shift << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
	}
	
	Bsymmetry		symC1;
	View*			views = asymmetric_unit_views(symC1, theta_step, phi_step, 1);
	
	long			nv = (long) (TWOPI/alpha_step)*count_list((char *)views);
	
	if ( verbose )
		cout << "Number of views:                " << nv << endl;
	
	int				i(0);
	View			best_view;
	View*			v;
	Vector3<double>	best_shift;
	double			a, CC, bestCC(0);
	Bimage*			prot;
	
	if ( verbose )
		cout << "#\tvx\tvy\tvz\tva\tCC" << endl;
	for ( v = views; v; v = v->next ) {
		for ( a = 0; a < TWOPI; a += alpha_step ) {
			i++;
			(*v)[3] = a;
			prot = rotate(size(), shift, *v);
			prot->symmetrize(sym, 1);
			CC = ptemp->correlate(prot);
			if ( bestCC < CC ) {
				bestCC = CC;
				best_view = *v;
				best_shift = shift;
			}
			if ( verbose & VERB_LABEL )
				cout << i << tab << v << tab << CC << endl;
			delete prot;
		}
	}
	
	kill_list((char *) views, sizeof(View));

	if ( verbose )
		cout << "Best view:\t" << best_view << tab << CC << endl;

	prot = rotate(size(), shift, best_view);

	return prot;
}

