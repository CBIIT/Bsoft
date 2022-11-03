/**
@file	Bimage_helix.cpp
@brief	Dealing with images with helical symmetry.
@author Bernard Heymann
@date	Created: 20021127
@date	Modified: 20160331
**/

#include "Bimage.h"
#include "mg_reconstruct.h"
#include "matrix_linear.h"
#include "symmetry.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

double		Bimage::helix_interpolate(long i, double helix_rise, double helix_angle,
				int zmin, int zmax, double radius, int norm_flag)
{
	long			xx, yy, zz, nn;
	coordinates(i, xx, yy, zz, nn);
	
	Vector3<double>	d = Vector3<double>(xx - image[nn].origin()[0], yy - image[nn].origin()[1], 0);
	double			r = d.length();
	
	if ( r > radius ) return 0;
	
	long			num(0), u;
	double			pix_rise(helix_rise/image->sampling()[2]);
	double			a, cosang, sinang, value(0);
	double			xs, ys, zs;
	long			umin = (long)((zmin-(float)zz)/pix_rise);
	if ( zz < zmin ) umin += 1;
	long			umax = (long)((zmax-(float)zz)/pix_rise);
	if ( zz > zmax ) umax -= 1;
	
	for ( u = umin; u <= umax; u++ ) {
		a = u*helix_angle;
		cosang = cos(a);
		sinang = sin(a);
		zs = zz + u*pix_rise;
		xs = d[0]*cosang - d[1]*sinang + image[nn].origin()[0];
		ys = d[0]*sinang + d[1]*cosang + image[nn].origin()[1];
		if ( within_boundaries((long)xs, (long)ys, (long)zs) ) {
			value += interpolate(xs, ys, zs, nn, 0);
			num += 1;
		}
	}
	
	if ( norm_flag && num ) value /= num;
	
	return value;
}

double		Bimage::dyad_interpolate(long i, int norm_flag)
{
	long			xx, yy, zz, nn;
	coordinates(i, xx, yy, zz, nn);
	
	double			value = (*this)[i];
	double			xd(xx);
	double			yd(2 * image[nn].origin()[1] - (double)yy);
	double			zd(2 * image[nn].origin()[2] - (double)zz);
	
	value += interpolate(xd, yd, zd, nn, background(nn));
	
	if ( norm_flag ) value /= 2;
	
	return value;
}

double		Bimage::tube_interpolate(long i, int h, int k,
				double latconst, int zmin, int zmax, double radius, int norm_flag)
{
	long			xx, yy, zz, nn;
	coordinates(i, xx, yy, zz, nn);
	
	Vector3<double>	d = Vector3<double>(xx - image[nn].origin()[0], yy - image[nn].origin()[1], 0);
	double			r = d.length();
	
	if ( r > radius ) return 0;
	
	long			num(0);
	double			da(atan2(d[1], d[0]));
	double			pix_rise(latconst/image->sampling()[2]);
	
	double			t(sqrt(h*h + h*k + k*k));
	double			a(atan2(k*sqrt(3.0)/2, h+k/2.0));
	Vector3<double>	u(cos(a), -sin(a), 0), v(cos(M_PI/3-a), sin(M_PI/3-a), 0);
	Vector3<double>	w(v - u);
	Vector3<double>	lat, start, loc;
	double			value(0), bottom(zmin/pix_rise), top(zmax/pix_rise);
	
	start[1] = zz/pix_rise;
	while ( start[1] < bottom ) start += w;
	while ( start[1] > bottom ) start -= w;
	
	while ( lat[1] < top ) {
		for ( lat = start; lat[0] < t + start[0]; lat += u ) {
			if ( lat[1] >= bottom && lat[1] <= top ) {
				a = da + TWOPI * lat[0]/t;
				loc[0] = r*cos(a) + image[nn].origin()[0];
				loc[1] = r*sin(a) + image[nn].origin()[1];
				loc[2] = lat[1] * pix_rise;
				value += interpolate(loc, nn, 0);
				num++;
			}
		}
		start += w;
	}
	
	if ( norm_flag && num ) value /= num;
	
	return value;
}

/**
@brief 	Symmetrizes an image given helical symmetry parameters.
@param 	helix_rise		rise per asymmetric unit (angstrom).
@param 	helix_angle		rotation angle per asymmetric unit (radians).
@param 	dyad_axis		dyad axis indicator: 2=dyad axis on x-axis, otherwise none.
@param 	zmin			mimimum z slice to include.
@param 	zmax			maximum z slice to include.
@param 	radius			radius to do symmetrizing over (pixels).
@param 	norm_flag		if 1, normalize
@return double	 		R factor.

	The data between the z limits are replicated along the helical axis
	according to the helical rise and rotation to fill the new volume.

**/
double 		Bimage::helix_symmetrize(double helix_rise, double helix_angle,
				int dyad_axis, int zmin, int zmax, double radius, int norm_flag)
{
	change_type(Float);
	
	if ( helix_rise < 0 ) {
		helix_rise = -helix_rise;
		helix_angle = -helix_angle;
	}

	if ( radius <= 0 ) radius = x/2;

	if ( zmin < 0 ) zmin = 0;
	if ( zmax >= z ) zmax = z - 1;

	double			pix_rise = helix_rise/image->sampling()[2];
	long			nun((long) ((zmax - zmin)/pix_rise));
	
	long			i;
	double			v, vr, d, nsum(0), sum1(0), sum2(0), sum1sq(0), sum2sq(0), sum12(0), R(1e37), CC(0);
	double			w, ndsum(0), dsum1(0), dsum2(0), dsum1sq(0), dsum2sq(0), dsum12(0), Rd, CCd;
	
	float*			symdata = new float[datasize];
	for ( i=0; i<datasize; i++ ) symdata[i] = 0;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Helical symmetrization:" << endl;
		cout << "Helical axis origin:            " << image->origin() << endl;
		cout << "Asymmetric unit rise:           " << helix_rise << " A" << endl;
		cout << "Asymmetric unit rotation angle: " << helix_angle*180/M_PI << " degrees" << endl;
		if ( dyad_axis == 2 ) cout << "Dyad axis along the x axis" << endl;
		cout << "Radius:                         " << radius << " pixels" << endl;
		cout << "Limits along helical axis:      " << zmin << " - " << zmax << " pixels" << endl;
		cout << "Units along helical axis:       " << nun << endl;
	}

	if ( norm_flag) w = 1;
	else w = nun;

	symmetry(symmetry_helical_label(helix_rise, helix_angle, dyad_axis, 1, 1));

#ifdef HAVE_GCD
	dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t k){
		symdata[k] = helix_interpolate(k, helix_rise, helix_angle, zmin, zmax, radius, norm_flag);
	});
#else
#pragma omp parallel for
	for ( i=0; i<datasize; i++ )
		symdata[i] = helix_interpolate(i, helix_rise, helix_angle, zmin, zmax, radius, norm_flag);
#endif

	for ( i=0; i<datasize; i++ ) {
		v = symdata[i];
		if ( v ) {
			vr = w * (*this)[i];
			sum1 += vr;
			sum2 += v;
			sum1sq += vr * vr;
			sum2sq += v * v;
			sum12 += vr * v;
			nsum += 1;
		}
		set(i, v);
	}

	if ( dyad_axis == 2 ) {

		if ( !norm_flag ) w = 2;

#ifdef HAVE_GCD
		dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t k){
			symdata[k] = dyad_interpolate(k, norm_flag);
		});
#else
#pragma omp parallel for
		for ( i=0; i<datasize; i++ )
			symdata[i] = dyad_interpolate(i, norm_flag);
#endif

		for ( i=0; i<datasize; i++ ) {
			v = symdata[i];
			if ( v ) {
				vr = w * (*this)[i];
				dsum1 += vr;
				dsum2 += v;
				dsum1sq += vr * vr;
				dsum2sq += v * v;
				dsum12 += vr * v;
				ndsum += 1;
			}
			set(i, v);
		}
	}

	if ( nsum ) {
		d = (sum1sq - sum1*sum1/nsum)*(sum2sq - sum2*sum2/nsum);
		if ( d > 0 ) {
			d = sqrt(d);
			R = sqrt((sum1sq - 2*sum12 + sum2sq)/d);
			CC = (sum12 - sum1*sum2/nsum)/d;
		}
	}
	
	if ( ndsum ) {
		d = (dsum1sq - dsum1*dsum1/ndsum)*(dsum2sq - dsum2*dsum2/ndsum);
		if ( d > 0 ) {
			d = sqrt(d);
			Rd = sqrt((dsum1sq - 2*dsum12 + dsum2sq)/d);
			CCd = (dsum12 - dsum1*dsum2/ndsum)/d;
		}
	}
	
	if ( verbose ) {
		cout << "Helical symmetry R factor:      " << R << endl;
		cout << "Helical symmetry correlation:   " << CC << endl;
		if ( dyad_axis == 2 ) {
			cout << "Dyad R factor:                  " << Rd << endl;
			cout << "Dyad correlation:               " << CCd << endl;
		}
	}
	
	if ( verbose >= VERB_LABEL )
		cout << endl;
	
	delete[] symdata;
	
	statistics();
	
	return R;
}

/**
@brief 	Symmetrizes an image given helical symmetry parameters and seam shift parameter.
@param 	helix_rise		rise per asymmetric unit (angstrom).
@param 	helix_angle		rotation angle per asymmetric unit (radians).
@param 	seam_shift		translation along the seam (subunit height units).
@param 	dyad_axis		dyad axis indicator: 2=dyad axis on x-axis, otherwise none.
@param 	zmin			mimimum z slice to include.
@param 	zmax			maximum z slice to include.
@param 	radius			radius to do symmetrizing over (pixels).
@param 	norm_flag		if 1, normalize
@return Bplot*			plot.

	The data between the z limits are replicated along the helical axis
	according to the helical rise and rotation and the seam to fill the new volume.

**/
Bplot* 		Bimage::seamed_helix_symmetrize(double helix_rise, double helix_angle, 
				double seam_shift, int dyad_axis, int zmin, int zmax, double radius, int norm_flag)
{
	change_type(Float);
	
	if ( helix_rise < 0 ) {
		helix_rise = -helix_rise;
		helix_angle = -helix_angle;
	}

	if ( zmin < 0 ) zmin = 0;
	if ( zmax >= z ) zmax = z - 1;
	if ( zmax < zmin ) swap(zmin, zmax);
		
	long			nproto((long) (TWOPI/fabs(helix_angle) + 0.5));
	double			seam_rise(nproto*helix_rise/seam_shift);
	double			seam_angle(angle_set_negPI_to_PI(nproto*helix_angle)/seam_shift);
	double			pix_helix_rise(helix_rise/image->sampling()[2]);
	double			pix_seam_rise(seam_rise/image->sampling()[2]);
	double			hr(helix_angle*pix_seam_rise - seam_angle*pix_helix_rise);
	double			daz(seam_angle/pix_seam_rise);
	
	if ( verbose & VERB_LABEL ) {
		cout << "Helical symmetrization with a seam with a shift of " << seam_shift << ":" << endl;
		cout << "Asymmetric unit rise:           " << helix_rise << " A" << endl;
		cout << "Asymmetric unit rotation angle: " << helix_angle*180/M_PI << " degrees" << endl;
		cout << "Number of protofilaments:       " << nproto << endl;
		cout << "Seam rise:                      " << seam_rise << " A" << endl;
		cout << "Seam rotation angle:            " << seam_angle*180/M_PI << " degrees" << endl;
		cout << "Origin:                         " << image->origin() << endl;
		cout << "Limits along helical axis:      " << zmin << " - " << zmax << " pixels" << endl;
		cout << "Units along helical axis:       " << (long) ((zmax - zmin)*nproto/pix_seam_rise) << endl << endl;
	}
	
	long   			i, xx, yy, zz, nn;
	long			u, h, k, kmin, kmax;
	double			xs, ys, zs, dx, dy, dz, d, d2;
	double			rad(0), radsq, a, sa, da, as, cosang, sinang;
	double			v, vr, CC, R(0);
	
	float*			symdata = new float[datasize];
	float*			num = new float[datasize];

	long			ia, na(360);
	double*			nsum = new double[na];
	double*			s1 = new double[na];
	double*			s2 = new double[na];
	double*			s1sq = new double[na];
	double*			s2sq = new double[na];
	double*			s12 = new double[na];
	
	for ( i=0; i<datasize; i++ ) symdata[i] = num[i] = 0;
	
	for ( u=0; u<na; u++ ) nsum[u] = s1[u] = s2[u] = s1sq[u] = s2sq[u] = s12[u] = 0;
	
	for ( i=nn=0; nn<n; nn++ ) {
		if ( radius > 0 ) rad = radius;
		if ( rad > image[nn].origin()[0] ) rad = image[nn].origin()[0];
		if ( rad > image[nn].origin()[1] ) rad = image[nn].origin()[1];
		if ( rad > x - image[nn].origin()[0] + 1 )
			rad = x - image[nn].origin()[0] + 1;
		if ( rad > y - image[nn].origin()[1] + 1 )
			rad = y - image[nn].origin()[1] + 1;
		if ( rad < 1 ) rad = x/2;
		radsq = rad*rad;
		for ( zz=0; zz<z; zz++ ) {
			dz = (double)zz - image[nn].origin()[2];
			sa = daz*dz;
			kmin = (long)((zmin-(double)zz)/pix_seam_rise);
			if ( zz < zmin ) kmin += 1;
			kmax = (long)((zmax-(double)zz)/pix_seam_rise);
			if ( zz > zmax ) kmax -= 1;
			for ( yy=0; yy<y; yy++ ) {
				dy = (double)yy - image[nn].origin()[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					dx = (double)xx - image[nn].origin()[0];
					d2 = dx*dx + dy*dy;
					if ( d2 <= radsq ) {
						a = atan2(dy, dx);
						if ( a*pix_seam_rise < seam_angle*dz ) a += TWOPI;
						da = a - sa;
						if ( da < 0 ) da += TWOPI;
						if ( da >= TWOPI ) da -= TWOPI;
						ia = (long) (na*da/TWOPI);
						if ( ia < 0 ) ia += na;
						if ( ia >= na ) ia -= na;
						u = (long) ((a*pix_seam_rise - seam_angle*dz)/hr);
						while ( u < 0 ) u += nproto;
						while ( u >= nproto ) u -= nproto;
						vr = (*this)[i];
						for ( k=kmin; k<=kmax; k++ ) {
							for ( h=-u; h<nproto-u; h++ ) {
								zs = zz + h*pix_helix_rise + k*pix_seam_rise;
								as = h*helix_angle + k*seam_angle;
								cosang = cos(as);
								sinang = sin(as);
								xs = dx*cosang - dy*sinang + image[nn].origin()[0];
								ys = dx*sinang + dy*cosang + image[nn].origin()[1];
								if ( within_boundaries((long)xs, (long)ys, (long)zs) ) {
									v = interpolate(xs, ys, zs, nn, 0.0);
									s1[ia] += vr;
									s2[ia] += v;
									s1sq[ia] += vr * vr;
									s2sq[ia] += v * v;
									s12[ia] += vr * v;
									nsum[ia] += 1;
									symdata[i] += v;
									num[i] += 1;
								}
							}
						}
					}
				}
			}
			if ( verbose & VERB_PROCESS )
				cout << "Completed: " << (nn*z+zz)*100.0/(n*z) << " %" << "\r" << flush;
		}
	}
	
	if ( norm_flag )
		for ( i=0; i<datasize; i++ )
			if ( num[i] ) symdata[i] /= num[i];

	data_assign((unsigned char *) symdata);
	
	statistics();
	
	symmetry(symmetry_helical_label(helix_rise, helix_angle, dyad_axis, 1, seam_shift));

	for ( ia=0; ia<na; ia++ ) {
		R = 1;
		CC = 0;
		if ( nsum[ia] ) {
			d = (s1sq[ia] - s1[ia]*s1[ia]/nsum[ia])*(s2sq[ia] - s2[ia]*s2[ia]/nsum[ia]);
			if ( d > 0 ) {
				R = sqrt(nsum[ia]*(s1sq[ia] - 2*s12[ia] + s2sq[ia])/d);
				CC = (s12[ia] - s1[ia]*s2[ia]/nsum[ia])/sqrt(d);
				s12[ia] = CC;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Helical symmetry R factor:      " << R << endl;
		cout << "Helical symmetry correlation:   " << CC << endl << endl;
	}
	
	
	long			umin, umax, ulen(na/(2*nproto) + 1);
	double			cc_min, cc_max, cc_cut, angle_best(0);
	
	Bstring			title("Seamed helix symmetrization");
	long			ncol(2);
	Bplot*			plot = new Bplot(1, na, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Angle");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).label("CC");
	plot->page(0).column(1).axis(3);
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(360);
	plot->page(0).axis(1).inc(30);
	
//	cout << na << endl;
	
	if ( verbose )
		cout << endl << "Angle\tCC" << endl;
	for ( ia=0, CC=0, cc_min=1, cc_max=0; ia<na; ia++ ) {
		CC += s12[ia];
		umin = ia - ulen;
		umax = ia + ulen;
		for ( u = umin, s1[ia] = 0; u < umax; u++ ) { //cout << u << tab << s1[ia] << endl;
			i = (u < 0)? u + na: u;
			if ( i >= na ) i -= na;
			s1[ia] += s12[i];
		}
		s1[ia] /= 2*ulen;
		if ( cc_min > s1[ia] ) cc_min = s1[ia];
		if ( cc_max < s1[ia] ) cc_max = s1[ia];
//		cout << ia << tab << s1[ia] << endl;
		(*plot)[ia] = ia*360.0/na;
		(*plot)[ia+na] = s1[ia];
		if ( verbose )
			cout << ia*360.0/na << tab << s1[ia] << endl;
	}
	
	CC /= na;
	cc_cut = (cc_max + cc_min)/2;
	
	for ( i=1; i<na/2; i++ ) {
		ia = na - i;
		a = i*TWOPI*1.0L/na;
		if ( s1[i] > cc_cut && s1[i-1] < cc_cut ) angle_best = a;
		if ( s1[i] < cc_cut && s1[i-1] > cc_cut ) angle_best = a;
		a = ia*TWOPI*1.0L/na;
		if ( s1[ia] > cc_cut && s1[ia-1] < cc_cut ) angle_best = a;
		if ( s1[ia] < cc_cut && s1[ia-1] > cc_cut ) angle_best = a;
	}
	
	if ( verbose ) {
		cout << endl;
		cout << "Correlation coefficient:        " << CC << endl;
		cout << "Best seam location angle:       " << angle_best*180.0/M_PI << " degrees" << endl;
		cout << "Min & max correlation:          " << cc_min << " " << cc_max << endl << endl;
	}
	
	Bstring		txt;
	
	txt = Bstring(helix_rise, "Helix rise:                   %g A");
	plot->page(0).add_text(txt);
	txt = Bstring(helix_angle*180.0/M_PI, "Helix rotation angle:    %g degrees");
	plot->page(0).add_text(txt);
	txt = Bstring(seam_rise, "Seam rise:                  %g A");
	plot->page(0).add_text(txt);
	txt = Bstring(seam_angle*180.0/M_PI, "Seam rotation angle:    %g degrees");
	plot->page(0).add_text(txt);
	txt = Bstring(angle_best*180.0/M_PI, "Seam location angle:    %g degrees");
	plot->page(0).add_text(txt);
	txt = Bstring(CC, "Correlation coefficient: %g");
	plot->page(0).add_text(txt);
	
	delete[] num;
	delete[] nsum;
	delete[] s1;
	delete[] s2;
	delete[] s1sq;
	delete[] s2sq;
	delete[] s12;
	
	return plot;
}

/**
@brief 	Calculates a cylindrically symmetrized map.
@param	flag		0: 
@return Bimage*		2D image with the cylindrical average.

	A 2D cylindrically symmetrized average is calculated.

**/
Bimage*		Bimage::symmetrize_cylinder(int flag)
{
	long			i, j, k, xx, yy, zz, nn, cc, ir;
	long			rmax = (long) sqrt(x*x + y*y);
	long			cylsize(n*z*rmax*c);
	double			dx, dy, r, f;
	
	if ( verbose & VERB_PROCESS )
		cout << "Symmetrizing cylindrically around " << image->origin() << endl << endl;
	
	Bimage*			pcyl = NULL;
	
	switch ( flag ) {
		case 1: pcyl = new Bimage(Float, compoundtype, z, rmax, 1, n); break;
		default: pcyl = new Bimage(Float, compoundtype, rmax, 1, z, n);
	}
	
	float*			num = new float[cylsize];
	for ( i=0; i<cylsize; i++ ) num[i] = 0;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				dy = yy - image[nn].origin()[0];
				for ( xx=0; xx<x; xx++ ) {
					dx = xx - image[nn].origin()[0];
					r = sqrt(dx*dx + dy*dy);
					ir = (long) r;
					f = r - ir;
					if ( flag < 1 ) {
						j = ((nn*z + zz)*rmax + ir)*c;
						k = j + c;
					} else {
						j = ((nn*rmax + ir)*z + zz)*c;
						k = j + z*c;
					}
					for ( cc=0; cc<c; cc++, i++, j++, k++ ) {
						pcyl->add(j, (1-f) * (*this)[i]);
						pcyl->add(k, f * (*this)[i]);
						num[j] += 1 - f;
						num[k] += f;
					}
				}
			}
		}
	}
	
	for ( i=0; i<cylsize; i++ ) if ( num[i] ) pcyl->set(i, (*pcyl)[i]/num[i]);
	
	delete[] num;
	
	return pcyl;
}

/**
@brief 	Calculates a cylindrically symmetrized map.
@return int			0.

	The image is replaced with its cylindrically symmetrized version.

**/
int			Bimage::symmetrize_cylinder()
{
	Bimage*			pcyl = symmetrize_cylinder(0);
	
	long			rmax(pcyl->x);
	long			i, j, k, xx, yy, zz, nn, cc, ir;
	double			dx, dy, r, f;
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				dy = yy - image[nn].origin()[0];
				for ( xx=0; xx<x; xx++ ) {
					dx = xx - image[nn].origin()[0];
					r = sqrt(dx*dx + dy*dy);
					ir = (long) r;
					f = r - ir;
					j = ((nn*z + zz)*rmax + ir)*c;
					k = j + c;
					for ( cc=0; cc<c; cc++, i++, j++, k++ )
						set(i, (1-f)*(*pcyl)[j] + f*(*pcyl)[k]);
				}
			}
		}
	}
	
	delete pcyl;
	
	statistics();
	
	return 0;
}

/**
@brief 	Symmetrizes an image given tubular lattice parameters.
@param 	h			units along u vector.
@param 	k			units along v vector.
@param 	latconst	lattice constant (angstrom).
@param 	zmin		mimimum z slice to include.
@param 	zmax		maximum z slice to include.
@param 	radius		radius to do symmetrizing over (pixels).
@param 	norm_flag	if 1, normalize
@return double	 	R factor.

	The data between the z limits are replicated along the helical axis
	according to the lattice parameters to fill the new volume.

**/
double 		Bimage::tube_symmetrize(int h, int k, double latconst,
				int zmin, int zmax, double radius, int norm_flag)
{
	change_type(Float);
	
	long			i;
	double			v, vr, d, nsum(0), sum1(0), sum2(0), sum1sq(0), sum2sq(0), sum12(0), R(1e37), CC(0);
	
	if ( zmin < 0 ) zmin = 0;
	if ( zmax >= z ) zmax = z - 1;
	
	if ( radius <= 0 ) radius = x/2;

	float*			symdata = new float[datasize];
	for ( i=0; i<datasize; i++ ) symdata[i] = 0;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Tubular symmetrization:" << endl;
		cout << "Tubular axis origin:            " << image->origin() << endl;
		cout << "Lattice parameters:             " << h << " x " << k << endl;
		cout << "Lattice constant:               " << latconst << " A" << endl;
		cout << "Radius:                         " << radius << " pixels" << endl;
		cout << "Limits along helical axis:      " << zmin << " - " << zmax << " pixels" << endl;
	}

	
#ifdef HAVE_GCD
	dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t k){
		symdata[k] = tube_interpolate(k, h, k, latconst, zmin, zmax, radius, norm_flag);
		cout << k*100.0/datasize << " %            \r" << flush;
	});
#else
#pragma omp parallel for
	for ( i=0; i<datasize; i++ )
		symdata[i] = tube_interpolate(i, h, k, latconst, zmin, zmax, radius, norm_flag);
#endif

	for ( i=0; i<datasize; i++ ) {
		v = symdata[i];
		if ( v ) {
			vr = (*this)[i];
			sum1 += vr;
			sum2 += v;
			sum1sq += vr * vr;
			sum2sq += v * v;
			sum12 += vr * v;
			nsum += 1;
		}
		set(i, v);
	}
	
	if ( nsum ) {
		d = (sum1sq - sum1*sum1/nsum)*(sum2sq - sum2*sum2/nsum);
		if ( d > 0 ) {
			d = sqrt(d);
			R = sqrt((sum1sq - 2*sum12 + sum2sq)/d);
			CC = (sum12 - sum1*sum2/nsum)/d;
		}
	}
	
	if ( verbose ) {
		cout << "Tubular symmetry R factor:      " << R << endl;
		cout << "Tubular symmetry correlation:   " << CC << endl;
	}
	
	if ( verbose >= VERB_LABEL )
		cout << endl;
	
	statistics();
	
	return R;
}

/**
@brief 	Converts a map to a helix.
@param 	helix_rise		rise per asymmetric unit
@param 	helix_angle		rotation angle per asymmetric unit.
@param 	offset			offset in the xy plane.
@return double			0.

	The data is offset from the central axis based on helical parameters along the axis.
	The data is rotated around the offset vector based on the helix rise 
	and helix angle by an angle:
		              distance * (1 - cos(helix_angle))
		angle2 = atan ---------------------------------
		                        helix_rise

**/
double 		Bimage::convert_to_helix(double helix_rise, double helix_angle, 
				Vector3<double> offset)
{
	change_type(Float);
	
	long			i, xx, yy, zz, nn;
	double			angle, dx, dy;
	double			distance(offset.length());
	double			angle0(atan2(offset[1], offset[0]));
	double			angle2(atan2(distance*sqrt(1-cos(helix_angle)), helix_rise*1.0));
	Vector3<double>	ori, old, rold, axis;
	Vector3<double>	u;
	Matrix3			mat;
	
	float*			nudata = new float[datasize];
	
	for ( i=0; i<datasize; i++ ) nudata[i] = 0;

	for ( i=nn=0; nn<n; nn++ ) {
		ori = image[nn].origin();
		u = image[nn].sampling();
		for ( zz=0; zz<z; zz++ ) {
			ori[2] = zz;
			angle = helix_angle*(zz - image[nn].origin()[2])/u[2] + angle0;
			dx = distance*cos(angle)/u[0];
			dy = distance*sin(angle)/u[1];
			axis[0] = dx;
			axis[1] = dy;
			axis.normalize();
			mat = Matrix3(axis, angle2);
			for ( yy=0; yy<y; yy++ ) {
				old[1] = yy - dy - ori[1];
				for ( xx=0; xx<x; xx++, i++ ) {
					old[0] = xx - dx - ori[0];
					rold = (mat * old) + ori;
					nudata[i] = interpolate(0, rold, nn, 0.0);
				}
			}
		}
	}
	
	data_assign((unsigned char *) nudata);
	
	statistics();
	
	return 0;
}

/**
@brief 	Transforms a tubular cylinder to an elliptical profile.
@param 	angle		angle of elliptical maximum.
@param 	shift		outward shift at elliptical maximum radius.
@return int			0.

	The density is shifted towards or away from the origin by a specific
	amount and in a direction determined by an elliptical profile.

**/
int			Bimage::distort_elliptically(double angle, double shift)
{
	change_type(Float);
	
	if ( verbose ) {
		cout << "Elliptical distortion: " << endl;
		cout << "Major axis angle:              " << angle*180.0/M_PI << " degrees" << endl;
		cout << "Shift at major axis:           " << shift << " pixels" << endl << endl;
	}
	
	long			xx, yy, zz, nn, i, j, inz, ix, iy;
	double			dx, dy, a, r, s, rx, ry, fx, fy;
	float*			nudata = new float[datasize];
	
	for ( i=0; i<datasize; i++ ) nudata[i] = 0;
	
	for ( j=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			inz = (nn*z + zz)*y*x;
			for ( yy=0; yy<y; yy++ ) {
				dy = yy - image[nn].origin()[0];
				for ( xx=0; xx<x; xx++, j++ ) {
					dx = xx - image[nn].origin()[0];
					r = sqrt(dx*dx + dy*dy);
					a = 2*(atan2(dy, dx) - angle);
					s = shift*cos(a)/r;
//					cout << x << tab << y << tab << z << tab << n << tab << s << endl;
					if ( s < 1 ) {
						rx = dx*(1 - s) + image[nn].origin()[0];
						ry = dy*(1 - s) + image[nn].origin()[1];
						ix = (long) rx;
						iy = (long) ry;
						if ( ix < x-1 && iy < y-1 ) {
							fx = rx - ix;
							fy = ry - iy;
							i = inz + iy*x + ix;
							nudata[j] = (1-fx)*(1-fy)*(*this)[i] + fx*(1-fy)*(*this)[i+1] +
								(1-fx)*fy*(*this)[i+x] + fx*fy*(*this)[i+x+1];
						}
					}
				}
			}
		}
	}
//	cout << "***" << endl;
	data_assign((unsigned char *) nudata);
	
	statistics();
	
	return 0;
}

Vector3<double>	Bimage::test_helix_parameters(double angle, double hires, double lores,
				Vector3<long> mask_size, Vector3<long> mask_start, long max_iter,
				fft_plan planf, fft_plan planb, double& cc)
{
	long			j;
	double			dd(1);
	Vector3<double>	shift(0,0,1), scale(1,1,1), translate, axis(0,0,1);
	Vector3<double>	origin(image->origin());
	Matrix3			mat(axis, angle);
	Bimage*			pt;
	
	for ( j=0; j<max_iter && fabs(shift[2]) > 0.001; ++j ) {
		pt = transform(size(), scale, origin, translate, mat, FILL_BACKGROUND, 0);
//		write_img("t.pif", pt);
		pt->edge(2, mask_size, mask_start, 2, FILL_BACKGROUND, 0);
		shift = find_shift(pt, NULL, hires, lores, z/4.0, 0, 1, planf, planb, cc);
		delete pt;
		if ( fabs(shift[2]) > 0.001 ) {
			translate[2] += dd*shift[2];
//			cout << angle*180.0/M_PI << tab << shift << tab << translate[2] << tab << cc << endl;
		}
		dd *= 0.9;
	}

	return translate;
}

/**
@brief 	Finds the best helix parameters for helical map.
@param 	angle_start		start value for angle.
@param 	angle_end		end value for angle.
@param 	angle_step		step size for angle.
@param 	bin				bin image for faster searching (limited to 1,2,3).
@param 	hires			high resolution limit in angstroms.
@param 	lores			low resolution limit in angstroms.
@param 	radius			radius for mask (voxels).
@return Bplot*			plot with search results.

	An incremental search is done for the rotation angle, with the rise
	inferred from the cross-correlation shift.

**/
Bplot*		Bimage::find_helix_parameters(double angle_start, double angle_end,
				double angle_step, int bin, double hires, double lores, double radius)
{
	angle_step = fabs(angle_step);
	if ( angle_step < 0.0001 ) angle_step = 0.0001;
	if ( angle_start > angle_end ) swap(angle_start, angle_end);
	 
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	
	if ( lores > 0 && hires > lores ) swap(hires, lores);
	check_resolution(hires);
	if ( lores < hires ) lores = x*image->sampling()[0];
	
	if ( radius < 1 ) radius = x/2;
	if ( radius > image->origin()[0] ) radius = image->origin()[0];
	if ( radius > x - image->origin()[0] ) radius = x - image->origin()[0];
	if ( radius > image->origin()[1] ) radius = image->origin()[1];
	if ( radius > y - image->origin()[1] ) radius = y - image->origin()[1];

	long			i, na, max_iter(10);
	double			rise, angle, offset, bestrise(0), bestangle(0);
	double			bestcc(0), mincc(1);

//	Bimage*			pt;
	
	if ( verbose ) {
		cout << "Searching for helical parameters by cross-correlation:" << endl;
		cout << "Angular range                   " << angle_start*180.0/M_PI << " - " 
			<< angle_end*180.0/M_PI << " degrees" << endl;
		cout << "Angular increment:              " << angle_step*180.0/M_PI << " degrees" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Binning:                        " << bin << endl << endl;
	}

	radius /= bin;
	angle_end += angle_step/2;	// To ensure the final angle is sampled
	
	Bimage*			pb = copy();
	pb->integer_interpolation(bin, 42);
	pb->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::find_helix_parameters: binned size = " << pb->size() << endl;

	Vector3<double>	bestorigin, scale(1,1,1), axis(0,0,1);
	Vector3<long>	mask_size((long) (2*radius), (long) (2*radius), (pb->sizeZ()*3)/4);
	Vector3<double>	mask_start(pb->image->origin()[0] - radius, pb->image->origin()[1] - radius, pb->sizeZ()/8);
	Matrix3			mat;
	Vector3<double>	shift;

	for ( na=0, angle = angle_start; angle <= angle_end; angle += angle_step ) na++;

	vector<double>		a(na);
	double*				cc = new double[na];
	Vector3<double>*	translate = new Vector3<double>[na];

	for ( i=0, angle = angle_start; angle <= angle_end; angle += angle_step, ++i )
		a[i] = angle;

	fft_plan		planf = pb->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = pb->fft_setup(FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(na, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bimage*		pt = pb->copy();
		translate[i] = pt->test_helix_parameters(a[i], hires, lores, mask_size, mask_start, max_iter, planf, planb, cc[i]);
		delete pt;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<na; i++ ) {
		Bimage*		pt = pb->copy();
		translate[i] = pt->test_helix_parameters(a[i], hires, lores, mask_size, mask_start, max_iter, planf, planb, cc[i]);
		delete pt;
	}
#endif

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	Bstring			title("Helical parameter search");
	int				ncol(3);
	Bplot*			plot = new Bplot(1, na, 3);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(3);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Rise (A)");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(4);
	plot->page(0).column(1).label("Rotation angle (degrees)");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).type(4);
	plot->page(0).column(2).label("CC");
	plot->page(0).column(2).axis(5);
	plot->page(0).axis(3).min(angle_start*180.0/M_PI);
	plot->page(0).axis(3).max(angle_end*180.0/M_PI);
	
	if ( verbose )
		cout << "Rise\tAngle\tOffset\tCC" << endl;
	for ( i=0; i<na; ++i ) {
		angle = a[i];
		rise = pb->sampling(0)[2]*translate[i][2];
		offset = pb->sampling(0)[2]*sqrt(translate[i][0]*translate[i][0] + translate[i][1]*translate[i][1]);
		(*plot)[i] = rise;
		(*plot)[i+na] = angle*180.0/M_PI;
		(*plot)[i+2*na] = cc[i];
		if ( bestcc < cc[i] ) {
			bestcc = cc[i];
			bestangle = angle;
			bestrise = rise;
		}
		if ( mincc > cc[i] ) mincc = cc[i];
		if ( verbose )
			cout << setprecision(6) << rise << tab << angle*180.0/M_PI
				<< tab << offset << tab << cc[i] << endl << flush;
	}
	
	delete pb;
	delete[] cc;
	delete[] translate;
	
	if ( bestrise < 0 ) {
		bestrise = -bestrise;
		bestangle = -bestangle;
	}
	
	if ( verbose ) {
		cout << "Best helical rise:              " << bestrise << " A" << endl;
		cout << "Best helical rotation angle:    " << bestangle*180.0/M_PI << " degrees" << endl;
		cout << "Best correlation coefficient:   " << bestcc << endl << endl;
	}

	Bstring		txt;
	
	txt = Bstring(bestrise, "Best helical rise:              %g A");
	plot->page(0).add_text(txt);
	txt = Bstring(bestangle*180.0/M_PI, "Best helical rotation angle:    %g degrees");
	plot->page(0).add_text(txt);
	txt = Bstring(bestcc,   "Best correlation coefficient:   %g");
	plot->page(0).add_text(txt);
	txt = Bstring(mincc,    "Minimum correlation coefficient: %g");
	plot->page(0).add_text(txt);
	
	return plot;
}

double		Bimage::test_helix_parameters(double rise, double angle,
				Vector3<long> mask_size, Vector3<long> mask_start)
{
	double			cc;
	Vector3<double>	scale(1,1,1), axis(0,0,1);
	Matrix3			mat = Matrix3(axis, angle);
	Vector3<double>	translate(0,0,rise/image->sampling()[2]);
	
	Bimage*			pt = transform(size(), scale, image->origin(), translate, mat, FILL_BACKGROUND, 0);

	pt->edge(2, mask_size, mask_start, 2, FILL_USER, 0);
	
	cc = correlate(pt);
	
	delete pt;
	
	return cc;
}

/**
@brief 	Finds the best helix parameters for helical map.
@param 	rise_start		start value for rise.
@param 	rise_end		end value for rise.
@param 	rise_step		step size for rise.
@param 	angle_start		start value for angle.
@param 	angle_end		end value for angle.
@param 	angle_step		step size for angle.
@param 	bin				bin image for faster searching (limited to 1,2,3).
@param 	radius			radius for mask (voxels).
@return Bplot*			plot with search results.

	An incremental search is done for the rise and rotation angle.

**/
Bplot*		Bimage::find_helix_parameters(double rise_start, double rise_end,
				double rise_step, double angle_start, double angle_end, double angle_step,
				int bin, double radius)
{
	rise_step = fabs(rise_step);
	if ( rise_step < 0.001 ) rise_step = 0.001;
	if ( rise_start > rise_end ) swap(rise_start, rise_end);
	
	angle_step = fabs(angle_step);
	if ( angle_step < 0.0001 ) angle_step = 0.0001;
	if ( angle_start > angle_end ) swap(angle_start, angle_end);
	 
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	
	if ( radius < 1 ) radius = x/2;
	if ( radius > image->origin()[0] ) radius = image->origin()[0];
	if ( radius > x - image->origin()[0] ) radius = x - image->origin()[0];
	if ( radius > image->origin()[1] ) radius = image->origin()[1];
	if ( radius > y - image->origin()[1] ) radius = y - image->origin()[1];

	long			i, na;
	double			rise, angle, bestrise(0), bestangle(0);
	double			bestcc(0), mincc(1);

	if ( verbose ) {
		cout << "Searching for helical parameters:" << endl;
		cout << "Rise range                      " << rise_start << " - "
			<< rise_end << " A" << endl;
		cout << "Rise increment:                 " << rise_step << " A" << endl;
		cout << "Angular range                   " << angle_start*180.0/M_PI << " - " 
			<< angle_end*180.0/M_PI << " degrees" << endl;
		cout << "Angular increment:              " << angle_step*180.0/M_PI << " degrees" << endl;
		cout << "Binning:                        " << bin << endl;
		cout << "Radius:                         " << radius << endl << endl;
	}

	radius /= bin;
	rise_end += rise_step/2;
	angle_end += angle_step/2;	// To ensure the final angle is sampled

	Bimage*			pb = copy();
	pb->integer_interpolation(bin, 42);
	pb->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::find_helix_parameters: binned size = " << pb->size() << endl;

	Vector3<double>	bestorigin, scale(1,1,1), translate, axis(0,0,1);
	Vector3<long>	mask_size((int) (2*radius), (int) (2*radius), (pb->sizeZ()*3)/4);
	Vector3<double>	mask_start(pb->image->origin()[0] - radius, pb->image->origin()[1] - radius, pb->sizeZ()/8);
	Matrix3			mat;

	for ( na=0, rise = rise_start; rise < rise_end; rise += rise_step )
		for ( angle = angle_start; angle <= angle_end; angle += angle_step ) na++;
	
	double*			cc = new double[na];
	vector<double>	r(na);
	vector<double>	a(na);
	
	for ( i=0, rise = rise_start; rise < rise_end; rise += rise_step )
		for ( angle = angle_start; angle <= angle_end; angle += angle_step, i++ ) {
			r[i] = rise;
			a[i] = angle;
		}

#ifdef HAVE_GCD
	dispatch_apply(na, dispatch_get_global_queue(0, 0), ^(size_t i){
		cc[i] = pb->test_helix_parameters(r[i], a[i], mask_size, mask_start);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<na; i++ )
		cc[i] = pb->test_helix_parameters(r[i], a[i], mask_size, mask_start);
#endif

	Bstring			title("Helical parameter search");
	int				ncol(3);
	Bplot*			plot = new Bplot(1, na, 3);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(3);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Rise (A)");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(4);
	plot->page(0).column(1).label("Rotation angle (degrees)");
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).type(4);
	plot->page(0).column(2).label("CC");
	plot->page(0).column(2).axis(5);
	plot->page(0).axis(3).min(angle_start*180.0/M_PI);
	plot->page(0).axis(3).max(angle_end*180.0/M_PI);
	
	if ( verbose )
		cout << "Rise\tAngle\tCC" << endl;
	for ( i=0; i<na; ++i ) {
		rise = r[i];
		angle = a[i];
		(*plot)[i] = rise;
		(*plot)[i+na] = angle*180.0/M_PI;
		(*plot)[i+2*na] = cc[i];
		if ( bestcc < cc[i] ) {
			bestcc = cc[i];
			bestangle = angle;
			bestrise = rise;
		}
		if ( mincc > cc[i] ) mincc = cc[i];
		if ( verbose )
			cout << setprecision(6) << rise << tab << angle*180.0/M_PI
				<< tab << cc[i] << endl << flush;
	}
	
	delete pb;
	delete[] cc;
	
	if ( bestrise < 0 ) {
		bestrise = -bestrise;
		bestangle = -bestangle;
	}
	
	if ( verbose ) {
		cout << "Best helical rise:              " << bestrise << " A" << endl;
		cout << "Best helical rotation angle:    " << bestangle*180.0/M_PI << " degrees" << endl;
		cout << "Best correlation coefficient:   " << bestcc << endl << endl;
	}

	Bstring		txt;
	
	txt = Bstring(bestrise, "Best helical rise:              %g A");
	plot->page(0).add_text(txt);
	txt = Bstring(bestangle*180.0/M_PI, "Best helical rotation angle:    %g degrees");
	plot->page(0).add_text(txt);
	txt = Bstring(bestcc,   "Best correlation coefficient:   %g");
	plot->page(0).add_text(txt);
	txt = Bstring(mincc,    "Minimum correlation coefficient: %g");
	plot->page(0).add_text(txt);
	
	return plot;
}


int			Bimage::helix_segment_correlation_one(long i,
				double angle_start, double angle_end, double angle_step,
				int bin, double hires, double lores, double radius,
				fft_plan planf, fft_plan planb, double* cc)
{
	long			k;
	long			m = (long) ((angle_end - angle_start)/angle_step) + 1;

	double			angle;
	Vector3<double>	axis(0,0,1);
	Bimage*			p1, *p2;

	for ( k = i*m, angle = angle_start; angle <= angle_end; angle += angle_step, k++ ) {
		p1 = extract(i);
		p2 = extract(i+1);
		p1->rotate(axis, angle);
		p1->find_shift(p2, NULL, hires, lores, p1->sizeX()/4, 0, 1, planf, planb, cc[k]);
		delete p2;
		delete p1;
	}
	
	return 0;
}

/**
@brief 	Calculates the correlation over rotation angles between helical segments.
@param 	thickness		segment thickness.
@param 	angle_start		start value for angle.
@param 	angle_end		end value for angle.
@param 	angle_step		step size for angle.
@param 	bin				bin image for faster searching (limited to 1,2,3).
@param 	hires			high resolution limit in angstroms.
@param 	lores			low resolution limit in angstroms.
@param 	radius			radius for mask (voxels).
@return Bplot*			plot with coefficients over all rotations results.

	The image is split up into segments in the z direction. Every adjacent pair
	of segments are cross-correlated over an angular range to find the best
	rotation. The cross-correlation can be resolution-limited. The image may
	also be masked beyond a radius and binned to speed up execution.

**/
Bplot*		Bimage::helix_segment_correlation(int thickness,
				double angle_start, double angle_end, double angle_step,
				int bin, double hires, double lores, double radius)
{
	if ( hires > lores ) swap(hires, lores);
	
	if ( radius <= 0 || radius > x/2 ) radius = x/2;
	
	Vector3<long> 	start, tile_size(x, y, thickness), step_size;
	Vector3<long>	mask_size((long) (2*radius), (long) (2*radius), z);
	Vector3<double>	mask_start(image->origin()[0] - radius, image->origin()[1] - radius, 0);

	calculate_background();
	edge(2, mask_size, mask_start, 2, FILL_BACKGROUND, 0);

//	cout << mask_start << tab << mask_size << endl;
//	cout << start << tab << p->size() << endl;
//	cout << tile_size << endl;
	
	Bimage*			ptiles = extract_tiles(0, start, size(), tile_size, step_size, 0);

//	write_img("ptiles.pif", ptiles);
	
	if ( bin > 1 ) {
		ptiles->integer_interpolation(bin, 42);
		ptiles->calculate_background();
	}

	long			i, k;
	long			nt(ptiles->images());
	long			m = (long) ((angle_end - angle_start)/angle_step) + 1;
	double			angle;

//	for ( i=0; i<nt; i++ ) ptiles->image[i].origin(ptiles->size()/2);
	ptiles->origin(ptiles->size()/2);

	if ( verbose )	{
		cout << "Calculating cross-correlations between helical segments:" << endl;
		cout << "Thickness:                      " << tile_size[2] << endl;
		cout << "Angular range                   " << angle_start*180.0/M_PI << " - "
			<< angle_end*180.0/M_PI << " degrees" << endl;
		cout << "Angular increment:              " << angle_step*180.0/M_PI << " degrees" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "Radius:                         " << radius << endl;
		cout << "Binning:                        " << bin << endl << endl;
	}

	double*			cc = new double[m*nt];
	for ( i=0; i<m*nt; i++ ) cc[i] = 0;

	fft_plan		planf = ptiles->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = ptiles->fft_setup(FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(nt-1, dispatch_get_global_queue(0, 0), ^(size_t i){
		ptiles->helix_segment_correlation_one(i,
				angle_start, angle_end, angle_step,
				bin, hires, lores, radius,
				planf, planb, cc);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nt; ++i )
		ptiles->helix_segment_correlation_one(i,
				angle_start, angle_end, angle_step,
				bin, hires, lores, radius,
				planf, planb, cc);
#endif

	Bstring			title("Helix segment correlation");
	int				ncol(nt);
	RGB<float>		color;
	Bplot*			plot = new Bplot(1, m, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Angle");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).label("CC");
		plot->page(0).column(i).axis(3);
		color.spectrum(i,1,ncol-1);
		plot->page(0).column(i).color(color.r(),color.g(),color.b());
	}
	plot->page(0).axis(1).label("Angle");
	plot->page(0).axis(1).min(angle_start*180.0/M_PI);
	plot->page(0).axis(1).max(angle_end*180.0/M_PI);
	plot->page(0).axis(3).label("CC");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	
	if ( verbose ) {
		cout << endl << "Angle";
		for ( i=0; i<nt-1; i++ ) cout << "\tCC(" << i+1 << ")";
		cout << endl;
		for ( k=0, angle = angle_start; angle <= angle_end; angle += angle_step, k++ ) {
			cout << setprecision(4) << angle*180.0/M_PI;
			for ( i=0; i<nt-1; i++ ) cout << tab << cc[i*m + k];
			cout << endl;
		}
	}
	
	for ( k=0, angle = angle_start; angle <= angle_end; angle += angle_step, k++ ) {
		(*plot)[k] = angle*180.0/M_PI;
		for ( i=0; i<nt-1; i++ ) (*plot)[(i+1)*m+k] = cc[i*m + k];
	}
	
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	delete[] cc;
	delete ptiles;
	
	return plot;
}

/*
@brief Transforms rows in a 2D image and shift the phases to the origin.
@return int			0, <0 on error.

	Each row of pixels is transformed and phase shifted.
	Only 2D images are transformed.
**/
int			Bimage::transform_lines()
{
	if ( z > 1 ) {
		cerr << "Error: Line transforms can only be calculated from 2D images!" << endl;
		return -1;
	}
	
	simple_to_complex();
	change_type(Float);

	fft_plan		plan = fft_setup_plan(x, 1, 1, FFTW_FORWARD, 1);

	if ( verbose & VERB_FULL )
		cout << "Transforming lines of size:     " << x << endl;
	
	long			i, xx, yy, nn, hx(x/2);
	double			dx, shift, phi, sx(x);
	Complex<float>*	data = (Complex<float> *) data_pointer();
	
	for ( i=nn=0; nn<n; nn++ ) {
		shift = -image[nn].origin()[0]/x;
		if ( verbose & VERB_FULL )
			cout << "Shift:                          " << shift << endl;
		for ( yy=0; yy<y; yy++ ) {
			data = (Complex<float> *) (data_pointer(c*i));
			fftw(plan, data);
			for ( xx=0; xx<x; xx++, i++ ) {
				dx = (xx < hx)? xx: xx - sx;
				phi = MIN2PI*dx*shift;
				data[xx].shift_phi(phi);
			}
		}
	}

	fft_destroy_plan(plan);

	if ( verbose & VERB_FULL )
		cout << "Finished transforming" << endl << endl;
	
	return 0;
}

/**
@brief 	Calculates a helical cross section from line transforms of a filament.
@param 	helix_rise	helical subunit rise in angstrom.
@param 	helix_angle	helical subunit rotation in radians.
@param 	scale		scale of cross section.
@param 	hires		high resolution limit.
@return Bimage* 	2D cross section image.

	The filament image must be oriented with the helical axis coinciding with 
	the y axis. Each row of pixels is transformed and written into a
	2D transform based on an orientation calculated from given helical
	parameters.

**/
Bimage*		Bimage::helical_cross_section(double helix_rise, double helix_angle,
				double scale, double hires)
{
	if ( z > 1 ) {
		cerr << "Error: Cross sections can only be calculated from 2D images!" << endl;
		return NULL;
	}
	
	if ( x > y ) {
		if ( verbose & VERB_PROCESS )
			cout << "Rotating image" << endl;
		Bstring		order("-yxz");
		reslice(order);
		order = 0;
	}
	
	check_resolution(hires);
	
	Vector3<long>	recsize((int) (x*scale), (int) (x*scale), 1);
	
	long			i, j, xx, yy, nn;
	long			ix, iy, pad_factor(2);
	double			a, cosang, sinang, dx;
	int				ft_size = part_ft_size(recsize[0], 1, pad_factor);
	double			pix_rise = helix_rise/image->sampling()[1];
	double			apz = helix_angle/pix_rise;

	Vector3<long>	pad_size(ft_size, y, z);
	Bimage*			pp = pad_copy(pad_size, FILL_BACKGROUND);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "2D reconstruction:" << endl;
		cout << "Size:                           " << recsize << endl;
		cout << "Helical rise and angle:         " << helix_rise << " " << helix_angle*180/M_PI << endl;
		cout << "Rotation per pixel along axis:  " << apz*180/M_PI << endl;
		cout << "Resolution limit:               " << hires << " Å" << endl;
		cout << "Fourier transform size:         " << ft_size << endl;
	}

	pp->transform_lines();
	
	Bimage* 		pcs = new Bimage(Float, TComplex, recsize, n);
	pcs->sampling(sampling(0)/scale);
	float*			w = new float[n*pcs->image_size()];
//	cout << "*** Channels = " << pcs->channels() << endl;
//	cout << "*** Compound type = " << pcs->compound_type() << endl;

	double			sx = pp->x;
	double			hx = sx/2;
	double			rscale = pcs->x/(scale*sx);
	double			mx = sx * pp->image->sampling()[0] / hires;
	if ( mx > sx ) mx = sx;
	
//	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::helical_cross_section: sampling=" << image->sampling()[0] << " mx=" << mx << endl;

	for ( i=nn=0; nn<pp->n; nn++ ) {
		for ( yy=0; yy<pp->y; yy++ ) {
			// Set view
			a = -(yy - pp->image[nn].origin()[1])*apz;
			cosang = cos(a);
			sinang = sin(a);
			// Calculate from helical parameters and write line into reconstruction
			for ( xx=0; xx<pp->x; xx++, i++ ) {
				dx = ( xx < hx )? xx: xx - sx;
				if ( fabs(dx) <= mx ) {
					dx *= rscale;
					ix = (long) floor(dx*cosang + 0.5);
					iy = (long) floor(dx*sinang + 0.5);
					if ( ix < 0 ) ix += pcs->x;
					if ( iy < 0 ) iy += pcs->y;
					if ( ix >= 0 && iy >= 0 && ix < pcs->x && iy < pcs->y ) {
						j = (nn*pcs->y + iy)*pcs->x + ix;
						pcs->add(j, pp->complex(i));
						w[j] += 1;
					}
				}
			}
		}
	}

	delete pp;
	
	// Weigh reconstructions
	if ( verbose & VERB_PROCESS )
		cout << "Weighing reconstruction" << endl;
	for ( i=0; i<pcs->size().volume() * pcs->n; i++ )
		if ( w[i] > 0 ) pcs->set(i, pcs->complex(i) / w[i]);

	delete[] w;
	
	pcs->phase_shift_to_center();
	
	if ( verbose & VERB_PROCESS )
		cout << "Backtransforming reconstruction" << endl << endl;
	pcs->fft_back();
	
	pcs->statistics();
	
	return pcs;
}

/**
@brief 	Calculates the ratio of variance inside and outside a given radius.
@param 	snradius	radius of signal.
@return double 		the SNR.

	The region outside the radial limit is considered noise and the region
	inside is considered signal plus noise.

**/
double		Bimage::snvariance(double snradius)
{
	if ( snradius < 1 || snradius > x/2 ) snradius = x/4;
	
	long			nn, i, xx, yy, fn, bn;
	double			x2, y2, r2(snradius*snradius), v, fa, fv, ba, bv, snr(0);
	
	if ( verbose )
		cout << "Image\tSNR" << endl;
	for ( i=nn=0; nn<n; nn++ ) {
		fn = bn = 0;
		fa = fv = ba = bv = 0;
		for ( yy=0; yy<y; yy++ ) {
			y2 = yy - image[nn].origin()[1];
			y2 *= y2;
			for ( xx=0; xx<x; xx++, i++ ) {
				x2 = xx - image[nn].origin()[0];
				x2 *= x2;
				v = (*this)[i];
				if ( x2 + y2 <= r2 ) {
					fa += v;
					fv += v*v;
					fn++;
				} else {
					ba += v;
					bv += v*v;
					bn++;
				}
			}
		}
		fa /= fn;
		ba /= bn;
		fv = fv/fn - fa*fa;
		bv = bv/bn - ba*ba;
		if ( fv < 0 ) fv = 0;
		if ( bv < 0 ) bv = 0;
		if ( bv ) {
			snr = fv/bv - 1;
		} else {
			snr = 0;
			cerr << "Error: The background variance is zero or less! (" << bv << ")" << endl;
		}
		if ( verbose )
			cout << nn+1 << tab << snr << endl;
	}
	
	if ( verbose )
		cout << endl;
	
	return snr;
}

/**
@brief 	Extrudes a 2D cross section into a 3D continuous helix.
@param 	length		length in z.
@param 	helix_rise	helical subunit rise in angstrom.
@param 	helix_angle	helical subunit rotation in radians.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value of edge voxels.
@return int			0.

	The helical axis is at the origin and alongthe a axis.

**/
int			Bimage::extrude_cross_section(long length, double helix_rise,
				double helix_angle, int fill_type, double fill)
{
	change_type(Float);

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	
	long			i, xx, yy, zz;
	double			a, ca, sa, dx, dy, cx, cy;
	
	if ( verbose ) {
		cout << "Extruding the image:" << endl;
		cout << "Helical rise and rotation:      " << helix_rise << " A  " << helix_angle*180.0/M_PI << "°" << endl;
		cout << "Length:                         " << length << endl;
	}
	
	long			nusize = (long) (size().volume() * length);
	float*			nudata = new float[nusize];
	for ( i=0; i<nusize; i++ ) nudata[i] = 0;
	
	for ( i=zz=0; zz<length; zz++ ) {
		a = zz*helix_angle*image->sampling()[0]/helix_rise;
//		cout << "z = " << z << "  angle = " << a*180.0/M_PI << endl;
		ca = cos(a);
		sa = sin(a);
		for ( yy=0; yy<y; yy++ ) {
			dy = yy - image->origin()[1];
			for ( xx=0; xx<x; xx++, i++ ) {
				dx = xx - image->origin()[0];
				cx = dx*ca - dy*sa + image->origin()[0];
				cy = dx*sa + dy*ca + image->origin()[1];
				nudata[i] = interpolate(cx, cy, 0, 0, fill);
			}
		}
	}
	
	sizeZ(length);
	data_assign((unsigned char *)nudata);
	
	return 0;
}

/**
@brief 	Estimates the width of a filament.
@param 	width		window size.
@param 	lim_lo		minimum filament width.
@param 	lim_hi		maximum filament width.
@return Bplot*		plot of filament widths.

	The filament axis must be along the y axis.

**/
Bplot*		Bimage::filament_width(long width, long lim_lo, long lim_hi)
{
	if ( lim_lo > lim_hi ) swap(lim_lo, lim_hi);
	if ( lim_lo < 1 ) lim_lo = 1;
	if ( lim_hi < 1 || lim_hi > x ) lim_hi = x;
	
	long			i, xx, yy, zz, nn, left, right;
	long			n1, n2, hw(width/2);
	double			ox(image->origin()[0]), a, b, slope_max, slope_min;
	if ( ox < hw || ox > x - hw ) ox = x/2;
	
	Vector3<double>	u(image->sampling());
	
	double*			vx = new double[x];
	double*			vy = new double[x];
	
	long			leftmin = (long) (ox - lim_hi/2);
	long			leftmax = (long) (ox - lim_lo/2);
	long			rightmin = (long) (ox + lim_lo/2);
	long			rightmax = (long) (ox + lim_hi/2);
	if ( leftmin < 0 ) leftmin = 0;
	if ( leftmax <= leftmin ) leftmax = leftmin + 1;
	if ( rightmax >= x ) rightmax = x - 1;
	if ( rightmin >= rightmax ) rightmin = rightmax - 1;
	
	for ( xx=0; xx<x; xx++ ) vx[xx] = xx;
	
	if ( verbose ) {
		cout << "Calculating filament width:" << endl;
		cout << "Origin:                         " << ox << endl;
		cout << "Moving line window width:       " << width << endl;
		cout << "Allowed filament width range:   " << lim_lo << " - " << lim_hi << endl;
		cout << "                                " <<
			lim_lo*u[1] << " - " << lim_hi*u[1] << " A" << endl << endl;
	}
	
	Bstring			title("Filament width");
	int				ncol(2);
	Bplot*			plot = new Bplot(1, y, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Distance (A)");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).label("Width (A)");
	plot->page(0).column(1).axis(3);
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(360);
//	plot->page(0).axis(1).inc(30);
	
	if ( verbose )
		cout << "y\ty(A)\tleft\tright\twidth\twidth(A)" << endl;
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				left = right = 0;
				for ( xx=0; xx<x; xx++, i++ ) vy[xx] = (*this)[i];
				slope_max = slope_min = 0;
				for ( xx=0; xx<x; xx++ ) {
					n1 = (long)xx - hw;
					if ( n1 < 0 ) n1 = 0;
					n2 = xx + hw;
					if ( n2 >= x ) n2 = x - 1;
					linear_least_squares(n1, n2, vx, vy, a, b);
//					cout << " " << b;
					if ( xx >= leftmin && xx <= leftmax && slope_max < b ) {
						slope_max = b;
						left = xx;
					}
					if ( xx >= rightmin && xx <= rightmax && slope_min > b ) {
						slope_min = b;
						right = xx;
					}
				}
				(*plot)[yy] = yy*u[1];
				(*plot)[yy+y] = (right - left)*u[0];
//				cout << endl;
				if ( verbose )
					cout << yy << tab << yy*u[1] << tab << left << tab
						<< right << tab << right - left << tab
						<< (right - left)*u[0] << endl;
			}
		}
	}
	
	delete[] vx;
	delete[] vy;
	
	return plot;
}

/**
@brief 	Estimates the density per pixel length of a filament.
@param 	width		filament width (pixels).
@return Bimage*		one dimensional image with line integrals.

	The filament axis must be along the long axis (x or y).
	The filament width must be about half of the image width.
	The background is calculated for each line from the regions 
	outside the width of the filament and subtracted from all
	values in the line.

**/
Bimage*		Bimage::filament_density(double width)
{
	if ( y < x ) {
		Bstring		order("y-xz");
		reslice(order);
		order = 0;
	}
	
	if ( width < 1 || width >= x ) width = x/2;
	
	if ( verbose & VERB_FULL ) {
		cout << "Calculating filament density per pixel length:" << endl;
		cout << "Filament width:                 " << width <<
			" (" << width*image->sampling()[0] << " A)" << endl;
		cout << "Image width:                    " << x << endl << endl;
	}

	long			i, xx, yy, zz, nn;
	double			t, tl, tr, hw((x-width)/2), hwm(x-hw);
	
	Bimage*			pd = new Bimage(Float, TSimple, y, 1, 1, n);
	pd->sampling(sampling(0));
	
	for ( nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				t = tl = tr = 0;
				i = index(0, 0, yy, zz, nn);
				for ( xx=0; xx<hw; xx++, i++ )
						tl += (*this)[i];
				for ( ; xx<hwm; xx++, i++ )
						t += (*this)[i];
				for ( ; xx<x; xx++, i++ )
						tr += (*this)[i];
				tl /= hw;
				tr /= hw;
				t /= width;
				if ( tl < tr ) t -= tl;
				else t -= tr;
				i = nn*pd->x + yy;
				pd->set(i, t);
//				cout << i << tab << t << endl;
			}
		}
	}
	
	pd->statistics();
	
	return pd;
}

/**
@brief 	Reconstructs a filament from a set of projections.
@param 	hi_res		high resolution limit.
@param 	flag		0=sequential, 1=random.
@return Bimage*		reconstructed filament.

	The views are set at equally spaced angles around the z-axis.

**/
Bimage*		Bimage::filament_from_projections(double hi_res, int flag)
{
	Bimage*			pf = new Bimage(Float, TComplex, x, x, x, 1);
	pf->sampling(image->sampling()[0], image->sampling()[0], image->sampling()[0]);

	long			nn;
	double			angle(0), dang(TWOPI/n), rang(TWOPI/get_rand_max());
	View			view;
	Vector3<double> scale(1,1,1);
	Bimage*			p = NULL;
	
	for ( nn=0; nn<n; ++nn ) {
		p = extract(nn);
		if ( flag ) angle = random()*rang;
		else angle = dang*nn;
		view[0] = cos(angle);
		view[1] = sin(angle);
		view[3] = angle;
		p->fft();
		p->phase_shift_to_origin();
		pf->fspace_pack_2D(p, view.matrix(), hi_res, 0, scale, 0, 1, 0);
		delete p;
	}

	pf->phase_shift_to_center();
	pf->fft_back();
	pf->origin(pf->size()/2);
	pf->statistics();
	pf->calculate_background();

	return pf;
}



