/**
@file	mg_ctf_fit.cpp
@brief	Functions for CTF (contrast transfer function) processing
@author Bernard Heymann
@date	Created: 19970715
@date	Modified: 20210817
**/

#include "rwimg.h"
#include "mg_ctf_fit.h"
#include "mg_processing.h"
#include "simplex.h"
#include "moving_average.h"
#include "utilities.h"
#include "timer.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
double		ctf_fit_baseline(Bimage* prad, double real_size, CTFparam& em_ctf, double lores, double hires);
double		ctf_fit_envelope(Bimage* prad, double real_size, CTFparam& em_ctf, double lores, double hires);

/**
@brief 	Calculates the power spectrum radial average corrected for astigmatism.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@return Bimage*		radial average, NULL on error.

	A power spectrum with its origin at the center.
	Functions:
		angle = atan(y/x) - astigmatism_angle
		s2 = reciprocal space distance squared
		defocus_min = defocus_avg - defocus_dev
		defocus_max = defocus_avg + defocus_dev
		smin2 = 1 - defocus_dev/defocus_avg
		smax2 = 1 + defocus_dev/defocus_avg
		radius = sqrt(2*s2*(smax2*cos(angle)*cos(angle)+
					smin2*sin(angle)*sin(angle))/(smax2+smin2))
	The radial average is returned as a new 1D image.

**/
Bimage*		img_ctf_radial_average(Bimage* p, long n, CTFparam& em_ctf)
{
	if ( em_ctf.check_defocus() || em_ctf.check_Cs() )
		cerr << "in img_ctf_radial_average" << endl;
	
	long 			nrad = p->sizeX()/2;
	Bimage*			prad = new Bimage(Double, TSimple, nrad, 1, 1, 1);
	prad->sampling(p->sampling(0)[0], 1, 1);
	
	long			i, x, y, z;
	long 			iradius;
	double			dx, dy, dz, dr, radius, angle, fraction, cosang, sinang;
	double			smin2 = 1-em_ctf.defocus_deviation()/em_ctf.defocus_average();
	double			smax2 = 1+em_ctf.defocus_deviation()/em_ctf.defocus_average();
	double			yscale(p->sizeX()*1.0L/p->sizeY());
	vector<double>	num(nrad);
	double*			rdata = (double *) prad->data_pointer();
	
	if ( verbose & VERB_FULL ) {
		cout << "Radial average for CTF fitting:" << endl;
		cout << "Size:                           " << nrad << endl;
		cout << "Power spectrum origin:          " << p->image->origin() << endl;
		cout << "Defocus, deviation, angle:      " <<
			em_ctf.defocus_average() << " " << em_ctf.defocus_deviation() << " " << em_ctf.astigmatism_angle() << endl;
	}
	
	if ( p->sizeZ() < 2 ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_ctf_radial_average: doing image " << n+1 << endl;
		for ( i=0; i<nrad; i++ ) num[i] = 0;
		for ( y=0; y<p->sizeY(); y++ ) {
			dy = ((double)y - p->image[n].origin()[1])*yscale;
			for ( x=0; x<p->sizeX(); x++ ) {
				dx = (double)x - p->image[n].origin()[0];
				radius = sqrt(dx*dx + dy*dy);
				iradius = (long) radius;
				if ( iradius < nrad ) {
					angle = atan2(dy, dx) - em_ctf.astigmatism_angle();
					cosang = cos(angle);
					sinang = sin(angle);
					dr = smax2*cosang*cosang + smin2*sinang*sinang;
					if ( dr >= 0) {
						radius *= sqrt(dr);
//					if ( !isfinite(radius) )
//						cerr << "Warning: radius for " << iradius << " not finite!" << endl;
						iradius = (long) radius;
						if ( iradius < nrad - 1 ) {
							fraction = radius - iradius;
							num[iradius] += 1.0 - fraction;
							i = p->index(0,x,y,0,n);
							rdata[iradius] += (1.0 - fraction)*(*p)[i];
							iradius++;
							if ( iradius < nrad ) {
								num[iradius] += fraction;
								rdata[iradius] += fraction*(*p)[i];
							}
						}
					}
				}
			}
		}
		for ( iradius=0; iradius<nrad; iradius++ )
			if ( num[iradius] ) rdata[iradius] /= num[iradius];
			else rdata[iradius] = rdata[iradius-1];
	} else {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_ctf_radial_average: calculating 3D RPS" << endl;
		for ( i=0; i<nrad; i++ ) num[i] = 0;
		for ( z=0; z<p->sizeZ(); z++ ) {
			dz = z - p->image[n].origin()[2];
			for ( y=0; y<p->sizeY(); y++ ) {
				dy = y - p->image[n].origin()[1];
				for ( x=0; x<p->sizeX(); x++ ) {
					dx = x - p->image[n].origin()[0];
					radius = sqrt(dx*dx + dy*dy + dz*dz);
					iradius = (long) radius;
					if ( iradius < nrad - 1 ) {
						fraction = radius - iradius;
						num[iradius] += 1.0 - fraction;
						i = p->index(0,x,y,z,n);
						rdata[n*nrad+iradius] += (1.0 - fraction)*(*p)[i];
						iradius++;
						if ( iradius < nrad ) {
							num[iradius] += fraction;
							rdata[iradius] += fraction*(*p)[i];
						}
					}
				}
			}
		}
		for ( iradius=0; iradius<nrad; iradius++ )
			if ( num[iradius] ) rdata[iradius] /= num[iradius];
	}
	
	prad->statistics();
	
	return prad;
}

double		img_ctf_fit_residual(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	long			rmin = (long) (p->real_size()[0]/lores);
	long			rmax = (long) (p->real_size()[0]/hires);

	long			i, m(0);
	double			s, s2, b, e, c, d, R(0);
	
	Bimage*			prad = img_ctf_radial_average(p, n, em_ctf);
	
	for ( i=rmin; i<=rmax; ++i ) {
		s = i/p->real_size()[0];
		s2 = s*s;
		b = em_ctf.calc_baseline(s);
		e = em_ctf.calc_envelope(s);
		c = em_ctf.calculate(s2, 0);
		c *= c;
		d = b + e*c - (*prad)[i];
		d /= b;
		R += d*d;
		m++;
		cout << i << tab << (*prad)[i] << tab << b + e*c << tab << d << endl;
	}
	
	delete prad;
	
	R = sqrt(R/m);
	
	return R;
}

/**
@brief 	Determines the CTF parameters from a power spectrum.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@param	def_start	defocus search start (default 1e3).
@param	def_end		defocus search end (default 2e5).
@param	def_inc		defocus search increment (default 1e3).
@param	flag		flag to determine astigmatism.
@return double		water ring index.

	Input: Power spectrum or its logarithm.
	A radial power spectrum is calculated.
	A range of defocus values is tested (100-200000 angstrom, 0.01-20 um), 
	defining the baseline as passing through the zeroes for each defocus 
	and fitting it to a 4th order polynomial.
	The envelope function is a simple gaussian on top of the baseline and
	fitted to minimize the RMSD between the calculated curve and the
	radial power spectrum logarithm.
	The fitting is limited to the spatial frequency range between the 
	first and last zeroes.
	Defocus values are positive for underfocus.
	Functions:
		angle = atan(y/x)
		s2 = reciprocal space distance squared
		defocus_average = (defocus_max + defocus_min)/2
		defocus_deviation = (defocus_max - defocus_min)/2
		defocus = defocus_average + defocus_deviation*cos(2*(angle - astigmatism_angle))
		phase = 0.5*PI*lambda*lambda*lambda*Cs*s2*s2 - PI*lambda*defocus*s2 - amp_shift;
		CTF = sin(phase)
	The new parameters are written into the CTPparam structure.

**/
double		img_ctf_fit(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires,
				double def_start, double def_end, double def_inc, int flag)
{
	p->check_resolution(hires);
	if ( lores < hires ) lores = p->real_size()[0];
	if ( lores > p->real_size()[0] ) lores = p->real_size()[0];
	
	long		best_tile_size(def_end*em_ctf.lambda()/
					(hires*p->sampling(0)[0]));
	if ( p->size()[0] < best_tile_size ) {
		cerr << "Warning: Tile size too small! (" << p->size()[0] << " < " << best_tile_size << ")" << endl;
		cerr << "Change the tile size, defocus maximum or high resolution limit." << endl << endl;
	}
	
	if ( hires > 3 && em_ctf.baseline_type() > 3 )
		em_ctf.baseline_type(em_ctf.baseline_type()-3);

	if ( verbose & VERB_PROCESS ) {
		cout << "Determining CTF parameters:" << endl;
		cout << "Sampling:                       " << p->sampling(0)[0] << " A" << endl;
		em_ctf.show();
		cout << "Resolution range:               " << hires << " - " << lores << " A" << endl;
		cout << "Baseline type:                  " << em_ctf.baseline_type() << endl;
		cout << "Envelope type:                  " << em_ctf.envelope_type() << endl;
		cout << "Defocus min, max, inc:          " << def_start << " - " << def_end << " ∆ " << def_inc << endl << endl;
	}
	
	if ( em_ctf.defocus_average() < 1 || em_ctf.defocus_average() > 2e5 )
		em_ctf.defocus_average(2e4);
	
	em_ctf.defocus_deviation(0);
	em_ctf.astigmatism_angle(0);
	em_ctf.fom(0);

	double			fom = img_ctf_find_defocus(p, n, em_ctf, lores, hires, def_start, def_end, def_inc);

	img_ctf_fit_baseline(p, n, em_ctf, lores, hires);
	
	img_ctf_fit_envelope(p, n, em_ctf, lores, hires);
	
	long			i(0), imax(10);
	double			pfom(0);
//	double			hires_astig = 2*hires*lores/(lores+hires);
//	double			hires_dec = (hires_astig - hires)/imax;
	
	if ( flag ) for ( i=0; i<imax && ( fabs(fom - pfom) > 1e-30 ); ++i ) {
//	if ( flag ) for ( i=0; i<imax; ++i ) {
		if ( verbose )
			cout << "Defocus [" << i << "]:                    " 
					<< em_ctf.defocus_average() 
					<< " +- " << em_ctf.defocus_deviation() 
					<< " @ " << em_ctf.astigmatism_angle()*180/M_PI 
					<< " (" << fom << ")" << endl;
		pfom = fom;
//		cout << i << tab << hires_astig << endl;
		def_start = em_ctf.defocus_average()/10;
		def_end = em_ctf.defocus_average()*10;
		def_inc = def_start;
//		cout << def_start << " - " << def_end << " ∆ " << def_inc << endl;
		fom = img_ctf_fit_astigmatism(p, n, em_ctf, lores, hires);
		img_ctf_find_defocus(p, n, em_ctf, lores, hires, def_start, def_end, def_inc);
		img_ctf_fit_baseline(p, n, em_ctf, lores, hires);
		img_ctf_fit_envelope(p, n, em_ctf, lores, hires);
//		fom = img_ctf_fit_residual(p, n, em_ctf, lores, hires);
//		hires_astig -= hires_dec;
	}

	double			wri = img_water_ring_index(p, n, em_ctf);
		
	if ( verbose ) {
		cout << "Best defocus:                   " 
			<< em_ctf.defocus_average() 
			<< " +- " << em_ctf.defocus_deviation() 
			<< " @ " << em_ctf.astigmatism_angle()*180/M_PI 
			<< " (" << fom << ")" << endl;
		cout << "Water ring index:               "
			<< wri << endl;
//			<< p->meta_data()["image"][n]["water_ring"].real() << endl;
	}

//	img_ctf_isotropy(p, n, em_ctf, lores, hires);

	return wri;
}
/*
double		isotropy_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*cos(2*(simp.parameter(2)+x[i])));
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}
*/
/*
	The power between the indicated resolution shells are averaged for
	each angle and fitted to an equation for anisotropy:
		P = Pavg + Pdev*cos(2(a-phi))
	where phi is the direction of maximum power.
	The power decays with shift inaccuracy ∆x as:
		cos^2(π*∆x*s)
	The residual shift in the minimum power direction is:
		∆x = acos(sqrt(2*Pdev/Pavg))/(π*s)
*/
/*double		img_ctf_isotropy(Bimage* p, long n, double lores, double hires)
{
	if ( lores <= 0 || lores > p->real_size()[0] ) lores = p->real_size()[0];
	if ( lores < hires ) swap(lores, hires);
	if ( hires < 2*p->sampling(0)[0] ) hires = 2*p->sampling(0)[0];
	
	long			i, j, xx, yy;
	long			kmin(p->real_size()[0]/lores);
	long			kmax(p->real_size()[0]/hires);
	long			na(TWOPI*kmax);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Power spectrum anisotropy:" << endl;
		cout << "Origin:                         " << p->image[n].origin() << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Pixel radii:                    " << kmin << " - " << kmax << endl;
		cout << "Angles:                         " << na << endl;
	}
	
	double			dx, dy, r, a(0), f, y_avg(0), y_var(0), xy_scale(p->sizeX()/p->sizeY());
	vector<double>	x(na,0), y(na,0);
		
	for ( i=n*p->size().volume(), yy=0; yy<p->sizeY(); ++yy ) {
		dy = (yy - p->image[n].origin()[1])*xy_scale;
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			dx = xx - p->image[n].origin()[0];
			r = sqrt(dx*dx + dy*dy);
			if ( r >= kmin && r <= kmax ) {
				a = (na/TWOPI)*atan2(dy, dx);
				if ( a < 0 ) a += na;
				j = long(a);
				f = a - j;
				x[j] += 1-f;
				y[j] += (1-f)*(*p)[i];
				j++;
				x[j] += f;
				y[j] += f*(*p)[i];
			}
		}
	}

	for ( i=0; i<na; ++i) {
//		cout << i << tab << x[i] << tab << y[i] << endl;
		y[i] /= x[i];
		x[i] = TWOPI*i*1.0/na;
		y_avg += y[i];
		y_var += y[i]*y[i];
	}
	y_avg /= na;
	y_var /= na;
	y_var -= y_avg*y_avg;
	
//	cout << "avg=" << y_avg << " var=" << y_var << endl;

	Bsimplex		simp(1, 3, 0, na, x, y);
	
	simp.parameter(0, y_avg);
	simp.parameter(1, sqrt(y_var));
	simp.parameter(2, 0);
	simp.limits(0, 0.5*y_avg, 2*y_avg);
	simp.limits(1, sqrt(y_var), 5*sqrt(y_var));
	simp.limits(2, -M_PI, M_PI);
	
	double			R = simp.run(10000, 1e-4, isotropy_R);
	double			aniso = simp.parameter(1)/simp.parameter(0);
	double			shift_residual = acos(sqrt(2*aniso))*p->real_size()[0]*2.0/(M_PI*(kmax+kmin));
	
//	cout << "limit_lo=" << simp.limit_low(1) << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Fitted function:                " << simp.parameter(0) << " + " << simp.parameter(1) <<
			" * cos(2*(a + " << simp.parameter(2)*180/M_PI << "))" << endl;
		cout << "Anisotropy:                     " << aniso << endl;
		cout << "Estimated residual movement:    " << shift_residual << " A" << endl;
		cout << "R:                              " << R << endl;
	}

	if ( verbose & VERB_FULL ) {
		cout << "Angle\tPower\tFit" << endl;
		for ( i=0; i<na; ++i)
			cout << x[i]*180.0/M_PI << tab << y[i] << tab << simp.parameter(0) + simp.parameter(1)*cos(2*(TWOPI*i*1.0/na+simp.parameter(2))) << endl;
		cout << endl;
	}
	
	return shift_residual;
}
*/
/**
@brief 	Fits only the baseline for a given CTF.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@return double		R factor.

	A radial power spectrum is calculated incorporating astigmatism.
	The baseline of the required type is fit.
	The new parameters are written into the CTFparam structure.

**/
double		img_ctf_fit_baseline(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	Bimage* 	prad = img_ctf_radial_average(p, n, em_ctf);
	
	if ( !prad ) {
		cerr << "Error in img_ctf_fit_baseline: No radial power spectrum calculated!" << endl;
		return -1;
	}
	
	double		R = ctf_fit_baseline(prad, p->real_size()[0], em_ctf, lores, hires);
	
	delete prad;
	
	if ( verbose & VERB_PROCESS )
		em_ctf.show_baseline();
	
	return R;
}
	
/**
@brief 	Fits only the envelope for a given CTF.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	*em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@return double		R factor.

	A radial power spectrum is calculated incorporating astigmatism.
	The envelope is fit.
	The new parameters are written into the CTFparam structure.

**/
double		img_ctf_fit_envelope(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	Bimage* 	prad = img_ctf_radial_average(p, n, em_ctf);
	
	if ( !prad ) {
		cerr << "Error in img_ctf_fit_envelope: No radial power spectrum calculated!" << endl;
		return -1;
	}
	
	double		R = ctf_fit_envelope(prad, p->real_size()[0], em_ctf, lores, hires);
	
	delete prad;
	
	if ( verbose & VERB_PROCESS )
		em_ctf.show_envelope();
	
	return R;
}

double		ctf_test_defocus(CTFparam& em_ctf, double def,
				double step_size, vector<double>& r, long rmin, long rmax)
{
	long			i, nx(r.size());
	long			nr(rmax-rmin+1);
	double			fom(0);
	
	em_ctf.defocus_average(def);
	
	vector<double>	c =	em_ctf.calculate(nx, 1, step_size);

	if ( !c.size() ) {
		cerr << "Error in ctf_test_defocus: CTF not calculated!" << endl;
		return fom;
	}
	
	for ( i=rmin; i<=rmax; i++ ) fom += (c[i]*c[i] - 0.5)*r[i];
	fom /= nr;
	
	return fom;
}

double		ctf_find_defocus(vector<double>& v, CTFparam& em_ctf,
				long rmin, long rmax, double step_size,
				double def_start, double def_end, double def_inc)
{
	if ( def_start < 100 ) def_start = 100;
	if ( def_end > 2e5 ) def_end = 2e5;
	if ( em_ctf.defocus_average() < def_start || em_ctf.defocus_average() > def_end )
		em_ctf.defocus_average((def_end - def_start)/2);
	
	long			i, nx(v.size());

	for ( i=0; i<rmin; i++ ) v[i] = v[rmin];

	vector<double>	r = moving_polynomial(2, v, nx/10);


	for ( i=0; i<nx; i++ )
		if ( r[i] ) r[i] = v[i]/r[i] - 1;
		else r[i] = 0;

	double			def, def_best(em_ctf.defocus_average());
	double			fom, fom_max(-2);
	
//	if ( verbose & VERB_PROCESS ) {
		cout << "Defocus range:                  " << def_start*1e-4 << " - " << def_end*1e-4 << " um" << endl;
		cout << "Defocus increment:              " << def_inc*1e-4 << " um" << endl;
		cout << "Frequency pixel range:          " << rmin << " - " << rmax << endl;
		cout << "Frequency step size:            " << step_size << " 1/A" << endl;
//	}
	
	double			ds(def_start), de(def_end);
	
	if ( verbose & VERB_FULL )
		cout << "Defocus\tCC" << endl;
	while ( def_inc >= 50 ) {
		for ( def=ds; def<=de; def+=def_inc ) {
			fom = ctf_test_defocus(em_ctf, def, step_size, r, rmin, rmax);
			if ( fom_max < fom ) {
				fom_max = fom;
				def_best = def;
			}
			if ( verbose & VERB_FULL )
				cout << def*1e-4 << tab << fom << endl;
		}
		ds = def_best - 2*def_inc;
		de = def_best + 2*def_inc;
		def_inc /= 1.6;
		if ( ds < def_start ) ds = def_start;
		if ( de > def_end ) de = def_end;
//		cout << "def_start=" << def_start << endl;
	}

	em_ctf.defocus_average(def_best);
	if ( fom_max > -1 ) em_ctf.fom(fom_max);
	
//	vector<double>	c =	em_ctf.calculate(nx, 1, step_size);
//	for ( i=0; i<nx; i++ )
//		cout << i*step_size << tab << v[i] << tab << r[i] << tab << c[i]*c[i]-0.5 << endl;

	return em_ctf.fom();
}

double		ctf_test_defocus(CTFparam& em_ctf, double def, long nx,
				double step_size, vector<double>& r, long rmin, long rmax)
{
	long			i;
	long			nr(rmax-rmin+1);
	double			fom(0);
	
	em_ctf.defocus_average(def);
	
	vector<double>	c =	em_ctf.calculate(nx, 1, step_size);

	if ( !c.size() ) {
		cerr << "Error in ctf_test_defocus: CTF not calculated!" << endl;
		return fom;
	}
	
	for ( i=rmin; i<=rmax; i++ ) fom += (c[i]*c[i] - 0.5)*r[i];
	fom /= nr;
	
	return fom;
}

/**
@brief 	Searches for the defocus based on correlation.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit.
@param	def_start	defocus search start (default 1e3).
@param	def_end		defocus search end (default 2e5).
@param	def_inc		defocus search increment (default 1e3).
@return double		R factor.

	A radial power spectrum is calculated.
	A range of defocus values is tested (1000-200000 angstrom, 0.1-20 um), 
	and the best fitting value accepted.
	The new parameters are written into the CTFparam structure.

**/
double		img_ctf_find_defocus(Bimage* p, long n, CTFparam& em_ctf, 
				double lores, double hires, 
				double def_start, double def_end, double def_inc)
{
	if ( def_start < 100 ) def_start = 100;
	if ( def_end > 2e5 ) def_end = 2e5;
	if ( em_ctf.defocus_average() < def_start || em_ctf.defocus_average() > def_end )
		em_ctf.defocus_average((def_end - def_start)/2);
	
	p->check_resolution(hires);
	if ( lores < hires ) lores = 100;
	if ( lores > 100 ) lores = 100;
	
	Bimage*			prad = img_ctf_radial_average(p, n, em_ctf);

	long			i, nx(p->sizeX()/2);
	double			step_size = 1/p->real_size()[0];
	long			rmin = (long) (p->real_size()[0]/lores);
	long			rmax = (long) (p->real_size()[0]/hires);

	for ( i=0; i<rmin; i++ ) prad->set(i, (*prad)[rmin]);

	vector<double>	r = moving_polynomial(2, nx, (double *) prad->data_pointer(), nx/5);

	for ( i=0; i<nx; i++ )
		if ( r[i] ) r[i] = (*prad)[i]/r[i] - 1;
		else r[i] = 0;

	double			def, def_best(em_ctf.defocus_average());
	double			fom, fom_max(-2);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Defocus range:                  " << def_start*1e-4 << " - " << def_end*1e-4 << " um" << endl;
		cout << "Defocus increment:              " << def_inc*1e-4 << " um" << endl;
		cout << "Frequency pixel range:          " << rmin << " - " << rmax << endl;
		cout << "Frequency step size:            " << step_size << " 1/A" << endl;
	}
	
	double			ds(def_start), de(def_end);
	
	if ( verbose & VERB_FULL )
		cout << "Defocus\tCC" << endl;
	while ( def_inc >= 50 ) {
		for ( def=ds; def<=de; def+=def_inc ) {
			fom = ctf_test_defocus(em_ctf, def, nx, step_size, r, rmin, rmax);
			if ( fom_max < fom ) {
				fom_max = fom;
				def_best = def;
			}
			if ( verbose & VERB_FULL )
				cout << def*1e-4 << tab << fom << endl;
		}
		ds = def_best - 2*def_inc;
		de = def_best + 2*def_inc;
		def_inc /= 1.6;
		if ( ds < def_start ) ds = def_start;
		if ( de > def_end ) de = def_end;
//		cout << "def_start=" << def_start << endl;
	}
	
	delete prad;
	
	em_ctf.defocus_average(def_best);
	if ( fom_max > -1 ) em_ctf.fom(fom_max);
	
	return em_ctf.fom();
}

//Expects a radial average
Bimage*		img_ctf_radial_subtract_baseline(Bimage* prad, double real_size, CTFparam& em_ctf)
{
	long		i;
	double		s, recip_interval = 1.0/real_size;
	Bimage*		prmb = prad->copy();
	
	for ( i=0; i<prmb->sizeX(); i++ ) {
		s = i*recip_interval;
		prmb->set(i, (*prmb)[i] - em_ctf.calc_baseline(s));
	}
	
	return prmb;
}

/*
	Expects a radial power spectrum
	Measures the astigmatism between the resolution limits
	The variable real_size is the box x-size in angstrom
	Zeroes are calculated, given in spatial frequency units
	For every two successive zeroes between the resolution limits:
		the minimum within a kernel is found at each zero: [ib1]:b1 and [ib2]:b2
		the maximum halfway between the zeroes within a kernel is found: [im]:m
	Measure: sum(im*(m - (b1+b2)/2))
		the intent is that this measure should be maximum for the correct astigmatism parameters
*/
double		img_ctf_astigmatism_measure(Bimage* p, long n, double real_size, CTFparam& em_ctf, double lores, double hires)
{
	Bimage* 		prad = img_ctf_radial_average(p, n, em_ctf);

	if ( !prad ) {
		cerr << "Error in img_ctf_astigmatism_measure: No radial power spectrum calculated!" << endl;
		return -1;
	}
	
	double			pixel_size = prad->sampling(0)[0];
	if ( pixel_size <= 0 ) pixel_size = 1;
	
	long			i, ib1, ib2, im, j, nn(0);
	long			k = prad->sizeX()/100;
	if ( k < 2 ) k = 2;
	long			ilo = (long) (real_size/lores + 0.5);
	long			ihi = (long) (real_size/hires + 0.5);
	if ( ihi >= prad->sizeX() ) ihi = prad->sizeX() - 1;
	double			b1, b2, m, OM(0);
	double			max_s = 0.45/pixel_size;
	
	Bimage*			prmb = img_ctf_radial_subtract_baseline(prad, real_size, em_ctf);
	
	vector<double>	zero = em_ctf.zeroes(max_s);
	long			nz(zero.size());

	if ( nz > prad->sizeX() ) nz = prad->sizeX();

	for ( i=1; i<nz; i++ ) {
		ib1 = (long) (zero[i-1]*real_size + 0.5);
		ib2 = (long) (zero[i]*real_size + 0.5);
		im = (ib1 + ib2)/2;
		if ( im >= ilo && im <= ihi ) {
			for ( b1=(*prmb)[ib1], j=ib1-k; j<=ib1+k; j++ )
				if ( b1 > (*prmb)[j] ) b1 = (*prmb)[j];
			for ( b2=(*prmb)[ib2], j=ib2-k; j<=ib2+k; j++ )
				if ( b2 > (*prmb)[j] ) b2 = (*prmb)[j];
			for ( m=(*prmb)[im], j=im-k; j<=im+k; j++ )
				if ( m < (*prmb)[j] ) m = (*prmb)[j];
			OM += im*(m - (b1 + b2)/2);
			nn += im;
		}
	}
	
	OM /= nn;
	
	delete prad;
	delete prmb;
	
	if ( verbose & VERB_FULL )
		cout << "Fitting measure:                " << OM << endl;
	
	return OM;
}

/**
@brief 	Fits the astigmatism with a given defocus, baseline and envelope.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@return double		objective measure (larger is better).

	A radial power spectrum is calculated and the baseline subtracted.
	The defocus deviation starts from a low value to get an estimate of
	the astigmatism angle. The defocus deviation is modified nased on
	the direction of improvements in the fit, at each iteration 
	narrowing the angular search for the astigmatism angle.
	The new parameters are written into the CTFparam structure.

**/
double		img_ctf_fit_astigmatism(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	long			i, j, na;
	double			real_size = p->real_size()[0];
	double			angle, angle_min(-M_PI_2), angle_max(M_PI_2), angle_range(M_PI), angle_inc(M_PI/9);
	double			def_dev(100), ddiff, prev_fom(0);
	double			best_dev(em_ctf.defocus_deviation()), best_ang(em_ctf.astigmatism_angle());
	if ( em_ctf.defocus_deviation() > 0 ) def_dev = em_ctf.defocus_deviation();	// Use the input deviation
	double*			fom;
	CTFparam*		ec;
	
	double			best_fom = img_ctf_astigmatism_measure(p, n, real_size, em_ctf, lores, hires);
	
	if ( verbose & VERB_FULL )
		cout << "Starting astigmatism: " << em_ctf.defocus_average() << " " 
			<< em_ctf.defocus_deviation() << " " << em_ctf.astigmatism_angle() << " " << best_fom << endl;
	
	long			maxiter(20);
	for ( i=0; i<maxiter && ( best_fom/prev_fom - 1 > 1e-10 ); i++ ) {
		for ( na=0, angle=angle_min; angle<=angle_max; angle+=angle_inc ) na++;
		ec = new CTFparam[na];
		fom = new double[na];
		em_ctf.defocus_deviation(def_dev);
		for ( angle=angle_min, j=0; angle<=angle_max; angle+=angle_inc, j++ ) {
			ec[j] = em_ctf;
			ec[j].astigmatism_angle(angle);
		}
		
#ifdef HAVE_GCD
		dispatch_apply(na, dispatch_get_global_queue(0, 0), ^(size_t k){
			fom[k] = img_ctf_astigmatism_measure(p, n, real_size, ec[k], lores, hires);
		});
#else
#pragma omp parallel for
		for ( j=0; j<na; j++ ) {
			fom[j] = img_ctf_astigmatism_measure(p, n, real_size, ec[j], lores, hires);
		}
#endif

		for ( j=0; j<na; j++ ) {
			if ( best_fom < fom[j] ) {
				prev_fom = best_fom;
				best_fom = fom[j];
				best_dev = def_dev;
				best_ang = ec[j].astigmatism_angle();
			}
			if ( verbose & VERB_FULL )
				cout << em_ctf.defocus_average() << tab << def_dev << tab << ec[j].astigmatism_angle() << tab << fom[j] << tab <<
					em_ctf.defocus_average() << tab << best_dev << tab << best_ang << tab << best_fom << endl;
		}
		ddiff = best_dev - def_dev;
		
		delete[] ec;
		delete[] fom;
		
		if ( ddiff == 0 ) {
			def_dev += 10;
			angle_range /= 1.5;
			angle_min = best_ang - angle_range/2;
			angle_max = best_ang + angle_range/2;
			angle_inc /= 1.5;
		} else def_dev += ddiff/1.5;
	}
	
	if ( best_ang < -M_PI_2 ) best_ang += M_PI;
	if ( best_ang > M_PI_2 ) best_ang -= M_PI;
	
	em_ctf.defocus_deviation(best_dev);
	em_ctf.astigmatism_angle(best_ang);
	
	if ( verbose & VERB_PROCESS )
		cout << "Astigmatism refined:            " << em_ctf.defocus_average() << " +- " 
			<< em_ctf.defocus_deviation() << " @ " << em_ctf.astigmatism_angle()*180.0/M_PI << " (OM = " << best_fom << ")" << endl;

	return best_fom;
}

/**
@brief 	Calculates the water ring index from a power spectrum.
@param 	*p			power spectrum.
@param	img_num		sub-image number.
@param 	&em_ctf		CTF parameter structure.
@return double		water ing index.

	A radial power spectrum is calculated.
	The the water ring index is defined as:
		wri = wp/b - 1
	where wp is the intensity at s=0.26 (3.8 Å) and b is the background.
	The background is estimated as the average of sections before and after 
	the water peak:
		b1 between s=0.1 and s=0.2
		b2 between s=0.3 and s=0.4

**/
double		img_water_ring_index(Bimage* p, long img_num, CTFparam& em_ctf)
{
	Bimage* 		prad = img_ctf_radial_average(p, img_num, em_ctf);

	double			wri = img_water_ring_index(prad);
	
	delete prad;
	
	if ( verbose & VERB_FULL )
		cout << "Water ring index for image " << img_num << ":   " << wri << endl;

	return wri;
}

/**
@brief 	Calculates the water ring index from a power spectrum.
@param 	*prad		radial power spectrum.
@return double		water ing index.

	The the water ring index is defined as:
		wri = wp/b - 1
	where wp is the average intensity between 0.2 and 0.3, covering the 
	peak at s=0.26 (3.8 Å).
	The background, b, is estimated as the average of sections before and after 
	the water peak:
		b1 between s=0.1 and s=0.2
		b2 between s=0.3 and s=0.4

**/
double		img_water_ring_index(Bimage* prad)
{
	if ( !prad ) {
		cerr << "Error in img_water_ring_index: No radial power spectrum calculated!" << endl;
		return -1;
	}
	
	long			i, n1(0), n2(0), nw(0);
	double			s, b1(0), b2(0), wp(0), wri(0);
	double			f = 0.5/prad->real_size()[0];

	for ( i=0, s=f/2; i<prad->sizeX() && s < 0.4; ++i, s+=f ) if ( s > 0.1 ) {
		if ( s < 0.2 ) {
			b1 += (*prad)[i];
			n1++;
		} else if ( s > 0.3 ) {
			b2 += (*prad)[i];
			n2++;
		} else {
			wp += (*prad)[i];
			nw++;
		}
	}

	n1 += n2;
	if ( n1 ) b1 = (b1 + b2)/n1;
	if ( nw ) wp /= nw;
	
	if ( b1 ) wri = wp/b1 - 1;
		
	return wri;
}

double		gaussian_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df, s2;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		s2 = x[i]*x[i];
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*exp(simp.parameter(2)*s2));
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}

double		gaussian1_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df, ds, w, wsum(0);
	double			cen(M_PI_2 / ( simp.points() - 1));
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		w = 1 - cos(cen*i);
		ds = simp.parameter(3) - x[i];
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*x[i] + 
			simp.parameter(2)*exp(simp.parameter(4)*ds*ds));
		R += w*df*df;
		wsum += w;
	}
	
	R = sqrt(R/wsum);
			
	return R;
}

double		double_gaussian_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df, s2;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		s2 = x[i]*x[i];
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*exp(simp.parameter(2)*s2) + 
			simp.parameter(3)*exp(simp.parameter(4)*s2));
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}

double		eman_baseline_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();

	for ( i=0; i<simp.points(); i++ ) {
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*
			exp(simp.parameter(2)*sqrt(x[i]) + simp.parameter(3)*x[i]*x[i]));
		R += df*df;
	}
		
	R = sqrt(R/i);
	
	return R;
}

/*
@brief 	Fits a double gaussian function.
@param 	n		number of data points.
@param 	*x		x array (at least order+1 values).
@param	*y		y array (at least order+1 values).
@param	coeff	array in which 5 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Double gaussian function:
		y = a + b*exp(c*s^2) + d*exp(f*s^2)
**/
double 		ctf_fit_polynomial(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long			nc(5);
	if ( n < 5 ) nc = n - 1;
	
	double			R = fit_polynomial(n, x, y, nc, coeff);
	
	return R;
}

/*
@brief 	Fits a double gaussian function.
@param 	n		number of data points.
@param 	*x		x array (at least order+1 values).
@param	*y		y array (at least order+1 values).
@param	coeff	array in which 5 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Double gaussian function:
		y = a + b*exp(c*s^2) + d*exp(f*s^2)
**/
double 		ctf_fit_double_gaussian_baseline(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;
	
	Bsimplex			simp(1, 5, 0, n, x, y);
	
	simp.parameters(5, coeff);
	simp.limits(0, 0, 1000);
	simp.limits(1, 0, 100);
	simp.limits(2, -1000, -0.1);
	
	if ( n > 5 ) {
		simp.limits(3, 0, 1000);
		simp.limits(4, -10000, -1);
	}
	
	double			R = simp.run(10000, 0.01, double_gaussian_R);
	
	for ( i=0; i<5; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

/*
@brief Fits an EMAN-style baseline function.
@param n		number of data points.
@param *x		x array (at least order+1 values).
@param *y		y array (at least order+1 values).
@param coeff	array in which 4 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	EMAN baseline function:
		y = a + b*exp( c*sqrt(s) + d*s^2 )
**/
double 		ctf_fit_eman_to_baseline(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;

	Bsimplex			simp(1, 4, 0, n, x, y);
	
	simp.parameters(4, coeff);
	simp.limits(0, 0, 1000);
	simp.limits(1, 0, 1000);
	simp.limits(2, -100, 100);
	simp.limits(3, -100, 0);
	
	double			R = simp.run(10000, 0.001, eman_baseline_R);
	
	for ( i=0; i<4; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

/*
@brief Fits a gaussian function.
@param n		number of data points.
@param *x		x array (at least order+1 values).
@param *y		y array (at least order+1 values).
@param coeff	array in which 3 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Gaussian function:
		y = a + b*exp(c*s^2)
**/
double 		ctf_fit_gaussian(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;
	
	Bsimplex			simp(1, 3, 0, n, x, y);
	
	simp.parameters(3, coeff);
	if ( coeff[0] ) simp.limits(0, 0, 1000);
	simp.limits(1, 0.001, 100);
	simp.limits(2, -1e4, -0.1);
	
	double			R = simp.run(10000, 0.01, gaussian_R);
	
	for ( i=0; i<3; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

double		ctf_fit_baseline_bump(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i, ilo(0), ns;
	double				slo(0.2), shi(0.3);
	
	coeff[5] = 0;		// Amplitude
	coeff[6] = 0.265;	// Location
	coeff[7] = -5e3;	// 1/2Width^2
	
	vector<double>		xb, yb;
	
	for ( i=ns=0; i<n; i++ ) if ( x[i] > slo && x[i] < shi ) {
		if ( !ilo ) ilo = i;
		if ( coeff[5] < y[i] ) {	 // Find the maximum
			coeff[5] = y[i];
			coeff[6] = x[i];
		}
		xb.push_back(x[i]);
		yb.push_back(y[i]);
		ns++;
//		cout << x[i] << tab << y[i] << endl;
	}
	
	double		slope((y[ilo+ns-1] - y[ilo])/(x[ilo+ns-1] - x[ilo]));
	double		inter(y[ilo] - slope*x[ilo]);
		
	coeff[5] -= 0.5*(y[ilo] + y[ilo+ns-1]);
	if ( coeff[5] < 1e-3 ) coeff[5] = 0.1*y[ilo];
	if ( coeff[6] < 0.25 || coeff[6] > 0.28 ) coeff[6] = 0.265;
		
	Bsimplex			simp(1, 5, 0, ns, xb, yb);
	
//	simp.parameters(5, &coeff[3]);
	simp.parameter(0, inter);
	simp.parameter(1, slope);
	simp.parameter(2, coeff[5]);
	simp.parameter(3, coeff[6]);
	simp.parameter(4, coeff[7]);
	simp.limits(0, 0.1*inter, 10*inter);
	simp.limits(1, -5*fabs(slope), 5*fabs(slope));
	simp.limits(2, 0.01*coeff[5], 10*coeff[5]);
	simp.limits(3, 0.25, 0.28);
	simp.limits(4, -2e4, -1e3);
	
//	for ( i=0; i<5; i++ ) cout << tab << simp.parameter(i);
//	cout << endl;
	
	double			R = simp.run(10000, 0.0001, gaussian1_R);
	
	for ( i=2; i<5; i++ ) coeff[i+3] = simp.parameter(i);
	
//	for ( i=0; i<5; i++ ) cout << tab << simp.parameter(i);
//	cout << tab << R << endl;
	
	return R;
}

/*
@brief 	Fits a double gaussian function.
@param	n		number of data points.
@param	*x		x array (at least order+1 values).
@param	*y		y array (at least order+1 values).
@param	coeff	array in which 5 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Double gaussian function:
		y = a + b*exp(c*s^2) + d*exp(f*s^2)
		with a = 0
**/
double 		ctf_fit_double_gaussian_envelope(long n, vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;
	
	Bsimplex			simp(1, 5, 0, n, x, y);
	
	simp.parameters(5, coeff);
	if ( coeff[0] ) simp.limits(0, 0, 1000);
	simp.limits(1, 0.001, 100);
	simp.limits(2, -500, -0.1);
	
	if ( n > 5 ) {
		simp.limits(3, 0.01, 100);
		simp.limits(4, -5000, -1);
	}
	
	double			R = simp.run(10000, 0.01, double_gaussian_R);
	
	for ( i=0; i<5; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

double		ctf_fit_baseline(Bimage* prad, double real_size, CTFparam& em_ctf, double lores, double hires)
{
	double			pixel_size = prad->sampling(0)[0];
	if ( pixel_size <= 0 ) pixel_size = 1;
	
	long			i, j, ib, nc(4);
	long			istart = (long) (real_size/lores);
	long			iend = (long) (real_size/hires);
	if ( istart < 1 ) istart = 1;
	if ( iend >= prad->sizeX() ) iend = prad->sizeX() - 1;
	long			w = (iend - istart)/10;
	
	double			recip_interval = 1.0/real_size;
	double			max_s = 0.45/pixel_size;
	long			k = prad->sizeX()/100;
	if ( k < 1 ) k = 1;
	
	double			R = 1e30;
	vector<double>	coeff(NCTFPARAM);
	for ( i=0; i<NCTFPARAM; i++ ) coeff[i] = em_ctf.baseline(i);
	
	vector<double>	s(prad->sizeX(), 0);
	vector<double>	d(prad->sizeX(), 0);
	
	vector<double>	zero = em_ctf.zeroes(max_s);
	long			n(zero.size());
	
	if ( n > prad->sizeX() ) n = prad->sizeX();

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_fit_baseline: zeroes=" << n << endl;
	
	double			min = 1e10;
	double*			mavg = NULL;
	
	if ( n > 3 ) {		// Enough zeroes to fit the baseline
		for ( i=0; i<n; i++ ) {
			ib = (long) (zero[i]*real_size + 0.5);
			s[i] = zero[i];
			d[i]=(*prad)[ib];
			for ( j=(ib>k)?ib-k:0; j<=ib+k && j<prad->sizeX(); j++ )	// minimum within kernel
				if ( d[i] > (*prad)[j] ) d[i] = (*prad)[j];
		}
	} else {			// Too few zeroes, use local minima
		mavg = moving_average(prad->sizeX(), (double *) prad->data_pointer(), prad->sizeX()/10);
		for ( i=istart, n=0; i<iend; i++ ) {
			if ( min > (*prad)[i] - mavg[i] ) {
				min = (*prad)[i] - mavg[i];
				s[n] = recip_interval*i;
				d[n] = (*prad)[i];
				if ( d[n] > mavg[i] ) d[n] = mavg[i];
			}
			if ( i%w == 0 ) {
				min = 1e10;
				n++;
			}
		}
		s[n] = 10;
		d[n] = 0;
		n++;
		delete[] mavg;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_fit_baseline: fitting baseline type " << em_ctf.baseline_type() << endl;
	
	if ( em_ctf.baseline_type() > 3 ) {
		ctf_fit_baseline_bump(n, s, d, coeff);
		for ( i=0; i<n; i++ )
			d[i] -= coeff[5] * exp(coeff[7] * (coeff[6] - s[i]) * (coeff[6] - s[i]));
	}
	
	switch ( em_ctf.baseline_type() ) {
		case 1:
		case 4:
			if ( n < 5 ) nc = n - 1;
			R = fit_polynomial(n, s, d, nc, coeff);
			break;
		case 2:
		case 5:
			R = ctf_fit_double_gaussian_baseline(n, s, d, coeff);
			break;
		case 3:
		case 6:
			R = ctf_fit_eman_to_baseline(n, s, d, coeff);
			break;
		default:
			cerr << "Error: Baseline type " << em_ctf.baseline_type() << " not supported" << endl;
	}
	
	if ( em_ctf.baseline_type() > 0 && em_ctf.baseline_type() < 7 )
		em_ctf.baseline(coeff);

	return R;
}

double		ctf_fit_envelope(Bimage* prad, double real_size, CTFparam& em_ctf, double lores, double hires)
{
	double		pixel_size = prad->sampling(0)[0];
	if ( pixel_size <= 0 ) pixel_size = 1;
	
	long		i, j, ib;
	long		k = prad->sizeX()/80;
	if ( k < 1 ) k = 1;
	
	double		max_s = 0.45/pixel_size;
	double		R = 1e30;
	
	vector<double>	zero = em_ctf.zeroes(max_s);
	long			nzero(zero.size());

	if ( nzero > prad->sizeX() ) nzero = prad->sizeX();
	
	double			b, si;
	vector<double>	coeff(NCTFPARAM);
	vector<double>	s(nzero);
	vector<double>	d(nzero);
	for ( i=0; i<NCTFPARAM; i++ ) coeff[i] = 0;
	for ( i=0; i<nzero; i++ ) s[i] = d[i] = 0;
	
	ib = (long) (zero[0]*real_size + 0.5);
	d[0] = 0;
	for ( i=ib/2; i<ib; i++ ) {
		si = i/real_size;
		b = (*prad)[i] - em_ctf.calc_baseline(si);
		if ( d[0] < b ) {
			d[0] = b;
			s[0] = si;
		}
	}
	
	for ( i=1; i<nzero; i++ ) {
		s[i] = (zero[i-1]+zero[i])*0.5;
		ib = (long) (s[i]*real_size + 0.5);
		b = em_ctf.calc_baseline(s[i]);
		d[i] = (*prad)[ib] - b;
		for ( j=(ib>k)?ib-k:0; j<=ib+k && j<prad->sizeX(); j++ )	// maximum within kernel
			if ( d[i] < (*prad)[j] - b ) d[i] = (*prad)[j] - b;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_envelope: type " << em_ctf.envelope_type() << endl;
		for ( i=0; i<5; i++ ) cout << tab << em_ctf.envelope(i);
		cout << endl;
	}

	switch ( em_ctf.envelope_type() ) {
		case 1:		// Single gaussian
			for ( i=0; i<2; i++ ) coeff[i+1] = em_ctf.envelope(i);
			R = ctf_fit_gaussian(nzero-1, s, d, coeff);
			for ( i=0; i<2; i++ ) coeff[i] = coeff[i+1];
			coeff[i] = 0;
			break;
		case 2:		// Single gaussian with constant
			for ( i=0; i<3; i++ ) coeff[i] = em_ctf.envelope(i);
			R = ctf_fit_gaussian(nzero-1, s, d, coeff);
			break;
		case 3:		// Double gaussian
			for ( i=0; i<4; i++ ) coeff[i+1] = em_ctf.envelope(i);
			R = ctf_fit_double_gaussian_envelope(nzero, s, d, coeff);
			for ( i=0; i<4; i++ ) coeff[i] = coeff[i+1];
			coeff[i] = 0;
			break;
		case 4:		// Double gaussian with constant
			for ( i=0; i<5; i++ ) coeff[i] = em_ctf.envelope(i);
			R = ctf_fit_double_gaussian_envelope(nzero, s, d, coeff);
			break;
		default: break;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_envelope:";
		for ( i=0; i<5; i++ ) cout << tab << coeff[i];
		cout << endl;
	}
	
	em_ctf.envelope(coeff);
	
	return R;
}

