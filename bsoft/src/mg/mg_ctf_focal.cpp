/**
@file	mg_ctf_focal.cpp
@brief	Processing focal series
@author 	Bernard Heymann
@date	Created: 20220808
@date	Modified: 20230127
**/

#include "rwimg.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "mg_ctf_focal.h"
#include "simplex.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates an aberration image.
@param 	cp				CTF & aberration parameters.
@param 	def_min			Minimum defocus.
@param 	def_max			Maximum defocus.
@param 	def_inc			Defocus increment.
@param 	size			new image size.
@param 	sam				new image pixel size.
@param 	lores			low resolution limit.
@param 	hires			high resolution limit.
@return Bimage*			new complex CTF function image.

	Functions:
		angle = atan(y/x)
		s2 = x*x + y*y
		defocus_average = (defocus_max + defocus_min)/2
		defocus_deviation = (defocus_max - defocus_min)/2
		defocus = defocus_average + defocus_deviation*cos(2*(angle - astigmatism_angle))
		phase = 0.5*PI*lambda*lambda*lambda*Cs*s2*s2 - PI*lambda*defocus*s2 - amp_shift;
		CTF = sin(phase)
	Note: Defocus is positive for underfocus and negative for overfocus.

**/
Bimage*		img_ctf_gradient(CTFparam& cp, double def_min, double def_max, double def_inc,
				Vector3<long> size, Vector3<double> sam, double lores, double hires)
{
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;
	
	cp.defocus_average((def_min+def_max)/2);
	
	double			shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	double			shi2(shi*shi), slo2(slo*slo);
	
	Bimage*			p = new Bimage(Float, TComplex, size, 1);
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	long 			i, x, y, z, n(0);
	double			d, sx, sy, sz, s2, a, dphi;
	Complex<double>	cv(1,0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF gradient function:" << endl;
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Defocus min, max, inc:          " << def_min << " - " << def_max << " ∆ " << def_inc << endl;
		cout << "First sinc node:                " << sqrt((def_max-def_min)*cp.lambda()/2) << " A" << endl;
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[1]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[1]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}
	
	for ( d = def_min, n=0; d <= def_max; d += def_inc, ++n ) {
		cp.defocus_average(d);
		for ( i=z=0; z<p->sizeZ(); ++z ) {
			sz = z;
			if ( z > h[2] ) sz -= p->sizeZ();
			sz *= freq_scale[2];
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); ++x, ++i ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					if ( s2 >= slo2 && s2 <= shi2 ) {
						a = atan2(sy,sx);
						dphi = cp.calculate_aberration_even(s2, a);
						cv = cp.aberration_odd_complex(s2, a);
						p->add(i, cv.conj() * sinl(dphi));
					}
				}
			}
		}
	}
	
	p->multiply(1.0L/n);
	
//	p->complex_to_real();
//	p->statistics();
	
//	write_img("c.grd", p, 0);
//	bexit(0);
	
	return p;
}

Bimage*		img_ctf_focal_series(CTFparam& cp, double def_start, double def_end, double def_inc,
				Vector3<long> size, Vector3<double> sam, double lores, double hires)
{
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;
	if ( def_start > def_end && def_inc > 0 ) def_inc = -def_inc;
	if ( def_start < def_end && def_inc < 0 ) def_inc = -def_inc;
	
	double			def(cp.defocus_average());
	cp.defocus_average(def_start+def_inc*size[2]/2.0);
	
//	cout << def_start << tab << def_end << tab << def_inc << endl;
	
	long			nimg((def_end-def_start)/def_inc);
	if ( size[2] > 1 ) nimg = 1;
	
	double			d, shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	double			shi2(shi*shi), slo2(slo*slo);
	
	Bimage*			p = new Bimage(Float, TComplex, size, nimg);
	if ( sam.volume() > 0 ) p->sampling(sam);
	p->origin(0,0,size[2]/2);
	
	long 			i, x, y, z, n;
	double			sx, sy, s, s2, a, dphi;
	Complex<double>	cv(1,0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF focal series:" << endl;
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Defocus start, end, increment:  " << def_start << " - " << def_end << " ∆ " << def_inc << endl;
		cout << "First sinc node:                " << sqrt(fabs(def_end-def_start)*cp.lambda()/2) << " A" << endl;
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		cp.show_envelope();
		cout << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[1]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[1]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}
	
	for ( d = def_start, i=n=0; n<nimg; ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			d += def_inc;
			cp.defocus_average(d);
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); ++x, ++i ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy;
					if ( s2 >= slo2 && s2 <= shi2 ) {
						s = sqrt(s2);
						a = atan2(sy,sx);
						dphi = cp.calculate_aberration_even(s2, a);
						cv = cp.aberration_odd_complex(s2, a);
						cv = cv.conj() * (sinl(dphi) * cp.calc_envelope(s));
						p->add(i, cv);
					}
				}
			}
		}
	}
	
//	p->multiply(1.0L/n);

	cp.defocus_average(def);

	return p;
}

/**
@brief 	Extracts a transverse section from focal series power spectra.
@param 	*p				Focal series.
@param 	which			0=x, 1=y.
@return Bimage*			transverse section.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_extract_section(Bimage* p, int which)
{
	long			i, j, xx, yy, zz;
	Vector3<long>	size(1,p->sizeZ(),1);
	if ( which ) size[0] = p->sizeY();
	else size[0] = p->sizeX();
	
	Bimage*			ps = new Bimage(Float, TSimple, size, 1);
	if ( which ) ps->sampling(p->sampling(0)[1], p->sampling(0)[2], 1);
	else ps->sampling(p->sampling(0)[0], p->sampling(0)[2], 1);
	ps->origin(ps->size()/2);

	if ( which ) {
		for ( j=zz=0; zz<p->sizeZ(); ++zz ) {
			xx = p->sizeX()/2;
			i = zz*p->sizeY()*p->sizeX() + xx;
			for ( yy=0; yy<p->sizeY(); ++yy, i+=p->sizeX(), ++j )
				ps->set(j, (*p)[i]);
		}
	} else {
		for ( j=zz=0; zz<p->sizeZ(); ++zz ) {
			yy = p->sizeY()/2;
			i = (2*zz+1)*yy*p->sizeX();
			for ( xx=0; xx<p->sizeX(); ++xx, ++i, ++j )
				ps->set(j, (*p)[i]);
		}
	}
	
	return ps;
}

/**
@brief 	Calculates a fit for a defocus value to a section from focal series modified power spectra.
@param 	*p				Section.
@param 	cp				CTF parameters.
@param 	def				defocus to test for.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_ctf_section_fit(Bimage* p, CTFparam cp, double def, double hires, double lores)
{
	double			df = p->sampling(0)[1];
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), v, c, cs(0), rs(0), CC(0);
//	long			h(p->sizeX()/2);
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);

	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		cp.defocus_average(def - df*(yy - p->sizeY()/2));
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				c = 2*c*c - 1;
				v = (*p)[i];
				CC += c*v;
				cs += c*c;
				rs += v*v;
			}
		}
	}

	CC /= sqrt(cs*rs);

	return CC;
}

/**
@brief 	Determines average defocus from a transverse section of focal series power spectra.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			Correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_find_section_defocus(Bimage* p, CTFparam& cp, double hires, double lores)
{
	double			df = p->sampling(0)[1];
	double			def, def_best(0), CC(0), CC_best(0);
	double			def_start(cp.defocus_average() - 1000);
	double			def_end(cp.defocus_average() + 1000);
	
	for ( def=def_start; def<=def_end; def+=df ) {
		CC = img_ctf_section_fit(p, cp, def, hires, lores);
		if ( verbose & VERB_FULL )
			cout << def << tab << CC << endl;
		if ( CC_best < CC ) {
			CC_best = CC;
			def_best = def;
		}
	}
	
	if ( verbose ) {
		cout << "Best defocus:                   " << def_best << endl;
		cout << "Correlation coefficient:        " << CC_best << endl << endl;
	}
	
	cp.defocus_average(def_best);

	return CC_best;
}

/**
@brief 	Calculates a transverse section of focal series power spectra from CTF parameters.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param 	res				High resolution limit.
@return Bimage*			Transverse section.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_ctf_section_calc(Bimage* p, CTFparam& cp, double res)
{
	double			def = cp.defocus_average();
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), c;
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smax = (res)? 1/res: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose )
		cout << "Calculating the section CTF for defocus " << def << " A" << endl << endl;

	Bimage*			pc = p->copy();
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		cp.defocus_average(def - p->sampling(0)[1]*(yy - p->sizeY()/2));
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				pc->set(i, c*c);
			} else {
				pc->set(i, 0);
			}
		}
	}
	
	cp.defocus_average(def);

	return pc;
}

double		focus_cs_amp_section_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			fac = -simp.constant(0)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i, ++j ) {
		v = simp.parameter(0);				// Amp
		v += simp.parameter(1) * x[j];		// Defocus
		v += simp.parameter(2)*x[j]*x[j];	// Cs
		w = exp(fac*x[j]);					// Envelope
		v += x[++j];						// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focus_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			fac = -simp.constant(2)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		v = simp.parameter(0) * x[j++];	// Average defocus
		w = exp(fac*x[j]);				// Envelope
		v += simp.constant(0) + simp.constant(1)*x[j]*x[j];	// Amp & Cs
		v += simp.parameter(1) * x[j++];	// Astigmatism
		v += simp.parameter(2) * x[j++];	// Astigmatism
		v += x[j++];					// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focal_aberration_fit_R(Bsimplex& simp)
{
	long			i, j, k;
	double			v, v2(0), f2(0), w(1), CC(0);
	long			iw = simp.constant(0);		// Index for s2
	double			fac = -simp.constant(1)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		w = exp(fac*x[j+iw]);			// Envelope
		for ( k=0, v=0; k<simp.parameters(); ++k )
			v += simp.parameter(k) * x[j++];
		v += x[j++];					// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

/**
@brief 	Fits 3 CTF parameters to a transverse section of focal series power spectra.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.
	
	The 3 parameters are the isotropic aberrations:
		constant phase shift (amplitude contrast)
		defocus
		spherical aberration
		

**/
double		img_ctf_fit_section(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting section:" << endl;

	cout << "origin = " << p->image->origin() << endl;
	
	long			i, xx, yy;
	double			wl(cp.lambda());
	double			dd, s2;
	vector<double>	t, terms, val;
		
	for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
		dd = M_PI*wl*p->image->sampling()[1]*(yy-p->sizeY()/2);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s2 = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s2 >= 0.5 ) s2 -= 1;
			s2 /= p->image->sampling()[0];
			s2 *= s2;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				terms.push_back(s2);
				terms.push_back(dd*s2);
				val.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 3, 1, val.size(), terms, val);
	simp.constant(0, Bfactor);
	simp.parameter(0, cp.aberration_weight(0,0));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, cp.aberration_weight(4,0));
	simp.limits(0, cp.aberration_weight(0,0)-0.2, cp.aberration_weight(0,0)+0.2);
	simp.limits(1, cp.aberration_weight(2,0)-50, cp.aberration_weight(2,0)+50);
	simp.limits(2, cp.aberration_weight(4,0)-50, cp.aberration_weight(4,0)+50);

	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_cs_amp_section_fit_R);
	CC = 1 - iR;
	cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_cs_amp_section_fit_R, 100);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(0,0,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(4,0,simp.parameter(2));

	if ( verbose ) {
		cout << "Amplitude phase:                " << cp.amp_shift()*180.0/M_PI << " degrees" << endl;
		cout << "Defocus average:                " << cp.defocus_average()*1e-4 << " um" << endl;
		cout << "Cs:                             " << cp.Cs()*1e-7 << " mm" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits defocus and astigmatism to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.

	The constant phase shift and spherical aberration are fixed at initial values.
	
**/
double		img_ctf_fit_astigmatism(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting astigmatism:" << endl;

	map<pair<long,long>,double>		wa;
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);

//	for ( auto w: wa )
//		cout << w.first.first << tab << w.first.second << tab << w.second << endl;
	cout << "origin = " << p->image->origin() << endl;
	
	long			nt(wa.size()), nc(3);
	
	long			i, j, xx, yy, zz;
	double			s2, dd;
	double			wl(cp.lambda());
	Vector3<double>	s;
	vector<double>	t, terms, val;
		
	for ( i=0, zz=0; zz<p->sizeZ(); ++zz ) {
//		dd = p->image->sampling()[2]*(zz-p->sizeZ()/2);
//		cout << zz << tab << dd << endl;
		dd = M_PI*wl*p->image->sampling()[2]*(zz-p->sizeZ()/2);
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
			if ( s[1] >= 0.5 ) s[1] -= 1;
			s[1] /= p->image->sampling()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
				if ( s[0] >= 0.5 ) s[0] -= 1;
				s[0] /= p->image->sampling()[0];
				s2 = s.length2();
				if ( s2 >= smin2 && s2 <= smax2 ) {
					t = aberration_terms(wa, s[0], s[1]);
					for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
					terms.push_back(dd*s2);
					val.push_back((*p)[i]);
				}
			}
		}
	}

//	cout << "independent variable size = " << terms.size() << endl;
//	cout << "dependent variable size =   " << val.size() << endl;
/*
	for ( i=j=0; i<val.size(); ++i ) {
		cout << i;
		for ( long k=0; k<4; ++k ) cout << tab << terms[j++];
		cout << tab << val[i] << endl;
	}
*/
	double				dlim(10);
	Bsimplex			simp(nt+1, nt, nc, val.size(), terms, val);
	simp.constant(0, cp.aberration_weight(0,0));
	simp.constant(1, cp.aberration_weight(4,0));
	simp.constant(2, Bfactor);
	simp.parameter(0, cp.aberration_weight(2,-2));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, cp.aberration_weight(2,2));
	simp.limits(0, -dlim, dlim);
	simp.limits(1, cp.aberration_weight(2,0)-10*dlim, cp.aberration_weight(2,0)+10*dlim);
	simp.limits(2, -dlim, dlim);

	if ( verbose ) {
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
		j = 0;
		for ( auto w: wa ) {
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
			j++;
		}
		cout << endl;
	}
	
	simp.show();
	
	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_fit_R);
	CC = 1 - iR;
	cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_fit_R, 100);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(2,-2,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(2,2,simp.parameter(2));

	if ( verbose ) {
		cout << "Defocus average:                " << cp.defocus_average() << " A" << endl;
		cout << "Defocus deviation:              " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180.0/M_PI << " °" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param 	wl				Electron wavelength.
@param 	&wa				Aberration weights (replaced).
@param 	wd				Aberration weight limits for fitting.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_aberration(Bimage* p, double wl, map<pair<long,long>,double>& wa,
				map<pair<long,long>,double> wd, double hires, double lores,
				double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting aberrations:" << endl;

//	for ( auto w: wa )
//		cout << w.first.first << tab << w.first.second << tab << w.second << endl;
//	cout << "origin = " << p->image->origin() << endl;
	
	long			nt(wa.size());
	
	long			i, j, xx, yy, zz;
	double			s2, dd;
	Vector3<double>	s;
	vector<double>	t, terms, val;
		
	for ( i=0, zz=0; zz<p->sizeZ(); ++zz ) {
//		dd = p->image->sampling()[2]*(zz-p->sizeZ()/2);
//		cout << zz << tab << dd << endl;
		dd = M_PI*wl*p->image->sampling()[2]*(zz-p->sizeZ()/2);
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
			if ( s[1] >= 0.5 ) s[1] -= 1;
			s[1] /= p->image->sampling()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
				if ( s[0] >= 0.5 ) s[0] -= 1;
				s[0] /= p->image->sampling()[0];
				s2 = s.length2();
				if ( s2 >= smin2 && s2 <= smax2 ) {
					t = aberration_terms(wa, s[0], s[1]);
					for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
					terms.push_back(dd*s2);
					val.push_back((*p)[i]);
				}
			}
		}
	}
/*
	cout << "independent variable size = " << terms.size() << endl;
	cout << "dependent variable size =   " << val.size() << endl;

	for ( i=j=0; i<val.size(); ++i ) {
		cout << i;
		for ( long k=0; k<=nt; ++k ) cout << tab << terms[j++];
		cout << tab << val[i] << endl;
	}
*/
	Bsimplex			simp(nt+1, nt, 2, val.size(), terms, val);
	simp.constant(0, 2);		// Index of s2 for weighting function
	simp.constant(1, Bfactor);	// B-factor of weighting function
	
	i = 0;
	for ( auto w: wa )
		simp.parameter(i++, w.second);

	i = 0;
	for ( auto w: wd ) {
		simp.limits(i, simp.parameter(i)-w.second, simp.parameter(i)+w.second);
		i++;
	}

	if ( verbose ) {
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
		j = 0;
		for ( auto w: wa ) {
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
			j++;
		}
		cout << endl;
	}
	
	simp.show();
	
	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focal_aberration_fit_R);
	CC = 1 - iR;
	cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focal_aberration_fit_R, 100);

	i = 0;
	for ( auto& w: wa )
		w.second = simp.parameter(i++);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits 5 even aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_even5(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	map<pair<long,long>,double>		wa, wd;
	wa[{0,0}] = cp.aberration_weight(0,0);
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);
	wa[{4,0}] = cp.aberration_weight(4,0);
	wd[{0,0}] = 0.1;
	wd[{2,-2}] = 100;
	wd[{2,0}] = 1000;
	wd[{2,2}] = 100;
	wd[{4,0}] = 100;

	double		CC = img_ctf_fit_aberration(p, cp.lambda(), wa, wd, hires, lores, Bfactor, maxiter);
	
	cp.aberration_weights(wa);
/*
	if ( verbose ) {
		cout << "Defocus average:                " << cp.defocus_average() << " A" << endl;
		cout << "Defocus deviation:              " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180.0/M_PI << " °" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}
*/
	return CC;
}

/**
@brief 	Fits 9 even aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_even9(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	map<pair<long,long>,double>		wa, wd;
	wa[{0,0}] = cp.aberration_weight(0,0);
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);
	wa[{4,-4}] = cp.aberration_weight(4,-4);
	wa[{4,-2}] = cp.aberration_weight(4,-2);
	wa[{4,0}] = cp.aberration_weight(4,0);
	wa[{4,2}] = cp.aberration_weight(4,2);
	wa[{4,4}] = cp.aberration_weight(4,4);
	wd[{0,0}] = 0.1;
	wd[{2,-2}] = 100;
	wd[{2,0}] = 1000;
	wd[{2,2}] = 100;
	wd[{4,-4}] = 10;
	wd[{4,-2}] = 10;
	wd[{4,0}] = 100;
	wd[{4,2}] = 10;
	wd[{4,4}] = 10;

	double		CC = img_ctf_fit_aberration(p, cp.lambda(), wa, wd, hires, lores, Bfactor, maxiter);
	
	cp.aberration_weights(wa);
/*
	if ( verbose ) {
		cout << "Defocus average:                " << cp.defocus_average() << " A" << endl;
		cout << "Defocus deviation:              " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180.0/M_PI << " °" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}
*/
	return CC;
}

/**
@brief 	Fits the CTF to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
	The image should be 3D with the third sampling interval the
	change in focus per 2D image.
**/
Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	p->set_hi_lo_resolution(hires, lores);
	
	Bimage*			pf = img_ctf_fit_prepare(p, 10);

	if ( verbose ) {
		cout << "Fitting focal series power spectra:" << endl;
		cout << "Focus change:                   " << p->sampling(0)[2] << " A" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "B-factor:                       " << Bfactor << " A2" << endl;
		cp.show();
	}

	Bimage*			ps = img_extract_section(pf, 0);

	double			Rx = img_find_section_defocus(ps, cp, hires, lores);

	Rx = img_ctf_fit_section(ps, cp, hires, lores, Bfactor, maxiter);

	double			dx = cp.defocus_average();

	Bimage*			psc = img_ctf_section_calc(ps, cp, hires);

	ps->replace_half(psc);
	
	write_img("px.grd", ps, 0);
	
	delete psc;
	delete ps;
	
	ps = img_extract_section(pf, 1);

	double			Ry = img_find_section_defocus(ps, cp, hires, lores);
	
	Ry = img_ctf_fit_section(ps, cp, hires, lores, Bfactor, maxiter);

	double			dy = cp.defocus_average();

	psc = img_ctf_section_calc(ps, cp, hires);

	ps->replace_half(psc);
	
	write_img("py.grd", ps, 0);
	
	delete psc;
	delete ps;
	
	if ( verbose ) {
		cout << "Best defocus in x:              " << dx << " A (" << Rx << ")" << endl;
		cout << "Best defocus in y:              " << dy << " A (" << Ry << ")" << endl;
	}
	
	cp.defocus_average(0.5*(dx+dy));

	img_ctf_fit_astigmatism(pf, cp, hires, lores, Bfactor, maxiter);

	img_ctf_fit_even5(pf, cp, hires, lores, Bfactor, maxiter);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	img_ctf_fit_even9(pf, cp, hires, lores, Bfactor, maxiter);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	double			def_inc = -p->sampling(0)[2];
	double			def_start = cp.defocus_average() - def_inc*p->sizeZ()/2;
	double			def_end = def_start + def_inc*(p->sizeZ() - 1);
	
	ps = img_ctf_focal_series(cp, def_start, def_end, def_inc,
				p->size(), p->image->sampling(), 0, hires/2);

	ps->complex_to_intensities();
	
	ps->shift_wrap(Vector3<double>(p->sizeX()/2, p->sizeY()/2, 0));
	ps->multiply(2);
	ps->add(-1);
	
	pf->replace_half(ps);
	
	delete ps;
	
	return pf;
}

/**
@brief 	Extracts a sphere corresponding to a given acceleration voltage.
@param 	*p				Focal series.
@param 	volt			Acceleration voltage (V).
@return Bimage*			Image with data from the sphere.
	
**/
Bimage*		img_fspace_extract_sphere(Bimage* p, double volt)
{
	long			i, j, xx, yy, zz, nn;
//	long			xy(p->sizeX()*p->sizeY());
	double			sx, sy, sz, s2;
	double			wl = electron_wavelength(volt);
	Vector3<double>	coor;
	
	if ( verbose )
		cout << "Extracting a sphere for a wavelength of " << wl << " A" << endl << endl;
	
	Bimage*			ps = new Bimage(Float, TSimple, p->sizeX(), p->sizeY(), 1, p->images());
	ps->sampling(p->image->sampling());
	ps->origin(p->image->origin());
	
	for ( nn=i=0; nn<p->images(); ++nn ) {
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			coor[1] = yy;
			sy = (yy - p->image->origin()[1])/p->real_size()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				coor[0] = xx;
				sx = (xx - p->image->origin()[0])/p->real_size()[0];
				s2 = sx*sx + sy*sy;
				sz = (wl/2)*s2;
//				zz = sz*p->real_size()[2] + p->image->origin()[2];
				coor[2] = sz*p->real_size()[2] + p->image->origin()[2];
				zz = coor[2] + 0.5;
				if ( zz >=0 && zz < p->sizeZ()-1 ) {
//					f = coor[2] - zz;
					j = p->index(xx, yy, zz, nn);
//					cout << xx << tab << yy << tab << zz << tab << j << endl;
					ps->set(i, (*p)[j]);
//					ps->set(i, sinc(f)*(*p)[j] + sinc(1-f)*(*p)[j+xy]);
//					ps->set(i, (1-f)*(*p)[j] + f*(*p)[j+xy]);
				}
			}
		}
	}
	
	return ps;
}
