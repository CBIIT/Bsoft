/**
@file	mg_ctf_fit.cpp
@brief	Functions for CTF (contrast transfer function) processing
@author 	Bernard Heymann
@date	Created: 19970715
@date	Modified: 20221213
**/

#include "rwimg.h"
#include "mg_ctf.h"
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

double		img_aberration_phase_fit_iter(Bimage* p, long nn, double lores, double hires, int flag, long iter, map<pair<long,long>,double>& wa);

int			complete_weights(map<pair<long,long>,double>& w, int flag)
{
	long					n(flag&1), m, t(2);
	if ( flag < 1 ) t = 1;
		
	for ( ; n<5; n+=t ) {
		for ( m=-n; m<=n; m+=2 ) {
			if ( w.find({n,m}) == w.end() )
				w[{n,m}] = 0;
		}
	}
	
	return 0;
}


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
	if ( em_ctf.check_defocus() )
		cerr << "in img_ctf_radial_average" << endl;
	
	long 			nrad = p->sizeX()/2;
	Bimage*			prad = new Bimage(Double, TSimple, nrad, 1, 1, 1);
	prad->sampling(p->sampling(0)[0], 1, 1);
	
	long			i, x, y, z;
	long 			iradius;
	double			dx, dy, dz, dr, radius, angle, fraction, cosang, sinang;
	double			yscale(p->sizeX()*1.0L/p->sizeY());
	vector<double>	num(nrad);
	double*			rdata = (double *) prad->data_pointer();

	double			smin2(1), smax2(1);
	if ( fabs(em_ctf.defocus_average()) > 1 ) {
		smin2 = 1-em_ctf.defocus_deviation()/em_ctf.defocus_average();
		smax2 = 1+em_ctf.defocus_deviation()/em_ctf.defocus_average();
	}
	
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

	long			i;
	double			s, s2, b, e, c, d, R(0), v(0);
	
	Bimage*			prad = img_ctf_radial_average(p, n, em_ctf);
	
	for ( i=rmin; i<=rmax; ++i ) {
		s = i/p->real_size()[0];
		s2 = s*s;
		b = em_ctf.calc_baseline(s);
		e = em_ctf.calc_envelope(s);
		c = em_ctf.calculate(s2, 0);
		c *= c;
		d = b + e;
		v += d*d;
		d = b + e*c - (*prad)[i];
		R += d*d;
		if ( verbose & VERB_FULL )
			cout << i << tab << (*prad)[i] << tab << b + e*c << tab << d << endl;
	}
	
	delete prad;
	
	R = sqrt(R/v);
	
	return R;
}

double		img_ctf_find_defocus2(Bimage* p, long n, CTFparam& em_ctf,
				double lores, double hires,
				double def_start, double def_end, double def_inc);

int			img_manage_power_spectrum_peak(Bimage* p, Bimage* pb, long radius)
{
	long		i, xx, yy, zz;
	double		dx, dy, f;
	
	for ( i=zz=0; zz<p->sizeZ(); ++zz ) {
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			dy = (double)yy - p->image->origin()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				dx = (double)xx - p->image->origin()[0];
				f = sqrt(dx*dx + dy*dy)/radius;
				if ( f <= 1 )
					pb->set(i, f*(*pb)[i] + (1.0-f)*(*p)[i]);
			}
		}
	}
	
	return 0;
}
/*
Bimage*		img_ctf_fit_prepare(Bimage* p, long n, double sigma)
{
	long		gauss_kernel(6*sigma);
	long		radius(p->real_size()[0]/30);
	
	Bimage*		ps = p->extract(n);
	Bimage*		pb = ps->copy();

	pb->filter_gaussian(gauss_kernel, sigma);
	img_manage_power_spectrum_peak(ps, pb, radius);

	Bimage*		pe = pb->copy();
	pe->add(p, -1, 0);
	pe->absolute();
	pe->filter_gaussian(gauss_kernel, sigma);
	pe->multiply(M_PI);
	
//	pb->add(pe, -0.5, 0);

	Bimage*		pf = ps->copy();
	pf->add(pb, -1, 0);
	pf->divide(pe, 1, 0);
	pf->truncate_to_min_max(-1,1);
//	pf->invert();
	
//	if ( verbose & VERB_DEBUG ) {
		write_img("pb.grd", pb, 0);
		write_img("pe.grd", pe, 0);
		write_img("pf.grd", pf, 0);
//	}
	
	delete ps;
	delete pb;
	delete pe;
	
	return pf;
}
*/
Bimage*		img_ctf_fit_prepare(Bimage* p, double sigma)
{
	long		gauss_kernel(6*sigma);
//	long		radius(p->real_size()[0]/30);
	long		radius(5);
	
	Bimage*		pb = p->copy();

	pb->filter_gaussian(gauss_kernel, sigma);		// Smooth to get B + E/2
	img_manage_power_spectrum_peak(p, pb, radius);	// Reconstruct central peak to counter excessive smoothing

	Bimage*		pe = pb->copy();
	pe->add(p, -1, 0);							// Subtract B + E/2 to get ~S/2
	pe->absolute();								// Convert to only positive values
	pe->filter_gaussian(gauss_kernel, sigma);	// Smooth to get ~E/pi
	pe->multiply(M_PI);							// Envelope
	
	Bimage*		pf = p->copy();
	pf->add(pb, -1, 0);							// Subtract B + E/2
	pf->divide(pe, 1, 0);						// Divide by E
	pf->truncate_to_min_max(-1,1);				// Ensure within [-1,1]
//	pf->invert();

	pb->add(pe, -0.5, 0);						// (B + E/2) - E/2
	pb->truncate_to_min_max(0,pb->maximum());	// Eliminate negative values

//	if ( verbose & VERB_DEBUG ) {
		write_img("pb.grd", pb, 0);
		write_img("pe.grd", pe, 0);
		write_img("pf.grd", pf, 0);
//	}
	
	pf->next = pb;
	pb->next = pe;
	
	return pf;
}

Bimage*		img_ctf_fit_prepare(Bimage* p, long n, double sigma)
{
	Bimage*		ps = p->extract(n);

	Bimage*		pf = img_ctf_fit_prepare(ps, sigma);
	
	delete ps;
	
	return pf;
}

/*
	CTFparam	em_ctf;
	em_ctf.volt(3e5);
	em_ctf.Cs(2e7);
	em_ctf.amp_shift(0.07);

	img_ctf_find_defocus2(pf, 0, em_ctf, 20, 5, 1e4, 5e4, 1e3);
	
	cout << "Best defocus: " << em_ctf.defocus_average() << endl;
	
	map<pair<long,long>,double>	wa;
	
	complete_weights(wa, 2);

	wa[{0,0}] = 0.07;
	wa[{2,0}] = M_PI * em_ctf.lambda() * em_ctf.defocus_average();
	wa[{4,0}] = M_PI_2 * em_ctf.Cs() * em_ctf.lambda() * em_ctf.lambda() * em_ctf.lambda();

		cout << "Aberration parameters:\nn\tm\tw" << endl;
	for ( auto w: wa )
		cout << w.first.first << tab << w.first.second << tab << w.second << endl;
	
	img_aberration_phase_fit_iter(pf, 0, 5, 4, 100, wa);

		cout << "Aberration parameters:\nn\tm\tw" << endl;
	for ( auto w: wa )
		cout << w.first.first << tab << w.first.second << tab << w.second << endl;

	Bimage*			pa = pf->copy();
	img_create_aberration(pa, wa, 2);
	pa->multiply(2);
	pa->cosine();
	pa->invert();

	write_img("ab.grd", pa, 0);
	
//	pa->subtract(pf);
	pa->replace_half(pf);

	write_img("abd.grd", pa, 0);

*/

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
double		img_ctf_fit_old(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires,
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
	
//	em_ctf.defocus_deviation(0);
//	em_ctf.astigmatism_angle(0);
	em_ctf.astigmatism(0,0);
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
		fom = img_ctf_fit_residual(p, n, em_ctf, lores, hires);
//		hires_astig -= hires_dec;
	}

	double			wri = img_water_ring_index(p, n, em_ctf);

//	img_ctf_fit_aberration(p, n, em_ctf, lores, hires);

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

	return wri;
}

double		img_ctf_fit(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires,
				double def_start, double def_end, double def_inc, int flag)
{
	if ( def_start < 1e3 || def_end < 1e4 || def_inc < 100 ) {
		cerr << "Warning: The defocus range parameters are too small:" << endl;
		cerr << tab << def_start << " A - " << def_end << " A ∆ " << def_inc << " A" << endl;
	}
	
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
		cout << "Defocus min, max, inc:          " << def_start << " - " << def_end << " ∆ " << def_inc << endl;
		cout << "Flag:                           " << flag << endl << endl;
	}
	
	
	if ( em_ctf.defocus_average() < 1 || em_ctf.defocus_average() > 2e5 )
		em_ctf.defocus_average(2e4);
	
	em_ctf.astigmatism(0,0);
	em_ctf.fom(0);
	
	double 			fom(0), sigma(10);
	Bimage*			pf = NULL;
	if ( flag & 32 ) {
		pf = img_ctf_fit_prepare(p, n, sigma);
		fom = img_ctf_find_defocus2(pf, 0, em_ctf, lores, hires, def_start, def_end, def_inc);
	} else {
		fom = img_ctf_find_defocus(p, n, em_ctf, lores, hires, def_start, def_end, def_inc);
	}
	
	if ( verbose )
		cout << "Best defocus: " << em_ctf.defocus_average() << tab << fom << endl;

	img_ctf_fit_baseline(p, n, em_ctf, lores, hires);
	
	img_ctf_fit_envelope(p, n, em_ctf, lores, hires);

	long			i(0), imax(10);
	double			pfom(0);
	
	if ( flag & 8 ) for ( i=0; i<imax && ( fabs(fom - pfom) > 1e-30 ); ++i ) {
		if ( verbose )
			cout << "Defocus [" << i << "]:                    "
					<< em_ctf.defocus_average()
					<< " +- " << em_ctf.defocus_deviation()
					<< " @ " << em_ctf.astigmatism_angle()*180/M_PI
					<< " (" << fom << ")" << endl;
		pfom = fom;
		def_start = em_ctf.defocus_average()/10;
		def_end = em_ctf.defocus_average()*10;
		def_inc = def_start;
//		cout << def_start << " - " << def_end << " ∆ " << def_inc << endl;
		fom = img_ctf_fit_astigmatism(p, n, em_ctf, lores, hires);
		img_ctf_find_defocus(p, n, em_ctf, lores, hires, def_start, def_end, def_inc);
		img_ctf_fit_baseline(p, n, em_ctf, lores, hires);
		img_ctf_fit_envelope(p, n, em_ctf, lores, hires);
		fom = img_ctf_fit_residual(p, n, em_ctf, lores, hires);
	} else if ( flag & 32 ) {
		map<pair<long,long>,double>	wa;
	
		complete_weights(wa, 2);

		wa[{0,0}] = em_ctf.aberration_weight(0,0);
		wa[{2,0}] = em_ctf.aberration_weight(2,0);
		wa[{4,0}] = em_ctf.aberration_weight(4,0);

		img_aberration_phase_fit_iter(pf, 0, lores, hires, 4, 1000, wa);
			
		em_ctf.aberration_weights(wa);
		
		delete pf;
	}
	
	if ( verbose )
		em_ctf.show_aberration();

//	if ( verbose ) {
//		cout << "Aberration parameters:\nn\tm\tw" << endl;
//		for ( auto w: wa )
//			cout << w.first.first << tab << w.first.second << tab << w.second << endl;
//	}

/*
	Bimage*			pa = pf->copy();
	img_create_aberration(pa, wa, 2);
	pa->multiply(2);
	pa->cosine();
	pa->invert();

	write_img("ab.grd", pa, 0);
	
//	pa->subtract(pf);
	pa->replace_half(pf);

	write_img("abd.grd", pa, 0);
*/
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

	return wri;
}


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
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_ctf_fit_baseline" << endl;

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
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_ctf_fit_envelope" << endl;

	Bimage* 	prad = img_ctf_radial_average(p, n, em_ctf);
	
	if ( !prad ) {
		cerr << "Error in img_ctf_fit_envelope: No radial power spectrum calculated!" << endl;
		return -1;
	}
	
	double		R = ctf_fit_envelope(prad, p->real_size()[0], em_ctf, lores, hires);
	
	delete prad;
	
	if ( verbose & VERB_PROCESS )
		em_ctf.show_envelope();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_ctf_fit_envelope: done" << endl;

	return R;
}

double		ctf_test_defocus(CTFparam& em_ctf, double def,
				double step_size, vector<double>& r, long rmin, long rmax)
{
//	long			i, nr(rmax-rmin+1);
	double			fom(0);
	
	em_ctf.defocus_average(def);
	
	vector<double>	c =	em_ctf.calculate(r.size(), 1, step_size);

	if ( !c.size() ) {
		cerr << "Error in ctf_test_defocus: CTF not calculated!" << endl;
		return fom;
	}
	
//	for ( i=rmin; i<=rmax; i++ ) fom += (c[i]*c[i] - 0.5)*r[i];
//	fom /= nr;

	double			c2, c2sum(0), r2sum(0);
	for ( long i=rmin; i<=rmax; ++i ) {
		c2 = c[i]*c[i] - 0.5;
		fom += c2*r[i];
		c2sum += c2*c2;
		r2sum += r[i]*r[i];
	}
	fom /= sqrt(c2sum*r2sum);

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
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Defocus range:                  " << def_start*1e-4 << " - " << def_end*1e-4 << " um" << endl;
		cout << "Defocus increment:              " << def_inc*1e-4 << " um" << endl;
		cout << "Frequency pixel range:          " << rmin << " - " << rmax << endl;
		cout << "Frequency step size:            " << step_size << " 1/A" << endl;
	}
	
	double			ds(def_start), de(def_end);
	
//	if ( verbose & VERB_FULL )
		cout << "Defocus\tCC" << endl;
	while ( def_inc >= 50 ) {
		for ( def=ds; def<=de; def+=def_inc ) {
			fom = ctf_test_defocus(em_ctf, def, step_size, r, rmin, rmax);
			if ( fom_max < fom ) {
				fom_max = fom;
				def_best = def;
			}
//			if ( verbose & VERB_FULL )
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
	
	return em_ctf.fom();
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
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_ctf_find_defocus" << endl;
		
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
	if ( rmax >= nx ) rmax = nx-1;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_ctf_find_defocus: lores=" << lores << tab << "hires=" << hires << endl;
		cout << "DEBUG img_ctf_find_defocus: nx=" << nx << tab << "rmin=" << rmin << tab << "rmax=" << rmax << endl;
	}
	
	for ( i=0; i<rmin; i++ ) prad->set(i, (*prad)[rmin]);

	vector<double>	r = moving_polynomial(2, nx, (double *) prad->data_pointer(), nx/5);

	for ( i=0; i<nx; i++ )
		if ( r[i] ) r[i] = (*prad)[i]/r[i] - 1;
		else r[i] = 0;

//	for ( i=0; i<nx; i++ )
//		cout << i << tab << r[i] << endl;
		
	double			def, def_best(em_ctf.defocus_average());
	double			fom, fom_max(-2);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Finding the defocus (old):" << endl;
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
//			fom = ctf_test_defocus(em_ctf, def, nx, step_size, r, rmin, rmax);
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
	
	delete prad;
	
	em_ctf.defocus_average(def_best);
	if ( fom_max > -1 ) em_ctf.fom(fom_max);
	
	return em_ctf.fom();
}

double		img_ctf_find_defocus2(Bimage* p, long n, CTFparam& em_ctf,
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

	vector<double>	r(nx,0);

	for ( i=0; i<nx; i++ )
		r[i] = (*prad)[i];

	double			def, def_best(em_ctf.defocus_average());
	double			fom, fom_max(-2);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Finding the defocus (new):" << endl;
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
//			fom = ctf_test_defocus(em_ctf, def, nx, step_size, r, rmin, rmax);
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
//		em_ctf.defocus_deviation(def_dev);
		for ( angle=angle_min, j=0; angle<=angle_max; angle+=angle_inc, j++ ) {
			ec[j] = em_ctf;
//			ec[j].astigmatism_angle(angle);
			ec[j].astigmatism(def_dev, angle);
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
	
	while ( best_ang < -M_PI_2 ) best_ang += M_PI;
	while ( best_ang > M_PI_2 ) best_ang -= M_PI;
	
	em_ctf.astigmatism(best_dev, best_ang);
	
	if ( verbose & VERB_PROCESS )
		cout << "Astigmatism refined:            " << em_ctf.defocus_average() << " +- " 
			<< em_ctf.defocus_deviation() << " @ " << em_ctf.astigmatism_angle()*180.0/M_PI << " (OM = " << best_fom << ")" << endl;

	return best_fom;
}

/**
@brief 	Fits a phase image with aberration functions.
@param 	*p			image structure.
@param	nn			sub-image number.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit.
@param	flag		0=full, 1=odd, 2=even.
@return map<pair<long,long>,double>	aberration weights.

	The aberration terms are calculated for each pixel and added to a symmetric matrix.
	The matrix is inverted and the aberration weights calculated.
	The origin of the image must be at {0,0}.

**/
map<pair<long,long>,double>	img_aberration_phase_fit(Bimage* p, long nn, double lores, double hires, int flag)
{
	if ( hires < p->image->sampling()[0] ) hires = p->image->sampling()[0];
	if ( lores && lores < hires ) lores = 1e5;
	
	if ( verbose & VERB_FULL )
		cout << "Image " << nn << ":" << endl;

	double			shi(1/hires);
	double			slo = (lores)? 1/lores: 0;
	long			nt(15);
	if ( flag == 1 ) nt = 6;
	else if ( flag == 2 ) nt = 9;
	
	long			i, j, k, xx, yy;
	double			slen;
	Vector3<double>	s;
	Matrix			A(nt,nt);
	vector<double>	b(nt,0);
	vector<double>	t(nt,0);
	
	for ( i=nn*p->image_size(), yy=0; yy<p->sizeY(); ++yy ) {
		s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
		if ( s[1] >= 0.5 ) s[1] -= 1;
		s[1] /= p->image->sampling()[1];
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s[0] >= 0.5 ) s[0] -= 1;
			s[0] /= p->image->sampling()[0];
			slen = s.length();
			if ( slen >= slo && slen <= shi ) {
				if ( flag == 1 ) t = aberration_odd_terms(nt,s[0],s[1]);
				else if ( flag == 2 ) t = aberration_even_terms(nt,s[0],s[1]);
				else t = aberration_terms(nt,s[0],s[1]);
				for ( j=0; j<nt; ++j ) {
					b[j] += t[j]*(*p)[i];
					for ( k=0; k<nt; ++k )
						A[j][k] += t[j]*t[k];
				}
			}
		}
	}

	A.singular_value_decomposition(b);

	if ( verbose & VERB_FULL ) {
		string		ws = "[" + concatenate(b) + "]";
		cout << "Image " << nn << ": " << ws << endl;
	}
	
	return	aberration_weights(b, flag);
}

double		aberration_fit_R(Bsimplex& simp)
{
	long			i, j, k, nt(simp.parameters());
	double			v, R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		for ( v=0, k=0; k<nt; ++k, ++j )
			v += simp.parameter(k) * x[j];
		df = angle_set_negPI_to_PI(f[i] - v);
		R += df*df;
	}
	
	R /= simp.points();
	
	return R;
}

double		sine_aberration_fit_R(Bsimplex& simp)
{
	long			i, j, k, nt(simp.parameters());
	double			v, R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		for ( v=0, k=0; k<nt; ++k, ++j )
			v += simp.parameter(k) * x[j];
		v = sin(v);
		df = f[i] - v*v;
		R += df*df;
	}
	
	R /= simp.points();
	
	return R;
}

double		cosine_aberration_fit_R(Bsimplex& simp)
{
	long			i, j, k, nt(simp.parameters());
	double			v, cc(0), R(0), w4(0);
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	if ( simp.constants() ) w4 = simp.constant(0);
	
	for ( i=j=0; i<simp.points(); ++i ) {
		if ( w4 ) v = w4*x[j+2]*x[j+2];
		else v = 0;
		for ( k=0; k<nt; ++k, ++j )
			v += simp.parameter(k) * x[j];
		v = -cos(2*v);
		cc += f[i]*v;
//		v2 += v*v;
//		f2 += f[i]*f[i];
	}
	
//	cc /= sqrt(v2*f2);
//	cc /= sqrt(v2);
	cc /= simp.points();
	
	R = 1 - cc;
	
	return R;
}

/*
	flag:
		0 = full/all, nt = 15
		1 = odd, nt = 6
		2 = even, nt = 9
		4 = even, nt = 9
		5 = even, nt = 4
*/
double		img_aberration_phase_fit_iter_old(Bimage* p, long nn, double lores, double hires, int flag, long iter, map<pair<long,long>,double>& wa)
{
	if ( hires < p->image->sampling()[0] ) hires = p->image->sampling()[0];
	if ( lores && lores < hires ) lores = 1e5;
	
	if ( verbose & VERB_FULL )
		cout << "Image " << nn << ":" << endl;

	double			shi(1/hires);
	double			slo = (lores)? 1/lores: 0;
	long			nt(15);
	if ( flag == 1 ) nt = 6;
	else if ( flag == 2 || flag == 4 ) nt = 9;
	
	long			i, j, xx, yy;
	double			slen;
	Vector3<double>	s;
	vector<double>	t, terms, phases;
		
	for ( i=nn*p->image_size(), yy=0; yy<p->sizeY(); ++yy ) {
		s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
		if ( s[1] >= 0.5 ) s[1] -= 1;
		s[1] /= p->image->sampling()[1];
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s[0] >= 0.5 ) s[0] -= 1;
			s[0] /= p->image->sampling()[0];
			slen = s.length();
			if ( slen >= slo && slen <= shi ) {
				phases.push_back((*p)[i]);
				if ( flag == 1 ) t = aberration_odd_terms(nt,s[0],s[1]);
				else if ( flag == 2 || flag == 4 ) t = aberration_even_terms(nt,s[0],s[1]);
				else t = aberration_terms(nt,s[0],s[1]);
				for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
			}
		}
	}

	Bsimplex			simp(nt, nt, 0, phases.size(), terms, phases);
	
	j = 0;
	for ( auto w: wa ) {
		simp.parameter(j, w.second);
		if ( w.second ) simp.limits(j++, w.second*0.8, w.second*1.2);
		else simp.limits(j++, -20, 20);
	}
	
	if ( flag < 2 ) {		// Odd
	}
	
	if ( flag%2 == 0 ) {	// Even
		simp.limits(0, wa[{0,0}]-0.1, wa[{0,0}]+0.1);	// amp
		simp.limits(1, -wa[{2,0}]*0.5, wa[{2,0}]*0.5);	// ∆∆f
		simp.limits(2, wa[{2,0}]*0.8, wa[{2,0}]*1.2);	// ∆f
		simp.limits(3, -wa[{2,0}]*0.5, wa[{2,0}]*0.5);	// ∆∆f
		simp.limits(4, -wa[{4,0}]*0.5, wa[{4,0}]*0.5);
		simp.limits(5, -wa[{4,0}]*0.5, wa[{4,0}]*0.5);
		simp.limits(6, wa[{4,0}]*0.8, wa[{4,0}]*1.2);	// Cs
		simp.limits(7, -wa[{4,0}]*0.5, wa[{4,0}]*0.5);
		simp.limits(8, -wa[{4,0}]*0.5, wa[{4,0}]*0.5);
	}
	
	double			R(0), iR(0);
	
	if ( flag < 4 ) {
		iR = simp.R(aberration_fit_R);
		iR /= simp.dependent_variance();
		R = simp.run(iter, 1e-5, aberration_fit_R, 0);
	} else {
//		R = simp.run(iter, 1e-5, sine_aberration_fit_R, 0);
		iR = simp.R(cosine_aberration_fit_R);
		iR /= simp.dependent_variance();
		R = simp.run(iter, 1e-5, cosine_aberration_fit_R, 50);
	}
	
	R /= simp.dependent_variance();
	p->image[nn].FOM(R);

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << "Image " << nn << ":" << endl;
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}
	
	map<pair<long,long>,double>	w = aberration_weights(simp.parameter_vector(), flag);

	for ( auto w1: w ) wa[w1.first] = w1.second;

	return R;
}

double		img_aberration_phase_fit_iter(Bimage* p, long nn, double lores, double hires, int flag, long iter, map<pair<long,long>,double>& wa)
{
	if ( hires < p->image->sampling()[0] ) hires = p->image->sampling()[0];
	if ( lores && lores < hires ) lores = 1e5;
	
	if ( verbose & VERB_FULL )
		cout << "Image " << nn << ":" << endl;

	double			shi(1/hires);
	double			slo = (lores)? 1/lores: 0;
	long			nt(wa.size()), nc(0);
	if ( flag == 5 ) {
		nt = 4;
		nc = 1;
	}
	
	long			i, j, xx, yy;
	double			slen, dwin(0.5);
	Vector3<double>	s;
	vector<double>	t, terms, phases;
		
	for ( i=nn*p->image_size(), yy=0; yy<p->sizeY(); ++yy ) {
		s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
		if ( s[1] >= 0.5 ) s[1] -= 1;
		s[1] /= p->image->sampling()[1];
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s[0] >= 0.5 ) s[0] -= 1;
			s[0] /= p->image->sampling()[0];
			slen = s.length();
			if ( slen >= slo && slen <= shi ) {
				phases.push_back((*p)[i]);
				t = aberration_terms(wa, s[0], s[1]);
				for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
			}
		}
	}

	Bsimplex			simp(nt, nt, nc, phases.size(), terms, phases);
	if ( nc ) simp.constant(0, wa[{4,0}]);
	
	if ( verbose )
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
	j = 0;
	for ( auto w: wa ) {
		if ( j >= nt ) break;	// flag = 5
		simp.parameter(j, w.second);
		if ( w.second ) {
			if ( w.first.second == 0 ) {
				simp.limits(j++, w.second*(1-dwin), w.second*(1+dwin));
			} else {
				simp.limits(j++, w.second - 5*dwin*fabs(w.second), w.second + 5*dwin*fabs(w.second));
			}
//			simp.limits(j++, -dwin*fabs(w.second), +dwin*fabs(w.second));
		} else {
			simp.limits(j++, -20, 20);
		}
		if ( verbose )
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j-1) << tab << simp.limit_low(j-1) << tab << simp.limit_high(j-1) << endl;
	}

	double			R(0), iR(0);
	
	if ( flag < 4 ) {
		iR = simp.R(aberration_fit_R);
		iR /= simp.dependent_variance();
		R = simp.run(iter, 1e-5, aberration_fit_R, 0);
	} else {
//		R = simp.run(iter, 1e-5, sine_aberration_fit_R, 0);
		iR = simp.R(cosine_aberration_fit_R);
		iR /= simp.dependent_variance();
		R = simp.run(iter, 1e-5, cosine_aberration_fit_R, 50);
	}
	
	R /= simp.dependent_variance();
	p->image[nn].FOM(R);

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << "Image " << nn << ":" << endl;
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}
	
	j = 0;
	for ( auto& w: wa ) {
		if ( j >= nt ) break;	// flag = 5
		w.second = simp.parameter(j++);
	}

	return R;
}

vector<map<pair<long,long>,double>>	img_aberration_phase_fit(Bimage* p, double lores, double hires, int flag,
		long iter)
{
	vector<map<pair<long,long>,double>>		win;
	return img_aberration_phase_fit(p, lores, hires, win, flag, iter);
}

/**
@brief 	Fits a phase image with aberration terms.
@param 	*p			phase image.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit.
@param	win			initial aberration parameters.
@param 	flag		0=full, 1=odd, 2=even.
@param	iter		maximum number of iterations.
@return vector<map<pair<long,long>,double>>	vector of aberration weights.

	

**/
vector<map<pair<long,long>,double>>	img_aberration_phase_fit(Bimage* p, double lores, double hires,
			vector<map<pair<long,long>,double>>& win, int flag, long iter)
{
	if ( hires < p->image->sampling()[0] ) hires = p->image->sampling()[0];
	if ( lores && lores < hires ) lores = 1e5;

	if ( verbose ) {
		cout << "Fitting a phase difference image:" << endl;
		if ( flag == 1 ) cout << "Odd terms" << endl;
		else if ( flag == 2 ) cout << "Even terms" << endl;
		else cout << "All terms" << endl;
		if ( win.size() ) {
			cout << "Input terms:\nn\tm\tw" << endl;
			for ( auto w: win[0] ) cout << w.first.first << tab << w.first.second << tab << w.second << endl;
		}
		cout << "Resolution limits:              " << hires << " - ";
		if ( lores ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Iterations:                     " << iter << endl;
	}

#ifdef HAVE_GCD
	__block	long	ndone(0);
	__block	vector<map<pair<long,long>,double>> wa(p->images());
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(p->images(), dispatch_get_global_queue(0, 0), ^(size_t nn){
		map<pair<long,long>,double>	w;
		if ( win.size() ) {
			if ( nn < win.size() ) w = win[nn];
			else w = win[0];
			complete_weights(w, flag);
		} else {
			w = img_aberration_phase_fit(p, nn, lores, hires, flag);
		}
		if ( iter > 0 ) img_aberration_phase_fit_iter(p, nn, lores, hires, flag, iter, w);
		wa[nn] = w;
		dispatch_sync(myq, ^{
			ndone++;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/p->images() << " %    \r" << flush;
		});
	});
#else
	long			ndone(0);
	vector<map<pair<long,long>,double>> wa(p->images());
#pragma omp parallel for
	for ( long nn=0; nn<p->images(); ++nn ) {
		map<pair<long,long>,double>	w;
		if ( win.size() ) {
			if ( nn < win.size() ) w = win[nn];
			else w = win[0];
			complete_weights(w, flag);
		} else {
			w = img_aberration_phase_fit(p, nn, lores, hires, flag);
		}
		if ( iter > 0 ) img_aberration_phase_fit_iter(p, nn, lores, hires, flag, iter, w);
		wa[nn] = w;
	#pragma omp critical
		{
			ndone++;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/p->images() << " %    \r" << flush;
		}
	}
#endif

	if ( verbose ) {
		string		ws;
		cout << endl << "Image\tR\tWeights" << endl;
		for ( long nn=0; nn<p->images(); ++nn ) {
			ws = aberration_weight_string(wa[nn]);
			cout << nn+1 << tab << p->image[nn].FOM() << tab << ws << endl;
		}
	}

	return wa;
}

double		img_ctf_fit_aberration(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	long		iter(1000);
	map<pair<long,long>,double>& wa(em_ctf.aberration_weights());
	
//	wa[{2,0}] = M_PI*em_ctf.lambda()*em_ctf.defocus_average();
//	wa[{4,0}] = M_PI_2*em_ctf.lambda()*em_ctf.lambda()*em_ctf.lambda()*em_ctf.Cs();
	
	double		R = img_aberration_phase_fit_iter(p, n, lores, hires, 4, iter, wa);
	
	em_ctf.update_aberration_weights(wa);
//	em_ctf.update_CTF_from_aberration_weights();
	
	em_ctf.show();
	
	return R;
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
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_water_ring_index" << endl;

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

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_water_ring_index: b1=" << b1 << " b2=" << b2 << " wp=" << wp << endl;

	n1 += n2;
	if ( n1 ) b1 = (b1 + b2)/n1;
	if ( nw ) wp /= nw;
	
	if ( b1 ) wri = wp/b1 - 1;

	return wri;
}

double		water_ring_R(Bsimplex& simp)
{
	long			i, j;
	double			R(0), df, s;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); i++ ) {
		s = x[j++];
		s -= simp.parameter(2) + simp.parameter(3)*x[j++];
		s /= simp.parameter(4);
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*exp(-0.5*s*s));
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}


vector<double>	img_fit_water_ring(Bimage* p)
{
	long			i, xx, yy;
	double			sx, sy, s;
	double			smin(0.2), smax(0.4);
	vector<double>	x, fx;
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		sy = (p->image->origin()[1] - yy)/p->real_size()[1];
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			sx = (p->image->origin()[0] - xx)/p->real_size()[0];
			s = sqrt(sx*sx + sy*sy);		// Spatial frequency
			if ( s > smin && s < smax ) {
				x.push_back(s);
				x.push_back(cos(2*atan2(sy, sx)));	// cos(2*angle)
				fx.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 5, 0, fx.size(), x, fx);

	simp.parameter(0, fx.back());	// Constant baseline
	simp.parameter(1, 1);			// Water ring power
	simp.parameter(2, 0.27);		// Spatial frequency average
	simp.parameter(3, 0.01);		// Spatial frequency deviation
	simp.parameter(4, 0.03);		// Spatial frequency sigma
	simp.limits(0, 0, 10*fx.back());
	simp.limits(1, 0, 100*fx.back());
	simp.limits(2, 0.25, 0.29);
	simp.limits(3, -0.1, 0.1);
	simp.limits(4, 0.01, 0.05);

	double			R = simp.run(10000, 0.01, water_ring_R);
	R /= sqrt(simp.dependent_variance());
	
	vector<double>	wrfit(6,0);
	for ( i=0; i<5; i++ ) wrfit[i] = simp.parameter(i);
	wrfit[i] = R;
	
	if ( verbose ) {
		cout << "Water ring fit:" << endl;
		cout << "Limits:                          " << smin << " - " << smax << " /A" << endl;
		cout << wrfit[0] << " + " << wrfit[1] << " * exp(-0.5*((s-" <<
			wrfit[2] << "-" << wrfit[3] << "*cos(2a))/" << wrfit[4] << ")^2)" << endl;
		cout << "R = " << wrfit[5] << endl << endl;
	}

	return wrfit;
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

double		gaussian_WR(Bsimplex& simp)
{
	long			i;
	double			R(0), df, s2;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		s2 = x[i]*x[i];
		df = 1 - (simp.parameter(0) + simp.parameter(1)*exp(simp.parameter(2)*s2))/f[i];
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

double		double_gaussian_WR(Bsimplex& simp)
{
	long			i;
	double			R(0), df, s2;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		s2 = x[i]*x[i];
		df = 1 - (simp.parameter(0) + simp.parameter(1)*exp(simp.parameter(2)*s2) +
			simp.parameter(3)*exp(simp.parameter(4)*s2))/f[i];
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
@param 	*x		x array (at least order+1 values).
@param	*y		y array (at least order+1 values).
@param	coeff	array in which 5 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Double gaussian function:
		y = a + b*exp(c*s^2) + d*exp(f*s^2)
**/
double 		ctf_fit_double_gaussian_baseline(vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;
	double				ymin(1e30);
	
	for ( i=0; i<y.size(); ++i ) if ( ymin > y[i] ) ymin = y[i];
	
//	cout << "ymin=" << ymin << endl;

	if ( coeff[0] < ymin/2 ) coeff[0] = ymin;
	if ( coeff[0] > 2*ymin ) coeff[0] = ymin;
	if ( coeff[1] < y[0] ) coeff[1] = y[0];
	if ( coeff[2] > -0.01 ) coeff[2] = -10;
	if ( coeff[3] < ymin ) coeff[3] = ymin;
	if ( coeff[4] > -0.01 ) coeff[4] = -100;

	Bsimplex			simp(1, 5, 0, y.size(), x, y);
	
	simp.parameters(5, coeff);
	simp.limits(0, ymin/2, 2*ymin);
	simp.limits(1, coeff[1]/5, 5*coeff[1]);
	simp.limits(2, -1e3, -0.1);
	
	if ( y.size() > 5 ) {
		simp.limits(3, 0.01, 2*coeff[3]);
		simp.limits(4, -1e4, -0.01);
	}

	double			R = simp.run(10000, 0.01, double_gaussian_R);
	
	for ( i=0; i<5; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

/*
@brief Fits an EMAN-style baseline function.
@param *x		x array (at least order+1 values).
@param *y		y array (at least order+1 values).
@param coeff	array in which 4 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	EMAN baseline function:
		y = a + b*exp( c*sqrt(s) + d*s^2 )
**/
double 		ctf_fit_eman_to_baseline(vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;

	Bsimplex			simp(1, 4, 0, y.size(), x, y);
	
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
@param *x		x array (at least order+1 values).
@param *y		y array (at least order+1 values).
@param coeff	array in which 3 coefficients are returned
	                (if NULL, no coefficients returned).
@return double 	R value.

	Gaussian function:
		y = a + b*exp(c*s^2)
**/
double 		ctf_fit_gaussian(vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;
	
	if ( coeff[0] ) if ( coeff[0] > y.back() ) coeff[0] = y.back();
	if ( coeff[1] < y[0] ) coeff[1] = y[0];
	if ( coeff[2] > -0.01 ) coeff[2] = -10;
	
	Bsimplex			simp(1, 3, 0, y.size(), x, y);
	
	simp.parameters(3, coeff);
	if ( coeff[0] ) simp.limits(0, 0, 2*coeff[0]);
	simp.limits(1, coeff[1]/2, 5*coeff[1]);
	simp.limits(2, -1e3, -0.1);
	
	double			R = simp.run(10000, 0.01, gaussian_WR);
	
	for ( i=0; i<3; i++ ) coeff[i] = simp.parameter(i);
	
	return R;
}

double		ctf_fit_baseline_bump(vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i, ilo(0), ns, n(y.size());
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
	simp.limits(4, -5000, -200);
	
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
double 		ctf_fit_double_gaussian_envelope(vector<double>& x, vector<double>& y, vector<double>& coeff)
{
	long				i;

	if ( coeff[0] ) if ( coeff[0] > y.back() ) coeff[0] = y.back();
	if ( coeff[1] < y[0] ) coeff[1] = y[0];
	if ( coeff[2] > -0.01 ) coeff[2] = -10;
	if ( coeff[3] < y.back() ) coeff[3] = y.back();
	if ( coeff[4] > -0.01 ) coeff[4] = -10;

	Bsimplex			simp(1, 5, 0, y.size(), x, y);
	
	simp.parameters(5, coeff);
	if ( coeff[0] ) simp.limits(0, 0, 2*coeff[0]);
	simp.limits(1, coeff[1]/2, 5*coeff[1]);
	simp.limits(2, -1e3, -0.1);
	
	if ( y.size() > 5 ) {
		simp.limits(3, 0.01, 2*coeff[3]);
		simp.limits(4, -1e4, -0.01);
	}
	
	double			R = simp.run(10000, 0.01, double_gaussian_WR);
	
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
	
	vector<double>	s, d;
	
	vector<double>	zero = em_ctf.zeroes(max_s);
	long			n(zero.size());
	
	if ( n > prad->sizeX() ) n = prad->sizeX();

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_baseline: zeroes=" << n << endl;
		for ( auto z: zero ) cout << z << endl;
	}
	
	double			min(1e10), dmin;
	
	if ( n > 3 ) {		// Enough zeroes to fit the baseline
		for ( i=0; i<n; i++ ) {
			ib = (long) (zero[i]*real_size + 0.5);
			dmin = (*prad)[ib];
			for ( j=(ib>k)?ib-k:0; j<=ib+k && j<prad->sizeX(); j++ )	// minimum within kernel
				if ( dmin > (*prad)[j] ) dmin = (*prad)[j];
			s.push_back(zero[i]);
			d.push_back(dmin);
		}
	} else {			// Too few zeroes, use local minima
		vector<double>	mavg = moving_average(prad->sizeX(), (double *) prad->data_pointer(), prad->sizeX()/10);
		for ( i=istart, n=0; i<iend; i++ ) {
			if ( min > (*prad)[i] - mavg[i] ) {
				min = (*prad)[i] - mavg[i];
				dmin = (*prad)[i];
				if ( dmin > mavg[i] ) dmin = mavg[i];
				s.push_back(recip_interval*i);
				d.push_back(dmin);
			}
			if ( i%w == 0 ) {
				min = 1e10;
				n++;
			}
		}
		s.push_back(10);
		d.push_back(0);
		n++;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_baseline: type " << em_ctf.baseline_type() << endl;
		cout << "Coefficients:" << endl;
		for ( i=0; i<5; i++ ) cout << tab << em_ctf.baseline(i);
		cout << endl << "Zeroes:" << endl;
		for ( i=0; i<s.size(); ++i ) cout << s[i] << tab << d[i] << endl;
	}

	if ( em_ctf.baseline_type() > 3 ) {
		ctf_fit_baseline_bump(s, d, coeff);
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
			R = ctf_fit_double_gaussian_baseline(s, d, coeff);
			break;
		case 3:
		case 6:
			R = ctf_fit_eman_to_baseline(s, d, coeff);
			break;
		default:
			cerr << "Error: Baseline type " << em_ctf.baseline_type() << " not supported" << endl;
	}

	if ( em_ctf.baseline_type() > 0 && em_ctf.baseline_type() < 7 )
		em_ctf.baseline(coeff);

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_baseline: type " << em_ctf.baseline_type() << endl;
		cout << "Coefficients fitted:" << endl;
		for ( i=0; i<5; i++ ) cout << tab << em_ctf.baseline(i);
		cout << endl << tab << "R=" << R << endl;
	}

	return R;
}

double		ctf_fit_envelope(Bimage* prad, double real_size, CTFparam& em_ctf, double lores, double hires)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_fit_envelope" << endl;

	double			pixel_size = prad->sampling(0)[0];
	if ( pixel_size <= 0 ) pixel_size = 1;
	
	double			smin = (lores)? 1/lores: 0;
	double			smax = (hires)? 1/hires: 1.0/(2*pixel_size);
	
	long			i, j, ib, ie;
	double			max_s(0.45/pixel_size);
	double			R(1e30);
	vector<double>	zero = em_ctf.zeroes(max_s);
	long			nzero(zero.size());

	if ( nzero > prad->sizeX() ) nzero = prad->sizeX();
	
	if ( nzero < 5 ) {
		cerr << "Error in ctf_fit_envelope: Cannot fit with fewer than 5 zeroes!" << endl;
		cerr << tab << "Defocus = " << em_ctf.defocus_average() << endl;
		return 0;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG ctf_fit_envelope: nzero=" << nzero << endl;
	
	double			b, si, v;
	vector<double>	coeff(NCTFPARAM, 0);
	vector<double>	s;
	vector<double>	d;

	for ( i=1; i<nzero; i++ ) {
		si = (zero[i-1]+zero[i])*0.5;
		if ( si >= smin && si <= smax ) {
			ib = (long) (zero[i-1]*real_size + 0.5);
			ie = (long) (zero[i]*real_size + 0.5);
			b = em_ctf.calc_baseline(si);
			v = (*prad)[ib] - b;
			for ( j=ib; j<=ie && j<prad->sizeX(); j++ )	// maximum within kernel
				if ( v < (*prad)[j] - b ) v = (*prad)[j] - b;
			s.push_back(si);
			d.push_back(v);
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_envelope: type " << em_ctf.envelope_type() << endl;
		cout << "Coefficients:" << endl;
		for ( i=0; i<5; i++ ) cout << tab << em_ctf.envelope(i);
		cout << endl << "Zeroes:" << endl;
		for ( i=0; i<s.size(); ++i ) cout << s[i] << tab << d[i] << endl;
	}

	switch ( em_ctf.envelope_type() ) {
		case 1:		// Single gaussian
			for ( i=0; i<2; i++ ) coeff[i+1] = em_ctf.envelope(i);
			R = ctf_fit_gaussian(s, d, coeff);
			for ( i=0; i<2; i++ ) coeff[i] = coeff[i+1];
			coeff[i] = 0;
			break;
		case 2:		// Single gaussian with constant
			for ( i=0; i<3; i++ ) coeff[i] = em_ctf.envelope(i);
			R = ctf_fit_gaussian(s, d, coeff);
			break;
		case 3:		// Double gaussian
			for ( i=0; i<4; i++ ) coeff[i+1] = em_ctf.envelope(i);
			R = ctf_fit_double_gaussian_envelope(s, d, coeff);
			for ( i=0; i<4; i++ ) coeff[i] = coeff[i+1];
			coeff[i] = 0;
			break;
		case 4:		// Double gaussian with constant
			for ( i=0; i<5; i++ ) coeff[i] = em_ctf.envelope(i);
			R = ctf_fit_double_gaussian_envelope(s, d, coeff);
			break;
		default: break;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG ctf_fit_envelope:";
		for ( i=0; i<5; i++ ) cout << tab << coeff[i];
		cout << endl << "R = " << R << endl;
	}
	
	em_ctf.envelope(coeff);
	
	return R;
}

