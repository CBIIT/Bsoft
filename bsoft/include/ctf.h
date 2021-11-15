/**
@file	ctf.h 
@brief	Header file for CTF (contrast transfer function) functions
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210531
**/

#include <cmath>
#include <vector>

#include "Bstring.h"
#include "Complex.h"
#include "json.h"
#include "utilities.h"

#define	NCTFPARAM	10

#ifndef _CTFparam_
/************************************************************************
@Object: class CTFparam
@Description:
	Contrast transfer function parameter class.
@Features:
	This defines the contrast transfer parameters for a microscope at
	a particular combination of settings.
*************************************************************************/
class CTFparam {
private:
	double			av;				// Acceleration voltage (V)
	double			fl;				// Focal length (angstrom)
	double			cs;				// Spherical aberration coefficient (angstrom)
	double			cc;				// Chromatic aberation (angstrom)
	double			a;				// Illumination half-angle (radians)
	double			de;				// Energy spread (eV)
	double			oa;				// Objective aperture radius (angstrom)
	double			sw;				// Energy filter slit width (eV), 0 if not used
	double			ac;				// Amplitude contrast phase shift (radian)
	double			def_avg;		// Defocus average (angstrom)
	double			def_dev;		// Defocus deviation (angstrom)
	double			ast_ang;		// Astigmatism angle (radian)
	double			firstzero;		// CTF first zero (angstrom)
	long			bt;				// Baseline type (1=poly, 2=double_gauss, 3=eman)
	vector<double>	base;			// Baseline (up to 10 coefficients)
	long			et;				// Envelope type (1=single_gauss, 2=c+single_gauss, 3=double_gauss, 4=c+double_gauss)
	vector<double>	env;			// Envelope (up to 10 coefficients)
	double			f;				// CTF fit figure-of-merit
	double			wl;				// Wave length
	long double		t1, t2;			// Terms in CTF
	void			initialize() {
		a = de = sw = def_avg = def_dev = ast_ang = firstzero = f = 0;
		av = 300000;
		fl = 3.5e7;
		oa = 2e6;
		cs = 2e7;
		cc = 2e7;
		ac = 0.07;
		for ( int i=0; i<NCTFPARAM; i++ ) {
			base.push_back(0);
			env.push_back(0);
		}
		bt = 1;
		base[0] = 1;
		et = 4;
		env[0] = 0.2;
		env[1] = -1000;
		wl = t1 = t2 = 0;
	}
	void			update_terms() {
		wl = lambda();
		t1 = M_PI_2*wl*wl*wl*cs;
		t2 = M_PI*wl;
	}
public:
	CTFparam() { initialize(); update_terms(); }
	CTFparam(double v, double sa, double ac) {
		initialize();
		volt(v); Cs(sa); amp_shift(ac);
	}
	int		update(CTFparam* ctf) {
		if ( !ctf ) return -1;
		return update(*ctf);
	}
	int		update(CTFparam& ctf);
	
	double	volt() { return av; }
	void	volt(double v) { av = v; update_terms(); }
	double	focal_length() { return fl; }
	void	focal_length(double v) { fl = v; }
	double	Cs() { return cs; }
	void	Cs(double v) { cs = v; update_terms(); }
	double	Cc() { return cc; }
	void	Cc(double v) { cc = v; }
	double	alpha() { return a; }
	void	alpha(double v) { a = v; }
	double	dE() { return de; }
	void	dE(double v) { de = v; }
	double	amp_shift() { return ac; }
	void	amp_shift(double v) { ac = v; }
//	double	phi_fac() { if ( pf <= 0 ) pf = sqrt(1-af*af); return pf; }
	double	objective_aperture() { return oa; }
	void	objective_aperture(double v) { oa = v; }
	double	slit_width() { return sw; }
	void	slit_width(double v) { sw = v; }
	double	defocus_average() { return def_avg; }
	void	defocus_average(double v) { def_avg = v; }
	double	defocus_deviation() { return def_dev; }
	void	defocus_deviation(double v) { def_dev = v; }
	double	astigmatism_angle() { return ast_ang; }
	void	astigmatism_angle(double v) { ast_ang = v; }
	long	baseline_type() { return bt; }
	void	baseline_type(long t) { bt = t; }
	vector<double>&	baseline() { return base; }
	double	baseline(int i) { if ( i>=0 && i < NCTFPARAM ) return base[i]; else return 0; }
	void	baseline(int i, double d) { if ( i>=0 && i < NCTFPARAM ) base[i] = d; }
	void	baseline(double* b) { for ( int i=0; i<NCTFPARAM; i++ ) base[i] = b[i]; }
	void	baseline(vector<double>& b) { for ( int i=0; i<NCTFPARAM; i++ ) base[i] = b[i]; }
	long	envelope_type() { return et; }
	void	envelope_type(long t) { et = t; }
	vector<double>&	envelope() { return env; }
	double	envelope(int i) { if ( i>=0 && i < NCTFPARAM ) return env[i]; else return 0; }
	void	envelope(int i, double d) { if ( i>=0 && i < NCTFPARAM ) env[i] = d; }
	void	envelope(double* v) { for ( int i=0; i<NCTFPARAM; i++ ) env[i] = v[i]; }
	void	envelope(vector<double>& v) { for ( int i=0; i<NCTFPARAM; i++ ) env[i] = v[i]; }
	double	fom() { return f; }
	void	fom(double v) { f = v; }
	
	bool	check_defocus() {
		if ( def_avg < 1 || def_avg > 2e5 ) {
			cerr << "Error: Defocus is out of range 0.0001 - 20 um! (" << def_avg*1e-4 << ")" << endl;
			if ( def_avg < 1 ) def_avg = 1;
			else if ( def_avg > 2e5 ) def_avg = 2e5;
			return 1;
		} else return 0;
	}
	bool	check_Cs() {
		if ( cs < 1e3 || cs > 1e8 ) {
			cerr << "Error: Cs is out of range 0.0001 - 10 mm! (" << cs*1e-7 << ")" << endl;
			if ( cs < 1e3 ) cs = 1e3;
			else if ( cs > 1e8 ) cs = 1e8;
			return 1;
		} else return 0;
	}
	double	lambda() {		// Electron wavelength
//		if ( av ) wl = 12.2643/sqrt(av*(1+av*0.97845e-6));
 		if ( av ) wl = PLANCK*1e10L/sqrt(2.0*EMASS*ECHARGE*av*(1.0+ECHARGE*av/(2.0*EMASS*LIGHTSPEED*LIGHTSPEED)));
		return wl;
	}
	double	frequency_cutoff() {		// Aperture cutoff frequency
		return oa/(2*fl*lambda());
	}
	double	term1() { return t1; }
	double	term2() { return t2; }
	long double	delta_phi(double s2, double angle) {
		double		defocus = def_avg + def_dev*cos(2*(angle - ast_ang));
//		long double	dphi = (t1*s2 - t2*defocus)*s2 - ac;
//		while ( dphi < -M_PI) dphi += TWOPI;
//		while ( dphi > M_PI) dphi -= TWOPI;
//		return dphi;
		return		(t1*s2 - t2*defocus)*s2 - ac;
	}
	long double	calculate(double s2, double angle) {
//		double		dphi = delta_phi(s2, angle);
//		return		pf*sin(dphi) - af*cos(dphi);
//		dphi -= asin(af);
		return		sinl(delta_phi(s2, angle));
	}
	vector<double>	calculate(int nrad, int npsi, double step_size);
	Complex<double>	calculate_complex(double s2, double angle) {
		double		dphi = delta_phi(s2, angle);
//		return Complex<double>(af*cos(dphi),-pf*sin(dphi));
		return Complex<double>(cos(dphi),sin(dphi));
	}
	double	calc_baseline(double s) {
		double		s2(s*s), b(0), ds;		
		switch ( bt ) {
			case 1:
			case 4:
				b = base[0] + base[1]*s + base[2]*s2 + base[3]*s*s2 + base[4]*s2*s2;
				break;
			case 2:
			case 5:
				b = base[0] + base[1]*exp(base[2]*s2) + base[3]*exp(base[4]*s2);
				break;
			case 3:
			case 6:
				b = base[0] + base[1]*exp(base[2]*sqrt(s) + base[3]*s2);
				break;
			default: ;
		}
		if ( bt > 3 ) {
//			for ( int i=5; i<8; ++i ) cout << "\t" << base[i];
//			cout << endl;
			ds = base[6] - s;
			b += base[5] * exp(base[7]*ds*ds);
		}
		if ( !isfinite(b) )
			cerr << "Error in CTFparam::calc_baseline: baseline value not finite!" << endl;
		else if ( b <= 0 ) b = 0;
		return b;
	}
	double	calc_envelope(double s) {
		double		s2(s*s), e(0);
		switch ( et ) {
			case 1:
				e = env[0]*exp(env[1]*s2);
				break;
			case 2:
				e = env[0] + env[1]*exp(env[2]*s2);
				break;
			case 3:
				e = env[0]*exp(env[1]*s2) + env[2]*exp(env[3]*s2);
				break;
			case 4:
				e = env[0] + env[1]*exp(env[2]*s2) + env[3]*exp(env[4]*s2);
				break;
			default: ;
		}
		if ( !isfinite(e) )
			cerr << "Error in CTFparam::calc_envelope: baseline value not finite!" << endl;
		else if ( e <= 0 ) e = 0;
		return e;
	}
	vector<double>	zeroes(double max_s);
	vector<double>	maxima(double max_s);
	Bstring	baseline_equation();
	Bstring	envelope_equation();
	int		parse_baseline_equation(Bstring base_eq);
	int		parse_envelope_equation(Bstring env_eq);
	vector<double>	envelope_partial_coherence(long n, double freq_step);
	vector<double>	envelope_energy_spread(long n, double freq_step);
	double	zero(int i) {
		double		l = lambda();
		double		t = i - ac/M_PI;
		double		s = t/(l*def_avg);
		if ( cs >= 1e3 ) {
			double		l3Cs = l*l*l*cs;
			double		ctf_z = def_avg/(l*l*cs);
			s = ctf_z - sqrt(ctf_z*ctf_z - t*2.0/l3Cs);
		}
		s = 1.0/sqrt(s);
		if ( !isfinite(s) ) s = -1;
		if ( i == 1 ) firstzero = s;
		return s;
	}
	double	defocus_for_first_zero(double s) {
		if ( s < 1e-6 ) return -1;
		double		l = lambda();
		return		0.5*l*l*cs*s*s + 1/(l*s*s);
	}
	void	show() {
		cout << "Defocus average:                " << def_avg*1e-4 << " um" << endl;
		cout << "Defocus deviation:              " << def_dev*1e-4 << " um" << endl;
		cout << "Astigmatism angle:              " << ast_ang*180/M_PI << " degrees" << endl;
		cout << "Amplitude phase shift:          " << ac*180/M_PI << " degrees" << endl;
		cout << "Voltage:                        " << av*1e-3 << " kV" << endl;
		cout << "Wavelength:                     " << lambda() << " A" << endl;
		cout << "Spherical aberration (Cs):      " << cs*1e-7 << " mm" << endl;
		cout << "Chromatic aberration (Cc):      " << cc*1e-7 << " mm" << endl;
		cout << "Illumination halfangle (alpha): " << a*1e3 << " mrad" << endl;
		cout << "Energy spread:                  " << de << " eV" << endl;
		cout << "Energy filter slit width:       " << sw << " eV" << endl;
		cout << "Objective aperture:             " << oa*1e-4 << " µm" << endl;
		cout << "Focal length:                   " << fl*1e-7 << " mm" << endl;
	}
	void	show_baseline() {
		cout << "Baseline coefficients:          " << 
			baseline(0) << " " << baseline(1) << " " << 
			baseline(2) << " " << baseline(3);
		switch ( bt ) {
			case 1:		// Polynomial
			case 4:		
				cout << " " << baseline(4) << " (polynomial";
				break;
			case 2:		// Double gaussian
			case 5:		
				cout << " " << baseline(4) << " (double gaussian";
				break;
			case 3:		// EMAN style
			case 6:		
				cout << " (EMAN";
				break;
			default:
				break;
		}
		if ( baseline_type() > 3 ) cout << " with gaussian bump";
		cout << ")" << endl;
	}
	void	show_envelope() {
		cout << "Envelope coefficients:         ";
		for ( int i=0; i<=envelope_type(); i++ )
			cout << " " << envelope(i);
		switch ( et ) {
			case 1:
				cout << " (single gaussian)";
				break;
			case 2:
				cout << " (single gaussian with constant)";
				break;
			case 3:
				cout << " (double gaussian)";
				break;
			case 4:
				cout << " (double gaussian with constant)";
				break;
			default:
				break;
		}
		cout << endl;
	}
} ;
#define _CTFparam_
#endif

// Function prototypes
JSvalue			ctf_to_json(CTFparam& cp);
CTFparam		ctf_from_json(JSvalue& js);
CTFparam		ctf_from_json(string filename);
int				ctf_update_from_json(CTFparam& cp, JSvalue& js);
double			electron_wavelength(double volt);
double			beta2(double volt);
double			beta(double volt);
vector<double>	C_curve(long n, double freq_step);
vector<double>	defocus_range_profile(CTFparam& ctf, double freq_step, double def_min, double def_max, double def_step);
double			defocus_factor(CTFparam& ctf, double def_min, double def_max, double def_step);

