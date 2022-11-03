/**
@file	ctf.h 
@brief	Header file for CTF (contrast transfer function) functions
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220606
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
	string			id;				// CTF group identifier
	long			sel;			// Selection integer
	double			av;				// Acceleration voltage (V)
	double			fl;				// Focal length (angstrom)
	double			cc;				// Chromatic aberation (angstrom)
	double			a;				// Illumination half-angle (radians)
	double			de;				// Energy spread (eV)
	double			oa;				// Objective aperture radius (angstrom)
	double			sw;				// Energy filter slit width (eV), 0 if not used
	double			tx, ty;			// Beam tilt (radian)
	map<pair<long,long>,double>	abw;	// Aberration weights
	long			bt;				// Baseline type (1=poly, 2=double_gauss, 3=eman)
	vector<double>	base;			// Baseline (up to 10 coefficients)
	long			et;				// Envelope type (1=single_gauss, 2=c+single_gauss, 3=double_gauss, 4=c+double_gauss)
	vector<double>	env;			// Envelope (up to 10 coefficients)
	double			f;				// CTF fit figure-of-merit
	double			wl;				// Wave length
	long double		t1;				// Term1: 0.5πl^3
	long double		t2;				// Term2: -πl
	void			initialize() {
		id = "1";
		sel = 1;
		a = de = sw = f = 0;
		av = 300000;
		fl = 3.5e7;
		oa = 2e6;
		cc = 2e7;
		tx = ty = 0;
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
		aberration_init();
	}
	void	aberration_init() {
		long			n(0), m;
		for ( n=0; n<5; n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw[{n,m}] = 0;
	}
	void			update_terms() {
		wl = lambda();
		t1 = M_PI_2*wl*wl*wl;
		t2 = -M_PI*wl;
	}
public:
	CTFparam() { initialize(); update_terms(); }
	CTFparam(string s) { initialize(); id=s; update_terms(); }
	CTFparam(double v, double sa, double ac) {
		initialize();
		volt(v); Cs(sa); amp_shift(ac);
		update_terms();
	}
	int		update(CTFparam* ctf) {
		if ( !ctf ) return -1;
		return update(*ctf);
	}
	int		update(CTFparam& ctf);
	
	void	identifier(string s) { id = s; }
	string&	identifier() { return id; }
	long	select() { return sel; }
	void	select(long i) { sel = i; }
	double	volt() { return av; }
	void	volt(double v) { av = v; update_terms(); }
	double	focal_length() { return fl; }
	void	focal_length(double v) { fl = v; }
	double	Cc() { return cc; }
	void	Cc(double v) { cc = v; }
	double	alpha() { return a; }
	void	alpha(double v) { a = v; }
	double	dE() { return de; }
	void	dE(double v) { de = v; }
	double	beam_tiltX() { return tx; }
	double	beam_tiltY() { return ty; }
	void	beam_tiltX(double v) { tx = v; }
	void	beam_tiltY(double v) { ty = v; }
	void	beam_tilt(double x, double y) { tx = x; ty = y; }
	double	aberration_weight(long n, long m) { return abw[{n,m}]; }
	void	aberration_weight(long n, long m, double v) { abw[{n,m}] = v; }
	void	add_aberration_weight(long n, long m, double v) { abw[{n,m}] += v; }
	map<pair<long,long>,double>&	aberration_weights() { return abw; }
	void	aberration_weights(map<pair<long,long>,double>& wa) {
		for ( auto w: wa ) abw[w.first] = w.second;
	}
	void	update_aberration_weights(map<pair<long,long>,double>& wa) {
		for ( auto w: wa ) abw[w.first] += w.second;
	}
	double	objective_aperture() { return oa; }
	void	objective_aperture(double v) { oa = v; }
	double	slit_width() { return sw; }
	void	slit_width(double v) { sw = v; }
	double	amp_shift() { return -abw[{0,0}]; }
	void	amp_shift(double v) { abw[{0,0}] = -v; }
	double	defocus_average() { return abw[{2,0}]/t2; }
	void	defocus_average(double v) { abw[{2,0}] = t2*v; }
	double	defocus_deviation() { return -sqrt(abw[{2,-2}]*abw[{2,-2}]+abw[{2,2}]*abw[{2,2}])/t2; }
	double	astigmatism_angle() { return angle_set_negPI_to_PI(atan2(-abw[{2,-2}], -abw[{2,2}])/2); }
	void	astigmatism(double dev, double ang) {
		dev = fabs(dev);
		abw[{2,2}] = t2*dev*cos(2*ang);
		abw[{2,-2}] = t2*dev*sin(2*ang);
	}
	double	Cs() { return abw[{4,0}]/t1; }
	void	Cs(double v) { abw[{4,0}] = t1*v; update_terms(); }
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
	string	aberration_weight_string() {
		string		ws;
		for ( auto w: abw ) ws += "," + to_string(w.second);
		ws[0] = '[';
		ws += ']';
		return ws;
	}
	vector<double>	aberration_even() {
		vector<double>	we;
		for ( auto w: abw ) if ( w.first.first%2 == 0 ) we.push_back(w.second);
		return we;
	}
	vector<double>	aberration_even_difference() {
		vector<double>	we(9,0);
		we[4] = abw[{4,-4}];
		we[5] = abw[{4,-2}];
		we[7] = abw[{4,2}];
		we[8] = abw[{4,4}];
		return we;
	}
	void	aberration_even(vector<double>& v) {
		long			i, n(0), m;
		for ( i=0; n<5 && i<v.size(); n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw[{n,m}] = v[i++];
	}
	void	aberration_even_update(vector<double>& v) {
		long			i, n(0), m;
		for ( i=0; n<5 && i<v.size(); n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw[{n,m}] += v[i++];
	}
//	vector<double>	aberration_odd() {
//		vector<double>	wo;
//		for ( auto w: abw ) if ( w.first.first%2 == 1 ) wo.push_back(w.second);
//		return wo;
//	}
	vector<double>	aberration_odd() {
		vector<double>	wo(6,0);
		long			i, n(1), m;
		for ( i=0; n<5 && i<wo.size(); n+=2 )
			for ( m=-n; m<=n; m+=2 )
				wo[i++] = abw[{n,m}] = 0;
		return wo;
	}
	void	aberration_odd(vector<double>& v) {
		long			i, n(1), m;
		for ( i=0; n<5 && i<v.size(); n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw[{n,m}] = v[i++];
	}
	void	aberration_odd_update(vector<double>& v) {
		long			i, n(1), m;
		for ( i=0; n<5 && i<v.size(); n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw[{n,m}] += v[i++];
	}
	void	convert_zernike() {
		abw[{0,0}] += -abw[{2,0}] + abw[{4,0}];
		abw[{1,-1}] -= abw[{3,-1}];
		abw[{1,1}] -= abw[{3,1}];
		abw[{2,-2}] -= 3*abw[{4,-2}];
		abw[{2,0}] = 2*abw[{2,0}] - 6*abw[{4,0}];
		abw[{2,2}] -= 3*abw[{4,2}];
		abw[{3,-1}] *= 3;
		abw[{3,1}] *= 3;
		abw[{4,-2}] *= 4;
		abw[{4,0}] *= 6;
		abw[{4,2}] *= 4;
	}
	void	add_zernike_even(vector<double>& v) {
		abw[{0,0}] += v[0] - v[2] + v[6];
		abw[{2,-2}] += v[1] - 3*v[5];
		abw[{2,0}] += 2*v[2] - 6*v[6];
		abw[{2,2}] += v[3] - 3*v[7];
		abw[{4,-4}] += v[4];
		abw[{4,-2}] += 4*v[5];
		abw[{4,0}] += 6*v[6];
		abw[{4,2}] += 4*v[7];
		abw[{4,4}] += v[8];
	}
	void	add_zernike_odd(vector<double>& v) {
		abw[{1,-1}] += v[0] - v[3];
		abw[{1,1}] += v[1] - v[4];
		abw[{3,-3}] += v[2];
		abw[{3,-1}] += 3*v[3];
		abw[{3,1}] += 3*v[4];
		abw[{3,3}] += v[5];
	}
	vector<double>	zernike_even() {
		vector<double>	ze(9,0);
		ze[0] = abw[{0,0}] + abw[{2,0}]/2 + abw[{4,0}]/3;
		ze[1] = abw[{2,-2}] + 0.75*abw[{4,-2}];
		ze[2] = (abw[{2,0}] + abw[{4,0}])/2;
		ze[3] = abw[{2,2}] + 0.75*abw[{4,2}];
		ze[4] = abw[{4,-4}];
		ze[5] = abw[{4,-2}]/4;
		ze[6] = abw[{4,0}]/6;
		ze[7] = abw[{4,2}]/4;
		ze[8] = abw[{4,4}];
		return ze;
	}
	vector<double>	zernike_odd() {
		vector<double>	zo(6,0);
		zo[0] = abw[{1,-1}] + abw[{3,-1}]/3;
		zo[1] = abw[{1,1}] + abw[{3,1}]/3;
		zo[2] = abw[{3,-3}];
		zo[3] = abw[{3,-1}]/3;
		zo[4] = abw[{3,1}]/3;
		zo[5] = abw[{3,3}];
		return zo;
	}
	void	delete_aberration(int which=3) {
		if ( which > 2 ) {
			abw.clear();
			return;
		}
		long					n(which%2), m;
		for ( ; n<5; n+=2 )
			for ( m=-n; m<=n; m+=2 )
				abw.erase({n,m});
	}
	bool	check_defocus() {
		double		def_avg(abw[{2,0}]/t2);
		if ( def_avg < 1 || def_avg > 2e5 ) {
			cerr << "Error: Defocus is out of range 0.0001 - 20 um! (" << def_avg*1e-4 << ")" << endl;
			if ( def_avg < 1 ) def_avg = 1;
			else if ( def_avg > 2e5 ) def_avg = 2e5;
			return 1;
		} else return 0;
	}
	double	lambda() {		// Electron wavelength
//		if ( av ) wl = 12.2643/sqrt(av*(1+av*0.97845e-6));
 		if ( av ) wl = PLANCK*1e10L/sqrt(2.0*EMASS*ECHARGE*av*(1.0+ECHARGE*av/(2.0*EMASS*LIGHTSPEED*LIGHTSPEED)));
 		if ( wl < 1e-10 ) {
 			cerr << "Warning: Wavelength is zero!" << endl;
 //			exit(-1);
		}
		return wl;
	}
	double	frequency_cutoff() {		// Aperture cutoff frequency
		return oa/(2*fl*lambda());
	}
	long double	calculate_aberration_even(double s2, double angle) {
		if (abw.size() < 5 ) return 0;	// At least the usual CTF parameters
		long		i(0), n(0), m(0);
		long double	dphi(abw[{0,0}]), s(s2);
		for ( i=1, n=2; n<5 && i<abw.size(); n+=2 ) {
			for ( m=-n; m<=n; m+=2, i++ ) {
				if ( m==0 ) dphi += abw[{n,m}]*s;
				else if ( m<0 ) dphi += abw[{n,m}]*s*sinl(-m*angle);
				else dphi += abw[{n,m}]*s*cosl(m*angle);
			}
			s *= s2;
		}
		return dphi;
	}
	long double	calculate_aberration_odd(double s2, double angle) {
		if (abw.size() < 2 ) return 0;	// At least beam tilt
		long		i(0), n(0), m(0);
		long double	dphi(0), s(sqrt(s2));
		for ( i=0, n=1; n<5 && i<abw.size(); n+=2 ) {
			for ( m=-n; m<=n; m+=2, i++ ) {
				if ( m<0 ) dphi += abw[{n,m}]*s*sinl(-m*angle);
				else dphi += abw[{n,m}]*s*cosl(m*angle);
			}
			s *= s2;
		}
		return dphi;
	}
	Complex<double>	aberration_odd_complex(double s2, double angle) {
		long double dphi = calculate_aberration_odd(s2, angle);
		return Complex<double>(cosl(dphi), sinl(dphi));
	}
	long double	delta_phi(double s2, double angle) {
		return calculate_aberration_even(s2, angle);
	}
/*	long double	delta_phi(double s2, double angle, int ewald=0) {
		double		defocus = def_avg + def_dev*cos(2*(angle - ast_ang));
		long double	dphi = (t1*s2 + t2*defocus)*s2 - ac;
		dphi += calculate_aberration_even(s2, angle);
//		while ( dphi < -M_PI) dphi += TWOPI;
//		while ( dphi > M_PI) dphi -= TWOPI;
		return dphi;
	}*/
	long double	calculate(double s2, double angle) {
		return		sinl(delta_phi(s2, angle));
	}
	vector<double>	calculate(int nrad, int npsi, double step_size);
	Complex<double>	calculate_complex(double s2, double angle) {
		double		dphi = delta_phi(s2, angle);
		return Complex<double>(cosl(dphi),sinl(dphi));
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
	int				parse_baseline_equation(Bstring base_eq);
	int				parse_envelope_equation(Bstring env_eq);
	double			partial_coherence(double s);
	vector<double>	envelope_partial_coherence(long n, double freq_step);
	double			energy_spread(double s2);
	vector<double>	envelope_energy_spread(long n, double freq_step);
	double			partial_coherence_and_energy_spread(double s2);
/*	double	zero(int i) {
		double		l = lambda();
		double		cs(Cs());
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
	}*/
	double	zero(int i) {
		double		s2;
//		abw[{0,0}] = ac;
//		abw[{2,0}] = t2*def_avg;
		double		t = (M_PI*i - abw[{0,0}]);
//		s2 = -t/abw[{2,0}];
		if ( abw[{4,0}] < 1 ) s2 = -t/abw[{2,0}];
		else s2 = (-abw[{2,0}] - sqrt(abw[{2,0}]*abw[{2,0}] - 4.0*abw[{4,0}]*t))/(2*abw[{4,0}]);
		double		r = 1.0/sqrt(s2);
		if ( !isfinite(r) ) r = -1;
		return r;
	}
/*	double	defocus_for_first_zero(double s) {
		if ( s < 1e-6 ) return -1;
		double		l = lambda();
		return		0.5*l*l*cs*s*s + 1/(l*s*s);
	}*/
	double	defocus_for_first_zero(double s) {
		if ( s < 1e-6 ) return -1;
		double		l = lambda();
		s *= s;
		return	(abw[{4,0}]*s + (abw[{0,0}] - M_PI)/s)/(M_PI*l);
	}
	void	show() {
		cout << "Optics group:                   " << id << endl;
		cout << "Number:                         " << sel << endl;
		cout << "Figure-of-merit:                " << f << endl;
		cout << "Defocus average:                " << defocus_average()*1e-4 << " um" << endl;
		cout << "Defocus deviation:              " << defocus_deviation()*1e-4 << " um" << endl;
		cout << "Astigmatism angle:              " << astigmatism_angle()*180/M_PI << " degrees" << endl;
		cout << "Amplitude phase shift:          " << -abw[{0,0}]*180/M_PI << " degrees" << endl;
		cout << "Voltage:                        " << av*1e-3 << " kV" << endl;
		cout << "Wavelength:                     " << lambda() << " A" << endl;
		cout << "Spherical aberration (Cs):      " << Cs()*1e-7 << " mm" << endl;
		cout << "Chromatic aberration (Cc):      " << cc*1e-7 << " mm" << endl;
		cout << "Beam tilt:                      " << tx << tab << ty << endl;
		cout << "Aberration parameters:\nn\tm\tw" << endl;
		for ( auto w: abw )
			cout << w.first.first << tab << w.first.second << tab << w.second << endl;
		cout << "Illumination halfangle (alpha): " << a*1e3 << " mrad" << endl;
		cout << "Energy spread:                  " << de << " eV" << endl;
		cout << "Energy filter slit width:       " << sw << " eV" << endl;
		cout << "Objective aperture:             " << oa*1e-4 << " µm" << endl;
		cout << "Focal length:                   " << fl*1e-7 << " mm" << endl;
		cout << "Frequency cutoff:               " << frequency_cutoff() << " (" << 1/frequency_cutoff() << " A)" << endl;
//		cout << "Ewald sphere application:       " << ew << endl;
		cout << endl;
	}
	void	show_aberration() {
		cout << "Optics group:                   " << id << endl;
		cout << "Number:                         " << sel << endl;
		cout << "Aberration parameters:\nn\tm\tw" << endl;
		for ( auto w: abw )
			cout << w.first.first << tab << w.first.second << tab << w.second << endl;
		cout << "Phase shift:                    " << -abw[{0,0}] << " radians" << endl;
		cout << "Image shift:                    " << abw[{1,1}]/TWOPI << tab << abw[{1,-1}]/TWOPI << " A" << endl;
		cout << "Defocus:                        " << defocus_average()*1e-4 << " um" << endl;
		cout << "Astigmatism:                    " << -abw[{2,2}]*1e-4/(M_PI*lambda()) << tab << -abw[{2,-2}]*1e-4/(M_PI*lambda()) << " um" << endl;
		cout << "Defocus deviation:              " << defocus_deviation()*1e-4 << " um" << endl;
		cout << "Astigmatism angle:              " << astigmatism_angle() << " degrees" << endl;
		cout << "Cs:                             " << Cs()*2e-7 << " mm" << endl;
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

