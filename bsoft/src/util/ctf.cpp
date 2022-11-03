/**
@file	ctf.cpp
@brief	Functions to manage CTF (contrast transfer function) parameters
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220322
**/

#include "ctf.h"
#include "string_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Updates a CTF structure from another.
@param 	ctf			CTF structure with new information.
@return int			0, <0 on error.
**/
int			CTFparam::update(CTFparam& ctf)
{
	int			i;
	
	if ( ctf.identifier().length() ) id = ctf.identifier();	// Optics group
	sel = ctf.sel;										// Optics group number
	if ( ctf.volt() ) volt(ctf.volt());					// Acceleration voltage (volts)
	if ( ctf.focal_length() ) fl = ctf.focal_length();	// Focal length (angstrom)
	if ( ctf.objective_aperture() ) oa = ctf.objective_aperture();		// Objective aperture (angstrom)
//	if ( ctf.Cs() ) Cs(ctf.Cs());						// Spherical aberation (angstrom)
	if ( ctf.Cc() ) cc = ctf.Cc();						// Chromatic aberation (angstrom)
	if ( ctf.beam_tiltX() ) tx = ctf.beam_tiltX();	// Beam tilt (radian)
	if ( ctf.beam_tiltY() ) ty = ctf.beam_tiltY();	// Beam tilt (radian)
	if ( ctf.alpha() ) a = ctf.alpha();					// Illumination half-angle (radians)
	if ( ctf.dE() ) de = ctf.dE();						// Energy spread (eV)
//	if ( ctf.amp_shift() ) amp_shift(ctf.amp_shift());	// CTF amplitude phase shift (radians)
//	if ( ctf.defocus_average() ) def_avg = ctf.defocus_average();		// CTF defocus average (angstrom)
//	if ( ctf.defocus_deviation() ) def_dev = ctf.defocus_deviation();	// CTF defocus deviation (angstrom)
//	if ( ctf.astigmatism_angle() ) ast_ang = ctf.astigmatism_angle();	// CTF astigmatism angle (radians)
	if ( ctf.baseline_type() ) {					// Baseline equation
		bt = ctf.baseline_type();
		for ( i=0; i<NCTFPARAM; i++ ) base[i] = ctf.baseline(i);
	}
	if ( ctf.envelope_type() ) {					// Envelope equation
		et = ctf.envelope_type();
		for ( i=0; i<NCTFPARAM; i++ ) env[i] = ctf.envelope(i);
	}
//	if ( ctf.aberration_even().size() ) we = ctf.aberration_even();
//	if ( ctf.aberration_odd().size() ) wo = ctf.aberration_odd();
	if ( ctf.abw.size() ) abw = ctf.abw;
	if ( ctf.f ) f = ctf.f;							// Figure-of-merit
//	ew = ctf.ew;									// Ewald sphere

	update_terms();
	
	return 0;
}

/**
@brief 	Calculates a CTF curve.
@param 	nrad			number of radii.
@param 	npsi			number of angular segments.
@param 	step_size 		reciprocal space step size (1/angstrom).
@return vector<double>	new CTF array.

	Defocus values are positive for underfocus.
	Functions:
		angle = atan(y/x)
		s2 = reciprocal space distance squared
		defocus_average = (defocus_max + defocus_min)/2
		defocus_deviation = (defocus_max - defocus_min)/2
		defocus = defocus_average + defocus_deviation*cos(2*(angle - astigmatism_angle))
		phase = 0.5*PI*lambda*lambda*lambda*Cs*s2*s2 - PI*lambda*defocus*s2 - amp_shift;
		CTF = sin(phase)
	The new CTF curve is returned.

**/
vector<double>	CTFparam::calculate(int nrad, int npsi, double step_size)
{
	vector<double>	ctf;

	if ( check_defocus() ) {
		cerr << "in CTFparam::calculate" << endl;
		return ctf;
	}
	
	if ( npsi < 1 ) npsi = 1;
	
	int 		i, j, k;
	double		s2, ang;
	double		dpsi(M_PI*2.0/npsi);
	
	ctf.resize(nrad*npsi);
	
//	defocus = def_avg;
	for ( j=k=0; j<npsi; j++ ) {
		ang = j*dpsi;
//		if ( npsi > 1 ) defocus = def_avg + def_dev*cos(2*(j*dpsi-ast_ang));
		for ( i=0; i<nrad; i++, k++ ) {
			s2 = i*step_size;
			s2 *= s2;
//			dphi = (t1*s2 + t2*defocus)*s2 - ac;
//			ctf[k] = sin(dphi);
			ctf[k] = calculate(s2, ang);
//			cout << k << tab << ctf[i] << endl;
		}
	}
	
	zero(1);
	
	return ctf;
}

/**
@brief 	Calculates the zeroes of a CTF curve on the spatial frequency scale.
@param 	max_s			maximum spatial frequency.
@return vector<double>		array of spatial frequencies for zeroes, NULL on error.

	The nth zero is given by the reciprocal space distance where the
	phase shift term is equal to -n*PI:
	phase = 0.5*PI*lambda^3*Cs*s^4 - PI*lambda*defocus*s^2 = -n*PI
	ctf_fz = defocus/(Cs*lambda^2)
	zero(n) = sqrt(ctf_fz - sqrt(ctf_fz^2 - 2.0*n/(Cs*lambda^3)))
	Defocus values are positive for underfocus.
	The array returned start with the first zero at index 0.

**/
/*vector<double>	CTFparam::zeroes(double max_s)
{
	vector<double>	zero;

	if ( check_defocus() ) {
		cerr << "in CTFparam::zeroes" << endl;
		return zero;
	}
	
	long		i;
	double		t, sm2(max_s*max_s);
    double		l = lambda();
	double		ps = abw[{0,0}]/M_PI;
	double		df = defocus_average();
    double		lf = l*df;
//	double		l3Cs = l*l*l*cs/2.0;
	double		l3Cs = abw[{4,0}]/M_PI;
	
//	cout << l*def_avg*sm2 << tab << l3Cs*sm2*sm2 << tab << ps << endl;
	
	long		nzero = (long) (lf*sm2 - l3Cs*sm2*sm2 + ps);
	
	if ( nzero < 1 ) nzero = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG CTFparam::zeroes: " << nzero << tab << max_s << tab << Cs() << endl;

	zero.resize(nzero);
	
	if ( abw[{4,0}] < 1 ) {
		for ( i=0; i<nzero; i++ ) {
			t = double(i)+1.0-ps;
			zero[i] = sqrt(t/lf);
		}
	} else {
		double		fz = df/(l*l*Cs()), fz2 = fz*fz;
		for ( i=0; i<nzero; i++ ) {
			t = (i+1.0-ps)/l3Cs;
			if ( fz2 > t )
				zero[i] = sqrt(fz - sqrt(fz2 - t));
			else
				zero[i] = -1;
			if ( !isfinite(zero[i]) ) zero[i] = 0.0001;
		}
	}

	if ( verbose & VERB_DEBUG )
		for ( i=0; i<nzero; i++ ) cout << i+1 << tab << zero[i] << tab << 1.0/zero[i] << endl;
	
	return zero;
}
*/
vector<double>	CTFparam::zeroes(double max_s)
{
	vector<double>	zero;

	if ( check_defocus() ) {
		cerr << "in CTFparam::zeroes" << endl;
		return zero;
	}
	
	long		i;
	double		t, z;

	if ( abw[{4,0}] < 1 ) {
		for ( i=0; i<100; i++ ) {
			t = abw[{0,0}] - M_PI*double(i+1);
			z = sqrt(t/abw[{2,0}]);
			if ( z <= max_s ) zero.push_back(z);
			else break;
		}
	} else {
		double		fz = -0.5*abw[{2,0}]/abw[{4,0}], fz2 = fz*fz;
//		cout << "fz = " << fz << endl;
		for ( i=0; i<100; i++ ) {
			t = (abw[{0,0}] - M_PI*double(i+1))/abw[{4,0}];
//			cout << i << tab << t << endl;
			if ( fz2 > t )
				z = sqrt(fz - sqrt(fz2 + t));
			else
				break;
			if ( !isfinite(z) )
				break;
			if ( z <= max_s )
				zero.push_back(z);
			else
				break;
		}
	}

	long		nzero(zero.size());

	if ( nzero < 1 ) {
		cerr << "Error in CTFparam::zeroes: No zeroes calculated!" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG CTFparam::zeroes: " << nzero << tab << max_s << endl;

	if ( verbose & VERB_DEBUG )
		for ( i=0; i<nzero; i++ ) cout << i+1 << tab << zero[i] << tab << 1.0/zero[i] << endl;
	
	return zero;
}
/*
vector<double>	CTFparam::zeroes(double max_s)
{
	vector<double>	zero;

	if ( check_defocus() ) {
		cerr << "in CTFparam::zeroes" << endl;
		return zero;
	}
	
	long		i;
	double		t, sm2(max_s*max_s);
	
	long		nzero = (long) ((abw[{0,0}] - abw[{2,0}]*sm2 - abw[{4,0}]*sm2*sm2)/M_PI);

	if ( nzero < 1 ) nzero = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG CTFparam::zeroes: " << nzero << tab << max_s << endl;

	zero.resize(nzero);
	
	if ( abw[{4,0}] < 1 ) {
		for ( i=0; i<nzero; i++ ) {
			t = abw[{0,0}] - M_PI*double(i+1);
			zero[i] = sqrt(t/abw[{2,0}]);
		}
	} else {
		double		fz = -0.5*abw[{2,0}]/abw[{4,0}], fz2 = fz*fz;
//		cout << "fz = " << fz << endl;
		for ( i=0; i<nzero; i++ ) {
			t = (abw[{0,0}] - M_PI*double(i+1))/abw[{4,0}];
//			cout << i << tab << t << endl;
			if ( fz2 > t )
				zero[i] = sqrt(fz - sqrt(fz2 + t));
			else
				zero[i] = -1;
			if ( !isfinite(zero[i]) ) zero[i] = 0.0001;
		}
	}

	if ( verbose & VERB_DEBUG )
		for ( i=0; i<nzero; i++ ) cout << i+1 << tab << zero[i] << tab << 1.0/zero[i] << endl;
	
	return zero;
}
*/

/**
@brief 	Calculates the maxima of a CTF curve on the spatial frequency scale.
@param 	max_s			maximum spatial frequency.
@return vector<double>		array of spatial frequencies for maxima, NULL on error.

	Thus uses the zeroes function to find the maxima.
	Defocus values are positive for underfocus.
	The array returned start with the first zero at index 0.

**/
vector<double>	CTFparam::maxima(double max_s)
{
	long			i;
	
	vector<double>	zero = zeroes(max_s);
	
	vector<double>	max(zero.size(),0);
	
	max[0] = zero[0]*0.75;
	
	for ( i=1; i<zero.size(); ++i )
		max[i] = (zero[i] + zero[i-1])/2;

	if ( verbose & VERB_DEBUG )
		for ( i=0; i<max.size(); ++i ) cout << i+1 << tab << max[i] << tab << 1.0/max[i] << endl;
	
	return max;
}

/**
@brief 	Generates a baseline equation string from 4 to 10 coefficients.
@return Bstring				baseline equation.

	The baseline string is constructed from the 4 to 10 coefficients
	according to the specified type.

**/
Bstring		CTFparam::baseline_equation()
{
	if ( base[0] <= 0 && base[1] == 0 ) base[0] = 1;		// Default baseline
	
	Bstring		bleq(base[0], "%.4g");
	
	switch ( bt ) {
		case 1:		// Polynomial
		case 4:		// Polynomial
			bleq += Bstring(base[1], " + %.4g*$s") + Bstring(base[2], " + %.4g*$s2")
				+ Bstring(base[3], " + %.4g*$s3") + Bstring(base[4], " + %.4g*$s4");
			break;
		case 2:		// Double gaussian
		case 5:		// Double gaussian
			bleq += Bstring(base[1], " + %.4g") + Bstring(base[2], "*exp(%.4g*$s2)") 
				+ Bstring(base[3], " + %.4g") + Bstring(base[4], "*exp(%.4g*$s2)");
			break;
		case 3:		// EMAN style
		case 6:		// EMAN style
			bleq += Bstring(base[1], " + %.4g") + Bstring(base[2], "*exp(%.4g*sqrt($s)") 
				+ Bstring(base[3], " + %.4g*$s2)");
			break;
		default: break;
	}

//	for ( int i=5; i<8; ++i ) cout << tab << base[i];
//	cout << endl;
	
	if ( bt > 3 )
		bleq += Bstring(base[5], " + %.4g") + Bstring(base[7], "*exp(%.4g")
			 + Bstring(base[6], "*(%.4g - $s)") + Bstring(base[6], "*(%.4g - $s))");
	
//	cout << bleq << endl;
	
	return bleq;
}

/**
@brief 	Generates an envelope equation string from 4 coefficients.
@return Bstring				envelope equation.

	The envelope string is constructed from the 4 double gaussian parameters.

**/
Bstring		CTFparam::envelope_equation()
{
	Bstring		enveq(env[0], "%.4g");

	switch ( et ) {
		case 1:		// Single gaussian
			enveq += Bstring(env[1], "*exp(%.4g*$s2)");
			break;
		case 2:		// Single gaussian with constant
			enveq += Bstring(env[1], " + %.4g") + Bstring(env[2], "*exp(%.4g*$s2)");
			break;
		case 3:		// Double gaussian
			enveq += Bstring(env[1], "*exp(%.4g*$s2)") +
				Bstring(env[2], " + %.4g") + Bstring(env[3], "*exp(%.4g*$s2)");
			break;
		case 4:		// Double gaussian with constant
			enveq += Bstring(env[1], " + %.4g") + Bstring(env[2], "*exp(%.4g*$s2)") +
				Bstring(env[3], " + %.4g") + Bstring(env[4], "*exp(%.4g*$s2)");
			break;
		default: break;
	}

	return enveq;
}

/**
@brief 	Extracts the coefficients from the baseline string.
@param 	base_eq			baseline equation string.
@return int				equation type.

	The baseline string is scanned in one of the three supported formats.

**/
int			CTFparam::parse_baseline_equation(Bstring base_eq)
{
//	cout << "Base eq = " << base_eq << endl;

	base_eq = base_eq.remove('\"');
		
	int 		n(0);
		
	if ( base_eq.contains("s4") ) {				// Polynomial
		bt = 1;	
		sscanf(base_eq.c_str(), "%lf + %lf*$s + %lf*$s2 + %lf*$s3 + %lf*$s4%n",
			   &base[0], &base[1], &base[2], &base[3], &base[4], &n);
	} else if ( base_eq.contains("sqrt") ) {	// EMAN style
		bt = 3;
		sscanf(base_eq.c_str(), "%lf + %lf*exp(%lf*sqrt($s) + %lf*$s2)%n",
			   &base[0], &base[1], &base[2], &base[3], &n);
	} else if ( base_eq.contains("exp") ) {  	// Double gaussian
		bt = 2;
		sscanf(base_eq.c_str(), "%lf + %lf*exp(%lf*$s2) + %lf*exp(%lf*$s2)%n",
			   &base[0], &base[1], &base[2], &base[3], &base[4], &n);
	} else {
		bt = 1;
		base[0] = 1;
		base[1] = base[2] = base[3] = base[4] = 0;
	}

	char*		ptr = &base_eq[n];
	if ( strstr(ptr, "exp") ) {
//		cout << ptr << endl;
		sscanf(ptr, " + %lf*exp(%lf*(%lf - $s)*(%lf - $s))",
			&base[5], &base[7], &base[6], &base[6]);
		bt += 3;
//		cout << base[5] << tab << base[6] << tab << base[7] << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG parse_baseline_equation: bt=" << bt << " baseline=" <<
			   base[0] << "," << base[1] << "," << base[2] << "," <<
			   base[3] << "," << base[4] << endl;
	
	return bt;
}

/**
@brief 	Extracts the coefficients from the envelope string.
@param 	env_eq			envelope equation string.
@return int				0.

	The envelope string is scanned to get the 4 double gaussian parameters.

**/
int			CTFparam::parse_envelope_equation(Bstring env_eq)
{
	env_eq = env_eq.remove('"');
	
	et = 2*env_eq.count('s') - 1;
	if ( env_eq.index('+') > 0 )
		if ( env_eq.index('+') < env_eq.index('*') ) et++;
	
//	cout << et << tab << env_eq << endl;

	switch ( et ) {
		case 1:		// Single gaussian
			sscanf(env_eq.c_str(), "%lf*exp(%lf*$s2)",
		   		&env[0], &env[1]);
			break;
		case 2:		// Single gaussian with constant
			sscanf(env_eq.c_str(), "%lf + %lf*exp(%lf*$s2)",
		   		&env[0], &env[1], &env[2]);
			break;
		case 3:		// Double gaussian
			sscanf(env_eq.c_str(), "%lf*exp(%lf*$s2) + %lf*exp(%lf*$s2)",
		   		&env[0], &env[1], &env[2], &env[3]);
			break;
		case 4:		// Double gaussian with constant
			sscanf(env_eq.c_str(), "%lf + %lf*exp(%lf*$s2) + %lf*exp(%lf*$s2)",
		   		&env[0], &env[1], &env[2], &env[3], &env[4]);
			break;
		default:
			et = 4;
			env[0] = 1;		// Default envelope
			env[1] = 0;
			env[2] = -10;
			env[3] = 0;
			env[4] = -10;
	}
	
	return 0;
}

/**
@brief 	Calculates the envelope curve based on partial coherence.
@param 	s			spatial frequency.
@return double		amplitude.

The curve is calculated at frequency s as:
	Epc(s) = exp(-(pi*alpha*(Cs*lamda^2*s^2 - def)*s)^2)
where
	Cs: Spherical aberation (~2e7 A)
	alpha: Beam spread/source size (~0.1 mrad)
	def: Defocus (~1e4 A)
	lamda: electron wavelength (~0.03 A)

References:
	Zhu et al. (1997) JSB 118, 197-219.
**/
double		CTFparam::partial_coherence(double s)
{
    double			Csl2(Cs()*wl*wl);
	double			pal(M_PI*a*wl);
	double			arg = pal*(Csl2*s*s - defocus_average())*s;
	
	return exp(-arg*arg);
}

/**
@brief 	Calculates the envelope curve based on partial coherence.
@param 	n				number of spatial frequency steps.
@param 	freq_step		size of spatial frequency step.
@return vector<double>		curve.

The curve is calculated at frequency s as:
	Epc(s) = exp(-(pi*alpha*(Cs*lamda^2*s^2 - def)*s)^2)
where
	Cs: Spherical aberation (~2e7 A)
	alpha: Beam spread/source size (~0.1 mrad)
	def: Defocus (~1e4 A)
	lamda: electron wavelength (~0.03 A)

References:
	Zhu et al. (1997) JSB 118, 197-219.
**/
vector<double>	CTFparam::envelope_partial_coherence(long n, double freq_step)
{
	long			i;
    double			l = lambda(), Csl2(Cs()*l*l);
	double			s, arg, pal(M_PI*a);
	vector<double>	curve(n);
	
	for ( i=0; i<n; i++ ) {
		s = i*freq_step;
		arg = pal*(Csl2*s*s - defocus_average())*s;
		curve[i] = exp(-arg*arg);
	}
	
	return curve;
}

/**
@brief 	Calculates the envelope curve based on energy spread.
@param 	s2			spatial frequency squared.
@return double		amplitude.

The curve is calculated at frequency s as:
	Ees(s) = exp(-0.5*(pi*lambda*Cc*(dE/V)*s^2)^2)
where
	Cc: Chromatic aberation (~2e7 A)
	dE: Energy spread (~1 eV)
	V: Acceleration voltage (~1e5 V)
	lamda: electron wavelength (~0.03 A)

References:
	Freitag et al. (2005) Ultramicroscopy 102, 209-14.
	Zhu et al. (1997) JSB 118, 197-219.
**/
double		CTFparam::energy_spread(double s2)
{
//	double			fac(1.0/(16*log(2.0)));
	double			fac(0.5);
    double			arg = M_PI*wl*cc*de*s2/av;
	
	return exp(-fac*arg*arg);
}

/**
@brief 	Calculates the envelope curve based on energy spread.
@param 	n				number of spatial frequency steps.
@param 	freq_step		size of spatial frequency step.
@return vector<double>	curve.

The curve is calculated at frequency s as:
	Ees(s) = exp(-0.5*(pi*lambda*Cc*(dE/V)*s^2)^2)
where
	Cc: Chromatic aberation (~2e7 A)
	dE: Energy spread (~1 eV)
	V: Acceleration voltage (~1e5 V)
	lamda: electron wavelength (~0.03 A)

References:
	Freitag et al. (2005) Ultramicroscopy 102, 209-14.
	Zhu et al. (1997) JSB 118, 197-219.
**/
vector<double>	CTFparam::envelope_energy_spread(long n, double freq_step)
{
	long			i;
	double			s, arg;
//	double			fac(1.0/(16*log(2.0)));
	double			fac(0.5);
    double			l = lambda(), deccl(M_PI*l*cc*de/av);
	vector<double>	curve(n);
	
	for ( i=0; i<n; i++ ) {
		s = i*freq_step;
		arg = deccl*s*s;
		curve[i] = exp(-fac*arg*arg);
	}
	
	return curve;
}

/**
@brief 	Calculates the envelope curve based on partial coherence and energy spread.
@param 	s2			spatial frequency squared.
@return double		amplitude.

References:
	Freitag et al. (2005) Ultramicroscopy 102, 209-14.
	Zhu et al. (1997) JSB 118, 197-219.
**/
double	CTFparam::partial_coherence_and_energy_spread(double s2)
{
	double			s(sqrt(s2));
	double			pc = M_PI*a*wl*(Cs()*wl*wl*s2 - defocus_average())*s;
//	double			fac(1.0/(16*log(2.0)));
	double			fac(0.5);
    double			es = M_PI*wl*cc*de*s2/av;
	
	return exp(-pc*pc-fac*es*es);
}

/**
@brief 	Converts microscope parameters to a JSON object list.
@param 	&cp			microscope parsmeter structure.
@return JSvalue		JSON object list.
**/
JSvalue		ctf_to_json(CTFparam& cp)
{
	JSvalue		js(JSobject);
	
	js["Optics_group"] = cp.identifier();
	js["Volt"] = cp.volt();
	js["Wavelength"] = cp.lambda();
	js["Focal_length"] = cp.focal_length();
	js["Objective_aperture"] = cp.objective_aperture();
	js["Slit_width"] = cp.slit_width();
	js["Cs"] = cp.Cs();
	js["Cc"] = cp.Cc();
	js["Beam_convergence"] = cp.alpha();
	js["Energy_spread"] = cp.dE();
	js["Amplitude_phase"] = cp.amp_shift();
	js["Defocus_average"] = cp.defocus_average();
	js["Defocus_deviation"] = cp.defocus_deviation();
	js["Astigmatism_angle"] = cp.astigmatism_angle()*180.0/M_PI;
	js["Aberration_even"] = "[" + concatenate(cp.aberration_even_difference()) + "]";
	js["Aberration_odd"] = "[" + concatenate(cp.aberration_odd()) + "]";
	
	return js;
}

/**
@brief 	Converts a JSON object list to microscope parameters.
@param 	filename	JSON file name.
@return CTFparam		microscope parsmeter structure.
**/
CTFparam	ctf_from_json(string filename)
{
	JSvalue			jsctf(JSobject);
	if ( filename.length() ) jsctf = JSparser(filename).parse();
	CTFparam		cp = ctf_from_json(jsctf);
	return cp;
}

/**
@brief 	Converts a JSON object list to microscope parameters.
@param 	&js			JSON object list.
@return CTFparam		microscope parameter structure.
**/
CTFparam	ctf_from_json(JSvalue& js)
{
	CTFparam		cp;
	ctf_update_from_json(cp, js);
	
	return cp;
}

/**
@brief 	Converts a JSON object list to microscope parameters.
@param	&cp			CTF parameter structure.
@param 	&js			JSON object list.
@return int			0.
**/
int			ctf_update_from_json(CTFparam &cp, JSvalue &js)
{
	if ( js.exists("Optics_group") ) cp.identifier(js["Optics_group"].value());
	if ( js.exists("Volt") ) cp.volt(js["Volt"].real());
	if ( js.exists("Focal_length") ) cp.focal_length(js["Focal_length"].real());
	if ( js.exists("Objective_aperture") ) cp.objective_aperture(js["Objective_aperture"].real());
	if ( js.exists("Slit_width") ) cp.slit_width(js["Slit_width"].real());
	if ( js.exists("Cs") ) cp.Cs(js["Cs"].real());
	if ( js.exists("Cc") ) cp.Cc(js["Cc"].real());
	if ( js.exists("Beam_convergence") ) cp.alpha(js["Beam_convergence"].real());
	if ( js.exists("Energy_spread") ) cp.dE(js["Energy_spread"].real());
	if ( js.exists("Amplitude_phase") ) cp.amp_shift(js["Amplitude_phase"].real());
	if ( js.exists("Defocus_average") ) cp.defocus_average(js["Defocus_average"].real());
	if ( js.exists("Defocus_deviation") && js.exists("Astigmatism_angle") )
		cp.astigmatism(js["Defocus_deviation"].real(), js["Astigmatism_angle"].real()*M_PI/180.0);
	if ( js.exists("Aberration_even") ) {
		vector<double>	v = parse_real_vector(js["Aberration_even"].value().substr(1));
		cp.aberration_even_update(v);
	}
	if ( js.exists("Aberration_odd") ) {
		vector<double>	v = parse_real_vector(js["Aberration_odd"].value().substr(1));
		cp.aberration_odd(v);
	}
//	cout << js["Cs"].real() << tab << cp.Cs() << endl;

	return 0;
}

/**
@brief 	Calculates the wavelength in angstrom from the acceleration voltage.
@param 	volt		acceleration voltage.
@return double		wavelength in angstrom, <0 on error.

	                   12.26
	lambda = ----------------------------
	         sqrt(volt*(1+volt*0.9788e-6)
	                       1e10 * h
	lambda = ---------------------------------------
	         sqrt(2*me*e*volt*(1+e*volt/(2*me*c^2)))

**/
double		electron_wavelength(double volt)
{
	if ( volt < 1 ) {
		error_show("Error in electron_wavelength", __FILE__, __LINE__);
		cerr << "The acceleration voltage is below 1 V! (" << volt << " V)" << endl;
		return -1;
	}
	
 //   return 12.2642596/sqrt(volt*(1+volt*0.978475598e-6));
 	
 	return 1e10*PLANCK/sqrt(2*EMASS*ECHARGE*volt*(1+ECHARGE*volt/(2*EMASS*LIGHTSPEED*LIGHTSPEED)));
}

/**
@brief 	Calculates the square of the relative electron velocity.
@param 	volt		acceleration voltage.
@return double		square of the relative electron velocity.

	                     1
	beta2 = 1 - ---------------------
	            (1+e*volt/(me*c^2))^2

**/
double			beta2(double volt)
{
	double		f, b2(0);
	
	f = 1 + ECHARGE*volt/(EMASS*LIGHTSPEED*LIGHTSPEED);
	
	b2 = 1 - 1/(f*f);
	
	return b2;
}

/**
@brief 	Calculates the relative electron velocity.
@param 	volt		acceleration voltage.
@return double		relative electron velocity.

	                          1
	beta = sqrt(1 - ---------------------)
	                (1+e*volt/(me*c^2))^2

**/
double			beta(double volt)
{
	return sqrt(beta2(volt));
}


/**
@brief 	Calculates the carbon scattering curve.
@param 	n				length of curve.
@param 	freq_step		frequency step size.
@return vector<double>	array with C-curve.

	The curve is based on the sum of five Gaussians.

**/
vector<double>	C_curve(long n, double freq_step)
{
	long			i;
	double			s2;
	vector<double>	curve(n);
	
	for ( i=0; i<n; i++ ) {
		s2 = i*freq_step;
		s2 *= s2;
		curve[i] = 0.893*exp(-0.2465*s2) + 0.2563*exp(-1.7100*s2) + 
					0.7570*exp(-6.4094*s2) + 1.0487*exp(-18.6113*s2) + 
					0.3575*exp(-50.2523*s2);
	}
	
	return curve;
}

/**
@brief 	Calculates a weighting profile associated with a range of defocus values.
@param 	&ctf			CTF parameters.
@param 	freq_step		frequency step size.
@param 	def_min			minimum defocus.
@param 	def_max			maximum defocus.
@param 	def_step		defocus step.
@return vector<double>	sum of all CTF curves.

	The CTF curves are calculated up to the frequency cutoff determined by the objective aperture.

**/
vector<double>	defocus_range_profile(CTFparam& ctf, double freq_step, double def_min, double def_max, double def_step)
{
	double			scut = ctf.frequency_cutoff();
	long			n(scut/freq_step+1);

	if ( def_min < 1 ) {
		vector<double>	c(n, 0.5);
		return c;
	}

	long			i, m(0);
	double			def, s2, v;
	vector<double>	c(n, 0);
	
	double			def_ori = ctf.defocus_average();
	
	for ( def = def_min; def <= def_max; def += def_step, ++m ) {
		ctf.defocus_average(def);
		for ( i=0; i<n; ++i ) {
			s2 = i*freq_step;
			s2 *= s2;
			v = ctf.calculate(s2, 0);
			c[i] += v*v;
		}
	}
	
	ctf.defocus_average(def_ori);	// Preserve the original defocus

	for ( i=0; i<n; ++i )
		c[i] /= m;
	
	return c;
}

/**
@brief 	Calculates an integrated defocus weighting factor for a range of defocus values.
@param 	&ctf			CTF parameters.
@param 	def_min			minimum defocus.
@param 	def_max			maximum defocus.
@param 	def_step		defocus step.
@return vector<double>	sum of all CTF curves.

	The CTF curves are calculated up to the frequency cutoff determined by the objective aperture.

**/
double		defocus_factor(CTFparam& ctf, double def_min, double def_max, double def_step)
{
	if ( def_min < 1 ) return 0.5;

	double				def_fac(0);
	vector<double>		c = defocus_range_profile(ctf, 0.001, def_min, def_max, def_step);
	
	for ( auto it = c.begin(); it != c.end(); ++it ) def_fac += *it;
	
	def_fac /= c.size();
	
	return def_fac;
}


