/**
@file	ps_ctf_plot.cpp
@brief	Postscript output functions for calculating CTF parameters.
@author Bernard Heymann
@date	Created: 20010515
@date	Modified: 20220221
**/
 
#include "ps_plot.h"
#include "ps_ctf_plot.h" 
#include "moving_average.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			ps_show_ctf_param(Bplot* plot, CTFparam& em_ctf)
{
	Bstring			txt;
	txt = Bstring(1e-3*em_ctf.volt(), "Accelerating voltage: %g kV");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-4*em_ctf.defocus_average(), "Defocus: %g um");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-7*em_ctf.Cs(), "Spherical aberration: %g mm");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-7*em_ctf.Cc(), "Chromatic aberration: %g mm");
	plot->page(0).add_text(txt);
	txt = Bstring(1e3*em_ctf.alpha(), "Beam spread/source size: %g mrad");
	plot->page(0).add_text(txt);
	txt = Bstring(em_ctf.dE(), "Energy spread: %g eV");
	plot->page(0).add_text(txt);
	
	return 0;
}

int			ps_show_ctf_fit(Bplot* plot, CTFparam& em_ctf)
{
	Bstring			txt;
	txt = Bstring(1e-3*em_ctf.volt(), "Accelerating voltage: %g kV");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-7*em_ctf.Cs(), "Spherical aberration: %g mm");
	plot->page(0).add_text(txt);
	txt = Bstring(em_ctf.amp_shift(), "Amplitude phase shift: %g");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-4*em_ctf.defocus_average(), "Defocus average: %g um");
	plot->page(0).add_text(txt);
	txt = Bstring(1e-4*em_ctf.defocus_deviation(), "Defocus deviation: %g um");
	plot->page(0).add_text(txt);
	txt = Bstring(em_ctf.astigmatism_angle()*180.0/M_PI, "Astigmatism angle: %g °");
	plot->page(0).add_text(txt);
	txt = "Baseline: " + em_ctf.baseline_equation();
	plot->page(0).add_text(txt);
	txt = "Envelope: " + em_ctf.envelope_equation();
	plot->page(0).add_text(txt);
	
	return 0;
}

/**
@brief 	Generates a postscript plot from CTF parameters for defocus and the envelope.
@param 	&filename	postscript file name.
@param 	&em_ctf		CTF parameter structure.
@param 	n			number of reciprocal space steps.
@param 	freq_step	spatial frequency increment per step.
@return int 		0.

	Postscript output is generated from microscope and defocus
	parameters.

**/
int 		ps_ctf_plot(Bstring& filename, CTFparam& em_ctf, size_t n, double freq_step)
{
	size_t			i, j, k, l, m, h, ncol(6);

	Bstring			title("Envelope curve");
	title = filename + ": " + title;
	
	Bplot*			plot = new Bplot(1, n, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("s");
	plot->page(0).column(1).label("Carbon");
	plot->page(0).column(2).label("PartialCoherence");
	plot->page(0).column(3).label("EnergySpread");
	plot->page(0).column(4).label("Envelope");
	plot->page(0).column(5).label("CTF");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(1).color(0,0,0);
	plot->page(0).column(2).color(0,0,1);
	plot->page(0).column(3).color(0,1,0);
	plot->page(0).column(4).color(1,0.5,0);
	plot->page(0).column(5).color(1,0,0);
	
	double			env;

	vector<double>	ec = C_curve(n, freq_step);
	vector<double>	epc = em_ctf.envelope_partial_coherence(n, freq_step);
	vector<double>	ees = em_ctf.envelope_energy_spread(n, freq_step);
	vector<double>	ctf = em_ctf.calculate(n, 1, freq_step);
	
	for ( i=0, j=n, k=2*n, l=3*n, m=4*n, h=5*n; i<n; i++, j++, k++, l++, m++, h++ ) {
		env = (ec[i]/ec[0]) * epc[i] * ees[i];
		(*plot)[i] = i*freq_step;
		(*plot)[j] = ec[i]/ec[0];
		(*plot)[k] = epc[i];
		(*plot)[l] = ees[i];
		(*plot)[m] = env;
		(*plot)[h] = env * ctf[i] * ctf[i];
//		cout << s << tab << env << endl;
	}
	
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(freq_step*n);
	plot->page(0).axis(1).label("Resolution (A)");
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).label("Power");

	ps_show_ctf_param(plot, em_ctf);
	
	ps_plot(filename, plot);
	
	delete plot;
	
	return 0;	
}

/**
@brief 	Generates a postscript plot of a contrast transfer function.
@param	n			number of elements in the radial power spectrum.
@param 	*rps		radial power spectrum.
@param	interval	spatial frequency step size.
@param 	*em_ctf		CTF parameter structure.
@param 	&filename	postscript file name.
@return int 		0.

	Postscript output is generated from fitted defocus, envelope and background
	parameters and compared to the radial average of the power spectrum image.

**/
int 		ps_ctf_plot(long n, double* rps, double interval,
				CTFparam* em_ctf, Bstring& filename)
{
	long			i, j, k, l, m, ncol(5);
	double			s, base, env;
	
//	double*			ctf = em_ctf->calculate(n, 1, interval);
	vector<double>	ctf = em_ctf->calculate(n, 1, interval);

	Bstring			title("Contrast transfer function");
	title = filename + ": " + title;
	
	Bplot*			plot = new Bplot(1, n, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("s");
	plot->page(0).column(1).label("RPS");
	plot->page(0).column(2).label("Baseline");
	plot->page(0).column(3).label("Envelope");
	plot->page(0).column(4).label("CTF");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(1).color(0,0,1);
	plot->page(0).column(2).color(0,1,0);
	plot->page(0).column(3).color(1,0.5,0);
	plot->page(0).column(4).color(1,0,0);
	
	double			min(0), max(0);
	
	for ( i=0, j=n, k=2*n, l=3*n, m=4*n; i<n; i++, j++, k++, l++, m++ ) {
		s = i*interval;
		base = em_ctf->calc_baseline(s);
		env = em_ctf->calc_envelope(s);
		(*plot)[i] = s;
		(*plot)[j] = rps[i] - base;
		(*plot)[k] = base - base;
		(*plot)[l] = env;
		(*plot)[m] = env * ctf[i] * ctf[i];
		if ( min > (*plot)[m] ) min = (*plot)[m];
		if ( max < (*plot)[m] ) max = (*plot)[m];
	}
	
//	max = em_ctf->calc_envelope(0);
	
//	delete[] ctf;

	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(interval*n);
	plot->page(0).axis(1).label("Resolution (A)");
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(3).min(min);
	plot->page(0).axis(3).max(max);
	plot->page(0).axis(3).label("Power");

	ps_show_ctf_fit(plot, *em_ctf);
	
	ps_plot(filename, plot);
	
	delete plot;
	
	return 0;	
}

/**
@brief 	Generates a postscript plot of defocus versus zeroes.
@param 	&filename	postscript output file name.
@param 	volts 		accelerating voltage (in volts).
@param 	Cs			spherical aberration (in angstrom).
@param 	amp_shift	amplitude contribution phase shift (radians).
@return int 		0.

	A postscript plot is generated from microscope and defocus parameters
	to show the corresponding defocus and zero values.
	Defocus is given in um.

**/
int 		ps_ctf_defocus_zeroes(Bstring& filename, double volts, double Cs, double amp_shift)
{
	CTFparam		cp(volts, Cs, amp_shift);
	int 			i, j, nzero(20);
	double			defocus, zero, x_max(20), x_inc(0.01), y_max(0);
	Bstring			title("Zeroes for different defocus values");
	Bstring			label;

	if ( verbose ) {
		cout << "Acceleration voltage:           " << volts/1000 << " kV" << endl;
		cout << "Electron wavelength:            " << cp.lambda() << " A" << endl;
		cout << "Spherical aberration:           " << Cs/1e7 << " mm" << endl;
		cout << "Amplitude phase shift:          " << amp_shift*180.0/M_PI << " degrees" << endl << endl;
	}
		
	long			nx = (long) (x_max/x_inc + 1);

	Bplot*			plot = new Bplot(1, nx, nzero+1);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(nzero+1);
	for ( i=0; i<=nzero; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Defocus");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<=nzero; i++ ) {
		label = Bstring(i, "Zero%d");
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).label(label);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(1).color(0,0,0);
	if ( nzero > 1 ) plot->page(0).column(2).color(1,0,0);
	if ( nzero > 2 ) plot->page(0).column(3).color(0,1,0);
	if ( nzero > 3 ) plot->page(0).column(4).color(0,0,1);
	if ( nzero > 4 ) plot->page(0).column(5).color(1,1,0);
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(x_max);
	plot->page(0).axis(1).label("Defocus (um)");
	plot->page(0).axis(3).label("Zero (A)");

	cp.show();

	for ( j=0, defocus=0; j<nx; defocus+=100, j++ ) {
		cp.defocus_average(defocus);
//		cp.show();
//		cout << defocus;
		(*plot)[j] = defocus/1e4;
		for ( i=1; i<=nzero; i++ ) {
			zero = cp.zero(i);
//			cout << tab << zero;
			if ( zero > y_max ) y_max = zero;
			(*plot)[i*nx+j] = zero;
		}
//		cout << endl;
	}

	y_max = 10*ceil(y_max/10);
	plot->page(0).axis(3).max(y_max);

	label = "Acceleration voltage:                " + Bstring(volts/1000, "%.0f") + " kV";
	plot->page(0).add_text(label);
	label = "Electron wavelength:                 " + Bstring(cp.lambda(), "%.4g") + " A";
	plot->page(0).add_text(label);
	label = "Spherical aberration coefficient: " + Bstring(Cs/1e7, "%.2g") + " mm";
	plot->page(0).add_text(label);
	label = "Amplitude phase shift:               " + Bstring(amp_shift*180.0/M_PI, "%.3g") + " degrees";
	plot->page(0).add_text(label);

	ps_plot(filename, plot);
	
	delete plot;	
	
	return 0;
}

int			point_spread(size_t n, vector<double>& a, double step)
{
	size_t			i, j;
	double			fac = M_PI*1.0L/n, v;
	vector<double>	t(n);
	
//	a[0] = 0;
	
	for ( i=0; i<n; i++ ) {
		t[i] = a[0]/2;
		v = fac*i;
		for ( j=1; j<n; j++ )
			t[i] += a[j]*cos(v*j);
		t[i] *= 2;
	}
	
	for ( i=0; i<n; i++ ) a[i] = t[i]*step;
	
	return 0;
}

/**
@brief 	Generates a point spread plot from CTF envelope parameters.
@param 	&filename	postscript file name.
@param 	&em_ctf		CTF parameter structure.
@param 	n			number of reciprocal space steps.
@param 	freq_step	spatial frequency increment per step.
@return int 		0.

	Postscript output is generated from microscope and defocus
	parameters.

**/
int 		ps_point_spread(Bstring& filename, CTFparam& em_ctf, size_t n, double freq_step)
{
	size_t			i, j, k, l, m, ncol(5);

	Bstring			title("Envelope curve");
	title = filename + ": " + title;
	
	Bplot*			plot = new Bplot(1, n, ncol);
	plot->title(title);
	
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Distance(A)");
	plot->page(0).column(1).label("Carbon");
	plot->page(0).column(2).label("PartialCoherence");
	plot->page(0).column(3).label("EnergySpread");
	plot->page(0).column(4).label("Envelope");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).axis(3);
		plot->page(0).column(i).element_size(0.5);
	}
	plot->page(0).column(1).color(0,0,0);
	plot->page(0).column(2).color(0,0,1);
	plot->page(0).column(3).color(0,1,0);
	plot->page(0).column(4).color(1,0.5,0);
//	plot->page(0).column(4).color(1,0,0);
	
	double			step(0.5/(n*freq_step));

	vector<double>	ec = C_curve(n, freq_step);
	vector<double>	epc = em_ctf.envelope_partial_coherence(n, freq_step);
	vector<double>	ees = em_ctf.envelope_energy_spread(n, freq_step);
	
	vector<double>	env(n);
	
	for ( i=0; i<n; i++ ) env[i] = ec[i] * epc[i] * ees[i];
	
	point_spread(n, ec, freq_step);
	point_spread(n, epc, freq_step);
	point_spread(n, ees, freq_step);
	point_spread(n, env, freq_step);
	
	for ( i=0, j=n, k=2*n, l=3*n, m=4*n; i<n; i++, j++, k++, l++, m++ ) {
		(*plot)[i] = i*step;
		(*plot)[j] = ec[i];
		(*plot)[k] = epc[i];
		(*plot)[l] = ees[i];
		(*plot)[m] = env[i];
//		cout << s << tab << env << endl;
	}
	
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(10);
	plot->page(0).axis(1).label("Distance (A)");
	plot->page(0).axis(3).label("Amplitude");
	
	ps_show_ctf_param(plot, em_ctf);
	
	ps_plot(filename, plot);
	
	delete plot;
	
	return 0;	
}

