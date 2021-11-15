/**
@file	mg_ctf.cpp
@brief	Functions for CTF (contrast transfer function) processing
@author Bernard Heymann
@date	Created: 19970715
@date	Modified: 20210817
**/

#include "rwimg.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "ps_plot.h" 
#include "ps_ctf_plot.h" 
#include "mg_processing.h"
//#include "matrix_linear.h"
#include "utilities.h"
#include "timer.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates a CTF complex image.
@param 	cp				CTF parameters.
@param 	action			type of CTF calculated (1-8).
@param 	wiener			Wiener factor (fraction).
@param 	size			new image size.
@param 	sam				new image pixel size.
@param 	lores			low resolution limit.
@param 	hires			high resolution limit.
@return Bimage*			new CTF function image.

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
Bimage*		img_ctf_calculate(CTFparam cp, int action, double wiener, Vector3<long> size, 
				Vector3<double> sam, double lores, double hires)
{
//	if ( cp.check_defocus() || cp.check_Cs() )
//		cerr << "in img_ctf_calculate" << endl;
	
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;

	double			shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	long 			i, x, y, z, sign;
	double			sx, sy, sz, s, s2, w, ctf_fac, base(0), env(0), ctf_env;
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
//	cout << h << endl;
//	cout << freq_scale << endl;

	Bstring			base_eq = cp.baseline_equation();
	Bstring			env_eq = cp.envelope_equation();
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF function:" << endl;
		switch ( action ) {
			case 1: cout << "Flipping phases" << endl; break;
			case 2: cout << "Applying a CTF" << endl; break;
			case 3: cout << "Correcting for the CTF" << endl; break;
			case 4: cout << "Correcting for the CTF with an optimal Wiener filter" << endl; break;
			case 5: cout << "Correcting for the CTF with baseline enhancement" << endl; break;
			case 6: cout << "Correcting for the CTF with baseline2 enhancement" << endl; break;
			case 7: cout << "Flipping phases with baseline enhancement" << endl; break;
			case 8: cout << "Correcting for the CTF with baseline3 enhancement" << endl; break;
//			case 9: cout << "Applying theoretical envelopes" << endl; break;
//			case 10: cout << "Applying a CTF with theoretical envelopes" << endl; break;
			default: break;
		}
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		if ( action == 3 || action == 5 || action == 6 )
			cout << "Wiener factor:                  " << wiener << endl;
		if ( action >= 4 ) {
			cout << "Envelope:                       " << env_eq << endl;
			cout << "Noise:                          " << base_eq << endl;
		}
		cout << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[0]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[0]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}
	
		for ( i=z=0; z<p->sizeZ(); z++ ) {
			sz = z;
			if ( z > h[2] ) sz -= p->sizeZ();
			sz *= freq_scale[2];
			for ( y=0; y<p->sizeY(); y++ ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					s = sqrt(s2);
					if ( s >= slo && s <= shi ) {
						ctf_fac = cp.calculate(s2, atan2(sy,sx));
//						cout << s << tab << ctf_fac << endl;
						w = 0;
						sign = ( ctf_fac < 0 )? -1: 1;
						if ( action > 3 ) base = cp.calc_baseline(s);
						switch ( action ) {			
							case 1: w = sign; break;			// Flip the phase
							case 2: w = ctf_fac; break;			// Apply the CTF
							case 3:
								w = ctf_fac/(ctf_fac*ctf_fac + wiener);
								break; // Spider method
							case 4:
								if ( base > 0 ) {
									env = cp.calc_envelope(s);
									ctf_env = ctf_fac*env;
									w = ctf_fac*sqrt(env)/(ctf_fac*ctf_env + base);
								}
								break; // Wiener filter method
							case 5:
								w = ctf_fac/(ctf_fac*ctf_fac*base + wiener);
								break; // Baseline enhancement
							case 6:
								w = 1/(ctf_fac*base + sign*wiener);
								break; // Baseline enhancement
							case 7:
								if ( base > 0 ) w = sign/sqrt(base);
								break;	// Flip phase with baseline enhancement
							case 8:
								env = cp.calc_envelope(s);
								w = sign/sqrt(ctf_fac*ctf_fac*env + base);
								break;	// Baseline enhancement
							default: break;
						}
						p->set(i, w);
					}
				}
			}
		}
	
	p->statistics();

	return p;
}

/**
@brief 	Calculates a wave aberration function.
@param 	cp				CTF parameters.
@param 	size			new image size.
@param 	sam				new image pixel size.
@return Bimage*			new wave aberration function image.

	Functions:
		angle = atan(y/x)
		s2 = x*x + y*y
		defocus_average = (defocus_max + defocus_min)/2
		defocus_deviation = (defocus_max - defocus_min)/2
		defocus = defocus_average + defocus_deviation*cos(2*(angle - astigmatism_angle))
		phase = 0.5*PI*lambda*lambda*lambda*Cs*s2*s2 - PI*lambda*defocus*s2;
	Note: Defocus is positive for underfocus and negative for overfocus.

**/
Bimage*		img_wave_aberration(CTFparam cp, Vector3<long> size, Vector3<double> sam)
{
	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	long 			i, n, x, y, z;
	double			sx, sy, sz, s2;
	Vector3<double>	freq_scale(1.0/p->real_size());
	Vector3<double>	h((p->sizeX() - 1)/2, (p->sizeY() - 1)/2, (p->sizeZ() - 1)/2);

	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a wave aberration function:" << endl;
		cout << "Defocus average & deviation:    " << cp.defocus_average() << " A  " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180/M_PI << " degrees" << endl;
		cout << "Voltage:                        " << cp.volt() << " V" << endl;
		cout << "Wavelength:                     " << cp.lambda() << " A" << endl;
		cout << "Cs:                             " << cp.Cs() << " A" << endl;
		cout << "Amplitude contrast:             " << sin(cp.amp_shift()) << endl;
		cout << endl;
	}
	
	for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			sz = z;
			if ( z > h[2] ) sz -= p->sizeZ();
			sz *= freq_scale[2];
			for ( y=0; y<p->sizeY(); y++ ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					p->set(i, cp.delta_phi(s2, atan2(sy,sx)));
				}
			}
		}
	}

	return p;
}

/**
@brief 	Applies or corrects for the contrast transfer function (CTF).
@param 	*p				image (modified).
@param 	*em_ctf			CTF parameter structure.
@param 	action			action to be taken.
@param 	wiener			Wiener factor (fraction).
@param 	lores			low resolution limit.
@param 	hires			high resolution limit.
@return int				0, <0 on error.

	The actions for this funtion are:
	1	flip phase (multiply by the sign of the CTF)
	2	apply a CTF (multiply with the CTF)
	3	correct for the CTF: ctf/(ctf^2 + wiener_factor)
	4	correct for the CTF: env*ctf/((env*ctf)^2 + noise^2)
	5	correct for the CTF with baseline: ctf/(ctf^2*noise^2 + wiener_factor)
	6	correct for the CTF with baseline: 1/(ctf*noise + sign*wiener_factor)

**/
int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires)
{
	if ( action < 1 || action > 10 ) return 0;
	
	if ( wiener < 0.01 ) wiener = 0.01;

	Bimage*			pctf = img_ctf_calculate(em_ctf, action, wiener, p->size(), 
						p->sampling(0), lores, hires);
	pctf->invert();
	
//	write_img("pctf.pif", pctf, 0);
//	bexit(-1);
	
	int				back_transform(0);
	if ( p->fourier_type() == NoTransform ) {
		p->fft();
		back_transform = 1;
	}

//	write_img("p.pif", p);
	
	long 			i, j, nn;
	
	for ( nn=j=0; nn<p->images(); nn++ )
		for ( i=0; i<p->image_size(); i++, j++ )
			p->set(j, p->complex(j) * (*pctf)[i]);

//	write_img("pf.pif", p);
	
	delete pctf;
	
	if ( back_transform ) p->fft(FFTW_BACKWARD);

//	p->statistics();	// fft includes statistics
	
	return 0;
}

int 		img_ctf_apply(Bimage* p, CTFparam em_ctf, int action, double wiener,
				double lores, double hires, fft_plan planf_2D, fft_plan planb_2D)
{
	if ( action < 1 || action > 10 ) return 0;

	if ( wiener < 0.01 ) wiener = 0.01;

	Bimage*			pctf = img_ctf_calculate(em_ctf, action, wiener, p->size(), 
						p->sampling(0), lores, hires);
	pctf->invert();
	
	int				back_transform(0);
	if ( p->fourier_type() == NoTransform ) {
		p->fft(planf_2D, 1);
		back_transform = 1;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_ctf_apply: image transformed" << endl;
	}

	long 			i, j, nn;
	
	for ( nn=j=0; nn<p->images(); nn++ )
		for ( i=0; i<p->image_size(); i++, j++ )
			p->set(j, p->complex(j) * (*pctf)[i]);
	
	delete pctf;
	
	if ( back_transform ) {
		p->fft(planb_2D, 1);
		p->complex_to_real();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_ctf_apply: image back-transformed" << endl;
	}

	p->statistics();
	
	return 0;
}

int			img_ttf_apply_one(Bimage* p, long nn, CTFparam ctf, int action, 
				double wiener, double def, double res_lo, double res_hi, 
				Vector3<long> psize, fft_plan planf, fft_plan planb)
{
	verbose = 1;
	
	ctf.defocus_average(def);
	
	Bimage*			pt = p->extract(nn);
	
	pt->calculate_background();
	
	Vector3<long>	oldsize(pt->size()), translate;

	pt->pad(psize);
	
	img_ctf_apply(pt, ctf, action, wiener, res_lo, res_hi, planf, planb);
	
	pt->resize(oldsize, translate);
	
	p->replace(nn, pt);
	
	delete pt;

	return 0;
}

/**
@brief 	Applies a CTF to a projection image.
@param 	*proj			projection image.
@param	em_ctf			CTF parameters.
@param 	defocus			defocus.
@param 	res_lo			low resolution limit (angstrom).
@param 	res_hi			high resolution limit (angstrom).
@param	planf_2D       	2D forward fourier transform plan.
@param 	planb_2D		2D backward fourier transform plan.
@return int				error code.

**/
int			img_ctf_apply_to_proj(Bimage* proj, CTFparam em_ctf, double defocus, double res_lo, double res_hi, fft_plan planf_2D, fft_plan planb_2D)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_ctf_apply_to_proj: def_avg=" << em_ctf.defocus_average() << endl;

	if ( defocus > 0 ) em_ctf.defocus_average(defocus);
	
	img_ctf_apply(proj, em_ctf, 2, 0, res_lo, res_hi, planf_2D, planb_2D);

	return 0;
}

/**
@brief 	Applies or corrects for the tilted contrast transfer function (TTF).
@param 	*p				image (modified).
@param 	ctf				CTF parameter structure.
@param 	action			action to be taken.
@param 	wiener			Wiener factor (fraction).
@param 	tile_size		tile size for tilted CTF operations.
@param 	tilt			tilt angle (radians).
@param 	axis			tilt axis angle (radians).
@param 	res_lo			high resolution limit.
@param 	res_hi			low resolution limit.
@return int				0, <0 on error.

	The actions for this funtion are:
	1	flip phase (multiply by the sign of the CTF)
	2	apply a CTF (multiply with the CTF)
	3	correct for the CTF: ctf/(ctf^2 + wiener_factor)
	4	correct for the CTF: env*ctf/((env*ctf)^2 + noise^2)
	5	correct for the CTF with baseline: ctf/(ctf^2*noise^2 + wiener_factor)
	6	correct for the CTF with baseline: 1/(ctf*noise + sign*wiener_factor)

**/
int			img_ttf_apply(Bimage* p, CTFparam ctf, int action, double wiener,
				Vector3<long> tile_size, double tilt, double axis, double res_lo, double res_hi)
{
	if ( res_lo > 0 && res_lo < res_hi ) swap(res_lo, res_hi);
	
	// Minimum and maximum tile sizes allowed
	long			min_size(64);
	long			max_size(2048);
	Vector3<long> 	start;
	Vector3<long> 	ext_size;

	tile_size = tile_size.min(p->size());
	
	tile_size = tile_size.max(min_size);
	tile_size = tile_size.min(max_size);

	tile_size = tile_size.min(p->size());
	
	p->change_type(Float);
	
//	if ( verbose )
//		cout << "Extracting tiles" << endl;
	Bimage* 		pt = p->extract_tiles(0, tile_size, 0.5);

	long			nn;
	double			ca(cos(axis)), sa(sin(axis)), tt(tan(tilt));
	Vector3<double>	coor;
	double*			d = new double[pt->images()];
	
//	if ( verbose )
//		cout << "Tile\tx\ty\t∆f" << endl;
	for ( nn=0; nn<pt->images(); nn++ ) {
		coor = (pt->image[nn].origin() + pt->size()/2 - p->image->origin()) * pt->sampling(nn);	// Center coordinates
//		d[nn] = ctf.defocus_average() + (coor[0]*sa - coor[1]*ca)*tt;
		d[nn] = ctf.defocus_average() + (coor[1]*ca - coor[0]*sa)*tt;
//		if ( verbose )
//			cout << nn << tab << pt->image[nn].origin() << tab << coor << tab << d[nn] << endl;
	}

	Bstring			base_eq = ctf.baseline_equation();
	Bstring			env_eq = ctf.envelope_equation();

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Applying a tilted image transfer function:" << endl;
		switch ( action ) {
			case 1: cout << "Flipping phases" << endl; break;
			case 2: cout << "Applying a CTF" << endl; break;
			case 3: cout << "Correcting for the CTF" << endl; break;
			case 4: cout << "Correcting for the CTF with an optimal Wiener filter" << endl; break;
			case 5: cout << "Correcting for the CTF with baseline enhancement" << endl; break;
			case 6: cout << "Correcting for the CTF with baseline2 enhancement" << endl; break;
			case 7: cout << "Flipping phases with baseline enhancement" << endl; break;
			case 8: cout << "Correcting for the CTF with baseline3 enhancement" << endl; break;
			default: break;
		}
	}
	if ( verbose & VERB_PROCESS ) {
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Defocus range:                  " << d[0] << " - " << d[nn-1] << " (" << nn << ")" << endl;
		cout << "Tilt angle and axis:            " << tilt*180.0/M_PI << " " << axis*180.0/M_PI << endl;
		ctf.show();
		cout << "Resolution range:               " << res_hi << " - ";
		if ( res_lo > 0 ) cout << res_lo << " A" << endl;
		else cout << "inf A" << endl;
		if ( action == 3 || action == 5 || action == 6 )
			cout << "Wiener factor:                  " << wiener << endl;
		if ( action >= 4 ) {
			cout << "Envelope:                       " << env_eq << endl;
			cout << "Noise:                          " << base_eq << endl;
		}
		cout << endl;
	}
	
	long			pad_factor(1);
	long			ft_size = findNextPowerOf((int)(pad_factor*pt->sizeX()), 2);
	Vector3<long>	psize(ft_size, ft_size, 1);
	
	fft_plan		planf = fft_setup_plan(psize, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(psize, FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(pt->images(), dispatch_get_global_queue(0, 0), ^(size_t i){
		img_ttf_apply_one(pt, i, ctf, action, wiener, d[i], res_lo, res_hi, psize, planf, planb);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<pt->images(); i++ )
		img_ttf_apply_one(pt, i, ctf, action, wiener, d[i], res_lo, res_hi, psize, planf, planb);
#endif

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);
	
	p->assemble_tiles(pt, 1);

	delete pt;
	delete[] d;

	return 0;
}


int			mg_ps_name(Bmicrograph* mg, int img_num, Bstring& path, Bstring& newname, Bstring insert)
{
	Bstring			fn, ext;

	mg->fps = newname;
	
	if ( mg->fps.length() < 3 ) {
		if ( img_num >= 0 ) insert = Bstring(img_num, "_%03d") + insert;
		fn = mg->fmg;
		if ( fn.length() < 1 ) fn = mg->fframe;
		if ( fn.length() < 1 ) fn = mg->fpart;
		if ( fn.length() ) {
			ext = fn.extension();
			if ( ( ext == "tif" ) || ( ext == "jpg" ) || ( ext == "jpeg" )
				|| ( ext == "png" ) || ext.contains("dm") || ext.contains("eer") )
							ext = "mrc";
			mg->fps = fn.base() + insert + ext;
		} else {
			mg->fps = "ps" + insert + ext;
		}
	}
	
	if ( path.length() > 1 )
		mg->fps = path + mg->fps.post_rev('/');
	
	return 0;
}


Bimage*		mg_ctf_prepare(Bmicrograph* mg, int action, double lores, double hires,
				Vector3<long> tile_size, double def_start, double def_end, double def_inc, int flags)
{
	int				img_type(0), img_num(-1);
	int				ps_flags(7);
	Vector3<double>	pixel_size(mg->pixel_size);
		
	if ( !mg->ctf ) mg->ctf = new CTFparam;
	
	Bstring			filename, insert, ext;
	if ( action == 12 && mg->fps.length() ) {
		filename = mg->fps;
		img_type = 3;
		img_num = mg->img_num;
		if ( verbose & VERB_PROCESS )
			cout << "Processing power spectrum " << mg->fps << " (" << img_num << ")" << endl;
	} else if ( flags & 1 ) {
		if ( mg->fframe.length() && ( flags & 16 ) ) {
			filename = mg->fframe;
			ps_flags |= 16;
		} else {
			filename = mg->fmg;
			img_num = mg->img_num;
		}
		img_type = 1;
		if ( verbose & VERB_PROCESS )
			cout << "Processing micrograph " << filename << " (" << img_num << ")" << endl;
	} else if ( mg->fpart.length() ) {
		filename = mg->fpart;
		img_type = 2;
		img_num = -1;
		ps_flags |= 16;
		if ( mg->part->pixel_size[0] > 0 )
			pixel_size = mg->part->pixel_size;
		if ( verbose & VERB_PROCESS )
			cout << "Processing particles " << filename << endl;
	} else if ( mg->part && mg->part->fpart.length() ) {
		img_type = 4;
		ps_flags |= 16;
		if ( verbose & VERB_PROCESS )
			cout << "Processing particles from micrograph " << mg->id << endl;
	}
	
	if ( img_type < 4 && !filename.length() ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: no file specified for micrograph " << mg->id << endl;
		mg->select = 0;
		return NULL;
	}
	
	if ( verbose && filename.length() )
		cout << "Processing file " << filename << " (" << img_num << ")" << endl;
		
	if ( action == 11 || action == 13 ) insert = "_ps.";
	else insert = "_ctf.";

	Bimage*			p = NULL;
	if ( filename.length() ) {
		p = read_img(filename, 1, img_num);
		if ( !p ) {
			cerr << "Warning: File " << filename << " not found or read!" << endl;
			mg->select = 0;
			return NULL;
		}
		p->sampling(pixel_size);
		p->change_type(Float);
		if ( ( flags & 16 ) && p->sizeZ() > 1 ) p->slices_to_images();
	}
	
	if ( !p ) return NULL;
	
	if ( img_type > 1 ) tile_size = p->size();

	Bimage*			ps = NULL;
	
	if ( flags & 2 )
		p->filter_extremes();
	
//	cout << "volt=" << mg->ctf->volt() << " Cs=" << mg->ctf->Cs() << endl;
//	cout << "p sampling = " << p->image->sampling() << endl;
	
	if ( action == 11 || action == 13 ) {
		double		lamda(electron_wavelength(mg->ctf->volt()));
		double		iCL2(1.0/(mg->ctf->Cs()*lamda*lamda));
		ps = p->powerspectrum_tiled_and_tilted(tile_size, mg->tilt_axis, 
				mg->tilt_angle, 0, mg->ctf->defocus_average(), iCL2, ps_flags);
//		ps = p->powerspectrum_tiled(0, tile_size, ps_flags);
//		cout << "size = " << tile_size << endl;
		delete p;
		p = ps;
	}
	
//	cout << "ps sampling = " << ps->image->sampling() << endl;
	
	if ( action == 12 || action == 13 )
		mg->wri = img_ctf_fit(p, 0, *mg->ctf, lores, hires, def_start, def_end, def_inc, (flags & 8));

	return p;
}

int 		rec_ctf_prepare(Breconstruction* rec, int action, double lores, double hires,
				Vector3<long> tile_size, double def_start, double def_end, double def_inc,
				Bstring& newname, int flags)
{
	int				img_type(0);
	
	if ( !rec->ctf ) rec->ctf = new CTFparam;
	
	if ( verbose & VERB_PROCESS )
		cout << "Processing reconstruction " << rec->id << " (" << rec->frec << ")" << endl;
	
	Bstring			filename, insert, ext;
	if ( action == 12 && rec->fps.length() ) {
		filename = rec->fps;
		img_type = 3;
	} else if ( ( flags & 1 ) && rec->frec.length() ) {
		filename = rec->frec;
		img_type = 1;
	} else if ( rec->fpart.length() ) {
		filename = rec->fpart;
		img_type = 2;
	} else if ( rec->part && rec->part->fpart.length() ) {
		img_type = 4;
	}
	
	if ( img_type < 4 && !filename.length() ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: no file specified for reconstruction " << rec->id << endl;
		rec->select = 0;
		return -1;
	}
	
	if ( verbose && filename.length() )
		cout << "Processing file " << filename << endl;
		
	if ( action == 11 || action == 13 ) insert = "_ps.";
	else insert = "_ctf.";

	Bimage*			p = NULL;
	if ( filename.length() ) {
		p = read_img(filename, 1, -1);
		if ( !p ) {
			cerr << "Warning: File " << filename << " not found or read!" << endl;
			rec->select = 0;
			return -1;
		}
		p->sampling(rec->voxel_size);
		p->change_type(Float);
		if ( verbose & VERB_LABEL )
			cout << "Processing image " << p->file_name() << endl;
	}

	Bimage*			pex = NULL;
	
	if ( p && ( flags & 2 ) )
		p->filter_extremes();
	
	if ( p && ( action == 11 || action == 13 ) ) {
		pex = p->powerspectrum_tiled(0, tile_size, 7);
		delete p;
		p = pex;
		if ( newname.length() )
			rec->fps = newname;
		else if ( rec->fps.length() < 3 ) {
			ext = filename.extension();
			if ( ( ext == "tif" ) || ( ext == "jpg" ) || ( ext == "jpeg" )
				|| ( ext == "png" ) || ext.contains("dm") )
							ext = "pif";
			rec->fps = filename.pre_rev('.') + insert + ext;
		}
		write_img(rec->fps, p, 0);
	}
	
	if ( action == 12 || action == 13 )
		rec->fom = img_ctf_fit(p, 0, *rec->ctf, lores, hires, def_start, def_end, def_inc, (flags & 8));

	delete p;

	return 0;
}

/**
@brief 	Calculates power spectra and optionally fits CTF curves.
@param 	*project	project parameter structure.
@param 	action		CTF processing action.
@param 	lores		low resolution limit for CTF operations.
@param 	hires		high resolution limit for CTF operations.
@param 	tile_size	tile size for power spectrum generation.
@param	def_start	defocus search start (default 1e3).
@param	def_end		defocus search end (default 2e5).
@param	def_inc		defocus search increment (default 1e3).
@param 	&path		new power spectrum directory for output.
@param 	&newname	new file name for output.
@param 	flags		1=use mg or rec, 2=filter, 4=background, 8=astigmatism, 16=frames
@return int			0, <0 on error.

	The default is to use the particle file. If the particle file is not
	specified, the micrograph is used.

**/
int 		project_ctf_prepare(Bproject* project, int action, double lores,
				double hires, Vector3<long> tile_size,
				double def_start, double def_end, double def_inc,
				Bstring& path, Bstring& newname, int flags)
{
	if ( action < 1 ) return 0;
	
	double			ti = timer_start();

	if ( lores < 0 ) lores = 0;
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	
	if ( path.length() > 1 ) {
		mkdir(path.c_str(), (mode_t)0755);
		if ( path[-1] != '/' ) path += "/";
	}
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	CTFparam*			em_ctf = project->field->mg->ctf;
	
	if ( verbose ) {
		if ( project->select < 1 )
			cout << "CTF operations on all micrographs:" << endl;
		else
			cout << "CTF operations on all reconstructions:" << endl;
		switch ( action ) {
			case 11: cout << "Calculating power spectra" << endl; break;
			case 12: cout << "Fitting power spectra" << endl; break;
			case 13: cout << "Calculating and fitting power spectra" << endl; break;
			default: break;
		}
		if ( action == 12 || action == 13 ) {
			em_ctf->show();
			cout << "Baseline type:                  " << em_ctf->baseline_type() << endl;
			cout << "Envelope type:                  " << em_ctf->envelope_type() << endl;
			cout << "Defocus min, max, inc:          " << def_start << " - " << def_end << " ∆ " << def_inc << endl << endl;
		}
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Flags:                          " << flags << endl << endl;
	}
	
	long				nmg, onefile(0);
	Bimage*				ps = NULL;;
	Bimage*				ps1;
	
	Bstring				filename, insert, ext;
	if ( action == 11 || action == 13 ) insert = "_ps.";
	else insert = "_ctf.";

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( nmg = 1, mg = field->mg; mg->next; mg = mg->next ) nmg++;
			if ( ( action == 11 || action == 13 ) && ( field->mg->fmg == mg->fmg ) ) {
//				ps = new Bimage(Float, TSimple, tile_size, nmg);
//				ps->sampling(field->mg->pixel_size);
				onefile = 1;
			} else onefile = 0;
			for ( mg = field->mg; mg; mg = mg->next ) {
				ps1 = mg_ctf_prepare(mg, action, lores, hires, tile_size, def_start, def_end, def_inc, flags);
				if ( ps1 ) { 
					if ( action == 11 || action == 13 ) {
						if ( onefile == 1 ) {
							ps = new Bimage(Float, TSimple, ps1->size(), nmg);
							onefile = 2;
						}
						if ( onefile == 2 ) {
							ps->replace(mg->img_num, ps1);
							ps->sampling(ps1->sampling(0));
							mg_ps_name(mg, -1, path, newname, insert);
						} else {
//							mg_ps_name(mg, mg->img_num, path, newname, insert);
							mg_ps_name(mg, -1, path, newname, insert);
							write_img(mg->fps, ps1, 0);
						}
					}
					delete ps1;
				} else {
					cerr << "Warning: No power spectrum calculated for micrograph "
						<< mg->id << endl;
				}
			}
			if ( ps ) {
				write_img(field->mg->fps, ps, 0);
				delete ps;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			rec_ctf_prepare(rec, action, lores, hires, tile_size,
				def_start, def_end, def_inc, newname, flags);
		}
	}

	timer_report(ti);

	return 0;
}


int 		part_ctf(Bparticle* partlist, int action, double lores, double hires,
				double wiener, DataType datatype, Bstring& partpath,
				Bstring& newname, int flags)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_ctf: partlist=" << partlist << endl;
	
	Bparticle*			part = partlist;
	if ( !part ) return -1;
	
	Bmicrograph*		mg = part->mg;
	Breconstruction*	rec = part->rec;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_ctf: mg=" << mg << " rec=" << rec << endl;
	
	int					fn_flag(0);	// particle file flag: 0=mg, 1=part
	long				npart(0);
	Vector3<double>		pixel_size;
	Bstring				filename, insert("_ctf."), ext;
	Bimage*				p = NULL;
	Bimage*				pnu = NULL;
	
	CTFparam*			ctf = NULL;
	if ( mg ) {
		filename = mg->fpart;
		pixel_size = mg->pixel_size;
		if ( mg->ctf ) ctf = mg->ctf;
		else ctf = new CTFparam;
	} else if ( rec ) {
		filename = rec->fpart;
		pixel_size = rec->voxel_size;
		if ( rec->ctf ) ctf = rec->ctf;
		else ctf = new CTFparam;
	}

	for ( part = partlist; part; part = part->next ) npart++;
	part = partlist;
	
	if ( part->fpart.length() ) {
		filename = part->fpart;
		fn_flag = 1;
	}
	
	if ( !filename.length() ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: no particle file specified!" << endl;
		return -1;
	}

	p = read_img(filename, 0, -1);
	if ( !p ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: particle file " << filename << " not read!" << endl;
		return -1;
	}

	if ( verbose )
		cout << p->file_name() << " (" << p->images() << ")" << endl;
	
	if ( pixel_size[0] < 0.01 ) pixel_size = p->sampling(0);
	
	DataType		nudatatype = datatype;
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_ctf: p->images()=" << p->images() << " npart=" << npart << endl;
	
	if ( fn_flag < 1 || p->images() > 1 ) {
		if ( p->images() != npart ) {
			cerr << "Error: Micrograph " << mg->id << " has " << npart << " particles," << endl;
			cerr << "    but particle image file " << p->file_name() << " has " << p->images() << " images!" << endl;
			delete p;
//			return -2;
			bexit(-1);
		}
		pnu = p->copy_header();
		pnu->data_type(nudatatype);
		pnu->data_alloc();
		pnu->sampling(pixel_size);
//		pnu->change_type(nudatatype);
	}
	
	delete p;

	CTFparam			part_ctf;
	part_ctf.update(ctf);
	
	double			cosax(1), sinax(0), sintilt(0);

	if ( mg ) {
		cosax = cos(mg->tilt_axis);
		sinax = sin(mg->tilt_axis);
		sintilt = sin(mg->tilt_angle);
	}
	
	for ( part = partlist; part; part = part->next ) {
		if ( part->def > 0 ) part_ctf.defocus_average(part->def);
		else if ( mg )
			part_ctf.defocus_average(mg->ctf->defocus_average() +
				fabs(cosax*part->loc[0] - sinax*part->loc[1])*sintilt);
		if ( part->dev ) part_ctf.defocus_deviation(part->dev);
		if ( part->ast ) part_ctf.astigmatism_angle(part->ast);
		part->pixel_size = pixel_size;
		
		if ( part->fpart.length() ) filename = part->fpart;
		else if ( mg ) filename = mg->fpart;
		else if ( rec ) filename = rec->fpart;
		if ( part->fpart.length() ) p = read_img(filename, 1, 0);
		else p = read_img(filename, 1, part->id - 1);
		p->sampling(pixel_size);
		p->change_type(Float);

		if ( verbose & VERB_LABEL )
			cout << "Processing image " << p->file_name() << " (" << part->id << ")" << endl;
		
		img_ctf_apply(p, part_ctf, action, wiener, lores, hires);
		
		if ( flags & 4 ) p->correct_background();
		p->change_type(nudatatype);
		
//		if ( nimg == 1 ) insert = Bstring(part->id, "_%04d.");
		if ( newname.length() ) filename = newname;
		else filename = filename.pre_rev('.') + insert + filename.post_rev('.');

		if ( partpath.length() > 1 ) {
			if ( filename.contains("/") ) filename = filename.post_rev('/');
			filename = partpath + filename;
		}
		
		if ( part->fpart.length() ) {
			write_img(filename, p, 0);
			part->fpart = filename;
		} else {
			pnu->replace(part->id - 1, p);
		}

		delete p;
		
		if ( fabs(sintilt) > 1e-6 ) part->def = part_ctf.defocus_average();
	}
	
	if ( pnu ) {
		write_img(filename, pnu, 0);
		if ( mg ) mg->fpart = filename;
		else if ( rec ) rec->fpart = filename;
	}

	delete pnu;

	return 0;
}

int 		mg_ctf(Bmicrograph* mg, Bimage* pmg, int action, double lores, double hires, double wiener,
				Vector3<long> tile_size, DataType datatype, Bstring& newname, int flags)
{
	if ( !mg->ctf ) mg->ctf = new CTFparam;
	
	Bstring			filename(mg->fmg), insert("_ctf."), ext("pif");
	
	if ( !filename.length() ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: no file specified for micrograph " << mg->id << endl;
		mg->select = 0;
		return -1;
	}
	
	Bimage*			p = NULL;
	if ( !pmg ) p = read_img(filename, 1, mg->img_num);
	else p = pmg->extract(mg->img_num);
	
	if ( !p ) {
		cerr << "Warning: File " << filename << " not found or read!" << endl;
		mg->select = 0;
		return -1;
	}
	
	p->origin(mg->origin);
	if ( p->image->origin().length() < 1 ) p->origin(p->size()/2);
	p->sampling(mg->pixel_size);
	p->change_type(Float);
	
	DataType		nudatatype = datatype;
	if ( pmg ) nudatatype = pmg->data_type();
	else if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	
	if ( p && ( flags & 2 ) )
		p->filter_extremes();
	
	if ( mg->tilt_angle ) {
		if ( verbose )
			cout << p->file_name() << " (" << mg->img_num << ") " << mg->tilt_angle*180.0/M_PI << endl;
		img_ttf_apply(p, *(mg->ctf), action, wiener,
				tile_size, mg->tilt_angle, mg->tilt_axis, lores, hires);
	} else {
		if ( verbose )
			cout << p->file_name() << " (" << mg->img_num << ")" << endl;
		img_ctf_apply(p, *(mg->ctf), action, wiener, lores, hires);
	}
	
	
	if ( flags & 4 ) p->correct_background();
	p->change_type(nudatatype);
	
	if ( pmg ) {
		pmg->replace(mg->img_num, p);
	} else {
		if ( newname.length() ) filename = newname;
		else {
			insert = Bstring(mg->img_num, "_%03d_ctf.");
			if ( filename.contains(".dm") )
				filename = filename.pre_rev('.') + insert + ext;
			else
				filename = filename.pre_rev('.') + insert + filename.post_rev('.');
		}
		mg->fmg = filename;
		mg->img_num = 0;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG mg_ctf: Writing file " << filename << endl;
		write_img(filename, p, 0);
	}

	delete p;

	return 0;
}

int 		rec_ctf(Breconstruction* rec, int action, double lores, double hires, double wiener,
				DataType datatype, Bstring& newname, int flags)
{
	if ( !rec->ctf ) rec->ctf = new CTFparam;
	
	Bstring			filename(rec->frec), insert("_ctf."), ext;
	
	if ( !filename.length() ) {
		if ( verbose & VERB_FULL )
			cerr << "Warning: no file specified for reconstruction " << rec->id << endl;
		rec->select = 0;
		return -1;
	}
	
	Bimage*			p = read_img(filename, 1, 0);
	if ( !p ) {
		cerr << "Warning: File " << filename << " not found or read!" << endl;
		rec->select = 0;
		return -1;
	}
	
	p->origin(rec->origin);
	p->sampling(rec->voxel_size);
	p->change_type(Float);

	if ( verbose )
		cout << p->file_name() << endl;
	
	DataType		nudatatype = datatype;
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	
	if ( flags & 2 )
		p->filter_extremes();
	
	img_ctf_apply(p, *(rec->ctf), action, wiener, lores, hires);

	if ( newname.length() ) filename = newname;
	else filename = filename.pre_rev('.') + insert + filename.post_rev('.');
	
	rec->frec = filename;
	
	if ( p ) {
		if ( flags & 4 ) p->correct_background();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG rec_ctf: Writing file " << filename << endl;
		p->change_type(nudatatype);
		write_img(filename, p, 0);
	}

	delete p;

	return 0;
}

/**
@brief 	Calculates or fits CTF curves to multiple power spectra.
@param 	*project		project parameter structure.
@param 	action			CTF processing action.
@param 	lores			low resolution limit for CTF operations.
@param 	hires			high resolution limit for CTF operations.
@param 	tile_size		tile size for tilted CTF operations.
@param 	wiener			Wiener factor.
@param 	datatype		corrected particle file data type.
@param 	&partpath		corrected particle file path.
@param 	&newname		new file name for output.
@param 	flags			1=use mg or rec, 2=filter, 4=background, 8=astigmatism, 16=use frmaes
@return int				0, <0 on error.

	The default is to use the particle file. If the particle file is not
	specified, the micrograph is used. The selection can also be done with
	the use_mg flag.

**/
int 		project_ctf(Bproject* project, int action, double lores, 
				double hires, Vector3<long> tile_size, double wiener, 
				DataType datatype, Bstring& partpath, Bstring& newname, int flags)
{
	if ( action < 1 ) return 0;
	
	if ( lores < 0 ) lores = 0;
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bimage*				p = NULL;
	Bstring				filename, insert("_ctf."), ext("pif");
	
	if ( verbose ) {
		if ( project->select < 1 )
			cout << "CTF operations on all micrographs:" << endl;
		else
			cout << "CTF operations on all reconstructions:" << endl;
		switch ( action ) {
			case 1: cout << "Flipping phases" << endl; break;
			case 2: cout << "Applying a CTF" << endl; break;
			case 3: cout << "Correcting for the CTF" << endl; break;
			case 4: cout << "Correcting for the CTF with an optimal Wiener filter" << endl; break;
			case 5: cout << "Correcting for the CTF with baseline enhancement" << endl; break;
			case 6: cout << "Correcting for the CTF with baseline2 enhancement" << endl; break;
			case 7: cout << "Flipping phases with baseline enhancement" << endl; break;
			case 8: cout << "Correcting for the CTF with baseline3 enhancement" << endl; break;
			default: break;
		}
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Flags:                          " << flags << endl << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			if ( flags & 1 ) {
				for ( mg = field->mg; mg->next; mg = mg->next ) ;
				if ( field->mg->fmg == mg->fmg ) {
					if ( verbose )
						cout << "Reading:                        " << mg->fmg << " (" << mg->img_num + 1 << ")" << endl;
					p = read_img(mg->fmg, 1, -1);
				}
			}
			if ( verbose )
				cout << "Image\tNumber\tTilt" << endl;
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( flags & 1 )
					mg_ctf(mg, p, action, lores, hires, wiener, tile_size,
						datatype, newname, flags);
				else
					part_ctf(mg->part, action, lores, hires, wiener,
						datatype, partpath, newname, flags);
			}
			if ( p ) {
				filename = p->file_name();
				if ( newname.length() ) filename = newname;
				else {
					if ( filename.contains(".dm") )
						filename = filename.pre_rev('.') + insert + ext;
					else
						filename = filename.pre_rev('.') + insert + filename.post_rev('.');
				}
				for ( mg = field->mg; mg; mg = mg->next ) mg->fmg = filename;
				if ( verbose )
					cout << "Writing:                        " << filename << endl;
				write_img(filename, p, 0);
				delete p;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( flags & 1 )
				rec_ctf(rec, action, lores, hires, wiener,
					datatype, newname, flags);
			else
				part_ctf(rec->part, action, lores, hires, wiener,
						datatype, partpath, newname, flags);
		}
	}

	return 0;
}

/**
@brief 	Calculates the isotropy of the power spectrum adjusted for astigmatism.
@param 	*p			image structure.
@param	n			sub-image number.
@param 	&em_ctf		CTF parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@return double		radial average std/avg ratio.

	A power spectrum with its origin at the center.

**/
double		img_ctf_isotropy(Bimage* p, long n, CTFparam& em_ctf, double lores, double hires)
{
	if ( hires < 2*p->sampling(0)[0] ) hires = 2*p->sampling(0)[0];
	
	vector<double>	maxima = em_ctf.maxima(1.0/hires);
	
	long			i, na(10), nt(180/na);
	double			k, r, a, aa, da(M_PI/180), ca, sa, tavg, tstd, ravg(0);
	double			ox(p->image[n].origin()[0]), oy(p->image[n].origin()[1]);
	double			smin2 = 1-em_ctf.defocus_deviation()/em_ctf.defocus_average();
	double			smax2 = 1+em_ctf.defocus_deviation()/em_ctf.defocus_average();
	
	if ( verbose ) {
		cout << "Power spectrum isotropy:" << endl;
		cout << "s(1/A)\tAvg\tStd\tRatio" << endl;
	}
//	cout << em_ctf.defocus_average() << " ± " << em_ctf.defocus_deviation() << " @ " << em_ctf.astigmatism_angle()*180.0/M_PI << endl;
	for ( auto s: maxima ) {
		k = p->real_size()[0]*s;
//		cout << s << tab << 1/s << tab << k << endl;
		vector<double>		avg(nt,0);
		for ( i = 0, a = 0; i < 180; ++i, a += da ) {
			aa = a - em_ctf.astigmatism_angle();
			ca = cos(aa);
			sa = sin(aa);
			r = k*sqrt(smin2*ca*ca + smax2*sa*sa);
//			cout << k << endl;
			avg[i/na] += p->interpolate(r*cos(a)+ox, r*sin(a)+oy, 0, n);
//			p->set(p->index(r*cos(a)+ox, r*sin(a)+oy), -1);
		}
		tavg = tstd = 0;
		for ( i = 0, a = na*da/2; i < nt; ++i, a += na*da ) {
			avg[i] /= na;
			tavg += avg[i];
			tstd += avg[i]*avg[i];
//			cout << a*180/M_PI << tab << avg[i] << endl;
		}
		tavg /= nt;
		tstd /= nt;
		tstd -= tavg*tavg;
		if ( tstd > 0 ) tstd = sqrt(tstd);
		else tstd = 0;
		ravg += tstd/tavg;
		if ( verbose )
			cout << s << tab << tavg << tab << tstd << tab << tstd/tavg << endl;
	}
	
	if ( maxima.size() ) ravg /= maxima.size();
	
	if ( verbose )
		cout << "Average ratio:                           " << ravg << " (" << exp(-ravg) << ")" << endl << endl;
	
//	write_img("test.mrc", p, 0);
	
	return ravg;
}
/*
double		img_ctf_isotropy(Bimage* p, long n, double lores, double hires)
{
	if ( lores <= 0 || lores > p->real_size()[0] ) lores = p->real_size()[0];
	if ( lores < hires ) swap(lores, hires);
	if ( hires < 2*p->sampling(0)[0] ) hires = 2*p->sampling(0)[0];
	
	long			i, j, xx, yy;
	long			kmin(p->real_size()[0]/lores);
	long			kmax(p->real_size()[0]/hires);
	long			na(TWOPI*kmax);
	
	if ( verbose ) {
		cout << "Power spectrum anisotropy:" << endl;
		cout << "Origin:                         " << p->image[n].origin() << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Pixel radii:                    " << kmin << " - " << kmax << endl;
		cout << "Angles:                         " << na << endl;
	}
	
	double			dx, dy, r, a(0), b(0), f;
	vector<double>	x(na,0), y(na,0);
		
	for ( i=n*p->size().volume(), yy=0; yy<p->sizeY(); ++yy ) {
		dy = yy - p->image[n].origin()[1];
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
		x[i] = cos(2.0*i*TWOPI/na);
	}

	double			R = linear_least_squares(0, 179, x, y, a, b);
	
	if ( verbose ) {
		cout << "Fitted coefficients:            " << a << tab << b << endl;
		cout << "Anisotropy:                     " << b/a << endl;
		cout << "R:                              " << R << endl;
	}

	for ( i=0; i<na; ++i)
		cout << 360.0*i*1.0/na << tab << x[i] << tab << y[i] << tab << a + b*x[i] << endl;

	return b/a;
}
*/

double		sinc_find_argument(double v)
{
	double			a(0);
	
	for ( a=0.001; a<M_PI; a+=0.001 )
		if ( sin(a) <= v*a ) return a;
	
	return a;
}

/**
@brief 	Calculates the isotropy at the CTF maxima.
@param 	*project	project parameter structure.
@param 	lores		low resolution limit.
@param 	hires		high resolution limit
@return int			0.
**/
int			project_powerspectrum_isotropy(Bproject* project, double lores, double hires)
{
	if ( lores && lores < hires ) swap(lores, hires);
	
	Bfield*			field;
	Bmicrograph*	mg = NULL;
	Bframe*			frame = NULL;
	Bimage*			p = NULL;
	long			nf(1);
	double			aniso(0), ratio(0), residual(0);
	vector<double>	param;

	if ( verbose ) {
		cout << "Power spectrum anisotropy:" << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Micrograph\tAniso\tAngle\tResidual\tRes/frame" << endl;
	}
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			if ( mg->fps.length() && mg->ctf ) {
				p = read_img(mg->fps, 1, mg->img_num);
				if ( !p ) {
					cerr << "Warning: File " << mg->fps << " not found or read!" << endl;
					mg->select = 0;
				} else {
//					img_ctf_isotropy(p, 0, *mg->ctf, lores, hires);
					if ( mg->frame ) for ( nf=0, frame = mg->frame; frame; frame = frame->next ) nf++;
					param = p->powerspectrum_isotropy(0, lores, hires);
					aniso = param[1]/param[0];
					ratio = (param[0] - param[1])/(param[0] + param[1]);
//					residual = acos(sqrt(ratio))*(lores+hires)/TWOPI;
					residual = sinc_find_argument(sqrt(ratio))*(lores+hires)/TWOPI;
					if ( verbose )
						cout << mg->id << tab << aniso << tab << param[2]*180.0/M_PI << tab << residual << tab << residual/nf << endl;
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Determines the minimum and maximum defocus values.
@param 	*project		project parameter structure.
@return JSvalue			JSON object with the minimum and maximum values.
**/
JSvalue		project_defocus_range(Bproject* project)
{
	long			ndef(0);
	double			def_min(1e37), def_max(0), def_avg(0), pdef(0);
	Bfield*			field;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;

	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
				if ( part->def ) {
					if ( def_min > part->def )
						def_min = part->def;
					if ( def_max < part->def )
						def_max = part->def;
					def_avg += part->def;
					ndef++;
					pdef = 1;
				}
			}
			if ( !pdef && mg->ctf ) {
				if ( def_min > mg->ctf->defocus_average() )
					def_min = mg->ctf->defocus_average();
				if ( def_max < mg->ctf->defocus_average() )
					def_max = mg->ctf->defocus_average();
				def_avg += mg->ctf->defocus_average();
				ndef++;
			}
		}
	}
	
	if ( ndef ) def_avg /= ndef;
	
	JSvalue			js(JSobject);
	js["defocus_minimum"] = def_min;
	js["defocus_maximum"] = def_max;
	js["defocus_average"] = def_avg;
	
	if ( verbose & VERB_PROCESS )
		cout << js << endl;	
	
	return js;
}

/**
@brief 	Calculates or fits CTF curves to multiple power spectra.
@param 	*project		project parameter structure.
@param 	&psname			postscript file name for output.
@return int						0, <0 on error.

	The default is to use the particle file. If the particle file is not
	specified, the micrograph is used. The selection can also be done with
	the use_mg flag.

**/
int 		project_ctf_average(Bproject* project, Bstring& psname)
{
	Bfield*			field;
	Bmicrograph*	mg = NULL;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mg->box_size[0] ) break;
	
	if ( !mg ) {
		cerr << "Error: No micrograph with a box size specified!" << endl;
		return -1;
	}
	
	int				i, j, n(0), maxrad = mg->box_size[0]/2;
	double			recip_interval = 1.0/(2*maxrad*mg->pixel_size[0]);
	vector<double>	ctf;

	int				ncol = 2;
	Bstring			title("Compound CTF curves"), txt;
	Bplot*			plot = new Bplot(1, maxrad, ncol);
	plot->title(title);

	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("s");
	plot->page(0).column(1).label("CTFavg");
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).type(2);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(1).label("Resolution (A)");
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(x_max);
	plot->page(0).axis(3).label("Average CTF");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->ctf ) {
				ctf = mg->ctf->calculate(maxrad, 1, recip_interval);
				for ( i=0, j=maxrad; i<maxrad; i++, j++ )
					(*plot)[j] += ctf[i]*ctf[i];
//				delete ctf;
				n++;
			}
		}
	}
	
	for ( i=0, j=maxrad; i<maxrad; i++, j++ ) {
		(*plot)[i] = i*recip_interval;
		(*plot)[j] /= n;
	}

	txt = "Number of micrographs:   " + Bstring(n, "%d");
	plot->page(0).add_text(txt);

	ps_plot(psname, plot);
	
	delete plot;
	
	return 0;
}

/**
@brief 	Averages multiple power spectra based on defocus estimates.
@param 	*project		project parameter structure.
@param 	deftarget		target defocus (angstrom).
@return Bimage*			average power spectrum.

**/
Bimage*		project_powerspectrum_average(Bproject* project, double deftarget)
{
	long			n(0);
	Bimage*			ps = NULL;
	Bimage*			pt = NULL;
	Bimage*			psavg = NULL;
	Bfield*			field;
	Bmicrograph*	mg = NULL;
	Vector3<double>	origin, scale(1,1,1), translate;
	Matrix3			mat(1);

	if ( verbose ) {
		cout << "Calculating an average power spectrum with target defocus: " << deftarget << endl;
		cout << "Micrograph\tDefocus\tScale" << endl;
	}
		
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->fps.length() )
				ps = read_img(mg->fps, 1, mg->img_num);
			if ( ps ) {
				if ( !psavg ) {
					psavg = ps->copy(1);
					psavg->origin(ps->size()/2);
					psavg->fill(0);
					psavg->sampling(mg->pixel_size);
				}
				ps->calculate_background();
				scale[0] = scale[1] = sqrt(mg->ctf->defocus_average()/deftarget);
				if ( verbose )
					cout << mg->id << tab << mg->ctf->defocus_average() << tab << scale << endl;
				pt = ps->transform(0, ps->size(), scale, psavg->image->origin(), translate, mat, FILL_BACKGROUND);
				psavg->add(pt);
				delete pt;
				delete ps;
				n++;
			}
		}
	}

	if ( n ) psavg->multiply(1.0/n);
	
	return psavg;
}

/**
@brief 	Puts CTF parameters from one project into another.
@param 	*project		project parameter structure with all parameters.
@param 	*ctfproject		project parameter structure with CTF parameters.
@return int				0.
**/
int			project_merge_CTF_parameters(Bproject* project, Bproject* ctfproject)
{
	if ( !project || ! ctfproject ) return 0;
	
	Bfield				*field, *ctffield;
	Bmicrograph			*mg, *ctfmg;
	Breconstruction		*rec, *ctfrec;

	if ( verbose & VERB_PROCESS )
		cout << "Merging CTF parameters into the main project" << endl << endl;
	
	for ( field=project->field, ctffield=ctfproject->field; field; field=field->next, ctffield=ctffield->next ) {
		for ( mg=field->mg, ctfmg=ctffield->mg; mg && ctfmg; mg=mg->next, ctfmg=ctfmg->next ) {
			if ( ctfmg->fps.length() > 0 ) mg->fps = ctfmg->fps;
			if ( ctfmg->magnification ) mg->magnification = ctfmg->magnification;
			if ( ctfmg->sampling ) mg->sampling = ctfmg->sampling;
			if ( ctfmg->ctf ) {
				if ( !mg->ctf ) mg->ctf = new CTFparam;
				mg->ctf->update(ctfmg->ctf);
			}
		}
	}
	
	for ( rec=project->rec, ctfrec=ctfproject->rec; rec && ctfrec; rec=rec->next, ctfrec=ctfrec->next ) {
		if ( ctfrec->fps.length() > 0 ) rec->fps = ctfrec->fps;
		if ( ctfrec->ctf ) {
			if ( !rec->ctf ) rec->ctf = new CTFparam;
			rec->ctf->update(ctfrec->ctf);
		}
	}
	
	return 0;
}

/**
@brief 	Transfers CTF parameters from micrographs to particles.
@param 	*project		project parameter structure with all parameters.
@return int						0.
**/
int			project_CTF_to_part(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( verbose & VERB_PROCESS )
		cout << "Transferring CTF parameters to particles" << endl << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			for ( part = mg->part; part; part = part->next ) {
				part->def = mg->ctf->defocus_average();
				part->dev = mg->ctf->defocus_deviation();
				part->ast = mg->ctf->astigmatism_angle();
			}
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->ctf ) {
		for ( part = rec->part; part; part = part->next ) {
			part->def = rec->ctf->defocus_average();
			part->dev = rec->ctf->defocus_deviation();
			part->ast = rec->ctf->astigmatism_angle();
		}
	}
	
	return 0;
}

/**
@brief 	Sets the defocus values of all the micrographs.
@param 	*project 		project parameter structure.
@param 	def_avg			defocus average.
@param 	def_dev			defocus deviation.
@param 	ast_angle 		astigmatism angle.
@return int						0.
**/
int			project_set_defocus(Bproject* project, double def_avg,
					double def_dev, double ast_angle)
{
	if ( !project ) return 0;
	
	long				n(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_defocus: avg=" << def_avg << " dev=" << def_dev << " ang=" << ast_angle*180.0/M_PI << endl;
	
	for ( n=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->defocus_average(def_avg);
			if ( def_dev ) {
				mg->ctf->defocus_deviation(def_dev);
				mg->ctf->astigmatism_angle(ast_angle);
			}
			n++;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_defocus: micrographs set = " << n << endl;
	
	for ( n=0, rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->defocus_average(def_avg);
		if ( def_dev ) {
			rec->ctf->defocus_deviation(def_dev);
			rec->ctf->astigmatism_angle(ast_angle);
		}
		n++;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_defocus: reconstructions set = " << n << endl;
	
	return 0;
}

/**
@brief 	Sets the defocus deviation and astigmatism angle of all the micrographs.
@param 	*project 		project parameter structure.
@param 	def_dev			defocus deviation.
@param 	ast_angle 		astigmatism angle.
@return int						0.
**/
int			project_set_astigmatism(Bproject* project, double def_dev, double ast_angle)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_astigmatism: dev=" << def_dev << " ang=" << ast_angle*180.0/M_PI << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->defocus_deviation(def_dev);
			mg->ctf->astigmatism_angle(ast_angle);
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->defocus_deviation(def_dev);
		rec->ctf->astigmatism_angle(ast_angle);
	}
	
	return 0;
}

/**
@brief 	Sets the acceleration voltage of all the micrographs.
@param 	*project 		project parameter structure.
@param 	jsctf			JSON parameters to be updated.
@return int				0.
**/
int			project_update_ctf(Bproject* project, JSvalue& jsctf)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose & VERB_FULL )
		cout << "CTF parameters:" << endl << jsctf << endl;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
//			mg->ctf[0] = ctf_from_json(jsctf);
			ctf_update_from_json(mg->ctf[0], jsctf);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
//		rec->ctf[0] = ctf_from_json(jsctf);
		ctf_update_from_json(rec->ctf[0], jsctf);
	}
	
	return 0;
}

/**
@brief 	Sets the acceleration voltage of all the micrographs.
@param 	*project 		project parameter structure.
@param 	volts			acceleration voltage.
@return int				0.
**/
int			project_set_volts(Bproject* project, double volts)
{
	if ( !project ) return 0;
	
	if ( volts < 1e4 ) volts *= 1000;	// Asume kV
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->volt(volts);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->volt(volts);
	}
	
	return 0;
}

/**
@brief 	Sets the spherical aberation constant of all the micrographs.
@param 	*project 		project parameter structure.
@param 	Cs				spherical aberation constant in angstrom.
@return int				0.
**/
int			project_set_Cs(Bproject* project, double Cs)
{
	if ( !project ) return 0;
	if ( Cs < 10 ) Cs *= 1e7; // Assume mm

	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->Cs(Cs);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->Cs(Cs);
	}
	
	return 0;
}

/**
@brief 	Sets the amplitude contribution of all the micrographs.
@param 	*project 		project parameter structure.
@param 	amp_shift		amplitude contribution phase shift.
@return int				0.
**/
int			project_set_amp_shift(Bproject* project, double amp_shift)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->amp_shift(amp_shift);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->amp_shift(amp_shift);
	}
	
	return 0;
}

/**
@brief 	Sets the focal length of all the micrographs.
@param 	*project 		project parameter structure.
@param 	focal_length	focal length in angstrom.
@return int				0.
**/
int			project_set_focal_length(Bproject* project, double focal_length)
{
	if ( !project ) return 0;
	if ( focal_length < 10 ) focal_length *= 1e7; // Assume mm
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->focal_length(focal_length);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->focal_length(focal_length);
	}
	
	return 0;
}

/**
@brief 	Sets the objective aperture of all the micrographs.
@param 	*project 		project parameter structure.
@param 	aperture		objective aperture in angstrom.
@return int				0.
**/
int			project_set_aperture(Bproject* project, double aperture)
{
	if ( !project ) return 0;
	if ( aperture < 1000 ) aperture *= 1e4; // Assume µm
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->objective_aperture(aperture);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->objective_aperture(aperture);
	}
	
	return 0;
}

/**
@brief 	Sets the energy filter slit width of all the micrographs.
@param 	*project 		project parameter structure.
@param 	slit			slit width in eV.
@return int				0.
**/
int			project_set_slit_width(Bproject* project, double slit)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->slit_width(slit);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->slit_width(slit);
	}
	
	return 0;
}

/**
@brief 	Sets the beam source size/divergence angle of all the micrographs.
@param 	*project 		project parameter structure.
@param 	alpha				beam source size/divergence angle (radians).
@return int						0.
**/
int			project_set_alpha(Bproject* project, double alpha)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->alpha(alpha);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->alpha(alpha);
	}
	
	return 0;
}

/**
@brief 	Sets the envelope type of all the micrographs.
@param 	*project 		project parameter structure.
@param 	type			envelope type.
@return int				0.
**/
int			project_set_envelope_type(Bproject* project, int type)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->envelope_type(type);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->envelope_type(type);
	}
	
	return 0;
}

/**
@brief 	Sets the envelope equations of all the micrographs.
@param 	*project 		project parameter structure.
@param 	type			envelope type.
@param 	*coeff			5 envelope coefficients.
@return int				0.
**/
int			project_set_envelope(Bproject* project, int type, double* coeff)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->envelope_type(type);
			mg->ctf->envelope(coeff);
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->envelope_type(type);
		rec->ctf->envelope(coeff);
	}
	
	return 0;
}

/**
@brief 	Replaces envelope equations with those based on partial coherence in all micrographs.
@param 	*project 		project parameter structure.
@return int				0.

	Partial coherence envelope:
		env = amp*exp(-(pi*alpha*defocus*s)^2)
	The amplitude, defocus and alpha values are taken from the fields 
	in each micrograph. The defocus must already be determined.

**/
int			project_set_coherence_envelope(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->envelope(1, -M_PI*M_PI*mg->ctf->alpha()*mg->ctf->alpha()*mg->ctf->defocus_average()*mg->ctf->defocus_average());
		}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->envelope(1, -M_PI*M_PI*rec->ctf->alpha()*rec->ctf->alpha()*rec->ctf->defocus_average()*rec->ctf->defocus_average());
	}
	
	return 0;
}

/**
@brief 	Sets the baseline type of all the micrographs.
@param 	*project 		project parameter structure.
@param 	type			baseline type.
@return int				0.
**/
int			project_set_baseline_type(Bproject* project, int type)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->baseline_type(type);
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->baseline_type(type);
	}
	
	return 0;
}

/**
@brief 	Sets the baseline equations of all the micrographs.
@param 	*project 		project parameter structure.
@param 	type				baseline type.
@param 	*coeff			5 baseline coefficients.
@return int						0.
**/
int			project_set_baseline(Bproject* project, int type, double* coeff)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->baseline_type(type);
			mg->ctf->baseline(coeff);
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( !rec->ctf ) rec->ctf = new CTFparam;
		rec->ctf->baseline_type(type);
		rec->ctf->baseline(coeff);
	}
	
	return 0;
}

/**
@brief 	Updates the first zero from the defocus average for all the micrographs.
@param 	*project 		project parameter structure.
@return int						0.
**/
int			project_update_first_zero(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_update_first_zero: start" << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			if ( mg->ctf && mg->ctf->defocus_average() > 0 && mg->ctf->volt() > 0 && mg->ctf->Cs() > 0 )
					mg->ctf->zero(1);
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		if ( rec->ctf && rec->ctf->defocus_average() > 0 && rec->ctf->volt() > 0 && rec->ctf->Cs() > 0 )
				rec->ctf->zero(1);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_update_first_zero: done" << endl;

	return 0;
}

/**
@brief 	Plots the CTF curves.
@param 	*project 		project parameter structure.
@param 	&filename		Postscript file name.
@return int				0.
**/
int			project_plot_ctf(Bproject* project, Bstring& filename)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Bimage*				p;
	Bimage*				prad;
	Bstring				basename = filename.pre_rev('.');
	Bstring				psname;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			p = read_img(mg->fps, 1, 0);
			if ( p ) {
				psname = basename + "_" + mg->fps.pre_rev('.') + ".ps";
				prad = img_ctf_radial_average(p, 0, *mg->ctf);
				ps_ctf_plot(prad->sizeX(), (double *)prad->data_pointer(), 
					1.0/p->real_size()[0], mg->ctf, psname);
				delete prad;
			}
		}
	}
	
	return 0;
}

