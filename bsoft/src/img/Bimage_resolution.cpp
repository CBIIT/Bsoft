/**
@file	Bimage_resolution.cpp
@brief	Library routines to estimate resolution 
@author 	Bernard Heymann
@date	Created: 20000611
@date	Modified: 20220324
**/

#include "Bimage.h"
#include "utilities.h"

#define	LORESLIM	200

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bimage*		Bimage::resolution_prepare(Bimage* p)
{
	if ( fouriertype != NoTransform ) {
		cerr << "Error: File " << file_name() << " must be a real space map!" << endl;
		return NULL; 
	}
	
	if ( p->fourier_type() != NoTransform ) {
		cerr << "Error: File " << p->file_name() << " must be a real space map!" << endl;
		return NULL; 
	}

	if ( verbose > VERB_RESULT )
		cout << "Determining resolution from images " << file_name() << " and " << p->file_name() << endl;
	
	// Pack the two images into one complex block
	Bimage*			pc = pack_two_in_complex(p);
	if ( pc == NULL ) return pc;
	
	pc->fft();

	pc->label(file_name() + " vs " + p->file_name());
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::resolution_prepare: FFT done" << endl;
	
	return pc;
}

Bimage*		Bimage::resolution_prepare(Bimage* p, fft_plan plan)
{
	if ( fouriertype != NoTransform ) {
		cerr << "Error: File " << file_name() << " must be a real space map!" << endl;
		return NULL; 
	}
	
	if ( p->fourier_type() != NoTransform ) {
		cerr << "Error: File " << p->file_name() << " must be a real space map!" << endl;
		return NULL; 
	}
	
	if ( verbose > VERB_RESULT )
		cout << "Determining resolution from images " << file_name() << " and " << p->file_name() << endl;

	// Pack the two images into one complex block
	Bimage*			pc = pack_two_in_complex(p);
	if ( pc == NULL ) return pc;
	
	pc->fft(plan, 0);
	
	pc->label(file_name() + " vs " + p->file_name());
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::resolution_prepare: FFT done" << endl;
	
	return pc;
}


/**
@brief 	Calculates an FSC and DPR curves from two images.
@param 	hi_res			high resolution limit.
@param 	sampling_ratio	radial sampling ratio (1 for per voxel sampling).
@param	flag			if 1 calculate only the FSC.
@return Bplot*			FSC curve.

	FRC: Fourier ring/shell correlation 
	-----------------------------------
	Saxton & Baumeister (1982) J. Microscopy 127, 127-138 
	de la Fraga et al. (1995) Ultramicroscopy 60, 385-391 
	           sum(|F1|*|F2|) 
	FRC/FSC = --------------------------------- 
	          sqrt( sum(|F1|^2) * sum(|F2|^2) ) 

**/
Bplot*		Bimage::fsc_dpr(double hi_res, double sampling_ratio, int flag)
{
	check_resolution(hi_res);
	if ( sampling_ratio <= 0 ) sampling_ratio = 1;
	
	// The frequency scaling is linked to the different dimensions of the
	// data set (important when the x, y and z dimensions are different)
	// The radius scaling is set to a value consistent with one pixel width
	// in reciprocal space for the longest dimension.
	Vector3<double>	freq_scale(1/real_size());
	double			rad_scale = real_size()[0]/sampling_ratio;
//	if ( hi_res < 2*sampling(0)[0] ) hi_res = 2*sampling(0)[0];
//	if ( hi_res > real_size().max()/2.0 ) hi_res = real_size().max()/2.0;
	
	long			maxrad = (long) (2 + rad_scale/hi_res);

	Complex<double>	sf1, sf2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fsc_dpr: FFT done" << endl;
	
	long			i, j, xx, yy, zz, ix, iy, iz, iradius, iradius2; 
	double			radius, fraction, fraction2, I1, I2, Isum, F1F2re, dphi, Id2;
	double			rx, ry, rz;
	
	if ( verbose > VERB_RESULT ) {
		cout << "Resolution limit:               " << hi_res << " A (" << maxrad << ")" << endl;
		cout << "Sampling ratio:                 " << sampling_ratio << endl << endl;
	}
	
	vector<double>	F1(maxrad,0);
	vector<double>	F2(maxrad,0);
	vector<double>	Fsum(maxrad,0);
	vector<double>	FSC(maxrad,0);
	vector<double>	DPR(maxrad,0);
	
	for ( zz=0; zz<z; zz++ ) {
		rz = zz; 
		if ( rz > (z - 1)/2 ) rz -= z;
		rz *= freq_scale[2];
		iz = -zz;
		if ( iz < 0 ) iz += z;
		for ( yy=0; yy<y; yy++ ) { 
			ry = yy; 
			if ( ry > (y - 1)/2 ) ry -= y;
			ry *= freq_scale[1];
			iy = -yy;
			if ( iy < 0 ) iy += y;
			for ( xx=0; xx<x/2; xx++ ) { 
				i = (zz*y + yy)*x + xx;
				rx = xx; 
				if ( xx > (x - 1)/2 ) rx -= x;
				rx *= freq_scale[0];
				ix = -xx;
				if ( ix < 0 ) ix += x; 
				radius = rad_scale*sqrt(rx*rx + ry*ry + rz*rz); 
				iradius = (long) radius;
				iradius2 = iradius + 1;
				if ( iradius2 < maxrad ) {
					fraction = radius - iradius;
					fraction2 = 1.0 - fraction;
					j = (iz*y + iy)*x + ix;
					sf1 = complex(i).unpack_first(complex(j));
					sf2 = complex(i).unpack_second(complex(j));
					I1 = sf1.power(); 
					I2 = sf2.power();
					dphi = angle_set_negPI_to_PI(sf1.phi() - sf2.phi());
					dphi *= dphi;
					F1[iradius] += fraction2*I1;
					F1[iradius2] += fraction*I1; 
					F2[iradius] += fraction2*I2; 
					F2[iradius2] += fraction*I2; 
					F1F2re = sf1.real()*sf2.real() + sf1.imag()*sf2.imag();
					FSC[iradius] += fraction2*F1F2re;
					FSC[iradius2] += fraction*F1F2re;
					Isum = I1 + I2;
					Fsum[iradius] += fraction2*Isum;
					Fsum[iradius2] += fraction*Isum;
					Id2 = Isum*dphi;
					DPR[iradius] += fraction2*Id2;
					DPR[iradius2] += fraction*Id2;
				}
			} 
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fsc_dpr: Creating plot" << endl;

	int			ncol(3);
	if ( flag ) ncol = 2;
	
	Bstring		title("Resolution"), txt;
	Bplot*		plot = new Bplot(1, maxrad, ncol);
	
	plot->title(title);
	plot->page(0).title(label());
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("Spatial Frequency A");
	plot->page(0).column(0).axis(1);

	if ( z > 1 ) plot->page(0).column(1).label("FSC");
	else plot->page(0).column(1).label("FRC");
	plot->page(0).column(1).type(2);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).color(1,0,0);
	
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(1/hi_res);
	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(1).label("Resolution (A)");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).inc(0.1);
	if ( z > 1 ) plot->page(0).axis(3).label("FSC");
	else plot->page(0).axis(3).label("FRC");
	
	if ( ncol > 2 ) {
		plot->page(0).column(2).label("DPR");
		plot->page(0).column(2).type(2);
		plot->page(0).column(2).axis(4);
		plot->page(0).column(2).color(0,0,1);

		plot->page(0).axis(4).min(0);
		plot->page(0).axis(4).max(120);
		plot->page(0).axis(4).inc(10);
		plot->page(0).axis(4).label("DPR");
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fsc_dpr: transferring data, maxrad=" << maxrad << endl;

	(*plot)[0] = 0;
	(*plot)[maxrad] = 1; 
	for ( i=1; i<maxrad; i++ ) {
//		cout << i << tab << F1[i] << tab << F2[i] << tab << FSC[i] << tab << endl;
		if ( F1[i]*F2[i] > SMALLFLOAT ) FSC[i] /= sqrt(F1[i]*F2[i]);
//		if ( FSC[i] < SMALLFLOAT ) FSC[i] = 0;
		if ( FSC[i] < -1 ) FSC[i] = -1;
		if ( FSC[i] > 1 ) FSC[i] = 1;
		if ( Fsum[i] > SMALLFLOAT ) DPR[i] = sqrt(DPR[i]/Fsum[i]);
		else DPR[i] = 0;
		(*plot)[i] = i/rad_scale;
		(*plot)[maxrad+i] = FSC[i];
		if ( ncol > 2 ) (*plot)[2*maxrad+i] = DPR[i]*180.0/M_PI;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::fsc_dpr: Plot done" << endl;

	return plot;
}

Bplot*		Bimage::fsc(double hi_res, double sampling_ratio, vector<double>& fsccut)
{
	check_resolution(hi_res);

	if ( verbose ) {
		cout << "Determining resolution for " << n << " images" << endl;
		cout << "Resolution limit:               " << hi_res << endl;
		cout << "Sampling ratio:                 " << sampling_ratio << endl << endl;
	}
	
	long			i, j, nn, ncol(n+1);
	double			r;
//	vector<double> 	fsccut{0.143, 0.3, 0.5, 0.9};
	Bstring			title("Resolution"), txt;
	Bimage*			p1 = NULL;
	Bplot*			plot = NULL;
	Bplot*			plot1 = NULL;
	
	if ( verbose ) {
		cout << "Image";
		for ( i=0; i<4; i++ ) if ( fsccut[i] ) cout << tab << "FSC(" << fsccut[i] << ")";
		cout << endl;
	}
	
	for ( nn=0; nn<n; nn++ ) {

		p1 = extract(nn);
		plot1 = p1->fsc_dpr(hi_res, sampling_ratio);
		delete p1;
		
		if ( !plot ) {
			plot = new Bplot(n, plot1->rows(), ncol);
			plot->title(title);
			for ( i=0; i<plot->rows(); i++ ) 
				(*plot)[i] = (*plot1)[i];
		}

		for ( i=plot->rows(), j=(nn+1)*plot->rows(); i<2*plot->rows(); i++, j++ )
			(*plot)[j] = (*plot1)[i];
	
		delete plot1;

//		title = label() + Bstring(nn+1, " image %d");
		title = label() + " image " + to_string(nn+1);
		plot->page(nn).title(title);
		plot->page(nn).columns(2);
		plot->page(nn).column(0).number(0);
		plot->page(nn).column(1).number(nn+1);
		plot->page(nn).column(0).label("Spatial Frequency (A)");
		plot->page(nn).column(0).axis(1);
		if ( z > 1 ) plot->page(nn).column(1).label("FSC");
		else plot->page(nn).column(1).label("FRC");
		plot->page(nn).column(1).type(2);
		plot->page(nn).column(1).axis(3);
		plot->page(nn).column(1).color(0,0,0);
		plot->page(nn).axis(1).min(0);
		plot->page(nn).axis(1).max(1/hi_res);
		plot->page(nn).axis(1).inc(0.1/hi_res);
		plot->page(nn).axis(1).flags(1);
		plot->page(nn).axis(1).label("Resolution (A)");
		plot->page(nn).axis(3).min(0);
		plot->page(nn).axis(3).max(1);
		plot->page(nn).axis(3).inc(0.1);
		if ( z > 1 ) plot->page(nn).axis(3).label("FSC");
		else plot->page(nn).axis(3).label("FRC");

		if ( verbose )
			cout << nn+1;
		for ( i=0; i<4; i++ ) if ( fsccut[i] ) {
			r = 1/plot->cut(nn+1,fsccut[i],-1);
			txt = plot->page(nn).column(1).label() + Bstring(fsccut[i], "(%g): ") 
				+ Bstring(r, "%g A");
			plot->page(nn).add_text(txt);
			if ( verbose )
				cout << tab << r;
		}
		if ( verbose )
			cout << endl;
	}

	return plot; 
}

Bplot*		Bimage::fsc(Bimage* p, double hi_res, double sampling_ratio)
{
	check_resolution(hi_res);
	
	if ( verbose ) {
		cout << "Determining resolution for " << n << " images" << endl;
		cout << "Resolution limit:               " << hi_res << endl;
		cout << "Sampling ratio:                 " << sampling_ratio << endl << endl;
	}
	
	long			i, j, nn, ncol(n+1);
	double			rad_scale(real_size()[0]/sampling_ratio);
	Bstring			title("Resolution"), txt;
	Bimage*			p1 = NULL;
	Bimage*			p2 = NULL;
	Bimage*			pr1 = NULL;
	Bimage*			pr2 = NULL;
	Bplot*			plot = NULL;
		
	for ( nn=0; nn<n; nn++ ) {

		p1 = extract(nn);
		p2 = p->extract(nn);
		pr1 = p1->fspace_radial_power(hi_res, sampling_ratio);
		pr2 = p2->fspace_radial_power(hi_res, sampling_ratio);
		p1->complex_conjugate_product(p2);
		p1->complex_to_real();
		delete p2;
		p1->origin(0,0,0);
//		p2 = p1->radial(0, maxrad, 1, 1);
		p2 = p1->fspace_radial_power(hi_res, sampling_ratio);
		p2->set(0, 1);
		for ( i=1; i<p2->sizeX(); ++i ) p2->set(i, (*p2)[i] / sqrt((*pr1)[i] * (*pr2)[i]));
		delete p1;
		delete pr1;
		delete pr2;
		
		if ( !plot ) {
			plot = new Bplot(n, p2->sizeX(), ncol);
			plot->title(title);
			for ( i=0; i<plot->rows(); ++i ) 
				(*plot)[i] = i/rad_scale;
		}

		for ( i=0, j=(nn+1)*plot->rows(); i<plot->rows(); ++i, ++j )
			if ( isfinite((*p2)[i]) ) (*plot)[j] = (*p2)[i];
			else (*plot)[j] = 0;
	
		delete p2;

//		title = label() + Bstring(nn+1, " image %d");
		title = "Image: " + to_string(nn+1);
		plot->page(nn).title(title);
		plot->page(nn).columns(2);
		plot->page(nn).column(0).number(0);
		plot->page(nn).column(1).number(nn+1);
		plot->page(nn).column(0).label("Spatial Frequency (A)");
		plot->page(nn).column(0).axis(1);
		if ( z > 1 ) plot->page(nn).column(1).label("FSC");
		else plot->page(nn).column(1).label("FRC");
		plot->page(nn).column(1).type(2);
		plot->page(nn).column(1).axis(3);
		plot->page(nn).column(1).color(0,0,0);
		plot->page(nn).axis(1).min(0);
		plot->page(nn).axis(1).max(1.0/hi_res);
		plot->page(nn).axis(1).inc(0.1/hi_res);
		plot->page(nn).axis(1).flags(1);
		plot->page(nn).axis(1).label("Resolution (A)");
		plot->page(nn).axis(3).min(0);
		plot->page(nn).axis(3).max(1);
		plot->page(nn).axis(3).inc(0.1);
		if ( z > 1 ) plot->page(nn).axis(3).label("FSC");
		else plot->page(nn).axis(3).label("FRC");
	}

	return plot; 
}

/**
@brief 	Determine the resolution for each concentric shell in a map.
@param 	*p			second image.
@param 	hi_res		high resolution limit.
@param 	*cutoff		correlation threshold(s).
@param 	thickness	shell thickness.
@param 	step		step size between shells.
@param 	minrad		minimum radius.
@param 	maxrad		maximum radius.
@param 	pad			padding factor.
@param 	smooth		flag for edge smooting.
@param 	fill		background fill value.
@return Bimage* 	a 1D image containing the shell resolutions.
**/
Bimage*		Bimage::fsc_shell(Bimage* p, double hi_res, double* cutoff, 
						int thickness, int step, int minrad, int maxrad, 
						int pad, int smooth, double fill)
{
	long   			i, j, r;
	double			r1, r2;
	int 			fill_type(FILL_AVERAGE);	// Fill type for resizing
	double			sampling_ratio(1);			// Sampling ratio in reciprocal space
	
	Bimage*			cp1 = NULL;
	Bimage*			cp2 = NULL;
	Bplot*			plot;
	long			ft_size(x), nmap(0), imgsize(image_size());
	double			res_est;					// Resolution estimate
	
    if ( pad > 0 ) {
		pad = x*(pad+1);
		ft_size = pad;
	}
    
	if ( minrad < 0 ) minrad = 0;
	if ( maxrad > x/2 ) maxrad = x/2;
	if ( maxrad < 1 ) maxrad = x/2;
	
	// Count the number of output maps
	for ( i=0; i<4; ++i ) if ( cutoff[i] > 0 ) nmap++;
    
	if ( verbose ) {
		cout << "Calculating shell resolution:" << endl;
		cout << "Shell thickness:                " << thickness << endl;
		cout << "Step size:                      " << step << endl;
		cout << "Radii:                          " << minrad << " " << maxrad << endl;
		cout << "Padding size:                   " << pad << endl;
		cout << "Smoothing:                      " << smooth << endl;
		cout << "Cutoff(s):                      ";
		for ( i=0; i<nmap-1; ++i ) cout << "," << cutoff[i];
		cout << cutoff[i] << endl;
		cout << endl;
	}
	
	change_type(Float) ;
	statistics();
	p->change_type(Float) ;
	p->statistics();
	p->rescale_to_avg_std(avg, std);
	
	Bimage*			pr = copy_header(nmap);
	pr->data_type(Float);
	pr->size((maxrad - minrad)/step + 1, 1, 1);
	pr->sampling(image->sampling()[0]*step, 1, 1);
	pr->image->origin(-minrad/step, 0, 0);
	
	pr->data_alloc();

	pr->fill(fill);

	fft_plan		plan = fft_setup_plan(ft_size, ft_size, ft_size, FFTW_FORWARD, 1);
	
	Bimage*			pc;
	
	if ( verbose )
		cout << "Radius\tResolution" << endl;
	
	for ( r=minrad, i=0; r<=maxrad; r+=step, i++ ) {
		r1 = r - thickness/2.0;
		r2 = r + thickness/2.0;
		if ( r1 < 0 ) r1 = 0;
		cp1 = extract_shell(0, r1, r2);
		cp2 = p->extract_shell(0, r1, r2);
		if ( smooth ) {
			cp1->filter_average(3);
			cp2->filter_average(3);
		}
		if ( pad > 0 ) {
			cp1->pad(pad, fill_type, avg);
			cp2->pad(pad, fill_type, p->average());
		}
		pc = cp1->resolution_prepare(cp2, plan);
		delete cp1;
		delete cp2;
		plot = pc->fsc_dpr(hi_res, sampling_ratio);
		delete pc;
		for ( j=0; j<4; ++j ) {
			if ( cutoff[j] ) {
				res_est = 1/plot->cut(1, cutoff[j], -1);
				if ( res_est < hi_res ) res_est = hi_res;
				if ( res_est < LORESLIM ) pr->set(i + j*imgsize, res_est);
			}
		}
		if ( verbose )
			cout << r << tab << res_est << endl;
		delete plot;
	}
	
	fft_destroy_plan(plan);
	
	return pr;
}

double*		Bimage::fsc_local_voxel(Bimage* p, double hi_res, int size,
					int pad, int taper, double* cutoff, fft_plan plan, long i)
{
	check_resolution(hi_res);

	long			nn(0), xx, yy, zz, cc;
	int 			fill_type(FILL_AVERAGE);   	// Fill type for resizing
	double			sampling_ratio(1);			// Sampling ratio in reciprocal space

	double			width = size/4.0;
	Vector3<long> 	boxsize(size,size,size);	// Resolution box
	Vector3<long> 	boxhalf(boxsize/2);			// Half resolution box
	
	Vector3<long> 	taper_size(boxhalf);
	Vector3<double> taper_start(boxsize/4);
	
	coordinates(i, cc, xx, yy, zz, nn);

	Vector3<double>	start = Vector3<double>(xx - boxhalf[0], yy - boxhalf[1], zz - boxhalf[2]);
	
	Bimage*			cp1 = extract(nn, start, boxsize);
	Bimage*			cp2 = p->extract(nn, start, boxsize);
	
	if ( taper == 2 ) {
		cp1->hanning_taper(cp1->average());
		cp2->hanning_taper(cp2->average());
	} else if ( taper == 1 ) {
		cp1->edge(1, taper_size, taper_start, width, fill_type, cp1->average());
		cp2->edge(1, taper_size, taper_start, width, fill_type, cp2->average());
	}
	
	if ( pad > 0 ) {
		cp1->pad(pad, fill_type, cp1->average());
		cp2->pad(pad, fill_type, cp2->average());
	}

	Bimage*			pc = cp1->resolution_prepare(cp2, plan);
	delete cp1;
	delete cp2;

	Bplot*			plot = pc->fsc_dpr(hi_res, sampling_ratio);
	delete pc;
	
	double*			res_est = new double[4];
	for ( long j=0; j<4; ++j ) {
		res_est[j] = 0;
		if ( cutoff[j] ) {
			res_est[j] = 1/plot->cut(1, cutoff[j], -1);
			if ( res_est[j] < hi_res ) res_est[j] = hi_res;
			if ( res_est[j] > LORESLIM ) res_est[j] = LORESLIM;
		}
	}
	
	delete plot;
	
	return res_est;
}

/**
@author Giovanni Cardone and Bernard Heymann
@brief 	Determine the local resolution at each masked voxel in a map.
@param 	*p			second image.
@param 	*pmask		mask.
@param 	hi_res		high resolution limit.
@param 	*cutoff		correlation threshold(s).
@param 	mask_level	mask level index.
@param 	size		kernel size.
@param 	pad			padding factor.
@param 	vedge		edge size.
@param 	step		voxel step size.
@param 	taper		kernel tapering function.
@param 	fill		background fill value.
@return Bimage* 	local resolution image.

	For each voxel specified in the mask, within the edge limits, and at
	the step size given, two kernels are extracted from the two input maps.
	These kernels are then compared to determine the resolution based FSC.
	The resultant resolution value is then written into a new map.

**/
Bimage*		Bimage::fsc_local(Bimage* p, Bimage* pmask, double hi_res, double* cutoff,
				int mask_level, int size, int pad, Vector3<long> vedge, 
				int step, int taper, double fill)
{
	check_resolution(hi_res);

	int				dovox;
	long   			xx, yy, zz, i, nmap(0), nvox(0);
	long			ft_size(size), slice_size(x*y), imgsize(image_size());
	
	if ( vedge[0] < 0 ) vedge = Vector3<long>(size/2, size/2 ,size/2);

    if ( pad > 0 ) {
		pad = size*(pad+1);
		ft_size = pad;
	}
	
	// Count the number of output maps
	for ( i=0; i<4; ++i ) if ( cutoff[i] > 0 ) nmap++;
    
	if ( verbose ) {
		cout << "Calculating local resolution:" << endl;
		if ( pmask )
			cout << "Mask:                           " << pmask->file_name() << " (" << mask_level << ")" << endl;
		cout << "Kernel size:                    " << size << endl;
		cout << "Step size:                      " << step << endl;
		cout << "Edge size:                      " << vedge << endl;
		cout << "Padding size:                   " << pad << endl;
		cout << "Smoothing/tapering:             " << taper << endl;
		cout << "Cutoff(s):                      " << setprecision(3);
		for ( i=0; i<nmap-1; ++i ) cout << cutoff[i] << ",";
		cout << cutoff[i] << endl;
		cout << endl ;
	}
	
	change_type(Float) ;
	statistics();
	p->change_type(Float) ;
	p->statistics();
	p->rescale_to_avg_std(avg, std);
	
	Bimage*			pr = copy_header(nmap);
	pr->data_alloc();
	pr->fill(fill);
	pr->background(fill);
	
	int				nomask(0);
	if ( !pmask ) {
		nomask = 1;
		pmask = copy_header();
		pmask->data_type(UCharacter);
		pmask->data_alloc();
		pmask->fill(1);
	}

	// Set up the mask for which voxels to calculate
	for ( i=zz=0; zz<pmask->sizeZ(); zz++ ) {
		for ( yy=0; yy<pmask->sizeY(); yy++ ) {
			for ( xx=0; xx<pmask->sizeX(); xx++, i++ ) {
				dovox = 1;
				if ( zz%step || yy%step || xx%step ) dovox = 0;
				if ( (*pmask)[i] < 0.5 ) dovox = 0;
				if ( zz < vedge[2] || zz > pmask->sizeZ()-vedge[2] ) dovox = 0;
				if ( yy < vedge[1] || yy > pmask->sizeY()-vedge[1] ) dovox = 0;
				if ( xx < vedge[0] || xx > pmask->sizeX()-vedge[0] ) dovox = 0;
				if ( mask_level > 0 && fabs((*pmask)[i] - mask_level) > 0.1 ) dovox = 0;
				pmask->set(i, dovox);
			}
		}
	}

//	write_img("locres_mask.map", pmask);

	fft_plan		plan = fft_setup_plan(ft_size, ft_size, ft_size, FFTW_FORWARD, 1);
	
	if ( verbose ) {
		for ( i=nvox=0; i<datasize; i++ ) if ( (*pmask)[i] ) nvox++;
		cout << "Boxes to calculate:             " << nvox << endl << endl;
	}
/*
#ifdef HAVE_GCD
	__block	long		ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(imgsize, dispatch_get_global_queue(0, 0), ^(size_t k){
		if ( (*pmask)[k] ) {
			double*		res_est = fsc_local_voxel(p, hi_res,
				size, pad, taper, cutoff, plan, k);
			for ( long j=0; j<nmap; ++j ) pr->set(k + j*imgsize, res_est[j]);
			delete[] res_est;
			if ( verbose )
				dispatch_sync(myq, ^{
					ndone++;
					cerr << "Complete:                       " << setprecision(3)
							<< ndone*100.0/nvox << " %    \r" << flush;
				});
		}
	});
#else
	long				ndone(0);
#pragma omp parallel for
	for ( long k=0; k<imgsize; k++ ) if ( (*pmask)[k] ) {
		double*		res_est = fsc_local_voxel(p, hi_res,
				size, pad, taper, cutoff, plan, k);
		for ( long j=0; j<nmap; ++j ) pr->set(k + j*imgsize, res_est[j]);
		delete[] res_est;
	#pragma omp critical
		{
			ndone++;
			if ( verbose )
					cerr << "Complete:                       " << setprecision(3)
							<< ndone*100.0/nvox << " %    \r" << flush;
		}
	}
#endif
*/

#ifdef HAVE_GCD
	__block	long		ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(z, dispatch_get_global_queue(0, 0), ^(size_t zz){
		long			i, k;
		for ( i=0, k=zz*slice_size; i<slice_size; ++i, ++k ) {
			if ( (*pmask)[k] ) {
				double*		res_est = fsc_local_voxel(p, hi_res,
								size, pad, taper, cutoff, plan, k);
				for ( long j=0; j<nmap; ++j ) pr->set(k + j*imgsize, res_est[j]);
				delete[] res_est;
				ndone++;
			}
		}
		if ( verbose ) {
			dispatch_sync(myq, ^{
				cerr << "Complete:                       " << setprecision(3)
						<< ndone*100.0/nvox << " %    \r" << flush;
			});
		}
	});
#else
	long				ndone(0);
#pragma omp parallel for
	for ( long zz=0; zz<z; zz++ ) {
		long			i, k;
		for ( i=0, k=zz*slice_size; i<slice_size; ++i, ++k ) {
			if ( (*pmask)[k] ) {
				double*		res_est = fsc_local_voxel(p, hi_res,
							size, pad, taper, cutoff, plan, k);
				for ( long j=0; j<nmap; ++j ) pr->set(k + j*imgsize, res_est[j]);
				delete[] res_est;
				ndone++;
			}
		}
	#pragma omp critical
		{
			if ( verbose )
					cerr << "Complete:                       " << setprecision(3)
							<< ndone*100.0/nvox << " %    \r" << flush;
		}
	}
#endif

	fft_destroy_plan(plan);
	
	for ( i=nvox=0; i<imgsize; i++ ) if ( fabs((*pr)[i] - fill) > 1e-10 ) nvox++;
	
	if ( verbose ) {
		cerr << endl;
		cout << "Boxes calculated:               " << nvox << endl << endl;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "x\ty\tz\tResolution" << endl;
		for ( zz=0; zz<pr->sizeZ(); zz+=step ) {
			for ( yy=0; yy<p->sizeY(); yy+=step ) {
				for ( xx=0; xx<p->sizeX(); xx+=step ) {
					pr->index(xx,yy,zz);
					if ( (*pmask)[i] )
						cout << xx << tab << yy << tab << zz << tab << (*pr)[i] << endl;
				}
			}
		}
	}

	if ( nomask ) {
		delete pmask;
		pmask = NULL;
	}
	
	pr->statistics();
	
	return pr;
}

double		Bimage::local_filter_voxel(Bimage* resmap, long size,
					fft_plan planf, fft_plan planb, long i)
{
	long			nn, xx, yy, zz, cc;
	Bimage*			cp = NULL;

	double			res_hi(0);					// High resolution limit
	double			res_lo(1e37); 				// Low resolution limit
	double			bandpass_width(0.002);		// Edge width for the bandpass filter

	Vector3<long> 	boxsize(size,size,size);	// Resolution box
	Vector3<long> 	boxhalf(boxsize/2);			// Half resolution box
    Vector3<double>	start;						// Origin of kernel to extract
	long			cenvox = (boxsize[2]/2*boxsize[1] + boxsize[1]/2)*boxsize[0] + boxsize[0]/2;
	
	coordinates(i, cc, xx, yy, zz, nn);

	start = Vector3<double>(xx - boxhalf[0], yy - boxhalf[1], zz - boxhalf[2]);
	res_hi = (*resmap)[i];
	cp = extract(nn, start, boxsize);
	cp->fspace_bandpass(res_hi, res_lo, bandpass_width, planf, planb);
	
	// read value of cp1 at the center of the box and assign to output
	double val =  (*cp)[cenvox];
	
	delete cp;
	
	return val;
}

/**
@author Giovanni Cardone and Bernard Heymann
@brief 	Applies a local resolution filter to a map.
@param 	*pmask		mask of areas to generate.
@param 	mask_level	mask level.
@param 	*resmap		local resolution map.
@param 	size		kernel size.
@param 	vedge		edge size.
@return Bimage* 	filtered image, NULL on error.

	At each voxel, a kernel/small box is extracted and lowpass filtered
	to the resolution limit for that voxel in the local resolution map.
	A mask may be used to limit the region of application.

**/
Bimage*		Bimage::local_filter(Bimage* pmask, int mask_level, Bimage* resmap, 
				int size, Vector3<long> vedge)
{
	if ( !resmap ) {
		cerr << "Error: A local resolution map must be provided!" << endl;
		return NULL;
	}

	change_type(Float);
	
	if ( vedge[0] < 0 ) vedge = Vector3<long>(size/2, size/2 ,size/2);

	long			xx, yy, zz, nn, i, nvox(0);
	int				dovox;

	if ( verbose ) {
		cout << "Filtering by local resolution:" << endl;
		cout << "Local resolution image:         " << resmap->file_name() << endl;
		if ( pmask )
			cout << "Mask:                           " << pmask->file_name() << " (" << mask_level << ")" << endl;
		cout << "Kernel size:                    " << size << endl;
		cout << "Edge size:                      " << vedge << endl;
		cout << endl;
	}
	
    Bimage*			p = copy_header();
	p->data_alloc();

	int				nomask(0);
	if ( !pmask ) {
		nomask = 1;
		pmask = copy_header();
		pmask->data_type(UCharacter);
		pmask->data_alloc();
		pmask->fill(1);
	}

	// Set up the mask for which voxels to calculate
	for ( i=zz=0; zz<pmask->sizeZ(); zz++ ) {
		for ( yy=0; yy<pmask->sizeY(); yy++ ) {
			for ( xx=0; xx<pmask->sizeX(); xx++, i++ ) {
				dovox = 1;
				if ( (*pmask)[i] < 0.99 ) dovox = 0;
				if ( zz < vedge[2] || zz > pmask->sizeZ()-vedge[2] ) dovox = 0;
				if ( yy < vedge[1] || yy > pmask->sizeY()-vedge[1] ) dovox = 0;
				if ( xx < vedge[0] || xx > pmask->sizeX()-vedge[0] ) dovox = 0;
				if ( mask_level > 0 && fabs((*pmask)[i] - mask_level) > 0.1 ) dovox = 0;
				pmask->set(i, dovox);
			}
		}
	}

	for ( i=0; i<datasize; i++ ) if ( (*pmask)[i] ) nvox++;
	
	fft_plan		planf = fft_setup_plan(size, size, size, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(size, size, size, FFTW_BACKWARD, 1);

//	write_img("locres_mask.map", pmask);

#ifdef HAVE_GCD
	__block	long		ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t k){
		if ( (*pmask)[k] ) {
			p->set(k, local_filter_voxel(resmap, size,  planf, planb, k));
			if ( verbose )
				dispatch_sync(myq, ^{
					ndone++;
					cerr << "Complete:                       " << setprecision(3)
							<< ndone*100.0/nvox << " %    \r" << flush;
				});
		}
	});
#else
	long				ndone(0);
#pragma omp parallel for
	for ( long j=0; j<datasize; j++ ) if ( (*pmask)[j] ) {
		p->set(j, local_filter_voxel(resmap, size,  planf, planb, j));
	#pragma omp critical
		{
			ndone++;
			if ( verbose )
					cerr << "Complete:                       " << setprecision(3)
							<< ndone*100.0/nvox << " %    \r" << flush;
		}
	}
#endif

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);

	if ( verbose & VERB_FULL ) {
		cout << "ImgNum\tx\ty\tz\tResolution" << endl;
		for ( i=nn=0; nn<p->images(); nn++ ) {
	    	for ( zz=0; zz<p->sizeZ(); zz++ ) {
	    		for ( yy=0; yy<p->sizeY(); yy++ ) {
	    	    	for ( xx=0; xx<p->sizeX(); xx++, i++ ) {
						if ( (*pmask)[i] )
							cout << nn << tab << xx << tab << yy << tab << zz << tab << (*p)[i] << endl;
					}
				}
			}
		}
	}
	
	if ( nomask ) {
		delete pmask;
		pmask = NULL;
	}
	
	p->statistics();
	
	return p;
}
