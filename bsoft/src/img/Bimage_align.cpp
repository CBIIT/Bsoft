/**
@file	Bimage_align.cpp
@brief	Functions to align images.
@author Bernard Heymann
@date	Created: 20000505
@date	Modified: 20210809
**/

#include "Bimage.h"
#include "simplex.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Aligns and sums a set of sub-images using a progressive algorithm.
@param 	nref		reference sub-image.
@param 	*pmask		cross-correlation mask.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit.
@param 	shift_limit	limit on extent of search for shift.
@param 	planf		plan for forward Fourier transform.
@param 	planb		plan for backward Fourier transform.
@return Bimage*		summed reference image.

	The images are aligned against the reference images with progressive
	addition to reference to decrease noise.

**/
Bimage*		Bimage::align_progressive(long nref, Bimage* pmask, 
				double hi_res, double lo_res, double shift_limit,
				fft_plan planf, fft_plan planb)
{
	if ( nref < 0 || nref >= n ) nref = 0;

	long			i, nn;
	double			cc(0), shift_avg(0), cc_avg(0);
	Vector3<double>	shift;
	Bimage*			pref = extract(nref);
	Bimage*			p1;
	
	if ( verbose ) {
		cout << "Progressive alignment:" << endl;
		cout << "Image\tShift\t\t\tCC" << endl;
	}
	
	for ( nn=nref+1, i=1; i<n; nn++, i++ ) {
		if ( nn >= n ) nn = 0;
		p1 = extract(nn);
//		cout << p1->average() << endl;
		shift = p1->find_shift(pref, pmask, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);
		image[nn].origin(pref->image->origin()+shift);
		image[nn].FOM(cc);
		shift = -shift;
		p1->shift_wrap(shift);
		pref->add(p1);
		delete p1;
//		shift_avg += shift.length();
		cc_avg += cc;
		if ( verbose )
			cout << nn+1 << tab << shift << tab << cc << endl;
	}

	for ( nn=1; nn<n; nn++ ) {
		shift = image[nn].origin() - image[nn-1].origin();
		shift_avg += shift.length();
	}
	
	shift_avg /= n-1;
	cc_avg /= n-1;
	pref->show_scale(shift_avg);
	pref->image->FOM(cc_avg);

	if ( verbose ) {
		cout << "Shift average:                 " << shift_avg << endl;
		cout << "CC average:                    " << cc_avg << endl;
	}
	
	return pref;
}

/**
@brief 	Aligns and sums a set of sub-images using a progressive algorithm.
@param 	nref		reference sub-image.
@param 	*pmask		cross-correlation mask.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit.
@param 	shift_limit	limit on extent of search for shift.
@param 	planf		plan for forward Fourier transform.
@param 	planb		plan for backward Fourier transform.
@return Bimage*		summed reference image.

	The images are aligned against the reference images with progressive
	addition to reference to decrease noise.

**/
Bimage*		Bimage::align_local(long nref, Bimage* pmask,
				double hi_res, double lo_res, double shift_limit,
				fft_plan planf, fft_plan planb)
{
	if ( nref < 0 || nref >= n ) nref = 0;

	long			nn;
	double			shift_avg(0), cc_avg(0);
	Vector3<double>	shift;
	Bimage*			p1;
	
	if ( verbose ) {
		cout << "Local alignment:" << endl;
		cout << "Image\tShift\t\t\tCC" << endl;
	}
	
#ifdef HAVE_GCD
	__block	vector<Vector3<double>>  sh(n);
	sh[0] = Vector3<double>(0,0,0);
	dispatch_apply(n-1, dispatch_get_global_queue(0, 0), ^(size_t nn){
		double			cc(0);
		Bimage*			p1 = extract(nn);
		Bimage*			p2 = extract(nn+1);
		sh[nn+1] = p2->find_shift(p1, pmask, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);
		image[nn+1].FOM(cc);
		delete p1;
		delete p2;
	});
#else
	vector<Vector3<double>>  	sh(n);
	sh[0] = Vector3<double>(0,0,0);
#pragma omp parallel for
	for ( nn=1; nn<n; ++nn ) {
		double			cc(0);
		Bimage*			p1 = extract(nn-1);
		Bimage*			p2 = extract(nn);
		sh[nn] = p2->find_shift(p1, pmask, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);
		image[nn].FOM(cc);
		delete p1;
		delete p2;
	}
#endif

	for ( nn=1; nn<n; ++nn )
		sh[nn] += sh[nn-1];

	for ( nn=0; nn<n; ++nn )
		sh[nn] -= sh[nref];

	Bimage*			pref = extract(nref);
	pref->origin(pref->size()/2);

	for ( nn=0; nn<n; ++nn ) {
		image[nn].origin(sh[nn] + pref->image->origin());
		shift = -sh[nn];
		if ( nn != nref ) {
			p1 = extract(nn);
			p1->shift_wrap(shift);
			pref->add(p1);
			delete p1;
		}
		shift_avg += shift.length();
		cc_avg += image[nn].FOM();
		if ( verbose )
			cout << nn+1 << tab << shift << tab << image[nn].FOM() << endl;
	}

	shift_avg /= n-1;
	cc_avg /= n-1;
	pref->show_scale(shift_avg);
	pref->image->FOM(cc_avg);

	if ( verbose ) {
		cout << "Shift average:                 " << shift_avg << endl;
		cout << "CC average:                    " << cc_avg << endl;
	}
	
	return pref;
}

Vector3<double>	Bimage::find_shift_in_transform(long nn, Bimage* pref, double shift_limit)
{
	Bimage*			p1 = extract(nn);
	cout << "origin = " << p1->image->origin() << endl;
	
//	long			h, k;
//	double			cc(0), ccmax(0);
//	Vector3<double>	shift, bestshift;
	
	p1->complex_conjugate_product(pref, 1);
/*
		for ( k=-1; k<2; ++k ) {
			shift[1] = k;
			for ( h=-1; h<2; ++h ) {
				shift[0] = h;
				cc = p1->correlation_coefficient(shift);
				if ( ccmax < cc ) {
					ccmax = cc;
					bestshift = shift;
				}
			}
		}
*/
	p1->fft(FFTW_BACKWARD, 0);

	p1->complex_to_real();	

	p1->find_peak(shift_limit, 0);
	cout << "origin = " << p1->image->origin() << endl;

	p1->refine_peak();
	
	Vector3<double>		shift = p1->image->origin();
	
	image[nn].FOM(p1->image->FOM());
	
	image[nn].origin(pref->image->origin()-shift);
//	image[nn].FOM(cc);

	cout << "shift = " << shift << endl;
	cout << "fom = " << p1->image->FOM() << endl;
	
	delete p1;
	
	return shift;
}

/**
@brief 	Aligns and sums a set of sub-images using a progressive algorithm.
@param 	nref		reference sub-image.
@param 	shift_limit	limit on extent of search for shift.
@return Bimage*		summed reference image.

	The images are aligned against the reference images with progressive
	addition to reference to decrease noise.

**/
Bimage*		Bimage::align_progressive_fast(long nref, double shift_limit)
{
	if ( nref < 0 || nref >= n ) nref = 0;

	long			i, nn;
	double			shift_avg(0), cc_avg(0);
	Vector3<double>	shift, shift_prev;
	Bimage*			pref = extract(nref);
	Bimage*			p1;
	
	if ( verbose ) {
		cout << "Progressive alignment:" << endl;
		cout << "Image\tShift\t\t\tCC" << endl;
	}
	
	for ( nn=nref+1, i=1; i<n; nn++, i++ ) {
		if ( nn >= n ) nn = 0;
		shift = find_shift_in_transform(nn, pref, shift_limit);
		p1 = extract(nn);
		p1->phase_shift(shift);
		pref->add(p1);
		delete p1;
		pref->fspace_normalize();
		cc_avg += image[nn].FOM();
		shift_prev = shift;
		if ( verbose )
			cout << nn+1 << tab << shift << tab << image[nn].FOM() << endl;
	}

	for ( nn=1; nn<n; nn++ ) {
		shift = image[nn].origin() - image[nn-1].origin();
		shift_avg += shift.length();
	}
	
	shift_avg /= n-1;
	cc_avg /= n-1;
	pref->show_scale(shift_avg);
	pref->image->FOM(cc_avg);

	if ( verbose ) {
		cout << "Shift average:                 " << shift_avg << endl;
		cout << "CC average:                    " << cc_avg << endl;
	}
	
	return pref;
}

vector<Vector3<double>>	interpolate_shifts(Vector3<double>* sh, long nimg, long window, long step)
{
	long					nsh = (nimg+step-1)/step;

	if ( nsh != nimg ) {
		cerr << "Error in interpolate_shifts: The number of shifts must equal the number of images!" << endl;
		bexit(-1);
	}
	
	vector<Vector3<double>>	shfull(nsh);
	
	if ( verbose )
		cout << "Interpolating shifts over window of " << window << " and step of " << step << endl;
	
	long			i, j, h((step-1)/2), w, ws, we;
	double			f;
	double			c = (window%2 || step==1)? 0: 0.5;		// Center of window
	
	for ( i = 0; i < nsh; ++i, c += step ) {
		ws = i*step - h;
		we = ws + step;
		if ( ws < 0 ) ws = 0;
		if ( we > nimg ) we = nimg;
		if ( i == nsh-1 ) we = nimg;
		for ( w = ws; w < we; ++w ) {
			f = (w-c)/step;
//			cout << c << tab << w << tab << ws << tab << f << endl;
			if ( i ) {
				j = i - 1;
				shfull[i] = sh[i] - (sh[i] - sh[j])*f;
			} else {
				shfull[0] = sh[0] - (sh[1] - sh[0])*f;
			}
//			if ( verbose & VERB_FULL )
				cout << i << tab << sh[i][0] << tab << sh[i][1] << tab
					<< f << tab << w << tab << shfull[i][0] << tab << shfull[i][1] << endl;
		}
	}
	
	return shfull;
}

/**
@brief 	Aligns and sums a set of sub-images, first progressively and then iteratively.
@param 	ref_num			reference sub-image.
@param 	window			moving sum window.
@param 	step			moving sum interval.
@param 	*pmask			cross-correlation mask.
@param 	hi_res			high resolution limit.
@param 	lo_res			low resolution limit.
@param 	shift_limit		limit on extent of search for shift.
@param 	edge_width		width of smoothing edge.
@param 	gauss_width		decay coefficient for smoothing edge.
@param 	aln_bin			3-value vector indicating binning level.
@param 	mode			flag to select initial alignment: 0=progressive, 1=local.
@return vector<Vector3<double>>	vector of shifts.

	The images are first aligned using a progressive algorith starting with 
	the indicated reference image.
	The images are then aligned iteratively using the average image to 
	improve the average.
	The image is extensively modified and should not be used for further processing.

**/
vector<Vector3<double>>	Bimage::align(long ref_num, long window, long step, Bimage* pmask,
				double hi_res, double lo_res, double shift_limit, 
				double edge_width, double gauss_width, Vector3<long> aln_bin, int mode)
{	
	if ( ref_num < 0 || ref_num >= n ) ref_num = 0;
	
	change_type(Float);

	origin(size()/2);

	calculate_background();
	
	bin(aln_bin);
	
	Bimage*			pt = this;
	if ( window > 1 ) pt = moving_sum(window, step, 0);
	
//	pt->information();

	Vector3<double>	edge_origin(edge_width, edge_width, 0);
	Vector3<long>	edge_size = pt->size() - (long) (2*edge_width);
	edge_size = edge_size.max(1);
	
//	cout << edge_size << tab << pt->image->background() << endl;
	
	if ( edge_width > 0 )
		pt->edge(0, edge_size, edge_origin, gauss_width, FILL_BACKGROUND);
	
	
//	write_img("t.pif", this);
//	exit(-1);
	
	fft_plan		planf = pt->fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = pt->fft_setup(FFTW_BACKWARD, 0);

	Bimage*			pref = NULL;
	if ( mode == 0 )
		pref = pt->align_progressive(ref_num, pmask, hi_res, lo_res, shift_limit, planf, planb);
	else
		pref = pt->align_local(ref_num, pmask, hi_res, lo_res, shift_limit, planf, planb);
	
	double			cc_avg(0), cc_avg_old(0), shift_avg, shift_avg_old(0);
	shift_avg = pref->show_scale() * aln_bin[0];
	cc_avg = pref->image->FOM();
	pref->origin(pref->size()/2);
	
//	write_img("progaln.pif", pref);
	
	long				i, nn;
	Bimage*				p1;
	Vector3<double>*	sh = new Vector3<double>[pt->images()];

	if ( verbose )
		cout << "Iterative refinement:" << endl;
	for ( i=1; i<=10 && fabs(shift_avg - shift_avg_old) > 0.01
			&& cc_avg > cc_avg_old; i++ ) {
		pt->origin(pref->size()/2);
		
		if ( verbose )
			cout << "Iteration " << i << endl;
#ifdef HAVE_GCD
		dispatch_apply(pt->images(), dispatch_get_global_queue(0, 0), ^(size_t nn){
			sh[nn] = pt->find_shift(nn, pref, pmask, hi_res, lo_res, shift_limit, planf, planb);
		});
#else
#pragma omp parallel for
		for ( nn=0; nn<pt->images(); nn++ ) {
			sh[nn] = pt->find_shift(nn, pref, pmask, hi_res, lo_res, shift_limit, planf, planb);
		}
#endif

		pref->clear();
		pref->origin(pref->size()/2);
		
		shift_avg_old = shift_avg;
		cc_avg_old = cc_avg;
	
		for ( nn=0; nn<n; nn++ ) sh[nn] = sh[ref_num] - sh[nn];
		
		if ( verbose )
			cout << "Image\tShift\t\t\tCC" << endl;
		for ( cc_avg=0, shift_avg=0, nn=0; nn<pt->images(); nn++ ) {
			p1 = pt->extract(nn);
			p1->shift_wrap(sh[nn]);
			pref->add(p1);
			delete p1;
			cc_avg += pt->image[nn].FOM();
			sh[nn] *= aln_bin;
			if ( nn ) shift_avg += sh[nn].distance(sh[nn-1]);
			if ( verbose )
				cout << nn+1 << tab << sh[nn] << tab << pt->image[nn].FOM() << endl;
		}
		
		shift_avg /= pt->images()-1;
		cc_avg /= pt->images();
		
		if ( verbose ) {
			cout << "Shift average per frame:       " << shift_avg << endl;
			cout << "CC average:                    " << cc_avg << endl;
		}
	}
	
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	delete pref;

	for ( nn=0; nn<pt->images(); nn++ ) sh[nn][2] = pt->image[nn].FOM();

	vector<Vector3<double>>		shfull = interpolate_shifts(sh, n, window, step);
	
/*
	vector<double>	vsh(3);
	JSvalue			job(JSobject);
	JSvalue			jsh(JSarray);
	
	for ( nn=0; nn<n; nn++ ) {
		image[nn].origin(sh[nn] + size()/2);
		image[nn].FOM(pt->image[nn].FOM());
		vsh = {sh[nn][0], sh[nn][1], pt->image[nn].FOM()};
		jsh.push_back(vsh);
	}

	job["shifts"] = jsh;
	job["shift_average"] = shift_avg;
	job["FOM"] = cc_avg;
*/
	if ( window > 1 ) delete pt;
//	delete[] sh;
	
//	return job;
	return shfull;
}
/*
JSvalue		Bimage::align(long ref_num, long window, long step, Bimage* pmask,
				double hi_res, double lo_res, double shift_limit,
				double edge_width, double gauss_width, Vector3<long> bin)
{
	if ( ref_num < 0 || ref_num >= n ) ref_num = 0;
	
	change_type(Float);

	origin(size()/2);

	calculate_background();
	
	Vector3<long>	obin(bin[0], bin[1], bin[2]);

	Bimage*			pt = bin_copy(bin);
	if ( window > 1 ) {
		Bimage*		pma = pt->moving_sum(window, step, 0);
		delete pt;
		pt = pma;
	}

	Vector3<double>	edge_origin(edge_width, edge_width, 0);
	Vector3<long>	edge_size = pt->size() - (long) (2*edge_width);
	edge_size = edge_size.max(1);
	
//	cout << edge_size << tab << pt->image->background() << endl;
	
	if ( edge_width > 0 )
		pt->edge(0, edge_size, edge_origin, gauss_width, FILL_BACKGROUND);
	
	
//	write_img("t.pif", this);
//	exit(-1);
	
	fft_plan		planf = pt->fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = pt->fft_setup(FFTW_BACKWARD, 0);

	Bimage*			pref = pt->align_progressive(ref_num,
						pmask, hi_res, lo_res, shift_limit, planf, planb);
	double			cc_avg(0), cc_avg_old(0), shift_avg, shift_avg_old(0);
	shift_avg = pref->show_scale() * bin[0];
	cc_avg = pref->image->FOM();
	pref->origin(pref->size()/2);
	
//	write_img("progaln.pif", pref);
	
	long				i, nn;
	Bimage*				p1;
	Vector3<double>*	sh = new Vector3<double>[pt->images()];

	if ( verbose )
		cout << "Iterative refinement:" << endl;
	for ( i=1; i<=10 && fabs(shift_avg - shift_avg_old) > 0.01
			&& cc_avg > cc_avg_old; i++ ) {
		pt->origin(pref->size()/2);
		
		if ( verbose )
			cout << "Iteration " << i << endl;
#ifdef HAVE_GCD
		dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
			sh[nn] = pt->find_shift(nn, pref, pmask, hi_res, lo_res, shift_limit, planf, planb);
		});
#else
#pragma omp parallel for
		for ( nn=0; nn<n; nn++ ) {
			sh[nn] = pt->find_shift(nn, pref, pmask, hi_res, lo_res, shift_limit, planf, planb);
		}
#endif

		pref->clear();
		pref->origin(pref->size()/2);
		
		shift_avg_old = shift_avg;
		cc_avg_old = cc_avg;
	
		for ( nn=0; nn<n; nn++ ) sh[nn] = sh[ref_num] - sh[nn];
		
		if ( verbose )
			cout << "Image\tShift\t\t\tCC" << endl;
		for ( cc_avg=0, shift_avg=0, nn=0; nn<n; nn++ ) {
			p1 = pt->extract(nn);
			p1->shift_wrap(sh[nn]);
			pref->add(p1);
			delete p1;
			cc_avg += pt->image[nn].FOM();
			sh[nn] *= obin;
			if ( nn ) shift_avg += sh[nn].distance(sh[nn-1]);
			if ( verbose )
				cout << nn+1 << tab << sh[nn] << tab << pt->image[nn].FOM() << endl;
		}
		
		shift_avg /= n-1;
		cc_avg /= n;
		
		if ( verbose ) {
			cout << "Shift average per frame:       " << shift_avg << endl;
			cout << "CC average:                    " << cc_avg << endl;
		}
	}
	
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	delete pref;

//	vector<Vector3<double>>		shfull = interpolate_shifts(sh, n, window, step);

	vector<double>	vsh(3);
	JSvalue			job(JSobject);
	JSvalue			jsh(JSarray);
	
	for ( nn=0; nn<n; nn++ ) {
		image[nn].origin(sh[nn] + size()/2);
		image[nn].FOM(pt->image[nn].FOM());
		vsh = {sh[nn][0], sh[nn][1], pt->image[nn].FOM()};
		jsh.push_back(vsh);
	}

	job["shifts"] = jsh;
	job["shift_average"] = shift_avg;
	job["FOM"] = cc_avg;

	delete pt;
	delete[] sh;
	
	return job;
}
*/

/**
@brief 	Aligns and sums a set of sub-images, first progressively and then iteratively.
@param 	ref_num		reference sub-image.
@param 	*pmask		cross-correlation mask.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit.
@param 	shift_limit	limit on extent of search for shift.
@param 	edge_width	width of smoothing edge.
@param 	gauss_width	decay coefficient for smoothing edge.
@return double		average shift.

	The images are first aligned using a progressive algorith starting with 
	the indicated reference image.
	The images are then aligned iteratively using the average image to 
	improve the average.

**/
JSvalue		Bimage::align_fast(long ref_num, Bimage* pmask,
				double hi_res, double lo_res, double shift_limit, 
				double edge_width, double gauss_width)
{	
	if ( ref_num < 0 || ref_num >= n ) ref_num = 0;
	check_resolution(hi_res);
	
	change_type(Float);

	origin(size()/2);

	calculate_background();
	
	Bimage*			pt = copy();

	Vector3<double>	edge_origin(edge_width, edge_width, 0);
	Vector3<long>	edge_size = pt->size() - (long) (2*edge_width);
	edge_size = edge_size.max(1);
	
//	cout << edge_size << tab << pt->image->background() << endl;
	
	if ( edge_width > 0 )
		pt->edge(0, edge_size, edge_origin, gauss_width, FILL_BACKGROUND);
	
	double			scale(hi_res/(2*sampling(0)[0]));
	Vector3<long>	nusize(size()/scale);
	nusize = nusize.max(1);
	if ( nusize[0]%2 ) nusize[0]++;
	if ( nusize[1]%2 ) nusize[1]++;
	if ( z > 1 && nusize[2]%2 ) nusize[2]++;
	scale = x/nusize[0];
	
	pt->fft();
	
	pt->zero_origin();
	
//	pt->change_transform_size(nusize);
	scale = 1;

	cout << "new size = " << nusize << endl;
	cout << "new sampling = " << pt->sampling(0) << endl;
	cout << "new origin = " << pt->image->origin() << endl;
	
	pt->fspace_bandpass(hi_res, lo_res, 0);
	
	pt->fspace_normalize();

	Bimage*			pref = pt->align_progressive_fast(ref_num, shift_limit);
	
	double			cc_avg(0), cc_avg_old(0), shift_avg, shift_avg_old(0);
	shift_avg = pref->show_scale() * scale;
	cc_avg = pref->image->FOM();
	pref->origin(pref->size()/2);
	
//	write_img("progaln.pif", pref, 0);
	
	long				i, nn;
	Bimage*				p1;
	Vector3<double>*	sh = new Vector3<double>[pt->images()];

	if ( verbose )
		cout << "Iterative refinement:" << endl;
	for ( i=1; i<=10 && fabs(shift_avg - shift_avg_old) > 0.01
			&& cc_avg > cc_avg_old; i++ ) {
		pt->origin(pref->size()/2);
		
		if ( verbose )
			cout << "Iteration " << i << endl;
#ifdef HAVE_GCD
		dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
			sh[nn] = pt->find_shift_in_transform(nn, pref, shift_limit);
		});
#else
#pragma omp parallel for
		for ( nn=0; nn<n; nn++ ) {
			sh[nn] = pt->find_shift_in_transform(nn, pref, shift_limit);
		}
#endif

		pref->clear();
		pref->origin(pref->size()/2);
		
		shift_avg_old = shift_avg;
		cc_avg_old = cc_avg;
	
		for ( nn=0; nn<n; nn++ ) sh[nn] = sh[ref_num] - sh[nn];
		
		if ( verbose )
			cout << "Image\tShift\t\t\tCC" << endl;
		for ( cc_avg=0, shift_avg=0, nn=0; nn<n; nn++ ) {
			p1 = pt->extract(nn);
			p1->phase_shift(sh[nn]);
			pref->add(p1);
			delete p1;
			cc_avg += pt->image[nn].FOM();
			sh[nn] *= scale;
			if ( nn ) shift_avg += sh[nn].distance(sh[nn-1]);
			if ( verbose )
				cout << nn+1 << tab << sh[nn] << tab << pt->image[nn].FOM() << endl;
		}
		
		pref->fspace_normalize();
		
		shift_avg /= n-1;
		cc_avg /= n;
		
		if ( verbose ) {
			cout << "Shift average per frame:       " << shift_avg << endl;
			cout << "CC average:                    " << cc_avg << endl;
		}
	}
	
	delete pref;
	
	vector<double>	vsh(3);
	JSvalue			job(JSobject);
	JSvalue			jsh(JSarray);
	
	for ( nn=0; nn<n; nn++ ) {
		vsh = {sh[nn][0], sh[nn][1], pt->image[nn].FOM()};
		jsh.push_back(vsh);
	}

	job["shifts"] = jsh;
	job["shift_average"] = shift_avg;
	job["FOM"] = cc_avg;
/*
	if ( verbose )
		cout << "Shifting original frames" << endl << endl;
	
	fft();
	
#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t nn){
		fspace_translate(nn, sh[nn]);
		image[nn].FOM(pt->image[nn].FOM());
	});
#else
#pragma omp parallel for
	for ( nn=0; nn<n; nn++ ) {
		fspace_translate(nn, sh[nn]);
		image[nn].FOM(pt->image[nn].FOM());
	}
#endif
	
	fft_back();
*/
	delete pt;
	delete[] sh;
	
	return job;
}

/*
@brief	Frequency space shift and sum of aligned images.
@param	shift		flag to shift to origin.
@return Bimage* 	frequency space sum of images linked to sum of power.

	The sub-image shifts must be encoded in the origins as translation from the center.
**/
Bimage*		Bimage::fspace_sum(int shift)
{
	if ( fouriertype == NoTransform )
		fft();

	long			i, j, nn, nimg(0), is(image_size());
	Complex<float>	cv;

	Bimage*			psum = copy_header(1);
	psum->origin(Vector3<double>(0,0,0));
	psum->data_alloc_and_clear();

	Bimage*			ppow = copy_header(1);
	ppow->origin(Vector3<double>(0,0,0));
	ppow->compound_type(TSimple);
	ppow->data_type(Float);
	ppow->data_alloc_and_clear();
	psum->next = ppow;

	for ( nn = 0; nn < n; ++nn ) {
		if ( image[nn].select() ) {
			nimg++;
			if ( shift ) fspace_translate(nn, image[nn].origin() - size()/2);
			for ( i=nn*is, j=0; j<is; ++j, ++i ) {
				cv = complex(i);
				psum->add(j, cv);
				ppow->add(j, cv.power());
			}
		}
	}

	return psum;
}

/*
@brief	Frequency space shift and sum of aligned images.
@param	subset		number of images to sum in each subset.
@param	flag		bit 1=shift origins, bit 2=progressive sum.
@return Bimage* 	frequency space sum of images linked to sum of power.

	The sub-image shifts must be encoded in the origins as translation from the center.
**/
Bimage*		Bimage::fspace_subset_sums(int subset, int flag)
{
	if ( fouriertype == NoTransform )
		fft();

	int				shift(flag&1);
	long			i, j, k, m, nn, is(image_size());
	long			nimg(n/subset);
	Complex<double>	cv;

	Bimage*			psum = copy_header(nimg);
	psum->origin(Vector3<double>(0,0,0));
	psum->data_alloc_and_clear();

	Bimage*			ppow = copy_header(nimg);
	ppow->origin(Vector3<double>(0,0,0));
	ppow->compound_type(TSimple);
	ppow->data_type(Float);
	ppow->data_alloc_and_clear();
	psum->next = ppow;

	for ( nn = 0; nn < n; ++nn ) {
		if ( shift ) fspace_translate(nn, image[nn].origin() - size()/2);
		for ( i=nn*is, k=is*(nn/subset), m=k-is, j=0; j<is; ++j, ++i, ++k, ++m ) {
			cv = complex(i);
			psum->add(k, cv);
			ppow->add(k, cv.power());
		}
	}
	
	return psum;
}

double		mtf_gaussian_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df, s2;
	double			invsig2(-0.5/(simp.parameter(2)*simp.parameter(2)));
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=1; i<simp.points(); i++ ) {
		s2 = x[i]*x[i];
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*exp(invsig2*s2));
		R += df*df;
	}
	
	R = sqrt(R/i);
	R /= simp.dependent_variance();
			
	return R;
}

double		coincidence_loss_R(Bsimplex& simp)
{
	long			i;
	double			R(0), sf, df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=1; i<simp.points(); ++i ) {
		sf = sin(simp.parameter(2)*x[i])/(simp.parameter(2)*x[i]);
		df = f[i] - simp.parameter(0)*(1 - simp.parameter(1)*sf*sf);
		R += df*df;
	}
	
	R /= i;
	R /= simp.dependent_variance();
	
	return R;
}

double		fit_mtf(Bimage* prad, long nimg, double& off, double& amp, double& sig)
{
	Bimage*			pradpow = prad->next;
	long			i;
	vector<double>	s(prad->sizeX(),0), f(prad->sizeX(),0);
	
//	cout << "sampling=" << prad->sampling(0)[0] << endl;
	
	for ( i=0; i<prad->sizeX(); ++i ) {
		s[i] = i*prad->sampling(0)[0];
		f[i] = ((*pradpow)[i] - (*prad)[i]/nimg);
//		cout << s[i] << tab << f[i] << endl;
	}
	
//	cout << f.back() << tab << f[1] << endl;
	off = f.back();
	amp = f[5];
	sig = s.back()/2;
	
	if ( verbose )
		cout << "Fitting mass transfer function:" << endl;
	if ( verbose & VERB_FULL ) {
		cout << "Starting values:" << endl;
		cout << "Offset:                        " << off << endl;
		cout << "Amplitude:                     " << amp << endl;
		cout << "Sigma:                         " << sig << " 1/Å" << endl;
	}
	
	Bsimplex			simp(1, 3, 0, s.size(), s, f);

	simp.parameter(0, off/2);
	simp.parameter(1, amp);
	simp.parameter(2, sig);
	simp.limits(0, 0, 2*off);
	simp.limits(1, -off, 5*amp);
	simp.limits(2, s[1], s.back());
	
	double			R = simp.run(10000, 1e-5, mtf_gaussian_R);
	
	off = simp.parameter(0);
	amp = simp.parameter(1);
	sig = simp.parameter(2);
	
	if ( verbose ) {
		cout << "Offset:                        " << off << endl;
		cout << "Amplitude:                     " << amp << endl;
		cout << "Sigma:                         " << sig << " 1/Å" << endl;
		cout << "R:                             " << R << endl;
	}

	return R;
}

double		fit_coincidence_loss(Bimage* prad, long nimg, double& power, double& supr, double& pxsig)
{
	Bimage*			pradpow = prad->next;
	long			i;
	double			pmin(1e30);
	vector<double>	s(prad->sizeX(),0), f(prad->sizeX(),0);
	
	for ( i=0; i<prad->sizeX(); ++i ) {
		s[i] = i*prad->sampling(0)[0];
		f[i] = ((*pradpow)[i] - (*prad)[i]/nimg);
		if ( pmin > f[i] ) pmin = f[i];
//		cout << s[i] << tab << f[i] << endl;
	}
	
//	cout << f.back() << tab << f[1] << endl;
	power = f.back();
	supr = 1 - pmin/power;
	pxsig = 1/prad->real_size()[0];
	
	if ( verbose )
		cout << "Fitting coincidence loss:" << endl;
	if ( verbose & VERB_FULL ) {
		cout << "Starting values:" << endl;
		cout << "Power:        " << power << endl;
		cout << "Suppression:  " << supr << endl;
		cout << "Pixel width:  " << pxsig << endl;
	}
	
	Bsimplex			simp(1, 3, 0, s.size(), s, f);

	simp.parameter(0, power);
	simp.parameter(1, supr);
	simp.parameter(2, pxsig);
	simp.limits(0, power/2, 2*power);
	simp.limits(1, 0, 0.9);
	simp.limits(2, 0.1*pxsig, 10*pxsig);
	
	double			R = simp.run(10000, 1e-5, coincidence_loss_R);
	
	power = simp.parameter(0);
	supr = simp.parameter(1);
	pxsig = simp.parameter(2);
	
	if ( verbose ) {
		cout << "Power:                         " << power << endl;
		cout << "Suppression:                   " << supr << endl;
		cout << "Pixel width:                   " << pxsig << " Å" << endl;
		cout << "R:                             " << R << endl;
	}

	return R;
}

/*
@brief	Calculates the spectral SNR from complex sum and summed power images.
@param 	nimg			number of images contributing to the summed image.
@param 	res_hi			high resolution limit.
@param 	sampling_ratio	radial sampling ratio (1 or larger).
@return Bplot* 			spectral SNR.

	The sub-image shifts must be encoded in the origins as translation from the center.
**/
Bplot*		Bimage::fspace_ssnr(long nimg, double res_hi, double sampling_ratio)
{
	if ( !next ) {
		cerr << "Error: A summed power image must be linked!" << endl;
		bexit(-1);
	}

	if ( res_hi < image->sampling()[0] * 2 ) res_hi = image->sampling()[0] * 2;

	if ( verbose ) {
		cout << "Calculating the SSNR:" << endl;
		cout << "Images:                        " << nimg << endl;
		cout << "High resolution limit:         " << res_hi << " A" << endl;
		cout << "Sampling ratio:                " << sampling_ratio << endl;
	}

	long			i;
	double			s, v, v2, snr, snr_min(0), snr_max(0), pow_min(0), pow_max(0);
	double			rad_scale(real_size()[0]/sampling_ratio);

	Bimage*			prad = fspace_radial_power(res_hi, sampling_ratio);
	Bimage*			pradpow = prad->next = next->fspace_radial_power(res_hi, sampling_ratio);

	double			off(0), amp(0), sig(0), power(0), supr(0), pxsig(0), sf;
	
	double			R = fit_mtf(prad, nimg, off, amp, sig);
	
	if ( amp < 0 ) R = fit_coincidence_loss(prad, nimg, power, supr, pxsig);

	long			maxrad(prad->sizeX());

	long			ncol(5);
	Bstring			title("Spectral signal-to-noise ratio");
	Bplot*			plot = new Bplot(1, maxrad, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("s(1/A)");
	plot->page(0).column(1).label("Signal");
	plot->page(0).column(2).label("Noise");
	plot->page(0).column(3).label("SSNR");
	plot->page(0).column(4).label("Nfit");
	plot->page(0).column(1).type(2);
	plot->page(0).column(2).type(2);
	plot->page(0).column(3).type(2);
	plot->page(0).column(4).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).axis(3);
	plot->page(0).column(3).axis(4);
	plot->page(0).column(4).axis(3);
	plot->page(0).column(1).color(0,0,1);
	plot->page(0).column(2).color(0.5,0.5,0.5);
	plot->page(0).column(3).color(1,0,0);
	plot->page(0).column(4).color(1,0.5,0);
	plot->page(0).axis(1).label("Spatial Frequency (1/A)");
	plot->page(0).axis(3).label("Power");
	plot->page(0).axis(4).label("SSNR");
	plot->page(0).axis(3).color(0,0,0);

	if ( verbose )
		cout << "Shell\ts\tSignal\tNoise\tSNR\tNfit" << setprecision(3) << endl;

	for ( i=1; i<maxrad; ++i ) {
		s = i/rad_scale;
		v = ((*prad)[i] - (*pradpow)[i]);
		v2 = ((*pradpow)[i] - (*prad)[i]/nimg);
		snr = v/v2;
		if ( snr < -1 || !isfinite(snr) ) snr = -1;
		sf = sin(pxsig*s)/(pxsig*s);
		(*plot)[i] = s;
		(*plot)[i+maxrad] = v;
		(*plot)[i+2*maxrad] = v2;
		(*plot)[i+3*maxrad] = snr;
		if ( amp < 0 ) (*plot)[i+4*maxrad] = power*(1-supr*sf*sf);
		else (*plot)[i+4*maxrad] = off + amp*exp(-0.5*s*s/(sig*sig));
		if ( s > 0.05 ) {
			if ( pow_min > v ) pow_min = v;
			if ( pow_max < v ) pow_max = v;
			if ( pow_min > v2 ) pow_min = v2;
			if ( pow_max < v2 ) pow_max = v2;
			if ( snr_min > snr ) snr_min = snr;
			if ( snr_max < snr ) snr_max = snr;
		}
		if ( verbose )
			cout << i << tab << s << tab << v << tab << v2 << tab << snr << tab << (*plot)[i+4*maxrad] << endl;
	}

	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(1.0/res_hi);
	plot->page(0).axis(3).min(pow_min);
	plot->page(0).axis(3).max(1.1*pow_max);
	plot->page(0).axis(4).min(snr_min);
	plot->page(0).axis(4).max(1.1*snr_max);

	Bstring		text(power, "Noise power: %g e/px");
	if ( amp < 0 ) {
		plot->page(0).add_text(text);
		text = Bstring(supr, "Suppression: %g");
		plot->page(0).add_text(text);
		text = Bstring(pxsig/(2*sampling(0)[0]), "Pixel kernel: %g px");
		plot->page(0).add_text(text);
	} else {
		text = Bstring(off, "Power minimum: %g e/px");
		plot->page(0).add_text(text);
		text = Bstring(amp, "Power amplitude: %g e/px");
		plot->page(0).add_text(text);
		text = Bstring(1/(TWOPI*sig*sampling(0)[0]), "Sigma: %g px");
		plot->page(0).add_text(text);
	}
	text = Bstring(R, "R: %g");
	plot->page(0).add_text(text);

	delete prad;

	return plot;
}

/*
@brief	Calculates the spectral SNR from a set of complex sum and summed power images.
@param	subset			number of images summed in each subset.
@param 	res_hi			high resolution limit.
@param 	sampling_ratio	radial sampling ratio (1 or larger).
@param	flag			0=subset sums, 1=progressive sums, 2=no noise.
@return Bplot* 			spectral SNR.

	The sub-image shifts must be encoded in the origins as translation from the center.
**/
Bplot*		Bimage::fspace_subset_ssnr(int subset, double res_hi, double sampling_ratio, int flag)
{
	if ( !next ) {
		cerr << "Error in fspace_subset_ssnr: A summed power image must be linked!" << endl;
		bexit(-1);
	}

	if ( res_hi < image->sampling()[0] * 2 ) res_hi = image->sampling()[0] * 2;

	long			i, j, nn;
	double			s, v, v2, snr, snr_min(0), snr_max(0);
	long			maxrad = fspace_maximum_radius(res_hi, sampling_ratio);
	double			rad_scale(real_size()[0]/sampling_ratio);

	double			dose_per_frame(1);
	if ( metadata.exists("dose") ) dose_per_frame = metadata["dose"].real()/n;

	Bimage*			prad = fspace_radial_power(res_hi, sampling_ratio);
	Bimage*			pradpow = next->fspace_radial_power(res_hi, sampling_ratio);

	Bstring			title, dose;
	if ( flag ) title = "Progressive spectral signal-to-noise ratio";
	else title = "Subset spectral signal-to-noise ratio";

	long			ncol(1+n);
	RGB<float>		color;
	Bplot*			plot = new Bplot(1, maxrad, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("s(1/A)");
	plot->page(0).column(0).axis(1);
	for ( i=1; i<ncol; i++ ) {
		dose = Bstring(dose_per_frame*i, "%g");
		plot->page(0).column(i).label(dose);
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).axis(3);
		color.spectrum(i,1,ncol-1);
		plot->page(0).column(i).color(color.r(),color.g(),color.b());
	}

	if ( verbose ) {
		if ( flag & 2 )
			cout << "Shell\ts\tSignal" << setprecision(3) << endl;
		else
			cout << "Shell\ts\tSignal\tNoise\tSNR" << setprecision(3) << endl;
	}
	for ( i=1; i<maxrad; ++i ) {
		s = i/rad_scale;
		(*plot)[i] = s;
		for ( nn=0, j=i+nn*maxrad; nn<n; ++nn, j+=maxrad ) {
			if ( flag & 2 ) {	// no noise
				if ( flag & 1 ) snr = (*prad)[j]/((nn+1)*subset);
				else snr = (*prad)[j]/subset;
			} else {
				v = (*prad)[j] - (*pradpow)[j];
				if ( flag & 1 ) v2 = (*pradpow)[j] - (*prad)[j]/((nn+1)*subset);
				else v2 = (*pradpow)[j] - (*prad)[j]/subset;
				snr = v/v2;
			}
//			if ( snr < 0 || !isfinite(snr) ) snr = 0;
			if ( !isfinite(snr) ) snr = 0;
			(*plot)[j+maxrad] = snr;
			if ( s > 0.05 ) {
				if ( snr_min > snr ) snr_min = snr;
				if ( snr_max < snr ) snr_max = snr;
			}
		}
		if ( verbose ) {
			if ( flag & 2 )
				cout << i << tab << s << tab << snr << endl;
			else
				cout << i << tab << s << tab << v << tab << v2 << tab << snr << endl;
		}
	}
	
	plot->page(0).axis(1).label("Spatial Frequency (1/A)");
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(1.0/res_hi);
	plot->page(0).axis(3).min(snr_min);
	plot->page(0).axis(3).max(1.1*snr_max);
	plot->page(0).axis(3).label("SSNR");
	
	delete prad;
	delete pradpow;

	return plot;
}

/*
@brief	Real space cross-correlation of two power spectrum annuli in two directions.
@param 	n			number of values in the annulus.
@param 	*x			first annulus.
@param 	*y			second annulus.
@return double* 	array of correlation coefficients.

	The first halves of the annuli is cross-correlated forward and the second
	halves backward, representing the correlations of the original annulus and
	its inverse.
**/
double*		cc_two_way(long nn, double* x, double* y)
{
	long		i, j, k, n2(nn/2);
	double		sx(0), sx2(0), sy(0), sy2(0), sxy, d;
	
	double*		cc = new double[nn];
	for ( i=0; i<nn; i++ ) cc[i] = 0;

	for ( i=0; i<n2; i++ ) {
		sx += x[i];
		sx2 += x[i]*x[i];
		sy += y[i];
		sy2 += y[i]*y[i];
	}
	
	d = (n2*sx2 - sx*sx)*(n2*sy2 - sy*sy);
	if ( d > 0 ) d = sqrt(d);
	else {
		cerr << "Error in cc_two_way: denominator is negative!" << endl;
		return cc;
	}
	
	for ( i=0; i<n2; i++ ) {
		sxy = 0;
		for ( j=0, k=i; j<n2; j++, k++ ) {
			if ( k >= nn ) k = 0;
			sxy += x[j]*y[k];
		}
		cc[i] = (n2*sxy - sx*sy)/d;
	}
	
	for ( i=n2; i<nn; i++ ) {
		sxy = 0;
		for ( j=n2, k=i; j<nn; j++, k-- ) {
			if ( k < 0 ) k = nn - 1;
			sxy += x[j]*y[k];
		}
		cc[i] = (n2*sxy - sx*sy)/d;
	}
	
	return cc;
}

/**
@brief 	Correlate annuli of a polar image.
@param 	*polref		reference image polar power spectrum transform.
@param 	ann_min 	low resolution cutoff (in inverse pixels).
@param 	ann_max 	high resolution cutoff (in inverse pixels).
@param 	ang_min 	minimum angle to test for (radians).
@param 	ang_max 	maximum angle to test for (radians).
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@param 	&cc_max		pointer to correlation coefficient for return.
@return double 		angle.

	The input images must be transformed polar images, with each annulus
	corresponding to a row in the image.
	The polar image annuli between the given limits are cross-correlated
	one by one, and a sum of all annuli calculated, weighted by annulus index.
	The maximum value in the annulus sum corresponds to the best fit angle.
	Both image and mirror image are assessed, yielding a negative angle
	when the mirror image gives a larger correlation.
	FFTW library (www.fftw.org).

**/
double		Bimage::correlate_annuli(Bimage* polref,
				int ann_min, int ann_max, double ang_min, double ang_max, 
				fft_plan planf, fft_plan planb, double& cc_max)
{
	if ( ann_min > ann_max ) swap(ann_min, ann_max);
	if ( ann_min < 0 ) ann_min = 0;
	if ( ann_max >= y ) ann_max = y - 1;
	
	while ( ang_min < -TWOPI ) ang_min += TWOPI;
	while ( ang_max > TWOPI ) ang_max -= TWOPI;
	
	if ( fabs(ang_min - ang_max) < 1e-10 ) {
		ang_min = 0;
		ang_max = TWOPI;
	}
	
	Bimage* 		pcc = polref->pack_two_in_complex(this);
	
	if ( verbose & VERB_PROCESS )
		cout << "Doing row cross-correlations from annulus " << ann_min << " to annulus " << ann_max << endl;
	
    // Do FFT's in place
    long     		i, j, ij, icc;
	long			offset, linesize(pcc->channels()*pcc->sizeX());
	double			scale, sum1, sum2;
	Complex<float>	d1, d2;
	
	float*			cc = new float[2*pcc->sizeX()];
	for ( i=0; i<2*pcc->sizeX(); i++ ) cc[i] = 0;
	
    // Specify FFTW arrays
    Complex<float>*	data;
    Complex<float>*	temp1 = new Complex<float>[pcc->sizeX()];
    Complex<float>*	temp2 = new Complex<float>[pcc->sizeX()];
	
	// Correlate polar images line by line
	for ( i=ann_min, offset=i*linesize; i<=ann_max; i++, offset+=linesize ) {
		data = (Complex<float> *) (pcc->data_pointer(offset));
		fftw(planf, data);
		sum1 = sum2 = 0;
		for ( j=0; j<pcc->sizeX(); j++ ) temp1[j] = 0;
		for ( j=0; j<pcc->sizeX()/2; j++ ) {
			ij = -j;
			if ( ij < 0 ) ij += pcc->sizeX();
			d1 = data[j].unpack_first(data[ij]);
			d2 = data[j].unpack_second(data[ij]);
			temp1[ij] = d1 * d2.conj();
			temp1[i] = temp1[ij].conj();
			sum1 += d1.power();
			sum2 += d2.power();
		}
		if ( sum1 && sum2 ) {
			scale = 0.5/sqrt(sum1*sum2);
			for ( j=0; j<pcc->sizeX(); j++ ) temp1[j] *= scale;
		}
		fftw(planb, temp1);
		for ( j=0; j<pcc->sizeX(); j++ ) cc[j] += temp1[j].real();
    }
	
	delete[] temp1;
	delete[] temp2;

	icc = 0;
	cc_max = cc[0];
	
	long			imin((ang_min/TWOPI)*pcc->sizeX());
	long			imax((ang_max/TWOPI)*pcc->sizeX());
	if ( imin < 0 && imax < 0 ) {
		imin += pcc->sizeX();
		imax += pcc->sizeX();
	}
	if ( imax > pcc->sizeX() ) imax = pcc->sizeX();
	
	if ( imin < 0 && imax > 0 ) {
		imin += pcc->sizeX();
		for ( j=0; j<imin; j++ ) {
			if ( cc_max < cc[j] ) {
				cc_max = cc[j];
				icc = j;
			}
		}
		for ( j=imax; j<pcc->sizeX(); j++ ) {
			if ( cc_max < cc[j] ) {
				cc_max = cc[j];
				icc = j;
			}
		}
	} else {
		if ( imin < 0 ) imin = 0;
		for ( j=imin; j<imax; j++ ) {
			if ( cc_max < cc[j] ) {
				cc_max = cc[j];
				icc = j;
			}
		}
	}
	double			angle = icc*M_PI*(2.0/pcc->sizeX());
//	if ( icc >= pcc->sizeX() ) angle = ((long)pcc->sizeX() - icc)*M_PI*(2.0/pcc->sizeX());
//	if ( icc >= pcc->sizeX() ) angle = (icc - 2.0*pcc->sizeX())*M_PI*(2.0/pcc->sizeX());
	
	if ( fabs(angle) < 1e-10 ) angle = 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::correlate_annuli: " << icc << tab << angle << tab << angle*180.0/M_PI << tab <<  cc[icc] << endl;
	
	delete[] cc;

	delete pcc;
	
	if ( verbose & VERB_PROCESS )
		cout << "Best in-plane rotation angle:   " << angle*180.0/M_PI << endl;
	
	return angle;
}

/*
@brief	Calculates the sums of the cross-correlation arrays for annuli of the power spectra.
@param 	*pref			reference 2D image.
@param 	res_hi			high resolution limit.
@param 	res_lo			low resolution limit.
@param 	nang			number of angles.
@param 	planf			FFT forward plan.
@return vector<double>	array of summed correlation coefficients.

	The images are packed into a complex form, Fourier transformed and unpacked.
	The two polar power spectra are calculated and the annuli cross-correlated.
	The input images must be equal-sized square 2D images.
**/
vector<double>	Bimage::pps_angular_correlation(Bimage* pref,
				double res_hi, double res_lo, long nang, fft_plan planf)
{
	Bimage*				pc = pack_two_in_complex(pref);

	pc->fft(planf, 1);
	
	pc->origin(0.0,0.0,0.0);
	
	Bimage*				pc2 = pc->unpack_combined_transform();
	
	pc->complex_to_amplitudes();
	pc2->complex_to_amplitudes();
	
	Bimage*				ppsref = pc2->polar_power_spectrum(res_hi, nang);
	Bimage*				pps = pc->polar_power_spectrum(res_hi, nang);

	delete pc;
	delete pc2;
	
	long				min_rad(1);
	if ( res_lo ) min_rad = (long) (real_size()[0]/res_lo);
	if ( min_rad < 1 ) min_rad = 1;
	
    long				i, j, na, r, m, offset;
	double				f;
	double*				cc;
	vector<double>		ccs(pps->sizeX(), 0);

	for ( r=min_rad, m=0, offset=r*pps->sizeX(); r<pps->sizeY(); r++, m++, offset+=pps->sizeX() ) {
		na = nang*r;
		cc = cc_two_way(na, (double*) (ppsref->data_pointer(offset)), (double*) (pps->data_pointer(offset)));
		f = na*1.0L/pps->sizeX();
		for ( i=0; i<pps->sizeX(); i++ ) {
			j = (long) (i*f + 0.5);
			if ( j >= na ) j -= na;
			ccs[i] += cc[j];
		}
		delete[] cc;
	}
	
	for ( i=0; i<pps->sizeX(); i++ ) ccs[i] /= m;

	delete pps;
	delete ppsref;
	
	return ccs;
}

/**
@brief 	Finds the best in-plane alignment for two 2D images using polar power spectra.
@param 	*pref		reference 2D image.
@param 	res_hi		high resolution limit.
@param 	res_lo		low resolution limit.
@param 	shift_limit	maximum shift from nominal origin of box.
@param 	angle_limit	maximum rotation from original in-plane rotation angle.
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@return double 		correlation coefficient.

	The rotation angle is found by cross-correlating the annuli
	of the polar transforms of the power spectra of the image and the
	reference image. If the full asu flag is not set, the annular cross-correlation
	is done both forwards and backwards to detect mirror images, with the sign 
	of the angle indicating which direction is best (>0 => forwards, <0 => backwards).
	The resolution limit specified in the image is used to low-pass filter 
	the cross-correlation. The shift is then determined by cross-correlation 
	after rotating the reference image, as well as its 180° rotation to deal with
	that ambiguity in the power spectrum.
	The correlation coefficient return is the cross-correlation peak.
	The input images must be equal-sized square 2D images.
	The origin and rotation angle is returned in the image structure.

**/
double		Bimage::align2D_pps(Bimage* pref, double res_hi, double res_lo, 
				double shift_limit, double angle_limit, 
				fft_plan planf, fft_plan planb)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::align2D_pps: res_hi=" << res_hi << " res_lo=" << res_lo << endl;

	long			nang(16);
	long			max_rad = (long) (real_size()[0]/res_hi);
	if ( max_rad > x/2 ) max_rad = x/2;
	long			max_ang(nang*max_rad);
    long			hang(max_ang/2);

	double			angle(image->view_angle() - atan2(image->view()[1], image->view()[0]));
//	double			angle(image->view_angle());
	long			imin(0), imax(hang);
	if ( angle_limit ) {
		imin = max_ang*(angle - angle_limit)/TWOPI;
		imax = max_ang*(angle + angle_limit)/TWOPI;
		if ( imax < 0 ) {
			imin += hang;
			imax += hang;
		}
	}
//	cout << "View: " << image->view() << endl;
//	cout << "Angle limits: " << (angle - angle_limit)*180.0/M_PI << tab << (angle + angle_limit)*180.0/M_PI << endl;
	
	long			i, j, shift(0);
	double			cc1, cc2, cc3, cc, f;
	
	// Forward cross-correlation
	vector<double>	ccs = pps_angular_correlation(pref, res_hi, res_lo, nang, planf);
	if ( imin >= 0 ) {
		for ( i=imin, cc2 = 0; i<imax; i++ ) if ( cc2 < ccs[i] ) {
			shift = i;
			cc2 = ccs[i];
		}
	} else {
		for ( i=0, cc2 = 0; i<imax; i++ ) if ( cc2 < ccs[i] ) {
			shift = i;
			cc2 = ccs[i];
		}
		for ( i=imin; i<hang; i++ ) if ( cc2 < ccs[i] ) {
			shift = i;
			cc2 = ccs[i];
		}
	}
	j = shift - 1;
	if ( j < 0 ) j += hang;
	cc1 = ccs[j];
	j = shift + 1;
	if ( j >= hang ) j -= hang;
	cc3 = ccs[j];
	f = (cc1 - cc3)/(2*(cc1 - 2*cc2 + cc3));
	angle = (shift + f)*M_PI*2.0L/max_ang;
	
	// Determine origins and CC's for both angles and select the best
	cc = cc1 = rotate_cross_correlate_two_way(pref, angle, res_hi, res_lo, shift_limit, planf, planb);

	return cc;
}


/**
@brief 	Finds the best in-plane alignment for two 2D images using polar power spectra.
@param 	*pref		reference 2D image.
@param 	res_polar	polar resolution limit.
@param 	ann_min 	minimum annulus (>=0).
@param 	ann_max 	maximum annulus (< image radius).
@param 	*prs_mask	dual mask.
@param 	shift_limit	maximum shift from nominal origin of box.
@param 	angle_limit	maximum rotation from original in-plane rotation angle.
@param 	planf_1D	FFT forward plan for polar images.
@param 	planb_1D	FFT backward plan for polar images.
@param 	planf_2D	FFT forward plan for 2D images.
@param 	planb_2D	FFT backward plan for 2D images.
@return double		correlation coefficient.

	Both the image and reference is converted to polar images using the 
	current origins. The annuli of the polar images are cross-correlated
	to find the rotation angle. The reference is rotated and cross-correlated
	with the image to determine a new origin for the image. This is iterated
	untill the rotation angle and origin of the image does not change any more,
	or a maximum number of iterations.
	The cross-correlation is done with the provided mask..
	The correlation coefficient return is the cross-correlation peak.
	The input images must be equal-sized square 2D images.

**/
double		Bimage::align2D(Bimage* pref, double res_polar,
				int ann_min, int ann_max, Bimage* prs_mask, double shift_limit, double angle_limit, 
				fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::align2D: ann_min=" << ann_min << " ann_max=" << ann_max << endl;

	double			angle(image->view_angle() - atan2(image->view()[1], image->view()[0]));
//	double			angle(image->view_angle());
	double			ang_min(0), ang_max(TWOPI);
	if ( angle_limit ) {
		ang_min = angle - angle_limit;
		ang_max = angle + angle_limit;
	}
//	cout << "View: " << image->view() << endl;
//	cout << "Angle limits: " << ang_min*180.0/M_PI << tab << ang_max*180.0/M_PI << endl;
	
	long 			nannuli(pref->sizeX()/2);	// Sampling same as image
	long 			nangles(NPOLANG);		// 0.5 degree step size
	double			pimgcc(0);				// Polar image correlation coefficient

    Bimage*         ppart = copy();
    Bimage*         pproj = pref->copy();
    
    if ( res_polar ) {
        ppart->fspace_bandpass(0, res_polar, 0, planf_2D, planb_2D);
        pproj->fspace_bandpass(0, res_polar, 0, planf_2D, planb_2D);
    }
	
//	cout << ppart->image->origin() << tab << pproj->image->origin() << endl;
	
	Bimage*			pol = NULL;
	Bimage*			polref = pproj->cartesian_to_cylindrical(nannuli, nangles, 1);
    delete pproj;
	
	long			i, done(0);
	double			best_angle(0);
	double			da, dx, dy, ap(1e37);
	double			cc(0), best_cc(-1e37);
	Vector3<double> translate, prev_translate(1e37, 1e37, 0), best_shift;
	Bimage*			prot;

	for ( i=0; i<5 && !done; i++ ) {
		pol = ppart->cartesian_to_cylindrical(nannuli, nangles, 1);
		angle = pol->correlate_annuli(polref, ann_min, ann_max, 
			ang_min, ang_max, planf_1D, planb_1D, pimgcc);
		delete pol;
		prot = pref->rotate(pref->size(), fabs(angle));
		translate = find_shift(prot, prs_mask, 0, 0, shift_limit, 0, 1, planf_2D, planb_2D, cc);
		delete prot;
		if ( best_cc < cc ) {
			best_cc = cc;
			best_shift = translate;
			best_angle = angle;
		}
		da = fabs(angle - ap);
		dx = fabs(translate[0] - prev_translate[0]);
		dy = fabs(translate[1] - prev_translate[1]);
		prev_translate = translate;
		ap = angle;
		if ( i > 0 && da < 0.005 && dx < 0.1 && dy < 0.1 ) done = 1;
	}
	
	delete polref;
    delete ppart;
	
	image->origin(pref->image->origin() + best_shift);
	image->view_angle(best_angle);
	image->FOM(best_cc);
//	image->view_angle(pref->image->view_angle() + best_angle);
	
	return best_cc;
}

/**
@brief 	Finds the best in-plane alignment for two 2D images using polar power spectra.
@param 	*pref		reference 2D image.
@param 	ann_min 	minimum annulus (>=0).
@param 	ann_max 	maximum annulus (< image radius).
@param 	res_lo		low resolution limit (angstrom).
@param 	res_hi		high resolution limit (angstrom).
@param 	shift_limit	maximum shift from nominal origin of box.
@param 	angle_limit	maximum rotation from original in-plane rotation angle.
@return double		correlation coefficient.

	This is a simplified interface to the other function Bimage::align2D.
	All the fourier transform and rotation functions are handled internally.

**/
int			Bimage::align2D(Bimage* pref, int ann_min, int ann_max,
				double res_lo, double res_hi, double shift_limit, double angle_limit)
{
 	// Generate the masks for cross-correlation and cross-validation
	Bimage*		prs_mask = new Bimage(Float, TSimple, size(), 1);
	prs_mask->sampling(sampling(0));
	
	vector<double>	band = fspace_default_bands(res_lo, res_hi);
	prs_mask->mask_fspace_banded(band);

   	// Specify FFTW arrays
	fft_plan		planf_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_FORWARD, 1);
	fft_plan		planb_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_BACKWARD, 1);	
	fft_plan		planf_2D = fft_setup_plan(size(), FFTW_FORWARD, 1);
	fft_plan		planb_2D = fft_setup_plan(size(), FFTW_BACKWARD, 1);


	double			cc = align2D(pref, 0, ann_min, ann_max, prs_mask, shift_limit,
							angle_limit, planf_1D, planb_1D, planf_2D, planb_2D);

	if ( verbose )
		cout << "Alignment correlation:          " << cc << endl << endl;
	
    fft_destroy_plan(planf_1D);
    fft_destroy_plan(planb_1D);
    fft_destroy_plan(planf_2D);
    fft_destroy_plan(planb_2D);
	
	delete prs_mask;
	
	rotate();
	
	return 0;
}

