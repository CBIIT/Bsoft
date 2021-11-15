/**
@file	mg_pick.cpp
@brief	Functions for picking single particle images from a micrograph.
@author Bernard Heymann
@date	Created: 20000505
@date	Modified: 20210510
**/

#include "mg_pick.h"
#include "mg_processing.h"
#include "mg_select.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* flags	1: quadric correct; 2: gaussian filter; 4: mask background */
Bimage*		img_prepare_for_picking(Bimage* p, int flags, long sigma, long bin)
{
	Bimage*			pb = p->bin_copy(bin);
	
	if ( flags & 1 ) pb->quadric_correct(pb->quadric_fit());

	long			ksize(9);
	if ( flags & 2 ) {
		sigma /= bin;
		if ( sigma < 1 ) sigma = 1;
		ksize = 2*( (long) (4*sigma + 0.9) ) + 1;
		pb->filter_gaussian(ksize, sigma);
	}
	
	pb->rescale_to_avg_std(0, 1);
	pb->statistics();

	double 			threshold(pb->average() - 3*pb->standard_deviation());
	Bimage*			pmask = NULL;
	
	if ( verbose & VERB_FULL )
		cout << "threshold = " << threshold << endl;
	
	if ( flags & 4 ) {
		pmask = pb->mask_by_threshold(threshold);
		pmask->mask_erode(5);
		pmask->change_type(Float);
		pmask->filter_average(11);
		pb->multiply(pmask, 1, 0);
		delete pmask;
	}
		
	return pb;
}


/**
@brief 	Finds the peaks in a cross-correlation map corresponding to particles.
@param 	*pcc			peak map (after binning).
@param	bin				binning to speed up calculations.
@param 	excl_dist		distance between peaks.
@param 	part_ori		particle origin.
@param 	fommin			minimum threshold to accept peaks.
@param 	fommax			maximum threshold to accept peaks.
@param 	maxnum			maximum number of peaks to pick.
@param 	pix_min			minimum peak width.
@param 	pix_max			maximum peak width.
@return Bparticle*		list of particles.

	The map is searched in increments of the particle radius to identify
	peaks above the threshold and within a box the size of the
	particle radius. The identified peaks are further examined to eliminate 
	ones that are too close to a higher scoring peak. The acceptable distance
	between peaks is set to 1.8 times the particle radius.

**/
Bparticle*	particles_from_peaks(Bimage* pcc, long bin, double excl_dist, 
				double part_ori, double& fommin, double fommax, 
				long maxnum, double pix_min, double pix_max)
{
	if ( part_ori < excl_dist ) part_ori = excl_dist;
	
	long			i, j, ncoor;
	double			fom;
	Vector3<double>	origin(part_ori, part_ori, part_ori);
	origin = origin.min(pcc->size()-1);
	Vector3<long>	min(origin/bin), max(pcc->size()-min);
	
	if ( verbose & VERB_PROCESS )
		cout << "Limits: " << min << tab << max << endl;
	
	Vector3<long>	obin(bin, bin, bin);
	obin = obin.min(pcc->size());

	if ( excl_dist < 10 ) excl_dist = 1.8*part_ori;
	
	excl_dist /= bin;
	
	Vector3<double>*	coor = pcc->find_peaks(excl_dist, ncoor, fommin, fommax, pix_min, pix_max);
	
	Bparticle*		partlist = NULL;
	Bparticle*		part = NULL;
	
	for ( i=j=0; i<ncoor && j<maxnum; i++ ) {
		fom = pcc->get(0, coor[i]);
		if ( fom >= fommin && fom <= fommax ) {
			if ( coor[i] >= min && coor[i] <= max ) {
				part = particle_add(&part, ++j);
				if ( !partlist ) partlist = part;
				part->loc = coor[i] * obin;
				part->ori = origin;
				part->fom[0] = fom;
			}
		}
	}
	
	delete[] coor;

	if ( verbose & VERB_PROCESS )
		cout << "Number of particles found:      " << j << endl << endl;
	
	return partlist;
}


/**
@brief 	Picks particles using cross-correlation.
@param 	*p				image to pick from.
@param 	*ptemp			template image.
@param 	*pmask			frequency space mask.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	fommin			minimum FOM cutoff.
@param 	fommax			maximum FOM cutoff.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return Bparticle*		list of particles.

	A template is cross-correlated with the input image including
	bandpass filtering to target the size of the particle.

**/
Bparticle*	particles_pick_cc(Bimage* p, Bimage* ptemp, Bimage* pmask, 
				double hires, double lores, double fommin, double fommax, double excl_dist, long bin)
{
	if ( !p ) return NULL;
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	if ( excl_dist < 1 ) excl_dist = 0.8*ptemp->sizeX();
	
	if ( verbose & VERB_FULL )
		cout << "Picking from " << p->file_name() << endl;

	long 			sigma(10);
	Bimage*			pb = img_prepare_for_picking(p, 3, sigma, bin);
	
	Bimage*			pt = ptemp->bin_copy(bin);
	
	fft_plan		planf = pb->fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = pb->fft_setup(FFTW_BACKWARD, 0);

	Bimage*			pcc = pb->find_template(pt, pmask, hires, lores, 1, planf, planb);
	
	fft_destroy_plan(planf);
	fft_destroy_plan(planb);
	delete pb;
	
	ptemp->image->FOM(pcc->image->FOM());
	
	if ( verbose & VERB_DEBUG )		
		write_img("pcc.pif", pcc, 0);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particles_pick: bin=" << bin << endl;
	
	Bparticle*		partlist = particles_from_peaks(pcc, bin, excl_dist, 
				ptemp->image->origin()[0], fommin, fommax);

	delete pcc;
	
	return partlist;
}

double		gauss_find_intersection(long ngauss, vector<double>& gauss)
{
	double		a1, a2, v, d(1), l((gauss[1]+gauss[2])/(gauss[0]+gauss[3]));
	double		cut((gauss[1] + gauss[4])/2);
	
	while ( d > 1e-6 ) {
		v = (cut - gauss[1])/gauss[2];
		a1 = gauss[0]*exp(-0.5*v*v);
		v = (cut - gauss[4])/gauss[5];
		a2 = gauss[3]*exp(-0.5*v*v);
		d = a2 - a1;
		cut -= l*d;
		l *= 0.99;
		if ( verbose & VERB_DEBUG )
			cout << a1 << tab << a2 << tab << cut << endl;
		d = fabs(d);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "u1=" << gauss[1] << tab << "u2=" << gauss[4] << tab << "cut=" << cut << endl;
	
	return cut;
}

/**
@brief 	Picks particles in variance map.
@param 	*p				image to pick from.
@param 	avg_kernel		averaging kernel size.
@param 	var_kernel		variance kernel size.
@param 	nsig			multiple of sigma above variance average to accept peaks.
@param 	part_ori		particle origin.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return Bparticle*		list of particles.

	A copy of the micrograph is filtered with an averaging kernel
	and a variance map calculated. The variance map is then used to
	find high variance peaks as candidate locations for particles.

**/
Bparticle*	particles_pick_var(Bimage* p, long avg_kernel, long var_kernel,
				double nsig, double part_ori, double excl_dist, long bin)
{
	if ( !p ) return NULL;
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	if ( excl_dist < 1 ) excl_dist = 1.6*part_ori;
	
	if ( verbose & VERB_FULL )
		cout << "Picking from " << p->file_name() << endl;
	
	avg_kernel /= bin;
	var_kernel /= bin;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particles_pick_var: nsig=" << nsig << endl;

	Bimage*			pc = p->bin_copy(bin);
	pc->filter_average(avg_kernel);
	pc->variance(var_kernel);
	
	// Two gaussians are fitted to the variance image histogram
	vector<double>	gauss = pc->histogram_gauss_fit(256, 2);
	double			fomcut, t;
	
	t = fabs(gauss[1]-gauss[4])/sqrt(gauss[2]*gauss[2]+gauss[5]*gauss[5]);
	
//	fomcut = gauss_find_intersection(2, gauss);
/*	
	if ( gauss[1] < gauss[4] ) {
		fomcut = gauss[1] + nsig*gauss[2];
		if ( t > 1 && fomcut > gauss[4] ) fomcut = gauss_find_intersection(2, gauss);
	} else {
		fomcut = gauss[4] + nsig*gauss[5];
		if ( t > 1 && fomcut > gauss[1] ) fomcut = gauss_find_intersection(2, gauss);
	}
*/	
	
	if ( gauss[1] < gauss[4] ) {
		fomcut = gauss[1] + nsig*gauss[2];
		if ( fomcut > gauss[4] ) fomcut = gauss[4];
	} else {
		fomcut = gauss[4] + nsig*gauss[5];
		if ( fomcut > gauss[1] ) fomcut = gauss[1];
	}	

	if ( verbose ) {
		cout << "Variance image min, max:        " << pc->minimum() << tab << pc->maximum() << endl;
		cout << "First peak:                     " << gauss[0] << tab << gauss[1] << tab << gauss[2] << endl;
		cout << "Second peak:                    " << gauss[3] << tab << gauss[4] << tab << gauss[5] << endl;
		cout << "T and FOM:                      " << t << tab << fomcut << endl;
	}
	
	if ( verbose & VERB_DEBUG )	{
		cout << "DEBUG particles_pick_var: fomcut=" << fomcut << endl;
		write_img("var.pif", pc, 0);
	}
	
	Bparticle*		partlist = particles_from_peaks(pc, bin, excl_dist, 
						part_ori, fomcut, 1e30, 1000000, 2, excl_dist);
	
	delete pc;
	
	return partlist;
}

Bparticle*	particles_pick_var(Bimage* p, long avg_kernel, long var_kernel,
				double cutmin, double cutmax, double part_ori, double excl_dist, long bin)
{
	if ( !p ) return NULL;
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	if ( excl_dist < 1 ) excl_dist = 1.6*part_ori;
	
	if ( verbose & VERB_FULL )
		cout << "Picking from " << p->file_name() << endl;
	
	avg_kernel /= bin;
	var_kernel /= bin;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particles_pick_var: cutmin=" << cutmin << " cutmax=" << cutmax << endl;

	Bimage*			pc = p->bin_copy(bin);
	pc->filter_average(avg_kernel);
	pc->variance(var_kernel);
	pc->truncate(cutmin, cutmax, 0, 0);
	
	Bparticle*		partlist = particles_from_peaks(pc, bin, excl_dist, 
						part_ori, cutmin, cutmax, 1000000, 2, excl_dist);
	
	delete pc;
	
	return partlist;
}


/**
@brief 	Picks particles using cross-correlation.
@param 	&filename		image to pick from.
@param 	img_num			sub-image number.
@param 	*ptemp			template image.
@param 	*pmask			frequency space mask.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	fommin			minimum FOM cutoff.
@param 	fommax			maximum FOM cutoff.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return Bparticle*		list of particles.

	A template is cross-correlated with the input image including
	bandpass filtering to target the size of the particle.
	The template must have the correct pixel size.

**/
Bparticle*	particles_pick_cc(Bstring& filename, long img_num, Bimage* ptemp, 
				Bimage* pmask, double hires, double lores, double fommin, double fommax, 
				double excl_dist, long bin)
{
	Bimage*				p = read_img(filename, 1, img_num);
	if ( !p ) return NULL;
	
	p->sampling(ptemp->sampling(0));
	
	Bparticle*			partlist = particles_pick_cc(p, ptemp, pmask, 
							hires, lores, fommin, fommax, excl_dist, bin);
	
	delete p;
	
	return partlist;
}

/**
@brief 	Picks particles in variance map.
@param 	&filename		image to pick from.
@param 	img_num			sub-image number.
@param 	avg_kernel		averaging kernel size.
@param 	var_kernel		variance kernel size.
@param 	nsig			multiple of sigma above variance average to accept peaks.
@param 	part_ori		particle origin.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return Bparticle*		list of particles.

	A copy of the micrograph is filtered with an averaging kernel
	and a variance map calculated. The variance map is then used to
	find high variance peaks as candidate locations for particles.

**/
Bparticle*	particles_pick_var(Bstring& filename, long img_num, 
				long avg_kernel, long var_kernel, double nsig, 
				double part_ori, double excl_dist, long bin)
{
	Bimage*				p = read_img(filename, 1, img_num);
	if ( !p ) return NULL;

	Bparticle*			partlist = particles_pick_var(p, avg_kernel, var_kernel,
							nsig, part_ori, excl_dist, bin);
	
	delete p;
	
	return partlist;
}

double		img_marker_set_difference_at_voxel(Bimage* p, long i, Bmarker* markin, Bmarker* markout)
{
	long				x, y, z, n;
	double				v, sumin(0), sumout(0);
	Bmarker*			mark = NULL;
	
	p->coordinates(i, x, y, z);
	Vector3<float>		c(x, y, z);
	
	for ( n=0, mark = markin; mark; mark = mark->next ) {
		v = p->get(0, mark->loc + c, 0);
		if ( v ) {
			sumin += v;
			n++;
		}
	}
	
	if ( n ) sumin /= n;
	
	for ( n=0, mark = markout; mark; mark = mark->next ) {
		v = p->get(0, mark->loc + c, 0);
		if ( v ) {
			sumout += v;
			n++;
		}
	}
	
	if ( n ) sumout /= n;
	
	return sumin - sumout;
}

/**
@brief 	Calculates the average difference between voxels tagged by two marker sets.
@param 	*p			image.
@param 	*markin		inner marker set (relates to particle).
@param 	*markout	outer marker set (relates to background).
@param 	contrast	contrast direction (foreground: white=1, black=0).
@return Bimage* 	map with resultant differences.

	The values in the image corresponding to the coordinates of the
	inner marker set is subtracted from those of the outer marker set.

**/
Bimage*		img_marker_set_difference(Bimage* p, Bmarker* markin, Bmarker* markout, int contrast)
{
	Bimage*				pt = p->copy_header(p->images());
	pt->data_type(Float);
	pt->data_alloc();
	
#ifdef HAVE_GCD
	dispatch_apply(p->data_size(), dispatch_get_global_queue(0, 0), ^(size_t i){
		pt->set(i, img_marker_set_difference_at_voxel(p, i, markin, markout));
	});
#else
#pragma omp parallel for
	for ( long i=0; i<p->data_size(); i++ )
		pt->set(i, img_marker_set_difference_at_voxel(p, i, markin, markout));
#endif

	pt->statistics();

	double			scale(1);
	if ( contrast > 0 ) scale = 1/pt->maximum();
	else scale = 1/pt->minimum();
	
	pt->rescale(scale, 0);
	
	return pt;
}

/**
@brief 	Locates particles using a fore/background difference measure.
@param 	*mg			micrograph parameter structure.
@param 	*p			micrograph image, if NULL, read from mg.
@param 	*markin		foreground marker set.
@param 	*markout	background marker set.
@param 	avg_kernel	averaging kernle to smooth difference map.
@param 	contrast	contrast direction (foreground: white=1, black=0).
@return double		threshold used to accept peaks.

	The two marker sets are used to calculate a difference
	between the foreground and background as an estimate of the
	presence of a particle.
	The resultant image is smoothed by an averaging kernel and rescaled to [0,1].
	This image is then scanned to identify particles.

**/
double		mg_pick_particles(Bmicrograph* mg, Bimage* p, Bmarker* markin,
				Bmarker* markout, int avg_kernel, int contrast)
{
	if ( !mg->select ) return 0;
	 
	if ( mg->pixel_size[0] > 0 ) p->sampling(mg->pixel_size);
	else mg->pixel_size = p->sampling(0);

	mg->origin = p->size()/2;

	p->filter_average(avg_kernel);
	
	Bimage*			pt = img_marker_set_difference(p, markin, markout, contrast);

	if ( verbose & VERB_DEBUG )		
		write_img("map.pif", pt, 0);

//	double			threshold = pt->average() + pt->standard_deviation();
	double			threshold(0);

//	if ( verbose )
//		cout << "Threshold:                      " << threshold << endl;
	
	if ( mg->part ) {
		particle_kill(mg->part);
		mg->part = NULL;
	}

	mg->part = particles_from_peaks(pt, 1, 0, mg->box_size[0]/2, threshold);
	
	delete pt;
	
	mg->fom = threshold;

	Bparticle*		part = NULL;
	for ( part = mg->part; part; part = part->next )
		part->ori = mg->box_size/2;

	micrograph_set_part_links(mg);
	
	return threshold;
}

/**
@brief 	Picks particles using cross-correlation.
@param 	*project		project parameter structure.
@param 	*ptemp			template image.
@param 	*pmask			frequency space mask.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	fommin			minimum FOM cutoff.
@param 	fommax			maximum FOM cutoff.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return double			minimum threshold used to accept peaks.

	Each micrograph is cross-correlated with the template image including
	bandpass filtering and frequency space masking.

**/
double		project_pick_particles(Bproject* project, Bimage* ptemp, Bimage* pmask, 
				double hires, double lores, double fommin, double fommax, double excl_dist, long bin)
{
	long				npart(0);
	double				min_thresh(1);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose ) {
		cout << "Finding particles by cross-correlation:" << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		if ( pmask )
			cout << "Mask:                           " << pmask->file_name() << endl;
		cout << "Resolution limits:              " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Thresholds:                     " << fommin << " - " << fommax << endl;
		cout << "Exclusion distance:             " << excl_dist << endl;
		cout << "Binning:                        " << bin << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
				if ( verbose )
					cout << "Picking from micrograph " << mg->id << flush;
				mg->box_size = ptemp->size();
				mg->part = particles_pick_cc(mg->fmg, mg->img_num, ptemp, pmask, hires, lores, fommin, fommax, excl_dist, bin);
				mg->fom = ptemp->image->FOM();
				npart = micrograph_set_part_links(mg);
				if ( min_thresh > mg->fom ) min_thresh = mg->fom;
				if ( verbose )
					cout << tab << npart << endl;
			}
		}
		npart = project_count_mg_particles(project);
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose )
				cout << "Picking from reconstruction " << rec->id << endl;
			rec->box_size = ptemp->size();
			rec->part = particles_pick_cc(rec->frec, 0, ptemp, pmask, hires, lores, fommin, fommax, excl_dist, bin);
			rec->fom = ptemp->image->FOM();
			if ( min_thresh > rec->fom ) min_thresh = rec->fom;
			npart = reconstruction_set_part_links(rec);
		}
		npart = project_count_rec_particles(project);
	}
	
	if ( verbose )
		cout << "Number of particles:            " << npart << endl << endl;
		
	return min_thresh;
}

/**
@brief 	Picks particles in variance map.
@param 	*project		project parameter structure.
@param 	avg_kernel		averaging kernel size.
@param 	var_kernel		variance kernel size.
@param 	nsig			multiple of sigma above variance average to accept peaks.
@param 	part_ori		particle origin.
@param 	excl_dist		minimum distance between particles.
@param 	bin				level of image binning.
@return double			minimum threshold used to accept peaks.

	A copy of the micrograph is filtered with an averaging kernel
	and a variance map calculated. The variance map is then used to
	find high variance peaks as candidate locations for particles.

**/
double		project_pick_particles(Bproject* project, long avg_kernel, long var_kernel,
				double nsig, double part_ori, double excl_dist, long bin)
{
	if ( part_ori < 1 ) part_ori = excl_dist;
	
	long				npart(0);
	double				min_thresh(1);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose ) {
		cout << "Finding particles from local variance map:" << endl;
		cout << "Averaging kernel size:          " << avg_kernel << endl;
		cout << "Variance kernel size:           " << var_kernel << endl;
		cout << "Multiple of sigma for cutoff:   " << nsig << endl;
		cout << "Exclusion distance:             " << excl_dist << endl;
		cout << "Binning:                        " << bin << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
				if ( verbose )
					cout << "Picking from micrograph " << mg->id << endl;
				if ( mg->box_size.length() < 1 )
					mg->box_size = Vector3<long>(2*part_ori, 2*part_ori, 1);
				mg->part = particles_pick_var(mg->fmg, mg->img_num, 
					avg_kernel, var_kernel, nsig, part_ori, excl_dist, bin);
				if ( min_thresh > mg->fom ) min_thresh = mg->fom;
				micrograph_set_part_links(mg);
			}
		}
		npart = project_count_mg_particles(project);
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose )
				cout << "Picking from reconstruction " << rec->id << endl;
			if ( rec->box_size.length() < 1 )
				rec->box_size = Vector3<long>(2*part_ori, 2*part_ori, 2*part_ori);
			rec->part = particles_pick_var(rec->frec, 0, 
				avg_kernel, var_kernel, nsig, part_ori, excl_dist, bin);
			if ( min_thresh > rec->fom ) min_thresh = rec->fom;
			reconstruction_set_part_links(rec);
		}
		npart = project_count_rec_particles(project);
	}
	
	if ( verbose )
		cout << "Number of particles:            " << npart << endl << endl;
		
	return min_thresh;
}

/**
@brief 	Locates particles using a fore/background difference measure.
@param 	*project		project parameter structure.
@param 	din				inner diameter corresponding to particle edge.
@param 	dout			outer diameter corresponding to background.
@param 	avg_kernel		averaging kernle to smooth difference map.
@param 	ainc			angular increment.
@param 	flags			flags: bit 1 = filter extremes.
@param 	contrast		contrast direction (foreground: white=1, black=0).
@return long			number of particles.

	Two marker sets are generated at the indicated diameters.
	At each pixel, the marker sets are used to calculate a difference
	between the foreground and background as an estimate of the
	presence of a particle.

**/
long		project_pick_particles(Bproject* project, double din, double dout,
				int avg_kernel, double ainc, int flags, int contrast)
{
	if ( din < 2 ) {
		cerr << "Error: The particle diameter must be specified!" << endl << endl;
		return -1;
	}

	int				w = (int) (dout + 10);
	Vector3<long>	size(w, w, 1);

	Bmarker*		markin = marker_set_at_radius(din/2, ainc);
	Bmarker*		markout = marker_set_at_radius(dout/2, ainc);
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bimage*			p;
	
	if ( verbose ) {
		cout << "Finding particles using a circular fore/background template:" << endl;
		cout << "Particle diameter:              " << din << endl;
		cout << "Background diameter:            " << dout << endl;
		cout << "Angular increment:              " << ainc*180.0/M_PI << " degrees" << endl;
		cout << "Averaging kernel size:          " << avg_kernel << endl << endl;
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			mg->box_size = size;
			p = read_img(mg->fmg, 1, mg->img_num);
			if ( p ) {
				if ( verbose )
					cout << "Picking from micrograph " << mg->id << endl;
				if ( flags & 1 ) p->filter_extremes(1);
				mg->fom = mg_pick_particles(mg, p, markin, markout, avg_kernel, contrast);
				micrograph_set_part_links(mg);
				delete p;
			}
		}
	}
	
	kill_list((char *) markin, sizeof(Bmarker));
	kill_list((char *) markout, sizeof(Bmarker));
	
	long			npart = project_count_mg_particles(project);

	if ( verbose )
		cout << "Number of particles:            " << npart << endl << endl;
			
	return npart;
}


int			part_mask_particles(Bparticle* part, Bimage* pmask, long box_size)
{
	if ( verbose ) {
		cout << "Masking the existing particles:" << endl;
		cout << "Box size:                       " << box_size << endl;
	}
	
	for ( ; part; part = part->next )
		pmask->sphere(part->loc, box_size, 0, FILL_USER, 0);
	
	return 0;
}

int			img_mask_above_threshold(Bimage* pvar, Bimage* p, double threshold)
{
	if ( verbose ) {
		cout << "Masking the variance map:" << endl;
		cout << "Threshold:                      " << threshold << endl;
	}
	
	for ( long i=0; i<p->data_size(); i++ ) if ( (*p)[i] > threshold ) pvar->set(i, 0);
	
	return 0;
}

double		mg_pick_background(Bmicrograph* mg, Bimage* p, long number, 
				long avg_kernel, long var_kernel, double excl_dist)
{
	if ( !mg->select ) return 0;

	p->filter_average(avg_kernel);

	if ( verbose )
		cout << "Filtered min, max, avg & std:   " << p->minimum() << " " <<
			p->maximum() << " " << p->average() << " " << p->standard_deviation() << endl;

	Bimage*			pvar = p->copy();
	
	pvar->variance(var_kernel);
	
	pvar->rescale(-1/pvar->standard_deviation(), pvar->average()/pvar->standard_deviation());
	
	if ( verbose )
		cout << "Variance min, max, avg & std:  " << pvar->minimum() << " " <<
			pvar->maximum() << " " << pvar->average() << " " << pvar->standard_deviation() << endl;

	part_mask_particles(mg->part, pvar, mg->box_size[0]);
	
	img_mask_above_threshold(pvar, p, p->average() + p->standard_deviation());
	
//	write_img("vt.pif", pvar);
	
	double			threshold(0);

	if ( mg->part ) {
		particle_kill(mg->part);
		mg->part = NULL;
	}

	mg->part = particles_from_peaks(pvar, 1, excl_dist, mg->box_size[0]/2, threshold, 1e30, number);
	
	delete pvar;

	micrograph_set_part_links(mg);
	
	return threshold;
}

double		rec_pick_background(Breconstruction* rec, Bimage* p, long number, 
				long avg_kernel, long var_kernel, double excl_dist)
{
	p->filter_average(avg_kernel);

	if ( verbose )
		cout << "Filtered min, max, avg & std:   " << p->minimum() << " " <<
			p->maximum() << " " << p->average() << " " << p->standard_deviation() << endl;

	Bimage*			pvar = p->copy();
	
	pvar->variance(var_kernel);
	
	pvar->rescale(-1/pvar->standard_deviation(), pvar->average()/pvar->standard_deviation());
	
	if ( verbose )
		cout << "Variance min, max, avg & std:  " << pvar->minimum() << " " <<
			pvar->maximum() << " " << pvar->average() << " " << pvar->standard_deviation() << endl;

	part_mask_particles(rec->part, pvar, rec->box_size[0]);
	
	img_mask_above_threshold(pvar, p, p->average() + p->standard_deviation());
	
//	write_img("vt.pif", pvar);
	
	double			threshold(0);

	if ( rec->part ) {
		particle_kill(rec->part);
		rec->part = NULL;
	}

	rec->part = particles_from_peaks(pvar, 1, excl_dist, rec->box_size[0]/2, threshold, 1e30, number);
	
	delete pvar;

	reconstruction_set_part_links(rec);
	
	return threshold;
}

/**
@brief 	Picks background areas not overlapping existing particles.
@param 	*project		project parameter structure.
@param 	number			maximum number of background images to pick.
@param 	avg_kernel		averaging kernel to smooth the image.
@param 	var_kernel		kernel to calculate a local variance image.
@param 	excl_dist		exclusion distance between areas.
@return long			number of background areas.

**/
long		project_pick_background(Bproject* project, long number, 
				long avg_kernel, long var_kernel, double excl_dist)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bimage*				p;
	long				npart(0);
	
	if ( verbose ) {
		cout << "Finding background areas:" << endl;
		cout << "Maximum number:                 " << number << endl;
		cout << "Averaging kernel size:          " << avg_kernel << endl;
		cout << "Variance kernel size:           " << var_kernel << endl;
		cout << "Exclusion distance:             " << excl_dist << endl << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) if ( !mg->select ) {
				p = read_img(mg->fmg, 1, mg->img_num);
				if ( p ) {
					if ( verbose )
						cout << "Picking from micrograph " << mg->id << endl;
					mg->origin = p->size()/2;
					mg->fom = mg_pick_background(mg, p, number, avg_kernel, var_kernel, excl_dist);
					delete p;
				}
			}
		}
		npart = project_count_mg_particles(project);
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			p = read_img(rec->frec, 1, 0);
			if ( p) {
				if ( verbose )
					cout << "Picking from reconstruction " << rec->id << endl;
				rec->origin = p->size()/2;
				rec->fom = rec_pick_background(rec, p, number, avg_kernel, var_kernel, excl_dist);
				delete p;
			}
		}
		npart = project_count_rec_particles(project);
	}
	
	if ( verbose )
		cout << "Background areas picked:        " << npart << endl << endl;

	return npart;
}

Bparticle*	part_pick_sym_axis(Bparticle* part, Bsymmetry& sym, View& refview, Vector3<double> refori, double axis_dist)
{
	long				j, g, s;
	Bparticle*			nupart_list = NULL;
	Bparticle*			nupart = NULL;
	View				vr;
	View*				views, *v;
	Quaternion			qr, qp;
	Vector3<double>		tr;

	qr = refview.quaternion();
		
	for ( j=0, g=1; part; part = part->next, ++g ) {
		tr = part->ori - refori + part->loc; 
		views = symmetry_get_all_views(sym, part->view);
		for ( s=1, v = views; v; v = v->next, ++s ) {
			qp = v->quaternion();
			vr = View(qr*qp);
			nupart = particle_add(&nupart, ++j);
			if ( !nupart_list ) nupart_list = nupart;
			nupart->group = g;
			nupart->sel = s;
			nupart->view = vr;
			nupart->loc = tr + vr.backward().vector3() * axis_dist;
			if ( part->loc[2] < 1 ) nupart->loc[2] = 0;	// for micrographs
			if ( verbose & VERB_FULL )
				cout << j << tab << nupart->view << endl;
		}
		kill_list((char *) views, sizeof(View));
	}

	return nupart_list;
}

/**
@brief 	Picks subregions in 3D particles on the given symmetry axis.
@param 	*project	parameter structure with all parameters.
@param 	sym			point group symmetry.
@param 	sym_axis	symmetry axis to pick subregions.
@param 	axis_dist	distance along symmetry axis.
@return long		number of new particles.

	The existing particles are replaced by the new particles.

**/
long		project_pick_sym_axis(Bproject* project, Bsymmetry& sym, int sym_axis, double axis_dist)
{
	long				n(0);
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	View				refview;
	Vector3<double>		refori;

	if ( sym_axis == 2 ) {
		if ( sym.point() > 200 && sym.point() < 300 ) {	// Dihedral
			refview = View(1, 0, 0, 0);
		} else if ( sym.point() == 432 ) {				// Octahedral
			refview = View(1, -1, 0, 0);
		}
	} else if ( sym_axis == 3 && sym.point() > 300 ) {	// Tetra-, Octa-, Icosahedral
		refview = View(1, 1, 1, 0);
	} else if ( sym_axis == 5 && sym.point() == 532 ) {	// Icosahedral
		refview = View(0, 1.0L/GOLDEN, 1, 0);
	}
	
	if ( verbose ) {
		cout << "Picking subregions based on symmetry " << sym.label() << ":" << endl;
		cout << "Symmetry axis:                   " << refview << " (" << sym_axis << ")" << endl;
		cout << "Distance along axis:             " << axis_dist << endl;
		if ( project->select )
			cout << "Box size:                        " << project->rec->box_size << endl;
		else
			cout << "Box size:                        " << project->field->mg->box_size << endl;
	}

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			refori = rec->box_size/2;
			part = part_pick_sym_axis(rec->part, sym, refview, refori, axis_dist);
			particle_kill(rec->part);
			rec->part = part;
		}
		n = project_count_rec_particles(project);
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				refori = mg->box_size/2;
				part = part_pick_sym_axis(mg->part, sym, refview, refori, axis_dist);
				particle_kill(mg->part);
				mg->part = part;
			}
		}
		n = project_count_mg_particles(project);
	}
	
	if ( verbose )
		cout << "Number of new particles:         " << n << endl << endl;

	return n;
}

Bimage*		img_prepare_projections(Bstring& tempfile, Bsymmetry& sym, double hires, long bin)
{
	if ( tempfile.length() < 1 ) return NULL;
	
	Bimage*			ptemp = read_img(tempfile, 1, 0);
	if ( !ptemp ) {
		cerr << "Error: The template/reference file was not read!" << endl;
		bexit(-1);
	}

	ptemp->check_resolution(hires);
	
	double			angle(hires/ptemp->real_size()[0]);
	View*			views = NULL;
	Bimage*			proj = NULL;
	FSI_Kernel*		kernel = NULL;

	if ( ptemp->sizeZ() > 1 ) {
		if ( verbose ) {
			cout << "Calculating projections:" << endl;
			cout << "Symmetry:                       " << sym.label() << endl;
			cout << "Angle step size:                " << angle*180.0/M_PI << endl;
		}
		kernel = new FSI_Kernel(8, 2);
		views = asymmetric_unit_views(sym, angle, angle, angle, 1);
		proj = ptemp->project(views, hires, kernel);
		delete ptemp;
		ptemp = proj;
		if ( views ) kill_list((char *) views, sizeof(View));
		if ( kernel ) delete kernel;
		if ( verbose ) 
			cout << "Number of projections:          " << ptemp->images() << endl;
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_projections: Calculating the background" << endl;
	ptemp->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_projections: Smoothing the edge" << endl;
	ptemp->edge(1, ptemp->size(), Vector3<double>(0,0,0), 2, FILL_BACKGROUND);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_projections: Doubling the size" << endl;
	double			scale(1.5);
	Vector3<long>	nusize((long)(scale*ptemp->sizeX()), (long)(scale*ptemp->sizeY()), 1);
	Vector3<long>	translate((nusize-ptemp->size())/2);
	ptemp->resize(nusize, translate, FILL_BACKGROUND);
	
	if ( bin > 1 ) ptemp->bin(bin);
	
	if ( verbose & VERB_DEBUG )
		write_img("temp.pif", ptemp, 0);
	
	return ptemp;
}

Bimage*		img_prepare_orientations(Bstring& tempfile, Bsymmetry& sym, double hires, long bin)
{
	if ( tempfile.length() < 1 ) return NULL;
	
	Bimage*			ptemp = read_img(tempfile, 1, 0);
	if ( !ptemp ) {
		cerr << "Error: The template/reference file was not read!" << endl;
		bexit(-1);
	}
	
	if ( ptemp->sizeZ() < 2 ) {
		cerr << "Error: The template must be a 3D volume!" << endl;
		return NULL;
	}

	ptemp->check_resolution(hires);
	
	double			angle(hires/ptemp->real_size()[0]);

	if ( verbose ) {
		cout << "Calculating orientations:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Angle step size:                " << angle*180.0/M_PI << endl;
	}
	
	View*			views = asymmetric_unit_views(sym, angle, angle, angle, 1);

	Bimage*			pmt = ptemp->orient(views);
	delete ptemp;
	
	ptemp = pmt;
	
	if ( views ) kill_list((char *) views, sizeof(View));

	if ( verbose ) 
			cout << "Number of orientations:         " << ptemp->images() << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_orientations: Calculating the background" << endl;
	ptemp->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_orientations: Smoothing the edge" << endl;
	ptemp->edge(1, ptemp->size(), Vector3<double>(0,0,0), 2, FILL_BACKGROUND);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_prepare_orientations: Doubling the size" << endl;
	double			scale(1.5);
	Vector3<long>	nusize((long)(scale*ptemp->sizeX()), (long)(scale*ptemp->sizeY()), (long)(scale*ptemp->sizeZ()));
	Vector3<long>	translate((nusize-ptemp->size())/2);
	ptemp->resize(nusize, translate, FILL_BACKGROUND);
	
	if ( bin > 1 ) ptemp->bin(bin);
	
	if ( verbose & VERB_DEBUG )
		write_img("temp.pif", ptemp, 0);
	
	return ptemp;
}


double		part_align(Bparticle* part, Bimage* p, Bimage* ptemp, 
				double hires, double lores, double shift_limit, long bin, 
				fft_plan planf, fft_plan planb)
{
	long			nn;
	Vector3<double>	shift;
	
	part->fom[0] = 0;

	if ( bin > 1 ) p->bin(bin);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Transforming" << endl;
	p->fft(planf, 0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Band-pass filtering" << endl;
	p->complex_bandpass(hires, lores);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Normalizing" << endl;
	p->complex_normalize();	
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Calculating the complex conjugate product" << endl;
	Bimage*			pcc = p->complex_conjugate_product_one2many(ptemp);
	if ( !pcc ) cerr << "Complex conjugate product failed!" << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Back-transforming" << endl;
	pcc->origin(shift);
	pcc->fft_back(planb, 0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_align: Finding peaks" << endl;
	pcc->find_peak(shift_limit);
	pcc->refine_peak();
	
	for ( nn=0; nn<pcc->images(); nn++ ) {
//		cout << nn << tab << pcc->image[nn].origin() * bin << tab << ptemp->image[nn].view() << tab << pcc->image[nn].FOM() << endl;
		if ( part->fom[0] < pcc->image[nn].FOM() ) {
			part->fom[0] = pcc->image[nn].FOM();
			part->view = ptemp->image[nn].view();
			shift = pcc->image[nn].origin();
		}
	}

	delete pcc;
	
	shift *= bin;

	if ( verbose )
		cout << right << setprecision(2) << part->id << tab << shift 
			<< setprecision(4) << tab << part->view << tab 
			<< setprecision(6) << part->fom[0] << endl;
	
//	cout << "box size = " << part->mg->box_size << endl;
	
	part->loc += shift;
	if ( part->mg ) part->ori = part->mg->box_size/2;
	else if ( part->rec ) part->ori = part->rec->box_size/2;
	
	return part->fom[0];
}


/**
@brief 	Finds the centers of picked particles within a micrograph.
@param 	*mg			micrograph parameter structure.
@param 	*ptemp 		template image.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	bin			level of image binning.
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@return long		number of particles.

	An image processing parameter structure loaded with micrograph
	information is used to extract particle images from the micrograph
	image using the particle coordinates in the parameter structure.
	The extracted particle images are each rotated by PI and the shift
	found by cross-correlation between the unrotated and rotated images.
	The particle coordinates in the parameter structure are updated with
	the shift.

**/
long		mg_extract_orient_particles(Bmicrograph* mg, Bimage* ptemp, 
				double hires, double lores, long bin, fft_plan planf, fft_plan planb)
{
	
	if ( verbose )
		cout << "Micrograph " << mg->id << " (" << mg->img_num << ")" << endl;

	Vector3<long>	extsize(ptemp->size()*bin);
	Vector3<double>	start, origin(extsize/2);
	Bparticle*		part;	
	Bimage*			pex = NULL;
	
	long			npart(particle_count(mg->part));
	if ( npart < 1 ) return 0;
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);

	if ( verbose & VERB_PROCESS )
		cout << "Finding particles from image " << p->file_name() << endl;
		
	double			shift_limit(ptemp->sizeX()/4);
	
	if ( verbose )
		cout << "ID\tdx\tdy\tdz\tvx\tvy\tvz\tva\tFOM" << endl;
	for ( part = mg->part; part; part = part->next ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG mg_extract_orient_particles: Extracting particle" << endl;
		start = part->loc - origin;
		start = start.max(0);
		start = start.min(p->size() - extsize);
		pex = p->extract(0, start, extsize);
		part_align(part, pex, ptemp, hires, lores, shift_limit, bin, planf, planb);
		delete pex;
	}
	
	delete p;
	
	return npart;
}

long		rec_extract_orient_particles(Breconstruction* rec, Bimage* ptemp, 
				double hires, double lores, long bin, fft_plan planf, fft_plan planb)
{
	
	if ( verbose )
		cout << "Reconstruction " << rec->id << endl;

	Vector3<long>	extsize(ptemp->size()*bin);
	Vector3<double>	start, origin(extsize/2);
	Bparticle*		part;	
	Bimage*			pex = NULL;
	
	long			npart(particle_count(rec->part));
	if ( npart < 1 ) return 0;
	
	Bimage*			p = read_img(rec->frec, 1, 0);

	if ( verbose & VERB_PROCESS )
		cout << "Finding particles from image " << p->file_name() << endl;
		
	double			shift_limit(ptemp->sizeX()/4);
	
	if ( verbose )
		cout << "ID\tdx\tdy\tdz\tvx\tvy\tvz\tva\tFOM" << endl;
	for ( part = rec->part; part; part = part->next ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG mg_extract_orient_particles: Extracting particle" << endl;
		start = part->loc - origin;
		start = start.max(0);
		start = start.min(p->size() - extsize);
		pex = p->extract(0, start, extsize);
		part_align(part, pex, ptemp, hires, lores, shift_limit, bin, planf, planb);
		delete pex;
	}
	
	delete p;
	
	return npart;
}


/**
@brief 	Picks particles using cross-correlation.
@param 	*project	project parameter structure.
@param 	&tempfile	template image.
@param 	sym			point group symmetry.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	bin			level of image binning.
@return double		minimum threshold used to accept peaks.

	Each micrograph is cross-correlated with the template image including
	bandpass filtering and frequency space masking.

**/
double		project_extract_orient_particles(Bproject* project, Bstring& tempfile, 
				Bsymmetry& sym, double hires, double lores, long bin)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bimage*				ptemp = NULL;
	
	if ( project->select < 1 ) 
		ptemp = img_prepare_projections(tempfile, sym, hires, bin);
	else 
		ptemp = img_prepare_orientations(tempfile, sym, hires, bin);
		
	if ( verbose ) {
		cout << "Extracting and orienting particles in micrographs:" << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		cout << "Template images:                " << ptemp->images() << endl;
		cout << "Template size:                  " << ptemp->size() << endl;
//		if ( pmask )
//			cout << "Mask:                           " << pmask->file_name() << endl;
		cout << "Resolution limits:              " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
//		cout << "Threshold:                      " << fomcut << endl;
//		cout << "Exclusion distance:             " << excl_dist << endl;
		cout << "Binning:                        " << bin << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Setting up FFT plans" << endl;
	
	fft_plan		planf = ptemp->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = ptemp->fft_setup(FFTW_BACKWARD, 1);
	
	if ( verbose & VERB_PROCESS )
		cout << "Transforming template" << endl;
	
	ptemp->fft(planf, 0);

	if ( verbose & VERB_PROCESS )
		cout << "Band-pass filtering template" << endl;
	
	ptemp->complex_bandpass(hires, lores);
	
	if ( verbose & VERB_PROCESS )
		cout << "Normalizing template" << endl;
	
	ptemp->complex_normalize();	
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				mg_extract_orient_particles(mg, ptemp, hires, lores, bin, planf, planb);
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			rec_extract_orient_particles(rec, ptemp, hires, lores, bin, planf, planb);
		}
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	return 0;
}

