/**
@file	mg_subtomo.cpp
@brief	Functions to find particles (subtomograms) in a tomographic reconstruction.
@author	Juha Huiskonen
@author Bernard Heymann
@date	Created:  20071010
@date	Modified: 20120124 (BH)
@date	Modified: 20120308
@date	Modified: 20120316
@date	Modified: 20120528
@date	Modified: 20121118 (fixed bug in img_divide_with_fom)
@date	Modified: 20150108 (BH) - incorporated into Bsoft
@date	Modified: 20150806 (BH)
	Based on the code from img_find.c
**/

#include "Bimage.h"
#include "mg_subtomo.h"
#include "mg_select.h"
#include "rwimg.h"
#include "mg_processing.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// It creates two temporary images the size of the input image - usually a tomogram
// A rotated template
// The cross correlation map
Bparticle*	img_search_view(Bimage* p, Bimage* ptemp, Bimage* pmask, View view,
				double hires, double lores, double shiftlimit, double mindist,
				double threshold, int maxhits, int refinepeaks,
				fft_plan planf, fft_plan planb)
{
	// View corresponds to a rotation matrix that takes the particle to the standard orientation
	// need to transpose this matrix to rotate the template from the standard orientation on the particle
	Matrix3	mat = view.matrix();
	mat = mat.transpose();

	// Transform the template to the box with size of the reconstruction
	Bimage*		ptemprot = ptemp->rotate(p->size(), mat);

	// Correlate the template and the map
	Bimage*		pcc = ptemprot->cross_correlate(p, hires, lores, pmask, planf, planb);

	delete ptemprot;
			
	Bparticle*	goodpeaks = img_find_refine_peaks(pcc, view, shiftlimit, 0, 0,
					mindist, threshold, maxhits, refinepeaks);

	delete pcc;

	return goodpeaks;
}

// It creates three temporary images the size of the template image
// A cropped input image
// An optional real space mask to apply to the cropped input image - deleted before generating more images
// A rotated template
// The cross correlation map
Bparticle*	img_refine_view(Bimage* pcrop, Bimage* ptemp, Bimage* pmask, Bimage* pmask2,
				View view, double hires, double lores, double shiftlimit,
				double shiftlimitz, double shiftlimitxy, double mindist, int refinepeaks,
				fft_plan planf, fft_plan planb)
{
//	Vector3<double>	scale(1,1,1), trans;

	// View corresponds to a rotation matrix that takes the particle to the standard orientation
	// need to transpose this matrix to rotate the template from the standard orientation on the particle
	Matrix3	mat = view.matrix();
	mat = mat.transpose();

	Bimage*			pcropmasked = pcrop->copy();
	
	// Rotate the real space mask to the current view
	if ( pmask2 ) {
		Bimage*		pmask2rot = pmask2->rotate(pmask2->size(), mat);
		pcropmasked->multiply(pmask2rot);
		delete pmask2rot;
	}
				
	// Rotate the template to the current view
	Bimage*		ptemprot = ptemp->rotate(ptemp->size(), mat);

	// Correlate the template and the map
	Bimage*		pcc = ptemprot->cross_correlate(pcropmasked, hires, lores, pmask, planf, planb);

	delete ptemprot;
	delete pcropmasked;
	
	Bparticle*	peak = NULL;
	long		i(0);
	
	// If no shift allowed, just return the CCC at the origin
	if ( shiftlimit == 0 && shiftlimitz == 0 && shiftlimitxy == 0) {
		peak = particle_add(&peak, 1);
		peak->view = view;
		peak->fom[0] = (*pcc)[i];
	}
	// Shifts allowed so the hit is looked around origin, with the radius of shiftlimit, or within a cylinder
	else {
		peak = img_find_refine_peaks(pcc, view,
					shiftlimit, shiftlimitz, shiftlimitxy,
					mindist, -1, -1, refinepeaks);
	}
	
	delete pcc;

	return peak;
}

/************************************************************************
@brief 	Searches a 3D density map for a template
@author	Juha Huiskonen
@param 	*rec			reconstruction parameters.
@param 	*p				the image.
@param 	*ptemp			the template to be searched for.
@param 	*pmask			reciprocal space mask for cross-correlation (ignored if NULL).
@param 	*pmask2			real space mask for cross-correlation (ignored if NULL).
@param 	alpha_step		angular step size around view vector (radians).
@param 	theta_step		angular step size around view vector (radians).
@param 	phi_step		angular step size around view vector (radians).
@param 	angle_limit		angular limit for refinement in alpha (radians).
@param 	thetaphi_limit	angular limit for refinement in theta & phi (radians).
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	shiftlimit		maximum shift from the original position (binned units).
@param 	mindist 		minimun distance for cc peaks (binned units).
@param 	maxhits			maximum number of hits
@param 	bin				binning for map, template and mask
@param 	*sym			symmetry to generate a list of views for search mode
@param 	refinepeaks		flag to run several iterations in refine
@param 	ccmax_file		file for cross-correlation map (max ccc for each position and rotation of the template)
@return double			the best correlation coefficient.

	The template is rotated and cross-correlated to find fits above the
        threshold.

**/
long	reconstruction_search_subtomo(Breconstruction* rec,
			Bimage* p, Bimage* ptemp, Bimage* pmask, Bimage* pmask2,
			double alpha_step, double theta_step, double phi_step,
			double alpha_limit, double thetaphi_limit,
			double hires, double lores, double shiftlimit, double mindist,
			double threshold, int maxhits, Vector3<long> bin, Bsymmetry& sym,
			int refinepeaks, Bstring ccmax_file)
{
	if ( !pmask ) { cerr << "Warning: No reciprocal space mask given." << endl; }
	if ( !pmask2 ) { cerr << "Warning: No real space mask given." << endl; }
	
	// Search whole rec by default
	if ( shiftlimit < 0 ) {
		if ( verbose & VERB_PROCESS )
			cout << endl << "Shift limit not specified, searching the whole reconstruction." << endl;
		shiftlimit = p->size().length()/2;
	}
		
	// Minimum distance between ccc peaks is the size of the template by default
	if ( mindist <= 0 ) { mindist = 2*ptemp->maximum_included_radius(); }

	// Clear any previous particle records
	if ( rec->part ) particle_kill(rec->part);
	rec->fpart = 0;
	rec->part = NULL;

	// Origin of the template (already binned units)
	Vector3<double>	origin = ptemp->image->origin();

	// Normalize and mask template
	if ( verbose & VERB_PROCESS ) { cout << endl << "Normalizing the template." << endl; }
	ptemp->rescale_to_avg_std(0, 1);

	if ( pmask2 ) {
		if ( verbose & VERB_PROCESS ) { cout << endl << "Masking the template." << endl; }
		ptemp->multiply(pmask2);
	}

	// Normalize map
	if ( verbose & VERB_PROCESS ) { cout << endl << "Normalizing the map." << endl; }
	p->rescale_to_avg_std(0, 1);

	// Search mode
	if ( verbose & VERB_PROCESS ) { cout << endl << "Running in search mode." << endl; }
	cout << endl;

	// Views to consider in Search
	View			refview;
	View*			allviews = NULL; // all views for search/refinement
	if ( sym.point() < 102 ) {
		allviews = views_within_limits(refview, theta_step, phi_step, alpha_step, thetaphi_limit, alpha_limit);
	} else {
		allviews = asymmetric_unit_views(sym, theta_step, phi_step, alpha_step, 1);
	}
	
	long			nviews = count_list((char *)allviews);

//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "Searching for the template in the map:" << endl;
		cout << "Map:                            " << p->file_name() << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		cout << "Map size:                       " << p->size() << endl;
		cout << "Template size:                  " << ptemp->size() << endl;
		cout << "Template origin:                " << setprecision(4) << origin << endl;
		cout << "Binning:                        " << bin << endl;
		cout << "Search radius:                  " << shiftlimit << endl;
		cout << "Minimum separation distance:    " << mindist << endl;
		if ( pmask2)
			cout << "Real space mask:                " << pmask2->file_name() << endl;
		if ( pmask )
			cout << "Reciprocal space mask:          " << pmask->file_name() << endl;
		cout << "Alpha step size:                " << alpha_step*180.0/M_PI << " degrees" << endl;
		cout << "Theta step size:                " << theta_step*180.0/M_PI << " degrees" << endl;
		cout << "Phi step size:                  " << phi_step*180.0/M_PI << " degrees" << endl;
		cout << "Limit in alpha:                 " << alpha_limit*180.0/M_PI << " degrees" << endl;
		if ( sym.point() > 101 )
			cout << "Symmetry:                       " << sym.label() << endl;
		else
			cout << "Limit in theta & phi:           " << thetaphi_limit*180.0/M_PI << " degrees" << endl;
		cout << "Number of views:                " << nviews << endl;
		if ( ccmax_file.length() )
			cout << "Cross correlation map:          " << ccmax_file << endl;
		cout << endl;
	}

	if ( pmask )
		pmask->mask_fspace_resize(p->size());

	long 			i, ipart(0);
	View*			v = NULL;
	View*			view_arr = new View[nviews];
	for ( i=0, v=allviews; v; v=v->next, i++ ) view_arr[i] = *v;

	Bparticle* 		part = NULL;
	Bparticle**		part_arr = new Bparticle*[nviews];

	fft_plan		planf = p->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = p->fft_setup(FFTW_BACKWARD, 1);
	
	if ( verbose & VERB_RESULT )
		cout << "#\tViewX\tViewY\tViewZ\tViewA\tNpart" << endl;
	
	// Loop through all views for the current particle and iteration
#ifdef HAVE_GCD
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t i){
		part_arr[i] = img_search_view(p, ptemp, pmask, view_arr[i], hires, lores,
				shiftlimit, mindist, threshold, maxhits, refinepeaks, planf, planb);
		long	np = count_list((char *) part_arr[i]);
		dispatch_sync(myq, ^{
			if ( verbose & VERB_RESULT )
				cout << i+1 << tab << view_arr[i] << tab << np << endl << flush;
		});
	});
#else
#pragma omp parallel for
	for ( long i=0; i < nviews; i++ ) {
		part_arr[i] = img_search_view(p, ptemp, pmask, view_arr[i], hires, lores,
				shiftlimit, mindist, threshold, maxhits, refinepeaks, planf, planb);
		long	np = count_list((char *) part_arr[i]);
	#pragma omp critical
		{
			if ( verbose & VERB_RESULT )
				cout << i+1 << tab << view_arr[i] << tab << np << endl << flush;
		}
	}
#endif

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	if ( verbose )
		cout << "Transferring the particles found to the reconstruction" << endl;
	for ( i=0; i < nviews; i++ ) if ( part_arr[i] ) {
		if ( !rec->part ) rec->part = part = part_arr[i];
		else part->next = part_arr[i];
		for ( ; part->next; part = part->next ) ;
	}
	
	Vector3<double>	translate;

	for ( ipart=0, part = rec->part; part; part = part->next ) {
		part->id = ++ipart;
		part->ori = origin * bin;
		translate = origin - part->loc;
		translate = vector3_set_PBC(translate, p->size());
		part->loc = translate * bin;
		if ( verbose & VERB_FULL )
			cout << ipart << tab << setprecision(2)
				<< part->loc << tab << setprecision(4) << part->view << tab << part->fom[0] << endl;
	}
	
	kill_list((char *) allviews, sizeof(View));
	delete[] view_arr;
	delete[] part_arr;

	if ( rec->box_size.volume() < 1 )
		rec->box_size = ptemp->size() * bin;
	
	return particle_count(rec->part);
}

/**
@brief 	Refines the view vectors for particles already in the project
@author	Juha Huiskonen
@param 	*rec			reconstruction parameters.
@param 	*p				the image.
@param 	*ptemp			the template to be searched for.
@param 	*pmask			reciprocal space mask for cross-correlation (ignored if NULL).
@param 	*pmask2			real space mask for cross-correlation (ignored if NULL).
@param 	alpha_step_orig		angular step size around view vector (radians).
@param 	theta_step_orig		angular step size around view vector (radians).
@param 	phi_step_orig		angular step size around view vector (radians).
@param 	alpha_limit_orig	angular limit for refinement in alpha (radians).
@param 	thetaphi_limit_orig	angular limit for refinement in theta & phi (radians).
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	shiftlimit_orig		maximum shift from the original position (binned units).
@param 	shiftlimitz_orig	maximum z-shift from the original position (binned units).
@param 	shiftlimitxy_orig	maximum xy-shift from the original position (binned units).
@param 	mindist 		minimun distance for cc peaks (binned units).
@param 	bin				binning for map, template and mask
@param 	*sym			symmetry to generate a list of views for search mode
@param 	iters			number of iterations in refine
@param 	refinepeaks		flag to run several iterations in refine
@param 	ccmax_file		file for cross-correlation map (max ccc for each position and rotation of the template)
@return double			the best correlation coefficient.

	The template is rotated and cross-correlated to find fits above the
        threshold.

**/
long	reconstruction_refine_subtomo(Breconstruction* rec,
			Bimage* p, Bimage* ptemp, Bimage* pmask, Bimage* pmask2,
			double alpha_step_orig, double theta_step_orig, double phi_step_orig,
			double alpha_limit_orig, double thetaphi_limit_orig,
			double hires, double lores, double shiftlimit_orig,
			double shiftlimitz_orig, double shiftlimitxy_orig, double mindist,
			Vector3<long> bin, Bsymmetry& sym, int iters, int refinepeaks, Bstring ccmax_file)
{
	if ( !pmask ) { cerr << "Warning: No reciprocal space mask given." << endl; }
	if ( !pmask2 ) { cerr << "Warning: No real space mask given." << endl; }

	// Refine shifts within 1/2 template by default
	if ( shiftlimit_orig < 0 ) {
		if ( verbose & VERB_PROCESS )
			cout << endl << "Shift limit not specified, using half the size of the template." << endl;
		shiftlimit_orig = ptemp->maximum_included_radius();
	}
	
	if ( mindist <= 0 ) { mindist = 2*ptemp->maximum_included_radius(); }

	// Normalize and mask template
	if ( verbose & VERB_PROCESS ) { cout << endl << "Normalizing the template." << endl; }
	ptemp->rescale_to_avg_std(0, 1);

	if ( pmask2 ) {
		if ( verbose & VERB_PROCESS ) { cout << endl << "Masking the template." << endl; }
		ptemp->multiply(pmask2);
	}

	// Normalize map
	if ( verbose & VERB_PROCESS ) { cout << endl << "Normalizing the map." << endl; }
	p->rescale_to_avg_std(0, 1);

	// Refine mode
	if ( verbose & VERB_PROCESS ) {	cout << "Running in refine mode" << endl; }

	if ( sym.point() > 101 && iters > 1) {
		iters = 1;
		cerr << "Warning: Symmetry selected, running only one iteration." << endl;
	}

	if ( verbose & VERB_RESULT ){  cout << endl << "Padding the template to a larger box." << endl;}

	double			maxshift = shiftlimit_orig + shiftlimitz_orig + shiftlimitxy_orig;

	// Origin of the template (already binned units)
	Vector3<double>	origin = ptemp->image->origin();

	Vector3<double>	scale(1,1,1);
	Vector3<double>	translate = Vector3<double>(maxshift,maxshift,maxshift);
	View			refview;
	Matrix3			mat = refview.matrix();
	Vector3<long>	ptempsize = ptemp->size() + Vector3<long>(2*(long)maxshift, 2*(long)maxshift, 2*(long)maxshift);
	ptemp->resize(ptempsize, translate, FILL_BACKGROUND, 0);
	ptemp->origin(origin+translate);
	Vector3<double>	neworigin = ptemp->image->origin();

	if ( pmask )
		pmask->mask_fspace_resize(p->size());

	if ( pmask2 ) {
		if ( verbose & VERB_RESULT ){  cout << endl << "Padding the real space mask to a larger box." << endl;}
		pmask2->resize(ptempsize, translate, FILL_BACKGROUND, 0);
	}

	if ( verbose ) {
		cout << endl << "Refining the template in the map around the starting view:" << endl;
		cout << "Map:                            " << p->file_name() << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		cout << "Map size:                       " << p->size() << endl;
		cout << "Template size:                  " << ptemp->size() << endl;
		cout << "Template origin:                " << setprecision(4) << ptemp->image->origin() << endl;
		cout << "Binning:                        " << bin << endl;
		cout << "Maximum translation:            " << shiftlimit_orig << " pixels" << endl;
		if ( shiftlimitz_orig > 0 )
		cout << "Additional shift along view:    " << shiftlimitz_orig << endl;
		if ( shiftlimitxy_orig > 0 )
		cout << "Additional shift orthogonal view: " << shiftlimitxy_orig << endl;
		if ( pmask2)
		cout << "Real space mask:                " << pmask2->file_name() << endl;
		if ( pmask )
		cout << "Reciprocal space mask:          " << pmask->file_name() << endl;
		cout << "Alpha step size:                " << alpha_step_orig*180.0/M_PI << " degrees" << endl;
		cout << "Theta step size:                " << theta_step_orig*180.0/M_PI << " degrees" << endl;
		cout << "Phi step size:                  " << phi_step_orig*180.0/M_PI << " degrees" << endl;
		cout << "Limit in alpha:                 " << alpha_limit_orig*180.0/M_PI << " degrees" << endl;
		if ( sym.point() < 102 )
		cout << "Limit in theta & phi:           " << thetaphi_limit_orig*180.0/M_PI << " degrees" << endl;;
		cout << endl;
	}

	long 			iter, nviews, i, ibest(0);
	long			npart = particle_count(rec->part);
	double			shiftlimit, shiftlimitz, shiftlimitxy;
	double			alpha_step,theta_step,phi_step;
	double			alpha_limit, thetaphi_limit;
	double			best_cc;
	Vector3<double>	bestshift;
	View			bestview;
	View*			v = NULL;
	View*			allviews = NULL;
	Bparticle* 		part = NULL;
	Bimage*			pcrop = NULL;

	fft_plan		planf = fft_setup_plan(ptempsize, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(ptempsize, FFTW_BACKWARD, 1);
	
	for ( part = rec->part; part; part = part->next ) if ( part->sel > 0 ) {

		if ( verbose )
			cout << "Particle: " << part->id << " of " << npart << endl;

		if ( verbose & VERB_PROCESS ) {
			cout << "Starting location:              " << part->loc << endl;
			cout << "Starting view:                  " << part->view << endl;
		}


		// Reset angle step sizes to originals for the current particle (halved after each iteration)
		alpha_step = alpha_step_orig;
		theta_step = theta_step_orig;
		phi_step = phi_step_orig;
		thetaphi_limit = thetaphi_limit_orig;
		alpha_limit = alpha_limit_orig;
		shiftlimit = shiftlimit_orig;
		shiftlimitxy = shiftlimitxy_orig;
		shiftlimitz = shiftlimitz_orig;

		// The initially best view is the given view
		bestview = part->view;

		for ( iter=1; iter <= iters; iter++) {

			// Views around the current best view to consider in refinement
			if ( sym.point() < 102 ) {
				allviews = views_within_limits(bestview, theta_step, phi_step, alpha_step, thetaphi_limit, alpha_limit);
			} else {
				allviews = asymmetric_unit_views(sym, theta_step, phi_step, alpha_step, 1);
			}

			nviews = count_list((char *)allviews);

			if ( verbose & VERB_PROCESS ) {
				cout << "Iteration:                      " << iter << " of " << iters << endl;
				cout << "Alpha step size:                " << alpha_step*180.0/M_PI << " degrees" << endl;
				cout << "Theta step size:                " << theta_step*180.0/M_PI << " degrees" << endl;
				cout << "Phi step size:                  " << phi_step*180.0/M_PI << " degrees" << endl;
				cout << "Limit in alpha:                 " << alpha_limit*180.0/M_PI << " degrees" << endl;
				if ( sym.point() < 102 )
				cout << "Limit in theta & phi:           " << thetaphi_limit*180.0/M_PI << " degrees" << endl;;
				cout << "Number of views:                " << nviews << endl;
				cout << endl;
			}

			// Crop reconstruction to the same size of a box as the template using the refined location
			if ( verbose & VERB_FULL ){ cout << endl << "Cropping the map." << endl << endl;}
			translate = -part->loc;
			if ( bin[0] ) { translate /= bin; }
			translate = translate + neworigin;
			mat = refview.matrix();
			pcrop = p->transform(ptempsize, scale, Vector3<double>(0,0,0), translate, mat, FILL_BACKGROUND, 0);
	
			// Normalize cropped map
			pcrop->rescale_to_avg_std(0, 1);
			pcrop->origin(Vector3<double>(0,0,0));

			// Since we are contracting around the previous best point, the best orientation is always found from the current views
			best_cc = -1;

			View*				view_arr = new View[nviews];
			for ( i=0, v=allviews; v; v=v->next, i++ ) view_arr[i] = *v;

			Bparticle**			part_arr = new Bparticle*[nviews];

#ifdef HAVE_GCD
			dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t i){
				part_arr[i] = img_refine_view(pcrop, ptemp, pmask, pmask2, view_arr[i],
								hires, lores, shiftlimit, shiftlimitz, shiftlimitxy,
								mindist, refinepeaks, planf, planb);
			});
#else
#pragma omp parallel for
			for ( i=0; i < nviews; i++ ) {
				part_arr[i] = img_refine_view(pcrop, ptemp, pmask, pmask2, view_arr[i],
								hires, lores, shiftlimit, shiftlimitz, shiftlimitxy,
								mindist, refinepeaks, planf, planb);
			}
#endif

			for ( i=0; i < nviews; i++ ) {
				if ( part_arr[i]->fom[0] > best_cc ) {
					best_cc = part_arr[i]->fom[0];
					ibest = i;
					bestview = part_arr[i]->view;
					bestshift = part_arr[i]->loc;
				}
				particle_kill(part_arr[i]);
			}
			
			delete[] view_arr;
			delete[] part_arr;
			
			translate = -bestshift;

			// wrap translation vectors (still working with binned units)
			translate = vector3_origin_to_shift(translate, pcrop->size());

			// Half the angle steps angle limits and shift limits for the next iteration
			alpha_step /= 2;
			theta_step /= 2;
			phi_step /= 2;
			alpha_limit /= 2;
			thetaphi_limit /= 2;
			shiftlimit /= 2;
			shiftlimitz /= 2;
			shiftlimitxy /= 2;
			
			kill_list((char *) allviews, sizeof(View));

			// Take binning into account in the translation
			if ( bin[0] ) { translate *= bin; }
		
			// Update refined orientation and origin + fom
			part->view = bestview;
			part->loc = part->loc + translate;
			part->fom[0] = best_cc;

			delete pcrop;

			if ( verbose )
				cout << "        \r" << flush << iter << tab << ibest+1 << tab << part->loc << tab << bestview << tab << best_cc << endl;
			
		} // for ( iter=1; iter <= iters; iter++)

//		if ( verbose )
//			cout << endl << "Final:\t" << ibest << tab << part->loc << tab <<
//						bestview << tab << best_cc << endl;

	} // if ( part->sel > 0 )

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	if ( rec->box_size.volume() < 1 )
		rec->box_size = ptemp->size() * bin;
	
	return particle_count(rec->part);
}

/**
@brief 	Finds peaks in an image to the nearest voxel iteratively. 
@author	Juha Huiskonen
@param 	*pcc			cross correlation map (not altered).
@param 	view			view of the particle to be refined
@param 	shift_limit		radius of spherical or cylindrical search space (if < 0, default 1e30).
@param 	shift_along		additional shift allowed in the direction of the view vector.
@param 	shift_orthogonal	additional shift allowed orthogonal to the view vector.
@param 	mindist			2 * template radius: used for a spherical mask
@param 	threshold		threshold. if value is <0, only the global maximum is returned
@param 	maxhits
@param 	refinepeaks
@return Bparticle*		list of peaks as particles.

	After a maximum value is found, it is masked with a spherical mask and the next largest value
	is found, until all the values are below the threshold .

**/
Bparticle*	img_find_refine_peaks(Bimage* pcc, View view, double shift_limit,
				double shift_along, double shift_orthogonal, double mindist,
				double threshold, int maxhits, int refinepeaks)
{
	double			peakcc;
	Bparticle*		goodpeaks = NULL;
	Bparticle*		peak = NULL;
	int				hits(0);

	do {
		// shift vectors are in actual pixel coordinates (no wrapping)
		peakcc = img_find_peak_subtomo(pcc, view, shift_limit, shift_along, shift_orthogonal);

		if ( refinepeaks ) { pcc->refine_peak(); }

		if ( (peakcc >= threshold) || (threshold < 0) ) {

			hits++;
			
			peak = particle_add(&goodpeaks, hits);
			peak->loc = pcc->image->origin();
			peak->view = view;
			peak->fom[0] = peakcc;
			
			if ( verbose && maxhits > 1 )
				cout << hits << tab << peak->loc << tab << view << tab << peakcc << endl;

			if ( verbose & VERB_PROCESS )
				cout << "Found a cross-correlation peak:\t" << peak->loc << tab << peakcc << endl;

			// return only the global maximum
			if (threshold < 0) { break; }

			// Enough particles found
			if (maxhits > 0 && hits >= maxhits) { break; }

		}
		
		pcc->shell_wrap(0, pcc->image->origin(), 0, mindist, 0, FILL_USER, 0);

	} while ( peakcc >= threshold );

	//write_img("temp_pcc.spi", pcc);

	if ( hits == 0 ) {
		if ( verbose & VERB_PROCESS )
			cout << "No cross-correlation peaks found above the threshold (the maximum was " << peakcc << ")" << endl;
	}
	
	return goodpeaks;
}

/**
@brief 	Finds the peak in an image to the nearest voxel.
@author	Juha Huiskonen
@param 	*p				image (not altered).
@param 	view			view of the particle to be refined
@param 	shift			radius of spherical or cylindrical search space (if < 0, default 1e30).
@param 	shift_along		additional shift allowed in the direction of the view vector.
@param 	shift_orthogonal	additional shift allowed orthogonal to the view vector.
@return double	 		peak maximum.

	An image is searched for the global maximum (typically used to find the shift vector in a cross-correlation map).
	The peak vector is returned in the image origin in actual pixel coordinates (no wrapping).
	The maximum is returned in the image FOM.

**/
double		img_find_peak_subtomo(Bimage* p, View view, double shift,
				double shift_along, double shift_orthogonal)
{
	if ( !p->data_pointer() ) return 0;
	
	if ( shift < 0 ) shift = 1e30;
	
	long		i(0), x, y, z;
	long				xh = p->sizeX()/2, yh = p->sizeY()/2, zh = p->sizeZ()/2;
	double				x2, y2, z2, r2;
	double				value(0), max;
	double				shift2 = shift*shift;
	Vector3<double>		peak;
	Vector3<double>		origin(0,0,0);
	Vector3<double>		point;
	Vector3<double>		endpoint1;
	Vector3<double>		endpoint2;
	Vector3<double>		closestpoint;
	
	view = view.backward();

	if ( shift_along > 0 ) {
		endpoint1 = Vector3<double>(view[0] * shift_along, view[1] * shift_along, view[2] * shift_along);
		endpoint2 = Vector3<double>(-view[0] * shift_along, -view[1] * shift_along, -view[2] * shift_along);
	}
	
		max = -1e37;
		for ( z=0; z<p->sizeZ(); z++ ) {
			z2 = origin[2] - z;
			if ( z2 < -zh ) z2 += p->sizeZ();
			if ( z2 >  zh ) z2 -= p->sizeZ();

			for ( y=0; y<p->sizeY(); y++ ) {
				y2 = origin[1] - y;
				if ( y2 < -yh ) y2 += p->sizeY();
				if ( y2 >  yh ) y2 -= p->sizeY();

				for ( x=0; x<p->sizeX(); x++, i++ ) {
					x2 = origin[0] - x;
					if ( x2 < -xh ) x2 += p->sizeX();
					if ( x2 >  xh ) x2 -= p->sizeX();

					value = (*p)[i];

					if ( shift_along > 0 ) {	
						point = Vector3<double>((double)x2, (double)y2, (double)z2);
						// squared distance to the closest point on the line segment
						r2 = closest_point_line_distance2(point, endpoint1, endpoint2);
					}
					else if ( shift_orthogonal > 0 ) {
						point = Vector3<double>((double)x2, (double)y2, (double)z2);
						// squared distance to the closest point on the disc
						r2 = closest_point_disc_distance2(point, Vector3<double>(0,0,0), view, shift_orthogonal);
					}
					else {
						// squared distance to the origin
						r2 = x2*x2 + y2*y2 + z2*z2;
					}
					if ( r2 > shift2 ) value = 0;
					if ( max < value ) {
						if ( r2 > shift2 ) value = 0;
						max = value;
						peak = Vector3<double>(x,y,z);
					}
				}
			}
		}
		p->image->origin(peak);
		p->image->FOM(max);
		if ( verbose & VERB_FULL )
			cout << "Peak: " << peak << tab << max << endl;
	
	return max;
}

/**
@brief 	Calculates the closest point on a line segment from a given point
@author	Juha Huiskonen
@param 	p	point
@param 	v	endpoint 1
@param 	w	endpoint 2
@return Vector3<double>		closest point

	Checks whether one of the end points is the closest point

**/

Vector3<double> closest_point_line( Vector3<double> p, Vector3<double> v, Vector3<double> w )
{
	double l2, t;
	Vector3<double> closest_point;
	Vector3<double> line;

	line  = w-v;
	l2 = line.length2();
 
	// parametric value of vector to the closest point is t = scalar(p - v, w - v) / l*l 
	t = (p-v).scalar(line) / l2;

	// the closest point is end point v
	if ( t < 0 ) { closest_point = v; }

	// the closest point is end point w
	else if ( t > 1 ) { closest_point = w; }

	// the closest point is the projection of p that is v + t*line
	else { closest_point = Vector3<double>( v[0] + t * line[0],
                                               v[1] + t * line[1],
                                               v[2] + t * line[2]);
        }
	return closest_point;
}

/**
@brief 	Calculates the squared distance to the closest point on a line segment from a given point
@author	Juha Huiskonen
@param 	p 	point
@param 	v	endpoint 1
@param 	w	endpoint 2
@return double				squared distance to the closest point

	Checks wether the closet point falls within the line segment. If not, one of the endpoints
	is the closest point.

**/

double closest_point_line_distance2( Vector3<double> p, Vector3<double> v, Vector3<double> w )
{
	double l2, t;
	Vector3<double> closest_point;
	Vector3<double> line;

	line  = w-v;
	l2 = line.length2();
 
	// parametric value of vector to the closest point is t = scalar(p - v, w - v) / l*l 
	t = (p-v).scalar(line) / l2;

	// the closest point is end point v
	if ( t < 0 ) { closest_point = v; }

	// the closest point is end point w
	else if ( t > 1 ) { closest_point = w; }

	// the closest point is the projection of p that is v + t*line
	else { closest_point = Vector3<double>( v[0] + t * line[0],
                                               v[1] + t * line[1],
                                               v[2] + t * line[2]);
        }
	return (closest_point-p).length2();
}


/**
@brief 	Calculates the squared distance to the closest point on a disc from a given point
@author	Juha Huiskonen
@param 	p					point
@param 	q					centre point of the disc
@param 	view
@param 	radius
@return double				squared distance to closest point

	Checks whether projection of the point (p') falls on the disc. If yes, this is the closest point
	Otherwise point on the circumference of the disc is the closest point.

**/

double closest_point_disc_distance2( Vector3<double> p, Vector3<double> q, View view, double radius )
{
	Vector3<double> pprime; // projection of p on to the plane

	double p_to_pprime2;  // squared distance from p' to the original point p

	double pprime_to_q2; // squared distance from p' to the center of the disc q
	double pprime_to_q; // distance from p' to the center of the disc q

	double pprime_to_c2; // squared distance from p' to the circumference of the disc c
	double p_to_c2; // squared distance from p to the circumference of the disc c

	double radius2 = radius * radius;
	double a, c;

	c = q[0]*view[0] + q[1]*view[1] + q[2]*view[2];

	a = (view[0] * p[0] + view[1] * p[1] + view[2] * p[2]) - c;
	pprime = Vector3<double>(p[0] - a * view[0], p[1] - a * view[1], p[2] - a * view[2]);

	pprime_to_q2 = (pprime-q).length2();
	p_to_pprime2 = (p-pprime).length2();

	if (pprime_to_q2 > radius2) {

		pprime_to_q = (pprime-q).length();
		pprime_to_c2 = (pprime_to_q - radius) * (pprime_to_q - radius);

		p_to_c2 = pprime_to_c2 + p_to_pprime2;

		return p_to_c2;
	}
	else {
		return p_to_pprime2;
	}
}

/**
@brief 	Least squares fit a sphere to 3D data (particle locations)
@author	Juha Huiskonen
@param 	*part		particle
@param 	N				iterations
@param 	Nstop		stopping condition: tolerance in change of sphere center
@return Sphere				fitted sphere struct

	Algorithm by ImaginaryZ
	From http://imaginaryz.blogspot.co.uk/2011/04/least-squares-fit-sphere-to-3d-data.html

	All you have to do is define:

	Error = Sum( |Position[n] - Center|^2 - Radius^2 )

	Then define the squared error:

	Squared Error = Sum( ( |Position[n] - Center|^2 - Radius^2 )^2 )

	And solve the summation using a iterative method (like newtons, below) after pulling out the summation terms.
	For example, if you do: Sum( (P.x[n] - Cx)^2 ) You get (after Expand):
	Sum( P.x[n]^2 - 2*P.x[n]*Cx + Cx^2 )
	And you can then split up the sum:
	Sum( P.x[n]^2 ) + Sum( P.x[n] ) * -2*Cx + Cx * Nelements
	Note you HAVE to ultimately divide the sums by Nelements
	
	Note that "Center" is A,B,C (3D) and I use Rsq as Radius^2.
	
	This method is not fast, but it converges, and the way the code is written it is independent of dataset size,
	but you do have to compute a number of sums and products before running the algorithm.
	
	Note this method is used to generate the equations used to compute linear and quadratic fits instantly, given you compute some sums first.

**/

Sphere		locations_fit_sphere(Bparticle* part, int N, double Nstop)
{
	Sphere sphere;
	long Pnpoints(0); 
	Bparticle* p = NULL;

	double PXsum=0, PXsumsq=0, PXsumcube=0;
	double PYsum=0, PYsumsq=0, PYsumcube=0;
	double PZsum=0, PZsumsq=0, PZsumcube=0;

	double PXYsum=0, PXZsum=0, PYZsum=0;
	double PX2Ysum=0, PX2Zsum=0, PY2Xsum=0, PY2Zsum=0, PZ2Xsum=0, PZ2Ysum=0;

	Pnpoints = count_list((char *)part);
	double x, y, z;

	for ( p = part; p; p = p-> next ) {

		x = p->loc[0];
		y = p->loc[1];
		z = p->loc[2];

		PXsum += x;
		PXsumsq += x * x;
		PXsumcube += x * x * x;

		PYsum += y;
		PYsumsq += y * y;
		PYsumcube += y * y * y;

		PZsum += z;
		PZsumsq += z * z;
		PZsumcube += z * z * z;

		PXYsum += x * y;
		PXZsum += x * z;
		PYZsum += y * z;

		PX2Ysum += x * x * y;
		PX2Zsum += x * x * z;
		PY2Xsum += y * y * x;
		PY2Zsum += y * y * z;
		PZ2Xsum += z * z * x;
		PZ2Ysum += z * z * y;
	}

	double Xn = PXsum/Pnpoints;        //sum( X[n] )
	double Xn2 = PXsumsq/Pnpoints;    //sum( X[n]^2 )
	double Xn3 = PXsumcube/Pnpoints;    //sum( X[n]^3 )
	double Yn = PYsum/Pnpoints;        //sum( Y[n] )
	double Yn2 = PYsumsq/Pnpoints;    //sum( Y[n]^2 )
	double Yn3 = PYsumcube/Pnpoints;    //sum( Y[n]^3 )
	double Zn = PZsum/Pnpoints;        //sum( Z[n] )
	double Zn2 = PZsumsq/Pnpoints;    //sum( Z[n]^2 )
	double Zn3 = PZsumcube/Pnpoints;    //sum( Z[n]^3 )

	double XY = PXYsum/Pnpoints;        //sum( X[n] * Y[n] )
	double XZ = PXZsum/Pnpoints;        //sum( X[n] * Z[n] )
	double YZ = PYZsum/Pnpoints;        //sum( Y[n] * Z[n] )
	double X2Y = PX2Ysum/Pnpoints;    //sum( X[n]^2 * Y[n] )
	double X2Z = PX2Zsum/Pnpoints;    //sum( X[n]^2 * Z[n] )
	double Y2X = PY2Xsum/Pnpoints;    //sum( Y[n]^2 * X[n] )
	double Y2Z = PY2Zsum/Pnpoints;    //sum( Y[n]^2 * Z[n] )
	double Z2X = PZ2Xsum/Pnpoints;    //sum( Z[n]^2 * X[n] )
	double Z2Y = PZ2Ysum/Pnpoints;    //sum( Z[n]^2 * Y[n] )

	//Reduction of multiplications
	double F0 = Xn2 + Yn2 + Zn2;
	double F1 = 0.5*F0;
	double F2 = -8.0*(Xn3 + Y2X + Z2X);
	double F3 = -8.0*(X2Y + Yn3 + Z2Y);
	double F4 = -8.0*(X2Z + Y2Z + Zn3);

	//Set initial conditions:
	double A = Xn;
	double B = Yn;
	double C = Zn;

	//First iteration computation:
	double A2 = A*A;
	double B2 = B*B;
	double C2 = C*C;
	double QS = A2 + B2 + C2;
	double QB = - 2*(A*Xn + B*Yn + C*Zn);

	//Set initial conditions:
	double Rsq = F0 + QB + QS;

	//First iteration computation:
	double Q0 = 0.5*(QS - Rsq);
	double Q1 = F1 + Q0;
	double Q2 = 8*( QS - Rsq + QB + F0 );
	double aA,aB,aC,nA,nB,nC,dA,dB,dC;

	//Iterate N times, ignore stop condition.
	int n(0);
	while( n != N ){
	    n++;
	
	    //Compute denominator:
	    aA = Q2 + 16*(A2 - 2*A*Xn + Xn2);
	    aB = Q2 + 16*(B2 - 2*B*Yn + Yn2);
	    aC = Q2 + 16*(C2 - 2*C*Zn + Zn2);
	    aA = (aA == 0) ? 1.0 : aA;
	    aB = (aB == 0) ? 1.0 : aB;
	    aC = (aC == 0) ? 1.0 : aC;


	    //Compute next iteration
	    nA = A - ((F2 + 16*( B*XY + C*XZ + Xn*(-A2 - Q0) + A*(Xn2 + Q1 - C*Zn - B*Yn) ) )/aA);
	    nB = B - ((F3 + 16*( A*XY + C*YZ + Yn*(-B2 - Q0) + B*(Yn2 + Q1 - A*Xn - C*Zn) ) )/aB);
	    nC = C - ((F4 + 16*( A*XZ + B*YZ + Zn*(-C2 - Q0) + C*(Zn2 + Q1 - A*Xn - B*Yn) ) )/aC);
	
	    //Check for stop condition
	    dA = (nA - A);
	    dB = (nB - B);
	    dC = (nC - C);
	    if( (dA*dA + dB*dB + dC*dC) <= Nstop ){ break; }

	    //Compute next iteration's values
	    A = nA;
	    B = nB;
	    C = nC;
	    A2 = A*A;
	    B2 = B*B;
	    C2 = C*C;
	    QS = A2 + B2 + C2;
	    QB = - 2*(A*Xn + B*Yn + C*Zn);
	    Rsq = F0 + QB + QS;
	    Q0 = 0.5*(QS - Rsq);
	    Q1 = F1 + Q0;
	    Q2 = 8*( QS - Rsq + QB + F0 );
	}

	double R = sqrt(Rsq);
	sphere.A = A;
	sphere.B = B;
	sphere.C = C;
	sphere.R = R;

	return (sphere);
}


