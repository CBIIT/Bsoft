/**
@file	Bimage_search.cpp
@brief	Searches orientation space for the best fit of a 3D map to a template.
@author Bernard Heymann
@date	Created: 20021027
@date	Modified: 20160608
**/

#include "Bimage.h"
//#include "rwimg.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Searches a 2D/3D density map for a template.
@param 	*ptemp			the template to be searched for.
@param 	*views			list of views to search.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	search_radius	radius for shift search.
@param 	*pmask			mask for cross-correlation (ignored if NULL).
@param 	&currview		best current view to return.
@param 	&currshift		best current shift to return.
@return double			the best correlation coefficient.

	The template is rotated and cross-correlated to find the best fit.
	The views must be calculated externally to allow for custom sets.

**/
double		Bimage::search_views(Bimage* ptemp, View* views,
				double hires, double lores, double search_radius,
				Bimage* pmask, View& currview, Vector3<double>& currshift)
{
	long 				i, nviews, ibest(0);
	double				best(-1e37);
	View				bestview;
	View*				v;

	for ( nviews=0, v=views; v; v=v->next ) nviews++;
		
	View*				view_arr = new View[nviews];
	double*				cc = new double[nviews];
	Vector3<double>*	translate = new Vector3<double>[nviews];
	
	for ( i=0, v=views; v; v=v->next, i++ ) view_arr[i] = *v;
			
	if ( search_radius < 0 ) search_radius = ptemp->sizeX()/4;
	
	change_type(Float);
	ptemp->change_type(Float);
	rescale_to_avg_std(0, 1);
	ptemp->rescale_to_avg_std(0, 1);
	calculate_background();
	ptemp->calculate_background();

	Vector3<double>	ori(ptemp->image->origin());
	
	if ( ori.length() == 0 ) {
		ori = ptemp->size()/2;
		ptemp->origin(ori);
	}
	
	if ( verbose ) {
		cout << "Finding a template in a map:" << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		cout << "Map:                            " << file_name() << endl;
		cout << "Template origin:                " << ptemp->image->origin() << endl;
		cout << "Map origin:                     " << image->origin() << endl;
		cout << "Number of views:                " << nviews << endl;
		cout << "Resolution range:               ";
		if ( lores <= 0 ) cout << "inf";
		else cout << lores;
		cout << " - " << hires << " A" << endl;
		cout << "Shift search radius:            " << search_radius << " pixels" << endl;
		if ( pmask )
			cout << "Frequency space mask:           " << pmask->file_name() << endl;
	}
	
	fft_plan		planf = fft_setup_plan(ptemp->size(), FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(ptemp->size(), FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t k){
		translate[k] = rotate_cross_correlate(ptemp, view_arr[k], hires, lores,
			search_radius, pmask, cc[k], planf, planb);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nviews; i++ ) {
		translate[i] = rotate_cross_correlate(ptemp, view_arr[i], hires, lores,
			search_radius, pmask, cc[i], planf, planb);
	}
#endif

	if ( verbose & VERB_RESULT )
		cout << "View\ttx\tty\ttz\tvx\tvy\tvz\tva\tCC" << endl;
	for ( i=0; i<nviews; i++ ) {
		if ( best < cc[i] ) {
			best = cc[i];
			ibest = i;
		}
		if ( verbose & VERB_RESULT )
			cout << i << tab << setprecision(2) << translate[i] << tab
					<< setprecision(4) << view_arr[i] << tab << cc[i] << endl;
	}

	bestview = view_arr[ibest];

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	image->origin(ptemp->image->origin() + translate[ibest]);
	view(bestview);
	image->FOM(best);
	
	if ( verbose )
		cout << "Best view:\t" << ibest << tab << setprecision(2) << translate[ibest] << tab
					<< setprecision(4) << bestview << tab << best << endl;

	delete[] cc;
	delete[] view_arr;
	
	currview = bestview;
	currshift = translate[ibest];
	
	delete[] translate;

	return best;
}

/**
@brief 	Searches a 2D/3D density map for a template using a specific view.
@param 	*ptemp		the template to be searched for.
@param 	view		view.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		mask for cross-correlation (ignored if NULL).
@param 	threshold	threshold value, if 0, threshold = FOMmax/2.
@param 	*pfit		view image with FOM block to hold results.
@return double		maximum FOM.

	The template is rotated to the view and cross-correlated.
	A view image is updated with the results where the correlation coefficients are higher.
	The views must be calculated externally to allow for custom sets.

**/
double		Bimage::search_volume_view(Bimage* ptemp, View view,
				double hires, double lores, Bimage* pmask, 
				double threshold, Bimage* pfit)
{
	Vector3<double> translate;
	Matrix3			mat = view.matrix();
	mat = mat.transpose();
	
	Bimage*			prot = NULL;

	if ( z > 1 ) prot = ptemp->rotate(size(), mat);
	else prot = ptemp->rotate_project(mat, translate, 0);
		
//	if ( z == 1 ) prot->project('z');
	
//	write_img("t.mrc", ptemp, 0);
//	write_img("t.mrc", prot, 0);
	
//	bexit(-1);

	check_if_same_size(prot);
	
	Bimage* 		pcc = cross_correlate(prot, hires, lores, pmask);
	
    long			i;
	double			fommax(0);
	
	for ( i=0; i<datasize; i++ ) {
		if ( (*pfit->next)[i] < (*pcc)[i] ) {
			pfit->next->set(i, (*pcc)[i]);
			pfit->set(i, view);
		}
		if ( fommax < (*pcc)[i] ) fommax = (*pcc)[i];
	}

	delete prot;
	delete pcc;
	
	return fommax;
}

/**
@brief 	Searches a 2D/3D density map for a template.
@param 	*ptemp		the template to be searched for.
@param 	*view		views.
@param 	alpha		rotation around view vector, <0 = use 2*PI (radians).
@param 	alpha_step	angular step size around view vector (radians).
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		mask for cross-correlation (ignored if NULL).
@param 	threshold	threshold value, if 0, threshold = FOMmax/2.
@return Bimage*		image with the best view at each voxel and FOM block.

	The template is rotated and cross-correlated to find the best fit.
	The views must be calculated externally to allow for custom sets.

**/
Bimage*		Bimage::search_volume(Bimage* ptemp, View* view, double alpha, 
				double alpha_step, double hires, double lores, 
				Bimage* pmask, double threshold)
{
	int				mode(0);		// Mode: 0=global, 1=refine, 2=symmetry
	long			i, nviews;
	double			alpha_min(0);
	double			alpha_max(TWOPI - alpha_step/2);
	double			a;
	double			fommin(1), fommax(0), fomavg(0), fomstd(0);
	View*			v;
	
	if ( alpha >= 0 ) {			// Case of refinement
		mode = 1;
		alpha_min = alpha - alpha_step;
		alpha_max = alpha + alpha_step;
	} else if ( alpha_step >= 1.999*M_PI ) mode = 2;	// Symmetry
	
	change_type(Float);
	ptemp->change_type(Float);
	rescale_to_avg_std(0, 1);
	ptemp->rescale_to_avg_std(0, 1);
	calculate_background();
	ptemp->calculate_background();

	Vector3<double>	ori(image->origin());
	
	if ( ori.sum() == 0 ) {
		ori = size()/2;
		origin(ori);
	}
	
	for ( nviews=0, v=view; v; v=v->next ) 
		for ( a=alpha_min; a<=alpha_max; a+=alpha_step ) nviews++;
		
	View*			view_arr = new View[nviews];
	
	for ( i=0, v=view; v; v=v->next ) 
		for ( a=alpha_min; a<=alpha_max; a+=alpha_step, i++ ) {
			view_arr[i] = *v;
			if ( mode < 2 ) view_arr[i][3] = a;
		}
		
	if ( verbose ) {
		cout << "Finding a template in a map:" << endl;
		cout << "Number of views:                " << nviews << endl;
		cout << "Alpha step size:                " << alpha_step*180.0/M_PI << endl;
		cout << "Origin:                         " << ori << endl;
	}

	Bimage*			pfit = new Bimage(Float, TView, size(), 1);
	pfit->next = new Bimage(Float, TSimple, size(), 1);
	pfit->sampling(sampling(0));

	if ( verbose & VERB_RESULT )
		cout << "View\tCCmax" << endl;

	for ( i=0; i<nviews; i++ ) {
		fommax = search_volume_view(ptemp, view_arr[i], hires, lores, pmask, threshold, pfit);
		if ( verbose )
			cout << view_arr[i] << tab << fommax << endl;
	}
	
	for ( i=0, fommax=-1, fommin=1, fomavg=fomstd=0; i<datasize; i++ ) {
		a = (*pfit->next)[i];
		if ( fommin > a ) fommin = a;
		if ( fommax < a ) fommax = a;
		fomavg += a;
		fomstd += a*a;
	}

	delete[] view_arr;
	
	if ( verbose ) {
		fomavg /= datasize;
		fomstd = fomstd/datasize - fomavg*fomavg;
		if ( fomstd > 0 ) fomstd = sqrt(fomstd);
		else fomstd = 0;
		cout << "FOM minimum and maximum:        " << fommin << " " << fommax << endl;
		cout << "FOM average:                    " << fomavg << endl;
		cout << "FOM standard deviation:         " << fomstd << endl << endl;
	}
	
	return pfit;
}


