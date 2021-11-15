/**
@file	mg_tomo_track.cpp
@brief	Functions to track fiducial markers in a tomographic series of images
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20210112
**/

#include "Bimage.h"
#include "mg_tomography.h"
#include "mg_tomo_track.h"
#include "mg_select.h"
#include "mg_img_proc.h"
#include "linked_list.h"
#include "qsort_functions.h"
#include "random_numbers.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			mg_marker_z_search(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref,
				Bimage* pgold, double thickness, int cc_type, int iz, double* z);
int			mg_z_matrix_update(Bmarker* model, int nmg, int nmark, double* z, double* w, int recenter);


double		fom_from_distance(double distance, double sigma)
{
	double		r = distance/sigma;
	
	double		fom = exp(-0.5*r*r);
	
	return fom;
}

int			marker_show(Bmarker* m1, Bmarker* m2)
{
	for ( ; m1 && m2; m1 = m1->next, m2 = m2->next )
		cout << m1->id << tab << m1->loc << tab << m2->loc << endl;
	
	return 0;
}

/**
@brief 	Generates a projection image of the global marker model.
@param 	size		size of image.
@param 	offset		offset to add to marker locations.
@param 	radius		gold fiducial marker radius.
@param 	*pmark		average marker image (can be NULL).
@param 	*mark		list of markers.
@return Bimage*		marker projection image.

	If the input marker image is NULL, synthetic particle images are 
	generated using the input radius.

**/
Bimage*		img_marker_projection(Vector3<long> size, Vector3<double> offset,
				double radius, Bimage* pmark, Bmarker* mark)
{
	double			outrad = sqrt(2.0)*radius;
	Bmarker*		m;
	Vector3<double>	loc;

	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	p->origin(p->size()/2);

	if ( pmark ) {
		p->fill(pmark->background(long(0)));
		for ( m = mark; m; m = m->next ) if ( m->sel ) {
			loc = m->loc - offset;
			p->place(0, pmark, loc, outrad, 1, 0, 1);
		}
	} else {
		for ( m = mark; m; m = m->next ) if ( m->sel ) {
			loc = m->loc - offset;
			p->sphere(loc, outrad, 2, FILL_USER, 1);
			p->sphere(loc, radius, 2, FILL_USER, -1);
		}
	}

	p->statistics();
	
//	write_img("mark_proj.map", p);
	
	return p;
}

double		mg_find_tilt_axis(Bmicrograph* mgp, Bmicrograph* mgn, Bmicrograph* mg_ref,
				double axis_angle, Bimage* pgold, double hi_res, double lo_res,
				double shift_limit, fft_plan planf, fft_plan planb)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_find_tilt_axis: axis_angle=" << axis_angle << endl;

	double			cc(0);
	Bmicrograph*	mg = micrograph_copy(mgp);
	
	mg->matrix = Matrix3(mg->tilt_angle, axis_angle);
//	mg->matrix = mg->matrix.transpose();
//	cout << mg->matrix << endl;
	
	mg_marker_update(mg, mg_ref->mark, mg_ref->origin, 1);
//	marker_show(mg_ref->mark, mg->mark);
	
	cc = mg_marker_shift(mg, pgold, hi_res, lo_res, shift_limit, planf, planb);
	
	micrograph_kill(mg);

	mg = micrograph_copy(mgn);
	
	mg->matrix = Matrix3(mg->tilt_angle, axis_angle);
//	mg->matrix = mg->matrix.transpose();
//	cout << mg->matrix << endl;
	
	mg_marker_update(mg, mg_ref->mark, mg_ref->origin, 1);
//	marker_show(mg_ref->mark, mg->mark);
	
	cc += mg_marker_shift(mg, pgold, hi_res, lo_res, shift_limit, planf, planb);
	cc /= 2;

	micrograph_kill(mg);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_find_tilt_axis: cc=" << cc << endl;
	
	return cc;
}

/**
@brief 	Finds the tilt axis for a tomographic tilt series.
@param 	*project		project parameter structure.
@param 	tilt			user-chosen tilt angle to use (radians).
@param 	axis_start		starting tilt axis angle.
@param 	axis_end		ending tilt axis angle.
@param 	axis_step		tilt axis angle step size (radians).
@param 	hi_res			high resolution limit for cross-correlation.
@param 	lo_res			low resolution limit for cross-correlation.
@param 	shift_limit		maximum micrograph shift to search for.
@return double			best correlation coefficient.

	The zero-tilt reference marker set must be defined.
	The micrographs closest to the positive and negative values of the
	given tilt angle is selected to find the tilt axis.
	The tilt axis for each of these micrographs is incremented from
	-PI to PI, markers generated from the reference seed, and 
	cross-correlated with the micrograph image.
	The tilt axis angle giving the best cross-correlation is chosen
	and assigned to all the micrographs.
	The markers for the chosen tilted micrographs are deleted.

**/
double		project_find_tilt_axis(Bproject* project, double tilt,
				double axis_start, double axis_end, double axis_step,
				double hi_res, double lo_res, double shift_limit)
{
	if ( !project->field ) return -1;
	
	double				a(0), abest(0), ccbest(-1);
	Bmicrograph*		mg;
	Bmicrograph*		mg_ref = field_find_low_tilt_mg_with_markers(project->field);
	Bmicrograph*		mgp = mg_ref;
	Bmicrograph*		mgn = mg_ref;
	
	if ( !mg_ref || !mg_ref->mark ) {
		cerr << "Error: No zero-tilt markers specified!" << endl;
		return 0;
	}

	for ( mg = project->field->mg; mg; mg = mg->next )
		if ( a < fabs(mg->tilt_angle) ) a = fabs(mg->tilt_angle);
	
	if ( a < tilt ) {
		cerr << "Error: The micrograph tilt angles are smaller than the requested tilt!" << endl;
		return 0;
	}
	
	// Choose micrographs close to the given tilt angle
	for ( mg = project->field->mg; mg; mg = mg->next ) {
		if ( fabs(mg->tilt_angle - tilt) < fabs(mgp->tilt_angle - tilt) ) mgp = mg;
		if ( fabs(mg->tilt_angle + tilt) < fabs(mgn->tilt_angle + tilt) ) mgn = mg;
	}
	
	if ( mgp && mgp->tilt_angle <= 0.05 ) {
		cerr << "Error: The postive test tilt is too small!" << endl;
		return 0;
	}
	
	if ( mgn && mgn->tilt_angle >= -0.05 ) {
		cerr << "Error: The negative test tilt is too large!" << endl;
		return 0;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_find_tilt_axis: mg_ref=" << mg_ref->id << " mgp=" << mgp->id << " mgn=" << mgn->id << endl;

	kill_list((char *) mgp->mark, sizeof(Bmarker));
	kill_list((char *) mgn->mark, sizeof(Bmarker));

	// Generate a composite gold particle image for cross-correlation
	int					cc_size = (int) (4*mg_ref->mark_radius);
	if ( cc_size < 30 ) cc_size = 30;
	Vector3<long>		gsize(cc_size, cc_size, 1);
	Bimage*				pgold = mg_composite_particle(mg_ref, gsize);
	
	if ( verbose ) {
		cout << "Finding the tilt axis angle:" << endl;
		cout << "Test tilt angle:                " << tilt*180/M_PI << endl;
		cout << "Micrograph tilt angles:         " << mgn->tilt_angle*180/M_PI
			<< " " << mgp->tilt_angle*180/M_PI << endl;
		cout << "Micrographs:                    " << mgn->id << " " << mgp->id << endl;
		cout << "Axis start, end and step:       " << axis_start*180/M_PI << " "
			<< axis_end*180/M_PI << " " << axis_step*180/M_PI << endl;
		cout << "Resolution limits:              " << hi_res << " - " << lo_res << endl << endl;
	}

	long				i, na(0);
	for ( a = axis_start; a < axis_end + 0.1*axis_step; a += axis_step ) na++;
	
	double*				ar = new double[na];
	double*				cc = new double[na];

	for ( i=0, a = axis_start; a < axis_end + 0.1*axis_step; a += axis_step, i++ ) {
		ar[i] = a;
		cc[i] = 0;
	}

	Vector3<long>	size = micrograph_get_size(mg_ref);
	
	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 1);
//	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 0);
//	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 0);

#ifdef HAVE_GCD
	dispatch_apply(na, dispatch_get_global_queue(0, 0), ^(size_t i){
		cc[i] = mg_find_tilt_axis(mgp, mgn, mg_ref, ar[i], pgold, hi_res, lo_res, shift_limit, planf, planb);
	});
#else
#pragma omp parallel for
	for ( i=0; i<na; i++ ) {
		cc[i] = mg_find_tilt_axis(mgp, mgn, mg_ref, ar[i], pgold, hi_res, lo_res, shift_limit, planf, planb);
	}
#endif

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	if ( verbose )
		cout << "Angle\tCC" << endl;
	for ( i=0; i<na; i++ ) {
		if ( ccbest < cc[i] ) {
			ccbest = cc[i];
			abest = ar[i];
		}
		if ( verbose )
			cout << ar[i]*180.0/M_PI << tab << cc[i] << endl;
	}
	
	delete[] ar;
	delete[] cc;
	delete pgold;
	
	if ( verbose ) {
		cout << "Best tilt axis angle:           " << abest*180.0/M_PI << " degrees" << endl;
		cout << "Alternative tilt axis angle:    " << angle_set_negPI_to_PI(abest-M_PI)*180.0/M_PI << " degrees" << endl;
		cout << "Best correlation coefficient:   " << ccbest << endl << endl;
	}
	
	kill_list((char *) mgp->mark, sizeof(Bmarker));
	kill_list((char *) mgn->mark, sizeof(Bmarker));
	mgp->mark = NULL;
	mgn->mark = NULL;
	
	for ( mg = project->field->mg; mg; mg = mg->next )
		mg->tilt_axis = abest;

	return ccbest;
}


/**
@brief 	Aligns markers in a tomographic series.
@param 	*project		project parameter structure.
@param 	hi_res			high resolution limit for cross-correlation.
@param 	lo_res			low resolution limit for cross-correlation.
@param 	shift_limit		maximum micrograph shift to search for.
@param 	thickness		estimated tomogram thickness (angstrom).
@param 	max_cycle		maximum number of iterations.
@param	target			target residual to terminate tracking.
@param	cc_type			indicates type of marker correlation: 0=real space, 1=cross correlation
@param	recenter		flag to recenter z coordinates.
@param 	paramfile		output parameter file name.
@return int				0, <0 on error.

	From the fiducial marker seed in the zero-degree tilt image, the 
	z coordinates of the markers as well as the image shift for 
	each micrograph is determined. The algorithm first attempts to 
	find the z-coordinate for each marker in an image by doing 
	real space correlations along a line determined by the tilt direction. 
	It then generates a projection image from the whole marker set at 
	the nominal tilt angle and cross-correlates it with the image to 
	find the shift. The process proceeds from the low-angle tilts to 
	higher tilts in both directions, using the lower dependence of the 
	low-tilt images on correct marker z-coordinates. This process is 
	iterated (typically 2-5 times) until the change in z-coordinates drops
	below one pixel on average or up to the maximum number of iterations. 
	The resolution limits are used in the cross-correlations.
	The shift limit prevents setting micrograph origin to far from the nominal origin.
	The thickness determines the extent of searching for the z coordinate of a marker.

**/
int			project_track_markers(Bproject* project, double hi_res, double lo_res, 
				double shift_limit, double thickness, int max_cycle,
				double target, int cc_type, int recenter, Bstring paramfile)
{
	if ( max_cycle > 20 ) max_cycle = 20;
	
	int					i, j, nmark, n, iref(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Bmicrograph*		mg_ref = field_find_low_tilt_mg_with_markers(field);
	Breconstruction*	rec = project->rec;
	Bmarker*			model = rec->mark;
	Bmarker*			mark = NULL;
	Bmarker*			mark2 = NULL;
	Bmarker*			model2 = NULL;
	
	if ( mg_ref->mark_radius < 1 )
		return error_show("Error: Cannot track markers with zero radius reference!", __FILE__, __LINE__);
	
	if ( target < 0.1 ) target = mg_ref->mark_radius/10;
	
	if ( !model )
		return error_show("No seed markers found!", __FILE__, __LINE__);
	
	if ( verbose ) {
		cout << "Tracking markers:" << endl;
		cout << "Maximum number of iterations:   " << max_cycle << endl;
		cout << "Marker radius:                  " << mg_ref->mark_radius << endl;
		cout << "Model origin:                   " << rec->origin << endl << endl;
	}
	
	if ( paramfile.length() < 1 ) paramfile = "btrack_param.star";
	Bstring				filename;
	
	// Generate a composite gold particle image for cross-correlation
	int					cc_size = (int) (4*mg_ref->mark_radius);
	if ( cc_size < 30 ) cc_size = 30;
	Vector3<long>		gsize(cc_size, cc_size, 1);
	Bimage*				pgold = mg_composite_particle(mg_ref, gsize);
	write_img("marker_track_ref.map", pgold, 0);
	
	for ( nmark=0, mark = model; mark; mark = mark->next, nmark++ ) ;
	model2 = markers_copy(model);

	int				nmg = field_count_micrographs(field);
	int				cycle, iz;
	double			R(100);
	
	double*			z = new double[nmark*nmg];
	double*			w = new double[nmg];
	Bmicrograph**	mglist = new Bmicrograph*[nmg];
	
	for ( i=0; i<nmg*nmark; i++ ) z[i] = -1e30;
	
	for ( i=0, mg = field->mg; mg; mg = mg->next, i++ ) {
		mglist[i] = mg;
		if ( mg == mg_ref ) iref = i;
		w[i] = 1/(1 + exp(0.2*(0.35 - fabs(mg->tilt_angle))));
	}

	int				imax = nmg - iref;
	if ( imax < iref ) imax = iref;

	Vector3<long>	size = micrograph_get_size(mg_ref);
	
//	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 1);
//	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 1);
	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 0);
	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 0);

	for ( cycle = 0; R > target && cycle < max_cycle; cycle++ ) {

		for ( mark=model, mark2=model2; mark; mark=mark->next, mark2=mark2->next )
			mark2->loc = mark->loc;
		
		for ( i=1; i<=imax; i++ ) {
			cout << "Processing images " << iref-i << " and " << iref+i << endl;
			iz = iref - i;
			if ( iz >= 0 ) {
				mg = mglist[iz];
				if ( mg->select ) {
					mg_marker_update(mg, model, rec->origin, 1);
					mg_marker_shift(mg, pgold, hi_res, lo_res, shift_limit, planf, planb);
					mg_marker_z_search(mg, model, rec->origin, pgold, thickness, cc_type, iz*nmark, z);
				}
			}
//			bexit(0);
			iz = iref + i;
			if ( iz < nmg ) {
				mg = mglist[iz];
				if ( mg->select ) {
					mg_marker_update(mg, model, rec->origin, 1);
					mg_marker_shift(mg, pgold, hi_res, lo_res, shift_limit, planf, planb);
					mg_marker_z_search(mg, model, rec->origin, pgold, thickness, cc_type, iz*nmark, z);
				}
			}
			mg_z_matrix_update(model, nmg, nmark, z, w, recenter);
		}
	
		for ( mg = field->mg; mg; mg = mg->next )
			mg_marker_update(mg, model, rec->origin, 1);

		R = 0;
		for ( n=0, mark=model, mark2=model2; mark; mark=mark->next, mark2=mark2->next ) if ( mark2->sel ) {
			R += (mark->loc - mark2->loc).length2();
			n++;
		}
		R = sqrt(R/n);
		cout << "Cycle " << cycle+1 << ":\tAverage change in positions = " << R << " (" << n << ")" << endl;
		
		if ( max_cycle > 1 ) {
			filename = paramfile.pre_rev('.') + Bstring(cycle+1, "_%02d.") + paramfile.post_rev('.');
			write_project(filename, project);
		}
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	kill_list((char *) model2, sizeof(Bmarker));
	
	cout << "Marker z coordinate matrix:\nMicrograph\tw";
	for ( j=0; j<nmark; j++ ) cout << tab << j+1;
	cout << endl;
	for ( i=0, mg=field->mg; i<nmg; i++, mg=mg->next ) if ( mg != mg_ref && mg->select ) {
		cout << mg->id << ":\t" << w[i];
		for ( j=0; j<nmark; j++ ) cout << tab << z[i*nmark+j];
		cout << endl;
	}
	
	delete[] z;
	delete[] w;
	delete[] mglist;
	delete pgold;

	if ( verbose ) {
		cout << "Marker\tx\ty\tz\tResidual" << endl;
		for ( mark = model; mark; mark = mark->next )
			cout << mark->id << tab << mark->loc[0] << tab << mark->loc[1] << tab << mark->loc[2] << tab << mark->res << endl;
		cout << endl;
	}
	
	return 0;
}

/**
@brief 	Aligns markers in a tomographic series.
@param 	*project		project parameter structure.
@param 	hi_res			high resolution limit for cross-correlation.
@param 	lo_res			low resolution limit for cross-correlation.
@param 	shift_limit		maximum micrograph shift to search for.
@param 	thickness		estimated tomogram thickness (angstrom).
@param 	max_cycle		maximum number of iterations.
@param	target			target residual to terminate tracking.
@param	cc_type			indicates type of marker correlation: 0=real space, 1=cross correlation
@param	recenter		flag to recenter z coordinates.
@param 	paramfile		output parameter file name.
@return int				0, <0 on error.

	From the fiducial marker seed in the zero-degree tilt image, the 
	z coordinates of the markers as well as the image shift for 
	each micrograph is determined. The algorithm first attempts to 
	find the z-coordinate for each marker in an image by doing 
	real space correlations along a line determined by the tilt direction. 
	It then generates a projection image from the whole marker set at 
	the nominal tilt angle and cross-correlates it with the image to 
	find the shift. The process proceeds from the low-angle tilts to 
	higher tilts in both directions, using the lower dependence of the 
	low-tilt images on correct marker z-coordinates. This process is 
	iterated (typically 2-5 times) until the change in z-coordinates drops
	below one pixel on average or up to the maximum number of iterations. 
	The resolution limits are used in the cross-correlations.
	The shift limit prevents setting micrograph origin to far from the nominal origin.
	The thickness determines the extent of searching for the z coordinate of a marker.

**/
int			project_track_markers_dual(Bproject* project, double hi_res, double lo_res, 
				double shift_limit, double thickness, int max_cycle,
				double target, int cc_type, int recenter, Bstring paramfile)
{
	if ( max_cycle > 20 ) max_cycle = 20;
	if ( thickness < 100 ) thickness = 1000;
	
	long				i, j, k, h, nmark, n, jref(0);
	Bstring				filename;
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Bmicrograph*		mgr = NULL;
	Bmicrograph*		mg_ref = field_find_low_tilt_mg_with_markers(field);
	Breconstruction*	rec = project->rec;
	Bmarker*			model = rec->mark;
	Bmarker*			mark = NULL;
	Bmarker*			mark2 = NULL;
	Bmarker*			model2 = NULL;
	
	if ( mg_ref->mark_radius < 1 )
		return error_show("Error: Cannot track markers with zero radius reference!", __FILE__, __LINE__);
	
	if ( target < 0.1 ) target = mg_ref->mark_radius/10;
	
	if ( !model )
		return error_show("No seed markers found!", __FILE__, __LINE__);
	
	if ( verbose ) {
		cout << "Tracking markers:" << endl;
		cout << "Maximum number of iterations:   " << max_cycle << endl;
		cout << "Marker radius:                  " << mg_ref->mark_radius << endl;
		cout << "Model origin:                   " << rec->origin << endl << endl;
	}
	
	if ( paramfile.length() < 1 ) paramfile = "btrack_param.star";
	
	// Generate a composite gold particle image for cross-correlation
	long				cc_size = (long) (4*mg_ref->mark_radius);
	if ( cc_size < 30 ) cc_size = 30;
	Vector3<long>		gsize(cc_size, cc_size, 1);
	Bimage*				pgold = mg_composite_particle(mg_ref, gsize);
	write_img("marker_track_ref.map", pgold, 0);
	
	for ( nmark=0, mark = model; mark; mark = mark->next, nmark++ ) ;
	model2 = markers_copy(model);

	long			nmg = project_count_micrographs(project);
	long			cycle;
	double			R(100);
	
	int*			seq = new int[nmg];
	double*			z = new double[nmark*nmg];
	double*			w = new double[nmg];
	Bmicrograph**	mglist = new Bmicrograph*[nmg];
	
	for ( i=0; i<nmg*nmark; i++ ) z[i] = -1e30;
	
	if ( verbose )
		cout << "Sequence:" << endl;
	for ( i=k=0, field = project->field; field; field = field->next ) {
		mgr = field_find_low_tilt_mg_with_markers(field);
		for ( n=0, mg = field->mg; mg; mg = mg->next, n++ )
			if ( mg == mgr ) jref = n;
		for ( j=1, h=0, mg = field->mg; mg; mg = mg->next, i++ ) {
			mglist[i] = mg;
			seq[i] = -1;
			if ( jref + j < n && h ) seq[i] = k + jref + j;
			else if ( jref - j >= 0 && h == 0 ) seq[i] = k + jref - j;
			if ( h ) j++;
			h = 1 - h;
			w[i] = 1/(1 + exp(0.2*(0.35 - fabs(mg->tilt_angle))));
			if ( verbose )
				cout << " " << seq[i];
		}
		k += n;
	}
	
	if ( verbose )
		cout << endl;

	Vector3<long>	size = micrograph_get_size(mg_ref);
	
//	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 1);
//	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 1);
	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 0);
	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 0);

	for ( cycle = 0; R > target && cycle < max_cycle; cycle++ ) {

		for ( mark=model, mark2=model2; mark; mark=mark->next, mark2=mark2->next )
			mark2->loc = mark->loc;
		
		for ( i=0; i<nmg; i++ ) if ( seq[i] >= 0 ) {
			cout << "Processing image " << seq[i] << endl;
			mg = mglist[seq[i]];
			if ( mg->select ) {
				mg_marker_update(mg, model, rec->origin, 1);
				mg_marker_shift(mg, pgold, hi_res, lo_res, shift_limit, planf, planb);
				mg_marker_z_search(mg, model, rec->origin, pgold, thickness, cc_type, seq[i]*nmark, z);
				mg_z_matrix_update(model, nmg, nmark, z, w, recenter);
			}
		}

		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				mg_marker_update(mg, model, rec->origin, 1);

		R = 0;
		for ( n=0, mark=model, mark2=model2; mark; mark=mark->next, mark2=mark2->next ) if ( mark2->sel ) {
			R += (mark->loc - mark2->loc).length2();
			n++;
		}
		R = sqrt(R/n);
		cout << "Cycle " << cycle+1 << ":\tAverage change in positions = " << R << " (" << n << ")" << endl;
		
		if ( max_cycle > 1 ) {
			filename = paramfile.pre_rev('.') + Bstring(cycle+1, "_%02d.") + paramfile.post_rev('.');
			write_project(filename, project);
		}
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	kill_list((char *) model2, sizeof(Bmarker));
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Marker z coordinate matrix:\nMicrograph\tw";
		for ( j=0; j<nmark; j++ ) cout << tab << j+1;
		cout << endl;
		for ( i=0, field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; i++, mg = mg->next ) if ( mg != mg_ref && mg->select ) {
				cout << mg->id << ":\t" << w[i];
				for ( j=0; j<nmark; j++ ) cout << tab << z[i*nmark+j];
				cout << endl;
			}
		}
	}
	
	delete[] seq;
	delete[] z;
	delete[] w;
	delete[] mglist;
	delete pgold;

	if ( verbose ) {
		cout << "Marker\tx\ty\tz\tResidual" << endl;
		for ( mark = model; mark; mark = mark->next )
			cout << mark->id << tab << mark->loc[0] << tab << mark->loc[1] << tab << mark->loc[2] << tab << mark->res << endl;
		cout << endl;
	}
	
	return 0;
}

double		marker_refine(Bmarker* mark, Bmarker* marklist, Bimage* p, Bimage* pgold,
				double radius, double hi_res, double lo_res)
{
	if ( !p->within_boundaries(mark->loc) ) return 0;

	Vector3<long>	extsize = pgold->size() * 2;
	if ( pgold->sizeZ() < 2 ) extsize[2] = 1;
	Vector3<double>	extori(extsize[0]/2, extsize[1]/2, extsize[2]/2);
	
	Bimage*			pmark = img_marker_projection(extsize, (mark->loc - extori),
								radius, pgold, marklist);
	
	Bimage*			pex = p->extract(0, (mark->loc - extori), extsize);
	pex->origin(extori);

	double			cc(0);
	mark->err = pex->find_shift(pmark, NULL, hi_res, lo_res, pex->sizeX()/4.0, 0.1, 1, cc);
	
	delete pmark;
	delete pex;
	
	mark->loc += mark->err;
	mark->res = mark->err.length();
	mark->fom = cc;
	if ( mark->fom < 0.0001 ) mark->fom = 0.0001;
	
	return cc;
}

/**
@brief 	Refines marker positions in a tomographic series.
@param 	*project	project parameter structure.
@param 	hi_res		high resolution limit for cross-correlation.
@param 	lo_res		low resolution limit for cross-correlation.
@return int			0, <0 on error.

	The area around a marker is extracted and cross-correlated with the
	corresponding projection from the 3D marker model.

**/
int			project_refine_markers(Bproject* project, double hi_res, double lo_res)
{
	int				i, nmg, nmark, nms(0), nmsi;

	Bfield*			field = project->field;
	Bmicrograph*	mg;
	Bmicrograph*	mg_ref = field_find_zero_tilt_mg(field);
	Bmarker*		mark = NULL;
	Bimage*			p = NULL;

	if ( !mg_ref->mark )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	int				cc_size = (int) (6*mg_ref->mark_radius);
	if ( cc_size < 40 ) cc_size = 40;
	Vector3<long>	gsize(cc_size, cc_size, 1);
	Bimage*			pgold = mg_composite_particle(mg_ref, gsize);
	write_img("marker_refine_ref.map", pgold, 0);
	
	for ( nmark=0, mark = mg_ref->mark; mark; mark = mark->next ) nmark++;
	double			mg_shift, avg_shift, mg_fom, avg_fom, cc;
	double*			mark_shift = new double[nmark];
	double*			mark_fom = new double[nmark];
	
	for ( i=0; i<nmark; i++ ) mark_shift[i] = mark_fom[i] = 0;
	
	if ( verbose )
		cout << "Image\tShift\tCC" << endl;
	for ( nmg=0, field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next, nmg++ ) if ( mg->select && mg->mark ) {
			p = read_img(mg->fmg, 1, mg->img_num);
			if ( !p )
				return error_show("project_refine_markers", __FILE__, __LINE__);
			nmsi = 0;
			mg_shift = mg_fom = 0;
			for ( i=0, mark = mg->mark; mark; mark = mark->next, i++ ) if ( mark->sel > 0 ) {
				cc = marker_refine(mark, mg->mark, p, pgold, mg->mark_radius, hi_res, lo_res);
				mark_shift[i] += mark->res;
				mg_shift += mark->res;
				mark_fom[i] += cc;
				mg_fom += cc;
				nmsi++;
			}
			delete p;
			nms += nmsi;
			if ( verbose ) 
				cout << mg->id << tab << mg_shift/nmsi << tab << mg_fom/nmsi << endl;
		}
	}
	
	if ( verbose )
		cout << "\nMarker\tShift\tCC" << endl;
	for ( avg_shift=avg_fom=0, i=0, mark=mg_ref->mark; i<nmark; i++, mark=mark->next ) {
		mark_shift[i] /= nmg;
		mark_fom[i] /= nmg;
		avg_shift += mark_shift[i];
		avg_fom += mark_fom[i];
		if ( verbose )
			cout << mark->id << tab << mark_shift[i] << tab << mark_fom[i] << endl;
	}
	
	if ( verbose ) {
		cout << endl << "Average shift:                  " << avg_shift/nms << endl;
		cout << "Average CC:                     " << avg_fom/nms << endl;
	}
	
	delete[] mark_shift;
	delete[] mark_fom;
	delete pgold;
	
	return 0;
}

/**
@brief 	Refines marker positions in a tomographic series.
@param 	*project	project parameter structure.
@param	id			marker identifier.
@param 	hi_res		high resolution limit for cross-correlation.
@param 	lo_res		low resolution limit for cross-correlation.
@return int			0, <0 on error.

	The area around a marker is extracted and cross-correlated with the
	corresponding projection from the 3D marker model.

**/
int			project_refine_one_marker(Bproject* project, int id, double hi_res, double lo_res)
{
	int				nmg, nms(0);

	Bfield*			field = project->field;
	Bmicrograph*	mg;
	Bmicrograph*	mg_ref = field_find_zero_tilt_mg(field);
	Bmarker*		mark = NULL;
	Bimage*			p = NULL;

	if ( !mg_ref->mark )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	int				cc_size = (int) (6*mg_ref->mark_radius);
	if ( cc_size < 40 ) cc_size = 40;
	Vector3<long>	gsize(cc_size, cc_size, 1);
	Bimage*			pgold = mg_composite_particle(mg_ref, gsize);
	write_img("marker_refine_ref.map", pgold, 0);
	
	double			cc;
	double			mark_shift(0);
	double			mark_fom(0);
	
	for ( nmg=0, field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next, nmg++ ) if ( mg->select && mg->mark ) {
			p = read_img(mg->fmg, 1, mg->img_num);
			if ( !p )
				return error_show("project_refine_one_marker", __FILE__, __LINE__);
			for ( mark = mg->mark; mark && mark->id != id; mark = mark->next ) ;
			if ( mark &&  mark->sel > 0 ) {
				cc = marker_refine(mark, mg->mark, p, pgold, mg->mark_radius, hi_res, lo_res);
				mark_shift += mark->res;
				mark_fom += cc;
				nms++;
			}
			delete p;
		}
	}
	
	if ( verbose ) {
		cout << "Marker " << id << ":" << endl;
		cout << "Average shift:                  " << mark_shift/nms << endl;
		cout << "Average CC:                     " << mark_fom/nms << endl << endl;
	}
	
	delete pgold;
	
	return 0;
}


/**
@brief 	Refines marker positions in a tomographic micrograph.
@param 	*mg			micrograph.
@param 	*pgold		marker reference image (can be NULL).
@param 	hi_res		high resolution limit for cross-correlation.
@param 	lo_res		low resolution limit for cross-correlation.
@return int			0.

	The area around a marker is extracted and cross-correlated with the
	corresponding projection from the 3D marker model.

**/
int			mg_refine_markers(Bmicrograph* mg, Bimage* pgold, double hi_res, double lo_res)
{
	if ( !mg->mark )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	int				nmark, nms(0);
	Bmarker*		mark = NULL;
	Bimage*			p = NULL;
	
	double			mg_shift, mg_fom, cc;
	
	p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p )
		return error_show("mg_refine_markers", __FILE__, __LINE__);
	mg_shift = mg_fom = 0;
	for ( nmark=0, mark = mg->mark; mark; mark = mark->next, nmark++ ) if ( mark->sel > 0 ) {
		cc = marker_refine(mark, mg->mark, p, pgold, mg->mark_radius, hi_res, lo_res);
		mg_shift += mark->res;
		mg_fom += cc;
		nms++;
	}
	delete p;
	
	if ( verbose & VERB_FULL )
		cout << mg->id << tab << mg_shift/nms << tab << mg_fom/nms << endl;
	
	return 0;
}

/**
@brief 	Refines all the alignment parameters.
@param 	*project	micrograph project.
@param 	iter		maximum number of iterations.
@param 	tol			tolerance for exit condition.
@param 	refop		string holding sequence of refinement operations.
@return double		marker RMSD.
**/
double		project_refine(Bproject* project, int iter, double tol, Bstring refop)
{
	long		i, j;
	double		R = 1e30, dR = 1e30;

	if ( verbose ) {
		cout << "Refining model and orientation parameters:" << endl;
		cout << "Model origin:                   " << project->rec->origin << endl;
		cout << "Sequence of refinements:        " << refop << endl << endl;
	}
	
	R = project_tomo_residuals(project, 0);

	if ( verbose )
		cout << "Cycle\tR\tdR" << endl;
	for ( i = 1; i <= iter && R > 0.1 && fabs(dR) > tol; i++ ) {
		dR = R;
		for ( j=0; j<refop.length(); j++ ) {
			switch ( refop[j] ) {
				case 'z': project_refine_z(project); break;
				case 'v': project_refine(project, 1, 0, 0); break;
				case 'o': project_refine(project, 0, 1, 0); break;
				case 's': project_refine(project, 0, 0, 1); break;
				case 'a': project_refine(project, 1, 1, 1); break;
			}
		}
		R = project_tomo_residuals(project, 0);
		dR -= R;
		if ( verbose )
			cout << i << tab << R << tab << dR << endl;
	}

	project_calculate_angles(project);

	if ( verbose )
		cout << endl;

	return R;
}

double		tomo_z_residual(Bproject* project, Bmarker* mark)
{
	Bfield*				field = project->field;
	Bmicrograph*		mg = field->mg;
	Bmarker*			m;
	Vector3<double>		oriref = project->rec->origin;
	Vector3<double>		loc, rloc;
	Matrix3				mat(1);
	long				n(0);
	double				d, R(0), w(0);

	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			mat = mg->matrix.transpose();	// Handedness
			for ( m = mg->mark; m && m->id != mark->id; m = m->next ) ;
			if ( m && m->id == mark->id && m->sel > 0 ) {
				loc = m->loc - mg->origin;
				rloc = mark->loc - oriref;
//				d = loc[0] - rloc[0]*mat[0] - rloc[1]*mat[1] - rloc[2]*mat[2];
				d = loc[0] - rloc.scalar(mat[0]);
				R += d*d*mark->fom;
//				d = loc[1] - rloc[0]*mat[3] - rloc[1]*mat[4] - rloc[2]*mat[5];
				d = loc[1] - rloc.scalar(mat[1]);
				R += d*d*mark->fom;
				n++;
				w += mark->fom;
			}
		}
	}
	
	if ( w ) R = sqrt(R/(2*w));

	return R;
}

/**
@brief 	Refines the z coordinates of the marker model.
@param 	*project	micrograph project.
@return double		overall RMS change.
**/
double		project_refine_z(Bproject* project)
{
	Bmarker*			mr;
	Bmarker*			model = project->rec->mark;
	long				n(0);
	double				oldz, zinc, d, Rz, Rp, Rm, R(0);

	if ( !model )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	if ( verbose & VERB_FULL )
		cout << "Marker\tOld z\tNew z\tR" << endl;
	for ( mr = model; mr; mr = mr->next ) if ( mr->sel ) {
		oldz = mr->loc[2];
		Rz = tomo_z_residual(project, mr);
		for ( zinc = 1; zinc > 0.001; ) {
			mr->loc[2] += zinc;
			Rp = tomo_z_residual(project, mr);
			if ( Rz > Rp ) {
				Rz = Rp;
			} else {
				mr->loc[2] -= 2*zinc;
				Rm = tomo_z_residual(project, mr);
				if ( Rz > Rm ) {
					Rz = Rm;
				} else {
					mr->loc[2] += zinc;
					zinc *= 0.7;
				}
			}
		}
		if ( verbose & VERB_FULL )
			cout << mr->id << tab << oldz << tab << mr->loc[2] << tab << Rz << endl;
		d = mr->loc[2] - oldz;
		R += d*d;
		n++;
	}
	
	R = sqrt(R/n);
	
	if ( verbose & VERB_PROCESS )
		cout << "RMS change in z:            " << R << endl;
	
	return R;
}

/**
@brief 	Refines the views from the marker positions and marker model.
@param 	*project	micrograph project.
@param 	do_view		refine micrograph views.
@param 	do_origin	refine micrograph origins.
@param 	do_scale	refine micrograph scales.
@return double		best residual.

	Requires the matrices in the micrograph structures to be defined.

**/
double		project_refine(Bproject* project, int do_view, int do_origin, int do_scale)
{
//	if ( ( do_origin < 1 ) && ( do_scale < 1 ) && ( do_view < 1 ) ) return 0;
	if ( !( do_origin || do_scale || do_view ) ) return 0;
	
//	cout << "refining " << do_view << do_origin << do_scale << endl;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec = project->rec;
	double				att(1e-6), beta(100);
	int					accept, ncyc(0), n(0), nacc;
	double				irm = 1.0/get_rand_max();
	double				R, pR, bR(0), d;
	Bmarker*			model = rec->mark;

	View				v, vt, bv;
	Vector3<double>		ori, oritemp, oribest;
	Vector3<double>		scale, scaletemp, scalebest;
	
	random_seed();
	
	if ( verbose & VERB_PROCESS )
		cout << "Image\tdx\tdy\tsx\tsy\tvx\tvy\tvz\tva\tR" << endl;
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			v = View(mg->matrix);
			ori = mg->origin;
			scale = mg->scale;
			v.normalize();
			bR = mg_tomo_residuals(mg, model, rec->origin);
			n = nacc = ncyc = 0;
			R = pR = bR = 1e30;
			while ( ncyc < 1000 && bR > 0.001 ) {
				oritemp = ori;
				scaletemp = scale;
				vt = v;
				if ( do_view ) {
					vt[0] = v[0] + att*(random()*irm*2.0 - 1);
					vt[1] = v[1] + att*(random()*irm*2.0 - 1);
					vt[2] = v[2] + att*(random()*irm*2.0 - 1);
					if ( vt[2] < 0 ) vt[2] = -vt[2];		// Is this good ?????
					vt[3] = v.angle() + att*(random()*irm*2.0 - 1);
					vt.normalize();
					mg->matrix = vt.matrix();
				}
				if ( do_origin ) {
					oritemp = ori + vector3_xy_random_gaussian(0.0, 1.0);
					mg->origin = oritemp;
				}
				if ( do_scale ) {
					scaletemp = scale + vector3_xy_random_gaussian(0.0, 0.001);
					mg->scale = scaletemp;
				}
				R = mg_tomo_residuals(mg, model, rec->origin);
				if ( bR > R ) {
					bR = R;
					bv = vt;
					oribest = oritemp;
					scalebest = scaletemp;
					ncyc = 0;
				}
				accept = 1;
				if ( R > pR ) {
					d = random()*irm;
					if ( exp(beta*(pR-R)) < d ) accept = 0;
				}
				if ( accept ) {
					pR = R;
					v = vt;
					ori = oritemp;
					scale = scaletemp;
					nacc++;
				} else {
					R = pR;
					vt = v;
					oritemp = ori;
					scaletemp = scale;
				}
				ncyc++;
				n++;
			}
			mg->matrix = bv.matrix();
			mg->origin = oribest;
			mg->scale = scalebest;
			if ( verbose & VERB_PROCESS )
				cout << mg->id << tab << bv[0] << tab << bv[1] << tab << bv[2] << tab << 
					bv.angle()*180/M_PI << tab << acos(bv[2])*180/M_PI << tab << bR << tab << nacc*100.0/n << endl;
		}
	}
	
	return bR;
}

/*double		project_tilt_axis_from_markers(Bproject* project)
{
	long				nmark(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Bmarker*			mark;
	
	long					h, i, j;
	Matrix3					a1(3,3);
	vector<Matrix3>			a;
	vector<double>			b1(3,0);
	vector<vector<double>>	b;
	Vector3<double>			loc, ori, n;

	cout << "Marker\tnx\tny\tnz" << endl;
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		mg = field_find_zero_tilt_mg(project->field);
		ori = mg->origin;
		for ( mark = mg->mark; mark; mark = mark->next ) {
			nmark++;
			a.push_back(a1);
			b.push_back(b1);
		}
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			for ( h=0, mark = mg->mark; mark; mark = mark->next, ++h ) {
				loc = mark->loc - ori;
				for ( i=0; i<3; i++ ) {
					for ( j=0; j<=i; j++ )
						a[h][j][i] += loc[i]*loc[j];
					b[h][i] += loc[i];
				}
			}
		}
		mg = field_find_zero_tilt_mg(field);
		for ( h=0, mark = mg->mark; mark; mark = mark->next, ++h ) {
			n = a[h].plane_normal(b[h]);
			cout << mark->id << tab << n << endl;
		}
		a.clear();
		b.clear();
	}
	
	return 0;
}
*/

Vector3<double>	fit_line(vector<Vector3<double>> d)
{
	long				i, j;
	Vector3<double>		v, avg;
//	Matrix				a(3,3);
	Matrix3				a;
	vector<double>		b(3,0);
//	Vector3<double>		b;

	for ( auto it=d.begin(); it != d.end(); ++it )
		avg += *it;
		
	avg /= d.size();

	for ( auto it=d.begin(); it != d.end(); ++it ) {
		v = *it - avg;
		cout << v << endl;
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<3; j++ )
				a[j][i] += v[i]*v[j];
			b[i] += v[i];
		}
	}

	a.singular_value_decomposition();

	cout << endl << a << endl;
	cout << b << endl;

//	return a.plane_normal(b);
	return Vector3<double>(b);
}

/**
@brief 	Calculates the tilt axis from the marker trajectories.
@param 	*project	micrograph project.
@return double		best residual.

	If the rotation is around a single axis, the trajectory of a marker should lie in a plane.

**/
double		project_tilt_axis_from_markers(Bproject* project)
{
	long				nmark(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Bmarker*			mark;
	
	long				i;
	Vector3<double>		loc;
	

	cout << "Marker\tnx\tny\tnz" << endl;
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		mg = field_find_zero_tilt_mg(project->field);
		for ( nmark=0, mark = mg->mark; mark; mark = mark->next )
			if ( nmark < mark->id ) nmark = mark->id;
		for ( i=1; i<=nmark; ++i ) {
			vector<Vector3<double>>	m;
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
				for ( mark = mg->mark; mark && mark->id != i; mark = mark->next ) ;
				if ( mark ) {
					loc = mark->loc - mg->origin;
					m.push_back(loc);
				}
			}
			Vector3<double>		v = fit_line(m);
			cout << i << tab << v << endl;
		}
	}
	
	return 0;
}

/**
@brief 	Finds the shift of a micrograph wrt the reference using markers.
@param 	*mg			micrograph.
@param 	*pgold		marker reference image (can be NULL).
@param 	hi_res		high resolution limit for cross-correlation.
@param 	lo_res		low resolution limit for cross-correlation.
@param 	shift_limit	maximum micrograph shift to search for.
@param 	planf		FFT forward plan.
@param 	planb		FFT backward plan.
@return double		correlation coefficient.

	The marker locations and marker radius must already be set.
	The correlation coefficient for the correlation between the micrograph
	and the 2D image generated from the markers is retuned.

**/
double		mg_marker_shift(Bmicrograph* mg, Bimage* pgold, double hi_res,
				double lo_res, double shift_limit, fft_plan planf, fft_plan planb)
{
	if ( !mg->mark )
		return error_show("No markers found!", __FILE__, __LINE__);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_marker_shift: mg=" << mg->id << endl;
	
	Bmarker*		mark = NULL;
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p )
		return error_show("mg_marker_shift", __FILE__, __LINE__);
	
	p->origin(mg->origin);
	
	Vector3<double>	offset;
	
	if ( shift_limit < 0 ) shift_limit = 0.2*p->sizeX();

	Bimage*			pmark = img_marker_projection(p->size(), offset, mg->mark_radius, pgold, mg->mark);

	double			cc(0);
	Vector3<double>	shift = p->find_shift(pmark, NULL, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);

//	write_img("pmark.pif", pmark);
	delete pmark;
	delete p;
	
	mg->origin += shift;
	for ( mark = mg->mark; mark; mark = mark->next ) if ( mark->sel )
		mark->loc += shift;
		
	if ( verbose & VERB_PROCESS )
		cout << "Image " << mg->img_num << " origin:\t" << mg->origin[0] << tab << mg->origin[1] << tab << "CC: " << cc << endl;
	
	return cc;
}

double		marker_fom(Bimage* p, Vector3<double> loc, double radius)
{
	double		sf, s;
	double		f = p->density(0, loc, radius, sf);
	double		d = p->density(0, loc, 2*radius, s);
	
//	return fabs(f/d - 1);
//	return sf/s - 1;
	return f - d;
}

double		mg_marker_z_search_one(Bmicrograph* mg, Bimage* p, Bmarker* modmark, Vector3<double> oriref,
				Bimage* pgold, double thickness, int cc_type, int i, double* z, fft_plan planf, fft_plan planb)
{
	Bmarker*		mark;
	
	for ( mark = mg->mark; mark && mark->id != modmark->id; mark = mark->next ) ;
	
	if ( !mark || !mark->sel ) return 0;
	
	Vector3<double>	mloc = modmark->loc - oriref;
	
	z[i] = mloc[2];
	mark->fom = 0;
	
	double			cc(-1);
	double			dz(1/sin(fabs(mg->tilt_angle)));		// Approximately one pixel in xy
	double			z_range(thickness/mg->pixel_size[0]);		// Thickness in voxels
	double			z_start = mloc[2] - z_range;
	double			z_end = mloc[2] + z_range;
	Matrix3			mat = mg->matrix;
	Vector3<long>	gsize(pgold->size());
	Vector3<double>	gori(pgold->image->origin());
	Vector3<double>	extori;
	Bimage*			pmark;
	Bimage*			pc;
	
	if ( z_range > 100*dz ) z_range = 100*dz;
	
	for ( mloc[2] = z_start; mloc[2] <= z_end; mloc[2] += dz ) {
		mark->loc = mg_location_from_3D_model(mloc, mat, mg->origin);
		if ( mark->loc[0] > gori[0] && mark->loc[0] < p->sizeX() - gori[0] &&
				mark->loc[1] > gori[1] && mark->loc[1] < p->sizeY() - gori[1] ) {
			extori = mark->loc - gori;
			pmark = p->extract(0, extori, gsize);
			pmark->origin(pmark->size()/2);
			if ( cc_type < 1 ) {
				cc = pmark->correlate(pgold, 0, pmark->sizeX()/2.0, NULL, 0);
//				cc = pmark->correlate(pgold);
			} else {
				pc = pmark->cross_correlate(pgold, 4*pmark->sampling(0)[0], 0, planf, planb);
				cc = pc->maximum();
				delete pc;
			}
//			write_img("m.map", pmark);
			delete pmark;
			if ( mark->fom < cc ) {
				mark->fom = cc;
				z[i] = mloc[2];
			}
		}
	}
	
//	mark->fom = marker_fom(p, mark->loc, mg->mark_radius)/pgold->image->FOM();
	if ( mark->fom < 0.0001 ) mark->fom = 0.0001;
	modmark->fom = mark->fom;

	if ( verbose & VERB_FULL )
		cout << " " << z[i];

	return cc;
}

/*
@brief 	Refines the z-coordinate of a global marker model.
@param 	*mg			micrograph.
@param 	*model		marker model.
@param 	oriref		reference origin.
@param 	*pgold		marker reference image.
@param 	thickness	estimated tomogram thickness (angstrom).
@param	cc_type		indicates type of correlation: 0=real space, 1=cross correlation
@param 	iz			starting index in the z-matrix.
@param 	*z			z-matrix.
@return int			0.
**/
int			mg_marker_z_search(Bmicrograph* mg, Bmarker* model, Vector3<double> oriref,
				Bimage* pgold, double thickness, int cc_type, int iz, double* z)
{
	if ( fabs(mg->tilt_angle) < 0.001 ) {
		cerr << "Warning: micrograph " << mg->id << " tilt angle is too small! (" << mg->tilt_angle << ")" << endl;
		return -1;
	}
	
	if ( thickness < 1000 ) thickness = 1000;		// Minimum thickness = 1000 angstrom (0.1 um)
	if ( thickness > 100000 ) thickness = 100000;	// Maximum thickness = 100000 angstrom (10 um)
	
	if ( !mg->mark )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	Bmarker*		modmark = NULL;

	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p )
		return error_show("mg_marker_z_search", __FILE__, __LINE__);

	long			i, cc_size = (long) (4*mg->mark_radius);
	if ( cc_size < 30 ) cc_size = 30;
	double			dz(1/sin(fabs(mg->tilt_angle)));		// Approximately one pixel in xy
	Vector3<long>	gsize;
	
	int				local_gold_particle(0);
	if ( pgold ) {
		gsize = pgold->size();
	} else {
		local_gold_particle = 1;
		gsize = Vector3<long>(cc_size, cc_size, 1);
		pgold = img_gold_particle(gsize, mg->mark_radius);
	}
	
	p->sampling(mg->pixel_size);
	pgold->sampling(p->sampling(0));
	
	long			nmark;
	for ( nmark=0, modmark = model; modmark; modmark = modmark->next, nmark++ ) ;
	
	Bmarker**		marr = new Bmarker*[nmark];
	for ( i=0, modmark = model; i<nmark && modmark; i++, modmark = modmark->next )
		marr[i] = modmark;

	fft_plan		planf = fft_setup_plan(gsize, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(gsize, FFTW_BACKWARD, 1);
	
	if ( verbose )
		cout << "Image " << mg->img_num << " z search: " << dz << endl;
#ifdef HAVE_GCD
	dispatch_apply(nmark, dispatch_get_global_queue(0, 0), ^(size_t i){
		mg_marker_z_search_one(mg, p, marr[i], oriref, pgold, thickness, cc_type, i+iz, z, planf, planb);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nmark; i++ ) {
		mg_marker_z_search_one(mg, p, marr[i], oriref, pgold, thickness, cc_type, i+iz, z, planf, planb);
	}
#endif

	if ( verbose & VERB_FULL )
		cout << endl;

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);

	delete p;
	delete[] marr;

	if ( local_gold_particle ) delete pgold;
	
	return 0;
}

/*
@brief	Updates the z-coordinate of a global marker model.
@param	model			reference marker list.
@param 	nmg				number of micrographs.
@param 	nmark			number of markers.
@param	*z				z-matrix.
@param	*w				array for micrograph weights.
@param	recenter		flag to recenter z coordinates.
@return	int				0.

	Calculates average and deviation of the z-coordinate of a global 
	marker model over all of the available z measurements.
	The ultimate z-value accepted is the median.
**/
int			mg_z_matrix_update(Bmarker* model, int nmg, int nmark, double* z, double* w, int recenter)
{
	if ( !model )
		return error_show("No markers found!", __FILE__, __LINE__);
	
	int			i, j, k, n, m;
	double		tw, avg_z, std_z, med_z, avg_all_z(0);
	double*		za = new double[nmg];
	Bmarker*	mark;
	
	for ( i=0; i<nmg; i++ ) za[i] = 0;
	
	if ( verbose )
		cout << "Marker\tZmed\tZavg\tZstd\tFOM\tCount" << endl;
	for ( i=m=0, mark = model; mark; i++, mark = mark->next ) if ( mark->sel ) {
		for ( avg_z=std_z=tw=0, j=n=0; j<nmg; j++ ) {
			k = j*nmark + i;
			if ( z[k] > -1e10 ) {
				avg_z += w[j]*z[k];
				std_z += w[j]*z[k]*z[k];
				za[n] = z[k];
				tw += w[j];
				n++;
			}
		}
		med_z = 0;
		if ( n ) {
			avg_z /= tw;
			std_z = std_z/tw - avg_z*avg_z;
			if ( std_z > 0 ) std_z = sqrt(std_z);
			else std_z = 0;
			qsort((void *) za, n, sizeof(double), QsortLargeToSmallDouble);
			med_z = za[n/2];
//			mark->loc.z = avg_z;
//			mark->loc[2] = med_z;
			mark->loc[2] = 0.5*(med_z + mark->loc[2]);
//			mark->loc[2] = 0.5*(med_z + avg_z);
			mark->res = std_z;
//			mark->fom = fom_from_distance(std_z, 10);
//			if ( mark->fom < 0.0001 ) mark->fom = 0.0001;
//			mark->fom = 1/(1 + std_z + fabs(avg_z - med_z));
		}
		avg_all_z += mark->loc[2];
		m++;
		if ( verbose )
			cout << setprecision(2) << mark->id << tab << med_z << tab << 
				avg_z << tab << std_z << tab << setprecision(5) << mark->fom << tab << n << endl;
	}
//	if ( verbose )
//		cout << endl;

	delete[] za;

	avg_all_z /= m;
	
	if ( recenter ) {
		if ( verbose )
			cout << "Recentering z coordinates: " << avg_all_z << endl;
		
		if ( verbose & VERB_FULL )
			cout << "Marker\tz\tFOM" << endl;
		for ( mark = model; mark; mark = mark->next ) {
			if ( mark->sel ) mark->loc[2] -= avg_all_z;
			else mark->loc[2] = mark->fom = 0;
			if ( verbose & VERB_FULL )
				cout << mark->id << tab << mark->loc[2] << tab << mark->fom << endl;
		}
	}
	
	if ( verbose )
		cout << endl;
	
	return 0;
}

/**
@brief 	Transfers the seed markers from the first to the second project.
@param 	*project	project with seed markers in first field.
@param 	rot_start	starting rotation angle.
@param 	rot_end		final rotation angle.
@param 	rot_step	angular search step size.
@param 	hi_res		high resolution limit for cross-correlation.
@param 	lo_res		low resolution limit for cross-correlation.
@param 	shift_limit	maximum micrograph shift to search for.
@return int			number of markers.

	The markers from the first series are rotated around the micrograph
	origin by the search angle. The markers are then used to generate
	an image with synthetic markers, and this image is cross-correlated 
	with the zero-tilt micrograph of the second series. The search angle
	giving the best correlation coefficient is selected and the seed marker
	locations for the second series are calculated.

**/
int			project_transfer_seed(Bproject* project, 
				double rot_start, double rot_end, double rot_step, 
				double hi_res, double lo_res, double shift_limit)
{
	Vector3<double>		origin, bestshift;
	double				cc(-1), bestcc(-1), bestangle(0);
	Bfield*				field1 = project->field;
	Bfield*				field2 = field1->next;
	
	Bmicrograph*		mg1 = field_find_zero_tilt_mg(field1);
	if ( !mg1->mark ) for ( mg1 = field1->mg; mg1 && !mg1->mark; mg1 = mg1->next ) ;
	
	if ( !mg1 ) {
		cerr << "Error: No markers defined for the first tilt series!" << endl;
		return -1;
	}

	if ( hi_res < 2*mg1->pixel_size[0] ) hi_res = 2*mg1->pixel_size[0];
	if ( lo_res < 2*hi_res ) lo_res = 100*hi_res;
	if ( shift_limit < 10 ) shift_limit = 10;
	
	Bmicrograph*		mg2 = field_find_zero_tilt_mg(field2);

	if ( !mg2 ) {
		cerr << "Error: No micrographs found in the second tilt series!" << endl;
		return -1;
	}
	
	if ( mg1->matrix.determinant() < 0.5 )
		project_mg_tilt_to_matrix(project);

	mg2->mark_radius = mg1->mark_radius;
	
	// Generate a composite gold particle image for cross-correlation
	int					cc_size = (int) (4*mg1->mark_radius);
	if ( cc_size < 30 ) cc_size = 30;
	Vector3<long>		gsize(cc_size, cc_size, 1);
	Bimage*				pgold = mg_composite_particle(mg1, gsize);
	
	Transform			t;
	t.origin = mg1->origin;
	
	if ( verbose ) {
		cout << "Finding the best angle to transform the first seed to the second" << endl;
		cout << "Micrograph 1:                   " << mg1->id << endl;
		cout << "Micrograph 2:                   " << mg2->id << endl;
		cout << "Resolution:                     " << hi_res << " - " << lo_res << endl;
		cout << "Shift limit:                    " << shift_limit << endl;
	}
	
	Vector3<long>	size = micrograph_get_size(mg1);
	
	fft_plan		planf = fft_setup_plan(size[0], size[1], 1, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(size[0], size[1], 1, FFTW_BACKWARD, 1);

	if ( verbose )
		cout << "Angle\tOriX\tOriY\tCC" << endl;
	for ( t.angle = rot_start; t.angle <= rot_end; t.angle += rot_step ) {
		kill_list((char *) mg2->mark, sizeof(Bmarker));
		mg2->mark = markers_copy(mg1->mark);
		marker_transform(mg2->mark, t);
		mg2->origin = mg1->origin;
		cc = mg_marker_shift(mg2, pgold, hi_res, lo_res, shift_limit, planf, planb);
		if ( bestcc < cc ) {
			bestcc = cc;
			bestangle = t.angle;
			bestshift = mg2->origin - mg1->origin;
		}
		if ( verbose )
			cout << t.angle*180.0/M_PI << tab << mg2->origin[0] << tab << mg2->origin[1] << tab << cc << endl;
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
	
	t.angle = bestangle;
	t.trans = bestshift;
	kill_list((char *) mg2->mark, sizeof(Bmarker));
	mg2->mark = markers_copy(mg1->mark);
	int				nmark = marker_transform(mg2->mark, t);
	
	if ( verbose ) {
		cout << "Best angle:                     " << bestangle*180.0/M_PI << " degrees" << endl;
		cout << "Best shift:                     " << bestshift[0] << " " << bestshift[1] << endl;
		cout << "Correlation coefficient:        " << bestcc << endl << endl;
	}

	field2->origin = t * mg1->origin;

//	t = marker_find_transform(mg2->mark, mg1->mark, mg1->origin);
//	cout << "axis=" << t.axis*180/M_PI << tab << "angle=" << t.angle*180/M_PI << endl;
	
	field2->matrix = Matrix3(t.axis, t.angle);
	cout << "Field origin: " << field2->origin << endl;
	cout << "Field matrix: " << field2->matrix << endl;
	for ( mg2 = field2->mg; mg2; mg2 = mg2->next ) {
		mg2->origin = field2->origin;
//		mg2->tilt_axis = mg1->tilt_axis;
//		mg2->matrix = mg2->matrix * field2->matrix.transpose();	
		mg2->matrix = field2->matrix.transpose() * mg2->matrix;	
	}
	
//	if ( project->rec && project->rec->mark )
//		marker_transform(project->rec->mark, t);
	
	delete pgold;
	
	return nmark;
}

/**
@brief 	Calculates the transformation of the second tilt series to fit the first.
@param 	*project	project with seed markers.
@return int			number of markers compared.

	The 3D marker coordinates from the second series are fitted to those
	of the first series to determine the rotation matrix and shifts.
	The 3D marker locations and the micrograph orientations and origins
	of the second series are then adkjusted to correspond to the first.
	Restrictions: The first two fields should contain the two tilt series
	and the first two reconstructions the corresponding 3D marker sets.

**/
int			project_transform_dual(Bproject* project)
{
	Breconstruction*	rec1 = project->rec;
	Breconstruction*	rec2 = project->rec->next;
	Bmarker*			mark1 = rec1->mark;
	Bmarker*			mark2 = rec2->mark;
		
	Transform			t = marker_find_transform(mark2, mark1, rec2->origin);

	Matrix3				mat = Matrix3(t.axis, t.angle);

	if ( verbose ) {
		cout << "Marker set 2 transformation:" << endl;
		cout << "Axis:                           " << t.axis << endl;
		cout << "Rotation angle:                 " << t.angle*180.0/M_PI << " degrees" << endl;
		cout << "Scale:                          " << t.scale << endl;
		cout << "Origin:                         " << t.origin << endl;
		cout << "Translation:                    " << t.trans << endl;
		matrix3_show_hp(mat);		
	}

	mat = t.scale * mat;
	
	long				n;
	double				d, R(0);
	Vector3<double>		loc, shift_avg;
	Vector3<double>		shift = t.origin*t.scale + t.trans;

	// Transform the markers to agree with the eventual reconstruction
	for ( mark2 = rec2->mark; mark2; mark2 = mark2->next ) {
		loc = mark2->loc - rec2->origin;
		loc = mat * loc;
		loc += shift;
		mark2->loc = loc;
	}

	for ( n=0, mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			loc = mark1->loc - mark2->loc;
			shift_avg += loc;
			n++;
		}
	}
	
	if ( n ) {
		shift_avg /= n;
		for ( mark2 = rec2->mark; mark2; mark2 = mark2->next ) mark2->loc += shift_avg;
	}
	
	if ( verbose )
		cout << "Marker\tResidual" << endl;
	for ( n=0, mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			loc = mark1->loc - mark2->loc;
			d = loc.length();
			R += d*d;
			n++;
			if ( verbose )
				cout << mark1->id << tab << d << endl;
		}
	}
	
	if ( n ) R = sqrt(R/n);
	
	if ( verbose ) {
		cout << "Markers selected:               " << n << endl;
		cout << "Additional shift:               " << shift_avg << endl;
		cout << "Residual:                       " << R << endl << endl;
	}
	
	Bmicrograph*		mg;
	Matrix3				tmat = mat.transpose();
	
	t.trans += shift_avg;
	cout << "Reconstruction origin = " << rec2->origin << endl;
	
	// Adjust the micrograph matrices and origins
	for ( mg = project->field->next->mg; mg; mg = mg->next ) {
		mg->matrix = mat * mg->matrix;
		tmat = mg->matrix.transpose();
		mg->origin -= tmat * t.trans;
//		cout << mg->id << tab << mg->origin << endl;
		mg->origin[2] = 0;
	}
	
	return n;
}


/**
@author Jessica Mavadia, Bernard Heymann
@brief 	Lists two sets of Z coordinates for comparison.
@param 	*project	project with two tilt series.
@return double		root-mean-square-difference.

	The root-mean-square-difference between the z coordinates is calculated.

**/
double		project_dual_zcompare (Bproject* project)
{
	int					n(0);
	double				d, R(0);
	Breconstruction*	rec1 = project->rec;
	Breconstruction*	rec2 = project->rec->next;
	Bmarker*			mark1 = NULL;
	Bmarker*			mark2 = NULL;
				
	cout << "Z Coordinates Compare" << endl << endl;
	cout << "Marker\tz1\tz2\tz1-z2" << endl;
	for ( mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			d = mark1->loc[2] - mark2->loc[2];
			R += d*d;
			n++;
			cout << mark1->id << tab << mark1->loc[2] << tab << mark2->loc[2] << tab << d << endl;
		}
	}
	
	if ( n ) R = sqrt(R/n);
	
	cout << "Z RMSD:                         " << R << " (" << n << ")" << endl << endl;
		
	cout << "Origin1:                        " << rec1->origin[0] << " " << rec1->origin[1] << endl;
	cout << "Origin2:                        " << rec2->origin[0] << " " << rec2->origin[1] << endl;
	cout << "Marker\tx1\ty1\tz1\tMarker\tx2\ty2\tz2" << endl;
	for ( mark1 = rec1->mark; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = rec2->mark; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			cout << mark1->id << tab << mark1->loc[0] << tab << mark1->loc[1] << tab << mark1->loc[2] << tab << 
				mark2->id << tab << mark2->loc[0] << tab << mark2->loc[1] << tab << mark2->loc[2] << endl;
		}
	}
	
	cout << endl;
	
	return R;
}

