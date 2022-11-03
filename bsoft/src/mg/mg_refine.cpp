/**
@file	mg_refine.cpp
@brief	Reciprocal space refinement of orientation parameters of particle images.
@author Bernard Heymann
@date	Created: 20070115
@date	Modified: 20220813
**/

#include "mg_processing.h"
#include "mg_refine.h"
#include "rwmg.h"
#include "mg_ctf.h"
#include "rwimg.h"
#include "Complex.h"
#include "Matrix3.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Local function prototypes
long		part_refine_orientation(Bparticle* part, Bstring& partfile,
				Bimage* pref, Bimage* pmask, Bsymmetry sym, int max_iter, double alpha_step, double accuracy,
				double shift_step, double shift_accuracy, double hi_res, double lo_res,
				int fom_type, vector<double> weight, double edge_radius,
				CTFparam* em_ctf, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags,
				FSI_Kernel* kernel, fft_plan planf);
double		img_refine_monte(Bimage* p, Bimage* pref, double hi_res, double lo_res,  int flags,
				int max_iter, int fom_type, vector<double>& weight, CTFparam* em_ctf, 
				double def_std, double shift_std, double view_std,
				double max_angle, double max_mag, FSI_Kernel* kernel, long &number);
double		img_refine_grid(Bimage* p, Bimage* pref, double hi_res, double lo_res, int flags,
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, double max_mag,
				int fom_type, vector<double>& weight, CTFparam* em_ctf, double def_std, FSI_Kernel* kernel, long &number);

/**
@brief 	Refine the orientation and origin with respect to a reference map.
@param 	*project		image processing parameter structure.
@param 	&reffile		file containing reference map.
@param 	&maskfile		file containing a real space mask (can be empty).
@param 	&sym_string		point group symmetry designator.
@param 	part_select		particle selection for processing.
@param 	max_iter		Monte Carlo maximum number of refining iterations.
@param 	alpha_step		grid search angular step size.
@param 	accuracy		grid search accuracy.
@param 	shift_step		grid shift size.
@param 	shift_accuracy	grid shift accuracy.
@param 	fom_type		type of resolution measure: 0=FRC, 1=DPR
@param 	weight			1D reciprocal space weight curve.
@param 	hi_res			high resolution limit (angstrom).
@param 	lo_res			low resolution limit (angstrom).
@param 	kernel_width	interpolation kernel width.
@param 	kernel_power	interpolation kernel power.
@param 	edge_radius		edge radius to mask background using previous origin.
@param 	def_std			random defocus standard deviation
@param 	shift_std		random origin shift standard deviation.
@param 	view_std		random view shift standard deviation.
@param 	max_angle		maximum random rotation angle adjustment.
@param 	max_mag			maximum magnification adjustment.
@param 	flags			option flags.
@return long				number of comparisons, <0 on error.

	The orientation, origin, magnitude and defocus are refined for each particle.
	The default method uses a grid search around the existing view and origin.
	Specifying the maximum number of iterations switches the algorithm to a
	Monte Carlo search for the best parameters.
	The FOM is either based on the FSC or the DPR.

**/
long		mg_refine_orientations(Bproject* project, Bstring& reffile, Bstring& maskfile,
				Bstring& sym_string, int part_select, int max_iter, 
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, int fom_type, vector<double> weight,
				double hi_res, double lo_res, int kernel_width, int kernel_power,
				double edge_radius, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags)
{
	if ( reffile.length() < 1 ) {
		cerr << "Error: A reference map must be provided!" << endl << endl;
		return -1;
	}

	if ( max_iter < 1 ) {
		if ( alpha_step < SMALLFLOAT || accuracy < SMALLFLOAT ) {
			if ( alpha_step < SMALLFLOAT )
				cerr << "Error: The angular step size is too small! (" << alpha_step*180.0/M_PI << ")" << endl;
			if ( accuracy < SMALLFLOAT )
				cerr << "Error: The accuracy is too small! (" << accuracy*180.0/M_PI << ")" << endl;
			return -1;
		}
		if ( alpha_step < accuracy ) {
			cerr << "Error: The angular step size is smaller than the accuracy! (" << alpha_step*180.0/M_PI << " < " << accuracy*180.0/M_PI << ")" << endl;
			return -1;
		}
	}
	
	random_seed();
	
	if ( sym_string.length() < 1 ) sym_string = "C1";

	Bsymmetry		sym(sym_string);

	long			npart(0), npix, number(0);
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	Bimage* 		pref = read_img(reffile, 1, -1);
	
	if ( !pref ) {
		cerr << "Error: The reference was not read!" << endl << endl;
		return -1;
	}

	pref->check_resolution(hi_res);
	
	if ( !weight.size() )
		weight = C_curve(pref->sizeX(), 1/pref->real_size()[0]);

	// Fourier transform size for the original map size
	long 			ft_size = pref->sizeX();

	if ( pref->compound_type() == TSimple ) {
		if ( verbose )
			cout << "Fourier transforming reference map " << pref->file_name() << endl;
		pref->fft();
	}
	
	// Reduce the reference map's transform size to save memory
	int				diam = (int) (2.2*pref->real_size()[0]/hi_res + 0.5);
	if ( diam%2 ) diam += 1;
	Vector3<long>	size(diam, diam, diam);
	Vector3<double>	scale(1,1,1);
	if ( diam < pref->sizeX() )
		scale = pref->change_transform_size(size);

	pref->set(0, 0);
	pref->set(1, 0);

	npix = pref->sizeX()*pref->sizeY();

	// Mask to apply to particles
	Bimage*		pmask = NULL;
	if ( maskfile.length() ) pmask = read_img(maskfile, 1, 0);

	pref->check_resolution(hi_res);
	
	// Interpolation kernel
	FSI_Kernel*		kernel = new FSI_Kernel(kernel_width, kernel_power);

	// FFT plan for the particle images
	fft_plan		planf = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	
	if ( verbose ) {
		cout << "Refining orientations:" << endl;
		if ( max_iter ) {
			cout << "Maximum iterations:             " << max_iter << endl;
			cout << "Defocus standard deviation:     " << def_std << " A" << endl;
			cout << "Shift standard deviation:       " << shift_std << " pixels" << endl;
			cout << "View standard deviation:        " << view_std << endl;
			cout << "Maximum angle variation:        " << max_angle*180.0/M_PI << " degrees" << endl;
			cout << "Maximum magnification variation: " << max_mag << endl;
		} else {
			cout << "Grid angular step and accuracy: " << alpha_step*180.0/M_PI 
					<< " " << accuracy*180.0/M_PI << " degrees" << endl;
			cout << "Grid shift step and accuracy:   " << shift_step
					<< " " << shift_accuracy << " pixels" << endl;
			cout << "Maximum magnification change:   " << max_mag << endl;
		}
		cout << "Reference:                      " << pref->file_name() << endl;
//		cout << "Reference origin:               " << pref->image->origin()/scale << endl;
		cout << "Reference origin:               " << pref->image->origin()/scale << endl;
    	if ( maskfile.length() ) cout << "Real space mask file:           " << maskfile << endl;
		cout << "Symmetry:                       " << sym_string << endl;
    	cout << "Resolution measurement type:    ";
		if ( fom_type < 1 ) cout << "FRC" << endl;
		else cout << "DPR" << endl;
    	cout << "Resolution limits:              " << hi_res << " - ";
		if ( lo_res > 0 ) cout << lo_res << " A" << endl;
		else cout << " inf A" << endl;
		cout << "Micrograph\tPID\tMag\tOriX\tOriY\tViewX\tViewY\tViewZ\tAngle\tCC\tCV\tTime" << endl;
	}

	pref->phase_shift_to_origin();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_refine_orientations: origin after phase shift: " << pref->image->origin() << endl;

	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	time_t			t = time(NULL);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->part ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) ) {
					number += part_refine_orientation(part, mg->fpart, 	pref, pmask, sym,
							max_iter, alpha_step, accuracy, shift_step, shift_accuracy, 
								hi_res, lo_res, fom_type, weight, edge_radius,
							mg->ctf, def_std, shift_std, view_std, max_angle, max_mag, 
								flags, kernel, planf);
					npart++;
					if ( verbose & VERB_RESULT )
						cout << setw(15) << mg->id << tab << part->id << tab << part->mag << tab << 
								setprecision(2) << part->ori[0] << tab << part->ori[1] << tab << 
								setprecision(4) << part->view[0] << tab << 
								part->view[1] << tab << part->view[2] << tab << 
								setprecision(2) << part->view.angle()*180.0/M_PI << tab << 
								setprecision(4) << part->fom[0] << tab << part->fom[1] << tab << 
								(long)(time(NULL) - t) << endl;
				} else {
					part->fom[0] = part->fom[1] = part->sel = 0;
				}
			}
		}
	}
	
	t = (time(NULL) - t);
	
	long		ts = (long) t;
	long		hr = ts/3600, min = (ts - 3600*hr)/60, sec = ts - 3600*hr - 60*min;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Number of comparisons:          " << number << endl;
		cout << "Comparisons per particle:       " << number*1.0/npart << endl;
		cout << "Time:                           " << ts << " s (" << hr << ":" << min << ":" << sec << ")" << endl;
		cout << "Time per particle:              " << ts*1.0/npart << " s/particle" << endl;
		cout << "Algorithm time:                 " << ts*1.0e6/(number*npix) << " us/pixel" << endl << endl;
	}
	
	delete kernel;

    fft_destroy_plan(planf);
	
	return number;
}

long		project_refine_orientations(Bproject* project, Bstring& reffile, Bstring& maskfile,
				Bstring& sym_string, int part_select, int max_iter, 
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, int fom_type, vector<double> weight,
				double hi_res, double lo_res, int kernel_width, int kernel_power,
				double edge_radius, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags)
{
	if ( reffile.length() < 1 ) {
		cerr << "Error: A reference map must be provided!" << endl << endl;
		return -1;
	}
	
	if ( max_iter < 1 ) {
		if ( alpha_step < SMALLFLOAT || accuracy < SMALLFLOAT ) {
			if ( alpha_step < SMALLFLOAT )
				cerr << "Error: The angular step size is too small! (" << alpha_step*180.0/M_PI << ")" << endl;
			if ( accuracy < SMALLFLOAT )
				cerr << "Error: The accuracy is too small! (" << accuracy*180.0/M_PI << ")" << endl;
			return -1;
		}
		if ( alpha_step < accuracy ) {
			cerr << "Error: The angular step size is smaller than the accuracy! (" << alpha_step*180.0/M_PI << " < " << accuracy*180.0/M_PI << ")" << endl;
			return -1;
		}
	}
	
	random_seed();
	
	if ( sym_string.length() < 1 ) sym_string = "C1";

	Bsymmetry		sym(sym_string);

	long			npix;
	Bimage* 		pref = read_img(reffile, 1, -1);
	
	if ( !pref ) {
		cerr << "Error: The reference was not read!" << endl << endl;
		return -1;
	}

	pref->check_resolution(hi_res);
//	pref->information();

	// Fourier transform size for the original map size
	long 			ft_size = pref->sizeX();

	if ( pref->compound_type() == TSimple ) {
		if ( verbose )
			cout << "Fourier transforming reference map " << pref->file_name() << endl;
		pref->fft();
	}
	
	// Reduce the reference map's transform size to save memory
	int				diam = (int) (2.2*pref->real_size()[0]/hi_res + 0.5);
	if ( diam%2 ) diam += 1;
	Vector3<long>	size(diam, diam, diam);
	Vector3<double>	scale(1,1,1);
	if ( diam < pref->sizeX() )
		scale = pref->change_transform_size(size);

	pref->set(0, 0);
	pref->set(1, 0);

	npix = pref->sizeX()*pref->sizeY();

	// Mask to apply to particles
	Bimage*		pmask = NULL;
	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		pmask->origin(pmask->size()/2);
	}

	pref->check_resolution(hi_res);
	
	// Interpolation kernel
	FSI_Kernel*		kernel = new FSI_Kernel(kernel_width, kernel_power);

	// FFT plan for the particle images
	fft_plan		planf = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	
	long			npart(0);
	Bparticle**		partarr = project_mg_particle_array(project, part_select, npart);
	
	if ( verbose ) {
		cout << "Refining orientations:" << endl;
		if ( max_iter ) {
			cout << "Maximum iterations:             " << max_iter << endl;
			cout << "Shift standard deviation:       " << shift_std << " pixels" << endl;
			cout << "View standard deviation:        " << view_std << endl;
			cout << "Maximum angle variation:        " << max_angle*180.0/M_PI << " degrees" << endl;
			cout << "Maximum magnification variation: " << max_mag << endl;
		} else {
			cout << "Grid step and accuracy:         " << alpha_step*180.0/M_PI 
					<< " " << accuracy*180.0/M_PI << " degrees" << endl;
			cout << "Maximum magnification change:   " << max_mag << endl;
		}
		cout << "Defocus standard deviation:     " << def_std << " A" << endl;
		cout << "Reference:                      " << pref->file_name() << endl;
		cout << "Reference origin:               " << pref->image->origin()/scale << endl;
    	if ( maskfile.length() ) cout << "Real space mask file:           " << maskfile << endl;
		cout << "Symmetry:                       " << sym_string << endl;
    	cout << "Resolution measurement type:    ";
		if ( fom_type < 1 ) cout << "FRC" << endl;
		else cout << "DPR" << endl;
    	cout << "Resolution limits:              " << hi_res << " - ";
		if ( lo_res > 0 ) cout << lo_res << " A" << endl;
		else cout << " inf A" << endl;
   		cout << "Weighting scheme:               ";
		if ( !weight.size() )
			cout << "C curve" << endl;
		else
			cout << "Input curve" << endl;
		cout << "Micrograph\tPID\tDef\tSamX\tSamY\tOriX\tOriY\tViewX\tViewY\tViewZ\tAngle\tCC\tCV\tTime" << endl;
	}

	pref->phase_shift_to_origin();
	
	if ( !weight.size() )
		weight = C_curve(pref->sizeX(), 1/pref->real_size()[0]);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_refine_orientations: origin after phase shift: " << pref->image->origin() << endl;

	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	time_t			t = time(NULL);
	long			number(0);
	int*			nar = new int[npart];
	
#ifdef HAVE_GCD
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(npart, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bparticle*		part = partarr[i];
		Bmicrograph*	mg = part->mg;
		nar[i] = part_refine_orientation(part,
					mg->fpart, pref, pmask, sym,
					max_iter, alpha_step, accuracy, shift_step, shift_accuracy,
					hi_res, lo_res, fom_type, weight, edge_radius,
					mg->ctf, def_std, shift_std, view_std, max_angle, max_mag, 
					flags, kernel, planf);
		dispatch_sync(myq, ^{
			if ( verbose & VERB_RESULT )
				cout << setw(15) << mg->id << tab << part->id << tab <<
					setprecision(0) << part->def << tab << 
					setprecision(2) << part->pixel_size[0] << tab << part->pixel_size[1] << tab << 
					setprecision(2) << part->ori[0] << tab << part->ori[1] << tab << 
					setprecision(4) << part->view[0] << tab << 
					part->view[1] << tab << part->view[2] << tab << 
					setprecision(2) << part->view.angle()*180.0/M_PI << tab << 
					setprecision(4) << part->fom[0] << tab << part->fom[1] << tab << 
					(long)(time(NULL) - t) << endl;
		});
	});
#else
#pragma omp parallel for
	for ( long i=0; i<npart; i++ ) {
		Bparticle*		part = partarr[i];
		Bmicrograph*	mg = part->mg;
		nar[i] = part_refine_orientation(part,
					mg->fpart, pref, pmask, sym,
					max_iter, alpha_step, accuracy, shift_step, shift_accuracy,
					hi_res, lo_res, fom_type, weight, edge_radius,
					mg->ctf, def_std, shift_std, view_std, max_angle, max_mag, 
					flags, kernel, planf);
	#pragma omp critical
		{
			if ( verbose & VERB_RESULT )
				cout << setw(15) << mg->id << tab << part->id << tab << 
					setprecision(0) << part->def << tab << 
					setprecision(2) << part->pixel_size[0] << tab << part->pixel_size[1] << tab << 
					setprecision(2) << part->ori[0] << tab << part->ori[1] << tab << 
					setprecision(4) << part->view[0] << tab << 
					part->view[1] << tab << part->view[2] << tab << 
					setprecision(2) << part->view.angle()*180.0/M_PI << tab << 
					setprecision(4) << part->fom[0] << tab << part->fom[1] << tab << 
					(long)(time(NULL) - t) << endl;
		}
	}
#endif

	for ( long i=0; i<npart; i++ ) number += nar[i];

	t = (time(NULL) - t);
	
	long		ts = (long) t;
	long		hr = ts/3600, min = (ts - 3600*hr)/60, sec = ts - 3600*hr - 60*min;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Number of comparisons:          " << number << endl;
		cout << "Comparisons per particle:       " << number*1.0/npart << endl;
		cout << "Time:                           " << ts << " s (" << hr << ":" << min << ":" << sec << ")" << endl;
		cout << "Time per particle:              " << ts*1.0/npart << " s/particle" << endl;
		cout << "Algorithm time:                 " << ts*1.0e6/(number*npix) << " us/pixel" << endl << endl;
	}
	
	delete[] partarr;
	delete[] nar;
	
	delete kernel;

    fft_destroy_plan(planf);
	
	return number;
}

/**
@brief 	Refine the orientation and origin of one particle with respect to a reference map.
@param 	*part			particle.
@param 	&partfile		particle filename.
@param 	*pref			reference map.
@param 	*pmask			mask.
@param 	*sym			symmetry structure.
@param 	max_iter		maximum number of refining iterations.
@param 	alpha_step		grid search angular step size.
@param 	accuracy		grid search accuracy.
@param 	shift_step		grid shift size.
@param 	shift_accuracy	grid shift accuracy.
@param 	hi_res			high resolution limit (angstrom).
@param 	lo_res			low resolution limit (angstrom).
@param 	fom_type		type of resolution measure: 0=FRC, 1=DPR
@param 	weight			1D reciprocal space weight curve.
@param 	edge_radius		edge radius to mask background using previous origin.
@param 	*em_ctf			CTF parameters (if NULL, not used).
@param 	def_std			random defocus standard deviation
@param 	shift_std		random origin shift standard deviation.
@param 	view_std		random view shift standard deviation.
@param 	max_angle		maximum random rotation angle adjustment.
@param 	max_mag			maximum magnification adjustment.
@param 	flags			option flags.
@param 	*kernel			interpolation kernel.
@param 	planf			FFT forward plan for 2D images.
@return long			number of comparisons.

	The orientation and origin are iteratively modified in small random steps,
	with selection based on the Fourier shell correlation.

**/
long		part_refine_orientation(Bparticle* part, Bstring& partfile,
				Bimage* pref, Bimage* pmask, Bsymmetry sym, int max_iter, double alpha_step, double accuracy,
				double shift_step, double shift_accuracy, double hi_res, double lo_res,
				int fom_type, vector<double> weight, double edge_radius,
				CTFparam* em_ctf, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags,
				FSI_Kernel* kernel, fft_plan planf)
{
	FOMType 	fom_tag[NFOM] = {FOM, FOM_CV};

	if ( flags & CHECK_PPX )
		if ( ppx_check(part, fom_tag) ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_refine_orientation: file=" << partfile << " image=" << part->id << endl;

	Bimage*		p = read_img(partfile, 1, part->id - 1);
	if ( !p )
		return error_show("Error in part_refine_orientation", __FILE__, __LINE__);
					
	if ( part->pixel_size[0] > 0 ) p->sampling(part->pixel_size);
	else if ( part->mg->pixel_size[0] > 0 ) p->sampling(part->mg->pixel_size);
	part->pixel_size = p->image->sampling();

	if ( part->ori[0] <= 0 ) {
		if ( p->image->origin()[0] > 0 ) part->ori[0] = p->image->origin()[0];
		else part->ori[0] = p->sizeX()/2;
	}
	if ( part->ori[1] <= 0 ) {
		if ( p->image->origin()[1] > 0 ) part->ori[1] = p->image->origin()[1];
		else part->ori[1] = p->sizeY()/2;
	}
	part->ori[2] = 0;
	p->origin(part->ori);
	
	if ( part->mag < 0.5 ) {
		if ( p->image->magnification() > 0.5 )
			part->mag = p->image->magnification();
		else part->mag = 1;
	}

	part->fom[0] = part->fom[1] = 0;
	
	p->change_type(Float);
	p->rescale_to_avg_std(0, 1);
	p->calculate_background();
	
	int				edge_size = (int) (2*edge_radius);
	Vector3<long>	size(edge_size, edge_size, 1);
	Vector3<double>	start(part->ori[0]-edge_radius, part->ori[1]-edge_radius, 0);
	if ( edge_size )
		p->edge(1, size, start, 1, FILL_BACKGROUND, 0);

	Vector3<double>	shift;
	Bimage*			pm;
	
	if ( pmask ) {
		pm = pmask->copy();
		shift = p->image->origin() - pm->image->origin();
		pm->shift(shift);
		p->multiply(pm);
		delete pm;
	}
	
	p->fft(planf, 0);

	p->set(0, 0);
	p->set(1, 0);
	
	// Reduce the particle transform size to match that of the reference
	Vector3<double>	scale(1,1,1);
	if ( pref->sizeX() < p->sizeX() ) {
		size = Vector3<int>(pref->sizeX(), pref->sizeY(), 1);
		scale = p->change_transform_size(size);
	}
	
	int				apply_ctf(0);
	double			bestdef(0);
	if ( em_ctf ) {
		bestdef = em_ctf->defocus_average();
		if ( def_std ) {
			if ( part->def < 1 ) part->def = em_ctf->defocus_average();
			bestdef = part->def;
			apply_ctf = 1;
			flags |= APPLY_CTF;
		}
	}

	long			i, number(0);
	double			fom, bestmag = part->mag;
	Vector3<double>	sam(p->image->sampling()), bestsam(sam);
	Vector3<double>	ori(p->image->origin()), bestorigin(ori);
	View			bestview = part->view;
	View*			view = symmetry_get_all_views(sym, part->view);
	View*			v;
	CTFparam*		ctf_copy = new CTFparam;
	
	for ( i=1, v = view; v; v = v->next, i++ ) {
		p->sampling(sam);
		p->origin(ori);
		p->view(*v);
		p->image->magnification(part->mag);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG part_refine_orientation: view=" << *v << endl;
	
		if ( em_ctf ) {
			ctf_copy->update(em_ctf);
			if ( apply_ctf && part->def > 1 ) ctf_copy->defocus_average(part->def);
		}
		
		if ( max_iter )
			fom = img_refine_monte(p, pref, hi_res, lo_res, flags, max_iter, fom_type,
				weight, ctf_copy, def_std,
				shift_std, view_std, max_angle, max_mag, kernel, number);
		else 
			fom = img_refine_grid(p, pref, hi_res, lo_res, flags, alpha_step, accuracy,
				shift_step, shift_accuracy, max_mag,
				fom_type, weight, ctf_copy, def_std, kernel, number);
		
		if ( part->fom[0] < fom ) {
			part->fom[0] = fom;
			part->fom[1] = p->image->FOM();
			part->sel = i;
			if ( ctf_copy && apply_ctf ) bestdef = ctf_copy->defocus_average();
			bestsam = p->image->sampling();
			bestorigin = p->image->origin();
			bestview = p->image->view();
			bestmag = p->image->magnification();
			if ( bestdef < 1000 ) {
				cerr << "Bad defocus! pid = " << part->id << tab << bestdef << endl;
			}
		}
		
	}

	if ( ctf_copy ) delete ctf_copy;

	if ( apply_ctf ) {
		if ( bestdef > 1 ) part->def = bestdef;
		part->dev = em_ctf->defocus_deviation();
		part->ast = em_ctf->astigmatism_angle();
	}
	part->pixel_size = bestsam*scale;
	part->ori = bestorigin/scale;
	part->view = bestview;
	part->mag = bestmag;
	
	kill_list((char *) view, sizeof(View));
	delete p;

	if ( flags & WRITE_PPX ) {
		Bstring		ppx_name = ppx_filename(part->mg->id, part->id);
		write_particle(ppx_name, part, 0, 0, fom_tag);
		ppx_name = 0;
	}
	
	return number;
}

/*
	FRC:
		     sum((FoFc*)re)
		-------------------------
		sqrt(sum(Fc^2)*sum(Fo^2))
	DPR:
		    (sum((|Fo|+|Fc|)*dPhi^2))
		sqrt|-----------------------|
		    (    sum(|Fo|+|Fc|)     )
*/
double		img_recip_space_fom(Bimage* p, Bimage* pref, double hi_res, double lo_res,
				int fom_type, vector<double>& weight, int flags, CTFparam* em_ctf)
{
	bool			apply_ctf(flags & APPLY_CTF);
	bool			invert(flags & INVERT);
	
	if ( pref->sampling(0)[0] < 0.01 ) pref->sampling(1,1,1);
	if ( lo_res > 0 && hi_res > lo_res ) swap(hi_res, lo_res);
	pref->check_resolution(hi_res);
	
	double			max_rad = pref->real_size()[0]/hi_res;
	double			max_rad_cv = 1.1*max_rad;
	if ( max_rad_cv <= max_rad ) max_rad_cv = max_rad + 1;
	
	double			max_rad_sq = max_rad*max_rad;
	double			max_rad_cv_sq = max_rad_cv*max_rad_cv;
	double			min_rad_sq(0);
	if ( lo_res > 0 ) {
		min_rad_sq = pref->real_size()[0]/lo_res;
		min_rad_sq *= min_rad_sq;
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_recip_space_fom: max_rad=" << max_rad << " (" << max_rad_sq << ")" << endl;
	
	double			dphi, ctf_fac;
	if ( !em_ctf ) apply_ctf = 0;
	else if ( em_ctf->defocus_average() < 1 ) apply_ctf = 0;
	
	long 			i, j, xx, yy;
	long			hx = (p->sizeX() - 1)/2, hy = (p->sizeY() - 1)/2;	
	double 			sx2, sy2, s2, d2, amp, denom, f, cv, wsum;
	Vector3<double>	m, d, iv;
	Vector3<double>	freq_scale(1.0/p->real_size()[0], 1.0/p->real_size()[1], 1.0);
//	cout << "size = " << p->size() << endl;
//	cout << "sampling = " << p->image->sampling() << endl;
//	cout << "real size = " << p->real_size() << endl;
//	bexit(-1);

	double*			w = new double[p->sizeX()];
	for ( j=0; j<p->sizeX(); j++ ) w[j] = 1;
	if ( weight.size() )
		for ( size_t j=0; j<weight.size(); j++ ) if ( weight[j] ) w[j] = weight[j];
		
	double*			fom = new double[p->sizeX()];
	double*			dsum = new double[p->sizeX()];
	double*			msum = new double[p->sizeX()];
	for ( i=0; i<p->sizeX(); i++ ) fom[i] = dsum[i] = msum[i] = 0;

	Complex<float>*	data = (Complex<float> *) p->data_pointer();
	Complex<float>*	refdata = (Complex<float> *) pref->data_pointer();
	Complex<float>	dtemp;
	
	for ( yy=i=0; yy<p->sizeY(); yy++ ) {
		iv[1] = yy;
		if ( iv[1] > hy ) iv[1] -= p->sizeY();
		sy2 = freq_scale[1]*iv[1];
		sy2 *= sy2;
		for ( xx=0; xx<p->sizeX(); xx++, i++ ) {
			iv[0] = xx;
			if ( iv[0] > hx ) iv[0] -= p->sizeX();
			sx2 = freq_scale[0]*iv[0];
			sx2 *= sx2;
			d2 = iv.length2();
			if ( d2 >= min_rad_sq && d2 <= max_rad_cv_sq ) {
				s2 = sx2 + sy2;
				j = (long) sqrt(d2);
				dtemp = refdata[i];
				if ( apply_ctf ) {
					ctf_fac = em_ctf->calculate(s2, atan2(iv[1],iv[0]));
					if ( invert ) ctf_fac = -ctf_fac;
					dtemp *= ctf_fac;
				}
				if ( fom_type < 1 ) {
					fom[j] += data[i].real()*dtemp.real() + data[i].imag()*dtemp.imag();
					msum[j] += dtemp.power();
					dsum[j] += data[i].power();
				} else {
					amp = data[i].amp() + dtemp.amp();
					dphi = angle_set_negPI_to_PI(data[i].phi() - dtemp.phi());
					fom[j] += amp*dphi*dphi;
					msum[j] += amp;
				}
			}
		}
	}

	if ( fom_type < 1 ) {
		for ( f=wsum=0, j=0; j<max_rad; j++ ) {
			wsum += w[j];
			denom = dsum[j]*msum[j];
			if ( denom > 0 ) f += w[j]*fom[j]/sqrt(denom);
		}
		f /= wsum;
		for ( cv=wsum=0, j=(long)max_rad; j<max_rad_cv; j++ ) {
			wsum += w[j];
			denom = dsum[j]*msum[j];
			if ( denom > 0 ) cv += w[j]*fom[j]/sqrt(denom);
		}
		cv /= wsum;
	} else {
		for ( f=wsum=0, j=0; j<max_rad; j++ ) {
			wsum += w[j];
			if ( msum[j] ) f += w[j]*cos(sqrt(fom[j]/msum[j]));
		}
		f /= wsum;
		for ( cv=wsum=0, j=(long)max_rad; j<max_rad_cv; j++ ) {
			wsum += w[j];
			if ( msum[j] ) cv += w[j]*cos(sqrt(fom[j]/msum[j]));
		}
		cv /= wsum;
	}
	
	delete[] fom;
	delete[] dsum;
	delete[] msum;	
	delete[] w;
	
	p->image->FOM(cv);
	
	return f;
}

double		img_fit_in_recip_space_monte(Bimage* p, Bimage* pref, double hi_res, double lo_res, int flags,
				int max_iter, int fom_type, vector<double>& weight, CTFparam* em_ctf, double def_std, double shift_std, 
				FSI_Kernel* kernel, long &number)
{
	// Defocus may be changed in this function
	bool			apply_ctf(flags&APPLY_CTF);
	double			predef(0), bestdef(0);
	if ( def_std && em_ctf && em_ctf->defocus_average() ) {
		predef = bestdef = em_ctf->defocus_average();
		apply_ctf = 1;
	}
	
	Vector3<double>	scale = pref->sampling(0)/p->sampling(0);
	View			view(p->image->view());
	Matrix3			mat = view.matrix();
//	mat = mat * p->image->magnification();
	mat[0] *= scale;
	mat[1] *= scale;
//	cout << mat << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_fit_in_recip_space_monte: view=" << view << endl;
	
	Bimage*			pcs = pref->central_section(mat, hi_res/1.11, kernel);
	
	long			iter;
	double 			fom, prefom(-1), bestfom(-1), bestcv(0);
	double			irm = 1.0/get_rand_max();
	Vector3<double>	origin(p->image->origin());
	Vector3<double>	r, preorigin = origin, bestorigin = origin;

	Bimage*			pcopy;
	
	for ( iter = 0; iter<max_iter; iter++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_fit_in_recip_space_monte: iter=" << iter << endl;
		if ( iter > 0 ) {
			if ( apply_ctf ) em_ctf->defocus_average(em_ctf->defocus_average() + random_gaussian(0, def_std));
			r = vector3_xy_random_gaussian(0, shift_std);
			origin += r;
		}	
			
		pcopy = p->copy();
		pcopy->origin(origin);
		pcopy->phase_shift_to_origin();
				
		fom = img_recip_space_fom(pcopy, pcs, hi_res, lo_res, fom_type, weight, flags, em_ctf);
		number++;

		if ( bestfom < fom ) {
			bestfom = fom;
			if ( apply_ctf ) bestdef = em_ctf->defocus_average();
			bestorigin = origin;
			bestcv = pcopy->image->FOM();
			if ( verbose & VERB_FULL )
				cout << iter << tab << fom << tab << bestcv << endl;
		}
		
		if ( fom/prefom > random()*irm ) {
			if ( apply_ctf ) predef = em_ctf->defocus_average();
			preorigin = origin;
			prefom = fom;
		} else {
			if ( apply_ctf ) em_ctf->defocus_average(predef);
			origin = preorigin;
			fom = prefom;
		}
		delete pcopy;
	}

	p->origin(bestorigin);
	p->image->FOM(bestcv);
	if ( apply_ctf ) em_ctf->defocus_average(bestdef);

	delete pcs;
	
	return bestfom;
}

double		img_fit_in_recip_space_grid(Bimage* p, Bimage* pref, double hi_res, double lo_res, int flags,
				double step, double accuracy, int fom_type, vector<double>& weight, CTFparam* em_ctf, double def_std, FSI_Kernel* kernel, long &number)
{
	// Defocus may be changed in this function
	bool			apply_ctf(flags&APPLY_CTF);
	double			bestdef(em_ctf->defocus_average());
	if ( def_std && em_ctf && bestdef ) {
		bestdef = em_ctf->defocus_average();
		apply_ctf = 1;
	}
	
	Vector3<double>	scale = pref->sampling(0)/p->sampling(0);
//	Vector3<double>	scale = p->sampling(0)/pref->sampling(0);
	scale[2] = 1;
	
	View			view(p->image->view());
	Matrix3			mat = view.matrix();
//	mat = mat * p->image->magnification();
//	mat = scale * mat;
	mat[0] *= scale;
	mat[1] *= scale;
//	cout << mat << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_fit_in_recip_space_grid: view=" << view << endl;
	
	Bimage*			pcs = pref->central_section(mat, hi_res/1.11, kernel);
	
	int				i, ibest;
	double			fom, prevfom, bestfom(-1);
	double			def, def_start(bestdef - 5*def_std), def_end(bestdef + 5*def_std), bestcv(0);
	Vector3<double>	ori, oldori(p->image->origin()), bestorigin(oldori);
	Bimage*			pcopy;

	if ( apply_ctf ) {
		pcopy = p->copy();
		pcopy->phase_shift_to_origin();
		for ( def = def_start; def <= def_end; def += def_std ) {
			em_ctf->defocus_average(def);
			fom = img_recip_space_fom(pcopy, pcs, hi_res, lo_res, fom_type, weight, flags, em_ctf);
			number++;
			if ( bestfom < fom ) {
				bestfom = fom;
				bestdef = def;
				bestcv = pcopy->image->FOM();
				if ( verbose & VERB_FULL )
					cout << def << tab << fom << tab << bestcv << endl;
			}			
		}
		delete pcopy;
		em_ctf->defocus_average(bestdef);
	}

	while ( step >= accuracy ) {
		prevfom = bestfom;
		bestfom = -1;
		for ( i=ibest=0, ori[1] = oldori[1] - step; ori[1] <= oldori[1] + step; ori[1] += step ) {
			for ( ori[0] = oldori[0] - step; ori[0] <= oldori[0] + step; ori[0] += step, ++i ) {
				pcopy = p->copy();
				pcopy->image->origin(ori);
				pcopy->phase_shift_to_origin();

				fom = img_recip_space_fom(pcopy, pcs, hi_res, lo_res, fom_type, weight, flags, em_ctf);
	
				number++;
	
				if ( bestfom < fom ) {
					bestfom = fom;
					ibest = i;
					bestorigin = ori;
					bestcv = pcopy->image->FOM();
					if ( verbose & VERB_FULL )
						cout << ori[0] << tab << ori[1] << tab << fom << tab << bestcv << endl;
				}
			
				delete pcopy;
			}
		}
		if ( ( ibest == i/2 ) || fabs(prevfom - bestfom) < 1e-6 )
			step /= 2;	// Contract around best case
	}

	p->origin(bestorigin);
	p->image->FOM(bestcv);

	delete pcs;
	
	return bestfom;
}

/**
@brief 	Refine the orientation and origin of one particle with respect to a reference map.
@param 	*p				particle image.
@param 	*pref			reference map.
@param 	hi_res			high resolution limit (angstrom).
@param 	lo_res			low resolution limit (angstrom).
@param 	flags			APPLY_CTF, INVERT.
@param 	max_iter		maximum number of refining iterations.
@param 	fom_type		type of resolution measure: 0=FRC, 1=DPR
@param 	weight			1D reciprocal space weight curve.
@param 	*em_ctf			CTF parameters (if NULL, not used).
@param 	def_std			random defocus standard deviation
@param 	shift_std		random origin shift standard deviation.
@param 	view_std		random view shift standard deviation.
@param 	max_angle		maximum random rotation angle adjustment.
@param 	max_mag			maximum magnification adjustment.
@param 	*kernel			interpolation kernel.
@param 	&number			number of comparisons.
@return double			best figure-of-merit.

	The orientation and origin are iteratively modified in small random steps,
	with selection based on the Fourier shell correlation.
	The best result is returned in the orientation parameters of the image structure.
	The FOM in the image structure is the first FOM determined, while the best
	FOM is returned.

**/
double		img_refine_monte(Bimage* p, Bimage* pref, double hi_res, double lo_res, int flags,
				int max_iter, int fom_type, vector<double>& weight, CTFparam* em_ctf, 
				double def_std, double shift_std, double view_std,
				double max_angle, double max_mag, FSI_Kernel* kernel, long &number)
{
	int				apply_ctf(0);
	double			def(0), predef(0), bestdef(0);
	if ( def_std && em_ctf && em_ctf->defocus_average() ) {
		def = predef = bestdef = em_ctf->defocus_average();
		apply_ctf = (flags & APPLY_CTF);
	}
	
	long			iter;
	double 			fom, prefom(-1), bestfom(-1), bestcv(0);
	double			a, irm(1.0/get_rand_max());
//	double			mag = p->image->magnification();
	double			magx(1), magy(1);
//	double			bestmag(mag), premag(mag);
	Vector3<double>	sam(p->image->sampling()), bestsam(sam), presam(sam);
	Vector3<double>	ori(p->image->origin()), r, bestorigin(ori);
	View			view(p->image->view());
	View			bestview = view, preview = view;

	for ( iter = 0; iter<max_iter; iter++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_refine_monte: iter=" << iter << endl;
		if ( iter > 0 ) {
			r = vector3_random_gaussian(0, view_std);
			a = max_angle*(random()*irm*2 - 1);
			view[0] += r[0];
			view[1] += r[1];
			view[2] += r[2];
			view[3] += a;
			view.normalize();
			if ( max_mag > 0 ) {
				magx = 1 + max_mag*(random()*irm*2 - 1);
				magy = 1 + max_mag*(random()*irm*2 - 1);
			}
		}
		p->view(view);
//		p->image->magnification(mag);
		p->image->sampling(magx * sam[0], magy * sam[1], 1);
		if ( apply_ctf ) em_ctf->defocus_average(def);	// Resets defocus to original
		p->origin(ori);									// Resets origin to original
		fom = img_fit_in_recip_space_monte(p, pref, hi_res, lo_res, flags, max_iter, fom_type, weight, em_ctf,
				def_std, shift_std, kernel, number);
		if ( bestfom < fom ) {
			bestfom = fom;
			bestview = view;
//			bestmag = mag;
			bestsam = p->image->sampling();
			bestorigin = p->image->origin();
			if ( apply_ctf ) bestdef = em_ctf->defocus_average();
			bestcv = p->image->FOM();
			if ( verbose & VERB_FULL )
				cout << iter << tab << fom << endl;
		}
		if ( fom/prefom > random()*irm ) {
			preview = view;
//			premag = mag;
			presam = sam;
			prefom = fom;
		} else {
			view = preview;
//			mag = premag;
			sam = presam;
			fom = prefom;
		}
	}

	p->view(bestview);
	p->image->origin(bestorigin);
//	p->image->magnification(bestmag);
	p->image->sampling(bestsam);
	p->image->FOM(bestcv);
	if ( apply_ctf ) em_ctf->defocus_average(bestdef);

	return bestfom;
}

/**
@brief 	Refine the orientation and origin of one particle with respect to a reference map.
@param 	*p				particle image.
@param 	*pref			reference map.
@param 	hi_res			high resolution limit (angstrom).
@param 	lo_res			low resolution limit (angstrom).
@param 	flags			APPLY_CTF, INVERT.
@param 	alpha_step		grid search angular step size.
@param 	accuracy		grid search accuracy.
@param 	shift_step		grid shift size.
@param 	shift_accuracy	grid shift accuracy.
@param 	max_mag			maximum magnification adjustment.
@param 	fom_type		type of resolution measure: 0=FRC, 1=DPR
@param 	weight			1D reciprocal space weight curve.
@param 	*em_ctf			CTF parameters (if NULL, not used).
@param 	def_std			defocus step size.
@param 	*kernel			interpolation kernel.
@param 	&number			number of comparisons.
@return double			best figure-of-merit.

	The orientation and origin are refined using a grid search
	with selection based on the FSC or DPR.
	The best result is returned in the orientation parameters of the image structure.
	The FOM in the image structure is the first FOM determined, while the best
	FOM is returned.

**/
double		img_refine_grid(Bimage* p, Bimage* pref, double hi_res, double lo_res, int flags,
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, double max_mag, int fom_type, vector<double>& weight, 
				CTFparam* em_ctf, double def_std, FSI_Kernel* kernel, long &number)
{
	if ( alpha_step < SMALLFLOAT || accuracy < SMALLFLOAT ) {
		if ( alpha_step < SMALLFLOAT )
			cerr << "Error: The angular step size is too small! (" << alpha_step*180.0/M_PI << ")" << endl;
		if ( accuracy < SMALLFLOAT )
			cerr << "Error: The accuracy is too small! (" << accuracy*180.0/M_PI << ")" << endl;
		return -1;
	}

	if ( alpha_step < accuracy ) {
		cerr << "Error: The angular step size is smaller than the accuracy! (" << alpha_step*180.0/M_PI << " < " << accuracy*180.0/M_PI << ")" << endl;
		return -1;
	}

	double			def(0), bestdef(0);
	if ( def_std && em_ctf && em_ctf->defocus_average() )
		def = bestdef = em_ctf->defocus_average();
	
	int				i, ibest;
	double 			fom, prefom, bestfom(-1), bestcv(0);
//	double			mag = p->image->magnification(), bestmag = mag;
	double			magx, magy, mmin(1.0-max_mag), mmax(1.0+max_mag), mstep(max_mag/10);
	Vector3<double>	sam(p->image->sampling()), bestsam(sam);
	Vector3<double>	ori(p->image->origin()), bestorigin(ori);
	View			view(p->image->view()), bestview(view);
	View			*views, *v;

	while ( alpha_step >= accuracy ) {
		prefom = bestfom;
		bestfom = -1;
		views = views_for_refinement(bestview, alpha_step);
		for ( i=ibest=0, v = views; v; v = v->next, i++ ) {
			p->view(*v);
			if ( def ) em_ctf->defocus_average(def);		// Resets defocus to original
			p->origin(ori);	// Resets origin to original
//			p->image->magnification(mag);		// Resets magnification to original
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG img_refine_grid: view=" << *v << endl;
			fom = img_fit_in_recip_space_grid(p, pref, hi_res, lo_res, flags, shift_step, shift_accuracy, fom_type, weight, em_ctf, def_std, kernel, number);
			if ( bestfom < fom ) {
				bestfom = fom;
				ibest = i;
				bestview = *v;
				bestorigin = p->image->origin();
				if ( def ) bestdef = em_ctf->defocus_average();
				bestcv = p->image->FOM();
			}
		}
		kill_list((char *) views, sizeof(View));
		if ( verbose & VERB_FULL )
			cout << alpha_step << tab << bestfom << tab << bestcv << endl;
		if ( ibest == i/2 || fabs(prefom - bestfom) < 1e-6 ) alpha_step /= 2;	// Contract around best case
	}

	p->view(bestview);
	ori = bestorigin;
//	if ( def ) {
//		em_ctf->defocus_average(bestdef);
//		def = bestdef;
//	}

	if ( max_mag > 0 ) {
		prefom = bestfom;
		bestfom = -1;
		for ( magy = mmin; magy <= mmax; magy += mstep ) {
			for ( magx = mmin; magx <= mmax; magx += mstep ) {
				p->image->sampling(sam[0] * magx, sam[1] * magy, 1);
//				cerr << p->sampling(0) << endl;
				p->origin(ori);
				if ( def ) em_ctf->defocus_average(def);		// Resets defocus to original
//		for ( mag = 1-max_mag; mag <= 1+max_mag; mag += 0.001 ) {
//			p->image->magnification(mag);
//			if ( verbose & VERB_DEBUG )
//				cout << "DEBUG img_refine_grid: mag=" << mag << endl;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG img_refine_grid: sam=" << p->image->sampling() << endl;
				fom = img_fit_in_recip_space_grid(p, pref, hi_res, lo_res, flags, 0.4, 0.1, fom_type, weight, em_ctf, def_std, kernel, number);
				if ( bestfom < fom ) {
					bestfom = fom;
//					bestmag = mag;
					bestsam = p->image->sampling();
					bestorigin = p->image->origin();
					if ( def ) bestdef = em_ctf->defocus_average();
					bestcv = p->image->FOM();
//					if ( verbose & VERB_FULL )
//						cout << mag << tab << fom << tab << bestcv << endl;
					if ( verbose & VERB_FULL )
						cout << magx << tab << magy << tab << fom << tab << bestcv << endl;
				}
			}
		}
	}

	p->sampling(bestsam);
	p->origin(bestorigin);
//	p->image->magnification(bestmag);
	p->image->FOM(bestcv);
	if ( def ) {	// Make sure the defocus is reasonable
		if ( bestdef > def - 5*def_std && bestdef < def + 5*def_std ) 
			em_ctf->defocus_average(bestdef);
		else
			em_ctf->defocus_average(def);
	}

	return bestfom;
}
