/**
@file	mg_orient.cpp
@brief	Determines orientation angles and x,y origins of single particle images 
@author	Bernard Heymann and David M. Belnap
@date	Created: 20010403
@date	Modified: 20200218 (BH)
**/

#include "mg_orient.h"
#include "mg_processing.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "mg_ctf.h"
#include "symmetry.h"
#include "linked_list.h"
#include "rwimg.h"
#include "ps_views.h"
#include "utilities.h"

#include <fstream>
#include <sys/stat.h>
#include <fcntl.h>


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int 		part_determine_orientation(Bparticle* part, Bimage* proj, Bimage* part_mask,
				Bimage* prs_mask, Bsymmetry& sym, int bin,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags,
				fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D);


double		correlation_coefficient_adjust(double cc, double std)
{
	int			i;
	double		ccc = cc, ccd = 1;
	double		fac = sqrt(2.0)*std;
	
	for ( i=0; i<100 && fabs(ccd) > 1e-10; i++ ) {
		ccd = cc - ccc*ccc*ccc - fac*ccc*(1-ccc);
		ccc += ccd;
	}
	
//	if ( verbose & VERB_DEBUG )
//		cout << "DEBUG correlation_coefficient_adjust: i=" << i << "  ccd=" << ccd << endl;
	
	return ccc;
}

/**
@brief 	Prepare reference images from particles.
@param 	project			micrograph project structure.
@param 	first			first particle image to use.
@param 	number			number of particle images to use.
@param 	bin				data compression by binning.
@param 	ctf_action		flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@return Bimage*			reference images, NULL on error.

	The particle images are shifted to center their origins.
	
**/
Bimage*		project_prepare_2D_references(Bproject* project, long first,
				long number, int bin, int ctf_action, double wiener)
{
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;

	Bparticle*		part = part_find_first(project);
	Bmicrograph*	mg = part->mg;
	Bfield*			field = project->field;
	Bimage*			p = NULL;
	Vector3<double>	pixel_size(part->pixel_size);

	if ( part->fpart.length() ) p = read_img(part->fpart, 1, 0);
	else p = read_img(mg->fpart, 1, part->id - 1);
	
	if ( pixel_size[0] < 0.01 ) pixel_size = p->sampling(0);
	if ( pixel_size[0] < 0.01 ) pixel_size = mg->pixel_size;
	
	Bimage*			pref = new Bimage(Float, TSimple, p->size(), number);
	pref->sampling(pixel_size);
	pref->origin(pref->size()/2);
	
	delete p;

	if ( verbose ) {
		cout << "Preparing particle images as references: " << endl;
		cout << "Size:                           " << pref->size() << endl;
		cout << "Pixel size:                     " << pixel_size << endl;
		cout << "Number:                         " << number << endl;
		cout << "Binning:                        " << bin << endl;
	}

	long			i(0), npart(0);
	
	for ( field = project->field; field && npart < number; field = field->next ) {
		for ( mg = field->mg; mg && npart < number; mg = mg->next ) {
			for ( part = mg->part; part && npart < number; part = part->next, ++i ) 
					if ( i > first ) {
				if ( part->fpart.length() ) p = read_img(part->fpart, 1, 0);
				else p = read_img(mg->fpart, 1, part->id - 1);
				if ( !p ) {
					error_show("project_prepare_2D_references", __FILE__, __LINE__);
					return NULL;
				}
				p->change_type(Float);
				p->sampling(pixel_size);
				if ( ctf_action && mg->ctf )
					img_ctf_apply(p, *mg->ctf, ctf_action, wiener, 0, 0);
				p->rescale_to_avg_std(0, 1);
				p->shift_background(0);
				if ( part->ori[0] > 0 && part->ori[1] > 0 ) p->origin(part->ori);
				else p->origin(p->size()/2);
				p->view(part->view);
				if ( verbose )
					cout << mg->id << tab << part->id << tab << p->image->origin() << tab << p->image->view() << endl;
				p->center();
				pref->replace(npart, p);
				npart++;
				delete p;
			}
		}
	}

	if ( bin > 1 ) pref->bin(bin);
	
	pref->statistics();
	pref->shift_background(0);
	
	if ( verbose & VERB_DEBUG )
		write_img("ref.pif", pref, 0);
	
	return pref;
}

/**
@brief 	Generate projections if a 3D file, otherwise clean projections.
@param 	&filename			file containing reference map or projections.
@param 	&mask_file			mask to apply to projections.
@param 	bin						data compression by binning.
@param 	*sym				point group symmetry structure.
@param 	theta_step			angular step size from primary symmetry axis (radians).
@param 	phi_step				angular step size around primary symmetry axis (radians).
@param 	side_ang				angular devaition from eqautor (radians).
@return Bimage*						projection images, NULL on error.

	If the input file is a 3D map, a set of projections are generated
	given the point group symmetry.
	Flags:
		FULL_ASU	projections for full asymmetric unit
		MULTI_FILE	projections in multiple files	

**/
Bimage*		img_prepare_projections(Bstring& filename, Bstring& mask_file,
				int bin, Bsymmetry& sym,
				double theta_step, double phi_step, double side_ang)
{
	if ( theta_step < SMALLFLOAT || phi_step < SMALLFLOAT ) {
		if ( theta_step < SMALLFLOAT )
			cerr << "Error: The theta step size is too small! (" << theta_step*180.0/M_PI << ")" << endl;
		if ( phi_step < SMALLFLOAT )
			cerr << "Error: The phi step size is too small! (" << phi_step*180.0/M_PI << ")" << endl;
		return NULL;
	}

	// Read the reference
	Bimage* 		pref = read_img(filename, 1, -1);
	if ( pref == NULL ) {
		error_show("Error in img_prepare_projections", __FILE__, __LINE__);
		return NULL;
	}
	
	pref->change_type(Float);
	pref->rescale_to_avg_std(0, 1);
	pref->statistics();
		
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	if ( bin > 1 ) pref->bin(bin);

	int				i;
	Vector3<double>	origin;
	for ( i=0; i<pref->images(); i++ ) {
		origin = pref->image[i].origin();
		if ( fabs(origin[0] - pref->sizeX()/2) > pref->sizeX()/10 )
			origin[0] = pref->sizeX()/2;
		if ( fabs(origin[1] - pref->sizeY()/2) > pref->sizeY()/10 )
			origin[1] = pref->sizeY()/2;
		if ( fabs(origin[2] - pref->sizeZ()/2) > pref->sizeZ()/10 )
			origin[2] = pref->sizeZ()/2;
		if ( pref->sizeZ() < 2 ) origin[2] = 0;
		pref->image[i].origin(origin);
	}

	if ( verbose ) {
		cout << "Preparing projections for the full asymmetric unit: " << endl;
		cout << "Reference:                      " << pref->file_name() << endl;
		if ( pref->images() > 1 ) cout << " (" << pref->images() << " projections)" << endl;
		cout << "Reference origin:               " << pref->image->origin() << endl;
		cout << "Theta and phi steps:            " << theta_step*180.0/M_PI << " " << phi_step*180.0/M_PI << endl;
		if ( side_ang < 0 )
			cout << "Symmetry:                       " << sym.label() << endl;
		else
			cout << "Side views within:              " << side_ang*180.0/M_PI << " degrees" << endl;
    	if ( mask_file.length() ) cout << "Mask file name:                 " << mask_file << endl;
	}

	Bimage*			proj = NULL;
	View*			views = NULL;
	int				nviews(0);
	
	time_t			t = time(NULL);

	FSI_Kernel*		kernel = new FSI_Kernel(10, 2);

	// Mask to apply to projections
	Bimage*		pmask = NULL;
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, 0);
		if ( bin > 1 ) pmask->bin(bin);
	}
	
	long	mem_req(0);
	long	mem_avail = system_memory();

	if ( pref->sizeZ() > 1 && pref->images() < 2 ) {
		if ( side_ang < 0 )
			views = asymmetric_unit_views(sym, theta_step, phi_step, 1);
		else
			views = side_views(sym, side_ang, theta_step, phi_step);
		nviews = count_list((char *)views);
		mem_req = pref->size()[0]*pref->size()[1]*sizeof(float)*nviews;
		if ( verbose ) {
			cout << "Number of projections:          " << nviews << endl;
			cout << "Memory required:                " << mem_req << endl;
			cout << "Memory available:               " << mem_avail << endl;
		}
		if ( mem_req > mem_avail ) {
			cerr << "Not enough memory to run!" << endl;
			bexit(-1);
		}
		proj = pref->project(views, pref->sampling(0)[0], kernel);
		if ( pmask ) proj->multiply(pmask);
		kill_list((char *) views, sizeof(View));
		delete pref;
		if ( verbose )
			cout << "Time to generate projections:   " << (long)(time(NULL) - t) << " seconds" << endl << endl;
	} else {
		proj = pref;
		proj->shift_background(0);
		mem_req = proj->size()[0]*proj->size()[1]*proj->images()*sizeof(float);
		if ( verbose ) {
			cout << "Number of projections:          " << proj->images() << endl;
			cout << "Memory required:                " << mem_req << endl;
			cout << "Memory available:               " << mem_avail << endl;
		}
		if ( mem_req > mem_avail ) {
			cerr << "Not enough memory to run!" << endl;
			bexit(-1);
		}
	}

	delete kernel;
	delete pmask;
	
	return proj;
}

Bimage*		img_read_projection(Bstring& proj_file, int multi_file, int proj_num)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_read_projection: proj_file=" << proj_file << " proj_num=" << proj_num << endl;
	
	Bimage*			proj = NULL;
	Bstring			proj_name = NULL;
	
	if ( multi_file ) {
		proj_name = proj_file.pre_rev('.') + Bstring(proj_num, "%04d.") + proj_file.post_rev('.');
		proj = read_img(proj_name, 1, 0);
	} else
		proj = read_img(proj_file, 1, proj_num);
	
	if ( !proj )
		error_show("Error in img_read_projection", __FILE__, __LINE__);
	
	return proj;
}

/**
@brief 	Find the orientation and origin of every image in a project.
@param 	*project		image processing parameter structure.
@param 	*proj			reference projections.
@param 	&mask_file		mask to apply to particles.
@param 	bin				data compression by binning.
@param 	*sym			point group symmetry structure.
@param 	part_select		particle selection for processing.
@param 	*band			array of reciprocal space bands.
@param 	res_lo			low resolution limit (angstrom).
@param 	res_hi			high resolution limit (angstrom).
@param 	res_polar       resolution limit for in-plane angular alignment (angstrom).
@param 	ann_min			minimum annulus (>=0).
@param 	ann_max			maximum annulus (< image radius).
@param 	shift_limit		maximum shift from nominal origin of box.
@param 	angle_limit		maximum rotation from original in-plane rotation angle.
@param 	edge_radius		edge radius to mask background using previous origin.
@param 	flags			option flags.
@return int				error code.

	The polar power spectrum (pps) of the reference projection is cross correlated
	with that of the image in order to find the angle of rotation.  The image
	is transformed using this angle and the shift found by cross correlation (cc).
	A  second iterative comparison is done with real space polar images
	to find the angle followed by cross correlation to find the shift.
	How much of this second comparison is done depends on the mode flag:
		mode=0	pps		projection selected only based on pps comparison
		mode=1	scc		several projections selected based on cutoff for pps cc's
		mode=2	ccc		all projections selected
	The angle and the x and y values are stored in the view_angle, and ox and oy 
	arrays of the micrograph parameter structure.
	The projections must already be binned.
	Flags:
		MODE		projection matching mode
		APPLY_CTF	apply CTF to projections
		PART_LOG	write log files in log directory

**/
int 		project_determine_orientations(Bproject* project, Bimage* proj, Bstring& mask_file,
				int bin, Bsymmetry& sym, int part_select, vector<double>& band,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags)
{
	int				mode = flags & MODE;
	int				ctf_apply = (flags & APPLY_CTF)? 1: 0;
	int				part_log = (flags & PART_LOG)? 1: 0;
	
//	Bfield*			field = project->field;

	Bparticle*		part = part_find_first(project);
	Bmicrograph*	mg = part->mg;

	if ( !part ) {
		cerr << "Error in project_determine_orientations: No particles found!" << endl;
		return -1;
	}
	
//	cout << proj->image->origin() << endl;
	
	// Note that the reference projections are already binned
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	edge_radius /= bin;
	
	Vector3<double>	pixel_size = mg->pixel_size*bin;
	
	// Read the reference (model) projections, convert to floating point and clean up
	long 			i, nproj(0), npix(0);
	
	nproj = proj->images();
	
	npix = proj->sizeX()*proj->sizeY();
	
//	cout << "Projection size:" << tab << proj->size() << tab << proj->images() << endl;
	
    // Specify FFTW arrays
	fft_plan		planf_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_FORWARD, 1);
	fft_plan		planb_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_BACKWARD, 1);	
	fft_plan		planf_2D = fft_setup_plan(proj->sizeX(), proj->sizeY(), 1, FFTW_FORWARD, 1);
	fft_plan		planb_2D = fft_setup_plan(proj->sizeX(), proj->sizeY(), 1, FFTW_BACKWARD, 1);

	// Set resolution limits
	if ( res_lo < res_hi ) swap(res_hi, res_lo);
	if ( res_lo > proj->sizeX()*mg->pixel_size[0] )
		res_lo = proj->sizeX()*mg->pixel_size[0];
	if ( res_hi > proj->sizeX()*mg->pixel_size[0]/2 )
		res_hi = proj->sizeX()*mg->pixel_size[0]/2;
	if ( res_lo < 4*pixel_size[0] ) res_lo = 4*pixel_size[0];
	if ( res_hi < 2*pixel_size[0] ) res_hi = 2*pixel_size[0];
	
//	cout << "Resolution:" << tab << res_hi << tab << res_lo << endl;

	// Check the annulus limits
	ann_min /= bin;
	ann_max /= bin;
	if ( ann_min > ann_max ) swap(ann_min, ann_max);
	if ( ann_min < 0 ) ann_min = 0;
	if ( ann_max > proj->sizeX()/2 - 1 ) ann_max = proj->sizeX()/2 - 1;
	
//	cout << "Annuli:" << tab << ann_min << tab << ann_max << endl;

	// Check the shift limit
	shift_limit /= bin;
	if ( shift_limit < 0 ) shift_limit = proj->sizeX()/10;	// Default shift limit
	
	// Generate the masks for cross-correlation and cross-validation
	Bimage*		prs_mask = new Bimage(SCharacter, TSimple, proj->size(), 1);
	prs_mask->sampling(pixel_size);
	if ( band.size() < 2 || band[0] < 1 )
		band = prs_mask->fspace_default_bands(res_lo, res_hi);

	prs_mask->mask_fspace_banded(band);
	write_img("mask.map", prs_mask, 0);
	
	if ( prs_mask->minimum() < 0 ) project->fom_tag[1] = FOM_CV;
	
	// Mask to apply to particles
	Bimage*		part_mask = NULL;
	if ( mask_file.length() ) {
		part_mask = read_img(mask_file, 1, 0);
		if ( bin > 1 ) part_mask->bin(bin);
	}

	long			npart(0);
	Bparticle**		partarr = project_mg_particle_array(project, part_select, npart);
	
	if ( verbose & VERB_RESULT ) {
		cout << "Determining orientations:" << endl;
		if ( mode == 0 ) cout << "Mode:                           PPS based followed by one CC (fast)" << endl;
		if ( mode == 1 ) cout << "Mode:                           PPS followed by selected CC's" << endl;
		if ( mode == 2 ) cout << "Mode:                           PPS followed by all CC's (slow)" << endl;
    	cout << "Projection file name:           " << proj->file_name() << endl;
    	cout << "Number of projections:          " << nproj << endl;
    	cout << "Binning:                        " << bin << endl;
		cout << "Pixel size:                     " << pixel_size << endl;
		if ( prs_mask->minimum() < 0 ) cout << "Reciprocal space mask file:     " << prs_mask->file_name() << endl;
		if ( band[0] > 0 ) {
			cout << "Reciprocal space bands:" << endl;
			for ( i=0; band[i+1]; i+=2 ) {
				cout << tab << setw(5) << band[i] << " - " << setw(5) << band[i+2];
				if ( band[i+1] > 0 ) cout << "\tcross-correlation" << endl;
				else cout << "\tcross-validation" << endl;
			}
		}
    	if ( ctf_apply ) cout << "CTF applied to projections." << endl;
    	cout << "Resolution limits (low, high):  " << res_lo << " - " << res_hi << " A" << endl;
    	if ( res_polar ) cout << "Polar resolution limit:         " << res_polar << " A" << endl;
    	cout << "Annular range:                  " << ann_min << " - " << ann_max << " (max = " << proj->sizeX()/2 - 1 << ")" << endl;
    	cout << "Shift limit:                    " << shift_limit << " pixels" << endl;
    	if ( angle_limit )
			cout << "Angle limit:                    " << angle_limit*180.0/M_PI << " degrees" << endl;
    	if ( edge_radius ) cout << "Edge mask radius:               " << edge_radius << " pixels" << endl;
    	if ( mask_file.length() ) cout << "Real space mask file:           " << mask_file << endl;
    	cout << endl;
		cout << "Micrograph\tPID\tOriX\tOriY\tViewX\tViewY\tViewZ\tAngle\tCC\tCV\tTime" << endl;
	}
	
	// Directory for individual particle log files
	if ( part_log ) {
		mkdir("log", O_CREAT );
		chmod("log", 0755);
	}
	
	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	time_t			t = time(NULL);

#ifdef HAVE_GCD
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(npart, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bparticle*		part = partarr[i];
		Bmicrograph*	mg = part->mg;
		part_determine_orientation(part, proj, part_mask,
				prs_mask, sym, bin,
				res_lo, res_hi, res_polar, ann_min, ann_max,
				shift_limit, angle_limit, edge_radius, flags,
				planf_1D, planb_1D, planf_2D, planb_2D);
		dispatch_sync(myq, ^{
			if ( verbose & VERB_RESULT ) {
				cout << setw(15) << mg->id << tab << part->id << tab << fixed <<
					setprecision(2) << part->ori[0] << tab << part->ori[1] << tab <<
					setprecision(4) << part->view[0] << tab <<
					part->view[1] << tab << part->view[2] << tab <<
					setprecision(2) << part->view.angle()*180.0/M_PI << tab <<
					setprecision(4) << part->fom[0] << tab << part->fom[1] << tab <<
					(long)(time(NULL) - t) << endl << flush;
			}
		});
	});
#else
#pragma omp parallel for
	for ( i=0; i<npart; i++ ) {
		Bparticle*		part = partarr[i];
		Bmicrograph*	mg = part->mg;
		part_determine_orientation(part, proj, part_mask,
				prs_mask, sym, bin,
				res_lo, res_hi, res_polar, ann_min, ann_max,
				shift_limit, angle_limit, edge_radius, flags,
				planf_1D, planb_1D, planf_2D, planb_2D);
	#pragma omp critical
		{
			if ( verbose & VERB_RESULT ) {
				cout << setw(15) << mg->id << tab << part->id << tab << fixed <<
					setprecision(2) << part->ori[0] << tab << part->ori[1] << tab <<
					setprecision(4) << part->view[0] << tab <<
					part->view[1] << tab << part->view[2] << tab <<
					setprecision(2) << part->view.angle()*180.0/M_PI << tab <<
					setprecision(4) << part->fom[0] << tab << part->fom[1] << tab <<
					(long)(time(NULL) - t) << endl << flush;
			}
		}
	}
#endif

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_determine_orientations: done, npart=" << npart << endl;

	t = (time(NULL) - t);

	long		ts = (long) t;
	long		hr = ts/3600, min = (ts - 3600*hr)/60, sec = ts - 3600*hr - 60*min;
	
	if ( verbose & VERB_RESULT && npart ) {
		cout << "Time:                           " << ts << " s (" << hr << ":" << min << ":" << sec << ")" << endl;
		cout << "Time per particle:              " << ts*1.0/npart << " s/particle" << endl;
		cout << "Algorithm time:                 " << ts*1.0e6/(npart*nproj*npix) << " us/pixel" << endl << endl;
	}

	delete[] partarr;
	delete prs_mask;
	delete part_mask;
	
    fft_destroy_plan(planf_1D);
    fft_destroy_plan(planb_1D);
    fft_destroy_plan(planf_2D);
    fft_destroy_plan(planb_2D);
	
	return 0;
}

int 		part_determine_orientation(Bparticle* part, Bimage* proj, Bimage* part_mask,
				Bimage* prs_mask, Bsymmetry& sym, int bin,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags,
				fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D)
{
	FOMType 	fom_tag[NFOM] = {FOM, FOM_CV};

	if ( flags & CHECK_PPX )
		if ( ppx_check(part, fom_tag) ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: id=" << part->id << " flags=" << flags << endl;
		
	int				mode = flags & MODE;
	int				ctf_apply = (flags & APPLY_CTF)? 1: 0;
	int				part_log = (flags & PART_LOG)? 1: 0;

	long 			k, imax(0), nproj(proj->images());
	double			cc_best(-1e37), cc_min, cc_max, cc_avg, cc_std, cc_cut;
	double*			angle = new double[nproj];
	double*			ox = new double[nproj];
	double*			oy = new double[nproj];
	View*			view = new View[nproj];
	double*			cc = new double[nproj];
	
	Bmicrograph*	mg = part->mg;
	if ( !mg )
		return error_show("part_determine_orientation", __FILE__, __LINE__);

	Vector3<double>	pixel_size = mg->pixel_size*bin;
	ofstream		flog;
	Bstring			log_name;
	long			edge_size = (long) (2*edge_radius);
	Vector3<long>	size(edge_size, edge_size, 1);
	Vector3<double>	start, shift;

	Bimage*			p = NULL;
	Bimage*			pm;
	Bimage*			proj1;
	
	if ( part_log ) {
		log_name = "log/" + mg->id + "_" + Bstring(part->id, "%04d") + "_orient.log";
		flog.open(log_name.c_str());
		if ( flog.fail() )
			return error_show("Log file cannot be written!", __FILE__, __LINE__);
		flog << mg->id << ": " << part->id << endl;
	}
	
//	cout << "Reading particle " << part->id << endl;
	if ( part->fpart.length() ) p = read_img(part->fpart, 1, 0);
	else if ( mg->fpart.length() ) p = read_img(mg->fpart, 1, part->id - 1);
	if ( !p )
		return error_show("part_determine_orientation", __FILE__, __LINE__);
	
//	cout << "Preparing particle " << part->id << endl;
	p->change_type(Float);
	if ( bin > 1 ) p->bin(bin);
	p->statistics();
	p->sampling(pixel_size);
//	if ( part->id == 1 ) cout << "Image sampling = " << p->sampling(0) << endl;
//	if ( part->ori[0] > 0 && part->ori[1] > 0 ) p->image->origin(part->ori/bin);
//	else p->image->origin(p->size());
	if ( part->ori[0] > 0 && part->ori[1] > 0 ) p->origin(part->ori/bin);
	else p->origin(p->size()/2);
	
	p->rescale_to_avg_std(0, 1);
	p->shift_background(0);
	
	if ( part_mask ) {
		pm = part_mask->copy();
		shift = p->image->origin() - pm->image->origin();
		pm->shift(shift);
		p->multiply(pm);
		delete pm;
		p->rescale_to_avg_std(0, 1);
		p->shift_background(0);
	}
	
	if ( edge_size ) {
//		cout << "Edge masking" << endl;
		start = Vector3<double>(p->image->origin()[0]-edge_radius, p->image->origin()[1]-edge_radius, 0);
		p->edge(1, size, start, 1, FILL_BACKGROUND, 0);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: preprocessing done for image " << part->id << endl;
	if ( part_log )
		flog << "PartID\tProjID\tox\toy\tvx\tvy\tvz\tva\tcc" << endl;

	for ( k=0; k<nproj; k++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG part_determine_orientation: pps with proj " << k << endl;
		proj1 = proj->extract(k);
		if ( !proj1 )
			return error_show("Error in part_determine_orientation", __FILE__, __LINE__);
		proj1->sampling(pixel_size);
//		cout << "Projection sampling = " << proj1->sampling(0) << endl;
		if ( ctf_apply && mg->ctf )
			img_ctf_apply_to_proj(proj1, *(mg->ctf), part->def, 1e6, res_hi, planf_2D, planb_2D);
		p->image->view(part->view);
		cc[k] = p->align2D_pps(proj1, res_hi, res_lo, shift_limit, angle_limit, planf_2D, planb_2D);
		angle[k] = p->image->view_angle();
		ox[k] = p->image->origin()[0];
		oy[k] = p->image->origin()[1];
		view[k] = proj1->image->view();
//		if ( part->id == 43 ) cout << ox[k] << tab << oy[k] << tab << cc[k] << endl;
		if ( angle[k] < 0 ) {
			view[k].negate();
			view[k][3] = angle_set_negPI_to_PI(M_PI - view[k].angle() + angle[k]);
		} else {
			view[k][3] = angle_set_negPI_to_PI(view[k].angle() - angle[k]);
		}
//		if ( verbose & VERB_DEBUG )
			if ( part_log )
				flog << fixed << setprecision(4) <<
					part->id << tab << k+1 << tab <<
					bin*ox[k] << tab << bin*oy[k] << tab << view[k][0] << tab << view[k][1] << tab <<
					view[k][2] << tab << view[k].angle()*180.0/M_PI << tab << cc[k] << endl;
//		cout << part->id << " " << k+1 << " " << ox[k] << " " << oy[k]
//			<< " " << view[k][0] << " " <<  view[k][1] << " " << view[k][2]
//			<< " " << view[k].angle()*180.0/M_PI << " " <<  cc[k] << endl;
		if ( cc_best < cc[k] ) {
			cc_best = cc[k];
			imax = k;
		}
		delete proj1;
	}
	
	cc_avg = cc_std = 0;
	cc_min = cc_max = cc[0];
	for ( k=imax=0; k<nproj; k++ ) {
		if ( cc_min > cc[k] ) cc_min = cc[k];
		if ( cc_max < cc[k] ) {
			cc_max = cc[k];
			imax = k;
		}
		cc_avg += cc[k];
		cc_std += cc[k]*cc[k];
	}
	cc_avg /= nproj;
	cc_std = cc_std/nproj - cc_avg*cc_avg;
	if ( cc_std > 0 ) cc_std = sqrt(cc_std);
	else cc_std = 0;
	cc_cut = cc_avg + cc_std;
	if ( cc_cut >= cc_max ) cc_cut = cc_max - 0.01;
	if ( mode < 1 ) cc_cut = cc_max - 0.01;
	else if ( mode > 1 ) cc_cut = cc_min - 1;
	if ( verbose & VERB_DEBUG )
		if ( part_log )
			flog << "Particle, min, max, avg, std, max/avg:" <<
				part->id << tab << cc_min << tab << cc_max << tab <<
				cc_avg << tab << cc_std << tab << cc_max/cc_avg << endl;
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG part_determine_orientation: cc_avg=" << cc_avg << " cc_std=" << cc_std << endl;
		cout << "DEBUG part_determine_orientation: cc_min=" << cc_min << " cc_max=" << cc_max << endl;
		cout << "DEBUG part_determine_orientation: cc_cut=" << cc_cut << endl;
	}
	
	cc_best = -1e37;
//	if ( mode > 0 )
	for ( k=0; k<nproj; k++ ) if ( cc[k] > cc_cut ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG part_determine_orientation: cc with proj " << k << endl;
		proj1 = proj->extract(k);
		if ( !proj1 )
			return error_show("Error in part_determine_orientation", __FILE__, __LINE__);
		proj1->sampling(pixel_size);
		if ( ctf_apply && mg->ctf )
			img_ctf_apply_to_proj(proj1, *(mg->ctf), part->def, 1e6, res_hi, planf_2D, planb_2D);
		p->image->view(part->view);
		p->image->view_angle(angle[k]);
		p->image->origin(ox[k], oy[k], 0);
//		cc[k] = img_align2D(p, proj1, 1, res_polar, ann_min, ann_max,
//			prs_mask, shift_limit, planf_1D, planb_1D, planf_2D, planb_2D);
		cc[k] = p->align2D(proj1, res_polar, ann_min, ann_max,
			prs_mask, shift_limit, angle_limit, planf_1D, planb_1D, planf_2D, planb_2D);
		angle[k] = p->image->view_angle();
		ox[k] = p->image->origin()[0];
		oy[k] = p->image->origin()[1];
		view[k] = proj1->image->view();
//		if ( part->id == 43 ) {
//			cout << k << tab << ox[k] << tab << oy[k] << tab << cc[k] << endl;
//			if ( cc[k] > 1 ) bexit(-1);
//		}
		if ( angle[k] < 0 ) {
			view[k].negate();
			view[k][3] = angle_set_negPI_to_PI(M_PI - view[k].angle() + angle[k]);
		} else {
			view[k][3] = angle_set_negPI_to_PI(view[k].angle() - angle[k]);
		}
		if ( part_log )
			flog << fixed << setprecision(4) <<
				part->id << tab << k+1 << tab <<
				bin*ox[k] << tab << bin*oy[k] << tab << view[k][0] << tab << view[k][1] << tab <<
				view[k][2] << tab << view[k].angle()*180.0/M_PI << tab << cc[k] << endl;
//		cout << part->id << " " << k+1 << " " << ox[k] << " " << oy[k]
//			<< " " << view[k][0] << " " <<  view[k][1] << " " << view[k][2]
//			<< " " << view[k].angle()*180.0/M_PI << " " <<  cc[k] << endl;
		if ( cc_best < cc[k] ) {
			cc_best = cc[k];
			imax = k;
		}
		delete proj1;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: projection matching done for image " << part->id << endl;
	part->ori = Vector3<double>(bin*ox[imax], bin*oy[imax], 0);
	part->view = view[imax];
	part->fom[0] = cc[imax];
	if ( prs_mask->minimum() < 0 ) {
		proj1 = proj->extract(imax);
		if ( !proj1 )
			return error_show("Error in part_determine_orientation", __FILE__, __LINE__);
		proj1->sampling(pixel_size);
		if ( ctf_apply && mg->ctf )
			img_ctf_apply_to_proj(proj1, *(mg->ctf), part->def, 1e6, res_hi, planf_2D, planb_2D);
		p->image->view_angle(angle[imax]);
		p->image->origin(ox[imax], oy[imax], 0);
		p->origin(0, ox[imax], oy[imax], 0.0);
		part->fom[1] = img_cross_validate(p, proj1, prs_mask, planf_2D);
		delete proj1;
	}
	part->view = find_asymmetric_unit_view(sym, part->view);
	part->view[3] = angle_set_negPI_to_PI(part->view.angle());
	part->sel = imax + 1;
	if ( part_log ) {
		flog << "Particle, min, max, avg, std, max/avg, cv:" << tab <<
			part->id << tab << cc_min << tab << cc_max << tab <<
			cc_avg << tab << cc_std << tab << cc_max/cc_avg << tab << part->fom[1] << endl;
		flog << fixed << setprecision(4) <<
			part->id << tab << part->ori[0] << tab << part->ori[1] << tab <<
			part->view[0] << tab << part->view[1] << tab <<
			part->view[2] << tab << part->view.angle()*180.0/M_PI << tab <<
			part->fom[0] << tab << part->fom[1] << endl;
		flog.close();
	}

	delete[] angle;
	delete[] ox;
	delete[] oy;
	delete[] view;
	delete[] cc;
	delete p;
	
	if ( flags & WRITE_PPX ) {
		log_name = ppx_filename(mg->id, part->id);
		write_particle(log_name, part, 0, 0, fom_tag);
	}
	
	return 1;
}

Bparticle	part_compare(Bimage* p, Bimage* proj, Bimage* prs_mask,
	long k, int bin, CTFparam* ctf,
	double res_lo, double res_hi, double res_polar, long ann_min, long ann_max,
	double shift_limit, double angle_limit, double edge_radius, int flags,
	fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D)
{
	Bimage*			pcopy = p->copy();
	Bparticle		part;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_compare: proj " << k << endl;
	
	Bimage*			proj1 = proj->extract(k);
	if ( !proj1 ) {
		error_show("Error in part_compare", __FILE__, __LINE__);
		return part;
	}
	
	proj1->sampling(p->sampling(0));
	
	if ( ctf )
		img_ctf_apply_to_proj(proj1, *ctf, ctf->defocus_average(), 1e6, res_hi, planf_2D, planb_2D);
	
	part.fom[0] = pcopy->align2D_pps(proj1, res_hi, res_lo, shift_limit, angle_limit, planf_2D, planb_2D);
	

	part.fom[0] = pcopy->align2D(proj1, res_polar, ann_min, ann_max,
				prs_mask, shift_limit, angle_limit, planf_1D, planb_1D, planf_2D, planb_2D);

	View			view = proj1->image->view();
	double			angle = pcopy->image->view_angle();
				
	if ( angle < 0 ) {
		view.negate();
		view[3] = angle_set_negPI_to_PI(M_PI - view.angle() + angle);
	} else {
		view[3] = angle_set_negPI_to_PI(view.angle() - angle);
	}

	part.ori = pcopy->image->origin() * bin;
	part.view = view;
	
	delete pcopy;
	delete proj1;

	return part;
}

int 		part_determine_orientation2(Bparticle* part, Bimage* proj, Bimage* part_mask,
				Bimage* prs_mask, Bsymmetry& sym, int bin,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags,
				fft_plan planf_1D, fft_plan planb_1D, fft_plan planf_2D, fft_plan planb_2D)
{
	FOMType 	fom_tag[NFOM] = {FOM, FOM_CV};

	if ( flags & CHECK_PPX )
		if ( ppx_check(part, fom_tag) ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: id=" << part->id << " flags=" << flags << endl;
		
	int				ctf_apply = (flags & APPLY_CTF)? 1: 0;
	int				part_log = (flags & PART_LOG)? 1: 0;

	long 			k, imax(0), nproj(proj->images());
	
	Bmicrograph*	mg = part->mg;
	if ( !mg )
		return error_show("part_determine_orientation", __FILE__, __LINE__);

	Vector3<double>	pixel_size = part->pixel_size*bin;
	ofstream		flog;
	Bstring			log_name;
	long			edge_size = (long) (2*edge_radius);
	Vector3<long>	size(edge_size, edge_size, 1);
	Vector3<double>	start, shift;

	Bimage*			p = NULL;
	Bimage*			pm;
	Bimage*			proj1;
	
	if ( part_log ) {
		log_name = "log/" + mg->id + "_" + Bstring(part->id, "%04d") + "_orient.log";
		flog.open(log_name.c_str());
		if ( flog.fail() )
			return error_show("Log file cannot be written!", __FILE__, __LINE__);
		flog << mg->id << ": " << part->id << endl;
	}
	
//	cout << "Reading particle " << part->id << endl;
	if ( part->fpart.length() ) p = read_img(part->fpart, 1, 0);
	else if ( mg->fpart.length() ) p = read_img(mg->fpart, 1, part->id - 1);
	if ( !p )
		return error_show("part_determine_orientation", __FILE__, __LINE__);
	
//	cout << "Preparing particle " << part->id << endl;
	p->change_type(Float);
	if ( bin > 1 ) p->bin(bin);
	p->statistics();
	p->sampling(pixel_size);
//	if ( part->id == 1 ) cout << "Image sampling = " << p->sampling(0) << endl;
//	if ( part->ori[0] > 0 && part->ori[1] > 0 ) p->image->origin(part->ori/bin);
//	else p->image->origin(p->size());
	if ( part->ori[0] > 0 && part->ori[1] > 0 ) p->origin(part->ori/bin);
	else p->origin(p->size()/2);
	p->image->view(part->view);
	
	p->rescale_to_avg_std(0, 1);
	p->shift_background(0);
	
	if ( part_mask ) {
		pm = part_mask->copy();
		shift = p->image->origin() - pm->image->origin();
		pm->shift(shift);
		p->multiply(pm);
		delete pm;
		p->rescale_to_avg_std(0, 1);
		p->shift_background(0);
	}
	
	if ( edge_size ) {
//		cout << "Edge masking" << endl;
		start = Vector3<double>(p->image->origin()[0]-edge_radius, p->image->origin()[1]-edge_radius, 0);
		p->edge(1, size, start, 1, FILL_BACKGROUND, 0);
	}
	
	Bparticle*		partarr = new Bparticle[nproj];
//	vector<Bparticle>	partarr(nproj);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: preprocessing done for image " << part->id << endl;
	if ( part_log )
		flog << "PartID\tProjID\tox\toy\tvx\tvy\tvz\tva\tcc" << endl;

#ifdef HAVE_GCD
	dispatch_apply(nproj, dispatch_get_global_queue(0, 0), ^(size_t k){
		Bmicrograph*	mg = part->mg;
		partarr[k] = part_compare(p, proj, prs_mask, k, bin, mg->ctf,
					res_lo, res_hi, res_polar, ann_min, ann_max,
					shift_limit, angle_limit, edge_radius, flags,
					planf_1D, planb_1D, planf_2D, planb_2D);
	});
#else
#pragma omp parallel for
	for ( k=0; k<nproj; ++k ) {
		Bmicrograph*	mg = part->mg;
		partarr[k] = part_compare(p, proj, prs_mask, k, bin, mg->ctf,
					res_lo, res_hi, res_polar, ann_min, ann_max,
					shift_limit, angle_limit, edge_radius, flags,
					planf_1D, planb_1D, planf_2D, planb_2D);
	}
#endif

	double			cc, cc_min(1), cc_avg(0), cc_std(0);

	part->fom[0] = -1;
	for ( k=0; k<nproj; ++k ) {
		cc = partarr[k].fom[0];
		cc_avg += cc;
		cc_std += cc * cc;
		if ( part->fom[0] < cc ) {
			part->fom[0] = cc;
			imax = k;
		}
		if ( cc_min > cc ) cc_min = cc;
	}
	cc_avg /= nproj;
	cc_std /= nproj;
	cc_std -= cc_avg * cc_avg;
	if ( cc_std > 0 ) cc_std = sqrt(cc_std);
	else cc_std = 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_determine_orientation: projection matching done for image " << part->id << endl;
	part->ori = partarr[imax].ori;
	part->view = partarr[imax].view;
	
	if ( prs_mask->minimum() < 0 ) {
		proj1 = proj->extract(imax);
		if ( !proj1 )
			return error_show("Error in part_determine_orientation", __FILE__, __LINE__);
		proj1->sampling(pixel_size);
		if ( ctf_apply && mg->ctf )
			img_ctf_apply_to_proj(proj1, *(mg->ctf), part->def, 1e6, res_hi, planf_2D, planb_2D);
		p->image->origin(part->ori);
		p->image->view(part->view);
//		p->image->view_angle(angle[imax]);
		part->fom[1] = img_cross_validate(p, proj1, prs_mask, planf_2D);
		delete proj1;
	}
	
	part->view = find_asymmetric_unit_view(sym, part->view);
	part->view[3] = angle_set_negPI_to_PI(part->view.angle());
	part->sel = imax + 1;
	
	if ( part_log ) {
		flog << "Particle, min, max, avg, std, max/avg, cv:" << tab <<
			part->id << tab << cc_min << tab << part->fom[0] << tab <<
			cc_avg << tab << cc_std << tab << part->fom[0]/cc_avg << tab << part->fom[1] << endl;
		flog << fixed << setprecision(4) <<
			part->id << tab << part->ori[0] << tab << part->ori[1] << tab <<
			part->view[0] << tab << part->view[1] << tab <<
			part->view[2] << tab << part->view.angle()*180.0/M_PI << tab <<
			part->fom[0] << tab << part->fom[1] << endl;
		flog.close();
	}

	delete p;
	delete[] partarr;
	
	if ( flags & WRITE_PPX ) {
		log_name = ppx_filename(mg->id, part->id);
		write_particle(log_name, part, 0, 0, fom_tag);
	}
	
	return 1;
}


/**
@brief 	Find the orientation and origin of every image in a project.
@param 	*project		image processing parameter structure.
@param 	*proj			reference projections.
@param 	&mask_file		mask to apply to particles.
@param 	bin				data compression by binning.
@param 	*sym			point group symmetry structure.
@param 	part_select		particle selection for processing.
@param 	*band			array of reciprocal space bands.
@param 	res_lo			low resolution limit (angstrom).
@param 	res_hi			high resolution limit (angstrom).
@param 	res_polar       resolution limit for in-plane angular alignment (angstrom).
@param 	ann_min			minimum annulus (>=0).
@param 	ann_max			maximum annulus (< image radius).
@param 	shift_limit		maximum shift from nominal origin of box.
@param 	angle_limit		maximum rotation from original in-plane rotation angle.
@param 	edge_radius		edge radius to mask background using previous origin.
@param 	flags			option flags.
@return int				error code.

	The polar power spectrum (pps) of the reference projection is cross correlated
	with that of the image in order to find the angle of rotation.  The image
	is transformed using this angle and the shift found by cross correlation (cc).
	A  second iterative comparison is done with real space polar images
	to find the angle followed by cross correlation to find the shift.
	How much of this second comparison is done depends on the mode flag:
		mode=0	pps		projection selected only based on pps comparison
		mode=1	scc		several projections selected based on cutoff for pps cc's
		mode=2	ccc		all projections selected
	The angle and the x and y values are stored in the view_angle, and ox and oy
	arrays of the micrograph parameter structure.
	The projections must already be binned.
	Flags:
		MODE		projection matching mode
		APPLY_CTF	apply CTF to projections
		PART_LOG	write log files in log directory

**/
int 		project_determine_orientations2(Bproject* project, Bimage* proj, Bstring& mask_file,
				int bin, Bsymmetry& sym, int part_select, vector<double>& band,
				double res_lo, double res_hi, double res_polar, int ann_min, int ann_max,
				double shift_limit, double angle_limit, double edge_radius, int flags)
{
	int				mode = flags & MODE;
	int				ctf_apply = (flags & APPLY_CTF)? 1: 0;
	int				part_log = (flags & PART_LOG)? 1: 0;
	
//	Bfield*			field = project->field;

	Bparticle*		part = part_find_first(project);
	Bmicrograph*	mg = part->mg;

	if ( !part ) {
		cerr << "Error in project_determine_orientations: No particles found!" << endl;
		return -1;
	}
	
//	cout << proj->image->origin() << endl;
	
	// Note that the reference projections are already binned
	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;
	edge_radius /= bin;
	
	Vector3<double>	pixel_size = part->pixel_size*bin;
	
	// Read the reference (model) projections, convert to floating point and clean up
	long 			i, nproj(0), npix(0);
	
	nproj = proj->images();
	
	npix = proj->sizeX()*proj->sizeY();
	
//	cout << "Projection size:" << tab << proj->size() << tab << proj->images() << endl;
	
    // Specify FFTW arrays
	fft_plan		planf_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_FORWARD, 1);
	fft_plan		planb_1D = fft_setup_plan(NPOLANG, 1, 1, FFTW_BACKWARD, 1);
	fft_plan		planf_2D = fft_setup_plan(proj->sizeX(), proj->sizeY(), 1, FFTW_FORWARD, 1);
	fft_plan		planb_2D = fft_setup_plan(proj->sizeX(), proj->sizeY(), 1, FFTW_BACKWARD, 1);

	// Set resolution limits
	if ( res_lo < res_hi ) swap(res_hi, res_lo);
	if ( res_lo > proj->sizeX()*pixel_size[0] )
		res_lo = proj->sizeX()*pixel_size[0];
	if ( res_hi > proj->sizeX()*pixel_size[0]/2 )
		res_hi = proj->sizeX()*pixel_size[0]/2;
	if ( res_lo < 4*pixel_size[0] ) res_lo = 4*pixel_size[0];
	if ( res_hi < 2*pixel_size[0] ) res_hi = 2*pixel_size[0];
	
//	cout << "Resolution:" << tab << res_hi << tab << res_lo << endl;

	// Check the annulus limits
	ann_min /= bin;
	ann_max /= bin;
	if ( ann_min > ann_max ) swap(ann_min, ann_max);
	if ( ann_min < 0 ) ann_min = 0;
	if ( ann_max > proj->sizeX()/2 - 1 ) ann_max = proj->sizeX()/2 - 1;
	
//	cout << "Annuli:" << tab << ann_min << tab << ann_max << endl;

	// Check the shift limit
	shift_limit /= bin;
	if ( shift_limit < 0 ) shift_limit = proj->sizeX()/10;	// Default shift limit
	
	// Generate the masks for cross-correlation and cross-validation
	Bimage*		prs_mask = new Bimage(SCharacter, TSimple, proj->size(), 1);
	prs_mask->sampling(pixel_size);
	if ( band.size() < 2 || band[0] < 1 )
		band = prs_mask->fspace_default_bands(res_lo, res_hi);

	prs_mask->mask_fspace_banded(band);
	write_img("mask.map", prs_mask, 0);
	
	if ( prs_mask->minimum() < 0 ) project->fom_tag[1] = FOM_CV;
	
	// Mask to apply to particles
	Bimage*		part_mask = NULL;
	if ( mask_file.length() ) {
		part_mask = read_img(mask_file, 1, 0);
		if ( bin > 1 ) part_mask->bin(bin);
	}

	long			npart(0);
	Bparticle**		partarr = project_mg_particle_array(project, part_select, npart);
	
	if ( verbose & VERB_RESULT ) {
		cout << "Determining orientations:" << endl;
		if ( mode == 0 ) cout << "Mode:                           PPS based followed by one CC (fast)" << endl;
		if ( mode == 1 ) cout << "Mode:                           PPS followed by selected CC's" << endl;
		if ( mode == 2 ) cout << "Mode:                           PPS followed by all CC's (slow)" << endl;
    	cout << "Projection file name:           " << proj->file_name() << endl;
    	cout << "Number of projections:          " << nproj << endl;
    	cout << "Binning:                        " << bin << endl;
		cout << "Pixel size:                     " << pixel_size << endl;
		if ( prs_mask->minimum() < 0 ) cout << "Reciprocal space mask file:     " << prs_mask->file_name() << endl;
		if ( band[0] > 0 ) {
			cout << "Reciprocal space bands:" << endl;
			for ( i=0; band[i+1]; i+=2 ) {
				cout << tab << setw(5) << band[i] << " - " << setw(5) << band[i+2];
				if ( band[i+1] > 0 ) cout << "\tcross-correlation" << endl;
				else cout << "\tcross-validation" << endl;
			}
		}
    	if ( ctf_apply ) cout << "CTF applied to projections." << endl;
    	cout << "Resolution limits (low, high):  " << res_lo << " - " << res_hi << " A" << endl;
    	if ( res_polar ) cout << "Polar resolution limit:         " << res_polar << " A" << endl;
    	cout << "Annular range:                  " << ann_min << " - " << ann_max << " (max = " << proj->sizeX()/2 - 1 << ")" << endl;
    	cout << "Shift limit:                    " << shift_limit << " pixels" << endl;
    	if ( angle_limit )
			cout << "Angle limit:                    " << angle_limit*180.0/M_PI << " degrees" << endl;
    	if ( edge_radius ) cout << "Edge mask radius:               " << edge_radius << " pixels" << endl;
    	if ( mask_file.length() ) cout << "Real space mask file:           " << mask_file << endl;
    	cout << endl;
		cout << "Micrograph\tPID\tOriX\tOriY\tViewX\tViewY\tViewZ\tAngle\tCC\tCV\tTime" << endl;
	}
	
	// Directory for individual particle log files
	if ( part_log ) {
		mkdir("log", O_CREAT );
		chmod("log", 0755);
	}
	
	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	time_t			t = time(NULL);

	for ( i=0; i<npart; ++i ) {
		part = partarr[i];
		mg = part->mg;
		part_determine_orientation2(part, proj, part_mask,
				prs_mask, sym, bin,
				res_lo, res_hi, res_polar, ann_min, ann_max,
				shift_limit, angle_limit, edge_radius, flags,
				planf_1D, planb_1D, planf_2D, planb_2D);
		if ( verbose & VERB_RESULT ) {
			cout << setw(15) << mg->id << tab << part->id << tab << fixed <<
				setprecision(2) << part->ori[0] << tab << part->ori[1] << tab <<
				setprecision(4) << part->view[0] << tab <<
				part->view[1] << tab << part->view[2] << tab <<
				setprecision(2) << part->view.angle()*180.0/M_PI << tab <<
				setprecision(4) << part->fom[0] << tab << part->fom[1] << tab <<
				(long)(time(NULL) - t) << endl << flush;
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_determine_orientations: done, npart=" << npart << endl;

	t = (time(NULL) - t);

	long		ts = (long) t;
	long		hr = ts/3600, min = (ts - 3600*hr)/60, sec = ts - 3600*hr - 60*min;
	
	if ( verbose & VERB_RESULT && npart ) {
		cout << "Time:                           " << ts << " s (" << hr << ":" << min << ":" << sec << ")" << endl;
		cout << "Time per particle:              " << ts*1.0/npart << " s/particle" << endl;
		cout << "Algorithm time:                 " << ts*1.0e6/(npart*nproj*npix) << " us/pixel" << endl << endl;
	}

	delete[] partarr;
	delete prs_mask;
	delete part_mask;
	
    fft_destroy_plan(planf_1D);
    fft_destroy_plan(planb_1D);
    fft_destroy_plan(planf_2D);
    fft_destroy_plan(planb_2D);
	
	return 0;
}

/**
@brief 	Find the origin of every image in a project.
@param 	*project			image processing parameter structure.
@param 	*proj				reference projections.
@param 	bin					data compression by binning.
@param 	sym					point group symmetry structure.
@param 	part_select			particle selection for processing.
@param 	res_lo				low resolution limit (angstrom).
@param 	res_hi				high resolution limit (angstrom).
@param 	shift_limit			maximum shift from nominal origin of box.
@param 	flags				option flags.
@return int					error code.

	The input view is used to find the corresponding reference projection.
	This projection is then rotated by the input view angle and cross-correlated 
	with the image to find the shift.
	Flags:
		PART_LOG	write log files in log directory

**/
int 		project_determine_origins(Bproject* project, Bimage* proj, int bin,
				Bsymmetry& sym, int part_select, double res_lo, double res_hi,
				double shift_limit, int flags)
{
	int				part_log = (flags & PART_LOG)? 1: 0;

	Bfield*			field = project->field;

	Bparticle*		part = part_find_first(project);
	Bmicrograph*	mg = part->mg;

	if ( !part ) {
		cerr << "Error in project_determine_origins: No particles found!" << endl;
		return -1;
	}
	
	Bimage*			p = NULL;
	Bimage*			prs_mask = NULL;

	if ( bin < 1 ) bin = 1;
	if ( bin > 4 ) bin = 4;

	Vector3<double>	pixel_size = mg->pixel_size*bin;
	
	// Read the reference (model) projections, convert to floating point and clean up
	long 			i, nproj(0), npart(0), npix(0), vmin(0);

	nproj = proj->images();
	
	npix = proj->sizeX()*proj->sizeY();

	// Specify FFTW arrays
	fft_plan		planf_2D = proj->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb_2D = proj->fft_setup(FFTW_BACKWARD, 1);

	// Set resolution limits
	if ( res_lo < res_hi ) swap(res_hi, res_lo);
	if ( res_lo < 4*pixel_size[0] ) res_lo = 4*pixel_size[0];
	if ( res_hi < 2*pixel_size[0] ) res_hi = 2*pixel_size[0];

	// Check the shift limit
	shift_limit /= bin;
	if ( shift_limit < 0 ) shift_limit = proj->sizeX()/10;	// Default shift limit

	double			cc, angle, da, da_min;
	View*			pv = new View[nproj];
	View			hasu_view;
	
	for ( i=0; i<nproj; i++ )
			pv[i] = proj->image->view();
	
	if ( verbose & VERB_RESULT ) {
		cout << "Determining origins:" << endl;
    	cout << "Projection file name:           " << proj->file_name() << endl;
    	cout << "Number of projections:          " << nproj << endl;
		cout << "Pixel size:                     " << pixel_size << endl;
    	cout << "Resolution limits (low, high):  " << res_lo << " - " << res_hi << " A" << endl;
    	cout << "Shift limit:                    " << shift_limit << " pixels" << endl << endl;
		cout << "Micrograph\tPID\tOriX\tOriY\tViewX\tViewY\tViewZ\tAngle\tCC\tTime" << endl;
	}
	
	
	ofstream		flog;
	Bstring			log_name;
	Bimage*			proj1;
	Vector3<double>	shift;
	
	// Directory for individual particle log files
	if ( part_log ) {
		mkdir("log", O_CREAT );
		chmod("log", 0755);
	}
	
	time_t			t = time(NULL);
	
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) if ( mg->part ) {
			for ( part=mg->part; part; part=part->next )
					if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) ) {
				if ( verbose & VERB_DEBUG ) {
					cout << "DEBUG project_determine_origins: file=";
					if ( part->fpart.length() ) cout << mg->fpart;
					else cout << part->fpart;
					cout << " image=" << part->id << endl;
				}
				if ( part_log ) {
					log_name = "log/" + mg->id + "_" + Bstring(part->id, "%04d") + "_orient.log";
					flog.open(log_name.c_str());
					if ( flog.fail() )
						return error_show("Error in project_determine_origins: Log file cannot be written!", __FILE__, __LINE__);
					flog << field->id << ": " << mg->id << ": " << part->id << endl;
				}
				if ( part->fpart.length() ) p = read_img(part->fpart, 1, 0);
				else p = read_img(mg->fpart, 1, part->id - 1);
				if ( !p )
					return error_show("Error in project_determine_origins", __FILE__, __LINE__);
				p->change_type(Float);
				if ( bin > 1 ) p->bin(bin);
				p->statistics();
				p->rescale_to_avg_std(0, 1);
				p->shift_background(0);
				p->sampling(pixel_size);
				hasu_view = part->view;
				hasu_view[3] += TWOPI;
				for ( i=vmin=0, da_min = 1e10; i<nproj; i++ ) {
					da = hasu_view.angle(pv[i]);
					da = fabs(angle_set_negPI_to_PI(da));
					if ( da < da_min ) {
						da_min = da;
						vmin = i;
					}
				}
				part->sel = vmin + 1;
				proj1 = proj->extract(vmin);
				if ( !proj ) {
					error_show("Error in project_determine_origins", __FILE__, __LINE__);
					return -1;
				}
				angle = hasu_view.angle() - pv[vmin].angle();
				proj1->rotate(angle);
				shift = p->find_shift(proj1, prs_mask, res_hi, res_lo, shift_limit, 0, 1, planf_2D, planb_2D, cc);
				part->ori = (proj1->image->origin() + shift) * bin;
				part->fom[0] = cc;
				delete proj1;
				if ( part_log ) {
					flog << part->id << tab << part->ori[0] << tab << part->ori[1] << tab << cc << endl;
					flog.close();
				}
				if ( verbose & VERB_RESULT )
					cout << setw(15) << mg->id << tab << part->id << tab << 
						setprecision(2) << part->ori[0] << tab << part->ori[1] << tab << 
						setprecision(4) << part->view[0] << tab << 
						part->view[1] << tab << part->view[2] << tab << 
						setprecision(2) << part->view.angle()*180.0/M_PI << tab << 
						setprecision(4) << part->fom[0] << tab << 
						(long)(time(NULL) - t) << endl;
				delete p;
				npart++;
			} else part->fom[0] = part->sel = 0;
		}
	}
	
	t = (time(NULL) - t);
	
	long		ts = (long) t;
	long		hr = ts/3600, min = (ts - 3600*hr)/60, sec = ts - 3600*hr - 60*min;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Time:                           " << ts << " s (" << hr << ":" << min << ":" << sec << ")" << endl;
		cout << "Time per particle:              " << ts*1.0/npart << " s/particle" << endl;
		cout << "Algorithm time:                 " << ts*1.0e6/(npart*nproj*npix) << " us/pixel" << endl << endl;
	}

    fft_destroy_plan(planf_2D);
    fft_destroy_plan(planb_2D);
	
	delete[] pv;
	
	return 0;
}


/**
@brief 	Rotates and shifts a reference image and calculates a cross-validation coefficient.
@param 	*p			2D image.
@param 	*pref		reference 2D image.
@param 	*prs_mask	dual mask.
@param 	planf		FFT forward plan.
@return double		cross-validation coefficient.

	The reference image is first rotated by the difference in the angles
	between the image and the reference, then shifted to the same origin
	as the image. Both are then Fourier transformed, the negative values in
	the mask used to select reciprocal space areas, and the complex product
	calculated. The sum of the zeroeth pixel then gives the cross-validation
	coefficient.
	The input images must be equal-sized square 2D images.

**/
double		img_cross_validate(Bimage* p, Bimage* pref, Bimage* prs_mask, fft_plan planf)
{
	Bimage*			prot = pref->rotate(pref->size(), fabs(p->image->view_angle() - pref->image->view_angle()));

	Bimage*			pcc = prot->pack_two_in_complex(p);

	delete prot;
	
	pcc->fft(planf, 0);
	pcc->complex_apply_negative_mask(prs_mask);
	pcc->combined_complex_product();
	
	double			cv(0);
	long			i, datasize = pcc->sizeX()*pcc->sizeY()*pcc->sizeZ()*pcc->images();

	for ( i=0; i<datasize; i++ )
		if ( (*prs_mask)[i] < 0 )
			cv += (pcc->complex(i)).real();

	delete pcc;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_cross_validate: cv=" << cv << endl;
	
	return cv;
}

