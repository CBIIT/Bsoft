/**
@file	mg_tomo_rec.cpp
@brief	Functions to do a tomographic reconstruction
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20181221
**/

#include "mg_tomo_rec.h"
#include "mg_reconstruct.h"
#include "mg_tomography.h"
#include "mg_select.h"
#include "mg_ctf.h"
#include "Complex.h"
#include "linked_list.h"
#include "utilities.h"
#include "timer.h"

#include <fstream>

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int			img_backtransform_z_lines(fstream* ftemp, Bimage* p);
int			img_write_data_block_with_type(fstream* ftemp, Bimage* p, 
				double avg, double std, double cutmin, double cutmax);

Bimage*		mg_tomo_rec_prepare(Bmicrograph* mg, int ft_size, Vector3<long> rec_size,
				double edge_width, double marker_radius,
				int fill_type, double fill,
				int action, double wiener, Vector3<long> tile_size, fft_plan plan)
{
	if ( verbose & VERB_FULL )
		cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
	
//	Bimage*			p = mg_erase_markers(mg, marker_radius);

	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);
	if ( !p ) {
		error_show("mg_erase_markers", __FILE__, __LINE__);
		return NULL;
	}
	
	p->change_type(Float);
	
	p->sampling(mg->pixel_size);
	
	p->origin(mg->origin[0], mg->origin[1], 0);

	Vector3<long>	size(rec_size[0], rec_size[1], 1);
	Vector3<long>	translate = (size - p->size())/2;
//	cout << size << tab << translate << endl;
	
	p->resize(size, translate, fill_type, fill);

	p->statistics();
	p->rescale_to_avg_std(0,1);

//	Vector3<long> 	tile_size(128,128,1);
//	Vector3<long> 	tile_size(64,64,1);
	if ( action ) {
		if ( verbose & VERB_PROCESS )
			cout << "Correcting for the CTF" << endl;
		if ( !mg->ctf ) {
			cerr << "Error: The CTF parameters are not specified! Abort!" << endl;
			bexit(-1);
		}
		img_ttf_apply(p, *(mg->ctf), action, wiener,
				tile_size, mg->tilt_angle, mg->tilt_axis, 0, 0);
	}
	
	if ( edge_width )
//		micrograph_clear_extraneous_areas(mg, p, rec_size[2], edge_width);
		img_clear_extraneous_areas(p, mg->tilt_axis, mg->tilt_angle, rec_size[2], edge_width);

	if ( marker_radius )
		img_erase_markers(p, mg->mark, marker_radius);
	
//	p->calculate_background();
//	p->statistics();
//	p->rescale_to_avg_std(0,1);
	
//	p->pad(ft_size, fill_type, fill);
	p->pad(ft_size);

	p->fft(plan, 1);
	
	p->phase_shift_to_origin();
	
	return p;
}


/**
@brief 	Reciprocal space reconstruction from the images in a multi-image file.  
@param 	*project 		image processing parameter structure.
@param 	hi_res			high resolution limit.
@param 	scale			scale of reconstruction.
@param 	size			size of reconstruction.
@param	interp_type		interpolation type.
@param 	pad_factor		factor that determines image padding.
@param 	edge_width		edge smoothing width for masks.
@param 	marker_radius	flag and radius to mask out markers.
@param 	fill_type		FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill			value to paint markers.
@param 	action			flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@return	Bimage*			reconstruction, NULL on failure.

	The orientation parameters, view vector, angle of rotation and origin,
	must all be set. Each image is padded to at least two times its size 
	and its Fourier transform packed into 3D reciprocal space.
	The figure-of-merit calculated for each reciprocal space voxel is:
		       sum(w*re)^2 + sum(w*im)^2
		FOM = ---------------------------
		      sum(w)*sum(w*(re^2 + im^2))
	where
		re	real part
		im	imaginary part
		w	weight (inverse distance of image pixel to closest grid point)
	For voxels with only one data pixel contributing to it, FOM = 0.
	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.

**/
Bimage*		project_tomo_reconstruct(Bproject* project, double hi_res,
				double scale, Vector3<long> size, int interp_type, int pad_factor,
				double edge_width, double marker_radius, 
				int fill_type, double fill, int action, double wiener)
{
	double			ti = timer_start();

	int				nmg = project_count_mg_selected(project);
	Bfield*			field;
	Bmicrograph*	mg = NULL;

	if ( nmg < 1 ) {
		cerr << "Error: No micrographs selected!" << endl;
		return NULL;
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select )
			if ( mg->fmg.length() ) break;
		if ( mg ) break;
	}
	
	if ( !mg ) {
		cerr << "No selected micrograph found!" << endl << endl;
		return NULL;
	}
	
	if ( !mg->fmg.length() ) {
		cerr << "No file name for micrograph " << mg->id << endl << endl;
		return NULL;
	}

	Bimage*			p = read_img(mg->fmg, 0, mg->img_num);
	
	if ( size.volume() < 1 )
		size = Vector3<long>(p->sizeX()*scale, p->sizeY()*scale, p->sizeX()*scale/10);

	int				ft_size = part_ft_size(p->sizeX(), scale, pad_factor);

	delete p;

	if ( hi_res < 2*mg->pixel_size[0]/scale )
		hi_res = 2*mg->pixel_size[0]/scale;
	double			rec_scale = 2*mg->pixel_size[0]/(scale*hi_res);
	if ( rec_scale > 1 ) rec_scale = 1;
	double			voxel_size(mg->pixel_size[0]/(scale*rec_scale));
	Vector3<double>	vscale(rec_scale*scale, rec_scale*scale, rec_scale*scale);
	Vector3<long>	rec_size = size*rec_scale;
	long   			ds = rec_size.volume();

	long			memreq1 = 2*size.volume()*sizeof(float);
	long			memreq2 = (ft_size*ft_size + 5*rec_size.volume())*sizeof(float);

//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "3D reciprocal space reconstruction:" << endl;
		cout << "Reconstruction size:            " << size << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Integration size:               " << rec_size << endl;
		cout << "Integration scale:              " << rec_scale << endl;
		cout << "Voxel size:                     " << voxel_size << " A" << endl;
		cout << "Resolution:                     " << hi_res << " A" << endl;
		cout << "Interpolation type:             " << interp_type << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Padding factor:                 " << pad_factor << endl;
		cout << "Edge smooting width:            " << edge_width << endl;
		if ( marker_radius )
			cout << "Erase markers using radius:     " << marker_radius << endl;
		if ( action )
			cout << "CTF correction:                 " << action << " (wiener=" << wiener << ")" << endl;
		cout << endl;
	}
	
	if ( memreq1 > memreq2) memory_check(memreq1);
	else memory_check(memreq2);

	Vector3<long> 	tile_size(128,128,1);
	fft_plan		plan = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 0);

	// Header structure for the final reconstruction
	Bimage* 		prec = new Bimage(Float, TComplex, rec_size, 1);
	prec->fourier_type(Standard);
	prec->sampling(voxel_size, voxel_size, voxel_size);
	prec->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	prec->next->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	prec->next->next->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	float*			fom = (float *) prec->next->data_pointer();
	float*			weight = (float *) prec->next->next->data_pointer();
	float*			weight2 = (float *) prec->next->next->next->data_pointer();

	double			pad_ratio(ft_size*1.0/prec->sizeX());
//	if ( pad_ratio > 1 ) pad_ratio *= pad_ratio*0.6;
//	if ( pad_ratio > 1 ) pad_ratio *= 0.25;
	
	if ( verbose )
		cout << "Transforming and packing " << nmg << " micrographs:" << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			p = mg_tomo_rec_prepare(mg, ft_size, size/scale, edge_width,
					marker_radius, fill_type, fill, action, wiener, tile_size, plan);
			if ( p ) {
				if ( verbose )
					cout << "Packing " << p->file_name() << " (" << mg->img_num << ")" << endl;
				prec->fspace_pack_2D(p, mg->matrix, hi_res, 0, vscale, 1/pad_ratio, interp_type);
//				prec->fspace_pack_2D(p, mg->matrix, hi_res, 0, vscale, 1, interp_type);
				delete p;
			}
		}
	}

	fft_destroy_plan(plan);
	
	timer_report(ti);

	if ( verbose )
		cout << "Weighing reconstruction" << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(ds, dispatch_get_global_queue(0, 0), ^(size_t i){
//		weight[i] += 1;
//		fom[i] += 1;
		if ( weight[i] > SMALLFLOAT ) {
			weight2[i] = weight[i] - weight2[i]/weight[i];
			prec->set(i, prec->complex(i) / weight[i]);
			if ( weight2[i] > 1e-3 ) fom[i] /= weight2[i];
			else if ( weight[i] > 1e-3 ) fom[i] /= weight[i];
//			fom[i] /= weight[i];
		} else fom[i] = 0;
		if ( fom[i] < 0 ) fom[i] = 0;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<ds; i++ ) {
		if ( weight[i] > SMALLFLOAT ) {
			weight2[i] = weight[i] - weight2[i]/weight[i];
			prec->set(i, prec->complex(i) / weight[i]);
			if ( weight2[i] > 1e-3 ) fom[i] /= weight2[i];
			else if ( weight[i] > 1e-3 ) fom[i] /= weight[i];
//			fom[i] /= weight[i];
		} else fom[i] = 0;
		if ( fom[i] < 0 ) fom[i] = 0;
	}
#endif

	long 			cov(0);
	for ( long i=0; i<ds; i++ ) if ( weight[i] > SMALLFLOAT ) cov++;

//	prec->fspace_reconstruction_stats(hi_res, 4);
	
	prec->image->origin(prec->size()/2);
	prec->phase_shift_to_origin();
	prec->image->origin(prec->size()/2);
	
	if ( verbose & VERB_RESULT )
		cout << "Coverage:                       " << cov << " (" << cov*100.0/ds << " %)" << endl << endl;
	
	delete prec->next;
	prec->next = NULL;

	timer_report(ti);

	return prec;
}

/**
@brief 	Reciprocal space reconstruction from the images in a multi-image file.  
@param 	*project 		image processing parameter structure.
@param 	hi_res			high resolution limit.
@param 	scale			scale of reconstruction.
@param 	size			size of reconstruction.
@param 	slab_start		start of reconstruction slab.
@param 	slab_end		end of reconstruction slab.
@param 	marker_radius	flag and radius to mask out markers.
@param 	fill_type		FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill			value to paint markers.
@param 	action			flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@return	Bimage*			reconstruction, NULL on failure.

	The orientation parameters, view vector, angle of rotation and origin,
	must all be set. Each image is padded to at least two times its size 
	and its Fourier transform packed into 3D reciprocal space.
	The figure-of-merit calculated for each reciprocal space voxel is:
		       sum(w*re)^2 + sum(w*im)^2
		FOM = ---------------------------
		      sum(w)*sum(w*(re^2 + im^2))
	where
		re	real part
		im	imaginary part
		w	weight (inverse distance of image pixel to closest grid point)
	For voxels with only one data pixel contributing to it, FOM = 0.
	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.

**/
Bimage*		project_fourier_reconstruction_slab(Bproject* project, double hi_res, 
				double scale, Vector3<long> size, int slab_start, int slab_end,
				double marker_radius, int fill_type, double fill, int action, double wiener) 
{
	double			ti;
//	if ( verbose & VERB_TIME )
		ti = timer_start();
		
	long			i;

	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
	Bimage* 		p = NULL;

	long			nmg = project_count_mg_selected(project);
	long			nfld = count_list((char *) project->field);
	long 			nrec(0);
	
	if ( nmg < 1 ) {
		cerr << "Error: No micrographs selected!" << endl;
		return NULL;
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
//			if ( mg->fft.length() ) break;
			if ( mg->fmg.length() ) break;
		}
		if ( mg ) break;
	}
	
	if ( !mg ) {
		cerr << "No selected micrograph with a micrograph file found!" << endl << endl;
		return NULL;
	}
	
	if ( !mg->fmg.length() ) {
		cerr << "No file name for micrograph " << mg->id << endl << endl;
		return NULL;
	}
	
	if ( mg->matrix.determinant() < 0.9 )
		project_mg_tilt_to_matrix(project);
	
	if ( mg->fmg.length() )
		p = read_img(mg->fmg, 0, 0);
	else {
		error_show("Error in project_fourier_reconstruction_slab ", __FILE__, __LINE__);
		cerr << "No micrograph image or transform file name given for micrograph " << mg->id << endl << endl;
		return NULL;
	}
	
	if ( size.volume() < 1 )
		size = {long(p->sizeX()*scale), long(p->sizeY()*scale), long(p->sizeX()*scale/10)};
		
/*	if ( size.volume() <= 0 ) {
		size[0] = (int) (p->sizeX()*scale);
		if ( size[0] < p->sizeY()*scale ) size[0] = (int) (p->sizeY()*scale);
		if ( size[0] < p->sizeZ()*scale ) size[0] = (int) (p->sizeZ()*scale);
		size[2] = size[1] = size[0];
	}
*/	
	if ( slab_start >= size[2] ) {
		cerr << "Error: Slab " << slab_start << " - " << slab_end << " is outside the volume limits!" << endl;
		return NULL;
	}
	
	if ( slab_start < 0 ) slab_start = 0;
	if ( slab_end < 0 || slab_end >= size[2] ) slab_end = size[2] - 1;
	if ( slab_end < slab_start ) slab_end = slab_start;
	
	long			slab_thickness = slab_end - slab_start + 1;

	int				pad = 2;
	long 			ft_size = findNextPowerOf((int)(pad*size[0]), 2);
	if ( ft_size < p->sizeX() ) ft_size = p->sizeX();

	delete p;
	
	// Header structure for the final reconstruction
	Bimage* 		prec = new Bimage(Float, TComplex, size[0], size[1], slab_thickness, 1);
	prec->sampling(mg->pixel_size[0]/scale, mg->pixel_size[0]/scale, mg->pixel_size[0]/scale);
	prec->origin(prec->sizeX()/2, prec->sizeY()/2, size[2]/2);
	prec->check_resolution(hi_res);
	prec->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	float*			fom = (float *) prec->next->data_pointer();

	long   			ds = prec->sizeX()*prec->sizeY()*prec->sizeZ();
	float*			weight = new float[ds];	
	float*			weight2 = new float[ds];	
	
//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "3D reciprocal space reconstruction on disk:" << endl;
		cout << "Map size:                       " << prec->sizeX() << " " << prec->sizeY() << " " << size[2] << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Slab:                           " << slab_start << " - " << slab_end << endl << endl;
	}


	Vector3<long> 	tile_size(128,128,1);
	fft_plan		plan = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 0);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			p = mg_tomo_rec_prepare(mg, ft_size, size/scale, 5.0,
					marker_radius, fill_type, fill, action, wiener, tile_size, plan);
			if ( !p ) {
				error_show("project_fourier_reconstruction_slab", __FILE__, __LINE__);
				return NULL;
			}
	
			img_pack_2D_in_recip_space_slab(p, prec, size[2], slab_start, weight, weight2,
					hi_res, mg->matrix, scale);
				
			delete p;
					
			nrec++;
					
			if ( verbose & ( VERB_TIME | VERB_PROCESS ) ) {
				cout << "Complete:                       " 
							<< setprecision(3) << nrec*100.0/(nfld*nmg) << " %\r";
				cout.flush();
			}
		}
	}

	fft_destroy_plan(plan);
	
	if ( nrec < 1 ) {
		delete prec;
		error_show("Error in project_fourier_reconstruction_slab: No micrographs used in this reconstruction!", __FILE__, __LINE__);
		return NULL;
	}

	if ( verbose & VERB_TIME )
		timer_report(ti);

	if ( verbose )
		cout << "Weighing reconstruction." << endl;
		
	long 	cov(0);

#ifdef HAVE_GCD
	dispatch_apply(ds, dispatch_get_global_queue(0, 0), ^(size_t i){
		if ( weight[i] > SMALLFLOAT ) {
			prec->set(i, prec->complex(i) / weight[i]);
			weight2[i] = weight[i]*weight[i] - weight2[i];	// Unbiased weighting
			if ( weight2[i] > 0.001 ) {
				double	rd2 = (prec->complex(i)).power();
				fom[i] = (fom[i] - weight[i]*rd2)/weight2[i];	// Noise variance estimator with noise reduction
				fom[i] = rd2/fom[i] - 1;	// SNR
			} else fom[i] = 0;
		}
	});
#else
#pragma omp parallel for
	for ( i=0; i<ds; i++ ) {
		if ( weight[i] > SMALLFLOAT ) {
			prec->set(i, prec->complex(i) / weight[i]);
			weight2[i] = weight[i]*weight[i] - weight2[i];	// Unbiased weighting
			if ( weight2[i] > 0.001 ) {
				double	rd2 = (prec->complex(i)).power();
				fom[i] = (fom[i] - weight[i]*rd2)/weight2[i];	// Noise variance estimator with noise reduction
				fom[i] = rd2/fom[i] - 1;	// SNR
			} else fom[i] = 0;
		}
	}
#endif

	for ( i=0; i<ds; i++ ) if ( weight[i] > SMALLFLOAT ) cov++;
	
	delete[] weight;
	delete[] weight2;

	img_phase_shift_slab_to_origin(prec, size[2], slab_start);
	
	long	vol = ds;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Coverage:" << endl;
		cout << "Total number of voxels:         " << vol << endl;
		cout << "Total coverage:                 " << cov << endl;
		cout << "Percentage coverage:            " << cov*100.0/vol << " %" << endl;
	}

//	if ( verbose & VERB_TIME )
		timer_report(ti);

	return prec;
}



/**
@brief 	Reconstructs individual particles from a tilt series.
@param 	*project		micrograph project.
@param 	*recpart		3D particle within the project.
@param	recsize			particle reconstructions ize
@param 	resolution		high resolution limit for reconstruction.
@param	interp_type		interpolation type.
@param 	ft_size			2D Fourier transform size.
@param	scale			reconstruction scale.
@param	planp			2D Fourier transform plan.
@param	planr			3D Fourier transform plan.
@param 	ctf_action		flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	*sym			always C1.
@param 	&partbase		particle base name for new particle reconstructions.
@param 	&partpath		directory for new particle reconstructions.
@param 	&partext		extension of new reconstructions.
@return long			number of particles.

	Requires the particles to be defined in all micrographs.
	The partbase, partpath and partext arguments can be left empty to
	use defaults.

**/
long		particle_tomo_reconstruct(Bproject* project, Bparticle* recpart,
				Vector3<long> recsize, double resolution, int interp_type,
				long ft_size, Vector3<double> scale, fft_plan planp, fft_plan planr, 
				int ctf_action, double wiener, Bsymmetry& sym,
				Bstring& partbase, Bstring& partpath, Bstring& partext)
{
	Bstring				filename("part");
	
	// Setting the base file name
	if ( partbase.length() ) filename = partbase + Bstring(recpart->id, "_part%03d.");
	else if ( recpart->fpart.length() ) filename = recpart->fpart.pre_rev('.') + "_r.";
//	else filename = rec->frec.pre_rev('.') + Bstring(recpart->id, "_part%03d.");
		
	// Setting the extension
	if ( partext.length() ) filename += partext;
	else if ( recpart->fpart.length() ) filename += recpart->fpart.post_rev('.');
//	else filename += rec->frec.post_rev('.');
		
	// Setting the path
	if ( partpath.length() ) {
		if ( filename.contains("/") ) filename = filename.post_rev('/');
		filename = partpath + filename;
	}
		
	if ( verbose & VERB_FULL )
		cout << "Reconstructing particle " << recpart->id << endl;

	Bfield*				field;
	Bmicrograph*		mg;
	Bparticle*			part;
	Bparticle*			partlist = NULL;
		
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			for ( part = mg->part; part; part = part->next )
				if ( recpart->id == part->id ) break;
			if ( part ) {
				particle_copy(&partlist, part);
			}
		}
	}

	Bimage*				prec = particle_reconstruct(partlist, sym, 0,
							resolution, scale, recsize, ft_size, planp, 
							interp_type, ctf_action, wiener, 0, 0);
	
	long				nr = prec->image->FOM();
	
	prec->fspace_reconstruction_weigh();
		
	prec->phase_shift_to_center();
	prec->origin(prec->size()/2);
	prec->fft_back(planr);
	prec->statistics();
	prec->correct_background();
//	prec->change_type(newdatatype);
	write_img(filename, prec, 0);
	recpart->fpart = filename;
	recpart->view = View(0,0,1,0);
	recpart->ori = prec->image->origin();
		
	delete prec;
	particle_kill(partlist);
		
	return nr;
}

/**
@brief 	Reconstructs particles from a tilt series.
@param 	*project		micrograph project.
@param 	resolution		high resolution limit for reconstruction.
@param	interp_type		interpolation type.
@param 	pad_factor		factor that determines image padding.
@param 	ctf_action		flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	*sym			point group symmetry.
@param 	&partbase		particle base name for new particle reconstructions.
@param 	&partpath		directory for new particle reconstructions.
@param 	&partext		extension of new reconstructions.
@return long			number of particles.

	Requires the particles to be defined in all micrographs.
	The partbase, partpath and partext arguments can be left empty to
	use defaults.

**/
long		project_tomo_reconstruct_particles(Bproject* project, 
				double resolution, int interp_type, int pad_factor, 
				int ctf_action, double wiener, Bsymmetry& sym,
				Bstring& partbase, Bstring& partpath, Bstring& partext)
{
	if ( !partpath.empty() ) {
		if ( partpath[-1] != '/' ) partpath += "/";
		mkdir(partpath.c_str(), (mode_t)0755);
	}

	if ( partext.length() < 1 ) partext = "pif";
	
	
	Bfield*				field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			partlist = NULL;
	Bparticle*			part = NULL;

	for ( rec = project->rec; rec; rec = rec->next )
		if ( rec->part ) {
			partlist = rec->part;
			break;
		}
	
	if ( !partlist ) {
		cerr << "No particle list found!" << endl << endl;
		return -1;
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mg ) break;
		if ( mg ) break;
	}
	
	if ( !mg ) {
		cerr << "No micrographs found!" << endl << endl;
		return -1;
	}
	
	if ( !mg->part ) {
		cerr << "No micrograph particles found!" << endl << endl;
		return -1;
	}
	
//	Vector3<long>		boxsize(mg->box_size);
	
	rec->box_size = mg->box_size.max(mg->box_size[0]);
//	Vector3<long>		recsize(boxsize);
//	recsize[2] = recsize[1] = recsize[0];

	rec->voxel_size = mg->pixel_size;
	rec->voxel_size[1] = rec->voxel_size[2] = rec->voxel_size[0];
	
	Vector3<double>		scale(rec->box_size[0]*1.0/mg->box_size[0],
							rec->box_size[1]*1.0/mg->box_size[1],
							rec->box_size[0]*1.0/mg->box_size[0]);
	
	long				ft_size = part_ft_size(mg->box_size[0], scale[0], pad_factor);
	
	fft_plan			planp = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	fft_plan			planr = fft_setup_plan(rec->box_size, FFTW_BACKWARD, 1);

	if ( verbose ) {
		cout << "Calculating individual particle reconstructions:" << endl;
		cout << "Reconstruction:                 " << rec->id << endl;
		cout << "Particle size:                  " << mg->box_size << endl;
		cout << "Reconstruction size:            " << rec->box_size << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Padding factor:                 " << pad_factor << endl;
		cout << "Reconstruction resolution:      " << resolution << " A" << endl;
		cout << "CTF application type:           " << ctf_action << endl;
		cout << endl;
	}

	long				i, npart = particle_count(partlist);
	
	Bparticle**			partarr = new Bparticle*[npart];
	
	for ( i=0, part = partlist; part; part = part->next, i++ ) partarr[i] = part;

	if ( verbose & VERB_RESULT )
		cout << "Part\t#" << endl;
	
#ifdef HAVE_GCD
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(npart, dispatch_get_global_queue(0, 0), ^(size_t i){
		long	nr = particle_tomo_reconstruct(project, partarr[i], rec->box_size, 
						resolution, interp_type, ft_size, scale, planp, planr, 
						ctf_action, wiener, sym, partbase, partpath, partext);
		dispatch_sync(myq, ^{
			if ( verbose & VERB_RESULT )
				cout << partarr[i]->id << tab << fixed << nr << endl << flush;
		});
	});
#else
#pragma omp parallel for
	for ( i=0; i<npart; i++ ) {
		long	nr = particle_tomo_reconstruct(project, partarr[i], rec->box_size, 
						resolution, interp_type, ft_size, scale, planp, planr, 
						ctf_action, wiener, sym, partbase, partpath, partext);
	#pragma omp critical
		{
			if ( verbose & VERB_RESULT ) {
				cout << partarr[i]->id << tab << fixed << nr << endl << flush;
			}
		}
	}
#endif

	delete[] partarr;
	
	fft_destroy_plan(planp);
	fft_destroy_plan(planr);
	
	return 0;
}


/**
@brief 	Backtransforms 2D slices in a 3D volume.  
@param 	*p				3D complex volume.
@return	int						0.

	Each slice is extracted, backtransformed and copied back into the
	original volume.
	Note: The phases are now not hermitian any more.

**/
int			img_backtransform_slices(Bimage* p)
{
	double			ti;
//	if ( verbose & VERB_TIME )
		ti = timer_start();
		
	long			i, j, z;
	long			slice_size = p->sizeX()*p->sizeY()*p->channels();
	Vector3<double>	extorigin;
	Vector3<long>	extsize(p->sizeX(), p->sizeY(), 1);
	Bimage*			ps = NULL;
	
	if ( verbose )
		cout << "Backtransforming reconstruction slices." << endl << endl;
	
	for ( z=0; z<p->sizeZ(); z++ ) {
		if ( verbose & VERB_FULL )
			cout << "Shifting and backtransforming slice " << z << endl;
		extorigin[2] = z;
		ps = p->extract(0, extorigin, extsize);
		ps->fft(FFTW_BACKWARD, 1);
		for ( i=z*slice_size, j=0; j<slice_size; i++, j++ ) p->set(i, (*ps)[j]);
		delete ps;
	}

//	if ( verbose & VERB_TIME )
		timer_report(ti);

	return 0;
}

/**
@brief 	Packs a 2D Fourier transform into a 3D reciprocal space volume.  
@param 	*p				2D Fourier transform
@param 	*prec			3D reciprocal space slab volume.
@param 	zsize			intentional z dimension of reconstruction.
@param 	slab_start		start of current slab.
@param 	*weight			weigth array.
@param 	*weight2		weight squared array.
@param 	hi_res			high resolution limit.
@param 	mat				affine matrix.
@param 	scale			scale of reconstruction and particle magnification.
@return	long			0.

	The rotation matrix is used to determine the plane in reciprocal space
	to which the 2D transform data is added. The map is assumed to be cubic
	and the 2D transform square. The real space 2D image must be supplied.
	This is then padded to more than twice its original size, fourier
	transformed, and packed into the 3D reciprocal space block.

**/
long		img_pack_2D_in_recip_space_slab(Bimage* p, Bimage* prec,
				long zsize, long slab_start, float* weight, float* weight2,
				double hi_res, Matrix3 mat, double scale)
{
	if ( !p->data_pointer() ) return 0;
	
	prec->check_resolution(hi_res);
	
	long 			i, j, x, y;
	long 			ix, iy, iz, izs, zs;
	double 			w, d2;
	Vector3<double>	m, d, iv;
	
	int				pad = 2;
	long 	ft_size = p->sizeX();
	Bimage*			pt;
	
	// Do transformation if not already done
	if ( p->compound_type() == TSimple ) {
		ft_size = findNextPowerOf((int)(pad*prec->sizeX()), 2);
		if ( ft_size < p->sizeX() ) ft_size = p->sizeX();
		pt = p->pad_copy(ft_size, FILL_AVERAGE, p->average());
		pt->fft();
		pt->phase_shift_to_origin();
	} else pt = p->copy();
	
	Vector3<double> 	vscale(prec->sizeX(), prec->sizeY(), zsize);
	vscale /= scale*ft_size;
	//	mat = mat.transpose();	Handedness
	mat = vscale * mat;	// Matrix includes scaling going from image to map
	
	if ( verbose )
		cout << "Packing " << p->file_name() << endl;
	
	if ( verbose & VERB_FULL ) {
		cout << "Resolution limit:             " << hi_res << " A" << endl;
		cout << "Slab start:                   " << slab_start << endl;
		cout << "Fourier transform size:       " << ft_size << endl;
		cout << "Transformation matrix:" << endl;
		cout << mat << endl;
		cout << endl;
	}
	
	float*			fom = (float *) prec->next->data_pointer();
	
	double			max_rad = pt->sampling(0)[0]/hi_res;
	if ( max_rad > 0.5 ) max_rad = 0.5;
	max_rad *= pt->sizeX();
	double			max_rad_sq = max_rad*max_rad;	// Maximum radius squared on ft scale
	double			xmin = floor(-max_rad);
	double			xmax = -xmin;
	double			ymin = floor(-max_rad);
	double			ymax = -ymin;
	
//	cout << "sam=" << pt->sampling(0)[0] << " hi_res=" << hi_res << " max_rad=" << max_rad << endl;

	for ( iv[1]=ymin; iv[1]<=ymax; iv[1]+=1 ) {
		y = (long) iv[1];
		if ( y < 0 ) y += (long)p->sizeY();
		for ( iv[0]=xmin; iv[0]<=xmax; iv[0]+=1 ) {
			x = (long) iv[0];
			if ( x < 0 ) x += (long)p->sizeX();
			d2 = iv.length2();
			if ( d2 < max_rad_sq ) {   // Nearest neighbour
				m = mat * iv;
				iz = izs = (long) floor(m[2] + 0.5);   // Nearest neighbour
				if ( izs < 0 ) izs += zsize;
				zs = izs - slab_start;
				if ( zs >= 0 && zs < prec->sizeZ() ) {
					ix = (long) floor(m[0] + 0.5);
					iy = (long) floor(m[1] + 0.5);
					d[0] = ix - m[0];
					d[1] = iy - m[1];
					d[2] = iz - m[2];
//					d /= vscale;
					w = 1 - d.length();
					if ( w > 1e-6 ) {
						if ( ix < 0 ) ix += prec->sizeX();
						if ( iy < 0 ) iy += prec->sizeY();
						i = y*pt->sizeX() + x;
						j = (zs*prec->sizeY() + iy)*prec->sizeX() + ix;
//						prec->set(j, prec->complex(j) + pt->complex(i) * w);
						prec->add(j, pt->complex(i) * w);
						fom[j] += w*(pt->complex(i)).power();
						weight[j] += w;
						weight2[j] += w*w;
					}
				}
			}
		}
	}
	
	delete pt;

	return 0;
}

/**
@brief 	Phase shifts a set of reflections to the image origin.
@param 	*p					complex image.
@param 	zsize				slab thickness.
@param 	slab_start			slab start.
@return int					0.

	A real space translation with wrapping is equivalent to phase shifting
	in reciprocal space. The phases are shifted based on the embedded
	sub-image origins.

**/
int			img_phase_shift_slab_to_origin(Bimage* p, int zsize, int slab_start)
{
	if ( p->compound_type() != TComplex ) return -1;
	
	long   i, x, y, z, zs;
	long 			h, k, l;
	double			skl, sl, phi;
	Vector3<double>	shift((p->sizeX()/2)*1.0/p->sizeX(), (p->sizeY()/2)*1.0/p->sizeY(), (zsize/2)*1.0/zsize);
	Vector3<double>	half((p->sizeX() - 1)/2, (p->sizeY() - 1)/2, (zsize - 1)/2);
	Complex<double>	temp;

	if ( verbose & VERB_PROCESS )
		cout << "Translate slice within unit cell to phase origin by " << shift << endl << endl;
	
	for ( i=zs=0, z=slab_start; zs<p->sizeZ(); zs++, z++ ) {
		l = z;
		if ( l > half[2] ) l -= zsize;
		sl = l*shift[2];
		for ( y=0; y<p->sizeY(); y++ ) {
			k = y;
			if ( k > half[1] ) k -= p->sizeY();
			skl = sl + k*shift[1];
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				h = x;
				if ( h > half[0] ) h -= p->sizeX();
				phi = MIN2PI*(h*shift[0] + skl);
				temp = Complex<double>(cos(phi), sin(phi));
				p->set(i, (p->complex(i))*temp);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Fourier transform micrographs and write to disk.  
@param 	*project 		image processing parameter structure.
@param 	size			intended reconstruction size.
@param 	scale			reconstruction scale.
@param 	pad_factor		factor that determines image padding.
@param 	datatype		datatype (default complex float).
@param 	marker_radius	flag to mask out markers.
@param 	fill_type		FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill			value to paint markers.
@return	Bimage*			reconstruction, NULL on failure.

	Each micrograph is padded to a square size that has power of 2
	dimensions. The micrograph is transformed and the phases shifted
	to the origin.
	A pad factor of zero indicates use of original size.

**/
int			mg_fft_write(Bproject* project, Vector3<int> size, double scale, int pad_factor,
				DataType datatype, double marker_radius, int fill_type, double fill)
{
	if ( pad_factor < 0 ) pad_factor = 0;	// 0 indicates using original size
	if ( pad_factor > 8 ) pad_factor = 8;
	
	if ( datatype < Short ) datatype = Float;
	
	int				n, nsel(0);
	
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bmarker*		mark = NULL;
	Bimage* 		p = read_img(mg->fmg, 0, 0);
	
	if ( !mg->fmg.length() ) {
		error_show("Error in mg_fft_write", __FILE__, __LINE__);
		cerr << "No micrograph image file name given for micrograph " << mg->id << endl << endl;
		return -1;
	}
	
	if ( size[2] < 1 ) size[2] = p->sizeX();
	if ( size.volume() < 1 ) size = Vector3<int>(p->size());
	if ( size[2] > p->sizeX() ) size[2] = p->sizeX();
	if ( size[2] > p->sizeY() ) size[2] = p->sizeY();
	
	long 	ft_size = findNextPowerOf((int)(pad_factor*size[0]*scale), 2);
	if ( ft_size < p->sizeX() ) ft_size = p->sizeX();
	if ( ft_size < p->sizeY() ) ft_size = p->sizeY();

	delete p;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Transforming micrographs:" << endl;
		cout << "Reconstruction size:            " << size << endl;
		cout << "Reconstruction scale:           " << scale << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl << endl;
	}

	for ( n=0, mg=field->mg; mg; mg=mg->next, n++ ) if ( mg->select ) {
		if ( verbose )
			cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
		if ( ( p = read_img(mg->fmg, 1, mg->img_num) ) == NULL ) {
			error_show("mg_fft_write", __FILE__, __LINE__);
			return -1;
		}
	
		p->image->origin(mg->origin);
		
		p->calculate_background();
		if ( size[2] < p->sizeX() && size[2] < p->sizeY() )
			img_clear_extraneous_areas(p, mg->tilt_axis, mg->tilt_angle, size[2], 5);
		
		if ( marker_radius > 0 ) {
			if ( marker_radius < mg->mark_radius ) marker_radius = 1.5*mg->mark_radius;
			p->statistics();
			for ( mark = mg->mark; mark; mark = mark->next )
				p->sphere(mark->loc, marker_radius, 2, fill_type, fill);
		}
		
		p->pad(ft_size, fill_type, fill);
		
		p->fft();
				
		p->phase_shift_to_origin();
		
		p->change_type(datatype);
		
		mg->fft = mg->fmg;
		if ( mg->fft.contains("/") ) mg->fft = mg->fft.post_rev('/');
		mg->fft = mg->fft.pre_rev('.') + "_ft" + Bstring(n, "%04d") + ".sup";
		write_img(mg->fft, p, 0);
//		mg->img_num = 0;
				
		delete p;
				
		nsel++;
	}
	
	return nsel;
}

/**
@brief 	Extracts a tile from the image with limits in y.  
@param 	*file_list	list of image file names.
@param 	ystart			first y index.
@param 	ysize			size in y.
@return  	Bimage*				tile, NULL on failure.

	From a series of z-slab images, tiles are extracted from a defined start in y,
	and with a defined size in y.
	A file indicating the origins of the tiles, "y.tiles", is written to
	be used with bpatch to assemble the tiles.

**/
Bimage*		img_extract_ytile(Bstring* file_list, int ystart, int ysize)
{
	Bimage*				p = NULL;
	Bimage*				pmap = NULL;
	int					err(0);
	long		i(0), j(0), y, z, yt, zt;
	Bstring*			thisfile = NULL;

	if ( verbose )
	
	// First read the file headers to identify problems
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		p = read_img(*thisfile, 0, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_DEBUG )
				cout << p->file_name() << ": " << 
						p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << " " << p->images() << " " << p->channels() << endl;
			if ( p->images() > 1 ) {
				cerr << p->file_name() << " must be a single image file!" << endl;
				err += -1;
			}
			if ( i == 0 ) {
				pmap = p;
			} else {
//				err += img_compatibility(pmap, p);
				if ( !p->compatible(pmap) ) err += -1;
				delete p;
			}
			i += p->sizeZ();
		} else err += -1;
	}
	
	if ( err < 0 || i == 0 ) {
		cerr << "Files not concatenated! Number of errors: " << -err << endl;
		return NULL;
	}

	pmap->sizeZ(i);		// Only z-slabs!

//	Bstring			filename = p->file_name().pre_rev('.') + ".tiles";
	Bstring			filename("y.tiles");
	ofstream		fd(filename.c_str());
	if ( fd.fail() ) {
		cerr << "Error: Tile parameter file %s not opened" << filename << endl;
		return NULL;
	}
	fd << pmap->sizeX() << " " << pmap->sizeY() << " " << pmap->sizeZ() << endl;
	fd << "0 0 0" << endl;
	y = ystart%ysize;
	if ( y ) fd << "0 0 0" << endl;;
	for ( ; y<pmap->sizeY(); y+=ysize ) fd << "0 " << y << " 0" << endl;
	fd.close();

	long		elementsize = pmap->channels()*pmap->data_type_size();
	long		linesize = elementsize*pmap->sizeX();
	
	if ( ystart + ysize > pmap->sizeY() ) ysize = pmap->sizeY() - ystart;

	pmap->sizeY(ysize);
	pmap->data_alloc();
	unsigned char*		tile = pmap->data_pointer();
	unsigned char*		data;
	
	if ( verbose ) {
		cout << "Extracting y tile from " << ystart << " to " << ystart+ysize-1 << endl;
		cout << "Tile size:                      " << pmap->size() << endl;
		cout << "Tile parameter file:            " << filename << endl;
	}

	for ( zt=0, thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		p = read_img(*thisfile, 1, 0);
		data = p->data_pointer();
		for ( z=0; z<p->sizeZ() && zt<pmap->sizeZ(); z++, zt++ ) {
			for ( y=ystart, yt=0; y<p->sizeY() && yt<pmap->sizeY(); y++, yt++ ) {
				i = (z*p->sizeY() + y)*linesize;
				j = (zt*pmap->sizeY() + yt)*linesize;
				memcpy(tile+j, data+i, linesize);
			}
		}
		delete p;
	}

	return pmap;
}

/**
@brief 	The lines along the z-dimensions of a disk-based block is Fourier back-transformed.  
@param 	*file_list	list of image file names.
@param 	&recfile	new reconstruction file name.
@param 	datatype	data type for new reconstruction file.
@param 	avg			target average.
@param 	std			target standard deviation.
@param 	cutmin		minimum for truncation.
@param 	cutmax		maximum for truncation.
@return	int					0.

	Each 2D xz plane is read from the raw complex data block.
	Each z-line in the plane is back-transformed.
	The transformed 2D plane is written back into the raw data block.

**/
Bimage*		img_backtransform_z_on_disk(Bstring* file_list, Bstring& recfile, 
				DataType datatype, double avg, double std, double cutmin, double cutmax)
{
	Bimage*				p = NULL;
	Bimage*				pmap = NULL;
	int					err(0);
	long		i(0), j(0);
	Bstring*			thisfile = NULL;
	Bstring				tempfile("temp.raw");
	
	// First read the file headers to identify problems
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		p = read_img(*thisfile, 0, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_DEBUG )
				cout << p->file_name() << ": " << 
						p->sizeX() << " " << p->sizeY() << " " << p->sizeZ() << " " << p->images() << " " << p->channels() << endl;
			if ( p->images() > 1 ) {
				cerr << p->file_name() << " must be a single image file!" << endl;
				err += -1;
			}
			if ( i == 0 ) {
				pmap = p;
			} else {
//				err += img_compatibility(pmap, p);
				if ( !p->compatible(pmap) ) err += -1;
				delete p;
			}
			i += p->sizeZ();
			j++;
		} else err += -1;
	}
	
	if ( err < 0 || j == 0 ) {
		cerr << "Files not concatenated! Number of errors: " << -err << endl;
		return NULL;
	}
	
	pmap->sizeZ(i);

	fstream*		ftemp = new fstream(tempfile.c_str());
	if ( ftemp->fail() ) {
		cerr << "Error: " << tempfile << " not opened!" << endl << endl;
		return NULL;
	}
	
	long	elementsize = pmap->channels()*pmap->data_type_size();
	long	datasize(0);
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		p = read_img(*thisfile, 1, -1);
		if ( !p ) {
			cerr << "Error: Image " << *thisfile << " not read!"  << endl;
			bexit(-1);
		} else {
			datasize = p->sizeX()*p->sizeY()*p->sizeZ()*p->images();
			ftemp->write((char *)p->data_pointer(), elementsize*datasize);
		}
		delete p;
	}
	
	img_backtransform_z_lines(ftemp, pmap);
	
	// Write the final real space map header
	pmap->file_name(recfile.str());
	if ( datatype == Unknown_Type ) pmap->data_type(Float);
	else pmap->data_type(datatype);

	img_write_data_block_with_type(ftemp, pmap, avg, std, cutmin, cutmax);
	
	ftemp->close();
	delete ftemp;
	
	remove(tempfile.c_str());

	return pmap;
}

int			img_backtransform_one_y_plane(Bimage* p, fft_plan plan, int y)
{
	long	i, x, z;
	long	slice_size = p->sizeX()*p->sizeY();
	long	zsize = p->sizeZ()*sizeof(Complex<float>);
	Complex<float>*	zline = new Complex<float>[zsize];

	if ( verbose & VERB_PROCESS )
		if ( y%10 == 0 ) cout << "Backtransforming z lines for y=" << y << endl;
	for ( x=0; x<p->sizeX(); x++ ) {
		for ( z=0, i=y*p->sizeX() + x; z<p->sizeZ(); z++, i+=slice_size ) zline[z] = p->complex(i);
		fftw(plan, zline);
		for ( z=0, i=y*p->sizeX() + x; z<p->sizeZ(); z++, i+=slice_size ) p->set(i, zline[z]);
	}

	delete[] zline;
	
	return 0;
}

int			img_backtransform_one_y_plane(fstream* ftemp, Bimage* p, fft_plan plan,
				int y, Complex<float>* zline, Complex<float>* cfdata)
{
	long	i, x, z, offset;
	long	cfsize = sizeof(Complex<float>);
	long	planesize = p->sizeX()*p->sizeZ();
	double			v;

	if ( verbose & VERB_PROCESS )
		if ( y%10 == 0 ) cout << "Backtransforming z lines for y=" << y << endl;
	for ( z=0; z<p->sizeZ(); z++ ) {
		offset = (z*p->sizeY() + y)*p->sizeX()*cfsize;
//		fseek(ftemp, offset, SEEK_SET);
//		fread(&cfdata[z*p->sizeX()], cfsize, p->sizeX(), ftemp);
		ftemp->seekg(offset, ios::beg);
		ftemp->read((char *)&cfdata[z*p->sizeX()], cfsize*p->sizeX());
	}
	for ( x=0; x<p->sizeX(); x++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) zline[z] = cfdata[z*p->sizeX() + x];
		fftw(plan, zline);
		for ( z=0; z<p->sizeZ(); z++ ) cfdata[z*p->sizeX() + x] = zline[z];
	}
	for ( z=0; z<p->sizeZ(); z++ ) {
		offset = (z*p->sizeY() + y)*p->sizeX()*cfsize;
//		fseek(ftemp, offset, SEEK_SET);
//		fwrite(&cfdata[z*p->sizeX()], cfsize, p->sizeX(), ftemp);
		ftemp->seekg(offset, ios::beg);
		ftemp->write((char *)&cfdata[z*p->sizeX()], cfsize*p->sizeX());
	}
	
	for ( i=0; i<planesize; i++ ) {
		v = cfdata[i].real();
		if ( p->minimum() > v ) p->minimum(v);
		if ( p->maximum() < v ) p->maximum(v);
		p->average(p->average() + v);
		p->standard_deviation(p->standard_deviation() + v*v);
	}
	
	return 0;
}

/**
@brief 	The lines along the z-dimensions of a disk-based block is Fourier back-transformed.  
@param 	*p			image header information (statistics updated).
@return	int					0.

	Each 2D xz plane is read from the raw complex data block.
	Each z-line in the plane is back-transformed.
	The transformed 2D plane is written back into the raw data block.

**/
int			img_backtransform_z_lines(Bimage* p)
{
	fft_plan		plan = fft_setup_plan(p->sizeZ(), 1, 1, FFTW_FORWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(p->sizeY(), dispatch_get_global_queue(0, 0), ^(size_t y){
		img_backtransform_one_y_plane(p, plan, y);
	});
#else
#pragma omp parallel for
	for ( long y=0; y<p->sizeY(); y++ )
		img_backtransform_one_y_plane(p, plan, y);
#endif

    fft_destroy_plan(plan);
	
	p->complex_to_real();
	p->statistics();
	
	if ( verbose ) {
		cout << "Reconstruction statistics:" << endl;
		cout << "Min and max:                    " << p->minimum() << " " << p->maximum() << endl;
		cout << "Avg and std:                    " << p->average() << " " << p->standard_deviation() << endl << endl;
	}
	
	return 0;
}

int			img_backtransform_z_lines(fstream* ftemp, Bimage* p)
{
	int				y;
	long	datasize = p->sizeX()*p->sizeY()*p->sizeZ()*p->images();

	fft_plan		plan = fft_setup_plan(p->sizeZ(), 1, 1, FFTW_FORWARD, 1);
	Complex<float>*	zline = new Complex<float>[p->sizeZ()];
	Complex<float>*	cfdata = new Complex<float>[p->sizeX()*p->sizeZ()];	

	p->minimum(1e37);
	p->maximum(-1e37);
	p->average(0);
	p->standard_deviation(0);
	
	for ( y=0; y<p->sizeY(); y++ )
		img_backtransform_one_y_plane(ftemp, p, plan, y, zline, cfdata);
	
    fft_destroy_plan(plan);
	delete[] zline;
	delete[] cfdata;	
	
	p->average(p->average() / datasize);
	p->standard_deviation(p->standard_deviation()/datasize - p->average()*p->average());
	if ( p->standard_deviation() > 0 ) p->standard_deviation(sqrt(p->standard_deviation()));
	else p->standard_deviation(0);
	
	if ( verbose ) {
		cout << "Reconstruction statistics:" << endl;
		cout << "Min and max:                    " << p->minimum() << " " << p->maximum() << endl;
		cout << "Avg and std:                    " << p->average() << " " << p->standard_deviation() << endl << endl;
	}
	
	return 0;
}


/**
@brief 	Get the size of a datatype.

	This function is used for calculating image data sizes for allocating
	memory and reading the data.

@param 	type	data type (defined in rwimg.h).
@return long   size of data type, if < 0 the data type is not supported.
**/
long   gettypesize(DataType type)
{
	long   size;
	
	switch ( type ) {
		case Bit:
		case UCharacter: case SCharacter:	size = sizeof(char); break;
		case UShort: case Short: size = sizeof(short); break;
		case UInteger: case Integer: 	size = sizeof(int); break;
		case ULong: case Long: 	size = sizeof(long); break;
		case Float:				size = sizeof(float); break;
		case Double:			size = sizeof(double); break;
		default: size = 0;
	}
	
	return size;
}


/**
@brief 	Set a pointer to a value with a given data type.

	A value is inserted into a given location with the given data type.
	The size of the allocated memory is channels*typesize.

@param 	*ptr		pointer to location.
@param 	value 	the value.
@param 	datatype	data type (defined in rwimg.h).
@return int				0, <0 on error.
**/
int			set_value_with_datatype(char* ptr, double value, DataType datatype)
{
	int 			typesize = gettypesize(datatype);
	if ( typesize < 1 ) return -1;
	
	switch ( datatype ) {
		case UCharacter:			*ptr = (unsigned char) value; break;
		case SCharacter: {
				signed char		cval = (signed char)    value;
				memcpy(ptr, &cval, typesize); break;
			}
		case UShort: {
				unsigned short	uval = (unsigned short) value;
				memcpy(ptr, &uval, typesize); break;
			}
		case Short: {
				short			sval = (short)          value;
				memcpy(ptr, &sval, typesize); break;
			}
		case UInteger: {
				unsigned int	ival = (unsigned int)	value;
				memcpy(ptr, &ival, typesize); break;
			}
		case Integer: {
				int 			ival = (int)            value;
				memcpy(ptr, &ival, typesize); break;
			}
		case ULong: {
				long	ival = (long)	value;
				memcpy(ptr, &ival, typesize); break;
			}
		case Long: {
				long 			ival = (long)			value;
				memcpy(ptr, &ival, typesize); break;
			}
		case Float: {
				float 			fval = (float)          value;
				memcpy(ptr, &fval, typesize); break;
			}
		case Double: {
				double 			dval = (double)         value;
				memcpy(ptr, &dval, typesize); break;
			}
		default: ;
	}
	
	return 0;
}

/*
@brief 	Rescales and writes a raw disk-based data block to a specific data type.  
@param 	fstream* ftemp		file pointer to the complex raw data block.
@param 	*p			image header information.
@return  	int					0.

	Each slice is read from the raw data block, rescaled, and written
	into the new file.
	Limitation: This function can only be used for an image file format
		with a contiguous data block (CCP4, MRC, PIF, etc.).

**/
int			img_write_data_block_with_type(fstream* ftemp, Bimage* p, 
				double avg, double std, double cutmin, double cutmax)
{

	long		i, z, offset;
	long		elementsize = p->data_type_size();
	long		cfsize = sizeof(Complex<float>);
	long		slicesize = p->sizeX()*p->sizeY();
    double				scale = 1;
	double				shift(0);

	switch ( p->data_type() ) {
		case UCharacter:			scale = 255/(p->maximum() - p->minimum());
							shift = -p->minimum()*scale; break;
		case SCharacter:			scale = 255/(p->maximum() - p->minimum()); 
							shift = -128 - p->minimum()*scale; break;
		case UShort:		scale = USHRT_MAX/(p->maximum() - p->minimum());
							shift = -p->minimum()*scale; break;
		case Short:			scale = SHRT_MAX/(p->maximum() - p->minimum()); break;
		case Integer:		scale = 2e9/(p->maximum() - p->minimum()); break;
		case Float:
		case Double:		scale = 1;  break;
		default: ;		
	}
		
	if ( std ) {
		scale = std/p->standard_deviation();
		shift = avg - p->average()*scale;
	}
	
	p->minimum(scale*p->minimum() + shift);
	p->maximum(scale*p->maximum() + shift);
	p->average(scale*p->average() + shift);
	p->standard_deviation(p->standard_deviation() * scale);

	if ( cutmin >= cutmax ) {
		cutmin = p->minimum();
		cutmax = p->maximum();
	}

	p->minimum(p->data_type_min());
	p->maximum(p->data_type_max());

	if ( cutmin < p->minimum() ) cutmin = p->minimum();
	if ( cutmax > p->maximum() ) cutmax = p->maximum();

	p->minimum(cutmin);
	p->maximum(cutmax);
	
	p->fourier_type(NoTransform);
	
	write_img(p->file_name(), p, 0);

	double			d;
	char*			data = new char[slicesize*elementsize];
	Complex<float>*	cfdata = new Complex<float>[slicesize];
	
	fstream			fmap(p->file_name().c_str());
	if ( fmap.fail() ) return -1;
	
	for ( z=0; z<p->sizeZ(); z++ ) {
		offset = z*slicesize*cfsize;
		ftemp->seekg(offset, ios::beg);
		ftemp->read((char *)cfdata, cfsize*slicesize);
		for ( i=0; i<slicesize; i++ ) {
			d = cfdata[i].real()*scale + shift;
			if ( d < cutmin ) d = cutmin;
			if ( d > cutmax ) d = cutmax;
			set_value_with_datatype(data+i*elementsize, d, p->data_type());
		}
		offset = p->data_offset() + z*slicesize*elementsize;
		fmap.seekg(offset, ios::beg);
		fmap.write((char *)data, elementsize*slicesize);
	}

	fmap.close();
	
	delete[] data;
	delete[] cfdata;
	
	return 0;
}

/**
@brief 	Creates a reciprocal space mask from the tilt series orientations.  
@param 	*project 	project parameter structure.
@param 	size		size of mask.
@param 	origin		origin of mask.
@param 	hi_res		high resolution limit.
@param 	scale		scale of mask.
@return	Bimage*		mask, NULL on failure.

	The orientation parameters of the micrographs must all be set. 
	A 2D mask is overlayed onto the 3D image for each micrograph.

**/
Bimage*		project_missing_mask(Bproject* project, Vector3<long> size,
				Vector3<double> origin, double hi_res, double scale)
{
	Bfield*			field = project->field;

	if ( !field ) {
		cerr << "Error in project_missing_mask: Field-of-view not defined!" << endl;
		return NULL;
	}
	
	Bmicrograph*	mg = field->mg;

	if ( !mg ) {
		cerr << "Error in project_missing_mask: Micrograph not defined!" << endl;
		return NULL;
	}
	
	Breconstruction*	rec = NULL;

	if ( size.volume() < 10 ) {
		rec = project->rec;
		if ( !rec ) {
			cerr << "Error in project_missing_mask: Reconstruction and mask size not defined!" << endl;
			return NULL;
		}
		size = rec->box_size;
	}
		
	if ( size.volume() < 10 ) {
		cerr << "Error in project_missing_mask: Mask size not defined!" << endl;
		return NULL;
	}
	
	if ( mg->matrix.determinant() < 0.9 )
		project_mg_tilt_to_matrix(project);
	
	// Image to pack
	Bimage* 		pimg = new Bimage(UCharacter, TSimple, size[0], size[1], 1, 1);
	pimg->fill(1);
	
	// Mask
	Bimage* 		pmask = new Bimage(UCharacter, TSimple, size, 1);
	pmask->sampling(mg->pixel_size[0]/scale, mg->pixel_size[0]/scale, mg->pixel_size[0]/scale);
	pmask->check_resolution(hi_res);

//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "Generating a 3D reciprocal space mask from a tilt series:" << endl;
		cout << "Mask size:                      " << pmask->sizeX() << " " << pmask->sizeY() << " " << size[2] << endl;
		cout << "Origin:                         " << origin << endl;
		cout << "Resolution limit:               " << hi_res << " A" << endl;
		cout << "Scale:                          " << scale << endl;
	}

	for ( mg=field->mg; mg; mg=mg->next ) if ( mg->select ) {
		
		pmask->mask_pack_plane(mg->matrix, hi_res, scale);
				
	}
	
	delete pimg;
	
	Vector3<long>	shift((long) -(origin[0]+0.5), (long) -(origin[1]+0.5), (long) -(origin[2]+0.5));
	pmask->shift_wrap(shift);
	
	pmask->statistics();
	
	return pmask;
}

