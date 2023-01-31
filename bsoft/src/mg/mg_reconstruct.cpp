/**
@file	mg_reconstruct.cpp
@brief	Functions for reconstruction
@author Bernard Heymann
@date	Created: 20010403
@date	Modified: 20220221
**/

#include "rwimg.h"
#include "rwmg.h"
#include "mg_processing.h"
#include "mg_select.h"
#include "mg_img_proc.h"
#include "mg_particle_select.h"
#include "mg_reconstruct.h"
#include "mg_ctf.h"
#include "Complex.h"
#include "symmetry.h"
#include "utilities.h"

#include <sys/stat.h>
#include <fcntl.h>


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Sets the Fourier transform size for 2D particle images for reconstruction.
@param 	xsize			size of x dimension.
@param 	scale			scale of reconstruction.
@param 	pad_factor		factor that determines image padding.
@return int				transform size.

	The reconstruction size must be set.

**/
int			part_ft_size(int xsize, double scale, int pad_factor) 
{
	if ( xsize <= 0 ) {
		error_show("Error in part_ft_size", __FILE__, __LINE__);
		cerr << "Size not set!" << endl << endl;
		return 0;		
	}

	int			ft_size(xsize);

	if ( scale > 1 ) ft_size = (int) (xsize*scale);

	if ( pad_factor < 1 ) return ft_size;

	if ( pad_factor > 8 ) pad_factor = 8;
	
	ft_size = findNextPowerOf((int)(pad_factor*xsize*scale), 2);
	
	if ( ft_size < xsize )
//		ft_size = findNextPowerOf(pad_factor*xsize, 2);
		ft_size = findNextPowerOf(xsize, 2);
	
//	if ( ft_size < xsize/scale )
//		ft_size = findNextPowerOf((int)(xsize/scale), 2);

	return ft_size;
}

/**
@brief 	Reciprocal space reconstruction from 2D particle images.  
@param 	*partlist		a list of 2D particle image parameters.
@param 	*sym			point group symmetry.
@param	sym_mode		0=apply symmetry, 1=C1, 2=random symmetry view
@param 	hi_res			high resolution limit.
@param 	scale			scale of reconstruction.
@param 	sam				sampling/voxel size of reconstruction.
@param 	size			size of reconstruction.
@param 	ft_size			Fourier transform size.
@param 	plan			Fourier transform plan.
@param	interp_type		interpolation type.
@param 	ctf_action		flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	flags			1=rescale particles, 2=2D reconstruction, 4=bootstrap, 8=Ewald.
@param 	first			flag to indicate the first thread.
@return	Bimage*			3D reconstructed map.

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
	A bootstrap reconstruction uses the particle selection to weigh each
	selected particle.

**/
Bimage*		particle_reconstruct(Bparticle* partlist, Bsymmetry sym, int sym_mode,
				double hi_res, Vector3<double> scale, Vector3<double> sam, Vector3<long> size,
				int ft_size, fft_plan plan, int interp_type,
				int ctf_action, double wiener, int flags, int first)
{
	random_seed();
	
	int				twoD_flag(flags & 2), bootstrap(flags & 4), ewald(flags & 8);
	Bmicrograph*	mg = partlist->mg;
	Bparticle*		part = partlist;
	
	if ( sam[0] < 0.1 ) sam = part->pixel_size;

//	Vector3<double>	pixel_size(mg->pixel_size/scale);
//	Vector3<double>	pixel_size(partlist->pixel_size/scale);
	Vector3<double>	pixel_size(sam/scale);
	pixel_size[1] = pixel_size[0];
	if ( twoD_flag ) pixel_size[2] = 1;			// Isotropic sampling in 2D
	else pixel_size[2] = pixel_size[0];			// Isotropic sampling in 3D

	Bimage*			p;
	Bimage* 		prec = new Bimage(Float, TComplex, size, 1);
	prec->sampling(pixel_size);

	prec->check_resolution(hi_res);
	
	long			nsel = particle_count(partlist);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particle_reconstruct: nsel = " << nsel << endl;

	CTFparam		em_ctf;
	
	long 			nrec(0), img_num;
	double			part_weight(1), pad_ratio(1);
//	Vector3<double>	part_scale(scale);
	Vector3<double>	part_scale(1,1,1);
	View			view;
	Bstring			partfile;
	
//	cout << "CTF action = " << ctf_action << endl;
//	cout << "volt=" << partlist->mg->ctf->volt() << endl;
	
	for ( part = partlist; part; part = part->next ) {
//		cout << "Particle: " << part->id << endl;
//		cout << "pixel_size=" << part->pixel_size << endl;
		mg = part->mg;
		if ( mg->ctf ) em_ctf.update(mg->ctf);
		partfile = mg->fpart;
		img_num = part->id - 1;
		if ( partfile.length() < 1 ) {
			partfile = part->fpart;
			img_num = 0;
		}
		if ( ( p = read_img(partfile, 1, img_num) ) == NULL ) {
			error_show("particle_reconstruct", __FILE__, __LINE__);
			return NULL;
		}
					
		p->change_type(Float);
		if ( flags & 1 ) p->rescale_to_avg_std(0, 1);
		p->calculate_background();
		p->sampling(part->pixel_size);
//		cout << "particle sampling = " << p->sampling(0) << endl;
		
		// set the origin
		if ( part->ori[0] <= 0 ) {
			if ( p->image->origin()[0] > 0 ) part->ori = p->image->origin();
			else part->ori = prec->image->origin();
		}
		part->ori[2] = 0;
		p->origin(part->ori);

		// set the view
		p->view(part->view);
		if ( sym_mode == 2 ) view = random_symmetric_view(part->view, sym);
		else view = part->view;
				
		if ( prec->sizeZ() < 2 ) {
			if ( p->image->view()[2] >= 0 ) p->image->view(0,0,1,p->image->view().angle());
			else p->image->view(0,0,-1,p->image->view().angle());
		}
		
		// set the particle scaling
		part_scale = Vector3<double>(1,1,1);
		if ( part->mag > 0 ) part_scale /= part->mag;

		if ( part->def > 0 ) em_ctf.defocus_average(part->def);

		part_weight = 1;
		if ( bootstrap ) part_weight = part->sel;
		
		pad_ratio = ft_size*1.0L/p->sizeX();
	
		if ( pad_ratio > 1 ) {
//			pad_ratio *= pad_ratio;
			pad_ratio *= pad_ratio*0.6;
//			pad_ratio *= pad_ratio*M_PI/4.0;
//			pad_ratio *= sqrt(pad_ratio);
			part_weight /= pad_ratio;
			p->pad(ft_size, FILL_BACKGROUND);
//			p->pad(ft_size, FILL_AVERAGE, p->average());
//			p->pad(ft_size, FILL_USER, 0);
		}
	
		p->fft(plan, 1);
		p->phase_shift_to_origin();
		
		//Ewald sphere offset
		double		ew_wl(0);
		if ( ewald ) ew_wl = em_ctf.lambda();

		if ( ctf_action )
	 		img_ctf_apply_complex(p, em_ctf, (ctf_action==1), wiener, 0, hi_res);
		
		if ( sym_mode )
			prec->fspace_pack_2D(p, view.matrix(), hi_res, 0, part_scale, ew_wl, part_weight, interp_type);
		else
			prec->fspace_pack_2D(p, view, sym, hi_res, 0, part_scale, ew_wl, part_weight, interp_type);

		if ( mg->ctf ) em_ctf.defocus_average(mg->ctf->defocus_average());
				
		delete p;
					
		nrec++;
					
		if ( first && ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) ) )
			cerr << "Complete:                       " << setprecision(3)
							<< nrec*100.0/nsel << " %    \r" << flush;
	}
	if ( first && ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) ) )
		cout << endl;
	
	prec->image->FOM(nrec);
	
	if ( verbose & VERB_FULL )
		cout << "Particles used:                 " << nrec << endl << endl;
	
	return prec;
}


/**
@brief 	Combines and weighs a map from several partial maps and weight sets.
@param 	**pacc			array of partial maps with linked weight maps.
@param 	imap			which output map to weigh.
@param 	nmaps			number of output maps (1,2,3).
@param 	maps_per_class	number of threads per map (= number of partial maps).
@param 	hi_res			high resolution limit.
@return Bimage*			weighed reconstruction with FOM block.

	The input is a set of partially integrated complex maps with associated 
	weigths as follows:
		voxel power sums			FOM block of map.
		voxel weight sums			linked image.
		voxel weight squared sums	linked image FOM block.
	The partial sums are completed into corresponding blocks in three possible
	ways based on the value of imap and nmap:
		nmap	imap	result
		1		0		one map from all input maps
		2		0,1		one map from half of the input maps
		3		0		one map from all input maps
		3		1,2		one map from half of the input maps
	The total number of maps in the array is nclasses*nmaps.

**/
Bimage*		img_reconstruction_sum_weigh(Bimage** pacc, int imap, int nmaps, int maps_per_class, double hi_res)
{
	int					iclass = imap/nmaps;

	Bimage*				prec = new Bimage(Float, TComplex, pacc[0]->size(), pacc[0]->images());
	prec->sampling(pacc[0]->sampling(0));
	prec->fourier_type(Standard);
	
	prec->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	prec->next->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	prec->next->next->next = new Bimage(Float, TSimple, prec->size(), prec->images());
	float*				fom = (float *) prec->next->data_pointer();
	float*				weight = (float *) prec->next->next->data_pointer();
	float*				weight2 = (float *) prec->next->next->next->data_pointer();
	
	long				i, j, k, nv(0);
	long			   	ds = prec->image_size();
	float*				fom_acc;
	float*				w_acc;
	float*				w2_acc;
	
	long				st(0), inc(1);
	if ( nmaps == 2 ) {
		st = imap%2;
		inc = 2;
	}
	if ( nmaps == 3 && imap%3 ) {
		st = imap%3 - 1;
		inc = 2;
	}

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_reconstruction_sum_weigh: maps_per_class=" << maps_per_class << endl;
		cout << "DEBUG img_reconstruction_sum_weigh: iclass=" << iclass << endl;
		cout << "DEBUG img_reconstruction_sum_weigh: st=" << st << endl;
		cout << "DEBUG img_reconstruction_sum_weigh: inc=" << inc << endl;
	}
	
//	double			ssw(0);
//	cout << "maps_per_class = " << maps_per_class << endl;
	for ( i=iclass*maps_per_class+st, j=0; j<maps_per_class; i+=inc, j+=inc ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG img_reconstruction_sum_weigh: Doing map " << i << endl;
		fom_acc = (float *) pacc[i]->next->data_pointer();
		w_acc = (float *) pacc[i]->next->next->data_pointer();
		w2_acc = (float *) pacc[i]->next->next->next->data_pointer();
//		double			sw(0);
		for ( k=0; k<ds; k++ ) {
			prec->add(k, pacc[i]->complex(k));
			fom[k] += fom_acc[k];
			weight[k] += w_acc[k];
			weight2[k] += w2_acc[k];
//			sw += w_acc[k];
		}
//		ssw += sw;
//		if ( imap == 0 ) cout << "total weight for " << i << " = " << sw << endl;
//		cout << "map " << imap << " block " << i << " = " << sw << endl;
		nv++;
	}
//	cout << "map " << imap << " total weight = " << ssw << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_reconstruction_sum_weigh: volumes accumulated: " << nv << endl;

	long			cov = prec->fspace_reconstruction_weigh();
	double			rad = prec->fspace_maximum_radius(hi_res) + 1;
	double			vol = (4.0/3.0)*M_PI*rad*rad*rad;
	
	double			Rf = prec->friedel_check();
//	if ( Rf > 0.001 ) prec->friedel_apply();
	
	if ( verbose ) {
		if ( inc == 1 )
			cout << "Full reconstruction:" << endl;
		else 
			cout << "Half reconstruction " << st+1 << ":" << endl;
		cout << "Friedel symmetry residual:      " << Rf << endl;
		cout << "Coverage:                       " << cov << " (" << cov*100.0/vol << " %)" << endl << endl;
	}
	
	if ( imap%nmaps == 0 )
		prec->fspace_reconstruction_stats(hi_res);

	prec->fspace_reconstruction_snr();
	
	if ( verbose & VERB_DEBUG )
		write_img("weight.map", prec->next->next, 0);
	
	delete prec->next->next;
	prec->next->next = NULL;

	prec->phase_shift_to_center();
	prec->origin(prec->size()/2);
	
	return prec;
}



double		part_single_reconstruction(Bparticle* part, Vector3<long> size,
				Bimage* pmask, Bsymmetry& sym, double hi_res,
				Vector3<double> scale, int ft_size,
				fft_plan img_plan, fft_plan map_plan, double part_weight,
				int interp_type, CTFparam* mg_ctf, int ctf_action,
				double wiener, int flags)
{
	FOMType 		fom_tag[NFOM] = {DENSITY, COVERAGE};

	if ( ppx_check(part, fom_tag) ) return part->fom[0];

	double			dens(0);
	Bparticle*		partone = NULL;
	particle_copy(&partone, part);

	Bimage*			prec = particle_reconstruct(partone, sym, 0, hi_res,
						scale, partone->pixel_size, size, ft_size, img_plan,
						interp_type, ctf_action, wiener, 0, 0);
	
	long			ncov = prec->fspace_reconstruction_weigh();
	double			sf = prec->real_size()[0]/hi_res;
	double			cov = ncov / ((4.0*M_PI/3.0)*sf*sf*sf);

	prec->phase_shift_to_center();
	prec->origin(prec->size()/2);

	prec->fft(map_plan, 1);
	prec->complex_to_real();

	prec->statistics();

//	cout << "avg=" << prec->average() << " std=" << prec->standard_deviation() << endl;
					
	if ( pmask ) dens = prec->relative_density(pmask);
	else dens = prec->standard_deviation();

//	Bstring			fn = "part_" + part->mg->id + Bstring(part->id, "_%04d.mrc");
//	write_img(fn, prec, 0);
	
	delete prec;
	
	part->fom[0] = dens;
	part->fom[1] = cov;
	
	if ( flags & WRITE_PPX ) {
		Bstring		ppx_name = ppx_filename(part->mg->id, part->id);
		write_particle(ppx_name, part, 0, 0, fom_tag);
		ppx_name = 0;
	}
	
	return dens;
}

/**
@brief 	Reciprocal space reconstruction from the images in a multi-image file.  
@param 	*project		image processing parameter structure.
@param 	&maskfile		mask to determine density statistics.
@param 	*sym			point group symmetry.
@param 	num_select		selection number from the selection column.
@param 	hi_res			high resolution limit.
@param 	scale			scale of reconstruction.
@param 	size			size of reconstruction.
@param 	pad_factor		image padding factor.
@param	interp_type		interpolation type (provisional).
@param 	ctf_action		flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	flags			option flags.
@return	long			particles processed.

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
	A bootstrap reconstruction uses the particle selection to weigh each
	selected particle.

**/
long		project_single_particle_reconstruction(Bproject* project, 
				Bstring& maskfile, Bsymmetry& sym,
				int num_select, double hi_res, Vector3<double> scale, Vector3<long> size, 
				int pad_factor, int interp_type, int ctf_action, double wiener, int flags)
{
	Vector3<double>		sam;
	
	if ( size.volume() <= 0 )
		size = project_set_reconstruction_size(project, sam, scale[0], 0);

	if ( size.volume() <= 0 ) {
		error_show("Error in project_single_particle_reconstruction", __FILE__, __LINE__);
		cerr << "No particles to determine the size from!" << endl << endl;
		return -1;
	}

	int				ft_size = part_ft_size(size[0], scale[0], pad_factor);

	fft_plan		img_plan = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	fft_plan		map_plan = fft_setup_plan(size[0], size[1], size[2], FFTW_BACKWARD, 1);
	
	long			i, j;

	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;

	for ( field = project->field; field && part == NULL; field = field->next )
		for ( mg = field->mg; mg && part == NULL; mg = mg->next )
			part = mg->part;
				
	if ( size.volume() <= 0 ) {
		error_show("Error in project_single_particle_reconstruction", __FILE__, __LINE__);
		cerr << "No reconstruction size set!" << endl << endl;
		return -1;
	}
	
	Bimage*			pmask = NULL;

	if ( hi_res < 2*part->pixel_size[0]/scale[0] ) hi_res = 2*part->pixel_size[0]/scale[0];
	
	if ( maskfile.length() )
		if ( ( pmask = read_img(maskfile, 1, 0) ) == NULL ) {
			error_show("project_single_particle_reconstruction", __FILE__, __LINE__);
			return -1;
		}
	
	long   datasize = (long) size.volume();
	long   maskedsize = datasize;
	
	if ( pmask ) for ( i=maskedsize=0; i<datasize; i++ ) if ( (*pmask)[i] == 1 ) maskedsize++;
	
	long			nsel = project_count_mg_part_selected(project, num_select);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_single_particle_reconstruction: num_select = " << num_select << endl;
	
	if ( verbose ) {
		cout << "3D reciprocal space reconstruction:" << endl;
		if ( pmask )
			cout << "Mask:                           " << pmask->file_name() << endl;
		cout << "Map size:                       " << size << endl;
		cout << "Resolution limit:               " << hi_res << " A" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Interpolation type:             " << interp_type << endl;
		cout << "CTF application type:           " << ctf_action << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Analysis volume:                " << maskedsize << " (" << maskedsize*100.0/datasize << " %)" << endl;
		cout << "Number of selected particles:   " << nsel << endl << endl;
	}

	long 			np(0);
	double			part_weight(1);
	Bparticle**		part_array;

	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( i=0, part = mg->part; part; part = part->next ) if ( part->sel ) {
				if ( num_select < 0 ) i++;
				else if ( part->sel == num_select ) i++;
			}
			part_array = new Bparticle*[i];
			for ( i=0, part = mg->part; part; part = part->next ) if ( part->sel ) {
				if ( num_select < 0 ) part_array[i++] = part;
				else if ( part->sel == num_select ) part_array[i++] = part;
			}

#ifdef HAVE_GCD
			dispatch_apply(i, dispatch_get_global_queue(0, 0), ^(size_t j){
				part_array[j]->fom[0] = part_single_reconstruction(part_array[j],
					size, pmask, sym, hi_res, scale,
					ft_size, img_plan, map_plan, part_weight,
					interp_type, mg->ctf, ctf_action, wiener, flags);
			});
#else
#pragma omp parallel for
			for ( j=0; j<i; j++ )
				part_array[j]->fom[0] = part_single_reconstruction(part_array[j],
					size, pmask, sym, hi_res, scale,
					ft_size, img_plan, map_plan, part_weight,
					interp_type, mg->ctf, ctf_action, wiener, flags);
#endif
			
			np += i;
			if ( verbose ) {
				cout << "Micrograph: " << mg->id << endl;
				cout << "PID\tFOM\tCov" << endl;
				for ( j=0; j<i; j++ ) cout << part_array[j]->id << tab << part_array[j]->fom[0]
					<< tab << part_array[j]->fom[1] << endl;
				cout << "Particles processed:    " << np << endl;
			}
			delete[] part_array;
		}
	}

	fft_destroy_plan(img_plan);
	fft_destroy_plan(map_plan);
	
	if ( verbose )
		cout << endl;
	
	return np;
}


/**
@brief 	Sets the reconstruction size for reconstruction from 2D particle images.
@param 	*project		project parameter structure.
@param 	&sam			sampling/voxel size of new map.
@param 	scale			scale of reconstruction.
@param 	twoD_flag		doing a 2D reconstruction rather than 3D.
@return Vector3<long>		reconstruction size.

	The reconstruction size is set from the first particle image size found.

**/
Vector3<long>	project_set_reconstruction_size(Bproject* project, Vector3<double>& sam, double scale, int twoD_flag)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_reconstruction_size: scale=" << scale << endl;
	
	Vector3<long>	size;
	Bfield*			field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bstring			partfile;

	for ( field = project->field; field && part == NULL; field = field->next )
		for ( mg = field->mg; mg && part == NULL; mg = (part)? mg: mg->next )
			part = mg->part;
				
	if ( !mg ) {
		error_show("Error in project_set_reconstruction_size", __FILE__, __LINE__);
		cerr << "No micrograph with particles!" << endl << endl;
		return size ;
	}
	
	partfile = mg->fpart;
	if ( !partfile.length() ) partfile = part->fpart;
	
	if ( !partfile.length() ) {
		error_show("Error in project_set_reconstruction_size", __FILE__, __LINE__);
		cerr << "No particle image file name given for micrograph " << mg->id << endl << endl;
		return size;
	}

	Bimage* 		p = read_img(partfile, 0, 0);
	
	if ( sam[0] < 0.1 ) sam = part->pixel_size;
	if ( sam[0] < 0.1 ) sam = mg->pixel_size;
	if ( sam[0] < 0.1 ) sam = p->sampling(0);
	sam[2] = sam[1] = sam[0];

	size[0] = (int) (p->sizeX()*scale);
	if ( size[0] < p->sizeY()/scale ) size[0] = (int) (p->sizeY()*scale);
	if ( size[0] < p->sizeZ()/scale ) size[0] = (int) (p->sizeZ()*scale);
	size[2] = size[1] = size[0];
	
	if ( twoD_flag ) {
		sam[2] = 1;
		size[2] = 1;
	}
	
	delete p;
	
	size = size.max(1);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_reconstruction_size: size=" << size << endl;
	
	return size;
}

/**
@brief 	Sets up the particle selection for reconstruction.
@param 	*project		project parameter structure.
@param 	classes			string specifying classes to use.
@param 	nmaps			number of maps per class (1,2,3).
@param 	nthreads		number of threads per class (must be even if nmaps > 1).
@return int				number of classes.

	The classes are specified in a string of comma-separated numbers,
	also allowing hyphened ranges (e.g., "2,5-7,9").
	The selection numbers for the particles in the project are set to
	calculate partial maps so that there are nthread maps per class.
	Each such partial map will be calculated in its own thread and
	integrated with others from the same class afterwards.

**/
int			project_configure_for_reconstruction(Bproject* project, Bstring classes, long& nmaps, int& nthreads)
{
	if ( !classes.length() ) classes = "selected";
	classes = classes.lower();
	if ( nmaps < 1 ) nmaps = 1;
	if ( nmaps > 3 ) nmaps = 3;
	if ( nthreads < 1 ) nthreads = 1;
	if ( nmaps&2 && nthreads%2 ) nthreads++;
	
	long				i, msel(1), nsel(0), nclass(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( classes[0] == 'a' ) {
		part_reset_selection(project, 3);
	} else if ( classes[0] == 's' ) {
		part_consolidate_selection(project, 1);
	} else {
		if ( project->select ) {
			for ( rec = project->rec; rec; rec = rec->next )
				for ( part = rec->part; part; part = part->next )
					if ( msel < part->sel ) msel = part->sel;
		} else {
			for ( field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					for ( part = mg->part; part; part = part->next )
						if ( msel < part->sel ) msel = part->sel;
		}
		if ( verbose )
			cout << "Maximum selection number:       " << msel << endl;
	}
	
	msel++;

	vector<int>			sel(msel, 0);
	vector<int>			cnt(msel, 0);
	vector<int>			num(msel, 0);

	if ( classes[0] == 'a' || classes[0] == 's' ) {
		sel[1] = 1;
		nclass = 1;
	} else {

		sel = select_numbers(classes, msel);
		for ( i=1; i<msel; i++ ) if ( sel[i] ) sel[i] = ++nclass;
	}

//	if ( verbose & VERB_FULL ) {
	if ( verbose ) {
		cout << "Selection:" << endl << "Class\tNumber" << endl;
		for ( i=0; i<msel; i++ ) cout << i << tab << sel[i] << endl;
	}
	
	// Set up the total number of threads
	int					maps_per_class = nthreads;
	if ( verbose ) {
		cout << "Classes:                        " << classes << " (" << nclass << ")" << endl;
		cout << "Maps per class:                 ";
		if ( nmaps == 1 ) cout << "Full set map" << endl;
		if ( nmaps == 2 ) cout << "Two half set maps" << endl;
		if ( nmaps == 3 ) cout << "Full set map and two half set maps" << endl;
		cout << "Threads per class:              " << nthreads << endl;
		cout << "Total number of threads:        " << nclass*nthreads << endl << endl;
	}
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				i = part->sel;
				if ( i>0 && i<msel && sel[i] ) {
					num[i]++;
					cnt[i]++;
					part->sel = (sel[i] - 1)*maps_per_class + cnt[i];
					if ( cnt[i] >= maps_per_class ) cnt[i] = 0;
					nsel++;
				} else {
					part->sel = 0;
					num[0]++;
				}
			}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					i = part->sel;
//					cout << mg->id << tab << part->id << tab << part->sel << endl;
					if ( i>0 && i<msel && sel[i] ) {
						num[i]++;
						cnt[i]++;
						part->sel = (sel[i] - 1)*maps_per_class + cnt[i];
//						cout << sel[i] << tab << cnt[i] << tab << part->sel << endl;
						if ( cnt[i] >= maps_per_class ) cnt[i] = 0;
						nsel++;
					} else {
						part->sel = 0;
						num[0]++;
					}
				}
	}
	
	if ( verbose ) {
		cout << "Particles selected:" << endl;
		cout << "Class\tNumber" << endl;
		for ( i=0; i<msel; i++ ) if ( sel[i] )
			cout << i << tab << num[i] << endl;
		cout << endl;
	}
	
	return nclass;
}

int			project_update_class_averages(Bproject* project, Bimage* prec, Bstring file_name)
{
	long			i;
	Bparticle*		part;
	
	project->class_avg = part = NULL;

	for ( i=0; i<prec->images(); ++i ) {
		part = particle_add(&part, i+1);
		if ( !project->class_avg ) project->class_avg = part;
		part->mag = 1;
		part->pixel_size = prec->image[i].sampling();
		part->ori = prec->image[i].origin();
		part->view = prec->image[i].view();
		part->fom[0] = prec->image[i].FOM();
		part->sel = prec->image[i].select();
		part->fpart = file_name;
	}

	project_show_class_averages(project);
	
	return 1;
}

/**
@brief 	Creates a 2D reconstruction from the images in a  multi-image file.  
@param 	*project 			image processing parameter structure.
@param	file_name			2D reconstruction file name.
@param 	transform_output	flag to output transformed images.
@return	Bimage* 				2D reconstruction image.

	The angle of rotation and the x,y origins must already have been found 
	and placed in the appropriate arrays within the Bproject structure.
	Each selected image is transformed, and then added to the reconstruction
	image corresponding to the original projection image chosen.
	If the transform_output flag is set, then the transformed images are written
	into a new image with a "_proj.spi" ending.

**/
Bimage* 	project_reconstruct_2D(Bproject* project, Bstring file_name, int transform_output)
{
	bool			invert(0);
	long			n, i, j, k, npart;

	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;

	Bimage* 		p = read_img(mg->fpart, 0, 0);

	Vector3<long>	size(p->size());
	double			px(part->pixel_size[0]), hires(2*px);

	delete p;
	
	// Determine the number of reconstructions to generate
	long 			nmap(0);
	for ( field=project->field; field; field=field->next )
		for ( mg=field->mg; mg; mg=mg->next )
			for ( part=mg->part; part; part=part->next )
				if ( part->sel > nmap ) nmap = part->sel;

	if ( nmap < 1 ) return NULL;
	
	Bimage* 		prec = new Bimage(Float, TSimple, size, nmap);	
	prec->origin(prec->size()/2);
	prec->sampling(px, px, 1);

	Bstring			filename;
	
	int* 			num = new int[nmap];
	float* 			fom = new float[nmap];
	for ( i=0; i<nmap; i++ ) {
		num[i] = 0;
		fom[i] = 0;
	}

	if ( verbose & VERB_RESULT ) {
		cout << "2D real space reconstruction:" << endl;
		cout << "Map size:                       " << prec->size() << endl;
		cout << "Map origin:                     " << prec->image->origin() << endl;
		cout << "Number of maps:                 " << prec->images() << endl << endl;
	}

	Vector3<double>	origin, translate;
	Vector3<double>	scale(1,1,1);
	Vector3<double> axis(0,0,1);
	Matrix3			mat(1);
	Euler			euler;
	
	Bimage* 		proj = NULL;
	
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			npart = micrograph_count_particles(mg);
			if ( verbose )
				cout << "Micrograph: " << mg->id << " with " << npart << " particles" << endl;
			if ( npart ) {
				if ( transform_output ) {
					proj = new Bimage(Float, TSimple, size, npart);
					proj->sampling(px, px, 1);
				}
				for ( part=mg->part; part; part=part->next ) {
					j = part->id - 1;
					n = part->sel - 1;
//					p = read_img(mg->fpart, 1, j);
					p = particle_read_img(part, 1);
					p->change_type(Float);
					p->rescale_to_avg_std(0, 1);
					p->shift_background(0);
					p->sampling(px, px, 1);
					if ( mg->ctf )
						img_ctf_apply(p, *(mg->ctf), 1, 0.2, 0, hires, invert);
					if ( n < 0 ) {
						translate = prec->image->origin() - part->ori;
					} else {
						translate = prec->image[n].origin() - part->ori;
					}
					origin = part->ori;
					euler = Euler(part->view);
					mat = Matrix3(axis, euler.psi());
//					mat = Matrix3(axis, part->view.angle());
					p->transform(scale, origin, translate, mat, FILL_BACKGROUND, 0);
					if ( part->sel ) {
						num[n]++;
						fom[n] += part->fom[0];
						prec->add(n, p);
					}
					if ( transform_output ) {
						proj->image[j].view(0,0,1,0);
						proj->replace(j, p);
					}
					delete p;
				}
				if ( transform_output ) {
					proj->origin(prec->image->origin());
					filename = mg->fpart.pre_rev('.') + "_proj.spi";
					if ( verbose & VERB_RESULT )
						cout << "Writing projections file:       " << filename << endl;
					write_img(filename, proj, 0);
					delete proj;
				}
			}
		}
	}
	
	npart = project_count_mg_particles(project);
	
	// Normalize the new reconstructions
	for ( i=n=0; n<prec->images(); n++ ) {
		if ( num[n] ) fom[n] /= num[n];
		prec->image[n].FOM(fom[n]);
		prec->image[n].select(num[n]);
		for ( k=0; k<prec->image_size(); i++, k++ )
			if ( num[n] ) prec->set(i, (*prec)[i] / num[n]);
	}
	
	prec->statistics();
	
	if ( verbose & VERB_RESULT ) {
		cout << "Particles used in reconstruction:" << endl;
		cout << "Map\tNumber\t%\tFOM" << endl;
		for ( n=0; n<nmap; n++ ) {
			cout << n+1 << tab << num[n] << tab << 
					num[n]*100.0/npart << tab << fom[n] << endl;
		}
		cout << endl;
	}
	
	delete[] num;
	delete[] fom;

	project_update_class_averages(project, prec, file_name);

	return prec;
}				

Bparticle*	img_reconstruct_2D(long i, Bproject* project, Bimage* prec, 
				double hires, fft_plan planf, fft_plan planb)
{
	bool			invert(0);
	long			n(0);
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;
	Bimage*			p = NULL;

	double			f(0);
	Vector3<double>	origin, translate;
	Vector3<double>	scale(1,1,1);
	Vector3<double> axis(0,0,1);
	Matrix3			mat(1);
	Euler			euler;
	
//	cout << "Map: " << i+1 << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) if ( part->sel == i+1 ) {
//				p = read_img(mg->fpart, 1, part->id - 1);
				p = particle_read_img(part, 1);
				if ( p ) {
					p->change_type(Float);
					p->rescale_to_avg_std(0, 1);
					p->shift_background(0);
					p->sampling(part->pixel_size);
					if ( mg->ctf )
						img_ctf_apply(p, *(mg->ctf), 1, 0.2, 0, hires, invert, planf, planb);
					translate = prec->image[i].origin() - part->ori;
					origin = part->ori;
					euler = Euler(part->view);
					mat = Matrix3(axis, euler.psi());
//					mat = Matrix3(axis, part->view.angle());
					p->transform(scale, origin, translate, mat, FILL_BACKGROUND, 0);
					n++;
					f += part->fom[0];
					prec->add(i, p);
					delete p;
				}
			}
		}
	}
	
	if ( n ) {
		f /= n;	
		prec->multiply(i, 1.0L/n);
	}
	
	prec->image[i].FOM(f);
	prec->image[i].select(n);
	
	Bparticle*		part_avg = new Bparticle;
	part_avg->id = i+1;
	part_avg->mag = 1;
	
	part_avg->ori = prec->image[i].origin();
	part_avg->fom[0] = f;
	part_avg->sel = n;
	
	return part_avg;
}

long		img_reconstruct_2D_fspace(long i, Bproject* project, Bimage* prec,
				double hires, fft_plan planf, fft_plan planb)
{
	bool			invert(0);

	long			maxrad = ( hires > 0 )? prec->real_size()[0]/hires: prec->sizeX()/2;
	if ( maxrad > prec->sizeX()/2 ) maxrad = prec->sizeX()/2;
	
	long			j, n(0);
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;
	Bimage*			p = NULL;
	Bimage*			psum = new Bimage(Float, TComplex, prec->size(), 1);
	Bimage*			pssq = new Bimage(Float, TSimple, prec->size(), 1);
	psum->origin(prec->image->origin());
	pssq->origin(prec->image->origin());
	psum->sampling(prec->sampling(0));
	pssq->sampling(prec->sampling(0));

	double			f(0), pwr, t(0.5);
	Vector3<double>	origin, translate;
	Vector3<double>	scale(1,1,1);
	Vector3<double> axis(0,0,1);
	Matrix3			mat(1);
	Euler			euler;
	
	long			nmg = project_count_mg_selected(project), m(0);
	
//	cout << "Map: " << i+1 << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			m++;
//			if ( mg->ctf ) cout << "CTF pointer: " << mg->ctf << endl << endl << endl;
			if ( verbose && i < 1 ) {
				cout << mg->id;
				if ( mg->ctf ) cout << tab << mg->ctf->defocus_average();
				cout << tab << m*100.0/nmg << "%" << endl;
			}
			for ( part = mg->part; part; part = part->next ) if ( part->sel == i+1 ) {
//				cout << "Reading particle " << part->id << endl;
				p = particle_read_img(part, 1);
//				cout << "Particle read " << p->file_name() << endl;
				if ( p ) {
					p->change_type(Float);
					p->rescale_to_avg_std(0, 1);
					p->shift_background(0);
					p->sampling(part->pixel_size);
					if ( mg->ctf )
						img_ctf_apply(p, *(mg->ctf), 1, 0.2, 0, hires, invert, planf, planb);
					translate = prec->image[i].origin() - part->ori;
					origin = part->ori;
					euler = Euler(part->view);
					mat = Matrix3(axis, euler.psi());
					p->transform(scale, origin, translate, mat, FILL_BACKGROUND, 0);
					n++;
					f += part->fom[0];
					p->fft(planf);
					psum->add(p);
					p->complex_to_intensities();
					if ( verbose & VERB_FULL )
						cout << "Particle " << part->id << endl;
					pssq->add(p);
//					cout << "data size = " << p->data_size() << endl;
					delete p;
				}
//				if ( verbose )
//					cout << "Particle " << part->id << " done " << endl;
			}
		}
	}

	if ( n ) {
		f /= n;	
		for ( j=0; j<psum->image_size(); ++j ) {
			pwr = psum->complex(j).power();
			psum->set(j, psum->complex(j)/n);
			pssq->set(j, (pwr - (*pssq)[j])/((*pssq)[j] - pwr/n));
		}
		psum->fft_back(planb);
		pssq->origin(Vector3<double>(0,0,0));
		Bimage*			pr = pssq->radial(0, pssq->sizeX()/2, 1, 1);
		if ( verbose ) {
			cout << "Radius\tSSNR" << endl;
			for ( j=0; j<maxrad; ++j )
				cout << j << tab << (*pr)[j] << endl;
			cout << endl;
		}
		for ( j=1; j<maxrad; ++j )
			if ( (*pr)[j] < t ) break;
		f = prec->real_size()[0]/j;
		if ( verbose )
			cout << j << tab << f << tab << (*pr)[j] << endl << endl;
		delete pr;
		prec->replace(i, psum);
		prec->image[i].FOM(f);
		prec->image[i].select(n);
	}
	
	delete psum;
	delete pssq;
		
	return n;
}

Bimage* 	project_reconstruct_2D_fast(Bproject* project, Bstring file_name)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = part_find_first(project);
	Bimage*			p = NULL;
	
	if ( !part ) return p;

	mg = part->mg;
	if ( !mg ) return p;
	
	p = particle_read_img(part, 0);
	if ( !p ) return p;

	Vector3<long>	size(p->size());
	double			px(part->pixel_size[0]), hires(2*px);

	delete p;
	
	// Determine the number of reconstructions to generate
	long 			nmap(0);
	for ( field=project->field; field; field=field->next )
		for ( mg=field->mg; mg; mg=mg->next )
			for ( part=mg->part; part; part=part->next )
				if ( part->sel > nmap ) nmap = part->sel;

	if ( nmap < 1 ) return NULL;
	
	Bimage* 		prec = new Bimage(Float, TSimple, size, nmap);
	prec->origin(prec->size()/2);
	prec->sampling(px, px, 1);

	if ( verbose & VERB_RESULT ) {
		cout << "2D frequency space averaging:" << endl;
		cout << "Size:                           " << prec->size() << endl;
		cout << "Origin:                         " << prec->image->origin() << endl;
		cout << "Sampling:                       " << prec->sampling(0) << " A" << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Number of averages:             " << prec->images() << endl << endl;
	}

	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(nmap, dispatch_get_global_queue(0, 0), ^(size_t i){
		img_reconstruct_2D_fspace(i, project, prec, hires, planf, planb);
		if ( verbose )
			cout << i+1 << tab << prec->image[i].select() << tab << prec->image[i].FOM() << endl;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nmap; i++ ) {
		img_reconstruct_2D_fspace(i, project, prec, hires, planf, planb);
		if ( verbose )
			cout << i+1 << tab << prec->image[i].select() << tab << prec->image[i].FOM() << endl;
	}
#endif

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);
	
	if ( verbose & VERB_DEBUG )
		cout << "done" << endl;
	
	prec->statistics();
	
	particle_kill(project->class_avg);

	project_update_class_averages(project, prec, file_name);
	
	return prec;
}				


/**
@brief 	Transforms 3D maps and calculates an average.
@param 	*project	parameter structure with all parameters.
@param 	selnum		selection number of reconstructions or particles.
@param 	calcfom		flag to calculate FOM block (1=var, 2=std).
@param 	size		size of particles to extract (only when extraction needed).
@return Bimage*		average map with FOM block defined.

	The orientations of the maps in the project must already be specified.
	Each map is transformed and added to an average map.
	A FOM block is optionally calculated with either the variance or
	standard deviation.

**/
Bimage*		project_reconstruct_3D(Bproject* project, long selnum, int calcfom, Vector3<long> size)
{
	int					n(0), img_num;
	Breconstruction*	rec = project->rec;
	Bparticle*			part;
	Bimage*				p = NULL;
	Bimage*				pmap = NULL;;
	
	if ( rec->part ) {
		if ( rec->part->fpart.length() ) p = read_img(rec->part->fpart, 0, 0);
		else if ( rec->fpart.length() ) p = read_img(rec->fpart, 0, 0);
		else pmap = read_img(rec->frec, 1, 0);
	} else p = read_img(rec->frec, 0, 0);
	
	if ( !p && !pmap ) {
		error_show("project_reconstruct_3D", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( p ) {
		size = p->size();
		delete p;
	}
		
	if ( size.volume() <= 0 ) {
		cerr << "Error: The particle box size must be specified!" << endl << endl;
		return NULL;
	}
	
	Bimage*				prec = new Bimage(Float, TSimple, size, 1);
	prec->sampling(rec->voxel_size);
	prec->origin(prec->default_origin());
	float*				fom = NULL;
	if ( calcfom ) {
		prec->next = new Bimage(Float, TSimple, size, 1);
		fom = (float *) prec->next->data_pointer();
	}
	
	long				i;
	Bstring				filename;
	Vector3<double>		start;
	double				std_sum(0);
	
	if ( verbose ) {
		cout << "Calculating an average reconstruction" << endl;
		cout << "Size:                           " << size << endl;
	}
	
	if ( rec->part ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			filename = rec->fpart;
			for ( part = rec->part; part; part = part->next )
					if ( ( selnum < 0 && part->sel ) || ( part->sel == selnum ) ) {
				img_num = part->id - 1;
				if ( part->fpart.length() ) {
					filename = part->fpart;
					img_num = 0;
				}
				if ( filename.length() ) {
					if ( ( p = read_img(filename, 1, img_num) ) == NULL ) {
						error_show("project_reconstruct_3D", __FILE__, __LINE__);
						return NULL;
					}
					if ( verbose ) 
						cout << "Adding particle: " << part->id << " (" << filename << ":" << img_num << ")" << endl;
				} else {
					start = part->loc - part->ori;
					p = pmap->extract(0, start, size);
					if ( verbose ) 
						cout << "Adding particle: " << part->id << endl;
				}
//				p->image->origin(part->ori);
				p->origin(part->ori);
				p->calculate_background();
				std_sum += p->standard_deviation();
				prec->rotate_and_add(p, part->ori, part->view);
				delete p;
				n++;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
				if ( ( selnum < 0 && rec->select ) || ( rec->select == selnum ) ) {
			if ( verbose ) 
				cout << "Adding reconstruction: " << rec->frec << endl;
			if ( ( p = read_img(rec->frec, 1, 0) ) == NULL ) {
				error_show("project_reconstruct_3D", __FILE__, __LINE__);
				return NULL;
			}
//			p->image->origin(rec->origin);
			p->origin(rec->origin);
			p->calculate_background();
			std_sum += p->standard_deviation();
			prec->rotate_and_add(p, rec->origin, rec->view);
			delete p;
			n++;
		}
	}

	if ( pmap ) delete pmap;
	
	if ( verbose )
		cout << "Number of particles:            " << n << endl;
	
	if ( n ) {
		for ( i=0; i<prec->data_size(); i++ ) {
			prec->set(i, (*prec)[i] / n);
			if ( fom ) {
				fom[i] = fom[i]/n - (*prec)[i]*(*prec)[i];
				if ( fom[i] < 0 ) fom[i] = 0;
				if ( calcfom > 1 ) fom[i] = sqrt(fom[i]);
			}
		}
		std_sum /= n;
	}
	
	prec->statistics();
	prec->calculate_background();
	
	if ( verbose )
		cout << "Signal enhancement:             " << 
			(sqrt(n)*prec->standard_deviation()/std_sum - 1)/(sqrt(n) - 1) << endl << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_reconstruct_3D: done!" << endl;
	
	return prec;
}

/**
@brief 	Transforms 3D maps and calculates an average.
@param 	*project	parameter structure with all parameters.
@param 	selnum		selection number of reconstructions or particles.
@param 	size		size of particles to extract (only when extraction needed).
@param 	resolution	maximum reconstruction resolution (in angstrom, default Nyquest).
@return Bimage*		average map with FOM block defined.

	The orientations of the maps in the project must already be specified.
	Each map is transformed and added to an average map.
	A FOM block is optionally calculated with either the variance or
	standard deviation.

**/
Bimage*		project_reconstruct_3D(Bproject* project, long selnum, Vector3<long> size, double resolution)
{
	long				n(0), img_num;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	Bimage*				p = NULL;
	Bimage*				pmap = NULL;
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		part = rec->part;
		if ( part ) break;
	}
	
	if ( !part ) {
		cerr << "Error in project_reconstruct_3D: No particles found!" << endl;
		return NULL;
	}

	Vector3<double>		voxel_size(rec->voxel_size);
	
	if ( part ) {
		if ( part->fpart.length() ) p = read_img(part->fpart, 0, 0);
		else if ( rec->fpart.length() ) p = read_img(rec->fpart, 0, 0);
		else pmap = read_img(rec->frec, 1, 0);
		if ( part->pixel_size[0] > 0 ) voxel_size = part->pixel_size;
	} else p = read_img(rec->frec, 0, 0);
	
	if ( !p && !pmap ) {
		error_show("project_reconstruct_3D", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( p ) {
		size = p->size();
		delete p;
	}
		
	if ( size.volume() <= 0 ) {
		cerr << "Error: The particle box size must be specified!" << endl << endl;
		return NULL;
	}
	
	Bimage*				prec = new Bimage(Float, TComplex, size, 1);
	prec->fourier_type(Standard);
	prec->sampling(voxel_size);
	prec->origin(prec->default_origin());
	prec->next = new Bimage(Float, TSimple, size, 1);
	
	double				threshold;
	Bstring				filename;
	Vector3<double>		start;
	
	if ( verbose ) {
		cout << "Calculating an average reconstruction" << endl;
		cout << "Size:                           " << size << endl;
		cout << "Resolution:                     " << resolution << endl;
	}
	
	if ( part ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			filename = rec->fpart;
			for ( part = rec->part; part; part = part->next )
					if ( ( selnum < 0 && part->sel ) || ( part->sel == selnum ) ) {
				img_num = part->id - 1;
				if ( part->fpart.length() ) {
					filename = part->fpart;
					img_num = 0;
				}
				if ( filename.length() ) {
					if ( ( p = read_img(filename, 1, img_num) ) == NULL ) {
						error_show("project_reconstruct_3D", __FILE__, __LINE__);
						return NULL;
					}
					if ( verbose ) 
						cout << "Adding particle: " << part->id << " (" << filename << ":" << img_num << ")" << endl;
				} else {
					start = part->loc - part->ori;
					p = pmap->extract(0, start, size);
					if ( verbose ) 
						cout << "Adding particle: " << part->id << endl;
				}
				threshold = p->standard_deviation()/1000;
				p->origin(part->ori);
				p->calculate_background();
				p->rotate(prec->image->origin() - part->ori, part->view);
				p->origin(prec->image->origin());
				p->fft();
				p->phase_shift_to_origin();
//				prec->add(p);
				prec->fspace_pack_3D(p, resolution, threshold);
				delete p;
				n++;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
				if ( ( selnum < 0 && rec->select ) || ( rec->select == selnum ) ) {
			if ( verbose ) 
				cout << "Adding reconstruction: " << rec->frec << endl;
			if ( ( p = read_img(rec->frec, 1, 0) ) == NULL ) {
				error_show("project_reconstruct_3D", __FILE__, __LINE__);
				return NULL;
			}
			threshold = p->standard_deviation()/1000;
			p->origin(rec->origin);
			p->calculate_background();
			p->rotate(prec->image->origin() - rec->origin, rec->view);
			p->fft();
			p->phase_shift_to_origin();
//			prec->add(p);
			prec->fspace_pack_3D(p, resolution, threshold);
			delete p;
			n++;
		}
	}

	if ( pmap ) delete pmap;
	
	if ( verbose ) {
		cout << "Number of particles:            " << n << endl;
		cout << "Rescaling average" << endl << endl;
	}
	
	prec->fspace_reconstruction_weigh();

	prec->fspace_reconstruction_stats(resolution);

	prec->phase_shift_to_center();
	prec->origin(prec->default_origin());
	
	prec->fft_back();
	
	prec->statistics();
	prec->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_reconstruct_3D: done!" << endl;
	
	return prec;
}

/**
@brief 	Back projects a set of 2D images into a 3D volume.
@param 	*project 	image processing parameter structure.
@param 	num_select	selection number from the selection column.
@param 	map_size	3-valued vector for the new map size.
@param 	sam			3-value vector for the voxel size in angstrom.
@param 	scale		reconstruction scale.
@param 	resolution	resolution for low-pass filtering.
@param 	planf		2D forward Fourier transform plan.
@param 	planb		2D backward Fourier transform plan.
@return Bimage* 	the new 3D reconstruction map.

	All the information needed to do a 3D reconstruction is passed in through
	an image processing structure. The new 3D volume is initialized.
	Each sub-image in each particle file is read individually and 
	back-projected within the new volume. The orientation parameters from
	the image processing structure is transferred to the sub-image
	structure in the 2D image before calling the function packing one
	image into the volume. The default origin is the center of the image.

**/
Bimage* 	project_back_projection(Bproject* project, long num_select,
				Vector3<long> map_size, Vector3<double> sam, double scale,
				double resolution, fft_plan planf, fft_plan planb)
{
	Bfield*			field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bimage*			p = NULL;

	if ( map_size.volume() <= 0 ) {
		error_show("Error in project_back_projection", __FILE__, __LINE__);
		cerr << "No reconstruction size set!" << endl << endl;
		return NULL;
	}

	for ( field = project->field; field && part == NULL; field = field->next )
		for ( mg = field->mg; mg && part == NULL; mg = mg->next )
			if ( ( part = mg->part ) ) break;
	
	if ( !mg || !part ) {
		cerr << "Error: No particles found!" << endl;
		return NULL;
	}
	
	// Initialize the new 3D volume
	Bimage*			pmap = new Bimage(Float, TSimple, map_size, 1);
	
	if ( sam.volume() > 0 ) pmap->sampling(sam);
	else if ( mg ) pmap->sampling(part->pixel_size);

	pmap->origin(pmap->default_origin());
	
	pmap->check();

	long			nsel = project_count_mg_part_selected(project, num_select);
	long 			psel, nrec(0);
	
	if ( verbose & VERB_FULL ) {
		cout << "3D backprojection:" << endl;
		cout << "Map size:                       " << pmap->size() << endl;
		cout << "Resolution:                     " << resolution << " A" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Sampling:                       " << pmap->sampling(0) << endl;
		cout << "Number of selected particles:   " << nsel << endl << endl;
	}

	// Do the back projection
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				psel = 0;
				if ( part->sel ) {
					if ( num_select < 0 ) {
						psel = 1;
					} else if ( part->sel == num_select ) {
						psel = 1;
					}
				}
				if ( psel ) {	// Use only selected images
					p = particle_read_img(part, 1);
//					if ( ( p = read_img(mg->fpart, 1, part->id-1) ) == NULL ) {
					if ( !p ) {
						error_show("project_back_projection", __FILE__, __LINE__);
						return NULL;
					}
				
					p->sampling(part->pixel_size);
					
					if ( part->ori[0] <= 0 ) {
						if ( p->image->origin()[0] > 0 ) part->ori[0] = p->image->origin()[0];
						else part->ori[0] = pmap->image->origin()[0];
					}
					
					if ( part->ori[1] <= 0 ) {
						if ( p->image->origin()[1] > 0 ) part->ori[1] = p->image->origin()[1];
						else part->ori[1] = pmap->image->origin()[1];
					}
					
					part->ori[2] = 0;
					
					p->origin(part->ori);
					
					p->image->view(part->view);
					
					p->image->magnification(part->mag / scale);
					
					pmap->back_project(p, resolution, 0, planf, planb);
					
					delete p;
					
					nrec++;
					
					if ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) )
						cout << "Complete:                       " << nrec*100.0/nsel << "%\r" << flush;
				}
			}
		}
	}
	
	if ( verbose ) {
		cout << endl << "Particles used:                 " << nrec << endl << endl;
	}

	return pmap;
}

/**
@brief 	Accumulates multiple backprojections.
@param 	**pacc		array of partial maps with linked weight maps.
@param 	imap		which output map to weigh.
@param 	nmaps		number of output maps (1,2,3).
@param 	nthreads	number of threads per map (= number of partial maps).
@return int 		0.

**/
Bimage*		img_backprojection_accumulate(Bimage** pacc, int imap,
				int nmaps, int nthreads)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_backprojection_accumulate:" << endl;

	int					maps_per_class = nthreads;
	int					iclass = imap/nmaps;

	Bimage*				prec = pacc[0]->copy_header();
	prec->data_alloc_and_clear();

	long				i, j, k;
	long   				datasize = (long) prec->size().volume();
	
	long				st(0), inc(1);
	if ( nmaps == 2 ) {
		st = imap%2;
		inc = 2;
	}
	if ( nmaps == 3 && imap%3 ) {
		st = imap%3 - 1;
		inc = 2;
	}

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_backprojection_weigh: nmaps=" << nmaps << " imap=" << imap << endl;
		cout << "DEBUG img_backprojection_weigh: maps_per_class=" << maps_per_class << " iclass=" << iclass << endl;
		cout << "DEBUG img_backprojection_weigh: st=" << st << " inc=" << inc << endl;
	}
	
	for ( i=iclass*maps_per_class+st, j=st; j<maps_per_class; i+=inc, j+=inc ) {
//		cout << "Doing map " << i << endl;
		for ( k=0; k<datasize; k++ )
			prec->add(k, (*pacc[i])[k]);
	}
	
	prec->statistics();
	prec->calculate_background();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_backprojection_weigh: Accumulation done" << endl;

	return prec;
}

