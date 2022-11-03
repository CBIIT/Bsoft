/**
@file	mg_tomo_resol.cpp
@brief	Functions to assess the resoltion of a tomographic tilt series
@author	Bernard Heymann
@date	Created: 20031205
@date	Modified: 20190221
**/

#include "mg_tomo_resol.h"
#include "mg_tomo_rec.h"
#include "mg_reconstruct.h"
#include "rwimg.h"
#include "mg_ctf.h"
#include "Complex.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bplot*		img_calculate_FRC_curve(vector<Bimage*> parr, double hi_res, double sampling_ratio);
double		plot_calculate_resolution(Bplot* plot, double cutoff);
int			plot_resolution(Bplot* plot, double hi_res, Bstring& psfile, Bstring& title);

int			img_pack_2D_into_central_section_old(Bimage* p, Bimage* prec, Bimage* prec2,
				long ft_size, int zsize, double scale, double hi_res, 
				Matrix3 matr, Matrix3 mat, int inplane);


/**
@brief 	Estimating the resolution of one micrograph in an aligned tilt series.
@param 	*project 		image processing parameter structure.
@param 	micrograph_id	micrograph number to use for resolution test.
@param 	hi_res			high resolution limit.
@param	sampling_ratio	ratio for averaging window.
@param 	scale			reconstruction scale.
@param 	size			reconstruction size.
@param 	fast_angle		angle to select micrographs for reconstruction.
@param 	action			flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@return	vector<Bimage*>	array of three images, length of zero on error.

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
	For voxels with only one data pixel contributing to it, FOM(0).
	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.
	Three images are returned: micrograph, full reconstruction, leave-one-out reconstruction.

**/
vector<Bimage*>	mg_tomo_res_reconstruct(Bproject* project, int micrograph_id, double hi_res,
				double sampling_ratio, double scale, Vector3<long> size,
				double fast_angle, int action, double wiener)
{
	vector<Bimage*>	parr;
	long   			i;
	int				pad_factor(0);
	Vector3<long> 	tile_size(128,128,1);

	Bfield*			field = project->field;
	Bmicrograph*	mg_res;
	Bmicrograph*	mg;
	
	for ( mg_res = field->mg, i=0; mg_res && i < micrograph_id; mg_res = mg_res->next, i++ ) ;

	if ( !mg_res ) {
		error_show("Error in mg_tomo_res_reconstruct", __FILE__, __LINE__);
		cerr << "Micrograph number " << micrograph_id << " not found!" << endl << endl;
		return parr;
	}

	if ( !mg_res->select ) {
		if ( verbose & VERB_FULL )
			cout << "Skipping micrograph number %d!" << micrograph_id << endl << endl;
		return parr;
	}
	
	if ( mg_res->fmg.length() < 0 ) {
		error_show("Error in mg_tomo_res_reconstruct", __FILE__, __LINE__);
		cerr << "No micrograph image file name given for micrograph " << mg_res->id << endl;
		return parr;
	}
	
	if ( hi_res < 2*mg_res->pixel_size[0] )
		hi_res = 2*mg_res->pixel_size[0];
	
	Bimage*			pt;
	Bimage* 		p = NULL;
	if ( pad_factor < 1 && mg_res->fft.length() > 0 ) p = read_img(mg_res->fft, 1, mg_res->img_num);
	else p = read_img(mg_res->fmg, 1, mg_res->img_num);
	
	p->sampling(mg_res->pixel_size);
	p->origin(mg_res->origin[0], mg_res->origin[1], 0.0);
	
	if ( size[2] < 1 ) size[2] = p->sizeX();
	if ( size.volume() < 1 ) size = Vector3<long>(p->size());
	if ( size[2] > p->sizeX() ) size[2] = p->sizeX();
	if ( size[2] > p->sizeY() ) size[2] = p->sizeY();
	
	long 	ft_size = findNextPowerOf((int)(pad_factor*size[0]*scale), 2);
	if ( ft_size < p->sizeX() ) ft_size = p->sizeX();
	if ( ft_size < p->sizeY() ) ft_size = p->sizeY();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_tomo_res_reconstruct: size=" << size << " ft_size=" << ft_size << endl;

	if ( action )
		img_ttf_apply(p, *(mg_res->ctf), action, wiener,
			tile_size, mg_res->tilt_angle, mg_res->tilt_axis, 0, 0, 0);
	
	// Do transformation if not already done
	if ( p->compound_type() == TSimple ) {
		if ( pad_factor ) {
			ft_size = findNextPowerOf((int)(pad_factor*size[0]*scale), 2);
			p->calculate_background();
		}
		p->pad(ft_size, FILL_AVERAGE, p->average());
		p->fft();
		p->phase_shift_to_origin();
	}
	
	// Header structure for the final reconstruction
	Bimage* 		prec = new Bimage(Float, TComplex, ft_size, ft_size, 1, 1);
	prec->sampling(mg_res->pixel_size);
	prec->origin(prec->size()/2);

	Bimage* 		prec2 = prec->copy();
	
	prec->next = new Bimage(Float, TSimple, prec->size(), 1);
	prec2->next = new Bimage(Float, TSimple, prec2->size(), 1);

	float*			fom = (float *) prec->next->data_pointer();
	float*			fom2 = (float *) prec2->next->data_pointer();

	long   			data_size = (long) prec->size().volume();
	
//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "NLOO resolution test:" << endl;
		cout << "Micrograph:                     " << mg_res->id << " (" << mg_res->img_num << ")" << endl;
		cout << "Tilt angle:                     " << mg_res->tilt_angle*180.0/M_PI << endl;
		if ( fast_angle )
			cout << "Maximum angular offset:         " << fast_angle*180.0/M_PI << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Reconstruction size:            " << size << endl;
		cout << "Resolution limit:               " << hi_res << " A" << endl;
		cout << "Sampling ratio:                 " << sampling_ratio << endl;
		if ( action )
			cout << "CTF correction:                 " << action << " (wiener=" << wiener << ")" << endl;
	}

	long 			nrec(0);
	long			nmg = count_list((char *) field->mg);
	
	for ( mg=field->mg; mg; mg=mg->next ) if ( mg->select )
			if ( fabs(mg->tilt_angle - mg_res->tilt_angle) <= fast_angle ) {
		if ( verbose )
			cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
		pt = NULL;
		if ( mg->fft.length() > 0 )
			pt = read_img(mg->fft, 1, 0);
		if ( !pt ) {
			if ( ( pt = read_img(mg->fmg, 1, mg->img_num) ) == NULL ) {
				error_show("mg_tomo_res_reconstruct", __FILE__, __LINE__);
				return parr;
			}
		}

		if ( mg == mg_res ) mg->fom = p->average();
	
		pt->origin(mg->origin[0], mg->origin[1], 0.0);

		if ( action )
			img_ttf_apply(pt, *(mg->ctf), action, wiener,
				tile_size, mg->tilt_angle, mg->tilt_axis, 0, 0, 0);
		
		if ( pt->compound_type() == TSimple ) {
			pt->pad(ft_size, FILL_AVERAGE, pt->average());
			pt->fft();
			pt->phase_shift_to_origin();
		}

		img_pack_2D_into_central_section_old(pt, prec, prec2, ft_size, size[2], 1, hi_res,
			mg_res->matrix, mg->matrix, (mg==mg_res));
//		img_pack_2D_into_central_section(pt, prec, prec2, ft_size, 1, hi_res,
//			mg_res->matrix, mg->matrix, (mg==mg_res));
				
		delete pt;
				
		nrec++;
				
		if ( verbose & ( VERB_TIME | VERB_PROCESS ) ) {
			cout << "Complete:                       " << nrec*100.0/nmg << "%\r";
			cout.flush();
		}
		
	}
	
	if ( nrec < 1 ) {
		delete prec;
		delete prec2;
		error_show("Error in mg_tomo_res_reconstruct: No micrographs used in this reconstruction!", __FILE__, __LINE__);
		return parr;
	}

	if ( verbose )
		cout << "Weighing reconstruction.          " << endl;
	
	for ( i=0; i<data_size; i++ ) {
		if ( fom[i] ) prec->set(i, prec->complex(i) / fom[i]);
		if ( fom2[i] ) prec2->set(i, prec2->complex(i) / fom2[i]);
		fom[i] = fom[i]/(fom[i] + 1);
		fom2[i] = fom2[i]/(fom2[i] + 1);
	}
	
	parr.push_back(p);
	parr.push_back(prec);
	parr.push_back(prec2);

	return parr;
}

/**
@brief 	Estimating the resolution of one micrograph in an aligned tilt series.  
@param 	*project 		image processing parameter structure.
@param 	micrograph_id	micrograph number to use for resolution test.
@param 	hi_res			high resolution limit.
@param	sampling_ratio	ratio for averaging window. 
@param 	scale			reconstruction scale.
@param 	size			reconstruction size.
@param 	fast_angle		angle to select micrographs for reconstruction.
@param 	action			flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	cutoff			FRC cutoff.
@param 	&psfile			postscript output file name.
@return	Bimage*			micrograph reconstruction, NULL on error.

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
	For voxels with only one data pixel contributing to it, FOM(0).
	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.

**/
Bimage*		mg_tomo_resolution(Bproject* project, int micrograph_id, double hi_res, 
				double sampling_ratio, double scale, Vector3<long> size,
				double fast_angle, int action, double wiener, double cutoff, Bstring& psfile)
{
	long   			i;
	Bfield*			field = project->field;
	Bmicrograph*	mg_res;

	for ( mg_res = field->mg, i=0; mg_res && i < micrograph_id; mg_res = mg_res->next, i++ ) ;

	if ( !mg_res ) {
		error_show("Error in mg_tomo_resolution", __FILE__, __LINE__);
		cerr << "Micrograph number " << micrograph_id << " not found!" << endl << endl;
		return NULL;
	}

	vector<Bimage*>	parr = mg_tomo_res_reconstruct(project, micrograph_id, hi_res,
				sampling_ratio, scale, size, fast_angle, action, wiener);

	Bplot*			plot = img_calculate_FRC_curve(parr, hi_res, sampling_ratio);
	double			resolution = plot_calculate_resolution(plot, cutoff);
	mg_res->fom = resolution;

	Bstring			filename(parr[0]->file_name());
	Bstring			title = "Resolution for " + filename + ": image " + Bstring(mg_res->img_num+1, "%d");
	if ( psfile.length() > 0 ) plot_resolution(plot, hi_res, psfile, title);
	
	if ( verbose )
		cout << "Resolution:                     " << setprecision(2) << resolution << endl << endl;
	

	delete parr[0];
	delete parr[2];
	delete plot;
	
	return parr[1];
}

int			plot_nloo3d(Bplot* plot)
{
	long		i, j, k, l, m, nr(plot->rows());
	
	for ( i=nr, j=2*nr, k=3*nr, l=4*nr, m=5*nr; i<2*nr; i++, j++, k++, l++, m++ )
		(*plot)[i] = ((*plot)[m]/(*plot)[l])*sqrt((*plot)[j]/(*plot)[k]);

	return 0;
}

/**
@brief 	Estimating the resolution of the tomogram from an aligned tilt series.
@param 	*project 		image processing parameter structure.
@param 	hi_res			high resolution limit.
@param	sampling_ratio	ratio for averaging window.
@param 	scale			reconstruction scale.
@param 	size			reconstruction size.
@param 	fast_angle		angle to select micrographs for reconstruction.
@param 	action			flag to apply CTF to projections.
@param 	wiener			Wiener factor.
@param 	cutoff			FRC cutoff.
@return	vector<Bplot*>		Two plots: Tilt-resolution and NLOO-3D.

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
	For voxels with only one data pixel contributing to it, FOM(0).
	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.

**/
vector<Bplot*>	project_tomo_resolution(Bproject* project, double hi_res,
				double sampling_ratio, double scale, Vector3<long> size,
				double fast_angle, int action, double wiener, double cutoff)
{
	long			nmg, i;
	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
//	Bstring			title("Micrograph resolution for a tilt series");
	Bstring			label;

	for ( nmg=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) nmg++;
	
/*	Bplot*			plot = new Bplot(1, nmg, 2);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(2);
	for ( i=0; i<2; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("TiltAngle");
	plot->page(0).column(1).label("NLOO");
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).color(0,0,0);
*/
	Bplot*			mgplot = NULL;
	Bplot*			tomoplot = NULL;
	
	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, i++ ) if ( mg->select ) {
			vector<Bimage*>	parr = mg_tomo_res_reconstruct(project, i, hi_res,
				sampling_ratio, scale, size, fast_angle, action, wiener);
			if ( verbose )
				cout << "Intensity:                      " << setprecision(2) << mg->fom << endl << endl;
			mgplot = img_calculate_FRC_curve(parr, hi_res, sampling_ratio);
			mg->fom = plot_calculate_resolution(mgplot, cutoff);
			if ( verbose )
				cout << "Resolution:                     " << setprecision(2) << mg->fom << endl << endl;
			delete parr[0];
			delete parr[1];
			delete parr[2];
			if ( tomoplot ) {
				tomoplot->add(mgplot, 1);
				delete mgplot;
			} else {
				tomoplot = mgplot;
			}
//			(*plot)[i] = mg->tilt_angle*180.0/M_PI;
//			(*plot)[i+nmg] = mg->fom;
		}
	}
	
	plot_nloo3d(tomoplot);

//	plot->page(0).axis(1).label("Tilt Angle (degrees)");
//	plot->page(0).axis(3).label("Resolution (A)");

	Bplot*			tiltplot = plot_tilt_resolution(project);
	label = Bstring(cutoff, "NLOO cutoff: %g");
	tiltplot->page(0).add_text(label);

//	Bstring			filename(parr[0]->file_name());
//	Bstring			title = "Resolution for " + filename + ": image " + Bstring(mg_res->img_num+1, "%d");
//	if ( psfile.length() > 0 ) plot_resolution(plot, hi_res, psfile, title);

	vector<Bplot*>	plot;
	plot.push_back(tiltplot);
	plot.push_back(tomoplot);

	return plot;
}

/**
@brief 	Estimates the resolution for each particle image in each micrograph.
@param 	*project		micrograph project.
@param 	hi_res			high resolution limit for resolution estimation.
@param	sampling_ratio	ratio for averaging window. 
@param 	fast_angle		angle to select micrographs for reconstruction.
@param 	cutoff			FRC cutoff to use.
@return Bplot*			plot with average particle resolutions.

	Requires the particles to be defined in all micrographs.
	The NLOO algorithm is used for each particle.

**/
Bplot*		project_tomo_particle_resolution(Bproject* project, double hi_res, 
				double sampling_ratio, double fast_angle, double cutoff)
{
	long				i, j, nmg(0);
	long				npart(0);
	int					pad_factor(0);
	double				scale(1), resolution;
	Bfield*				field, *field_res;
	Bmicrograph*		mg, *mg_res;
	Breconstruction*	rec;
	Bparticle*			partlist = NULL;
	Bparticle*			recpart = NULL;
	Bparticle*			part = NULL, *part_res = NULL;
	Bimage*				prec, *prec2;
	Bimage*				pmg;
	Bimage*				ppart, *ppart_res = NULL;
	Bplot*				plot = NULL;
	Bstring				psfile;

	for ( rec = project->rec; rec; rec = rec->next )
		if ( rec->part ) {
			partlist = rec->part;
			break;
		}
	
	if ( !partlist ) {
		cerr << "No particle list found!" << endl << endl;
		return NULL;
	}
	
	for ( nmg=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) nmg++;
	
	if ( !nmg ) {
		cerr << "No micrographs found!" << endl << endl;
		return NULL;
	}
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mg && mg->part ) break;
		if ( mg && mg->part ) break;
	}
	
	if ( !mg || !mg->part ) {
		cerr << "No micrograph particles found!" << endl << endl;
		return NULL;
	}
	
	double*				res_avg = new double[nmg];
	double*				res_std = new double[nmg];
	for ( i=0; i<nmg; i++ ) res_avg[i] = res_std[i] = 0;
	
	Vector3<long>		boxsize(mg->box_size);

	if ( hi_res < 2*mg->pixel_size[0] ) hi_res = 2*mg->pixel_size[0];

	Vector3<double>		realsize(boxsize[0]*mg->pixel_size[0], boxsize[1]*mg->pixel_size[1], 1);
	double				rad_scale = realsize[0];
	if ( rad_scale < realsize[1] ) rad_scale = realsize[1];
	if ( rad_scale < realsize[2] ) rad_scale = realsize[2];	
	if ( boxsize[0] > 256 ) rad_scale *= 256.0/boxsize[0];
	long		maxrad = (long) (rad_scale/hi_res);
	
	int					ft_size = boxsize[0];
	if ( pad_factor ) ft_size = part_ft_size(boxsize[0], scale, pad_factor);
	
	fft_plan			planp = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);

	if ( verbose ) {
		cout << "Calculating individual particle resolutions:" << endl;
		cout << "Number of micrographs:          " << nmg << endl;
		cout << "Particle size:                  " << boxsize << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "High resolution limit:          " << hi_res << " A" << endl;
		cout << "FRC cutoff:                     " << cutoff << endl;
		cout << endl;
	}

	double			s, sp, FRC, FRCp, fraction;
	double*			sum_r_2 = new double[maxrad] ;
	double*			sum_r2_2 = new double[maxrad] ;
	double*			sum_i_r = new double[maxrad] ;
	double*			sum_i_r2 = new double[maxrad] ;
	
	for ( npart=0, recpart = partlist; recpart; recpart = recpart->next ) {
		if ( verbose ) {
			cout << "Resolution for particle " << recpart->id << ":" << endl;
			cout << "#mg\tAngle\tNLOO(" << cutoff << ")" << endl;
		}
		
		for ( j=0; j<maxrad; j++ ) sum_r_2[j] = sum_r2_2[j] = sum_i_r[j] = sum_i_r2[j] = 0;
		
		for ( i=0, field_res = project->field; field_res; field_res = field_res->next ) if ( field_res->select ) {
			for ( mg_res = field_res->mg; mg_res; mg_res = mg_res->next ) if ( mg_res->select ) {
				for ( part_res = mg_res->part; part_res; part_res = part_res->next )
					if ( recpart->id == part_res->id ) break;

				if ( part_res ) {

					prec = new Bimage(Float, TComplex, ft_size, ft_size, 1, 1);
					prec->sampling(mg_res->pixel_size);
					prec->origin(prec->size()/2);

					prec2 = prec->copy();

					prec->next = new Bimage(Float, TSimple, prec->size(), 1);
					prec2->next = new Bimage(Float, TSimple, prec2->size(), 1);

					for ( field = project->field; field; field = field->next ) if ( field->select ) {
						for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select )
								if ( fabs(mg->tilt_angle - mg_res->tilt_angle) <= fast_angle ) {
							for ( part = mg->part; part; part = part->next )
								if ( recpart->id == part->id ) break;
							
							if ( part ) {
								if ( ( pmg = read_img(mg->fmg, 1, mg->img_num) ) == NULL ) {
									error_show("project_tomo_particle_resolution", __FILE__, __LINE__);
									return NULL;
								}
								ppart = pmg->extract(0, part->loc - boxsize/2, boxsize);
								delete pmg;
								ppart->sampling(part->pixel_size);
//								ppart->image->origin(part->ori);
								ppart->origin(part->ori);
								ppart->image->view(part->view);
								ppart->fft(planp);
								ppart->phase_shift_to_origin();
								img_pack_2D_into_central_section_old(ppart, prec, prec2,
									ft_size, ppart->sizeX(), scale, hi_res,
									mg_res->matrix, mg->matrix, (mg==mg_res));
								if ( mg == mg_res ) ppart_res = ppart;
								else delete ppart;
							}
						}
					}
					vector<Bimage*>		parr;
					parr.push_back(ppart_res);
					parr.push_back(prec);
					parr.push_back(prec2);
					plot = img_calculate_FRC_curve(parr, hi_res, sampling_ratio);
					resolution = plot_calculate_resolution(plot, cutoff);
					delete ppart_res;
					delete prec;
					delete prec2;
					for ( j=0; j<maxrad; j++ ) {
						sum_r_2[j] += (*plot)[2*maxrad+j];
						sum_r2_2[j] += (*plot)[3*maxrad+j];
						sum_i_r[j] += (*plot)[4*maxrad+j];
						sum_i_r2[j] += (*plot)[5*maxrad+j];
					}
					delete plot;
					part_res->fom[0] = resolution;
					res_avg[i] += resolution;
					res_std[i] += resolution*resolution;
					i++;
					if ( verbose )
						cout << i << tab << setprecision(2) << mg_res->tilt_angle*180.0/M_PI << tab << resolution << endl;
				}
			}
		}
		npart++;
		if ( verbose ) {
			cout << endl;
			cout << "3D resolution: " << endl;
			cout << "s(1/A)\tRes(A)\tNLOO3D" << endl;
			sp = 0;
			FRCp = 1;
			resolution = 1/rad_scale;
			for ( j=1; j<maxrad; j++ ) {
				s = j/rad_scale;
				if ( fabs(sum_i_r2[j]) < SMALLFLOAT || fabs(sum_i_r[j]) < SMALLFLOAT )
					FRC = 0;
				else
					FRC = (sum_i_r2[j]*sqrt(sum_r_2[j]))/(sum_i_r[j]*sqrt(sum_r2_2[j]));
				if ( FRC < 0 ) FRC = 0;
				if ( FRC > cutoff ) resolution = 1/s;
				if ( FRCp > cutoff && FRC <= cutoff ) {
					fraction = (FRCp - cutoff)/(FRCp - FRC);
					if ( sp ) resolution = (1-fraction)/sp + fraction/s;
				}
				sp = s;
				FRCp = FRC;
				cout << setprecision(4) << s << tab << setprecision(1) << 1/s << tab << setprecision(4) << FRC << endl;
			}
			cout << "NLOO3D resolution: " << resolution << endl << endl;
		}
	}

	delete[] sum_r_2;
	delete[] sum_r2_2;
	delete[] sum_i_r;
	delete[] sum_i_r2;

	fft_destroy_plan(planp);
	
	Bstring			title("Particle resolution for a tilt series");
	plot = new Bplot(1, nmg, 3);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(3);
	for ( i=0; i<3; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("TiltAngle");
	plot->page(0).column(1).label("NLOOavg");
	plot->page(0).column(2).label("NLOOstd");
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(1/hi_res);
//	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).column(1).type(2);
	plot->page(0).column(2).type(3);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(2).axis(3);
	plot->page(0).column(1).color(0,0,0);
	plot->page(0).column(2).color(1,0,0);
	
	Bstring			label = Bstring(cutoff, "NLOO cutoff: %g");
	plot->page(0).add_text(label);
	label = Bstring(boxsize[0], "Box size: %g") + Bstring(boxsize[1], " x %g");
	plot->page(0).add_text(label);
	
	double			min(1e37), max(0);
	
	if ( verbose ) {
		cout << "Average resolution for " << npart << " particles:" << endl;
		cout << "#mg\tAngle\tFRCavg\tFRCstd" << endl;
	}
	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, i++ ) {
			res_avg[i] /= npart;
			res_std[i] = res_std[i]/npart - res_avg[i]*res_avg[i];
			if ( res_std[i] > 0 ) res_std[i] = sqrt(res_std[i]);
			else res_std[i] = 0;
			if ( min > res_avg[i] - res_std[i] ) min = res_avg[i] - res_std[i];
			if ( max < res_avg[i] + res_std[i] ) max = res_avg[i] + res_std[i];
			(*plot)[i] = mg->tilt_angle*180.0/M_PI;
			(*plot)[i+nmg] = res_avg[i];
			(*plot)[i+2*nmg] = res_std[i];
			if ( verbose )
				cout << i+1 << tab << setprecision(2) << mg->tilt_angle*180.0/M_PI << tab << 
					res_avg[i] << tab << res_std[i] << endl;
		}
	}
	if ( verbose )
		cout << endl;
	
	delete[] res_avg;
	delete[] res_std;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_tomo_particle_resolution: min=" << min << " max=" << max << endl;

	plot->page(0).axis(1).label("Tilt Angle (degrees)");
	plot->page(0).axis(3).label("Resolution (A)");
	plot->page(0).axis(3).min(min);
	plot->page(0).axis(3).max(max);
	
	return plot;
}

/**
@brief 	Packs a 2D Fourier transform into a 3D reciprocal space volume.  
@param 	*p			2D Fourier transform.
@param 	*prec		3D central section.
@param 	*prec2		3D central section with in-plane micrograph omitted.
@param 	ft_size		Fourier transform size.
@param 	scale		reconstruction scale.
@param 	hi_res		high resolution limit.
@param 	matr		in plane or reference matrix.
@param 	mat			matrix of image being packed.
@param 	inplane		flag to indicate an in plane image.
@return	long		number voxels packed.

	The rotation matrix is used to determine the plane in reciprocal space
	to which the 2D transform data is added. The map is assumed to be cubic
	and the 2D transform square. The orientation parameters must be written
	into the image structure. The real space 2D image must be supplied.
	This is then padded to more than twice its original size, fourier
	transformed, and packed into the 3D reciprocal space block.

**/
long		img_pack_2D_into_central_section(Bimage* p, Bimage* prec, Bimage* prec2,
				long ft_size, double scale, double hi_res,
				Matrix3 matr, Matrix3 mat, int inplane)
{
	if ( !p->data_pointer() ) return 0;
	
	prec->check_resolution(hi_res);
	
	if ( verbose & VERB_FULL )
		cout << "Packing an image into reciprocal space up to " << hi_res << " A resolution" << endl << endl;
	
	long 			i, j, x, y, kx, ky;
	long 			ix, iy, ixx, iyy, nover(0);
	double 			w, dx, dy, d2;
	Vector3<double>	m, d(0,0,1), iv;
	
//	Vector3<double> vscale(prec->sizeX()/(scale*ft_size), prec->sizeY()/(scale*ft_size), 1/scale);

	float*			fom = (float *) prec->next->data_pointer();
	float*			fom2 = (float *) prec2->next->data_pointer();
	
	int				notinplane = 1 - inplane;		// Flag to indicate if this image corresponds to testing image

	if ( verbose & VERB_FULL & inplane ) cout << "In plane image" << endl;
	if ( verbose & VERB_FULL )
		cout << mat << endl;
	
	double			maxrad = p->sampling(0)[0]/hi_res;
	double			maxrad_sq = maxrad*maxrad;	// Maximum radius squared on map scale
	double			ymin = floor(-maxrad*p->sizeY());
	double			ymax = -ymin;
	double			xmin = floor(-maxrad*p->sizeX());
	double			xmax = -xmin;

	// Vector to calculate z-shift
	Vector3<double>	vn(matr[2][0]/matr[2][2], matr[2][1]/matr[2][2], 0);

	matr = matr.transpose();
//	matr = vscale * matr;
	
	for ( iv[1]=ymin; iv[1]<=ymax; iv[1]+=1 ) {
		if ( iv[1] >= 0 ) y = (long)iv[1];
		else y = (long)(iv[1] + p->sizeY());
		dy = iv[1]/p->sizeY();
		for ( iv[0]=xmin; iv[0]<=xmax; iv[0]+=1 ) {
			if ( iv[0] >= 0 ) x = (long)iv[0];
			else x = (long)(iv[0] + p->sizeX());
			dx = iv[0]/p->sizeX();
			d2 = dx*dx + dy*dy;
			if ( d2 <= maxrad_sq ) {
				m = mat * iv;
				m[2] = m.scalar(vn);
				m = matr * m;
				d[2] = fabs(m[2]);
				if ( d[2] <= 1 ) {
					ix = (long) floor(m[0]);
					iy = (long) floor(m[1]);
					d[0] = m[0] - ix;
					d[1] = m[1] - iy;
					i = y*p->sizeX() + x;
					for ( ky=0; ky<2; ky++ ) {
						iyy = iy + ky;
						if ( iyy < 0 ) iyy += prec->sizeY();
						d[1] = 1 - d[1];
						for ( kx=0; kx<2; kx++ ) {
							ixx = ix + kx;
							if ( ixx < 0 ) ixx += prec->sizeX();
							d[0] = 1 - d[0];
							j = iyy*prec->sizeX() + ixx;
							w = d.volume();
							fom[j] += w;
							prec->add(j, p->complex(i) * w);
							if ( notinplane ) {
								fom2[j] += w;
								prec2->add(j, p->complex(i) * w);
							}
						}
					}
					nover++;
				}
			}
		}
	}

	return nover;
}

int			img_pack_2D_into_central_section_old(Bimage* p, Bimage* prec, Bimage* prec2,
				long ft_size, int zsize, double scale, double hi_res, 
				Matrix3 matr, Matrix3 mat, int inplane)
{
	if ( !p->data_pointer() ) return 0;
	
	prec->check_resolution(hi_res);
	
	if ( verbose & VERB_FULL )
		cout << "Packing an image into reciprocal space up to " << hi_res << " A resolution" << endl << endl;
	
	long 			i, j, x, y, kx, ky;
	long 			ix, iy, ixx, iyy;
	double 			w, d2;
	Vector3<double>	m, d, iv;
	
	Vector3<double>	invsize(1.0/prec->sizeX(), 1.0/prec->sizeY(), 1.0/zsize);
	Vector3<double> vscale(prec->sizeX()/(scale*ft_size), prec->sizeY()/(scale*ft_size), zsize/(scale*ft_size));

	float*			fom = (float *) prec->next->data_pointer();
	float*			fom2 = (float *) prec2->next->data_pointer();
	
	int				notinplane = 1 - inplane;		// Flag to indicate if this image corresponds to testing image

	if ( verbose & VERB_FULL & inplane ) cout << "In plane image" << endl;
	if ( verbose & VERB_FULL )
		cout << mat << endl;
	
	double			maxrad = prec->sampling(0)[0]/hi_res;
	double			maxrad_sq = maxrad*maxrad;	// Maximum radius squared on map scale
	double			ymin = floor(-maxrad*p->sizeY());
	double			ymax = -ymin;
	double			xmin = floor(-maxrad*p->sizeX());
	double			xmax = -xmin;

	Vector3<double>	vt(mat[0][2], mat[1][2], 0);
	if ( vt[0] || vt[1] ) vt.normalize();
	vt *= 1 - mat[2][2]/matr[2][2];
//	if ( fabs(mat[2]) > fabs(mat[5]) ) { if ( mat[2] > 0 ) vt *= -1; }
//	else { if ( mat[5] > 0 ) vt *= -1; }
//	if ( tilt_angle < 0 ) vt *= -1;
	if ( mat[0][2] > 0 ) vt *= -1;
//	vt = 0;

//	cout << mat << endl;
	
	mat = mat.transpose();
	mat = vscale * mat;
	mat = mat * matr;
	
	for ( iv[1]=ymin; iv[1]<=ymax; iv[1]+=1 ) {
		if ( iv[1] >= 0 ) y = (long)iv[1];
		else y = (long)(p->sizeY() + iv[1]);
		for ( iv[0]=xmin; iv[0]<=xmax; iv[0]+=1 ) {
			if ( iv[0] >= 0 ) x = (long)iv[0];
			else x = (long)(p->sizeX() + iv[0]);
			m = mat * iv;
			m += iv * vt;
			d[2] = fabs(m[2]);
			if ( d[2] <= 1 ) {
//				d[2] = 1 - d[2];
				d[2] = 1;
				d2 = (m * invsize).length2();
				if ( d2 <= maxrad_sq ) {
					ix = (long) floor(m[0]);
					iy = (long) floor(m[1]);
					d[0] = m[0] - ix;
					d[1] = m[1] - iy;
					i = y*p->sizeX() + x;
					for ( ky=0; ky<2; ky++ ) {
						iyy = iy + ky;
						if ( iyy < 0 ) iyy += prec->sizeY();
						d[1] = 1 - d[1];
						for ( kx=0; kx<2; kx++ ) {
							ixx = ix + kx;
							if ( ixx < 0 ) ixx += prec->sizeX();
							d[0] = 1 - d[0];
							j = iyy*prec->sizeX() + ixx;
							w = d.volume();
							fom[j] += w;
							prec->add(j, p->complex(i) * w);
							if ( notinplane ) {
								fom2[j] += w;
								prec2->add(j, p->complex(i) * w);
							}
						}
					}
				}
			}
		}
	}

	return 0;
}


/**
@brief 	Calculates the FRC/FSC for a micrograph or particle image in a series.  
@param 	parr			image array with 2D Fourier transform, 3D central section, 3D LOO central section
@param 	hi_res			high resolution limit.
@param	sampling_ratio	ratio for averaging window. 
@return	Bplot*			plot structure with FRC/FSC curve.
**/
Bplot*		img_calculate_FRC_curve(vector<Bimage*> parr, double hi_res, double sampling_ratio)
{
	if ( parr.size() < 3 ) return NULL;

	Bimage*			p = parr[0];
	Bimage*			prec = parr[1];
	Bimage*			prec2 = parr[2];

	if ( !p->data_pointer() ) return NULL;
	
	prec->check_resolution(hi_res);
	
	if ( sampling_ratio < 1 ) sampling_ratio = 1;

	long			i, j, x, y, z, iradius, iradius2, nt(0), nf(0);
	long			xx, yy, zz;
	double			rx, ry, rz, radius, fraction, fraction2, v;
	
	Vector3<double>	realsize(p->real_size());
	Vector3<double>	freq_scale(1.0/realsize);
	double			rad_scale_max(realsize[0]/sampling_ratio);
	double			rad_scale = realsize[0];
	if ( rad_scale < realsize[1] ) rad_scale = realsize[1];
	if ( rad_scale < realsize[2] ) rad_scale = realsize[2];
	
//	if ( p->sizeX() > rad_scale_max ) rad_scale *= rad_scale_max/p->sizeX();
	if ( rad_scale > rad_scale_max ) rad_scale = rad_scale_max;
	long	maxrad = (long) (rad_scale/hi_res);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_calculate_FRC_curve: maxrad=" << maxrad << endl;
	
	vector<double>	n(maxrad);
	vector<double>	sum_i_2(maxrad);
	vector<double>	sum_r_2(maxrad);
	vector<double>	sum_r2_2(maxrad);
	vector<double>	sum_i_r(maxrad);
	vector<double>	sum_i_r2(maxrad);
	
	for ( i=0; i<maxrad; i++ )
		n[i] = sum_i_2[i] = sum_r_2[i] = sum_r2_2[i] = sum_i_r[i] = sum_i_r2[i] = 0;
	
	float*			fom = (float *) prec2->next->data_pointer();
	
	for ( i=z=0; z<p->sizeZ(); z++ ) {
		zz = z; 
		if ( z > (p->sizeZ() - 1)/2 ) zz -= p->sizeZ();
		rz = zz*freq_scale[2];
		for ( y=0; y<p->sizeY(); y++ ) { 
			yy = y; 
			if ( y > (p->sizeY() - 1)/2 ) yy -= p->sizeY();
			ry = yy*freq_scale[1];
			for ( x=0; x<p->sizeX(); x++, i++ ) {
				xx = x; 
				if ( x > (p->sizeX() - 1)/2 ) xx -= p->sizeX();
				rx = xx*freq_scale[0];
				radius = rad_scale*sqrt(rx*rx + ry*ry + rz*rz); 
				iradius = (long) radius;
				iradius2 = iradius + 1;
//				if ( fom[i] && ( iradius < maxrad ) ) {
				if ( iradius2 < maxrad ) {
					fraction = radius - iradius;
					fraction2 = 1.0 - fraction;
					v = (p->complex(i)).power();
					sum_i_2[iradius] += fraction2*v;
					sum_i_2[iradius2] += fraction*v;
					v = (prec->complex(i)).power();
					sum_r_2[iradius] += fraction2*v;
					sum_r_2[iradius2] += fraction*v;
					v = (prec2->complex(i)).power();
					sum_r2_2[iradius] += fraction2*v;
					sum_r2_2[iradius2] += fraction*v;
					v = (p->complex(i)).real()*(prec->complex(i)).real() + (p->complex(i)).imag()*(prec->complex(i)).imag();
					sum_i_r[iradius] += fraction2*v;
					sum_i_r[iradius2] += fraction*v;
					v = (p->complex(i)).real()*(prec2->complex(i)).real() + (p->complex(i)).imag()*(prec2->complex(i)).imag();
					sum_i_r2[iradius] += fraction2*v;
					sum_i_r2[iradius2] += fraction*v;
					n[iradius] += fraction2;
					n[iradius2] += fraction;
					nt++;
					if ( fom[i] ) nf++;
				}
			}
		}
	}

	Bplot*		plot = new Bplot(1, maxrad, 6);
	(*plot)[0] = 0;
	(*plot)[maxrad] = 1;
	
	double	s, FRC;

	if ( verbose & VERB_PROCESS )
		cout << "Pixel\ts\tResol\tsumI2\tsumR2\tsumRo2\tsumIR\tsumIRo\tFRC\tNum" << endl;
	else if ( verbose )
		cout << "Pixel\ts\tResol\tFRC" << endl;
	for ( i=1; i<maxrad; i++ ) {
		s = i/rad_scale;
		FRC = 1.;
		if ( sum_i_r[i] < SMALLFLOAT || sum_r2_2[i] < SMALLFLOAT )
			FRC = 0;
		else
//			FRC = (sum_i_r2[i]*sqrt(sum_r_2[i]))/(sum_i_r[i]*sqrt(sum_r2_2[i]));
			FRC = (sum_i_r2[i]/sum_i_r[i])*sqrt(sum_r_2[i]/sum_r2_2[i]);
		if ( FRC > 1 ) FRC = 1;
		if ( FRC < 0 ) FRC = 0;
		if ( verbose & VERB_PROCESS )
			cout << i << tab
				<< setprecision(4) << s << tab 
				<< setprecision(2) << 1/s << tab 
				<< setprecision(0) << sum_i_2[i] << tab 
				<< sum_r_2[i] << tab << sum_r2_2[i] << tab 
				<< sum_i_r[i] << tab << sum_i_r2[i] << tab 
				<< setprecision(4) << FRC << tab 
				<< setprecision(0) << n[i] << endl;
		else if ( verbose )
			cout << i << tab
				<< setprecision(4) << s << tab 
				<< setprecision(2) << 1/s << tab 
				<< setprecision(4) << FRC << endl;
		(*plot)[i] = s;
		j = i + maxrad;
		(*plot)[j] = FRC;
		j += maxrad;
		(*plot)[j] = sum_r_2[i];
		j += maxrad;
		(*plot)[j] = sum_r2_2[i];
		j += maxrad;
		(*plot)[j] = sum_i_r[i];
		j += maxrad;
		(*plot)[j] = sum_i_r2[i];
	}

	if ( verbose & VERB_PROCESS ) {
		cout << "Total voxels:                   " << nt << endl;
		cout << "Coverage:                       " << nf << " (" << nf*100.0/nt << " %)" << endl;
	}
	
	return plot;
}

/**
@brief 	Calculates the resolution estimate from a FRC/FSC curve.  
@param 	*plot		plot structure.
@param 	cutoff		FRC/FSC cutoff to use.
@return  double			estimated resolution.
**/
double		plot_calculate_resolution(Bplot* plot, double cutoff)
{
	double			s, sp(0), FRC, FRCp(1), fraction, resolution(1000);
	long	i, maxrad = (long) plot->rows();

	for ( i=1; i<maxrad; i++ ) {
		s = (*plot)[i];
		FRC = (*plot)[i+maxrad];
		if ( FRC > cutoff ) resolution = 1/s;
		if ( FRCp > cutoff && FRC <= cutoff ) {
			fraction = (FRCp - cutoff)/(FRCp - FRC);
			if ( sp ) resolution = (1-fraction)/sp + fraction/s;
		}
		sp = s;
		FRCp = FRC;
	}

	Bstring			txt = Bstring(cutoff, "FSC(%g): ") + Bstring(resolution, "%g A");
	plot->page(0).add_text(txt);

	return resolution;
}

/**
@brief 	Plots the resolution curve in a plot structure.  
@param 	*plot			plot structure.
@param 	hi_res		high resolution limit.
@param 	&psfile		postscript output file name.
@param 	&title		plot title.
@return	int					0.
**/
int			plot_resolution(Bplot* plot, double hi_res, Bstring& psfile, Bstring& title)
{
	plot->title("Resolution");
	plot->page(0).title(title);
	plot->page(0).columns(2);
	plot->page(0).column(1).number(1);
	plot->page(0).column(0).label("Resolution(A)");
	plot->page(0).column(1).label("FRC");
	plot->page(0).axis(1).min(0);
	plot->page(0).axis(1).max(1/hi_res);
	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).inc(0.1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);

	if ( psfile.length() > 0 ) ps_plot(psfile, plot);

	return 0;
}

/**
@brief 	Plots the estimated resolution against the tilt angle.
@param 	*project		project structure.
@return	*plot			new plot.

	The resolution estimates must be encoded in the micrograph FOM's.
**/
Bplot*	plot_tilt_resolution(Bproject* project)
{
	long			nmg(0), i;
	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
	Bstring			title("Micrograph resolution for a tilt series");

	for ( nmg=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) nmg++;
	
	Bplot*			plot = new Bplot(1, nmg, 2);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(2);
	for ( i=0; i<2; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("TiltAngle");
	plot->page(0).column(1).label("NLOO");
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);
	plot->page(0).column(1).color(0,0,0);

	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, i++ ) if ( mg->select ) {
			(*plot)[i] = mg->tilt_angle*180.0/M_PI;
			(*plot)[i+nmg] = mg->fom;
		}
	}

	plot->page(0).axis(1).label("Tilt Angle (degrees)");
	plot->page(0).axis(3).label("Resolution (A)");

	return plot;
}
