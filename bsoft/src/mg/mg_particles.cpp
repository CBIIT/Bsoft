/**
@file	mg_particles.cpp
@brief	Calculates centers of single particle images.
@author Bernard Heymann
@date	Created: 20080424
@date	Modified: 20200401
**/

#include "mg_processing.h"
#include "mg_particles.h"
#include "mg_select.h"
#include "mg_img_proc.h"
#include "rwimg.h"
#include "matrix_linear.h"
#include "file_util.h"
#include "utilities.h"
#include "timer.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Aligns single particle images based on orientation parameters.
@param 	*project 	micrograph processing parameter structure.
@param 	part_select	selection number from the selection column.
@param 	nuavg		rescale to new average.
@param 	nustd		rescale to new standard deviation.
@return	int					0.

	Each selected particle image is rotated and shifted.

**/
int			project_align_particles(Bproject* project, int part_select,
				double nuavg, double nustd)
{
	long				psel, n, nc(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part ;
	Bimage*				p;
	Bimage*				p1;
	
	if ( verbose ) {
		cout << "Aligning particle images" << endl;
		cout << "Selection:                      " << part_select << endl;
		if ( nustd > 0 )
			cout << "Rescale avg & stdev:            " << nuavg << " " << nustd << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( ( p = read_img(mg->fpart, 1, -1) ) == NULL ) {
					error_show("project_align_particles", __FILE__, __LINE__);
					return -1;
				}
				for ( n=0, part = mg->part; part && n<p->images(); part = part->next, n++ ) {
					psel = 0;
					if ( part->sel ) {
						if ( part_select < 0 ) {
							psel = 1;
						} else if ( part->sel == part_select ) {
							psel = 1;
						}
					}
					if ( psel ) {	// Use only selected images
						p1 = p->extract(n);
						p1->rotate(p1->size()/2 - part->ori, part->view);
						if ( nustd > 0 ) p1->rescale_to_avg_std(nuavg, nustd);
//						p1->image->origin(p1->size()/2);
						p1->origin(p1->size()/2);
						p1->image->view(0,0,1,0);
						part->ori = p1->image->origin();
						part->view = p1->image->view();
						p->replace(n, p1);
						delete p1;
						nc++;
					}
				}
				mg->fpart = mg->fpart.pre_rev('.') + "_aln." + mg->fpart.post_rev('.');
				write_img(mg->fpart, p, 0);
				delete p;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( n=0, part = rec->part; part; part = part->next, n++ ) {
				psel = 0;
				if ( part->sel ) {
					if ( part_select < 0 ) {
						psel = 1;
					} else if ( part->sel == part_select ) {
						psel = 1;
					}
				}
				if ( psel ) {	// Use only selected images
					if ( ( p1 = read_img(part->fpart, 1, 0) ) == NULL ) {
						error_show("project_align_particles", __FILE__, __LINE__);
						return -1;
					}
					p1->rotate(p1->size()/2 - part->ori, part->view);
					if ( nustd > 0 ) p1->rescale_to_avg_std(nuavg, nustd);
//					p1->image->origin(p1->size()/2);
					p1->origin(p1->size()/2);
					p1->image->view(0,0,1,0);
					part->ori = p1->image->origin();
					part->view = p1->image->view();
					part->fpart = part->fpart.pre_rev('.') + "_aln." + part->fpart.post_rev('.');
					write_img(part->fpart, p1, 0);
					delete p1;
					nc++;
				}
			}
		}
	}
	
	if ( verbose )
		cout << "Particles oriented:             " << nc << endl << endl;
	
	return 0;
}

/**
@brief 	Centers single particle images based on the parametric center.  
@param 	*project 	micrograph processing parameter structure.
@param 	part_select	selection number from the selection column.
@param 	nuavg		rescale to new average.
@param 	nustd		rescale to new standard deviation.
@return	int					0.

	Each particle image is shifted to center the origin.

**/
int			project_center_particles(Bproject* project, int part_select,
				double nuavg, double nustd)
{
	long				psel, n, np, nc(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Bimage*				p = NULL;
	Bimage*				pm = NULL;
	
	if ( verbose ) {
		cout << "Centering particle images" << endl;
		cout << "Selection:                      ";
		if ( part_select < 0 ) cout << "all" << endl;
		else cout  << part_select << endl;
		if ( nustd > 0 )
			cout << "Rescale avg & stdev:            " << nuavg << " " << nustd << endl;
	}

	if ( project->select < 1 ) {
		for ( field=project->field; field; field=field->next ) {
			for ( mg=field->mg; mg; mg=mg->next ) {
				for ( np=0, part=mg->part; part; part=part->next, np++ ) ;
				for ( n=0, part=mg->part; part; part=part->next, n++ ) {
					if ( ( p = particle_read_img(part, 1) ) == NULL ) {
						error_show("project_center_particles", __FILE__, __LINE__);
						return -1;
					}
					p->calculate_background();
					psel = 0;
					if ( part->sel ) {
						if ( part_select < 0 ) {
							psel = 1;
						} else if ( part->sel == part_select ) {
							psel = 1;
						}
					}
					if ( psel ) {	// Use only selected images
						p->origin(part->ori);
						part->loc -= p->size()/2 - part->ori;
						part->ori = p->size()/2;
					} else {
						p->origin(p->size()/2);
					}
					p->center(FILL_BACKGROUND);
					if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
					if ( verbose & VERB_FULL )
						cout << part->id << tab << p->image->origin() << tab 
							<< p->image->background() << endl;
					nc++;
					if ( part->fpart.length() ) {
						part->fpart = part->fpart.pre_rev('.') + "_cen." + part->fpart.post_rev('.');
						write_img(part->fpart, p, 0);
					} else {
						if ( n == 0 ) {
							pm = new Bimage(p->data_type(), p->compound_type(), p->size(), np);
							pm->origin(pm->default_origin());
							pm->sampling(p->sampling(0));
						}
						pm->replace(n, p);
					}
					delete p;
				}
				if ( pm && mg->fpart.length() ) {
					mg->fpart = mg->fpart.pre_rev('.') + "_cen." + mg->fpart.post_rev('.');
					write_img(mg->fpart, pm, 0);
					delete pm;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( n=0, part = rec->part; part && n<p->images(); part = part->next, n++ ) {
				psel = 0;
				if ( part->sel ) {
					if ( part_select < 0 ) {
						psel = 1;
					} else if ( part->sel == part_select ) {
						psel = 1;
					}
				}
				if ( psel ) {	// Use only selected images
					if ( ( p = read_img(part->fpart, 1, 0) ) == NULL ) {
						error_show("project_align_particles", __FILE__, __LINE__);
						return -1;
					}
					p->image[n].origin(part->ori);
					part->loc -= p->size()/2 - part->ori;
					part->ori = p->size()/2;
					p->calculate_background();
					if ( verbose & VERB_FULL )
						cout << part->id << tab << p->image->origin() << tab 
							<< p->background(long(0)) << endl;
					p->center(FILL_BACKGROUND);
					if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
					rec->fpart = rec->fpart.pre_rev('.') + "_cen." + rec->fpart.post_rev('.');
					write_img(rec->fpart, p, 0);
					delete p;
					nc++;
				}
			}
		}
	}

	if ( verbose )
		cout << "Particles oriented:             " << nc << endl << endl;
	
	return 0;
}

/**
@brief 	Calculates centers of single particle images by cross-correlation with a reference.  
@param 	*project 	micrograph processing parameter structure.
@param 	*pref		reference image.
@param 	part_select	selection number from the selection column.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return	int					0.

	Each particle image shift is determined by cross correlation with
	the reference image.
	The reference must be the same size as the particle images.

**/
int			project_set_particle_centers(Bproject* project, Bimage* pref, 
				int part_select, double hires, double lores)
{
	int				psel;
//	Vector3<float>*	shift;
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;
	Bimage*			p;
	
	if ( verbose )
		cout << "Setting the centers of particle images" << endl << endl;

	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			for ( part=mg->part; part; part=part->next ) {
				psel = 0;
				if ( part->sel ) {
					if ( part_select < 0 ) {
						psel = 1;
					} else if ( part->sel == part_select ) {
						psel = 1;
					}
				}
				if ( psel ) {	// Use only selected images
					if ( ( p = read_img(mg->fpart, 1, part->id - 1) ) == NULL ) {
						error_show("mg_set_particle_centers", __FILE__, __LINE__);
						return -1;
					}
					p->origin(p->size()/2);
//					img_find_shift(pref, p, NULL, hires, lores, p->sizeX()/4, 0, 1);
					pref->find_shift(p, NULL, hires, lores, p->sizeX()/4, 0, 1);
					part->ori = p->image->origin();
					delete p;
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Finds the centers of picked particles within a micrograph.
@param 	*project		project parameter structure.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	filter_flag		flag to filter extremes in particles.
@return int				0.

	An image processing parameter structure loaded with micrograph
	information is used to extract particle images from the micrograph
	image using the particle coordinates in the parameter structure.
	The extracted particle images are each rotated by PI and the shift
	found by cross-correlation between the unrotated and rotated images.
	The particle coordinates in the parameter structure are updated with
	the shift.

**/
long		project_find_part_centers_in_mgs(Bproject* project, 
				double hires, double lores, int filter_flag)
{
	long			npart(0);
	
	if ( verbose & VERB_LABEL )
		cout << "Finding particle centers with resolution limits " << hires << " - " << lores << " A" << endl;
	
	Bfield*			field;
	Bmicrograph*	mg;
	
	for ( field=project->field; field; field=field->next )
		for ( mg=field->mg; mg; mg=mg->next )
			npart += mg_find_part_centers(mg, hires, lores, filter_flag);
	
	return npart;
}

/**
@brief 	Finds the centers of picked particles within a micrograph.
@param 	*mg				micrograph parameter structure.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	filter_flag		flag to filter extremes in particles.
@return long			number of particles.

	An image processing parameter structure loaded with micrograph
	information is used to extract particle images from the micrograph
	image using the particle coordinates in the parameter structure.
	The extracted particle images are each rotated by PI and the shift
	found by cross-correlation between the unrotated and rotated images.
	The particle coordinates in the parameter structure are updated with
	the shift.

**/
long		mg_find_part_centers(Bmicrograph* mg, 
				double hires, double lores, int filter_flag)
{
	Vector3<long> 		box_size(mg->box_size);
	
	if ( verbose & VERB_LABEL )
		cout << "Finding particle centers with resolution limits " << hires << " - " << lores << " A" << endl;
	
	Vector3<double>	start;
	Bparticle*		part;	
	Bimage*			pex = NULL;
	Bimage*			pexrot = NULL;
	
	long			npart(particle_count(mg->part));
	if ( npart < 1 ) return 0;
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);

	if ( filter_flag ) p->filter_extremes(1);
	
	if ( verbose & VERB_PROCESS )
		cout << "Finding particle origins from image " << p->file_name() << endl;
	
	for ( part=mg->part; part; part=part->next ) {
		start = part->loc - box_size/2;
		pex = p->extract(0, start, box_size);
		pexrot = pex->rotate(pex->size(), M_PI);
		pex->find_shift(pexrot, NULL, hires, lores, pex->sizeX()/8.0, 0, 1);
		part->loc += (pexrot->image->origin() - pex->image->origin())/2;
		part->fom[0] = pexrot->image->FOM();
		delete pex;
		delete pexrot;
	}
	
	delete p;
	
	return npart;
}

/**
@brief 	Calculates centers of single particle images.  
@param 	*project 	image processing parameter structure.
@param 	max_iter	maximum number of iterations.
@param 	part_select	selection number from the selection column.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@return	Bimage*		final image composite reference.

	A composite image is generated from the selected particles and radially
	symmetrized. Each image shift is then determined by cross correlation.
	This is repeated until the chnage in image shifts decreases below
	a threshold.

**/
Bimage*		project_find_particle_centers(Bproject* project, int max_iter, 
				int part_select, double hires, double lores)
{
	int				iter, psel, n;
	double			d, R(1e30), change(1e30);
	Vector3<float>	origin, shift1;
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bparticle*		part = mg->part;
	View			ref_view;
	Bimage*			p;
	Bimage* 		pref = NULL;
	
	if ( verbose )
		cout << "Finding the centers of particle images" << endl << endl;
	
	if ( verbose & VERB_RESULT )
		cout << "Iter\tR" << endl;

	for ( iter=0; iter<max_iter && R > 0.1 && change > 0.001; iter++ ) {
		if ( pref ) pref->clear();
		for ( field=project->field; field; field=field->next ) {	// Generating the reference image
			for ( mg=field->mg; mg; mg=mg->next ) {
				for ( part=mg->part; part; part=part->next ) {
					psel = 0;
					if ( part->sel ) {
						if ( part_select < 0 ) {
							psel = 1;
						} else if ( part->sel == part_select ) {
							psel = 1;
						}
					}
					if ( psel ) {	// Use only selected images
						p = particle_read_img(part, 1);
						p->change_type(Float);
						if ( part->ori.length() < 1 ) part->ori = p->size()/2;
						if ( !pref ) {
							pref = p->copy();
//							pref->image->origin(pref->size()/2);
							pref->origin(pref->size()/2);
							pref->clear();
						}
						shift1 = pref->image->origin() - part->ori;
						if ( verbose & VERB_FULL )
							cout << shift1 << endl;
						p->shift_wrap(shift1);
						pref->add(p);
						delete p;
					}
				}
			}
		}
		pref->calculate_background();
		pref->symmetrize_cylinder();
		pref->calculate_background();
		pref->shift_background(0);
		write_img("part_avg.pif", pref, 0);
		n = 0;
		change = R;
		R = 0;
		for ( field=project->field; field; field=field->next ) {	// Determining the shift for centering
			for ( mg=field->mg; mg; mg=mg->next ) {
				for ( part=mg->part; part; part=part->next ) {
					psel = 0;
					if ( part->sel ) {
						if ( part_select < 0 ) {
							psel = 1;
						} else if ( part->sel == part_select ) {
							psel = 1;
						}
					}
					if ( psel ) {	// Use only selected images
						p = particle_read_img(part, 1);
						p->change_type(Float);
						p->origin(part->ori);
						p->find_shift(pref, NULL, hires, lores, p->sizeX()/8, 0, 1);
//						origin = - p->image->origin() + p->size();
						origin = p->image->origin();
						origin[2] = 0;
						if ( verbose & VERB_FULL )
							cout << origin << endl;
						d = (origin - part->ori).length2();
						R += d;
						part->ori = origin;
						delete p;
						n++;
					}
				}
			}
		}
		R = sqrt(R/n);
		change -= R;
		if ( verbose & VERB_RESULT )
			cout << iter+1 << tab << R << endl;
	}
	
	if ( verbose & VERB_RESULT )
		cout << endl;
	
	return pref;
}

long		particles_mask(Bmicrograph* mg, Bimage* pmask, Bstring& partpath)
{
	long			n(0);
	Matrix3			mat;
	Vector3<double>	translate;
	Bparticle*		part = mg->part;
	Bimage*			p;
	Bimage* 		ppart = NULL;
	Bimage*			pm;
	
	if ( mg->part && mg->fpart.length() && mg->part->fpart.length() < 1 ) {
		ppart = read_img(mg->fpart, 1, -1);	
		ppart->sampling(part->pixel_size);
		ppart->change_type(Float);
	}
	
	for ( part = mg->part; part; part = part->next ) {
		if ( part->fpart.length() ) {
			p = read_img(part->fpart, 1, 0);
		} else {
			p = ppart->extract(part->id-1);
		}
		p->change_type(Float);
		p->rescale_to_avg_std(0,1);
		p->image->view(part->view);
//		p->image->origin(part->ori);
		p->origin(part->ori);
		translate = part->ori - p->size()/2;
		mat = p->image->view().matrix();
		pm = pmask->rotate_project(mat, translate, p->sizeX()/2.0);
		pm->truncate(0, 0.001*pm->standard_deviation(), 0, 1);
		pm->filter_average(7);
		p->multiply(pm);
		delete pm;
		if ( part->fpart.length() ) {
			part->fpart = part->fpart.pre_rev('.') + "_m." + part->fpart.post_rev('.');
			if ( partpath.length() > 1 ) 
				part->fpart = partpath + part->fpart.post_rev('/');
			write_img(part->fpart, p, 0);
		} else {
			ppart->replace(part->id - 1, p);
		}
		delete p;
		n++;
	}
	
	if ( ppart ) {
		mg->fpart = mg->fpart.pre_rev('.') + "_m." + mg->fpart.post_rev('.');
		if ( partpath.length() > 1 ) 
			mg->fpart = partpath + mg->fpart.post_rev('/');
		write_img(mg->fpart, ppart, 0);
		delete ppart;
	}
	
	return n;
}

/**
@brief 	Calculates centers of single particle images.  
@param 	*project 	image processing parameter structure.
@param 	*pmask	 	3D volume mask to be projected.
@param 	&partpath	new path to particle files.
@return	long		number of particles masked.

	A 3D mask is projected into each particle view and the particle image masked.

**/
long		project_mask_particles(Bproject* project, Bimage* pmask, Bstring& partpath)
{
	pmask->change_type(Float);
	
	if ( !partpath.empty() ) if ( partpath[-1] != '/' ) partpath += "/";
	
	long			nmg(0);
	Matrix3			mat;
	Vector3<double>	translate;
	
//	Bfield*			field = project->field;

	if ( partpath.length() > 1 )
		mkdir(partpath.c_str(), (mode_t)0755);
	
	if ( verbose )
		cout << "Masking particle images" << endl << endl;

	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);

#ifdef HAVE_GCD
	__block	long	n(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nmg, dispatch_get_global_queue(0, 0), ^(size_t i){
		long	n1 = particles_mask(mgarr[i], pmask, partpath);
		dispatch_sync(myq, ^{
			n += n1;
			if ( verbose & VERB_RESULT )
				cout << setw(15) << mgarr[i]->id << tab << n << endl;
		});
	});
#else
	long			n(0);
#pragma omp parallel for
	for ( long i=0; i<nmg; i++ ) {
		long	n1 = particles_mask(mgarr[i], pmask, partpath);
	#pragma omp critical
		{
			n += n1;
			if ( verbose & VERB_RESULT )
				cout << setw(15) << mgarr[i]->id << tab << n << endl;
		}
	}
#endif
	
	delete[] mgarr;
	
/*	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose ) cout << "Micrograph: " << mg->id << endl;
			particles_mask(mg, pmask, partpath);
		}
	}
*/	
	if ( verbose & VERB_RESULT )
		cout << endl;
	
	return n;
}

/**
@brief 	Compares the coordinates of particles between two files.  
@param 	*project 	project parameter structure.
@param 	*projcomp 	comparable project parameter structure.
@return	long		number of common particles.

	The coordinates of particles in one parameter file is compared to
	that of a reference parameter file.
	The two parameter files must have the same field and micrograph ID's.

**/
long		project_compare_particles(Bproject* project, Bproject* projcomp)
{
	long			nf, tp(0), tn(0), fp(0);
	double			d, r;
	
	Bfield*			field, *fieldcomp;
	Bmicrograph*	mg, *mgcomp;
	Bparticle*		part, *partcomp;
	
	for ( field = project->field; field; field = field->next ) {
		for ( fieldcomp = projcomp->field; fieldcomp && fieldcomp->id != field->id; fieldcomp = fieldcomp->next ) ;
		if ( fieldcomp ) {
			for ( mg = field->mg; mg ; mg = mg->next ) {
				for ( mgcomp = fieldcomp->mg; mgcomp && mgcomp->id != mg->id; mgcomp = mgcomp->next ) ;
				if ( mgcomp ) {
					r = mg->box_size[0]/2;
					for ( partcomp = mgcomp->part; partcomp; partcomp = partcomp->next ) {
						partcomp->sel = 0;
						partcomp->fom[0] = 2*r;
					}
					for ( part = mg->part; part; part = part->next ) {
						part->sel = 0;
						nf = 0;
						for ( partcomp = mgcomp->part; partcomp; partcomp = partcomp->next ) {
							d = part->loc.distance(partcomp->loc);
							if ( d < r && partcomp->fom[0] > d ) {
								partcomp->fom[0] = d;
								partcomp->sel = part->id;
								part->sel = partcomp->id;
								nf++;
							}
						}
						if ( nf ) tp++;
						else tn++;
					}
					for ( partcomp = mgcomp->part; partcomp; partcomp = partcomp->next )
						if ( partcomp->sel < 1 ) fp++;
				}
			}
		}
	}
	
	if ( verbose ) {
		cout << "Comparing " << project->filename << " against " << projcomp->filename << endl;
		cout << "True positive:                  " << tp << endl;
		cout << "True negative:                  " << tn << endl;
		cout << "False positive:                 " << fp << endl << endl;
	}
	
	return tp;
}

/**
@brief 	Finds the tilt axis from refined particle defocus values.  
@param 	*project 	project parameter structure.
@return	double		average tilt axis.


**/
double		project_tilt_from_particle_defocus(Bproject* project)
{
	long			i, j, nmg(0), nmgsel(0), npart(0), nptot(0), npsel(0);
	double			axis, axis_avg(0), axis_std(0), a1(0), da;
	double			tilt, tilt_avg(0), tilt_std(0);
	double			ddef, ddef_rmsd(0);
	Vector3<double>	loc, normal;
	vector<double>	b(3);
	Matrix			a(3,3);
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	nptot = project_count_mg_particles(project);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg ; mg = mg->next, nmg++ ) if ( mg->select ) {
			for ( i=0; i<3; i++ ) {
				for ( j=0; j<3; j++ )
					a[j][i] = 0;
				b[i] = 0;
			}
			npart = 0;
			ddef_rmsd = 0;
			for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
				loc = part->loc - mg->origin;
				ddef = part->def - mg->ctf->defocus_average();
				ddef_rmsd += ddef*ddef;
				loc[2] = ddef/mg->pixel_size[0];
				for ( i=0; i<3; i++ ) {
					for ( j=0; j<3; j++ )
						a[j][i] += loc[i]*loc[j];
					b[i] += loc[i];
				}
				npart++;
			}
			npsel += npart;
			normal = fit_plane(a, b);
			if ( normal[2] < 0 ) normal = -normal;
			axis = atan2(-normal[0], normal[1]);
			tilt = acos(normal[2]);
			if ( nmg > 0 ) {
				da = fabs(angle_set_negPI_to_PI(axis - a1));
				if ( da > M_PI_2 ) {
					axis = angle_set_negPI_to_PI(axis + M_PI);
					tilt = -tilt;
				}
			} else {
				a1 = axis;
			}
			if ( npart )
				ddef_rmsd = sqrt(ddef_rmsd/npart);
			if ( npart > 5 ) {
				axis_avg += axis;
				axis_std += axis*axis;
				tilt_avg += tilt;
				tilt_std += tilt*tilt;
				nmgsel++;
			}
			if ( verbose )
				cout << mg->id << tab << npart << tab << ddef_rmsd << tab << axis*180.0/M_PI << tab << tilt*180.0/M_PI << endl;
		}
	}
	
	if ( nmgsel ) {
		axis_avg /= nmgsel;
		axis_std = axis_std/nmgsel - axis_avg*axis_avg;
		if ( axis_std > 0 ) axis_std = sqrt(axis_std);
		tilt_avg /= nmgsel;
		tilt_std = tilt_std/nmgsel - tilt_avg*tilt_avg;
		if ( tilt_std > 0 ) tilt_std = sqrt(tilt_std);
	}
	
	if ( verbose ) {
		cout << "Number of micrographs:    " << nmgsel << " (" << nmgsel*100.0/nmg << "%)" << endl;
		cout << "Number of particles:      " << npsel << " (" << npsel*100.0/nptot << "%)" << endl;
		cout << "Axis average & stdev:     " << axis_avg*180.0/M_PI << tab << axis_std*180.0/M_PI << endl;
		cout << "Tilt average & stdev:     " << tilt_avg*180.0/M_PI << tab << tilt_std*180.0/M_PI << endl;
	}
	
	return axis_avg;
}

/**
@brief 	Finds the tilt axis from refined particle defocus values.  
@param 	*project 	project parameter structure.
@param 	axis	 	tilt axis angle (radians).
@param 	tilt	 	tilt angle (radians).
@return	long		number of selected particles.


**/
long		project_set_particle_defocus_from_tilt(Bproject* project, double axis, double tilt)
{
	long			npart(0);
	double			ca(cos(axis)), sa(sin(axis)), tt(tan(tilt));
	Vector3<double>	coor;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg ; mg = mg->next ) if ( mg->select ) {
			for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
				coor = (part->loc - mg->origin)*mg->pixel_size;			
				part->def = mg->ctf->defocus_average() + (coor[0]*sa - coor[1]*ca)*tt;
				if ( verbose )
					cout << part->loc << tab << part->def - mg->ctf->defocus_average() << endl;
				npart++;
			}
		}
	}
	
	return npart;
}

