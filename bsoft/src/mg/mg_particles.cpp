/**
@file	mg_particles.cpp
@brief	Analyzes and manipulates single particle images.
@author Bernard Heymann
@date	Created: 20080424
@date	Modified: 20220531
**/

#include "mg_processing.h"
#include "mg_particles.h"
#include "mg_particle_select.h"
#include "mg_select.h"
#include "mg_img_proc.h"
#include "mg_ctf.h"
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
@param 	*project 		micrograph processing parameter structure.
@param 	part_select		selection number from the selection column.
@param 	nuavg			rescale to new average.
@param 	nustd			rescale to new standard deviation.
@return	int				0.

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
@param 	*project 		micrograph processing parameter structure.
@param 	part_select		selection number from the selection column.
@param 	nuavg			rescale to new average.
@param 	nustd			rescale to new standard deviation.
@return	int				0.

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
@param 	*project 		micrograph processing parameter structure.
@param 	*pref			reference image.
@param 	part_select		selection number from the selection column.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@return	int				0.

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
@return long				number of particles.

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
@param 	*project 		image processing parameter structure.
@param 	max_iter		maximum number of iterations.
@param 	part_select		selection number from the selection column.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@return	Bimage*			final image composite reference.

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
@return	long			number of particles masked.

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
@return	long			number of common particles.

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
@return	long			number of selected particles.


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

Bimage*		particle_correlation_sum(Bparticle* part, Bimage* pref, double hires, FSI_Kernel* kernel, fft_plan planf)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( pref->fourier_type() != Standard ) {
		pref->fft();
		pref->phase_shift_to_origin();
	}

	bool			invert(0);
	long			ndone(0);
	double			wl(0);
	Matrix3			mat;
	Bmicrograph*	mg = part->mg;
	Bimage*			psec;
	Bimage*			psec2;
	Bimage*			ppart;

	//			if ( mg->ctf ) wl = mg->ctf->lambda();

	Bimage*			psum = new Bimage(Float, TComplex, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(part->pixel_size);
	psum->origin(psum->size()/2);

	psum->next = new Bimage(Float, TSimple, pref->sizeX(), pref->sizeY(), 1, 1);

	for ( ; part; part = part->next ) if ( part->sel > 0 ) {
		ppart = read_img(mg->fpart, 1, part->id-1);
		ppart->sampling(part->pixel_size);
		ppart->origin(part->ori);
		ppart->view(part->view);
		ppart->fft(planf);
		ppart->phase_shift_to_origin();
		mat = part->view.matrix();
		part->ori[2] = 0;
//		cout << mat << part->ori << endl;
		psec = pref->central_section(mat, hires, kernel, wl);
		if ( mg->ctf ) {
			if ( part->def > 0 ) mg->ctf->defocus_average(part->def);
			img_ctf_apply(psec, *mg->ctf, 2, 0.1, 0, hires, invert);
		}
		psec2 = psec->copy();
		psec2->complex_to_intensities();
		psec->complex_conjugate_product(ppart, 0);
		psum->add(psec);
		psum->next->add(psec2);
		delete ppart;
		delete psec;
		delete psec2;
		ndone++;
	}
	
	psum->image->select(ndone);

	return psum;
}

/**
@brief 	Correlates each particle with the CTF applied reference projection and sum the correlations.
@param 	*project 	project parameter structure.
@param	*pref		reference map.
@param 	hires	 	high resolution limit (angstrom).
@param 	*kernel	 	frequency space interpolation kernel lookup table.
@return	Bimage*		sum of cross-correlations.


**/
Bimage*		project_correlation_sum(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	long			nsel = project_count_mg_part_selected(project);
	
	if ( verbose ) {
		cout << "Cross-correlating particles with CTF-applied reference projections:" << endl;
		cout << "Reference map:                  " << pref->file_name() << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Selected particles:             " << nsel << endl;
	}
	
	if ( pref->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the reference map" << endl;
		pref->fft();
	}

	pref->phase_shift_to_origin();

	long			nmg(0), nogrp(0);
	Matrix3			mat;
	Bparticle*		part = part_find_first(project);
	
	if ( verbose )
		cout << "Setting up the micrograph array" << endl;
	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);
	
	map<string,long>	optgrp;
	for ( long i=0; i<nmg; ++i )
		if ( mgarr[i]->ctf )
			if ( optgrp.find(mgarr[i]->ctf->identifier()) == optgrp.end() )
				optgrp[mgarr[i]->ctf->identifier()] = 0;

	for ( auto& g: optgrp ) g.second = nogrp++;	// Ensure the order is encoded properly
	
	if ( nogrp < 1 ) nogrp = 1;
	
	Bimage*			psum = new Bimage(Float, TComplex, pref->sizeX(), pref->sizeY(), 1, nogrp);
	psum->sampling(part->pixel_size);
	psum->origin(psum->size()/2);
	
	psum->next = new Bimage(Float, TSimple, pref->sizeX(), pref->sizeY(), 1, nogrp);
	psum->next->fill(1);

	fft_plan		planf = fft_setup_plan(psum->size(), FFTW_FORWARD, 1);
	
#ifdef HAVE_GCD
	__block	long	ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nmg, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp.at(mgarr[i]->ctf->identifier());
		Bimage*		p1 = particle_correlation_sum(mgarr[i]->part, pref, hires, kernel, planf);
		dispatch_sync(myq, ^{
			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->next->add(n, p1->next);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		});
	});
#else
	long			ndone(0);
#pragma omp parallel for
	for ( long i=0; i<nmg; ++i ) {
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp[mgarr[i]->ctf->identifier()];
		Bimage*		p1 = particle_correlation_sum(mgarr[i]->part, pref, hires, kernel, planf);
	#pragma omp critical
		{
			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->next->add(n, p1->next);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		}
	}
#endif
	
	fft_destroy_plan(planf);

	if ( verbose & VERB_RESULT )
		cerr << endl;

	if ( verbose ) {
		cout << "Group\tParticles" << endl;
//		for ( long n=0; n<optgrp.size(); ++n )
//			cout << optgrp[n].first << tab << psum->image[n].select() << endl;
		for ( auto g: optgrp )
			cout << g.first << tab << psum->image[g.second].select() << endl;
		cout << endl;
	}

	psum->divide(psum->next);

	delete[] mgarr;
	delete psum->next;
	psum->next = NULL;

	return psum;
}

int			img_add_aberration_terms(Bimage* psum, Bimage* psec, Bimage* ppart, CTFparam& cp)
{
	long 			i, j, n, x, y, z;
	double			sx, sy, sz, s2, phi, pwr;
	Complex<double>	cv, ca;
	Vector3<double>	freq_scale(1.0/psum->real_size());
	Vector3<double>	h((psum->sizeX() - 1)/2, (psum->sizeY() - 1)/2, (psum->sizeZ() - 1)/2);

	if ( verbose & VERB_DEBUG ) {
		cout << "Calculating aberration terms" << endl;
		cout << endl;
	}
	
	for ( i=j=n=0; n<psum->images(); n++ ) {
		for ( z=0; z<psum->sizeZ(); z++ ) {
			sz = z;
			if ( z > h[2] ) sz -= psum->sizeZ();
			sz *= freq_scale[2];
			for ( y=0; y<psum->sizeY(); y++ ) {
				sy = y;
				if ( y > h[1] ) sy -= psum->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<psum->sizeX(); x++, i++ ) {
					sx = x;
					if ( x > h[0] ) sx -= psum->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					phi = atan2(sy,sx);
					ca = cp.calculate_complex(s2, phi);
					cv = psec->complex(i);
					pwr = cv.power();
					cv *= ppart->complex(i).conj();
					psum->add(j++, cv.real() * ca.imag());			// (S*P')re * sin g
					psum->add(j++, cv.imag() * ca.imag());			// (S*P')im * sin g
					psum->add(j++, cv.real() * ca.real());			// (S*P')re * cos g
					psum->add(j++, pwr * ca.real() * ca.real());	// S2 * cos2 g
					psum->add(j++, pwr * ca.real() * ca.imag());	// S2 * cos g * sin g
					psum->add(j++, pwr * ca.imag() * ca.imag());	// S2 * sin2 g
				}
			}
		}
	}

	return 0;
}


Bimage*		particle_aberration_sum(Bparticle* part, Bimage* pref, double hires, FSI_Kernel* kernel, fft_plan planf)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( pref->fourier_type() != Standard ) {
		pref->fft();
		pref->phase_shift_to_origin();
	}

	long			ndone(0);
	double			wl(0);
	Matrix3			mat;
	Bmicrograph*	mg = part->mg;
	Bimage*			psec;
	Bimage*			ppart;
	
	if ( !mg->ctf ) {
		cerr << "Error: The CTF parameters must be defined for micrograph " << mg->id << endl;
		bexit(-1);
	}

	wl = mg->ctf->lambda();

	Bimage*			psum = new Bimage(Float, 6, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(part->pixel_size);

	for ( ; part; part = part->next ) if ( part->sel > 0 ) {
		if ( part->def > 0 ) mg->ctf->defocus_average(part->def);
		ppart = read_img(mg->fpart, 1, part->id-1);
		ppart->sampling(part->pixel_size);
		ppart->origin(part->ori);
		ppart->view(part->view);
		ppart->fft(planf);
		ppart->phase_shift_to_origin();
		mat = part->view.matrix();
		part->ori[2] = 0;
//		cout << mat << part->ori << endl;
		psec = pref->central_section(mat, hires, kernel, wl);
		img_add_aberration_terms(psum, psec, ppart, *mg->ctf);
		delete ppart;
		delete psec;
		ndone++;
	}
	
	psum->image->select(ndone);

	return psum;
}

int			img_check_contrast(Bimage* psum)
{
	int				contrast(0);
	long			j, xx, yy, zz;
	
	for ( zz=0; zz<3 && zz<psum->sizeZ(); ++zz ) {
		for ( yy=0; yy<3 && yy<psum->sizeY(); ++yy ) {
			for ( xx=0; xx<3 && xx<psum->sizeX(); ++xx ) {
				j = psum->index(0, xx, yy, zz, 0);
				if ( fabs(atan2((*psum)[j+2],(*psum)[j])) < M_PI_2 ) contrast += 1;
				else contrast -= 1;
			}
		}
	}
	
	if ( contrast < 0 ) contrast = -1;
	else contrast = 1;
	
	if ( verbose )
		cout << "Contrast check:                " << contrast << endl;

	return contrast;
}

Bimage*		img_calculate_phase_differences(Bimage* psum)
{
	long			i, j, k, n;
	double			ae, ao;
	vector<double>	t(2);
	Matrix			r(2,2);

	int				contrast = img_check_contrast(psum);
	
	if ( verbose )
		cout << "Calculating phase difference maps" << endl;

	Bimage*			pphi = new Bimage(Float, TComplex, psum->size(), psum->images());
	pphi->sampling(psum->sampling(0));

	for ( n=i=j=0; n<pphi->images(); ++n ) {
		for ( k=0; k<pphi->image_size(); ++i, ++k, j+=psum->channels() ) {
			t[0] = (*psum)[j+2];
			t[1] = (*psum)[j];
			r[0][0] = (*psum)[j+3];
			r[0][1] = r[1][0] = (*psum)[j+4];
			r[1][1] = (*psum)[j+5];
			r.singular_value_decomposition(t);
			ae = atan2(contrast*t[0],contrast*t[1]);
			ao = atan2(contrast*(*psum)[j+1],contrast*(*psum)[j]);
			pphi->set(i, Complex<double>(ae, ao));
		}
	}
	
	pphi->statistics();

	return pphi;
}

/**
@brief 	Calculates the aberration phase difference from a set of particles.
@param 	*project 	project parameter structure.
@param	*pref		reference map.
@param 	hires	 	high resolution limit (angstrom).
@param 	*kernel	 	frequency space interpolation kernel lookup table.
@return	Bimage*		aberration phase differences as a complex image.

	Uses the method as described in Zivanov et al (2020).
	The return image is complex with the even and odd phase differences
	in real and imaginary parts, respectively.

**/
Bimage*		project_aberration_phase_difference(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	long			nsel = project_count_mg_part_selected(project);
	
	if ( verbose ) {
		cout << "Calculating phase difference between particles and CTF-applied reference projections:" << endl;
		cout << "Reference map:                  " << pref->file_name() << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Selected particles:             " << nsel << endl;
	}
	
	if ( pref->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the reference map" << endl;
		pref->fft();
	}

	pref->phase_shift_to_origin();

	long			nmg(0), nogrp(0);
	Matrix3			mat;
	Bparticle*		part = part_find_first(project);
	
	part_select_micrographs_with_selected_particles(project);
	
	if ( verbose )
		cout << "Setting up the micrograph array" << endl;
	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);
	
	map<string,long>	optgrp;
	for ( long i=0; i<nmg; ++i )
		if ( mgarr[i]->ctf )
			if ( optgrp.find(mgarr[i]->ctf->identifier()) == optgrp.end() )
				optgrp[mgarr[i]->ctf->identifier()] = 0;

	for ( auto& g: optgrp ) g.second = nogrp++;	// Ensure the order is encoded properly
	
	if ( nogrp < 1 ) nogrp = 1;
	
	Bimage*			psum = new Bimage(Float, 6, pref->sizeX(), pref->sizeY(), 1, nogrp);
	psum->sampling(part->pixel_size);

	fft_plan		planf = fft_setup_plan(psum->size(), FFTW_FORWARD, 1);
	
#ifdef HAVE_GCD
	__block	long	ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nmg, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp.at(mgarr[i]->ctf->identifier());
		Bimage*		p1 = particle_aberration_sum(mgarr[i]->part, pref, hires, kernel, planf);
		dispatch_sync(myq, ^{
//			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		});
	});
#else
	long			ndone(0);
#pragma omp parallel for
	for ( long i=0; i<nmg; ++i ) {
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp[mgarr[i]->ctf->identifier()];
		Bimage*		p1 = particle_aberration_sum(mgarr[i]->part, pref, hires, kernel, planf);
	#pragma omp critical
		{
//			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		}
	}
#endif
	
	fft_destroy_plan(planf);

	if ( verbose & VERB_RESULT )
		cerr << endl;

	if ( verbose ) {
		cout << "#\tGroup\tParticles" << endl;
		for ( auto g: optgrp )
			cout << g.second << tab << g.first << tab << psum->image[g.second].select() << endl;
		cout << endl;
	}

	Bimage*		pphi = img_calculate_phase_differences(psum);

	delete[] mgarr;
	delete psum;

	return pphi;
}

int			img_add_ewald_terms(Bimage* psum, Bimage* psec, Bimage* ppart, CTFparam& cp)
{
	long 			i, j, n, x, y, z, ix, iy, iz;
	double			sx, sy, sz, s2, phi;
	Complex<double>	cv, cup, clo, ca, csum, cdif;
	Vector3<double>	freq_scale(1.0/psum->real_size());
	Vector3<double>	h((psum->sizeX() - 1)/2, (psum->sizeY() - 1)/2, (psum->sizeZ() - 1)/2);

	if ( verbose & VERB_DEBUG ) {
		cout << "Calculating aberration terms" << endl;
		cout << endl;
	}
	
	for ( i=j=n=0; n<psum->images(); n++ ) {
		for ( z=iz=0; z<psum->sizeZ(); z++ ) {
			if ( z ) iz = psum->sizeZ() - z;
			sz = z;
			if ( z > h[2] ) sz -= psum->sizeZ();
			sz *= freq_scale[2];
			for ( y=iy=0; y<psum->sizeY(); y++ ) {
				if ( y ) iy = psum->sizeY() - y;
				sy = y;
				if ( y > h[1] ) sy -= psum->sizeY();
				sy *= freq_scale[1];
				for ( x=ix=0; x<psum->sizeX(); x++, i++ ) {
					if ( x ) ix = psum->sizeX() - x;
					sx = x;
					if ( x > h[0] ) sx -= psum->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					phi = atan2(sy,sx);
					cv = ppart->complex(i).conj();
					ca = cp.calculate_complex(s2, phi);
					cup = psec->complex(i);
					clo = psec->complex(psec->index(ix,iy,iz,n)).conj();
					csum = (cup + clo)*ca.imag();
					cdif = (cup - clo)*ca.real();
					psum->add(j++, csum.real() * cv.real() + csum.imag() * cv.imag());	// (Fs*F')re * sin g
					psum->add(j++, cdif.real() * cv.real() + cdif.imag() * cv.imag());	// (Fd*F')re * cos g
					psum->add(j++, csum.power());										// |Fs * sin g|2
					psum->add(j++, csum.real() * cdif.real() + csum.imag() * cdif.imag());	// (Fs*Fd')re * sin g * cos g
					psum->add(j++, cdif.power());										// |Fd * cos g|2
				}
			}
		}
	}

	return 0;
}

int			img_add_ewald_terms_new(Bimage* psum, Bimage* psec, Bimage* ppart, CTFparam& cp)
{
	long 			i, j, k, m, n, x, y, z, ix, iy, iz;
	double			sx, sy, sz, s2, phi;
	Complex<double>	cv, cup, clo, ca, csum, cdif;
	vector<Complex<double>>	ct(4);
	Vector3<double>	freq_scale(1.0/psum->real_size());
	Vector3<double>	h((psum->sizeX() - 1)/2, (psum->sizeY() - 1)/2, (psum->sizeZ() - 1)/2);

	if ( verbose & VERB_DEBUG ) {
		cout << "Calculating aberration terms" << endl;
		cout << endl;
	}
	
	for ( i=j=n=0; n<psum->images(); n++ ) {
		for ( z=iz=0; z<psum->sizeZ(); z++ ) {
			if ( z ) iz = psum->sizeZ() - z;
			sz = z;
			if ( z > h[2] ) sz -= psum->sizeZ();
			sz *= freq_scale[2];
			for ( y=iy=0; y<psum->sizeY(); y++ ) {
				if ( y ) iy = psum->sizeY() - y;
				sy = y;
				if ( y > h[1] ) sy -= psum->sizeY();
				sy *= freq_scale[1];
				for ( x=ix=0; x<psum->sizeX(); x++, i++ ) {
					if ( x ) ix = psum->sizeX() - x;
					sx = x;
					if ( x > h[0] ) sx -= psum->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					phi = atan2(sy,sx);
					cv = ppart->complex(i);
					ca = cp.calculate_complex(s2, phi);
					cup = psec->complex(i);
//					clo = psec->complex(psec->index(ix,iy,iz,n));
					clo = psec->complex(psec->index(ix,iy,iz,n)).conj();
					csum = cup + clo;
					cdif = cup - clo;
					ct[0] = csum*ca.imag();
					ct[1] = csum*ca.real();
					ct[2] = cdif*ca.real();
					ct[3] = cdif*ca.imag();
					for ( k=0; k<4; ++k )
						psum->add(j++, ct[k].real() * cv.real() + ct[k].imag() * cv.imag());
					for ( k=0; k<4; ++k )
						for ( m=k; m<4; ++m )
							psum->add(j++, ct[k].real() * ct[m].real() + ct[k].imag() * ct[m].imag());
				}
			}
		}
	}

	return 0;
}

Bimage*		particle_ewald_sum(Bparticle* part, Bimage* pref, double hires, FSI_Kernel* kernel, fft_plan planf)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( pref->fourier_type() != Standard ) {
		pref->fft();
		pref->phase_shift_to_origin();
	}

	long			ndone(0);
	Matrix3			mat;
	Bmicrograph*	mg = part->mg;
	Bimage*			psec;
	Bimage*			ppart;
	
	if ( !mg->ctf ) {
		cerr << "Error: The CTF parameters must be defined for micrograph " << mg->id << endl;
		bexit(-1);
	}

	double			wl = mg->ctf->lambda();

	Bimage*			psum = new Bimage(Float, 5, pref->sizeX(), pref->sizeY(), 1, 1);
//	Bimage*			psum = new Bimage(Float, 14, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(part->pixel_size);

	for ( ; part; part = part->next ) if ( part->sel > 0 ) {
		if ( part->def > 0 ) mg->ctf->defocus_average(part->def);
		ppart = read_img(mg->fpart, 1, part->id-1);
		ppart->sampling(part->pixel_size);
		ppart->origin(part->ori);
		ppart->view(part->view);
		ppart->fft(planf);
		ppart->phase_shift_to_origin();
		mat = part->view.matrix();
		part->ori[2] = 0;
//		cout << mat << part->ori << endl;
		psec = pref->central_section(mat, hires, kernel, wl);
		img_add_ewald_terms(psum, psec, ppart, *mg->ctf);
		delete ppart;
		delete psec;
		ndone++;
	}
	
	psum->image->select(ndone);

	return psum;
}

Bimage*		img_calculate_ewald_phase(Bimage* psum)
{
	long			i, j, k, n;
	vector<double>	t(2);
	Matrix			r(2,2);
	
//	Bimage*		pphi = new Bimage(Float, TComplex, psum->size(), psum->images());
	Bimage*		pphi = new Bimage(Float, TSimple, psum->size(), psum->images());
	pphi->sampling(psum->sampling(0));

	for ( n=i=j=0; n<pphi->images(); ++n ) {
		for ( k=0; k<pphi->image_size(); ++i, ++k ) {
			t[0] = 2 * (*psum)[j++];
			t[1] = 2 * (*psum)[j++];
			r[0][0] = (*psum)[j++];
			r[0][1] = r[1][0] = (*psum)[j++];
			r[1][1] = (*psum)[j++];
			r.singular_value_decomposition(t);
//			pphi->set(i, Complex<double>(-t[0],-t[1]));
			pphi->set(i, angle_set_negPI_to_PI(atan2(t[0],t[1])));
		}
	}
	
	pphi->statistics();
	
	return pphi;
}

Bimage*		img_calculate_ewald_phase_old2(Bimage* psum)
{
	long			i, j, ii, jj, k, n;
//	double			geven, gewald;
	vector<double>	t(4);
	Matrix			r(4,4);
	
//	Bimage*		pphi = new Bimage(Float, TComplex, psum->size(), psum->images());
	Bimage*		pphi = new Bimage(Float, 4, psum->size(), psum->images());
	pphi->sampling(psum->sampling(0));

	for ( n=i=j=0; n<pphi->images(); ++n ) {
		for ( k=0; k<pphi->image_size(); ++k ) {
			for ( ii=0; ii<4; ++ii )
				t[ii] = 2 * (*psum)[j++];
			for ( ii=0; ii<4; ++ii )
				for ( jj=ii; jj<4; ++jj )
					r[ii][jj] = r[jj][ii] = (*psum)[j++];
			r.singular_value_decomposition(t);
//			geven = atan2(t[1],t[0]);
//			gewald = atan2(t[2],t[0]);
//			pphi->set(i, Complex<double>(geven,gewald));
			t[3] = -t[3];
			for ( ii=0; ii<4; ++ii )
				pphi->set(i++, -t[ii]);
		}
	}
	
	pphi->statistics();
	
	return pphi;
}

Bimage*		img_calculate_ewald_phase_new(Bimage* psum)
{
	long			i, j, ii, jj, k, n;
	double			geven, gewald, v;
	vector<double>	t(4);
	Matrix			r(4,4);
	
	Bimage*		pphi = new Bimage(Float, TComplex, psum->size(), psum->images());
//	Bimage*		pphi = new Bimage(Float, 4, psum->size(), psum->images());
	pphi->sampling(psum->sampling(0));

	for ( n=i=j=0; n<pphi->images(); ++n ) {
		for ( k=0; k<pphi->image_size(); ++k ) {
			for ( ii=0; ii<4; ++ii )
				t[ii] = 2 * (*psum)[j++];
			for ( ii=0; ii<4; ++ii )
				for ( jj=ii; jj<4; ++jj )
					r[ii][jj] = r[jj][ii] = (*psum)[j++];
			r.singular_value_decomposition(t);
			geven = angle_set_negPI_to_PI(atan2(-t[1],-t[0]));
			v = cos(geven);
/*			if ( fabs(v) > 0.5 ) {
				gewald = angle_set_negPI_to_PI(asin(-t[2]/v));
			} else {
				v = sin(geven);
				gewald = angle_set_negPI_to_PI(asin(t[3]/v));
			}*/
			if ( fabs(v) > 0.5 ) {
				gewald = angle_set_negPI_to_PI(acos(-t[0]/v));
			} else {
				v = sin(geven);
				gewald = angle_set_negPI_to_PI(acos(-t[1]/v));
			}
			pphi->set(i++, Complex<double>(geven,gewald));
		}
	}
	
	pphi->statistics();
	
	return pphi;
}

/**
@brief 	Calculates the ewald phase from a set of particles.
@param 	*project 	project parameter structure.
@param	*pref		reference map.
@param 	hires	 	high resolution limit (angstrom).
@param 	*kernel	 	frequency space interpolation kernel lookup table.
@return	Bimage*		aberration phase differences as a complex image.

	The return image is complex with the even and odd phase differences
	in real and imaginary parts, respectively.

**/
Bimage*		project_ewald_phase(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	long			nsel = project_count_mg_part_selected(project);
	
	if ( verbose ) {
		cout << "Calculating the Ewald phase from the difference between particles and CTF-applied reference projections:" << endl;
		cout << "Reference map:                  " << pref->file_name() << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Selected particles:             " << nsel << endl;
	}
	
	if ( pref->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the reference map" << endl;
		pref->fft();
	}

	pref->phase_shift_to_origin();

	long			nmg(0), nogrp(0);
	Matrix3			mat;
	Bparticle*		part = part_find_first(project);
	
	part_select_micrographs_with_selected_particles(project);
	
	if ( verbose )
		cout << "Setting up the micrograph array" << endl;
	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);
	
	map<string,long>	optgrp;
	for ( long i=0; i<nmg; ++i )
		if ( mgarr[i]->ctf )
			if ( optgrp.find(mgarr[i]->ctf->identifier()) == optgrp.end() )
				optgrp[mgarr[i]->ctf->identifier()] = 0;

	for ( auto& g: optgrp ) g.second = nogrp++;	// Ensure the order is encoded properly
	
	if ( nogrp < 1 ) nogrp = 1;
	
	Bimage*			psum = new Bimage(Float, 5, pref->sizeX(), pref->sizeY(), 1, nogrp);
//	Bimage*			psum = new Bimage(Float, 14, pref->sizeX(), pref->sizeY(), 1, nogrp);
	psum->sampling(part->pixel_size);

	fft_plan		planf = fft_setup_plan(psum->size(), FFTW_FORWARD, 1);
	
#ifdef HAVE_GCD
	__block	long	ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nmg, dispatch_get_global_queue(0, 0), ^(size_t i){
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp.at(mgarr[i]->ctf->identifier());
		Bimage*		p1 = particle_ewald_sum(mgarr[i]->part, pref, hires, kernel, planf);
		dispatch_sync(myq, ^{
//			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		});
	});
#else
	long			ndone(0);
#pragma omp parallel for
	for ( long i=0; i<nmg; ++i ) {
		long		n(0);
		if ( mgarr[i]->ctf ) n = optgrp[mgarr[i]->ctf->identifier()];
		Bimage*		p1 = particle_ewald_sum(mgarr[i]->part, pref, hires, kernel, planf);
	#pragma omp critical
		{
//			cout << mgarr[i]->id << tab << mgarr[i]->ctf->identifier() << tab << n << tab << p1->image->select() << endl;
			ndone += p1->image->select();
			psum->add(n, p1);
			psum->image[n].select(psum->image[n].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		}
	}
#endif
	
	fft_destroy_plan(planf);

	if ( verbose & VERB_RESULT )
		cerr << endl;

	if ( verbose ) {
		cout << "#\tGroup\tParticles" << endl;
		for ( auto g: optgrp )
			cout << g.second << tab << g.first << tab << psum->image[g.second].select() << endl;
		cout << endl;
	}
	
//	write_img("t.grd", psum, 0);

	Bimage*		pphi = img_calculate_ewald_phase(psum);

	delete[] mgarr;
	delete psum;

	return pphi;
}

int			img_add_correlation_terms(Bimage* psum, Bimage* psec, Bimage* ppart, double hires)
{
	if ( hires <= psec->image->sampling()[0] ) hires = 2*psec->image->sampling()[0];
	
	long 			i, j, n, x, y, z;
	double			sx, sy, sz, s2;
	double			s2max(1/(hires*hires));
	Complex<double>	cv, cr, ca;
	Vector3<double>	freq_scale(1.0/psum->real_size());
	Vector3<double>	h((psum->sizeX() - 1)/2, (psum->sizeY() - 1)/2, (psum->sizeZ() - 1)/2);

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG img_add_correlation_terms: Calculating correlation terms" << endl;
		cout << endl;
	}
	
	for ( i=j=n=0; n<psum->images(); n++ ) {
		for ( z=0; z<psum->sizeZ(); z++ ) {
			sz = z;
			if ( z > h[2] ) sz -= psum->sizeZ();
			sz *= freq_scale[2];
			for ( y=0; y<psum->sizeY(); y++ ) {
				sy = y;
				if ( y > h[1] ) sy -= psum->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<psum->sizeX(); x++, i++ ) {
					sx = x;
					if ( x > h[0] ) sx -= psum->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy + sz*sz;
					if ( s2 <= s2max ) {
						cv = ppart->complex(i);
						cr = psec->complex(i);
						psum->add(j++, cr.real() * cv.real() + cr.imag() * cv.imag());	// (Fs*F')re
						psum->add(j++, cv.power());										// |F|2
						psum->add(j++, cr.power());										// |Fs|2
					} else {
						j += 3;
					}
				}
			}
		}
	}

	return 0;
}

Bimage*		particle_ewald_correlation(Bparticle* part, Bimage* pref, double hires, FSI_Kernel* kernel, fft_plan planf)
{
	if ( !kernel ) {
		error_show("Error in project_correlation_sum: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( pref->fourier_type() != Standard ) {
		pref->fft();
		pref->phase_shift_to_origin();
	}

	bool			invert(0);
	long			ndone(0);
	Matrix3			mat;
	Bmicrograph*	mg = part->mg;
	Bimage*			psec;
	Bimage*			ppart;
/*
	if ( !mg->ctf ) {
		cerr << "Error: The CTF parameters must be defined for micrograph " << mg->id << endl;
		bexit(-1);
	}
*/
	Bimage*			psum = new Bimage(Float, 3, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(part->pixel_size);

	for ( ; part; part = part->next ) if ( part->sel > 0 ) {
		ppart = read_img(mg->fpart, 1, part->id-1);
		ppart->sampling(part->pixel_size);
		ppart->origin(part->ori);
		ppart->view(part->view);
		ppart->fft(planf);
		ppart->phase_shift_to_origin();
		mat = part->view.matrix();
//		part->ori[2] = 0;
//		cout << mat << part->ori << endl;
		psec = pref->central_section(mat, hires, kernel, 0);
		if ( mg->ctf ) {
			if ( part->def > 0 )
				mg->ctf->defocus_average(part->def);
			if ( mg->ctf->defocus_average() > 0 )
				img_ctf_apply(psec, *mg->ctf, 2, 0.1, 0, hires, invert);
		}
		img_add_correlation_terms(psum, psec, ppart, hires);
		delete ppart;
		delete psec;
		ndone++;
	}
	
	psum->image->select(ndone);

	return psum;
}

Bimage*		img_calculate_ewald_correlation(Bimage* psum)
{
	long			i, j, k, n;
	double			cor, div;
	
	Bimage*			pcor = new Bimage(Float, TSimple, psum->size(), psum->images());
	pcor->sampling(psum->sampling(0));

	for ( n=i=j=0; n<pcor->images(); ++n ) {
		for ( k=0; k<pcor->image_size(); ++i, ++k ) {
			cor = (*psum)[j++];
			div = (*psum)[j++];
			div *= (*psum)[j++];
			if ( div )
				pcor->set(i, cor/sqrt(div));
		}
	}
	
	pcor->statistics();
	
	return pcor;
}


/**
@brief 	Calculates the correlation between particle images  and the corresponding central sections.
@param 	*project 	project parameter structure.
@param	*pref		reference map.
@param 	hires	 	high resolution limit (angstrom).
@param 	*kernel	 	frequency space interpolation kernel lookup table.
@return	Bimage*		correlation images.

	The return image is the equivalent of FRC.

**/
Bimage*		project_ewald_correlation(Bproject* project, Bimage* pref, double hires, FSI_Kernel* kernel)
{
	if ( !kernel ) {
		error_show("Error in project_ewald_correlation: No kernel defined for projection!", __FILE__, __LINE__);
		return NULL;
	}
	
	if ( hires < pref->image->sampling()[0] ) hires = 2*pref->image->sampling()[0];
	
	long			nsel = project_count_mg_part_selected(project);
	
	long			set_size(nsel/system_processors()+1);
	if ( set_size < 1 ) set_size = 1;

	long			nset = part_select_sets(project, set_size, 0);

	if ( verbose ) {
		cout << "Calculating the correlation between particles and CTF-applied reference projections:" << endl;
		cout << "Reference map:                  " << pref->file_name() << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Selected particles:             " << nsel << endl;
		cout << "Number of sets:                 " << nset << endl;
		cout << "Set size:                       " << set_size << endl;
	}
	
	if ( pref->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the reference map" << endl;
		pref->fft();
	}

	pref->phase_shift_to_origin();

	Bimage*			psum = new Bimage(Float, 3, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(pref->image->sampling()[0], pref->image->sampling()[1], 1);

	fft_plan		planf = fft_setup_plan(psum->size(), FFTW_FORWARD, 1);

	if ( verbose )
		cout << "Starting analysis" << endl;
		
#ifdef HAVE_GCD
	__block	long	ndone(0);
	dispatch_queue_t 	myq = dispatch_queue_create(NULL, NULL);
	dispatch_apply(nset, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bparticle*	partlist = project_selected_partlist(project, i+1, 0);
		Bimage*		p1 = particle_ewald_correlation(partlist, pref, hires, kernel, planf);
		particle_kill(partlist);
		dispatch_sync(myq, ^{
			ndone += p1->image->select();
			psum->add(0, p1);
			psum->image[0].select(psum->image[0].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		});
	});
#else
	long			ndone(0);
#pragma omp parallel for
	for ( long i=0; i<nset; i++ ) {
		Bparticle*	partlist = project_selected_partlist(project, i+1, 0);
		Bimage*		p1 = particle_ewald_correlation(partlist, pref, hires, kernel, planf);
		particle_kill(partlist);
	#pragma omp critical
		{
			ndone += p1->image->select();
			psum->add(0, p1);
			psum->image[0].select(psum->image[0].select() + p1->image->select());
			delete p1;
			if ( verbose & VERB_RESULT )
				cerr << "Complete:                       " << setprecision(3)
					<< ndone*100.0/nsel << " %    \r" << flush;
		}
	}
#endif

	fft_destroy_plan(planf);

	if ( verbose & VERB_RESULT )
		cerr << endl;

	Bimage*		pphi = img_calculate_ewald_correlation(psum);

	delete psum;

	return pphi;
}
