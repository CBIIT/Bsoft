/**
@file	mg_extract.cpp
@brief	Functions to extract particles from micrographs
@author	Bernard Heymann
@date	Created: 20040406
@date	Modified: 20200616
**/

#include "mg_extract.h"
#include "mg_select.h"
#include "linked_list.h"
#include "spline.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			mask_filament_particles(Bstring& filename, Bparticle* partlist, int mask_width);
int			rotate_mask_filament_particles(Bstring& filename, Bparticle* partlist, 
				int rotation_axis, int back_flag, int mask_width);

/**
@brief 	Sets up individual particle file names for extraction.
@param 	*part		particle list.
@param 	filename	particle file name base.
@param 	partpath	particle path with optional wild card characters.
@return long		number of particles.

	The particle path can have wildcard characters ('?') that specify
	the number of digits to use to insert numbers in the path and
	the file name.

**/
int			particle_setup_filenames(Bparticle* part, Bstring filename, Bstring partpath)
{
	long				n(0);
	char				format[32] = "%03d";
	Bstring				path;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG particle_setup_filenames: base=" << filename << endl;

	if ( partpath.length() ) {
		if ( !partpath.empty() )
			if ( partpath[-1] == '/' ) partpath = partpath.truncate(1);
		n = partpath.count('?');
		if ( n > 0 ) {
			partpath = partpath.remove('?');
			snprintf(format, 32, "%%0%ldd", n);
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG particle_setup_filenames: path=" << partpath << endl;
	}
		
	if ( partpath.length() > 1 ) if ( filename.contains("/") ) filename = filename.post_rev('/');
	for ( ; part; part = part->next ) {
		if ( partpath.length() > 1 ) {
			if ( n < 1 ) path = partpath;
			else path = partpath + Bstring(part->id, format);
			part->fpart = path + "/" + filename;
		} else part->fpart = filename;
		mkdir(path.c_str(), (mode_t)0755);
		part->fpart = part->fpart.pre_rev('.') + 
				Bstring(part->id, format) + "." + part->fpart.post_rev('.');
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG particle_setup_filenames: id=" << part->id << " path=" << part->fpart << endl;
	}
	
	return 0;
}

/**
@brief 	Sets up individual filament file names for extraction.
@param 	*fil		filament list.
@param 	filename	filament file name base.
@param 	filpath	filament path with optional wild card characters.
@return long				number of filaments.

	The filament path can have wildcard characters ('?') that specify
	the number of digits to use to insert numbers in the path and
	the file name.

**/
int			filament_setup_filenames(Bfilament* fil, Bstring filename, Bstring filpath)
{
	long				n(0);
	char				format[32] = "%03d";
	Bstring				path;

	if ( filpath.length() ) {
		if ( !filpath.empty() )
			if ( filpath[-1] == '/' ) filpath = filpath.truncate(1);
		n = filpath.count('?');
		if ( n > 0 ) {
			filpath = filpath.remove('?');
			snprintf(format, 32, "%%0%ldd", n);
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG filament_setup_filenames: path=" << filpath << endl;
	}
		
	if ( filpath.length() > 1 ) if ( filename.contains("/") ) filename = filename.post_rev('/');
	for ( ; fil; fil = fil->next ) {
		if ( filpath.length() > 1 ) {
			if ( n < 1 ) path = filpath;
			else path = filpath + Bstring(fil->id, format);
			fil->ffil = path + "/" + filename;
		} else fil->ffil = filename;
		mkdir(path.c_str(), (mode_t)0755);
		fil->ffil = fil->ffil.pre_rev('.') + 
				Bstring(fil->id, format) + "." + fil->ffil.post_rev('.');
	}
	
	return 0;
}

/**
@brief 	Calculates a spline curve from a set of filament nodes.
@param 	*fnode			node list.
@param 	&nspline		number of elements in the spline array.
@return Vector3<double>*	array of spline coordinates.

	The node coordinates are copied into a new array used to calculate the spline array.

**/
Vector3<double>*	vector3_spline_from_nodes(Bfilnode* fnode, long& nspline)
{
	int			n;
	Bfilnode*	fn;
	
	for ( n=0, fn = fnode; fn; fn=fn->next ) n++;
	
	Vector3<double>*	coords = new Vector3<double>[n];
	
	for ( n=0, fn = fnode; fn; fn=fn->next, n++ ) coords[n] = fn->loc;

	Vector3<double>*	spline = vector3_catmull_rom_spline(n, coords, nspline);
	
	delete[] coords;
	
	return spline;
}

/**
@brief 	Extracts particle images from micrographs defined in a project.
@param 	*project	micrograph project.
@param 	scale		scale to extract (usually 1).
@param 	back_flag	background correction flag.
@param 	norm_flag	normalization flag.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		value to fill in new regions.
@param 	mask_width	filament mask width, if 0, don't apply.
@param 	split		flag to split images into separate files.
@param 	&partbase	file name base.
@param 	&partpath	path to particle file.
@param 	&partext	particle file extension.
@return long		number of particles, <0 on error.
**/
long		project_extract_particles(Bproject* project, double scale,
				int back_flag, int norm_flag, int fill_type, double fill, int mask_width,
				int split, Bstring& partbase, Bstring& partpath, Bstring& partext)
{
	long				i(1), npart(0);
	Bstring				filename;
	Bfield*				field = project->field;
	Bmicrograph*		mg;
	Breconstruction*	rec = project->rec;
	Bimage*				p = NULL;
	
	if ( !partpath.empty() ) if ( partpath[-1] != '/' ) partpath += "/";
	
	if ( partext.length() < 1 ) partext = "mrc";
	
	long				nmg = project_count_micrographs(project);
	long				nrec = project_count_reconstructions(project);
	
	long				multi_img(0);

	Vector3<long>		size;
	if ( project->select && rec ) size = rec->box_size;
	else if ( field && field->mg ) size = field->mg->box_size;
	if ( size.volume() < 1 ) {
		cerr << "Error in project_extract_particles: Box size = " << size << endl;
		return -1;
	}

	if ( verbose ) {
		cout << "Extracting particles:" << endl;
		cout << "Particle size:                   " << size << endl;
		cout << "Particle scale:                  " << scale << endl;
		if ( mask_width )
			cout << "Filament mask width:             " << mask_width << endl;
		if ( back_flag ) cout << "With background correction" << endl;
		if ( norm_flag ) cout << "With normalization" << endl;
		if ( partbase.length() )
			cout << "Particle base name:              " << partbase << endl;
		if ( partext.length() )
			cout << "Particle extension:              " << partext << endl;
		if ( partpath.length() )
			cout << "Particle path:                   " << partpath << endl;
		if ( project->select > 1 )
			cout << "Projecting" << endl;
		cout << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg=field->mg; mg; mg=mg->next, i++ ) {
				if ( mg->next ) {
					if ( mg->fmg == mg->next->fmg ) multi_img = 1;
					else multi_img = 0;
				}
//				cout << mg->fmg << endl;
				npart = particle_count(mg->part);
				if ( npart > 0 ) {
					if ( partbase.length() ) {
						if ( nmg == 1 )
							filename = partbase;
						else
							filename = partbase + Bstring(i, "_%03d");
					} else {
						if ( mg->fpart.length() < 1 ) mg->fpart = mg->fmg;
						if ( mg->fmg.length() < 1 ) mg->fpart = mg->ffil;
						filename = mg->fpart.pre_rev('.');
						if ( mg->fmg == mg->fpart ) {
							if ( multi_img ) filename = filename + Bstring(mg->img_num+1, "_part%03d");
							else filename = filename + "_part";
						}
						if ( mg->ffil == mg->fpart ) {
							filename = filename + "_fil";
						}
					}
					if ( filename.length() < 1 )
						cerr << "Error in project_extract_particles: No particle file name!" << endl;
					filename = filename + "." + partext;
					if ( partpath.length() > 1 ) {
						if ( filename.contains("/") ) filename = filename.post_rev('/');
						filename = partpath + filename;
					}
					if ( split ) {
						particle_setup_filenames(mg->part, filename, partpath);
						mg->fpart = 0;
					} else {
						mg->fpart = filename;
						if ( partpath.length() > 1 ) mkdir(partpath.c_str(), (mode_t)0755);
//						cout << mg->fmg << " ---> " << mg->fpart << endl;
					}
					if ( mg->fmg.length() ) p = read_img(mg->fmg, 1, mg->img_num);
					else if ( mg->ffil.length() ) p = read_img(mg->ffil, 1, 0);
					if ( !p ) {
						error_show(mg->fmg.c_str() , __FILE__, __LINE__);
						return -1;
					}
					if ( fill_type == FILL_AVERAGE ) fill = p->average();
					if ( fill_type == FILL_BACKGROUND ) {
						if ( fabs(p->background(long(0))) < 1e-6 )
							p->calculate_background();
						fill = p->background(long(0));
					}
					micrograph_extract_particles(mg, p, scale, back_flag, norm_flag, fill, mask_width);
					delete p;
				}
			}
		}
	} else if ( project->select == 1 ) {
		for ( rec=project->rec; rec; rec=rec->next, i++ ) {
			npart = particle_count(rec->part);
			if ( npart > 0 ) {
				if ( partbase.length() ) {
					if ( nrec == 1 )
						filename = partbase;
					else
						filename = partbase + Bstring(i, "_%03d");
				} else {
					if ( rec->fpart.length() < 1 ) rec->fpart = rec->frec;
					if ( rec->frec.length() < 1 ) rec->fpart = rec->ffil;
					filename = rec->fpart.pre_rev('.');
					if ( rec->frec == rec->fpart ) {
//						if ( rec->img_num ) filename = filename + Bstring(rec->img_num, "_part%03d");
//						else filename = filename + "_part";
						filename = filename + "_part";
					}
					if ( rec->ffil == rec->fpart ) {
						filename = filename + "_fil";
					}
				}
				if ( filename.length() < 1 )
					cerr << "Error in project_extract_particles: No particle file name!" << endl;
				filename = filename + "." + partext;
				if ( partpath.length() > 1 ) {
					if ( filename.contains("/") ) filename = filename.post_rev('/');
					filename = partpath + filename;
				}
				if ( split ) {
					particle_setup_filenames(rec->part, filename, partpath);
					rec->fpart = 0;
				} else {
					rec->fpart = filename;
					if ( partpath.length() > 1 ) mkdir(partpath.c_str(), (mode_t)0755);
				}
				if ( rec->frec.length() ) p = read_img(rec->frec, 1, 0);
				else if ( rec->ffil.length() ) p = read_img(rec->ffil, 1, 0);
				if ( !p ) {
					error_show(rec->frec.c_str() , __FILE__, __LINE__);
					return -1;
				}
				if ( fill_type == FILL_AVERAGE ) fill = p->average();
				if ( fill_type == FILL_BACKGROUND ) {
					if ( fabs(p->background(long(0))) < 1e-6 )
						p->calculate_background();
					fill = p->background(long(0));
				}
				reconstruction_extract_particles(rec, p, scale, back_flag, norm_flag, fill, mask_width);
				delete p;
			}
		}
	} else {
		if ( rec->frec.length() ) p = read_img(rec->frec, 1, 0);
		if ( !p ) {
			error_show(rec->frec.c_str() , __FILE__, __LINE__);
			return -1;
		}
		if ( fill_type == FILL_AVERAGE ) fill = p->average();
		if ( fill_type == FILL_BACKGROUND ) {
			if ( fabs(p->background(long(0))) < 1e-6 )
				p->calculate_background();
			fill = p->background(long(0));
		}
		mg = field_find_zero_tilt_mg(project->field);
		if ( mg->part ) particle_kill(mg->part);
		mg->part = reconstruction_project_extract_particles(rec, p, scale, back_flag, norm_flag, fill, mask_width);
		if ( p->file_name().length() ) mg->fpart = p->file_name();
		delete p;
	}
	
	return npart;
}

/**
@brief 	Extracts particle images from a micrograph.
@param 	*mg			micrograph parameters.
@param 	*p			micrograph image.
@param 	scale		scale to extract (usually 1).
@param 	back_flag	background correction flag.
@param 	norm_flag	normalization flag.
@param 	fill		value to fill in new regions.
@param 	mask_width	filament mask width, if 0, don't apply.
@return long		number of particles.
**/
long		micrograph_extract_particles(Bmicrograph* mg, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width)
{
	long		npart = particle_count(mg->part);
	
	if ( npart < 1 ) return 0;
	
	if ( verbose )
		cout << mg->id << tab << npart << endl;
	
	if ( mg->pixel_size[0] > 0 ) p->sampling(mg->pixel_size);

	Bimage*		ppart = particle_extract(mg->part, mg->bad, p, mg->box_size,
					scale, mg->bad_radius, back_flag, norm_flag, fill, mask_width);

	if ( ppart ) {
		write_img(mg->fpart, ppart, 0);
		delete ppart;
	}
	
	return npart;
}

/**
@brief 	Extracts particle images from a reconstruction.
@param 	*rec		reconstruction parameters.
@param 	*p			micrograph image.
@param 	scale		scale to extract (usually 1).
@param 	back_flag	background correction flag.
@param 	norm_flag	normalization flag.
@param 	fill		value to fill in new regions.
@param 	mask_width	filament mask width, if 0, don't apply.
@return long		number of particles.
**/
long		reconstruction_extract_particles(Breconstruction* rec, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width)
{
	long		npart = particle_count(rec->part);
	
	if ( npart < 1 ) return 0;
	
	if ( verbose )
		cout << rec->id << tab << npart << endl;

	if ( rec->voxel_size[0] > 0 ) p->sampling(rec->voxel_size);
	
	Bimage*		ppart = particle_extract(rec->part, rec->bad, p, rec->box_size,
					scale, rec->bad_radius, back_flag, norm_flag, fill, mask_width);

	if ( ppart ) {
		write_img(rec->fpart, ppart, 0);
		delete ppart;
	}
	
	return npart;
}

/**
@brief 	Extracts and projects 3D particle images from a reconstruction and return a 2D list.
@param 	*rec		reconstruction parameters.
@param 	*p			micrograph image.
@param 	scale		scale to extract (usually 1).
@param 	back_flag	background correction flag.
@param 	norm_flag	normalization flag.
@param 	fill		value to fill in new regions.
@param 	mask_width	filament mask width, if 0, don't apply.
@return Bparticle*		list of 2D particles.

	The file name of the new particles is returned in the image file name.
	
**/
Bparticle*		reconstruction_project_extract_particles(Breconstruction* rec, Bimage* p, double scale,
				int back_flag, int norm_flag, double fill, int mask_width)
{
	long		npart = particle_count(rec->part);
	
	if ( npart < 1 ) return 0;
	
	if ( verbose )
		cout << rec->id << tab << npart << endl;

	if ( rec->voxel_size[0] > 0 ) p->sampling(rec->voxel_size);
	
	Bparticle*	part = rec->part;

	Bstring		filename = rec->frec.pre_rev('.') + "_part.mrc";
	if ( part->fpart.length() && part->fpart.contains("/") )
		filename = part->fpart.pre_rev('/') + "/" + filename;

	// Set up particle list for projections
	Bparticle*	part2D = particle_copy(rec->part);

	// Delete filenames to avoid splitting image file
	for ( part = part2D; part; part = part->next )
		part->fpart = 0;

	Bimage*		ppart = particle_extract(part2D, rec->bad, p, rec->box_size,
					scale, rec->bad_radius, back_flag, norm_flag, fill, mask_width);

	// Set the z-related parameters
	for ( part = part2D; part; part = part->next ) {
		part->pixel_size[2] = 1;
		part->loc[2] = 0;
		part->ori[2] = 0;
	}

	if ( ppart ) {
		Bimage* 	pproj = ppart->project('z');
		write_img(filename, pproj, 0);
		delete ppart;
		delete pproj;
		p->file_name(filename.c_str());
	} else {
		p->file_name("");
	}
	
	return part2D;
}

/**
@brief 	Extracts particle images from an image.
@param 	*particles	particle parameters.
@param 	*bad_areas	bad area parameters.
@param 	*p			image.
@param 	size		size of box to extract.
@param 	scale		scale to extract (usually 1).
@param 	bad_radius	radius of bad area.
@param 	back_flag	background correction flag.
@param 	norm_flag	normalization flag.
@param 	fill		value to fill in new regions.
@param 	mask_width	filament mask width, if 0, don't apply.
@return Bimage*		multi-particle image.

	If the background flag is specified, the particle mask is set in the 
	background, defined as outside the inscribing circle.
	The mask is set within every bad area in the micrograph and
	transferred to the mask for a particle where it overlaps.
	
**/
Bimage*		particle_extract(Bparticle* particles, Bbadarea* bad_areas, Bimage* p, 
				Vector3<long> size, double scale, double bad_radius,
				int back_flag, int norm_flag, double fill, int mask_width)
{
	int					split(0);
	if ( particles->fpart.length() ) split = 1;		// Flag to write separate particle images
	
	if ( scale < 0.1 ) scale = 1;
	
//	if ( mask_width < 1 ) mask_width = (p->sizeX() > p->sizeY())? p->sizeX(): p->sizeY();
	
	long				id, i, j, k, c, n, x, y, z;
	long				oldx, oldy, oldz, oldn(0);
	double				d(0), dh(0), halfwidth(mask_width/2.0), ori_tol;
	double				iscale(1/scale);
	Bparticle*			part = particles;
	Bbadarea*			bad;
	Bimage*				pmask = NULL;
	Bimage*				ppart = NULL;
	Bimage*				partmask = NULL;
	Vector3<long>		box, radius;
	Vector3<double>		u, v;
	Euler				euler;

	long				npart = particle_count(particles);

	if ( npart < 1 ) return NULL;
	
	Vector3<double>		pixel_size;
	
	if ( part->mg ) pixel_size = part->mg->pixel_size;
	else if ( part->rec ) pixel_size = part->rec->voxel_size;
	if ( pixel_size[0] < 0.01 ) pixel_size = p->image->sampling();
	
	if ( pixel_size[0] < 0.01 ) {
		cerr << "Error in particle_extract: Pixel size is too small! " << pixel_size << endl;
		return NULL;
	}
	
	pixel_size *= iscale;
	box = size; // Box size
	box = box.min(p->size()*scale);
	box = box.max(1);
	radius = box/2;
	ori_tol = radius.length()/4;
	if ( ori_tol > 2 ) ori_tol = 2;
	
	if ( box.volume() < 1 ) {
		cerr << "Error: Particle box size is zero!" << endl;
		return NULL;
	}
	
	if ( box[2] == 1 ) pixel_size[2] = 1;
	
	if ( verbose )
		cout << "Extracting " << npart << " particles from " << p->file_name() << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of particles:            " << npart << endl;
		cout << "Particle box size:              " << box << endl;
		cout << "Particle scale:                 " << scale << endl;
		cout << "Particle radii:                 " << radius << endl;
		cout << "Bad region radius:              " << bad_radius << endl;
		cout << "Mask width:                     " << mask_width << endl;
		cout << "Fill value:                     " << fill << endl;
		if ( back_flag ) cout << "Background corrected" << endl;
		if ( norm_flag ) cout << "Normalized" << endl;
		cout << endl;
	}

	if ( bad_areas && bad_radius > 1 ) {
		pmask = new Bimage(UCharacter, p->compound_type(), p->size(), p->images());
		for ( bad = bad_areas; bad; bad = bad->next )
			pmask->sphere(bad->loc, bad_radius, 0, FILL_USER, 1);
	}
	
	if ( split ) {
		partmask = new Bimage(UCharacter, TSimple, box, 1);
	} else {
		ppart = new Bimage(p->data_type(), p->compound_type(), box, npart);
		ppart->fill(fill);
		partmask = new Bimage(UCharacter, TSimple, box, npart);
	}
	
	for ( id=1, i=n=0, part = particles; part; part = part->next, n++, id++ ) {
		if ( split ) {
			i = n = 0;
			ppart = new Bimage(p->data_type(), p->compound_type(), box, 1);
			ppart->fill(fill);
		}
//		cout << "extracting " << id << endl;
		part->id = id;
		if ( part->ori.distance(radius) > ori_tol ) part->ori = radius;
		part->pixel_size = pixel_size;
		ppart->sampling(pixel_size);
		ppart->origin(n, part->ori);
		ppart->image[n].view(part->view);
		euler = Euler(part->view);
		u[0] = cos(euler.psi());
		u[1] = -sin(euler.psi());
		for ( z=0; z<ppart->sizeZ(); z++ ) {
			oldz = (z - radius[2])*iscale + (long) part->loc[2];
			v[2] = z - ppart->image[n].origin()[2];
			for ( y=0; y<ppart->sizeY(); y++ ) {
				oldy = (y - radius[1])*iscale + (long) part->loc[1];
				v[1] = y - ppart->image[n].origin()[1];
				for ( x=0; x<ppart->sizeX(); x++, i++ ) {
					oldx = (x - radius[0])*iscale + (long) part->loc[0];
					v[0] = x - ppart->image[n].origin()[0];
					if ( oldx >= 0 && oldx < p->sizeX() && oldy >= 0 && oldy < p->sizeY() && oldz >= 0 && oldz < p->sizeZ()) {
						j = (((oldn*p->sizeZ() + oldz)*p->sizeY() + oldy)*p->sizeX() + oldx)*p->channels();
						if ( halfwidth > 0 ) dh = (u.cross(v)).length() - halfwidth;
						if ( back_flag ) d = (v/ppart->size()).length() - 0.5;
//						cout << "d=" << d << endl;
						if ( ( pmask && (*pmask)[j] > 0 ) || ( d > 0 ) || ( dh > 0 ) ) {
							partmask->set(i, 0);
							if ( back_flag ) for ( c=0, k=i*p->channels(); c<p->channels(); c++ )
//								ppart->set(k++, (*p)[j++]);
								ppart->set(k++, p->average(c, oldx, oldy, oldz, oldn, iscale));
						} else {
							partmask->set(i, 1);
							for ( c=0, k=i*p->channels(); c<p->channels(); c++ )
//								ppart->set(k++, (*p)[j++]);
								ppart->set(k++, p->average(c, oldx, oldy, oldz, oldn, iscale));
						}
					}
				}
			}
		}
		if ( split ) {
			if ( back_flag ) ppart->correct_background(partmask, 1);	// The mask is 1 for the background
			else ppart->calculate_background();
			ppart->statistics();
			if ( ppart->standard_deviation() > 0 ) {
				if ( norm_flag )
					ppart->normalize(ppart->average(), ppart->standard_deviation(), 0);
			} else {
				part->sel = 0;
			}
			write_img(part->fpart, ppart, 0);
			delete ppart;
			ppart = NULL;
		}
	}
		
	delete pmask;
	
//	cout << "pixel size = " << ppart->sampling(0) << endl;
	
//	write_img("part.mrc", ppart);
//	write_img("mask.mrc", partmask);
	
	if ( !split ) {
		if ( back_flag ) ppart->correct_background(partmask, 1);
		else ppart->calculate_background();	
		ppart->statistics();
		if ( norm_flag ) ppart->normalize(ppart->average(), ppart->standard_deviation(), 0);
		for ( n=0, part = particles; part; part = part->next, n++ ) {
			if ( ppart->image[n].standard_deviation() <= 0 ) part->sel = 0;
		}
	}
	
	delete partmask;	
	
	return ppart;
}

/**
@brief 	Extracts gold particle images from a micrograph.
@param 	*mg		micrograph parameters.
@param 	*p			micrograph image.
@param 	radius		radius of gold particle.
@return Bimage*				multi-particle image.
**/
Bimage*		micrograph_extract_gold(Bmicrograph* mg, Bimage* p, double radius)
{
	if ( radius < 1 ) radius = mg->mark_radius;
	
	Bimage*				pmark = marker_extract_gold(mg->mark, p, mg->img_num, radius);
	
	if ( mg->pixel_size[0] > 0 ) pmark->sampling(mg->pixel_size);

	return pmark;
}

/**
@brief 	Extracts gold particle images from a reconstruction.
@param 	*rec	reconstruction parameters.
@param 	*p				reconstruction image.
@param 	radius			radius of gold particle.
@return Bimage*					multi-particle image.
**/
Bimage*		reconstruction_extract_gold(Breconstruction* rec, Bimage* p, double radius)
{
	if ( radius < 1 ) radius = rec->mark_radius;
	
	Bimage*				pmark = marker_extract_gold(rec->mark, p, 0, radius);
	
	if ( rec->voxel_size[0] > 0 ) pmark->sampling(rec->voxel_size);

	return pmark;
}

/**
@brief 	Extracts gold particle images from an image.
@param 	*marker_list	marker parameters.
@param 	*p				image.
@param 	img_num				sub-image number.
@param 	radius			radius of gold particle.
@return Bimage*					multi-particle image.
**/
Bimage*		marker_extract_gold(Bmarker* marker_list, Bimage* p, int img_num, double radius)
{
	long				i, j, c, n, x, y, z;
	long				oldx, oldy, oldz;
	Bmarker*			mark;
	long				nmark;
	Vector3<long>		box, iradius;
	
	for ( nmark=0, mark=marker_list; mark; mark=mark->next ) nmark++;
	
	if ( nmark < 1 ) return NULL;
		
	iradius[0] = iradius[1] = iradius[2] = (int)(2*radius); // Half box size
	box[0] = box[1] = box[2] = (int)(4*radius); // Box size
	if ( p->sizeZ() < 2 ) {
		iradius[2] = 0;
		box[2] = 1;
	}

	if ( verbose & VERB_PROCESS )
		cout << "Extracting " << nmark << " markers from image " << p->file_name()
			<< " (" << img_num << ") with box size " << box << endl;

	Bimage*				pmark = new Bimage(p->data_type(), p->compound_type(), box[0], box[1], box[2], nmark);
	pmark->sampling(p->sampling(0));
	
	for ( n=0, mark=marker_list; mark && n<pmark->images(); mark=mark->next, n++ ) {
		pmark->origin(n, iradius);	// Set the origin to the nominal center
		for ( z=0; z<pmark->sizeZ(); z++ ) {
			oldz = z - iradius[2] + (int) mark->loc[2];
			for ( y=0; y<pmark->sizeY(); y++ ) {
				oldy = y - iradius[1] + (int) mark->loc[1];
				for ( x=0; x<pmark->sizeX(); x++ ) {
					oldx = x - iradius[0] + (int) mark->loc[0];
					if ( oldx >= 0 && oldx < p->sizeX() && oldy >= 0 && oldy < p->sizeY() && oldz >= 0 && oldz < p->sizeZ()) {
						j = (((img_num*p->sizeZ() + oldz)*p->sizeY() + oldy)*p->sizeX() + oldx)*p->channels();
						i = (((n*pmark->sizeZ() + z)*pmark->sizeY() + y)*pmark->sizeX() + x)*pmark->channels();
						for ( c=0; c<p->channels(); c++, i++, j++ ) pmark->set(i, (*p)[j]);
					}
				}
			}
		}
	}
	
	pmark->statistics();
	
	pmark->correct_background();
	
	pmark->normalize(pmark->average(), pmark->standard_deviation(), 0);
	
	return pmark;
}

/**
@brief 	Extracts filament images from micrographs defined in a project.
@param 	*project		micrograph project.
@param 	filament_width	extracted filament width.
@param 	axis			helical axis alignment: x=1, y=2, z=3.
@param 	&base			file name base.
@param 	&path			path to filament file.
@param 	&ext			filament file extension.
@param 	split			flag to split filaments into individual files.
@return long				number of filaments, <0 on error.
**/
long		project_extract_filaments(Bproject* project, int filament_width, int axis, 
				Bstring& base, Bstring& path, Bstring& ext, int split)
{
	long				i(1), nfil(0);
	Bstring				filename;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bimage*				p;
	
	if ( path.length() && access(path.c_str(), F_OK) ) {
		cerr << "Error: The directory " << path << " must exist to extract filaments!" << endl;
		return -1;
	}
	
	if ( !path.empty() ) if ( path[-1] != '/' ) path += "/";
	
	if ( ext.length() < 1 ) ext = "mrc";
	
	long				nmg = project_count_micrographs(project);
	long				nrec = project_count_reconstructions(project);

	if ( project->select < 1 ) {
		if ( verbose )
			cout << "Extracting filaments from micrographs" << endl << endl;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next, i++ ) {
				nfil = filament_count(mg->fil);
				if ( nfil > 0 ) {
					if ( base.length() ) {
						if ( nmg == 1 )
							filename = base;
						else
							filename = base + Bstring(i, "_%03d");
					} else {
						if ( mg->ffil.length() < 1 ) mg->ffil = mg->fmg;
						filename = mg->ffil.pre_rev('.');
						if ( mg->fmg == mg->ffil )
							filename = filename + Bstring(mg->img_num, "_fil%03d");
					}
					filename = filename + "." + ext;
					if ( filename.length() < 1 )
						cerr << "Error in project_extract_filaments: No filament file name!" << endl;
					if ( split ) {
						filament_setup_filenames(mg->fil, filename, path);
						mg->ffil = 0;
					} else {
						if ( path.length() ) {
							if ( filename.contains("/") ) filename = filename.post_rev('/');
							filename = path + "/" + filename;
						}
						mg->ffil = filename;
					}
					p = read_img(mg->fmg, 1, mg->img_num);
					if ( !p ) return -1;
					micrograph_extract_filaments(mg, p, filament_width, axis);
					delete p;
				}
			}
		}
	} else {
		if ( verbose )
			cout << "Extracting filaments from reconstructions" << endl << endl;
		for ( rec = project->rec; rec; rec = rec->next, i++ ) {
			nfil = filament_count(rec->fil);
			if ( nfil > 0 ) {
				if ( base.length() ) {
					if ( nrec == 1 )
						filename = base;
					else
						filename = base + Bstring(i, "_%03d");
				} else {
					if ( rec->ffil.length() < 1 ) rec->ffil = rec->frec;
					filename = rec->ffil.pre_rev('.');
					if ( rec->frec == rec->ffil )
						filename = filename + "_fil.";
				}
				filename = filename + "." + ext;
				if ( filename.length() < 1 )
					cerr << "Error in project_extract_filaments: No filament file name!" << endl;
				if ( split ) {
					filament_setup_filenames(rec->fil, filename, path);
					rec->ffil = 0;
				} else {
					if ( path.length() ) {
						if ( filename.contains("/") ) filename = filename.post_rev('/');
						filename = path + "/" + filename;
					}
					rec->ffil = filename;
				}
				p = read_img(rec->frec, 1, 0);
				if ( !p ) return -1;
				reconstruction_extract_filaments(rec, p, filament_width, axis);
				delete p;
			}
		}
	}
	
	return nfil;
}

/**
@brief 	Extracts filament images from a micrograph.
@param 	*mg		micrograph parameters.
@param 	*p			micrograph image.
@param 	width			width of box around filament to extract.
@param 	axis				helical axis alignment: x=1, y=2, z=3.
@return long				number of filaments.
**/
long		micrograph_extract_filaments(Bmicrograph* mg, Bimage* p, double width, int axis)
{
	if ( width < 1 ) width = mg->filament_width;

	long			nfil = count_list((char *) mg->fil);

	if ( mg->pixel_size[0] > 0 ) p->sampling(mg->pixel_size);
	
	Bimage*			pfil = filament_extract(mg->fil, p, width, axis);
	
	if ( pfil ) {
		write_img(mg->ffil, pfil, 0);
		delete pfil;
	}
	
	return nfil;
}

/**
@brief 	Extracts filament images from a micrograph.
@param 	*rec	reconstruction parameters.
@param 	*p				reconstruction image.
@param 	width				width of box around filament to extract.
@param 	axis				helical axis alignment: x=1, y=2, z=3.
@return long					number of filaments.
**/
long		reconstruction_extract_filaments(Breconstruction* rec, Bimage* p, double width, int axis)
{
	if ( width < 1 ) width = rec->filament_width;
	
	long			nfil = count_list((char *) rec->fil);

	if ( rec->voxel_size[0] > 0 ) p->sampling(rec->voxel_size);
	
	Bimage*			pfil = filament_extract(rec->fil, p, width, axis);
	
	if ( pfil ) {
		write_img(rec->ffil, pfil, 0);
		delete pfil;
	}
	
	return nfil;
}

/**
@brief 	Extracts filament images from a micrograph.
@param 	*filaments	filament parameters.
@param 	*p			image.
@param 	width		width of box around filament to extract.
@param 	axis		helical axis alignment: x=1, y=2, z=3.
@return Bimage*		multi-filament image.

	If the file names for individual filaments are specified, they are
	used to write individual filaments to files.
	Default axis alignment (for axis==0):
		2D images:	x-axis (1)
		3D images:	z-axis (3)

**/
Bimage*		filament_extract(Bfilament* filaments, Bimage* p, double width, int axis)
{
	if ( axis < 1 || axis > 3 ) {	// Default axis selection
		if ( p->sizeZ() > 1 ) axis = 3;
		else axis = 1;
	}
	
	int						split(0);
	if ( filaments->ffil.length() ) split = 1;
	
	long			n;
	long					nfil(1), nspline;
	double					length, maxlen(0);
	Bfilament*				fil;
	Vector3<double>*		spline;
	Vector3<long>			size((int)width, (int)width, (int)width);
	Bimage*					pfil = NULL;
	Bimage*					pone = NULL;

	if ( p->sizeX() < width ) size[0] = p->sizeX();
	if ( p->sizeY() < width ) size[1] = p->sizeY();
	if ( p->sizeZ() < width ) size[2] = p->sizeZ();
	
	if ( !split ) {
		for ( nfil=0, fil=filaments; fil; fil=fil->next ) {
			length = filament_length(fil) + 1;
			if ( maxlen < length ) maxlen = length;
			nfil++;
		}
		size[axis-1] = (int) maxlen;
		pfil = new Bimage(Float, p->compound_type(), size[0], size[1], size[2], nfil);
		pfil->sampling(p->sampling(0));
		pfil->label(p->label());
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG filament_extract: size=" << size << endl;
		cout << "DEBUG filament_extract: axis=" << axis << " split=" << split << endl;
	}

	for ( n=0, fil=filaments; fil; fil=fil->next, n++ ) {
		spline = vector3_spline_from_nodes(fil->node, nspline);
//		pone = img_extract_filament(p, 0, width, axis, nspline, spline);
		pone = p->extract_filament(0, width, axis, nspline, spline);
		pone->calculate_background();
		if ( split ) {
			if ( fil->ffil.length() ) write_img(fil->ffil, pone, 0);
			else cerr << "Error: No filename for filament " << fil->id << endl;
		} else {
			pfil->replace(n, pone, pone->background(long(0)));
		}
		delete pone;
		delete [] spline;
	}
	
	return pfil;
}

/**
@brief 	Converts filaments to sets of particles.
@param 	*project		project parameters.
@param 	box_size		size of particle box.
@param 	boxing_interval	step size between boxes.
@param 	rise			rise per asymmetric unit in angstrom.
@param 	angle			angular rotation per asymmetric unit in radians. 
@return int				number of filaments converted.

	Particle coordinates are calculated along a spline curve through
	the filament nodes, separated by half the given width.
	Different filaments in a micrograph are indicated by the selection number.

**/
int			project_filaments_to_particles(Bproject* project, Vector3<long> box_size,
				double boxing_interval, double rise, double angle)
{
	int					nfil(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( project->select < 1 ) {
		for ( field=project->field; field; field=field->next )
			for ( mg=field->mg; mg; mg=mg->next )
				nfil += micrograph_filaments_to_particles(mg, box_size, boxing_interval, rise, angle);
	} else {
		for ( rec=project->rec; rec; rec=rec->next )
			nfil += reconstruction_filaments_to_particles(rec, box_size, boxing_interval, rise, angle);
	}
	
	return nfil;
}

/**
@brief 	Converts filaments to sets of particles.
@param 	*mg				micrograph parameters.
@param 	box_size		size of particle box.
@param 	boxing_interval	step size between boxes.
@param 	rise			rise per asymmetric unit in angstrom.
@param 	angle			angular rotation per asymmetric unit in radians. 
@return int				number of filaments converted.

	Particle coordinates are calculated along a spline curve through
	the filament nodes, separated by half the given width.
	Different filaments in a micrograph are indicated by the selection number.

**/
int			micrograph_filaments_to_particles(Bmicrograph* mg, Vector3<long> box_size,
				double boxing_interval, double rise, double angle)
{
	if ( box_size[0] < 1 ) box_size[0] = 40;
	if ( box_size[1] < 1 ) box_size[1] = box_size[0];
	box_size[2] = 1;
	
	mg->box_size = box_size;
	if ( boxing_interval < 1 ) boxing_interval = box_size[0]/2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_filaments_to_particles: box_size=" << mg->box_size << endl;
	
	particle_kill(mg->part);
	
	mg->part = filaments_to_particles(mg->fil, mg->pixel_size, mg->box_size, boxing_interval, rise, angle);

	mg_part_links(mg);
	
	return filament_count(mg->fil);
}
	
/**
@brief 	Converts filaments to sets of particles.
@param 	*rec			reconstruction parameters.
@param 	box_size		size of particle box.
@param 	boxing_interval	step size between boxes.
@param 	rise			rise per asymmetric unit in angstrom.
@param 	angle			angular rotation per asymmetric unit in radians. 
@return int				number of filaments converted.

	Particle coordinates are calculated along a spline curve through
	the filament nodes, separated by half the given width.
	Different filaments in a micrograph are indicated by the selection number.

**/
int			reconstruction_filaments_to_particles(Breconstruction* rec, Vector3<long> box_size,
				double boxing_interval, double rise, double angle)
{
	if ( box_size[0] < 1 ) box_size[0] = 40;
	if ( box_size[1] < 1 ) box_size[1] = box_size[0];
	if ( box_size[2] < 1 ) box_size[2] = box_size[0];
	
	rec->box_size = box_size;
	if ( boxing_interval < 1 ) boxing_interval = box_size[0]/2;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_filaments_to_particles: box_size=" << rec->box_size << endl;
	
	particle_kill(rec->part);
	
	rec->part = filaments_to_particles(rec->fil, rec->voxel_size, rec->box_size, boxing_interval, rise, angle);

	rec_part_links(rec);
	
	return filament_count(rec->fil);
}
	
/**
@brief 	Converts filaments to sets of particles.
@param 	*filaments		filament parameters.
@param 	pixel_size		sampling.
@param 	box_size		size of particle box.
@param 	boxing_interval	step size between boxes.
@param 	rise			rise per asymmetric unit in angstrom.
@param 	angle			angular rotation per asymmetric unit in radians.
@return Bparticle*		pointer to new list of particles.

	Particle coordinates are calculated along a spline curve through
	the filament nodes, separated by half the given width.
	Different filaments in a micrograph are indicated by the selection number.

**/
Bparticle*	filaments_to_particles(Bfilament* filaments,
				Vector3<double> pixel_size, Vector3<long> box_size,
				double boxing_interval, double rise, double angle)
{
	long				i, n(1), nspline, it, nfil(0);
	Bparticle*			part_list = NULL;
	Bparticle*			part = NULL;
	Bfilament*			fil;
	Vector3<double>*	spline;
	Vector3<double>		dir;
	
	double				bi, dr(rise/pixel_size[0]), phi(0), st, t, od(0);

	if ( verbose & VERB_FULL ) {
		cout << "Converting filaments to particles:" << endl;
		cout << "Box size:                       " << box_size << endl;
		cout << "Boxing interval:                " << boxing_interval << endl;
		cout << "Helical rise and angle:         " << rise << " A " << angle << " degrees" << endl << endl;
	}
	
	for ( fil=filaments; fil; fil=fil->next ) {
		nfil++;
		spline = vector3_spline_from_nodes(fil->node, nspline);
		for ( bi=0; bi<nspline; bi += boxing_interval, n++ ) {
			i = (int) (bi + 0.5);
			part = particle_add(&part, n);
			if ( !part_list ) part_list = part;
			part->group = fil->id;
			part->loc = spline[i];
			part->ori = box_size/2;
			if ( i == 0 ) dir = spline[1] - spline[0];
			else dir = spline[i] - spline[i-1];
			st = sqrt(1 - dir[2]*dir[2]);
			if ( angle ) {
				t = i/dr;
				it = (int) (t + 0.5);
				phi = angle * it;
				od = (t - it)*dr;
				part->ori -= dir*od;
			}
			part->view = View(cos(phi)*st, sin(phi)*st, dir[2], angle_set_negPI_to_PI(phi - atan2(dir[1],dir[0])));
//			if ( spline[i][2] < 1 ) part->ori[2] = 0;
			part->sel = nfil;
		}
		delete[] spline;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG filaments_to_particles: npart=" << n << endl;
	
	return part_list;
}

/**
@brief 	Masks particles extracted from filaments.
@param 	*project		project parameters.
@param 	mask_width		width of mask to apply.
@return int				number of particles masked.

	The orientation of the filament is inferred from the view angle and
	the adjacent regions are set to the background.

**/
int			project_mask_filament_particles(Bproject* project, int mask_width)
{
	int					npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;

	if ( project->select < 1 ) {
		for ( field=project->field; field; field=field->next )
			for ( mg=field->mg; mg; mg=mg->next )
				npart += mask_filament_particles(mg->fpart, mg->part, mask_width);
	} else {
		for ( rec=project->rec; rec; rec=rec->next )
			npart += mask_filament_particles(rec->fpart, rec->part, mask_width);
	}
	
	return npart;
}

/**
@brief 	Rotates and masks particles extracted from filaments.
@param 	*project		project parameters.
@param 	rotation_axis	axis to rotate to: 1=x, 2=y, 3=z.
@param 	back_flag		background correction flag.
@param 	mask_width		width of mask to apply.
@return int				number of particles masked.

	The orientation of the filament is inferred from the view angle and
	the particle is rotated to orient the filament axis along a cartesian axis. 
	The mask is applied and adjacent regions are set to the background.

**/
int			project_rotate_mask_filament_particles(Bproject* project, 
				int rotation_axis, int back_flag, int mask_width)
{
	if ( project->select ) return 0;	// No 3D yet!
	
	if ( rotation_axis < 1 || rotation_axis > 3 ) return 0;

	if ( project->select < 1 && rotation_axis > 2 ) return 0;
	
	int					npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose ) {
		cout << "Rotating and masking filament particles:" << endl;
		cout << "Rotation axis:                  " << (char) (119+rotation_axis) << endl;
		cout << "Mask width:                     " << mask_width << endl;
		if ( back_flag ) cout << "With background correction" << endl;
		cout << endl;
	}

	if ( project->select < 1 ) {
		for ( field=project->field; field; field=field->next )
			for ( mg=field->mg; mg; mg=mg->next )
				npart += rotate_mask_filament_particles(mg->fpart, mg->part, 
							rotation_axis, back_flag, mask_width);
	} else {
		for ( rec=project->rec; rec; rec=rec->next )
			npart += rotate_mask_filament_particles(rec->fpart, rec->part, 
						rotation_axis, back_flag, mask_width);
	}
	
	return npart;
}

/*
@brief 	Masks particles extracted from filaments.
@param 	&filename		particle file name.
@param 	*partlist		particle parameters.
@param 	mask_width		width of mask to apply.
@return int				number of particles masked.

	The orientation of the filament is inferred from the view angle and
	the adjacent regions are set to the background.

**/
int			mask_filament_particles(Bstring& filename, Bparticle* partlist, int mask_width)
{
	int					npart(0);
	long		i, n, x, y, z;
	double				d, halfwidth = mask_width/2.0;
	Vector3<double>		u, v;
	Euler				euler;
	Bparticle*			part;
	
	Bimage*				p = read_img(filename, 1, -1);

	long		imagesize = p->sizeX()*p->sizeY()*p->sizeZ()*p->channels();
	
	for ( n=0, part = partlist; part; part = part->next, n++ ) {
		euler = Euler(part->view);
		u[0] = cos(euler.psi());
		u[1] = -sin(euler.psi());
		for ( i=n*imagesize, z=0; z<p->sizeZ(); z++ ) {
			v[2] = z - p->image[n].origin()[2];
			for ( y=0; y<p->sizeY(); y++ ) {
				v[1] = y - p->image[n].origin()[1];
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					v[0] = x - p->image[n].origin()[0];
					d = (u.cross(v)).length();
					if ( d > halfwidth ) p->set(i, p->background(n));
				}
			}
		}
		npart++;
	}
	
	write_img(filename, p, 0);
	
	delete p;
	
	return npart;
}

/*
@brief	Rotates and masks particles extracted from filaments.
@param 	&filename		particle file name.
@param 	*partlist		particle parameters.
@param 	rotation_axis	axis to rotate to: 1=x, 2=y, 3=z.
@param 	back_flag		background correction flag.
@param 	mask_width		width of mask to apply.
@return	int				number of particles masked.

	The orientation of the filament is inferred from the view angle and
	the particle is rotated to orient the filament axis along a cartesian axis. 
	The mask is applied and adjacent regions are set to the background.
**/
int			rotate_mask_filament_particles(Bstring& filename, Bparticle* partlist, 
				int rotation_axis, int back_flag, int mask_width)
{
	if ( !partlist ) return 0;
	
	int					npart(0);
	long				i, j, n, x, y, z;
	double				d(0), dx(0), dy(0), dz(0), halfwidth;
	Euler				euler;
	Bparticle*			part;
	Bimage*				pex;
	Bimage*				p = read_img(filename, 1, -1);
	
	if ( !p ) {
		error_show("Error in rotate_mask_filament_particles: No image!", __FILE__, __LINE__);
		return -1;
	}
	
	p->change_type(Float);
	
	if ( mask_width < 1 ) mask_width = (p->sizeX() > p->sizeY())? p->sizeX(): p->sizeY();
	halfwidth = mask_width/2.0;
	
	long		imagesize = p->sizeX()*p->sizeY()*p->sizeZ()*p->channels();
	
	for ( n=0, part = partlist; part; part = part->next, n++ ) {
		pex = p->extract(n);
		pex->image[0] = p->image[n];
		euler = Euler(part->view);
		euler[1] = euler[2] = 0;
		if ( rotation_axis == 2 ) euler[0] += M_PI_2;
		pex->rotate(euler.matrix());
//		img_rotate(pex, euler);
//		pex->rotate(euler.view());
		for ( i=n*imagesize, j=0, z=0; z<p->sizeZ(); z++ ) {
			if ( p->sizeZ() > 1 ) {
				dz = fabs(z - p->image[n].origin()[2]);
				dz *= dz;
			}
			for ( y=0; y<p->sizeY(); y++ ) {
				dy = fabs(y - p->image[n].origin()[1]);
				if ( p->sizeZ() > 1 ) dy *= dy;
				for ( x=0; x<p->sizeX(); x++, i++, j++ ) {
					dx = fabs(x - p->image[n].origin()[0]);
					if ( p->sizeZ() < 2 ) {
						if ( rotation_axis == 2 ) d = dx;
						else d = dy;
					} else {
						dx *= dx;
						d = sqrt(dx + dy + dz);
					}
					if ( d > halfwidth ) p->set(i, p->background(n));
					else p->set(i, (*pex)[j]);
				}
			}
		}
		part->view[3] -= euler.psi();
		delete pex;
		npart++;
	}

	if ( back_flag )
		p->correct_background();
	
	write_img(filename, p, 0);
	
	delete p;
	
	return npart;
}

