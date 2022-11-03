/**
@file	rwmg.cpp
@brief	Library routines to read and write micrograph parameters
@author Bernard Heymann
@date	Created: 20030418
@date	Modified: 20220406
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_ctf.h"
#include "mg_tags.h"
#include "rwmgSTAR.h"
#include "rwmgXML.h"
#include "rwmgEMX.h"
#include "rwmgRELION.h"
#include "rwmgIMOD.h"
#include "rwmgSerialEM.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string	command;		// Command line

// Internal prototype
int			project_check(Bproject* project, int flags=0);

int			project_display_counts(Bproject* project)
{
	long			nfield = project_count_fields(project);
	long			nmg = project_count_micrographs(project);
	long			nmgs = project_count_mg_selected(project);
	long			npart = project_count_mg_particles(project);
	long			nsel = project_count_mg_part_selected(project);
	long			ngrp = project_count_mg_groups(project);
	long			ngrpsel = project_count_mg_groups_selected(project);
	long			nfil = project_count_mg_filaments(project);
	long			nfn = project_count_mg_filament_nodes(project);
	
	cout << endl;
	cout << "Fields:                         " << nfield << endl;
	cout << "Micrographs:                    " << nmg << " ( " << nmgs << " selected)" << endl;
	cout << "Particles:                      " << npart << " ( " << nsel << " selected)" << endl;
	cout << "Particle groups:                " << ngrp << " ( " << ngrpsel << " selected)" << endl;
	cout << "Filaments and nodes:            " << nfil << " " << nfn << endl;
	cout << endl;
	
	long			nrec = project_count_reconstructions(project);
	long			nrecs = project_count_rec_selected(project);
	npart = project_count_rec_particles(project);
	nsel = project_count_rec_part_selected(project);
	nfil = project_count_rec_filaments(project);
	nfn = project_count_rec_filament_nodes(project);	
	
	cout << "Reconstructions:                " << nrec << " ( " << nrecs << " selected)" << endl;
	cout << "Particles:                      " << npart << " ( " << nsel << " selected)" << endl;
	cout << "Filaments and nodes:            " << nfil << " " << nfn << endl;
	cout << endl;
	
	return 0;
}

/* The filename is for the current input parameter file */
int			field_resolve_file_access(Bfield* field, Bmicrograph* mg, Bstring filename, int flags)
{
	if ( !field ) return 0;
	if ( filename.length() < 1 ) return 0;
	
	Bstring		path;
	Bparticle*	part;

	if ( filename.contains("/") ) path = filename.pre_rev('/');
	
	if ( verbose )
		cout << "Checking for micrographs in path: " << path << endl;
	
	for ( ; field; field = field->next ) {
		if ( !mg ) mg = field->mg;
		for ( ; mg; mg = mg->next ) {
			mg->fmg = find_file(mg->fmg, path, flags);
			mg->fframe = find_file(mg->fframe, path, flags);
			mg->fpart = find_file(mg->fpart, path, flags);
			mg->ffil = find_file(mg->ffil, path, flags);
			mg->fft = find_file(mg->fft, path, flags);
			mg->fps = find_file(mg->fps, path, flags);
			if ( mg->part && mg->part->fpart.length() )
				for ( part = mg->part; part; part = part->next )
					part->fpart = find_file(part->fpart, path, flags);
		}
	}
	
	return 0;
}

/* The filename is for the current input parameter file */
int			reconstruction_resolve_file_access(Breconstruction* rec, Bstring filename, int flags)
{
	if ( !rec ) return 0;
	if ( filename.length() < 1 ) return 0;
	
	Bstring		path;
	Bparticle*	part;
	
	if ( filename.contains("/") ) path = filename.pre_rev('/');
	
	if ( verbose )
		cout << "Checking for reconstructions in path: " << path << endl;
	
	for ( ; rec; rec = rec->next ) {
		rec->frec = find_file(rec->frec, path, flags);
		rec->fft = find_file(rec->fft, path, flags);
		rec->fps = find_file(rec->fps, path, flags);
		rec->fpart = find_file(rec->fpart, path, flags);
		rec->ffil = find_file(rec->ffil, path, flags);
		for ( part = rec->part; part; part = part->next )
			part->fpart = find_file(part->fpart, path, flags);
	}
	
	return 0;
}

int			project_resolve_file_access(Bproject* project, Bstring filename, int flags)
{
	if ( filename.length() < 1 ) return 0;
	
	field_resolve_file_access(project->field, NULL, filename, flags);

	reconstruction_resolve_file_access(project->rec, filename, flags);

	return 0;
}

/**
@brief 	Reading micrograph parameters from files.
@param 	&filename		file name (or comma-delimited list).
@param	flags			flags to modify behavior.
@return Bproject*			project structure, NULL if reading failed.

	Flags (bits):
	1	check particle images against images.
	8	warn if files not found.
	16	delete file names of files not found.
	32	update micrograph intensities.
	64	update STAR tags.

**/
Bproject*	read_project(const char* filename, int flags)
{
	Bstring		thefile(filename);
	return read_project(thefile, flags);
}

Bproject*	read_project(Bstring& filename, int flags)
{
	Bstring		xsdfile;
	return read_project(filename, xsdfile, flags);
}

Bproject*	read_project(Bstring& filename, Bstring& xsdfile, int flags)
{
	if ( filename.empty() ) {
		error_show("No micrograph parameter filename!", __FILE__, __LINE__);
		return NULL;
	}
	
	Bstring				ext = filename.extension();
	
	if ( verbose & VERB_LABEL )
		cout << "Parameter filename: " << filename << endl;

	int					err(0);
	FileType			type;
	Bproject*			project = new Bproject;
	project->filename = filename;
	
	if ( ext.contains("star") ) {
		type = file_type(filename);
		if ( type == Micrograph )
			err = read_project_star(filename, project, flags&64);
		else if ( type == MgRelion )
			err = read_project_relion(filename, project);
		else
			cerr << "Error: The STAR file is not a valid form Bsoft can read!" << endl;
	} else if ( ext.contains("xml") ) {
		err = read_project_xml(filename, project);
	} else if ( ext.contains("emx") ) {
		err = read_project_emx(filename, project, xsdfile);
	} else if ( ext.contains("imod") ) {
		err = read_project_imod(&filename, project, 0);
	} else if ( ext.contains("mdoc") ) {
		err = read_project_serialem(filename, project, 0);
	} else {
		cerr << "Error: Extension \"" << ext << "\" not valid for parameter files!" << endl;
		cerr << "	Parameter file \"" << filename << "\" not read." << endl;
		err = -1;
	}
	
	if ( err < 0 ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		project_kill(project);
		return NULL;
	}

	field_resolve_file_access(project->field, NULL, filename, flags & 8);

	reconstruction_resolve_file_access(project->rec, filename, flags & 8);
	
	project_check(project, flags);
	
	if ( verbose & VERB_PROCESS )
		project_display_counts(project);
	
	return project;
}

Bproject*	read_project(Bstring* file_list, int flags)
{
	Bstring		xsdfile;
	return read_project(file_list, xsdfile, flags);
}

Bproject*	read_project(Bstring* file_list, Bstring& xsdfile, int flags)
{
	if ( !file_list ) {
		error_show("No micrograph parameter filename!", __FILE__, __LINE__);
		return NULL;
	}
	
	Bstring				ext = file_list->extension();
	Bstring*			thisfile;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Reading parameter files:";
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	int					err(0);
	FileType			type;
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	project->filename = *file_list;
	
	if ( ext.contains("imod") ) {
		err = read_project_imod(file_list, project, 0);
	} else if ( ext.contains("mdoc") ) {
		err = read_project_serialem(*file_list, project, 0);
	} else for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		for ( field = project->field; field && field->next; field = field->next ) ;
		if ( field ) for ( mg = field->mg; mg && mg->next; mg = mg->next ) ;
		for ( rec = project->rec; rec && rec->next; rec = rec->next ) ;
		ext = thisfile->extension();
		if ( ext.contains("star") ) {
			type = file_type(*thisfile);
			if ( type == Micrograph )
				err = read_project_star(*thisfile, project, flags&64);
			else if ( type == MgRelion )
				err = read_project_relion(*thisfile, project);
			else
				cerr << "Error: The STAR file is not a valid form Bsoft can read!" << endl;
		} else if ( ext.contains("xml") ) {
			err = read_project_xml(*thisfile, project);
		} else if ( ext.contains("emx") ) {
			err = read_project_emx(*thisfile, project, xsdfile);
		} else {
			cerr << "Error: Extension \"" << ext << "\" not valid for parameter files!" << endl;
			cerr << "	Parameter file \"" << *thisfile << "\" not read." << endl;
			err = -1;
		}
		if ( !err ) {
			if ( !field ) field = project->field;
//			else field = field->next;
//			cout << "field " << field->id << endl;
			if ( field ) {
				if ( !mg ) mg = field->mg;
//				else mg = mg->next;
//				cout << "mg " << mg->id << endl;
				if ( mg ) field_resolve_file_access(field, mg, *thisfile, flags & 8);
			}
			if ( !rec ) rec = project->rec;
//			else rec = rec->next;
			if ( rec ) reconstruction_resolve_file_access(rec, *thisfile, flags & 8);
		}
	}
	
	if ( err < 0 ) {
		error_show(file_list->c_str(), __FILE__, __LINE__);
		project_kill(project);
		return NULL;
	}
	
	project_check(project, flags);
	
	if ( verbose & VERB_PROCESS )
		project_display_counts(project);
	
	return project;
}

long		append_project(Bproject* project, Bstring* file_list, Bstring& xsdfile, int flags)
{
	long			nmg(0);
	
	Bproject*		project_append = read_project(file_list, xsdfile, flags);
	
	Bfield*			field;
	Bmicrograph*	mg;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, nmg++ ) ;
		if ( project_append && !field->next ) {
			field->next = project_append->field;
			project_append->field = NULL;
			project_kill(project_append);
			project_append = NULL;
		}
	}
	
	return nmg;
}

Bparticle*	read_particle(Bstring& filename, Bparticle** partlist, FOMType fom_tag[NFOM])
{
#if defined(HAVE_XML)
	return read_particle_xml(filename, partlist, fom_tag);
#else
	return NULL;
#endif
}


/**
@brief 	Writing micrograph parameters to a file.
@param	*filename		file name (or comma-delimited list).
@param 	*project		project structure.
@return int				error code (<0 means failure).
**/
int			write_project(const char* filename, Bproject* project)
{
	Bstring		thefile(filename);
	return write_project(thefile, project, 0);
}

/**
@brief 	Writing micrograph parameters to a file.
@param	*filename		file name (or comma-delimited list).
@param 	*project		project structure.
@param 	mg_select		flag to only write selected micrographs.
@param 	rec_select		flag to only write selected reconstructions.
@return int				error code (<0 means failure).

**/
int			write_project(Bstring& filename, Bproject* project, 
				int mg_select, int rec_select)
{
	int		flags(0);
	if ( mg_select ) flags |= 2;
	if ( rec_select ) flags |= 4;
	return write_project(filename, project, flags);
}

/**
@brief 	Writing micrograph parameters to a file.
@param	&filename		file name (or comma-delimited list).
@param 	*project		project structure.
@param	flags			flags to modify behavior.
@return int				error code (<0 means failure).

	Flags (bits):
	1.	check particle images.
	2.	write only selected micrographs.
	4.	write only selected reconstructions.
	8.	warn if files not found.
	16.	delete file names of files not found.
	32.	update micrograph intensities.
	
**/
int			write_project(Bstring& filename, Bproject* project, int flags)
{
	int				err(0);
	int 			mg_select(flags & 2), rec_select(flags & 4);
	
	project_resolve_file_access(project, filename, flags);
	
	project_check(project, flags);
	
	project->comment += command_line_time();

	if ( verbose )
		project_display_counts(project);
	
	if ( project->split == 9 )
		return write_project_star(filename, project, mg_select, rec_select);
	
	if ( filename.empty() )
		return error_show("No micrograph parameter filename!", __FILE__, __LINE__);
	
	Bstring			ext = filename.extension();

	if ( verbose )
		cout << "Writing parameter file: " << filename << endl;

    if ( ext.contains("star") ) {
		err = write_project_star(filename, project, mg_select, rec_select);
	} else if ( ext.contains("xml") ) {
		err = write_project_xml(filename, project, mg_select, rec_select);
	} else if ( ext.contains("emx") ) {
		err = write_project_emx(filename, project, mg_select, rec_select);
	} else if ( ext.contains("imod") ) {
		err = write_project_imod(filename, project);
	} else {
		cerr << "Error: Extension \"" << ext << "\" not valid for parameter files!" << endl;
		cerr << "	Parameter file \"" << filename << "\" not written." << endl;
		err = -1;
	}
	
	if ( err < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	else
		project->filename = filename;
	
	return err;
}

/**
@brief 	Writing particle coordinates to a file.
@param	&filename		file name (or comma-delimited list).
@param 	*project		project structure.
@param 	flags			bits: 1=in angstrom
@return int				error code (<0 means failure).

**/
long		write_particle_list(Bstring& filename, Bproject* project, int flags)
{
	ofstream		ftxt(filename.c_str());
	if ( ftxt.fail() ) {
		cerr << "Error: Not able to write " << filename << endl;
		return -1;
	}

	long				i(0), npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Bstring				fn;
	Vector3<double>		loc;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) if ( field->select > 0 ) {
			for ( mg = field->mg; mg; mg = mg->next, ++i ) if ( mg->select > 0 ) {
				for ( part = mg->part; part; part = part->next ) if ( part->sel > 0 ) {
					if ( mg->fmg.length() ) fn = mg->fmg;
					else if ( mg->fpart.length() ) fn = mg->fpart;
					else if ( part->fpart.length() ) fn = part->fpart;
					else fn = Bstring(i, "%d");
					if ( fn.contains("/") ) fn = fn.post_rev('/');
					if ( flags & 1 ) loc = part->loc * mg->pixel_size;
					else loc = part->loc;
					ftxt << fn << " " << loc[0] << " " << loc[1] << endl;
					npart++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next, ++i ) if ( rec->select > 0 ) {
			for ( part = rec->part; part; part = part->next ) if ( part->sel > 0 ) {
				if ( rec->frec.length() ) fn = rec->frec;
				else if ( rec->fpart.length() ) fn = rec->fpart;
				else if ( part->fpart.length() ) fn = part->fpart;
				else fn = Bstring(i, "%d");
				if ( fn.contains("/") ) fn = fn.post_rev('/');
				if ( flags & 1 ) loc = part->loc * rec->voxel_size;
				else loc = part->loc;
				ftxt << fn << " " << loc[0] << " " << loc[1] << endl;
				npart++;
			}
		}
	}
	
	ftxt.close();
	
	if ( verbose )
		cout << npart << " particles exported to " << filename << endl;
	
	return npart;
}


int			write_particle(Bstring& filename, Bparticle* part, int euler_flag, int omega_flag, FOMType fom_tag[NFOM])
{
#if defined(HAVE_XML)
	return write_particle_xml(filename, part, euler_flag, omega_flag, fom_tag);
#else
	return 0;
#endif
}

Bstring		ppx_filename(Bstring& id, int part_id)
{
	Bstring		filename = "ppx/" + id + "_" + Bstring(part_id, "%04d") + ".ppx";
	return filename;
}

int			ppx_exists(Bparticle* part, int flag)
{
	int			exist(0);
	Bstring		filename;

	if ( part->mg )
		filename = ppx_filename(part->mg->id, part->id);
	else if ( part->rec )
		filename = ppx_filename(part->rec->id, part->id);

	if ( access(filename.c_str(), R_OK) == 0 )
		exist = 1;
	
	// flag=1 for absent, flag=2 for present
	if ( flag == exist + 1 )
		cout << part->sel << tab << filename << endl;
	
	return exist;
}

int			ppx_check(Bparticle* part, FOMType fom_tag[NFOM])
{
	Bstring		filename;

	if ( part->mg )
		filename = ppx_filename(part->mg->id, part->id);
	else if ( part->rec )
		filename = ppx_filename(part->rec->id, part->id);

	if ( access(filename.c_str(), R_OK) )
		return 0;

	Bparticle*	part_list = NULL;

	Bparticle*	part_nu = read_particle(filename, &part_list, fom_tag);
	
	if ( part_nu ) {
		particle_update(part, part_nu);
		particle_kill(part_nu);
		return 1;
	}
	
	return 0;
}

int			ppx_check(Bparticle* part)
{
	FOMType 	fom_tag[NFOM] = {FOM, FOM_CV};
	
	return ppx_check(part, fom_tag);
}

int			project_check(Bproject* project, int flags)
{
	long				i, j, num_img_files(0), nwarn(0);
//	Vector3<double>		d, ori;
	double				d;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Bbadarea*			bad;
	Bstring				part_file;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_check: " << endl;
	 
	if ( project->field == NULL && project->rec ) project->select = 1;
	
	for ( field=project->field; field; field=field->next ) {
		i = 1;
		for ( mg=field->mg; mg; mg=mg->next ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_check: mg=" << mg->id << endl;
			if ( mg->next ) {
				if ( mg->id == mg->next->id ) {
					mg->id = mg->id + Bstring(i++, "_%03d");
					if ( verbose )
						cout << "Resetting micrograph ID: " << mg->id << endl;
					if ( !mg->next->next ) {
						mg->next->id = mg->next->id + Bstring(i, "_%03d");
						if ( verbose )
							cout << "Resetting micrograph ID: " << mg->next->id << endl;
					}
				}
			}
			if ( mg->fmg.length() ) num_img_files++;
			if ( mg->fframe.length() ) num_img_files++;
			if ( mg->fpart.length() ) num_img_files++;
			if ( mg->ffil.length() ) num_img_files++;
			if ( mg->fps.length() ) num_img_files++;
			if ( num_img_files < 1 )
				cerr << "Warning: No image files referenced for micrograph " << mg->id << endl;
			if ( mg->sampling < 1e3 ) mg->sampling *= 1e4;
			if ( mg->magnification < 1e3 ) mg->magnification *= 1e3;
			mg->frame_pixel_size[2] = 1;
			if ( mg->frame_pixel_size[0] < 0.01 || mg->frame_pixel_size[1] < 0.01 ) {
				mg->frame_pixel_size = Vector3<double>(1,1,1);
				if ( mg->magnification && mg->sampling ) {
					mg->frame_pixel_size[0] = mg->sampling/mg->magnification;
					if ( mg->frame_pixel_size[0] < 0.01 ) mg->frame_pixel_size[0] *= 1e4; // Assume um
					mg->frame_pixel_size[1] = mg->frame_pixel_size[0];
				}
				if ( verbose && nwarn < 10 ) {
					cerr << "Warning: Resetting frame pixel size for micrograph " <<
						mg->id << " to " << mg->frame_pixel_size << endl;
					nwarn++;
				} else if ( nwarn == 10 ) {
					cerr << "Warning: ..." << endl;
					nwarn++;
				}
			}
			mg->pixel_size[2] = 1;
			if ( mg->pixel_size[0] < 0.01 || mg->pixel_size[1] < 0.01 ) {
				if ( mg->frame_pixel_size.volume() ) {
					mg->pixel_size = mg->frame_pixel_size;
				} else {
					mg->pixel_size = Vector3<double>(1,1,1);
					if ( mg->magnification && mg->sampling ) {
						mg->pixel_size[0] = mg->sampling/mg->magnification;
						if ( mg->pixel_size[0] < 0.01 ) mg->pixel_size[0] *= 1e4; // Assume um
						mg->pixel_size[1] = mg->pixel_size[0];
					}
					cerr << "Warning: Resetting pixel size for micrograph " <<
						mg->id << " to " << mg->pixel_size << endl;
				}
			}
			if ( flags & 32 && mg->intensity <= 0 ) micrograph_intensity(mg, flags);
			if ( mg->scale[0] <= 0 ) mg->scale[0] = 1;
			if ( mg->scale[1] <= 0 ) mg->scale[1] = 1;
			if ( mg->scale[2] <= 0 ) mg->scale[2] = 1;
			if ( mg->ctf ) {
				if ( mg->ctf->volt() < 1e3 ) mg->ctf->volt(1e3*mg->ctf->volt());
				if ( mg->ctf->volt() < 1e3 )
					cerr << "Warning: The acceleration voltage must be specified!" << endl;
				mg->ctf->zero(1);
			}
			if ( mg->part ) part_file = mg->part->fpart;
			for ( part = mg->part; part; part = part->next ) {
				if ( part->id < 1 )
					cerr << "Error: Particle has an ID less than 1! (" << mg->id << ")" << endl;
				part->mg = mg;
				if ( part_file.length() && part_file != part->fpart ) part_file = 0;
				if ( part->pixel_size[0] < 0.01 ) part->pixel_size = mg->pixel_size;
				if ( part->pixel_size[1] < 0.01 ) part->pixel_size[1] = part->pixel_size[0];
				part->pixel_size[2] = 1;
				d = (part->ori - mg->box_size/2).length();
				if ( d > mg->box_size[0]/4 )
					part->ori = mg->box_size/2;
				part->ori[2] = 0;
			}
			if ( part_file.length() ) {
				mg->fpart = part_file;
				for ( part = mg->part; part; part = part->next ) part->fpart = 0;
				part_file = 0;
			}
			for ( j=0, bad = mg->bad; bad; bad = bad->next ) bad->id = --j;
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_check: micrographs=" << project_count_micrographs(project) << endl;
	 
	i = 1;
	for ( rec=project->rec; rec; rec=rec->next ) {
		if ( rec->next ) {
			if ( rec->id == rec->next->id ) {
				rec->id = rec->id + Bstring(i++, "_%03d");
				if ( verbose )
					cout << "Resetting reconstruction ID: " << rec->id << endl;
				if ( !rec->next->next ) {
					rec->next->id = rec->next->id + Bstring(i, "_%03d");
					if ( verbose )
						cout << "Resetting reconstruction ID: " << rec->next->id << endl;
				}
			}
		}
		if ( rec->voxel_size[0] < 0.01 )
			rec->voxel_size = Vector3<double>(1,1,1);
		if ( rec->scale[0] <= 0 ) rec->scale[0] = 1;
		if ( rec->scale[1] <= 0 ) rec->scale[1] = 1;
		if ( rec->scale[2] <= 0 ) rec->scale[2] = 1;
		for ( part = rec->part; part; part = part->next ) {
			if ( part->id < 1 )
				cerr << "Error: Particle has an ID less than 1! (" << mg->id << ")" << endl;
			part->rec = rec;
			if ( part->pixel_size[0] < 0.01 )
				part->pixel_size = rec->voxel_size;
			d = (part->ori - rec->box_size/2).length();
			if ( d > rec->box_size[0]/4 )
				part->ori = rec->box_size/2;
		}
		for ( j=0, bad = rec->bad; bad; bad = bad->next ) bad->id = --j;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_check: reconstructions=" << project_count_reconstructions(project) << endl;

	field = project->field;
	if ( ( flags & 1 ) && field ) {
		if ( project_check_particles(project) < 0 ) bexit(-1);
	}

	project_update_first_zero(project);
		
	return 0;
}

/**
@brief 	Finds a particle associated with a file name.
@param 	*project 	micrograph processing parameter structure.
@param 	fn			filename to search for.
@return	Bparticle*		particle structure if found, otherwise NULL.

	The file name to search for is altered to represent an accurate path.

**/
Bparticle*	project_find_particle(Bproject* project, Bstring& fn)
{
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	Bstring				path;
	if ( project->filename.contains("/") )
		path = project->filename.pre_rev('/');

	fn = find_file(fn, path);
//	if ( verbose & VERB_FULL )
		cout << "Searching for file " << fn << endl;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				part = mg->part;
				if ( mg->fpart == fn ) return part;
				for ( ; part; part = part->next )
					if ( part->fpart == fn ) return part;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			part = rec->part;
			if ( rec->fpart == fn ) return part;
			for ( ; part; part = part->next )
				if ( part->fpart == fn ) return part;
		}
	}
	
	if ( verbose )
		cerr << "Error: File " << fn << " not found!" << endl;
	
	return NULL;
}

/**
@brief 	Updates a project with ppx files.
@param 	*project 	project parameter structure.
@return	long			number of particles updated.

**/
long		project_update_from_ppx(Bproject* project)
{
	long				npart(0);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
		
	FOMType 	fom_tag[NFOM] = {FOM, FOM_CV};

	if ( verbose ) {
		cout << "Updating particles from ppx files" << endl;
		cout << "FOM tags:";
		for ( int f=0; f<NFOM && fom_tag[f]; ++f )
			cout << tab << fom_tag[f];
		cout << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					npart += ppx_check(part, fom_tag);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				npart += ppx_check(part, fom_tag);
	}
	
	if ( verbose )
		cout << "Number of particles updated: " << npart << endl << endl;

	return 0;
}

/**
@brief 	Lists existing ppx files for a project.
@param 	*project 	project parameter structure.
@param	flag		select: 0 = all, 1 = existing, 2 = absent.
@return	long			number of particles updated.
**/
long		project_list_ppx(Bproject* project, int flag)
{
	long				npart(0), nexist(0);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
		
	if ( verbose ) {
		if ( flag == 2 )
			cout << "Existing ppx files" << endl;
		else if ( flag == 1 )
			cout << "Absent ppx files" << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next, npart++ )
					nexist += ppx_exists(part, flag);
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next, npart++ )
				nexist += ppx_exists(part, flag);
	}
	
	if ( verbose ) {
		cout << "Number of particles:        " << npart << endl;
		cout << "Number of ppx files:        " << nexist << endl;
		cout << "Number of absent ppx files: " << npart - nexist << endl << endl;
	}
	
	return 0;
}


/**
@brief 	Splits a project into sets of particle selections and write the files.
@param 	&filename	base filename modified to include the particle selection number.
@param 	*project 	the project.
@return int 				number of sets generated.

	The relevant micrographs are selected for each particle selection number
	and written into a file.

**/
int			project_split_write(Bstring& filename, Bproject* project)
{
	long				n(1), c, max(0);
	Bstring				fn;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( max < part->sel ) max = part->sel;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( max < part->sel ) max = part->sel;
	}
	
	if ( verbose )
		cout << "Splitting into " << max << " sets of particles" << endl;

	for ( n=1; n<=max; n++ ) {
		if ( project->select < 1 ) {
			for ( c=0, field = project->field; field; field = field->next ) {
				for ( mg = field->mg; mg; mg = mg->next ) {
					mg->select = 0;
					for ( part = mg->part; part; part = part->next ) if ( part->sel == n ) {
						mg->select = 1;
						c++;
					}
				}
			}
			for ( rec = project->rec; rec; rec = rec->next )
				rec->select = 0;
		} else {
			for ( field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					mg->select = 0;
			for ( c=0, rec = project->rec; rec; rec = rec->next ) {
				rec->select = 0;
				for ( part = rec->part; part; part = part->next ) if ( part->sel == n ) {
					rec->select = 1;
					c++;
				}
			}
		}
		if ( c ) {
			fn = filename.pre_rev('.') + Bstring(n, "_%04d.") + filename.post_rev('.');
			if ( verbose )
				cout << "Writing " << fn << " with " << c << " particles" << endl;
			write_project(fn, project, 1, 1);
		}
	}
	
	return n;
}

/**
@brief 	Splits a project into sets of particle selections and write the files.
@param 	*project 	the project.
@return int 				number of field files generated.

	The relevant micrographs are selected for each particle selection number
	and written into a file.

**/
int			project_split_field_write(Bproject* project)
{
	Bproject*			project_temp = new Bproject;
	project_temp->comment = project->comment;

	int					n(0);
	Bfield*				field, *field_next;
	Bstring				fn;

	for ( field = project->field; field; field = field->next ) {
		field_next = field->next;
		field->next = NULL;
		project_temp->field = field;
		fn = field->id + ".star";
		write_project(fn, project_temp, 1, 1);
		field->next = field_next;
	}
		
	return n;
}

