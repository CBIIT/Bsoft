/**
@file	mg_img_proc.cpp
@brief	Functions for image processing from micrograph structures
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20210806
**/

#include "mg_img_proc.h"
#include "mg_processing.h"
#include "mg_select.h"
#include "mg_ctf.h"
#include "img_combine.h"
#include "qsort_functions.h"
#include "rwimg.h"
#include "linked_list.h"
#include "math_util.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bmicrograph*		micrograph_create_from_image(Bimage* p, Bstring type);
Breconstruction*	reconstruction_create_from_image(Bimage* p, Bstring type);

/**
@brief 	Converts an image type string to a standard form.
@param 	type		a string with a type definition.
@return Bstring		standard type definition.

	The input string is converted to the following strings:
	mg		micrograph
	frame	micrograph frames
	rec		reconstruction
	part	particle
	fil		filament
	ps		power spectrum
	ft		Fourier transform

**/
Bstring		image_type(Bstring type)
{
	type = type.lower();
	
	if ( type.contains("micro") ) type = "mg";
	else if ( type.contains("frame") ) type = "frame";
	else if ( type.contains("rec") ) type = "rec";
	else if ( type.contains("part") ) type = "part";
	else if ( type.contains("fil") ) type = "fil";
	else if ( type.contains("power") ) type = "ps";
	else if ( type.contains("trans") ) type = "ft";
	
	return type;
}

/**
@brief 	Creates a project structure using image file names.
@param 	*p			list of images.
@param 	type		type of images: mg, frame, rec, part, fil.
@return Bproject*		new project parameter structure.

	The function sets up the project hierarchy from a list of file names.
	Each file may represent a micrograph or a picked particle file.
	If the image is equal or larger than 1024x1024, it is assumed to be
	a micrograph, its name will be assigned as a micrograph
	file name, and no particle tags will be added.
	If the image is smaller than 1024x1024 or the make_part flag is set,
	it is taken to be picked particles and the file name will 
	be assigned as a particle file name.

**/
Bproject*	project_create_from_image_old(Bimage* p, Bstring type)
{
	long				j, k, pid, nid, flen, finc(1);
	Vector3<double>		pixel_size(p->sampling(0));
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bframe*				frame = NULL;
	Bparticle*			part = NULL;
	Bparticle**			partlist = NULL;
	Bfilament*			fil = NULL;
	Bfilament**			fillist = NULL;
	Bfilnode*			fnode = NULL;
	Vector3<double>		start, fildir;
	Bstring				id, base(p->file_name());
	base = base.alnum();
	
	type = image_type(type);

	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			if ( p->sizeZ() > 1 ) {
				type = "rec";									// Assume reconstructions
				project->select = 1;
			} else type = "mg";									// Assume micrographs
		} else {
			if ( p->sizeZ() > 1 ) {
				project->select = 1;
				if ( p->sizeX() != p->sizeY() ) type = "fil";	// Assume filaments
				else type = "part";								// Assume particles
			} else type = "frame";								// Assume frames
		}
	}
	
	if ( project->select < 1 ) {
		if ( type == "frame" ) if ( p->sizeZ() > 1 ) p->slices_to_images();
		field = field_add(&project->field, base);
		if ( verbose & VERB_FULL )
			cout << "Creating field \"" << field->id << endl;
		k = 1;
		if ( type == "mg" ) k = p->images();
		for ( j=0; j<k; j++ ) {
			id = base;
			if ( k > 1 ) id += Bstring(j+1, "_%03d");
			mg = micrograph_add(&field->mg, id);
			if ( verbose & VERB_FULL )
				cout << "Micrograph \"" << mg->id << endl;
			if ( type == "mg" ) mg->fmg = p->file_name();
			else if ( type == "frame" ) mg->fframe = p->file_name();
			mg->img_num = j;
			mg->block = j;
			mg->pixel_size = pixel_size;
			mg->origin = p->size()/2;
			mg->intensity = p->image[j].average();
			if ( type == "mg" ) mg->origin = p->image[j].origin();
			mg->matrix = p->image[j].view().matrix();
			mg->select = 1;
		}
		if ( type == "frame" ) {
			k = p->images();
			mg->frame_pixel_size = pixel_size;
			mg->origin = p->size()/2;
			mg->intensity = 0;
			for ( j=0; j<k; j++ ) {
				frame = frame_add(&mg->frame, j+1);
				frame->shift = p->image[j].origin() - mg->origin;
				mg->intensity += p->image[j].average();
			}
		}
		mg = field_find_zero_tilt_mg(field);
		field->origin = mg->origin;
		if ( type == "part" ) {
			mg->fpart = p->file_name();
			partlist = &mg->part;
			pid = 0;
			if ( verbose & VERB_FULL )
				cout << " with " << p->images() << " particles";
		} else if ( type == "fil" ) {
			mg->ffil = p->file_name();
			fillist = &mg->fil;
			fildir = Vector3<double>(1,0,0);
			start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, 0);
			mg->filament_width = p->sizeY();
			if ( p->sizeX() < p->sizeY() ) {
				mg->filament_width = p->sizeX();
				fildir = Vector3<double>(0,1,0);
				start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, 0);
			}
			finc = (long) mg->filament_width/2;
			pid = 0;
			if ( verbose & VERB_FULL )
				cout << " with " << p->images() << " filaments";
		}
		if ( verbose & VERB_FULL )
			cout << endl;
	} else {
		rec = reconstruction_add(&project->rec, base);
		if ( verbose & VERB_FULL )
			cout << "Creating reconstruction \"" << rec->id << "\"" << endl;
		if ( type == "rec" ) rec->frec = p->file_name();
		else if ( type == "part" ) rec->fpart = p->file_name();
		else if ( type == "fil" ) rec->ffil = p->file_name();
		rec->voxel_size = pixel_size;
		if ( type == "rec" ) rec->origin = p->image->origin();
		rec->select = 1;
		if ( type == "part" ) partlist = &rec->part;
		if ( type == "fil" ) {
			fillist = &rec->fil;
			fildir = Vector3<double>(1,0,0);
			start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, p->sizeY()/2);
			rec->filament_width = p->sizeY();
			if ( p->sizeX() < p->sizeY() ) {
				rec->filament_width = p->sizeX();
				fildir = Vector3<double>(0,1,0);
				start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, p->sizeX()/2);
			}
			finc = (long) rec->filament_width/2;
		}
	}

	if ( partlist ) {
		for ( j=pid=0; j<p->images(); j++ ) {
			part = particle_add(partlist, ++pid);
			part->pixel_size = pixel_size;
			part->ori = p->image[j].origin();
			if ( fabs(part->ori[0] - p->sizeX()/2) > p->sizeX()/4 ) part->ori[0] = p->sizeX()/2;
			if ( fabs(part->ori[1] - p->sizeY()/2) > p->sizeY()/4 ) part->ori[1] = p->sizeY()/2;
			if ( fabs(part->ori[2] - p->sizeZ()/2) > p->sizeZ()/4 ) part->ori[2] = p->sizeZ()/2;
			part->view = p->image[j].view();
			part->sel = j+1;
		}
	}
		
	if ( fillist ) {
		flen = p->sizeX() - 2*finc;
		if ( p->sizeY() > p->sizeX() ) flen = p->sizeY() - 2*finc;
		for ( j=0; j<p->images(); j++ ) {
			fil = filament_add(fillist, ++pid);
			for ( nid = 0, k = 0; k < flen; k += finc ) {
				fnode = filament_node_add(&fil->node, ++nid);
				fnode->loc = start + fildir * k;
			}
		}
	}
	
	return project;
}

/**
@brief 	Creates a project structure using image file names.
@param 	*file_list	list of file names.
@param 	type		type of images: mg, frame, rec, part, fil.
@return Bproject*	new project parameter structure.

	The function sets up the project hierarchy from a list of file names.
	Each file may represent a micrograph or a picked particle file.
	If the image is equal or larger than 1024x1024, it is assumed to be
	a micrograph, its name will be assigned as a micrograph
	file name, and no particle tags will be added.
	If the image is smaller than 1024x1024 or the make_part flag is set,
	it is taken to be picked particles and the file name will 
	be assigned as a particle file name.

**/
Bproject*	project_create_from_images_old(Bstring* file_list, Bstring type)
{
	if ( file_list->empty() ) {
		error_show("Error in project_create_from_images: No filenames given to create a project!", __FILE__, __LINE__);
		return NULL;
	}
	
	type = image_type(type);

	Bimage*				p = read_img(*file_list, 0, -1);
	
	if ( !p ) {
		cerr << "Error: Input image file " << *file_list << " not found!" << endl;
		return NULL;
	}
	
	int					multi(0);	// Flag for multiple files
	
	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			if ( p->sizeZ() > 1 ) type = "rec";				// Assume reconstructions
			else type = "mg";								// Assume micrographs
		} else {
			if ( p->sizeX() != p->sizeY() ) type = "fil";	// Assume filaments
			else type = "part";								// Assume particles
		}
	}
	
//	if ( type == "part" && p->images() == 1 ) multi = 1;	// Single particle files - assume one mg or rec only
	
	if ( verbose )
		cout << "Creating a project from images of type \"" << type << "\"" << endl << endl;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_create_from_images: type=" << type << " multi=" << multi << endl;
		
	long				i, j, k, pid(0), nid, flen, finc(1);
	Vector3<double>		pixel_size(p->sampling(0));
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bframe*				frame = NULL;
	Bparticle*			part = NULL;
	Bparticle**			partlist = NULL;
	Bfilament*			fil = NULL;
	Bfilament**			fillist = NULL;
	Bfilnode*			fnode = NULL;
	Bstring*			thisfile = NULL;
	Bstring				base("1"), id;
	Vector3<double>		start, fildir;

	if ( ( type == "part" || type == "mg" || type == "frame" ) && multi ) {
		if ( type == "frame" ) if ( p->sizeZ() > 1 ) p->slices_to_images();
		if ( p->sizeZ() == 1 ) {
			field = field_add(&project->field, base);
			mg = micrograph_add(&field->mg, base);
			if ( verbose & VERB_FULL )
				cout << "Creating field \"" << field->id << "\", micrograph \"" << mg->id << "\"" << endl;
			mg->select = 1;
			if ( type == "part" ) partlist = &mg->part;
			if ( type == "fil" ) fillist = &mg->fil;
		} else {
			rec = reconstruction_add(&project->rec, base);
			if ( verbose & VERB_FULL )
				cout << "Creating reconstruction \"" << rec->id << "\"" << endl;
			rec->select = 1;
			if ( type == "part" ) partlist = &rec->part;
			if ( type == "fil" ) fillist = &rec->fil;
		}
	}
	
	delete p;
	
	for ( i=0, thisfile = file_list; thisfile; thisfile = thisfile->next, i++ ) {
		base = thisfile->base().alnum();
		p = read_img(*thisfile, 0, -1);
		if ( !multi ) {
			partlist = NULL;
			fillist = NULL;
			if ( type == "frame" || type == "mg" ) if ( p->sizeZ() > 1 ) p->slices_to_images();
			if ( p->image->origin().length() < 1 ) p->origin(p->size()/2);
			if ( p->sizeZ() == 1 ) {		// 2D files are assumed to be micrographs
				field = field_add(&project->field, base);
				if ( verbose & VERB_FULL )
					cout << "Creating field \"" << field->id << endl;
				k = 1;
				if ( type == "mg" ) k = p->images();
				for ( j=0; j<k; j++ ) {
					id = base;
					if ( k > 1 ) id += Bstring(j+1, "_%03d");
					mg = micrograph_add(&field->mg, id);
					if ( verbose & VERB_FULL )
						cout << "Micrograph \"" << mg->id << endl;
					if ( type == "mg" ) mg->fmg = *thisfile;
					else if ( type == "frame" ) mg->fframe = *thisfile;
					mg->img_num = j;
					mg->block = i + j;
//					mg->pixel_size = p->sampling(0);
					mg->pixel_size = pixel_size;
					mg->intensity = p->image[j].average();
					if ( type == "mg" ) mg->origin = p->image[j].origin();
					mg->origin[2] = 0;
					if ( mg->origin.length() < 1 ) mg->origin = p->size()/2;
					mg->matrix = p->image[j].view().matrix();
					mg->select = 1;
				}
				if ( type == "frame" ) {
					k = p->images();
//					mg->frame_pixel_size = p->sampling(0);
					mg->frame_pixel_size = pixel_size;
					mg->origin = p->size()/2;
					mg->intensity = 0;
					for ( j=0; j<k; j++ ) {
						frame = frame_add(&mg->frame, j+1);
						frame->shift = p->image[j].origin() - mg->origin;
						mg->intensity += p->image[j].average();
					}
				}
				mg = field_find_zero_tilt_mg(field);
				field->origin = mg->origin;
				if ( type == "part" ) {
					mg->fpart = *thisfile;
					partlist = &mg->part;
					pid = 0;
					if ( verbose & VERB_FULL )
						cout << " with " << p->images() << " particles";
				} else if ( type == "fil" ) {
					mg->ffil = *thisfile;
					fillist = &mg->fil;
					pid = 0;
					if ( verbose & VERB_FULL )
						cout << " with " << p->images() << " filaments";
				}
				if ( verbose & VERB_FULL )
					cout << endl;
			} else {						// 3D files are assumed to be reconstructions
				rec = reconstruction_add(&project->rec, base);
				if ( verbose & VERB_FULL )
					cout << "Creating reconstruction \"" << rec->id << "\"";
				rec->block = i;
				rec->select = 1;
				if ( type == "part" ) {
					rec->fpart = *thisfile;
					partlist = &rec->part;
					pid = 0;
					if ( verbose & VERB_FULL )
						cout << " with " << p->images() << " particles";
				} else if ( type == "fil" ) {
					rec->ffil = *thisfile;
					fillist = &rec->fil;
					pid = 0;
					if ( verbose & VERB_FULL )
						cout << " with " << p->images() << " filaments";
				} else {
					rec->frec = *thisfile;
					rec->origin = p->image->origin();
				}
				if ( verbose & VERB_FULL )
					cout << endl;
			}
		}

		if ( mg ) {
			if ( type == "part" ) {
				mg->box_size = p->size();
			} else if ( type == "fil" ) {
				fildir = Vector3<double>(1,0,0);
				start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, 0);
				mg->filament_width = p->sizeY();
				if ( p->sizeX() < p->sizeY() ) {
					mg->filament_width = p->sizeX();
					fildir = Vector3<double>(0,1,0);
					start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, 0);
				}
				finc = (long) mg->filament_width/2;
			} else if ( type == "ft" || type == "ps" ) {
				mg->filament_width = p->sizeY();
				mg->helix_axis = 0;
				if ( p->sizeX() < p->sizeY() ) {
					mg->filament_width = p->sizeX();
					mg->helix_axis = M_PI_2;
				}
			}
		} else if ( rec ) {
			rec->voxel_size = p->sampling(0);
			if ( type == "part" ) {
				rec->box_size = p->size();
			} else if ( type == "fil" ) {
				fildir = Vector3<double>(1,0,0);
				start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, p->sizeY()/2);
				rec->filament_width = p->sizeY();
				if ( p->sizeX() < p->sizeY() ) {
					rec->filament_width = p->sizeX();
					fildir = Vector3<double>(0,1,0);
					start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, p->sizeX()/2);
				}
				finc = (long) rec->filament_width/2;
			} else if ( type == "ft" || type == "ps" ) {
				rec->filament_width = p->sizeX();
				if ( p->sizeY() < rec->filament_width ) rec->filament_width = p->sizeY();
				if ( p->sizeZ() < rec->filament_width ) rec->filament_width = p->sizeZ();
			}
		}
		
		if ( partlist ) {
			for ( j=0; j<p->images(); j++ ) {
				part = particle_add(partlist, ++pid);
				part->ori = p->image[j].origin();
				if ( fabs(part->ori[0] - p->sizeX()/2) > p->sizeX()/4 ) part->ori[0] = p->sizeX()/2;
				if ( fabs(part->ori[1] - p->sizeY()/2) > p->sizeY()/4 ) part->ori[1] = p->sizeY()/2;
				if ( fabs(part->ori[2] - p->sizeZ()/2) > p->sizeZ()/4 ) part->ori[2] = p->sizeZ()/2;
				part->view = p->image[j].view();
//				cout << part->id << tab << part->view << endl;
//				part->sel = 1;
				part->sel = j+1;
			}
			if ( multi ) part->fpart = *thisfile;
		}
		
		if ( fillist ) {
			flen = p->sizeX() - 2*finc;
			if ( p->sizeY() > p->sizeX() ) flen = p->sizeY() - 2*finc;
			for ( j=0; j<p->images(); j++ ) {
				fil = filament_add(fillist, ++pid);
				for ( nid = 0, k = 0; k < flen; k += finc ) {
					fnode = filament_node_add(&fil->node, ++nid);
					fnode->loc = start + fildir * k;
				}
			}
			if ( multi ) fil->ffil = *thisfile;
		}
		
		delete p;
	}
	
	return project;
}

/**
@brief 	Creates a project structure using image file names.
@param 	*p			image.
@param 	type		type of images: mg, frame, rec, part, fil.
@return Bproject*		new project parameter structure.

	The function sets up the project hierarchy from one image base on the type.
	If the type is not specified, it is guessed based on the following rules:
		#sub-images = 1
			z=1 => micrograph
			z>1 => reconstruction
		#sub-images > 1
			z=1 => frames
			z>1
				x=y => particles
				x≠y => filaments

**/
Bproject*	project_create_from_image(Bimage* p, Bstring type)
{
	type = image_type(type);

	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			if ( p->sizeZ() > 1 ) type = "rec";					// Assume reconstructions
			else type = "mg";									// Assume micrographs
		} else {
			if ( p->sizeZ() > 1 )  {
				if ( p->sizeX() != p->sizeY() ) type = "fil";	// Assume filaments
				else type = "part";								// Assume particles
			} else type = "frame";								// Assume frames
		}
	}
	
	if ( verbose )
		cout << "Creating a project from images of type \"" << type << "\"" << endl << endl;
		
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Vector3<double>		start, fildir;
	Bstring				id, base(p->file_name());
	base = base.alnum();

	if ( type == "frame" || type == "mg" ) if ( p->sizeZ() > 1 ) p->slices_to_images();
		
	if ( p->image->origin().length() < 1 ) p->origin(p->size()/2);
		
	if ( p->sizeZ() == 1 ) {		// 2D files are assumed to be micrographs
		field = field_add(&project->field, base);
		if ( verbose & VERB_FULL )
			cout << "Creating field \"" << field->id << endl;
		field->mg = micrograph_create_from_image(p, type);
		mg = field_find_zero_tilt_mg(field);
		field->origin = mg->origin;
		if ( verbose & VERB_FULL )
			cout << endl;
	} else {						// 3D files are assumed to be reconstructions
		project->rec = rec = reconstruction_create_from_image(p, type);
		if ( verbose & VERB_FULL )
			cout << "Creating reconstruction \"" << rec->id << "\"";
		if ( verbose & VERB_FULL )
			cout << endl;
	}
	
	return project;
}

/**
@brief 	Creates a project structure using image file names.
@param 	*file_list	list of file names.
@param 	type		type of images: mg, frame, rec, part, fil.
@return Bproject*	new project parameter structure.

	The function sets up the project hierarchy from a list of file names.
	Each file may represent a micrograph or a picked particle file.
	If the image is equal or larger than 1024x1024, it is assumed to be
	a micrograph, its name will be assigned as a micrograph
	file name, and no particle tags will be added.
	If the image is smaller than 1024x1024 or the make_part flag is set,
	it is taken to be picked particles and the file name will
	be assigned as a particle file name.

**/
Bproject*	project_create_from_images(Bstring* file_list, Bstring type)
{
	if ( file_list->empty() ) {
		error_show("Error in project_create_from_images: No filenames given to create a project!", __FILE__, __LINE__);
		return NULL;
	}
	
	type = image_type(type);

	Bimage*				p = read_img(*file_list, 0, -1);
	
	if ( !p ) {
		cerr << "Error: Input image file " << *file_list << " not found!" << endl;
		return NULL;
	}
	
	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			if ( p->sizeZ() > 1 ) type = "rec";					// Assume reconstructions
			else type = "mg";									// Assume micrographs
		} else {
			if ( p->sizeZ() > 1 )  {
				if ( p->sizeX() != p->sizeY() ) type = "fil";	// Assume filaments
				else type = "part";								// Assume particles
			} else type = "frame";								// Assume frames
		}
	}
	
	if ( verbose )
		cout << "Creating a project from images of type \"" << type << "\"" << endl << endl;
		
	long				i;
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bstring*			thisfile = NULL;
	Bstring				base("1"), id;
	Vector3<double>		start, fildir;

	delete p;
	
	for ( i=0, thisfile = file_list; thisfile; thisfile = thisfile->next, i++ ) {
		base = thisfile->base().alnum();
		p = read_img(*thisfile, 0, -1);
//		cerr << "file=" << *thisfile << tab << p->file_name() << endl;
		
		if ( type == "frame" || type == "mg" ) if ( p->sizeZ() > 1 ) p->slices_to_images();
		
		if ( p->image->origin().length() < 1 ) p->origin(p->size()/2);
		
		if ( p->sizeZ() == 1 ) {		// 2D files are assumed to be micrographs
			field = field_add(&project->field, base);
			if ( verbose & VERB_FULL )
				cout << "Creating field \"" << field->id << endl;
			field->mg = micrograph_create_from_image(p, type);
			mg = field_find_zero_tilt_mg(field);
			field->origin = mg->origin;
			if ( verbose & VERB_FULL )
				cout << endl;
		} else {						// 3D files are assumed to be reconstructions
			if ( rec ) rec = rec->next = reconstruction_create_from_image(p, type);
			else project->rec = rec = reconstruction_create_from_image(p, type);
			if ( verbose & VERB_FULL )
				cout << "Creating reconstruction \"" << rec->id << "\"";
			rec->block = i;
			if ( verbose & VERB_FULL )
				cout << endl;
		}

		delete p;
	}
	
	return project;
}

Bparticle*		particles_from_image(Bimage* p)
{
	long			j, pid;
	Bparticle*		partlist = NULL;
	Bparticle*		part;
	
	for ( j=pid=0; j<p->images(); j++ ) {
		part = particle_add(&partlist, ++pid);
		part->pixel_size = p->sampling(j);
		part->ori = p->image[j].origin();
		if ( fabs(part->ori[0] - p->sizeX()/2) > p->sizeX()/4 ) part->ori[0] = p->sizeX()/2;
		if ( fabs(part->ori[1] - p->sizeY()/2) > p->sizeY()/4 ) part->ori[1] = p->sizeY()/2;
		if ( fabs(part->ori[2] - p->sizeZ()/2) > p->sizeZ()/4 ) part->ori[2] = p->sizeZ()/2;
		part->view = p->image[j].view();
		part->sel = pid;
	}
	
	return partlist;
}

Bfilament*		filaments_from_image(Bimage* p, Vector3<double> fildir, Vector3<double> start, long finc)
{
	long			j, k, pid, nid, flen;
	Bfilament*		fillist = NULL;
	Bfilament*		fil = NULL;
	Bfilnode*		fnode = NULL;

	if ( p->sizeX() > p->sizeY() ) flen = p->sizeX() - 2*finc;
	else flen = p->sizeY() - 2*finc;
	
	for ( j=pid=0; j<p->images(); j++ ) {
		fil = filament_add(&fillist, ++pid);
		for ( nid=k=0; k < flen; k += finc ) {
			fnode = filament_node_add(&fil->node, ++nid);
			fnode->loc = start + fildir * k;
		}
	}
	
	return fillist;
}

Bmicrograph*	micrograph_create_from_image(Bimage* p, Bstring type)
{
	type = image_type(type);

	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			type = "mg";								// Assume micrographs
		} else {
			type = "frame";								// Assume frames
		}
	}

	if ( type == "frame" ) if ( p->sizeZ() > 1 ) p->slices_to_images();

	if ( p->sizeZ() > 1 ) {
		cerr << "Error: The image is 3D and a micrograph record cannot be created!" << endl;
		return NULL;
	}
	
	long				j, k(1), finc(1);
	Vector3<double>		pixel_size(p->sampling(0));
	Bmicrograph*		mglist = NULL;
	Bmicrograph*		mg = NULL;
	Bframe*				frame = NULL;
	Vector3<double>		start, fildir;
	Bstring				id, base(p->file_name());
	base = base.pre_rev('.').alnum();
	
	if ( type == "mg" ) k = p->images();
	
	for ( j=0; j<k; j++ ) {
		id = base;
		if ( k > 1 ) id += Bstring(j+1, "_%03d");
		mg = micrograph_add(&mglist, id);
		if ( verbose & VERB_FULL )
			cout << "Micrograph \"" << mg->id << endl;
		if ( type == "mg" ) mg->fmg = p->file_name();
		else if ( type == "frame" ) mg->fframe = p->file_name();
		mg->img_num = j;
		mg->block = j;
		mg->pixel_size = pixel_size;
		mg->origin = p->size()/2;
//		mg->intensity = p->image[j].average();
		if ( type == "mg" ) mg->origin = p->image[j].origin();
		mg->matrix = p->image[j].view().matrix();
		mg->select = 1;
		micrograph_intensity(mg, p, 1);
		if ( p->meta_data().exists("dose") )
			mg->dose = (*p)["dose"].real();
//		else if ( type == "frame" )
//			mg->dose = p->average()*p->images()/(pixel_size[0]*pixel_size[1]);
//		else
//			mg->dose = p->average()/(pixel_size[0]*pixel_size[1]);
		if ( p->meta_data().exists("exposure") )
			mg->exposure = (*p)["exposure"].real();
	}
	
//	cerr << type << tab << k << tab << p->file_name() << tab << mg->fframe << endl;
			
	if ( type == "frame" ) {
		k = p->images();
		mg->frame_pixel_size = pixel_size;
		mg->origin = p->size()/2;
//		mg->intensity = 0;
		for ( j=0; j<k; j++ ) {
			frame = frame_add(&mg->frame, j+1);
			frame->shift = p->image[j].origin() - mg->origin;
//			mg->intensity += p->image[j].average();
		}
		if ( mg->fframe.length() < 1 )
			cerr << "Error: No frame file name assigned!" << endl;
	}

	if ( type == "ft" || type == "ps" ) {
		mg->filament_width = p->sizeY();
		mg->helix_axis = 0;
		if ( p->sizeX() < p->sizeY() ) {
			mg->filament_width = p->sizeX();
			mg->helix_axis = M_PI_2;
		}
	}

	if ( type == "part" ) {
		mg->fpart = p->file_name();
		mg->box_size = p->size();
		mg->part = particles_from_image(p);
		if ( verbose & VERB_FULL )
			cout << " with " << p->images() << " particles";
		if ( mg->fpart.length() < 1 )
			cerr << "Error: No particle file name assigned!" << endl;
	} else if ( type == "fil" ) {
		mg->ffil = p->file_name();
		fildir = Vector3<double>(1,0,0);
		start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, 0);
		mg->filament_width = p->sizeY();
		if ( p->sizeX() < p->sizeY() ) {
			mg->filament_width = p->sizeX();
			fildir = Vector3<double>(0,1,0);
			start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, 0);
		}
		finc = (long) mg->filament_width/2;
		mg->fil = filaments_from_image(p, fildir, start, finc);
		if ( verbose & VERB_FULL )
			cout << " with " << p->images() << " filaments";
		if ( mg->ffil.length() < 1 )
			cerr << "Error: No filament file name assigned!" << endl;
	}
	if ( verbose & VERB_FULL )
		cout << endl;

	return mglist;
}

Breconstruction*	reconstruction_create_from_image(Bimage* p, Bstring type)
{
	type = image_type(type);

	if ( type.length() < 1 ) {
		if ( p->images() == 1 ) {
			type = "mg";								// Assume micrographs
		} else {
			type = "frame";								// Assume frames
		}
	}

	long				finc(1);
	Vector3<double>		pixel_size(p->sampling(0));
	Vector3<double>		start, fildir;
	Breconstruction*	rec = NULL;
	Bstring				id, base(p->file_name());
	base = base.alnum();
	
	rec = reconstruction_add(&rec, base);
	
	if ( verbose & VERB_FULL )
		cout << "Creating reconstruction \"" << rec->id << "\"";
	rec->select = 1;
	rec->voxel_size = pixel_size;
	if ( type == "part" ) {
		rec->fpart = p->file_name();
		rec->box_size = p->size();
		rec->part = particles_from_image(p);
		if ( verbose & VERB_FULL )
			cout << " with " << p->images() << " particles";
	} else if ( type == "fil" ) {
		rec->ffil = p->file_name();
		fildir = Vector3<double>(1,0,0);
		start = Vector3<double>(p->sizeY()/2, p->sizeY()/2, p->sizeY()/2);
		rec->filament_width = p->sizeY();
		if ( p->sizeX() < p->sizeY() ) {
			rec->filament_width = p->sizeX();
			fildir = Vector3<double>(0,1,0);
			start = Vector3<double>(p->sizeX()/2, p->sizeX()/2, p->sizeX()/2);
		}
		finc = (long) rec->filament_width/2;
		rec->fil = filaments_from_image(p, fildir, start, finc);
		if ( verbose & VERB_FULL )
			cout << " with " << p->images() << " filaments";
	} else if ( type == "ft" || type == "ps" ) {
		rec->filament_width = p->sizeX();
		if ( p->sizeY() < rec->filament_width ) rec->filament_width = p->sizeY();
		if ( p->sizeZ() < rec->filament_width ) rec->filament_width = p->sizeZ();
	} else {
		rec->frec = p->file_name();
		rec->origin = p->image->origin();
	}

	if ( verbose & VERB_FULL )
		cout << endl;

	return rec;
}

/**
@brief 	Gets the average of a micrograph and sets the dose.
@param 	*mg			micrograph parameter structure.
@param	*p			image.
@param 	flag		1=force calculation of statistics of not available; 2=check Poisson
@return double		intensity.

	The micrograph image header is read.

**/
double		micrograph_intensity(Bmicrograph* mg, Bimage* p, int flag)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_intensity: " << mg->id << tab << flag << endl;

	int				reread(0);
	
	if ( flag ) {
		if ( p->standard_deviation() < 1e-10 ) {
			reread = 1;
			if ( mg->fmg.length() ) p = read_img(mg->fmg, 1, mg->img_num);
			else if ( mg->fframe.length() ) p = read_img(mg->fframe, 1, -1);
			p->statistics();
		}
	}
	
	mg->intensity = p->average();
	
	if ( flag&2 ) p->poisson_statistics_check();
	
	if ( mg->dose <= 0 )
		mg->dose = p->images()*p->average()*p->average()/(p->sampling(0)[0]*p->sampling(0)[1]*p->standard_deviation()*p->standard_deviation());
	
	if ( reread ) delete p;
	
	return mg->intensity;
}

/**
@brief 	Gets the average of a micrograph and sets the dose.
@param 	*mg			micrograph parameter structure.
@param 	flag		flag to force calculation of statistics of not available.
@return double		intensity.

	The micrograph image header is read.

**/
double		micrograph_intensity(Bmicrograph* mg, int flag)
{
	Bimage*				p = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_intensity: " << mg->id << tab << flag << endl;

	if ( mg->fmg.length() ) p = read_img(mg->fmg, 0, mg->img_num);
	else if ( mg->fframe.length() ) p = read_img(mg->fframe, 0, -1);

	if ( !p ) {
		error_show("micrograph_intensity", __FILE__, __LINE__);
		return 0;
	}
	
	micrograph_intensity(mg, p, flag);

	delete p;
	
	return mg->intensity;
}

/**
@brief 	Retrieves the average micrograph intensities.
@param 	*project	project parameter structure.
@return int 			0.

	For each micrograph the FOM is set to the micrograph average.

**/
int			project_mg_avg_intensities(Bproject* project)
{
	Bfield*				field;
	Bmicrograph*		mg;

	if ( verbose & VERB_FULL )
		cout << "Getting the average micrograph intensities" << endl << endl;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			micrograph_intensity(mg, 1);
	
	return 0;
}

/**
@brief 	Concatenates all micrograph image files into one file.
@param 	*project		project structure.
@return int				0.

	The new file name is the common part of the original
	micrograph file names.

**/
int			project_catenate_micrographs(Bproject* project)
{
	long			i;
	Bfield*			field;
	Bmicrograph*	mg, *mg1;
	Bstring*		file_list = NULL;
	Bimage*			pcat = NULL;
	
	Bstring			rawstring, filename;
	Vector3<long>	nusize;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = mg1 = field->mg; mg; mg = mg->next )
			if ( mg->fmg != mg1->fmg )
				string_add(&file_list, mg->fmg);
		if ( file_list ) {
			if ( verbose )
				cout << "Catenating micrographs in field " << field->id << endl;
//			filename = file_list->common(*(file_list->next));
			filename = field->id;
			filename += "." + file_list->extension();
			cout << "New file name = " << filename << endl;
			pcat = img_catenate(file_list, rawstring, Unknown_Type, 
				nusize, 0, FILL_USER, 0, 0, 0);
			if ( verbose )
				cout << "Writing concatenated file " << filename << endl;
			write_img(filename, pcat, 0);
			for ( i=0, mg = field->mg; mg; mg = mg->next, ++i ) {
				mg->fmg = filename;
				mg->img_num = i;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Reads a particle image file.
@param 	*part		particle.
@param 	readflag	flag to indicate reading the data.
@return Bimage*		image (NULL means failure).

	The file name is taken from the particle record by preference,
	otherwise from the micrograph record.

**/
Bimage*		particle_read_img(Bparticle* part, int readflag)
{
	if ( !part ) return NULL;
	
	Bstring			filename(part->fpart);
	long			img_num(0);
	
	if ( !part->mg ) cerr << "Error: no micrograph for particle " << part->id << endl;
	
	if ( filename.length() < 2 ) {
		if ( part->mg ) {
			filename = part->mg->fpart;
			img_num = part->id - 1;
		} else if ( part->rec ) {
			filename = part->rec->fpart;
			img_num = part->id - 1;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particle image filename:" << tab << filename << " (" << img_num << ")" << endl;
	
	Bimage*			p = NULL;
	
	if ( filename.length() ) {
		p = read_img(filename, readflag, img_num);
		if ( !p ) {
			error_show("Particle file not read", __FILE__, __LINE__);
			return  NULL;
		}
	}
	
	return p;
}

/**
@brief 	Writes a particle image file.
@param 	*part		particle.
@param	*p			image to be written.
@param 	compression	flag to indicate compression.
@return int			images written.

	The file name is taken from the particle record by preference,
	otherwise from the micrograph record.

**/
int			particle_write_img(Bparticle* part, Bimage* p, int compression)
{
	if ( !part ) return 0;
	
	Bstring			filename(part->fpart);
	long			img_num(p->images());
	
	if ( !part->mg ) cerr << "Error: no micrograph for particle " << part->id << endl;
	
	if ( filename.length() < 2 ) {
		if ( part->mg )
			filename = part->mg->fpart;
		else if ( part->rec )
			filename = part->rec->fpart;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particle image filename:" << tab << filename << " (" << img_num << ")" << endl;
	
	if ( filename.length() ) {
    	write_img(filename, p, compression);
		if ( !p ) {
			error_show("Particle file not written", __FILE__, __LINE__);
			return 0;
		}
	}
	
	return 1;
}

/**
@brief 	Checks particle numbers and box sizes.
@param 	*project	project.
@return int			error code (<0 means failure).

	Reads each particle file header, checks the number of particles
	and sets the box size for the micrograph.

**/
int			project_check_particles(Bproject* project)
{
	if ( !project ) return -1;

	int				oldverbose = verbose;
	verbose = 0;
	
	int				err(0);
	long			n;
	Bimage*       	p = NULL;      // image currently being processed
	Bfield*       	field = NULL;  // field-of-view currently being processed
	Bmicrograph*  	mg = NULL;     // micrograph currently being processed
	Bparticle*    	part = NULL;   // particle currently being processed

	cout << "Checking particle numbers:" << endl;
	
	for ( field = project->field; field; field = field->next) {
		for ( mg = field->mg; mg; mg = mg->next) if ( mg->fpart.length() ) {
			p = read_img(mg->fpart, 0, -1);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_check_particles: " << mg->fpart << ":" << p->images() << endl;
			if ( p ) {
				mg->box_size = p->size();
				for ( n=0, part = mg->part; part; part = part->next, n++ ) {
					if ( part->pixel_size[0] < 0.01 )
						part->pixel_size = mg->pixel_size;
					if ( part->ori[0]*part->ori[1] < mg->box_size.volume()/16 )
						part->ori += mg->box_size/2;
				}
				if ( n != p->images() ) {
					cerr << "Error: Micrograph " << mg->id << " has " << n << " particles," << endl;
					cerr << "    but particle image file " << p->file_name() << " has " << p->images() << " images!" << endl;
					err--;
				} else {
					cout << mg->id << ": " << n << " = " << p->images() << endl;
				}
				delete p;
				p = NULL;
			}
		}
	}
	
	verbose = oldverbose;

	return err;
}

/**
@brief 	Gets the size of a micrograph.
@param 	*mg			micrograph parameter structure.
@return Vector3<long> 	size.

	The micrograph image header is read.

**/
Vector3<long>	micrograph_get_size(Bmicrograph* mg)
{
	Vector3<long>		size;

	Bimage*				p = NULL;
	
	if ( mg->fmg.length() )
		p = read_img(mg->fmg, 0, 0);
	else if ( mg->fframe.length() )
		p = read_img(mg->fframe, 0, 0);
	
	if ( !p )
		error_show("micrograph_get_size", __FILE__, __LINE__);
	else {
		size = p->size();
		delete p;
	}
	
	return size;
}

/**
@brief 	Gets the nominal origin for a micrograph.
@param 	*mg				micrograph parameter structure.
@return Vector3<double> nominal origin.

	The nominal origin is defined as the center of the micrograph.

**/
Vector3<double>	micrograph_get_nominal_origin(Bmicrograph* mg)
{
	Vector3<double>		origin;

	Bimage*				p = NULL;
	
	if ( mg->fmg.length() )
		p = read_img(mg->fmg, 0, 0);
	else if ( mg->fframe.length() )
		p = read_img(mg->fframe, 0, 0);
	
	if ( !p )
		error_show("micrograph_get_nominal_origin", __FILE__, __LINE__);
	else {
		origin = Vector3<double>(p->size()/2);
		delete p;
	}
	
	return origin;
}

/**
@brief 	Sets micrograph origins to the centers of the micrographs.
@param 	*project	project parameter structure.
@return int 		0.

	For each micrograph the micrograph origin is set to the center.

**/
int			project_set_nominal_mg_origins(Bproject* project)
{
	Bmicrograph*		mg;
	Bimage*				p = NULL;

	if ( verbose & VERB_FULL )
		cout << "Setting nominal micrograph origins" << endl << endl;

	for ( mg = project->field->mg; mg; mg = mg->next ) {
		p = read_img(mg->fmg, 0, mg->img_num);
		if ( !p )
			error_show("project_set_nominal_mg_origins", __FILE__, __LINE__);
		else {
			mg->origin = p->size()/2;
			delete p;
		}
	}
	
	return 0;
}

/**
@brief 	Resets particle origins to nominal values.
@param 	*project	project parameter structure.
@return int 		0.

	For each micrograph a particle image header is read and the origin
	is set to to the middle of the image.

**/
int			project_reset_origins(Bproject* project)
{
	if ( !project ) return 0;
	if ( !project->field ) return 0;
	if ( !project->field->mg ) return 0;
	if ( !project->field->mg->part ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;
	Bimage*			p;
	
	if ( verbose & VERB_FULL )
		cout << "Resetting origins" << endl << endl;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->fpart.length() ) p = read_img(mg->fpart, 0, -1);
			else p = read_img(mg->fmg, 0, -1);
			if ( !p )
				error_show("project_reset_origins", __FILE__, __LINE__);
			for ( part = mg->part; part; part = part->next )
				part->ori = p->size()/2;
			delete p;
		}

	return 0;
}

/**
@author  D. Belnap
@brief 	Writes particle origins into particle image files.
@param 	*project	project.
@return int			error code (<0 means failure).

	Sets the origins (offsets from the first voxel) in an image header to
	values set within a project.  Rewrites image to file. 

**/
int			project_set_part_img_origins(Bproject* project)
{
	if ( !project ) return -1;

	Bimage*       	p = NULL;      // image currently being processed
	Bfield*       	field = NULL;  // field-of-view currently being processed
	Bmicrograph*  	mg = NULL;     // micrograph currently being processed
	Bparticle*    	part = NULL;   // particle currently being processed

	for ( field = project->field; field; field = field->next) {
		for ( mg = field->mg; mg; mg = mg->next) {
			p = read_img(mg->fpart, 1, -1);
			for ( part = mg->part; part; part = part->next )
				p->image[part->id-1].origin(part->ori);
			write_img(mg->fpart, p, 0);
			delete p;
		}
	}

	return 0;
}

/**
@brief 	Deselects particles that are too close to the image edges.
@param 	*project	project.
@return long		number of particles selected.

**/
long		project_delesect_edge_particles(Bproject* project)
{
	if ( !project ) return -1;

	long			nsel(0);
	Bimage*       	p = NULL;      // image currently being processed
	Bfield*       	field = NULL;  // field-of-view currently being processed
	Bmicrograph*  	mg = NULL;     // micrograph currently being processed
	Bparticle*    	part = NULL;   // particle currently being processed
	Vector3<double>	min, max;

	for ( field = project->field; field; field = field->next) {
		for ( mg = field->mg; mg; mg = mg->next) {
			p = read_img(mg->fmg, 0, 0);
			min = mg->box_size/2;
			max = p->size() - min;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_delesect_edge_particles: limits:" << tab << min << tab << max << endl;
			for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
				if ( part->loc >= min || part->loc <= max ) nsel++;
				else part->sel = 0;
			}
			delete p;
		}
	}

	return 0;
}

/**
@brief 	Gets the box size from particle image files.
@param 	*project	project.
@return int			error code (<0 means failure).

	Reads each particle file header and sets the box size for the micrograph. 

**/
int			project_get_part_box_size(Bproject* project)
{
	if ( !project ) return -1;

	Bimage*       	p = NULL;      // image currently being processed
	Bfield*       	field = NULL;  // field-of-view currently being processed
	Bmicrograph*  	mg = NULL;     // micrograph currently being processed
	Bparticle*    	part = NULL;   // particle currently being processed

	for ( field = project->field; field; field = field->next) {
		for ( mg = field->mg; mg; mg = mg->next) if ( mg->fpart.length() ) {
			p = particle_read_img(mg->part, 0);
			if ( p ) {
				mg->box_size = p->size();
				delete p;
				for ( part = mg->part; part; part = part->next )
					if ( part->ori[0]*part->ori[1] < mg->box_size.volume()/16 )
						part->ori += mg->box_size/2;
			}
		}
	}

	return 0;
}

/**
@brief 	Gets the box size from particle image file.
@param 	*part			particle.
@return Vector3<long>	box size (0 if image does not exist).

	Reads the particle file header and returns the box size.

**/
Vector3<long>	particle_get_box_size(Bparticle* part)
{
	Vector3<long>	size;
	
	Bimage*			p = particle_read_img(part, 0);
	
	if ( p ) {
		size = p->size();
		delete p;
	}
	
	return size;
}

/**
@brief 	Writes particles to different stacks based on class.
@param 	*project	project.
@return int			error code (<0 means failure).

	Writes new particle image files, numbered by the selection number. 

**/
int			project_write_particle_classes(Bproject* project)
{
	if ( !project ) return -1;

	long	i, n(0);
	Bfield*       	field = NULL;
	Bmicrograph*  	mg = NULL;
	Bparticle*    	part = NULL;
	Bimage*      	p = NULL;
	Bimage*      	pt = NULL;
	Bstring			filename;

	for ( field = project->field; field; field = field->next)
		for ( mg = field->mg; mg; mg = mg->next)
			for ( part = mg->part; part; part = part->next )
				if ( part->sel > n ) n = part->sel;

	long*	cn = new long[n];
	Bimage**		pn = new Bimage*[n];
	
	if ( verbose )
		cout << "Writing " << n << " particle classes to separate files:" << endl;

	for ( field = project->field; field; field = field->next) {
		for ( mg = field->mg; mg; mg = mg->next) {
			if ( n ) {
				p = read_img(mg->fpart, 1, -1);
				for ( i=0; i<n; i++ ) cn[i] = 0;
				for ( part = mg->part; part; part = part->next )
					if ( part->sel > 0 ) cn[part->sel-1]++;
				for ( i=0; i<n; i++ ) if ( cn[i] ) {
					if ( verbose & VERB_FULL )
						cout << "Creating image with " << cn[i] << " particles" << endl;
					pn[i] = p->copy_header(cn[i]);
					pn[i]->data_alloc();
				}
				for ( i=0; i<n; i++ ) cn[i] = 0;
				for ( part = mg->part; part; part = part->next )
					if ( part->sel > 0 ) {
						pt = p->extract(part->id-1);
						pn[part->sel-1]->replace(cn[part->sel-1]++, pt);
						delete pt;
					}
				delete p;
				for ( i=0; i<n; i++ ) if ( cn[i] ) {
					filename = mg->fpart.pre_rev('.') + Bstring(i+1, "_%03d.") + mg->fpart.post_rev('.');
					if ( verbose )
						cout << "Writing " << filename << " with " << cn[i] << " particles" << endl;
					write_img(filename, pn[i], 0);
					delete pn[i];
				}
			}
		}
	}

	delete[] cn;
	delete[] pn;
	
	return 0;
}

/**
@brief 	Deletes selected class averages from the project..
@param 	*project	project parameter structure.
@param 	&list		selection list.
@return int 		0.

	The new class average file name has an insert of "_del".

**/
int			project_trim_class_averages(Bproject* project, Bstring& list)
{
	if ( ! project->class_avg ) return 0;
	if ( list.length() < 1 ) return 0;
	
	long			navg(0), i, id;
	Bparticle*		part, *part_list, *part_nu;

	for ( part = project->class_avg; part; part = part->next ) navg++;

//	int*				numsel = new int[navg];
//	select_numbers(list, navg, numsel);
	vector<int>		numsel = select_numbers(list, navg);
	
	Bstring			filename(project->class_avg->fpart);
	Bimage*			pavg = read_img(filename, 1, -1);
	pavg->delete_images(list);

	filename = filename.pre_rev('.') + "_del." + filename.post_rev('.');
	write_img(filename, pavg, 0);
	
	for ( i=id=0, part_list = NULL, part = project->class_avg; part; part = part->next, ++i ) {
		if ( !numsel[i] ) {
			part_nu = particle_copy(&part_list, part);
			part_nu->id = ++id;
			part_nu->fpart = filename;
		}
	}
	
	particle_kill(project->class_avg);
	project->class_avg = part_list;
	
//	delete[]			numsel;
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of class averages:       " << id << endl << endl;
	
	return 0;
}

/**
@brief 	Reverses one or more particle coordinates.
@param 	*project	project parameter structure.
@param 	flip		axes to flip.
@return int 		0.

	The specification of axes to flip is embedded in the flip number:
		first bit  - x
		second bit - y
		third bit  - z

**/
int			project_flip_particle_coordinates(Bproject* project, int flip)
{
	Bfield* 		field;
	Bmicrograph*	mg;
	Breconstruction*	rec;
	Bparticle*		part;
	Bbadarea*		bad;
	Bimage*			p = NULL;
	int				fz = flip/4;
	int				fy = (flip - 4*fz)/2;
	int				fx = flip - 4*fz - 2*fy;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Flipping locations along ";
		if ( fx ) cout << "x";
		if ( fy ) cout << "y";
		if ( fz ) cout << "z";
		cout << endl;
	}

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( rec->frec.length() ) p = read_img(rec->frec, 0, -1);
			for ( part = rec->part; part; part = part->next ) {
				if ( fx ) part->loc[0] = p->sizeX() - part->loc[0];
				if ( fy ) part->loc[1] = p->sizeY() - part->loc[1];
				if ( fz ) part->loc[2] = p->sizeZ() - part->loc[2];
			}
			for ( bad = rec->bad; bad; bad = bad->next ) {
				if ( fx ) bad->loc[0] = p->sizeX() - bad->loc[0];
				if ( fy ) bad->loc[1] = p->sizeY() - bad->loc[1];
				if ( fz ) bad->loc[2] = p->sizeZ() - bad->loc[2];
			}
			delete p;
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->fmg.length() ) p = read_img(mg->fmg, 0, -1);
				else if ( mg->fframe.length() ) p = read_img(mg->fframe, 0, -1);
				for ( part = mg->part; part; part = part->next ) {
					if ( fx ) part->loc[0] = p->sizeX() - part->loc[0];
					if ( fy ) part->loc[1] = p->sizeY() - part->loc[1];
					if ( fz ) part->loc[2] = p->sizeZ() - part->loc[2];
				}
				for ( bad=mg->bad; bad; bad=bad->next ) {
					if ( fx ) bad->loc[0] = p->sizeX() - bad->loc[0];
					if ( fy ) bad->loc[1] = p->sizeY() - bad->loc[1];
					if ( fz ) bad->loc[2] = p->sizeZ() - bad->loc[2];
				}
				delete p;
			}
		}
	}

	return 0;
}

/**
@brief 	Gets views from particle images.
@param 	*project		project parameter structure.
@return int				0.
**/
int			project_set_views_from_images(Bproject* project)
{	
	if ( !project ) return 0;
	
	long				n;
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;
	Bimage* 		p = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Getting parameters from picked particle files" << endl << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->fpart.length() && ( p = read_img(mg->fpart, 0, -1) ) ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG project_set_views_from_images: File=" << mg->fpart << " particles=" << p->images() << endl;
				if ( p->images() < 2 ) mg->fmg = mg->fpart;
				if ( !mg->part ) 
					for ( n=0; n<p->images(); n++ ) particle_add(&mg->part, n+1);
				if ( mg->pixel_size[0] < 0.01 && p->sampling(0)[0] > 0 )
					mg->pixel_size = p->sampling(0);
				for ( part = mg->part; part; part = part->next ) {
					n = part->id - 1;
					if ( part->pixel_size[0] < 0.01 && p->sampling(0)[0] > 0 )
						part->pixel_size = p->sampling(0);
					part->ori = p->size()/2;
					if ( p->image[n].origin()[0] >= 0 && p->image[n].origin()[0] < p->sizeX() )
						part->ori[0] = p->image[n].origin()[0];
					if ( p->image[n].origin()[1] >= 0 && p->image[n].origin()[1] < p->sizeY() )
						part->ori[1] = p->image[n].origin()[1];
					if ( p->image[n].origin()[2] >= 0 && p->image[n].origin()[2] < p->sizeZ() )
						part->ori[2] = p->image[n].origin()[2];
					part->view = p->image[n].view();
				}
				delete p;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Sets views in particle images from the project parameter structure.
@param 	*project		project parameter structure.
@return int				0.
**/
int			project_set_views_in_images(Bproject* project)
{	
	if ( !project ) return 0;
	
	long				n;
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;
	Bimage* 		p = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Setting view parameters in picked particle files" << endl << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			p = read_img(mg->fpart, 1, -1);
			if ( p ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG project_set_views_in_images: File=" << mg->fpart << " particles=" << p->images() << endl;
				if ( p->images() < 2 ) mg->fmg = mg->fpart;
				if ( mg->pixel_size[0] > 0 )
					p->sampling(mg->pixel_size);
				if ( mg->part ) {
					for ( part = mg->part; part; part = part->next ) {
						n = part->id - 1;
						if ( part->pixel_size[0] > 0 )
							p->image[n].sampling(part->pixel_size);
						p->image[n].origin(p->size()/2);
						if ( part->ori[0] < 0 || part->ori[0] >= p->sizeX() )
							part->ori[0] = p->sizeX()/2;
						if ( part->ori[1] < 0 || part->ori[1] >= p->sizeY() )
							part->ori[1] = p->sizeY()/2;
						if ( part->ori[2] < 0 || part->ori[2] >= p->sizeZ() )
							part->ori[2] = p->sizeZ()/2;
						p->image[n].origin(part->ori);
						p->image[n].view(part->view);
					}
					write_img(mg->fpart, p, 0);
				}
				delete p;
			}
		}
	}
	
	return 0;
}

Bstring			mg_img_scale(Bstring& filename, Bstring& path, Bstring& insert, 
					double scale, Vector3<double> origin, Vector3<double>  pixel_size)
{
	if ( filename.length() < 1 ) return filename;
	
	Bimage*		p = read_img(filename, 1, -1);
	
	if ( !p ) {
		delete p;
		return filename;
	}

	if ( verbose )
		cout << "Binning " << filename << endl;
	
	p->sampling(pixel_size);
	
	Bstring		nufile(filename);
	
	if ( nufile.contains("/") ) nufile = nufile.post_rev('/');
	nufile = nufile.pre_rev('.') + insert + nufile.post_rev('.');
	if ( path.length() ) nufile = path + nufile;
	
	Vector3<long>		size(p->size() * scale);
	size = size.max(1);

	Vector3<double>		svec(scale, scale, 1);

	Vector3<double>		translate = (origin*scale - origin).round(0);

	Matrix3 			mat(1);
	
	Bimage*				ps = p->transform(size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
	
	write_img(nufile, ps, 0);
	
	delete p;
	delete ps;

	return nufile;
}

/**
@brief 	Changes all the micrograph fields linked to the box size and write new particle files.
@param 	*project		project parameter structure.
@param 	scale			scaling factor.
@param 	mgpath			scaled micrograph path.
@param 	partpath		scaled particle path.
@param 	insert			string to insert into new file names.
@return int 			0.

	The fields linked to pixel size are:
		pixel_size
		shift
		box, bad and marker radii
		particle origins and locations
		bad area and marker locations

**/
int			project_scale_box(Bproject* project, double scale, Bstring mgpath, Bstring partpath, Bstring insert)
{
	if ( !project ) return 0;
	
	if ( scale <= 0 ) return 0;
	
	Vector3<long>		size;
	Vector3<double>		shift, origin, translate;
	Vector3<double>		svec(scale, scale, scale);
	Matrix3 			mat(1);
	Bstring				newfile;
	Bimage* 			p = NULL;
	Bimage* 			ps = NULL;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Bfilament*			fil;
	Bfilnode*			fnode;
	Bbadarea*			bad;
	Bmarker*			mark;
	
	if ( verbose )
		cout << "Scaling by " << scale << endl;

	if ( mgpath == "." || mgpath == "./" ) mgpath = 0;
	if ( partpath == "." || partpath == "./" ) partpath = 0;
	
	if ( mgpath.length() ) {
		if ( mgpath[-1] != '/' ) mgpath += "/";
		if ( verbose )
			cout << "Micrograph path:                " << mgpath << endl;
	}
	
	if ( partpath.length() ) {
		if ( partpath[-1] != '/' ) partpath += "/";
		if ( verbose )
			cout << "Particle path:                  " << partpath << endl;
	}

	if ( verbose )
		cout << endl;
	
	if ( !project->select ) {
		svec[2] = 1;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG project_scale_box: scale=" << scale << endl;
/*				mg->pixel_size /= scale;
				origin = mg->origin;
				mg->origin = mg->origin * scale;
				p = NULL;
				if ( mg->fmg.length() ) p = read_img(mg->fmg, 1, -1);
				if ( p ) {
					newfile = mg->fmg;
					if ( newfile.contains("/") ) newfile = newfile.post_rev('/');
					newfile = newfile.pre_rev('.') + insert + newfile.post_rev('.');
					if ( mgpath.length() ) newfile = mgpath + newfile;
					mg->fmg = newfile;
					size = p->size() * scale;
					size = size.max(1);
					translate = (mg->origin - origin).round(0);
					ps = p->transform(size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
					ps->sampling(mg->pixel_size);
					write_img(mg->fmg, ps, 0);
					delete p;
					delete ps;
				}*/
				if ( mg->fmg.length() )
					mg->fmg = mg_img_scale(mg->fmg, mgpath, insert, scale, mg->origin, mg->pixel_size);
				if ( mg->fframe.length() )
					mg->fframe = mg_img_scale(mg->fframe, mgpath, insert, scale, mg->origin, mg->pixel_size);
				mg->origin = mg->origin * scale;
				mg->pixel_size /= scale;
				origin = mg->box_size/2;
				mg->box_size = mg->box_size * scale;
				mg->box_size = mg->box_size.max(1);
				mg->bad_radius = bround(mg->bad_radius*scale, 0);
				mg->mark_radius = bround(mg->mark_radius*scale, 0);
				for ( part = mg->part; part; part = part->next ) {
					part->loc *= scale;
					part->pixel_size /= scale;
					shift = part->loc.remainder(1);
					part->loc = part->loc.floor(0);
					part->ori = (part->ori * scale) + shift;
				}
				for ( fil=mg->fil; fil; fil=fil->next ) {
					for ( fnode=fil->node; fnode; fnode=fnode->next )
						fnode->loc = (fnode->loc * scale).round(0);
				}
				for ( bad=mg->bad; bad; bad=bad->next )
					bad->loc = (bad->loc * scale).round(0);
				for ( mark=mg->mark; mark; mark=mark->next )
					mark->loc = (mark->loc * scale).round(0);
/*				if ( mg->fpart.length() ) {
					p = read_img(mg->fpart, 1, -1);
					if ( p ) {
						newfile = mg->fpart;
						if ( newfile.contains("/") ) newfile = newfile.post_rev('/');
						newfile = newfile.pre_rev('.') + insert + newfile.post_rev('.');
						if ( partpath.length() ) newfile = partpath + newfile;
						mg->fpart = newfile;
						translate = mg->box_size/2 - origin;
						ps = p->transform(mg->box_size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
						ps->sampling(mg->pixel_size);
						write_img(mg->fpart, ps, 0);
						delete p;
						delete ps;
					}
				}*/
				if ( mg->fpart.length() )
					mg->fpart = mg_img_scale(mg->fpart, partpath, insert, scale, origin, mg->pixel_size);
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_scale_box: scale=" << scale << endl;
			rec->voxel_size /= scale;
			origin = rec->origin;
			rec->origin = rec->origin * scale;
			p = NULL;
			if ( rec->frec.length() ) p = read_img(rec->frec, 1, -1);
			if ( p ) {
				newfile = rec->frec;
				if ( newfile.contains("/") ) newfile = newfile.post_rev('/');
				newfile = newfile.pre_rev('.') + insert + newfile.post_rev('.');
				if ( mgpath.length() ) newfile = mgpath + newfile;
				rec->frec = newfile;
				size = p->size() * scale;
				translate = (rec->origin - origin).round(0);
//				img_transform(p, size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
				ps = p->transform(size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
				ps->sampling(rec->voxel_size);
				write_img(rec->frec, ps, 0);
				delete p;
				delete ps;
			}
			rec->box_size = rec->box_size * scale;
			rec->bad_radius = bround(rec->bad_radius*scale, 0);
			rec->mark_radius = bround(rec->mark_radius*scale, 0);
			for ( part=rec->part; part; part=part->next ) {
				part->loc *= scale;
				part->pixel_size /= scale;
				shift = part->loc.remainder(1);
				part->loc = part->loc.floor(0);
				part->ori = (part->ori * scale) + shift;
			}
			for ( fil=rec->fil; fil; fil=fil->next ) {
				for ( fnode=fil->node; fnode; fnode=fnode->next )
					fnode->loc = (fnode->loc * scale).round(0);
			}
			for ( bad=rec->bad; bad; bad=bad->next )
				bad->loc = (bad->loc * scale).round(0);
			for ( mark=rec->mark; mark; mark=mark->next )
				mark->loc = (mark->loc * scale).round(0);
			if ( rec->fpart.length() ) {
				p = read_img(rec->fpart, 1, -1);
				if ( p ) {
					newfile = rec->fpart;
					if ( newfile.contains("/") ) newfile = newfile.post_rev('/');
					newfile = newfile.pre_rev('.') + insert + newfile.post_rev('.');
					if ( partpath.length() ) newfile = partpath + newfile;
					rec->fpart = newfile;
//					translate = (rec->box_size/2 - origin).round(0);
					translate = rec->box_size/2 - origin;
//					img_transform(p, rec->box_size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
					ps = p->transform(rec->box_size, svec, origin, translate, mat, FILL_BACKGROUND, 0);
					ps->sampling(rec->voxel_size);
					write_img(rec->fpart, ps, 0);
					delete p;
					delete ps;
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Bins all micrographs a project.
@param 	*project	project parameter structure.
@param 	bin			binning value.
@param 	mgpath		binned micrograph path (must be allocated).
@param 	partpath	binned particle path (must be allocated).
@return int			0.

	All micrographs in a project are binned by the indicated value.
	New micrograph file names are generated with a "_b<n>" insert,
	where the n indicates the bin value.
	The path to the binned micrograph can be specified.

**/
int			project_bin_micrographs(Bproject* project, int bin, Bstring mgpath, Bstring partpath)
{
	if ( bin <= 1 ) return -1;
	
	if ( verbose )
		cout << "Binning micrographs by " << bin << endl;

	double		scale = 1.0/bin;
	Bstring		insert(bin, "_b%d.");
	
	project_scale_box(project, scale, mgpath, partpath, insert);
	
	return 0;
}

/**
@brief 	Changes all the micrograph fields linked to the pixel size.
@param 	*project		project parameter structure.
@param 	new_pixel_size	new pixel size.
@param 	mgpath			binned micrograph path (must be allocated).
@param 	partpath		binned particle path (must be allocated).
@return int 			0.

	The fields linked to pixel size are:
		pixel_size
		shift
		box, bad and marker radii
		particle origins and locations
		bad area and marker locations

**/
int			project_change_pixel_size(Bproject* project,
				Vector3<double> new_pixel_size,
				Bstring mgpath, Bstring partpath)
{
	if ( !project ) return 0;
	
	if ( new_pixel_size[0] <= 0 ) return 0;

	double			scale(0);
	
	if ( project->select ) {
		if ( project->rec )
			scale = project->rec->voxel_size[0]/new_pixel_size[0];
	} else {
		if ( project->field && project->field->mg )
			scale = project->field->mg->pixel_size[0]/new_pixel_size[0];
	}

	Bstring		insert(scale, "_s%3.1lf.");

	if ( scale ) project_scale_box(project, scale, mgpath, partpath, insert);
	
	return 0;
}


