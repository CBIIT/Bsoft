/**
@file	mg_processing.cpp
@brief	Functions for micrograph processing
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20210722
**/

#include "mg_processing.h"
#include "mg_select.h"
#include "mg_tags.h"
#include "qsort_functions.h"
#include "linked_list.h"
#include "math_util.h"
#include "random_numbers.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Creates a micrograph processing project parameter structure.
@param 	nmg				number of micrograph structures.
@param 	nrec			number of reconstruction structures.
@return Bproject*			project parameter structure.

	The function allocates memory for the project structure, a field-of-view structure,
	and the requested number of micrographs and reconstructions.

**/
Bproject*	project_create(int nmg, int nrec)
{
	int					i, j=0;
	Bstring				id("1");
	Bproject*			project = new Bproject;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	
	if ( nmg > 0 ) {
		field_add(&project->field, id);
		for ( i=0; i<nmg; i++, j++ ) {
			id = Bstring(i, "%03d");
			mg = micrograph_add(&project->field->mg, id);
			mg->block = j;
			mg->img_num = i;
		}
	}
	
	if ( nrec > 0 ) {
		for ( i=0; i<nrec; i++, j++ ) {
			id = Bstring(i, "%03d");
			rec = reconstruction_add(&project->rec, id);
			rec->block = j;
		}
	}
	
	if ( nmg < 1 && nrec > 0 ) project->select = 1;
	
	return project;
}

/**
@brief 	Sets micrograph and particle files to be the same.
@param 	*project	project parameter structure.
@return int 				0.

	If the micrograph file name exists, that is used, otherwise the
	particle file name is used.

**/
int			project_equal_mg_part_files(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	
	if ( verbose & VERB_PROCESS )
		cout << "Setting micrograph and particle file names to be the same" << endl << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
//			cout << "mg=" << mg->fmg << " part=" << mg->fpart << endl;
			if ( mg->fmg.length() ) mg->fpart = mg->fmg;
			else if ( mg->fpart.length() ) mg->fmg = mg->fpart;
			else cerr << "No file name for micrograph " << mg->id << endl;
//			cout << "mg=" << mg->fmg << " part=" << mg->fpart << endl;
		}
	}
	
	return 0;
}

/**
@brief 	Adds a field-of-view parameter structure to a linked list.
@param 	**field		pointer to any field-of-view in the list.
@param 	field_id	field identifier.
@return Bfield* 			new field-of-view.

	The function allocates memory for a new field-of-view structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bfield* 	field_add(Bfield** field, Bstring& field_id)
{
	Bfield* 		this_field = *field;
	Bfield* 		nu_field = new Bfield;
	
	if ( field_id.length() ) nu_field->id = field_id.c_str();
	
	if ( !this_field )
		*field = nu_field;
	else {
		while ( this_field->next ) this_field = this_field->next;
		this_field->next = nu_field;
	}
	
	return nu_field;
}

/**
@brief 	Adds a micrograph parameter structure to a linked list.
@param 	**mg	pointer to any micrograph in the list.
@param 	mg_id		micrograph identifier.
@return Bmicrograph* 		new micrograph.

	The function allocates memory for a new micrograph structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bmicrograph*	micrograph_add(Bmicrograph** mg, Bstring& mg_id)
{
	Bmicrograph*	this_mg = *mg;
	Bmicrograph*	nu_mg = new Bmicrograph;
	
	if ( mg_id.length() ) {
		if ( mg_id.contains("/") ) mg_id = mg_id.post_rev('/');
		nu_mg->id = mg_id;
	}
	
	if ( !this_mg )
		*mg = nu_mg;
	else {
		while ( this_mg->next ) this_mg = this_mg->next;
		this_mg->next = nu_mg;
	}
	
	return nu_mg;
}

/**
@brief 	Adds a movie frame parameter structure to a linked list.
@param 	**frame		pointer to any frame in the list.
@param 	pid 		particle number in file (starts at 1).
@return Bframe* 	new frame.

	The function allocates memory for a new frame structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bframe* 	frame_add(Bframe** frame, int pid)
{
	Bframe*			this_frame = *frame;
	Bframe*			nu_frame = new Bframe;
	
	nu_frame->id = pid;
	
	if ( !this_frame )
		*frame = nu_frame;
	else {
		while ( this_frame->next ) this_frame = this_frame->next;
		this_frame->next = nu_frame;
	}
	
	return nu_frame;
}

/**
@brief 	Adds a particle parameter structure to a linked list.
@param 	**part		pointer to any particle in the list.
@param 	pid 		particle number in file (starts at 1).
@return Bparticle* 	new particle.

	The function allocates memory for a new particle structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bparticle* 	particle_add(Bparticle** part, int pid)
{
	Bparticle*		this_part = *part;
	Bparticle*		nu_part = new Bparticle;
	
	nu_part->id = pid;
	
	if ( !this_part )
		*part = nu_part;
	else {
		while ( this_part->next ) this_part = this_part->next;
		this_part->next = nu_part;
	}
	
	return nu_part;
}

/**
@brief 	Sets the links from particles to micrographs.
@param 	*mg			pointer to micrograph.
@return int			0.
**/
int			mg_part_links(Bmicrograph* mg)
{
	Bparticle*			part;
	
//	cout << "Setting links from particles to micrographs" << endl;
	
	for ( part = mg->part; part; part = part->next ) part->mg = mg;
	
	return 0;
}
	
/**
@brief 	Sets the links from particles to reconstructions.
@param 	*rec		pointer to reconstruction.
@return int			0.
**/
int			rec_part_links(Breconstruction* rec)
{
	Bparticle*			part;
	
	for ( part = rec->part; part; part = part->next ) part->rec = rec;
	
	return 0;
}
	

/**
@brief 	Finds the first particle in a project.
@param 	*project	pointer to project.
@return Bparticle* 	first particle with a filename.

	The function searches for the first particle with a filename.

**/
Bparticle*	part_find_first(Bproject* project)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part = NULL;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				part = mg->part;
				if ( part ) return part;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			part = rec->part;
			if ( part ) return part;
		}
	}

	return part;
}

/**
@brief 	Adds a filament parameter structure to a linked list.
@param 	**fil			pointer to any filament in the list.
@param 	pid 				filament number in file (starts at 1).
@return Bfilament* 				new filament.

	The function allocates memory for a new filament structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bfilament* 	filament_add(Bfilament** fil, int pid)
{
	Bfilament*		this_fil = *fil;
	Bfilament*		nu_fil = new Bfilament;
	
	nu_fil->id = pid;
	
	if ( !this_fil )
		*fil = nu_fil;
	else {
		while ( this_fil->next ) this_fil = this_fil->next;
		this_fil->next = nu_fil;
	}
	
	return nu_fil;
}

/**
@brief 	Adds a filament node parameter structure to a linked list.
@param 	**fnode		pointer to any filament node in the list.
@param 	pid 				filament node number in file (starts at 1).
@return Bfilnode* 				new filament node.

	The function allocates memory for a new filament node structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bfilnode* 	filament_node_add(Bfilnode** fnode, int pid)
{
	Bfilnode*		this_fnode = *fnode;
	Bfilnode*		nu_fnode = new Bfilnode;
	
	nu_fnode->id = pid;
	
	if ( !this_fnode )
		*fnode = nu_fnode;
	else {
		while ( this_fnode->next ) this_fnode = this_fnode->next;
		this_fnode->next = nu_fnode;
	}
	
	return nu_fnode;
}

/**
@brief 	Calculates the length of a filament.
@param 	*fil		pointer to a filament.
@return double 				length in coordinate units.

	The length is defined as the sum of the link lengths connecting the nodes.

**/
double 		filament_length(Bfilament* fil)
{
	double 		length(0);
	Bfilnode*	fnode;
	
	for ( fnode=fil->node; fnode->next; fnode=fnode->next )
		length += (fnode->loc - fnode->next->loc).length();
	
	return length;
}

/**
@brief 	Adds a reconstruction parameter structure to a linked list.
@param 	**rec			pointer to any reconstruction in the list.
@param 	&rec_id			reconstruction identifier.
@return Breconstruction* 		new reconstruction.

	The function allocates memory for a new micrograph structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Breconstruction*	reconstruction_add(Breconstruction** rec, Bstring& rec_id)
{
	Breconstruction*	this_rec = *rec;
	Breconstruction*	nu_rec = new Breconstruction;
	
	if ( rec_id.length() ) nu_rec->id = rec_id;
	
	nu_rec->select = 1;

	if ( !this_rec )
		*rec = nu_rec;
	else {
		while ( this_rec->next ) this_rec = this_rec->next;
		this_rec->next = nu_rec;
	}
	
	return nu_rec;
}

/**
@brief 	Frees a whole project.

	Frees all the structures down the hierarchy.

@param 	*project	project parameter structure.
@return int 				0.
**/
int 		project_kill(Bproject* project)
{
	if ( !project ) return 0;
	
	project->comment = 0;
	project->filename = 0;
	
	string_kill(project->reference);
	
	Bfield			*field, *field2;
	Breconstruction	*rec, *rec2;
	
	for ( field = project->field; field; ) {
		field2 = field->next;
		field_kill(field);
		field = field2;
	}
		
	for ( rec = project->rec; rec; ) {
		rec2 = rec->next;
		reconstruction_kill(rec);
		rec = rec2;
	}
		
	delete project;
	
	return 0;
}

/**
@brief 	Frees a whole field-of-view structure.

	Frees all the structures down the hierarchy.

@param 	*field		field-of-view parameter structure.
@return int 				0.
**/
int 		field_kill(Bfield* field)
{
	if ( !field ) return 0;
	
	field->id = 0;
	
	Bmicrograph *mg, *mg2;
	
	for ( mg = field->mg; mg; ) {
		mg2 = mg->next;
		micrograph_kill(mg);
		mg = mg2;
	}
		
	delete field;
	
	return 0;
}

/**
@brief 	Frees a whole micrograph parameter structure.

	Frees all the structures down the hierarchy.

@param 	*mg		micrograph parameter structure.
@return int 				0.
**/
int 		micrograph_kill(Bmicrograph* mg)
{
	if ( !mg ) return 0;
	
	mg->id = 0;
	mg->fmg = 0;
	mg->fframe = 0;
	mg->fpart = 0;
	mg->ffil = 0;
	mg->fft = 0;
	mg->fps = 0;
	
	delete mg->ctf;
	
	particle_kill(mg->part);
	filament_kill(mg->fil);
	
	kill_list((char *) mg->bad, sizeof(Bbadarea));
	kill_list((char *) mg->sf, sizeof(Bstrucfac));
	kill_list((char *) mg->layer, sizeof(Blayerline));
	kill_list((char *) mg->mark, sizeof(Bmarker));
	
	delete mg;
	
	return 0;
}

/**
@brief 	Frees a list of particles.

	Frees all the structures down the hierarchy.

@param 	*part		particle parameter structure.
@return int 				0.
**/
int 		particle_kill(Bparticle* part)
{
	if ( !part ) return 0;
	
	Bparticle	*part2;
	
	for ( ; part; ) {
		part2 = part->next;
		part->fpart = 0;
		delete part;
		part = part2;
	}
	
	return 0;
}

/**
@brief 	Frees a list of filaments.

	Frees all the structures down the hierarchy.

@param 	*fil		filament parameter structure.
@return int 				0.
**/
int 		filament_kill(Bfilament* fil)
{
	if ( !fil ) return 0;
	
	Bfilament	*fil2;
	Bfilnode	*fnode, *fnode2;
	
	for ( ; fil; ) {
		for ( fnode=fil->node; fnode; ) {
			fnode2 = fnode->next;
			delete fnode;
			fnode = fnode2;
		}
		fil2 = fil->next;
		fil->ffil = 0;
		delete fil;
		fil = fil2;
	}
	
	return 0;
}

/**
@brief 	Frees a whole reconstruction parameter structure.

	Frees all the structures down the hierarchy.

@param 	*rec	reconstruction parameter structure.
@return int						0.
**/
int 		reconstruction_kill(Breconstruction* rec)
{
	if ( !rec ) return 0;
	
	rec->id = 0;
	rec->symmetry = 0;
	rec->frec = 0;
	rec->fpart = 0;
	rec->ffil = 0;
	rec->fft = 0;
	rec->fps = 0;
	
	particle_kill(rec->part);
	filament_kill(rec->fil);
	
	kill_list((char *) rec->bad, sizeof(Bbadarea));
	kill_list((char *) rec->sf, sizeof(Bstrucfac));
	kill_list((char *) rec->mark, sizeof(Bmarker));

	string_kill(rec->model);
	
	delete rec;
	
	return 0;
}

/**
@brief 	Adds one particle onto a list of particle parameters.
@param 	**partlist		destination particle list.
@param 	*part			particle.
@return Bparticle*		new particle.
**/
Bparticle*	particle_copy(Bparticle** partlist, Bparticle* part)
{
	int			i;
	Bparticle*	partnu = NULL;
	
	partnu = particle_add(partlist, part->id);
	partnu->fpart = part->fpart;		// Particle image file name
	partnu->group = part->group;		// Group membership (such as filaments or crystals)
	partnu->mag = part->mag;			// Magnification
	partnu->def = part->def;			// Defocus (angstrom)
	partnu->dev = part->dev;			// Defocus deviation (angstrom)
	partnu->ast = part->ast;			// Astigmatism angle (radians)
	partnu->loc = part->loc;			// Coordinates in the micrograph or tomogram
	partnu->pixel_size = part->pixel_size;	// Sampling in angstrom/pixel
	partnu->ori = part->ori;			// Origin of particle in voxel units
	partnu->view = part->view; 			// View: 3-value vector and angle (radians)
	for ( i=0; i<NFOM; i++ ) partnu->fom[i] = part->fom[i];		// Figure-of-merit, types defined in project structure
	partnu->sel = part->sel;			// Selection flag
	partnu->mg = part->mg;
	partnu->rec = part->rec;
	
	return partnu;
}
	
/**
@brief 	Copies a list of particle parameters to a new list.
@param 	*partlist		source particle list.
@return Bparticle*		new list.
**/
Bparticle*	particle_copy(Bparticle* partlist)
{
	Bparticle*	part;
	Bparticle*	partlistnu = NULL;
	
	for ( part = partlist; part; part = part->next )
		particle_copy(&partlistnu, part);
	
	return partlistnu;
}
	
Bbadarea*	particle_bad_copy(Bbadarea* badlist)
{
	Bbadarea*	bad;
	Bbadarea*	badnu = NULL;
	Bbadarea*	badlistnu = NULL;
	
	for ( bad = badlist; bad; bad = bad->next ) {
		badnu = (Bbadarea *) add_item((char **) &badlistnu, sizeof(Bbadarea));
		badnu->id = bad->id;
		badnu->loc = bad->loc;
	}
	
	return badlistnu;
}
	
Bfilament*	filament_copy(Bfilament* fillist)
{
	Bfilament*	fil;
	Bfilament*	filnu = NULL;
	Bfilament*	fillistnu = NULL;
	Bfilnode*	fnode = NULL;
	Bfilnode*	fnodenu = NULL;
	
	for ( fil = fillist; fil; fil = fil->next ) {
		filnu = filament_add(&fillistnu, fil->id);
		filnu->ffil = fil->ffil;
		for ( fnode = fil->node; fnode; fnode = fnode->next ) {
			fnodenu = filament_node_add(&filnu->node, fnode->id);
			fnodenu->loc = fnode->loc;
		}
	}
	
	return fillistnu;
}
	
Bstrucfac*	structurefactor_copy(Bstrucfac* sflist)
{
	Bstrucfac*	sf;
	Bstrucfac*	sfnu = NULL;
	Bstrucfac*	sflistnu = NULL;
	
	for ( sf = sflist; sf; sf = sf->next ) {
		sfnu = (Bstrucfac *) add_item((char **) &sflistnu, sizeof(Bstrucfac));
		sfnu->loc = sf->loc;			// Coordinates in the transform
		sfnu->index = sf->index;		// Miller indices
		sfnu->amp = sf->amp;			// Amplitude
		sfnu->sigamp = sf->sigamp;		// Amplitude deviation
		sfnu->phi = sf->phi;			// Phase
		sfnu->sigphi = sf->sigphi;		// Phase deviation
		sfnu->fom = sf->fom;			// Figure-of-merit
		sfnu->sel = sf->sel;			// Selection flag
	}
	
	return sflistnu;
}
	
Blayerline*	layerline_copy(Blayerline* layerlist)
{
	Blayerline*	line;
	Blayerline*	linenu = NULL;
	Blayerline*	linelistnu = NULL;
	
	for ( line = layerlist; line; line = line->next ) {
		linenu = (Blayerline *) add_item((char **) &linelistnu, sizeof(Blayerline));
		linenu->number = line->number;		// Layer line number
		linenu->order = line->order;		// Bessel order
		linenu->distance = line->distance;	// Distance along helical axis from origin
		linenu->freq = line->freq;			// Spatial frequency
		linenu->amp = line->amp;			// Amplitude
		linenu->fom = line->fom;			// Figure-of-merit
		linenu->sel = line->sel;			// Selection flag
	}
	
	return linelistnu;
}

/**
@brief 	Copies a micrograph.
@param 	*mg			micrograph structure to be copied.
@return Bmicrograph* 			new micrograph.
**/
Bmicrograph*	micrograph_copy(Bmicrograph* mg)
{
	Bmicrograph*	mgnu = NULL;
	
	if ( verbose & VERB_FULL )
		cout << "Copying micrograph " << mg->id << endl;

	micrograph_add(&mgnu, mg->id);

	mgnu->select = mg->select;
	mgnu->fmg = mg->fmg;					// Micrograph image file
	mgnu->fframe = mg->fframe;				// Micrograph frames image file
	mgnu->fpart = mg->fpart;				// Image file with picked particles
	mgnu->ffil = mg->ffil;					// Image file with filaments
	mgnu->fft = mg->fft;					// Image file with Fourier transform
	mgnu->fps = mg->fps;					// Image file with power spectrum
	mgnu->img_num = mg->img_num;			// Image number in file
	mgnu->magnification = mg->magnification;	// Microscope magnification
	mgnu->sampling = mg->sampling;			// Scanner sampling (angstrom)
	mgnu->pixel_size = mg->pixel_size; 		// Nominal micrograph pixel size
	mgnu->exposure = mg->exposure; 			// Acquisition time
	mgnu->dose = mg->dose; 					// Direct beam/Dose
	mgnu->intensity = mg->intensity; 		// Micrograph average
	mgnu->tilt_axis = mg->tilt_axis;		// Tilt axis angle, origin at x-axis (radians)
	mgnu->tilt_angle = mg->tilt_angle; 		// Tilt angle, right-handed around tilt axis (radians)
	mgnu->level_angle = mg->level_angle; 	// Level angle (radians)
	mgnu->rot_angle = mg->rot_angle;		// In-plane rotation angle of micrograph or specimen (radians)
	mgnu->origin = mg->origin;				// Origin within field-of-view
	mgnu->scale = mg->scale;				// Scale with respect to field-of-view
	mgnu->matrix = mg->matrix;				// Rotation matrix
	mgnu->box_size = mg->box_size;			// Particle size
	mgnu->bad_radius = mg->bad_radius;		// Radius of bad area around coordinates
	mgnu->filament_width = mg->filament_width;	// Filament width
	mgnu->fil_node_radius = mg->fil_node_radius;// Radius of filament node
	mgnu->sf_radius = mg->sf_radius;		// Radius of marker
	mgnu->mark_radius = mg->mark_radius;	// Radius of marker
	mgnu->fom = mg->fom;					// Micrograph figure-of-merit

	if ( mg->ctf ) mgnu->ctf->update(mg->ctf);
	
	mgnu->part = particle_copy(mg->part);
	
	mgnu->bad = particle_bad_copy(mg->bad);
	
	mgnu->fil = filament_copy(mg->fil);
	
	mgnu->sf = structurefactor_copy(mg->sf);
	
	mgnu->layer = layerline_copy(mg->layer);
	
	mgnu->mark = markers_copy(mg->mark);
	
	return mgnu;
}


/**
@brief 	Updates an existing project with new information.
@param 	*project		project parameter structure to be updated.
@param 	*proj_new		project parameter structure with new information.
@param 	fom_index			index of FOM to select on.
@return int 				0.


	The new information is encoded as a project hierarchy.
	The original project hierarchy is searched for the micrograph ID
	tags that correspond to the original project hierarchy and those
	micrograph parameters are updated. Because the micrograph ID tags
	should be unique in themselves, the field ID tags are ignored and
	the two project hierarchies may have different field designations.
	If a micrograph ID from the new project hierarchy is not found in
	the original project hierarchy, it is added.
	The new project hierarchy is not modified.

**/
int			project_update(Bproject* project, Bproject* proj_new, int fom_index)
{
	int				notfound;
	Bfield 			*field, *nu_field;
	Bmicrograph		*mg = NULL, *nu_mg;
	Breconstruction	*rec = NULL, *nu_rec;
	
	for ( nu_field = proj_new->field; nu_field; nu_field = nu_field->next ) {
		for ( nu_mg = nu_field->mg; nu_mg; nu_mg = nu_mg->next ) {
			notfound = 1;
			for ( field = project->field; field && notfound; field = field->next ) {
				for ( mg = field->mg; mg && notfound; ) {
					if ( mg->id == nu_mg->id )
						notfound = 0;
					else
						mg = mg->next;
				}
			}
			if ( notfound )
				cerr << "Error: Micrograph " << nu_mg->id << " not found!" << endl;
			else
				micrograph_update(mg, nu_mg, fom_index, 63);
		}
	}
	
	for ( nu_rec = proj_new->rec; nu_rec; nu_rec = nu_rec->next ) {
		notfound = 1;
		for ( rec = project->rec; rec && notfound; ) {
			if ( rec->id == nu_rec->id )
				notfound = 0;
			else
				rec = rec->next;
		}
		if ( notfound )
			cerr << "Error: Reconstruction " << nu_rec->id << " not found!" << endl;
		else
			reconstruction_update(rec, nu_rec, fom_index, 63);
	}
	
	return 0;
}

/**
@brief 	Updates an existing micrograph with new information.
@param 	*mg		micrograph structure to be updated.
@param 	*nu_mg	micrograph structure with new information.
@param 	fom_index		index of FOM to select on.
@param 	flags			flags indicating which parts to update.
@return int 				0.

	The new information is encoded as a micrograph structure.
	All fields are updated if the new fields contain non-default data.
	Defaults are defined here as zero-length strings or zeroes.
	The particles in the original micrograph are matched by ID with 
	those in the new micrograph and updated, with any new particles 
	added from the new micrograph.
	The bad area and marker coordinates are replaced if they exist 
	in the new structure.
	The new micrograph structure is not modified, except for deletion 
	of bad areas and markers.

	Flags:
		1	particles
		2	bad areas
		4	filaments
		8	structure factors
		16	layer lines
		32	markers
**/
int			micrograph_update(Bmicrograph* mg, Bmicrograph* nu_mg, int fom_index, int flags)
{
	if ( nu_mg->select < 1 ) return 0;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG micrograph_update: " << nu_mg->id << tab << nu_mg->select << tab << flags << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Updating micrograph " << mg->id << endl;
	
	mg->select =  nu_mg->select;
	if ( nu_mg->fmg.length() ) mg->fmg = nu_mg->fmg;			// Micrograph image file
	if ( nu_mg->fframe.length() ) mg->fframe = nu_mg->fframe;	// Micrograph frames image file
	if ( nu_mg->fpart.length() ) mg->fpart = nu_mg->fpart;		// Image file with picked particles
	if ( nu_mg->ffil.length() ) mg->ffil = nu_mg->ffil;			// Image file with filaments
	if ( nu_mg->fft.length() ) mg->fft = nu_mg->fft;			// Image file with Fourier transform
	if ( nu_mg->fps.length() ) mg->fps = nu_mg->fps;			// Image file with power spectrum
	if ( nu_mg->img_num ) mg->img_num = nu_mg->img_num;			// Image number in file
	if ( nu_mg->magnification ) mg->magnification = nu_mg->magnification;	// Microscope magnification
	if ( nu_mg->sampling ) mg->sampling = nu_mg->sampling;		// Scanner sampling (angstrom)
	if ( nu_mg->pixel_size[0] ) mg->pixel_size = nu_mg->pixel_size; // Nominal micrograph pixel size
	if ( nu_mg->exposure ) mg->dose = nu_mg->exposure; 			// Acquisition time
	if ( nu_mg->dose ) mg->dose = nu_mg->dose; 					// Direct beam/Dose
	if ( nu_mg->intensity ) mg->intensity = nu_mg->intensity; 	// Micrograph average
	if ( nu_mg->tilt_axis ) mg->tilt_axis = nu_mg->tilt_axis;	// Tilt axis angle, origin at x-axis (radians)
	if ( nu_mg->tilt_angle ) mg->tilt_angle = nu_mg->tilt_angle; // Tilt angle, right-handed around tilt axis (radians)
	if ( nu_mg->rot_angle ) mg->rot_angle = nu_mg->rot_angle;	// In-plane rotation angle of micrograph or specimen (radians)
	if ( nu_mg->origin.length2() ) mg->origin = nu_mg->origin;	// Origin within field-of-view
	if ( nu_mg->scale.length2() ) mg->scale = nu_mg->scale;		// Scale with respect to field-of-view
	if ( nu_mg->box_size[0] ) mg->box_size = nu_mg->box_size;	// Particle size
	if ( nu_mg->bad_radius ) mg->bad_radius = nu_mg->bad_radius;	// Radius of bad area around coordinates
	if ( nu_mg->mark_radius ) mg->mark_radius = nu_mg->mark_radius;	// Radius of marker
	
	if ( nu_mg->ctf ) {
		if ( !mg->ctf ) mg->ctf = new CTFparam;
		mg->ctf->update(nu_mg->ctf);
	}
	
	if ( flags & 1)
		particle_update(&mg->part, nu_mg->part, fom_index);
	
	if ( (flags & 2) && nu_mg->bad ) {
		if ( verbose & VERB_PROCESS )
			cout << "Replacing bad areas" << endl;
		kill_list((char *) mg->bad, sizeof(Bbadarea));
		mg->bad = nu_mg->bad;
		nu_mg->bad = NULL;
	}
	
	if ( (flags & 4) && nu_mg->fil ) {
		if ( verbose & VERB_PROCESS )
			cout << "Replacing filaments" << endl;
		kill_list((char *) mg->fil, sizeof(Bfilament));
		mg->fil = nu_mg->fil;
		nu_mg->fil = NULL;
	}
	
	if ( (flags & 8) && nu_mg->sf ) {
		if ( verbose & VERB_PROCESS )
			cout << "Replacing structure factors" << endl;
		kill_list((char *) mg->sf, sizeof(Bstrucfac));
		mg->sf = nu_mg->sf;
		nu_mg->sf = NULL;
	}
	
	if ( (flags & 16) && nu_mg->layer ) {
		if ( verbose & VERB_PROCESS )
			cout << "Replacing layer lines" << endl;
		kill_list((char *) mg->layer, sizeof(Blayerline));
		mg->layer = nu_mg->layer;
		nu_mg->layer = NULL;
	}
	
	if ( (flags & 32) && nu_mg->mark ) {
		if ( verbose & VERB_FULL )
			cout << "Replacing markers" << endl;
		kill_list((char *) mg->mark, sizeof(Bmarker));
		mg->mark = nu_mg->mark;
		nu_mg->mark = NULL;
	}
	
	return 0;
}

/**
@brief 	Updates an existing reconstruction with new information.
@param 	*rec			reconstruction structure to be updated.
@param 	*nu_rec			reconstruction structure with new information.
@param 	fom_index		index of FOM to select on.
@param 	flags			flags indicating which parts to update.
@return int 			0.

	The new information is encoded as a reconstruction structure.
	The update is only done if the FOM of the new reconstruction is better.
	All fields are updated if the new fields contain non-default data.
	Defaults are defined here as zero-length strings or zeroes.
	The particles in the original reconstruction are matched by ID with 
	those in the new reconstruction and updated, with any new particles 
	added from the new reconstruction.
	The bad area and marker coordinates are replaced if they exist 
	in the new structure.
	The new reconstruction structure is not modified, except for deletion 
	of bad areas and markers.

**/
int			reconstruction_update(Breconstruction* rec, Breconstruction* nu_rec, int fom_index, int flags)
{
	if ( nu_rec->select < 1 ) return 0;
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reconstruction_update: " << nu_rec->id << tab << nu_rec->select << tab << flags << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Updating reconstruction " << rec->id << endl;
	
	rec->select = nu_rec->select;
	rec->fom = nu_rec->fom;
	if ( nu_rec->symmetry.length() ) rec->symmetry = nu_rec->symmetry;	// Reconstruction symmetry
	if ( nu_rec->frec.length() ) rec->frec = nu_rec->frec;				// Reconstruction image file
	if ( nu_rec->fpart.length() ) rec->fpart = nu_rec->fpart;				// Image file with picked particles
	if ( nu_rec->ffil.length() ) rec->ffil = nu_rec->ffil;				// Image file with filaments
	if ( nu_rec->fft.length() ) rec->fft = nu_rec->fft;					// Image file with Fourier transform
	if ( nu_rec->fps.length() ) rec->fps = nu_rec->fps;					// Image file with power spectrum
	if ( nu_rec->voxel_size[0] ) rec->voxel_size = nu_rec->voxel_size;		// Nominal voxel size
	if ( nu_rec->origin.length2() ) rec->origin = nu_rec->origin;			// Origin within field-of-view
	if ( nu_rec->scale.length2() ) rec->scale = nu_rec->scale;			// Scale with respect to field-of-view
	if ( nu_rec->box_size[0] ) rec->box_size = nu_rec->box_size;			// Particle size
	if ( nu_rec->bad_radius ) rec->bad_radius = nu_rec->bad_radius;		// Radius of bad area around coordinates
	if ( nu_rec->mark_radius ) rec->mark_radius = nu_rec->mark_radius;	// Radius of marker
	rec->view = nu_rec->view;												// Reconstruction orientation
	
	if ( flags & 1)
		particle_update(&rec->part, nu_rec->part, fom_index);
	
	return 0;
}

/**
@brief 	Updates one existing particle with new information.
@param 	*part		particle structure to be updated.
@param 	*nu_part	particle structure with new information.
@return int			0.

	The new information is encoded as a particle structure.
	All fields are updated if the new fields contain non-default data.
	Defaults are defined here as zeroes.
	The new particle structure is not modified.

**/
int			particle_update(Bparticle* part, Bparticle* nu_part)
{
	if ( nu_part->def > 0 ) {
		part->def = nu_part->def;
		part->dev = nu_part->dev;
		part->ast = nu_part->ast;
	}
	
	if ( nu_part->mag > 0 && nu_part->mag != 1 )
		part->mag = nu_part->mag;
		
	if ( nu_part->loc.length2() > 0 )
		part->loc = nu_part->loc;
		
	if ( nu_part->pixel_size.volume() > 0 )
		part->pixel_size = nu_part->pixel_size;
		
	if ( nu_part->ori.length2() > 0 )
		part->ori = nu_part->ori;
		
	if ( nu_part->view.vector_size() > 0.1 )
		part->view = nu_part->view;
		
	for ( int f=0; f<NFOM; f++ )
		if ( nu_part->fom[f] ) part->fom[f] = nu_part->fom[f];
	
	part->sel = nu_part->sel;
	
	if ( nu_part->fpart.length() )
		part->fpart = nu_part->fpart;
	
	return 0;
}

/**
@brief 	Updates an existing particle list with new information.
@param 	**pnt_part	particle structure list to be updated.
@param 	*nu_part	particle structure list with new information.
@param 	fom_index	index of FOM to select on.
@return int			0.

	The new information is encoded as a particle structure.
	All fields are updated if the new fields contain non-default data.
	Defaults are defined here as zeroes.
	The new particle structure is not modified.

**/
int			particle_update(Bparticle** pnt_part, Bparticle* nu_part, int fom_index)
{
	if ( verbose & VERB_FULL )
		cout << "Updating particles" << endl;
	
	int				notfound(1), set_new(0);
	Bparticle*		part;

	for ( ; nu_part; nu_part = nu_part->next ) {
		notfound = 1;
		set_new = 0;
		for ( part = *pnt_part; part && notfound; ) {
			if ( part->id == nu_part->id ) notfound = 0;
			else part = part->next;
		}
		if ( notfound ) {
			cerr << "Particle " << nu_part->id << " not found!" << endl;
			part = particle_add(pnt_part, nu_part->id);
			part->group = 0;
			part->sel = 0;
			part->fom[fom_index] = 0;
			set_new = 1;
		}
		if ( nu_part->sel > 0 ) {
			set_new = 1;
			if ( part->sel > 0 ) {
				if ( nu_part->fom[fom_index] < part->fom[fom_index] ) set_new = 0;
//				if ( nu_part->fom[fom_index] > 0.999 && nu_part->view[2] > 0.999 ) set_new = 0;
			}
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG particle_update: " << part->id << tab << part->sel << tab << nu_part->sel << tab 
				<< part->fom[fom_index] << tab << nu_part->fom[fom_index] << tab << set_new << endl;
		if ( set_new ) particle_update(part, nu_part);
	}
	
	return 0;
}

/**
@brief 	Puts particle parameters from one project into another.
@param 	*project		project parameter structure with all parameters.
@param 	*partproject	project parameter structure with particle parameters.
@return int				0.
**/
int			project_merge_part_parameters(Bproject* project, Bproject* partproject)
{
	if ( !project || ! partproject ) return 0;
	
	Bfield				*field, *partfield;
	Bmicrograph			*mg, *partmg;
	Breconstruction		*rec, *partrec;

	if ( verbose & VERB_PROCESS )
		cout << "Merging particle parameters into the main project" << endl << endl;
	
	for ( field=project->field, partfield=partproject->field; field; field=field->next, partfield=partfield->next ) {
		for ( mg=field->mg, partmg=partfield->mg; mg && partmg; mg=mg->next, partmg=partmg->next ) {
			particle_update(&mg->part, partmg->part, 0);
		}
	}
	
	for ( rec=project->rec, partrec=partproject->rec; rec && partrec; rec=rec->next, partrec=partrec->next ) {
		particle_update(&rec->part, partrec->part, 0);
	}
	
	return 0;
}


/**
@brief 	Sets the links from particles back to the micrographs and reconstructions.
@param 	*project	project parameter structure.
@return long		number of particles.
**/
long		project_set_part_links(Bproject* project)
{
	long				npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			npart += micrograph_set_part_links(mg);
	
	for ( rec = project->rec; rec; rec = rec->next )
		npart += reconstruction_set_part_links(rec);
	
	return npart;
}

/**
@brief 	Sets the links from particles back to the micrograph.
@param 	*mg			micrograph parameter structure.
@return long		number of particles.
**/
long		micrograph_set_part_links(Bmicrograph* mg)
{
	long			npart(0);
	Bparticle*		part;
	
	for ( part = mg->part; part; part = part->next, npart++ ) part->mg = mg;
	
	return npart;
}

/**
@brief 	Sets the links from particles back to the reconstruction.
@param 	*rec		reconstruction parameter structure.
@return long		number of particles.
**/
long		reconstruction_set_part_links(Breconstruction* rec)
{
	long			npart(0);
	Bparticle*		part;
	
	for ( part = rec->part; part; part = part->next, npart++ ) part->rec = rec;
	
	return npart;
}

/**
@brief 	Divides a project into a number of micrograph subsets.
@param 	*project	project parameter structure.
@param 	n			number of subsets.
@return long		number of new projects created.

	Only selected micrographs are considered in the subdivision.

**/
long		project_divide(Bproject* project, long n)
{
	long				i, j, k, t, nmg(0), nrec(0), ndiv(0);
	Bproject*			project2 = project;
	Bfield*				field, *field_prev;
	Bmicrograph*		mg;
	Breconstruction*	rec, *rec_prev;

	for ( nmg=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select)
			nmg++;

	for ( nrec=0, rec = project->rec; rec; rec = rec->next )
		if ( rec->select ) nrec++;
	
	field = field_prev = project->field;
	rec = rec_prev = project->rec;

	for ( i=j=k=0; i<n; i++ ) {
		t = (nmg*(i+1))/n;
		for ( ; field && j<t; field_prev = field, field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				if ( mg->select) j++;
		t = (nrec*(i+1))/n;
		for ( ; rec && k<t; rec_prev = rec, rec = rec->next )
			if ( rec->select) k++;
		if ( field || rec ) {
			if ( field_prev ) field_prev->next = NULL;
			if ( rec_prev ) rec_prev->next = NULL;
			project2->next = new Bproject;
			project2 = project2->next;
			project2->field = field;
			project2->rec = rec;
			for ( i=0; i<NFOM; i++ ) project2->fom_tag[i] = project->fom_tag[i];
			ndiv++;
		}
	}
	
	return ndiv;
}

/**
@brief 	Sets up an array of selected micrograph pointers.
@param 	*project		project parameter structure.
@param 	&nmg			array length (modified).
@return Bmicrograph**	micrograph pointer array.

	FOM's of non-selected micrographs are zeroed.

**/
Bmicrograph**	project_micrograph_array(Bproject* project, long &nmg)
{
	Bfield*			field;
	Bmicrograph*	mg;

	for ( nmg=0, field=project->field; field; field=field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select )
			nmg++;

	Bmicrograph**	mgarr = new Bmicrograph*[nmg];
	
	for ( nmg=0, field=project->field; field; field=field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select )
			mgarr[nmg++] = mg;

	return mgarr;
}

/**
@brief 	Sets up an array of selected particle pointers.
@param 	*project		project parameter structure.
@param 	part_select		particle selection (-1 if none).
@param 	&npart			array length (modified).
@return Bparticle**		particle pointer array.

	FOM's of non-selected particles are zeroed.

**/
Bparticle**	project_mg_particle_array(Bproject* project, int part_select, long &npart)
{
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;

	for ( npart=0, field=project->field; field; field=field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				part->mg = mg;
				if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
					npart++;
				else
					part->fom[0] = part->fom[1] = part->sel = 0;
			}
		}
	}

	Bparticle**		partarr = new Bparticle*[npart];
	
	for ( npart=0, field=project->field; field; field=field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
					partarr[npart++] = part;
			}
		}
	}

	return partarr;
}

/**
@brief 	Sets up an array of selected particle pointers.
@param 	*project		project parameter structure.
@param 	part_select		particle selection (-1 if none).
@param 	&npart			array length (modified).
@return Bparticle**		particle pointer array.

	FOM's of non-selected particles are zeroed.

**/
Bparticle**	project_rec_particle_array(Bproject* project, int part_select, long &npart)
{
	Breconstruction*	rec;
	Bparticle*			part;

	for ( npart=0, rec=project->rec; rec; rec=rec->next ) {
		for ( part = rec->part; part; part = part->next ) {
			if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
				npart++;
			else
				part->fom[0] = part->fom[1] = part->sel = 0;
		}
	}

	Bparticle**		partarr = new Bparticle*[npart];
	
	for ( npart=0, rec=project->rec; rec; rec=rec->next ) {
		for ( part = rec->part; part; part = part->next ) {
			if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
				partarr[npart++] = part;
		}
	}

	return partarr;
}

/**
@brief 	Sets up an array of selected particle pointers.
@param 	*partlist		particle linked list.
@param 	part_select		particle selection (-1 if none).
@param 	&npart			array length (modified).
@return Bparticle**		particle pointer array.

	FOM's of non-selected particles are zeroed.

**/
Bparticle**	particle_array(Bparticle* partlist, int part_select, long &npart)
{
	Bparticle*			part;

	for ( npart=0, part = partlist; part; part = part->next ) {
		if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
			npart++;
		else
			part->fom[0] = part->fom[1] = part->sel = 0;
	}

	Bparticle**		partarr = new Bparticle*[npart];
	
	for ( npart=0, part = partlist; part; part = part->next )
		if ( ( part_select<0 && part->sel ) || ( part->sel==part_select ) )
			partarr[npart++] = part;

	return partarr;
}


/**
@brief 	Reverts file names to what it was in an old project.
@param 	*project		project parameter structure.
@param 	*project_old	project parameter structure with old file names.
@param	flag			1=mg, 2=ps, 4=rec, 8=part, 16=frames.
@return int				0.
**/
int			project_revert_filenames(Bproject* project, Bproject* project_old, int flag)
{
	if ( !project ) return 0;
	if ( !project_old ) return 0;
	
	Bfield*				field, *field_old;
	Bmicrograph*		mg, *mg_old;
	Breconstruction*	rec, *rec_old;
	Bparticle*			part, *part_old;

	if ( verbose & VERB_PROCESS )
		cout << "Reverting to old particle file names in: " << project_old->filename << endl << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( field_old = project_old->field; field_old && field_old->id != field->id; field_old = field_old->next ) ;
		if ( field_old ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( mg_old = field_old->mg; mg_old && mg_old->id != mg->id; mg_old = mg_old->next ) ;
				if ( mg_old ) {
					if ( flag & 1 ) mg->fmg = mg_old->fmg;
					if ( flag & 16 ) mg->fframe = mg_old->fframe;
					if ( flag & 2 ) mg->fps = mg_old->fps;
					if ( flag & 8 ) {
						mg->fpart = mg_old->fpart;
						part_old = mg_old->part;
						for ( part = mg->part; part; part = part->next ) {
							for ( ; part_old && part_old->id != part->id; part_old = part_old->next ) ;
							if ( part_old ) {
								part->fpart = part_old->fpart;
							}
						}
					}
				}
			}
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		for ( rec_old = project_old->rec; rec_old && rec_old->id != rec->id; rec_old = rec_old->next ) ;
		if ( rec_old ) {
			if ( flag & 4 ) rec->frec = rec_old->frec;
			if ( flag & 8 ) {
				part_old = rec_old->part;
				for ( part = rec->part; part; part = part->next ) {
					for ( ; part_old && part_old->id != part->id; part_old = part_old->next ) ;
					if ( part_old ) {
						part->fpart = part_old->fpart;
					}
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Sets the field id's according to a regular series.

	The micrograph data blocks must be arranged in a regular order with
	every nseries micrographs from the same field-of-view.
	The field id's are taken from the first micrograph id in a series.

@param 	*project	project parameter structure.
@param 	nseries			number of micrographs per field-of-view.
@param 	&field_id	a user-specified field ID.
@return int 				0.
**/
int			project_set_field_id(Bproject* project, int nseries, Bstring& field_id)
{
	int				i(0), nf(0);
	Bfield* 		field;
	Bfield* 		nu_field = NULL;
	Bfield* 		nu_field_start = NULL;
	Bmicrograph*	mg = NULL;
	Bmicrograph*	nu_mg = NULL;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( i%nseries == 0 ) {
				nu_field = field_add(&nu_field, mg->id);
				if ( i == 0 ) nu_field_start = nu_field;
				else nu_mg->next = NULL;
				nu_field->mg = nu_mg = mg;
				nf++;
			} else {
				nu_mg->next = mg;
				nu_mg = mg;
			}
			i++;
		}
	}
	
	for ( field=project->field; field; ) {
		nu_field = field->next;
		delete field;
		field = nu_field;
	}
	
	project->field = nu_field_start;
	
	if ( field_id.length() > 0 ) {
		if ( nf < 2 ) {
			project->field->id = field_id;
		} else {
			for ( i=1, field = project->field; field; field = field->next, i++ )
				field->id = field_id + Bstring(i, "_%04d");
		}
	}
	
	return 0;
}

/**
@brief 	Shows the project hierarchy.
@param 	*project	project parameter structure.
@return long 		number of particles selected.
**/
long		project_show_hierarchy(Bproject* project)
{
	if ( !project ) return 0;
	
	if ( !verbose ) return 0;
	
	long				i, j, k, l, nf(0), nmg(0), nmgsel(0);
	long				npart(0), nsel(0), nfil(0), nfn(0);
	double				nmgpercent(0), npercent(0);
	int					g;
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;

	long				gmax = 10000000;
	int*				ng = new int[gmax];
	for ( i=0; i<gmax; i++ ) ng[i] = 0;

	if ( verbose > VERB_RESULT )
		cout << "Project hierarchy:" << endl;
	if ( field ) {
		for ( field = project->field; field; field = field->next, nf++ ) {
			if ( verbose > VERB_RESULT )
				cout << "Field: " << field->id << " (" << field->select << ")" << endl;
			for ( mg = field->mg; mg; mg = mg->next, nmg++ ) {
				if ( mg->select ) nmgsel++;
				for ( part = mg->part; part; part = part->next ) {
					g = part->group;
					if ( g >= gmax ) g = gmax-1;
					ng[g]++;
				}
				i = particle_count(mg->part);
				j = particle_count_selected(mg->part);
				k = filament_count(mg->fil);
				l = filament_node_count(mg->fil);
				npart += i;
				nsel += j;
				nfil += k;
				nfn += l;
				if ( verbose > VERB_RESULT )
					cout << "\tMicrograph: " << mg->id << " ("
						<< mg->select << ")\tParticles: " << i << " ("
						<< j << ")\tFilaments: " << k << " (" << l << ")" << endl;
			}
		}
		if ( nmg ) nmgpercent = nmgsel*100.0/nmg;
		if ( npart ) npercent = nsel*100.0/npart;
		cout << "Totals:" << endl; 
		cout << nf << " fields\t" << nmg << " micrographs\t" 
			<< nmgsel << " selected\t(" << nmgpercent << "%)" << endl;
		if ( npart )
			cout << npart << " particles\t" << nsel << " selected\t(" 
				<< npercent << "%)" << endl;
		if ( nfil )
			cout << nfil << " filaments\t" << nfn << " filament nodes" << endl;
		if ( npart ) {
			cout << "Group\tCount" << endl;
			for ( g=0; g<gmax; g++ ) if ( ng[g] )
				cout << g << tab << ng[g] << endl;
		}
	}
	
	for ( g=0; g<gmax; g++ ) ng[g] = 0;
	if ( rec ) {
		npart = nsel = nfil = nfn = 0;
		npercent = 0;
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				g = part->group;
				if ( g >= gmax ) g = gmax-1;
				ng[g]++;
			}
			if ( verbose > VERB_RESULT ) {
				cout << "Reconstruction: " << rec->id;
				if ( rec->frec.length() ) cout << ", " << rec->frec;
				cout << endl;
			}
			i = particle_count(rec->part);
			j = particle_count_selected(rec->part);
			k = filament_count(rec->fil);
			l = filament_node_count(rec->fil);
			npart += i;
			nsel += j;
			nfil += k;
			nfn += l;
			if ( verbose > VERB_RESULT )
				cout << "\tParticles: " << i << " (" << j << ")\tFilaments: " << k << " (" << l << ")" << endl;
		}
		if ( npart ) npercent = nsel*100.0/npart;
		cout << "Totals:\t" << npart << " particles\t" << nsel << " selected\t(" 
			<< npercent << "%)\t" << nfil << " filaments\t" << nfn << " filament nodes" << endl;
		if ( npart ) {
			cout << "Group\tCount" << endl;
			for ( g=0; g<gmax; g++ ) if ( ng[g] )
				cout << g << tab << ng[g] << endl;
		}
	}
	
	cout << endl;

	delete[] ng;
	
	return nsel;
}

/**
@brief 	Shows the class averages in a project.
@param 	*project	project parameter structure.
@return long 		number of class averages.
**/
long		project_show_class_averages(Bproject* project)
{
	if ( !verbose ) return 0;
	
	if ( !project->class_avg ) return 0;
	
	long				n(0);
	long				npart = project_count_mg_part_selected(project);
	Bparticle*			part;
	
	cout << "Class averages:" << endl;
	cout << "Class\tNumber\t   %    \tFOM" << endl;
	for ( part = project->class_avg; part; part = part->next ) {
		cout << part->id << tab << part->sel << tab << 
			setw(8) << part->sel*100.0/npart << tab << 
					part->fom[0] << endl;
		n++;
	}
	cout << endl;
	
	return n;
}

/**
@brief 	Dumps particle info in the project hierarchy.
@param 	*project	project parameter structure.
@param 	&filename	output file name.
@return long 				number of particles dumped.
**/
long		project_dump(Bproject* project, Bstring& filename)
{
	if ( !project ) return 0;
	
	int					f;
	long				i, j, npart(0), nsel(0);
	double				npercent(0);
	Euler				euler;
	Bfield*				field;
	Bmicrograph*		mg;
	Bparticle*			part;
	Breconstruction*	rec;
	
	ofstream			fd(filename.c_str());
	if ( fd.fail() ) return -1;

	fd << "Project hierarchy:" << endl;
	for ( field = project->field; field; field = field->next ) {
		fd << "Field: " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next ) {
			i = particle_count(mg->part);
			j = particle_count_selected(mg->part);
			npart += i;
			nsel += j;
			fd << "Micrograph: " << mg->id << tab << i << " particles\t" << j << " selected" << endl;
			if ( !project->euler_flag ) fd << "ID\tox\toy\toz\tvx\tvy\tvz\tva";
			else fd << "ID\tox\toy\toz\tphi\ttheta\tpsi";
			for ( f=0; f<NFOM; f++ ) fd << "\tFOM" << f;
			fd << endl;
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel ) {
					if ( !project->euler_flag ) {
						fd << part->id << tab << part->ori[0] << tab 
							<< part->ori[1] << tab << part->ori[2] << tab
							<< part->view[0] << tab << part->view[1] << tab 
							<< part->view[2] << tab << part->view.angle()*180/M_PI;
					} else {
						euler = Euler(part->view);
						fd << part->id << tab << part->ori[0] << tab 
							<< part->ori[1] << tab << part->ori[2] << tab 
							<< euler.phi()*180/M_PI << tab << euler.theta()*180/M_PI 
							<< tab << euler.psi()*180/M_PI;
					}
					for ( f=0; f<NFOM; f++ ) fd << tab << part->fom[f];
					fd << endl;
				}
			}
			fd << endl;
		}
	}
	if ( npart ) npercent = nsel*100.0/npart;
	fd << "Totals:\t" << npart << " particles\t" << nsel << " selected\t(" << npercent << "%)" << endl << endl;
	
	if ( project->rec && project->rec->part ) {
		for ( npart=nsel=0, rec = project->rec; rec; rec = rec->next ) {
			i = particle_count(rec->part);
			j = particle_count_selected(rec->part);
			npart += i;
			nsel += j;
			fd << "Reconstruction: " << rec->id << tab << i << " particles\t" << j << " selected" << endl;
			if ( !project->euler_flag ) fd << "ID\tox\toy\toz\tvx\tvy\tvz\tva";
			else fd << "ID\tox\toy\toz\tphi\ttheta\tpsi";
			for ( f=0; f<NFOM; f++ ) fd << "\tFOM" << f;
			fd << endl;
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel ) {
					if ( !project->euler_flag ) {
						fd << part->id << tab << part->ori[0] << tab 
							<< part->ori[1] << tab << part->ori[2] << tab
							<< part->view[0] << tab << part->view[1] << tab 
							<< part->view[2] << tab << part->view.angle()*180/M_PI;
					} else {
						euler = Euler(part->view);
						fd << part->id << tab << part->ori[0] << tab 
							<< part->ori[1] << tab << part->ori[2] << tab 
							<< euler.phi()*180/M_PI << tab << euler.theta()*180/M_PI 
							<< tab << euler.psi()*180/M_PI;
					}
					for ( f=0; f<NFOM; f++ ) fd << tab << part->fom[f];
					fd << endl;
				}
			}
			fd << endl;
		}
		if ( npart ) npercent = nsel*100.0/npart;
		if ( npart ) npercent = nsel*100.0/npart;
		fd << "Totals:\t" << npart << " particles\t" << nsel << " selected\t(" << npercent << "%)" << endl << endl;
	}

	fd.close();
	
	return nsel;
}

Bstring		path_update(Bstring& filename, Bstring& path)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG path_update: Old name = " << filename << endl;

	if ( path[-1] != '/' ) path += "/";
//	if ( path == "./" ) path = 0;
	
	Bstring			newname = filename;
	if ( newname.contains("/") ) newname = newname.post_rev('/');
	if ( path.length() && path != "./" ) newname = path + newname;
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG path_update: New name = " << newname << endl;
		
	return newname;
}

/**
@brief 	Sets the path to micrograph files for all the micrographs.
@param 	*project 	project parameter structure.
@param 	&path		path to micrograph files.
@return int			0.

	If the requested path is "." or "./", the path is completely removed.

**/
int			project_set_micrograph_path(Bproject* project, Bstring& path)
{
	if ( !project ) return 0;
	if ( path.length() < 1 ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	
	if ( verbose )
		cout << "Setting micrograph path to " << path << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->fmg.length() ) {
			mg->fmg = path_update(mg->fmg, path);
		}
	}
	
	return 0;
}

/**
@brief 	Sets the path to micrograph frames files for all the micrographs.
@param 	*project 	project parameter structure.
@param 	&path		path to micrograph frames files.
@return int			0.

	If the requested path is "." or "./", the path is completely removed.

**/
int			project_set_frame_path(Bproject* project, Bstring& path)
{
	if ( !project ) return 0;
	if ( path.length() < 1 ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	
	if ( verbose )
		cout << "Setting micrograph frames path to " << path << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->fframe.length() ) {
			mg->fframe = path_update(mg->fframe, path);
		}
	}
	
	return 0;
}

/**
@brief 	Sets the path to powerspectrum files for all the micrographs.
@param 	*project 	project parameter structure.
@param 	&path		path to powerspectrum files.
@return int			0.

	If the requested path is "." or "./", the path is completely removed.

**/
int			project_set_powerspectrum_path(Bproject* project, Bstring& path)
{
	if ( !project ) return 0;
	if ( path.length() < 1 ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	if ( verbose )
		cout << "Setting power spectrum path to " << path << endl;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->fps.length() ) {
			mg->fps = path_update(mg->fps, path);
		}
	}
	
	return 0;
}

/**
@brief 	Sets the path to particle files for all the micrographs.
@param 	*project 	project parameter structure.
@param 	&path		path to particle files.
@return int			0.

	If the requested path is "." or "./", the path is completely removed.

**/
int			project_set_particle_path(Bproject* project, Bstring& path)
{
	if ( !project ) return 0;
	if ( path.length() < 1 ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part = NULL;

	if ( verbose )
		cout << "Setting particle path to " << path << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( rec->part && rec->part->fpart.length() ) {
				for ( part = rec->part; part; part = part->next )
					part->fpart = path_update(part->fpart, path);
			} else if ( rec->fpart.length() ) {
				rec->fpart = path_update(rec->fpart, path);
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( mg->part && mg->part->fpart.length() ) {
					for ( part = mg->part; part; part = part->next )
						part->fpart = path_update(part->fpart, path);
				} else if ( mg->fpart.length() ) {
					mg->fpart = path_update(mg->fpart, path);
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_particle_path: done" << endl;

	return 0;
}

/**
@brief 	Sets the path to filament files for all the micrographs.
@param 	*project 	project parameter structure.
@param 	&path		path to filament files.
@return int			0.

	If the requested path is "." or "./", the path is completely removed.

**/
int			project_set_filament_path(Bproject* project, Bstring& path)
{
	if ( !project ) return 0;
	if ( path.length() < 1 ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bstring				newname, insert = "_fil.";

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_set_particle_path: Old name = " << rec->ffil << endl;
			newname = rec->ffil;
			if ( rec->ffil.length() < 1 )
				newname = rec->frec.pre_rev('.') + insert + rec->frec.post_rev('.');
//			if ( newname.contains("/") ) newname = newname.post_rev('/');
//			rec->ffil = newname;
//			if ( path.length() ) rec->ffil = path + rec->ffil;
			rec->ffil = path_update(rec->ffil, path);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG project_set_filament_path: New name = " << rec->ffil << endl;
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG project_set_particle_path: Old name = " << mg->ffil << endl;
				newname = mg->ffil;
				if ( mg->ffil.length() < 1 )
					newname = mg->fmg.pre_rev('.') + insert + mg->fmg.post_rev('.');
//				if ( newname.contains("/") ) newname = newname.post_rev('/');
//				mg->ffil = newname;
//				if ( path.length() ) mg->ffil = path + mg->ffil;
				mg->ffil = path_update(mg->ffil, path);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG project_set_filament_path: New name = " << mg->ffil << endl;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Sets the magnification for all the micrographs.
@param 	*project 	project parameter structure.
@param 	mag			magnification.
@return int					0.
**/
int			project_set_magnification(Bproject* project, double mag)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->magnification = mag;
	
	return 0;
}

/**
@brief 	Sets the scan sampling for all the micrographs.
@param 	*project 	project parameter structure.
@param 	sampling	sampling.
@return int			0.

	Sampling is isotropic.

**/
int			project_set_scan_sampling(Bproject* project, double sampling)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->sampling = sampling;
	
	return 0;
}

/**
@brief 	Sets the pixel size for all the micrographs.
@param 	*project 	project parameter structure.
@param 	pixel_size	pixel size.
@return int			0.
**/
int			project_set_mg_pixel_size(Bproject* project, Vector3<double> pixel_size)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_mg_pixel_size: pixel_size = " << pixel_size << endl;

	pixel_size[2] = 1;
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
				mg->pixel_size = pixel_size;
	
	return 0;
}

/**
@brief 	Sets the pixel size for all the micrograph frames.
@param 	*project 	project parameter structure.
@param 	pixel_size	pixel size.
@return int			0.
**/
int			project_set_frame_pixel_size(Bproject* project, Vector3<double> pixel_size)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_mg_pixel_size: pixel_size = " << pixel_size << endl;

	pixel_size[2] = 1;
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			mg->frame_pixel_size = pixel_size;
			if ( mg->pixel_size.volume() < 1e-3 || mg->pixel_size.volume() > 10*mg->frame_pixel_size.volume() )
				mg->pixel_size = pixel_size;
		}
	}
	
	return 0;
}

/**
@brief 	Sets the pixel size for all the reconstructions.
@param 	*project 	project parameter structure.
@param 	pixel_size	pixel size.
@return int			0.
**/
int			project_set_rec_voxel_size(Bproject* project, Vector3<double> pixel_size)
{
	if ( !project ) return 0;
	
	Breconstruction*	rec;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_rec_voxel_size: voxel_size = " << pixel_size << endl;

	for ( rec = project->rec; rec; rec = rec->next )
			rec->voxel_size = pixel_size;
	
	return 0;
}

/**
@brief 	Sets the pixel size for all the particles.
@param 	*project 	project parameter structure.
@param 	pixel_size	pixel size.
@return int			0.
**/
int			project_set_part_pixel_size(Bproject* project, Vector3<double> pixel_size)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_part_pixel_size: pixel_size = " << pixel_size << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				part->pixel_size = pixel_size;
	} else {
		pixel_size[2] = 1;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					part->pixel_size = pixel_size;
	}
	
	return 0;
}

/**
@brief 	Sets the tilt parameters for all the micrographs.
@param 	*project 			project parameter structure.
@param 	tilt_axis			micrograph tilt axis angle.
@param 	tilt_angle			micrograph tilt angle.
@return int					0.
**/
int			project_set_tilt(Bproject* project, double tilt_axis, double tilt_angle)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			mg->tilt_axis = tilt_axis;
			mg->tilt_angle = tilt_angle;
		}
	}
	
	return 0;
}

/**
@brief 	Sets the aquisition time of all the micrographs.
@param 	*project	project parameter structure.
@param 	exposure	aquisition time  (seconds).
@return int			0.
**/
int			project_set_exposure(Bproject* project, double exposure)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->exposure = exposure;
	
	return 0;
}

/**
@brief 	Sets the electron dose of all the micrographs.
@param 	*project	project parameter structure.
@param 	dose		electron dose (e/A2).
@return int			0.
**/
int			project_set_dose(Bproject* project, double dose)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg->dose = dose;
	
	return 0;
}

/**
@brief 	Sets the dose fractionation scheme of all the micrographs.
@param 	*project	project parameter structure.
@param 	dose_frac	electron dose fractionation scheme.
@return int			0.
**/
int			project_set_dose(Bproject* project, JSvalue& dose_frac)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	Bframe*			frame;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_set_dose: " << dose_frac << endl;

	double			dose(0), exposure(0), dose_per_frame(0), frame_rate(0);
	if ( dose_frac.exists("dose_per_frame") )
		dose_per_frame = dose_frac["dose_per_frame"].real();
	if ( dose_frac.exists("frame_rate") )
		frame_rate = dose_frac["frame_rate"].real();
	else if ( dose_frac.exists("exposure") )
		exposure = dose_frac["exposure"].real();

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->frame ) {
				if ( dose_per_frame )
					for ( dose = 0, frame = mg->frame; frame; frame = frame->next )
						dose += dose_per_frame;
				if ( frame_rate )
					for ( exposure = 0, frame = mg->frame; frame; frame = frame->next )
						exposure += 1/frame_rate;
			} else {
				dose = dose_per_frame;	// Assume a single image
				if ( frame_rate )
					exposure = 1/frame_rate;
			}
			if ( dose ) mg->dose = dose;
			if ( exposure ) mg->exposure = exposure;
		}
	}
	
	return 0;
}

/**
@brief 	Sets micrograph origins to the given origin.
@param 	*project	project parameter structure.
@param 	origin		particle origin.
@return int 			0.

	For each micrograph the micrograph origin is set to the given origin.

**/
int			project_set_micrograph_origins(Bproject* project, Vector3<double> origin)
{
	if ( !project ) return 0;
	if ( !project->field ) return 0;
	if ( !project->field->mg ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	
	if ( verbose & VERB_FULL )
		cout << "Setting micrograph origins to " << origin << endl;

	for ( field = project->field; field; field = field->next ) {
		field->origin = origin;
		for ( mg = field->mg; mg; mg = mg->next )
			mg->origin = origin;
	}

	return 0;
}

/**
@brief 	Adds the particle origins to the micrograph coordinates.
@param 	*project	project parameter structure.
@return int 		0.

	The coordinates of particles from a micrograph and particle origins
	in an image processing parameter structure are added. The old
	micrograph coordinates are overwritten with the results.
	Requirement: The box radius for picked particles must be specified.

**/
int			project_add_origins_to_coords(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				part->loc += part->ori - rec->box_size/2;
				part->ori = rec->box_size/2;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					part->loc += part->ori - mg->box_size/2;
					part->ori = mg->box_size/2;
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Flip origin coordinates.
@param 	*project	project parameter structure.
@param 	flip		axes to flip.
@return int 		0.

	The particle origins are reversed.
	The specification of axes to flip is embedded in the flip number:
		first bit  - x
		second bit - y
		third bit  - z
**/
int			project_flip_origins(Bproject* project, int flip)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	int					fx(flip & 1);
	int					fy(flip & 2);
	int					fz(flip & 4);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Flipping origins along ";
		if ( fx ) cout << "x";
		if ( fy ) cout << "y";
		if ( fz ) cout << "z";
		cout << endl;
	}

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( fx ) part->ori[0] = rec->box_size[0] - part->ori[0];
				if ( fy ) part->ori[1] = rec->box_size[1] - part->ori[1];
				if ( fz ) part->ori[2] = rec->box_size[2] - part->ori[2];
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( fx ) part->ori[0] = mg->box_size[0] - part->ori[0];
					if ( fy ) part->ori[1] = mg->box_size[1] - part->ori[1];
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Renumbers particles.
@param 	*project	project parameter structure.
@return long 		total number of particles.

	For each micrograph the particles are renumbered starting from 1.

**/
long		project_renumber_particles(Bproject* project)
{
	if ( !project ) return 0;
	
	long				n, np(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_FULL )
		cout << "Renumbering particles" << endl << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( n=1, part = rec->part; part; part = part->next, n++, np++ )
				part->id = n;
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( n=1, part = mg->part; part; part = part->next, n++, np++ )
					part->id = n;
	}

	return np;
}

/**
@brief 	Sets particle box size to the given value.
@param 	*project	project parameter structure.
@param 	box_size	particle box size.
@return int 		0.

	For each micrograph or reconstruction the particle box size is set 
	to the given radius and the particle origins are adjusted as well.

**/
int			project_set_particle_box_size(Bproject* project, Vector3<long> box_size)
{
	if ( !project ) return 0;
	
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_FULL )
		cout << "Setting the box size to " << box_size << endl << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next )
				part->ori -= rec->box_size/2;
			rec->box_size = box_size;
			for ( part = rec->part; part; part = part->next )
				part->ori += rec->box_size/2;
		}
	} else {
		box_size[2] = 1;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next )
					part->ori -= mg->box_size/2;
				mg->box_size = box_size;
				for ( part = mg->part; part; part = part->next )
					part->ori += mg->box_size/2;
			}
	}

	return 0;
}

/**
@brief 	Sets particle box size to the given value.
@param 	*project	project parameter structure.
@param 	box_size	particle box size.
@return int 		0.

	For each micrograph or reconstruction the particle box size is set 
	to the given radius and the particle origins are adjusted as well.

**/
int			project_set_particle_box_size(Bproject* project, long box_size)
{
	if ( !project ) return 0;
	
	Vector3<long>		bsize(box_size, box_size, 1);
	if ( project->select ) bsize[2] = box_size;

	return project_set_particle_box_size(project, bsize);
}	

/**
@brief 	Sets particle origins to given values.
@param 	*project		project parameter structure.
@param 	origin	particle origin.
@return int 					0.

	For each micrograph or reconstruction the particle origin is set to the given origin.

**/
int			project_set_particle_origins(Bproject* project, Vector3<double> origin)
{
	if ( !project ) return 0;
	
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_FULL )
		cout << "Setting origins to " << origin << endl << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				part->ori = origin;
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					part->ori = origin;
	}

	return 0;
}

/**
@brief 	Sets particle views within the asymmetric unit.
@param 	*project		project parameter structure.
@param 	&symmetry_asu	point group string.
@return int						0.

	For each micrograph or reconstruction the particle view is set to within the asymmetric unit.

**/
int			project_set_particle_asu_views(Bproject* project, Bstring& symmetry_asu)
{
	if ( !project ) return 0;
	
	Bsymmetry			sym(symmetry_asu);

	project_set_particle_asu_views(project, sym);

	return 0;
}

int			project_set_particle_asu_views(Bproject* project, Bsymmetry& sym)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose )
		cout << "Setting particle views to the asymmetric unit for point group " << sym.label() << endl << endl;

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				part->view = find_asymmetric_unit_view(sym, part->view);
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					part->view = find_asymmetric_unit_view(sym, part->view);
	}

	return 0;
}

/**
@author Eduardo Sanz-Garcia
@brief 	Rotates particle views with respect to a reference view.
@param 	*project		project parameter structure.
@param 	view				reference view.
@return int						0.
**/
int			project_rotate_particle_views(Bproject* project, View view)
{
	if ( !project ) return 0;
	
	if ( view[2] >= 0.999999 ) return 0;	// Standard view
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Quaternion			q;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Rotating particle views relative to the reference view:" << endl;
		cout << view << endl << endl;
	}

	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				q = view.quaternion() * part->view.quaternion();
				part->view = View(q);
			}
	} else {
		for ( field=project->field; field; field=field->next )
			for ( mg=field->mg; mg ; mg=mg->next )
				for ( part=mg->part; part; part=part->next) {
					q = view.quaternion() * part->view.quaternion();
					part->view = View(q);
				}
	}

	return 0;
}

/**
@brief 	Change particle pixel sizes for particles based on map magnifications.
@param 	*project		project parameter structure.
@param 	mag_num				number of maps;
@param 	*mag				array of map magnifications.
@return int						0.
**/
int			project_apply_map_magnifications(Bproject* project, int mag_num, float* mag)
{
	if ( !project ) return 0;
	
	int					i;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( verbose & VERB_PROCESS ) {
		cout << "Applying map magnifications to particles:" << endl;
		cout << "Number of reference maps:       " << mag_num << endl;
		for ( i=0; i<mag_num; i++ ) cout << i+1 << tab << mag[i] << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel > 0 && part->sel <= mag_num ) {
						part->mag = mag[part->sel-1];
					}
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel > 0 && part->sel <= mag_num ) {
					part->mag = mag[part->sel-1];
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Reset a particle parameter from its micrograph.
@param 	*project		project parameter structure.
@param 	&reset			a string specifying the parameter.
@return long			number of particles changed.
**/
long			project_reset(Bproject* project, Bstring& reset)
{
	if ( !project ) return 0;
	
	if ( reset[0] == 'd' ) reset = "defocus";
	else return 0;
	
	long				np(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	if ( verbose )
		cout << "Resetting " << reset << " from the micrographs" << endl;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel && mg->ctf ) {
						part->def = mg->ctf->defocus_average();
						np++;
					}
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel && rec->ctf ) {
					part->def = rec->ctf->defocus_average();
					np++;
				}
			}
		}
	}
	
	if ( verbose )
		cout << "Number of particles changed: " << np << endl << endl;

	return np;
}

/**
@brief 	Retrieves the particle views from a project.
@param 	*project	project parameter structure.
@param 	selection		selection number (-1 selects positives, 0 selects all).
@return View* 				linked list of views.
**/
View*		views_from_project(Bproject* project, int selection)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	View*				view = NULL;
	View*				v = NULL;
	
	if ( project->select < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Retrieving the particle views from the micrographs" << endl << endl;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					if ( selection==0 || ( selection==part->sel ) || ( selection<0 && part->sel>0 ) ) {
						v = (View *) add_item((char **) &v, sizeof(View));
						if ( !view ) view = v;
						*v = part->view;
						v->next = NULL;
					}
				}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "Retrieving the particle views from the reconstructions" << endl << endl;
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				if ( selection==0 || ( selection==part->sel ) || ( selection<0 && part->sel>0 ) ) {
					v = (View *) add_item((char **) &v, sizeof(View));
					if ( !view ) view = v;
					*v = part->view;
					v->next = NULL;
				}
			}
	}
	
	return view;
}

/**
@brief 	Returns the tag associated with a particular FOM.
@param 	fom_type	type of FOM.
@return Bstring 		FOM tag.
**/
Bstring		get_fom_tag(FOMType fom_type)
{
	Bstring		fom_tag;
	
	switch ( fom_type ) {
		case FOM: fom_tag = PARTICLE_FOM; break;
		case FOM_CC: fom_tag = PARTICLE_CC; break;
		case FOM_CV: fom_tag = PARTICLE_FOM_CV; break;
		case FOM_SNR: fom_tag = PARTICLE_FOM_SNR; break;
		case FOM_CC_AVG: fom_tag = PARTICLE_FOM_AVG; break;
		case FOM_CC_STD: fom_tag = PARTICLE_FOM_STD; break;
		case FOM_HAND_A: fom_tag = PARTICLE_HANDA_FOM; break;
		case FOM_HAND_B: fom_tag = PARTICLE_HANDB_FOM; break;
		case FOM_PFT_CC: fom_tag = PARTICLE_PFT_CC; break;
		case FOM_PFT_PRJ: fom_tag = PARTICLE_PRJ_CC; break;
		case FOM_PFT_CMP: fom_tag = PARTICLE_CMP_CC; break;
		case FOM_RFACTORAB: fom_tag = PARTICLE_RFACTORAB; break;
		case COVERAGE: fom_tag = PARTICLE_COVERAGE; break;
		case DENSITY: fom_tag = PARTICLE_DENSITY; break;
		case FOM1: fom_tag = PARTICLE_FOM1; break;
		case FOM2: fom_tag = PARTICLE_FOM2; break;
		case FOM3: fom_tag = PARTICLE_FOM3; break;
		case FOM4: fom_tag = PARTICLE_FOM4; break;
		case FOM5: fom_tag = PARTICLE_FOM5; break;
		case FOM6: fom_tag = PARTICLE_FOM6; break;
		case FOM7: fom_tag = PARTICLE_FOM7; break;
		case FOM8: fom_tag = PARTICLE_FOM8; break;
		case FOM9: fom_tag = PARTICLE_FOM9; break;
		default: break;
	}
	
	return fom_tag;
}
