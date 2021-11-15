/**
@file	mg_select.cpp
@brief	Functions for micrograph processing
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20210515
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
@brief 	Finds a field-of-view based on its identifier.
@param 	*field		pointer to first field-of-view in the list.
@param 	&field_id	field-of-view identifier.
@return Bfield* 			field-of-view or NULL if not found.

	The function searches a linked list for the field-of-view
	identifier and returns a pointer to that structure or NULL if
	it cannot find it.

**/
Bfield* 	field_find_id(Bfield* field, Bstring& field_id)
{
	while ( field && field->id != field_id )
		field = field->next;
	
	return field;
}

/**
@author  David Belnap
@brief 	Returns a pointer to a specific micrograph in a field-of-view 
	structure based on a unique characteristic of that micrograph.
@param 	*field			field-of-view.
@param 	mg_select			selection criterion.
@param 	mg_index         	Reference by its index in field (for mg_ref_select=0|1|2).
@param 	mg_ang			Reference by micrograph rotation or tilt angle (for mg_ref_select=3|4).
@return Bmicrograph*			Pointer to the selected micrograph

	Calls functions that find a micrograph by index number within the
	field, by focus level, in-plane rotational angle, tilt angle, or
	tilt angle plus focus level or index number:

	criterion                        mg_select   mg_index  mg_ang
	--------------------------------------------------------------
	nth closest-to-focus                 0           n        -
	nth farthest-from-focus              1           n        -
	index number in series (1...n)       2         index      -
	in-plane rotation angle              3           -      angle
	tilt angle                           4           -      angle

**/
Bmicrograph*  field_find_micrograph(Bfield* field, int mg_select, int mg_index, double mg_ang)
{
	Bmicrograph*  mg = NULL;    // Selected micrograph

	switch ( mg_select )  {
		case 0:
		case 1:
			mg = field_find_micrograph_by_focus(field, mg_select, mg_index);
			break;
		case 2:
			mg = field_find_micrograph_n(field, mg_index);
			break;
		case 3:
			mg = field_find_micrograph_by_rotang(field, mg_ang);
			break;
		case 4:
			mg = field_find_micrograph_by_tiltang(field, mg_ang);
			break;
	}
  
	// Error checking (primarily to catch programming mistakes)
	if ( !mg )
		error_show("Error: No appropriate micrograph has been selected!", __FILE__, __LINE__);

	return mg;
}

/**
@author  David Belnap and Bernard Heymann
@brief 	Returns a pointer to the nth micrograph in a field-of-view 
	structure
@param 	*field	a field-of-view
@param 	n        	nth micrograph index (first index number = 1)
@return Bmicrograph*	Pointer to the nth micrograph

	Loop through micrographs in the field-of-view until the nth 
	micrograph is reached.  Return pointer to that micrograph.
	Tests whether n is within the appropriate range.

**/
Bmicrograph*  field_find_micrograph_n(Bfield* field, int n)
{
	int           i;				// index
	Bmicrograph*  mg = NULL;		// micrograph
  
	for ( i=1, mg = field->mg; mg && i < n; mg = mg->next, i++ ) ;
  
	if ( !mg )  {  // Make sure n is in appropriate range
		cerr << "Error:  The micrograph index you wish to find (" << n << ") is greater" << endl;
		cerr << "than the number of micrographs (" << i << ") in the field (" << field->id << ")" << endl;
	}

	return mg;
}

/**
@author  David Belnap and Bernard Heymann
@brief 	Returns a pointer to the micrograph in a field-of-view structure
	based on the focus level and iteration index
@param 	*field       a field-of-view structure
@param	focus_opt   0=closest-to-focus, 1=farthest-from-focus
@param	iselect       0,1,...,n; select (index+1)th-closest or (index+1)th-farthest focus
@return Bmicrograph*   mg          Micrograph at specified focus level
	Loops through micrographs in the field-of-view.  Orders focus 
	values in an array from smallest to largest.  Selects the 
	(index+1)th closest-to-focus or farthest-from-focus value.
	(focus_opt determines whether the closest or farthest value is
	selected.)  Loops through micrographs again to find the 
	micrograph with the selected value.  Tests if there are zero or 
	multiple micrographs with the same focus value and exits program 
	if so.  Otherwise, the pointer to the micrograph with the 
	specified focus level is returned.
	  An example, if index=1 and focus_opt=0, a pointer to the 
	2nd closest-to-focus micrograph will be returned.  To get the
	closest-to-focus or farthest-from-focus value, the index must 
	be zero.

**/
Bmicrograph*  field_find_micrograph_by_focus(Bfield* field, int focus_opt, int iselect)
{
	int				i;				// index
	int				nmgs;			// number of micrographs in field
	int				count;			// used to test if multiple micrographs have same focus level
	double			focus(0);		// designated focus level, depends on focus_opt and index
	Bmicrograph*  	mg = NULL;		// pointer to micrograph at desired focal level

	nmgs = field_count_micrographs(field);
  
	float*  focal_series = new float[nmgs];

	// Put defocal values into an array, then put in order min to max
	for ( i = 0, mg = field->mg; mg; mg = mg->next, i++ ) {
		focal_series[i] = 0;
		if ( mg->ctf )	focal_series[i] = mg->ctf->defocus_average();
		else cerr << "Warning in field_find_micrograph_by_focus: CTF structure not allocated!" << endl;
	}

	qsort((void *) focal_series, nmgs, sizeof(float), (int (*)(const void *, const void *)) QsortSmallToLargeFloat);

	// Select (index+1)th closest-to-focus or farthest-from-focus value
	if ( focus_opt == 0 )  focus = focal_series[0+iselect];
	if ( focus_opt == 1 )  focus = focal_series[nmgs-1-iselect];
	for ( count = 0, i=0; i<nmgs; i++ ) if ( focus == focal_series[i] ) count++;

	delete[] focal_series;

	// Find micrograph with the selected focus level
	for ( mg = field->mg; mg && mg->ctf && mg->ctf->defocus_average() != focus; mg = mg->next ) ;

	if (verbose & VERB_PROCESS)
		cout << "Focus of selected micrograph: " << focus << endl;

	if ( (count > 1) || (count == 0) )  {
		cerr << "Error:  There are " << count << " micrographs with focus levels of " << focus << endl;
		cerr << "Please check data or select with a different criterion" << endl;
	}

	return mg;
}

/**
@author  David Belnap and Bernard Heymann
@brief 	Returns a pointer to the micrograph in a field-of-view structure
	with the specified in-plane rotational angle
@param 	*field    a field-of-view
@param	rotang   in-plane rotational angle of micrograph (in radians)
@return Bmicrograph*   mg       Pointer to micrograph with rotang

	Loop through all micrographs in the field-of-view to find
	micrograph with the specified angle.  Test if there are zero or
	multiple micrographs with the same rotational angle, halt program 
	if true.  If not true, then return pointer to the micrograph with 
	the specified angle.

**/
Bmicrograph*  field_find_micrograph_by_rotang(Bfield* field, double rotang)
{
	int             count;              // used to test if multiple micrographs have same rotang
	double           closest_angle(0);	// Closest angle to requested angle
	double           mindiff = 1e37;     // Minimum difference between requested and available angle
	Bmicrograph*    mg = NULL;          // pointer to micrograph with specified rotang

	// Find micrograph with the input rot. angle
	for ( mindiff = 1e37, mg = field->mg; mg; mg = mg->next ) {
		if ( fabs(rotang - mg->rot_angle) < mindiff ) {
			mindiff = fabs(rotang - mg->rot_angle);
			closest_angle = mg->rot_angle;
		}
	}
	
	for ( count = 0, mg = field->mg; mg; mg = mg->next)
		if ( fabs(closest_angle - mg->rot_angle) < 1e-10 ) count++;

	for ( mg = field->mg; mg && closest_angle != mg->rot_angle; mg = mg->next ) ;

	if ( (count > 1) || (count == 0) )  {
		cerr << "Error:  There are " << count << " micrographs with the rotational angle of " 
			<< closest_angle*180/M_PI << endl;
		cerr << "Please check data or select with a different criterion" << endl;
	}

	return mg;
}

/**
@author  David Belnap
@brief 	Returns a pointer to the micrograph in a field-of-view structure
	with the specified tilt angle
@param 	*field    	a field-of-view
@param 	tiltang  	tilt angle of micrograph (in radians)
@return Bmicrograph*		micrograph with rotang

	Loop through all micrographs in the field-of-view to find
	micrograph with the specified tilt angle.  Test if there are zero
	or multiple micrographs with the same angle, halt program if true.
	If not true, then return pointer to the micrograph with the 
	specified angle.

**/
Bmicrograph*  field_find_micrograph_by_tiltang(Bfield* field, double tiltang)
{
	int             count;              // used to test if multiple micrographs have same rotang
	double    		closest_angle(0);	// Closest angle to requested angle
	double          mindiff = 1e37;     // Minimum difference between requested and available angle
	Bmicrograph*    mg = NULL;          // pointer to micrograph with specified rotang

	// Find micrograph with the input tilt angle
	for ( mindiff = 1e37, mg = field->mg; mg; mg = mg->next ) {
		if ( fabs(tiltang - mg->tilt_angle) < mindiff ) {
			mindiff = fabs(tiltang - mg->tilt_angle);
			closest_angle = mg->tilt_angle;
		}
	}

	for ( count = 0, mg = field->mg; mg; mg = mg->next)
		if ( fabs(closest_angle - mg->tilt_angle) < 1e-10 ) count++;

	for ( mg = field->mg; mg && (closest_angle != mg->tilt_angle); mg = mg->next ) ;

	if ( (count > 1) || (count == 0) )  {
		cerr << "Error:  There are " << count << " micrographs with the tilt angle of " 
			<< closest_angle*180/M_PI << endl;
		cerr << "Please check data or select with a different criterion" << endl;
	}

	return mg;
}

/**
@brief 	Finds the micrograph closest to a zero degree tilt in a series.
@param 	*field		field-of-view.
@return Bmicrograph* 		micrograph closest to zero degrees tilt.

	The first micrograph with the smallest deviation from zero degree
	tilt is returned.

**/
Bmicrograph*	field_find_zero_tilt_mg(Bfield* field)
{
	if ( !field ) return NULL;
	
	double			min_angle(M_PI);
	Bmicrograph*	mg;
	Bmicrograph*	mg_ref = field->mg;
	
	for ( mg = field->mg; mg; mg = mg->next ) {
		if ( min_angle > fabs(mg->tilt_angle) ) {
			min_angle = fabs(mg->tilt_angle);
			mg_ref = mg;
		}
//		cout << mg->tilt_angle << tab << min_angle << endl;
	}
	
	if ( !mg_ref )
		cerr << "Warning: No zero tilt micrograph found!" << endl;
	
	return mg_ref;
}

/**
@brief 	Finds the micrograph closest to a zero degree tilt in a series with markers.
@param 	*field		field-of-view.
@return Bmicrograph* 		micrograph closest to zero degrees tilt with markers.

	The first micrograph with the smallest deviation from zero degree
	tilt with defined markers is returned.

**/
Bmicrograph*	field_find_low_tilt_mg_with_markers(Bfield* field)
{
	if ( !field ) return NULL;
	
	double			min_angle = M_PI;
	Bmicrograph*	mg;
	Bmicrograph*	mg_ref = NULL;
	
	for ( mg = field->mg; mg; mg = mg->next ) {
		if ( mg->mark && min_angle > fabs(mg->tilt_angle) ) {
			min_angle = fabs(mg->tilt_angle);
			mg_ref = mg;
		}
	}
	
	return mg_ref;
}

/**
@brief 	Counts all the fields in the project structure.
@param 	*project	project parameter structure.
@return long				number of fields in project.
**/
long		project_count_fields(Bproject* project)
{
	Bfield*			field;
	long			nfield(0);

	for ( field = project->field; field; field = field->next ) nfield++;
	
	return nfield;
}

/**
@brief 	Counts all the micrographs in the project structure.
@param 	*project	project parameter structure.
@return long				number of micrographs in project.
**/
long		project_count_micrographs(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	long			nmg(0);

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) nmg++;
	
	return nmg;
}

/**
@brief 	Counts all the selected micrographs in the project structure.
@param 	*project	project parameter structure.
@return long				number of micrographs selected in project.
**/
long		project_count_mg_selected(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	long			nmg(0);

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) nmg++;
	
	return nmg;
}

/**
@brief 	Counts all the reconstructions in the project structure.
@param 	*project	project parameter structure.
@return long				number of reconstructions in project.
**/
long		project_count_reconstructions(Bproject* project)
{
	Breconstruction*	rec;
	long				nrec(0);

	for ( rec = project->rec; rec; rec = rec->next ) nrec++;
	
	return nrec;
}

/**
@brief 	Counts all the selected reconstructions in the project structure.
@param 	*project	project parameter structure.
@return long				number of reconstructions selected in project.
**/
long		project_count_rec_selected(Bproject* project)
{
	Breconstruction*	rec;
	long				nrec(0);

	for ( rec = project->rec; rec; rec = rec->next )
		if ( rec->select ) nrec++;
	
	return nrec;
}

/**
@brief 	Counts all the particles in the project structure.
@param 	*project	project parameter structure.
@return long				number of particles in project.
**/
long		project_count_mg_particles(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	long			npart(0);

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next ) npart++;
	
	return npart;
}

/**
@brief 	Counts the number of particles in a project.
@param 	*project	project.
@return long 				number of particles selected.

	The function counts all the selected particles in a project.

**/
long 		project_count_mg_part_selected(Bproject* project)
{
	long 			nsel(0);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next )
				if ( part->sel > 0 ) nsel++;
	
	return nsel;
}

/**
@brief 	Counts the number of particles in a project.
@param 	*project	project.
@param 	num_select	selection number.
@return long 			number of particles selected.

	The function counts all the selected particles in a project.

**/
long 		project_count_mg_part_selected(Bproject* project, int num_select)
{
	if ( num_select < 0 ) return project_count_mg_part_selected(project);
	
	long 			nsel(0);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next )
				if ( part->sel == num_select ) nsel++;
	
	return nsel;
}

/**
@brief 	Counts the number of groups of particles in the project structure.
@param 	*project	project parameter structure.
@return long			number of groups in project.
**/
long		project_count_mg_groups(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	long			ngrp(0);

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part && part->next; part = part->next )
				if ( part->group != part->next->group ) ngrp++;
			if ( part ) ngrp++;
		}
	
	return ngrp;
}

/**
@brief 	Counts the number of groups of particles in the project structure.
@param 	*project	project parameter structure.
@return long			number of groups selected in project.
**/
long		project_count_mg_groups_selected(Bproject* project)
{
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	long			ngrp(0);

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part && part->next; part = part->next ) if ( part->sel )
				if ( part->group != part->next->group ) ngrp++;
			if ( part ) ngrp++;
		}
	
	return ngrp;
}


/**
@brief 	Counts all the particles in the project structure.
@param 	*project	project parameter structure.
@return long			number of particles in project.
**/
long		project_count_rec_particles(Bproject* project)
{
	Breconstruction*	rec;
	Bparticle*			part;
	long				npart(0);

	for ( rec = project->rec; rec; rec = rec->next )
		for ( part = rec->part; part; part = part->next ) npart++;
	
	return npart;
}

/**
@brief 	Counts all the particles in the project structure.
@param 	*project	project parameter structure.
@return long				number of particles in project.
**/
long		project_count_rec_part_selected(Bproject* project)
{
	Breconstruction*	rec;
	Bparticle*			part;
	long				npart(0);

	for ( rec = project->rec; rec; rec = rec->next )
		for ( part = rec->part; part; part = part->next )
			if ( part->sel ) npart++;
	
	return npart;
}

/**
@brief 	Counts all the particles in the project structure.
@param 	*project	project parameter structure.
@return long				number of particles in project.
**/
long		project_count_rec_groups(Bproject* project)
{
	Breconstruction*	rec;
	Bparticle*			part;
	long				ngrp(0);

	for ( rec = project->rec; rec; rec = rec->next ) {
		for ( part = rec->part; part && part->next; part = part->next )
			if ( part->group != part->next->group ) ngrp++;
		if ( part ) ngrp++;
	}
	
	return ngrp;
}

/**
@brief 	Counts the number of filaments in a project.
@param 	*project	project.
@return long 			number of filaments.

	The function counts all the filaments in a project.

**/
long 		project_count_mg_filaments(Bproject* project)
{
	long 			nfil(0);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bfilament*		fil = NULL;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( fil=mg->fil; fil; fil=fil->next )
				nfil++;
	
	return nfil;
}

/**
@brief 	Counts the number of filament nodes in a project.
@param 	*project	project.
@return long 			number of filament nodes.

	The function counts all the filament nodes in a project.

**/
long 		project_count_mg_filament_nodes(Bproject* project)
{
	long 			nfn(0);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bfilament*		fil = NULL;
	Bfilnode*		fnode = NULL;
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( fil=mg->fil; fil; fil=fil->next )
				for ( fnode=fil->node; fnode; fnode=fnode->next )
					nfn++;
	
	return nfn;
}

/**
@brief 	Counts the number of filaments in a project.
@param 	*project	project.
@return long 			number of filaments.

	The function counts all the filaments in a project.

**/
long 		project_count_rec_filaments(Bproject* project)
{
	long				nfil(0);
	Breconstruction*	rec;
	Bfilament*			fil = NULL;
	
	for ( rec = project->rec; rec; rec = rec->next )
		for ( fil = rec->fil; fil; fil = fil->next )
			nfil++;
	
	return nfil;
}

/**
@brief 	Counts the number of filaments in a project.

	The function counts all the filaments in a project.

@param 	*project	project.
@return long 				number of filaments.
**/
long 		project_count_rec_filament_nodes(Bproject* project)
{
	long				nfil(0);
	Breconstruction*	rec;
	Bfilament*			fil = NULL;
	Bfilnode*			fnode = NULL;
	
	for ( rec = project->rec; rec; rec = rec->next )
		for ( fil = rec->fil; fil; fil = fil->next )
			for ( fnode = fil->node; fnode; fnode = fnode->next )
				nfil++;
	
	return nfil;
}

/**
@brief 	Counts the number of fields-of-view in a linked list.

	The function counts fields-of-view in the list from the given pointer,
	and does not count fields-of-view prior to this one.

@param 	*field 		pointer to any field-of-view in the list.
@return long 				number of fields-of-view.
**/
long 		field_count(Bfield* field)
{
	long 		nfield;
	
	for ( nfield=0; field; nfield++ ) field = field->next;
	
	return nfield;
}

/**
@brief 	Counts all the micrographs in the field-of-view structure.
@param 	*field		field-of-view parameter structure.
@return long				number of micrographs of field.
**/
long		field_count_micrographs(Bfield* field)
{
	Bmicrograph*	mg;
	long			nmg(0);

	for ( mg = field->mg; mg; mg = mg->next ) nmg++;
	
	return nmg;
}

/**
@brief 	Counts all the micrographs in the field-of-view structure.
@param 	*field		field-of-view parameter structure.
@return long				number of micrographs of field.
**/
long		field_count_mg_selected(Bfield* field)
{
	Bmicrograph*	mg;
	long			nmg(0);

	for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) nmg++;
	
	return nmg;
}

/**
@brief 	Counts all the particles in the field-of-view structure.
@param 	*field		field-of-view parameter structure.
@return long				number of particles in field.
**/
long		field_count_particles(Bfield* field)
{
	Bmicrograph*	mg;
	Bparticle*		part;
	long			npart(0);

	for ( mg = field->mg; mg; mg = mg->next )
		for ( part = mg->part; part; part = part->next ) npart++;
	
	return npart;
}

/**
@brief 	Counts the number of micrographs in a linked list.

	The function counts micrographs in the list from the given pointer,
	and does not count micrographs prior to this one.

@param 	*mg 	pointer to any micrograph in the list.
@return long 				number of micrographs.
**/
long 		micrograph_count(Bmicrograph* mg)
{
	long 		nmg;
	
	for ( nmg=0; mg; nmg++ ) mg = mg->next;
	
	return nmg;
}

/**
@brief 	Counts all the particles in the micrograph structure.
@param 	*mg		micrograph parameter structure.
@return long				number of particles in micrograph.
**/
long		micrograph_count_particles(Bmicrograph* mg)
{
	Bparticle*		part;
	long			npart(0);

	for ( part = mg->part; part; part = part->next ) npart++;
	
	return npart;
}

/**
@brief 	Counts the number of particles in a linked list.

	The function counts particles in the list from the given pointer,
	and does not count particles prior to this one.

@param 	*part		pointer to any particle in the list.
@return long 				number of particles.
**/
long 		particle_count(Bparticle* part)
{
	long 		npart;
	
	for ( npart=0; part; npart++ ) part = part->next;
	
	return npart;
}

/**
@brief 	Counts the number of particles in a linked list.

	The function counts particles in the list from the given pointer,
	and does not count particles prior to this one.

@param 	*part		pointer to any particle in the list.
@return long 				number of particles selected.
**/
long 		particle_count_selected(Bparticle* part)
{
	long 			nsel;
	
	for ( nsel=0; part; part=part->next ) if ( part->sel > 0 ) nsel++;
	
	return nsel;
}

/**
@brief 	Counts the number of filaments in a linked list of filaments.

	The function counts filaments in the list from the given pointer,
	and does not count filaments prior to this one.

@param 	*fil		pointer to any filament in the list.
@return long 				number of filaments.
**/
long 		filament_count(Bfilament* fil)
{
	long 		nfil;
	
	for ( nfil=0; fil; fil=fil->next ) nfil++;
	
	return nfil;
}

/**
@brief 	Counts the number of filament nodes in a linked list of filaments.

	The function counts filament nodes in the list from the given pointer,
	and does not count filaments prior to this one.

@param 	*fil		pointer to any filament in the list.
@return long 				number of filament nodes.
**/
long 		filament_node_count(Bfilament* fil)
{
	long 		nfil;
	Bfilnode*	fnode;
	
	for ( nfil=0; fil; fil=fil->next )
		for ( fnode=fil->node; fnode; fnode=fnode->next ) nfil++;
	
	return nfil;
}

/**
@brief 	Finds the maximum selection number for a project.
@param 	*project	project parameter structure.
@return long 			maximum selection number.
**/
long		project_maximum_selection(Bproject* project)
{
	long				smax(0);
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
 	Bmicrograph*		mg;
	Bparticle*			part;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( smax < part->sel ) smax = part->sel;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( smax < part->sel ) smax = part->sel;
	}
	
	return smax;
}

/**
@brief 	Shows the selection numbers for the particles in a project.
@param 	*project	project parameter structure.
@return long 		number of particles selected.
**/
long		project_show_selection_numbers(Bproject* project)
{
	if ( !project ) return 0;
	
	long				i, maxsel, nsel(0);
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
 	Bmicrograph*		mg;
	Bparticle*			part;
	int*				sel = NULL;

	if ( field && field->mg && field->mg->part ) {
		for ( maxsel=0, field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( maxsel < part->sel ) maxsel = part->sel;
		maxsel++;
		sel = new int[maxsel];
		for ( i=0; i<maxsel; i++ ) sel[i] = 0;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( part->sel > 0 ) sel[part->sel]++;	
		cout << "Micrograph selection numbers:              ";
		for ( i=1; i<maxsel; i++ )
			if ( sel[i] ) cout << tab << i;
		cout << endl;
		delete[] sel;
	}
	
	if ( rec && rec->part ) {
		for ( maxsel=0, rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( maxsel < part->sel ) maxsel = part->sel;
		maxsel++;
		sel = new int[maxsel];
		for ( i=0; i<maxsel; i++ ) sel[i] = 0;
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( part->sel > 0 ) sel[part->sel]++;
		cout << "Reconstruction selection numbers:              ";
		for ( i=1; i<maxsel; i++ )
			if ( sel[i] ) cout << tab << i;
		cout << endl;
		delete[] sel;
	}
	
	return nsel;
}

/**
@brief 	Shows the selected particles in a project.
@param 	*project	project parameter structure.
@return long 				number of particles selected.
**/
long		project_show_selected(Bproject* project)
{
	if ( !project ) return 0;
	
	long				i, f, smax(0);
	long				nmg(0), nrec(0), ntotsel(0);
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;

	double				pp(0), fomavg[NFOM], fomstd[NFOM], fommin[NFOM], fommax[NFOM];

	long				npart = project_count_mg_particles(project);

	if ( npart && field && field->mg ) {
		nmg = project_count_mg_selected(project);
		for ( smax=0, field = project->field; field; field = field->next )
			if ( field->select )
				for ( mg = field->mg; mg; mg = mg->next )
					if ( mg->select )
						for ( part = mg->part; part; part = part->next )
							if ( smax < part->sel ) smax = part->sel;
		if ( verbose ) {
			cout << "Micrographs selected:           " << nmg << endl;
			cout << "Maximum selection number:       " << smax << endl;
		}
		vector<long>		nsel(smax+1, 0);
		vector<long>		ntot(smax+1, 0);
		vector<long>		npercent(smax+1, 0);
		vector<double>		fom(smax+1, 0);
		for ( f=0; f<NFOM; f++ ) {
			fomavg[f] = fomstd[f] = 0;
			fommin[f] = 1;
			fommax[f] = -1;
		}
		if ( verbose & VERB_PROCESS ) {
			cout << "Field\tMgraph";
			for ( i=0; i<=smax; i++ ) cout << tab << i;
			cout << endl;
		}
		for ( field = project->field; field; field = field->next ) if ( field->select ) {
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
				for ( i=0; i<=smax; i++ ) nsel[i] = 0;
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel > -1 ) {
						nsel[part->sel]++;
						fom[part->sel] += part->fom[0];
						if ( part->sel > 0 ) {
							ntotsel++;
							for ( f=0; f<NFOM; f++ ) {
								fomavg[f] += part->fom[f];
								fomstd[f] += part->fom[f]*part->fom[f];
								if ( fommin[f] > part->fom[f] ) fommin[f] = part->fom[f];
								if ( fommax[f] < part->fom[f] ) fommax[f] = part->fom[f];
							}
						}
					}
				}
				for ( i=0; i<=smax; i++ ) ntot[i] += nsel[i];
				if ( verbose & VERB_PROCESS ) {
					cout << field->id << tab << mg->id;
					for ( i=0; i<=smax; i++ ) cout << tab << nsel[i];
					cout << endl;
				}
			}
		}
		for ( i=0; i<=smax; i++ ) {
			if ( ntot[i] ) fom[i] /= ntot[i];
			if ( npart ) npercent[i] = ntot[i]*100.0/npart;
		}
		if ( ntotsel ) {
			for ( f=0; f<NFOM; f++ ) {
				fomavg[f] /= ntotsel;
				fomstd[f] = fomstd[f]/ntotsel - fomavg[f]*fomavg[f];
				if ( fomstd[f] > 0 ) fomstd[f] = sqrt(fomstd[f]);
				else fomstd[f] = 0;
			}
		}
		if ( verbose & VERB_PROCESS ) {
			cout << "Totals\t";
			for ( i=0; i<=smax; i++ ) cout << tab << ntot[i];
			cout << endl << "Percent\t";
			for ( i=0; i<=smax; i++ ) cout << tab << npercent[i];
			cout << endl << "FOM\t";
			for ( i=0; i<=smax; i++ ) cout << tab << fom[i];
			cout << endl << endl;
		}
		if ( npart ) pp = ntotsel*100.0/npart;
		cout << "Overall selected FOM average and standard deviation:" << endl;
		cout << "Selected particle images:       " << ntotsel << " (" << pp << "%)" << endl;
		for ( f=0; f<NFOM; f++ ) if ( project->fom_tag[f] )
			cout << "FOM" << f << ":                           " << 
				fommin[f] << tab << fommax[f] << tab << fomavg[f] << tab << fomstd[f] << endl;
		cout << endl;
	}

	npart = project_count_rec_particles(project);

	if ( npart && rec && rec->part ) {
		nrec = project_count_rec_selected(project);
		for ( smax=0, rec = project->rec; rec; rec = rec->next )
			if ( rec->select )
				for ( part = rec->part; part; part = part->next )
					if ( smax < part->sel ) smax = part->sel;
		if ( verbose ) {
			cout << "Reconstructions selected:         " << nrec << endl;
			cout << "Maximum particle selection number:" << smax << endl;
		}
		vector<int>		nsel(smax+1, 0);
		vector<int>		ntot(smax+1, 0);
		vector<int>		npercent(smax+1, 0);
		vector<int>		fom(smax+1, 0);
		for ( f=0; f<NFOM; f++ ) {
			fomavg[f] = fomstd[f] = 0;
			fommin[f] = 1;
			fommax[f] = -1;
		}
		if ( verbose & VERB_PROCESS ) {
			cout << "Recons";
			for ( i=0; i<=smax; i++ ) cout << tab << i;
			cout << endl;
		}
		for ( rec = project->rec; rec; rec = rec->next ) if ( rec->select ) {
			for ( i=0; i<=smax; i++ ) nsel[i] = 0;
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel > -1 ) {
					nsel[part->sel]++;
					fom[part->sel] += part->fom[0];
					if ( part->sel > 0 ) {
						ntotsel++;
						for ( f=0; f<NFOM; f++ ) {
							fomavg[f] += part->fom[f];
							fomstd[f] += part->fom[f]*part->fom[f];
							if ( fommin[f] > part->fom[f] ) fommin[f] = part->fom[f];
							if ( fommax[f] < part->fom[f] ) fommax[f] = part->fom[f];
						}
					}
				}
			}
			for ( i=0; i<=smax; i++ ) ntot[i] += nsel[i];
			if ( verbose & VERB_PROCESS ) {
				cout << rec->id;
				for ( i=0; i<=smax; i++ ) cout << tab << nsel[i];
				cout << endl;
			}
		}
		for ( i=0; i<=smax; i++ ) {
			if ( ntot[i] ) fom[i] /= ntot[i];
			if ( npart ) npercent[i] = ntot[i]*100.0/npart;
		}
		if ( ntotsel ) {
			for ( f=0; f<NFOM; f++ ) {
				fomavg[f] /= ntotsel;
				fomstd[f] = fomstd[f]/ntotsel - fomavg[f]*fomavg[f];
				if ( fomstd[f] > 0 ) fomstd[f] = sqrt(fomstd[f]);
				else fomstd[f] = 0;
			}
		}
		if ( verbose & VERB_PROCESS ) {
			cout << "Totals";
			for ( i=0; i<=smax; i++ ) cout << tab << ntot[i];
			cout << endl << "Percent";
			for ( i=0; i<=smax; i++ ) cout << tab << npercent[i];
			cout << endl << "FOM";
			for ( i=0; i<=smax; i++ ) cout << tab << fom[i];
			cout << endl << endl;
		}
		pp = 0;
		if ( npart ) pp = ntotsel*100.0/npart;
		cout << "Overall selected FOM average and standard deviation:" << endl;
		cout << "Selected particle images:       " << ntotsel << " (" << pp << " %)" << endl;
		for ( f=0; f<NFOM; f++ ) if ( project->fom_tag[f] )
			cout << "FOM" << f << ":                           " << 
				fommin[f] << tab << fommax[f] << tab << fomavg[f] << tab << fomstd[f] << endl;
		cout <<endl;
	}
	
	if ( project->class_avg ) {
		cout << "Class averages:                 " << particle_count(project->class_avg) << endl << endl;
	}
	
	return ntotsel;
}

/**
@brief 	Shows the selected particle parameters in a project.
@param 	*project	project parameter structure.
@param 	show		selection number to show.
@return long 		number of particles selected.
**/
long		project_show_selected_parameters(Bproject* project, int show)
{
	if ( !project ) return 0;
	
	if ( show < 0 ) return 0;
	
	int					f;
	long				nsel(0), npart(0);
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;

	if ( field && field->mg && field->mg->part ) {
		for ( field = project->field; field; field = field->next ) if ( field->select ) {
			cout << "Field: " << field->id << endl;
			for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
				cout << "Micrograph: " << mg->id << endl;
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel == show ) {
						cout << part->id << tab << part->ori[0] << tab
							<< part->ori[1] << tab << part->view[0] << tab
							<< part->view[1] << tab << part->view[2] << tab
							<< part->view.angle()*180/M_PI;
						for ( f=0; f<NFOM; f++ ) cout << tab << part->fom[f];
						cout << endl;
						nsel++;
					}
					npart++;
				}
				cout << endl;
			}
		}
		cout << "Particles selected:             " << nsel << " (" << nsel*100.0/npart << " %)" << endl << endl;
	}
	
	if ( rec && rec->part ) {
		for ( npart=nsel=0, rec = project->rec; rec; rec = rec->next ) if ( rec->select ) {
			cout << "Reconstruction: " << rec->id << endl;
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel == show ) {
					cout << part->id << tab << part->ori[0] << tab
						<< part->ori[1] << tab << part->view[0] << tab
						<< part->view[1] << tab << part->view[2] << tab
						<<  part->view.angle()*180/M_PI;
					for ( f=0; f<NFOM; f++ ) cout << tab << part->fom[f];
					cout << endl;
					nsel++;
				}
				npart++;
			}
			cout << endl;
		}
		cout << "Particles selected:             " << nsel << " (" << nsel*100.0/npart << " %)" << endl << endl;
	}
	
	return nsel;
}

int			mg_show_parameter(Bmicrograph* mg, Bstring& tag)
{
	cout << mg->id << tab;
	
	if ( tag == MICROGRAPH_FILE ) cout << mg->fmg;
	else if ( tag == MICROGRAPH_FRAMES_FILE ) cout << mg->fframe;
	else if ( tag == PARTICLE_FILE ) cout << mg->fpart;
	else if ( tag == MICROGRAPH_FILAMENT_FILE ) cout << mg->ffil;
	else if ( tag == MICROGRAPH_TRANSFORM_FILE ) cout << mg->fft;
	else if ( tag == MICROGRAPH_POWERSPEC_FILE ) cout << mg->fps;
	else if ( tag == MICROGRAPH_NUMBER ) cout << mg->img_num;
	else if ( tag == MICROGRAPH_SELECT ) cout << mg->select;
	else if ( tag == MICROGRAPH_FOM ) cout << mg->fom;
	else if ( tag == MICROGRAPH_PIXEL_X ) cout << mg->pixel_size[0];
	else if ( tag == MICROGRAPH_PIXEL_Y ) cout << mg->pixel_size[1];
	else if ( tag == MICROGRAPH_DOSE ) cout << mg->dose;
	else if ( tag == MICROGRAPH_WATER_RING ) cout << mg->wri;
	else if ( tag == MICROGRAPH_ORIGIN_X ) cout << mg->origin[0];
	else if ( tag == MICROGRAPH_ORIGIN_Y ) cout << mg->origin[1];
	else if ( tag == MICROGRAPH_SCALE_X ) cout << mg->scale[0];
	else if ( tag == MICROGRAPH_SCALE_Y ) cout << mg->scale[1];
	else if ( tag == MICROGRAPH_TILT_AXIS ) cout << mg->tilt_axis;
	else if ( tag == MICROGRAPH_TILT_ANGLE) cout << mg->tilt_angle;
	else if ( tag == MICROGRAPH_LEVEL_ANGLE) cout << mg->level_angle;
	else if ( tag == MICROGRAPH_ROT_ANGLE ) cout << mg->rot_angle;
	cout << endl;

	return 0;
}

int			part_show_parameter(Bparticle* part, Bstring& tag)
{
	int				float_flag(1), idata(0);
	float			fdata(0);
	Euler			euler;
	
	euler = Euler(part->view);
	if ( tag == PARTICLE_X ) fdata = part->loc[0];
	else if ( tag == PARTICLE_Y ) fdata = part->loc[1];
	else if ( tag == PARTICLE_Z ) fdata = part->loc[2];
	else if ( tag == PARTICLE_PIXEL ) fdata = part->pixel_size[0];
	else if ( tag == PARTICLE_PIXEL_X ) fdata = part->pixel_size[0];
	else if ( tag == PARTICLE_PIXEL_Y ) fdata = part->pixel_size[1];
	else if ( tag == PARTICLE_PIXEL_Z ) fdata = part->pixel_size[2];
	else if ( tag == PARTICLE_ORIGIN_X ) fdata = part->ori[0];
	else if ( tag == PARTICLE_ORIGIN_Y ) fdata = part->ori[1];
	else if ( tag == PARTICLE_ORIGIN_Z ) fdata = part->ori[2];
	else if ( tag == PARTICLE_DEFOCUS ) fdata = part->def;
	else if ( tag == PARTICLE_DEF_DEV ) fdata = part->dev;
	else if ( tag == PARTICLE_AST_ANG ) fdata = part->ast;
	else if ( tag == PARTICLE_MAGNIF ) fdata = part->mag;
	else if ( tag == PARTICLE_VIEW_X ) fdata = part->view[0];
	else if ( tag == PARTICLE_VIEW_Y ) fdata = part->view[1];
	else if ( tag == PARTICLE_VIEW_Z ) fdata = part->view[2];
	else if ( tag == PARTICLE_VIEW_ANGLE ) fdata = part->view.angle();
	else if ( tag == PARTICLE_PSI ) fdata = euler.psi()*180/M_PI;
	else if ( tag == PARTICLE_THETA ) fdata = euler.theta()*180/M_PI;
	else if ( tag == PARTICLE_PHI ) fdata = euler.phi()*180/M_PI;
	else if ( tag == PARTICLE_OMEGA ) fdata = -euler.psi()*180/M_PI;
	else if ( tag == PARTICLE_FOM ) fdata = part->fom[0];
	else if ( tag == PARTICLE_FOM_CV ) fdata = part->fom[1];
	else if ( tag == PARTICLE_HANDA_FOM ) fdata = part->fom[1];
	else if ( tag == PARTICLE_HANDB_FOM ) fdata = part->fom[2];
	else if ( tag == PARTICLE_SELECT ) {
		idata = part->sel;
		float_flag = 0;
	} else float_flag = -1;
	
	if ( float_flag == 0 ) {
		cout << part->id << tab << idata << endl;
	} else {
		cout << part->id << tab << fdata << endl;
	}

	return 0;
}

/**
@brief 	Shows the micrograph parameter indicated by the tag in a project.
@param 	*project	project parameter structure.
@param 	&tag		parameter tag.
@return long 		number of micrographs selected.
**/
long		project_show_mg_parameter(Bproject* project, Bstring& tag)
{
	if ( !project ) return 0;
	
	if ( tag.length() < 1 ) return 0;
	
	long				nsel(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg;

	cout << "MGID" << tab << tag << endl;
	
	for ( field = project->field; field; field = field->next ) if ( field->select ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			mg_show_parameter(mg, tag);
			nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Shows the particle parameter indicated by the tag in a project.
@param 	*project	project parameter structure.
@param 	&tag		parameter tag.
@return long 		number of particles selected.
**/
long		project_show_part_parameter(Bproject* project, Bstring& tag)
{
	if ( !project ) return 0;
	
	if ( tag.length() < 1 ) return 0;
	
	long				nsel(0), npart(0);
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;

	cout << "PID" << tab << tag << endl;
	
	if ( field && field->mg && field->mg->part ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel ) {
						part_show_parameter(part, tag);
						nsel++;
					}
					npart++;
				}
			}
		}
	}
	
	if ( rec && rec->part ) {
		for ( npart=nsel=0, rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel ) {
					part_show_parameter(part, tag);
					nsel++;
				}
				npart++;
			}
		}
	}
	
	return nsel;
}

/**
@brief 	Shows a histogram of the FOM values.
@param 	*project	project parameter structure.
@param	bins		number of bins, if 0, set to 100.
@param	min			FOM for the first bin.
@param	max			FOM for the last bin.
@return long 			number of particles selected.
**/
long		project_show_fom_histogram(Bproject* project, long bins,
				double min, double max)
{
	if ( !project ) return 0;
	if ( bins < 2 ) bins = 100;
	if ( max < min ) swap(min, max);
	if ( max == 0 ) max = 1;
	
	int					f;
	long				i, nh(bins-1), sum, ssum, nsel(0), npart;
	double				scale(nh*1.0/(max - min));
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;

	vector<long>		hist(bins);
	vector<long>		hsel(bins);

	if ( field && field->mg && field->mg->part )
		cout << "Micrograph figures-of-merit:" << endl;
	for ( f=0; f<NFOM; f++ ) if ( project->fom_tag[f] ) {
		nsel = npart = 0;
		for ( i=0; i<=nh; i++ ) hist[i] = hsel[i] = 0;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					i = (long) ((part->fom[f]-min)*scale + 0.5);
					if ( i < 0 ) i = 0;
					if ( i > nh ) i = nh;
					hist[i]++;
					npart++;
					if ( part->sel ) {
						hsel[i]++;
						nsel++;
					}
				}
			}
		}
		if ( nsel ) {
			cout << endl << "FOM" << f << " histogram:" << endl;
			cout << "  Bin\tFOM\tCount\t%%\tSum\tSel\t%%\tSelSum" << endl;
			sum = npart;
			ssum = nsel;
			for ( i=0; i<=nh; i++ ) {
				cout << i << tab << i*1.0/scale + min
					<< tab << hist[i] << tab
					<< hist[i]*100.0/npart << tab
					<< sum << tab << hsel[i] << tab
					<< hsel[i]*100.0/nsel << tab << ssum << endl;
				sum -= hist[i];
				ssum -= hsel[i];
			}
			cout << endl;
		}
	}
	
	if ( rec && rec->part )
		cout << "Reconstruction figures-of-merit:" << endl;
	for ( f=0; f<NFOM; f++ ) if ( project->fom_tag[f] ) {
		nsel = npart = 0;
		for ( i=0; i<=nh; i++ ) hist[i] = hsel[i] = 0;
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				i = (long) ((part->fom[f]-min)*scale + 0.5);
				if ( i < 0 ) i = 0;
				if ( i > nh ) i = nh;
				hist[i]++;
				npart++;
				if ( part->sel ) {
					hsel[i]++;
					nsel++;
				}
			}
		}
		if ( nsel ) {
			cout << endl << "FOM" << f << " histogram:" << endl;
			cout << "  Bin\tFOM\tCount\t%%\tSum\tSel\t%%\tSelSum" << endl;
			sum = npart;
			ssum = nsel;
			for ( i=0; i<=nh; i++ ) {
				cout << i << tab << i*1.0/scale + min
					<< tab << hist[i] << tab
					<< hist[i]*100.0/npart << tab
					<< sum << tab << hsel[i] << tab
					<< hsel[i]*100.0/nsel << tab << ssum << endl;
				sum -= hist[i];
				ssum -= hsel[i];
			}
			cout << endl;
		}
	}
	
	return nsel;
}

/**
@brief 	Shows the histogram for magnification values for selected particles.
@param 	*project		project parameter structure.
@param 	bins			number of bins.
@param 	increment		increment between bins.
@return long 				number of particles selected.
**/
long		project_show_mag_histogram(Bproject* project, long bins, double increment)
{
	if ( !project ) return 0;
	
	int				i, nsel(0);
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	double			min = 1 - (bins - 1)*increment/2;
//	int*			hist = new int[bins];
//	for ( i=0; i<bins; i++ ) hist[i] = 0;
	vector<long>	hist(bins,0);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel ) {
					i = (int) ((part->mag - min)/increment + 0.5);
					if ( i < 0 ) i = 0;
					if ( i >=bins ) i = bins - 1;
					hist[i]++;
					nsel++;
				}
			}
		}
	}

	cout << "Mag\tCount" << endl;
	
	for ( i=0; i<bins; i++ )
		cout << min + i*increment << tab << hist[i] << endl;
	
	cout << endl << "Total\t" << nsel << endl << endl;
	
//	delete[] hist;
	
	return nsel;
}

/**
@brief 	Selects a field from the project and deletes the rest.
@param 	*project 	project parameter structure.
@param 	&field_id	field to select.
@return int					0.
**/
int			project_select_field(Bproject* project, Bstring& field_id)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bfield* 		selected_field = NULL;

	if ( project->field->id == field_id ) {
		selected_field = project->field;
	} else {
		for ( field=project->field; field && !selected_field; field=field->next ) {
			if ( field->next->id == field_id ) {
				selected_field = field->next;
				field->next = NULL;
			}
		}
		field_kill(project->field);
	}
	
	field_kill(selected_field->next);
	selected_field->next = NULL;
	project->field = selected_field;
	
	return 0;
}

/**
@brief 	Selects a micrograph from the project and deletes the rest.
@param 	*project 	project parameter structure.
@param 	&mg_id		micrograph to select.
@return int			0.
**/
int			project_select_micrograph(Bproject* project, Bstring& mg_id)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	Bmicrograph*	selected_mg = NULL;

	for ( field=project->field; field && !selected_mg; field=field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mg->id == mg_id ) {
				selected_mg = mg;
				project_select_field(project, field->id);
			}
	
	field=project->field;
	if ( field->mg->id == mg_id ) {
		selected_mg = field->mg;
	} else {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg->next->id == mg_id ) {
				selected_mg = mg->next;
				mg->next = NULL;
			}
		}
		micrograph_kill(field->mg);
	}

	micrograph_kill(selected_mg->next);
	selected_mg->next = NULL;
	field->mg = selected_mg;
	
	return 0;
}

/**
@brief 	Selects micrographs with particles.
@param 	*project 		project parameter structure.
@param	part_sel		flag to select only with this selection number
@return long			number of particles.
**/
long		project_select_with_particles(Bproject* project, long part_sel)
{
	if ( !project ) return 0;
	
	long			np(project_count_mg_particles(project));
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;

	for ( field=project->field; field; field=field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			mg->select = 0;
			if ( part_sel < -1 ) {
				if ( mg->part ) mg->select = 1;
			} else if ( part_sel < 0 ) {
				for ( part = mg->part; part; part = part->next )
					if ( part->sel ) mg->select = 1;
			} else {
				for ( part = mg->part; part; part = part->next )
					if ( part->sel == part_sel ) mg->select = 1;
			}
		}
	
	return np;
}

/**
@brief 	Sorts the micrographs in a field by tilt angle.
@param 	*field 				field parameter structure.
@return int					0.
**/
int			field_sort_by_tilt(Bfield* field)
{
	if ( !field ) return 0;
	
	long				i(0), nmg = field_count_micrographs(field);
	double				minang;
	
	vector<Bmicrograph*>	mgarr(nmg);
	
	Bmicrograph*		mg, *mg2;
	Bmicrograph*		mgsel = NULL;

	for ( mg = field->mg; mg; mg = mg->next ) mg->select = 1;

	for ( mg = field->mg; mg; mg = mg->next ) {
		minang = 100;
		for ( mg2 = field->mg; mg2; mg2 = mg2->next ) if ( mg2->select ) {
			if ( minang > mg2->tilt_angle ) {
				minang = mg2->tilt_angle;
				mgsel = mg2;
			}
		}
		mgarr[i++] = mgsel;
		mgsel->select = 0;
	}
	
	field->mg = mg = mgarr[0];
	
	for ( i=1; i<nmg; ++i, mg = mg->next ) mg->next = mgarr[i];
	mg->next = NULL;

	for ( mg = field->mg; mg; mg = mg->next ) mg->select = 1;
	
	return 0;
}

/**
@brief 	Sorts the micrographs in each field by tilt angle.
@param 	*project 			project parameter structure.
@return int					0.
**/
int			project_sort_by_tilt(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield* 		field;

	for ( field = project->field; field; field = field->next )
		field_sort_by_tilt(field);
	
	return 0;
}

/**
@brief 	Sorts the micrographs by a selected parameter.
@param 	*project 		project parameter structure.
@param 	tag 				parameter tag.
@return vector<pair<Bmicrograph*,double>>	array of micrograph links and values.
**/
vector<pair<Bmicrograph*,double>>	project_mg_sort(Bproject* project, Bstring tag)
{
	long				i(0);
	Bfield*				field = project->field;
	Bmicrograph*		mg = field->mg;
	vector<std::pair<Bmicrograph*,double>> 	mg_value;
	pair<Bmicrograph*,double>				mgvp;

	if ( !mg ) {
		cerr << "Error: Micrograph parameters not defined!" << endl;
		return mg_value;
	}
	
	if ( tag == CTF_DEF_AVG && !mg->ctf ) {
		cerr << "Error: CTF parameters not defined!" << endl;
		return mg_value;
	}
	
	if ( verbose )
		cout << "Ordering micrographs based on tag " << tag << endl;
	
	
	for ( field = project->field; field; field = field->next ) {
		 for ( mg = field->mg; mg; mg = mg->next ) {
			i++;
			mgvp.first = mg;
			if ( tag == MICROGRAPH_FOM ) mgvp.second = mg->fom;
			else if ( tag == MICROGRAPH_INTENSITY )
				mgvp.second = mg->intensity;
			else if ( tag == CTF_DEF_AVG )
				mgvp.second = mg->ctf->defocus_average();
			else if ( tag == MICROGRAPH_TILT_ANGLE )
				mgvp.second = mg->tilt_angle;
			else if ( tag == MICROGRAPH_WATER_RING )
				mgvp.second = mg->wri;
			else if ( tag == PARTICLE )
				mgvp.second = particle_count(mg->part);
			else
				mgvp.second = i;
			mg_value.push_back(mgvp);
		}
	}
	
	if ( tag.length() ) {
		sort(mg_value.begin(), mg_value.end(),
			[=](pair<Bmicrograph*,double> a, pair<Bmicrograph*,double> b) {
				return a.second < b.second;
			}
		);
	}

	return mg_value;
}
