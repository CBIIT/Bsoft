/**
@file	mg_marker.cpp
@brief	Functions for dealing with markers and linkers
@author Bernard Heymann
@date	Created: 20150212
@date	Modified: 20150212
**/

#include "mg_marker.h"
#include "rwimg.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Lists micrograph markers.
@param 	*project	micrograph project.
@return long		number of markers (<0 means failure).
**/
long		project_marker_lists(Bproject* project)
{
	long		ntot(0), ntotsel(0), n, nsel;
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bmarker*			mark;

	if ( field && field->mg && field->mg->mark ) {
		cout << "Field\tMicrograph\tMarkers\tSelected" << endl;
		for ( field = project->field; field; field = field->next ) {
			n = nsel = 0;
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( mark = mg->mark; mark; mark = mark->next ) {
					n++;
					if ( mark->sel ) nsel++;
				}
				cout << field->id << tab << mg->id << tab << n << tab << nsel << endl;
			}
			ntot += n;
			ntotsel += nsel;
		}
	}
	
	if ( rec && rec->mark ) {
		cout << "Reconstruction\tMarkers\tSelected" << endl;
		n = nsel = 0;
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( mark = rec->mark; mark; mark = mark->next ) {
				n++;
				if ( mark->sel ) nsel++;
			}
			cout << rec->id << tab << n << tab << nsel << endl;
		}
		ntot += n;
		ntotsel += nsel;
	}
	
	cout << "Total number of markers:      " << ntot << tab << ntotsel << endl;

	return ntot;
}


/**
@brief 	Enumerates markers within particles.
@param 	*project	micrograph project.
@return long		number of markers (<0 means failure).
**/
long		project_marker_in_particle(Bproject* project)
{
	long		i, n(0);
	double				maxrad, d;
	int					h[100];
	for ( i=0; i<100; i++ ) h[i] = 0;
	
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;
	Bmarker*			mark;

	if ( !project->select && field && field->mg && field->mg->mark ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				maxrad = mg->box_size[0]/2;
				for ( part = mg->part; part; part = part->next ) {
					part->sel = 0;
					for ( mark = mg->mark; mark; mark = mark->next ) {
						d = part->loc.distance(mark->loc);
						if ( d < maxrad ) {
							part->sel++;
						}
					}
					h[part->sel]++;
					n++;
				}
			}
		}
	}
	
	if ( project->select && rec && rec->mark ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			maxrad = rec->box_size[0]/2;
			for ( part = rec->part; part; part = part->next ) {
				part->sel = 0;
				for ( mark = rec->mark; mark; mark = mark->next ) {
					d = part->loc.distance(mark->loc);
					if ( d < maxrad ) {
						part->sel++;
					}
				}
				h[part->sel]++;
				n++;
			}
		}
	}
	
	cout << "Number\tCount\t%" << endl;
	for ( i=0; i<100; i++ ) if ( h[i] )
		cout << i << tab << h[i] << tab << fixed << setprecision(2) << h[i]*100.0/n << endl;
	cout << "Total" << tab << n << tab << 100 << endl;

	return n;
}

/**
@brief 	Enumerates markers within particles.
@param 	*project			micrograph project.
@param 	marker_radius		search radius to associate with marker.
@return long		number of markers (<0 means failure).
**/
long		project_marker_in_particle_image(Bproject* project, double marker_radius)
{
	long		i, n(0), np, fill_type(FILL_USER);
	double				maxrad, d;
	Vector3<float>		loc, box_rad;
	int					h[100];
	for ( i=0; i<100; i++ ) h[i] = 0;
	
	Bfield*				field = project->field;
	Breconstruction*	rec = project->rec;
	Bmicrograph*		mg;
	Bparticle*			part;
	Bmarker*			mark;
	Bimage*				p;
	Bimage*				pt;

	if ( !project->select && field && field->mg && field->mg->mark ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				maxrad = mg->box_size[0]/2;
				box_rad = mg->box_size/2;
				mg->mark_radius = marker_radius;
				for ( np=0, part = mg->part; part; part = part->next ) np++;
				if ( np ) {
					p = new Bimage(Float, TSimple, mg->box_size, np);
					pt = new Bimage(Float, TSimple, mg->box_size, 1);
					for ( i=0, part = mg->part; part; part = part->next, i++ ) {
						pt->clear();
						part->sel = 0;
						for ( mark = mg->mark; mark; mark = mark->next ) {
							d = part->loc.distance(mark->loc);
							if ( d < maxrad ) {
								part->sel++;
								loc = mark->loc - part->loc + part->ori;
								pt->sphere(loc, marker_radius, 2, fill_type, 1);
							}
						}
						h[part->sel]++;
						n++;
						p->replace(i, pt);
					}
					mg->fpart = mg->fpart.pre_rev('.') + "_mark." + mg->fpart.post_rev('.');
					cout << mg->fpart << endl;
					write_img(mg->fpart, p, 0);
					delete p;
				}
			}
		}
	}
	
	if ( project->select && rec && rec->mark ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			maxrad = rec->box_size[0]/2;
			for ( part = rec->part; part; part = part->next ) {
				part->sel = 0;
				for ( mark = rec->mark; mark; mark = mark->next ) {
					d = part->loc.distance(mark->loc);
					if ( d < maxrad ) {
						part->sel++;
					}
				}
				h[part->sel]++;
				n++;
			}
		}
	}
	
	cout << "Number\tCount\t%" << endl;
	for ( i=0; i<100; i++ ) if ( h[i] )
		cout << i << tab << h[i] << tab << fixed << setprecision(2) << h[i]*100.0/n << endl;
	cout << "Total" << tab << n << tab << 100 << endl;

	return n;
}

