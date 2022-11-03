/**
@file	mg_particle_select.cpp
@brief	Select particles
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220623
**/

#include "mg_processing.h"
#include "mg_tomography.h"
#include "mg_select.h"
#include "mg_particle_select.h"
#include "mg_tags.h"
#include "symmetry.h"
#include "matrix_linear.h"
#include "Matrix.h"
#include "random_numbers.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int_float* 	project_fom_order(Bproject* project, long& npart, int fom_index, int defocus_fit)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;
	
	long				nn(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select < 1 ) npart = project_count_mg_particles(project);
	else npart = project_count_rec_particles(project);

	if ( npart < 1 ) {
		cerr << "Error: No particles found!" << endl;
		return NULL;
	}
	
	double				slope, intercept, fom_adjust(0);
		
	if ( verbose & VERB_FULL )
		cout << "Ordering " << npart << " particles based on FOM:" << endl;
	
	int_float*			rank = new int_float[npart];
	
	if ( project->select < 1 ) {
		if ( defocus_fit ) part_fom_defocus_fit(project, fom_index, intercept, slope);
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( defocus_fit && mg->ctf ) fom_adjust = slope*(mg->ctf->defocus_average() - 1e4);
				for ( part = mg->part; part; part = part->next ) {
					rank[nn].i = nn;
					if ( part->sel > 0 )
//						rank[nn].f = fom_adjust + part->fom[fom_index];
						rank[nn].f = part->fom[fom_index] - fom_adjust;
					else
						rank[nn].f = 0;
					nn++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				rank[nn].i = nn;
				if ( part->sel > 0 )
					rank[nn].f = part->fom[fom_index];
				else
					rank[nn].f = 0;
				nn++;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_fom_order: Ranking the FOMs:" << endl;

	// Rank according to the FOMs
	qsort((void *) rank, npart, sizeof(int_float), 
			(int (*)(const void *, const void *)) QsortLargeToSmallIntFloat);
	
	return rank;
}

int_float* 	part_fom_order(Bparticle* partlist, long& npart, int fom_index)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;
	
	long				nn(0);
	
	Bparticle*			part;
	
	npart = particle_count(partlist);

	if ( npart < 1 ) {
		cerr << "Error: No particles found!" << endl;
		return NULL;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Ordering " << npart << " particles based on FOM:" << endl;
	
	int_float*			rank = new int_float[npart];
	
	for ( part = partlist; part; part = part->next ) if ( part->sel ) {
		rank[nn].i = nn;
		rank[nn].f = part->fom[fom_index];
		nn++;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_fom_order: Ranking the FOMs:" << endl;

	// Rank according to the FOMs
	qsort((void *) rank, npart, sizeof(int_float), 
			(int (*)(const void *, const void *)) QsortLargeToSmallIntFloat);
	
	return rank;
}

/**
@brief 	Resets selection to all particles.
@param 	*project		project parameter structure with all parameters.
@param	flag			flag to limit resetting: 1=all, 2=part.
@return long			number of particles.
**/
long		part_reset_selection(Bproject* project, int flag)
{
	long				npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		if ( verbose & VERB_PROCESS )
			cout << "Resetting reconstruction particle selections" << endl << endl;
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( flag & 1 ) rec->select = 1;
			if ( flag & 2 ) for ( part = rec->part; part; part = part->next, npart++ )
				part->sel = 1;
		}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "Resetting micrograph particle selections" << endl << endl;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( flag & 1 ) mg->select = 1;
				if ( flag & 2 ) for ( part = mg->part; part; part = part->next, npart++ )
					part->sel = 1;
			}
	}
	
	return npart;
}

/**
@brief 	Unsets selection to no particles.
@param 	*project		project parameter structure with all parameters.
@return long			number of particles.
**/
long		part_unset_selection(Bproject* project)
{
	long				npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		if ( verbose & VERB_PROCESS )
			cout << "Unsetting reconstruction particle selections" << endl << endl;
		for ( rec = project->rec; rec; rec = rec->next ) {
//			rec->select = 1;
			for ( part = rec->part; part; part = part->next, npart++ )
				part->sel = 0;
		}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "Unsetting micrograph particle selections" << endl << endl;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) {
//				mg->select = 1;
				for ( part = mg->part; part; part = part->next, npart++ )
					part->sel = 0;
			}
	}
	
	return npart;
}


/**
@brief 	Inverts the selection.
@param 	*project		project parameter structure with all parameters.
@return long			number of particles.
**/
long		part_invert_selection(Bproject* project)
{
	long				npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		if ( verbose & VERB_PROCESS )
			cout << "Inverts reconstruction particle selections" << endl << endl;
		for ( rec = project->rec; rec; rec = rec->next ) {
//			rec->select = 1;
			for ( part = rec->part; part; part = part->next, npart++ ) {
				if ( part->sel ) part->sel = 0;
				else part->sel = 1;
			}
		}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "Inverts micrograph particle selections" << endl << endl;
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) {
//				mg->select = 1;
				for ( part = mg->part; part; part = part->next, npart++ ) {
					if ( part->sel ) part->sel = 0;
					else part->sel = 1;
				}
			}
	}
	
	return npart;
}


/**
@brief 	Sets selection to the micrographs indicated.
@param 	*project		project parameter structure with all parameters.
@param 	&mgselect		string with selection.
@return long				number of particles.

	Only the micrograph selection fields are modified.
**/
long		part_select_micrograph(Bproject* project, Bstring& mgselect)
{
	long			i, nsel(0);
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	long			nmg = project_count_micrographs(project);
	
	vector<int>		numsel = select_numbers(mgselect, nmg+1);
	
	for ( i=1; i<=nmg; ++i ) if ( numsel[i] ) nsel++;
	
	if ( verbose )
		cout << "Micrographs selected:           " << nsel << " (" << nsel*100.0/nmg << " %)" << endl;
	
	for ( i=1, nsel=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, ++i ) {
			if ( numsel[i] ) {
				mg->select = 1;
				for ( part = mg->part; part; part = part->next )
					if ( part->sel ) nsel++;
			} else mg->select = 0;
		}
	}
		
	return nsel;
}

/**
@brief 	Sets selection to the micrographs with selected particles.
@param 	*project		project parameter structure with all parameters.
@return long				number of particles.

	Only the micrograph selection fields are modified.
**/
long		part_select_micrographs_with_selected_particles(Bproject* project)
{
	long			nsel(0);
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	long			nmg = project_count_micrographs(project);
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->select ) {
			mg->select = 0;
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel ) {
					mg->select = 1;
					break;
				}
			}
			nsel += mg->select;
		}
	}
		
	if ( verbose )
		cout << "Micrographs selected:           " << nsel << " (" << nsel*100.0/nmg << " %)" << endl;
	
	return nsel;
}

/**
@brief 	Consolidates selected particles under one selection number.
@param 	*project		project parameter structure with all parameters.
@param 	number				new selection number (1 if < 1).
@return long					number of particles.
**/
long		part_consolidate_selection(Bproject* project, int number)
{
	long				nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( part->sel > 0 ) {
					part->sel = number;
					nsel++;
				}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( part->sel > 0 ) {
						part->sel = number;
						nsel++;
					}
	}
	
	return nsel;
}

/**
@brief 	Selects one selection number and sets all others to zero.
@param 	*project	project parameter structure with all parameters.
@param 	number		selection number to keep.
@return long			number of particles selected.
**/
long		part_set_selection(Bproject* project, int number)
{
	long				nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( part->sel == number ) {
					nsel++;
				} else {
					part->sel = 0;
				}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( part->sel == number ) {
						nsel++;
					} else {
						part->sel = 0;
					}
	}
	
	return nsel;
}

/**
@brief 	Sets selection numbers sequentially.
@param 	*project	project parameter structure with all parameters.
@return long			number of particles selected.
**/
long		part_set_sequential(Bproject* project)
{
	long				cnt(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( part->sel ) part->sel = ++cnt;
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( part->sel ) part->sel = ++cnt;
	}
	
	return cnt;
}


/**
@brief 	Sets the FOM for all particles.
@param 	*project	project parameter structure with all parameters.
@param 	fom_index		index of FOM value to set.
@param 	fom			new FOM.
@return long				number of particles.
**/
long		part_set_FOM(Bproject* project, int fom_index, double fom)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;

	long				npart(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select )
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next, npart++ )
				part->fom[fom_index] = fom;
	else
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next, npart++ )
					part->fom[fom_index] = fom;
	
	return npart;
}

/**
@brief 	Deselects particles with selection numbers from a list.
@param 	*project	parameter structure with all parameters.
@param 	list		comma-separated list of selection numbers.
@return long		number of particles selected.
**/
long		part_deselect_from_list(Bproject* project, Bstring list)
{
	long				nsel = project_count_mg_part_selected(project);
	
	if ( list.length() < 1 ) return nsel;
	
	long				smax = project_maximum_selection(project);

	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	long				n(smax+1);
//	int*				numsel = new int[n];
//	select_numbers(list, n, numsel);
	vector<int>		numsel = select_numbers(list, n);
	
	nsel = 0;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					if ( numsel[part->sel] ) part->sel = 0;
					else nsel++;
				}
			}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel > 0 ) {
						if ( numsel[part->sel-1] ) part->sel = 0;
						else nsel++;
					}
				}
	}
	
//	delete[]			numsel;
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}


/**
@brief 	Deselects particles below a given FOM cutoff.
@param 	*project	parameter structure with all parameters.
@param 	fom_index	index of FOM value to test for.
@param 	fommin		minimum threshold for deselection.
@param 	fommax		maximum threshold for deselection.
@return long		number of particles selected.
**/
long		part_deselect(Bproject* project, int fom_index, double fommin, double fommax)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;

	long				nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				if ( part->fom[fom_index] < fommin || part->fom[fom_index] > fommax )
					part->sel = 0;
				if ( part->sel ) nsel++;
			}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					if ( part->fom[fom_index] < fommin || part->fom[fom_index] > fommax ) 
						part->sel = 0;
					if ( part->sel ) nsel++;
				}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "FOM range:                      " << fommin << " - " << fommax << endl;
		cout << "Particles selected:             " << nsel << endl << endl;
	}
	
	return nsel;
}

int			QsortLargeToSmallParticle(const void *x, const void *y)
{
	Bparticle** 		part1 = (Bparticle**) x;
	Bparticle** 		part2 = (Bparticle**) y;
	if( (*part1)->fom[0] > (*part2)->fom[0] ) return -1;
	else return 1;
}

/**
@brief 	Deselects particles overlapping with better ones.
@param 	*partlist	particle linked list.
@param 	excl_dist	minimum distance between particles.
@param 	part_select	initial selection number (-1 means all >0).
@param 	fom_index	index of FOM value to test for.
@return long		number of particles selected.
**/
long		part_deselect_redundant(Bparticle* partlist, double excl_dist, int part_select, int fom_index)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;

	long		i, j, npart(0), nsel(0);
	
	Bparticle**	parr = particle_array(partlist, part_select, npart);

	// Rank according to the FOMs
	qsort((void *) parr, npart, sizeof(Bparticle*), 
			(int (*)(const void *, const void *)) QsortLargeToSmallParticle);
	
	for ( i=0; i<npart; i++ ) {
		for ( j=i+1; j<npart; j++ ) if ( parr[j]->sel ) {
			if ( parr[i]->loc.distance(parr[j]->loc) < excl_dist ) 
				parr[j]->sel = 0;
		}
		if ( parr[i]->sel ) nsel++;
	}
	
	return nsel;
}

/**
@brief 	Deselects particles overlapping with better ones.
@param 	*project	parameter structure with all parameters.
@param 	excl_dist	minimum distance between particles.
@param 	part_select	initial selection number (-1 means all >0).
@param 	fom_index	index of FOM value to test for.
@return long			number of particles selected.
**/
long		part_deselect_redundant(Bproject* project, double excl_dist, int part_select, int fom_index)
{
	if ( fom_index < 0 ) fom_index = 0;
	if ( fom_index >= NFOM ) fom_index = NFOM - 1;

	long				nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			nsel +=	part_deselect_redundant(rec->part, excl_dist, part_select, fom_index);
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				nsel +=	part_deselect_redundant(mg->part, excl_dist, part_select, fom_index);
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Exclusion distance:             " << excl_dist << endl;
		cout << "Particles selected:             " << nsel << endl << endl;
	}
	
	return nsel;
}

/**
@brief 	Sets selection for generating multiple maps.
@param 	*project		parameter structure with all parameters.
@param 	part_select		initial selection number (-1 means all >0).
@param 	nmaps			desired number of maps.
@return long				number of particles selected.

	Selected particles are sequentially assigned increasing integers
	up to desired number of maps (nmaps).

**/
long		part_set_multi_maps(Bproject* project, int part_select, int nmaps)
{
	long				nsel(0);
	int					i = 1;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( nmaps < 1 ) nmaps = 1;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				if ( part_select < 0 ) {
					if ( part->sel > 0 ) part->sel = i;
				} else {
					if ( part->sel == part_select ) part->sel = i;
					else part->sel = 0;
				}
				if ( part->sel ) {
					nsel++;
					i++;
					if ( i > nmaps ) i = 1;
				}
			}
	} else {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					if ( part_select < 0 ) {
						if ( part->sel > 0 ) part->sel = i;
					} else {
						if ( part->sel == part_select ) part->sel = i;
						else part->sel = 0;
					}
					if ( part->sel ) {
						nsel++;
						i++;
						if ( i > nmaps ) i = 1;
					}
				}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Generates a selection number from the group number for each filament.
@param 	*project		parameter structure with all parameters.
@return long				number of maps selected.
**/
long		part_set_filament_maps(Bproject* project)
{
	long				nsel(0);
	int					prev_group(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					if ( part->group != prev_group ) {
						nsel++;
						prev_group = part->group;
					}
					part->sel = nsel;
				}
			}
			prev_group = 0;
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel > 0 ) {
						if ( part->group != prev_group ) {
							nsel++;
							prev_group = part->group;
						}
						part->sel = nsel;
					}
				}
				prev_group = 0;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Maps selected:                  " << nsel << endl << endl;
	
	return nsel;
}


int			part_reselect(Bparticle* part, Bstring& tag, double reselect_min, double reselect_max)
{
	int				float_flag(1), idata(0);
	double			fdata;
	Euler			euler;
	
	fdata = reselect_min - 1;
	euler = Euler(part->view);
	if ( tag == PARTICLE_ID ) {
		idata = part->id;
		float_flag = 0;
	}
	if ( tag == PARTICLE_X ) fdata = part->loc[0];
	if ( tag == PARTICLE_Y ) fdata = part->loc[1];
	if ( tag == PARTICLE_Z ) fdata = part->loc[2];
	if ( tag == PARTICLE_ORIGIN_X ) fdata = part->ori[0];
	if ( tag == PARTICLE_ORIGIN_Y ) fdata = part->ori[1];
	if ( tag == PARTICLE_ORIGIN_Z ) fdata = part->ori[2];
	if ( tag == PARTICLE_DEFOCUS ) fdata = part->def;
	if ( tag == PARTICLE_DEF_DEV ) fdata = part->dev;
	if ( tag == PARTICLE_AST_ANG ) fdata = part->ast;
	if ( tag == PARTICLE_MAGNIF ) fdata = part->mag;
	if ( tag == PARTICLE_VIEW_X ) fdata = part->view[0];
	if ( tag == PARTICLE_VIEW_Y ) fdata = part->view[1];
	if ( tag == PARTICLE_VIEW_Z ) fdata = part->view[2];
	if ( tag == PARTICLE_VIEW_ANGLE ) fdata = part->view.angle();
	if ( tag == PARTICLE_PSI ) fdata = euler.psi()*180/M_PI;
	if ( tag == PARTICLE_THETA ) fdata = euler.theta()*180/M_PI;
	if ( tag == PARTICLE_PHI ) fdata = euler.phi()*180/M_PI;
	if ( tag == PARTICLE_OMEGA ) fdata = -euler.psi()*180/M_PI;
	if ( tag == PARTICLE_FOM ) fdata = part->fom[0];
	if ( tag == PARTICLE_FOM_CV ) fdata = part->fom[1];
	if ( tag == PARTICLE_HANDA_FOM ) fdata = part->fom[1];
	if ( tag == PARTICLE_HANDB_FOM ) fdata = part->fom[2];
	if ( tag == PARTICLE_SELECT ) {
		idata = part->sel;
		float_flag = 0;
	}
		
	if ( !float_flag ) {
		if ( idata < reselect_min || idata > reselect_max ) part->sel = 0;
	} else {
		if ( fdata < reselect_min || fdata > reselect_max ) part->sel = 0;
	}

	return 0;
}

/**
@brief 	Reselects particles identified by an input STAR tag and between a minimum and maximum.
@param 	*project		parameter structure with all parameters.
@param 	&tag			tag indicating which data to select on.
@param 	reselect_min		minimum value.
@param 	reselect_max		maximum value.
@return long					number of particles selected.

	Only particles already selected are subject to reselection.

**/
long		part_reselect(Bproject* project, Bstring& tag, double reselect_min, double reselect_max)
{	
	long				npart(0), nmg(0), nrec(0), npsel(0), nmsel(0), nrsel(0);
	double				fdata, mpct(0), rpct(0), ppct(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( tag[0] == '_' ) tag = tag.erase(0);
	
	if ( verbose ) {
		cout << "Reselecting on tag:             " << tag << endl;
		cout << "Range:                          " << reselect_min << " - " << reselect_max << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				fdata = reselect_min - 1;
				if ( mg->ctf && tag.contains(CTF) ) {
					if ( tag == CTF_DEF_AVG ) fdata = mg->ctf->defocus_average();
					if ( tag == CTF_DEF_DEV ) fdata = mg->ctf->defocus_deviation();
					if ( tag == CTF_AST_ANG ) fdata = mg->ctf->astigmatism_angle();
					if ( fdata < reselect_min || fdata > reselect_max ) mg->select = 0;
				}
				if ( tag.contains(PARTICLE) )
					for ( part = mg->part; part; part = part->next )
						part_reselect(part, tag, reselect_min, reselect_max);
				if ( mg->select ) nmsel++;
				else for ( part = mg->part; part; part = part->next ) part->sel = 0;
				for ( part = mg->part; part; part = part->next, npart++ ) if ( part->sel ) npsel++;
				nmg++;
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			fdata = reselect_min - 1;
			if ( rec->ctf && tag.contains(CTF) ) {
				if ( tag == CTF_DEF_AVG ) fdata = rec->ctf->defocus_average();
				if ( tag == CTF_DEF_DEV ) fdata = rec->ctf->defocus_deviation();
				if ( tag == CTF_AST_ANG ) fdata = rec->ctf->astigmatism_angle();
				if ( fdata < reselect_min || fdata > reselect_max ) rec->select = 0;
			}
			if ( tag.contains(PARTICLE) )
				for ( part = rec->part; part; part = part->next )
					part_reselect(part, tag, reselect_min, reselect_max);
			if ( rec->select ) nrsel++;
			else for ( part = rec->part; part; part = part->next ) part->sel = 0;
			for ( part = rec->part; part; part = part->next, npart++ ) if ( part->sel ) npsel++;
			nrec++;
		}
	}
	
	if ( npart ) ppct = npsel*100.0/npart;
	if ( nmg ) mpct = nmsel*100.0/nmg;
	if ( nrec ) rpct = nrsel*100.0/nrec;
	
	if ( verbose & VERB_PROCESS ) {
		if ( project->select < 1 )
			cout << "Micrographs reselected:         " << nmsel << " (" << mpct << " %)" << endl;
		else
			cout << "Reconstructions reselected:     " << nrsel << " (" << rpct << " %)" << endl;
		cout << "Particles reselected:           " << npsel << " (" << ppct << " %)" << endl << endl;
	}
	
	return npsel;
}

/**
@brief 	A function to select particles based on comparison of the orientations within a series.
@param 	*project		parameter structure with all parameters.
@param 	*sym			symmetry.
@param 	angle_cutoff		angle cutoff value (radians).
@return double					average angular deviation.

	Only particles already selected are subject to the test.

**/
double		part_series_comparison(Bproject* project, Bsymmetry& sym, double angle_cutoff)
{	
	int 			i, max(0), nsel(0);
	double			angle, min_angle, avg_angle(0);
	View			view_pred;
	View*			symview;
	Vector3<double> axis(1,0,0);
	Quaternion		qt, q;
	
	Bfield*			field;
	Bmicrograph*	mg, *mg0;
	Bparticle*		part, *part0;

	int				hist[180];
	for ( i=0; i<180; i++ ) hist[i] = 0;
	
	int 			nsymviews(sym.order());
	
	if ( verbose & VERB_LABEL )
		cout << "Comparison of particles in micrograph series with " << angle_cutoff*180.0/M_PI << " degrees cutoff" << endl;
	
	for ( field = project->field; field; field = field->next ) {
		mg0 = field->mg;
		axis[0] = cos(mg0->tilt_axis);
		axis[1] = sin(mg0->tilt_axis);
		for ( mg=field->mg->next; mg; mg=mg->next ) {
//			qt = quaternion_from_angle_and_axis3(mg0->tilt_angle - mg->tilt_angle, axis);
			qt = Quaternion(axis, mg0->tilt_angle - mg->tilt_angle);
			for ( part=mg->part, part0=mg0->part; part; part=part->next, part0=part0->next ) {
				q = part->view.quaternion();
				q = qt * q;
				view_pred = View(q);
				symview = symmetry_get_all_views(sym, part0->view);
				min_angle = M_PI;
				for ( i=0; i<nsymviews; i++ ) {
					angle = view_pred.angle(symview[i]);
					if ( min_angle > angle ) min_angle = angle;
				}
				if ( min_angle > angle_cutoff ) {
					part0->sel = 0;
				} else {
					avg_angle += min_angle;
					nsel++;
				}
				i = (int) (min_angle*180.0/M_PI);
				hist[i]++;
				if ( max < i ) max = i;
				delete[] symview;
			}
		}
		for ( mg=field->mg->next; mg; mg=mg->next )
			for ( part=mg->part, part0=mg0->part; part; part=part->next, part0=part0->next )
				part->sel = part0->sel;		// Set all micrographs to have the same selection
	}
	if ( nsel ) avg_angle /= nsel;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Angular deviation histogram:\nAngle\tCount" << endl;
		for ( i=0; i<max; i++ ) cout << i << tab << hist[i] << endl;
		cout << endl;
	}
	
	if ( verbose & VERB_LABEL ) {
		cout << "Average angle deviation:        " << avg_angle*180/M_PI << endl;
		cout << "Particle sets selected:         " << nsel << endl << endl;
	}
	
	return avg_angle;
}

long		part_convert_selection(Bproject* project, int* sel)
{
	long				n(0), nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					part->sel = sel[n];
					if ( part->sel > 0 ) nsel++;
					n++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				part->sel = sel[n];
				if ( part->sel > 0 ) nsel++;
				n++;
			}
		}
	}

	return nsel;
}

/**
@brief 	Selects groups of particles based on the FOM.
@param 	*project	parameter structure with all parameters.
@param 	ngroups		number of groups of particles.
@param 	fom_index	index of FOM value to select on.
@param 	defocus_fit	flag to compensate for defocus.
@return long		number of particles selected in the first group.

	Particles are ranked according to the figure-of-merit and equal 
	numbers are distributed to the requested number of particle groups.
	The particle group selections are then written into a selection 
	array in the STAR data base, with the group with the best FOM's first.
	Only particles already selected are subject to the test.

**/
long		part_select_FOM_groups(Bproject* project, int ngroups, int fom_index, int defocus_fit)
{
	long				i, n(0), npart(0), nsel(0);
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting " << ngroups << " particle groups based on FOM ranking:" << endl;
	
	int_float*			rank = project_fom_order(project, npart, fom_index, defocus_fit);
	
	if ( npart < 1 ) {
		cerr << "Error: No particles found!" << endl;
		return -1;
	}

	if ( project->select < 1 ) nsel = project_count_mg_part_selected(project);
	else nsel = project_count_rec_part_selected(project);
	
	int*				sel = new int[npart];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_FOM_groups: Writing the selection array:" << endl;
	// Write the ranking groups into the selection array
	for ( i=0; i<nsel; i++ ) {
		n = ngroups - (i*ngroups)/nsel;
		sel[rank[i].i] = n;
	}
	for ( ; i<npart; i++ ) sel[rank[i].i] = 0;
	
	// Write the selection array back into the parameter hierarchy
	nsel = part_convert_selection(project, sel);

	delete[] sel;
	delete[] rank;
	
	return nsel;
}

/**
@author David Belnap & Bernard Heymann
@brief 	Selects a percentage of particles based on the FOM.
@param 	*project	parameter structure with all parameters.
@param 	percentage	percentage of particles (0 - 100).
@param 	fom_index	which FOM value to select on.
@param 	defocus_fit	flag to compensate for defocus.
@return long		number of particles selected.

	Particles are ranked according to the figure-of-merit and the
	desired percentage selected.

**/
long		part_select_percentage(Bproject* project, double percentage, int fom_index, int defocus_fit)
{
	if ( percentage < 0 ) percentage = 0;
	if ( percentage > 100 ) percentage = 100;
	
	long			npart(0), nsel(0), n;
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting " << percentage << "% particles based on FOM[" << fom_index << "] ranking:" << endl;
	
	int_float*		rank = project_fom_order(project, npart, fom_index, defocus_fit);
	
	if ( npart < 1 ) {
		cerr << "Error: No particles found!" << endl;
		return -1;
	}

	if ( project->select < 1 ) nsel = project_count_mg_part_selected(project);
	else nsel = project_count_rec_part_selected(project);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_percentage: npart=" << npart << " nsel=" << nsel << endl;
	
	// Get the FOM cutoff
	long			num_to_select = (long) (nsel*percentage/100.0);
	if ( num_to_select < 0 ) num_to_select = 0;
	if ( num_to_select >= npart ) num_to_select = npart - 1;
	
	int*				sel = new int[npart];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_percentage: Writing the selection array:" << endl;
	// Write the top number into the selection array
	for ( n=0; n<num_to_select; n++ ) sel[rank[n].i] = 1;
	for ( ; n<npart; n++ ) sel[rank[n].i] = 0;

	// Write the selection array back into the parameter hierarchy
	nsel = part_convert_selection(project, sel);
	
	delete[] sel;
	delete[] rank;
	
	return nsel;
}

/**
@brief 	Selects the best particles based on the FOM.
@param 	*project	parameter structure with all parameters.
@param 	number		number of particles.
@param 	fom_index	which FOM value to select on.
@param 	defocus_fit	flag to compensate for defocus.
@return long		number of particles selected.

	Particles are ranked according to the figure-of-merit and the
	desired percentage selected.

**/
long		part_select_best(Bproject* project, long number, int fom_index, int defocus_fit)
{
	if ( number < 1 ) return 0;
	
	long			npart(0), nsel(0), n;
	
	if ( project->select < 1 ) nsel = project_count_mg_part_selected(project);
	else nsel = project_count_rec_part_selected(project);
	
	if ( number > nsel ) number = nsel;
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting " << number << " particles based on FOM[" << fom_index << "] ranking:" << endl;
	
	int_float*		rank = project_fom_order(project, npart, fom_index, defocus_fit);
	
	if ( npart < 1 ) {
		cerr << "Error: No particles found!" << endl;
		return -1;
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_best: npart=" << npart << " nsel=" << nsel << endl;
	
	int*				sel = new int[npart];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_best: Writing the selection array:" << endl;

	// Write the top number into the selection array
	for ( n=0; n<number; n++ ) sel[rank[n].i] = 1;
	for ( ; n<npart; n++ ) sel[rank[n].i] = 0;

	// Write the selection array back into the parameter hierarchy
	nsel = part_convert_selection(project, sel);
	
	delete[] sel;
	delete[] rank;
	
	return nsel;
}

/**
@brief 	Selects particles based on FOM average and standard deviation.
@param 	*project		parameter structure with all parameters.
@param 	factor 			factor to multiply standard deviation with.
@param 	fom_index			which FOM value to select on.
@return long					number of particles selected.

	Particles are selected using the average and standard deviation,
	with the level set as a function of difference from the average:
		selected > average + factor*std_dev
	Note the multiplying factor can be negative.
	Only particles already selected are subject to the test.

**/
long		part_select_FOM_avg_std(Bproject* project, double factor, int fom_index)
{
	long				n(0);
	double				avg(0), std(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & (VERB_LABEL | VERB_PROCESS) )
		cout << "Selecting particles based on the FOM average and standard deviation" << endl;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					avg += part->fom[fom_index];
					std += part->fom[fom_index]*part->fom[fom_index];
					n++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				avg += part->fom[fom_index];
				std += part->fom[fom_index]*part->fom[fom_index];
				n++;
			}
		}
	}

	if ( n ) avg /= n;
	if ( n > 1 ) std = sqrt(std/n - avg*avg);
	else std = 0;
	
	if ( verbose & VERB_PROCESS )
		cout << "Average and standard deviation: " << avg << " " << std << endl << endl;

	long			nsel = part_deselect(project, fom_index, avg + factor*std);
	
	return nsel;
}

/**
@brief 	Selects a specific number of particles randomly.
@param 	*project		parameter structure with all parameters.
@param 	number 			number to select.
@return long					number of particles selected.

	An array is set up for all selected particles. The given number of random
	elements is selected in the array and transfered as particle selections.

**/
long		part_select_random(Bproject* project, long number)
{
	random_seed();
	
	long				i, n, nsel;

	if ( project->select < 1 ) {
		nsel = project_count_mg_part_selected(project);
	} else {
		nsel = project_count_rec_part_selected(project);
	}
	
	if ( number >= nsel ) return nsel;

	long				d = (long) (get_rand_max()/nsel);
	if ( d < 1 ) d = 1;
	
	int					v = (number < nsel/2)? 1: 0;
	int*				narr = new int[nsel];
	for ( i=0; i<nsel; i++ ) narr[i] = 1-v;
	
	for ( i=0; i<number; i++ ) {
		while ( ( n = random()/d ) >= nsel ) ;
		while ( n < nsel && narr[n] == v ) n++;
		if ( n >= nsel ) n = nsel - 1;
		while ( n >= 0 && narr[n] == v ) n--;
		if ( n < 0 ) n = 0;
		if ( narr[n] == v )
			cerr << "Warning: A random selection has been missed!" << endl;
		else
			narr[n] = v;
	}
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_PROCESS )
		cout << "Selecting " << number << " particles randomly" << endl;
	
	n = nsel = 0;
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
					if ( narr[n] < 1 ) part->sel = 0;
					if ( part->sel ) nsel++;
					n++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) if ( part->sel ) {
				if ( narr[n] < 1 ) part->sel = 0;
				if ( part->sel ) nsel++;
				n++;
			}
		}
	}
	
	delete[] narr;
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects particles randomly based on a fraction of the total.
@param 	*project		parameter structure with all parameters.
@param 	fraction 		fraction of total to select.
@return long					number of particles selected.

	A random number between 0 and 1 is generated for each particle and if 
	it is smaller than the given fraction, the particle is selected.

**/
long		part_select_random_fraction(Bproject* project, double fraction)
{
	random_seed();
	
	long				nsel(0);
	long				cut = (long) (fraction*get_rand_max());
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_PROCESS )
		cout << "Selecting " << fraction*100 << "% of the particles randomly" << endl;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( random() > cut ) part->sel = 0;
					if ( part->sel ) nsel++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( random() > cut ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects groups randomly to give at least a number of particles.
@param 	*project		parameter structure with all parameters.
@param 	number 			number to select.
@return long				number of particles selected.

	A random number between 0 and 1 is generated for each particle and if 
	it is smaller than the given fraction, the particle is selected.

**/
long		part_select_random_group(Bproject* project, long number)
{
	random_seed();
	
	long				i, max(0), nsel(0);
	double				irm;
//	int*				sel;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_PROCESS )
		cout << "Picking groups randomly to select at least " << number << " particles" << endl;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					if ( max < part->group ) max = part->group;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				if ( max < part->group ) max = part->group;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG part_select_random_group: max=" << max << endl;
	
	max++;
	irm = max*1.0L/(get_rand_max() + 1.0L);
//	sel = new int[max];
//	for ( i=0; i<max; i++ ) sel[i] = 0;
	
	vector<long>		sel(max,0);
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next )
					sel[part->group]--;
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next )
				sel[part->group]--;
	}

	while ( nsel < number ) {
		i = (long) (irm * random());
		if ( sel[i] < 0 ) {
			sel[i] = -sel[i];
			nsel += sel[i];
		}
	}

	nsel = 0;
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next )
				for ( part = mg->part; part; part = part->next ) {
					if ( sel[part->group] <= 0 ) part->sel = 0;
					if ( part->sel ) nsel++;
				}
	} else {
		for ( rec = project->rec; rec; rec = rec->next )
			for ( part = rec->part; part; part = part->next ) {
				if ( sel[part->group] <= 0 ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
	}
	
//	delete[] sel;
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Generates a selection number from the group number for each filament.
@param 	*project		parameter structure with all parameters.
@param 	nmaps			number of maps to select for.
@return long				number of particles selected.
**/
long		part_select_random_filaments(Bproject* project, int nmaps)
{
	random_seed();

	long				r = get_rand_max()/nmaps;

	long				nsel(0);
	int					prev_group(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	vector<long>		sel(nmaps+1,0);
	
	if ( verbose )
		cout << "Selecting filaments for " << nmaps << " maps:" << endl;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					if ( part->group != prev_group ) {
						nsel = random()/r + 1;
						if ( nsel > nmaps ) nsel = nmaps;
						prev_group = part->group;
					}
					part->sel = nsel;
				}
				sel[part->sel]++;
			}
			prev_group = 0;
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel > 0 ) {
						if ( part->group != prev_group ) {
							nsel = random()/r + 1;
							if ( nsel > nmaps ) nsel = nmaps;
							prev_group = part->group;
						}
						part->sel = nsel;
					}
					sel[part->sel]++;
				}
				prev_group = 0;
			}
		}
	}
	
	if ( verbose ) {
		cout << "Map\tCount" << endl;
		for ( long i=0; i<sel.size(); ++i )
			cout << i << tab << sel[i] << endl;
		cout << endl;
	}
	
	return nsel;
}

/**
@brief 	Selects a given number of particles randomly with replacement.
@param 	*project		parameter structure with all parameters.
@param 	number				number to select.
@return long					number of particles selected.

	A random number between 1 and the number of selected particles is 
	generated the given number of times. The selection value for the
	selected particle is incremented each time. A particle may therefore
	be selected more than once.

**/
long		part_select_bootstrap(Bproject* project, int number)
{
	random_seed();

	long			i, rn, n = project_count_mg_part_selected(project), nsel(0);
	
	if ( verbose & VERB_PROCESS )
		cout << "Selecting " << number << " from " << n << " particles randomly with replacement for bootstrapping" << endl;
	
	double			rnf = n*1.0L/get_rand_max();
	
	int*			sel = new int[n];
	for ( i=0; i<n; i++ ) sel[i] = 0;
	
	for ( i=0; i<number; i++ ) {
		rn = n;
		while ( rn >= n ) rn = (long) (rnf*random());
		sel[rn]++;
	}
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	for ( i=nsel=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel ) {
					part->sel = sel[i++];
					if ( part->sel ) nsel++;
				}
			}
		}
	}
	
	delete[] sel;

	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects a given number of particles within each view.
@param 	*project		parameter structure with all parameters.
@param 	*sym			symmetry structure.
@param 	theta_step		angular step size from primary symmetry axis (radians).
@param 	phi_step		angular step size around primary symmetry axis (radians).
@param 	number			number within view to select.
@return long			number of particles selected.

	A random number between 1 and the number of selected particles is 
	generated the given number of times. The selection value for the
	selected particle is incremented each time. A particle may therefore
	be selected more than once.
	This selection is meant to be used with the bootstrap reconstruction.

**/
long		part_select_random_within_view(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, int number)
{
	random_seed();

	View*			viewlist = asymmetric_unit_views(sym, theta_step, phi_step, 1);
	
	long	nview = count_list((char *) viewlist);
	long			i, j, k, n0;
	double			d, dmin;
	double			rnf = 1.0L/get_rand_max();
	
	int*			np = new int[nview];
	int*			nr = new int[number];
	int*			nrn = new int[number];
	for ( i=0; i<nview; i++ ) np[i] = 0;
	for ( i=0; i<number; i++ ) nr[i] = nrn[i] = 0;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	View*			view;
	
	// Assign views and count the number of particles per view
	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					for ( j=0, dmin=10, view = viewlist; view; view = view->next, j++ ) {
						d = view->distance(part->view);
						if ( d < dmin ) {
							dmin = d;
							part->sel = -(j+1);	// Negative to not be confused with positive selection later
						}
					}
					np[-part->sel-1]++;
				}
			}
		}
	}
	
	kill_list((char *) viewlist, sizeof(View));

	for ( j=n0=0; j<nview; j++ ) if ( np[j] < 1 ) n0++;
	
	if ( verbose ) {
		cout << "Selecting a maximum of " << number << " particles from each view:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Number of views:                " << nview << endl;
		cout << "View coverage:                  " << 100*(1 - n0*1.0/nview) << " %" << endl;
	}
		
	// For each view, select the desired number of particles randomly
	for ( j=0; j<nview; j++ ) if ( np[j] ) {
		for ( k=0; k<number; k++ ) {
			nr[k] = (int) (random()*rnf*np[j]);
			nrn[k] = 1;
			for ( i=0; i<k; i++ ) if ( nr[i] == nr[k] ) nrn[i]++;
		}
		if ( verbose & VERB_DEBUG ) {
			cout << "DEBUG part_select_random_within_view: " << j+1 << "(" << np[j] << ") ";
			for ( k=0; k<number; k++ ) cout << nr[k] << "(" << nrn[k] << ") ";
			cout << endl;
		}
		for ( i=0, field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel == -(j+1) ) {
						for ( k=0; k<number && nr[k] != i; k++ ) ;
						if ( k >= number ) part->sel = 0;
						else part->sel = nrn[k];	// Selection is set to be used with bootstrapping reconstruction
						i++;
					}
				}
			}
		}
	}
	
	delete[] np;
	delete[] nr;
	delete[] nrn;

	long			nsel = project_count_mg_part_selected(project);

	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects particles based on a smoothed maximum surface.
@param 	*project		parameter structure with all parameters.
@param 	*sym			symmetry structure.
@param 	theta_step		angular step size from primary symmetry axis (radians).
@param 	phi_step		angular step size around primary symmetry axis (radians).
@param 	threshfrac		fraction of maximum threshold.
@param 	sigma			smoothing parameter: gaussian sigma.
@param 	fom_index		which FOM value to select on.
@return long			number of particles selected.

	The best particles within each view is selected.

**/
long		part_select_maxsmooth(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, double threshfrac, double sigma, int fom_index)
{
	View*			viewlist = asymmetric_unit_views(sym, theta_step, phi_step, 1);
	
	long			nview = count_list((char *) viewlist);
	long			i, j, n0;
	double			isig2(-0.5/(sigma*sigma)), d, dmin, f;
	
	long*			np = new long[nview];
	double*			fommax = new double[nview];
	for ( i=0; i<nview; i++ ) fommax[i] = np[i] = 0;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	View			*view, *view2, symview;
	
	// Assign views and count the number of particles per view
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					for ( j=0, dmin=10, view = viewlist; view; view = view->next, j++ ) {
						d = view->angle(part->view);
						if ( d < dmin ) {
							dmin = d;
							part->sel = j+1;
						}
					}
					np[part->sel-1]++;
					if ( fommax[part->sel-1] < part->fom[fom_index] )
						fommax[part->sel-1] = part->fom[fom_index];
				}
			}
		}
	}

	for ( i=n0=0; i<nview; i++ ) if ( np[i] < 1 ) n0++;
	
	if ( verbose ) {
		cout << "Selecting based on smoothed maxima surface:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Number of views:                " << nview << endl;
		cout << "View coverage:                  " << 100*(1 - n0*1.0/nview) << " %" << endl;
		cout << "Threshold fraction:             " << threshfrac << endl;
		cout << "Smoothing sigma:                " << sigma*180.0/M_PI << " degrees" << endl << endl;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Maxima:" << endl << "#\tViewX\tViewY\tCount\tFOMmax" << endl;
		for ( i=0, view = viewlist; view; view = view->next, i++ )
			cout << i+1 << tab << view->x() << tab << view->y() << tab << np[i] << tab << fommax[i] << endl;
	}
	
	// Smooth maxima surface
	for ( i=0, view = viewlist; view; view = view->next, i++ ) {
		for ( j=0, view2 = viewlist; view2; view2 = view2->next, j++ ) 
				if ( i != j ) {
			symview = find_closest_symmetric_view(sym, *view, *view2);
			d = view->angle(symview);
			f = exp(d*d*isig2);
			if ( fommax[i] > fommax[j] ) {
				f *= fommax[i];
				if ( fommax[j] < f ) fommax[j] = f;
			} else {
				f *= fommax[j];
				if ( fommax[i] < f ) fommax[i] = f;
			}
		}
	}

	if ( verbose & VERB_FULL ) {
		cout << "Smoothed maxima:" << endl << "#\tViewX\tViewY\tCount\tFOMmax" << endl;
		for ( i=0, view = viewlist; view; view = view->next, i++ )
			cout << i+1 << tab << view->x() << tab << view->y() << tab << np[i] << tab << fommax[i] << endl;
	}
	
	kill_list((char *) viewlist, sizeof(View));

	// Set the thresholds
	for ( i=0; i<nview; i++ ) fommax[i] *= threshfrac;

	// Select particles based on smoothed maxima surface
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next )
				if ( part->sel > 0 )
					if ( part->fom[fom_index] < fommax[part->sel-1] )
						part->sel = 0;

	long			nsel = project_count_mg_part_selected(project);

	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}


/**
@brief 	Selects the best number of particles within each view.
@param 	*project		parameter structure with all parameters.
@param 	*sym			symmetry structure.
@param 	theta_step		angular step size from primary symmetry axis (radians).
@param 	phi_step			angular step size around primary symmetry axis (radians).
@param 	number				number within view to select.
@param 	fom_index			which FOM value to select on.
@return long					number of particles selected.

	The best particles within each view is selected.

**/
long		part_select_best_within_view(Bproject* project, Bsymmetry& sym, 
				double theta_step, double phi_step, int number, int fom_index)
{
	View*			viewlist = asymmetric_unit_views(sym, theta_step, phi_step, 1);
	
	long	nview = count_list((char *) viewlist);
	long			j, k, n, n0, imin;
	double			d, dmin, fommin;
	
	int*			np = new int[nview];
	Bparticle**		bestpart = new Bparticle*[number];
	for ( j=0; j<nview; j++ ) np[j] = 0;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	View*			view;
	
	// Assign views and count the number of particles per view
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( part->sel > 0 ) {
					for ( j=0, dmin=10, view = viewlist; view; view = view->next, j++ ) {
						d = view->distance(part->view);
						if ( d < dmin ) {
							dmin = d;
							part->sel = j+1;
						}
					}
					np[part->sel-1]++;
				}
			}
		}
	}
	
	kill_list((char *) viewlist, sizeof(View));

	for ( j=n0=0; j<nview; j++ ) if ( np[j] < 1 ) n0++;
	
	if ( verbose ) {
		cout << "Selecting a maximum of " << number << " best particles from each view:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Number of views:                " << nview << endl;
		cout << "View coverage:                  " << 100*(1 - n0*1.0/nview) << " %" << endl;
	}
		
	// For each view, select the best desired number of particles
	for ( j=0; j<nview; j++ ) if ( np[j] ) {
		for ( k=n=0; k<number; k++ ) bestpart[k] = NULL;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->sel == j+1 ) {
						if ( n < number ) bestpart[n++] = part;	// Fill up the array
						else {
							for ( k=0, imin=0, fommin=1e30; k<number; k++ ) {	// Find the minimum
								if ( fommin > bestpart[k]->fom[fom_index] ) {
									fommin = bestpart[k]->fom[fom_index];
									imin = k;
								}
							}
							if ( part->fom[fom_index] > fommin ) {	// Replace the minimum
								bestpart[imin]->sel = 0;
								bestpart[imin] = part;
							} else part->sel = 0;
						}
					}
				}
			}
		}
	}
	
	delete[] np;
	delete[] bestpart;

	long			nsel = project_count_mg_part_selected(project);

	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Transfer selection numbers to the group identifiers.
@param 	*project	parameter structure with all parameters.
@return long		number of particles selected.
**/
long		part_select_to_group(Bproject* project)
{
	if ( !project ) return 0;
	
	long				nsel(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_LABEL )
		cout << "Transferring selections to groups:" << endl;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				part->group = part->sel;
				if ( part->sel ) ++nsel;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					part->group = part->sel;
					if ( part->sel ) ++nsel;
				}
			}
		}
	}
	
	if ( verbose & VERB_LABEL )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects particles with a given group identifier.
@param 	*project	parameter structure with all parameters.
@param 	group		group identifier.
@return long		number of particles selected.
**/
long		part_select_group(Bproject* project, int group)
{
	if ( !project ) return 0;
	
	long				nsel(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting particles of group " << group << endl;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) if ( part->sel ) {
				if ( part->group == group ) part->sel = ++nsel;
				else part->sel = 0;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
					if ( part->group == group ) part->sel = ++nsel;
					else part->sel = 0;
				}
			}
		}
	}
	
	if ( verbose & VERB_LABEL )
		cout << "Particles selected:             " << nsel << endl << endl;

	return nsel;
}

/**
@brief 	Selects sets of particles, each set with the same size.
@param 	*project	parameter structure with all parameters.
@param 	size		number of particles in each set.
@param 	flag		flag to not count across mg or rec boundaries.
@return long			number of sets selected.

	Sets up sets of particles, each set identified as a number in the
	selection array.

**/
long		part_select_sets(Bproject* project, int size, int flag)
{
	if ( !project || size <= 0 ) return 0;
	
	long				nsel(0), count(0), np, number(0), ntot(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & VERB_LABEL )
		cout << "Generating sets of " << size << " particles each" << endl;
	
	if ( project->select ) {
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) mg->select = 0;
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( flag && count ) {
				count = 0;
				number++;
			}
			for ( np=0, part = rec->part; part; part = part->next, ntot++ ) {
				if ( part->sel ) {
					if ( count >= size ) {
						number++;
						count = 0;
					}
					if ( number < 1 ) number = 1;
					part->sel = number;
					count++;
					nsel++;
					np++;
				}
			}
			if ( np < 1 ) rec->select = 0;
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) rec->select = 0;
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( flag && count ) {
					count = 0;
					number++;
				}
				for ( np=0, part = mg->part; part; part = part->next, ntot++ ) {
					if ( part->sel ) {
						if ( count >= size ) {
							number++;
							count = 0;
						}
						if ( number < 1 ) number = 1;
						part->sel = number;
						count++;
						nsel++;
						np++;
					}
				}
				if ( np < 1 ) mg->select = 0;
			}
		}
	}
	
	if ( verbose ) {
		if ( ntot < 1 ) ntot = 1;
		cout << "Number of sets generated:       " << number << endl;
		cout << "Number of particles selected:   " << nsel << " (" << nsel*100.0/ntot << " %)" << endl << endl;
	}
	
	return number;
}

/**
@brief 	Select particles from frames in a series.
@param 	*project	parameter structure with all parameters.
@param 	frame_start	first frame (starts at 1).
@param 	frame_end	last frame.
@return long		number of particles selected.

	Particles in frames from a series in the same field-of-view are 
	selected. Only particles that are already selected are considered.

**/
long		part_select_frames(Bproject* project, int frame_start, int frame_end)
{
	int				n, sel;
	long			nsel(0);
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	if ( verbose & (VERB_LABEL | VERB_PROCESS) )
		cout << "Selecting particle series." << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( n=1, mg = field->mg; mg; mg = mg->next, n++ ) {
			sel = 0;
			if ( n >= frame_start && n <= frame_end ) sel = 1;
			for ( part = mg->part; part; part = part->next ) {
				if ( !sel ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
		}
	}

	return nsel;
}

/**
@brief 	List of particle images for reciprocal space reconstruction.
@param 	*project		image processing parameter structure.
@param 	num_select		selection number from the selection column.
@param 	bootstrap		flag to indicate a bootstrap reconstruction.
@return	Bimage*			3D reconstructed map.

	An image is used in the reconstruction if its selection flag has been set.
	If the selection number is less than zero, all particles with selection flags
	greater than zero are used. If the selection number is zero or above, all
	particles with the selection flag set to the same number are used.
	A bootstrap reconstruction uses the particle selection to weigh each
	selected particle.

**/
Bparticle* 	project_selected_partlist(Bproject* project,
				int num_select, int bootstrap)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Breconstruction*	rec = project->rec;
	Bparticle*		part = mg->part;
	Bparticle*		partlist = NULL;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG project_selected_partlist: num_select = " << num_select << endl;

	long 			psel;

	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					psel = 0;
					if ( part->sel ) {
						if ( bootstrap ) {
							psel = 1;
						} else if ( num_select < 0 ) {
							psel = 1;
						} else if ( part->sel == num_select ) {
							psel = 1;
						}
					}
					if ( psel )	// Use only selected images
						particle_copy(&partlist, part);
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				psel = 0;
				if ( part->sel ) {
					if ( bootstrap ) {
						psel = 1;
					} else if ( num_select < 0 ) {
						psel = 1;
					} else if ( part->sel == num_select ) {
						psel = 1;
					}
				}
				if ( psel )	// Use only selected images
					particle_copy(&partlist, part);
			}
		}
	}
	
/*	if ( num_select == 4) {
		cout << num_select;
		for ( part = partlist; part; part = part->next ) cout << tab << part->id;
		cout << endl;
	}
*/
	return partlist;
}

long		part_fix_fil_direction(Bparticle* partlist)
{
	long				nchanged(0);
	int					i, group(-1), n[2], dir(0);
	double				da, a[2];
	View				view, pview;
	Bparticle*			part;
	Bparticle*			fil_start = NULL;
	Bparticle*			fil_part = NULL;

	for ( part = partlist; part; part = part->next ) {
		view = part->view.backward();
		if ( part->group != group ) {
			fil_start = part;
			dir = 0;
			n[0] = 1;
			n[1] = 0;
			a[0] = a[1] = 0;
			group = part->group;
		} else {
			da = view.angle(pview);
			if ( da > M_PI_2 ) {
				dir = 1 - dir;
				da -= M_PI;
				if ( verbose & VERB_FULL )
					cout << group << tab << part->id << tab << da << tab << 1-2*dir << endl;
			}
			a[dir] += da*da;
			n[dir]++;
			if ( verbose & VERB_FULL )
				cout << group << tab << part->id << tab << da << tab << 1-2*dir << endl;
		}
		if ( part->next == NULL || part->next->group != group ) {
			for ( i=0; i<2; i++ ) if ( n[i] ) a[i] = sqrt(a[i]/n[i]);
			if ( verbose & VERB_PROCESS )
				cout << part->group << tab << n[0] << tab << n[1] << tab << 
					a[0]*180.0/M_PI << tab << a[1]*180.0/M_PI << endl;
			pview = fil_start->view.backward();
			if ( n[1] > n[0] ) { pview.negate(); pview[3] = -pview.angle(); }
			for ( fil_part = fil_start; fil_part && fil_part->group == group; fil_part = fil_part->next ) {
				view = fil_part->view.backward();
				da = view.angle(pview);
				if ( da > M_PI_2 ) {
					view.negate(); view[3] = -view.angle();
					fil_part->view = view.backward();
					nchanged++;
				}
				pview = view;
			}
		}
		pview = part->view.backward();
	}
	
	return nchanged;
}

long		part_fix_fil_direction_old(Bparticle* partlist, double minpct)
{
	if ( minpct < 1 ) minpct *= 100;
	
	long				n(0), up(0), dn(0), nt(0);
	int					group(-1), fin(0), dir;
	double				tol(20*M_PI/180), pup, pdn;
	Euler				euler;
	Bparticle*			part;
	Bparticle*			fil_part = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Fil#\t#part\tup\t%up\tdown\t%down\toutliers" << endl;

	for ( part = partlist; part; part = part->next ) {
		euler = Euler(part->view);
		
		if ( part->group != group ) {
			n = up = dn = fin = 0;
			pup = pdn = 0;
			group = part->group;
			fil_part = part;
		}
		
		n++;
		if ( euler.psi() <= -M_PI_2 + tol && euler.psi() >= -M_PI_2 - tol ) up++;
		if ( euler.psi() <= M_PI_2 + tol && euler.psi() >= M_PI_2 - tol ) dn++;

		if ( part->next ) {
			if ( part->next->group != group ) fin = 1;
		} else fin = 1;
		
		if ( fin ) {
			if ( n ) {
				pup = up*100.0/n;
				pdn = dn*100.0/n;
				nt += up + dn;
			}
			
			if ( verbose & VERB_PROCESS )
				cout << group << tab << n << tab << up << tab << pup << tab << dn << tab << pdn << tab << n-up-dn << endl;
				
			if ( pup >= minpct ) {
				dir = 1;
			} else if ( pdn >= minpct ) {
				dir = -1;
			} else {
				dir = 0;
			}
			
			for ( ; fil_part && fil_part->group == group; fil_part = fil_part->next ) if ( fil_part->sel > 0 ) {
				euler = Euler(fil_part->view);
				fil_part->sel = 0;
				if ( euler.psi() <= -M_PI_2 + tol && euler.psi() >= -M_PI_2 - tol && dir == 1 ) fil_part->sel = 1;
				if ( euler.psi() <= M_PI_2 + tol && euler.psi() >= M_PI_2 - tol && dir == -1 ) fil_part->sel = 1;
			}
		}
	}
	
	return nt;
}

vector<long>	part_fix_fil_direction(Bparticle* partlist, double minpct)
{
	if ( minpct < 1 ) minpct *= 100;
	
	long				n(0), up(0), dn(0), nt(0), nf(0), fup(0), fdn(0);
	long				group(-1), fin(0), dir;
	double				tol(M_PI/4.0), pup, pdn, theta(0), da, ccmax;
//	Matrix3				mat1, mat;
	Euler				euler;
	Bparticle*			part;
	Bparticle*			fil_part = NULL;
	Bparticle*			part_max = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Fil#\t#part\t%up\t%down\tdirection" << endl;
//		cout << "Fil#\t#part\tup\t%up\tdown\t%down\toutliers" << endl;

	for ( part = partlist; part; part = part->next ) {

		if ( part->group != group ) {	// start of filament
			group = part->group;
			ccmax = part->fom[0];
			// Find the best particle in the filament to set the reference view
			for ( fil_part = part_max = part; fil_part && fil_part->group == group; fil_part = fil_part->next )
				if ( ccmax < fil_part->fom[0] ) {
					ccmax = fil_part->fom[0];
					part_max = fil_part;
				}
//			mat1 = part_max->view.matrix();
			euler = Euler(part_max->view);
			theta = euler[1];
			fil_part = part;
			n = up = dn = fin = 0;
			pup = pdn = 0;
//			cout << part_max->id << tab << part->group << tab << part_max->view << tab << ccmax << endl;
		}
		
		n++;
//		mat = part->view.matrix();
		euler = Euler(part->view);
//		da = fabs(angle_set_negPI_to_PI(mat.angle(mat1)));
		da = fabs(angle_set_negPI_to_PI(theta - euler[1]));
//		cout << part->id << tab << part->group << tab << part->view << tab << da << endl;
		if ( da < tol ) up++;
		else if ( da > M_PI - tol ) dn++;

		if ( part->next ) {
			if ( part->next->group != group ) fin = 1;
		} else fin = 1;
		
		if ( fin ) {
			nf++;
			if ( n ) {
				pup = up*100.0/n;
				pdn = dn*100.0/n;
				nt += up + dn;
			}
			
			if ( pup >= minpct || pdn >= minpct ) {
				if ( pup >= pdn ) {
					dir = 1;
					fup++;
				} else {
					dir = -1;
					fdn++;
				}
			} else {
				dir = 0;
			}
			
			if ( verbose & VERB_PROCESS )
				cout << group << tab << n << tab << pup << tab << pdn << tab << dir << endl;
//				cout << group << tab << n << tab << up << tab << pup << tab << dn << tab << pdn << tab << n-up-dn << endl;
			
			for ( ; fil_part && fil_part->group == group; fil_part = fil_part->next ) if ( fil_part->sel > 0 ) {
				fil_part->sel = 0;
				if ( dir != 0 ) {
					fil_part->sel = 1;
//					mat = fil_part->view.matrix();
					euler = Euler(fil_part->view);
//					da = fabs(angle_set_negPI_to_PI(mat.angle(mat1)));
					da = fabs(angle_set_negPI_to_PI(theta - euler[1]));
					if ( ( da < M_PI_2 && dir == -1 ) || ( da > M_PI_2 && dir == 1 ) ) {
						euler[1] = angle_set_negPI_to_PI(euler[1]+M_PI);
						fil_part->view = euler.view();
					}
				}
			}
		}
	}
	
	vector<long>	nfstats(3);
	nfstats[0] = nf;
	nfstats[1] = fup;
	nfstats[2] = fdn;
	
	return nfstats;
}

vector<long>	part_fix_fil_direction(Bparticle* partlist, Bparticle* partlist2, double minpct)
{
	if ( minpct < 1 ) minpct *= 100;
	
	long				n(0), up(0), dn(0), nt(0), nf(0), fup(0), fdn(0);
	long				group(-1), fin(0), dir;
	double				tol(M_PI/4.0), pup, pdn, theta(0), da, ccmax;
	Euler				euler;
	Bparticle*			part;
	Bparticle*			part2;
	Bparticle*			fil_part = NULL;
	Bparticle*			fil_part2 = NULL;
	Bparticle*			part_max = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Fil#\t#part\t%up\t%down\tdirection" << endl;
//		cout << "Fil#\t#part\tup\t%up\tdown\t%down\toutliers" << endl;

	for ( part = partlist, part2 = partlist2; part && part2; part = part->next, part2 = part2->next ) {

		if ( part->group != group ) {	// start of filament
			group = part->group;
			ccmax = part->fom[0];
			// Find the best particle in the filament to set the reference view
			for ( fil_part = part_max = part; fil_part && fil_part->group == group; fil_part = fil_part->next )
				if ( ccmax < fil_part->fom[0] ) {
					ccmax = fil_part->fom[0];
					part_max = fil_part;
				}
			euler = Euler(part_max->view);
			theta = euler[1];
			fil_part = part;
			fil_part2 = part2;
			n = up = dn = fin = 0;
			pup = pdn = 0;
//			cout << part_max->id << tab << part->group << tab << part_max->view << tab << ccmax << endl;
		}
		
		n++;
		euler = Euler(part->view);
		da = fabs(angle_set_negPI_to_PI(theta - euler[1]));
//		cout << part->id << tab << part->group << tab << part->view << tab << da << endl;
		if ( da < tol ) up++;
		else if ( da > M_PI - tol ) dn++;

		if ( part->next ) {
			if ( part->next->group != group ) fin = 1;
		} else fin = 1;
		
		if ( fin ) {
			nf++;
			if ( n ) {
				pup = up*100.0/n;
				pdn = dn*100.0/n;
				nt += up + dn;
			}
			
			if ( pup >= minpct || pdn >= minpct ) {
				if ( pup >= pdn ) {
					dir = 1;
					fup++;
				} else {
					dir = -1;
					fdn++;
				}
			} else {
				dir = 0;
			}
			
			if ( verbose & VERB_PROCESS )
				cout << group << tab << n << tab << pup << tab << pdn << tab << dir << endl;
			
			for ( ; fil_part && fil_part2 && fil_part->group == group; fil_part = fil_part->next, fil_part2 = fil_part2->next ) {
				if ( fil_part->fom[0] < fil_part2->fom[0] ) {
					particle_update(fil_part, fil_part2);
					euler = Euler(fil_part->view);
					euler[1] = angle_set_negPI_to_PI(euler[1]+M_PI);
					fil_part->view = euler.view();
				}
				fil_part->sel = 0;
				if ( dir != 0 ) {
					fil_part->sel = 1;
					euler = Euler(fil_part->view);
					da = fabs(angle_set_negPI_to_PI(theta - euler[1]));
					if ( ( da < M_PI_2 && dir == -1 ) || ( da > M_PI_2 && dir == 1 ) ) {
						euler[1] = angle_set_negPI_to_PI(euler[1]+M_PI);
						fil_part->view = euler.view();
					}
				}
			}
		}
	}
	
	vector<long>	nfstats(3);
	nfstats[0] = nf;
	nfstats[1] = fup;
	nfstats[2] = fdn;
	
	return nfstats;
}

/**
@brief 	Checks the direction of particles associated with filaments.
@param 	*project	project structure with all parameters.
@param 	minpct		minimum percentage.
@return long			number of particles with clear direction.

	The selected particles associated with a filament is checked for their
	in-plane direction. If the percentage of particles exceed the
	minimum desired, the particles with the corresponding direction
	are selected, and the rest deselected.

**/
long 		part_filament_direction(Bproject* project, double minpct)
{
	if ( !project ) return 0;
	if ( minpct < 1 ) minpct *= 100;
	
	long				npart(0);
	vector<long>		nfs(3,0), nfsall(3,0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Checking and fixing filament directions:" << endl;
		cout << "Minimum percentage directed particles: " << minpct << endl << endl;
	}
	
	if ( project->select ) {
		npart = project_count_rec_particles(project);
		for ( rec = project->rec; rec; rec = rec->next ) {
			if ( verbose )
				cout << "Reconstruction: " << rec->id << endl;
			nfs = part_fix_fil_direction(rec->part, minpct);
			transform(nfsall.begin(), nfsall.end(), nfs.begin(), nfs.begin(), plus<long>());
		}
	} else {
		npart = project_count_mg_particles(project);
		for ( field = project->field; field; field = field->next )
			for ( mg = field->mg; mg; mg = mg->next ) {
				if ( verbose )
					cout << "Micrograph: " << mg->id << endl;
				nfs = part_fix_fil_direction(mg->part, minpct);
				transform(nfsall.begin(), nfsall.end(), nfs.begin(), nfs.begin(), plus<long>());
			}
	}
	
//	if ( verbose & VERB_LABEL )
//		cout << endl << "Particles with direction:       " << ndir << " (" << ndir*100.0/npart << "%)" << endl << endl;

	if ( verbose ) {
		cout << "Filaments:                      " << nfsall[0] << endl;
		cout << "Filaments up:                   " << nfsall[1] << endl;
		cout << "Filaments down:                 " << nfsall[2] << endl;
	}

	return npart;
}

/**
@brief 	Checks the direction of particles associated with filaments and compares with opposite filament direction.
@param 	*project	project structure with all parameters.
@param 	*project2	second project structure with opposite filamant directions.
@param 	minpct		minimum percentage.
@return long			number of particles with clear direction.

	The selected particles associated with a filament is checked for their
	in-plane direction. If the percentage of particles exceed the
	minimum desired, the particles with the corresponding direction
	are selected, and the rest deselected.

**/
long 		part_filament_direction(Bproject* project, Bproject* project2, double minpct)
{
	if ( !project ) return 0;
	if ( minpct < 1 ) minpct *= 100;
	
	long				npart(0);
	vector<long>		nfs(3,0), nfsall(3,0);
	
	Bfield*				field, *field2;
	Bmicrograph*		mg, *mg2;
	Breconstruction*	rec, *rec2;
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Checking and fixing filament directions:" << endl;
		cout << "Minimum percentage directed particles: " << minpct << endl << endl;
	}
	
	if ( project->select ) {
		npart = project_count_rec_particles(project);
		for ( rec = project->rec, rec2 = project2->rec; rec && rec2; rec = rec->next, rec2 = rec2->next ) {
			if ( verbose )
				cout << "Reconstruction: " << rec->id << endl;
			nfs = part_fix_fil_direction(rec->part, rec2->part, minpct);
			transform(nfsall.begin(), nfsall.end(), nfs.begin(), nfs.begin(), plus<long>());
		}
	} else {
		npart = project_count_mg_particles(project);
		for ( field = project->field, field2 = project2->field; field && field2; field = field->next, field2 = field2->next )
			for ( mg = field->mg, mg2 = field2->mg; mg && mg2; mg = mg->next, mg2 = mg2->next ) {
				if ( verbose )
					cout << "Micrograph: " << mg->id << endl;
				nfs = part_fix_fil_direction(mg->part, mg2->part, minpct);
				transform(nfsall.begin(), nfsall.end(), nfs.begin(), nfs.begin(), plus<long>());
			}
	}
	
//	if ( verbose & VERB_LABEL )
//		cout << endl << "Particles with direction:       " << ndir << " (" << ndir*100.0/npart << "%)" << endl << endl;

	if ( verbose ) {
		cout << "Filaments:                      " << nfsall[0] << endl;
		cout << "Filaments up:                   " << nfsall[1] << endl;
		cout << "Filaments down:                 " << nfsall[2] << endl;
	}
	
	return npart;
}

/**
@brief 	Selects particles based on a distance from a reference view.
@param 	*project	project structure with all parameters.
@param 	view		central view to search for.
@param 	angle		angular distance from the view.
@return long		number of particles selected.

	Only particles already selected are subject to the test.

**/
long 		part_view_select(Bproject* project, View view, double angle)
{	
	long				nsel(0);
	double				da;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & (VERB_LABEL|VERB_PROCESS) )
		cout << "Selecting particles based on a view" << endl;
	if ( verbose & VERB_PROCESS ) {
		cout << "View:                           " << view << endl;
		cout << "Deviation angle:                " << angle*180/M_PI << " degrees" << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					da = part->view.angle(view);
					if ( da > angle ) part->sel = 0;
					if ( part->sel ) nsel++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				da = part->view.angle(view);
				if ( da > angle ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Selects particles based on a distance from a side view.
@param 	*project	project structure with all parameters.
@param 	angle		angular distance from the side view.
@return long		number of particles selected.

	Only particles already selected are subject to the test.

**/
long 		part_side_view_select(Bproject* project, double angle)
{	
	long				nsel(0);
	double				dz(sin(angle));
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & (VERB_LABEL|VERB_PROCESS) )
		cout << "Selecting particles with side views" << endl;
	if ( verbose & VERB_PROCESS ) {
		cout << "Deviation angle:                " << angle*180/M_PI << " degrees" << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( part->view[2] > dz ) part->sel = 0;
					if ( part->sel ) nsel++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( part->view[2] > dz ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Selects particles based on Euler angle limits.
@param 	*project	project structure with all parameters.
@param 	*euler6		6-value vector of Euler angle limits.
@return long		number of particles selected.

	Only particles already selected are subject to the test.

**/
long 		part_euler_angle_select(Bproject* project, double* euler6)
{	
	long				i, nsel(0);
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Euler				euler;
	
	if ( verbose & (VERB_LABEL|VERB_PROCESS) )
		cout << "Selecting particles based on Euler angle limits" << endl;
	if ( verbose & VERB_PROCESS ) {
		cout << "Psi limits:                     " << euler6[0]*180/M_PI << " - " << euler6[1]*180/M_PI << " degrees" << endl;
		cout << "Theta limits:                   " << euler6[2]*180/M_PI << " - " << euler6[3]*180/M_PI << " degrees" << endl;
		cout << "Phi limits:                     " << euler6[4]*180/M_PI << " - " << euler6[5]*180/M_PI << " degrees" << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) if ( part->sel ) {
					euler = Euler(part->view);
					for ( i=0; i<3; i++ ) {
						if ( euler[i] < euler6[2*i] ) euler[i] += TWOPI;
						if ( euler[i] > euler6[2*i+1] ) euler[i] -= TWOPI;
						if ( euler[i] < euler6[2*i] || euler[i] > euler6[2*i+1] ) part->sel = 0;
					}
					if ( part->sel ) nsel++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) if ( part->sel ) {
				euler = Euler(part->view);
				for ( i=0; i<3; i++ ) {
					if ( euler[i] < euler6[2*i] ) euler[i] += TWOPI;
					if ( euler[i] > euler6[2*i+1] ) euler[i] -= TWOPI;
					if ( euler[i] < euler6[2*i] || euler[i] > euler6[2*i+1] ) part->sel = 0;
				}
				if ( part->sel ) nsel++;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Selects particles based on the distance from the nominal origin.
@param 	*project	project structure with all parameters.
@param 	origin		nominal origin.
@param 	distance	maximum distance to accept.
@return long		number of particles selected.

	Only particles already selected are subject to the test.

**/
long 		part_origin_select(Bproject* project, Vector3<double> origin, double distance)
{	
	long				nsel(0);
	double				dd;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	if ( verbose & (VERB_LABEL|VERB_PROCESS) )
		cout << "Selecting particles based on the origin" << endl;
	if ( verbose & VERB_PROCESS ) {
		cout << "Origin:                         " << origin << endl;
		cout << "Maximum distance:               " << distance << " pixels" << endl;
	}
	
	if ( project->select < 1 ) {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					dd = part->ori.distance(origin);
					if ( dd > distance ) part->sel = 0;
					if ( part->sel ) nsel++;
				}
			}
		}
	} else {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				dd = part->ori.distance(origin);
				if ( dd > distance ) part->sel = 0;
				if ( part->sel ) nsel++;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}


/**
@brief 	Sets all orientation parameters for a particle to the first in the set.
@param 	*project	parameter structure with all parameters.
@return int			0.

	This function sets the views of all the particles in the field-of-view to 
	the views for the first micrograph.

**/
long		part_set_first_view_in_series(Bproject* project)
{
	int				f;
	Bfield*			field;
	Bmicrograph*	mg, *mg0;
	Bparticle*		part, *part0;
	
	if ( verbose & (VERB_LABEL | VERB_PROCESS) )
		cout << "Setting all to the first orientations in the first micrograph in the field-of-view" << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg0=field->mg, mg=field->mg->next; mg; mg=mg->next ) {
			for ( part=mg->part, part0=mg0->part; part; part=part->next, part0=part0->next ) {
				for ( f=0; f<NFOM; f++ ) part->fom[f] = part0->fom[f];
				part->sel = part0->sel;
				part->view = part0->view;
			}
		}
	}

	return 0;
}

/**
@brief 	Sets all orientation parameters for a particle to the best in the set.
@param 	*project	parameter structure with all parameters.
@param 	fom_index	index of FOM value to select on.
@return int			0.

	This function uses the FOM to select the best orientation parameters
	for a particle in a set (typically a focal series) and sets all
	orientation parameters to the same values.

**/
long		part_set_best_view_in_series(Bproject* project, int fom_index)
{
	int				f;
	Bfield*			field;
	Bmicrograph*	mg, *mg0;
	Bparticle*		part, *part0;
	
	if ( verbose & (VERB_LABEL | VERB_PROCESS) )
		cout << "Setting all to the best orientations in the field-of-view" << endl;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg0=field->mg, mg=field->mg->next; mg; mg=mg->next ) {
			for ( part=mg->part, part0=mg0->part; part; part=part->next, part0=part0->next ) {
				if ( part->fom[fom_index] > part0->fom[fom_index] ) {
					for ( f=0; f<NFOM; f++ ) part0->fom[f] = part->fom[f];
					part0->sel = part->sel;
					part0->view = part->view;
				}
			}
		}
		for ( mg0=field->mg, mg=field->mg->next; mg; mg=mg->next ) {
			for ( part=mg->part, part0=mg0->part; part; part=part->next, part0=part0->next ) {
				for ( f=0; f<NFOM; f++ ) part->fom[f] = part0->fom[f];
				part->sel = part0->sel;
				part->view = part0->view;
			}
		}
	}

	return 0;
}

/**
@brief 	Select series where all are already selected.
@param 	*project	parameter structure with all parameters.
@param 	size		number of particles in each set.
@param 	flag		flag to not count across mg or rec boundaries.
@return long		number of particles selected.

	This function selects those series where all the particles
	are already selected, and deselects the rest. Sets are selected
	with each series part of the same set.

**/
long		part_select_series(Bproject* project, int size, int flag)
{
	if ( size < 1 ) size = 1000000000;
	
	int				n;
	long			count(0), nser(0), nsel(0);
	Bfield*			field;
	Bmicrograph*	mg, *mg0;
	Bparticle*		part, *part0;
	
	if ( verbose & (VERB_LABEL | VERB_PROCESS) )
		cout << "Selecting particle series." << endl;
	
	for ( field = project->field; field; field = field->next ) {
		if ( flag ) nser = 0;
		for ( n=0, mg0=field->mg, mg = field->mg; mg; mg = mg->next, n++ ) {
			for ( part0=mg0->part, part = mg->part; part; part = part->next, part0=part0->next ) {
				if ( part->sel ) {
					if ( mg0 == mg ) part0->sel = 1;
					else part0->sel++;
				}
			}
		}
		for ( part0=field->mg->part; part0; part0=part0->next ) {
			if ( part0->sel == n ) {
				if ( nser%size == 0 ) count++;
				part0->sel = count;
				nser++;
			} else {
				part0->sel = 0;
			}
		}
		for ( mg0=field->mg, mg = field->mg; mg; mg = mg->next ) {
			for ( part0=mg0->part, part = mg->part; part; part = part->next, part0=part0->next ) {
				if ( part0->sel ) {
					part->sel = part0->sel;
					nsel++;
				} else {
					part->sel = 0;
				}
			}
		}
	}

	return nsel;
}

/**
@brief 	Calculates the locations of particles in a series from a seed.
@param 	*project	micrograph project.
@param 	flags		bit 1=invert z, bit 2=use original particle view.
@return long		number of particles created.

	The seed particle locations is either in the zero-degree tilt 2D micrograph,
	or a 3D reconstruction. The latter also gives the z-coordinate, resulting
	in a better definition of the particle location.
	The project selection flag indicates micrograph or reconstruction.

**/
long		part_series_from_seed(Bproject* project, int flags)
{
	int					zflip(flags & 1), use_view(flags & 2);
	long				i, n, npart(0);
	Matrix3				mat(1);
	Bfield*				field;
	Bmicrograph*		mg;
	Bmicrograph*		mg_ref = NULL;
	Breconstruction*	rec;
	Bparticle*			part_list = NULL;
	Bparticle			*part_ref, *part;
	Vector3<long>		part_size;
	Vector3<double>		origin, scale(1,1,1), loc, part_ori;

	if ( !project->rec ) project->select = 0;
	
	// Find the reference particle list
	if ( project->select ) {
		if ( verbose )
			cout << "Calculating micrograph particle locations from 3D locations" << endl;
		for ( rec = project->rec; rec && !part_list; rec = rec->next )
			if ( rec->part ) break;
		if ( rec ) {
			part_list = rec->part;
			origin = rec->origin;
			scale /= rec->scale;
			part_size = scale*rec->box_size;
		}
	} else {
		if ( verbose )
			cout << "Calculating micrograph particle locations from zero-tilt micrograph" << endl;
		for ( field = project->field; field && !part_list; field = field->next ) {
			for ( mg = field->mg; mg && !part_list; mg = mg->next ) if ( mg->part ) {
				mg_ref = mg;
				part_list = mg->part;
				origin = mg->origin;
				part_size = mg->box_size;
			}
		}
	}
	
	if ( !part_list ) {
		cerr << "Error: No particle records found!" << endl;
		return -1;
	}

	part_size[2] = 1;
	part_ori = part_size/2;

	// Set the selection number
	for ( n=0, part_ref=part_list; part_ref; part_ref=part_ref->next )
		part_ref->sel = ++n;

	double			ca(1), sa(0), tt(0);
	
	// Generate new particle records
	if ( verbose )
		cout << "Number of particles:            " << n << endl;
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg != mg_ref ) {
			if ( mg->part ) particle_kill(mg->part);
			mg->part = part = NULL;
			mg->box_size = part_size;
			ca = cos(mg->tilt_axis);
			sa = sin(mg->tilt_axis);
			tt = tan(mg->tilt_angle);
			for ( part_ref=part_list; part_ref; part_ref=part_ref->next ) {
				part = particle_add(&part, part_ref->id);
				if ( !mg->part ) mg->part = part;
				loc = (part_ref->loc - origin)*scale;
				if ( zflip ) loc[2] = -loc[2];
				part->loc = mg_location_from_3D_model(loc, mg->matrix, mg->origin, mg->scale);
				part->ori = part_ori;
				if ( use_view )
					mat = part_ref->view.matrix();
				part->view = View(mat * mg->matrix);
				for ( i=0; i<NFOM; i++ ) part->fom[i] = part_ref->fom[i];
				part->sel = part_ref->sel;
				if ( mg->ctf ) {
					loc = (part->loc - mg->origin)*mg->pixel_size;
					part->def = mg->ctf->defocus_average() +
						(loc[0]*sa - loc[1]*ca)*tt;
					part->dev = mg->ctf->defocus_deviation();
					part->ast = mg->ctf->astigmatism_angle();
				}
				npart++;
			}
		}
	}

	return npart;
}

/**
@brief 	Fits the particle FOMs as a function of the defocus.
@param 	*project	parameter structure with all parameters.
@param 	fom_index	index of FOM value to test for.
@param 	&intercept	fit intercept.
@param 	&slope		fit slope.
@return double		correlation index of fit.

	The FOM is fit as a linear function of the defocus average.

**/
double		part_fom_defocus_fit(Bproject* project, int fom_index, double& intercept, double& slope)
{
	slope = 0;
	intercept = 1;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	long			i, npart, nmg = project_count_micrographs(project);
	double			sum, CI(0);
	double*			x = new double[nmg];
	double*			y = new double[nmg];
	for ( i=0; i<nmg; i++ ) x[i] = y[i] = 0;

	if ( verbose & VERB_FULL )
		cout << "Defocus\tFOMavg" << endl;
	for ( i=0, field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			for ( npart=0, sum=0, part = mg->part; part; part = part->next )
					if ( part->sel ) {
						sum += part->fom[fom_index];
						npart++;
					}
			if ( npart ) {
				x[i] = mg->ctf->defocus_average();
				y[i] = sum/npart;
				if ( verbose & VERB_FULL )
					cout << x[i] << tab << y[i] << endl;
				i++;
			}
		}
	if ( verbose & VERB_FULL )
		cout << endl;

	if ( i ) CI = linear_least_squares(0, i-1, x, y, intercept, slope);
	
	delete[] x;
	delete[] y;

	if ( verbose & VERB_PROCESS ) {
		cout << "Linear fit:                     FOM = " << slope << " * defocus + " << intercept << endl;
		cout << "Correlation index:              " << CI << endl << endl;
	}
	
	return CI;
}

/**
@brief 	Deselects particles below a given FOM cutoff, adjusted for the defocus.
@param 	*project	parameter structure with all parameters.
@param 	fom_index	index of FOM value to test for.
@param 	cutoff		threshold for deselection.
@return double		correlation index of fit.

	The FOM is fit as a linear function of the defocus average.
	Particles are deselected based on the adjusted FOM cutoff:
		adj_cut = slope * defocus + cut

**/
double		part_fom_defocus_fit_deselect(Bproject* project, int fom_index, double cutoff)
{
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	
	double			slope, intercept, this_cutoff;
	double			CI = part_fom_defocus_fit(project, fom_index, intercept, slope);
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
//			this_cutoff = slope*mg->ctf->defocus_average() + cutoff;
			this_cutoff = slope*(mg->ctf->defocus_average() - 1e4) + cutoff;
			for ( part = mg->part; part; part = part->next )
				if ( part->fom[fom_index] < this_cutoff) part->sel = 0;
		}
	
	return CI;
}

/**
@brief 	Deselects all particles except those from micrographs closest to focus.
@param 	*project	parameter structure with all parameters.
@return long		number of particles selected.
**/
long		part_select_closest_to_focus(Bproject* project)
{
	long			nsel(0);
	Bfield*			field;
	Bmicrograph*	mg, *mg_ctf;
	Bparticle*		part;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = mg_ctf = field->mg; mg; mg = mg->next ) {
			if ( mg_ctf->ctf->defocus_average() > mg->ctf->defocus_average() )
				mg_ctf = mg;
		}
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg != mg_ctf )
				for ( part = mg->part; part; part = part->next )
					part->sel = 0;
			else
				for ( part = mg->part; part; part = part->next )
					if ( part->sel ) nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Deselects all particles except those from micrographs furthest from focus.
@param 	*project	parameter structure with all parameters.
@return long		number of particles selected.
**/
long		part_select_furthest_from_focus(Bproject* project)
{
	long			nsel(0);
	Bfield*			field;
	Bmicrograph*	mg, *mg_ctf;
	Bparticle*		part;
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = mg_ctf = field->mg; mg; mg = mg->next ) {
			if ( mg_ctf->ctf->defocus_average() < mg->ctf->defocus_average() )
				mg_ctf = mg;
		}
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( mg != mg_ctf )
				for ( part = mg->part; part; part = part->next )
					part->sel = 0;
			else
				for ( part = mg->part; part; part = part->next )
					if ( part->sel ) nsel++;
		}
	}
	
	return nsel;
}

/**
@brief 	Deletes deselected particles and renumber the remaining.
@param 	*project	parameter structure with all parameters.
@return long			number of particles remaining.
**/
long		part_delete_deselected(Bproject* project)
{
	long				id, nsel(0);
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part, *part_nu, *part_list;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( id=0, part_list = NULL, part = rec->part; part; part = part->next ) {
				if ( part->sel ) {
					part_nu = particle_copy(&part_list, part);
					part_nu->id = ++id;
				}
			}
			particle_kill(rec->part);
			rec->part = part_list;
			nsel += id;
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( id=0, part_list = NULL, part = mg->part; part; part = part->next ) {
					if ( part->sel ) {
						part_nu = particle_copy(&part_list, part);
						part_nu->id = ++id;
					}
				}
				particle_kill(mg->part);
				mg->part = part_list;
				nsel += id;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Particles selected:             " << nsel << endl << endl;
	
	return nsel;
}

/**
@brief 	Deletes deselected particles and renumber the remaining.
@param 	**partlist	particle parameter structure.
@return long			number of particles remaining.

	The old list is deallocated.
	
**/
long		part_delete_deselected(Bparticle** partlist)
{
	long				id(0);
	Bparticle*			part, *part_nu, *nulist;
	
	for ( nulist = NULL, part = *partlist; part; part = part->next ) {
		if ( part->sel ) {
			part_nu = particle_copy(&nulist, part);
			part_nu->id = ++id;
		}
	}
	
	particle_kill(*partlist);
	
	*partlist = nulist;
	
	return id;
}

/**
@brief 	Resets defocus values for particles thata are too far off.
@param 	*project	parameter structure with all parameters.
@param 	max_dev		maximum defocus difference allowed.
@return long			number of defocus values reset.
**/
long		part_fix_defocus(Bproject* project, double max_dev)
{
	long				nd(0);
	double				dmin, dmax;
	Bfield*				field;
	Bmicrograph*		mg;
	Bparticle*			part;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Resetting particle defocus:" << endl;
		cout << "Maximum difference allowed:     " << max_dev << " A" << endl;
	}
		
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg->ctf ) {
			dmin = mg->ctf->defocus_average() - max_dev;
			dmax = mg->ctf->defocus_average() + max_dev;
			for ( part = mg->part; part; part = part->next ) {
				if ( part->def < dmin || part->def > dmax ) {
					part->def = mg->ctf->defocus_average();
					nd++;
				}
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of defocus values reset: " << nd << endl << endl;
	
	return nd;
}

/**
@brief 	Sorts the Particles by a selected parameter.
@param 	*project 			project parameter structure.
@param 	tag 				parameter tag.
@return vector<pair<Bparticle*,double>>	array of particle links and values.
**/
vector<pair<Bparticle*,double>>	project_part_sort(Bproject* project, Bstring tag)
{
	Bfield*				field;
	Bmicrograph*		mg;
	Bparticle*			part;
	
	vector<std::pair<Bparticle*,double>> 	part_value;
	pair<Bparticle*,double>					pv;
	
	for ( field = project->field; field; field = field->next ) {
		 for ( mg = field->mg; mg; mg = mg->next ) {
		 	for ( part = mg->part; part; part = part->next ) {
				pv.first = part;
				if ( tag == PARTICLE_FOM ) pv.second = part->fom[0];
				if ( tag == PARTICLE_FOM_CV ) pv.second = part->fom[1];
				part_value.push_back(pv);
			}
		}
	}
				
	sort(part_value.begin(), part_value.end(), [=](pair<Bparticle*,double> a, std::pair<Bparticle*,double> b) {
			return a.second < b.second;
		}
	);

	return part_value;
}
