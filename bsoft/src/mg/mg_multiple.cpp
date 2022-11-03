/**
@file	mg_multiple.cpp
@brief	Selection of single particle parameters from multiple files for classification
@author Bernard Heymann
@date	Created: 20010319
@date	Modified: 20220831
**/

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_select.h"
#include "mg_particle_select.h"
#include "rwmg.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

long		project_check_same_number_particles(Bproject* project_list)
{
	long				i, npart(0), npart2(0), err(0);
	Bproject*			project = NULL;

	if ( project_list->select < 1 ) {
		npart = project_count_mg_particles(project_list);
		if ( verbose )
			cout << "Number of particles:            " << npart << endl;
		for ( i=2, project=project_list->next; project; project=project->next, ++i ) {
			npart2 = project_count_mg_particles(project);
			if ( npart2 != npart ) {
				err++;
				cerr << "Error: Incorrect number of particles in project " << i << ": " << npart2 << endl;
			}
		}
	} else {
		npart = project_count_rec_particles(project_list);
		if ( verbose )
			cout << "Number of particles:            " << npart << endl;
		for ( i=2, project=project_list->next; project; project=project->next, ++i ) {
			npart2 = project_count_rec_particles(project);
			if ( npart2 != npart ) {
				err++;
				cerr << "Error: Incorrect number of particles in project " << i << ": " << npart2 << endl;
			}
		}
	}

	if ( err ) bexit(-1);
	
	return npart;
}

/**
@brief 	Merges multiple projects.
@param 	*file_list		linked list of parameter file names.
@param 	fom_index		index of FOM to select on.
@param 	flags			flags: 1=rec, 2=template, 4=reset.
@return Bproject*			merged project.

	The selected particles are merged.
	Requirement: The project component ID's must correspond.
	If bit 2 of the flag is set (2), it signifies a template and the FOM's of the template
	are set to zero.

**/
Bproject*	project_multi_merge(Bstring* file_list, int fom_index, int flags)
{
	int					use_rec(flags & 1);
	int					temp_flag(flags & 2);
	int					reset(flags & 4);
	
	long				nmg(0), nrec(0);
	Bstring*			filename;
	Bproject*			project = NULL;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bproject*			project2 = NULL;
	Bfield*				field2 = NULL;
	Bmicrograph*		mg2 = NULL;
	Breconstruction*	rec2 = NULL;
	
	if ( verbose ) {
		cout << "Merging parameter files:" << endl;
		cout << "FOM index:                  " << fom_index << endl;
		if ( temp_flag ) cout << "First project is a template" << endl;
		cout << endl;
	}

	for ( filename = file_list; filename; filename = filename->next ) {
		project2 = read_project(*filename);
		if ( use_rec ) project2->select = 1;
		if ( reset ) part_reset_selection(project2, 3);
		if ( !project ) {
			project = project2;
			if ( temp_flag ) {
				if ( verbose )
					cout << "Resetting the template selections" << endl << endl;
				part_set_selection(project, -1);
				part_set_FOM(project, 0, 0);
				part_set_FOM(project, 1, 0);
			}
		} else {
			for ( field2 = project2->field; field2;  ) {	// Field to be merged
				// Find the same field
				for ( field = project->field; field && field->id != field2->id; field = field->next ) ;
				if ( !field ) {
					if ( project->field ) {
						for ( field = project->field; field->next; field = field->next ) ;
						field->next = field2;
					} else {
						project->field = field2;
					}
					field = field2;
					field2 = field2->next;
					field->next = NULL;
				} else {
					for ( mg2 = field2->mg; mg2;  ) {
						// Find the same micrograph
						for ( mg = field->mg; mg && mg->id != mg2->id; mg = mg->next ) ;
						if ( !mg ) {
							if ( field->mg ) {
								for ( mg = field->mg; mg->next; mg = mg->next ) ;
								mg->next = mg2;
							} else {
								field->mg = mg2;
							}
							mg = mg2;
							mg2 = mg2->next;
							mg->next = NULL;
						} else {
							mg2->select = 1;
							micrograph_update(mg, mg2, fom_index, 63);
							mg = mg2;
							mg2 = mg2->next;
							micrograph_kill(mg);
							nmg++;
						}
					}
					field2->mg = NULL;
					field = field2;
					field2 = field2->next;
					field_kill(field);
				}
			}
			project2->field = NULL;
			for ( rec2 = project2->rec; rec2;  ) {
				// Find the same reconstruction
				for ( rec = project->rec; rec && rec2->id != rec->id; rec = rec->next ) ;
				if ( rec ) {
					if ( rec2->fom >= rec->fom ) {
						rec2->select = 1;
						reconstruction_update(rec, rec2, fom_index, 63);
					}
					rec = rec2;
					rec2 = rec2->next;
					reconstruction_kill(rec);
					nrec++;
				} else {
					if ( project->rec ) {
						for ( rec = project->rec; rec->next; rec = rec->next ) ;
						rec->next = rec2;
					} else {
						project->rec = rec2;
					}
					rec = rec2;
					rec2 = rec2->next;
					rec->next = NULL;
				}
			}
			project2->rec = NULL;
			project_kill(project2);
		}
	}

	long			npart(0), nsel(0);
	double			percentage(0);
	npart += project_count_mg_particles(project);
	nsel += project_count_mg_part_selected(project);
	npart += project_count_rec_particles(project);
	nsel += project_count_rec_part_selected(project);
	
	if ( npart ) percentage = nsel*100.0/npart;

	if ( verbose ) {
		cout << "Micrographs updated:            " << nmg << endl;
		cout << "Reconstructions updated:        " << nrec << endl;
		cout << "Number of particles selected:   " << nsel << " (" << percentage << " %)" << endl << endl;
	}
	
	return project;
}

/**
@brief 	Merges multiple projects with an existing one.
@param 	*project			project to add/merge to.
@param 	*file_list			linked list of parameter file names.
@param 	fom_index			index of FOM to select on.
@return long				selected particles.

	The selected particles are merged with existing ones.
	Requirement: The project component ID's must correspond.

**/
long		project_multi_add(Bproject* project, Bstring* file_list, int fom_index)
{
	long				nmg(0), nrec(0);
	Bstring*			filename;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bproject*			project2 = NULL;
	Bfield*				field2 = NULL;
	Bmicrograph*		mg2 = NULL;
	Breconstruction*	rec2 = NULL;
	
	if ( verbose ) {
		cout << "Merging parameter files:" << endl;
		cout << "FOM index:                  " << fom_index << endl;
		cout << endl;
	}

	for ( filename = file_list; filename; filename = filename->next ) {
		project2 = read_project(*filename);
		if ( !project ) {
			project = project2;
		} else {
			for ( field2 = project2->field; field2;  ) {	// Field to be merged
				// Find the same field
				for ( field = project->field; field && field->id != field2->id; field = field->next ) ;
				if ( !field ) {
					if ( project->field ) {
						for ( field = project->field; field->next; field = field->next ) ;
						field->next = field2;
					} else {
						project->field = field2;
					}
					field = field2;
					field2 = field2->next;
					field->next = NULL;
				} else {
					for ( mg2 = field2->mg; mg2;  ) {
						// Find the same micrograph
						for ( mg = field->mg; mg && mg->id != mg2->id; mg = mg->next ) ;
						if ( !mg ) {
							if ( field->mg ) {
								for ( mg = field->mg; mg->next; mg = mg->next ) ;
								mg->next = mg2;
							} else {
								field->mg = mg2;
							}
							mg = mg2;
							mg2 = mg2->next;
							mg->next = NULL;
						} else {
							mg2->select = 1;
							micrograph_update(mg, mg2, fom_index, 63);
							mg = mg2;
							mg2 = mg2->next;
							micrograph_kill(mg);
							nmg++;
						}
					}
					field2->mg = NULL;
					field = field2;
					field2 = field2->next;
					field_kill(field);
				}
			}
			project2->field = NULL;
			for ( rec2 = project2->rec; rec2;  ) {
				// Find the same reconstruction
				for ( rec = project->rec; rec && rec2->id != rec->id; rec = rec->next ) ;
				if ( rec ) {
					if ( rec2->fom >= rec->fom ) {
						rec2->select = 1;
						reconstruction_update(rec, rec2, fom_index, 63);
					}
					rec = rec2;
					rec2 = rec2->next;
					reconstruction_kill(rec);
					nrec++;
				} else {
					if ( project->rec ) {
						for ( rec = project->rec; rec->next; rec = rec->next ) ;
						rec->next = rec2;
					} else {
						project->rec = rec2;
					}
					rec = rec2;
					rec2 = rec2->next;
					rec->next = NULL;
				}
			}
			project2->rec = NULL;
			project_kill(project2);
		}
	}

	long			npart(0), nsel(0);
	double			percentage(0);
	npart += project_count_mg_particles(project);
	nsel += project_count_mg_part_selected(project);
	npart += project_count_rec_particles(project);
	nsel += project_count_rec_part_selected(project);
	
	if ( npart ) percentage = nsel*100.0/npart;

	if ( verbose ) {
		cout << "Micrographs updated:            " << nmg << endl;
		cout << "Reconstructions updated:        " << nrec << endl;
		cout << "Number of particles selected:   " << nsel << " (" << percentage << " %)" << endl << endl;
	}
	
	return nsel;
}

/**
@brief 	Adds particles from a list of projects to the first one.
@param 	*file_list			linked list of parameter file names.
@return Bproject*				project with renumbered particles.

	The particles from other files are added to existing ones in the first file.
	Where project component ID's correspond the particles are added,
	otherwise new components are added.
	Particles are renumbered and the references to old particle files removed.

**/
Bproject*	project_multi_add_particles(Bstring* file_list)
{
	long				nmg(0), nrec(0);
	Bproject*			project = NULL;
	Bstring*			filename;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	Bproject*			project2 = NULL;
	Bfield*				field2 = NULL;
	Bmicrograph*		mg2 = NULL;
	Breconstruction*	rec2 = NULL;
	
	if ( verbose ) {
		cout << "Adding particles from other parameter files" << endl;
		cout << endl;
	}

	for ( filename = file_list; filename; filename = filename->next ) {
		project2 = read_project(*filename);
		if ( !project ) {
			project = project2;
		} else {
			for ( field2 = project2->field; field2;  ) {	// Field to be merged
				// Find the same field
				for ( field = project->field; field && field->id != field2->id; field = field->next ) ;
				if ( !field ) {
					if ( project->field ) {
						for ( field = project->field; field->next; field = field->next ) ;
						field->next = field2;
					} else {
						project->field = field2;
					}
					field = field2;
					field2 = field2->next;
					field->next = NULL;
				} else {
					for ( mg2 = field2->mg; mg2;  ) {
						// Find the same micrograph
						for ( mg = field->mg; mg && mg->id != mg2->id; mg = mg->next ) ;
						if ( !mg ) {
							mg2->fpart = 0;
							if ( field->mg ) {
								for ( mg = field->mg; mg->next; mg = mg->next ) ;
								mg->next = mg2;
							} else {
								field->mg = mg2;
							}
							mg = mg2;
							mg2 = mg2->next;
							mg->next = NULL;
						} else {
							mg->fpart = 0;
							if ( mg->part ) {
								for ( part = mg->part; part->next; part = part->next ) ;
								part->next = mg2->part;
							} else {
								mg->part = mg2->part;
							}
							mg2->part = NULL;
							mg = mg2;
							mg2 = mg2->next;
							micrograph_kill(mg);
							nmg++;
						}
					}
					field2->mg = NULL;
					field = field2;
					field2 = field2->next;
					field_kill(field);
				}
			}
			project2->field = NULL;
			for ( rec2 = project2->rec; rec2;  ) {
				// Find the same reconstruction
				for ( rec = project->rec; rec && rec2->id != rec->id; rec = rec->next ) ;
				if ( rec ) {
					rec->fpart = 0;
					if ( rec->part ) {
						for ( part = rec->part; part->next; part = part->next ) ;
						part->next = rec2->part;
					} else {
						rec->part = rec2->part;
					}
					rec2->part = NULL;
					rec = rec2;
					rec2 = rec2->next;
					reconstruction_kill(rec);
					nrec++;
				} else {
					rec2->fpart = 0;
					if ( project->rec ) {
						for ( rec = project->rec; rec->next; rec = rec->next ) ;
						rec->next = rec2;
					} else {
						project->rec = rec2;
					}
					rec = rec2;
					rec2 = rec2->next;
					rec->next = NULL;
				}
			}
			project2->rec = NULL;
			project_kill(project2);
		}
	}

	project_renumber_particles(project);

	long			npart(0), nsel(0);
	double			percentage(0);
	npart += project_count_mg_particles(project);
	nsel += project_count_mg_part_selected(project);
	npart += project_count_rec_particles(project);
	nsel += project_count_rec_part_selected(project);
	
	if ( npart ) percentage = nsel*100.0/npart;

	if ( verbose ) {
		cout << "Micrographs updated:            " << nmg << endl;
		cout << "Reconstructions updated:        " << nrec << endl;
		cout << "Number of particles selected:   " << nsel << " (" << percentage << " %)" << endl << endl;
	}
	
	return project;
}

/**
@brief 	Adjusts FOM to the avergae of the first project.
@param 	*project_list	linked list of project structures with all parameters.
@param 	fom_index		index of FOM to select on.
@return long				number of particles selected.

	Requirement: The project structures must be of identical form.

**/
long 		project_multi_adjust_FOM(Bproject* project_list, int fom_index)
{
	long				np = count_list((char *)project_list);
	long				npart = project_check_same_number_particles(project_list);
	
	Bproject*			project;
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	
	long				i;
	vector<double>		fom(np,0);
	
	if ( verbose ) {
		cout << "Adjusting the FOM of index " << fom_index << endl;
		cout << "Project\tAdjustment" << endl;
	}
	
	if ( project_list->select < 1 ) {
		for ( i=0, project = project_list; project; project = project->next, ++i )
			for ( field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					for ( part = mg->part; part; part = part->next )
						fom[i] += part->fom[fom_index];
		for ( i=1; i<np; ++i ) {
			fom[i] = fom[0]/fom[i];
			if ( verbose )
				cout << i+1 << tab << fom[i] << endl;
		}
		fom[0] = 1;
		for ( i=0, project = project_list; project; project = project->next, ++i )
			for ( field = project->field; field; field = field->next )
				for ( mg = field->mg; mg; mg = mg->next )
					for ( part = mg->part; part; part = part->next )
						part->fom[fom_index] *= fom[i];
	} else {
		for ( i=0, project = project_list; project; project = project->next, ++i )
			for ( rec = project->rec; rec; rec = rec->next )
				for ( part = rec->part; part; part = part->next )
					fom[i] += part->fom[fom_index];
		for ( i=1; i<np; ++i ) {
			fom[i] = fom[0]/fom[i];
			if ( verbose )
				cout << i+1 << tab << fom[i] << endl;
		}
		fom[0] = 1;
		for ( i=0, project = project_list; project; project = project->next, ++i )
			for ( rec = project->rec; rec; rec = rec->next )
				for ( part = rec->part; part; part = part->next )
					part->fom[fom_index] *= fom[i];
	}
	
	return npart;
}

/**
@brief 	Selects the best particle orientations from multiple parameter files.
@param 	*project_list	linked list of project structures with all parameters.
@param 	fom_cut			FOM cutoff to eliminate bad particles.
@param 	fom_index		index of FOM to select on.
@param 	fom_def_flag	flag to adjust the cutoff for defocus.
@return long				number of particles selected.

	The orientation parameters for each particle is selected from
	the file with the highest FOM. The selection flag gets the number 
	of the file from which the particle was selected, with the first
	file getting the number 1. If the best FOM is under the cutoff,
	the selection flag is set to zero. The selected orientation
	parameters are written into the first micrograph parameter
	structure.
	Requirement: The project structures must be of identical form.

**/
long 		project_multi_select_best_FOM(Bproject* project_list, double fom_cut, int fom_index, int fom_def_flag)
{
	long			np = count_list((char *)project_list);
	
	if ( np < 1 ) {
		cout << "Selecting best FOM: Too few parameter files provided!" << endl;
		return -1;
	}
	
	if ( verbose ) {
		cout << "Selecting from " << np << " data sets:" << endl;
		cout << "FOM cutoff:                     " << fom_cut << endl;
		cout << "FOM index:                      " << fom_index << endl;
	}

	long				i, f, nsel(0), npart(0);
	double				intercept, intercept_avg(0), slope, slope_avg(0), cutoff = fom_cut;
	
	Bproject*			project = project_list;
	Bfield**			field = NULL;
	Bmicrograph**		mg = NULL;
	Breconstruction**	rec = NULL;
	Bparticle**			part = NULL;
	
	if ( project->select < 1 ) {
		npart = project_check_same_number_particles(project_list);
		field = new Bfield*[np];
		mg = new Bmicrograph*[np];
		part = new Bparticle*[np];
		for ( i=0, project=project_list; i<np; i++, project=project->next ) {
			field[i] = project->field;
			if ( fom_def_flag ) {
				part_fom_defocus_fit(project, fom_index, intercept, slope);
				intercept_avg += intercept;
				slope_avg += slope;
			}
		}
		if ( fom_def_flag ) {
			intercept_avg /= np;
			slope_avg /= np;
			if ( verbose )
				cout << "FOM vs defocus slope:           " << slope_avg << endl;
		}
		for ( ; field[0]; ) {
			if ( verbose & VERB_FULL )
				cout << "Field: " << field[0]->id << endl;
			for ( i=0; i<np; i++ ) mg[i] = field[i]->mg;
			for ( ; mg[0]; ) {
				if ( fom_def_flag && mg[0]->ctf )
					cutoff = slope_avg*mg[0]->ctf->defocus_average() + fom_cut;
				for ( i=0; i<np; i++ ) part[i] = mg[i]->part;
				for ( ; part[0]; ) {
					if ( part[0]->sel ) {
						for ( i=1; i<np; i++ ) {
							if ( part[0]->fom[fom_index] < part[i]->fom[fom_index] ) {
								part[0]->sel = i + 1;
								part[0]->view = part[i]->view;
								part[0]->ori = part[i]->ori;
								part[0]->mag = part[i]->mag;
								for ( f=0; f<NFOM; f++ )
									part[0]->fom[f] = part[i]->fom[f];
							}
						}
					}
					if ( part[0]->fom[fom_index] < cutoff ) part[0]->sel = 0;
					if ( part[0]->sel > 0 ) nsel++;
					for ( i=0; i<np; i++ ) part[i] = part[i]->next;
				}
				for ( i=0; i<np; i++ ) mg[i] = mg[i]->next;
			}
			for ( i=0; i<np; i++ ) field[i] = field[i]->next;
		}
		delete[] field;
		delete[] mg;
		delete[] part;
	} else {
		npart = project_check_same_number_particles(project_list);
		rec = new Breconstruction*[np];
		part = new Bparticle*[np];
		for ( i=0, project=project_list; i<np; i++, project=project->next ) rec[i] = project->rec;
		for ( ; rec[0]; ) {
			if ( verbose & VERB_FULL )
				cout << "Reconstruction: " << rec[0]->id << endl;
			for ( i=0; i<np; i++ ) part[i] = rec[i]->part;
			for ( ; part[0]; ) {
				if ( part[0]->sel ) {
					for ( i=1; i<np; i++ ) {
						if ( part[0]->fom[fom_index] < part[i]->fom[fom_index] ) {
							part[0]->sel = i + 1;
							part[0]->view = part[i]->view;
							part[0]->ori = part[i]->ori;
							part[0]->mag = part[i]->mag;
							for ( f=0; f<NFOM; f++ )
								part[0]->fom[f] = part[i]->fom[f];
						}
					}
				}
				if ( part[0]->fom[fom_index] < fom_cut ) part[0]->sel = 0;
				if ( part[0]->sel > 0 ) nsel++;
				for ( i=0; i<np; i++ ) part[i] = part[i]->next;
			}
			for ( i=0; i<np; i++ ) rec[i] = rec[i]->next;
		}
		delete[] rec;
		delete[] part;
	}
	
	if ( verbose )
		cout << "Number of particles selected:   " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
	
	
	return nsel;
}

/**
@brief 	Counts particle assignment distributions using selection arrays.
@param 	*project_list	linked list of project structures with all parameters.
@return long					number of particles selected.

	Requirement: The project structures must be of identical form.

**/
long 		project_multi_selection_stats(Bproject* project_list)
{
	long			np = count_list((char *)project_list);
	
	if ( np < 1 ) {
		cerr << "Error: Selection statistics: Too few parameter files provided!" << endl;
		return -1;
	}
	
	long				i, j, m, n, f, nsel(0), npart(0);
	long				comb = 1<<np;
	int*				sel = new int[comb];
	for ( i=0; i<comb; i++ ) sel[i] = 0;
	
	Bproject*			project = project_list;
	Bfield**			field = NULL;
	Bmicrograph**		mg = NULL;
	Breconstruction**	rec = NULL;
	Bparticle**			part = NULL;
	
	if ( project->select < 1 ) {
		npart = project_check_same_number_particles(project_list);
		field = new Bfield*[np];
		mg = new Bmicrograph*[np];
		part = new Bparticle*[np];
		for ( i=0, project=project_list; i<np; i++, project=project->next ) field[i] = project->field;
		while ( field[0] ) {
			for ( i=0; i<np; i++ ) mg[i] = field[i]->mg;
			while ( mg[0] ) {
				for ( i=0; i<np; i++ ) part[i] = mg[i]->part;
				while ( part[0] ) {
					m = n = 0;
					for ( i=0; i<np; i++ ) {
						if ( part[i]->sel ) {
							m += 1<<i;
							if ( n == 0 ) n = i + 1;
							else n = -1;
						}
					}
					sel[m] += 1;
					part[0]->sel = 0;
					if ( n == 1 ) part[0]->sel = 1;
					if ( n > 1 ) {
						part[0]->sel = n;
						i = n - 1;
						part[0]->sel = i + 1;
						part[0]->view = part[i]->view;
						part[0]->ori = part[i]->ori;
						part[0]->mag = part[i]->mag;
						for ( f=0; f<NFOM; f++ )
							part[0]->fom[f] = part[i]->fom[f];
					}
					if ( n > 0 ) nsel++;
					for ( i=0; i<np; i++ ) part[i] = part[i]->next;
				}
				for ( i=0; i<np; i++ ) mg[i] = mg[i]->next;
			}
			for ( i=0; i<np; i++ ) field[i] = field[i]->next;
		}
		delete[] field;
		delete[] mg;
		delete[] part;
	} else {
		npart = project_check_same_number_particles(project_list);
		rec = new Breconstruction*[np];
		part = new Bparticle*[np];
		for ( i=0, project=project_list; i<np; i++, project=project->next ) rec[i] = project->rec;
		while ( rec[0] ) {
			for ( i=0; i<np; i++ ) part[i] = rec[i]->part;
			while ( part[0] ) {
				m = n = 0;
				for ( i=0; i<np; i++ ) {
					if ( part[i]->sel ) {
						m += 1<<i;
						if ( n == 0 ) n = i + 1;
						else n = -1;
					}
				}
				sel[m] += 1;
				part[0]->sel = 0;
				if ( n == 1 ) part[0]->sel = 1;
				if ( n > 1 ) {
					part[0]->sel = n;
					i = n - 1;
					part[0]->sel = i + 1;
					part[0]->view = part[i]->view;
					part[0]->ori = part[i]->ori;
					part[0]->mag = part[i]->mag;
					for ( f=0; f<NFOM; f++ )
						part[0]->fom[f] = part[i]->fom[f];
				}
				if ( n > 0 ) nsel++;
				for ( i=0; i<np; i++ ) part[i] = part[i]->next;
			}
			for ( i=0; i<np; i++ ) rec[i] = rec[i]->next;
		}
		delete[] rec;
		delete[] part;
	}
	
//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "Selection distribution:" << endl;
		for ( i=0; i<comb; i++ ) {
			for ( j=0; j<np; j++ ) {
				m = (i>>j) & 1;
				cout << m;
			}
			cout << tab << sel[i] << endl; 
		}
		cout << "Number of particles selected with no overlap: " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
	}
	
	delete[] sel;
	
	return nsel;
}

/**
@brief 	Selects the best particle orientations from multiple parameter files.
@param 	*project_list	linked list of project structures with all parameters.
@param 	*sym			symmetry.
@param 	origin_dev		cutoff to accept origins (pixels).
@param 	view_dev			cutoff to accept views (radians).
@param 	angle_dev		cutoff to accept rotation angles (radians).
@param 	mag_dev			cutoff to accept magnifications (fraction).
@return long					number of particles selected.

	The standard deviations of particle parameters are calculated as 
	follows:
		origin_std = sqrt(var(origin_x) + var(origin_y))
		view_std = sqrt(var(view_x) + var(view_y) + var(view_z))
		angle_std = sqrt(var(angle))
		size_std = sqrt(var(magnification))
	Requirement: The origin, view and rotation angle (or Euler angle), 
		and magnification arrays in the micrographs must be defined.

**/
long 		project_multi_select_low_variance(Bproject* project_list, Bsymmetry& sym,
				double origin_dev, double view_dev, double angle_dev, double mag_dev)
{
	long			np = count_list((char *)project_list);
	
	if ( np < 2 ) {
		cout << "Selecting lowest variance: Too few parameter files provided!" << endl;
		return -1;
	}
	
	if ( verbose )
		cout << "Selecting based on the lowest variance for symmetry " << sym.label() << endl;
	
	int				i, j, h, f, nsel(0), psel, flag3D(0);
	double			var, ox_sum, oy_sum, oz_sum, ox_sqsum, oy_sqsum, oz_sqsum, origin_std, view_std, angle_std;
	View			vsum, vsqsum;
	double			size_sum, size_sqsum, size_std, defsum, defsqsum, def_std;
	double			def_avg_dev(0), origin_avg_dev(0), view_avg_dev(0), angle_avg_dev(0), size_avg_dev(0);
	
	int				hist_ang[181], hist_view[181], hist_ori[101];
	for ( i=0; i<181; i++ ) hist_ang[i] = hist_view[i] = 0;
	for ( i=0; i<101; i++ ) hist_ori[i] = 0;
	
	long   nmg = project_count_micrographs(project_list);
	long   npart = project_count_mg_particles(project_list);

	Bproject*		project;
	Bfield**		field = new Bfield*[np];
	Bmicrograph**	mg = new Bmicrograph*[np];
	Bparticle**		part = new Bparticle*[np];
//	View*			view, *v;
	View			theview;
//	double			da, da_min;
	
	for ( i=0, project=project_list; i<np; i++, project=project->next ) field[i] = project->field;
	while ( field[0] ) {
		if ( verbose & VERB_FULL )
			cout << field[0]->id << endl;
		for ( i=0; i<np; i++ ) mg[i] = field[i]->mg;
		while ( mg[0] ) {
			if ( verbose & VERB_FULL )
				cout << mg[0]->id << endl;
			defsum = defsqsum = def_std = 0;
			if ( mg[0]->ctf ) for ( i=0; i<np; i++ ) {
				defsum += mg[i]->ctf->defocus_average();
				defsqsum += mg[i]->ctf->defocus_average()*mg[i]->ctf->defocus_average();
			}
			var = defsqsum - defsum*defsum/np;
			if ( var > 0 ) def_std = sqrt(var/(np-1));
			def_avg_dev += def_std;
			for ( i=0; i<np; i++ ) part[i] = mg[i]->part;
			while ( part[0] ) {
				for ( i=psel=0; i<np; i++ ) if ( part[i]->sel ) psel++;
				if ( psel ) {
					ox_sum = oy_sum = oz_sum = ox_sqsum = oy_sqsum = oz_sqsum = 0;
					for ( i=0; i<4; i++ ) vsum[i] = vsqsum[i] = 0;
					size_sum = size_sqsum = 0;
					if ( verbose & VERB_FULL )
						cout << part[0]->id << tab << part[0]->ori << endl;
					for ( i=0; i<np; i++ ) {
						theview = find_closest_symmetric_view(sym, part[0]->view, part[i]->view);
						ox_sum += part[i]->ori[0];
						ox_sqsum += part[i]->ori[0]*part[i]->ori[0];
						oy_sum += part[i]->ori[1];
						oy_sqsum += part[i]->ori[1]*part[i]->ori[1];
						if ( flag3D ) {
							oz_sum += part[i]->ori[2];
							oz_sqsum += part[i]->ori[2]*part[i]->ori[2];
						}
						if ( theview[3] - part[0]->view[3] < -M_PI ) theview[3] += TWOPI;
						if ( theview[3] - part[0]->view[3] >  M_PI ) theview[3] -= TWOPI;
						for ( j=0; j<4; j++ ) vsum[j] += theview[j];
						for ( j=0; j<4; j++ ) vsqsum[j] += theview[j]*theview[j];
						size_sum += part[i]->mag;
						size_sqsum += part[i]->mag*part[i]->mag;
					}
					if ( verbose & VERB_FULL ) cout << endl;
					origin_std = view_std = angle_std = size_std = 0;
					var = ox_sqsum + oy_sqsum + oz_sqsum - 
						(ox_sum*ox_sum + oy_sum*oy_sum + oz_sum*oz_sum)/np;
					if ( var > 0 ) origin_std = sqrt(var/(np-1));
					if ( origin_std > origin_dev ) part[0]->sel = 0;
					h = (int) (origin_std*10 + 0.5);
					if ( h > 100 ) h = 100;
					hist_ori[h]++;
					var = vsqsum[0] + vsqsum[1] + vsqsum[2] - 
							(vsum[0]*vsum[0] + vsum[1]*vsum[1] + vsum[2]*vsum[2])/np;
					if ( var > 0 ) view_std = sqrt(var/(np-1));
					if ( view_std > view_dev ) part[0]->sel = 0;
					h = (int) (view_std*180.0/M_PI + 0.5);
					if ( h > 180 ) h = 180;
					hist_view[h]++;
					var = vsqsum.angle() - vsum.angle()*vsum.angle()/np;
					if ( var > 0 ) angle_std = sqrt(var/(np-1));
					if ( angle_std > angle_dev ) part[0]->sel = 0;
					h = (int) (angle_std*180.0/M_PI + 0.5);
					if ( h > 180 ) h = 180;
					hist_ang[h]++;
					var = size_sqsum - size_sum*size_sum/np;
					if ( var > 0 ) size_std = sqrt(var/(np-1));
					if ( size_std > mag_dev ) part[0]->sel = 0;
					for ( i=1; i<np; i++ ) part[i]->sel = part[0]->sel;
					if ( part[0]->sel > 0 ) {
						nsel++;
						origin_avg_dev += origin_std;
						view_avg_dev += view_std;
						angle_avg_dev += angle_std;
						size_avg_dev += size_std;
					}
					for ( i=1; i<np; i++ ) if ( part[0]->fom[0] < part[i]->fom[0] ) 
						for ( f=0; f<NFOM; f++ ) part[0]->fom[f] = part[i]->fom[f];
				}
				for ( i=0; i<np; i++ ) part[i] = part[i]->next;
			}
			for ( i=0; i<np; i++ ) mg[i] = mg[i]->next;
		}
		for ( i=0; i<np; i++ ) field[i] = field[i]->next;
	}
	
	def_avg_dev /= nmg;
	if ( nsel ) {
		origin_avg_dev /= nsel;
		view_avg_dev /= nsel;
		angle_avg_dev /= nsel;
		size_avg_dev /= nsel;
	}
	
	if ( verbose & VERB_RESULT ) {
		cout << "Average defocus deviation:      " << def_avg_dev << " A" << endl;
		cout << "Number of particles selected with low variance: " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
		cout << "Average origin deviation:       " << origin_avg_dev << " pixels" << endl;
		cout << "Average view deviation:         " << view_avg_dev*180.0/M_PI << " degrees" << endl;
		cout << "Average angle deviation:        " << angle_avg_dev*180.0/M_PI << " degrees" << endl;
		cout << "Average magnification deviation:" << size_avg_dev << endl << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Origin deviation histogram:" << endl;
		cout << "Dist\tCount" << endl;
		for ( h=0; h<101; h++ ) cout << h*0.1 << tab << hist_ori[h] << endl;
		cout << endl;
	
		cout << "View and in-plane rotation angle deviation histograms:" << endl;
		cout << "Angle\tView\tAngle" << endl;
		for ( h=0; h<181; h++ ) cout << h << tab << hist_view[h] << tab << hist_ang[h] << endl;
		cout << endl;
	}
	
	delete[] field;
	delete[] mg;
	delete[] part;
	
	return nsel;
}

/**
@brief 	Selects the best particle orientations from two parameter files.
@param 	*project1	first project structure with all parameters.
@param 	*project2	second project structure with all parameters.
@param 	*sym		symmetry.
@param 	origin_err	cutoff to accept origins (pixels).
@param 	view_err	cutoff to accept views (radians).
@param 	angle_err 	cutoff to accept rotation angles (radians).
@param 	mag_err		cutoff to accept magnifications (fraction).
@return long		number of particles selected.

	The error is calculated between parameters. 
	Requirement: The origin, view and rotation angle (or Euler angle), 
		and magnification arrays in the micrographs must be defined.

**/
long 		project_multi_select_low_difference(Bproject* project1, Bproject* project2, Bsymmetry& sym,
				double origin_err, double view_err, double angle_err, double mag_err)
{
	long   			i(0), nmgsel(0), npart(0), nsel(0), ndsel(0);
	double			dori, dview, dang, dsize;
	double			ddef, def_avg_err(0), pdef_avg_err(0);
	double			origin_avg_err(0), view_avg_err(0), angle_avg_err(0), size_avg_err(0);
	
	if ( verbose )
		cout << "Selecting based on the lowest difference for symmetry " << sym.label() << endl;
	
	int				h, hist_ang[361], hist_view[361], hist_ori[101];
	for ( h=0; h<361; h++ ) hist_ang[h] = hist_view[h] = 0;
	for ( h=0; h<101; h++ ) hist_ori[h] = 0;	

	Bfield			*field1, *field2;
	Bmicrograph		*mg1, *mg2;
//	Breconstruction	*rec1, *rec2;
	Bparticle		*part1, *part2;
	Bparticle**		partarr1;
	Bparticle**		partarr2;
	View			view1, view2;
	
	if ( project1->select < 1 ) {
//		nmg = project_count_micrographs(project1);
		partarr1 = project_mg_particle_array(project1, -1, npart);
		partarr2 = project_mg_particle_array(project2, -1, npart);
		for ( field1=project1->field, field2=project2->field; field1 && field2; field1=field1->next, field2=field2->next ) {
			if ( verbose & VERB_FULL )
				cout << field1->id << tab << field2->id << endl;
			for ( mg1=field1->mg, mg2=field2->mg; mg1 && mg2; mg1=mg1->next, mg2=mg2->next ) {
				if ( verbose & VERB_FULL )
					cout << tab << mg1->id << tab << mg2->id << endl;
				if ( mg1->ctf && mg2->ctf ) {
					if ( verbose & VERB_DEBUG )
						cout << "DBEUG project_multi_select_low_difference: Calculating defocus difference" << endl;
					ddef = fabs(mg1->ctf->defocus_average() - mg2->ctf->defocus_average());
					def_avg_err += ddef;
					nmgsel++;
				}
			}
		}
	} else {
		partarr1 = project_rec_particle_array(project1, -1, npart);
		partarr2 = project_rec_particle_array(project2, -1, npart);
	}
	
	if ( verbose )
		cout << "Number of particles:            " << npart << endl;
	
	for ( i=0; i<npart; i++ ) {
		part1 = partarr1[i];
		part2 = partarr2[i];
		if ( verbose & VERB_FULL )
			cout << "particle " << i << endl;
//		cout << part1 << tab << part2 << endl;
		if ( part1->sel && part2->sel ) {
			if ( part1->def > 1 && part2->def > 1 ) {
				ddef = fabs(part1->def - part2->def);
				pdef_avg_err += ddef;
				ndsel++;
			}
			dori = (part1->ori - part2->ori).length();
			view1 = part1->view;
			view2 = find_closest_symmetric_view(sym, view1, part2->view);
//			cout << part1->id << tab << view1 << tab << view2 << endl;
//			cout << view1.vector3().scalar(view2.vector3()) << endl;
			dview = view1.angle(view2);
			dang = fabs(angle_set_negPI_to_PI(view1.angle() - view2.angle()));
//			if ( dang > 0.1 )
//				cout << setprecision(3) << part1->id << tab << view1 << tab << view2 << tab << dang*180.0/M_PI << endl;
			dsize = fabs(part1->mag - part2->mag);
//			cout << part1->id << tab << dori << tab << dview << tab << dang << endl;
			h = (int) (dori*10 + 0.5);
			if ( h > 100 ) h = 100;
			hist_ori[h]++;
			h = (int) (dview*360.0/M_PI + 0.5);
			if ( h > 360 ) h = 360;
			hist_view[h]++;
			h = (int) (dang*360.0/M_PI + 0.5);
			if ( h > 360 ) h = 360;
			hist_ang[h]++;
			if ( dori > origin_err ) part1->sel = part2->sel = 0;
			if ( dview > view_err ) part1->sel = part2->sel = 0;
			if ( dang > angle_err ) part1->sel = part2->sel = 0;
			if ( dsize > mag_err ) part1->sel = part2->sel = 0;
			if ( part1->sel > 0 && part2->sel > 0 ) {
				nsel++;
				origin_avg_err += dori;
				view_avg_err += dview;
				angle_avg_err += dang;
			}
		}
	}
	
	delete[] partarr1;
	delete[] partarr2;
	
	if ( nmgsel ) def_avg_err /= nmgsel;
	if ( ndsel ) pdef_avg_err /= ndsel;
	if ( nsel ) {
		origin_avg_err /= nsel;
		view_avg_err /= nsel;
		angle_avg_err /= nsel;
		size_avg_err /= nsel;
	}
	
	if ( verbose & VERB_RESULT ) {
		if ( nmgsel )
			cout << "Average defocus error:          " << def_avg_err << " A" << endl;
		if ( ndsel )
			cout << "Average particle defocus error: " << pdef_avg_err << " A" << endl;
		cout << "Number of particles selected with low variance: " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
		cout << "Average origin error:           " << origin_avg_err << " pixels" << endl;
		cout << "Average view error:             " << view_avg_err*180.0/M_PI << " degrees" << endl;
		cout << "Average angle error:            " << angle_avg_err*180.0/M_PI << " degrees" << endl;
		cout << "Average magnification error:    " << size_avg_err << endl << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Origin error histogram:" << endl;
		cout << "Dist\tCount" << endl;
		for ( h=0; h<101; h++ ) cout << h*0.1 << tab << hist_ori[h] << endl;
		cout << endl;
	
		cout << "View and in-plane rotation angle error histograms:" << endl;
		cout << "Angle\tView\tAngle" << endl;
		for ( h=0; h<361; h++ ) cout << h*0.5 << tab << hist_view[h] << tab << hist_ang[h] << endl;
		cout << endl;
	}
	
	return nsel;
}

/**
@brief 	Selects the best particle orientations from two parameter files.
@param 	*project1	first project structure with all parameters.
@param 	*project2	second project structure with all parameters.
@param 	*sym		symmetry.
@param 	origin_rmsd	cutoff to accept origins (pixels).
@param 	view_rmsd	cutoff to accept views (radians).
@param 	angle_rmsd 	cutoff to accept rotation angles (radians).
@param 	mag_rmsd	cutoff to accept magnifications (fraction).
@param	flag		0=RMSD, 1=MSD
@return long			number of particles selected.

	The error is calculated between parameters. 
	Requirement: The origin, view and rotation angle (or Euler angle), 
		and magnification arrays in the micrographs must be defined.

**/
long 		project_multi_select_low_rmsd(Bproject* project1, Bproject* project2,
				Bsymmetry& sym, double origin_rmsd, double view_rmsd,
				double angle_rmsd, double mag_rmsd, int flag)
{
	int				nsel(0);
	double			ddef, dori, dview, dang, dsize;
	double			def_avg_rmsd(0), origin_avg_rmsd(0), view_avg_rmsd(0), angle_avg_rmsd(0), size_avg_rmsd(0);
	
	if ( verbose )
		cout << "Selecting based on the lowest RMSD for symmetry " << sym.label() << endl;
	
	long   			i(0), nmg(0), npart(0);

	Bfield			*field1, *field2;
	Bmicrograph		*mg1, *mg2;
//	Breconstruction	*rec1, *rec2;
	Bparticle		*part1, *part2;
	Bparticle**		partarr1;
	Bparticle**		partarr2;
	View			view1, view2;
	
	if ( project1->select < 1 ) {
		nmg = project_count_micrographs(project1);
		partarr1 = project_mg_particle_array(project1, -1, npart);
		partarr2 = project_mg_particle_array(project2, -1, npart);
		for ( field1=project1->field, field2=project2->field; field1 && field2; field1=field1->next, field2=field2->next ) {
			if ( verbose & VERB_FULL )
				cout << field1->id << endl;
			for ( mg1=field1->mg, mg2=field2->mg; mg1 && mg2; mg1=mg1->next, mg2=mg2->next ) {
				if ( verbose & VERB_FULL )
					cout << mg1->id << endl;
				if ( mg1->ctf && mg2->ctf ) {
					ddef = mg1->ctf->defocus_average() - mg2->ctf->defocus_average();
					def_avg_rmsd += ddef*ddef;
				}
			}
		}
	} else {
		partarr1 = project_rec_particle_array(project1, -1, npart);
		partarr2 = project_rec_particle_array(project2, -1, npart);
	}
	
	for ( i=0; i<npart; i++ ) {
		part1 = partarr1[i];
		part2 = partarr2[i];
		if ( part1->sel && part2->sel ) {
			dori = (part1->ori - part2->ori).length();
			view1 = part1->view;
			view2 = find_closest_symmetric_view(sym, view1, part2->view);
			dview = view1.angle(view2);
			dang = fabs(angle_set_negPI_to_PI(view1.angle() - view2.angle()));
			dsize = part1->mag - part2->mag;
			if ( dori > origin_rmsd ) part1->sel = part2->sel = 0;
			if ( dview > view_rmsd ) part1->sel = part2->sel = 0;
			if ( dang > angle_rmsd ) part1->sel = part2->sel = 0;
			if ( dsize > mag_rmsd ) part1->sel = part2->sel = 0;
			if ( part1->sel > 0 ) {
				nsel++;
				origin_avg_rmsd += dori*dori;
				view_avg_rmsd += dview*dview;
				angle_avg_rmsd += dang*dang;
				size_avg_rmsd += dsize*dsize;
			}
		}
	}
	
	delete[] partarr1;
	delete[] partarr2;
	
	if ( nmg ) def_avg_rmsd = sqrt(def_avg_rmsd/nmg);
	
	if ( nsel ) {
		origin_avg_rmsd /= nsel;
		view_avg_rmsd /= nsel;
		angle_avg_rmsd /= nsel;
		size_avg_rmsd /= nsel;
		if ( !flag ) {
			origin_avg_rmsd = sqrt(origin_avg_rmsd);
			view_avg_rmsd = sqrt(view_avg_rmsd);
			angle_avg_rmsd = sqrt(angle_avg_rmsd);
			size_avg_rmsd = sqrt(size_avg_rmsd);
		}
	}
	
	if ( verbose & VERB_RESULT ) {
		if ( flag ) {
			cout << "Average defocus MSD:            " << def_avg_rmsd << " A" << endl;
			cout << "Number of particles selected with low MSD: " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
			cout << "Average origin MSD:             " << origin_avg_rmsd << " pixels" << endl;
			cout << "Average view MSD:               " << view_avg_rmsd*180.0/M_PI << " degrees" << endl;
			cout << "Average angle MSD:              " << angle_avg_rmsd*180.0/M_PI << " degrees" << endl;
			cout << "Average magnification MSD:      " << size_avg_rmsd << endl << endl;
		} else {
			cout << "Average defocus RMSD:           " << def_avg_rmsd << " A" << endl;
			cout << "Number of particles selected with low RMSD: " << nsel << " (" << nsel*100.0/npart << " %)" << endl;
			cout << "Average origin RMSD:            " << origin_avg_rmsd << " pixels" << endl;
			cout << "Average view RMSD:              " << view_avg_rmsd*180.0/M_PI << " degrees" << endl;
			cout << "Average angle RMSD:             " << angle_avg_rmsd*180.0/M_PI << " degrees" << endl;
			cout << "Average magnification RMSD:     " << size_avg_rmsd << endl << endl;
		}
	}
	
	return nsel;
}

/**
@brief 	Calculates the number of selections that are the same.
@param 	*project1	first project structure with all parameters.
@param 	*project2	second project structure with all parameters.
@return long		number of particles selections the same.

**/
long 		project_multi_selection_compare(Bproject* project1, Bproject* project2)
{
	if ( verbose )
		cout << "Comparing the selection for two projects: " << endl;
	
	long   			i(0), npart(0), npart2(0), nsame(0);

	Bparticle		*part1, *part2;
	Bparticle**		partarr1;
	Bparticle**		partarr2;
	
	if ( project1->select < 1 ) {
		partarr1 = project_mg_particle_array(project1, -1, npart);
		partarr2 = project_mg_particle_array(project2, -1, npart2);
	} else {
		partarr1 = project_rec_particle_array(project1, -1, npart);
		partarr2 = project_rec_particle_array(project2, -1, npart2);
	}
	
	if ( npart2 != npart ) {
		cerr << "Error: The number of particles must be the same in both projects:" << endl;
		cerr << "\tnumber[1] = " << npart << "\tnumber[2] = " << npart2 << endl;
		return 0;
	}
	
	for ( i=0; i<npart; i++ ) {
		part1 = partarr1[i];
		part2 = partarr2[i];
		if ( part1->sel == part2->sel ) nsame++;
		else if ( verbose )
			cout << part1->rec->id << tab << part1->id << tab << part1->sel << tab << part2->sel << endl;
	}
	
	delete[] partarr1;
	delete[] partarr2;
	
	if ( verbose & VERB_RESULT )
		cout << "Selections the same:            " << nsame << " (" << nsame*100.0/npart << "%)" << endl << endl;
	
	return nsame;
}

