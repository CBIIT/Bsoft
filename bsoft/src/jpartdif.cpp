/**
@file	jpartdif.cpp
@brief	Analysis of particle differences between two starfiles
@author Juha Huiskonen
@author	Bernard Heymann
@date	Created:  20071214
@date	Modified: 20111506
@date	Modified: 20120106
@date	Modified: 20150108 (BH) - incorporated into Bsoft
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_multiple.h"
#include "mg_particle_select.h"
#include "Matrix.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"
#include "Vector3.h"
#include "View.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: jpartdif [options] input1.star input2.star",
"-------------------------------------------------",
"Calculates the difference in the particle parameters,",
"or selects particles based on the fom difference.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-consolidate             Only consolidate particles based on distance.",
"-locdiff                 Calculate the difference in the locations.",
"-fomdiff                 Calculate the difference in the foms.",
"-printlocdiff            Print the difference in the locations.",
"-printanglediff          Print the difference in the view angle.",
"-printalphadiff          Print the difference in the alpha angles.",
"-thrlocdiff 10.0         Analyse only particles below threshold (pixels).",
"-thranglediff 20.0       Analyse only particles below threshold (degrees).",
"-remove                  Remove all particles with selection number 0.",
" ",
"Output:",
"-output file.star        Output parameter file.",

" ",
NULL
};

int			main(int argc, char** argv)
{
	// Initializing variables
	Bstring			outfile;					// Output parameter file
	int				optind;
	int				consolidate(0);
	int				locdiff(0), fomdiff(0);
	int				printlocdiff(0), printanglediff(0), printalphadiff(0);
	double			anglediffthr(-1), locdiffthr(-1);

	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "consolidate" )
		       	consolidate = 1;
		if ( curropt->tag == "locdiff" )
		       	locdiff = 1;
		if ( curropt->tag == "printlocdiff" )
		       	printlocdiff = 1;
		if ( curropt->tag == "printanglediff" )
		       	printanglediff = 1;
		if ( curropt->tag == "printalphadiff" )
		       	printalphadiff = 1;
		if ( curropt->tag == "fomdiff" )
		       	fomdiff = 1;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "thrlocdiff" )
			if ( ( locdiffthr = curropt->value.real() ) < 0.001 )
				cerr << "-thrlocdiff: A threshold must be specified!" << endl;
		if ( curropt->tag == "thranglediff" )
			if ( ( anglediffthr = curropt->value.real() ) < 0.001 )
				cerr << "-thranglediff: A threshold must be specified!" << endl;

	}
	
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read all the parameter files
	Bproject*		project1 = NULL;
	Bproject*		project2 = NULL;
	Bproject*		project = NULL;
	Breconstruction*	rec1 = NULL;
	Breconstruction*	rec2 = NULL;
	Bparticle*		part1 = NULL;
	Bparticle*		part2 = NULL;
	Bparticle*		part3 = NULL;
	Bparticle*		part = NULL;
	View			view1;
	View			view2;
	Vector3<double>	part1_view;
	Vector3<double>	part2_view;

	double			locresid(0), angleresid(0), alpharesid(0), sumlocresid2(0), sumangleresid2(0), sumalpharesid2(0), alpha1, alpha2;

	project1 = read_project(argv[optind++]);
	project2 = read_project(argv[optind++]);
	project = project_create(0,1);

	double d, dmin;
	int closest_part_id(-1), found(0), npart(0);

	if ( consolidate ) {
		for ( rec1 = project1->rec, rec2=project2->rec ; rec1 && rec2; rec1=rec1->next, rec2=rec2->next ) {
		for ( part1 = rec1->part; part1 ; part1 = part1->next ) {
			// find the closest particle (part2) for the current hit (part1)
			dmin=10000000;
			for ( part2 = rec2->part; part2 ; part2 = part2->next ) {
				d = part2->loc.distance(part1->loc);
				if ( d < dmin ) {
					dmin = d;
					closest_part_id = part2->id;
				}
			}
			// check if a particle was assigned to this hit already in the output structure (project)
			for ( part3 = project->rec->part ; part3 ; part3 = part3->next ) {
				if ( part3->id == closest_part_id ) {
					found = 1;
				}
				if ( found ) { break ; }
			}
			// if the hit corresponds to a new found particle, add to the outputs
			if ( !found ) {
				part=particle_add(&project->rec->part, closest_part_id);
				part->loc = part1->loc;
				part->ori = part1->ori;
				part->view = part1->view;
				part->fom[0] = part1->fom[0];
			}
			// if the hit corresponds to a particle that was assigned already, check if this hit is better and then update
			else {
				if (part1->fom[0] > part3->fom[0] ) {
					part3->loc = part1->loc;
					part3->ori = part1->ori;
					part3->view = part1->view;
					part3->fom[0] = part1->fom[0];
				}
			}
			found = 0;	
		}
		}
	}
	else {
		for ( rec1 = project1->rec, rec2=project2->rec ; rec1 && rec2; rec1=rec1->next, rec2=rec2->next ) {
		for ( part1 = rec1->part; part1 ; part1 = part1->next ) {
			// find the same particle in the other file (based on ID)
			for ( part2 = rec2->part; part2 ; part2 = part2->next ) {
				if ( part1->id == part2->id ) {
					part=particle_add(&project->rec->part, part1->id);
					//default values from particle 1
					part->loc = part1->loc;
					part->ori = part1->ori;
					part->view = part1->view;
					part->fom[0] = part1->fom[0];

					// calculate differences for selected particles

					if ( part1->sel > 0 ) {
						if ( locdiff ) { part->loc = part1->loc - part2->loc; }		
						if ( fomdiff ) { part->fom[0] = part1->fom[0] - part2->fom[0]; }

						view1 = View(part1->view[0],part1->view[1],part1->view[2],part1->view.angle());
						view2 = View(part2->view[0],part2->view[1],part2->view[2],part2->view.angle());

						view1 = view1.backward();
						view2 = view2.backward();
	
						locresid = part2->loc.distance(part1->loc);

						part1_view = Vector3<double>(view1[0], view1[1], view1[2]);
						part2_view = Vector3<double>(view2[0], view2[1], view2[2]);

						angleresid = part1_view.angle(part2_view);
						angleresid = angleresid *180.0/M_PI;

						//angleresid = view_difference((part1->view).backward(),(part2->view).backward())*180.0/M_PI;

						alpha1 = (view1.angle()*180/M_PI);
						alpha2 = (view2.angle()*180/M_PI);


						// warning C3 symmetry assumed here					
						while ( alpha1 < 0 ) { alpha1 += 120; }
						while ( alpha2 < 0 ) { alpha2 += 120; }

						while ( alpha1 >= 120 ) { alpha1 -= 120; }
						while ( alpha2 >= 120 ) { alpha2 -= 120; }

						//while ( alpha1 < 0 ) { alpha1 += 360; }
						//while ( alpha2 < 0 ) { alpha2 += 360; }

						alpharesid =  fabs(alpha1 - alpha2);
						if ( (120 - fabs(alpha1 - alpha2)) < alpharesid ) { alpharesid = 120 - fabs(alpha1 - alpha2); }
						//if ( (360 - fabs(alpha1 - alpha2)) < alpharesid ) { alpharesid = 360 - fabs(alpha1 - alpha2); }
		 
						if ( verbose & VERB_PROCESS ) { cout << "Part1 ID: " << part1->id << "\tPart2 ID: " << part2->id << tab << locresid << tab << angleresid << tab << alpharesid << tab << alpha1 << tab << alpha2 << endl; }

						if ( ( locdiffthr < 0 || locresid <= locdiffthr ) && ( anglediffthr < 0 || angleresid <= anglediffthr )) {
							npart++;
							sumlocresid2 += locresid * locresid;
							sumangleresid2 += angleresid * angleresid;
							sumalpharesid2 += alpharesid * alpharesid;
						}
						else {
							if ( verbose & VERB_PROCESS ) { cout << "Skipping particle " << part1->id << endl; }
						}
					}
				}
			}
		}
		}
	}

	if ( printlocdiff ) { cout << "Location RMSD\t" << sqrt(sumlocresid2/npart) << "\tParticles: " << npart << endl; }
	if ( printanglediff ) { cout << "Angle RMSD:\t" << sqrt(sumangleresid2/npart) << "\tParticles: " << npart << endl; }
	if ( printalphadiff ) { cout << "Alpha RMSD:\t" << sqrt(sumalpharesid2/npart) << "\tParticles: " << npart << endl; }
	
	if ( verbose & VERB_PROCESS )
		project_show_selected(project);
	
	// Write an output parameter file if a name is given
        if ( outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}

	project_kill(project1);
	project_kill(project2);
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

