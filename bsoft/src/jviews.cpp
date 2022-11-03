/**
@file	jviews.cpp
@brief	Generates views covering an asymmetric unit or all symmetry related views.
@author Juha Huiskonen
@author	Bernard Heymann
@date	Created:  20071203
@date 	Modified: 20080925
@date	Modified: 20081126 added shift option
@date	Modified: 20090204 added asymmetric unit views
@date	Modified: 20090205 added random views
@date	Modified: 20090206 added PDB output
@date	Modified: 20090408 added generate views option
@date	Modified: 20100527 added -Grid option
@date	Modified: 20110510 changed the orientation of the views when creating a Chimera marker file
@date	Modified: 20111208 added -Top option
@date	Modified: 20120123 (BH)
@date	Modified: 20120316 added -randomizealpha option
@date	Modified: 20150108 (BH) - incorporated into Bsoft

	Writes particle views into a Chimera markerfile for visualization
	Writes particle orientations into a PDB file as BIOIMT matrices 
	(these can be combined with a pseudoatomic model of the template, and all orientations 
	can be visualized in Chimera)
**/

#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"
#include "mg_processing.h"
#include "mg_extract.h"
#include "mg_subtomo.h"
#include "spline.h"
#include "random_numbers.h"

#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates PDB file matrices from particle views and locations.
@param 	*part     		list of particles.
@param	sampling		spatial sampling.
@param	&filename       PDB file with matrices.
@return	int  			0.
	Writes a PDB file.
	Uses the BIOMT records.
**/
int			pdb_matrices_from_particles(Bparticle* part, Vector3<double> sampling, Bstring& filename) {

	ofstream		fpdb;
	Matrix3			mat;
	Vector3<double>	shift;

	fpdb.open(filename.c_str());
	if ( fpdb.fail() ) return 0;
    
	for ( ; part; part = part->next ) {

		mat = (part->view.backward()).matrix();
		shift = (part->loc - part->ori) * sampling;	

		for ( int i=0; i<3; i++ )
			fpdb << fixed << "REMARK 350   BIOMT" << i+1 << setw(4) << part->id
				<< setprecision(6) << setw(10) << mat[i][0]
				<< setw(10) << mat[i][1]
				<< setw(10) << mat[i][2]
				<< setprecision(5) << setw (15) << shift[i] << endl;
	}
	
	fpdb.close();
	
	return 1;
}

int			cmm_views(Bstring& filename, Bparticle* part, Vector3<double> sampling,
				long part_sel, double fom_cutoff, double colorfom_min,
				double colorfom_max, RGB<double> color, double diameter, double length, double thickness)
{
	ofstream		fcmm;

	fcmm.open(filename.c_str());
	if ( fcmm.fail() ) return 0;
    
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Creating a Chimera marker file. Only the first datablock is used." << endl;
	}

	Bstring			id = filename.base();
	fcmm << "<marker_set name=\"" << id << "\">" << endl;

	int				n(0);
	double			fom;
	View			view;
	Vector3<double>	loc;

//		cout << part_sel << tab << fom_cutoff << endl;
	
	for ( ; part; part = part->next ) {

		fom = part->fom[0];
		view = part->view;
		view = view.backward();
		view.normalize();
		loc = part->loc * sampling;
		
//		cout << part->id << tab << part->sel << tab << fom << endl;

		if ( ((part_sel < 0) && (fom_cutoff < 0 )) || (part->sel == part_sel) || ((fom_cutoff > 0) && (fom >= fom_cutoff)) ) {

			if ( colorfom_min >= 0 ) {
				// truncate between min and max
				if ( fom < colorfom_min ) { fom = colorfom_min; }
				if ( fom > colorfom_max ) { fom = colorfom_max; }
	
				// scale between 0 and 1
				fom = fom - colorfom_min;
				fom = fom / (colorfom_max - colorfom_min);
	
				// calculate rgb from the fom value
				color[0] = fom;
				color[1] = 0;
				color[2] = 1 - fom;
			}

			// put marker at the origin
			fcmm << "<marker id=\"" << part->id << "99999\" x=\"" << loc[0]
					<< "\" y=\"" << loc[1]
					<< "\" z=\"" << loc[2]
					<< "\" r=\"" << color.r() << "\" g=\"" << color.g() << "\" b=\"" << color.b()
					<< "\" radius=\"" << thickness/2 << "\"/>" << endl;

			// marker at the end of the current view vector
			fcmm << "<marker id=\"" << part->id << "\" x=\"" << loc[0] + view[0] * length
					<< "\" y=\"" << loc[1] + view[1] * length
					<< "\" z=\"" << loc[2] + view[2] * length
					<< "\" r=\"" << color.r() << "\" g=\"" << color.g() << "\" b=\"" << color.b()
					<< "\" radius=\"" << diameter/2 << "\"/>" << endl;

			//link this point to the origin
			fcmm << "<link id1=\"" << part->id << "\" id2=\"" << part->id
					<< "99999\" r=\"" << color.r() << "\" g=\"" << color.g() << "\" b=\"" << color.b()
					<< "\" radius=\"" << thickness/2 << "\"/>" << endl;

			n++;
		}

	}
	
	fcmm << "</marker_set>" << endl;

	fcmm.close();
		
	cout << n << " markers written to " << filename << endl;

	return n;
}

// Usage assistance
const char* use[] = {
" ",
"Usage: jviews [options] [input.star]",
"------------------------------------",
"Creates and modifies particle views.",
" ",
"Parameters:",
"-Even                      Output evenly distributed views.",
"-Asymmetricunit            Output all views within asymmetric unit.",
"-Grid 10                   Output all views within angular distance from a given reference view.",
"-Related                   Output all symmetry related views for a given reference view and symmetry.",
"-Radial 75,75,75           Generate a view vector for each particle extending from the given point.",
"-Top 200                   Assing one of the two top views (up/down) for each particle, depending if the particle is above or below the chosen Z slice.",
"-Segment                   Generate a view vector for each particle extending from the closest point on a line segment (defined as a filament in the input.star).",
"-Filament 100,10,5         Generate views covering a filament with a given radius, helical rise & twist (use with option -symmetry)",
"                           and input.star defining the ends of a filament.",
"-Random 100                Generate a number of random views.",
"-randomizealpha            Randomize the angle around view vector (alpha).",
"-symmetry C10              Symmetry to use (default: C1 symmetry).",
"-recfile 3d.map            Associate views with this reconstruction filename.",
"-angle 2                   Angular sampling stepsize for theta and phi (default: 45 degrees).",
"-radius 75                 Place particles views at given radius (default 0).",
"-shift 5                   Shift particles by given amount along the view vectors (default 0).",
"-verbose 7                 Verbosity of output.",
"-reference 0,0.618,1,0     Reference view (dafault: 0,0,1,0)",
"-sampling 1.5,1.5,1.5      Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44          View origin (default: 0,0,0).",
"-discardz 20               Discard views with Z locations further from the Z origin  (pixels).",
"-select 1                  Process only particles with this selection number (default: all).",
"-invert                    Invert views in outputfile.",
" ",
"Chimera parameters:",
"-diameter 5                Chimera marker diameter (default: 2)",
"-length 100                Chimera link length (default: 10).",
"-thickness 5               Chimera link thickness (default: 2)",
"-color 0,0,1               Chimera marker & link color in RGB (default: 1,1,1).",
"-fomcolor 0.2,0.3          Chimera marker & link coloring range from blue to red based on particle FOM.",
" ",
"Output:",
"-output output.star        Output parameter file",
"-cmm output.cmm            Output Chimera marker file.",
"-pdb output.pdb            Output PDB file with BIOMT matrixes.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				mode(0);
	int 			set_origin(0);
	int				part_sel(-1);
	int				randomize_alpha(0);
	int				nviews(0);
	int 			invert(0);
	int				shift(0);
	int				locz(0);
	double			filament_ang(0);
	double			filament_rise(0);
	Vector3<double> sampling(1,1,1);
	Bstring			symmetry_string = "C1";
	double			angle(M_PI_4);
	Bstring			outcmm;
	Bstring			outstar;
	Bstring			outpdb;
	Bstring			recfile;

	double			length(10);
	double			angle_limit(0);
	double			diameter(2);
	double			thickness(2);
	double			fom_cutoff(-1);	
	double			colorfom_min(-1), colorfom_max(-1);
	double			rad(0);
	double			shiftparticle(0);
	double			discardz(-1);
	Vector3<double>	origin(0,0,0);
	Vector3<double>	vieworigin(0,0,0);
	Bproject*		project = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*		part = NULL;
	Bparticle*		part2 = NULL;
	Bparticle*		part_original = NULL;
	RGB<double>		color(1,1,1);
	Vector3<double>	loc;
	Vector3<double>	loc1;
	Vector3<double>	loc2;
	View			view;
	View			view_inv;
	View			std_view;
	View			ref_view;
	View*			v;
	Quaternion		q1, q2;

	int	           	i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {

		if ( curropt->tag == "Even") { mode = 1; }

		if ( curropt->tag == "Asymmetricunit") { mode = 1; }

		if ( curropt->tag == "Related") { mode = 2; }
		if ( curropt->tag == "symmetry") {
			symmetry_string = curropt->symmetry_string(); 
		}
		if ( curropt->tag == "Random") {
			if ( ( nviews = curropt->value.integer() ) < 1 ) {
				cerr << "-Random: A number of views must be specified!" << endl;
				bexit(-1);
			}
			else {
				mode = 3;
			}
		}
		if ( curropt->tag == "Radial") {
			vieworigin = curropt->origin();
			mode = 4;
		}
		if ( curropt->tag == "Grid") {
			if ( ( angle_limit = curropt->value.real() ) < 0.001 ) {
				cerr << "-Grid: An angular distance must be specified!" << endl;
				bexit(-1);
			}
			else {
				angle_limit *= M_PI/180.0;
				mode = 5;
			}
		}
		if ( curropt->tag == "Filament") {
			if ( curropt->values(shiftparticle, filament_rise, filament_ang) < 3 ) {
				cerr << "-Filament: Filament radius, angular spacing and angular raise must be specified!" << endl << endl;
				bexit(-1);
				shift = 1;
			}
			else {
			 	filament_ang *= M_PI/180.0;
				mode = 6;
			}
		}
		if ( curropt->tag == "Segment") {
			mode = 7;
		}
		if ( curropt->tag == "Top") {
			if ( ( locz = curropt->value.integer() ) < 0 ) {
				cerr << "-Top: A value for the Z slice must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				mode = 8;
			}
		}
		if ( curropt->tag == "recfile")
			recfile = curropt->filename();
		if ( curropt->tag == "radius") {
			if ( ( rad = curropt->value.real() ) < 1 ) {
				cerr << "-radius: A value for radius must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				if ( rad < 0 ) { rad = 0; }
			}
		}
		if ( curropt->tag == "shift") {
			if ( ( shiftparticle = curropt->value.real() ) < 0.001 ) {
				cerr << "-shift: A value for shift must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				shift = 1;
			}
		}
		if ( curropt->tag == "sampling")
			sampling = curropt->scale();
		if ( curropt->tag == "origin") {
			origin = curropt->origin();
			set_origin = 1;
		}
		if ( curropt->tag == "reference")
			ref_view = curropt->view();

		if ( curropt->tag == "length") {
			if ( ( length = curropt->value.real() ) < 0.001 ) {
				cerr << "-length: A value for length must be specified!" << endl << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "angle") {
			if ( ( angle = curropt->value.real() ) < 0.001 ) {
				cerr << "-angle: An angular stepsize must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				angle *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "diameter") {
			if ( ( diameter = curropt->value.real() ) < 0.001 ) {
				cerr << "-diameter: A value for diameter must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				if ( diameter < 1 ) { diameter = 1; }
			}
		}
		if ( curropt->tag == "thickness") {
			if ( ( thickness = curropt->value.real() ) < 0.001 ) {
				cerr << "-thickness: A value for thickness must be specified!" << endl << endl;
				bexit(-1);
			}
			else {
				if ( thickness < 1 ) { thickness = 1; }
			}
		}
		if ( curropt->tag == "color") {
			if ( curropt->values(color[0], color[1], color[2]) < 3 ) {
				cerr << "-color: Three values for color must be specified!" << endl << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "fomcolor") {
			if ( curropt->values(colorfom_min, colorfom_max) < 2 ) {
				cerr << "-fomcolor: Two values for FOM range must be specified!" << endl << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "output")
			outstar = curropt->filename();
		if ( curropt->tag == "cmm")
			outcmm = curropt->filename();
		if ( curropt->tag == "pdb")
			outpdb = curropt->filename();
		if ( curropt->tag == "select") {
			if ( ( part_sel = curropt->value.integer() ) < 1 ) {
				cerr << "-select: A selection number must be specified!" << endl << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "fom") {
			if ( ( fom_cutoff = curropt->value.real() ) < 0.001 ) {
				cerr << "-fom: A threshold value must be specified!" << endl << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "invert")
			invert = 1;
		if ( curropt->tag == "randomizealpha")
			randomize_alpha = 1;
	}
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read a parameter file or create a new project
	Bstring* file = NULL;
	if (argv[optind]) {
		string_add(&file, argv[optind]);
		if ( file->contains(".star") || file->contains(".STAR") ) {
			project = read_project(file);
		}
	} else {
		project = project_create(0,1);
		rec=project->rec;
		part = particle_add(&rec->part, 1);
		part->view = std_view;
		part->ori = origin;
		part->loc = origin;
	}
	string_kill(file);

	ref_view.normalize();

	Bsymmetry		sym(symmetry_string);
	Matrix3			mat;
	Matrix3			mat2;
	Quaternion		quat;
	Quaternion		quat_inv;
	Quaternion		quat_ref;

	// Views withing an asymmetric unit
	if ( mode == 1 ) {

		View*	allviews;

		if ( verbose & VERB_PROCESS ) {
			cout << endl << "Creating all views within asymmetric unit:" << endl;
		}

		allviews = asymmetric_unit_views(sym, angle, angle, 1);

		rec=project->rec;
		i=1;
		part_original = rec->part;
		rec->part = NULL;

		for ( v=allviews; v; v=v->next ) {
			view = *v;

			part = particle_add(&rec->part, i);
			part->loc = part_original->loc;
			part->ori = part_original->ori;
			part->view = view.backward();
			part->sel = part_original->sel;
			part->fom[0] = part_original->fom[0];
			part->fpart = part_original->fpart;
	
			i++;
		}
	}

	// All symmetry related views
	if ( mode == 2 ) {

		View* allviews;

		if ( verbose & VERB_PROCESS ) {
			cout << endl << "Creating all symmetry related views:" << endl;
		}

		for ( rec=project->rec; rec; rec=rec->next ) {
			i=1;
			part_original= rec->part;
			rec->part = NULL;

			for ( part = part_original; part ; part=part->next ) {

				mat2 = (ref_view.matrix()).transpose();

				view = part->view;
				view.normalize();

				if ( verbose & VERB_PROCESS ) {
					cout << "Using the reference view: " << view << endl;
				}

				allviews = symmetry_get_all_views(sym, view);

				for ( v=allviews; v; v=v->next ) {	
					mat = (v->matrix()).transpose();
					mat = mat * mat2;
					mat = mat.transpose();
					view = View(mat);

					part2 = particle_add(&rec->part, i);
					part2->loc = part->loc;
					part2->ori = part->ori;
					part2->view = view;
					part2->sel = part->sel;
					part2->fom[0] = part->fom[0];
					part2->fpart = part->fpart;
			
					if ( verbose & VERB_PROCESS ) {
						cout << tab << view << endl;
					}
					
					i++;
				} 
			}
		}
	}

	// Random views
	if ( mode == 3 ) {

		View*	allviews;

		rec=project->rec;
		i=1;
		part_original = rec->part;
		rec->part = NULL;

		if ( verbose & VERB_PROCESS ) {
			cout << endl << "Creating " << nviews << " random views:" << endl;
		}

		allviews = random_views(nviews);

		for ( v=allviews; v; v=v->next ) {
			view = *v;

			if ( verbose & VERB_PROCESS ) {
				cout << view << endl;
			}

			part = particle_add(&rec->part, i);
			part->loc = part_original->loc;
			part->ori = part_original->ori;
			part->view = view;
			part->sel = part_original->sel;
			part->fom[0] = part_original->fom[0];
			part->fpart = part_original->fpart;

			i++;	
		}
	} 

	// Views extending radially from a given point
	if ( mode == 4 ) {

		if ( verbose & VERB_PROCESS ) {
			cout << endl << "Assigning views extending radially from the given point." << endl;
		}

		for ( rec=project->rec; rec; rec=rec->next ) {
			i=1;

			for ( part=rec->part; part; part=part->next ) {

				view = View(part->loc[0] - vieworigin[0], part->loc[1] - vieworigin[1], 
					part->loc[2] - vieworigin[2], 0);

				view.normalize();
				view = view.backward();

				if ( verbose & VERB_PROCESS ) {
					cout << view << endl;
				}

				part->view = view;

				i++;
			}
		}
	} 

	// Views within specified angular distance
	if ( mode == 5 ) {

                View*      allviews;

		if ( verbose & VERB_PROCESS ) {
			cout << endl << "Creating all views within specified angular distance:" << endl << endl;
		}

                for ( rec=project->rec; rec; rec=rec->next ) {
	                i=1;
			part_original = rec->part;
			rec->part = NULL;
        
                        for ( part=part_original; part; part=part->next ) {
				view = part->view;
				view.normalize();

				if ( verbose & VERB_PROCESS ) {
					cout << "Using the reference view: " << view << endl;
				}

				allviews = views_within_limits(view, angle, angle, 1, angle_limit, 0);

				for ( v=allviews; v; v=v->next ) {
					view = *v;

					part2 = particle_add(&rec->part, i);
					part2->loc = part->loc;
					part2->ori = part->ori;
					part2->view = view;
					part2->sel = part->sel;
					part2->fom[0] = part->fom[0];
					part2->fpart = part->fpart;

					if ( verbose & VERB_PROCESS ) {
						cout << view << endl;
					}
					
					i++;
				}
			} 
                }
        }

	// Views for a filamentous object
        if ( mode == 6 ) {

                View*      allviews;
                Bfilament  fil;
                Bfilnode   node;

                if ( verbose & VERB_PROCESS ) {
                        cout << endl << "Creating all views for a filamentous object:" << endl;
                }

                mat2 = View(0,1,0,0).matrix();

                for ( rec=project->rec; rec; rec=rec->next ) {
	                i=1;
					part_original = filaments_to_particles(rec->fil, sampling,
						Vector3<double>(filament_rise,filament_rise,filament_rise),
						filament_rise, filament_rise, filament_ang);

					rec->part = NULL;

					for ( part=part_original; part; part=part->next ) {
                                allviews = symmetry_get_all_views(sym, part->view);

                                for ( v=allviews; v; v=v->next ) {      

                                        mat = v->matrix();
                                        mat = mat.transpose() * mat2.transpose();
										mat = mat.transpose();
                                        view = View(mat);

                                        part2 = particle_add(&rec->part, i);
                                        part2->loc = part->loc;
                                        part2->ori = part->ori;
                                        part2->view = view;
                
                                        if ( verbose & VERB_PROCESS ) {	cout << view << endl; }
                                        i++;
                                }
                        }
                }

        }

	// Views extending from the closest point on a line segment
        if ( mode == 7 ) {

              Bfilnode   node;
		Vector3<double> point1;
		Vector3<double> point2;

		for ( rec=project->rec; rec; rec=rec->next ) {
			i=1;

			point1 = rec->fil->node->loc;
			point2 = rec->fil->node->next->loc;

	                if ( verbose & VERB_PROCESS ) {
	                        cout << endl << "Assigning views extending from the closest point on a line segment from " << point1 << " to " << point2 << endl;
	                }

			for ( part=rec->part; part; part=part->next ) {

				// calculate the closest point on the line segment
				vieworigin = closest_point_line(part->loc, point1, point2);

				cout << endl << "Closest point:\t" << vieworigin << endl;

				view = View(part->loc[0] - vieworigin[0], part->loc[1] - vieworigin[1], part->loc[2] - vieworigin[2], 0);

				view.normalize();
				view = view.backward();

				if ( verbose & VERB_PROCESS ) {
					cout << view << endl;
				}

				part->view = view;

				i++;
			}
		}

        }

	// Top views (up / down)
        if ( mode == 8 ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {	
				if ( part->loc[2] < locz ) { part->view = View(0,0,-1,0); }
				else { part->view = View(0,0,1,0); }
			}
		}
	}

	// Randomize alpha
	if ( randomize_alpha ) {
		double irm = 2.0/get_rand_max();
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {
				if ( verbose & VERB_PROCESS ) {
					cout << endl << "Original view" << endl;
					cout << part->view << endl;
				}
				part->view = View(part->view[0],part->view[1],part->view[2],M_PI*(random()*irm - 1));
				if ( verbose & VERB_PROCESS ) {
					cout << endl << "View with random alpha" << endl;
					cout << part->view << endl;
				}
			}
		}
	}

	// Invert views
	if ( invert ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {
				if ( verbose & VERB_PROCESS ) {
					cout << endl << "Original view" << endl;
					cout << part->view << endl;
				}
				part->view = part->view.backward();
				if ( verbose & VERB_PROCESS ) {
					cout << endl << "Inverted view" << endl;
					cout << part->view << endl;
				}
			}
		}
	}

	// Shift particles along view vectors to place them one given radius
	if ( rad > 0 ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {
				view = part->view.backward();
				view.normalize();
				origin = part->ori;
				part->loc = Vector3<double>((origin[0] + rad * view[0]),
											(origin[1] + rad * view[1]),
											(origin[2] + rad * view[2]));
			}
		}
	}

	// Shift particles along view vectors
	if ( shift > 0 ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {
				view=part->view;
				view.normalize();
				view = view.backward();
				loc = part->loc;
				part->loc = Vector3<double>((loc[0] + shiftparticle * view[0]),
											(loc[1] + shiftparticle * view[1]),
											(loc[2] + shiftparticle * view[2]));
			}
		}
	}


	// Set origins
	if ( set_origin == 1 ) {
		if ( verbose & VERB_PROCESS ) {	cout << endl << "Setting particle origins." << endl; }		
		for ( rec=project->rec; rec; rec=rec->next ) {
			for (part = rec->part; part; part=part->next ) { part->ori = origin; }
		}
	}

	// Dicard particles which have Z locations too far from the central slize
	if ( discardz > 0 ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			for ( part=rec->part; part; part=part->next ) {
				view=part->view;
				origin = part->ori;
				loc = part->loc;
                                if (loc[2] < (origin[2] - discardz)) {part->sel = 0;}
                                if (loc[2] > (origin[2] + discardz)) {part->sel = 0;}
			}
		}
	}

	if ( outcmm.length() )
		if ( !cmm_views(outcmm, project->rec->part, sampling,
				part_sel, fom_cutoff, colorfom_min,
				colorfom_max, color, diameter, length, thickness) )
			cerr << "Cannot open cmm file for writing!" << endl << endl;

	if ( outpdb.length() )
		if ( !pdb_matrices_from_particles(project->rec->part, sampling, outpdb) )
			cerr << "Cannot open PDB file for writing!" << endl << endl;

	if ( recfile.length() ) {
		for ( rec=project->rec; rec; rec=rec->next ) {
			rec->frec = recfile;
		}
	}

	if ( outstar.length() ) {
		write_project(outstar, project, 0, 0);
	}

	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
