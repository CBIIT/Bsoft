/**
@file	bheight.cpp
@brief	Projecting a 3D map and calculating comparison statistics of the projections.
@author Bernard Heymann
@date	Created: 20170613
@date	Modified: 20170613
**/

#include "rwimg.h"
#include "symmetry.h"
#include "ps_views.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bheight [options] input.img output.img",
"---------------------------------------------",
"Calculates 2D height images from a 3D map.",
" ",
"Actions:",
"-axis y                  Project down a major axis (x, y, z or 0, 1, 2).",
"-View 0.3,-0.5,0.8,33    View to generate symmetry-related projections.",
"-Tilt -65,72,2.5,55      Tilt series of projections with angle min, max, step and axis angle.",
"-side 15                 Generate side view projections within the given angle from the equator.",
"-random 24               Generate a number of random projections.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 10,-10,20        Origin for rotation (voxels, default 0,0,0).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-symmetry D6             Symmetry: Point group identifier.",
"-angles 8,6              Step sizes for theta and phi in the asymmetric unit, one value sets both.",
"-threshold 1.2           Density threshold to indicate boundary of object (default 0).",
" ",
"Output:",
"-Plotviews plot.ps       Output postscript file with a plot of projection vectors.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType		nudatatype(Unknown_Type);		// Conversion to new type
	char			axis('n');						// Projection axis
	Vector3<double>	origin;				// Origin
	int				set_origin(0);					// Flag to set origin
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Bstring			symmetry_string("C1");			// Default: asymmetric or C1 point group
	double			theta_step(0);					// Angular step size for theta
	double			phi_step(0);					// Angular step size for phi
	double			ang_min(0), ang_max(0), ang_step(0), ang_axis(0); // Tilt series
	double			side_ang(-1);					// Side view variation angle
	int				nviews(0);						// Number of views
	int				view_flag(1);					// Flag for projection generation
	View			theview;						// View to generate symmetry-related projections
	int				symviews(0);					// Flag to generate symmetry-related projections
	double			threshold(0);					// Density threshold
	int				ps_flag(0);						// Flag for postscript output
	Bstring			ps_file;
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( curropt->tag == "axis" )
			if ( ( axis = curropt->value[0] ) < 'a' )
				cerr << "-axis: An axis for projection must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "random" )
			if ( ( nviews = curropt->value.integer() ) < 1 )
				cerr << "-random: A number of views must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				if ( i < 2 ) phi_step = theta_step;
				else phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "View" ) {
			theview = curropt->view();
			symviews = ps_flag = 1;
		}
		if ( curropt->tag == "Tilt" ) {
			if ( curropt->values(ang_min, ang_max, ang_step, ang_axis) < 3 )
				cerr << "-Tilt: Three angles must be specified!" << endl;
			else {
				ang_min *= M_PI/180.0;
				ang_max *= M_PI/180.0;
				ang_step *= M_PI/180.0;
				ang_axis *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "side" ) {
			if ( ( side_ang = curropt->value.real() ) < 1 )
				cerr << "-side: One angle must be specified!" << endl;
			else
				side_ang *= M_PI/180.0;
		}
		if ( curropt->tag == "threshold" )
			if ( ( threshold = curropt->value.real() ) < 0 )
				cerr << "-threshold: A positive value must be specified!" << endl;
		if ( curropt->tag == "Plotviews" )
			ps_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	Bsymmetry 	sym(symmetry_string);
	View*		views = NULL;
	
	if ( axis == 'n' ) {
		if ( nviews ) {
			views = random_views(nviews);
		} else if ( side_ang > -1 ) {
			views = side_views(sym, side_ang, theta_step, phi_step);
		} else if ( ang_step ) {
			views = tilt_views(ang_min, ang_max, ang_step, ang_axis);
		} else if ( symviews ) {
			views = symmetry_get_all_views(sym, theview);
		} else if ( theta_step && phi_step ) {
			views = asymmetric_unit_views(sym, theta_step, phi_step, view_flag);
			ps_flag = 1;
		}
		if ( views ) {
			nviews = count_list((char *) views);
			if ( verbose )
				cout << "Generating " << nviews << " projections" <<endl << endl;
			if ( ps_file.length() )
				ps_views(ps_file, symmetry_string, views, ps_flag);
		}
	} else if ( axis == 'z' ) {
		views = new View(0,0,1,0);
	} else if ( axis == 'y' ) {
		views = new View(0,1,0,0);
	} else if ( axis == 'x' ) {
		views = new View(1,0,0,0);
	}
	
	if ( optind >= argc ) {
		if ( views ) kill_list((char *) views, sizeof(View));
		bexit(0);
	}

	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();		// Preserve the old type
	
	if ( sam.volume() ) p->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
	
	p->calculate_background();
	
	Bimage* 	pheight = p->height(views, threshold);
	
	if ( views ) kill_list((char *) views, sizeof(View));
	
	if ( optind < argc ) {
		pheight->change_type(nudatatype);
		write_img(argv[optind], pheight, 0);
	}
	
	delete p;	
	delete pheight;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


