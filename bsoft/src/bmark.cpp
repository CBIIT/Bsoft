/**
@file	bmark.cpp
@brief	Program to generate symmetry axes for mark group symmetries
@author Bernard Heymann
@date	Created: 20020619
@date	Modified: 20150212
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_tomography.h"
#include "mg_tomo_track.h"
#include "mg_marker.h"
#include "marker.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bmark [options] in.star [in2.star...]",
"--------------------------------------------",
"Operates on sets of micrograph markers.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-list                    List all micrograph and reconstruction markers.",
"-rotate 0.3,0.1,-0.5,81  Rotate around an axis by an angle.",
"-scale 0.5,3,1           Scale.",
"-translate 12.2,0.5,-50  Shift after rotation and scaling.",
"-plane                   Calculate a plane through the markers.",
"-fom 0.55                Select markers with a figure-of-merit above this cutoff.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-origin 0,22.5,30        Set the origin for rotation.",
//"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
" ",
"Output:",
"-output out.star         Output micrograph parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	int				use_rec(0);					// Flag to process reconstructions
	int				list(0);
	int 			settransform(0);
	int				setorigin(0);
	int				calcplane(0);				// Flag to calculate a plane
	double			fom(0);						// FOM cutoff to select markers
	Transform		t;							// All transformation parameters
//	Vector3<double>	samp = {1,1,1};			// Pixel size
	Bstring			output;						// Output micrograph parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "list" ) list = 1;
		if ( curropt->tag == "rotate" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else {
				t.axis = Vector3<double>(d[0], d[1], d[2]);
				t.angle = d[3] * M_PI/180.0;
				settransform = 1;
			}
		}
		if ( curropt->tag == "scale" ) {
			t.scale = curropt->vector3();
			if ( t.scale.volume() < 1 )
				cerr << "-scale: At least one value must be specified!" << endl;
			else {
				if ( t.scale[1] <= 0 ) t.scale[2] = t.scale[1] = t.scale[0];
				settransform = 1;
			}
		}
		if ( curropt->tag == "translate" ) {
			t.trans = curropt->vector3();
			if ( t.trans.length() < 0.1 )
				cerr << "-translate: At least one translation value must be specified!" << endl;
			else
				settransform = 1;
		}
		if ( curropt->tag == "plane" ) calcplane = 1;
		if ( curropt->tag == "fom" )
			if ( ( fom = curropt->value.real() ) < 0.0001 )
				cerr << "-fom: A figure-of-merit must be specified!" << endl;
		if ( curropt->tag == "origin" ) {
			t.origin = curropt->origin();
			setorigin = 1;
		}
//		if ( curropt->tag == "sampling" )
//			sam = curropt->scale();
		if ( curropt->tag == "output" )
			output = curropt->filename();
    }
	option_kill(option);	
	
	double		ti = timer_start();
	
	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);

	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( use_rec ) project->select = 1;

	if ( setorigin ) project->rec->origin = t.origin;
	
	if ( settransform )
		marker_transform(project->rec->mark, t);
	
	if ( calcplane )
		marker_plane(project->rec->mark, project->rec->origin);
//		project_tilt_axis_from_markers(project);
	
	if ( fom > 0 )
		project_mg_marker_select(project, fom);
	
	if ( verbose ) 
		marker_stats(project->rec->mark);

	if ( list ) {
		project_marker_lists(project);
//		project_marker_in_particle(project);
		double marker_radius(10);
		project_marker_in_particle_image(project, marker_radius);
	}
	
	// write an output STAR format file if a name is given
    if ( output.length() ) {
		write_project(output, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


