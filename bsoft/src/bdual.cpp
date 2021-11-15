/**
@file	bdual.cpp
@brief	Program to handle dual tomographic tilt series
@author  Bernard Heymann, Jessica Mavadia, Charith Jayasinghe
@date	Created: 20071114
@date	Modified: 20160617
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "mg_tomo_track.h"
#include "ps_marker.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bdual [options] input1.star [input2.star]",
"------------------------------------------------",
"Handles dual tomographic tilt series of micrographs.",
"Expects exactly two input parameter files corresponding to the two tilt series.",
"Outputs only the second series to a parameter file.",
" ",
"Actions:",
"-seed 1.3,-95,75         Transfers seed markers from first to second series,",
"                         using the angular increment to test for rotation angle,",
"                         within the angular limits (default -180,180).",
"-refine                  Refines markers from the second series after transfer.",
"-transform               Calculates the correct micrograph transformations for the 2nd series.",
"-deselect 5-9,33         Deselect markers by id from all micrographs and models.",
"-zcompare                Compare Z coordinates of the first and second series.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 900,300      High and low resolution limits for cross-correlation (default 0.1,1000 angstrom).",
"-shiftlimit 155.5        Limit on micrograph shift search in pixels (default 20% of micrograph width).",
" ",
"Output:",
"-output file.star        Output parameter file for 2nd series.",
"-Postscript out.ps       Postscript output with Z coordinate error plots.",
"-bild out.bld            BILD output with fitted marker models.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	double			rot_start(-180), rot_end(180);	// Rotation search limits
	double			rot_step(0);					// Rotation search step size for seed transfer
	int				refine(0);						// Flag to refine seed markers
	int				transform(0);					// Flag to transform micrograph parameters
	Bstring			deselect;						// String specifying markers to deselect
	int				zcompare(0);					// Flag to enable Z compare option
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	double			hi_res(0), lo_res(1000);		// Default resolution range
	double			shift_limit(-1);				// Micrograph shift search limit
	Bstring			paramfile;						// Output parameter file
	Bstring			psfile;							// Output postscript file
	Bstring			bildfile;						// Output bild file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "seed" ) {
    	    if ( curropt->values(rot_step, rot_start, rot_end) < 1 )
				cerr << "-rotation: A rotation angle step size must be specified." << endl;
			else {
				rot_step *= M_PI/180.0;
				rot_start *= M_PI/180.0;
				rot_end *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "refine" ) refine = 1;
		if ( curropt->tag == "transform" ) transform = 1;
		if ( curropt->tag == "deselect" ) deselect = curropt->value;
		if ( curropt->tag == "zcompare" ) zcompare = 1 ;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( hi_res > lo_res )
				swap(hi_res, lo_res);
        }
		if ( curropt->tag == "shiftlimit" )
    	    if ( ( shift_limit = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A shift limit in pixels must be specified." << endl;
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
		if ( curropt->tag == "bild" )
			bildfile = curropt->filename();
	}
	option_kill(option);

	double		ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	if ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}
	
	Bproject*		project = read_project(file_list);
	if ( !project || !project->field || !project->field->mg ) {
		cerr << "Error: No first tilt series!" << endl;
		bexit(-1);
	}

	if ( !project->field->next ) {
		cerr << "Error: Two tilt series (fields) must be specified!" << endl;
		bexit(-1);
	}

	project_sort_markers_by_id(project);
	
	if ( deselect.length() )
		project_deselect_markers(project, deselect);

	if ( project->field->mg->origin[0] < 1 || project->field->mg->origin[1] < 1 )
		project_set_nominal_mg_origins(project);

	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	if ( rot_step > 0 ) {
		project_transfer_seed(project, rot_start, rot_end, rot_step, hi_res, lo_res, shift_limit);
		project->field->select = 0;
		if ( refine )
			project_refine_markers(project, hi_res, lo_res);
		project->field->select = 1;
	}
	
	if ( transform )
		project_transform_dual(project);

	if ( bildfile.length() )
		marker_sets_to_bild(bildfile, project->rec->mark, project->rec->next->mark);
	 
	if ( zcompare )
		project_dual_zcompare(project);
	
	if ( psfile.length() )
		ps_dual_zcompare(psfile, project);

//	project_tomo_errors(project2);

//	if ( transform )
//		project_merge_rec_markers(project);
	
//	if ( project->rec ) project_tomo_residuals(project, 1);

    if ( paramfile.length() )
		write_project(paramfile, project, 1, 1);
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

