/**
@file	dmerge.cpp
@brief	Merge images in a defocal or other series, align images before adding together
@author David Belnap & Bernard Heymann
@date	Created: 20021004
@date	Modified:  20041124 (DB)
**/

#include "mg_merge.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Global variables
extern int   verbose;     // Level of output to the screen

const char* use[] = {
" ",
"Usage: dmerge [options] input1.star [input2.star...]",
"----------------------------------------------------",
"Merges corresponding single particle images in a series.",
"The images can be aligned according to their origin values and",
"the in-plane rotational angle (of the micrograph).",
" ",
"You may use the index number of a micrograph in a series, a relative focus",
"level, or an in-plane rotation angle of the micrograph to select a reference",
"micrograph--the one to which the others are aligned.  The default is to use",
"the closest-to-focus micrograph in each field-of-view.",
" ",
"The default is to use origins from parameter file(s), while origins may also",
"be determined by cross-correlation or obtained from particle image files.",
" ",
"Default image output is to one file for every field named",
"'[reference_particle_file_prefix]_merged.img'.",
" ",
"Actions:",
"-unmerge merged.star     Un-merge, corresponding particles in field get orientations,",
"                         figure-of-merit (FOM), and selection from 'merged.star',",
"                         micrograph rotation angle is accounted for",
" ",
"Parameters:",
"-verbose 7               Verbosity of screen output, default = 0 (none)",
"-sampling 1.5            Global sampling (angstrom/pixel)",
" ",
"Parameters for merging:",
"-origins cross           Get particle origins from parameter files (param, default),",
"                         by cross-correlation (cross), from image files (image),",
"                         or no alignment (none)",
"-reference near          Reference micrograph: near, nearest-to-focus; far, farthest-from-focus;",
"                         index,1, micrograph index; or angle,0, micrograph rotation angle",
"                         Note: the index option assumes micrographs were entered",
"                         into the parameter file(s) in consistent order",
" ",
"Parameters for unmerging:",
"-keeporient 0.1          Keep orientation, FOM, & selection if FOM difference >= this value",
" ",
"Output:",
"-output file.star        Output parameter file",
"-summed output.img       Write summed images to output*.img (* = _micrograph_ids), but only",
"                         if micrograph_ids are <= 6 characters, if not default name is used ",
" ",
NULL
};



int   main(int argc, char** argv)
{
	Vector3<double>	sam;    		// For resetting the grid size
	Bproject*		orientations = NULL;	// Data for orientations to be applied to all in field
	double			fom_diff(-1.0);     	// FOM-difference threshold to not change (for un-merging option)
	int				mg_ori_select(0);		// Flag for getting origins (default use origins in parameter file)
	int				mg_ref_select(0);		// Reference micrograph flag (default is closest-to-focus)
	int				mg_index(0);			// Reference micrograph index
	double			mg_rot_ang(0);			// Reference micrograph rotation angle
	Bstring			orientfile;				// Parameter file with orientations to be applied to all in field
	Bstring			outimg;					// Filename prefix and suffix for output images
	Bstring			outfile;				// Output parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next )  {
		if ( curropt->tag == "origins" )  {  	// Source of origins
			if ( curropt->value.contains("no") )   mg_ori_select = -1;
			if ( curropt->value.contains("par") )  mg_ori_select = 0;
			if ( curropt->value.contains("cro") )  mg_ori_select = 1;
			if ( curropt->value.contains("im") )   mg_ori_select = 2;
		}
		if ( curropt->tag == "keeporient" )  {  // do not change orient. in certain cases (only with unmerge option)
			if ( ( fom_diff = curropt->value.real() ) < 0.0001 )
				cerr << "-keeporient: A FOM-difference threshold must be specified." << endl;
			if ( fom_diff <= 0 )
				cerr << "-keeporient: The FOM-difference threshold (" << fom_diff << ") must be greater than zero." << endl;
		}
		if ( curropt->tag == "sampling" )       // set grid size of images
			sam = curropt->scale();
		if ( curropt->tag == "reference" )  {      // reference micrograph selection
			if ( curropt->value.contains("near") ) mg_ref_select = 0;
			else if ( curropt->value.contains("far") ) mg_ref_select = 1;
			else if ( curropt->value.contains("ind") ) {
				mg_ref_select = 2;
				mg_index = curropt->value.post(',').integer();
			} else if ( curropt->value.contains("ang") ) {
				mg_ref_select = 3;
				mg_rot_ang = curropt->value.post(',').real();
				mg_rot_ang *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "summed" )         // output image file
			outimg = curropt->filename();
		if ( curropt->tag == "unmerge" )        // input file with orientations for all members of field
			orientfile = curropt->filename();
 		if ( curropt->tag == "output" )      // output parameter file
			outfile = curropt->filename();
	}
	option_kill(option);

	double       ti = timer_start();

	if ( orientfile.length() && (fom_diff > 0) )  {  // option to not change orientation, etc. values only works with un-merge option
		cerr << "ERROR:  You selected the 'unchange' option.  It only works if the unmerge option also is used." << endl;
		bexit(-1);
	}

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

	if ( sam[0] > 0.1 )			// Reset sampling size for project
		project_set_mg_pixel_size(project, sam);


	// Align and sum images and write to a file or "un-merge" parameter data
	if ( orientfile.length() )  {
		orientations = read_project(orientfile);
		mg_particle_unmerge(project, orientations, fom_diff);
	} else
		mg_particle_merge_series(project, mg_ref_select, mg_index, mg_rot_ang, mg_ori_select, outimg);


	// Output parameters
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}

	project_kill(project);

	if ( verbose & VERB_TIME )  timer_report(ti);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

