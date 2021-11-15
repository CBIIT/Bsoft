/**
@file	dhand.cpp
@brief	Compare projections of two oppositely handed reconstructions to images of tilted specimens
@author  David Belnap
@date	Created: 20011003
@date	Modified:  20121013 (BH)
**/

#include "mg_hand.h"
#include "mg_tomography.h"
#include "mg_particle_select.h"
#include "rwmg.h"
#include "options.h"
#include "timer.h"
#include "utilities.h"

// Global variables
extern int   verbose;   // Level of output to the screen

// Internal function prototype
int   hand_get_options_mgselect(Bstring& tag, Bstring& tagvalue, int* mg_select, int* mg_index, double* mg_ang);


// Usage information
const char* use[] = {
" ",
"Usage: dhand [options] input1.star [input2.star...]",
"----------------------------------------------------------------------------",
"Compares projections of two oppositely handed reconstructions to images",
"of tilted specimens.  The orientation of the second view is predicted",
"from the first view's orientation and from the tilt-axis direction and the",
"rotation angle.  For particles from the two views, output selection values",
"are set to 1 (hand A), 2 (hand B), or 0 (unassigned, FOM < min. value or ",
"|FOMA-FOMB| < min. acceptable difference).  Particles with input selection",
"values of zero are ignored unless option '-all' is used.",
" ",
"Actions:",
"-difference 1            Compute |difference| between measured (input) and predicted",
"                         orientations, 0=do not (default), 1=pred. orient. and avg",
"                         |difference|, 2=#1 + individual values.",
"-org2 cross              Get particle origins of second view by cross-correlation",
"                         with projections of predicted view (cross), from image",
"                         files (image), or from parameter files (param, default).",
"-selconsist              Sets selection values in unused micrographs to be consistent",
"                         with those in views 1 and 2.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-all                     Include all particles (default sel>0).",
"-numberperfield 3        Number of micrographs per field-of-view (default 1).",
"-fieldID dbgd            Set the field ID (multiple micrographs will get numbered ID's).",
"-axis 74.7               Tilt axis angle relative to x-axis (default from parameter file).",
//"-angle -15               Tilt angle (default from parameter file).
"-origin 50,5.3,50.4      Set the x,y,z origin of 3D map, in pixel units.",
"-mindiff 0.05            minimum acceptable |difference| between hand A & B ",
"                         correlation coefficients, default = 0.050.",
"-minfom 0.1              min. acceptable value for FOM of correct hand, default=0.",
"-radii 5,65              Radial limits (in pixels), default = all radii.",
"-resolution 100,17       Resolution limits in angs. (low,high), default = 400,20.",
"-sampling 2.54           x,y,z pixel sampling (input map and all micrographs).",
" ",
"Parameters for selecting first and second views (use mgsel or both mg1 & mg2):",
"-mgsel angle,0,-5.0      Criterion, value for first view, value for second view",
"-mg1 index,1             Micrograph in each set that is the first view",
"-mg2 angle,-5.0          Micrograph in each set that is the second view",
"                         Options:",
"                           angle,a (tilt angle)",
"                           index,n (index within a series, beginning with 1)",
"                           near,n  (nth nearest-to-focus, beginning with 1)",
"                           far,n   (nth farthest-from-focus, beginning with 1)",
" ",
"Input:",
"-map input.img           Reference map.",
" ",
"Output:",
"-output output.star      Outputs data in STAR format to output.star",
"-outprj output.img       Outputs projections to output_field#_img#_imgID.img",
" ",
NULL
};



int   main(int argc, char* argv[])
{

	int             all(0);					// Flag to reselect all particles
	long 			nseries(0);				// Number of micrographs in a series from a field
	Bstring			field_id;				// User-specified field ID
	double			tilt_axis(-1000);		// Tilt axis angle
	double			tilt_angle(0);			// Tilt angle start
	double			AB_min(0);				// FOM for A or B must be this value or greater
	double			AmB_min(0.05);			// |FOMA - FOMB| must be this value or greater
	long            count,mcount,fcount;	// Placeholder for number of fields-of-view, micrographs, and particles
	int             diff_out(0);			// Flag to output difference between measured and predicted orientations
	int             handedness;				// Output to shell:  1=hand A, 2=hand B, 0=undetermined
	Bstring         inmap;					// Input 3D map (image-format) file
	int             MapOriginFlag(0);		// Flag to set origin of input 3D map
	double			mg_ang[2];				// values and flags for selection of micrographs corresponding to views 1 and 2
	int             mg_index[2];
	int             mg_select[2];
	int             mgsel(0);				// Flag to make sure proper number of mg selections is made
	int             origins2(0);			// Flag to determine origins for second view
	Vector3<double>	origind;				// X,Y,Z origin for 3D map
	Bstring         outimg;					// Prefix and suffix for output image files
	Bstring         outfile;				// Output parameter file
	Vector3<double>	sam;    		// User-entered pixel size for map and particle images
	Bproject*       project=NULL;			// Data for all micrographs
	Bimage          *penantiomer=NULL;		// 3D maps, pmirror is mirror image of penantiomer
	double			rad_min(0);				// minimum and maximum radius, given in pixels
	double			rad_max(-1.0);
	double			res_min(400);			// minimum and maximum resolution, given in angstroms
	double			res_max(20);
	int             sel_consist(0);			// Flag to set consistent selection values in all micrographs


	// Get and check command-line input, program terminates if there is a problem
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next )  {
		if ( curropt->tag == "all" )            // all particles will be processed
			all = 1;
		if ( curropt->tag == "numberperfield" )
			if ( ( nseries = curropt->value.integer() ) < 1 ) {
				cerr << "-numberperfield: The number of micrographs in a field-of-view must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "fieldID" ) {
			field_id = curropt->value;
			if ( field_id.length() < 1 )
				cerr << "-fieldID: A field ID must be specified!" << endl;
		}
		if ( curropt->tag == "axis" ) {
    	    if ( ( tilt_axis = curropt->value.real() ) < 0.00001 )
				cerr << "-axis: The tilt axis angle must be specified." << endl;
			else
				tilt_axis = angle_set_negPI_to_PI(tilt_axis*M_PI/180.0);
		}
		if ( curropt->tag == "angle" ) {
    	    if ( ( tilt_angle = curropt->value.real() ) < 0.00001 )
				cerr << "-angle: The tilt angle must be specified." << endl;
			else
				tilt_angle = angle_set_negPI_to_PI(tilt_angle*M_PI/180.0);
        }
		if ( curropt->tag == "difference" )  {    // difference meas. & pred.
			diff_out = curropt->value.integer();
			if ( (diff_out < 0) || (diff_out > 2) )  {
				cerr << "Error:  -difference, a difference output option (0-2) must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "map" )              // input 3D map
			inmap = curropt->filename();
		if ( curropt->tag == "origin" )  {     // set origin of input map
			origind = curropt->origin();
			MapOriginFlag = 1;
		}
		if ( curropt->tag == "mgsel" ||		// selection of 1st and 2nd micrographs (views)
				curropt->tag == "mg1" || 
				curropt->tag == "mg2" )  {
			if ( curropt->tag == "mgsel" )  mgsel += 2;
			else  mgsel += 1;
			if ( mgsel > 2 )  {
				cerr << "Error:  too many micrograph selection choices were entered!" << endl;
				bexit(-1);
			}
			hand_get_options_mgselect(curropt->tag, curropt->value, mg_select, mg_index, mg_ang);
		}
		if ( curropt->tag == "minfom" )  {         // min. allowable CC value
			if ( ( AB_min = curropt->value.real() ) < 0.0001 )  {
				cerr << "Error:  -minfom, threshold value must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "mindiff" )  {         // min. allowable difference between hand A & B
			if ( ( AmB_min = curropt->value.real() ) < 0.0001 )  {
				cerr << "Error:  -mindiff, threshold value must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "org2" )  {      // Source of origins for images of second view
			if ( curropt->value.contains("par") )       origins2 = 0;
			else if ( curropt->value.contains("cro") )  origins2 = 1;
			else if ( curropt->value.contains("im") )   origins2 = 2;
			else  {
				cerr << "Error:  -org2, " << curropt->value << " is an improper value for this option" << endl;
				bexit(-1);
			}
		}
 		if ( curropt->tag == "outprj" )  {        // Output projections
			outimg = curropt->filename();
		}
 		if ( curropt->tag == "output" )           // Output parameter file
			outfile = curropt->filename();
		if ( curropt->tag == "radii" )  {         // Radial parameters
			if ( curropt->values(rad_min, rad_max) < 2 )  {
				cerr << "Error:  -radii, radial limits (min,max) in pixels must be specified!" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "resolution" )  {    // Resolution parameters
			if ( curropt->values(res_min, res_max) < 2 )  {
				cerr << "Error:  -resolution, two resolution parameters must be specified!" << endl;
				bexit(-1);
			}
			else  {
				if ( res_max > res_min)  {
					swap(res_min, res_max);
					if (verbose)  cout << "Resolution max > resolution min, so the two were swapped" << endl;
				}
			}
		}
		if ( curropt->tag == "sampling" )          // set pixel size
			sam = curropt->scale();
		if ( curropt->tag == "selconsist" )     // set selection values in unused micrographs to those in views 1 and 2
			sel_consist = 1;
	}
	option_kill(option);
	
    double			ti = timer_start();

	if ( mgsel != 2 )  {
		cerr << "Error:  no or not enough micrograph selection input was entered." << endl;
		bexit(-1);
	}

	// Read all the parameter files
	Bstring*   file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}
	project = read_project(file_list);
	string_kill(file_list);
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	// Select all particles
	if ( all )   part_reset_selection(project, 3);

	if ( nseries )
		project_set_field_id(project, nseries, field_id);

	if ( tilt_axis )
		project_set_tilt_axis(project, tilt_axis);

	// Data count, there must be at least two micrographs in each field-of-view
	fcount = project_count_fields(project);
	mcount = project_count_micrographs(project);
	if ( mcount < 2*fcount )  {
		cerr << "Error:  there must be at least two micrographs in each field-of-view:" << endl;
		cerr << "Number of micrographs = " << mcount << "     Number of fields-of-view = " << fcount << endl;
		bexit(-1);
	}
	if ( verbose < VERB_PROCESS )  {
		cout << "Fields-of-view:                 " << fcount << endl;
		cout << "Micrographs:                    " << mcount << endl;
		count = project_count_mg_particles(project);
		cout << "Particles in all micrographs:   " << count << endl;
		count = project_count_mg_part_selected(project);
		cout << "Selected particles:             " << count << endl << endl;
	}


	// Read input enantiomer (hand A map),
	if (verbose)  cout << "Map file name (hand A):         \"" << inmap << "\"" << endl;
	penantiomer = read_img(inmap, 1, -1);


	// Set map origin if entered
	if ( MapOriginFlag )		// Map origin set to command-line input value
		penantiomer->origin(origind);


	// Set pixel size for input map and all micrographs, if entered
	if ( sam.volume() > 0.0 ) {
		penantiomer->sampling(sam);
		project_set_mg_pixel_size(project, sam);
	}


	// Set radial max. cut-off to max. radius of map, if it wasn't entered
	if (rad_max < 0.0)
		rad_max = penantiomer->maximum_included_radius();


	// Determine handedness
	handedness = project_get_handedness(penantiomer, project, mg_ang, mg_index, mg_select, rad_min, rad_max,
										res_min, res_max, AmB_min, AB_min, diff_out, origins2, outimg);


	// Set consistent selection values for all micrographs in a field-of-view
	if ( sel_consist )
		hand_select_consist(project, mg_ang, mg_index, mg_select, sel_consist);


	// Output parameters
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}


	// Memory cleanup
	delete penantiomer;
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(handedness);

}   // End of main program

/**************************************************************************
@author  David Belnap
@brief 	Gets criteria for selecting micrographs (from a field-of-view) for a
	handedness-determination calculation.
@param 	&tag		"mgsel", "mg1", or "mg2".
@param 	&tagvalue	string with one of options listed above.
@param 	*mg_select	micrograph selection criteria, views 1 & 2.
@param 	*mg_index	indices for micrograph selection, views 1 & 2.
@param 	*mg_ang		angles for micrograph selection, views 1 & 2.
@return int 		0.

	User has four options: 
    1.  angle,a (tilt angle)
    2.  index,n (index within a series, beginning with 1)
    3.  near,n  (nth nearest-to-focus, beginning with 1)
    4.  far,n   (nth farthest-from-focus, beginning with 1)

**/
int   hand_get_options_mgselect(Bstring& tag, Bstring& tagvalue, int* mg_select, int* mg_index, double* mg_ang)
{
	Bstring*	strarr = tagvalue.split(",");
	
	if ( tag == "mgsel" ) {
		if ( tagvalue.contains("ang") ) {
			mg_select[0] = mg_select[1] = 4;
			mg_ang[0] = strarr->next->real() * M_PI/180.0;
			mg_ang[1] = strarr->next->next->real() * M_PI/180.0;
		} else {
			mg_index[0] = strarr->next->integer();
			mg_index[1] = strarr->next->next->integer();
			if ( tagvalue.contains("near") ) mg_select[0] = mg_select[1] = 0;
			else if ( tagvalue.contains("far") ) mg_select[0] = mg_select[1] = 1;
			else if ( tagvalue.contains("ind") ) mg_select[0] = mg_select[1] = 2;
		}
	} else if ( tag == "mg1" )  {
		if ( tagvalue.contains("ang") ) {
			mg_select[0] = 4;
			mg_ang[0] = strarr->next->real() * M_PI/180.0;
		} else {
			mg_index[0] = strarr->next->integer();
			if ( tagvalue.contains("near") ) mg_select[0] = 0;
			else if ( tagvalue.contains("far") ) mg_select[0] = 1;
			else if ( tagvalue.contains("ind") ) mg_select[0] = 2;
		}
	} else if ( tag == "mg2" )  {
		if ( tagvalue.contains("ang") ) {
			mg_select[1] = 4;
			mg_ang[1] = strarr->next->real() * M_PI/180.0;
		} else {
			mg_index[1] = strarr->next->integer();
			if ( tagvalue.contains("near") ) mg_select[1] = 0;
			else if ( tagvalue.contains("far") ) mg_select[1] = 1;
			else if ( tagvalue.contains("ind") ) mg_select[1] = 2;
		}
	}
	
	string_kill(strarr);
	
	return 0;
}

