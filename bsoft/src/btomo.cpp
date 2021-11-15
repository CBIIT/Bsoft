/**
@file	btomo.cpp
@brief	Program to process tomographic series of images
@author	Bernard Heymann
@date	Created: 20020416
@date	Modified: 20210311
**/

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "mg_align.h"
#include "rwmg.h"
#include "rwmgIMOD.h"
#include "matrix_linear.h"
#include "rwimg.h"
#include "qsort_functions.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototype
int			read_imod_tlt(Bstring& filename, Bproject* project);

// Usage assistance
const char* use[] = {
" ",
"Usage: btomo [options] input.star/input.img [input.star/input.img]",
"------------------------------------------------------------------",
"Manipulates parameter files for tilt series and aligns them by cross-correlation.",
" ",
"Actions:",
"-all                     Select all markers.",
"-updatematrix            Update micrograph matrices from tilt and axis angles.",
"-resetmodel              Reset the x and y coordinates to that of the zero-degree tilt image.",
"-checkmarkers            Determine if markers all fall within micrograph boundaries.",
"-deselect 5-9,33         Deselect markers by id from all micrographs and models.",
"-deletemarkers 5-9,33    Delete markers by id from all micrographs and models.",
"-rotationaxis            Determine the best rotation axis from refined markers.",
"-clear 240,7.5           Clear extraneous areas on micrographs with this intended reconstruction thickness",
"                         and smoothing edge width (default 5).",
"-removemarkers 14        Mask out markers with this radius (pixels) before transformation.",
"-thickness 234           Calculate mean free path from a given thickness (in angstrom).",
"-emfp 2650               Calculate thickness from the effective MFP (in angstrom).",
"-catenate                Concatenate separate micrograph image files into one file.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for tilt axis (default 0,0,0).",
"-axis 74.7               Tilt axis angle relative to x-axis (default 0 or from parameter file).",
"-tilt -65,1.5            Tilt series start and step size (default from parameter file).",
"-reference 48            Reference image for alignment (default closest to zero tilt).",
"-goldradius 12           Gold particle radius (default from parameter file).",
" ",
"Input:",
"-rawtlt series.rawtlt    Text file with tilt angles.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-image file.pif          Output multi-image file with aligned images.",
"-Postscript inten.ps     Postscript file with fit of tilt series intensities.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int				select_all(0);			// Flag to select all markers
	int				update_matrix(0);		// Flag to update micrograph matrix
	int				reset_model(0);			// Flag to reset the x and y coordinates of the model
	int				checkmarkers(0);		// Flag to determine if markers all fall within micrograph boundaries
	Bstring			deselect;				// String specifying markers to deselect
	Bstring			deletemarkers;			// String specifying markers to delete
	int				rotation_axis(0);		// Flag to determine the best rotation axis
//	int				crosscorrelate(0);		// Flag to align by cross-correlation
	double			thickness(0);			// Calculates the mean free path
	double			emfp(0);				// Calculates the thickness
	int				catenate(0);			// Falg to concatenate micrographs into one file
	int				clearareas(0);			// Thickness for clearing extraneous areas
	double			edge_width(5);			// Smoothing edge width
	double			marker_erase_radius(0);	// Radius to mask out markers
	double			marker_radius(0);		// Gold marker radius
	DataType 		datatype(Unknown_Type);	// Conversion to new type
	int				ref_img(-1);			// Reference image for alignment, <0 means find zero tilt
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Vector3<double>	origin;					// Tilt axis origin
	double			tilt_axis(0);			// Tilt axis angle
	double			tilt_start(0);			// Tilt angle start
	double			tilt_step(0);			// Tilt_angle step size
//	int 			fill_type(FILL_USER);	// Type of fill value
//	double			fill(0);				// Fill value for new areas
	int				read_flags(32);			// Flags to pass to the parameter file reading function
    Bstring			tltfile;				// Input tilt angle file
    Bstring			maskfile;				// Mask to use for cross-correlation
	Bstring			paramfile;				// Output parameter file
	Bstring			imgfile;				// Output multi-image file
	Bstring			psfile;					// Output Postscript intensities plot file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) select_all = 1;
		if ( curropt->tag == "resetmodel" ) reset_model = 1;
		if ( curropt->tag == "updatematrix" ) update_matrix = 1;
		if ( curropt->tag == "checkmarkers" ) checkmarkers = 1;
		if ( curropt->tag == "deselect" ) deselect = curropt->value;
		if ( curropt->tag == "deletemarkers" ) deletemarkers = curropt->value;
		if ( curropt->tag == "rotationaxis" ) rotation_axis = 1;
//		if ( curropt->tag == "crosscorrelate" ) crosscorrelate = 1;
		if ( curropt->tag == "thickness" )
			if ( ( thickness = curropt->value.real() ) < 1 )
				cerr << "-thickness: A thickness must be specified!" << endl;
		if ( curropt->tag == "emfp" )
			if ( ( emfp = curropt->value.real() ) < 1 )
				cerr << "-emfp: A mean free path must be specified!" << endl;
		if ( curropt->tag == "catenate" ) catenate = 1;
		if ( curropt->tag == "clear" )
    	    if ( curropt->values(clearareas, edge_width) < 1 )
				cerr << "-clear: A thickness must be specified." << endl;
		if ( curropt->tag == "removemarkers" )
			if ( ( marker_erase_radius = curropt->value.real() ) < 1 )
				cerr << "-removemarkers: A marker radius must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			datatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "reference" )
    	    if ( ( ref_img = curropt->value.integer() ) < 0 )
				cerr << "-reference: An image number must be specified." << endl;
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "axis" ) {
    	   	tilt_axis = angle_set_negPI_to_PI(curropt->value.real()*M_PI/180.0);
			update_matrix = 1;
		}
		if ( curropt->tag == "tilt" ) {
    	    if ( curropt->values(tilt_start, tilt_step) < 2 )
				cerr << "-tilt: Two angles must be specified." << endl;
			else {
				tilt_start = angle_set_negPI_to_PI(tilt_start*M_PI/180.0);
				tilt_step = angle_set_negPI_to_PI(tilt_step*M_PI/180.0);
				update_matrix = 1;
			}
        }
		if ( curropt->tag == "goldradius" )
    	    if ( ( marker_radius = curropt->value.real() ) < 1 )
				cerr << "-goldradius: A gold marker radius must be specified." << endl;
//		if ( curropt->tag == "fill" )
//			fill = curropt->fill(fill_type);
		if ( curropt->tag == "rawtlt" )
			tltfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();

	Bproject*		project = NULL;
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	if ( file_type(*file_list) == Micrograph ) {
		project = read_project(file_list, read_flags);
	} else {
		project = project_create_from_images(file_list, "mg");
		project_calculate_angles(project);
	}
		
	string_kill(file_list);

	if ( !project ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( verbose & VERB_PROCESS )
		project_show_markers(project);

	if ( catenate )
		project_catenate_micrographs(project);

	if ( origin[0] )
		project_set_micrograph_origins(project, origin);
	
	if ( project->field && project->field->mg->origin.length() < 1 )
		project_set_nominal_mg_origins(project);
	
	if ( tilt_axis )
		project_set_tilt_axis(project, tilt_axis);

	if ( tltfile.length() )
		read_imod_tlt(tltfile, project);
	else if ( tilt_step )
		project_set_tilt_angles(project, tilt_start, tilt_step);
	
	if ( project->field && project->field->mg->matrix.determinant() < 0.5 )
		update_matrix = 1;

	if ( update_matrix )
		project_mg_tilt_to_matrix(project);
	
	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	if ( clearareas )
		project_clear_extraneous_areas(project, clearareas, edge_width);

	double			tlr(0);
	Bstring			txt;
	Bplot*			plot = NULL;
	if ( thickness || emfp ) {
		plot = project_intensity_plot(project);
		tlr = project_fit_intensities(project, plot);
		if ( thickness ) emfp = thickness/tlr;
		if ( emfp ) {
			thickness = emfp*tlr;
			txt = Bstring(emfp, "Effective MFP: %lg A");
			plot->page(0).add_text(txt);
			txt = Bstring(thickness, "Thickness: %lg A");
			plot->page(0).add_text(txt);
		}
		if ( thickness && emfp ) {
			if ( verbose ) {
				cout << "Thickness:                      " << thickness << " A" << endl;
				cout << "Effective mean free path:       " << emfp << " A" << endl << endl;
			}
		}
		if ( plot && psfile.length() )
			ps_plot(psfile, plot);
	}
	
	int				i;
	Bmicrograph*	mg_ref = NULL;
	
	if ( marker_radius )
		project_set_marker_radius(project, marker_radius);
	
	if ( ref_img > -1 && project->field )
		for ( i=0, mg_ref=project->field->mg; mg_ref && ( i<ref_img ); mg_ref=mg_ref->next, i++ ) ;
	
	if ( !mg_ref )
		mg_ref = field_find_zero_tilt_mg(project->field);

	if ( checkmarkers )
		project_check_markers(project, 1);

	if ( select_all )
		project_mg_marker_select(project, -1);
	
	if ( deselect.length() )
		project_deselect_markers(project, deselect);

	if ( deletemarkers.length() )
		project_delete_markers(project, deletemarkers);

	if ( reset_model && project->rec && project->rec->mark )
		mg_reset_model(mg_ref, project->rec->mark);

	if ( rotation_axis )
		project_marker_rotation_axis(project);
	
	if ( project->rec && project->rec->mark )
		project_tomo_residuals(project, 1);
	
	if ( marker_erase_radius )
		project_erase_markers(project, marker_erase_radius);
	
	if ( imgfile.length() )
		project_write_aligned_images(project, NULL, imgfile, datatype);

    if ( paramfile.length() )
		write_project(paramfile, project, 0, 0);
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

