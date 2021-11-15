/**
@file	btrack.cpp
@brief	Program to track fiducial markers in a tomographic series of images
@author  Bernard Heymann
@date	Created: 20020416
@date	Modified: 20210311
**/

#include "rwmg.h"
#include "rwimg.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_align.h"
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
"Usage: btrack [options] input.star [input.star]",
"-----------------------------------------------",
"Track fiducial markers in a tomographic tilt series of micrographs.",
" ",
"Actions:",
"-all                     Select all markers.",
"-updatematrix            Update micrograph matrices from tilt and axis angles.",
"-emfp 2650               Calculate thickness from the effective MFP (in angstrom).",
"-adjusttilt              Apply tilt adjustments based on intensity variation.",
"-pick 1,20,1             Pick seed markers in reference image (1) or all images (2),",
"                         excluding the distance from the edge (in pixels).",
"                         and a flag to add to rather than replace existing markers.",
"-generatemarkers         Generates markers for all micrographs from a seed.",
"-findaxis 10,1.5,-30,50  Finds the tilt axis using seed markers and micrographs ",
"                         near the given tilt angle and with the axis angle step size,",
"                         start and end (default axis parameters: 1,-180,180).",
"-track 5,1.5             Track markers: number of iterations and residual target (default 10% of marker radius).",
"                         using zero-tilt reference markers as seed.",
"-crossmarker             Cross-correlate markers during z-search (default real space correlate).",
"-recenter                Recenter z coordinates during tracking (default not).",
"-refine shifts           Refine marker parameters:",
"                         markers: Refine marker micrograph positions.",
"                         z:       Refine marker z positions.",
"                         origins: Refine micrograph origins.",
"                         scales:  Refine micrograph scales.",
"                         views:   Refine view transforms.",
"                         all:     Refine origins, scales and views simultaneously.",
"                         8,v,o,z: Refine iteratively in sequence views, origins and z's.",
"-resetmodel              Reset the x and y marker coordinates to that of the zero-degree tilt image.",
"-angles                  Calculate micrograph orientation angles.",
"-check                   Determine if markers all fall within micrograph boundaries.",
"-fix                     Fix markers falling outside micrograph boundaries.",
"-errors 1.6              Show errors above this threshold (distance in pixels).",
"-exclude 0,1,134-140     Exclude these micrographs (can be \"none\").",
"-deselect 5-9,33         Deselect markers by id from all micrographs and models.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 900,300      High and low resolution limits for cross-correlation (default 0.1,1000 angstrom).",
"-origin 110,50,44        Origin for tilt axis (default 0,0,0).",
"-axis 74.7               Tilt axis angle relative to x-axis (default 0 or from parameter file).",
"-tilt -65,1.5            Tilt series start and step size (default from parameter file).",
"-reference 48            Reference image for alignment (default closest to zero tilt).",
"-edge 23,12              Smooth the edge to a given width, with gaussian decay of a given width.",
"-goldradius 12           Gold particle radius (default from parameter file).",
"-shiftlimit 155.5        Limit on micrograph shift search in pixels (default 25% of micrograph width).",
"-thickness 2156.5        Estimated tomogram thickness (default and minimum 1000 angstrom).",
" ",
"Output:",
"-output file.star        Output parameter file (produces numbered files for tracking).",
"-image file.pif          Output multi-image file with aligned images.",
"-Postscript out.ps       Postscript output with marker error plots.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int				select_all(0);			// Flag to select all markers
//	int				update_matrix(0);		// Flag to update micrograph matrix
	int				adjust_tilt(0);			// Flag to apply tilt adjustments based on intensity variation
	int				reset_model(0);			// Flag to reset the x and y coordinates of the model
	int				genmarkers(0);			// Generates markers for all micrographs from a seed
	int				pickmarkers(0);			// Flag to pick markers in reference image (1) or all images (2)
	int				addmarkers(0);			// Flag to add picked markers rather than replace old markers
	double			emfp(0);				// Calculates the thickness
	double			findaxis(0);			// Flag and tilt angle to find tilt axis
	double			axis_step(1);			// Tilt axis angle step size (in degrees)
	double			axis_start(-180);		// Tilt axis angle start (in degrees)
	double			axis_end(180);			// Tilt axis angle end (in degrees)
	int				trackmarkers(0);		// Flag to align using markers
	double			target(0);				// Target resdual to terminate tracking
	int				cc_type(0);				// Type of marker correlation: 0=real space, 1=cross correlation
	int				recenter(0);			// Flag to recenter z coordinates
	int				edge(0);				// Edge to exclude from picking markers
	int				ref_iter(0);			// Refine geometric transformations: iterations
	double			ref_tol(1e-4);			// Refine geometric transformations: tolerance
	Bstring			refop;					// String to define refine operations
	int				fixmarkers(0);			// Flag to fix markers falling outside micrograph boundaries
	int				checkmarkers(0);		// Flag to determine if markers all fall within micrograph boundaries
	int				calc_angles(0);			// Flag to calculate micrograph orientation angles
	double			error_cutoff(0);		// Threshold to show errors
	Bstring			mg_exclude;				// Micrographs to exclude
	Bstring			deselect;				// String specifying markers to deselect
	double			gold_radius(0);			// Gold marker radius
	DataType 		datatype(Unknown_Type);	// Conversion to new type
	int				ref_img(-1);			// Reference image for alignment, <0 means find zero tilt
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	double			hi_res(0),lo_res(1000);	// Default resolution range
	Vector3<double>	origin;					// Tilt axis origin
	double			tilt_axis(0);			// Tilt axis angle
	double			tilt_start(0);			// Tilt angle start
	double			tilt_step(0);			// Tilt_angle step size
	double			edge_width(0), gauss_width(0);	// Edge parameters
	double			shift_limit(-1);		// Micrograph shift search limit
	double			thickness(0);			// Estimated tomogram thickness in angstrom
	int				read_flags(32);			// Flags to pass to the parameter file reading function
	Bstring			paramfile;				// Output parameter file
	Bstring			imgfile;				// Output multi-image file
	Bstring			psfile;					// Output postscript file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) select_all = 1;
//		if ( curropt->tag == "updatematrix" ) update_matrix = 1;
		if ( curropt->tag == "adjusttilt" ) adjust_tilt = 1;
		if ( curropt->tag == "emfp" )
			if ( ( emfp = curropt->value.real() ) < 1 )
				cerr << "-emfp: A mean free path must be specified!" << endl;
		if ( curropt->tag == "pick" )
    	    if ( curropt->values(pickmarkers, edge, addmarkers) < 1 )
				cerr << "-pick: A 1 or 2 must be specified." << endl;
		if ( curropt->tag == "generatemarkers" ) genmarkers = 1;
		if ( curropt->tag == "findaxis" ) {
    	    if ( curropt->values(findaxis,
					axis_step, axis_start, axis_end) < 1 )
				cerr << "-findaxis: A tilt angle must be specified." << endl;
			else {
				findaxis *= M_PI/180.0;
				axis_step *= M_PI/180.0;
				axis_start *= M_PI/180.0;
				axis_end *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "track" )
    	    if ( curropt->values(trackmarkers, target) < 1 )
				cerr << "-track: A maximum number of iterations must be specified." << endl;
		if ( curropt->tag == "crossmarker" ) cc_type = 1;
		if ( curropt->tag == "recenter" ) recenter = 1;
		if ( curropt->tag == "refine" ) {
			ref_iter = curropt->value.integer();
			if ( ref_iter < 1 ) {
				ref_iter = 1;
				refop = curropt->value.substr(0,1).lower();
			} else {
				refop = curropt->value.post(',').remove(',');
 			}
			if ( refop.length() < 1 ) {
				cerr << "-refine: The option is misformed! (" << curropt->value << ")" << endl;
				bexit(-1);
			}
		}
		if ( curropt->tag == "resetmodel" ) reset_model = 1;
		if ( curropt->tag == "angles" ) calc_angles = 1;
		if ( curropt->tag == "check" ) checkmarkers = 1;
		if ( curropt->tag == "fix" ) fixmarkers = 1;
		if ( curropt->tag == "errors" )
    	    if ( ( error_cutoff = curropt->value.real() ) < 0.001 )
				cerr << "-errors: An error cutoff value must be specified." << endl;
		if ( curropt->tag == "exclude" ) mg_exclude = curropt->value;
		if ( curropt->tag == "deselect" ) deselect = curropt->value;
		if ( curropt->tag == "datatype" )
			datatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "reference" )
    	    if ( ( ref_img = curropt->value.integer() ) < 0 )
				cerr << "-reference: An image number must be specified." << endl;
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( hi_res > lo_res )
				swap(hi_res, lo_res);
        }
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "axis" ) {
			tilt_axis = curropt->value.real();
			tilt_axis = angle_set_negPI_to_PI(tilt_axis*M_PI/180.0);
//			update_matrix = 1;
 /*   	    if ( ( tilt_axis = curropt->value.real() ) == 0 )
				cerr << "-axis: The tilt axis angle must be specified." << endl;
			else {
				tilt_axis = angle_set_negPI_to_PI(tilt_axis*M_PI/180.0);
				update_matrix = 1;
			}*/
		}
		if ( curropt->tag == "tilt" ) {
    	    if ( curropt->values(tilt_start, tilt_step) < 2 )
				cerr << "-tilt: Three angles must be specified." << endl;
			else {
				tilt_start = angle_set_negPI_to_PI(tilt_start*M_PI/180.0);
				tilt_step = angle_set_negPI_to_PI(tilt_step*M_PI/180.0);
//				update_matrix = 1;
			}
        }
		if ( curropt->tag == "edge" )
    	    if ( curropt->values(edge_width, gauss_width) < 1 )
				cerr << "-edge: An edge width must be specified." << endl;
		if ( curropt->tag == "goldradius" )
    	    if ( ( gold_radius = curropt->value.real() ) < 1 )
				cerr << "-goldradius: A gold marker radius must be specified." << endl;
		if ( curropt->tag == "shiftlimit" )
    	    if ( ( shift_limit = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A shift limit in pixels must be specified." << endl;
		if ( curropt->tag == "thickness" )
    	    if ( ( thickness = curropt->value.real() ) < 1 )
				cerr << "-thickness: An estimated tomogram tickness must be specified." << endl;
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}
	
	Bproject*		project = read_project(file_list, read_flags);		
	string_kill(file_list);
	
	if ( !project || !project->field || !project->field->mg ) {
		cerr << "Error: No micrographs in input file!" << endl;
		bexit(-1);
	}
	
	if ( verbose & VERB_PROCESS )
		project_show_markers(project);

	project_sort_markers_by_id(project);
	
	if ( origin[0] )
		project_set_micrograph_origins(project, origin);

	if ( project->field->mg->origin[0] < 1 || project->field->mg->origin[1] < 1 )
		project_set_nominal_mg_origins(project);

	if ( mg_exclude.length() )
		project_mg_exclude(project, mg_exclude);

	if ( deselect.length() )
		project_deselect_markers(project, deselect);

	if ( tilt_axis )
		project_set_tilt_axis(project, tilt_axis);

	if ( tilt_step )
		project_set_tilt_angles(project, tilt_start, tilt_step);
	
//	if ( project->field->mg->matrix.determinant() < 0.5 )
//		update_matrix = 1;

//	if ( update_matrix )
//		project_mg_tilt_to_matrix(project);
	
	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	double			tlr(0);
	Bplot*			plot = NULL;
	if ( adjust_tilt || emfp ) {
		plot = project_intensity_plot(project);
		tlr = project_fit_intensities(project, plot, adjust_tilt);
		if ( emfp ) {
			if ( thickness < 100 ) thickness = emfp*tlr;
			if ( verbose ) {
				cout << "Thickness:                      " << emfp*tlr << " A" << endl;
				cout << "Effective mean free path:       " << emfp << " A" << endl << endl;
			}
		}
	}
	
	int				i;
	Bmicrograph*	mg_ref = NULL;
	if ( ref_img > -1 )
		for ( i=0, mg_ref=project->field->mg; mg_ref && ( i<ref_img ); mg_ref=mg_ref->next, i++ ) ;
	
	if ( !mg_ref ) 
		mg_ref = field_find_zero_tilt_mg(project->field);
//		cout << "Reference micrograph:           " << mg_ref->id << " (" << mg_ref->img_num << ")" << endl;

	if ( !mg_ref->mark && !pickmarkers )
		for ( mg_ref = project->field->mg; mg_ref && !mg_ref->mark; mg_ref = mg_ref->next ) ;
//		cout << "Reference micrograph:           " << mg_ref->id << " (" << mg_ref->img_num << ")" << endl;
	
	if ( !mg_ref ) {
		mg_ref = field_find_zero_tilt_mg(project->field);
		if ( !pickmarkers ) pickmarkers = 1;
	}
	
	if ( verbose )
		cout << "Reference micrograph:           " << mg_ref->id << " (" << mg_ref->img_num << ")" << endl;
	
	if ( pickmarkers == 1 )
		mg_find_markers(mg_ref, edge, addmarkers);
	else if ( pickmarkers == 2 )
		project_find_markers(project, edge, addmarkers);
	
	if ( !project->rec )
		reconstruction_add(&project->rec, project->field->id);

	// Set a uniform voxel size
	sam = mg_ref->pixel_size;
	sam[2] = sam[1] = sam[0];
	project_set_rec_voxel_size(project, sam);
	
	if ( project->rec->origin.length() < 1 )
		project->rec->origin = micrograph_get_nominal_origin(mg_ref);

	if ( gold_radius ) {
		project_set_marker_radius(project, gold_radius);
	} else {
		project->rec->mark_radius = mg_ref->mark_radius;
	}
	
	if ( select_all )
		project_mg_marker_select(project, -1);
	
	if ( findaxis )
		project_find_tilt_axis(project, findaxis, axis_start, axis_end, axis_step, hi_res, lo_res, shift_limit);
	
	if ( !project->rec->mark )
		project_calculate_model(project);
	
	if ( genmarkers )
		project_generate_markers(project);

	if ( reset_model && project->rec->mark )
		mg_reset_model(mg_ref, project->rec->mark);
	
	if ( checkmarkers )
		project_check_markers(project, 1);
	
	if ( fixmarkers && project->rec->mark )
		project_fix_markers(project);

	if ( trackmarkers > 0 )
//		project_track_markers(project, hi_res, lo_res, shift_limit,
//			thickness, trackmarkers, target, cc_type, recenter, paramfile);
		project_track_markers_dual(project, hi_res, lo_res, shift_limit,
			thickness, trackmarkers, target, cc_type, recenter, paramfile);
	
	if ( ref_iter && refop[0] == 'm' )
		project_refine_markers(project, hi_res, lo_res);
	
	if ( ref_iter && project->rec->mark )
		project_refine(project, ref_iter, ref_tol, refop);

	if ( calc_angles )
		project_calculate_angles(project);
	
	if ( checkmarkers || fixmarkers || ref_iter ) {
		project_tomo_errors(project);
		project_tomo_residuals(project, 1);
	}

	if ( error_cutoff ) 
		project_show_errors(project, error_cutoff);
	
	if ( psfile.length() ) {
		if ( verbose )
			cout << "Writing " << psfile << endl;
		ofstream*	fps = NULL;
		if ( plot ) {
			fps = ps_open_and_init(psfile, plot);
			ps_graph(fps, plot, 1);
		} else {
			fps = ps_open_and_init(psfile, psfile, 3, 600, 800);
		}
		ps_marker_plots(fps, psfile, project);
		ps_close(fps);
	}
	
    if ( paramfile.length() ) {
		if ( verbose )
			cout << "Writing " << paramfile << endl;
		write_project(paramfile, project, 0, 0);
	}
	
	if ( imgfile.length() ) {
		if ( verbose )
			cout << "Writing " << imgfile << endl;
		project_write_aligned_images(project, NULL, imgfile, datatype);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

