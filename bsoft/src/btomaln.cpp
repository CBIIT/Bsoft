/**
@file	btomaln.cpp
@brief	Program to do fiducialless alignment of a tilt series
@author	Bernard Heymann
@date	Created: 20040407
@date	Modified: 20210706
**/

#include "mg_processing.h"
#include "mg_align.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: btomaln [options] input.star [input.star]",
"------------------------------------------------",
"Aligns a tilt series of micrographs by progressive reconstruction.",
" ",
"Actions:",
"-align 5,1.5,8           Align micrographs: number of iterations, stopping change,",
"                         and number of adjacent images for limited reconstructions",
"                         (defaults 1 and ± 5 images).",
"-emfp 2650               Calculate thickness from the effective MFP (in angstrom).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-axis 74.7               Tilt axis angle relative to x-axis (default 0 or from parameter file).",
"-tilt -65,1.5            Tilt series start and step size (default from parameter file).",
"-thickness 235           Reconstruction thickness (default 20% of x size).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 900,300      High and low resolution limits for cross-correlation (default 0.1,1000 angstrom).",
"-shiftlimit 155.5        Limit on micrograph shift search in pixels (default 25% of micrograph width).",
"-origin 110,50,44        Origin for tilt axis (default 0,0,0).",
"-edge 23,12              Smooth the edge to a given width, with gaussian decay of a given width.",
"-fill 127                Fill value for particles extending beyond image (default background).",
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
	DataType 		datatype(Unknown_Type);	// Conversion to new type
	long			iter(0);				// Number of alignment iterations
	double			dchange(1);				// Threshold for change in origins
	long			dimg(5);				// Number of adjacent images in reconstructions
	Vector3<double>	sampling;				// Units for the three axes (A/pixel)
	double			tilt_axis(0);			// Tilt axis angle
	double			tilt_start(0);			// Tilt angle start
	double			tilt_step(0);			// Tilt_angle step size
	long			thickness(0);			// Reconstruction thickness
	double			emfp(0);				// Calculates the thickness
	double			hi_res(0), lo_res(1e10);	// Default resolution range
	double			shift_limit(-1);		// Micrograph shift search limit
	Vector3<double>	origin;					// Tilt axis origin
	double			edge_width(0), gauss_width(0);	// Edge parameters
	int 			fill_type(FILL_BACKGROUND);
	double			fill(0);
	int				read_flags(32);			// Flags to pass to the parameter file reading function
	Bstring			paramfile;				// Output parameter file
	Bstring			imgfile;				// Output multi-image file
	Bstring			psfile;					// Output Postscript intensities plot file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "align" )
    	    if ( curropt->values(iter, dchange, dimg) < 1 )
				cerr << "-align: A number of iterations must be specified." << endl;
		if ( curropt->tag == "emfp" )
			if ( ( emfp = curropt->value.real() ) < 1 )
				cerr << "-emfp: A mean free path must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			datatype = curropt->datatype();
		if ( curropt->tag == "axis" ) {
			tilt_axis = curropt->value.real();
			tilt_axis = angle_set_negPI_to_PI(tilt_axis*M_PI/180.0);
		}
		if ( curropt->tag == "tilt" ) {
    	    if ( curropt->values(tilt_start, tilt_step) < 2 )
				cerr << "-tilt: Three angles must be specified." << endl;
			else {
				tilt_start = angle_set_negPI_to_PI(tilt_start*M_PI/180.0);
				tilt_step = angle_set_negPI_to_PI(tilt_step*M_PI/180.0);
			}
        }
		if ( curropt->tag == "thickness" ) thickness = curropt->value.integer();
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( hi_res > lo_res )
				swap(hi_res, lo_res);
        }
		if ( curropt->tag == "shiftlimit" )
    	    if ( ( shift_limit = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A shift limit in pixels must be specified." << endl;
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "edge" )
    	    if ( curropt->values(edge_width, gauss_width) < 1 )
				cerr << "-edge: An edge width must be specified." << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
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

	if ( !project ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( sampling[0] > 0 )
		project_set_mg_pixel_size(project, sampling);

	if ( tilt_axis )
		project_set_tilt_axis(project, tilt_axis);

	if ( tilt_step )
		project_set_tilt_angles(project, tilt_start, tilt_step);
	
	double			tlr(0);
	Bstring			txt;
	Bplot*			plot = NULL;
	if ( thickness || emfp ) {
		plot = project_intensity_plot(project);
		tlr = project_fit_intensities(project, plot);
		if ( thickness ) emfp = thickness/tlr;
		if ( emfp ) {
			thickness = emfp*tlr;
			txt = Bstring(emfp, "Effective MFP: %lg");
			plot->page(0).add_text(txt);
			txt = Bstring(thickness, "Thickness: %lg");
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
	
		
	if ( iter > 0 ) {
		thickness *= 1.5;
		project_tomo_align(project, thickness, iter, dchange, dimg, hi_res, shift_limit, edge_width, gauss_width);
	}
	
	if ( imgfile.length() )
		project_write_aligned_images(project, NULL, imgfile, datatype);
	
    if ( paramfile.length() )
		write_project(paramfile, project, 0, 0);
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}


