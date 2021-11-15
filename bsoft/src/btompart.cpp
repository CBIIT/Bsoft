/**
@file	btompart.cpp
@brief	Program to extract single particle images from a tilt series and generate individual reconstructions
@author Bernard Heymann
@date	Created: 20120307
@date	Modified: 20160114
**/

#include "rwmg.h"
#include "mg_reconstruct.h"
#include "mg_tomo_rec.h"
#include "mg_tomo_resol.h"
#include "mg_tomography.h"
#include "mg_particle_select.h"
#include "symmetry.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: btompart [options] input.star [input.star]",
"-------------------------------------------------",
"Averages reconstructions or subvolumes of reconstructions based on",
"the selection in parameter files.",
" ",
"Actions:",
"-add                     Adjust micrograph coordinates to particle origins.",
"-zflip                   Invert the particle z coordinate before transfer.",
"-transfer mg             Calculate micrograph particle coordinates from micrograph (mg) or reconstruction (rec).",
"-reconstruct             Reconstruct individual particles from micrographs.",
"-removemarkers 14        Mask out markers with this radius (pixels) before reconstruction.",
"-nloo 34.7,0.3           Determine particle resolution up to the given limit and with the given cutoff.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from parameter file; a single value can be given).",
"-views                   Use particle views during transfer.",
"-boxsize 55              Set the particle box size (default from parameter file).",
//"-resolution 4.5,130      Resolution range for correlation (default 0 - 1e6 angstrom).",
" ",
"Parameters for particle reconstruction:",
"-maxresolution 20        Reconstruction resolution limit (angstrom, default Nyquist).",
"-interpolation weighted  Interpolation type: nearest (default), weighted, trilinear.",
"-pad 3                   Image padding factor (default 2).",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
"-symmetry C5             Point group symmetry.",
"-base subvol             Set the particle base file name (no extension).",
"-path dir/subdir         Set the particle file paths.",
"-extension pif           Set the particle image file format.",
" ",
"Parameters for particle resolution estimation:",
"-ratio 2.5               Radial sampling to Cartesian sampling ratio (default 1).",
"-fast 10                 Fast mode, only use micrographs within the given angle from the selected one (default all).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-Postscript file.ps      Output plot file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				add_coordinates(1);			// Add coordinates and origins before other operations
	int				transfer(0);				// Flag to transfer particle coordinates
	int				zflip(0);					// Flag to invert z coordinate before transfer
	int				use_view(0);				// Flag to use original particle views during transfer
	int				part_recon(0);				// Flag to reconstruct individual particles
	double			marker_radius(0);			// Radius to mask out markers
	double			nloo(0);					// Flag and resolution for NLOO resolution estimation
	double			sampling_ratio(1);			// Averaging ratio for FRC calculation
	double			fast_angle(M_PI);			// Select all micrographs to do a reconstruction
	double			cutoff(0.3);				// NLOO cutoff
//	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	int				box_size(0);				// Particle box size
	double			hires(0), lores(1e6);		// Limiting resolution range (hires must be > 0 to be set)
	double 			resolution(0); 				// Must be set > 0 to limit resolution
	int				interp_type(0);				// Interpolation type
	int				pad_factor(2);				// Image padding factor
	int				ctf_action(0);				// Default no CTF operation
	double			wiener(0.2);				// Wiener for CTF correction
	Bstring			symmetry_string("C1");		// Point group
	Bstring			partbase;					// Particle base file name
	Bstring			partpath;					// Particle file path
	Bstring			partext;					// Particle image file format
	Bstring			outfile;					// Output parameter file name
	Bstring			psfile;						// Output Postscript file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "add" ) add_coordinates = 1;
		if ( curropt->tag == "zflip" ) zflip = 1;
		if ( curropt->tag == "transfer" ) {
			if ( curropt->value[0] == 'm' ) transfer = 1;
			if ( curropt->value[0] == 'r' ) transfer = 2;
		}
		if ( curropt->tag == "reconstruct" ) part_recon = 1;
		if ( curropt->tag == "removemarkers" )
			if ( ( marker_radius = curropt->value.real() ) < 1 )
				cerr << "-removemarkers: A marker radius must be specified!" << endl;
		if ( curropt->tag == "nloo" )
			if ( curropt->values(nloo, cutoff) < 1 )
				cerr << "-nloo: A resolution limit must be specified!" << endl;
//		if ( curropt->tag == "datatype" )
//			nudatatype = curropt->datatype();
		if ( curropt->tag == "views" ) use_view = 2;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "boxsize" )
        	if ( ( box_size = curropt->value.integer() ) < 1 )
				cerr << "-boxsize: A particle box size in pixels must be specified" << endl;
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
		if ( curropt->tag == "maxresolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-maxresolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "interpolation" ) {
			if ( curropt->value[0] == 'n' ) interp_type = 0;
			if ( curropt->value[0] == 'w' ) interp_type = 1;
			if ( curropt->value[0] == 't' ) interp_type = 2;
		}
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: An integer factor must be specified!" << endl;
		if ( curropt->tag == "CTF" )
			ctf_action = curropt->ctf_action();
		if ( curropt->tag == "wiener" ) {
			if ( ( wiener = curropt->value.real() ) < 0.000001 )
				cerr << "-wiener: A Wiener factor must be specified!" << endl;
			else {
				if ( wiener < 0.01 ) wiener = 0.01;
//				if ( wiener > 1 ) wiener = 1;
			}
		}
		if ( curropt->tag == "ratio" )
			if ( ( sampling_ratio = curropt->value.real() ) < 1 )
				cerr << "-ratio: A ratio must be specified!" << endl;
		if ( curropt->tag == "fast" ) {
			if ( ( fast_angle = curropt->value.real() ) < 1 )
				cerr << "-fast: An angle must be specified!" << endl;
			else
				fast_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "base" ) {
			partbase = curropt->value;
			if ( partbase.length() < 1 )
				cerr << "-base: The particle base file name must be specified!" << endl;
		}
		if ( curropt->tag == "path" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-path: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "extension" ) {
			partext = curropt->value;
			if ( partext.length() < 1 )
				cerr << "-extension: The particle file extension must be specified!" << endl;
		}
		if ( curropt->tag == "output" )
        	outfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
        	psfile = curropt->filename();
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

	Bproject*			project = read_project(file_list);
	string_kill(file_list);

	if ( !project )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( add_coordinates ) {
		project->select = transfer - 1;
		project_add_origins_to_coords(project);
	}
	
	if ( sam[0] ) {
		project_set_mg_pixel_size(project, sam);
		project_set_rec_voxel_size(project, sam);
	}
	
	if ( box_size ) {
		project->select = 0;
		project_set_particle_box_size(project, box_size);
		project->select = 1;
		project_set_particle_box_size(project, box_size);
	}
	
	if ( transfer ) {
		project->select = transfer - 1;
		if ( !part_series_from_seed(project, zflip | use_view) ) {
			cerr << "Error: No particles transferred!" << endl;
			bexit(-1);
		}
	}

	Bsymmetry			sym(symmetry_string);
	
	if ( part_recon )
		project_tomo_reconstruct_particles(project, resolution,
			interp_type, pad_factor, ctf_action, wiener, sym, 
			partbase, partpath, partext);

	Bplot*				plot = NULL;
	if ( nloo ) {
		plot = project_tomo_particle_resolution(project, nloo, sampling_ratio, 
				fast_angle, cutoff);
		if ( psfile.length() ) ps_plot(psfile, plot);
		delete plot;
	}
	
	if ( outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

