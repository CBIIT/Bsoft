/**
@file	bpick.cpp
@brief	Pick particles in a micrograph given the particle coordinates, or pick particles in a focal series of micrographs with alignment of the micrographs.
@author Bernard Heymann and Samuel Payne
@date	Created: 20000505
@date	Modified: 20200616 (BH)
**/

#include "mg_pick.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_extract.h"
#include "mg_particles.h"
#include "mg_particle_select.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpick [options] input.star [input2.star ...]",
"---------------------------------------------------",
"Picks and extracts single particles from 2D (micrograph) or 3D (tomogram) images.",
" ",
"Actions:",
"-pick cc                 Picking type: cc, var, ring, background (default none).",
"-orient ref.pif          Orient using a reference map.",
"                         High resolution and symmetry determine number of projections.",
"-center                  Center each particle image (specify resolution limits).",
"-axis 5,13.6             Pick on a symmetry axis and distance along the axis (voxels, specify -symmetry).",
"-extract 0.5             Extract particles into new images at the given scale.",
"-split                   Split particles into individual files on extraction.",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-project                 Project 3D particles to 2D and assign to zero-tilt micrograph.",
"-add                     Adjust micrograph coordinates to particle origins.",
"-filterextremes          Filter extremes before picking or centering.",
"-background              Correct background of extracted particles.",
"-normalize               Normalize extracted particles.",
"-delete                  Delete all non-selected particles.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-box 123,145,130         Particle box size (voxels, one sets all).",
"-resolution 23.5,200     Resolution limits for picking and centering (default 0 - 1e6 angstrom).",
"-fill 127                Fill value for particles extending beyond image (default background).",
"-mgpath dir/subdir       Set the micrograph file paths.",
"-pspath dir/subdir       Set the power spectrum file paths.",
"-partpath dir/subdir     Set the particle file paths.",
"-base subvol             Set the particle base file name (no extension).",
"-extension pif           Set the particle image file format.",
" ",
"Parameters for picking particles:",
"-average 11              Averaging/smoothing kernel size.",
"-exclusion 210           Minimum distance between particles (default 0.8*x-size).",
"-bin 2                   Binning to apply (default 1).",
"-fom 0.1,0.8             Selecting picked particles based on an FOM range (default 0,1e30).",
" ",
"Parameters for picking particles in variance map:",
"-variance 51             Calculate a local variance image using the given size kernel.",
"-sigma 2.5               Sigma units above average to select particles (default 5).",
"-maximum 241             Maximum number of background areas to pick (default 1000000).",
" ",
"Parameters for picking particles by ring pattern:",
"-diameters 115,130       Particle and background ring diameters.",
"-contrast black          Foreground (particles) are white or black (default).",
"-angle 2.5               Angular increment (default 1 degree).",
" ",
"Parameter for orienting and selecting a symmetry axis:",
"-symmetry C5             Point group symmetry (used with -orient and -axis).",
" ",
"Input:",
"-Template file.map       Input template file for picking by cross-correlation.",
"-Mask file.map           Input frequency space mask file.",
"-compare file.star       Input parameter file for comparison.",
" ",
"Output:",
"-output file.star        Output parameter file.",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize optional variables
	int				pick_type(0);			// 1=cc, 2=var, 3=ring, 4=bkg
	int				use_rec(0);				// Flag to process reconstructions
	int				center(0);				// Flag to center particles
	int				sym_axis(0);			// Symmetry axis
	int				delete_deselected(0);	// Flag to delete deselected particles
	double			ext_scale(0);			// Flag and scale to extract particles
	double			axis_dist(0);			// Distance along icosahedral axis
	Vector3<double> sampling;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;					// Origin
	int				set_origin(0);			// Flag to set origin
	int 			add_coordinates(0);		// For adding mg coordinates and particle origins
	int				flags(0);				// Bit 1=filter
	Bsymmetry		sym;				// Point group
	double			hires(0), lores(0);		// Resolution limits for centering
	double			excl_dist(0);			// Minimum distance between particles picked
	long			bin(1);					// Binning to apply
	double			fommin(0), fommax(1e30);	// FOM range for particle selection
	long			average_kernel(0);		// Average filter kernel size
	long 			var_kernel(0);			// Variance filter kernel size
	double			sigma(5);				// Sigma units above average to select particles
	long			maxnum(1000000);		// Maximum number of background areas to pick
	double			rp_din(0), rp_dout(0);	// Diameters for finding particles
	double			rp_anginc(M_PI/180);	// Angular increment for finding particles
	int				contrast(0);			// Foreground: white=1, black=0
	int				split(0);				// Flag to split output into individual files
	Vector3<long>	box_size;				// Box size for picking and extracting particles
	int				back_flag(0);			// Background correction for extracting particles
	int				norm_flag(0);			// Normalization for extracting particles
	int 			fill_type(FILL_BACKGROUND);
	double			fill(0);
	Bstring			mgpath;					// Micrograph file path
	Bstring			pspath;					// Power spectrum file path
	Bstring			partpath;				// Particle file path
	Bstring			partbase;				// Particle base file name
	Bstring			partext;				// Particle image file format
	Bstring			tempfile;				// Input template file name for picking
	Bstring			orientfile;				// Input template file name for orientation
	Bstring			maskfile;				// Input frequency space mask file name
	Bstring			compfile;				// Input comparison parameter file name
	Bstring			outfile;				// Output parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "pick" ) {
			if ( curropt->value[0] == 'c' ) pick_type = 1;
			if ( curropt->value[0] == 'v' ) pick_type = 2;
			if ( curropt->value[0] == 'r' ) pick_type = 3;
			if ( curropt->value[0] == 'b' ) pick_type = 4;
		}
		if ( curropt->tag == "axis" )
			if ( curropt->values(sym_axis, axis_dist) < 2 )
				cerr << "-axis: Both an axis and distance must be specified!" << endl;
		if ( curropt->tag == "orient" )
			orientfile = curropt->filename();
		if ( curropt->tag == "delete" )
			delete_deselected = 1;
		if ( curropt->tag == "extract" )
			if ( ( ext_scale = curropt->value.real() ) < 0.01 )
				cerr << "-extract: An extraction scale must be specified! (usually 1)" << endl;
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "project" ) use_rec = 2;
		if ( curropt->tag == "add" ) add_coordinates = 1;
		if ( curropt->tag == "diameters" )
			if ( curropt->values(rp_din, rp_dout) < 2 )
				cerr << "-diameters: Two diameters must be specified!" << endl;
		if ( curropt->tag == "contrast" ) if ( curropt->value[0] == 'w' ) contrast = 1;
		if ( curropt->tag == "filterextremes" ) flags |= 1;
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "background" ) back_flag = 1;
		if ( curropt->tag == "normalize" ) norm_flag = 1;
		if ( curropt->tag == "split" ) split = 1;
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "box" ) box_size = curropt->size();
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "exclusion" )
			if ( ( excl_dist = curropt->value.real() ) <= 0 )
				cerr << "-exclusion: An separation distance must be specified!" << endl;
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: An integer must be specified!" << endl;
		if ( curropt->tag == "fom" )
			if ( curropt->values(fommin, fommax) < 1 )
				cerr << "-fom: A FOM cutoff must be specified!" << endl;
		if ( curropt->tag == "average" )
			if ( ( average_kernel = curropt->integer() ) < 1 )
				cerr << "-average: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "variance" )
			if ( ( var_kernel = curropt->integer() ) < 1 )
				cerr << "-variance: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "sigma" )
			if ( ( sigma = curropt->value.real() ) <= 0 )
				cerr << "-sigma: A value must be specified!" << endl;
		if ( curropt->tag == "maximum" )
			if ( ( maxnum = curropt->integer() ) < 1 )
				cerr << "-maximum: The maximum number of background areas must be specified!" << endl;
		if ( curropt->tag == "angle" ) {
			if ( ( rp_anginc = curropt->value.real() ) < 0.01 )
				cerr << "-angle: An angular increment must be specified!" << endl;
			else
				rp_anginc *= M_PI/180.0;
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "mgpath" ) {
			mgpath = curropt->value;
			if ( mgpath.length() < 1 )
				cerr << "-mgpath: The micrograph file path must be specified!" << endl;
			else
				if ( mgpath[-1] != '/' ) mgpath += "/";
		}
		if ( curropt->tag == "pspath" ) {
			pspath = curropt->value;
			if ( pspath.length() < 1 )
				cerr << "-pspath: The power spectrum file path must be specified!" << endl;
			else
				if ( pspath[-1] != '/' ) pspath += "/";
		}
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "base" ) {
			partbase = curropt->value;
			if ( partbase.length() < 1 )
				cerr << "-base: The particle base file name must be specified!" << endl;
		}
		if ( curropt->tag == "extension" ) {
			partext = curropt->value;
			if ( partext.length() < 1 )
				cerr << "-extension: The particle file extension must be specified!" << endl;
		}
		if ( curropt->tag == "Template" )
			tempfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "compare" )
			compfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
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

	project->select = use_rec;
	
	Bproject*		projcomp = NULL;
	if ( compfile.length() ) {
		projcomp = read_project(compfile);
		if ( use_rec ) projcomp->select = 1;
	}

//	cout << "mgpath = " << mgpath << endl;
	if ( mgpath.length() )
			project_set_micrograph_path(project, mgpath);
	if ( pspath.length() )
			project_set_powerspectrum_path(project, pspath);
	if ( partpath.length() )
			project_set_particle_path(project, partpath);

	if ( sampling[0] > 0 )
		project_set_mg_pixel_size(project, sampling);

	if ( box_size[0] > 0 )
		project_set_particle_box_size(project, box_size);
		
	if ( add_coordinates )
		project_add_origins_to_coords(project);

	Bimage*			ptemp = NULL;
	if ( tempfile.length() ) {
		ptemp = read_img(tempfile, 1, 0);
		if ( set_origin ) {
			if ( set_origin == 2 ) ptemp->origin(ptemp->size()/2);
			else ptemp->origin(origin);
		} else
			origin = ptemp->image->origin();
	}
	
	if ( origin.length() < 1 ) origin = box_size/2;
		
	Bimage*			pmask = NULL;
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, 0);
	
	if ( ptemp && excl_dist < 1 ) excl_dist = 0.8*ptemp->sizeX();
	
	if ( pick_type == 1 ) {
		if ( ptemp )
			project_pick_particles(project, ptemp, pmask, hires, lores, fommin, fommax, excl_dist, bin);
		else
			cerr << "Error: Picking particles by cross-correlation requires a template!" << endl;
	} else if ( pick_type == 2 ) {
		project_pick_particles(project, average_kernel, var_kernel, sigma, origin[0], excl_dist, bin);
	} else if ( pick_type == 3 ) {
		project_pick_particles(project, rp_din, rp_dout, average_kernel, rp_anginc, flags, contrast);
	} else if ( pick_type == 4 ) {
		project_pick_background(project, maxnum, average_kernel, var_kernel, excl_dist);
	}

	project_set_part_links(project);

	if ( sym.point() > 101 && sym_axis > 1 )
		project_pick_sym_axis(project, sym, sym_axis, axis_dist);
		
	if ( center )
		project_find_part_centers_in_mgs(project, hires, lores, flags);

	if ( orientfile.length() ) {
		project_extract_orient_particles(project, orientfile, sym, hires, lores, bin);
		project_delesect_edge_particles(project);
	}

	if ( fommin > 0 )
		part_deselect(project, 0, fommin, fommax);
	
	if ( excl_dist ) {
		part_deselect_redundant(project, excl_dist, -1, 0);
		project_delesect_edge_particles(project);
	}

	if ( delete_deselected ) part_delete_deselected(project);

	if ( projcomp )
		project_compare_particles(project, projcomp);

	if ( ext_scale )
		project_extract_particles(project, ext_scale, back_flag, norm_flag,
			fill_type, fill, 0, split, partbase, partpath, partext);
	
	if ( verbose ) project_show_selected(project);
	
	if ( project && outfile.length() )
		write_project(outfile, project, 0, 0);
	
	project_kill(project);
	if ( projcomp ) project_kill(projcomp);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}
