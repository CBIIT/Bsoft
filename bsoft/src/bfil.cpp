/**
@file	bfil.cpp
@brief	Extract filaments from a micrograph as straightened filaments or boxed particles.
@author Bernard Heymann
@date	Created: 20080513
@date	Modified: 20210414
**/

#include "mg_processing.h"
#include "mg_helix.h"
#include "mg_extract.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bfil [options] input.star [input2.star ...]",
"--------------------------------------------------",
"Extracts filaments from a micrograph as straightened filaments or boxed particles.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-center 85               Center filament nodes with the given window width.",
"-density 25              Estimate filament density per length with the given filament width.",
"-extract 150,y           Extract filaments with this width (voxels) along an axis (optional).",
"-split                   Split filaments into individual files.",
" ",
"Actions for boxing:",
"-box 125,125,1           Box particles along filaments with this size (voxels, one sets all).",
"-interval 25.5           Interval of boxes (default: box radius).",
"-maskwidth 50            Width to mask boxed images (voxels).",
"-background              Correct background of extracted particles.",
"-normalize               Normalize extracted particles.",
"-rotate y                Rotate particles from filaments to align with an axis.",
"-powerspectrum 200       Calculates average power spectra from extracted particles with optional padding.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-base fil_001            Set the filament/particle base file name (no extension).",
"-path dir/subdir         Set the filament/particle file paths.",
"-extension pif           Set the filament/particle image file format.",
" ",
"Parameters for boxing:",
"-helix 16.2,39.5         Helical parameters for filament boxing: helical rise ",
"                         and angle per asymmetric unit (angstroms and degrees).",
"-fill 127                Fill value for particles extending beyond image (default background).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-Densities img.mrc       One dimensional image with filament densities.",
"-Postscript plot.ps      Histogram of filament densities.",
"-export particles.txt,1  Export selected particle coordinates to a text file.",
"                         Flag: convert coordinates to physical location in angstroms.",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize optional variables
	int				use_rec(0);					// Flag to process reconstructions
	int				split(0);					// Flag to split output into individual files
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Bstring			base;						// Particle base file name
	Bstring			path;						// Particle file path
	Bstring			ext;						// Particle image file format
	int				center_width(0);			// Filament width for centering nodes
	int				filament_density(0);		// Filament width for density estimates
	int				filament_width(0);			// Filament width for extraction
	int				filament_axis(0);			// Filament axis for extraction
	char			axis_char;					// Filament axis character
	double			boxing_interval(0);			// Boxing interval for defining particles along a filament
	int				mask_width(0);				// Filament mask width
	double			helix_rise(0), helix_angle(0);	// Helical parameters for filaments
	Vector3<long>	box_size;					// Box size for extracting particles
	int				back_flag(0);				// Background correction for extracting particles
	int				norm_flag(0);				// Normalization for extracting particles
	int				powerspectrum(0);			// Flag to calculate average power spectra
	int				pad(0);						// Padding for power spectrum
	int 			fill_type(FILL_BACKGROUND);
	double			fill(0);
	int				rotation_axis(0);			// Axis to rotate particles from filaments, 0=no rotation
	Bstring			outfile;					// Output parameter file name
	Bstring			densfile;					// Output filament densities file
	Bstring			psfile;						// Postscript output file
	Bstring			expfile;					// Particle coordinate export file
	int				exp_flags(0);				// Export flags
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "split" ) split = 1;
		if ( curropt->tag == "center" )
    	    if ( ( center_width = curropt->value.integer() ) < 1 )
				cerr << "-center: A filament width must be specified." << endl;
		if ( curropt->tag == "density" )
			if ( ( filament_density = curropt->value.integer() ) < 1 )
				cerr << "-density: A filament width must be specified." << endl;
		if ( curropt->tag == "extract" ) {
			filament_width = curropt->value.integer();
			if ( filament_width < 1 )
 				cerr << "-extract: A filament width must be specified." << endl;
			else if ( curropt->value.contains(",") ) {
				axis_char = curropt->value.post(',')[0];
				if ( axis_char == '1' || axis_char == 'x' ) filament_axis = 1;
				if ( axis_char == '2' || axis_char == 'y' ) filament_axis = 2;
				if ( axis_char == '3' || axis_char == 'z' ) filament_axis = 3;
			}
		}
		if ( curropt->tag == "box" ) {
			box_size = curropt->size();
    	    if ( box_size.volume() < 1 )
				cerr << "-box: A boxing size must be specified." << endl;
		}
		if ( curropt->tag == "interval" )
    	    if ( ( boxing_interval = curropt->value.real() ) < 1 )
				cerr << "-interval: A boxing interval must be specified." << endl;
		if ( curropt->tag == "maskwidth" )
    	    if ( ( mask_width = curropt->value.integer() ) < 1 )
				cerr << "-maskwidth: A mask width must be specified." << endl;
		if ( curropt->tag == "background" ) back_flag = 1;
		if ( curropt->tag == "normalize" ) norm_flag = 1;
		if ( curropt->tag == "powerspectrum" ) {
			powerspectrum = 1;
			if ( ( pad = curropt->value.integer() ) < 1 )
				cerr << "-powerspectrum: A padding value must be specified!" << endl;
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "base" ) {
			base = curropt->value;
			if ( base.length() < 1 )
				cerr << "-base: The particle base file name must be specified!" << endl;
		}
		if ( curropt->tag == "path" ) {
			path = curropt->value;
			if ( path.length() < 1 )
				cerr << "-path: The particle file path must be specified!" << endl;
			else
				if ( path[-1] != '/' ) path += "/";
		}
		if ( curropt->tag == "extension" ) {
			ext = curropt->value;
			if ( ext.length() < 1 )
				cerr << "-extension: The particle file extension must be specified!" << endl;
		}
		if ( curropt->tag == "helix" ) {
    	    if ( curropt->values(helix_rise, helix_angle) < 2 )
				cerr << "-helix: Both helical rise and angle per unit must be specified." << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "rotate" ) {
			if ( ( axis_char = curropt->value[0] ) < 1 )
				cerr << "-rotate: An axis must be specified!" << endl;
			else {
				if ( axis_char == '1' || axis_char == 'x' ) rotation_axis = 1;
				if ( axis_char == '2' || axis_char == 'y' ) rotation_axis = 2;
				if ( axis_char == '3' || axis_char == 'z' ) rotation_axis = 3;
			}
		}
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "Densities" )
			densfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
		if ( curropt->tag == "export" ) {
			if ( curropt->value.contains(",") ) {
				expfile = curropt->value.pre(',');
				exp_flags = curropt->value.post(',').integer();
			} else {
				expfile = curropt->value;
			}
		}
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
	else box_size[2] = 1;
	
	Bimage*			pdens = NULL;
	Bplot*			plot = NULL;
	
	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	if ( box_size[0] > 0 )
		project_set_particle_box_size(project, box_size);
		
	if ( center_width > 0 )
		project_center_filaments(project, center_width);

	if ( filament_density ) {
		pdens = project_filament_density(project, filament_density);
		plot = pdens->histogram_gauss_plot(100, 2);
//		if ( verbose )
//			cout << pdens->average() << tab << pdens->standard_deviation() << endl;
		pdens->statistics();
		if ( densfile.length() ) write_img(densfile, pdens, 0);
		delete pdens;
	}
	
	if ( filament_width > 0 )
		project_extract_filaments(project, filament_width, filament_axis, base, path, ext, split);
		
	if ( path.length() && box_size[0] )
		project_set_particle_path(project, path);
	
	if ( boxing_interval > 0 )
		project_filaments_to_particles(project, box_size, boxing_interval, helix_rise, helix_angle);

	if ( box_size[0] > 0 ) {
		if ( rotation_axis > 0 ) {
			project_extract_particles(project, 1, back_flag, norm_flag,
				fill_type, fill, 0, project->select, base, path, ext);
			project_rotate_mask_filament_particles(project, rotation_axis, back_flag, mask_width);
		} else {
			project_extract_particles(project, 1, back_flag, norm_flag,
				fill_type, fill, mask_width, project->select, base, path, ext);
		}
	}
	
	if ( powerspectrum )
		project_filament_powerspectrum(project, pad, rotation_axis, path);

	if ( expfile.length() )
		write_particle_list(expfile, project, exp_flags);

	if ( outfile.length() )
		write_project(outfile, project, 0, 0);
	
	if ( psfile.length() && plot )
		ps_plot(psfile, plot);

	project_kill(project);
	delete plot;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	return 0;
}

