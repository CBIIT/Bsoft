/**
@file	bspr.cpp
@brief	3D reconstruction from single particle images
@author	Bernard Heymann
@date	Created: 20120212
@date	Modified: 20151031
**/

#include "mg_processing.h"
#include "mg_reconstruct.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bspr [options] input.star [input.star]",
"---------------------------------------------",
"Analyzing single particle reconstructions.",
"The mask should have 4 levels.",
"The relative density is calculated as (level1 - level3)/(level2 - level3).",
" ",
"Actions:",
"-all                     Reset selection to all particles.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-postsymmetrize          Symmetrization after reconstruction (default during).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-select 3                Selection number from STAR file selection column (default 1).",
"-sampling 1.5,1.5,1.5    Sampling of input images (A/pixel; a single value can be given).",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
"-scale 1.2               Scale or magnification of reconstruction compared to original images (default 1).",
"-size 100,120,90         Size of reconstruction (default from images).",
"-interpolation weighted  Interpolation type: nearest (default), weighted, trilinear.",
"-pad 3                   Image padding factor (default 2).",
"-symmetry D6             Symmetry: Point group identifier",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
" ",
"Input:",
"-Mask mask.pif           Mask with 4 levels: 0=ignore, 1=analyze, 2=reference density, 3=background.",
" ",
"Output:",
"-output file.star        Output parameter file name.",
"-ppx                     Write temporary particle parameter files to directory \"ppx\".",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
//	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			nuavg, nustd(0); 		// Values for rescaling
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Vector3<long>	map_size;				// Size of reconstruction
	int 			reset(0);				// Keep selection as read from file
	int				interp_type(0);			// Interpolation type
	Vector3<double>	scale(1,1,1);			// Scale or magnification of reconstruction compared to original images
	int				ctf_action(0);			// Default no CTF operation
	double			wiener(0.2);			// Wiener for CTF correction
	int				pad_factor(2);			// Image padding factor
	double 			resolution(0); 			// Must be set > 0 to limit resolution
	Bsymmetry		sym;				// No symmetry specified
//	int 			sym_mode(0);			// Symmetrization during reconstruction
	int 			part_select(-1000); 	// Selection number from parameter file selection column
	int				flags(0);				// Flags for processing options
	Bstring			maskfile;				// Mask file
	Bstring			outfile;				// Output parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" )
			reset = 1;
		if ( curropt->tag == "interpolation" ) {
			if ( curropt->value[0] == 'n' ) interp_type = 0;
			if ( curropt->value[0] == 't' ) interp_type = 1;
			if ( curropt->value[0] == 'k' ) interp_type = 2;
		}
//		if ( curropt->tag == "datatype" )
//			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 1 )
				cerr << "-select: A selection number must be specified!" << endl;
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
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "size" )
			map_size = curropt->size();
		if ( curropt->tag == "scale" )
			scale = curropt->scale();
		if ( curropt->tag == "interpolation" ) {
			if ( curropt->value[0] == 'n' ) interp_type = 0;
			if ( curropt->value[0] == 'w' ) interp_type = 1;
			if ( curropt->value[0] == 't' ) interp_type = 2;
		}
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: An integer factor must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
//		if ( curropt->tag == "postsymmetrize" )
//			sym_mode = 1;
		if ( curropt->tag == "ppx" ) flags |= WRITE_PPX | CHECK_PPX;
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
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

	if ( reset ) part_reset_selection(project, 3);

	if ( sam[0] > 0 ) {
		if ( sam[0] < 0.1 ) sam[0] = 1;
		project_set_mg_pixel_size(project, sam);
	}

	View				ref_view;
	Bstring				id;

	if ( maskfile.length() )
		project_single_particle_reconstruction(project, maskfile, sym, part_select,
				resolution, scale, map_size, pad_factor, 
				interp_type, ctf_action, wiener, flags);
	
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

