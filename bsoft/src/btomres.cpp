/**
@file	btomres.cpp
@brief	Program to determine the resolution of a tilt series
@author  Bernard Heymann
@date	Created: 20031205
@date	Modified: 20200723
**/

#include "mg_tomo_resol.h"
#include "mg_tomo_rec.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: btomres [options] input.star [input2.star]",
"-------------------------------------------------",
"Program to determine the resolution of micrographs in a tilt series.",
"If the -micrograph option is used, only one micrograph will be processed",
"and the postscript file will contain the FSC plot."
"If the -micrograph option is not used, the whole tilt series will be",
"processed and the postrscript file will contain a plot of the tilt angle",
"versus the resolution, as well as a second plot, the NLOO3D FSC curve."
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling of input images (A/pixel; a single value can be given).",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
"-ratio 2.5               Radial sampling to Cartesian sampling ratio (default 1).",
"-scale 1.2               Scale or magnification of reconstruction compared to original images (default 1).",
"-size 100,120,90         Size of reconstruction (default from images).",
"-cutoff 0.3              Resolution cutoff for FRC (default 0.5).",
"-micrograph 56           Micrograph number to use for the resolution test (first = 0).",
"-fast 10                 Fast mode, only use micrographs within the given angle from the selected one (default all).",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
" ",
"Output:",
"-output file.star        Output parameter file name.",
"-reconstruction file.mrc Output 2D reconstruction file name.",
"-Postscript file.ps      Postscript output filename.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	Vector3<long>	size;						// Size of reconstruction
	double			scale(1);					// Scale or magnification of reconstruction compared to original images
	double 			resolution(0);				// Must be set > 0 to limit resolution
	double			sampling_ratio(1);			// Averaging ratio for FRC calculation
	double			cutoff(0.5);				// FRC cutoff
	int				mg_num(-1);					// Default use all micrographs
	double			fast_angle(M_PI);			// Select all micrographs to do a reconstruction
	int				ctf_action(0);				// Default no CTF operation
	double			wiener(0.2);				// Wiener for CTF correction
	Bstring			outfile;					// Output parameter file
	Bstring			recfile;					// Output reconstruction file
	Bstring			ps_file;					// Output postscript file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.1 )
				cerr << "-resolution: A positive resolution limit must be specified!" << endl;
		if ( curropt->tag == "ratio" )
			if ( ( sampling_ratio = curropt->value.real() ) < 1 )
				cerr << "-ratio: A ratio must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.1 )
				cerr << "-resolution: A positive resolution limit must be specified!" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.1 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( cutoff = curropt->value.real() ) < 0.01 )
				cerr << "-cutoff: A cuttoff must be specified!" << endl;
		if ( curropt->tag == "micrograph" )
			if ( ( mg_num = curropt->value.integer() ) < 0 )
				cerr << "-micrograph: A micrograph number must be specified!" << endl;
		if ( curropt->tag == "fast" ) {
			if ( ( fast_angle = curropt->value.real() ) < 0.1 )
				cerr << "-fast: An angle must be specified!" << endl;
			else
				fast_angle *= M_PI/180.0;
		}
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
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "reconstruction" )
			recfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
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

	if ( sam[0] > 0 ) {
		if ( sam[0] < 0.1 ) sam = Vector3<double>(1,1,1);
		project_set_mg_pixel_size(project, sam);
	}

	long   			i;
	Bfield*			field = project->field;
	Bmicrograph*	mg = NULL;
	Bmicrograph*	mg_res = NULL;
	Bimage*			prec = NULL;
	Bstring			title("Micrograph resolution for a tilt series");
	
	if ( mg_num >= 0 ) { 
		prec = mg_tomo_resolution(project, mg_num, resolution, sampling_ratio,
						scale, size, fast_angle, ctf_action, wiener, cutoff, ps_file);
	
		for ( mg = field->mg, i=0; mg; mg = mg->next, i++ )
			if ( i == mg_num ) mg_res = mg;
			else mg->select = 0;

		if ( prec && recfile.length() ) {
			prec->origin(0, 0, 0);
			prec->phase_shift_to_center();
			prec->fft_back();
			prec->statistics();
			prec->change_type(nudatatype);
			write_img(recfile, prec, 0);
			if ( mg_res ) {
				mg_res->fmg = recfile;
				mg_res->img_num = 0;
				mg_res->origin = prec->image->origin();
			}
		}

		delete prec;
		
	} else {
		vector<Bplot*>	plot = project_tomo_resolution(project,
				resolution, sampling_ratio, scale, size,
				fast_angle, ctf_action, wiener, cutoff);

		if ( ps_file.length() ) {
			ps_plot(ps_file, plot[0]);
			ps_file = ps_file.pre_rev('.') + "_nloo3d." + ps_file.post_rev('.');
			ps_plot(ps_file, plot[1]);
		}
		
		delete plot[0];
		delete plot[1];
	}

	if ( project && outfile.length() ) {
		write_project(outfile, project, 1, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}


