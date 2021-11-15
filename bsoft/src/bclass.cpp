/**
@file	bclass.cpp
@brief	Classifies raw single particle images with respect to multiple models
@author Bernard Heymann
@date	Created: 20010222
@date	Modified: 20151008
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "mg_class.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bclass [options] input.star [input.star]",
"-----------------------------------------------",
"Classifies raw single particle images with respect to multiple reference maps.",
"The reference maps can be specified in the parameter file or on the command line.",
"FOM choices:",
"Correlation coefficient (FOM = CC)",
"R-factor (FOM = 1 - R)",
"Phase difference (FOM = cos(dPhi))",
" ",
"Actions:",
"-kernel 10,2             Reciprocal space projection: kernel size and power.",
"-CTF                     Apply CTF to projections.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-resolution 15,200       Resolution limits for comparisons (angstrom).",
"-fomtype phase           Type FOM to use for classification (default 0=CC=correlation,).",
"                         (other: 1=R=R-factor, 2=PD=phasedifference).",
"-cutoff 0.2              FOM cutoff to accept a particle (default 0).",
" ",
"Input:",
"-reference m1.map,m2.map List of reference maps (replaces any in the parameter file).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-projections             Output projection files (with insertion _map<n>).",
"-differences             Output difference files (with insertion _map<n>_diff).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				kernel_width(0);		// Reciprocal space kernel width
	int				kernel_power(2);		// Reciprocal space kernel power
	double			resolution_lo(1e30);	// Lower resolution limit (angstrom)
	double			resolution_hi(0);		// Upper resolution limit (angstrom)
	int				fom_type(0);			// Default correlation coefficient
	double			fom_cut(0); 			// Classify all particles
	int				ctf_apply(0);			// Flag to apply CTF to projections
	int				img_out(0);				// Don't output projection or difference images
	Bstring			ref_file;				// Comma-delimited list of reference maps
	Bstring			outfile;				// Output STAR file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "kernel" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-kernel: At least the kernel size must be specified!" << endl;
		if ( curropt->tag == "CTF" ) ctf_apply = 1;
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(resolution_hi, resolution_lo) < 1 )
				cerr << "-resolution: At least the high resolution limit must be specified!" << endl;
			else
				if ( resolution_hi > resolution_lo )
					swap(resolution_hi, resolution_lo);
		}
		if ( curropt->tag == "fomtype" ) {
			if ( isdigit(curropt->value[0]) ) {
				fom_type = curropt->value.integer();
			} else {
				curropt->value = curropt->value.lower();
				if ( curropt->value[0] == 'c' )
					fom_type = 0;
				else if ( curropt->value[0] == 'r' || curropt->value.contains("fac") )
					fom_type = 1;
				else if ( curropt->value[0] == 'p' )
					fom_type = 2;
			}
			if ( fom_type < 0 || fom_type > 2 ) fom_type = 0;
		}
		if ( curropt->tag == "cutoff" )
			if ( ( fom_cut = curropt->value.real() ) < 0.00001 )
				cerr << "-cutoff: A FOM cutoff must be specified!" << endl;
		if ( curropt->tag == "projections" )
			img_out += 1;
		if ( curropt->tag == "differences" )
			img_out += 2;
		if ( curropt->tag == "reference" )
			ref_file = curropt->value;
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

	if ( ref_file.length() ) {
		string_kill(project->reference);
		project->reference = NULL;
		project->reference = ref_file.split(",");
	}
	
	if ( !project->reference ) {
		cerr << "Error: No reference files found!\n" << endl;
		bexit(-1);
	}

	FSI_Kernel*		kernel = NULL;
	if ( kernel_width > 1 )
		kernel = new FSI_Kernel(kernel_width, kernel_power);
	
	mg_classify(project, resolution_hi, resolution_lo, fom_type, fom_cut, kernel, ctf_apply, img_out);
	
	// Write an output parameter file if a name is given
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);
	if ( kernel ) delete kernel;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

