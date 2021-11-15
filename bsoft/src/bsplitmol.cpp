/**
@file	bsplitmol.cpp
@brief	A program to concatenate coordinates from different files
@author Bernard Heymann
@date	Created: 20050119
@date 	Modified: 20050119
**/

#include "rwmolecule.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bsplitmol [options] input.pdb output.pdb",
"-----------------------------------------------",
"Splits sets of coordinates to different numbered files.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-first 5                 Number given to the first file (default 0).",
"-digits 3                Number of digits inserted before the last period in the output file name.",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	Bstring			atom_select("all");			// Selection
	int 			first_number(0);			// Number given to first file
	int 			digits(3);					// File number size
	Bstring			paramfile = NULL;			// Use default parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "first" )
			if ( ( first_number = curropt->value.integer() ) < 0 )
				cerr << "-first: A number must be specified!" << endl;
		if ( curropt->tag == "digits" )
			if ( ( digits = curropt->value.integer() ) < 1 )
				cerr << "-digits: A number of digits must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();

	Bstring			inputfile = argv[optind++];
	Bmolgroup*		molgroup = read_molecule(inputfile, atom_select, paramfile);
	if ( molgroup == NULL ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
    
	if ( optind >= argc ) bexit(0);
	
	long			i, nmol = 0;
	Bmolecule*		mol, *mol2;
	for ( mol=molgroup->mol; mol; mol=mol->next ) nmol++;
	
	if ( verbose & VERB_LABEL )
		cout << "Splitting " << inputfile << " up into " << nmol << " files" << endl << endl;
	
	Bstring		outputfile = argv[optind];
	Bstring		filename;
	char		format[32];
	
	if ( digits < 0.44*log(1.0*nmol) + 1 ) digits = (int) (0.44*log(1.0*nmol) + 1);
	snprintf(format, 32, "_%%0%dd.", digits);
	
	Bmolgroup*		molgroup_out = molgroup_init();
	molecule_update_comment(molgroup_out, argc, argv);
	for ( i=0, mol=molgroup->mol; mol; mol=mol2, i++ ) {
		mol2 = mol->next;
		mol->next = NULL;
		molgroup_out->mol = mol;
		filename = outputfile.pre_rev('.') + Bstring(i+first_number, format) + outputfile.post_rev('.');
		write_molecule(filename, molgroup_out);
		molecule_kill(mol);
		molgroup_out->mol = NULL;
	}
	
	molgroup->mol = NULL;
	molgroup_kill(molgroup);
	molgroup_kill(molgroup_out);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

