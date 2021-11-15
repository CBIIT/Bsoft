/**
@file	bcatmol.cpp
@brief	A program to concatenate coordinates from different files
@author Bernard Heymann
@date	Created: 20031018
@date	Modified: 20170123
**/

#include "rwmolecule.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bcatmol [options] input.pdb [input2.pdb ...]",
"---------------------------------------------------",
"Concatenates sets of coordinates with shifts from sets of files.",
" ",
"Actions:",
"-select CA               Atom selection (default all).",
"-rename D                Rename molecules from the given letter.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-box 10,10,10            Size of enclosing box (angstrom, default from coordinates).",
"-pbc                     Resolve periodic boundary bonds before catenation.",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
" ",
"Output:",
"-output output.pdb       Output coordinate file (default temp.pdb).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	Bstring			atom_select("all");		// Selection
	Vector3<double>	box;					// Size of enclosing box, default from coordinates
    int 			set_pbc(0);				// Flag for resolving wrapped coordinates
	char			first_name(0);			// First molecule name for renaming
	Bstring			catfile("temp.pdb");	// Default output file
	Bstring			paramfile;				// Use default parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*			curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
        	if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "box" ) {
			box = curropt->vector3();
			if ( box.volume() < 1 )
				cerr << "-box: All three dimensions must be specified" << endl;
		}
		if ( curropt->tag == "pbc" )
			set_pbc = 1; 		// Resolve periodic boundaries
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			catfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();

	if ( !catfile.length() ) {
		cerr << "Error: No output file name given!" << endl;
		bexit(-1);
	}

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No molecule files specified!" << endl;
		bexit(-1);
	}

	Bmolgroup*	molgroup = read_molecule(file_list, set_pbc, box, atom_select, paramfile);
	
	if ( first_name )
		molgroup_rename(molgroup, first_name);
	
	molecule_update_comment(molgroup, argc, argv);
	write_molecule(catfile, molgroup);

	molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

