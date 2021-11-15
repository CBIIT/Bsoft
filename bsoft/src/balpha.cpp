/**
@file	balpha.cpp
@brief	Program to make and analyze alpha helices.
@author Bernard Heymann
@date	Created: 20050315
@date 	Modified: 20050322
**/

#include "rwmolecule.h"
#include "mol_alpha.h"
#include "mol_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: balpha [options] [input.pdb] output.pdb",
"----------------------------------------------",
"Creates and manipulates alpha helices.",
" ",
"Actions:",
"-create 12               Create an alpha helix of a given length.",
"-helix 213,236           Set residues to a helix (all molecules).",
"-find                    Find a helix center and orientation (only for a single helix).",
"-set                     Set a helix to the standard orientation (only for a single helix).",
"-consolidate             Consolidate a set of helices to eliminate redundancy.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	int				create_number(0);				// Length of helix to create
	int				helix_start(0), helix_end(0);	// Setting to a helix
	int				find_orientation(0);			// Flag to find the helix orientation
	int				set_orientation(0);				// Flag to set standard orientation
	int				consolidate(0);					// Flag to consolidate a set of helices
	Bstring			atom_select("all");				// Selection
	Bstring			paramfile;
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( curropt->tag == "create" )
			if ( ( create_number = curropt->value.integer() ) < 1 )
				cerr << "-create: The desired helix length must be specified!" << endl;
 		if ( curropt->tag == "helix" )
			if ( curropt->values(helix_start, helix_end) < 2 )
				cerr << "-helix: Both the start and end of the helix must be specified!" << endl;
 		if ( curropt->tag == "find" )
        	find_orientation = 1;
 		if ( curropt->tag == "set" )
        	find_orientation = set_orientation = 1;
 		if ( curropt->tag == "consolidate" )
        	consolidate = 1;
    }
	option_kill(option);

	double		ti = timer_start();
	
    // Read the molecule file
    Bmolgroup*	molgroup = NULL;
	Bstring		filename(argv[optind++]);
	
	if ( create_number > 0 )
		molgroup = molgroup_generate_alpha_helix(create_number);
	else
		molgroup = read_molecule(filename, atom_select, paramfile);
	
	if ( !molgroup ) {
		cerr << "Error: No molecules created or read!" << endl;
		bexit(-1);
	}
	
	if ( helix_end )
		molgroup_set_alpha_helix(molgroup, helix_start, helix_end);
	
	molgroup_find_helical_axes(molgroup);

	if ( find_orientation )
		mol_find_alpha_orientation(molgroup->mol, set_orientation);
	
	if ( consolidate )
		molgroup = molgroup_consolidate_alpha(molgroup);
	
    // Write an output file if a file name is given
	if ( optind < argc ) {
		molecule_update_comment(molgroup, argc, argv);
		write_molecule(argv[optind], molgroup);
	}
	
	molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

