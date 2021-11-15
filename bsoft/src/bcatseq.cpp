/**
@file	bcatseq.cpp
@brief	A program to concatenate protein sequences from different files
@author Bernard Heymann
@date	Created: 20000814
@date	Modified: 20060722
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
"Usage: bcatseq [options] input.seq input2.seq",
"---------------------------------------------",
"Concatenates protein sequences from different files.",
"Multiple sequences in the different files must be in the same order.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Output:",
"-output output.embl      Output sequence file (default temp.star).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	Bstring			atom_select("all");			// Selection
	Bstring			catfile("temp.star");		// Default output file
	Bstring			paramfile;					// Use default parameter file
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*			curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "output" )
			catfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();
	
	// Read all the sequence files
	Bmolgroup*		catgroup = NULL;
	Bmolgroup*		molgroup = NULL;
	Bmolecule*		catmol = NULL;
	Bmolecule*		mol = NULL;
	Bstring			filename;
	
	i = 0;
	while ( optind < argc ) {
		filename = argv[optind++];
		molgroup = read_molecule(filename, atom_select, paramfile);
		if ( i == 0 ) catgroup = molgroup;
		else {
			for ( mol = molgroup->mol, catmol = catgroup->mol; mol; mol = mol->next, catmol = catmol->next ) {
				catmol->seq += mol->seq;
				catmol->id += mol->id;
				catmol->nres += mol->nres;
			}
			molgroup_kill(molgroup);
		}
		i++;
	}
	
	write_molecule(catfile, catgroup);
	
	molgroup_kill(catgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

