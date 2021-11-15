/**
@file	bparam.cpp
@brief	A tool to extract parameters from coordinate files
@author Bernard Heymann
@date	Created: 20050304
@date	Modified: 20150512
**/

#include "rwmolecule.h"
#include "rwmd.h"
#include "mol_param.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bparam [options] param.star",
"----------------------------------",
"Manipulate molecular parameter files.",
" ",
"Actions:",
"-select CA               Atom selection (default all).",
"-bonds                   List all bond lengths.",
"-angles                  List all angles.",
"-length 1.85             Length to define bonds.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-elements                Use elements as identifiers rather than atom names.",
" ",
"Input:",
"-from input.pdb          Extract parameters from an atomic coordinate file.",
"-parameters parm.star    Molecular parameter file (default atom_prop.star)",
" ",
"Output:",
"-output param.star       Output parameter file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
    Bstring    		atom_select("ALL");
	int				show_bonds(0);
	int				show_angles(0);
	double			length_cutoff(0);
	int				elements(0);				// Flag to use element names
	Bstring			paramfile;					// Use default parameter file
	Bstring			atomfile;					// Atom coordinate file
	Bstring			outfile;					// Output parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
			if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "bonds" ) show_bonds = 1;
		if ( curropt->tag == "angles" ) show_angles = 1;
		if ( curropt->tag == "length" )
			length_cutoff = curropt->value.real();
		if ( curropt->tag == "elements" ) elements = 1;
		if ( curropt->tag == "from" )
			atomfile = curropt->filename();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bmd*		md = NULL;
	Bstring		filename;
	
	if ( optind < argc ) {
		filename = argv[optind];
		md = read_md_parameters(filename);
		if ( !md ) {
			cerr << "Error: Parameter file " << argv[optind] << " not read!" << endl;
			bexit(-1);
		}
	}

    Bmolgroup*	molgroup = NULL;
	if ( atomfile.length() ) {
		molgroup = read_molecule(atomfile, atom_select, paramfile);
		if ( !molgroup ) {
			cerr << "Error: Coordinate file not read!" << endl;
			bexit(-1);
		}
		md_kill(md);
		if ( length_cutoff ) {
			bond_kill(molgroup->bond);
			molgroup->bond = NULL;
			molgroup_bond_list_generate(molgroup, length_cutoff, 0);
		}
		md = md_calculate_parameters(molgroup, elements, show_bonds, show_angles);
		molgroup_kill(molgroup);
	}
	
	if ( outfile.length() ) {
//		molecule_update_comment(molgroup, argc, argv);
//		if ( write_molecule(argv[optind], molgroup) < 1 )
//			cerr << "Error: %s not written!\n", argv[optind]);
		write_md_parameters(outfile, md);
	}

	md_kill(md);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

