/**
@file	bcrystal.cpp
@brief	Calculates coordinates for multiple unit cells.
@author Bernard Heymann
@date	Created: 20001015
@date	Modified: 20080225
**/

#include "rwmolecule.h"
#include "mol_symmetry.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bcrystal [options] in.pdb out.pdb",
"----------------------------------------",
"Calculates coordinates for multiple unit cells.",
" ",
"Actions:",
"-translate 0,-50,22      Translate (angstrom)",
"-cells 2,5,3             Number of unit cells in each lattice direction",
" ",
"Parameters:",
"-verbose 7               Verbose output",
"-rename D                Rename molecules from the given letter.",
"-symmetry 1              Space group",
"-unitcell 50,50,50,90,90,90 Unit cell parameters (angstrom & degrees)",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star)",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    /* Initialize variables */
	char			first_name(0);				// First molecule name for renaming
    Vector3<double> 	t;							// Shift
	Vector3<long>	number(1,1,1);				// Number of unit cells in the 3 directions
	int 			spacegroup(1);
	UnitCell		uc(0,0,0,M_PI_2,M_PI_2,M_PI_2);
    Bstring    		atom_select("all");
	Bstring			paramfile;					// Use default parameter file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "translate" ) {
        	t = curropt->vector3();
       		if ( t.length() < 1 )
				cerr << "-translate: Three values must be specified!" << endl;
		}
		if ( curropt->tag == "cells" ) {
			number = curropt->vector3();
			if ( number.volume() < 1 )
				cerr << "-cells: Three values must be specified" << endl;
		}
		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
		if ( curropt->tag == "symmetry" )
			if ( ( spacegroup = curropt->value.integer() ) < 1 )
				cerr << "-symmetry: The space group number must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bstring		filename = argv[optind++];
	Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( molgroup == NULL ) {
		cerr << "Error: No coordinate file given!" << endl;
		bexit(-1);
	}
	
	if ( t[0] != 0 || t[1] != 0 || t[2] != 0 )
		molgroup_coor_shift(molgroup, t);
	
	if ( number.volume() > 1 )
		molgroup_generate_crystal(molgroup, uc, number);
		
    if ( first_name )
		molgroup_rename(molgroup, first_name);
    
	if ( optind < argc ) {
		molecule_update_comment(molgroup, argc, argv);
		write_molecule(argv[optind], molgroup);
	}

    molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

	

