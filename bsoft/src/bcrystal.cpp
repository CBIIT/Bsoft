/**
@file	bcrystal.cpp
@brief	Calculates coordinates for multiple unit cells.
@author Bernard Heymann
@date	Created: 20001015
@date	Modified: 20230120
**/

#include "rwmodel.h"
#include "model_symmetry.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bcrystal [options] in.pdb",
"--------------------------------",
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
"Output:",
"-output newmod.star      New model file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    /* Initialize variables */
	char			first_name(0);				// First molecule name for renaming
    Vector3<double> t;							// Shift
	Vector3<long>	lattice(1,1,1);				// Number of unit cells in the 3 directions
	int 			spacegroup(1);
	UnitCell		uc;
    Bstring    		atom_select("all");
	Bstring			paramfile;					// Use default parameter file
	Bstring			outfile;					// Output parameter file name
    
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
			lattice = curropt->vector3();
			if ( lattice.volume() < 1 )
				cerr << "-cells: Three values must be specified" << endl;
		}
		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
		if ( curropt->tag == "symmetry" )
			if ( ( spacegroup = curropt->value.integer() ) < 1 )
				cerr << "-symmetry: The space group number must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	Bmodel*			model = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( file_list ) {
		model = read_model(file_list, paramfile);
		string_kill(file_list);
	}
	
	if ( !model ) {
		cerr << "Error: No model file read!" << endl;
		bexit(-1);
	}
	
	if ( t.length() )
		models_shift(model, t);
	
	if ( lattice.volume() > 1 ) {
		model_generate_lattice(model, uc, lattice);
		model_merge(model);
	}
		
    if ( first_name )
		model_rename(model, first_name);
    
	model_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model )
		write_model(outfile, model);

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

	

