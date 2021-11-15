/**
@file	bmolhel.cpp
@brief	A program to perform helical symmetry operations on coordinate files
@author Bernard Heymann
@date	Created: 20090121
@date 	Modified: 20130426
**/

#include "rwmolecule.h"
#include "mol_symmetry.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmolhel [options] in.pdb out.pdb",
"---------------------------------------",
"Performs helical symmetry operations on molecules.",
" ",
"Actions:",
"-generate 5,7            Generate new molecules downwards and upwards.",
"-translate 0,-50,22      Translate (angstrom)",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-origin 10,-10,20        Origin for rotation (default 0,0,0, can be \"center\" for COM).",
"-helix 25.3,67.1         Helical rise and angle.",
"-reference 0,1.5,-0.2,35 Reference symmetry axis and rotation angle (default 0,0,1,0).",
"-rename D                Rename molecules from the given letter.",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int				gen_down(0), gen_up(0);	// Number of molecules to generate
    Vector3<double> 	t;						// Shift
	char			first_name(0);			// First molecule name for renaming
	int				com_origin(0);			// Flag to set the center-of-mass as origin of rotation
	Vector3<double>	origin;					// Rotation origin
	double			helix_rise(0);			// Rise per asymmetric unit
	double			helix_angle(0);			// Rotation angle per asymmetric unit
    Bstring    		atom_select("ALL");
	View			ref_view;				// Reference view for symmetry application
	Bstring			paramfile;				// Use default parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "generate" )
			if ( curropt->values(gen_down, gen_up) < 1 )
				cerr << "-generate: At least one number must be specified!" << endl;
		if ( curropt->tag == "translate" ) {
			t = curropt->vector3();
        	if ( t.length() < 0.0001 )
				cerr << "-translate: Three values must be specified!" << endl;
		}
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value.contains("cen") )
				com_origin = 1;
			else
				origin = curropt->origin();
		}
		if ( curropt->tag == "helix" ) {
			if ( curropt->values(helix_rise, helix_angle) < 2 )
				cerr << "-helix: Both rise and angle must be specified" << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "reference" )
			ref_view = curropt->view();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bstring		filename(argv[optind++]);
    Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: Coordinate file not read!" << endl;
		bexit(-1);
	}
	
	if ( t[0] != 0 || t[1] != 0 || t[2] != 0 )
		molgroup_coor_shift(molgroup, t);
	
	if ( com_origin )
		origin = molgroup_center_of_mass(molgroup);
	
	if ( origin[0] || origin[1] || origin[2] )
		molgroup_coor_shift(molgroup, -origin);
	
	molgroup_stats(molgroup);
	
	if ( gen_down + gen_up > 0 && helix_rise && helix_angle )
		molgroup_generate_helix(molgroup, ref_view, helix_rise, helix_angle, gen_down, gen_up);
	
	molgroup_from_molgroup_list(molgroup);
	
    if ( first_name )
		molgroup_rename(molgroup, first_name);

	if ( optind < argc ) {
		molecule_update_comment(molgroup, argc, argv);
		if ( write_molecule(argv[optind], molgroup) < 1 )
			cerr << "Error: " << argv[optind] << " not written!" << endl;
	}

	molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

