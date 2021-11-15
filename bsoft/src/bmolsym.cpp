/**
@file	bmolsym.cpp
@brief	A tool to perform symmetry operations on coordinate files
@author Bernard Heymann
@date	Created: 19980214
@date 	Modified: 20150813
**/

#include "rwmolecule.h"
#include "mol_symmetry.h"
#include "mol_transform.h"
#include "mol_edit.h"
#include "mol_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmolsym [options] in.pdb out.pdb",
"---------------------------------------",
"Performs symmetry operations on molecules.",
" ",
"Actions:",
"-translate 0,-50,22      Translate (angstrom).",
"-rotate 0.5,0,5.8,45     Rotate around a vector (x,y,z) by an angle.",
"-select CA               Atom selection (default all).",
"-apply C5                Apply point group symmetry.",
"-find C5                 Find standard orientation for this point group symmetry.",
"-Bfactor D5              Calculate B factors for this point group symmetry.",
"-pdbsymmetry             Read SMTRY matrices from a PDB file and apply them.",
"-pdbbiomt                Read BIOMT matrices from a PDB file and apply them.",
"-rename D                Rename molecules from the given letter.",
"-show                    Show operational symmetry matrices.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-origin 10,-10,20        Origin for rotation (default 0,0,0, can be \"center\" for COM).",
"-reference 0,1.5,-0.2,35 Reference symmetry axis and rotation angle (default 0,0,1,0).",
"-distance 3.5            Distance between overlapping atoms allowed (default 3 angstrom, use with -apply).",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
"-similarity sim.star     Residue similarity matrix for finding orientation (default blosum62.star).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	char			first_name(0);			// First molecule name for renaming
	int				com_origin(0);			// Flag to set the center-of-mass as origin of rotation
    Transform		t;						// Rotation origin, axis and angle
    Bstring    		atom_select("ALL");
	View			ref_view;				// Reference view for symmetry application
	Bsymmetry		sym;					// Point group
	Bsymmetry		sym_B;					// Symmetry for B factor calculation
	int				sym_flag(0);
	int				pdbsym_flag(0);			// Flag to get SMTRY matrices from a PDB file
	int				pdbbio_flag(0);			// Flag to get BIOMT matrices from a PDB file
	double			distance(0);			// Distance allowed between overlapping atoms in angstrom
	int				show(0);				// Flag to show operational matrices
	Bstring			paramfile;				// Use default parameter file
	Bstring			simfile;				// Use default similarity matrix file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "translate" )
			t.trans = curropt->vector3();
		if ( curropt->tag == "rotate" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else {
				t.axis = Vector3<double>(d[0], d[1], d[2]);
				t.angle = d[3] * M_PI/180.0;
			}
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value.contains("cen") )
				com_origin = 1;
			else
				t.origin = curropt->origin();
		}
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
        	if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "apply" ) {
 			sym = curropt->symmetry();
			if ( sym.point() < 102 )
				cerr << "-apply: A point group must be specified!" << endl;
			else
				sym_flag = 1;
		}
		if ( curropt->tag == "find" ) {
 			sym = curropt->symmetry();
			if ( sym.point() < 102 )
				cerr << "-find: A point group must be specified!" << endl;
			else
				sym_flag = 2;
		}
		if ( curropt->tag == "Bfactor" )
 			sym_B = curropt->symmetry();
		if ( curropt->tag == "pdbsymmetry" ) pdbsym_flag = 1;
		if ( curropt->tag == "pdbbiomt" ) pdbbio_flag = 1;
		if ( curropt->tag == "reference" )
			ref_view = curropt->view();
		if ( curropt->tag == "distance" )
			if ( ( distance = curropt->value.real() ) < 0.1 )
				cerr << "-distance: A separation distance must be specified!" << endl;
		if ( curropt->tag == "show" )
			show = 1;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "similarity" )
			simfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	if ( show ) sym_show_operational_matrices(sym);
	
	if ( optind >= argc ) bexit(-1);
    
	Bstring		filename(argv[optind++]);
    Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: Coordinate file not read!" << endl;
		bexit(-1);
	}
	
	if ( com_origin )
		t.origin = molgroup_center_of_mass(molgroup);
	
	if ( t.angle == 0 )
		molgroup_coor_shift(molgroup, t.trans);
	else
		molgroup_coor_rotate(molgroup, t);

	molgroup_stats(molgroup);
	
	Bresidue_matrix*	simat;
	if ( sym_flag == 2 && sym.point() > 101 ) {
//		molgroup_find_standard_view(molgroup, sym, ref_view);
		simat = get_residue_matrix(simfile);
		molgroup_orient_to_standard_view(molgroup, sym, ref_view, simat);
		residue_matrix_kill(simat);
	}
	
	if ( sym_flag == 1 && sym.point() > 101 ) {
		molgroup_apply_point_group(molgroup, sym, ref_view);
		if ( distance ) molgroup_remove_overlapping_atoms(molgroup, distance);
	}
	
	if ( pdbsym_flag )
		molgroup_apply_symmetry_from_pdb(molgroup, filename);
	
	if ( pdbbio_flag )
		molgroup_apply_matrices_from_pdb(molgroup, filename);

	if ( molgroup && sym_B.point() > 101 )
		molgroup_symmetry_B(molgroup, sym_B);
	
	if ( molgroup && sym.point() > 101 )
		molgroup_symmetry_RMSD(molgroup, sym);
	
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
