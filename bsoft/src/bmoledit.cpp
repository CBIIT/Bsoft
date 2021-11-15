/**
@file	bmoledit.cpp
@brief	A tool to edit coordinate files
@author Bernard Heymann
@date	Created: 19980214
@date 	Modified: 20211029
**/

#include "rwmolecule.h"
#include "mol_symmetry.h"
#include "mol_transform.h"
#include "mol_select.h"
#include "mol_edit.h"
#include "mol_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

long		molgroup_trim(Bmolgroup* molgroup, Vector3<double> box);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmoledit [options] in.pdb out.pdb",
"----------------------------------------",
"Manipulates atomic coordinate files.",
" ",
"Actions:",
"-rename D                Rename molecules from the given letter.",
"-renumber 52             Renumber residues from the given number.",
"-select CA               Atom selection (default all).",
"-elements                Use elements as identifiers for atom types.",
"-split                   Split molecule groups, output multiple numbered files.",
"-slices 23.7             Split coordinates into slices of this thickness, output multiple numbered files.",
"-translate 0,-50,22      Translate (angstrom).",
"-trim 200,230,180        Trim to new box size (use with -translate).",
"-noH                     Remove all hydrogens.",
"-SS 2.2                  Add disulfide bonds within the given separation.",
"-nobonds                 Remove all bond specifications.",
"-bonds 2.1               Add bonds between atoms closer than the given distance.",
"-pseudobonds 3           Add pseudo-atoms on bonds given the number per bond.",
"-random 1.2              Add a random displacement to each coordinate within a maximum (angstrom).",
"-B 23.4                  Add a random displacement to each coordinate based on a B-factor (angstrom^2).",
"-displace 120,2.7        Displacement a given number of atoms given a standard deviation (angstrom).",
"-prune                   Remove overlapping molecules.",
"-untangle 3.5,0.2        Eliminate overlaps by moving molecules apart (sampling and damping factor).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-distance 3.5            Distance between overlapping atoms allowed (default 3 angstrom, use with -insert).",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
"-insert in.pdb           Insert this into the input file, deleting overlapping atoms in the input file.",
"                         (use with -distance to indicate cutoff to eliminate atoms).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int				elements(0);		// Flag to use element names
	char			first_name(0);		// First molecule name for renaming
    int     		first_res(0);		// Don't renumber
	int				split(0);			// Flag to split molecule groups
    double	 		slice_thickness(0);	// Slice thickness for splitting up coordinates
    Transform		t;					// Transformation details
	Vector3<double>	trim;				// New enclosing box
    int				removeH(0);			// Flag to remove hydrogens
	double			SSdistance(0);		// Maximum distance between sulfur atoms
    int				removeBonds(0);		// Flag to remove bonds
    double			addBonds(0);		// Cutoff distance to add bonds, 0=don't add
    int				pseudobonds(0);		// Pseudo-atoms to add
	Bstring			insert;				// Molecules to insert
	double			random_max(0);		// Random displacement maximum in angstrom
	double			B(0);				// B-factor in angstrom^2
	double			random_number(0);	// Random displacement number of atoms
	double			random_stdev(0);	// Random displacement standard deviation in angstrom
	double			distance(0);		// Distance allowed between overlapping atoms in angstrom
    int 			set_wrap(0);		// Flag for wrapping coordinates
	int				prune(0);			// Flag to remove overlapping molecules
	double			untangle(0);		// Sampling for moving overlapping molecules
	double			lambda(0.1);		// Damping factor for moving overlapping molecules
	Vector3<double>	box;				// Enclosing box
    Bstring    		atom_select("ALL");
	Bstring			paramfile;			// Use default parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
        	if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "elements" ) elements = 1;
		if ( curropt->tag == "split" ) split = 1;
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "renumber" )
        	if ( ( first_res = curropt->value.integer() ) < 1 )
				cerr << "-renumber: The first residue number must be specified!" << endl;
		if ( curropt->tag == "slices" ) {
        	if ( ( slice_thickness = curropt->value.real() ) < 1 )
				cerr << "-slices: A slice thickness must be specified!" << endl;
			else
				if ( slice_thickness <= 0 ) slice_thickness = 0;
		}
		if ( curropt->tag == "translate" )
			t.trans = curropt->vector3();
		if ( curropt->tag == "trim" )
			trim = curropt->size();
		if ( curropt->tag == "noH" )
        	removeH = 1;
		if ( curropt->tag == "SS" )
			if ( ( SSdistance = curropt->value.real() ) < 0.1 )
				cerr << "-SS: A maximum separation distance must be specified!" << endl;
		if ( curropt->tag == "nobonds" )
        	removeBonds = 1;
		if ( curropt->tag == "bonds" )
			if ( ( addBonds = curropt->value.real() ) < 0.1 )
				cerr << "-bonds: A maximum bond length must be specified!" << endl;
		if ( curropt->tag == "pseudobonds" )
			if ( ( pseudobonds = curropt->value.integer() ) < 1 )
				cerr << "-pseudobonds: A number of pseudo-atoms must be specified!" << endl;
		if ( curropt->tag == "insert" )
        	insert = curropt->filename();
		if ( curropt->tag == "random" )
			if ( ( random_max = curropt->value.real() ) < 0.1 )
				cerr << "-random: A random displacement maximum must be specified!" << endl;
		if ( curropt->tag == "B" )
			if ( ( B = curropt->value.real() ) < 0.1 )
				cerr << "-B: A B-factor must be specified!" << endl;
		if ( curropt->tag == "displace" )
			if ( curropt->values(random_number, random_stdev) < 2 )
				cerr << "-displace: A number and standard deviation must be specified!" << endl;
		if ( curropt->tag == "prune" ) prune = 1;
		if ( curropt->tag == "untangle" )
			if ( curropt->values(untangle, lambda) < 1 )
				cerr << "-untangle: Grid sampling in angstrom must be specified!" << endl;
		if ( curropt->tag == "distance" )
			if ( ( distance = curropt->value.real() ) < 1 )
				cerr << "-distance: A separation distance must be specified!" << endl;
		if ( curropt->tag == "wrap" ) {
			box = curropt->size();
        	set_wrap = 1; 		// Wrap within periodic boundaries
		}
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
    }
	option_kill(option);
	
	if ( optind >= argc ) bexit(-1);
    
	double		ti = timer_start();
	
	Bstring		filename = argv[optind++];
    Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: Coordinate file not read!" << endl;
		bexit(-1);
	}
	
	if ( elements )
		molgroup_set_atom_types_to_elements(molgroup);
	
    if ( first_name )
		molgroup_rename(molgroup, first_name);
    
    if ( first_res > 0 )
		molgroup_residue_renumber(molgroup, first_res);
    
	if ( t.trans.length() ) molgroup_coor_shift(molgroup, t.trans);

	if ( trim.volume() ) molgroup_trim(molgroup, trim);
	else if ( set_wrap ) molgroup->box = box;

	if ( removeH )
		molgroup_remove_hydrogens(molgroup);

	if ( SSdistance )
		molgroup_add_disulfides(molgroup, SSdistance);
	
	if ( removeBonds ) {
		bond_kill(molgroup->bond);
		molgroup->bond = NULL;
	}
	
	if ( addBonds )
		molgroup_bond_list_generate(molgroup, addBonds, set_wrap);
	
	if ( pseudobonds )
		molgroup_bond_pseudo_atoms(molgroup, pseudobonds, set_wrap);
	
    Bmolgroup*	molinsert = NULL;
	if ( insert.length() ) {
    	molinsert = read_molecule(insert, atom_select, paramfile);
		molgroup_insert(molgroup, molinsert, distance);	// Remember that the group inserted is deallocated
		molgroup_kill(molinsert);
	}
	
	if ( random_max )
		molgroup_randomize(molgroup, random_max);
	
	if ( B )
		molgroup_randomize_B(molgroup, B);
		
	if ( random_stdev )
		molgroup_random_displace_number(molgroup, random_number, random_stdev);
		
	if ( prune )
		molgroup_prune_molecules(molgroup);

	if ( untangle )
		molgroup_untangle_molecules(molgroup, untangle, lambda);
	
	int				i, num(0);
	Bmolgroup**		molgroup_set;
	if ( optind < argc ) {
		filename = argv[optind];
		molecule_update_comment(molgroup, argc, argv);
		if ( split ) {
			molecules_to_molgroups(molgroup);
			molgroup_list_write(filename, molgroup);
		} else if ( slice_thickness > 0 ) {
			molgroup_set = molgroup_split_into_slices(molgroup, slice_thickness, num);
			for ( i=0; i<num; i++ ) {
				filename = filename.pre_rev('.') + Bstring(i+1, "_%03d.") + filename.post_rev('.');
				write_molecule(filename, molgroup_set[i]);
				molgroup_kill(molgroup_set[i]);
			}
			delete[] molgroup_set;
		} else {
			i = write_molecule(argv[optind], molgroup);
			if ( i < 1 )  {
				cerr << "Error: No output file written!" << endl;
				bexit(-1);
			}
		}
	}

	molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief  Deletes atoms outside the box.
@param 	*molgroup 	molecule group structure to be modified.
@param 	box			box size.
@return long		number of remaining atoms.
**/
long		molgroup_trim(Bmolgroup* molgroup, Vector3<double> box)
{
	if ( box.volume() < 1 ) return 0;
	
	long			n(0), nsel(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Vector3<double>	start, end(box);
	
	if ( verbose & VERB_PROCESS )
		cout << "Trimming to a box of size:      " << box << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next, n++ ) {
				if ( atom->coord.within(start, end) ) {
					atom->sel = 1;
					nsel++;
				} else {
					atom->sel = 0;
				}
			}
		}
    }

	long	ndel = molgroup_delete_deselected_atoms(molgroup);
	
	if ( n != nsel + ndel )
		cerr << "Error: There is a discrepancy between the the total number of atoms ("
			<< n << ") and the sum of selected (" << nsel
			<< ") and deleted (" << ndel << ") atoms!" << endl;
	
	molgroup->box = box;

	if ( verbose & VERB_PROCESS )
		cout << "Number of atoms remaining:      " << nsel << endl;
	
	return nsel;
}
