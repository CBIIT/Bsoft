/**
@file	bmd.cpp
@brief	Molecular dynamics - humble beginnings
@author Bernard Heymann
@date	Created: 20001014
@date 	Modified: 20060412
**/

#include "rwmd.h"
#include "mol_md.h"
#include "rwmolecule.h"
#include "mol_bonds.h"
#include "mol_util.h"
#include "Matrix.h"
#include "options.h"
#include "linked_list.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmd [options] in.pdb out.pdb",
"-----------------------------------",
"Calculates forces on an atomic structure and adjust coordinates.",
" ",
"Actions:",
"-dynamics 4000           Number of iterations for molecular dynamics.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
//"-origin 0,0,0            Origin (angstrom).",
"-size 10,10,10           Simulation box size (angstrom).",
"-bonds intra             Generate all or intramolecular bonds only.",
//"-length 1.6              Bond length (angstrom).",
"-wrap                    Wrap around (periodic boundaries).",
"-timestep 0.01           Integration time step (default 0.001).",
"-velocitylimit 0.01      Limit on the velocity (default 0.1 per time step).",
"-friction 0.2            Friction constant (default 1 = no friction).",
"-Kbond 150               Bond strength (default 1).",
"-Kangle 4                Angle strength (default 1).",
"-Kvdw 0.1                Van der Waals strength (default 0).",
"-Kelectrostatic 0.4      Electrostatic strength (default 0).",
"-cutoff 7.8              Distance cutoff for non-bonded forces (default 5 A).",
" ",
"Input:",
"-parameters md.star      Molecular dynamics parameter file (default md_param.star).",
" ",
"Output:",
"-output param.star       Output parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    /* Initialize variables */
//	int				calcrdf(0);				// Flag to calculate RDF
	int 			wrap(0);				// No wrapping as default
	Vector3<double>	size;					// Size for generating water
	int				bondtype(0);			// Bond types to consider: 0=all, 1=intramolecular only
	double			bondlength(1);			// Bond length
	int 			max_iter(0);			// Number of iterations/cycles for minimization
	double			timestep(0.001);		// Integration time step
	double			velocitylimit(0.1);		// Limit on velocity per time step
	double			Kfriction(1);			// Friction constant, 1=no friction
	double			Kbond(1);				// Bond strength
	double			Kangle(1);				// Angle strength
	double			Kelec(0);				// Electrostatic strength
	double			Kvdw(0);				// Van der Waals strength
	double			cutoff(5);				// Distance cutoff for non-bonded forces
    Bstring    		atom_select("all");
	Bstring			paramfile("md_param.star");	// Default parameter file
	Bstring			paramout;				// Output parameter file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "size" ) {
			size = curropt->size();
			if ( size.volume() < 3 )
				cerr << "-size: All three dimensions must be specified" << endl;
		}
		if ( curropt->tag == "wrap" ) wrap = 1;
		if ( curropt->tag == "bonds" ) {
			if ( curropt->value.contains("int") ) bondtype = 1;
			if ( curropt->value.contains("val") ) bondtype = 2;
		}
		if ( curropt->tag == "length" )
			if ( ( bondlength = curropt->value.real() ) < 1e-6 )
				cerr << "-length: A bond length must be specified!" << endl;
		if ( curropt->tag == "dynamics" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-dynamics: The number of iterations must be specified!" << endl;
		if ( curropt->tag == "timestep" )
			if ( ( timestep = curropt->value.real() ) < 1e-6 )
				cerr << "-timestep: The time step must be specified!" << endl;
		if ( curropt->tag == "velocitylimit" )
			if ( ( velocitylimit = curropt->value.real() ) < 1e-6 )
				cerr << "-velocitylimit: The velocity limit must be specified!" << endl;
		if ( curropt->tag == "friction" )
			if ( ( Kfriction = curropt->value.real() ) < 1e-6 )
				cerr << "-friction: The friction constant must be specified!" << endl;
		if ( curropt->tag == "Kbond" )
			if ( ( Kbond = curropt->value.real() ) < 1e-6 )
				cerr << "-Kbond: The bond strength must be specified!" << endl;
		if ( curropt->tag == "Kangle" )
			if ( ( Kangle = curropt->value.real() ) < 1e-6 )
				cerr << "-Kangle: The angle strength must be specified!" << endl;
		if ( curropt->tag == "Kelectrostatic" )
			if ( ( Kelec = curropt->value.real() ) < 1e-6 )
				cerr << "-Kelectrostatic: The electrostatic strength must be specified!" << endl;
		if ( curropt->tag == "Kvdw" )
			if ( ( Kvdw = curropt->value.real() ) < 1e-6 )
				cerr << "-Kvdw: The Van der Waals strength must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( cutoff = curropt->value.real() ) < 1e-6 )
				cerr << "-cutoff: The cutoff distance must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			paramout = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bstring		filename(argv[optind++]);
	Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: No molecule to work with!" << endl;
		bexit(-1);
	}
	
	if ( size[0] > 0 && size[1] > 0 && size[2] > 0 )
		molgroup->box = size;
	else
		molgroup->box = molgroup->max - molgroup->min;
	
	Bmd*		md = read_md_parameters(paramfile);
	if ( !md ) {
		cerr << "Error: File " << paramfile << " not read!" << endl;
		bexit(-1);
	}
	
	md->timestep = timestep;
	md->Kfriction = Kfriction;
	md->Kbond = Kbond;
	md->Kangle = Kangle;
	md->Kelec = Kelec;
	md->Kvdw = Kvdw;
	md->cutoff = cutoff;
	md->wrap = wrap;

	md_show_bonds(molgroup);
	
	Bbond*		bondlist = NULL;
	if ( bondtype == 1 )
		bondlist = md_generate_molecular_bond_list(molgroup, md);
	else if ( bondtype == 2 )
		bondlist = md_generate_bond_list_with_valence(molgroup, md, 3);
	else
		bondlist = md_generate_bond_list(molgroup, md);

	md_bond_list_set_parameters(bondlist, md->bond);

	md_show_bond_stats(molgroup);

//	md_generate_angle_list(molgroup, md);

	if ( verbose & VERB_PROCESS )
		molgroup_density(molgroup);
	
	if ( verbose ) {
		md_calculate_deviations(molgroup, md->wrap);
		md_calculate_radial_deviation(molgroup);		
		md_show_bond_types(molgroup, md->bond);
	}
	
	if ( molgroup->mol && max_iter ) {
		md_leapfrog(molgroup, md, max_iter, velocitylimit);
		if ( verbose & VERB_PROCESS )
			molgroup_density(molgroup);
		if ( verbose ) {
			md_calculate_deviations(molgroup, md->wrap);
			md_calculate_radial_deviation(molgroup);
			md_show_bond_types(molgroup, md->bond);
		}
	}
	
	if ( paramout.length() )
		write_md_parameters(paramout, md);
	
	if ( optind < argc ) {
		molecule_update_comment(molgroup, argc, argv);
		write_molecule(argv[optind], molgroup);
	}

    molgroup_kill(molgroup);
	md_kill(md);
		
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

