/**
@file	bwater.cpp
@brief	Molecular dynamics - humble beginnings
@author Bernard Heymann
@date	Created: 20001014
@date 	Modified: 20060412
**/

#include "rwmd.h"
#include "mol_md.h"
#include "rwmolecule.h"
#include "mol_water.h"
#include "mol_bonds.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bwater [options] in.pdb out.pdb",
"--------------------------------------",
"Generates and does molecular dynamics on blocks of water.",
" ",
"Actions:",
"-create rect             Generate water, 1=random, 2=rectangular, 3=tetrahedral (use with -size option).",
"-rdf 0.05                Calculate the radial distribution function with the given sampling (angstrom).",
"-dynamics 4000           Number of iterations for molecular dynamics.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
//"-origin 0,0,0            Origin (angstrom).",
"-size 10,10,10           Size for generating water (angstrom).",
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
	int				genwater(0);			// Flag to generate water
	double			rdf(0);					// Calculate RDF at this sampling (0=not)
	int 			wrap(0);				// No wrapping as default
//	int 			set_origin(0); 			// Flag to indicate if origin is set
	Vector3<double>	size;					// Size for generating water
//	Vector3<double>	origin;					// Coordinate origin in angstrom
	int 			max_iter(0);			// Number of iterations/cycles for minimization
	double			timestep(0.001);		// Integration time step
	double			velocitylimit(0.1);		// Limit on velocity per time step
	double			Kfriction(1);			// Friction constant, 1=no friction
	double			Kbond(1);				// Bond strength
	double			Kangle(1);				// Angle strength
	double			Kelec(0);				// Electrostatic strength
	double			Kvdw(0);				// Van der Waals strength
	double			cutoff(5);				// Distance cutoff for non-bonded forces and rdf
    Bstring    		atom_select("all");
	Bstring			paramfile("water_param.star");	// Default parameter file
	Bstring			paramout;				// Output parameter file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" ) {
			if ( curropt->value.contains("ran") ) genwater = 1;
			else if ( curropt->value.contains("rec") ) genwater = 2;
			else if ( curropt->value.contains("tet") ) genwater = 3;
			else genwater = curropt->value.integer();
		}
		if ( curropt->tag == "rdf" )
			if ( ( rdf = curropt->value.real() ) < 0.01 )
				cerr << "-rdf: A sampling interval must be specified!" << endl;
		if ( curropt->tag == "size" ) {
			size = curropt->vector3();
			if ( size.volume() < 1 )
				cerr << "-size: All three dimensions must be specified" << endl;
		}
//		if ( curropt->tag == "origin" ) {
//			origin = curropt->origin();
//			set_origin = 1;
//		}
		if ( curropt->tag == "wrap" )
			wrap = 1;
		if ( curropt->tag == "dynamics" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-dynamics: The number of iterations must be specified!" << endl;
		if ( curropt->tag == "timestep" )
			if ( ( timestep = curropt->value.real() ) < 1e-30 )
				cerr << "-timestep: The time step must be specified!" << endl;
		if ( curropt->tag == "velocitylimit" )
			if ( ( velocitylimit = curropt->value.real() ) < 1e-30 )
				cerr << "-velocitylimit: The velocity limit must be specified!" << endl;
		if ( curropt->tag == "friction" )
			if ( ( Kfriction = curropt->value.real() ) < 1e-30 )
				cerr << "-friction: The friction constant must be specified!" << endl;
		if ( curropt->tag == "Kbond" )
			if ( ( Kbond = curropt->value.real() ) < 1e-30 )
				cerr << "-Kbond: The bond strength must be specified!" << endl;
		if ( curropt->tag == "Kangle" )
			if ( ( Kangle = curropt->value.real() ) < 1e-30 )
				cerr << "-Kangle: The angle strength must be specified!" << endl;
		if ( curropt->tag == "Kelectrostatic" )
			if ( ( Kelec = curropt->value.real() ) < 1e-30 )
				cerr << "-Kelectrostatic: The electrostatic strength must be specified!" << endl;
		if ( curropt->tag == "Kvdw" )
			if ( ( Kvdw = curropt->value.real() ) < 1e-30 )
				cerr << "-Kvdw: The Van der Waals strength must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( cutoff = curropt->value.real() ) < 1e-30 )
				cerr << "-cutoff: The cutoff distance must be specified!" << endl;
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			paramout = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bmolgroup*	molgroup = NULL;
	Bstring		filename(argv[optind++]);
	
	if ( genwater == 1 )
		molgroup = molgroup_generate_random_water(size);
	else if ( genwater > 1 )
		molgroup = molgroup_generate_regular_water(size, genwater);
	else {
		molgroup = read_molecule(filename, atom_select, paramfile);
		if ( size[0] > 0 && size[1] > 0 && size[2] > 0 )
			molgroup->box = size;
		else
			molgroup->box = molgroup->max - molgroup->min;
	}
	
	if ( !molgroup ) {
		cerr << "Error: No molecule to work with!" << endl;
		bexit(-1);
	}
	
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
	
	water_bond_list(molgroup);
	
	water_angle_list(molgroup);
	
	if ( verbose & VERB_PROCESS )
		molgroup_density(molgroup);
	
	if ( verbose )
		md_calculate_deviations(molgroup, md->wrap);

	if ( molgroup->mol && max_iter ) {
		md_leapfrog(molgroup, md, max_iter, velocitylimit);
		if ( verbose & VERB_PROCESS )
			molgroup_density(molgroup);
		if ( verbose )
			md_calculate_deviations(molgroup, md->wrap);
	}
	
	if ( molgroup->mol && rdf > 0 ) molgroup_calc_water_rdf(molgroup, rdf, cutoff);
	
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

