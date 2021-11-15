/**
@file	bmonte.cpp
@brief	Program to use a monte carlo metroplis algorithm to energy minimize molecular positions.
@author Bernard Heymann
@date	Created: 20041230
@date 	Modified: 20071223
**/

#include "rwmolecule.h"
#include "rwimg.h"
#include "rwmd.h"
#include "mol_monte.h"
#include "mol_bonds.h"
#include "mol_transform.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmonte [options] input.pdb output.pdb",
"--------------------------------------------",
"Minimizes the energy of a set of molecular structures using a Monte Carlo Metropolis algorithm.",
"The -grid and -orientation options only work with rigid bodies, not atoms,",
"	and their output is to a model parameter file (.star).",
" ",
"Actions:",
"-rigid mol               Rigidity: all/whole = single rigid body (default),",
"                                   group = each file with a molecule group is a rigid body,",
"                                   molecule = each molecule is a rigid body,",
"                                   atom = each atom can move.",
"-concurrent              Run all models simultaneously (default each model separately),",
"-grid 3.5,5,4.2          Search on a grid with this sampling and within the bounding box (angstroms).",
"-orientations 12.5       Search all orientations with the given angular step size (degrees).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-bbox 0,82.5,50,50,80,50 Bounding box center and size (default from coordinates or map).",
"-wrap                    Wrap around (periodic boundaries).",
"-resolution 20           Resolution limit (angstrom).",
"-iterations 60           Maximum number of iterations (default 1).",
"-Kbond 150               Bond strength (default 0).",
"-Kangle 4                Angle strength (default 0).",
"-Kvdw 0.1                Van der Waals strength (default 0).",
"-Kelectrostatic 0.4      Electrostatic strength (default 0).",
"-Kseparation 0.02        Separation energy constant (default 0).",
"-Kmap 1.5                Map energy constant (default 1 if -Map option is used).",
"-separation 4.5          Separation distance for overlap calculation (default 4 A).",
"-cutoff 7.8              Distance cutoff for non-bonded forces (default 5 A).",
"-beta 12                 Inverse of mean energy per atom and degree of freedom (default 10).",
"-angle 1.5               Maximum angular increment per iteration (default 1 degree).",
"-shift 3.1               Maximum shift per iteration (default 1 angstrom).",
"-bondsteps 5             Number of steps along a bond for map fitting (default none).",
"-location 12,5.3,6       Location of harmonic force to apply to center-of-mass.",
"-Klocation 16.8,0.2      Magnitude of harmonic force and decay constant.",
"-SS 2.05                 Set disulphide reference bond length.",
" ",
"Input:",
"-parameters md.star      Molecular dynamics parameter file.",
"-Model locations.star    Model parameter file with initial rigid body locations.",
"-Map file.map            Map to use as an additional restraint.",
"-Mask mask.mrc           Mask to limit grid-searches (must be byte data type).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	Vector3<double>	sam;    			// Map sampling
	Vector3<double>	grid_sampling;				// Grid sampling
	Vector3<double>	bbox_center;				// Bounding box center
	Vector3<double>	bbox_size;					// Bounding box size
	double			angle_step(0);				// Angle step size for orientation search
	int 			wrap(0);					// No wrapping as default
	double	 			resolution(0); 				// Must be set > 0 to limit resolution
	unsigned long	max_iter(1);				// Maximum number of iterations
	double			Kbond(0);					// Bond strength
	double			Kangle(0);					// Angle strength
	double			Kelec(0);					// Electrostatic strength
	double			Kvdw(0);					// Van der Waals strength
	double			Ksep(0);					// Separation energy constant
	double			Kmap(0);					// Map energy constant
	double			separation(4);				// Separation distance for overlap calculation
	double			cutoff(5);					// Distance cutoff for non-bonded forces
	double			beta(10);					// Equivalent of 1/kT
	double			max_angle(M_PI/180.0);		// Maximum angular deviation
	double			max_shift(1);				// Maximum allowed shift
	int				bond_steps(0);				// Number of steps along bond
	double			Kpoint(0);					// Harmonic location force constant
	double			pointdecay(0.1);			// Point force decay constant
	Vector3<double>	location;					// Point for harmonic force
	int				rigid(0);					// Flag to treat the whole ensemble as a rigid body
	int				concurrent(0);				// Flag to run all models simultaneously
	double			ss(0);						// SS bond length
	Bstring			atom_select("all");			// Selection
	Bstring			modelfile;					// Input model file
	Bstring			mapfile;					// Input map file
	Bstring			maskfile;					// Input map file
	Bstring			paramfile;					// Parameter file
    
	random_seed();
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rigid" ) {
			if ( curropt->value[0] == 'g' ) rigid = 1;
			if ( curropt->value[0] == 'm' ) rigid = 2;
			if ( curropt->value[0] == 'a' ) rigid = 3;
		}
		if ( curropt->tag == "concurrent" ) concurrent = 1;
 		if ( curropt->tag == "grid" )
        	grid_sampling = curropt->scale();
 		if ( curropt->tag == "orientations" ) {
			if ( ( angle_step = curropt->value.real() ) < 0.1 )
				cerr << "-orientations: An angle step size must be specified!" << endl;
			else
				angle_step *= M_PI/180.0;
		}
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
		if ( curropt->tag == "bbox" )
			if ( curropt->box(bbox_center, bbox_size) < 6 )
				cerr << "-bbox: All 6 values must be specified!" << endl;
		if ( curropt->tag == "wrap" )
			wrap = 1;
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "iterations" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number must be specified!" << endl;
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
		if ( curropt->tag == "Kseparation" )
			if ( ( Ksep = curropt->value.real() ) < 1e-30 )
				cerr << "-Kseparation: The separation energy constant strength must be specified!" << endl;
		if ( curropt->tag == "Kmap" )
			if ( ( Kmap = curropt->value.real() ) < 1e-30 )
				cerr << "-Kmap: The map energy constant strength must be specified!" << endl;
		if ( curropt->tag == "separation" )
			if ( ( separation = curropt->value.real() ) < 1e-30 )
				cerr << "-separation: The separation distance must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( cutoff = curropt->value.real() ) < 1e-30 )
				cerr << "-cutoff: The cutoff distance must be specified!" << endl;
		if ( curropt->tag == "beta" )
			if ( ( beta = curropt->value.real() ) < 1e-30 )
				cerr << "-beta: A value must be specified!" << endl;
		if ( curropt->tag == "angle" ) {
			if ( ( max_angle = curropt->value.real() ) < 1e-30 )
				cerr << "-angle: An angle must be specified!" << endl;
			else
				max_angle = angle_set_negPI_to_PI(max_angle*M_PI/180.0);
		}
		if ( curropt->tag == "shift" )
			if ( ( max_shift = curropt->value.real() ) < 1e-30 )
				cerr << "-shift: A distance must be specified!" << endl;
		if ( curropt->tag == "bondsteps" )
			if ( ( bond_steps = curropt->value.integer() ) < 1 )
				cerr << "-bondsteps: The bond step length must be specified!" << endl;
		if ( curropt->tag == "location" )
			location = curropt->vector3();
		if ( curropt->tag == "Klocation" )
			if ( curropt->values(Kpoint, pointdecay) < 1 )
				cerr << "-Klocation: The location force constant must be specified!" << endl;
		if ( curropt->tag == "SS" )
			if ( ( ss = curropt->value.real() ) < 1e-30 )
				cerr << "-SS: A bond length must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "Model" )
			modelfile = curropt->filename();
		if ( curropt->tag == "Map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();
	
	Bmd*		md = NULL;
	if ( paramfile.length() ) {
		md = read_md_parameters(paramfile);
	} else {
		md = md_init();
	}
	md->Kbond = Kbond;
	md->Kangle = Kangle;
	md->Kelec = Kelec;
	md->Kvdw = Kvdw;
	md->Ksep = Ksep;
	md->Kpoint = Kpoint;
	md->point = location;
	md->pointdecay = pointdecay;
	md->sepdist = separation;
	md->bondsteps = bond_steps;
	md->cutoff = cutoff;
	md->wrap = wrap;
	
    // Read the molecule file
	Bstring		filename(argv[optind++]);
    Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: No molecules read!" << endl;
		bexit(-1);
	}

	molecule_update_comment(molgroup, argc, argv);
	
	if ( rigid == 3 ) {
		md_generate_bond_list(molgroup, md);
		md_bond_list_set_parameters(molgroup->bond, md->bond);
//		md_generate_angle_list(molgroup, md);
//		md_angle_list_set_parameters(molgroup->angle, md->angle);
	}
	
	if ( ss ) {
		Bbondtype*		bt = new Bbondtype;
		strcpy(bt->type1, "S");
		strcpy(bt->type2, "S");
		bt->covlength = ss;
		md_bond_list_set_parameters(molgroup->bond, bt);
		delete bt;
	}
	
	// Read the map file
	long		natom_in(0);
	Bimage*		map = NULL;
	Bimage*		pmask = NULL;
	Bmodel*		model = NULL;
	
	filename = 0;
	if ( optind < argc )
		filename = argv[optind];
	
	if ( mapfile.length() ) {
		map = read_img(mapfile, 1, 0);
		if ( sam.volume() ) map->sampling(sam);
		molgroup_set_box_to_map_boundaries(molgroup, map);
		natom_in = molgroup_test_if_within_box(molgroup, molgroup->min, molgroup->max);
		if ( natom_in < 2 ) {
			molgroup_place_at_coordinates(molgroup, ((molgroup->max + molgroup->min) * 0.5));
			molgroup_set_box_to_map_boundaries(molgroup, map);
		}
		md->Kmap = 1;
		if ( Kmap ) md->Kmap = Kmap;
	}

	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		if ( sam.volume() ) pmask->sampling(sam);
	}
	
	if ( bbox_size[0] > 0 && bbox_size[1] > 0 && bbox_size[2] > 0 ) {
		molgroup->min = bbox_center - (bbox_size * 0.5);
		molgroup->max = bbox_center + (bbox_size * 0.5);
		molgroup->box = bbox_size;
	}
	
	if ( modelfile.length() ) {
		model = read_model(&modelfile);
	} else if ( grid_sampling[0] > 0 ) {
		model = molgroup_generate_masked_grid_list(molgroup, grid_sampling, pmask);
	} else if ( angle_step > 0 ) {
		model = molgroup_generate_orientation_list(molgroup, angle_step);
	}
	
	if ( model ) {
		if ( concurrent ) 
			mcm_molecule_groups(molgroup, model, md, map, beta, max_angle, max_shift, max_iter, rigid);
		else
			mcm_molecule_list(molgroup, model, md, map, beta, max_angle, max_shift, max_iter, rigid);
		if ( filename.length() )
			write_model(filename, model);
	} else {
		if ( rigid < 3 ) {
			molgroup = monte_carlo_metropolis(molgroup, md, map, beta, max_angle, max_shift, 
				max_iter, rigid, monte_rigid_body_fit_energy, molgroup_rigid_body_transform);
		} else {
			if ( bond_steps > 0 )
				molgroup = monte_carlo_metropolis(molgroup, md, map, beta, max_angle, max_shift, 
					max_iter, rigid, monte_bond_fit_energy, molgroup_move_atoms_down_energy);
			else
				molgroup = monte_carlo_metropolis(molgroup, md, map, beta, max_angle, max_shift, 
					max_iter, rigid, monte_atom_fit_energy, molgroup_move_atoms_down_energy);
		}
		if ( filename.length() )
			write_molecule(filename, molgroup);
	}
	
	delete map;
	delete pmask;
	molgroup_list_kill(molgroup);
	model_kill(model);
	md_kill(md);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

