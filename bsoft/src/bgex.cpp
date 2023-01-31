/**
@file	bgex.cpp
@brief	Calculating and comparing atomic models and maps
@author Bernard Heymann
@date	Created: 19970914
@date 	Modified: 20210528
**/

#include "molecule_to_map.h"
#include "mol_transform.h"
#include "mol_select.h"
#include "mol_util.h"
#include "rwmolecule.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bgex [options] out.pdb out.map",
"-------------------------------------",
"Calculates and compares atomic models and maps.",
"A coordinate file is required.",
"If a map file is given, the occupancy of the atoms are calculated.",
" ",
"Selections:",
"-select CA               Atom selection (default all).",
"-chains B,G,O,Z          Select chains or molecules by ID (default all).",
"-limits -50.5,10,0,40,-100,-20 Coordinate limits (default -10000,-10000,-10000,10000,10000,10000).",
" ",
"Actions:",
"-center                  Center coordinates before calculations.",
"-potential               Calculate an atomic potential map (default simple gaussian expansion).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new image data type.",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-size 10,10,10           Map size, input map size otherwise (voxels).",
"-sampling 1              Sampling (default 1 angstrom/voxel).",
"-resolution 3            Resolution (default 2 angstrom, only simple gaussian expansion).",
"-wrap                    Wrap around periodic boundaries (default off).",
"-Bfactor 30              Global B-factor (default 0 angstrom squared).",
"-symmetry 1              Space group (default 1).",
"-unitcell 50,50,50,90,90,90 Unit cell parameters (angstrom & degrees).",
" ",
"Input:",
"-coordinates file.pdb    Input coordinate file.",
"-map file.map            Input map file.",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
    Bstring    		atom_select("ALL");
	Bstring			chain_select;
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	origin;					// Coordinate origin placement
	int				set_origin(0);			// Flag to set origin
	Vector3<long>	size;					// New map size
	Vector3<double>	sam;    				// Sampling in angstrom/voxel side
	double			resolution(2); 			// Resolution = 2 angstrom
	int 			wrap(0);				// No wrapping as default
	int				center(0);				// Flag to center coordinates
	int				gextype(0);				// Type of gaussian used: 0 = single, 1 = atomic potential
	double			Bfactor(0);				// B-factor = 30 angstrom squared
	int 			spacegroup(1);
	UnitCell		uc;
	Bstring			coorfile;				// Atomic coordinates
	Bstring			mapfile;				// Map to compare with
	Bstring			paramfile;				// Use default parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*			curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
			if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "chains" ) {
			chain_select = curropt->value;
			if ( chain_select.length() < 1 )
				cerr << "-chains: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "potential" )
			gextype = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->filename();
		if ( curropt->tag == "map" )
        	mapfile = curropt->filename();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "size" )
        	size = curropt->size();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
        	if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: Resolution must be specified!" << endl;
		if ( curropt->tag == "wrap" )
        	wrap = 1;
		if ( curropt->tag == "Bfactor" )
        	if ( ( Bfactor = curropt->value.real() ) < 0.001 )
				cerr << "-Bfactor: The B-factor must be specified!" << endl;
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

	Bmolgroup*		molgroup = NULL;
	Bimage* 		p = NULL;
	Bstring			ext;
    
	if ( coorfile.length() ) {
		molgroup = read_molecule(coorfile, atom_select, paramfile);
		if ( !molgroup->mol->res ) {
			cerr << "Error: Problem with coordinate file " << coorfile << ", exiting!" << endl;
			bexit(-1);
		}
	}
	
	if ( mapfile.length() ) {
		p = read_img(mapfile, 1, -1);
		if ( !p ) {
			cerr << "Error: Map file " << mapfile << " not read!" << endl;
			bexit(-1);
		}
	}
	
	if ( molgroup == NULL ) {
		cerr << "Error: No coordinate file given!" << endl;
		bexit(-1);
	}
	
	// Set some parameters from input
	if ( p ) {
		if ( set_origin ) {
			if ( set_origin == 2 ) origin = p->default_origin();
		} else origin = p->image->origin();
		if ( sam.volume() > 0 ) p->sampling(sam);
		if ( size[0] <= 0 || size[1] <= 0 || size[2] <= 0 ) {
			size = {p->sizeX(), p->sizeY(), p->sizeZ()};
			if ( !uc.a() ) uc[0] = p->sizeX()*p->sampling(0)[0];
			if ( !uc.b() ) uc[1] = p->sizeY()*p->sampling(0)[1];
			if ( !uc.c() ) uc[2] = p->sizeZ()*p->sampling(0)[2];
		}
	}

	if ( set_origin == 2 ) origin = {double(size[0]/2), double(size[1]/2), double(size[2]/2)};
		
	if ( center )
		molgroup_shift_to_center_of_mass(molgroup);

	if ( chain_select.length() )
		molgroup_select_chains(molgroup, chain_select);

	// Calculate the density
	Bimage* 	pcalc = img_from_molecule(molgroup, origin, 
					size, sam, resolution, Bfactor, wrap, gextype, spacegroup, uc);

	double		mw = molgroup_weight_from_atoms(molgroup);

	if ( pcalc )
		pcalc->mass_threshold(0, mw, RHO);
	
	// If a map is also given, compare it to the generated density
	if ( p )
		compare_mol_map(molgroup, pcalc, p);
	
    // Write output file(s)
	while ( optind < argc ) {
		ext = argv[optind];
		ext = ext.extension();
		if ( ext.contains("pdb") || ext.contains("gro") ) {
    		if ( molgroup ) write_molecule(argv[optind++], molgroup);
		} else if ( pcalc ) {
			pcalc->change_type(nudatatype);
			write_img(argv[optind++], pcalc, 0);
		}
	}

    if ( molgroup ) molgroup_kill(molgroup);
	if ( p ) delete p;
	if ( pcalc ) delete pcalc;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

