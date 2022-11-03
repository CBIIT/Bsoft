/**
@file	bsf.cpp
@brief	Calculating structure factors from atomic models
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
#include "ctf.h"
#include "scatter.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bsf [options] out.pdb out.hkl",
"------------------------------------",
"Calculates structure factors from atomic models.",
"A coordinate file is required.",
"Note: The occupancy for each atom is used - make sure it is set properly.",
"If a structure factor file is given, it is compared to the calculated structure factors.",
" ",
"Selections:",
"-select CA               Atom selection (default all).",
"-chains B,G,O,Z          Select chains or molecules by ID (default all).",
"-occupancy 27,99,0.6     Set an occupancy range (res1,res2) to a value.",
" ",
"Actions:",
"-center                  Center coordinates before calculations.",
"-realspace               Transform the output back to real space.",
"-scatter 300             Correct scattering intensities for the given acceleration voltage.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-size 100,80,70          Size (default automatic voxels).",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-sampling 2.5,2.5,2.5    Sampling (default 1 angstrom/voxel).",
"-resolution 3            Resolution (default 0.1 angstrom).",
"-wrap                    Wrap around periodic boundaries (default off).",
"-Bfactor 30              Global B-factor (default 0 angstrom squared).",
"-symmetry 1              Space group (default 1).",
"-unitcell 50,50,50,90,90,90 Unit cell parameters (angstrom & degrees).",
"-curves C,Au,Mn,H        Selection of elements for scattering curve output.",
" ",
"Input:",
"-coordinates file.pdb    Input coordinate file.",
//"-hkl file.hkl            Input structure factor file.",
"-parameters scat.star    Parameter file with scattering coefficients.",
" ",
"Output:",
"-output outfile.txt      File for scattering curve output (use with -curves option).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    /* Initialize variables */
    Bstring    		atom_select("ALL");
	Bstring			chain_select;
    double			occupancy(0);
    int 			range_first(0);			// For resetting occupancy
    int 			range_last(0);
    int 			set_occupancy(0);		// Flag for setting occupancy
	int				center(0);				// Flag to center coordinates
	int 			set_backtransform(0);	// Flag for back transformation
	double			volt(0);				// In volts
	Vector3<double>	origin;					// Coordinate origin placement
	int				set_origin(0);			// Flag to set origin
	Vector3<long>	size;					// New map size
	Vector3<double>	sam;    				// Sampling in angstrom/voxel side
	double			resolution(0.1);		// Resolution = 0.1 angstrom
	int 			wrap(0);				// No wrapping as default
	double			Bfactor(0);				// B-factor = 0 angstrom squared
	int 			spacegroup(1);
	UnitCell		uc(0,0,0,M_PI_2,M_PI_2,M_PI_2);
	Bstring			curves;					// String with elements for curve output
	Bstring			coorfile;				// Atomic coordinates
	Bstring			mapfile;				// Map to compare with
	Bstring			paramfile;				// Use default parameter file
	Bstring			outfile;				// File name for curve output

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
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
		if ( curropt->tag == "occupancy" ) {
        	if ( curropt->values(range_first, range_last, occupancy) < 3 )
				cerr << "-occupancy: A range and value must be specified!" << endl;
			else
				set_occupancy = 1;
		}
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "realspace" )
			set_backtransform = 1;
		if ( curropt->tag == "scatter" ) {
			if ( ( volt = curropt->value.real() ) < 1 )
				cerr << "-scatter: A voltage must be specified!" << endl;
			else
				if ( volt < 1e3 ) volt *= 1e3;	// Assume kilovolts
		}
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
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
		if ( curropt->tag == "curves" )
			curves = curropt->value;
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->value;
		if ( curropt->tag == "map" )
        	mapfile = curropt->filename();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
		
	if ( outfile.length() )
		write_scattering_curves(paramfile, outfile, curves, resolution);
	
	Bmolgroup*		molgroup = NULL;
	Bimage* 		p = NULL;
	Bstring			ext;
	Bstring*		file_list = NULL;
	Vector3<double> box;
    
	if ( coorfile.length() ) {
//		molgroup = read_molecule(coorfile, atom_select, paramfile);
		file_list = coorfile.split(",");
		molgroup = read_molecule(file_list, 0, box, atom_select, paramfile);
		
		if ( !molgroup->mol || !molgroup->mol->res ) {
			cerr << "Error: Problem with coordinate file " << coorfile << ", exiting!" << endl;
			bexit(-1);
		}
	}
	
	if ( molgroup == NULL ) {
		cerr << "No coordinates read!" << endl;
		bexit(-1);
	}

    if ( set_occupancy )
		molgroup_coor_reset_occupancy(molgroup, range_first,
				range_last, occupancy);
    
	if ( center ) 
		molgroup_shift_to_center_of_mass(molgroup);
	
	if ( mapfile.length() ) {
		p = read_img(mapfile, 1, -1);
		if ( !p ) {
			cerr << "Error: Map file " << mapfile << " not read!" << endl;
			bexit(-1);
		}
	}
	
	// Get default parameters from the input image if possible
	if ( p ) {
		if ( set_origin ) {
			if ( set_origin == 2 ) origin = p->default_origin();
		} else origin = p->image->origin();
		if ( sam.volume() > 0 ) p->sampling(sam);
		if ( size[0] <= 0 || size[1] <= 0 || size[2] <= 0 ) {
			size = {p->sizeX(), p->sizeY(), p->sizeZ()};
			if ( !uc.a() ) uc[0] = p->real_size()[0];
			if ( !uc.b() ) uc[1] = p->real_size()[1];
			if ( !uc.c() ) uc[2] = p->real_size()[2];
		}
	}

	if ( set_origin == 2 ) origin = {double(size[0]/2), double(size[1]/2), double(size[2]/2)};
	
	if ( chain_select.length() )
		molgroup_select_chains(molgroup, chain_select);

	// Calculate the structure factors
	Bimage* 	pcalc = img_sf_from_molecule(molgroup, origin, size, 
					sam, resolution, spacegroup, uc, wrap, Bfactor, paramfile);
	
	double		mw = molgroup_weight_from_atoms(molgroup);

#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	
	if ( pcalc ) {
		if ( volt ) {
			double		lambda = electron_wavelength(volt);
			double		b2 = beta2(volt);
			double		scale = 2*lambda/sqrt(1 - b2);
			if ( verbose ) {
				cout << "Scaling to correct scattering amplitudes:" << endl;
				cout << "Voltage:                        " << volt*1e-3 << " kV" << endl;
				cout << "Wavelength:                     " << lambda << " A" << endl;
				cout << "Beta:                           " << sqrt(b2) << endl;
				cout << "Scale:                          " << scale << endl;
			}
//			scale /= sqrt(POTPREFAC);	// Remove old scaling factor
			scale /= POTPREFAC;	// Remove old scaling factor
			pcalc->multiply(scale);
		}
		if ( set_backtransform ) {
			pcalc->fft(FFTW_BACKWARD, 0);
			pcalc->complex_to_real();
			if ( pcalc->sizeZ() > 1 ) pcalc->mass_threshold(0, mw, RHO);
		}
	}
	
#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif

    // Write output file(s)
	while ( optind < argc ) {
		ext = argv[optind];
		ext = ext.extension();
		if ( ext.contains("pdb") || ext.contains("gro") ) {
    		if ( molgroup ) write_molecule(argv[optind++], molgroup);
		} else if ( pcalc ) {
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

