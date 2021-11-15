/**
@file	bmol.cpp
@brief	A tool to edit and analyze coordinate files
@author Bernard Heymann
@date	Created: 19980214
@date 	Modified: 20211029
**/

#include "rwmolecule.h"
#include "mol_transform.h"
#include "mol_edit.h"
#include "mol_select.h"
#include "mol_compare.h"
#include "mol_util.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmol [options] in.pdb out.pdb",
"------------------------------------",
"Edits, analyzes and transforms atomic coordinate files.",
" ",
"Selections:",
"-select CA               Atom selection (default all).",
"-random 100              Select a number of atoms randomly.",
"-chains B,G,O,Z          Select chains or molecules by ID (default all).",
"-limits -50.5,10,0,40,-100,-20 Coordinate limits (default -10000,-10000,-10000,10000,10000,10000).",
" ",
"Actions:",
"-invert 10,-5.3,2.6      Invert through a point.",
"-translate 0,-50,22      Translate (angstrom).",
"-scale 2,1.5,2.2         Scale the coordinates.",
"-rotate 0.5,0,5.8,45     Rotate around a vector (x,y,z) by an angle.",
"-toview 0.3,0.5,0.8,33   Rotate to view: vector {xyz} and angle (default 0,0,1,0).",
"-fromview 0.3,0.5,0.8,33 Rotate from view: vector {xyz} and angle (default 0,0,1,0).",
"-radius 41.8             Force all coordinates to this radius (default not).",
"-rename D                Rename molecules from the given letter.",
"-renumber 52             Renumber residues from the given number.",
"-occupancy 27,99,0.6     Set an occupancy range (res1,res2) to a value.",
"-box 15,66,54            Consolidate structure within box (0,0,0 means no change).",
"-wrap 15,66,54           Pack coordinates within periodic box.",
"-Bfactor                 Analyze B-factors.",
"-Sequence                Output sequence as string of letters.",
"-Mass                    Calculate the total molecular weight and volume.",
"-Volume 1                Calculate the Van der Waals volume, wrapping flag (default 0).",
"-Composition             Calculate the composition.",
"-Center                  Calculate the center of mass.",
"-Translate               Translate the coordinate set to the center of mass.",
"-Place -43.2,33,0        Place the center-of-mass at the given coordinates.",
"-Axes                    Calculate the orthogonal axes.",
"-Sphericity 1.5          Calculate the sphericity with the given angular sampling.",
"-Radialdensity 0.05,6,1  Calculate the radial density function at a given sampling",
"                         and distance cutoff (angstrom) with optional wrapping.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-origin 10,-10,20        Origin for rotation (default 0,0,0, can be \"center\" for COM)",
" ",
"Input:",
"-parameters parm.star    Atomic properties parameter file (default atom_prop.star).",
"-compare mol.pdb         Reference to calculate RMSD.",
"-distance mol.pdb        Reference to calculate distance matrix.",
" ",
"Output:",
"-matrix file.mat         Distance matrix.",
"-image file.map          Distance matrix as an image.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
    Bstring    		atom_select("ALL");
	Bstring			chain_select;
	long			select_random(0);	// Number of atoms to select
    int 			set_minmax(0); 		// Flag for setting min and max
    double	   		rmin(0);   	    	// Inner radius
    double	   		rmax(10000);  		// Outer radius
	int 			set_invert(0);
	int				com_origin(0);		// Flag to set the center-of-mass as origin of rotation
    Vector3<double>	point;				// Inversion point
    Transform		t;					// Transformation details
	View			view;				// View to rotate to
	int 			rot_from(0); 		// Rotation from view flag
	Vector3<double>	min(-1e6,-1e6,-1e6); // Coordinate minima
	Vector3<double>	max(1e6,1e6,1e6);	// Coordinate maxima
	Vector3<double>	box = max - min;	// Enclosing box
	Vector3<double>	place;				// Coordinates for placing a molecule
	char			first_name(0);		// First molecule name for renaming
    int     		first_res(0);		// Don't renumber
	double			set_radius(0);		// Spherical radius to enforce
    double			occupancy(0);
    int 			range_first(0);		// For resetting occupancy
    int 			range_last(0);
    int 			set_occupancy(0);	// Flag for setting occupancy
    int 			set_pbc(0);			// Flag for resolving wrapped coordinates
    int 			set_wrap(0);		// Flag for wrapping coordinates
	int 			show_seq(0);		// Flag to show the sequence(s)
	int 			show_MW(0);			// Flag for calculating the MW(s)
	int 			show_Volume(0);		// Flag for calculating the Van der Waals volume
    int				show_composition(0); // Flag for calculating the composition
	int 			show_COM(0);		// Flag for calculating the center of mass
	int 			translate_COM(0);	// Flag for translating to the center of mass
	int				orth_axes(0);		// Flag for calculating orthogonal axes
	double			dang(0);			// Angular step size for sphericity
	int 			setBfac(0);
	double	 		rdf_sampling(0);	// Sampling for radial density function
	double			rdf_cutoff(5);		// Distance cutoff for radial density function
	int 			wrap(0);			// No wrapping as default
	Bstring			paramfile;			// Use default parameter file
	Bstring			compfile;			// Reference to calculate RMSD
	Bstring			distfile;			// Reference to calculate distance matrix
	Bstring			matfile;			// Distance matrix
	Bstring			imgfile;			// Distance matrix as an image
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			atom_select = curropt->value;
			if ( atom_select.length() < 1 )
				cerr << "-select: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "random" )
        	if ( ( select_random = curropt->value.integer() ) < 1 )
				cerr << "-random: A number of atoms must be specified!" << endl;
		if ( curropt->tag == "chains" ) {
			chain_select = curropt->value;
			if ( chain_select.length() < 1 )
				cerr << "-chains: A selection must be specified!" << endl;
		}
		if ( curropt->tag == "invert" ) {
			point = curropt->vector3();
            set_invert = 1;
		}
		if ( curropt->tag == "translate" )
			t.trans = curropt->vector3();
		if ( curropt->tag == "scale" )
			t.scale = curropt->vector3();
		if ( curropt->tag == "rotate" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else {
				t.axis = Vector3<double>(d[0], d[1], d[2]);
				t.angle = d[3] * M_PI/180.0;
			}
		}
		if ( curropt->tag == "toview" )
			view = curropt->view();
		if ( curropt->tag == "fromview" ) {
			view = curropt->view();
			rot_from = 1;
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value.contains("cen") )
				com_origin = 1;
			else
				t.origin = curropt->origin();
		}
		if ( curropt->tag == "limits" ) {
        	if ( curropt->box(min, max) < 6 )
				cerr << "-limits: All 6 limits must be specified!" << endl;
			else
				set_minmax = 1;
		}
		if ( curropt->tag == "radius" )
        	if ( ( set_radius = curropt->value.real() ) < 1 )
				cerr << "-radius: The radius must be specified!" << endl;
		if ( curropt->tag == "rename" )
        	if ( ( first_name = curropt->value[0] ) < 1 )
				cerr << "-rename: The first molecule name must be specified!" << endl;
		if ( curropt->tag == "renumber" )
        	if ( ( first_res = curropt->value.integer() ) < 1 )
				cerr << "-renumber: The first residue number must be specified!" << endl;
		if ( curropt->tag == "occupancy" ) {
        	if ( curropt->values(range_first, range_last, occupancy) < 3 )
				cerr << "-occupancy: A range and value must be specified!" << endl;
			else
				set_occupancy = 1;
		}
		if ( curropt->tag == "box" ) {
			box = curropt->size();
        	set_pbc = 1; 		// Resolve periodic boundaries
		}
		if ( curropt->tag == "wrap" ) {
			box = curropt->size();
			set_wrap = 1; 		// Wrap within periodic boundaries
		}
 		if ( curropt->tag == "Bfactor" )
       		setBfac = 1;
		if ( curropt->tag == "Sequence" )
        	show_seq = 1;
		if ( curropt->tag == "Mass" )
        	show_MW = 1;
		if ( curropt->tag == "Volume" ) {
        	if ( ( show_Volume = curropt->value.integer() ) < 1 )
				cerr << "-Volume: A wrapping flag must be specified!" << endl;
			else
				show_Volume++;
		}
		if ( curropt->tag == "Composition" )
        	show_composition = 1;
 		if ( curropt->tag == "Center" )
       		show_COM = 1;
		if ( curropt->tag == "Translate" )
        	translate_COM = 1;
		if ( curropt->tag == "Place" )
			place = curropt->vector3();
		if ( curropt->tag == "Axes" )
        	orth_axes = 1;
		if ( curropt->tag == "Sphericity" ) {
         	if ( ( dang = curropt->value.real() ) < 0.01 )
				cerr << "-Sphericity: The angular step size must be specified!" << endl;
			else
				dang *= M_PI/180.0;
		}
		if ( curropt->tag == "Radialdensity" )
        	if ( curropt->values(rdf_sampling, rdf_cutoff, wrap) < 1 )
				cerr << "-Radialdensity: The sampling must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "compare" )
			compfile = curropt->filename();
		if ( curropt->tag == "distance" )
			distfile = curropt->filename();
		if ( curropt->tag == "matrix" )
			matfile = curropt->filename();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
    }
	option_kill(option);
	
	if ( optind >= argc ) bexit(-1);
    
	double		ti = timer_start();
	
	Bstring		filename(argv[optind++]);
	Bmolgroup*	mgtemp = NULL;
    Bmolgroup*	molgroup = read_molecule(filename, atom_select, paramfile);
	if ( !molgroup ) {
		cerr << "Error: Coordinate file not read!" << endl;
		bexit(-1);
	}
	
	if ( chain_select.length() )
		molgroup_select_chains(molgroup, chain_select);

    /* Select some atoms in a ring between rmin and rmax */
    if ( ( rmin > 0 ) || ( rmax < 1000 ) )
		molgroup_coor_select_ring(molgroup, rmin, rmax);
    
    /* Select atoms between desired limits */
	if ( set_minmax )
		molgroup_coor_select(molgroup, min, max);
    
    if ( select_random )
		molgroup_select(molgroup, select_random);
	
	if ( com_origin )
		t.origin = molgroup_center_of_mass(molgroup);
	else if ( show_COM )
		molgroup_show_center_of_mass(molgroup);
	
	if ( translate_COM )
		molgroup_shift_to_center_of_mass(molgroup);
		
	if ( place.length() > 0 )
		molgroup_place_at_coordinates(molgroup, place);
		
	if ( orth_axes )
		molgroup_principal_axes(molgroup, NULL);
	
	if ( dang > 0 )
		molgroup_sphericity(molgroup, dang);
	
	if ( set_radius )
		molgroup_set_radius(molgroup, set_radius);

	if ( set_invert )
		molgroup_coor_invert(molgroup, point);
	
	if ( fabs(t.scale.volume() - 1) < 1e-10 ) {
		if ( view[2] < 1 ) {
			if ( rot_from )
				mgtemp = molgroup_rotate_from_view(molgroup, view, t.origin, t.trans);
			else
				mgtemp = molgroup_rotate_to_view(molgroup, view, t.origin, t.trans);
			molgroup_kill(molgroup);
			molgroup = mgtemp;
		} else if ( t.angle == 0 ) {
			molgroup_coor_shift(molgroup, t.trans);
		} else {
			molgroup_coor_rotate(molgroup, t);
		}
	} else {
		molgroup_coor_transform(molgroup, t);
	}

	if ( set_pbc || set_wrap ) molgroup->box = box;
    
    if ( set_pbc ) molgroup_resolve_pbc(molgroup);
	else if ( set_wrap ) molgroup_pack_in_periodic_box(molgroup);
    
	molgroup_stats(molgroup);
	
	if ( setBfac )
		molgroup_Bfactors(molgroup);

	if ( show_seq )
		molgroup_print_sequence(molgroup);
	
	if ( show_MW ) {
		molgroup_weight_from_sequence(molgroup);
		molgroup_weight_from_atoms(molgroup);
	}
	
	if ( show_Volume )
		molgroup_volume(molgroup, paramfile, show_Volume-1);
	
	if ( show_composition )
		molgroup_composition(molgroup, paramfile);
//		molgroup_elements(molgroup);
	
    if ( first_name )
		molgroup_rename(molgroup, first_name);
    
    if ( first_res > 0 )
		molgroup_residue_renumber(molgroup, first_res);
    
    if ( set_occupancy )
		molgroup_coor_reset_occupancy(molgroup, range_first, 
				range_last, occupancy);
    
	if ( rdf_sampling > 0 )
		molgroup_radial_density(molgroup, rdf_sampling, rdf_cutoff, wrap);
	
	Bmolgroup*	molgroup2 = NULL;
	if ( compfile.length() ) {
		molgroup2 = read_molecule(compfile, atom_select, paramfile);
		if ( !molgroup2 ) {
			cerr << "Error: Template file not read!" << endl;
		} else {
			molgroup_calculate_rmsd(molgroup, molgroup2);
			molgroup_kill(molgroup2);
		}
	}

	if ( distfile.length() ) {
		molgroup2 = read_molecule(distfile, atom_select, paramfile);
		if ( !molgroup2 ) {
			cerr << "Error: Template file not read!" << endl;
		} else {
			Matrix	mat = mol_distance_matrix(molgroup->mol, molgroup2->mol);
			molgroup_kill(molgroup2);
			if ( matfile.length() && mat.rows() )
				mat.write(matfile);

			if ( imgfile.length() && mat.rows() ) {
				Bimage*		pimg = new Bimage(mat, 1);
//				pimg->change_type(nudatatype);
				write_img(imgfile, pimg, 0);
				delete pimg;
			}
		}
	}
			
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
