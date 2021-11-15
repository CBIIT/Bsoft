/**
@file	bpoly.cpp
@brief	A tool to manipulate polyhedral coordinate files
@author Bernard Heymann
@date	Created: 20010828
@date 	Modified: 20120920
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_mechanics.h"
#include "model_transform.h"
#include "model_shell.h"
#include "model_path.h"
#include "model_select.h"
#include "model_views.h"
#include "model_links.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bpoly [options] in1.star [in2.star...]",
"---------------------------------------------",
"Manipulate model polyhedral coordinate files.",
" ",
"Actions:",
"-center                  Center before all other operations.",
"-rename                  Rename components based on the number of links.",
"-vertextypes 1           Determine vertex types from adjacent polygons (1=3-digit, 2=6-digit).",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-analyze                 Analyze the polyhedrons.",
"-angles                  Analyze the angles in polyhedrons.",
"-energy 112.4            Calculate the energy components of models using the given reference angle.",
"                         If the angle is 0, the inner polygon angle will be used as reference.",
"-pentagon                Calculate the pentagon adjacency.",
"-sphericity              Estimate the spherical nature of the polyhedron.",
"-ellipsoidicity          Estimate the ellipsoid nature of the polyhedron.",
"-wiener                  Calculate the Wiener index for each polyhedron.",
"-eigen                   Do an eigenanalysis for each polyhedron.",
"-coordinates             Generate new coordinates based on spherical harmonics.",
"-dual                    Calculate the dual of a polyhedron.",
"-faces 5                 Generate a dual using only the faces of the given order.",
"-find 34.7               Find symmetry axes with an RMSD below the given value.",
"-regularize 150          Number of iterations for regularizing a polyhedron.",
" ",
"Selections:",
"-all                     Reset selection to all components before other selections.",
"-closed order,3          Select based on valency (valency,<n>) or polygon order (order,<n>).",
"-select #232@14          Select models and components.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
" ",
"Regularization parameters:",
"-distance 212            Interaction reference distance.",
"-linklength 150          Vertex link reference length.",
"-Kdistance 0.1           Distance force constant (default 0).",
"-Klink 0.01              Link force constant (default 0.1).",
"-Kpolyangle 0.1          Angle force constant (default 0.01).",
"-Kpolygon 0.001          Polygon regularity force constant (default 0).",
"-Kpolyplane 0.1          Polygon planarity force constant (default 0).",
"-Kpoint 0.1,0.02         Point force constant and decay (default 0,0.01).",
" ",
"Input:",
"-parameters parm.star    Molecular parameter file (default atom_prop.star).",
" ",
"Output:",
"-output file.star        Output model parameter file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int 			reset(0);				// Keep selection as read from file
	int				closure_rule(0);		// Closure rule: 1=valency, 2=order
	int				val_order(0);			// Valency or order - depending on rule
	Bstring			mod_select;				// Model and component selection
	int				center(0);				// Flag to center the structure
	int				rename(0);				// Flag to rename components
	int				vertextypes(0);			// Flag to determine vertex types
	Bstring			calc_views;				// Mode to calculate component views
	int				angles(0);				// Flag to analyze the polygon angles
	int				analyze(0);				// Flag to analyze the structure
	double			angle_ref(0);			// Reference angle and flag to calculate the energy components of models
	int				pentagon(0);			// Flag to calculate the pentagon adjacency
	int				sphere(0);				// Flag to estimate the sphericity
	int				ellips(0);				// Flag to estimate the ellipsoidicity
	int				wiener(0);				// Flag to calculate the Wiener index
	int				eigen(0);				// Flag to do an eigenanalysis
	int				coor(0);				// Flag to generate new coordinates
	int				dual(0);				// Flag to calculate dual of polyhedron
	int				faces(0);				// Flag to calculate dual of polyhedron restricted to specific faces
	double			find(0);				// RMSD cutoff to find symmetry axes
	long 			reg_iter(0);			// Number of iterations/cycles for regularization
	double			distance(0);			// Unit length separating vertices
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	double			linklength(0);			// Unit length linking vertices
	double			Kdistance(0);			// Distance interaction strength
	double			Klink(0.1);				// Link strength
	double			Kpolyangle(0.01);		// Angle strength
	double			Kpolygon(0);			// Polygon regularity strength
	double			Kpolyplane(0);			// Polygon planarity strength
	double			Kpoint(0);				// Point force strength
	double			decay(0.01);			// Point force decay with distance
	Bstring			map_path;				// Path to density map
    Bstring    		atom_select("ALL");
	Bstring			paramfile;				// Use default parameter file
	Bstring			outfile;				// Output parameter file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*			curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "closed" ) {
			if ( curropt->value[0] == 'v' ) closure_rule = 1;
			if ( curropt->value[0] == 'o' ) closure_rule = 2;
			if ( curropt->value.contains(",") )
				val_order = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "rename" ) rename = 1;
		if ( curropt->tag == "vertextypes" )
			if ( ( vertextypes = curropt->value.integer() ) < 1 )
				cerr << "-vertextypes: A selection number must be specified!" << endl;
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "analyze" ) analyze = 1;
		if ( curropt->tag == "angles" ) angles = 1;
		if ( curropt->tag == "energy" ) {
			if ( ( angle_ref = curropt->value.real() ) < 1e-30 )
				cerr << "-energy: A reference angle must be specified!" << endl;
			else {
				angle_ref *= M_PI/180.0;
				if ( angle_ref <= 0 ) angle_ref = -1;
			}
		}
		if ( curropt->tag == "pentagon" ) pentagon = 1;
		if ( curropt->tag == "sphericity" ) sphere = 1;
		if ( curropt->tag == "ellipsoidicity" ) ellips = 1;
		if ( curropt->tag == "wiener" ) wiener = 1;
		if ( curropt->tag == "eigen" ) eigen = 1;
		if ( curropt->tag == "coordinates" ) coor = 1;
		if ( curropt->tag == "dual" ) dual = 1;
		if ( curropt->tag == "faces" )
			if ( ( faces = curropt->value.integer() ) < 1 )
				cerr << "-faces: A polygon order must be specified!" << endl;
		if ( curropt->tag == "find" )
			if ( ( find = curropt->value.real() ) < 1e-30 )
				cerr << "-find: A RMSD cutoff must be specified!" << endl;
		if ( curropt->tag == "regularize" )
			if ( ( reg_iter = curropt->value.integer() ) < 1 )
				cerr << "-regularize: The number of iterations must be specified!" << endl;
		if ( curropt->tag == "distance" )
			if ( ( distance = curropt->value.real() ) < 1e-30 )
				cerr << "-distance: The reference distance separating vertices must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 1e-30 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 1e-30 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 1e-30 )
				cerr << "-linklength: The link length must be specified!" << endl;
		if ( curropt->tag == "Kdistance" )
			if ( ( Kdistance = curropt->value.real() ) < 1e-30 )
				cerr << "-Kdistance: The distance interaction strength must be specified!" << endl;
		if ( curropt->tag == "Klink" )
			if ( ( Klink = curropt->value.real() ) < 1e-30 )
				cerr << "-Klink: The link strength must be specified!" << endl;
		if ( curropt->tag == "Kpolyangle" )
			if ( ( Kpolyangle = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolyangle: The angle strength must be specified!" << endl;
		if ( curropt->tag == "Kpolygon" )
			if ( ( Kpolygon = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolygon: The polygon regularity strength must be specified!" << endl;
		if ( curropt->tag == "Kpolyplane" )
			if ( ( Kpolyplane = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolyplane: The polygon planarity strength must be specified!" << endl;
		if ( curropt->tag == "Kpoint" )
			if ( curropt->values(Kpoint, decay) < 1 )
				cerr << "-Kpoint: The point force strength must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No model files specified!" << endl;
		bexit(-1);
	}

	Bmodel*		model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( reset ) models_process(model, model_reset_selection);
	
	if ( !model->poly ) model_poly_generate(model);
	
	if ( closure_rule ) model_select_closed(model, closure_rule, val_order);
	
	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( rename ) models_process(model, model_rename_components);
	
	if ( coor ) model_poly_sphere_coor(model);		// Pointer to eigenvalue array not freed! Single model
	
	if ( linklength > 0 ) models_process(model, linklength, model_set_link_length);
	
	if ( center ) models_process(model, model_center);

	if ( vertextypes == 1 ) model_vertex_types(model);
	else if ( vertextypes == 2 ) model_extended_vertex_types(model);

	if ( calc_views.length() ) model_calculate_views(model, calc_views);
	
	if ( analyze ) model_poly_analyze(model);

	if ( angles ) model_poly_angles(model);
	
	if ( angle_ref ) model_poly_energy(model, angle_ref);
	
	if ( pentagon ) model_poly_pentagon_adjacency(model);

	Bmodel*			mp;
	double			sph, ell;
	if ( sphere && ellips ) {
		cout << "Model\tSphericity\tEllipsoidicity" << endl; 
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			sph = model_sphericity(mp);
			ell = model_ellipsoidicity(mp);
			cout << mp->identifier() << tab << sph << tab << ell << endl;
		}
	} else if ( sphere ) {
		cout << "Model\tSphericity" << endl; 
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			sph = model_sphericity(mp);
			cout << mp->identifier() << tab << sph << endl;
		}
	} else if ( ellips ) {
		cout << "Model\tEllipsoidicity" << endl; 
		for ( mp = model; mp; mp = mp->next ) if ( mp->select() ) {
			ell = model_ellipsoidicity(mp);
			cout << mp->identifier() << tab << ell << endl;
		}
	}
	
	if ( wiener ) model_wiener_index(model);

	if ( eigen ) model_poly_eigenvalues(model, 1);

	if ( dual || faces ) {
		Bmodel*		new_model = model_poly_dual(model, faces);
		if ( new_model ) {
			model_kill(model);
			model = new_model;
		}
	}
	
	model_check(model, map_path);

	model_selection_stats(model);

	if ( find ) {
		if ( verbose )
			cout << "Model\tPgroup" << endl;
		for ( mp = model; mp; mp = mp->next ) {
			model_poly_find_symmetry(mp, find);
			if ( verbose )
				cout << mp->identifier() << tab << mp->symmetry() << endl;
		}
		if ( verbose )
			cout << endl;
	}
	
	if ( reg_iter )
		model_regularize(model, reg_iter, distance, Kdistance, 
				Klink, Kpolyangle, Kpolygon, Kpolyplane, Kpoint, decay);

	if ( compradius > 0 ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius > 0 ) models_process(model, linkradius, model_set_link_radius);
	
	if ( outfile.length() ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


