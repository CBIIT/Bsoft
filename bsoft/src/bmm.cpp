/**
@file	bmm.cpp
@brief	A program to do model mechanics
@author Bernard Heymann
@date	Created: 20100223
@date 	Modified: 20210224
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "model_mechanics.h"
#include "model_transform.h"
#include "model_color.h"
#include "model_select.h"
#include "model_links.h"
#include "model_neighbors.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmm [options] in1.star [in2.star ...]",
"--------------------------------------------",
"Mechanical modeling.",
" ",
"Actions:",
"-center                  Center before all other operations.",
"-minimize 150            Number of iterations for minimizing a model.",
"-dynamics 500            Number of iterations for running dynamics on a model.",
" ",
"Selections:",
"-all                     Reset selection to all components before other selections.",
"-select #Mod1@14         Select models and components.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
"-map image.pif,2         Map and image number associated with model.",
" ",
"Parameters for minimization:",
"-shift 6.8               Maximum shift per iteration (angstrom, default infinite).",
" ",
"Parameters for dynamics:",
"-timestep 0.01           Integration time step (default 0.001).",
"-velocitylimit 0.01      Limit on the velocity (default 0.1 per time step).",
"-friction 0.2            Friction constant (default 1 = no friction).",
" ",
"Parameters for mechanics:",
"-neighbors 7,45          Neighbor designation: number and neighborhood scale (default 6).",
"-cutoff 8.5              Cutoff for distance interactions.",
"-Kdistance 0.1           Distance force constant (default 0).",
"-Klink 0.01              Link force constant (default 0).",
"-Kangle 0.1              Angle force constant (default 0).",
"-Kpolyangle 0.1          Polygon angle force constant (default 0).",
"-Kpolygon 0.001          Polygon regularity force constant (default 0).",
"-Kpolyplane 0.1          Polygon planarity force constant (default 0).",
"-Kpoint 0.1,0.02         Point force constant and decay (default 0,0.01).",
"-Kradial 0.1             Radial force constant (default 0).",
"-radius 212              Reference radial distance (default 0).",
"-Kplane 0.2              Neighbor plane force constant (default 0).",
"-Kguide 0.5              Polyhedron guide force constant (default 0).",
"-Kmap 0.5                Map force constant (default 0).",
"-sigma 17                Sigma determining kernel size for fitting a map (default map pixel size).",
" ",
"Input:",
"-parameters parm.star    Model parameter file (default atom_prop.star).",
"-guide file.star         Polyhedron guide model file.",
" ",
"Output:",
"-output file.star        Output model file.",
"-writeparam dm.star      Write model parameter file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int 			reset(0);				// Keep selection as read from file
	int				center(0);				// Flag to center the structure
	long 			max_iter(0);			// Number of iterations/cycles for minimization
	int 			mm_type(0);				// Type of mechanics: 0=minimization, 1=dynamics
	Bstring			map_name;				// Density map reference
	long			img_num(0);				// Image number in density map file
	Bstring			map_path;				// Density map path
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	double			max_shift(0);			// Maximum shift per iteration
	double			velocitylimit(0.1);		// Limit on velocity per time step
	long			neighbors(6);			// Number of neighbors to generate
	double			neighbor_scale(0);		// Neighborhood scale
	Bstring			mod_select;				// Model and component selection
    Bstring    		atom_select("ALL");
	Bstring			paramfile;				// Use default parameter file
	Bstring			guidefile;				// Polyhedron guide file
	Bstring			outfile;				// Output model file name
	Bstring			mdfile;					// Distance matrix file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	Bmodparam		md;

	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "minimize" ) {
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-minimize: The number of iterations must be specified!" << endl;
			else mm_type = 0;
		}
		if ( curropt->tag == "dynamics" ) {
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-dynamics: The number of iterations must be specified!" << endl;
			else mm_type = 1;
		}
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 0.001 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 0.001 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "map" ) {
			map_name = curropt->value;
			img_num = (map_name.post(',')).integer();
			map_name = map_name.pre(',');
			if ( map_name.length() < 1 )
				cerr << "-map: A file name must be specified!" << endl;
		}
		if ( curropt->tag == "path" ) map_path = curropt->value;
		if ( curropt->tag == "shift" )
			if ( ( max_shift = curropt->value.real() ) < 0.001 )
				cerr << "-shift: The maximum shift distance must be specified!" << endl;
		if ( curropt->tag == "timestep" )
			if ( ( md.timestep = curropt->value.real() ) < 0.001 )
				cerr << "-timestep: The time step must be specified!" << endl;
		if ( curropt->tag == "velocitylimit" )
			if ( ( velocitylimit = curropt->value.real() ) < 0.001 )
				cerr << "-velocitylimit: The velocity limit must be specified!" << endl;
		if ( curropt->tag == "friction" )
			if ( ( md.Kfriction = curropt->value.real() ) < 1e-30 )
				cerr << "-friction: The friction constant must be specified!" << endl;
		if ( curropt->tag == "neighbors" )
			if ( curropt->values(neighbors, neighbor_scale) < 1 )
				cerr << "-neighbors: The number of neighbors must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( md.cutoff = curropt->value.real() ) < 1e-30 )
				cerr << "-cutoff: A distance must be specified!" << endl;
		if ( curropt->tag == "Kdistance" )
			if ( ( md.Kdistance = curropt->value.real() ) < 1e-30 )
				cerr << "-Kdistance: The distance force constant must be specified!" << endl;
		if ( curropt->tag == "Klink" )
			if ( ( md.Klink = curropt->value.real() ) < 1e-30 )
				cerr << "-Klink: The link force constant must be specified!" << endl;
		if ( curropt->tag == "Kangle" )
			if ( ( md.Kangle = curropt->value.real() ) < 1e-30 )
				cerr << "-Kangle: The angle force constant must be specified!" << endl;
		if ( curropt->tag == "Kpolyangle" )
			if ( ( md.Kpolyangle = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolyangle: The angle force constant for polygons must be specified!" << endl;
		if ( curropt->tag == "Kpolygon" )
			if ( ( md.Kpolygon = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolygon: The polygon regularity force constant must be specified!" << endl;
		if ( curropt->tag == "Kpolyplane" )
			if ( ( md.Kpolyplane = curropt->value.real() ) < 1e-30 )
				cerr << "-Kpolyplane: The polygon planarity force constant must be specified!" << endl;
		if ( curropt->tag == "Kpoint" )
			if ( curropt->values(md.Kpoint, md.pointdecay) < 1 )
				cerr << "-Kpoint: The point force constant must be specified!" << endl;
		if ( curropt->tag == "Kradial" )
			if ( ( md.Kradial = curropt->value.real() ) < 1e-30 )
				cerr << "-Kradial: The radial force constant must be specified!" << endl;
		if ( curropt->tag == "radius" )
			if ( ( md.radius = curropt->value.real() ) < 1e-30 )
				cerr << "-radius: The reference radial distance must be specified!" << endl;
		if ( curropt->tag == "Kplane" )
			if ( (md.Kplane = curropt->value.real() ) < 1e-30 )
				cerr << "-Kplane: The neighbor plane force constant must be specified!" << endl;
		if ( curropt->tag == "Kguide" )
			if ( ( md.Kguide = curropt->value.real() ) < 1e-30 )
				cerr << "-Kguide: The polyhedron guide force constant must be specified!" << endl;
		if ( curropt->tag == "Kmap" )
			if ( ( md.Kmap = curropt->value.real() ) < 1e-30 )
				cerr << "-Kmap: The map force constant must be specified!" << endl;
		if ( curropt->tag == "sigma" )
			if ( ( md.sigma = curropt->value.real() ) < 1e-30 )
				cerr << "-sigma: The sigma value must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "guide" )
			guidefile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "writeparam" )
			mdfile = curropt->filename();
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

	if ( mod_select.length() ) model_select(model, mod_select);

	if ( map_name.length() ) {
		model->mapfile(map_name.str());
		model->image_number(img_num);
		model_check(model, map_path);
	} else if ( map_path.length() ) {
		model_check(model, map_path);
	}

	if ( paramfile.length() )
		update_dynamics_parameters(md, paramfile);
	else
		model_param_generate(md, model);

	if ( md.Kguide > 0 && guidefile.length() )
		md.guide = read_model(guidefile, paramfile);
	
//	model_param_set_type_indices(model, md);
	model->update_component_types(md.comptype);
	
//	model_param_display(md);
	
	if ( center ) models_process(model, model_center);

	if ( neighbor_scale ) model_set_neighbors(model, neighbors, neighbor_scale);
	else model_set_neighbors(model, neighbors);

	model_check(model, map_path);

	model_selection_stats(model);

	if ( max_iter )
		model_mechanics(model, md, mm_type, max_iter, max_shift, velocitylimit);

	model_color_by_fom(model);

	if ( compradius > 0 ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius > 0 ) models_process(model, linkradius, model_set_link_radius);
	
	if ( outfile.length() ) {
		write_model(outfile, model);
	}
	
	if ( mdfile.length() )
		write_dynamics_parameters(mdfile, md);

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

