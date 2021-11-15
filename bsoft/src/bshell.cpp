/**
@file	bshell.cpp
@brief	Generates a shell point component model
@author Bernard Heymann
@date	Created: 20071210
@date 	Modified: 20160802
**/

#include "rwmodel.h"
#include "model_create.h"
#include "model_shell.h"
#include "model_map.h"
#include "model_select.h"
#include "model_transform.h"
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
"Usage: bshell [options] in1.star [in2.star ...]",
"-----------------------------------------------",
"Generates and operates on shell models.",
" ",
"Actions for creation action:",
"-shell 134.2             Generate a shell with the given radius.",
"-dodecahedron 215,3,0.3  Generate an dodecahedral shell with the given radius and divisions,",
"                         and a fraction to make it spherical.",
"-icosahedron 238.5,3,0.5 Generate an icosahedral shell with the given radius and divisions,",
"                         and a fraction to make it spherical.",
"-circle 527,243          Generate a circle with the given radius and z offset.",
"-ellipse 620.5,314.9,5.2 Generate an ellipse with the given x and y semi-axes lengths and z offset.",
"-ellipsoid 24.3,51.2,68  Generate an ellipsoid with the given x,y,z semi-axes.",
"-cylinder 0,1,0,88.3,155 Generate a cylinder with the given direction, radius and length.",
"-add -14.7               Add a shell offset by a radial distance from selected components.",
"-fromcomponents 18.6,1   Generate shells from components with given radius and flag for 2D.",
" ",
"Actions:",
"-merge                   Merge models into one after all other operations.",
"-guide poly.star         Input guide polyhedron to adjust the shell model.",
"-linklength 56.3         Generate links with this maximum length for the sphere model.",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-rdf 1.2                 Calculate the radial distribution function with the given sampling (angstrom).",
"-fit 15.6,200,1          Fits a shell model to a map using the given resolution limits,",
"                         and a flag to indicate negative contrast.",
"-radial                  Calculate a radial profile based on a shell model.",
"-powerspectrum out.map   Calculate 2D power spectra around components in the shell.",
" ",
"Selections:",
"-all                     Reset selection to all models and components before other selections.",
"-select #Mod1@14         Select models and components.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-id model1               Set model identifier.",
"-newtype MEM             Type for all new components.",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
" ",
"Parameters for generating a model:",
"-points 152              Number of points.",
"-separation 2.5          Separation distance between points.",
" ",
"Parameter for adjusting to a guide model:",
"-fraction 0.6            Fraction of adjustment to a guide model (default 1).",
"-curved                  Flag to indicate curved hull or shell surface (default not).",
" ",
"Parameters for power spectra:",
"-size 45,40,20           Size of densities to extract (default 20,20,20 voxels).",
"-origin 0.8,-10,15.7     Set the origin for extracted densities.",
"-fftsize 256             Size of power spectrum (default from images).",
"-resolution 20,360       High and low resolution limits for cross-correlation (angstrom).",
"-annuli 2,95             Real space annular limits (default 0,inf pixels).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
"-map file.map            Input map file to associate with shell model.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			reset(0);					// Keep selection as read from file
	int				merge(0);					// Flag to merge models rather than concatenate
	double			shell_radius(-1);			// Radius of sphere
	int				dod_divisions(0);			// Number of divisions for dodecahedral sphere
	int				ico_divisions(0);			// Number of divisions for icosahedral sphere
	double			sphere_fraction(0);			// Fraction for spherical shape
	double			circle_radius(0);			// Circle radius
	double			circle_offset(0);			// Offset of circle in z
	Vector3<double>	ellipse_semi_z;				// Ellipse semi-axes and offset in z
	Vector3<double> ellipsoid_semi;				// Ellipsoid semi-axes
	Vector3<double>	cyl_dir;					// Cylinder direction
	double			cyl_rad(0), cyl_len(0);		// Cylinder radius and length
	double			add_distance(0);			// Distance for adding a new shell
	double			comp2shell(0), twod(0);		// Distance for converting components to shells
	double			linklength(0);				// Link length for generating links
	Bstring			calc_views;					// Mode to calculate component views
	double			rdf_interval(0);			// Calculate RDF at this sampling (0=not)
	int				fit(0);						// Flag to indicate fitting a shell to a map
	int				radial(0);					// Flag to calculate radial profile
	Bstring			mod_select;					// Model and component selection
	double			hires(0), lores(1e6);		// Limiting resolution for fitting to a map
	long 			ann_min(0);					// Minimum annulus
	long 			ann_max(1000);	 			// Maximum annulus, reset to maximum radius in image
	int				neg_flag(0);				// Flag to indicate negative contrast
	long			points(0);					// Number of points on the sphere
	Bstring			model_id;					// Model identifier
	Bstring			nutype;						// New component type to set
	double			compradius(0);				// Component display radius
	double			linkradius(0);				// Link display radius
	double			separation(0);				// Distance between points
	double			guide_fraction(1);			// Fraction to adjust to the guide polyhedron
	int				curv_flag(0);				// Adjust to guide with curved surfaces
	Vector3<long>	size(20,20,20);				// Size of densities to extract
	Vector3<double> origin;						// Origin for extracted densities
	int				ft_size(0);					// Power spectrum size
    Bstring    		atom_select("all");
	Bstring			paramfile;					// Input parameter file name
	Bstring			mapfile;					// Input map file name
	Bstring			guidefile;					// Guide model file
	Bstring			outfile;					// Output model parameter file name
	Bstring			outps;						// Output power spectrum file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "shell" )
			if ( ( shell_radius = curropt->value.real() ) < 0.001 )
				cerr << "-shell: The radius must be specified!" << endl;
		if ( curropt->tag == "dodecahedron" ) {
			if ( curropt->values(shell_radius, dod_divisions, sphere_fraction) < 1 )
				cerr << "-dodecahedron: The radius and divisions must be specified!" << endl;
			else if ( dod_divisions < 1 ) dod_divisions = 1;
		}
		if ( curropt->tag == "icosahedron" ) {
			if ( curropt->values(shell_radius, ico_divisions, sphere_fraction) < 1 )
				cerr << "-icosahedron: The radius and divisions must be specified!" << endl;
			else if ( ico_divisions < 1 ) ico_divisions = 1;
		}
		if ( curropt->tag == "circle" )
			if ( curropt->values(circle_radius, circle_offset) < 1 )
				cerr << "-circle: The circle radius must be specified!" << endl;
		if ( curropt->tag == "ellipse" ) {
			ellipse_semi_z = curropt->vector3();
			if ( ellipse_semi_z.length() < 1 )
				cerr << "-ellipse: The major and minor semiaxes must be specified!" << endl;
		}
		if ( curropt->tag == "ellipsoid" ) {
			ellipsoid_semi = curropt->vector3();
			if ( ellipsoid_semi.volume() < 1 )
				cerr << "-ellipsoid: The 3 semiaxes must be specified!" << endl;
		}
		if ( curropt->tag == "cylinder" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 5 )
				cerr << "-cylinder: All 5 parameters must be specified!" << endl;
			else {
				cyl_dir = Vector3<double>(d[0], d[1], d[2]);
				cyl_rad = d[3];
				cyl_len = d[4];
			}
		}
		if ( curropt->tag == "add" )
			if ( ( add_distance = curropt->value.real() ) < 0.001 )
				cerr << "-add: A distance must be specified!" << endl;
		if ( curropt->tag == "fromcomponents" )
			if ( curropt->values(comp2shell, twod) < 1 )
				cerr << "-fromcomponents: A distance must be specified!" << endl;
		if ( curropt->tag == "guide" )
			guidefile = curropt->filename();
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 1 )
				cerr << "-linklength: A link length must be specified!" << endl;
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "rdf" )
			if ( ( rdf_interval = curropt->value.real() ) < 0.1 )
				cerr << "-rdf: A sampling interval must be specified!" << endl;
		if ( curropt->tag == "fit" ) {
			if ( curropt->values(hires, lores, neg_flag) < 1 )
				cerr << "-fit: A high resolution limit must be specified!" << endl;
			else
				fit = 1;
		}
		if ( curropt->tag == "radial" ) radial = 1;
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "id" ) {
			model_id = curropt->value;
			if ( model_id.length() < 1 )
				cerr << "-id: An identifier must be specified!" << endl;
		}
		if ( curropt->tag == "newtype" ) {
			nutype = curropt->value;
			if ( nutype.length() < 1 )
				cerr << "-newtype: An type identifier must be specified!" << endl;
		}
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 1 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 1 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "points" )
			if ( ( points = curropt->value.integer() ) < 1 )
				cerr << "-points: The number of points must be specified!" << endl;
		if ( curropt->tag == "separation" )
			if ( ( separation = curropt->value.real() ) < 0.001 )
				cerr << "-separation: The separation distance must be specified!" << endl;
		if ( curropt->tag == "fraction" )
			if ( ( guide_fraction = curropt->value.real() ) < 0.001 )
				cerr << "-fraction: A fraction must be specified!" << endl;
		if ( curropt->tag == "curved" ) curv_flag = 1;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "fftsize" )
			if ( ( ft_size = curropt->value.integer() ) < 1 )
				cerr << "-fftsize: A Fourier transform size must be specified!" << endl;
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
			else
				if ( lores < hires ) swap(hires, lores);
		}
		if ( curropt->tag == "annuli" )
			if ( curropt->values(ann_min, ann_max) < 1 )
				cerr << "-annuli: At least a minimum annulus must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "powerspectrum" )
			outps = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();

	Bstring*		file_list = NULL;
	Bmodel*			model = NULL;		
	Bmodel*			mp = NULL;		
	Bmodel*			sph = NULL;		
	Bmodel*			guide = NULL;
	
	// Read all the model parameter files
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	
	if ( file_list ) {
		model = read_model(file_list);		
		string_kill(file_list);
		if ( !model ) {
			cerr << "Error: Input file not read!" << endl;
			bexit(-1);
		}
		if ( reset )
			for ( mp = model; mp; mp = mp->next )
				model_reset_selection(mp);
	}

	if ( guidefile.length() ) guide = read_model(guidefile, paramfile);

	sph = model;
	
	if ( shell_radius >= 0 ) {
		if ( dod_divisions > 0 )
			sph = model_create_dodecahedron(shell_radius, dod_divisions, sphere_fraction);
		else if ( ico_divisions > 0 )
			sph = model_create_icosahedron(shell_radius, ico_divisions, sphere_fraction);
		else
			sph = model_create_shell(points, shell_radius, separation);
	} else if ( circle_radius > 0) {
		sph = model_create_circle(circle_radius, circle_offset, separation);
	} else if ( ellipse_semi_z[0] > 0 && ellipse_semi_z[1] > 0) {
		sph = model_create_ellipse(ellipse_semi_z, separation);
	} else if ( ellipsoid_semi.volume() > 0) {
		sph = model_create_ellipsoid(ellipsoid_semi, separation);
	} else if ( cyl_len > 0 ) {
		sph = model_create_cylinder(cyl_dir, cyl_rad, cyl_len, separation);
	}
	
	if ( !sph ) {
		cerr << "Error: Ellipsoid model not generated!" << endl;
		bexit(-1);
	}
	
	if ( !model ) model = sph;
	else if ( model != sph ) {
		for ( mp = model; mp->next; mp = mp->next ) ;
		mp->next = sph;
	}
	
	if ( reset ) model_reset_selection(sph);
	
	if ( mod_select.length() ) model_select(sph, mod_select);

	if ( nutype.length() ) model_set_type(sph, nutype);
	
	if ( guide ) model_align_to_guide(sph, guide);
	
	if ( model_id.length() ) sph->identifier(model_id);
	
	if ( add_distance ) model_add_shell(sph, add_distance, nutype);

	if ( comp2shell ) {
		sph = model_components_to_shells(model, comp2shell, nutype, twod);
		model_kill(model);
		model = sph;
	}
	
	if ( mapfile.length() ) sph->mapfile(mapfile.str());

	if ( linklength > 0 ) model_link_list_generate(sph, linklength);
	
	if ( calc_views.length() ) model_calculate_views(sph, calc_views);

	if ( rdf_interval > 0 ) model_radial_distribution(sph, rdf_interval);

	if ( fit ) {
		if ( sph->mapfile().length() < 2 ) sph->mapfile(model->mapfile());
		model_shell_fit(sph, hires, lores, neg_flag);
	}

	if ( guide && guide_fraction > 0 )
		model_adjust_shell_to_guide(sph, guide, guide_fraction, curv_flag);
	
	if ( radial ) model_shell_radial_profile(sph);

	Bimage*		p = NULL;
	if ( outps.length() ) {
		p = model_shell_power_spectrum(sph, size, origin, ft_size, ann_min, ann_max, hires, lores);
		write_img(outps, p, 0);
		delete p;
	}

	if ( compradius ) model_set_component_radius(sph, compradius);

	if ( linkradius ) model_set_link_radius(sph, linkradius);

	if ( merge ) model_merge(model);
	
	model_selection_stats(model);
	
	// Write an output parameter format file if a name is given
    if ( outfile.length() ) {
		write_model(outfile, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

