/**
@file	bpartsel.cpp
@brief	Selection of single particles for 3D reconstruction
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20210517
**/

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_ctf.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "ps_micrograph.h" 
#include "symmetry.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpartsel [options] input.star [input.star]",
"-------------------------------------------------",
"Selects single particles for reconstruction and further processing.",
"Input particle selections are honoured, so that multiple restrictive selections",
"    can be made. To reset for a new selection, use the -all option.",
" ",
"Actions:",
"-check                   Check particle records.",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-renumber                Renumber particles starting from 1 for each micrograph or reconstruction.",
"-all                     Reset selection to all particles before other selections.",
"-none                    Unset selection to all particles before other selections.",
"-mgselect 3-5,8,12-15    Select micrographs:",
"                         Argument: 4-11,34: Select by the indicated micrographs (first = 1).",
"                         Argument: part: Select only micrographs with selected particles.",
"-divide 3                Divide the project into a number of subsets of fields.",
"-consolidate 5           Consolidate selection under the given selection number.",
"-setnumber 3             Select given selection number and set others to zero.",
"-setfom 0.3              Set the FOM (must be in the range 0-1).",
"-setasu D3               Set the views to the asymmetric unit.",
"-setgroups               Transfer selection numbers to groups.",
"-fixdefocus 0.15         Reset particle defocus differing from micrograph defocus.",
"-group 3022              Select particles belonging to this group.",
"-sets 20,1               Generate sets of this size of selected particles, ",
"                         with a flag not to select across micrograph or reconstruction boundaries.",
"-series 25,1             Select: Particles that are all selected in a series and generating sets",
"                         with a flag not to select across micrograph or reconstruction boundaries.",
"-frames 3,11             Select: Particles from a series of frames of a field-of-view.",
"-alternate 3             Select: Assign successive particles to different numbers.",
"-reselect psi,0.1,0.9    Reselect based on a range of values associated with this tag.",
"-origin 150,150,0,8      Select: Within a distance from the nominal origin: x,y,z,distance.",
"-View 0.2,-0.5,0.3,15    Select: Within an angular distance from a view: x,y,z,angle.",
"-sideview 17.8           Select: Side views within an angular distance from the equator.",
"-Euler 2,50,-30,2,0,40   Select: Euler angle ranges: psi, theta, phi.",
"-rank 4,1                Select a number of groups from FOM ranking and flag to adjust for defocus.",
"-fom 0.25,1              Select: FOM cutoff and flag to adjust for defocus.",
"-cv 0.55,1               Select: Cross validation cutoff and flag to adjust for defocus.",
"-exclusion 210           Select: Deselect overlapping particles (default 0).",
"-top 12.0,1              Select: number or % particles and flag to adjust for defocus, ranked by FOM[index].",
"-deviation -2.0          Select: FOM >= average + factor*standard_deviation, using FOM[index].",
"-angle 2.4               Select: difference angle between view vectors in focal or tilt series.",
"-closesttofocus          Select: Particles from micrographs closest to focus.",
"-furthestfromfocus       Select: Particles from micrographs furthest from focus.",
"-random 23.6%            Select: random selection to get the given number, ",
"                         number in groups (g), filaments (f), fraction (<1) or percentage (%).",
"-bootstrap 527           Select: random selection with replacement up to the given number.",
"-views 4.5,7.2,3         Select: random selection within views defined on a grid with the ",
"                         given theta and phi step sizes (degrees) and up to the given number.",
"-bestviews 4.5,7.2,3     Select: best particles within views defined on a grid with the ",
"                         given theta and phi step sizes (degrees), up to the given number.",
"-maxsmooth 0.7,3,5,4.4   Select: above fraction specified by smoothed maxima surface.",
"                         given theta and phi step sizes and smoothing sigma (degrees).",
"-Fomfirst                Set orientations to first in field-of-view series.",
"-Fombest                 Set orientations to best FOM[index] in field-of-view series.",
"-filament 60,opp.star    Selects filament particles with clear direction based on percentage.",
"                         Optionally compare with project with opposite direction.",
"-deselect 4,6,22         Deselect particles with the listed selection numbers.",
"-getselection            Prints all the selection numbers.",
"-show particle.fom       Show particle parameters indicated by the tag or selection number.",
"-Fomhistogram 100,0.1.2  Show FOM histogram with parameters: bins, minimum and maximum.",
"-focus                   Show defocus minimum and maximum.",
"-reset defocus           Reset particle defocus to micrograph defocus.",
"-delete                  Delete all non-selected particles from the parameter file.",
"-remove                  Do not write non-selected micrographs into the parameter file.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-symmetry C5             Point group symmetry (required with the -angle and -views options).",
"-index 1                 FOM index to select on (default 0, options -setfom, -top, -deviation, -bestviews, -Fombest).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-split 3                 Split micrograph or reconstruction data blocks into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": micrograph ID's are used as file names.",
"                         Argument: \"set\": particle sets based on selection are used as output.",
"-Particle views.ps       Postscript output with particle view and origin distributions.",
"-Micrograph out.ps       Postscript output with particle micrograph distributions.",
"-FOM histogram.ps        Postscript output with particle FOM distributions.",
"-classes                 Write stacks of particle images based on classes.",
"-export particles.txt,1  Export selected particle coordinates to a text file.",
"                         Flag: convert coordinates to physical location in angstroms.",
" ",
NULL
};

Bstring		read_list(Bstring filename)
{
	ifstream		f(filename.c_str());
	if ( f.fail() ) {
		cerr << "File " << filename << " not read!" << endl;
		bexit(-1);
	}
	
	string			s;
	Bstring			list;
	
	f >> s;
	list = s;
	
	while ( !f.eof() ) {
		f >> s;
		list += "," + s;
	}
	
	f.close();
	
	cout << list << endl;
	
	return list;
}


int			main(int argc, char** argv)
{
	// Initializing variables
	int				check_flag(0);				// Flag to check particle records
	int				use_rec(0);					// Flag to process reconstructions
	int				split(0);					// Output one big STAR file
	int				symmetry_asu(0);			// Set views to the asymmetric unit
	Bstring			symmetry_string("C1");		// Point group string
	int 			renumber(0);				// Flag to renumber particles
	double			max_def_dev(0);				// Maximum defocus deviation allowed
	int 			all(0);						// Flag to reset selection
	int 			unset(0);					// Flag to unset selection
	Bstring			mgselect;					// Micrograph selection
	int				divide(0);					// Number of output project subsets
	int 			consolidate(0);				// Flag to combine selections
	int				set_number(0);				// Selection number to keep
	int				set_groups(0);				// Flag to transfer selection numbers to groups
	int				group(0);					// Group to select
	int				set_size(0);				// Size of sets
	int				set_flag(0);				// Flag to keep set within micrograph
	int				select_series(0);			// Size of sets and to select particles where all are selected in a series
	int				alternate(0);				// Assign successive particles to different numbers
	int				frame_start(0), frame_end(1000);	// Frame selection from a field-of-view series
	double			setfom(-1);					// Set the FOM
	Bstring			tag_reselect;				// STAR tag for reselection
	double			reselect_min(-1e37);		// Reselection minimum
	double			reselect_max(1e37);			// Reselection maximum
	Vector3<double>	origin;						// Nominal origin
	double			dev_origin(0);				// Maximum distance to accept
	View			view;						// View to select around
	double			dev_angle(0);				// Angle from view
	double			side_angle(0);				// Angle from side view
	double			euler[6] = {0,0,0,0,0,0};	// Euler angle ranges
	int				ngroups(0);					// Group selection turned off
	double			angle_cutoff(0);			// No selection based on angle deviation between focal pairs
	int				fom_index(0);                       	// Index of FOM
	double			fom_cutoff[5] = {0,0,0,0,0};			// No selection based on FOM
	int				fom_defocus_adjust[5] = {0,0,0,0,0};	// Flag to adjust the FOM cutoff by defocus
	double			excl_dist(0);				// Minimum distance between particles
	double			fom_percentage(-1);			// Percentage particles to accept
	int				fom_best(0);				// Best number of particles to select
	double			fom_std_factor(-1e37);		// FOM standard deviation multiplier for cutoff
	long			random_number(0);			// Number of particles for random selection
	double			random_fraction(0);			// Fraction of particles for random selection
	long			random_group(0);			// Number of particles for random group selection
	long			random_filament(0);			// Number of maps for random filament selection
	int				bootstrap(0);				// Number of particles for bootstrap selection
	int				nviews(0);					// Number within views to select randomly
	double			threshfrac(0), maxsmooth(0);// Maxima smoothed surface for selection: fraction and sigma
	int				bestviews(0);				// Best number within views to select
	double			theta_step(0), phi_step(0);	// Anglular step sizes for views selection
	int				closest_to_focus(0);		// Flag to select closest to focus
	int				furthest_from_focus(0);		// Flag to select furthest from focus
	int				select_first_in_series(0);	// Flag to select first in series
	int				select_best_in_series(0);	// Flag to select best in series
	double			fil_dir_pct(0);				// Flag and percentage to check filament directions
	Bstring			fil_opp_name;				// Project with opposite filament directions
	Bstring			reset;                      // Reset a particle parameter from its micrograph
	Bstring			deselect;					// List of particles to deselect
	int				delete_deselected(0);		// Flag to delete deselected particles
	int				getselection(0);			// Flag to show all selection numbers
	int				show(-1);					// Show selected particles
	Bstring			tag_show;					// STAR tag for parameter selection
	long			FOMbins(0);					// Flag and bins to show the FOM histogram
	double			FOMmin(0), FOMmax(0);		// FOM histogram minimum and maximum
	int				focus(0);					// Flag to show minumum and maximum defocus
	Bstring			outfile;					// Output parameter file name
	Bstring			viewps;						// Postscript file name for views
	Bstring			mgps;						// Postscript file name for micrograph locations
	Bstring			FOMps;						// Postscript file name for FOM histogram
	int				write_classes(0);			// Flag to write stacks of particle classes
	int				write_flags(0);				// Flags to pass to the parameter file writing function
	Bstring			expfile;					// Particle coordinate export file
	int				exp_flags(0);				// Export flags
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "check" ) check_flag = 1;
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "renumber" ) renumber = 1;
		if ( curropt->tag == "all" ) all = 2;
		if ( curropt->tag == "none" ) unset = 1;
		if ( curropt->tag == "mgselect" ) mgselect = curropt->value;
		if ( curropt->tag == "setasu" ) {
			symmetry_asu = 1;
			symmetry_string = curropt->symmetry_string();
		}
		if ( curropt->tag == "fixdefocus" ) {
			if ( ( max_def_dev = curropt->value.real() ) < 0.001 )
				cerr << "-fixdefocus: A defocus difference must be specified!" << endl;
			else
				if ( max_def_dev < 10 ) max_def_dev *= 1e4;	// assume in um
		}
		if ( curropt->tag == "divide" )
			if ( ( divide = curropt->value.integer() ) < 2 )
				cerr << "-divide: A number of subsets must be specified!" << endl;
		if ( curropt->tag == "consolidate" )
			if ( ( consolidate = curropt->value.integer() ) < 1 )
				cerr << "-consolidate: A number must be specified!" << endl;
		if ( curropt->tag == "setnumber" )
			if ( ( set_number = curropt->value.integer() ) < 1 )
				cerr << "-setnumber: A selection number must be specified!" << endl;
		if ( curropt->tag == "sets" )
			if ( curropt->values(set_size, set_flag) < 1 )
				cerr << "-sets: A set size must be specified!" << endl;
		if ( curropt->tag == "setgroups" ) set_groups = 1;
		if ( curropt->tag == "group" )
			if ( ( group = curropt->value.integer() ) < 1 )
				cerr << "-group: A group must be specified!" << endl;
		if ( curropt->tag == "series" )
			if ( curropt->values(select_series, set_flag) < 1 )
				cerr << "-series: A set size must be specified!" << endl;
		if ( curropt->tag == "frames" )
			if ( curropt->values(frame_start, frame_end) < 1 )
				cerr << "-frames: A starting frame must be specified!" << endl;
		if ( curropt->tag == "alternate" )
			if ( ( alternate = curropt->value.integer() ) < 2 )
				cerr << "-alternate: A number must be specified!" << endl;
		if ( curropt->tag == "setfom" ) {
			if ( ( setfom = curropt->value.real() ) <= 0 )
				cerr << "-setfom: A FOM value must be specified!" << endl;
			else {
				if ( setfom < 0 || setfom > 1 ) {
					cerr << "-setfom: The FOM value must be in the range 0-1! The FOM will not be changed." << endl;
					setfom = -1;
				}
			}
		}
		if ( curropt->tag == "reselect" ) {
			tag_reselect = curropt->value;
			Bstring*		slist = tag_reselect.split(",");
			tag_reselect = *slist;
			reselect_min = slist->next->real();
			reselect_max = slist->next->next->real();
			string_kill(slist);
		}
		if ( curropt->tag == "origin" )
			if ( curropt->values(origin[0], origin[1], origin[2], dev_origin) < 4 )
				cerr << "-origin: A 3-value origin and a distance must be specified!" << endl;
		if ( curropt->tag == "View" ) {
			view = curropt->view();
			if ( view.angle() ) dev_angle = view.angle();
		}
		if ( curropt->tag == "sideview" ) {
			if ( ( side_angle = curropt->value.real() ) < 0.01 )
				cerr << "-sideview: An angle must be specified!" << endl;
			else
				side_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "Euler" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			for ( size_t i=0; i<d.size(); i++ ) euler[i] = d[i] * M_PI/180;
			if ( d.size() < 6 )
				cerr << "-Euler: All 6 values for Euler angle ranges must be specified!" << endl;
		}
		if ( curropt->tag == "rank" ) 
			if ( curropt->values(ngroups, fom_defocus_adjust[0]) < 1 )
				cerr << "-rank: The number of groups must be specified!" << endl;
		if ( curropt->tag == "fom" ) {
			if ( curropt->values(fom_cutoff[0], fom_defocus_adjust[0]) < 1 )
				cerr << "-fom: A FOM cutoff must be specified!" << endl;
			else
				if ( fom_defocus_adjust[0] > 0 ) fom_defocus_adjust[0] = 1;
		}
		if ( curropt->tag == "cv" ) {
			if ( curropt->values(fom_cutoff[1], fom_defocus_adjust[1]) < 1 )
				cerr << "-cv: A cross validation cutoff must be specified!" << endl;
			else
				if ( fom_defocus_adjust[1] > 0 ) fom_defocus_adjust[1] = 1;
		}
		if ( curropt->tag == "exclusion" )
			if ( ( excl_dist = curropt->value.real() ) <= 0 )
				cerr << "-exclusion: An separation distance must be specified." << endl;
		if ( curropt->tag == "top" ) {
			if ( curropt->values(fom_percentage, fom_defocus_adjust[fom_index]) < 1 )
				cerr << "-top: A percentage must be specified!" << endl;
			if ( !curropt->value.contains("%") && !curropt->value.contains(".") ) {
				fom_best = (int) fom_percentage;
				fom_percentage = -1;
			}
		}
		if ( curropt->tag == "deviation" )
			if ( ( fom_std_factor = curropt->value.real() ) < 0.001 )
				cerr << "-deviation: A standard deviation multiplying factor must be specified!" << endl;
		if ( curropt->tag == "angle" ) {
			if ( ( angle_cutoff = curropt->value.real() ) < 1 )
				cerr << "-angle: An angle cutoff must be specified!" << endl;
			else {
				if ( angle_cutoff < 1 ) angle_cutoff = 1;
				if ( angle_cutoff > 180 ) angle_cutoff = 180;
				angle_cutoff *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "closesttofocus" )
			closest_to_focus = 1;
		if ( curropt->tag == "furthestfromfocus" )
			furthest_from_focus = 1;
		if ( curropt->tag == "Fombest" )
			select_best_in_series = 1;
		if ( curropt->tag == "Fomfirst" )
			select_first_in_series = 1;
		if ( curropt->tag == "filament" ) {
			if ( ( fil_dir_pct = curropt->value.real() ) < 0.01 )
				cerr << "-filament: A percentage must be specified!" << endl;
			else if ( curropt->value.contains(",") )
				fil_opp_name = curropt->value.post(',');
		}
		if ( curropt->tag == "deselect" )
			deselect = curropt->value;
		if ( curropt->tag == "getselection" )
			getselection = 1;
		if ( curropt->tag == "show" ) {
			if ( ( show = curropt->value.integer() ) < 1 ) {
				tag_show = curropt->value;
				if ( tag_show.length() < 1 )
					cerr << "-show: A selection number or tag must be specified!" << endl;
			}
		}
		if ( curropt->tag == "Fomhistogram" )
			if ( curropt->values(FOMbins, FOMmin, FOMmax) < 1 )
				cerr << "-Fomhistogram: A number of bins must be specified!" << endl;
		if ( curropt->tag == "focus" )
			focus = 1;
		if ( curropt->tag == "random" ) {
			if ( curropt->value.contains("g") ) {
				if ( ( random_group = curropt->value.integer() ) < 1 )
					cerr << "-random: A group must be specified!" << endl;
			} else if ( curropt->value.contains("f") ) {
				if ( ( random_filament = curropt->value.integer() ) < 1 )
					cerr << "-random: A number of maps must be specified!" << endl;
			} else if ( ( random_fraction = curropt->value.real() ) < 0.001 ) {
				cerr << "-random: A number, fraction or percentage must be specified!" << endl;
			} else {
				if ( curropt->value.contains("%") ) {
					random_fraction /= 100;
				} else if ( random_fraction >= 1.0 ) {
					random_number = (long) random_fraction;
					random_fraction = 0;
				}
			}
		}
		if ( curropt->tag == "bootstrap" )
			if ( ( bootstrap = curropt->value.integer() ) < 1 )
				cerr << "-bootstrap: A number of particles must be specified!" << endl;
		if ( curropt->tag == "views" ) {
			if ( curropt->values(theta_step, phi_step, nviews) < 3 )
				cerr << "-views: Two angles and a number of particles must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "bestviews" ) {
			if ( curropt->values(theta_step, phi_step, bestviews) < 3 )
				cerr << "-bestviews: Two angles and a number of particles must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "maxsmooth" ) {
			if ( curropt->values(threshfrac, theta_step, phi_step, maxsmooth) < 4 )
				cerr << "-maxsmooth: A fraction, two angles and a sigma value must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
				maxsmooth *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "index" ) {
			if ( ( fom_index = curropt->value.integer() ) < 0 )
				cerr << "-index: A FOM index must be specified!" << endl;
			else if ( fom_index < 0 || fom_index >= NFOM ) fom_index = 0;
		}
		if ( curropt->tag == "reset" )
			reset = curropt->value;
		if ( curropt->tag == "delete" )
			delete_deselected = 1;
		if ( curropt->tag == "remove" ) write_flags |= 2;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( curropt->value.contains("set") || curropt->value.contains("SET") ) split = -1;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
		}
		if ( curropt->tag == "Particle" )
			viewps = curropt->filename();
		if ( curropt->tag == "Micrograph" )
			mgps = curropt->filename();
		if ( curropt->tag == "FOM" )
			FOMps = curropt->filename();
		if ( curropt->tag == "classes" ) write_classes = 1;
		if ( curropt->tag == "export" ) {
			if ( curropt->value.contains(",") ) {
				expfile = curropt->value.pre(',');
				exp_flags = curropt->value.post(',').integer();
			} else {
				expfile = curropt->value;
			}
		}
    }
	option_kill(option);

	double		ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list, check_flag);
	string_kill(file_list);

	if ( !project ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	Bproject*		project_fil = NULL;
	if ( fil_opp_name.length() ) {
		project_fil = read_project(fil_opp_name, check_flag);
		if ( !project_fil ) {
			cerr << "Error: Opposite filament input file not read!" << endl;
			bexit(-1);
		}
	}
	
	if ( use_rec ) project->select = 1;
	
	Bsymmetry		sym(symmetry_string);

	if ( renumber ) project_renumber_particles(project);
	
	if ( all ) part_reset_selection(project, all);
	
	if ( unset ) part_unset_selection(project);
	
	if ( mgselect.length() ) {
		if ( mgselect[0] == 'p' ) project_select_with_particles(project, -1);
		else part_select_micrograph(project, mgselect);
	}
	
	if ( group > 0 ) part_select_group(project, group);

	if ( consolidate ) part_consolidate_selection(project, consolidate);

	if ( alternate > 0 ) part_set_multi_maps(project, -1, alternate);

	if ( set_number ) part_set_selection(project, set_number);

	if ( reset.length() ) project_reset(project, reset);
	
	if ( tag_reselect.length() )
		part_reselect(project, tag_reselect, reselect_min, reselect_max);
	
	if ( symmetry_asu )
		project_set_particle_asu_views(project, symmetry_string);

	if ( max_def_dev )
		part_fix_defocus(project, max_def_dev);

	if ( dev_origin )
		part_origin_select(project, origin, dev_origin);
	
	if ( dev_angle )
		part_view_select(project, view, dev_angle);
	
	if ( side_angle )
		part_side_view_select(project, side_angle);

	if ( euler[0]+euler[2]+euler[4] != 0 || euler[1]+euler[3]+euler[5] != 0 )
		part_euler_angle_select(project, euler);
	
	for ( i=0; i<5; i++ ) {
		if ( fom_cutoff[i] ) {
			if ( fom_defocus_adjust[i] )
				part_fom_defocus_fit_deselect(project, i, fom_cutoff[i]);
			else
				part_deselect(project, i, fom_cutoff[i]);
		}
	}
	
	if ( excl_dist )
		part_deselect_redundant(project, excl_dist, -1, fom_index);

	
	if ( fom_percentage >= 0 )
		part_select_percentage(project, fom_percentage, fom_index, fom_defocus_adjust[fom_index]);
	
	if ( fom_best > 0 )
		part_select_best(project, fom_best, fom_index, fom_defocus_adjust[fom_index]);
	
	if ( fom_std_factor > -1e36 )
		part_select_FOM_avg_std(project, fom_std_factor, fom_index);
	
	if ( angle_cutoff )
		part_series_comparison(project, sym, angle_cutoff);
	
	if ( closest_to_focus )
		part_select_closest_to_focus(project);
		
	if ( furthest_from_focus )
		part_select_furthest_from_focus(project);
		
	if ( select_first_in_series )
		part_set_first_view_in_series(project);
		
	if ( select_best_in_series )
		part_set_best_view_in_series(project, fom_index);

	if ( ngroups > 1 )
		part_select_FOM_groups(project, ngroups, fom_index, fom_defocus_adjust[fom_index]);
	
	if ( random_number )
		part_select_random(project, random_number);
	else if ( random_fraction )
		part_select_random_fraction(project, random_fraction);
	else if ( random_group )
		part_select_random_group(project, random_group);
	else if ( random_filament )
		part_select_random_filaments(project, random_filament);
		
	if ( bootstrap )
		part_select_bootstrap(project, bootstrap);

	if ( nviews )
		part_select_random_within_view(project, sym, theta_step, phi_step, nviews);

	if ( maxsmooth )
		part_select_maxsmooth(project, sym, theta_step, phi_step, threshfrac, maxsmooth, fom_index);

	if ( bestviews )
		part_select_best_within_view(project, sym, theta_step, phi_step, bestviews, fom_index);

	if ( setfom > -1 )
		part_set_FOM(project, fom_index, setfom);

	if ( select_series ) part_select_series(project, select_series, set_flag);
	
	if ( set_size ) part_select_sets(project, set_size, set_flag);

	if ( frame_start ) part_select_frames(project, frame_start, frame_end);

	if ( fil_dir_pct > 0 ) {
		if ( project_fil ) part_filament_direction(project, project_fil, fil_dir_pct);
		else part_filament_direction(project, fil_dir_pct);
	}

	if ( set_groups ) part_select_to_group(project);
	
	if ( deselect.length() ) {
		if ( deselect.contains(".") ) deselect = read_list(deselect);
		part_deselect_from_list(project, deselect);
		project_trim_class_averages(project, deselect);
		deselect = 0;
	}

	if ( delete_deselected ) part_delete_deselected(project);
	
	if ( verbose ) project_show_selected(project);
	if ( verbose > VERB_RESULT ) project_show_class_averages(project);
	
	if ( getselection ) project_show_selection_numbers(project);
	
	if ( show > -1 ) project_show_selected_parameters(project, show);
	if ( tag_show.length() > 0 ) {
		if ( tag_show.contains("micro") ) project_show_mg_parameter(project, tag_show);
		else project_show_part_parameter(project, tag_show);
	}
	
	Bstring			title;
	if ( FOMbins > 1 )
		project_show_fom_histogram(project, FOMbins, FOMmin, FOMmax);
	
	if ( focus ) project_defocus_range(project);
	
	if ( FOMps.length() )
		ps_part_fom_histogram(FOMps, project);
	
	if ( viewps.length() ) {
		title = "Particle view distributions";
		ps_particle_views_origins(viewps, title, symmetry_string, project, -1);
//		ps_particle_phi_theta(viewps, title, project, -1);
	}
	
	if ( mgps.length() ) {
		title = "Particle micrograph distributions";
		ps_mg_particle_positions(mgps, title, project);
	}
	
	if ( write_classes )
		project_write_particle_classes(project);
	
	if ( expfile.length() )
		write_particle_list(expfile, project, exp_flags);
		
	// Write an output parameter format file if a name is given
	if ( project && ( outfile.length() || split == 9 ) ) {
		if ( divide > 1 ) {
			Bstring			filename;
			Bproject*		project2;
			project_divide(project, divide);
			for ( i=1, project2 = project; project2; project2 = project2->next, ++i ) {
				filename = outfile.pre_rev('.') + Bstring(i, "_%02d.") + outfile.post_rev('.');
				write_project(filename, project2, write_flags);
			}
    	} else if ( split < 0 ) {
			project_split_write(outfile, project);
		} else {
			project->split = split;
			write_project(outfile, project, write_flags);
		}
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

