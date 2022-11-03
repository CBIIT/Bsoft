/**
@file	bpartmulti.cpp
@brief	Selection of single particle parameters from multiple files for classification
@author Bernard Heymann
@date	Created: 20010319
@date	Modified: 20220831
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_multiple.h"
#include "mg_particle_select.h"
#include "Matrix.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpartmulti [options] input.star [input.star]",
"---------------------------------------------------",
"Selection and analysis of single particles from multiple parameter files.",
"Selects particles based on correlation coefficient values",
"returning the number of the input file with the best value.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-all                     Reset selection to all particles before other selections (default not).",
"-setasu D3               Set the views to the asymmetric unit.",
"-select var              Select based on: fom (FOM), err (error between first 2),",
"                         rmsd (RMSD between first 2), msd or var (variance).",
"-adjust                  Adjust FOM averages to that of first project.",
"-merge                   Merge multiple files, keeping selected particle parameters.",
"-add                     Add particles from multiple files.",
"-divide 3                Divide the project into a number of subsets of fields.",
"-Compare                 Compare selections.",
"-Statistics              Do selection statistics.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-symmetry D6             Symmetry: Point group identifier.",
"-fom 0.2,2,1             FOM cutoff, index and flag to weigh with defocus (default 0,0,0).",
"-Origin 1.2              Origin deviation to accept (default 10 pixels).",
"-View 13.5               View deviation to accept (default 180 degrees).",
"-Angle 5.8               Rotation angle deviation to accept (default 180 degrees).",
"-Magnification 0.03      Magnification deviation to accept (default 2).",
" ",
"Input:",
"-template input.star     Template parameter file to merge in data blocks from other input files.",
" ",
"Output:",
"-output file.star        Output parameter file.",
" ",
NULL
};

int			main(int argc, char** argv)
{
	// Initializing variables
	int				flags(0);				// Flags to process reconstructions (bit 1), use template (bit 2), reset selections (bit 3)
	Bstring			selection_type;
	bool			merge(0);				// Flag to merge files
	bool			add(0);					// Flag to add particles from multiple files
	int				divide(0);				// Number of output project subsets
	Bstring			symmetry_string;		// No symmetry specified
	Bstring			symmetry_asu;			// Point group to set views to the asymmetric unit
	double			fom_cut(0);				// Minimum FOM cutoff
	int				fom_index(0);			// Index of FOM to use
	int				fom_def_flag(0);		// Flag to adjust the cutoff for defocus
	int				fom_adjust(0);			// Flag to adjust FOM averages to that of first project
	double			origin_dev(10);			// Origin deviation
	double			view_dev(M_PI);			// View vector deviation
	double			angle_dev(M_PI);		// View rotation angle deviation
	double			mag_dev(2);				// Magnification deviation
	int				selection_compare(0);	// Flag to compare selections
	int 			selection_stats(0);		// Flag to do selection statistics
	Bstring			masterfile;				// Template parameter file
	Bstring			outfile;				// Output parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) flags |= 1;
		if ( curropt->tag == "all" ) flags |= 4;
		if ( curropt->tag == "setasu" )
			symmetry_asu = curropt->symmetry_string();
		if ( curropt->tag == "merge" )
			merge = true;
		if ( curropt->tag == "add" )
			add = true;
		if ( curropt->tag == "divide" )
			if ( ( divide = curropt->value.integer() ) < 2 )
				cerr << "-divide: A number of subsets must be specified!" << endl;
		if ( curropt->tag == "select" ) {
			selection_type = curropt->value.lower();
			if ( selection_type.length() < 1 )
				cerr << "-select: A selection type (fom/var/err) must be specified!" << endl;
		}
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "adjust" )
			fom_adjust = 1;
		if ( curropt->tag == "fom" ) {
			if ( curropt->values(fom_cut, fom_index, fom_def_flag) < 1 )
				cerr << "-fom: A cutoff value must be specified" << endl;
			else
				selection_type = "fom";
		}
		if ( curropt->tag == "Origin" ) {
			if ( ( origin_dev = curropt->value.real() ) <= 0 )
				cerr << "-Origin: An acceptable origin deviation must be specified!" << endl;
			else
				if ( origin_dev < 0.01 ) origin_dev = 0.01;
		}
		if ( curropt->tag == "View" ) {
			if ( ( view_dev = curropt->value.real() ) <= 0 )
				cerr << "-View: An acceptable view deviation must be specified!" << endl;
			else {
				view_dev *= M_PI/180.0;
				if ( view_dev < 0.001 ) view_dev = 0.001;
			}
		}
		if ( curropt->tag == "Angle" ) {
			if ( ( angle_dev = curropt->value.real() ) <= 0 )
				cerr << "-Angle: An acceptable angle deviation must be specified!" << endl;
			else {
				angle_dev *= M_PI/180.0;
				if ( angle_dev < 0.001 ) angle_dev = 0.001;
			}
		}
		if ( curropt->tag == "Magnification" ) {
			if ( ( mag_dev = curropt->value.real() ) <= 0 )
				cerr << "-Magnification: An acceptable magnification deviation must be specified!" << endl;
			else
				if ( mag_dev < 0.0001 ) mag_dev = 0.0001;
		}
		if ( curropt->tag == "Compare" )
			selection_compare = 1;
		if ( curropt->tag == "Statistics" )
			selection_stats = 1;
		if ( curropt->tag == "template" )
			masterfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bsymmetry		sym(symmetry_string);
	
	// Read all the parameter files
	Bproject*		first_project = NULL;
	Bproject*		project = NULL;
	Bstring*		file_list = NULL;
	int 			np = 0;
	
	if ( merge || add ) {
		if ( masterfile.length() ) {
			flags |= 2;
			string_add(&file_list, masterfile);
		}
		while ( optind < argc ) string_add(&file_list, argv[optind++]);
		if ( !file_list ) {
			cerr << "Error: No parameter or image files specified!" << endl;
			bexit(-1);
		}
		if ( merge )
			first_project = project_multi_merge(file_list, fom_index, flags);
		else
			first_project = project_multi_add_particles(file_list);
		string_kill(file_list);
		if ( flags & 4 ) part_reset_selection(first_project);
	} else {
		for ( np=0; optind < argc; np++ ) {
			if ( first_project ) {
				project->next = read_project(argv[optind++]);
				project = project->next;
			} else {
				first_project = project = read_project(argv[optind++]);
			}
			if ( flags & 1 ) project->select = 1;
//			project_show_selected(project);
		}
		if ( verbose & VERB_PROCESS )
			cout << np << " parameter files read" << endl << endl;
		if ( np < 1 ) {
			cerr << "Error: No parameter or image files specified!" << endl;
			bexit(-1);
		}
		
		for ( project = first_project; project; project = project->next ) {
			if ( flags & 4 ) part_reset_selection(project, 3);
			if ( symmetry_asu.length() )
				project_set_particle_asu_views(project, symmetry_asu);
		}

		if ( selection_type.contains("err") )
			project_multi_select_low_difference(first_project, first_project->next, sym,
					origin_dev, view_dev, angle_dev, mag_dev);
		
		if ( selection_type.contains("var") )
			project_multi_select_low_variance(first_project, sym,
					origin_dev, view_dev, angle_dev, mag_dev);

		if ( selection_type.contains("msd") )
			project_multi_select_low_rmsd(first_project, first_project->next,
					sym, origin_dev, view_dev, angle_dev, mag_dev,
					1-selection_type.contains("rmsd"));
	
		if ( selection_type.contains("fom") ) {
			if ( fom_adjust )
				project_multi_adjust_FOM(first_project, fom_index);
			project_multi_select_best_FOM(first_project, fom_cut, fom_index, fom_def_flag);
		}
		
		if ( selection_compare )
			project_multi_selection_compare(first_project, first_project->next);
		
		if ( selection_stats )
			project_multi_selection_stats(first_project);
	
		// Delete all but the first project which should have all the results
		project_kill(first_project->next);
		first_project->next = NULL;
	}

	if ( verbose & VERB_PROCESS )
		project_show_selected(first_project);
	
	// Write an output parameter file if a name is given
	Bstring			filename;
    if ( outfile.length() ) {
		if ( divide > 1 ) {
			project_divide(first_project, divide);
			for ( np=1, project = first_project; project; project = project->next, np++ ) {
				filename = outfile.pre_rev('.') + Bstring(np, "_%02d.") + outfile.post_rev('.');
				write_project(filename, project, 0, 0);
			}
    	} else {
			write_project(outfile, first_project, 0, 0);
		}
	}
	
	project_kill(first_project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

