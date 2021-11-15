/**
@file	bmodsel.cpp
@brief	Manipulates selections from sets of components as solutions of fits
@author Bernard Heymann
@date	Created: 20060908
@date 	Modified: 20210318
**/

#include "rwmolecule.h"
#include "rwmodel.h"
#include "rwimg.h"
#include "model_select.h"
#include "model_util.h"
#include "ps_model.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodsel [options] in1.star [in2.star...]",
"-----------------------------------------------",
"Selects and manipulates models and components.",
"Selection syntax:",
"	#model_id",
"	^model_type_id",
"	@component_id",
"	%component_type_id",
"	#model_id@component_id",
"	#model_id%component_type_id",
"	^model_type_id@component_id",
"	^model_type_id%component_type_id",
"	wild cards allowed: .",
" ",
"Actions for preparation:",
"-all                     Reset selection to all components before other selections.",
"-merge                   Merge models in different files rather than concatenate.",
"-invert                  Inverts selection.",
" ",
"Selections:",
"-select #Mod1@14         Select models and components.",
"-group 3                 Select components with this selection flag number.",
"-number 24,36            Select models based on the number of components.",
"-sets 20,1               Generate sets of this size of selected components, ",
"                         with a flag not to select across model boundaries.",
"-closed order,3          Select models based on valency (valency,<n>) or polygon order (order,<n>).",
"-fullerene               Select fullerene type models.",
"-nonfullerene            Select non-fullerene type models.",
"-fom 0.25                Select: FOM cutoff.",
"-fraction 0.6            Select: FOM/FOMmax ratio cutoff.",
"-rank 5                  Select: Rank components into a number of groups.",
"-prune fom               Type of pruning: simple, fom, similar, large.",
"-shell 7,18,23,5.4,-3.5  Select within a shell within given radii from the center (default {0,0,0}).",
"-delete                  Delete non-selected models and components.",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components after other operations.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-distance 12.8           Distance criterion for pruning or overlap (default 10 angstrom).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
"-overlap file.star       Deselect components that overlap with this model (use with -distance).",
"-mask mask.pif           Select components only within the mask.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-FOM histogram.ps        Postscript output with component FOM distributions.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			all(0);					 	// Keep selection as read from file
	int 			reset(0);					// Keep selection as ouput
	int				merge(0);					// Flag to merge models rather than concatenate
	int 			invert(0);					// Flag to invert selection
	Bstring			mod_select;					// Model and component selection
	long			group(0);					// Group selection
	long			ncomp_min(0), ncomp_max(0);	// Minimum and maximum number of components to select for
	long			set_size(0);				// Size of sets
	int				set_flag(0);				// Flag to keep set within model
	int				closure_rule(0);			// Closure rule: 1=valency, 2=order
	int				val_order(0);				// Valency or order - depending on rule
	int				fullerene(0);				// Flag to select fullerenes
	int				nonfullerene(0);			// Flag to select non-fullerenes
	double			minrad(0), maxrad(0);		// Shell radii
	Vector3<double> shell_center;				// Shell center
	double			fom_cutoff(0);				// No selection based on FOM
	double			fom_fraction(0);			// No selection based on FOM/FOMmax fraction
	int				rank(0);					// Ranking selection into a number of groups
	int				prune_type(0);				// Prune type: 1=simple, 2=fom, 3=fit, 4=large
	Bstring			prune_type_string("none");	// Prune type string
	double			distance(10);				// Prune distance between components criterion
	int				mod_comp_del(0);			// Flag to delete non-selected models and components
	Bstring    		atom_select("all");
	Bstring			paramfile;					// Input parameter file name
	Bstring			ovlapfile;					// Reference file name for deselecting overlaps
	Bstring			maskfile;					// Mask for selection
	Bstring			outfile;					// Output parameter file name
	Bstring			FOMps = NULL;				// Output postscript file name for FOM histogram
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "invert" ) invert = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "group" )
			if ( ( group = curropt->integer() ) < 1 )
				cerr << "-group: A selection number must be specified!" << endl;
		if ( curropt->tag == "number" )
			if ( curropt->values(ncomp_min, ncomp_max) < 1 )
				cerr << "-number: A number of components must be specified!" << endl;
		if ( curropt->tag == "delete" ) mod_comp_del = 1;
		if ( curropt->tag == "sets" )
			if ( curropt->values(set_size, set_flag) < 1 )
				cerr << "-sets: A set size must be specified!" << endl;
		if ( curropt->tag == "closed" ) {
			if ( curropt->value[0] == 'v' ) closure_rule = 1;
			if ( curropt->value[0] == 'o' ) closure_rule = 2;
			if ( curropt->value.contains(",") )
				val_order = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "fullerene" ) fullerene = 1;
		if ( curropt->tag == "nonfullerene" ) nonfullerene = 1;
		if ( curropt->tag == "shell" )
			if ( curropt->values(minrad, maxrad, shell_center[0], shell_center[1], shell_center[2]) < 1 )
				cerr << "-shell: A minimum radius must be specified!" << endl;
		if ( curropt->tag == "fom" )
			if ( ( fom_cutoff = curropt->value.real() ) < 0.0001 )
				cerr << "-fom: A FOM cutoff must be specified!" << endl;
		if ( curropt->tag == "fraction" )
			if ( ( fom_fraction = curropt->value.real() ) < 0.0001 )
				cerr << "-fraction: A fraction must be specified!" << endl;
		if ( curropt->tag == "rank" )
			if ( ( rank = curropt->value.integer() ) < 1 )
				cerr << "-rank: A number must be specified!" << endl;
		if ( curropt->tag == "prune" ) {
			prune_type_string = curropt->value;
			if ( prune_type_string.length() < 1 )
				cerr << "-prune: A prune type must be specified!" << endl;
			else {
				if ( prune_type_string.contains("simp") ) prune_type = 1;
				if ( prune_type_string.contains("fom") ) prune_type = 2;
				if ( prune_type_string.contains("fit") ) prune_type = 3;
				if ( prune_type_string.contains("simi") ) prune_type = 4;
				if ( prune_type_string.contains("larg") ) prune_type = 5;
			}
		}
		if ( curropt->tag == "distance" )
			if ( ( distance = curropt->value.real() ) < 1 )
				cerr << "-distance: The overlap distance must be specified!" << endl;
/*		if ( curropt->tag == "minradius" )
			if ( ( minrad = curropt->value.real() ) < 1 )
				cerr << "-minradius: The minimum allowed distance from the origin must be specified!" << endl;
		if ( curropt->tag == "maxradius" )
			if ( ( maxrad = curropt->value.real() ) < 1 )
				cerr << "-maxradius: The maximum allowed distance from the origin must be specified!" << endl;*/
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "overlap" )
			ovlapfile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "FOM" )
			FOMps = curropt->filename();
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bmodel*		model = read_model(file_list);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( all ) models_process(model, model_reset_selection);
	
	if ( merge ) model_merge(model);

	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( group ) model_select(model, group);
	
	if ( ncomp_min ) model_select_number_of_components(model, ncomp_min, ncomp_max);

	if ( minrad > 0 || maxrad > 0 )
		model_select_within_shell(model, shell_center, minrad, maxrad);
	
	if ( closure_rule ) model_select_closed(model, closure_rule, val_order);

	if ( fullerene ) model_select_fullerene(model);
	
	if ( nonfullerene ) model_select_non_fullerene(model);
	
	if ( set_size ) model_select_sets(model, set_size, set_flag);
	
	if ( fom_cutoff > 0 ) model_fom_deselect(model, fom_cutoff);

	if ( fom_fraction > 0 ) model_fom_max_fraction_deselect(model, fom_fraction);
	
	if ( rank > 1 ) model_fom_ranking(model, rank);

//	if ( minrad || maxrad ) models_radius_deselect(model, minrad, maxrad);

	if ( prune_type > 4 ) models_process(model, distance, model_prune_large);
	else if ( prune_type == 4 ) models_process(model, model_prune_similar);
	else if ( prune_type == 3 ) models_process(model, distance, model_prune_fit);
	else if ( prune_type == 2 ) models_process(model, distance, model_prune_fom);
	else if ( prune_type == 1 ) models_process(model, distance, model_prune_simple);
	
	if ( ovlapfile.length() )
		model_find_overlap(model, ovlapfile, distance);
	
	if ( maskfile.length() ) {
		Bimage*		pmask = read_img(maskfile, 1, 0);
		model_select_in_mask(model, pmask);
		delete pmask;
	}
	
	if ( invert )
		models_process(model, model_invert_selection);
	
	if ( FOMps.length() )
		ps_model_fom_histogram(FOMps, model);
	
	if ( mod_comp_del ) model_delete_non_selected(&model);
	
	if ( reset ) models_process(model, model_reset_selection);
	
	model_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	model_kill(model);
		
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

