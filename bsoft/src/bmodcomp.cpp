/**
@file	bmodcomp.cpp
@brief	A tool to compare polyhedra.
@author Bernard Heymann
@date	Created: 20080102
@date 	Modified: 20190125
**/

#include "rwmodel.h"
#include "model_select.h"
#include "model_poly.h"
#include "model_compare.h"
#include "model_util.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodcomp [options] in.star",
"---------------------------------",
"Compares models internally or to reference models.",
" ",
"Actions:",
"-all                     Select all models for comparison.",
"-distance                Calculate a distance matrix.",
"-unknown                 Select unknown models for comparison.",
"-closed order,3          Select based on valency (valency,<n>) or polygon order (order,<n>).",
"-consensus 12.5          Calculate a consensus model: components are considered the same if within the given distance.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Input:",
"-parameters parm.star    Molecular parameter file (default atom_prop.star).",
"-reference ref.star      File with reference models.",
" ",
"Output:",
"-output new.star         Output model file.",
"-writereference ref.star Reference output model file.",
"-matrix file.mat         Distance matrix.",
"-image file.map          Distance matrix as an image.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
 	int				all(0);					// Flag to select all models
 	int				calc_dist(0);			// Flag to calculate a distance matrix
 	int				unknown(0);				// Flag to select unknown models
	int				closure_rule(0);		// Closure rule: 1=valency, 2=order
	int				val_order(0);			// Valency or order - depending on rule
	double			distance(0);			// Cutoff distance for consensus model
	Bstring    		atom_select("ALL");
	Bstring			paramfile;				// Use default parameter file
	Bstring			reffile;				// Reference model file
	Bstring			outfile;				// Output model file
	Bstring			refoutfile;				// Reference output model file
	Bstring			matfile;				// Output matrix file
	Bstring			imgfile;				// Output image file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "distance" ) calc_dist = 1;
		if ( curropt->tag == "unknown" ) unknown = 1;
		if ( curropt->tag == "closed" ) {
			if ( curropt->value[0] == 'v' ) closure_rule = 1;
			if ( curropt->value[0] == 'o' ) closure_rule = 2;
			if ( curropt->value.contains(",") )
				val_order = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "consensus" )
			if ( ( distance = curropt->real() ) < 1 )
				cerr << "-consensus: A cutoff distance must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "writereference" )
			refoutfile = curropt->filename();
		if ( curropt->tag == "matrix" )
			matfile = curropt->filename();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
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

	if ( all ) models_process(model, model_reset_selection);
	
	Bmodel*		refmodel = NULL;		
	if ( reffile.length() )
		refmodel = read_model(reffile, paramfile);		

	if ( !model->poly ) model_poly_generate(model);

	if ( calc_dist ) {
		Matrix		mat;
		if ( refmodel )
			mat = model_distance_matrix(model, refmodel);
		else
			mat = model_distance_matrix(model, model);
//		for ( mp = model; mp; mp = mp->next ) {
//			mat = model_distance_matrix(mp, 0);
//			cout << mp->identifier() << endl << mat << endl;
//			cout << mp->identifier() << tab << mp->mapfile() << tab << mat[0][1] << endl;
//		}
		if ( imgfile.length() && mat.rows() ) {
			Bimage*		pimg = new Bimage(mat, 1);
//			pimg->change_type(nudatatype);
			write_img(imgfile, pimg, 0);
			delete pimg;
		}
	}
	
	if ( unknown ) {
		model_select_unknowns(model);
		refmodel = model;
	}
	
	if ( closure_rule ) model_select_closed(model, closure_rule, val_order);
	
	if ( refmodel ) {
		if ( model->poly )
			model_poly_compare(model, refmodel);
		else
			model_compare(model, refmodel);
	}
	
	Bmodel*		numod = NULL;
	if ( distance > 0 ) {
		numod = models_consensus(model, distance);
		model_kill(model);
		model = numod;
	}
	
	if ( outfile.length() ) {
		write_model(outfile, model);
	}

	if ( refoutfile.length() ) {
		write_model(refoutfile, refmodel);
	}

	if ( model != refmodel ) model_kill(model);
	model_kill(refmodel);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

