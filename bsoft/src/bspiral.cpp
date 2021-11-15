/**
@file	bspiral.cpp
@brief	A tool to generate polyhedra using the spiral algorithm.
@author Bernard Heymann
@date	Created: 20071127
@date 	Modified: 20080408
**/

#include "rwmodel.h"
#include "model_poly.h"
#include "model_poly_spiral.h"
#include "model_transform.h"
#include "model_links.h"
#include "model_util.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bspiral [options] out.star",
"---------------------------------",
"Generates a polyhedron using the spiral algorithm.",
" ",
"Actions:",
"-analyze                 Analyze the polyhedron.",
"-enantiomorph            Generate the enantiomorph of the polyhedron.",
"-sequence 5566665566565  Sequence of polygons to use in stead of permutations.",
"-cone 14,30,25           Polygons to generate for the tip, body and base of a cone.",
"-lozenge 2,7             Unit lengths for the top and body of a lozenge.",
"-coffin 3,5,2            Unit lengths for the top, body and base of a coffin.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-vertices 36             Number of vertices (default 20).",
"-linklength 150          Vertex separation (default 10 A).",
"-requirements 0          Polyhedron acceptance requirements (default 1 = strict).",
" ",
"Input:",
"-parameters parm.star    Molecular parameter file (default atom_prop.star).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int				analyze(0);				// Flag to analyze the structure
	int				enantiomorph(0);		// Flag to generate enantiomorphs
	int				vertices(20);			// Number of vertices
	int				valence(3);				// Vertex valence
	double			linklength(1);			// Unit length separating vertices/triskelions
	int				use_one(0);				// Flag to indicate using only one sequence
	int				conetip(0);
	int				conebody(0);
	int				conebase(0);
	int				loztop(0);
	int				lozbody(0);
	int				coftop(0);
	int				cofbody(0);
	int				cofbase(0);
	int				requirements(1);		// Polyhedron acceptance requirements
	Bstring			seq("555555555555");	// Initial sequence
	Bstring			paramfile;				// Use default parameter file
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "analyze" ) analyze = 1;
		if ( curropt->tag == "enantiomorph" ) enantiomorph = 1;
		if ( curropt->tag == "vertices" ) {
			if ( ( vertices = curropt->value.integer() ) < 1 )
				cerr << "-vertices: The number of vertices must be specified!" << endl;
			else {
				use_one = 0;
			}
		}
		if ( curropt->tag == "sequence" ) {
			seq = curropt->value;
			vertices = 20 + 2*(seq.length() - 12);
			use_one = 1;
		}
		if ( curropt->tag == "cone" )
			if ( curropt->values(conetip, conebody, conebase) < 3 )
				cerr << "-cone: Three integers must be specified!" << endl;
		if ( curropt->tag == "lozenge" )
			if ( curropt->values(loztop, lozbody) < 2 )
				cerr << "-lozenge: Two integers must be specified!" << endl;
		if ( curropt->tag == "coffin" )
			if ( curropt->values(coftop, cofbody, cofbase) < 3 )
				cerr << "-coffin: Three integers must be specified!" << endl;
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 1 )
				cerr << "-linklength: The link length separating vertices must be specified!" << endl;
		if ( curropt->tag == "requirements" )
			if ( ( requirements = curropt->value.integer() ) < 1 )
				cerr << "-requirements: An integer must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	mkdir("temp", O_CREAT );
	chmod("temp", 0755);
	
	Bmodel*		model = NULL;
	
	if ( use_one ) {
		model = model_poly_gen_sequence(seq, valence, enantiomorph, requirements, 0);
		if ( model )
			cout << model->identifier() << ": " << seq << ": " << model->symmetry() << endl;
	} else if ( conetip > 0 && conebody > 0 && conebase > 0 ) {
		model = model_poly_gen_cone(conetip, conebody, conebase, valence, enantiomorph, requirements);
	} else if ( loztop > 0 ) {
		model = model_poly_gen_lozenge(loztop, lozbody, valence, enantiomorph, requirements);
	} else if ( coftop > 0 && cofbody > 0 && cofbase > 0 ) {
//		model = model_poly_gen_coffin(coftop, cofbody, cofbase, valence, enantiomorph, requirements);
//		model = model_poly_gen_coffin_loose(coftop, cofbody, cofbase, valence, enantiomorph, requirements);
		model = model_poly_gen_coffin_jiggle(coftop, cofbody, cofbase, valence, enantiomorph, requirements);
	} else {
		model = model_poly_gen_permutations(vertices, valence, enantiomorph);
	}

	i = count_list((char *)model);
	
	if ( !i )
		if ( verbose ) cout << "No model generated!" << endl;
	
	Vector3<double>	scale(linklength,linklength,linklength);
	Vector3<double>	origin;
	model_scale(model, scale, origin);

	models_process(model, 0.1*linklength, model_set_component_radius);

	models_process(model, 0.1*linklength, model_set_link_radius);
	
	if ( analyze ) model_poly_analyze(model);

	Bstring			filename;
	if ( optind < argc ) {
		filename = argv[optind];
		write_model(filename, model);
	}

	model_kill(model);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(i);
}

