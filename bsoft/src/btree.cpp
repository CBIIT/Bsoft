/**
@file	btree.cpp
@brief	Program to process trees.
@author Bernard Heymann
@date	Created: 20010722
@date	Modified: 20030608
**/

#include "rwtree.h"
#include "ps_tree.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: btree [options] tree.ph tree.ps",
"--------------------------------------",
"Draws phylogenetic trees as postscript.",
"Handles phylip format trees as input.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-scale 350               Scale for postscript output (default 200).",
"-rotate 110              Rotate for postscript output (default 0 degrees).",
" ",
NULL
};

int 		main(int argc, char** argv)
{
	// Initialize variables
	double			scale = 200;	// Postscript output scale
	double			angle(0);		// Angle to rotate the tree
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.0001 )
				cerr << "-scale: A scale must be specified" << endl;
		if ( curropt->tag == "rotate" ) {
			if ( ( angle = curropt->value.real() ) < 1 )
				cerr << "-rotate: An angle must be specified" << endl;
			else
				if ( angle > M_PI || angle < -M_PI ) angle *= M_PI/180;
        }
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bnode*		root = read_tree(argv[optind++]);
	
	if ( angle != 0 ) tree_rotate(root, angle);
	
	Bstring		output;
	if ( optind < argc ) {
		output = argv[optind];
		ps_draw_tree(output, root, scale);
	}
	
	tree_delete_node(root);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

