/**
@file	bar.cpp
@brief	Program to do simple arithmetic on images
@author Bernard Heymann
@date	Created: 20040727
@date	Modified: 20110730
**/

/* test*/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#define	MAXOP	16

// Usage assistance
const char* use[] = {
" ",
"Usage: bar [options] input.img output.img",
"-----------------------------------------",
"Does simple arithmetic on images.",
"Operations are done in the order given on the command line.",
" ",
"Actions:",
"-add -5.78               Add a constant to each voxel.",
"-multiply 6.123          Multiply each voxel with a constant.",
"-power 0.75              Raise to the power (must be positive).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			add[MAXOP];					// Additions
	double			multiply[MAXOP];			// Multipliers
	double			power[MAXOP];				// Powers
	
	int				i;
	for ( i=0; i<MAXOP; i++ ) add[i] = multiply[i] = power[i] = 0;
	
	int				n(0), optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "add" ) {
			if ( ( add[n] = curropt->value.real() ) == 0 )
				cerr << "-add: A constant must be specified!" << endl;
			else
				n++;
		}
		if ( curropt->tag == "multiply" ) {
			if ( ( multiply[n] = curropt->value.real() ) ==0 )
				cerr << "-multiply: A constant must be specified!" << endl;
			else
				n++;
		}
		if ( curropt->tag == "power" ) {
			if ( ( power[n] = curropt->value.real() ) == 0 )
				cerr << "-power: A constant must be specified!" << endl;
			else
				n++;
		}
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);
    
	if ( optind >= argc ) bexit(0);
	
	if ( nudatatype > p->data_type() ) p->change_type(nudatatype);
	
	for ( i=0; i<MAXOP && i<n; i++ ) {
		if ( fabs(add[i]) > SMALLFLOAT ) p->add(add[i]);
		if ( fabs(multiply[i]) > SMALLFLOAT ) p->multiply(multiply[i]);
		if ( fabs(power[i]) > SMALLFLOAT ) p->power(power[i]);
	}
	
	p->change_type(nudatatype);
	write_img(argv[optind], p, 0);

	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

