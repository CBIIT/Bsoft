/**
@file	bar.cpp
@brief	Program to do simple arithmetic on images
@author Bernard Heymann
@date	Created: 20040727
@date	Modified: 20220209
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bar [options] input.img output.img",
"-----------------------------------------",
"Does simple arithmetic on images.",
"Operations are done in the order given on the command line.",
" ",
"Actions:",
"-invert                  Invert the image values.",
"-absolute                Convert to absolute values.",
"-add -5.78               Add a constant to each voxel.",
"-multiply 6.123          Multiply each voxel with a constant.",
"-power 0.75              Raise to the power (must be positive).",
"-phaseadd 2.3            Add an angle to a phase image.",
"-sine                    Calculate the sine of a phase image.",
"-cosine                  Calculate the cosine of a phase image.",
"-tangent                 Calculate the tangent of a phase image.",
"-arcsine                 Calculate the arcsine of an image.",
"-arccosine               Calculate the arccosine of an image.",
"-arctangent              Calculate the arctangent of an image.",
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
	vector<int>		ops;			// Operations
	vector<double>	op_param;		// Operational parameter
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "absolute" ) {
			ops.push_back(1);
			op_param.push_back(0);
		}
		if ( curropt->tag == "add" ) {
			ops.push_back(2);
			op_param.push_back(curropt->value.real());
		}
		if ( curropt->tag == "multiply" ) {
			ops.push_back(3);
			op_param.push_back(curropt->value.real());
		}
		if ( curropt->tag == "power" ) {
			ops.push_back(4);
			op_param.push_back(curropt->value.real());
		}
		if ( curropt->tag == "phaseadd" ) {
			ops.push_back(5);
			op_param.push_back(curropt->angle());
//			cout << "Angle = " << op_param.back() << endl;
		}
		if ( curropt->tag == "sine" ) {
			ops.push_back(6);
			op_param.push_back(0);
		}
		if ( curropt->tag == "cosine" ) {
			ops.push_back(7);
			op_param.push_back(0);
		}
		if ( curropt->tag == "tangent" ) {
			ops.push_back(8);
			op_param.push_back(0);
		}
		if ( curropt->tag == "arcsine" ) {
			ops.push_back(9);
			op_param.push_back(0);
		}
		if ( curropt->tag == "arccosine" ) {
			ops.push_back(10);
			op_param.push_back(0);
		}
		if ( curropt->tag == "arctangent" ) {
			ops.push_back(11);
			op_param.push_back(0);
		}
		if ( curropt->tag == "invert" ) {
			ops.push_back(12);
			op_param.push_back(0);
		}
   }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);
    
	if ( optind >= argc ) bexit(0);
	
	if ( nudatatype > p->data_type() ) p->change_type(nudatatype);
	
	for ( long i=0; i<ops.size(); ++i ) {
		switch ( ops[i] ) {
			case 1: p->absolute(); break;
			case 2: p->add(op_param[i]); break;
			case 3: p->multiply(op_param[i]); break;
			case 4: p->power(op_param[i]); break;
			case 5: p->phase_add(op_param[i]); break;
			case 6: p->sine(); break;
			case 7: p->cosine(); break;
			case 8: p->tangent(); break;
			case 9: p->arcsine(); break;
			case 10: p->arccosine(); break;
			case 11: p->arctangent(); break;
			case 12: p->invert(); break;
		}
	}
	
	p->change_type(nudatatype);
	write_img(argv[optind], p, 0);

	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

