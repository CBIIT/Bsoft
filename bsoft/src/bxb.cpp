/**
@file	bxb.cpp
@brief	Extracting a part of a density map and building a new map
@author Bernard Heymann
@date	Created: 20060411
@date 	Modified: 20190821
**/

#include "rwmodel.h"
#include "model_extract_build.h"
#include "model_select.h"
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
"Usage: bxb [options] in1.star [in2.star ...]",
"--------------------------------------------",
"Extracting a part of a density map and building a new map based on a model.",
"Component orientational and positional refinement can be done as:",
"	1. random modifications over the specified iterations,",
"	2. a contracting grid search when the steps and accuracy are specified.",
" ",
"Actions:",
"-setviews local          Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-invert                  Invert views for selected components.",
"-linklength 58.4         Reference link length to generate links (angstrom, default 1).",
"-type component          Type of element: component (default) or link.",
"-extract 80,50,30        Extract a part of the density with the given size.",
"-build 200,250,150       Build a new map of the given size.",
"-refine template.map     Refine the positions by cross-correlation.",
"-average 2               Average extracted densities: averages per type (use with -extract).",
"-componenttype A,Ver     Names to assign to component types in the template file.",
" ",
"Selections:",
"-all                     Reset selection to all models and components before other selections.",
"-select #Mod1@14         Select models and components.",
"-linkselect 2            Link selection for building a new map (first = 1, default all).",
"-fom 0.25                Select: FOM cutoff.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-datatype u              Force writing of a new data type.",
" ",
"Parameters for building:",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin for extracted or built map",
"                         (origins for extract and build cannot be used together).",
"-weigh                   Weigh by contributions when building a new map (default not).",
"-separate                Separate maps for each component type (default one map).",
" ",
"Common refinement parameters:",
"-resolution 4.5,130      Resolution range for cross-correlation (default 0 - 1e6 angstrom).",
"-shift 2.3               Maximum shift allowed per iteration (angstrom, default 10% map width).",
"-viewangle 1.4           Maximum rotation angle of view vector allowed per iteration (default 0 degrees).",
"-rotationangle 3.5       Maximum rotation angle around view vector allowed per iteration (default 0 degrees).",
"-normal                  Flag to shift only along the component view vector.",
"-lateral                 Flag to shift only perpendicular to the component view vector.",
" ",
"Parameters for random refinement:",
"-iterations 80           Maximum number of refining iterations (default 0 - only positional refinement).",
" ",
"Parameters for grid refinement:",
"-viewstep 2.5            Angular steps in 2 view directions (default 0 degrees).",
"-rotationstep 1.3        Angular step around view (default 0 degrees).",
"-accuracy 0.7            Angular refinement accuracy (default 0.5 degrees).",
" ",
"Output image parameters:",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
"-map file.map            Input map file.",
"-mask mask.mrc           Real space mask to be applied to template and subvolumes.",
"-Mask mask.mrc           Mask for cross-correlation during refinement.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-construct file.pif      Constructed map file.",
"-std file.pif            Standard deviation map file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			reset(0);					// Keep selection as read from file
	Bstring			calc_views;					// Mode to calculate component views
	int				inv_views(0);				// Flag to invert views for selected components
	int				type(1);					// Component = 1, link = 2
	Vector3<long>	extract;					// Size for extracting a part
	Vector3<long>	build;						// Size to build a new map
	int				extract_average(0);			// Average extracted densities
	double			fom_cutoff(0);				// No selection based on FOM
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    			// Map sampling
	Vector3<double>	origin;						// New image origin
    double			hires(0), lores(1e6);		// Limiting resolution range (hires must be > 0 to be set)
	double			max_shift(0);				// Maximum shift allowed, default infinite
	double			max_view_angle(0);			// Maximum rotation angle of view vector allowed
	double			max_rot_angle(0);			// Maximum rotation angle around view vector allowed
	int				shift_flag(0);				// Flag to shift only along normal (1) or perpendicular to it (2)
	double			bias(1);					// Selection bias for first CC
	int				max_iter(0);				// Maximum refining iterations, 0 means only positional
	double			viewstep(0), rotstep(0);	// Angular refinement steps in degrees
	double			accuracy(0.5);				// Accuracy in degrees
	double			linklength(1);				// Default link length
	int				link_select(0);				// Template selection for building
	Bstring			mod_select;					// Model and component selection
	int				flags(0);					// Flags to weigh during building (1) and build separate maps (2)
	double			nuavg(0), nustd(0);			// Rescaling to average and stdev
	Bstring			astring;
	Bstring*		ct_names = NULL;			// List of component names
	Bstring			refine_template;			// Template for refinement
//	Bstring    		comp_select("all");
	Bstring			paramfile;					// Input parameter file
	Bstring			mapfile;					// Input map file
	Bstring			maskfile;					// Input real space mask file
	Bstring			fsmaskfile;					// Input fourier space mask file
	Bstring			modelout;					// Output model file
	Bstring			mapout;						// Output map file
	Bstring			stdfile;					// Output standard deviation map file
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "setviews" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "invert" ) inv_views = 1;
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'c' || curropt->value[0] == 'C' ) type = 1;
			if ( curropt->value[0] == 'l' || curropt->value[0] == 'L' ) type = 2;
		}
		if ( curropt->tag == "extract" )
			extract = curropt->size();
		if ( curropt->tag == "build" )
			build = curropt->size();
		if ( curropt->tag == "average" )
			if ( ( extract_average = curropt->value.integer() ) < 1 )
				cerr << "-average: An integer must be specified" << endl;
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
 		if ( curropt->tag == "linkselect" )
			if ( ( link_select = curropt->value.integer() ) < 0 )
				cerr << "-linkselect: A number must be specified" << endl;
		if ( curropt->tag == "fom" )
			if ( ( fom_cutoff = curropt->value.real() ) < 0.001 )
				cerr << "-fom: A FOM cutoff must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
 		if ( curropt->tag == "shift" )
			if ( ( max_shift = curropt->value.real() ) < 0.001 )
				cerr << "-shift: A maximum shift distance must be specified" << endl;
 		if ( curropt->tag == "viewangle" ) {
			if ( ( max_view_angle = curropt->value.real() ) < 0.001 )
				cerr << "-viewangle: A maximum rotation angle must be specified" << endl;
			else
				max_view_angle *= M_PI/180.0;
		}
 		if ( curropt->tag == "rotationangle" ) {
			if ( ( max_rot_angle = curropt->value.real() ) < 0.001 )
				cerr << "-rotationangle: A maximum rotation angle must be specified" << endl;
			else
				max_rot_angle *= M_PI/180.0;
		}
 		if ( curropt->tag == "bias" )
			if ( ( bias = curropt->value.real() ) < 0.001 )
				cerr << "-bias: A bias must be specified" << endl;
 		if ( curropt->tag == "normal" ) shift_flag = 1;
 		if ( curropt->tag == "lateral" ) shift_flag = 2;
		if ( curropt->tag == "iterations" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "viewstep" ) {
			if ( ( viewstep = curropt->value.real() ) < 0.001 )
				cerr << "-viewstep: An angle must be specified" << endl;
			else
				viewstep *= M_PI/180.0;
		}
		if ( curropt->tag == "rotationstep" ) {
			if ( ( rotstep = curropt->value.real() ) < 0.001 )
				cerr << "-rotationstep: An angle must be specified" << endl;
			else
				rotstep *= M_PI/180.0;
		}
 		if ( curropt->tag == "accuracy" ) {
			if ( ( accuracy = curropt->value.real() ) < 0.001 )
				cerr << "-accuracy: An angle must be specified" << endl;
			else
				accuracy *= M_PI/180.0;
		}
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 1 )
				cerr << "-linklength: A link length must be specified" << endl;
 		if ( curropt->tag == "weigh" ) flags |= 1;
 		if ( curropt->tag == "separate" ) flags |= 2;
		if ( curropt->tag == "componenttype" ) {
			astring = curropt->value;
			if ( astring.length() ) ct_names = astring.split(",");
		}
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "refine" )
			refine_template = curropt->filename();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			fsmaskfile = curropt->filename();
		if ( curropt->tag == "output" )
			modelout = curropt->filename();
		if ( curropt->tag == "construct" )
			mapout = curropt->filename();
		if ( curropt->tag == "std" )
			stdfile = curropt->filename();
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

	Bmodel*			model = read_model(file_list, paramfile);		
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}

	if ( reset ) models_process(model, model_reset_selection);

	if ( mod_select.length() ) model_select(model, mod_select);

	if ( fom_cutoff > 0 ) model_fom_deselect(model, fom_cutoff);

	if ( calc_views.length() ) model_calculate_views(model, calc_views);

	if ( inv_views ) model_invert_views(model);
	
	Bimage*			p = NULL;
	Bimage*			pn = NULL;
	Bimage*			pt = NULL;
	Bimage*			pm = NULL;
	Bimage*			pfsm = NULL;

	if ( mapfile.length() ) {
		if ( type == 2 ) {
			p = read_img(mapfile, 1, -1);
			if ( p == NULL ) {
				cerr << "Error: No image file read!" << endl;
				bexit(-1);
			}
		}
	}

	if ( type == 2 && !model->link )
		model_link_list_generate(model, linklength);

	if ( refine_template.length() ) {
		if ( mapfile.length() ) model->mapfile(mapfile.str());
		pt = read_img(refine_template, 1, -1);
		if ( pt == NULL ) {
			cerr << "Error: No image file read!" << endl;
			bexit(-1);
		}
		if ( maskfile.length() ) pm = read_img(maskfile, 1, -1);
		if ( fsmaskfile.length() ) pfsm = read_img(fsmaskfile, 1, -1);
		if ( type == 1 ) {
			model_refine_components(model, ct_names, pt, pm, pfsm, max_iter,
				viewstep, rotstep, hires, lores, accuracy,
				max_shift, max_view_angle, max_rot_angle, shift_flag);
		} else if ( type == 2 ) {
			model_refine_link_positions(model, pt, pm, pfsm, hires, lores,
				max_shift, shift_flag, bias);
		}
//		if ( pt->n > 1 ) model_type_from_selection(model, ct_names, pt->file_name());
		delete pt;
		delete pm;
		delete pfsm;
	} else if ( ct_names ) {
		astring = 0;
//		model_type_from_selection(model, ct_names, astring);
	}
	
	if ( extract.volume() ) {
		if ( mapfile.length() ) model->mapfile(mapfile.str());
		if ( type == 1 ) {
			if ( extract_average > 0 ) {
				pn = model_average_component_density(model, extract, origin, extract_average);
				if ( mapout.length() ) model_set_comptype_filenames(model, mapout);
			} else {
				pn = model_extract_component_densities(model, extract, origin);
			}
		} else if ( type == 2 ) {
			pn = model_average_link_density(model, extract, origin);
		}
		origin = 0;
	}
	
	if ( build.volume() ) {
		if ( pn ) {
			delete p;
			p = pn;
		}
		if ( type == 1 ) {
			if ( mapfile.length() ) model_set_comptype_filenames(model, mapfile);
			pn = model_build_from_component_density(model, build, origin, flags);
		} else if ( type == 2 ) {
			pn = model_build_from_link_density(model, mapfile, build, origin, link_select, flags);
		}
		if ( mapout.length() ) model->mapfile(mapout.str());
	}

	if ( modelout.length() ) {
		write_model(modelout, model);
	}

	if ( pn ) pn->calculate_background();


    // Write an output file if a file name is given
    if ( mapout.length() ) {
		if ( sam.volume() > 0 ) pn->sampling(sam);
		if ( nustd > 0 ) pn->rescale_to_avg_std(nuavg, nustd);
		pn->change_type(nudatatype);
    	write_img(mapout, pn, 0);
	}
	
    // Write an standard deviation output file if a file name is given
	if ( stdfile.length() && pn->next )
		write_img(stdfile, pn->next, 0);

	model_kill(model);
	delete p;
	delete pn;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

