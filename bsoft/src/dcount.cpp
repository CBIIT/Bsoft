/**
@file	dcount.cpp 
@brief	A program to determine the occupancy of components in models. 
@author Daniel Nemecek and Bernard Heymann
@date	Created: 20091202 (DN)
@date	Modified: 20110729 (BH)
@date	Updated 20110507 (DN) 
      - determine density threshold from unmasked particles
      - mask particles before refinement of locations and determination of occupancy 
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "model_occupancy.h"
#include "model_select.h"
#include "model_util.h"
#include "ps_model.h"
#include "timer.h"
#include "options.h"
#include "utilities.h"


// Declaration of global variables
extern int      verbose;                // Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: dcount [options] input.star [input2.star...]",
"---------------------------------------------------",
"Determines the occupancy of components in the maps associated with models.",
"The maps must have positive density for occupancy determination (use -invert option).",
" ",
"Actions:",
"-invert                  Invert density for determination of the threshold.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-componentradius 8.4     Set radius for all components.",
"-shift 40                Maximum shift allowed (angstrom, default infinite).",
"-resolution 50,500       Resolution range for cross-correlation (default 0 - 1e6 angstrom).",
"-rho 0.8                 Protein density (default 0.81 Da/A3, use with -mass).",
"-mass 12.645M            Molecular mass of the particle (determines the threshold.",
"-cutoff 0.3              Cutoff for voxels above threshold for determining occupancy (default 0.5),",
"-bins 50                 Number of bins for histogram of average intensities at the marker locations in all images.",
"-fit 2                   Number of fitted binomial distributions to the distribution of occupancy.",
" ",
"Input:",
"-model model.star        Input reference model of expected locations in the particle.",
"-template templ.pif      Template map for cross-correlation refinement of component locations.",
"-Mask shell.map          Binary mask to hide denisity in some regions of the map before determination of occupancy.",
" ",
"Output:",
"-output model.star       Model output file.",
"-split 3                 Split models into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": model ID's are used as file names.",
"-Postscript plots.ps     Postscript output file.",
" ",
"Examples:",
"dcount -v 7 -model model_template.star -comp 15.6 -template map_template.pif -shift 55 ",
"	-resol 50,500 -invert -mass 12.645M -cutoff 50 -bins 20 -fit 2 ",
"	-output model_occ.star -Postscript model_occ.ps model.star",
" ",
NULL
};


int             main(int argc, char **argv) 
{
	// Initialize variables
	double			comprad(0);			// Component radius
	int				bins(256);          // Number of bins in output histogram of average intensities
	int				invert_flag(0);     // Invertion of intensity for threshold determination
	int				nfit(1);            // Number of fitted binomial distributions
	double			mol_weight(0);      // Molecular weight for threshold finding
	double			hires(0);			// Limiting resolution range (hires must be > 0 to be set)
	double			lores(1e6);         // Limiting resolution range (hires must be > 0 to be set)
	double			max_shift(0);       // Maximum shift allowed, default infinite
	double			cutoff(0.5);		// Limits for cutoff determination
	double			rho(RHO);           // Protein density in Da/A3
	Bstring			modelin;            // Input STAR file with a model of expected locations
	Bstring			templfile;          // Filename of the template file for cross correlation
	Bstring			maskfile;           // Filename of the input mask
	Bstring			PSname;             // Name of the output PostScript file
	Bstring			outmod;             // Rootname of output generized models
	int				split(0);			// Sets output of multiple single-model files
        
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	for ( curropt = option; curropt; curropt = curropt->next ) {
        if ( curropt->tag == "invert" ) invert_flag = 1;
        if ( curropt->tag == "fit" )
            if ( ( nfit = curropt->value.integer() ) < 1 )
				cerr << "-fit: Number of fitted distributions must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
        if ( curropt->tag == "model" )	
			modelin = curropt->filename();
       if ( curropt->tag == "Mask" )	
			maskfile = curropt->filename();
        if ( curropt->tag == "template" )	
			templfile = curropt->filename();
        if ( curropt->tag == "shift" )
			if ( ( max_shift = curropt->value.real() ) < 1 )
				cerr << "-shift: A maximum shift distance must be specified" << endl;
        if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
        if ( curropt->tag == "rho" )
			if ( ( rho = curropt->value.real() ) < 1 )
				cerr << "-rho: A protein density value must be specified!" << endl;
        if ( curropt->tag == "mass" ) {
			if ( ( mol_weight = curropt->value.real() ) < 1 )
				cerr << "-mass: A molecular weight must be specified!" << endl;
            else {
				if ( curropt->value.contains("k") || curropt->value.contains("K") )	mol_weight *= 1e3;
                if ( curropt->value.contains("m") || curropt->value.contains("M") )	mol_weight *= 1e6;
                if ( curropt->value.contains("g") || curropt->value.contains("G") )	mol_weight *= 1e9;
			}
		}
        if ( curropt->tag == "cutoff" ) {
			if ( ( cutoff = curropt->value.real() ) < 1 )
				cerr << "-cuttof: A cutoff fraction must be specified!" << endl;
            else
				if ( cutoff >= 1 ) cutoff /= 100;
		}
		if ( curropt->tag == "bins" )
			if ( ( bins = curropt->value.integer() ) < 1 )
				cerr << "-bins: Number of bins must be specified!" << endl;
		if ( curropt->tag == "output" ) {
			outmod = curropt->value;
			if ( outmod.length() < 1 )
				cerr << "-output: A filename must be specified!" << endl;
		}
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
		}
		if ( curropt->tag == "Postscript" ) {
			PSname = curropt->value;
			if ( PSname.length() < 1 )
				cerr << "-Postscript: A filename for the output PostScript file must be specified!" << endl;
		}
	}
    option_kill(option);
        
    double			ti = timer_start();
        
	// Read the input models
	Bstring*		file_list = NULL;
        
    while ( optind < argc ) string_add(&file_list, argv[optind++]);
    if ( !file_list ) {
		cerr << "Error: No model or image files specified!" << endl;
		bexit(-1);
	}

    Bmodel*         model = read_model(file_list);
	string_kill(file_list);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
    Bmodel*         modref = NULL;
	
	// Read a reference model and replace all components in the input models
	if ( modelin.length() ) {
		modref = read_model(modelin);
		if ( !modref ) {
			cerr << "Error: No model was read!\n" << endl;
			bexit(-1);
		}
		model_replace_components(model, modref);
		model_kill(modref);
	}

	if ( comprad > 0 ) models_process(model, comprad, model_set_component_radius);

//    Bimage*         p = NULL;
    Bimage*         pt = NULL;
	Bimage*         pmask = NULL;
	Bmodel*			mp;
	
	// Read a common mask for individual images
    if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, -1);
		if ( pmask == NULL ) {
			cerr << "Error: No mask image was read!\n" << endl;
			bexit(-1);
		}
	}

	// Refine component locations if the template map is defined
    if ( templfile.length() ) {
		pt = read_img(templfile, 1, -1);
		if ( pt == NULL ) {
			cerr << "Error: No template image was read!\n" << endl;
			bexit(-1);
		}
		for ( i=0, mp=model; mp; mp = mp->next, i++ )
			model_refine_comp_for_occupancy(mp, pmask, pt, NULL, hires, lores, max_shift);
	}

	model_occupancy(model, pmask, mol_weight, rho, cutoff, invert_flag);
	
	// Occupancy analysis and binomial fit
	long			ncomp;
	double			R(0);
    vector<double>	prob(2*nfit);         // array for parameters of fitted B(n,p)

	vector<double>	distrib = model_occupancy_distribution(model, cutoff, nfit, ncomp, prob, R);
	
	if ( PSname.length() )
		ps_model_occupancy(model, cutoff, bins, nfit, distrib, prob, R, PSname);
	
	model_selection_stats(model);
	
	// Write the model 
    if ( model && ( outmod.length() || split == 9 ) )
		write_model(outmod, model);

    delete pt;
    delete pmask;
    model_kill(model);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

