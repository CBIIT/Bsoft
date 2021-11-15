/**
@file	bzfft.cpp
@brief	Disk-based z-line transforms for big reconstructions
@author	Bernard Heymann
@date	Created: 20051221
@date	Modified: 20110804
**/

#include "rwmg.h"
#include "mg_tomo_rec.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bzfft [options] input.img/input.star [input2.img/input2.star...]",
"-----------------------------------------------------------------------",
"Assembles slabs and back-transform z-lines for tomogram reconstruction.",
" ",
"Actions:",
"-ytile 50,25             Tile start and size along y dimension (default 0 and y size).",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-truncate -0.5,1.2       Truncate data to minimum and maximum after rescaling.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling of input images (A/pixel; a single value can be given).",
" ",
"Output:",
"-output new.star         Output parameter file name.",
"-reconstruction new.map  Output image file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			avg(0), std(0);				// Values for rescaling
	double			cutmin(0), cutmax(0);		// Truncation
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	int				ystart(0), ysize(0);		// Y tile start and size
	Bstring			recfile;					// Output image file
	Bstring			outfile;					// Output parameter file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(avg, std) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( std <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "truncate" ) {
			if ( curropt->values(cutmin, cutmax) < 2 )
				cerr << "-truncate: Both min and max must be specified!" << endl;
			if ( cutmin > cutmax ) swap(cutmin, cutmax);
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "ytile" )
			if ( curropt->values(ystart, ysize) < 2 )
				cerr << "-ytile: Both start and size must be specified!" << endl;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "reconstruction" )
			recfile = curropt->filename();
    }
	option_kill(option);
	
	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl;
		bexit(-1);
	}
	
	double				ti = timer_start();
	
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	Bstring				paramfile(*file_list);
	Bproject*			project = NULL;
	Breconstruction*	rec;

	if ( file_type(paramfile) == Micrograph ) {
		project = read_project(file_list);
		string_kill(file_list);
		file_list = NULL;
		for ( rec = project->rec; rec; rec = rec->next )
			string_add(&file_list, rec->frec);
		project_kill(project);
	}

	Bimage*			pmap = NULL;
	
	if ( ysize > 0 ) {
		pmap = img_extract_ytile(file_list, ystart, ysize);
		if ( !pmap ) bexit(-1);
		img_backtransform_z_lines(pmap);
		if ( std > 0 ) pmap->rescale_to_avg_std(avg, std);
		if ( cutmin < cutmax ) pmap->truncate_to_min_max(cutmin, cutmax);
		pmap->change_type(nudatatype);
		if ( sam.volume() > 0 ) pmap->sampling(sam);
		write_img(recfile, pmap, 0);
	} else {
		pmap = img_backtransform_z_on_disk(file_list, recfile, nudatatype, 
					avg, std, cutmin, cutmax);
		if ( !pmap ) bexit(-1);
		if ( sam.volume() > 0 ) pmap->sampling(sam);
	}
	
	Bstring			id("1");
	if ( outfile.length() ) {
		project = read_project(paramfile);
		if ( !project->rec ) rec = reconstruction_add(&project->rec, id);
		else rec = project->rec;
		rec->frec = recfile;
		write_project(outfile, project, 0, 0);
		project_kill(project);
	}

	delete pmap;
	string_kill(file_list);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

