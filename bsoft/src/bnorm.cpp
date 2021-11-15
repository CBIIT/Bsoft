/**
@file	bnorm.cpp
@brief	Program to normalize a set of images based on their histograms.
@author Bernard Heymann
@date	Created: 20030411
@date	Modified: 20160617
**/

#include "rwimg.h"
#include "rwmg.h"
#include "mg_processing.h"
#include "mg_tomography.h"
#include "mg_img_proc.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bnorm [options] input.star/input.img [input2.star/output.img]",
"--------------------------------------------------------------------",
"Normalizes a set of images based on their histograms.",
" ",
"Actions:",
"-replacemaxima 57.5      Replace maxima above the given threshold with surrounding average.",
"-type simple             Normalization type: simple, Gaussian or Poisson.",
"-local 25                Normalization using local average and standard deviation.",
"-truncate -0.5,1.2       Truncate data to minimum and maximum before mass normalization.",
"-rescale 0.5,1.5         Average and standard deviation for output.",
"-images                  Input slices of a single 3D image as 2D images.",
"-slices                  Output 2D images as z-slices of a single 3D image.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
" ",
"Output:",
"-output file.star        Output parameter file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	int 			norm_type(-1);				// Normalization type
	int 			local(0);					// Local normalization kernel size
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			replace_threshold(0);		// Threshold for replacing maxima
	int				replace_flag(0);			// Flag for replacing maxima
	int				setinputZslices(0);			// Interpret 3D slices as separate 2D images
	int				setoutputZslices(0);		// Interpret 3D slices as separate 2D images
	double			cutmin(0), cutmax(0);		// Truncation
	double			avg(0), std(0);				// Average and standard deviation for output
	Bstring			paramfile;					// Output parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "replacemaxima" ) {
			if ( ( replace_threshold = curropt->value.real() ) < 0.0000001 )
				cerr << "-replacemaxima: A threshold must be specified!" << endl;
			else
				replace_flag = 1;
		}
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 's' || curropt->value[0] == 'S' ) norm_type = 0;
			if ( curropt->value[0] == 'g' || curropt->value[0] == 'G' ) norm_type = 1;
			if ( curropt->value[0] == 'p' || curropt->value[0] == 'P' ) norm_type = 2;
		}
		if ( curropt->tag == "local" )
			if ( ( local = curropt->value.integer() ) < 3 )
				cerr << "-local: A kernel size must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "truncate" ) {
			if ( curropt->values(cutmin, cutmax) < 2 )
				cerr << "-truncate: Both min and max must be specified!" << endl;
			if ( cutmin > cutmax ) swap(cutmin, cutmax);
		}
		if ( curropt->tag == "rescale" )
			if ( curropt->values(avg, std) < 2 ) {
				cerr << "-rescale: Both the average and standard deviation must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "images" ) setinputZslices = 1;
		if ( curropt->tag == "slices" ) setoutputZslices = 1;
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	Bproject*		project = NULL;
	Bimage*			p = NULL;
	Bstring			filename = argv[optind];
	Bstring*		file_list = NULL;

	if ( file_type(filename) == Micrograph ) {
		while ( optind < argc ) string_add(&file_list, argv[optind++]);
		project = read_project(file_list);
		string_kill(file_list);
		project_mass_normalize(project, avg, std, norm_type,
				nudatatype, setinputZslices, setoutputZslices,
				cutmin, cutmax, replace_threshold);
	} else {
		p = read_img(filename, 1, -1);
		if ( p == NULL ) bexit(-1);
	
		if ( nudatatype == Unknown_Type )
			nudatatype = p->data_type();
		else if ( nudatatype > p->data_type() )
			p->change_type(nudatatype);
	
		if ( setinputZslices ) 
			if ( p->slices_to_images() < 0 ) bexit(-1);

		if ( cutmin < cutmax ) p->truncate_to_min_max(cutmin, cutmax);
	
		if ( replace_flag )
			p->replace_maxima(replace_threshold);
		
		if ( local > 0 )
			p->normalize_local(local);
		else if ( norm_type >= 0 )
			p->normalize(avg, std, norm_type);

		if ( setoutputZslices ) 
			if ( p->images_to_slices() < 0 ) bexit(-1);
	
		if ( optind < argc ) {
			filename = argv[++optind];
			if ( p->data_type() != nudatatype ) {
				if ( nudatatype == UCharacter ) p->truncate_to_min_max(0, 255);
				p->change_type(nudatatype);
			}
			write_img(filename, p, 0);
		}
		delete p;
	
		file_list = &filename;
		if ( paramfile.length() && filename.length() )
			project = project_create_from_images(file_list, "mg");
	}
	
	if ( project && paramfile.length() ) {
		if ( verbose )
			cout << "Writing " << paramfile << endl;
		write_project(paramfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

