/**
@file	bmgalign.cpp
@brief	Aligns micrographs specified in parameter files.
@author Bernard Heymann and Samuel Payne
@date	Created: 20000505
@date	Modified: 20150730 (BH)
**/

#include "rwimg.h"
#include "mg_processing.h"
#include "rwmg.h"
#include "mg_pick.h"
#include "mg_align.h"
#include "Matrix.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmgalign [options] input.star [input2.star...]",
"-----------------------------------------------------",
"Aligns micrographs within sets such as focal series specified in parameter files.",
" ",
"Actions:",
"-align mic               Alignment type: par = particle coordinates, mic = micrograph images, fea = micrograph features.",
"-add                     Add micrograph coordinates and particle origins (default not).",
"-filterextremes          Filter micrograph extremes before aligning.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-number 3                Number of micrographs per field-of-view (default 1).",
"-fieldID dbgd            Set the field ID (multiple micrographs will get numbered ID's).",
"-reference 2             Reference micrograph in a series (default 1=first micrograph).",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-correlate 1024,1024,5   Size of subimages for correlation (default 512,512,1, use with -align mic).",
"-resolution 900,300      High and low resolution limits for bandpass (default 0.1,1000 angstrom).",
"-maxshift 55             Maximum allowed shift (default 1/4 of tile).",
"-threshold -3.0          Threshold used to pick features (in std dev units; default -3, use with -align fea).",
"-features 1,55           Features type (1=mass, 2=area, 3=darkest pixel) and maximum number (default 1,unlimited).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-image file.pif          Output multi-image file with aligned images.",
" ",
"Examples:",
"The set of input files reference a series of micrographs arranged in ",
"fields-of-view where each field contains a focal pair or set. The program then",
"calculates for each set the alignment of all other micrographs with respect",
"to the reference micrograph (typically the further-from-focus micrograph),",
"and outputs a new STAR format file with the alignment parameters and particle",
"coordinates corresponding to the alignment. The alignment can be done in",
"three different ways:",
"1. Alignment by cross-correlation (recommended):",
"		bmgalign -v 7 -reference 2 -align mic -correlate 300,300,1 -out out.star input.star",
"2.	Alignment by feature extraction:",
"		bmgalign -v 7 -ref 2 -resol 300,1000 -align fea -thresh -2.8 -out out.star input.star",
"3.	Alignment using picked particle coordinates:",
"		bmgalign -v 7 -ref 2 -align part -output out.star input.star",
" ",
"All input is derived from a STAR format parameter file.",
"A series (such as a focal pair) is assumed to be consecutive data blocks in the STAR file.",
"Orientational parameters are calculated so that applying them to the first",
"data set or image produces the others in the focal series.",
"Angles are given with positive counter-clockwise.",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize optional variables
	DataType 		datatype(Unknown_Type);		// Conversion to new type
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	int 			nseries(0);					// Number of micrographs in a series from a field
	Bstring			field_id;					// User-specified field ID
	int 			ref_img(0); 				// Reference image in series (first is 1)
	int 			add_coordinates(0); 		// For adding mg coordinates and particle origins
	int 			set_align(0);				// For alignment options
	int				filter_flag(0);
//	int				set_center(0);				// Centering images before output
	Vector3<long>	tile_size(512,512,1);		// Size of subimage to correlate
	int 			max_features(1000000000);	// Maximum number of features to extract
	double	 		resolution_low(1000);		// Low resolution limit
	double	 		resolution_high(0.1);		// High resolution limit
	double			max_shift(0);
	double			thresh(-3.0);				// Threshold to find features
	int				extract_method(0); 			// Method for defining a feature
	Bstring			outfile;					// Output parameter file name
	Bstring			imgfile;					// Output multi-image file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "align" ) {
			if ( curropt->value.contains("par") ) set_align = 1;
			else if ( curropt->value.contains("mic") ) set_align = 3;
			else if ( curropt->value.contains("fea") ) set_align = 4;
			else
				cerr << "-align: An alignment method must be specified." << endl;
		}
		if ( curropt->tag == "filterextremes" )
			filter_flag = 1;
		if ( curropt->tag == "number" )
			if ( ( nseries = curropt->value.integer() ) < 1 ) {
				cerr << "-number: The number of micrographs in a field-of-view must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "datatype" )
			datatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "fieldID" ) {
			field_id = curropt->value;
			if ( field_id.length() < 1 )
				cerr << "-fieldID: A field ID must be specified!" << endl;
		}
		if ( curropt->tag == "reference" )
			if ( ( ref_img = curropt->value.integer() ) < 1 )
				cerr << "-reference: The reference micrograph number in a series must be specified." << endl;
		if ( curropt->tag == "add" )
        		add_coordinates = 1;
		if ( curropt->tag == "correlate" ) {
			tile_size = curropt->size();
			if ( tile_size.volume() < 1 )
				cerr << "-correlate: All three dimensions must be specified." << endl;
		}
		if ( curropt->tag == "features" )
			if ( curropt->values(extract_method, max_features) < 1 )
				cerr << "-features: A method for extraction must be specified." << endl;
		if ( curropt->tag == "threshold" )
			if ( ( thresh = curropt->value.real() ) < -3 )
				cerr << "-threshold: A threshold must be specified." << endl;
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(resolution_high, resolution_low) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( resolution_high > resolution_low )
				swap(resolution_high, resolution_low);
        }
		if ( curropt->tag == "maxshift" )
    	    if ( ( max_shift = curropt->value.real() ) < 0.1 )
				cerr << "-maxshift: A maximum shift in pixels must be specified." << endl;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);

	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	if ( nseries )
		project_set_field_id(project, nseries, field_id);

	if ( add_coordinates )
		project_add_origins_to_coords(project);

	switch ( set_align ) {
		case 1:
			mg_align_coordinates(project, ref_img);
			break;
		case 2:
		case 3:
			mg_align_micrographs(project, ref_img, tile_size, 
					resolution_low, resolution_high, max_shift, filter_flag, set_align - 2);
			break;
		case 4:
			mg_align_feature_extraction(project, max_features,
					resolution_low, resolution_high, thresh, extract_method);
	}
	
	if ( imgfile.length() )
		project_write_aligned_images(project, NULL, imgfile, datatype);
	
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
