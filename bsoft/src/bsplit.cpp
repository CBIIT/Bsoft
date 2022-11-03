/**
@file	bsplit.cpp
@brief	Program to seperate a multi-image file into individual image files
@author Bernard Heymann
@date	Created: 20001026
@date	Modified: 20220203
**/

#include "rwimg.h"
#include "Bstring.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bsplit [options] input.img output.img",
"--------------------------------------------",
"Splits multi-image files of the following formats:",
"Digital Instruments (di), Imagic (img), Image Magick (miff), PIF (pif)", 
"SPIDER (spi), TIFF (tif), RAW (raw).",
"Files of single image 3D image formats can be split with the -images option.",
" ",
"Actions:",
"-images                  Interpret slices of a single 3D image as 2D images.",
"-channels                Splits channels into a series of images.",
"-invert                  Invert density in the image.",
"-center                  Center images.",
"-reslice -z+xy           Reslice = switch axes (default xyz).",
"-select 1,5-13,122       Sub-images to select (first image = 0, default all).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-first 5                 Number given to the first file (default 0).",
"-digits 3                Number of digits inserted before the last period in the output file name.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setZslices(0);             	// Interpret 3D slices as separate 2D images
	int				setchan(0);					// Split channels into images
	int 			setinvert(0);				// Flag to invert density
	int				setcenter(0);				// Flag to center images with wrapping
	int 			setreslice(0);
	Bstring			order("xyz");
	int 			first_number(0);			// Number given to first file
	int 			digits(3);					// File number size
	Bstring			select_list;				// List of sub-images to select, default all
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "images" )
			setZslices = 1;
		if ( curropt->tag == "channels" )
			setchan = 1;
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "center" )
			setcenter = 1;
		if ( curropt->tag == "reslice" ) {
			if ( curropt->value.length() < 3 )
				cerr << "-reslice: The order of x, y and z must be specified!" << endl;
			else {
				order = curropt->value;
				setreslice = 1;
			}
		}
		if ( curropt->tag == "select" ) {
			select_list = curropt->value;
			if ( select_list.length() < 1 )
				cerr << "-select: Image numbers must be specified!" << endl;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "first" )
			if ( ( first_number = curropt->value.integer() ) < 0 )
				cerr << "-first: A number must be specified!" << endl;
		if ( curropt->tag == "digits" )
			if ( ( digits = curropt->value.integer() ) < 1 )
				cerr << "-digits: A number of digits must be specified!" << endl;
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read image file
	Bstring		inputfile = argv[optind++];
	Bimage*		p = read_img(inputfile, 1, -1);
	if ( p == NULL ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
    
	if ( optind >= argc ) bexit(0);
	
	Bstring		filename;
	Bstring		outputfile = argv[optind];
	
	char			format[32];
	if ( digits < 0.44*log(1.0*p->images()) + 1 ) digits = (int) (0.44*log(1.0*p->images()) + 1);
	snprintf(format, 32, "_%%0%dd.", digits);
	
    if ( setZslices ) 
		if ( p->slices_to_images() < 0 ) bexit(-1);

	if ( setinvert ) p->invert();
	
	if ( setcenter ) p->center_wrap();

	Bimage*			pone = NULL;

	if ( setchan ) {
		if ( verbose & VERB_LABEL )
			cout << "Splitting " << inputfile << " into " << p->channels() << " files" << endl;
	
		Bimage*		pchan = p->split_channels_to_images();
		pone = pchan;
		
		for ( i=0; pone; pone = pone->next, ++i ) {
			pone->change_type(nudatatype);
			if ( setreslice ) pone->reslice(order);
			filename = outputfile.pre_rev('.') + Bstring(i+first_number, format) + outputfile.post_rev('.');
			if ( verbose )
				cout << "Image " << i << ":  " << filename << endl;
			write_img(filename, pone, 0);
		}
		
		delete pchan;
		
	} else {
		if ( verbose & VERB_LABEL )
			cout << "Splitting " << inputfile << " into " << p->images() << " files" << endl;

		vector<int>		imgnum = select_numbers(select_list, p->images());

		for ( i=0; i<p->images(); i++ ) if ( imgnum[i] ) {
			pone = p->extract(i);
			pone->image[0] = p->image[i];
			pone->label(p->label());
			pone->change_type(nudatatype);
			if ( setreslice ) pone->reslice(order);
			filename = outputfile.pre_rev('.') + Bstring(i+first_number, format) + outputfile.post_rev('.');
			if ( verbose )
				cout << "Image " << i << ":  " << filename << endl;
			write_img(filename, pone, 0);
			delete pone;
		}
	
		if ( verbose & VERB_PROCESS )
			cout << endl;
	}
	
	delete p;
	order = 0;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
