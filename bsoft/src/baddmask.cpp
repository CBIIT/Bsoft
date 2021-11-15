/**
@file	baddmask.cpp
@brief	Program to catenate image files
@author Bernard Heymann
@date	Created: 20091117
@date	Modified: 20170117
**/

#include "rwimg.h"
#include "img_combine.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bimage*		img_composite_mask(Bstring* file_list, int mask_type);

// Usage assistance
const char* use[] = {
" ",
"Usage: baddmask [options] input.images",
"--------------------------------------",
"Adds multiple masks to generate a new composite mask.",
"Any number of input images may be given and may include the wild card \"*\".",
" ",
"Actions:",
"-type bit                Type of composite mask: uniform/bit/level.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
" ",
"Output:",
"-output output.img       Output image (default sum.map).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				mask_type(0);				// Default: simply add images
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;			// New image origin
	int				set_origin(0);				// Flag to set origin
	Bstring			maskfile("mask.tiff");		// Default output image file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'u' ) mask_type = 1;
			if ( curropt->value[0] == 'b' ) mask_type = 2;
			if ( curropt->value[0] == 'l' ) mask_type = 3;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "output" )
        	maskfile = curropt->filename();
    }
	option_kill(option);
	
	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();

	if ( optind >= argc ) {
		cerr << "No input files!" << endl;
		bexit(-1);
	}
	
	// Set up image file names and weights
	int				nfiles(0);
	Bstring*		file_list = NULL;
	while ( optind < argc ) {
		string_add(&file_list, argv[optind++]);
		nfiles++;
	}

	Bimage*			pmask = NULL;
	vector<double>	weight;
	
	switch ( mask_type ) {
		case 1:
			pmask = img_add_weighed(file_list, weight);
			pmask->truncate_to_min_max(0, 1);
			break;
		case 2:
		case 3:
			pmask = img_composite_mask(file_list, mask_type);
			break;
		default:
			pmask = img_add_weighed(file_list, weight);
	}
	
	if ( sam.volume() > 0 ) pmask->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) pmask->origin(pmask->default_origin());
		else pmask->origin(origin);
	}
	pmask->mask_stats();

	if ( maskfile.length() ) {
		pmask->change_type(nudatatype);
		write_img(maskfile, pmask, 0);
	}
	
	delete pmask;
	string_kill(file_list);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Generates a composite of multiple masks.
@param 	*file_list 	list of file names.
@param 	mask_type	type of mask to generate (2/3).
@return Bimage* 	composite mask (4-byte integer).

	The composite mask can be either a multi-bit mask (mask_type = 2),
	or a multi-level mask (mask_type = 3).

**/
Bimage*		img_composite_mask(Bstring* file_list, int mask_type)
{
	long		nimg(0);
	Bimage*	 	pmask = img_setup_combined(file_list, nimg);
	
	if ( !pmask ) return NULL;
	
	pmask->change_type(Integer);
		
	Bstring*		filename;
	Bimage*			p = NULL;
	int				val, bit, base_bit;
	long			i, j, n, nf, nfiles(0), base_level, max_level;

	for ( nfiles=0, filename = file_list; filename; filename = filename->next ) nfiles++;
	
	if ( verbose ) {
		if ( mask_type == 2 )
			cout << "Generating a multi-bit mask from " << nfiles << " images:" << endl;
		else
			cout << "Generating a multi-level mask from " << nfiles << " images:" << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "New mask size:                  " << pmask->size() << endl;
	
	long			imgsize(pmask->image_size());
	
	for ( nf=0, base_bit=0, base_level=max_level=0, filename = file_list; 
			filename && nf<nfiles; filename = filename->next, nf++ ) {
		p = read_img(*filename, 1, -1);
		if ( p != NULL ) {
			if ( verbose & VERB_LABEL )
				cout << "Adding mask " << nf << endl;
			p->change_type(Float) ;
			for ( n=j=0; n<p->images(); n++ ) {
				pmask->origin(n, p->image[n].origin());
				for ( i=0; i<imgsize; i++, j++ ) {
					if ( (*p)[j] > 0.5 ) {		// What is zero stay whatever it currently is
						val = (int) ((*p)[j] + 0.5);
						if ( mask_type == 2 ) {
							bit = val;
							if ( base_bit > 0 ) bit <<= base_bit;
							pmask->set(j, (char)(*pmask)[j] | bit);
						} else {
							pmask->set(j, base_level + val);
						}
						if ( max_level < (*pmask)[j] ) max_level = (int) (*pmask)[j];
					}
				}
			}
			delete p;
			for ( base_bit = 0, val = max_level; val > 0; val >>= 1, base_bit++ ) ; 
			base_level = max_level;
		}
	}
	
	pmask->statistics();
	
	return pmask; 
}

