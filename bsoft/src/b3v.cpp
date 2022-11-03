/**
@file	b3v.cpp
@brief	Program to extract 3 orthogonal views from a map
@author Bernard Heymann
@date	Created: 20100428
@date	Modified: 20150718
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
"Usage: b3v [options] input.img output.img",
"-----------------------------------------",
"Generates 3 orthogonal views of a map.",
" ",
"Actions:",
"-orthogonal 24,67,45     Extract orthogonal slices around this voxel.",
"                         The voxel can be specified as center or origin.",
"-montage 2               Montage the slices into a number of columns.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-select 5                Image number.",
"-size 100,80,120         Size of orthogonal slices to extract.",
"-pad 5                   Padding to add around slices or images (default 0).",
"-fill 127                Fill value for padding (default average).",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int				orth_type(0);				// Type of voxel specification: 1=center, 2=origin, 3=coor
	long			montage_cols(0);			// Number of montage columns
	Vector3<long>	voxel;						// Intersection of slices
	long			img_num(0);					// Image number
	Vector3<long>	size;						// Size of slices to extract
	int 			pad(0);						// Padding around images
	double			fill(0);	 				// Fill value for resizing
	int 			fill_type(FILL_BACKGROUND);	// Fill type for resizing
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "montage" )
			if ( ( montage_cols = curropt->value.integer() ) < 1 )
				cerr << "-montage: A number of columns must be specified!" << endl;
		if ( curropt->tag == "orthogonal" ) {
			if ( curropt->value[0] == 'c' ) orth_type = 1;
			else if ( curropt->value[0] == 'o' ) orth_type = 2;
			else {
				if ( curropt->values(voxel[0], voxel[1], voxel[2]) < 3 )
					cerr << "-orthogonal: A voxel must be specified!" << endl;
				else orth_type = 3;
			}
		}
		if ( curropt->tag == "select" )
			if ( ( img_num = curropt->value.integer() ) < 0 )
				cerr << "-select: Image numbers start from 0!" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "pad" )
			if ( ( pad = curropt->value.integer() ) < 1 )
				cerr << "-pad: A number of pixels must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
	}
	option_kill(option);
	
	double		ti = timer_start();
	
	Bimage*		p = read_img(argv[optind++], 1, -1);
	if ( p == NULL ) bexit(-1);
    
	if ( optind >= argc ) {
		delete p;
		bexit(0);
	}
	
	Bimage*			p3 = NULL;
	
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) {
		if ( fabs(p->background(long(0))) < 1e-6 )
			p->calculate_background();
		fill = p->background(long(0));
	}
	if ( fill_type == FILL_MIN ) fill = p->minimum();
	if ( fill_type == FILL_MAX ) fill = p->maximum();
	
	p->background(fill);
	fill_type = FILL_BACKGROUND;
	
	if ( orth_type ) {
		if ( orth_type < 3 ) {
			Vector3<double>	vori;
			if ( orth_type == 1 )
				vori = p->default_origin();
			else if ( orth_type == 2 )
				vori = p->image[img_num].origin();
			voxel = Vector3<long>(vori[0], vori[1], vori[2]);
		}
		if ( montage_cols ) {
			p3 = p->orthogonal_slices(img_num, voxel, size);
		} else {
			p3 = p->orthogonal_montage(voxel, size, pad, fill_type, fill);
			pad = 0;
		}
		delete p;
		p = p3;
	
		if ( pad ) {
			Vector3<long>	newsize = {p->sizeX()+2*pad, p->sizeY()+2*pad, p->sizeZ()};
			Vector3<long>	translate = {pad, pad, 0};
			p->resize(newsize, translate, fill_type, fill);
		}

		if ( montage_cols ) {
			p3 = p->montage(0, montage_cols, 4 - montage_cols, 0, 1);
			delete p;
			p = p3;
		}
	}
	
	p->change_type(nudatatype);
	write_img(argv[optind], p, 0);

    delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

