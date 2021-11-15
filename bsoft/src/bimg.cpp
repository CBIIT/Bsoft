/**
@file	bimg.cpp
@brief	General image processing program
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210624
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
"Usage: bimg [options] input.img output.img",
"------------------------------------------",
"Manipulates images in a few general ways.",
"Image format extensions supported:",
"BioRad (pic), Brix (brx), CCP4 (map), Digital Instruments (di)",
"EM (em), Goodford (pot), GRD (grd), Imagic (img), JPEG, (jpg), MFF (mff)", 
"Image Magick (miff), MRC (mrc), PIC (bp), PIF (pif), Spider (spi)", 
"Suprim (spm), TIFF (tif), RAW (raw)",
" ",
"Actions:",
"-delete 1,5-13,122       Delete the indicated sub-images (first image = 0).",
"-slices                  Interpret multiple 2D images as z-slices of a single 3D image",
"-images                  Interpret slices of a single 3D image as 2D images",
"-invert                  Invert density in the image.",
"-project                 Project along z axis (after reslicing).",
"-sum                     Sum images in a multi-image file.",
"-movingsum 4,3,1         Calculate moving sum of multi-image with step size and averaging flag.",
"-average                 Average images in a multi-image file.",
"-reslice -z+xy           Reslice = switch axes (default xyz).",
"-size 10,50,8            New image size (pixels/voxels).",
"-translate -5,12,50      Translate (pixels, with wrapping if option -wrap is used).",
"-Wrap 10,50,8            Shrinking and wrapping an image (pixels).",
"-truncate -0.5,1.2       Truncate data to minimum and maximum.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation (after -truncate).",
"-minmax -0.5,1.2         Rescale data to minimum and maximum (after -truncate).",
"-range -2.3,5.4          Fill with fill value outside this range.",
"-logarithm               Calculate the logarithm of an image.",
"-exponential             Calculate the exponential of an image.",
"-Levels 8                Restrict to a number of levels.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-select 1                Select an image (first image = 0; default: all images).",
"-fill 127                Fill value for resizing: average (default), background, or value.",
"-setbackground 0.5       Set background.",
"-wrap                    Turn wrapping on (default off).",
" ",
"Output:",
"-std stdev.mrc           Standard deviation map.",
"-compression 2           Compression type: 5=LZW (TIFF only).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	int 			setinvert(0);				// Flag to invert density
	int 			setreslice(0);
	Bstring			order("xyz");				// Reslice order
	int 			setproject(0);
	int 			setsum(0);
	int 			setaverage(0);
	int 			setrescale(0);
	int 			setlog(0);
	int 			setexp(0);
	long			mov_win(0), mov_step(1), mov_flag(0);	// Moving sum window and step
	double			nuavg(0), nustd(0);			// Rescaling to average and stdev
	double			numin(0), numax(0);			// Rescaling to new min and max
	double			rangemin(0), rangemax(0);	// Range to keep
	int 			settruncate(0);
	double			cutmin(0), cutmax(0); 		// Truncation
	int 			setshrinkwrap(0); 			// Shrinking and wrapping flag
	int 			settranslate(0);			// Translation flag
	int 			setwrap(0);					// Wrapping flag
	Vector3<int> 	nusize;						// New image size
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<int> 	shift;						// Shift for integral translation
	double			fill(0);	 				// Fill value for resizing
	int 			fill_type(FILL_AVERAGE);	// Fill type for resizing
	int 			nlevels(0);
	int 			setimg(-1);					// Select all images
	Bstring			select_list;				// List of sub-images to select
	Bstring			delete_list;				// List of sub-images to delete
	int 			setbackground(0);			// Background
	double	 		nubackground(-1);			// New background
	int				znswitch(0);				// 0=not, 1=n2z, 2=z2n
	Bstring			stdfile;					// Standard deviation output map
//	int				write_flags(0);				// Flags for image properties during saving
	int				compression(0);				// Output compression type
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "select" ) {
			if ( curropt->value.contains(",") ||
					curropt->value.contains("-") )
				select_list = curropt->value;
			else
				setimg = curropt->value.integer();
		}
		if ( curropt->tag == "delete" ) {
			delete_list = curropt->value;
			if ( delete_list.length() < 1 )
				cerr << "-delete: Image numbers must be specified!" << endl;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "project" )
			setproject = 1;
		if ( curropt->tag == "average" )
			setaverage = 1;
		if ( curropt->tag == "sum" )
			setsum = 1;
		if ( curropt->tag == "wrap" )
			setwrap = 1;
		if ( curropt->tag == "Wrap" ) {
			nusize = curropt->size();
			setshrinkwrap = 1;
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "reslice" ) {
			if ( curropt->value.length() < 3 )
				cerr << "-reslice: The order of x, y and z must be specified!" << endl;
			else {
				order = curropt->value;
				setreslice = 1;
			}
		}
		if ( curropt->tag == "movingsum" )
			if ( curropt->values(mov_win, mov_step, mov_flag) < 1 )
				cerr << "-movingsum: A window size must be specified!" << endl;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else
				setrescale = 1;
			if ( nustd <= 0 ) setrescale = 0;
		}
		if ( curropt->tag == "minmax" ) {
			if ( curropt->values(numin, numax) < 2 )
				cerr << "-minmax: Both min and max must be specified!" << endl;
			else
				setrescale = 2;
			if ( numin == numax ) setrescale = 0;
		}
		if ( curropt->tag == "range" )
			if ( curropt->values(rangemin, rangemax) < 2 )
				cerr << "-range: Both range min and max must be specified!" << endl;
		if ( curropt->tag == "logarithm" )
			setlog = 1;
		if ( curropt->tag == "exponential" )
			setexp = 1;
		if ( curropt->tag == "truncate" ) {
			if ( curropt->values(cutmin, cutmax) < 2 )
				cerr << "-truncate: Both min and max must be specified!" << endl;
			else
				settruncate = 1;
			if ( cutmin > cutmax ) swap(cutmin, cutmax);
			if ( cutmin == cutmax ) settruncate = 0;
		}
		if ( curropt->tag == "size" )
			nusize = curropt->size();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "translate" ) {
			shift = curropt->vector3();
			if ( shift.length() ) settranslate = 1;
		}
		if ( curropt->tag == "Levels" ) nlevels = curropt->value.integer();
		if ( curropt->tag == "setbackground" ) {
			nubackground = curropt->value.real();
			setbackground = 1;
        }
		if ( curropt->tag == "slices" )
			znswitch = 1;
		if ( curropt->tag == "images" )
			znswitch = 2;
		if ( curropt->tag == "std" )
			stdfile = curropt->filename();
		if ( curropt->tag == "compression" )
			if ( ( compression = curropt->integer() ) < 1 )
				cerr << "-compression: A number must be specified!" << endl;
    }
	option_kill(option);
	
	double 		ti = timer_start();
	
	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, setimg);
	if ( p == NULL ) bexit(-1);
	
	if ( znswitch == 1 )
		if ( p->images_to_slices() < 0 ) bexit(-1);
	
	if ( znswitch == 2 ) 
		if ( p->slices_to_images() < 0 ) bexit(-1);

	if ( select_list.length() )
		p->select_images(select_list);
	
	if ( delete_list.length() )
		p->delete_images(delete_list);

	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	else if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( setbackground ) p->background(nubackground);
	
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) {
		if ( fabs(p->background(long(0))) < 1e-6 )
			p->calculate_background();
		fill = p->background(long(0));
	}
	
	if ( nusize.volume() > 0 ) {
		if ( setshrinkwrap )
			p->shrink_wrap(nusize, shift);
		else {
			if ( verbose & VERB_PROCESS ) {
				if ( setwrap ) cout << "Resizing with wrapping:" << endl;
				else cout << "Resizing:" << endl;
				cout << "Shift:                          " << shift << endl;
				cout << "New size:                       " << nusize << endl;
				if ( fill_type != FILL_BACKGROUND )
					cout << "Fill value:                     " << fill << endl;
				cout << endl;
			} else if ( verbose & VERB_LABEL )
				cout << "Resizing" << endl << endl;

			if ( setwrap ) p->resize_wrap(nusize, shift);
			else p->resize(nusize, shift, fill_type, fill);
		}
	} else if ( settranslate ) {
		if ( setwrap ) p->shift_wrap(shift);
		else p->shift(shift, fill_type, fill);
	}
	
	if ( setinvert ) p->invert();
	
	if ( settruncate ) {
		if ( fill_type == FILL_USER || fill_type == FILL_BACKGROUND )
			p->truncate(cutmin, cutmax, fill, fill);
		else p->truncate_to_min_max(cutmin, cutmax);
	}
	
	if ( setrescale == 1 ) p->rescale_to_avg_std(nuavg, nustd);
	
	if ( setrescale == 2 ) p->rescale_to_min_max(numin, numax);
	
	if ( rangemin < rangemax ) {
		if ( fill_type == FILL_AVERAGE ) p->truncate_to_avg(rangemin, rangemax);
		else if ( fill_type == FILL_BACKGROUND ) p->truncate_to_background(rangemin, rangemax);
		else p->truncate(rangemin, rangemax, fill, fill);
	}
	
	if ( setlog ) p->logarithm();
	
	if ( setexp ) p->exponential();
	
	if ( setreslice ) p->reslice(order);
	
	if ( nlevels ) p->limit_levels(nlevels);

	Bimage*		p2 = NULL;
	if ( setproject ) {
		p2 = p->project('z', 1);
		delete p;
		p = p2;
	}
	
	if ( mov_win > 1 ) {
		p2 = p->moving_sum(mov_win, mov_step, mov_flag);
		delete p;
		p = p2;
	}
	
	if ( setsum )
		p->sum_images();
	
	if ( setaverage ) {
		if ( stdfile.length() ) {
			p2 = p->average_images(1);
			delete p;
			p = p2;
		} else {
			p->average_images();
		}
	}
	
	if ( sam.volume() > 0 ) p->sampling(sam);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}

    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, compression);
	}

    // Write a standard deviation output file if a file name is given
    p2 = p->next;
	if ( stdfile.length() && p2 )
		write_img(stdfile, p2, 0);

	delete p;
	order = 0;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
