/**
@file	bcolour.cpp
@brief	Program to generate nice colour images
@author Bernard Heymann
@date	Created: 20001004
@date	Modified: 20180412
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
"Usage: bcolour [options] input.img output.img",
"---------------------------------------------",
"Generates nice colour images.",
"Only image file formats capable of handling RGB can be used for output",
"	(TIFF, PNG, JPEG).",
" ",
"Actions:",
"-invert                  Invert image.",
"-red -4.5,88.3           Convert a grayscale range to red.",
"-green -4.5,88.3         Convert a grayscale range to green.",
"-blue -4.5,88.3          Convert a grayscale range to blue.",
"-spectrum 1,22.4         Convert a grayscale range to a spectral colour scale.",
"-rwb 0,56,167,245        Convert a grayscale range to a red, white and blue image.",
"-phases                  Colour the phases in a polar map or complex transform.",
"-intensities             Calculates the intensities without colours.",
"-pure                    Calculates pure colours without intensity.",
"-alpha                   Add an alpha (transparency) channel.",
"-noalpha                 Remove the alpha (transparency) channel.",
"-torgb                   Convert grayscale or CMYK to RGB.",
"-tocmyk                  Convert grayscale or RGB to CMYK.",
"-split                   Split the channels into individual images.",
"-combine 3               Combine channels from individual images.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-subtractive             Flag to calculate subtractive one-color conversions.",
"-grayextremes            Flag to keep gray extremes beyond colorization limits.",
" ",
"Input:",
"-add file.tif            Add this file to the input file (must be the same size).",
" ",
"Output:",
"-scale file.tif          Output spectral scale image.",
"-phasewheel pcw.png      Output a phase color wheel (size 512,512,1)."
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setimg(-1);					// Selects all images
	int 			setinvert(0);
	int				flag(0);					// Flag: 1=subtractive, 2=gray scale extremes
	int				setred(0);
	int				setgreen(0);
	int				setblue(0);
	int 			setspectrum(0);
	int 			setrwb(0);
	int				setalpha(0);				// No alpha channel change
	int				to_cmyk(0);
	int				to_rgb(0);
	int 			colour_phases(0);
	int				calcintensity(0);
	int				calcpure(0);
	int				split(0);
	int				combine(0);
	double			cmin(0), cmax(0);			// Min/max for conversion to spectrum
	double			red_min(0), white_min(0), white_max(0), blue_max(0);
	Bstring			add_image;					// Image to add to input file
	Bstring			scale_image;				// Spectral scale image
	Bstring			pcw;						// Phase color wheel image
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "red" ) {
			if ( curropt->values(cmin, cmax) < 2 )
				cerr << "-red: A grayscale range must be specified!" << endl;
			else
				setred = 1;
		}
		if ( curropt->tag == "green" ) {
			if ( curropt->values(cmin, cmax) < 2 )
				cerr << "-green: A grayscale range must be specified!" << endl;
			else
				setgreen = 1;
		}
		if ( curropt->tag == "blue" ) {
			if ( curropt->values(cmin, cmax) < 2 )
				cerr << "-blue: A grayscale range must be specified!" << endl;
			else
				setblue = 1;
		}
		if ( curropt->tag == "spectrum" ) {
			if ( curropt->values(cmin, cmax) < 2 )
				cerr << "-spectrum: A grayscale range must be specified!" << endl;
			else
				setspectrum = 1;
		}
		if ( curropt->tag == "rwb" ) {
        	 if ( curropt->values(red_min, white_min,
						white_max, blue_max) < 4 )
				cerr << "-rwb: 4 grayscale values must be specified!" << endl;
			else
				setrwb = 1;
		}
		if ( curropt->tag == "phases" )
        	colour_phases = 1;
		if ( curropt->tag == "intensities" )
        	calcintensity = 1;
		if ( curropt->tag == "pure" )
        	calcpure = 1;
		if ( curropt->tag == "alpha" )
        	setalpha = 1;
		if ( curropt->tag == "noalpha" )
        	setalpha = -1;
		if ( curropt->tag == "torgb" )
        	to_rgb = 1;
		if ( curropt->tag == "tocmyk" )
        	to_cmyk = 1;
		if ( curropt->tag == "split" )
        	split = 1;
		if ( curropt->tag == "combine" )
 			if ( ( combine = curropt->value.integer() ) < 1 )
				cerr << "-combine: The number of channels must be specified!" << endl;
		if ( curropt->tag == "subtractive" )
        	flag |= 1;
		if ( curropt->tag == "grayextremes" )
        	flag |= 2;
		if ( curropt->tag == "add" )
			add_image = curropt->filename();
		if ( curropt->tag == "scale" )
			scale_image = curropt->filename();
		if ( curropt->tag == "phasewheel" )
			pcw = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	if ( pcw.length() ) {
		Bimage*		pnu = new Bimage(UCharacter, TRGB, 512, 512, 1, 1);
		pnu->phase_colour_wheel();
		write_img(pcw, pnu, 0);
		delete pnu;
	}
	
	if ( optind >= argc ) bexit(0);
	
	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, setimg);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	Bimage*		p2 = NULL;
	Bimage*		pscale = NULL;
	
	if ( optind < argc ) {
		if ( setinvert ) p->invert();
		
		if ( to_cmyk ) {
			p->simple_to_rgb();
			p->rgb_to_cmyk();
		}
		
		if ( to_rgb ) {
			p->simple_to_rgb();
			p->cmyk_to_rgb();
		}
	
		if ( setred ) p->color_red(cmin, cmax, flag);

		if ( setgreen ) p->color_green(cmin, cmax, flag);

		if ( setblue ) p->color_blue(cmin, cmax, flag);

		if ( setspectrum ) pscale = p->color_spectrum(cmin, cmax);
	
		if ( setrwb )
			pscale = p->red_white_blue(red_min, white_min, white_max, blue_max);
	
		if ( colour_phases ) {
			p2 = p->intensities_phase_colored(0);
			if ( p2 ) {
				delete p;
				p = p2;
			}
		}
		
		if ( calcintensity ) p->color_to_simple();
		
		if ( calcpure ) p->pure_color();
		
		if ( setalpha > 0 ) p->rgb_to_rgba();
		else if ( setalpha < 0 ) p->rgba_to_rgb();
		
		if ( add_image.length() ) {
			p2 = read_img(add_image, 1, -1);
			p->color_combine(p2);
			delete p2;
		}
		
		if ( split ) {
			Bimage*		pnu = p->split_channels();
			delete p;
			p = pnu;
		}
    	
		if ( combine ) p->combine_channels(combine);
    	
		p->change_type(nudatatype);
	
			
		write_img(argv[optind], p, 0);
	}
	
	if ( scale_image.length() && pscale ) {
		write_img(scale_image, pscale, 0);
		delete pscale;
	}
	
	delete p;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
