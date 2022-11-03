/**
@file	bcomplex.cpp
@brief	Program for handing complex images.
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20220620
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
"Usage: bcomplex [options] input.img output.img",
"----------------------------------------------",
"Converts complex image forms.",
" ",
"Actions:",
"-size 120,102,200        Resize the image as a Fourier transform.",
"-tocomplex               Convert a phase image to a complex image (phases in radians).",
"-from2                   Convert two sub-images to a complex image.",
"-invert                  Inverts a complex image.",
"-logarithm               Calculate the logarithm of the image.",
"-amplitudes              Convert a complex image to amplitudes.",
"-intensities             Convert a complex image to intensities.",
"-real                    Convert a complex image to real values.",
"-imaginary               Convert a complex image to imaginary values.",
"-phases                  Convert a complex image to phases.",
"-signed                  Convert a complex image to amplitudes, with sign based on phase.",
"-positive                Convert a complex image to be postive definite.",
"-friedel                 Check and apply Friedel symmetry.",
"-ewald opp               Convert an Ewald projection to the opposite or combined.",
"-shift -3,4.5,-0.3       Shift phases (pixels, can be half).",
"-correlate other.sup     Calculate correlation.",
"-color 1.5               Generate a phase-coloured power spectrum: amplitude scaling.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<long>	nusize;				// New size if > 0
	Vector3<double>	sam;    			// Units for the three axes (A/pixel)
	bool			tocomplex(0);		// Flag to convert a phase image to complex
	bool			from2(0);			// Flag to convert a phase image to complex
	bool			invert(0);			// Flag to invert
	bool 			setlogarithm(0);	// Flag for logarithmic power spectrum
	int 			conv_comp(0);		// Flag for complex to real conversions
	bool			positive(0);		// Flag to convert to positive definite
	bool			friedel(0);			// Flag for applying Fiedel symmetry
	int				ewald(0);			// Flag for converting an Ewald projection
	bool			half(0);			// Flag to shift half of the image size
	Vector3<double>	shift;				// Shift vector
	double	 		colour_phase_scale(0);		// Scale for coloured phases
	Bstring			file_for_corr;		// Complex image to calculate difference
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "size" )
			nusize = curropt->size();
		if ( curropt->tag == "tocomplex" )
        	tocomplex = 1;
		if ( curropt->tag == "from2" )
        	from2 = 1;
		if ( curropt->tag == "invert" )
        	invert = 1;
		if ( curropt->tag == "logarithm" )
			setlogarithm = 1;
		if ( curropt->tag == "amplitudes" )
			conv_comp = 1;
		if ( curropt->tag == "intensities" )
        	conv_comp = 2;
		if ( curropt->tag == "real" )
        	conv_comp = 3;
		if ( curropt->tag == "imaginary" )
        	conv_comp = 4;
		if ( curropt->tag == "phases" )
        	conv_comp = 5;
		if ( curropt->tag == "signed" )
        	conv_comp = 6;
		if ( curropt->tag == "positive" )
        	positive = 1;
		if ( curropt->tag == "friedel" )
        	friedel = 1;
		if ( curropt->tag == "ewald" ) {
			if ( curropt->value[0] == 'c' ) ewald = 2;
			else ewald = 1;
		}
		if ( curropt->tag == "shift" ) {
			if ( curropt->value[0] == 'h' ) half = 1;
			else shift = curropt->vector3();
		}
 		if ( curropt->tag == "color" )
			if ( ( colour_phase_scale = curropt->value.real() ) < 0.001 )
		    		cerr << "-color: A scale for the phase colours must be specified" << endl;
		if ( curropt->tag == "correlate" )
			file_for_corr = curropt->filename();
    }
	option_kill(option);
	
 	double		ti = timer_start();
	
	int 		dataflag(friedel);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( !dataflag ) bexit(0);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	else if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	if ( tocomplex ) p->phase_to_complex();
	
	if ( from2 ) p->two_to_complex();
	
	if ( nusize.volume() > 0 ) p->change_transform_size(nusize);
	
	if ( invert ) p->complex_invert();
	
	if ( half ) shift = p->size()/2;
	if ( shift.length() ) {
		p->phase_shift(shift);
	}
	
	if ( friedel ) {
		p->friedel_check();
		p->friedel_difference();
//		p->friedel_apply();
	}
	
	if ( ewald == 1 ) p->opposite_ewald();
	if ( ewald == 2 ) p->combine_ewald();
	
	if ( file_for_corr.length() ) {
		Bimage*		pref = read_img(file_for_corr, 1, -1);
		p->complex_conjugate_product(pref, 1);
		delete pref;
	}
	
	if ( positive )
		p->fspace_positive();
	
	if ( conv_comp && p->compound_type() == TComplex ) {
		switch ( conv_comp ) {
			case 1: p->complex_to_amplitudes(); break;
			case 2: p->complex_to_intensities(); break;
			case 3: p->complex_to_real(); break;
			case 4: p->complex_to_imaginary(); break;
			case 5: p->complex_to_phases(); break;
			case 6: p->complex_to_signed_amplitudes(); break;
		}
		p->statistics();
	} else if ( colour_phase_scale > 0 ) {
		Bimage*		pnu = p->intensities_phase_colored(colour_phase_scale);
		if ( pnu ) {
			delete p;
			p = pnu;
			nudatatype = p->data_type();
		}
		Bstring		filename("ps_color_wheel.png");
		pnu = new Bimage(UCharacter, TRGB, 512, 512, 1, 1);
		pnu->phase_colour_wheel();
		write_img(filename, pnu, 0);
		delete pnu;
	}
	
	if ( setlogarithm ) p->logarithm();
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}

	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
