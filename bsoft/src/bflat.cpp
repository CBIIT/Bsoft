/**
@file	bflat.cpp
@brief	Solvent flattening
@author Bernard Heymann
@date	Created: 20041113
@date	Modified: 20150717

	An implementation based on the method proposed in:
	van Heel, M. (2001). "Do single (ribosome) molecules phase themselves?" Cold Spring Harb Symp Quant Biol 66: 77-86.
**/

#include "rwimg.h"
#include "symmetry.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bflat [options] input.img output.img",
"-------------------------------------------",
"Does phase extension by solvent flattening.",
"The amplitudes are defined in three possible mutually exclusive ways:",
"1. Derived from the input map (default).",
"2. Derived from a reference map.",
"3. Specified in an external amplitude map.",
" ",
"Actions:",
"-symmetrize C5           Apply point group symmetry.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-resolution 35,15.5,1.5  Resolution limit start, end and step (default 30,20,1).",
"-iterations 80           Iterations (default 1).",
"-exit 0.1                Exit condition based on amplitude RMSD (default 0.01).",
"-output 5                Number of iterations between writing output maps (default 100).",
" ",
"Input:",
"-Reference file.map      File to generate amplitudes from (default amplitudes from input map).",
"-Amplitudes file.map     File with amplitudes (default amplitudes from input map).",
"-Mask mask.tif           Mask file.",
" ",
NULL
};

int main (int argc, char **argv)
{
	// Initialize variables
	Bsymmetry		sym;					// Point group
	DataType        nudatatype(Unknown_Type);	// Conversion to new type
	int   			i;							// loop variables
	int   			imax(1);					// largest iteration number
	int   			iout(100);					// iterations between writing output maps
	double			hi_res(20), lo_res(30);		// High and low resolution limits
	double			step_res(1);				// Resolution step size
	double			Rtarget(0.01);				// Target for amplitude RMSD
    Bstring			maskfile;					// Mask to use
    Bstring			reffile;					// Reference file
    Bstring			ampfile;					// Amplitude file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "symmetrize" )
			sym = curropt->symmetry();
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(lo_res, hi_res, step_res) < 1 )
				cerr << "-resolution: A high resolution limit must be specified!" << endl;
		if ( curropt->tag == "iterations" )
        	if ( ( imax = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "output" )
        	if ( ( iout = curropt->value.integer() ) < 1 )
				cerr << "-output: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "exit" )
        	if ( ( Rtarget = curropt->value.real() ) < 1 )
				cerr << "-exit: A target RMSD must be specified!" << endl;
 		if ( curropt->tag == "Reference" )
			reffile = curropt->filename();
 		if ( curropt->tag == "Amplitude" )
			ampfile = curropt->filename();
 		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, 0);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( optind >= argc ) {
		cerr << "Error: No output file specified!" << endl;
		bexit(-1);
	}
	

	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	
	p->change_type(Float);
	p->rescale_to_avg_std(0, 1);
	
	Bimage*		pmask = NULL;
	if ( maskfile.length() ) {
		pmask = read_img(maskfile, dataflag, 0);
		if ( pmask == NULL ) {
			cerr << "Error: Mask file " << maskfile << " not read!" << endl;
			bexit(-1);
		}
	}
	
	Bimage*		pamp = NULL;
	Bimage*		pref = NULL;
	
	if ( reffile.length() ) {
		pref = read_img(reffile, dataflag, 0);
	} else if ( ampfile.length() ) {
		pamp = read_img(ampfile, dataflag, 0);
	} else {
		pref = p->copy();
	}

	if ( pref ) pref->fft();
	
	if ( hi_res ) p->fspace_bandpass(hi_res, 1e10, 0);
		
	Bstring		filename, outname(argv[optind]);
	Bimage*		pint = NULL;
	double		res = lo_res, R = 1e37;
	
	if ( verbose )
		cout << "Iter\tRes\tRMSD" << endl;
	for ( i=0; i<imax && R>Rtarget; i++ ) {
		if ( res > hi_res ) res -= step_res;

		p->subtract_background();
		p->multiply(pmask);
		
		p->fft();
		
		if ( pref ) {
			R = p->merge_amplitudes_and_phases(pref, res, lo_res);
		} else {
			R = p->merge_amplitudes_and_phases(pamp);
			p->fspace_bandpass(res, 1e10, 0);
		}

		if ( verbose )
			cout << i+1 << tab << res << tab << R << endl;
		
		p->fft_back();
		
		if ( sym.point() > 101 ) p->symmetrize(sym, 1);
	
		p->statistics();
		
		if ( (i+1)%iout == 0 && (i+1) != imax ) {
			filename = outname.pre_rev('.') + Bstring(i+1, "_%04d.") + outname.post_rev('.');
			pint = p->copy();
			pint->change_type(nudatatype);
			write_img(filename, pint, 0);
			delete pint;
		}
	} 

	if ( optind < argc ) {
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete pamp;
	delete pmask;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

