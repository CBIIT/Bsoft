/**
@file	bcc.cpp
@brief	Program for calculating cross-correlations
@author Bernard Heymann
@date	Created: 19980805
@date 	Modified: 20220801
**/

#include "rwimg.h"
#include "Matrix.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bcc [options] file.spi file.mrc",
"--------------------------------------",
"Calculates various comparisons using fast Fourier transforms.",
"Actions:",
"	auto-correlation.",
"	centering by correlation with a phase-inverted image.",
"	cross-correlation with a reference image.",
"	phase-correlation with a reference image.",
"	cross-correlation and cross-validation with a reference image.",
" ",
"Actions:",
"-action phase            Type: auto, center, cross, phase, valid.",
//"-Autocorrelate           Auto-correlate, output auto-correlation map.",
//"-Crosscorrelate          Cross-correlate with reference file.mrc, output shifted image.",
//"-Validate                Cross-correlate and cross-validate with reference file.mrc, output shifted image.",
//"-Center                  Find center by cross-correlation with an inverted copy.",
"-halfshift               Shift correlation map by half the size to put the origin in the middle.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-resolution 4.5,130      Resolution range for correlation (default 0 - 1e6 angstrom).",
"-origin 0,-10,30         Set origin, used with -size option (default 0,0,0).",
"-size 22,33,14           Size for extraction or hkl input (voxels, default from data).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-limit 25                Shift search limit in cross-correlation map (voxels, default: quarter of x-dimension).",
//"-wrap                    Turn wrapping on (default off).",
" ",
"Input:",
"-Reference ref.mrc       Reference file for cross-correlation.",
"-Mask mask.tif           Reciprocal space mask to use for cross-correlation.",
" ",
"Output:",
"-Map ccmap.tif           Cross-correlation map.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
    int				type(0);				// 0=none; 1=AC; 2=center; 3=CC; 4=PC; 5=Validate
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
//    int				findcenter(0);			// Flag to find image center
//    int 			setauto(0);				// Flag to do auto-correlation
	int				halfshift(0);			// Flag to shift origin to middle
	double			search_limit(-1);		// Search limit in cross-correlation map
    Bstring			reffile;				// File for cross-correlation
//    Bstring			valfile;				// File for cross-validation
    Bstring			ccfile;					// Cross-correlation map output file
    Bstring			maskfile;				// Mask to use for cross-correlation
    
	Vector3<double>	origin;					// Origin
	int				set_origin(0);			// Flag for setting the origin
	Vector3<double>	sam;					// Sampling
//	int 			setwrap(0);				// Wrapping flag
    
	double			hires(0), lores(0);		// Limiting resolution range (hires must be > 0 to be set)

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "action" ) {
			if ( curropt->value[0] == 'a' ) type = 1;
			if ( curropt->value[0] == 'c' ) {
				if ( curropt->value[1] == 'e' )type = 2;
			 	else type = 3;
			}
			if ( curropt->value[0] == 'p' ) type = 4;
			if ( curropt->value[0] == 'v' ) type = 5;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
// 		if ( curropt->tag == "Autocorrelate" )
//       		setauto = 1;
// 		if ( curropt->tag == "Center" )
//       		findcenter = 1;
 		if ( curropt->tag == "halfshift" ) halfshift = 1;
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
  		if ( curropt->tag == "limit" )
			if ( ( search_limit = curropt->value.real() ) < 0 )
				cerr << "-limit: A search limit must be specified" << endl;
//		if ( curropt->tag == "wrap" )
//			setwrap = 1;
// 		if ( curropt->tag == "Crosscorrelate" ) type = 0;
//			reffile = curropt->filename();
 //		if ( curropt->tag == "Validate" )
//			valfile = curropt->filename();
 		if ( curropt->tag == "Reference" )
			reffile = curropt->filename();
 		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "Map" )
			ccfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();
	
	if ( type > 2 && reffile.length() < 1 ) {
		cerr << "Error: A reference file must be specified for this action!" << endl;
		bexit(-1);
	}
	
	// Allocate memory for the structure factors and image parameters
    Bimage*	 		p = NULL;
    Bimage*	 		pref = NULL;
    Bimage*	 		pmask = NULL;
    Bimage*	 		pcc = NULL;
    
    // Read the file
	if ( ( p = read_img(argv[optind++], 1, -1) ) == NULL ) {
		cerr << "File not read!" << endl;
		bexit(-1);
	}
	
    if ( reffile.length() ) {
    	pref = read_img(reffile, 1, -1);
		if ( !pref ) bexit(-1);
	}
	
	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, -1);
		if ( !pmask ) bexit(-1);
	}

	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	else if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( sam.volume() > 0 ) p->sampling(sam);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}

	
	if ( search_limit < 0 ) search_limit = p->sizeX()/4.0;
	
    // Do cross correlation or a simple Fourier transform
	long   			n;
	double			v, ccavg(0), ccstd(0);
    if ( type == 3 || type == 4 ) {
		if ( type == 3 )
			pcc = p->cross_correlate(pref, hires, lores, pmask);
		else
			pcc = p->phase_correlate(pref, hires, lores, pmask);
//		pcc = p->cross_correlate_fspace(pref, hires, lores, search_limit);
//		origin = pref->image->origin();
//		cout << "origin=" << pcc->image->origin() << endl;
		pcc->find_peak(search_limit, 0);
//		cout << "origin=" << pcc->image->origin() << endl;
		pcc->refine_peak();
//		cout << "origin=" << pcc->image->origin() << endl;
		for ( n=0; n<pcc->images(); n++ ) {
			p->origin(n, pref->image[n].origin() + pcc->image[n].origin());
			p->image[n].FOM(pcc->image[n].FOM());
		}
		if ( verbose ) {
			cout << "Image\t   dx\t   dy\t   dz\t  CC\t   P" << endl;
			for ( n=0; n<p->images(); n++ ) {
				v = p->image[n].FOM();
				ccavg += v;
				ccstd += v*v;
				cout << n+1 << tab << fixed << setprecision(4) << 
					pref->image[n].origin() - p->image[n].origin() << tab <<
					p->image[n].FOM() << tab << pcc->ccmap_confidence(n) << endl;
			}
			cout << endl;
			ccavg /= n;
			ccstd /= n;
			ccstd -= ccavg*ccavg;
			if ( ccstd > 1 ) ccstd = sqrt(ccstd);
			cout << "Correlation avg & stdev:         " << ccavg << tab << ccstd << endl << endl;
		}
		delete pref;
		if ( optind < argc ) p->center_wrap();
    } else if ( type == 5 ) {
		pcc = p->cross_correlate_validate(pref, pmask);
//		origin = pref->image->origin();
		pcc->find_peak(search_limit, 0);
		pcc->refine_peak();
		for ( n=0; n<pcc->images(); n++ )
			p->origin(n, pref->image[n].origin() + pcc->image[n].origin());
		delete pref;
		if ( verbose ) {
			cout << "Image\t   dx\t   dy\t   dz\t  CC\t   P" << endl;
			for ( n=0; n<p->images(); n++ ) {
				cout << n+1 << tab << fixed << setprecision(4) << 
					pref->image[n].origin() - p->image[n].origin() << tab <<
					p->image[n].FOM() << tab << pcc->ccmap_confidence(n) << endl;
			}
			cout << endl;
		}
		if ( optind < argc ) p->center_wrap();
	} else if ( type == 1 ) {
		p->auto_correlate(hires, lores);
		if ( halfshift ) {
			p->origin(0,0,0);
			p->center_wrap();
		}
		pcc = p;
	} else if ( type == 2 ) {
		p->find_center(pmask, hires, lores, search_limit, 0, 1);
		if ( halfshift ) p->center_wrap();
		pcc = p;
	}

	if ( pcc && ccfile.length() ) {
		if ( !pcc->next ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG bcc: writing pcc to " << ccfile << endl;
			if ( halfshift ) pcc->center_wrap();
			write_img(ccfile, pcc, 0);
		} else {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG bcc: writing pcc->next to " << ccfile << endl;
			if ( halfshift ) pcc->next->center_wrap();
			write_img(ccfile, pcc->next, 0);
		}
	}
	
	if ( optind < argc ) {
		p->change_type(nudatatype);
	    write_img(argv[optind], p, 0);
	}
	
	if ( p != pcc ) delete pcc;
	if ( pmask ) delete pmask;
	delete p;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

