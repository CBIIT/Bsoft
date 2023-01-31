/**
@file	bfocal.cpp
@brief	Processing focal series
@author Bernard Heymann
@date	Created: 20220808
@date	Modified: 20230127
**/

#include "rwimg.h"
#include "mg_ctf_focal.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/* Usage assistance */
const char* use[] = {
" ",
"Usage: bfocal [options] img.mrc out.pif",
"---------------------------------------",
"Processes micrographs taken as a focal series.",
" ",
"Actions:",
"-series 1.2,1.4,0.05     Calculate a focal series function: start, end, increment.",
"-fit 1000                Fit a focal series CTF, using the given maximum iterations.",
"-zft                     Back transform z columns.",
"-center                  Center the simulated image.",
"-power                   Convert to intensities.",
"-sphere                  Extract sphere from 3D frequency space.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-size 200,300,50         Size of focal series (must be specified for focal series).",
"-resolution 2.3,120      Resolution limits for fitting (A).",
"-Bfactor 44              B-factor (default 0 A^2).",
"-Defocus 1.2,1.0,47      Defocus average & deviation, and astigmatism angle (default 2 um, 0, 0).",
"-Astigmatism 0.3,-34     Set defocus deviation and astigmatism angle.",
"-basetype 2              Baseline type: (default 1)",
"                         Type 1: Polynomial with 5 coefficients.",
"                         Type 2: Double Gaussian with 5 coefficients.",
"                         Type 3: EMAN style with 4 coefficients.",
"                         Types 4-6: 1-3 with gaussian fit of water ring (only for high resolution limit < 3 A).",
"-baseline 1,1.5,-2.6,12.9,30,118 Baseline type and 4 or 5 coefficients: (default 1,0,0,0,0)",
"-envtype 3               Envelope type: (default 4)",
"                         Type 1: Single gaussian (2 coefficients).",
"                         Type 2: Single gaussian with constant (3 coefficients).",
"                         Type 3: Double gaussian (4 coefficients).",
"                         Type 4: Double gaussian with constant (5 coefficients).",
"-envelope 28,-563        Envelope coefficients: default 1,-10.",
" ",
#include "use_ctf.inc"
" ",
"Input:",
"-json file.json          Input JSON file with CTF parameters.",
" ",
"Output:",
"-jsonout file.json       Output CTF parameters to a JSON file.",
//"-envelope env.mrc        The Ewald sphere envelope map.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
 	int				fit(0);						// Flag to fit the CTF and number of iterations
	int				ctf_flag(0);				// Flag for CTF image: bit 1 = center, bit 2 = power
	int				extract_sphere(0);			// Flag to extract a sphere based on the voltage
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<long>	size;						// Size for new focal series image
	double			hires(0), lores(0);			// High and low resolution limits
	double			Bfactor(0); 				// B-factor
	double			def_avg(0);					// In angstrom
	double			def_dev(0);					// In angstrom
	double			ast_angle(0);	 			// Used to limit astigmatism
	int				basetype(1);				// Baseline type: 1=poly, 2=double_gauss, 3=EMAN
	int				setbase(0);
	vector<double>	base = {1,0,0,0,0};			// Baseline coefficients
	int				envtype(1);					// Envelope type: 1=gauss, 2=gauss+, 3=double_gauss, 4=double_gauss+
	int				setenv(0);
	vector<double>	env = {1,-1,0,0,0};			// Envelope coefficients
 	double			def_ser_start(0), def_ser_end(0), def_ser_inc(0);	// In angstrom
	Bstring			jsin, jsout;				// JSON files

	double			v;
	JSvalue			jsctf(JSobject);

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "series" )
			if ( curropt->real_units(def_ser_start, def_ser_end, def_ser_inc) < 3 )
				cerr << "-series: All three values must be specified!" << endl;
		if ( curropt->tag == "fit" )
			if ( ( fit = curropt->integer() ) < 1 )
				cerr << "-fit: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "zft" ) ctf_flag |= 4;
		if ( curropt->tag == "center" ) ctf_flag |= 1;
		if ( curropt->tag == "power" ) ctf_flag |= 2;
		if ( curropt->tag == "sphere" ) extract_sphere = 1;
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
		if ( curropt->tag == "size" )
			size = curropt->size();
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
		if ( curropt->tag == "Bfactor" )
			if ( ( Bfactor = curropt->real() ) < 0.1 )
				cerr << "-Bfactor: A B-factor must be specified!" << endl;
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->real_units(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: At least the defocus average must be specified!" << endl;
			else {
				ast_angle *= M_PI/180;				// Assume degrees
			}
		}
		if ( curropt->tag == "Astigmatism" ) {
			if ( curropt->real_units(def_dev, ast_angle) < 1 )
				cerr << "-Astigmatism: A defocus value must be specified!" << endl;
			else {
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
		if ( curropt->tag == "basetype" ) {
			basetype = curropt->value.integer();
			if ( basetype < 1 || basetype > 6 ) {
				basetype = 1;
				cerr << "Warning: The baseline type must be 1, 2 or 3. Reset to 1." << endl;
			} else
				setbase = 1;
		}
		if ( curropt->tag == "baseline" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			for ( size_t i=0; i<d.size(); i++ ) base[i] = d[i];
			if ( d.size() < 1 )
				cerr << "-baseline: At least one coefficient must be specified!" << endl;
			else
				setbase = 2;
		}
		if ( curropt->tag == "envtype" ) {
			envtype = curropt->value.integer();
			if ( envtype < 1 || envtype > 4 ) {
				envtype = 4;
				cerr << "Warning: The envelope type must be 1, 2, 3 or 4. Reset to 4." << endl;
			} else
				setenv = 1;
		}
		if ( curropt->tag == "envelope" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			for ( size_t i=0; i<d.size(); i++ ) env[i] = d[i];
			if ( d.size() < 1 )
				cerr << "-envelope: At least an envelope amplitude must be specified!" << endl;
			else
				setenv = 2;
		}
#include "ctf.inc"
		if ( curropt->tag == "json" )
			jsin = curropt->filename();
		if ( curropt->tag == "jsonout" )
			jsout = curropt->filename();
    }
	option_kill(option);
	
	double			ti;
	if ( verbose & VERB_TIME )
		ti = timer_start();

	if ( jsin.length() ) jsctf = JSparser(jsin.c_str()).parse();
	CTFparam		cp = ctf_from_json(jsctf);
	cp.defocus_average(def_avg);
	cp.astigmatism(def_dev, ast_angle);
	if ( setbase ) cp.baseline(basetype, base);
	if ( setenv ) cp.envelope(envtype, env);

    Bimage* 		p = NULL;
    
    if ( size.volume() < 1 && optind < argc )
    	p = read_img(argv[optind++], 1, -1);
    
	if ( !p && size.volume() < 1 ) {
		cerr << "Error: No input file read or size specified!" << endl;
		bexit(-1);
	}
	
	if ( p ) {
		if ( sam.volume() ) p->sampling(sam);
		else sam = p->image->sampling();
		size = p->size();
	}
	
	if ( sam.volume() < 1e-3 ) sam = Vector3<double>(1,1,1);

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}

	if ( def_ser_inc ) {
		Bimage*		pctf = img_ctf_focal_series(cp, def_ser_start, def_ser_end, def_ser_inc,
						size, sam, 0, hires);
		delete p;
		p = pctf;
	}

	if ( fit ) {
		Bimage*	pfit = img_ctf_focal_fit(p, cp, hires, lores, Bfactor, fit);
		delete p;
		p = pfit;
	}

	if ( jsout.length() ) {
		jsctf = ctf_to_json(cp);
		jsctf.write(jsout.str());
	}

	if ( ctf_flag & 4 ) p->fftz(FFTW_BACKWARD, 1);
	if ( ctf_flag & 1 ) p->center_wrap();
	if ( ctf_flag & 2 ) p->complex_to_intensities();
	
	if ( extract_sphere ) {
		Bimage*		ps = img_fspace_extract_sphere(p, cp.volt());
		delete p;
		p = ps;
	}

    // Write an output file if a file name is given
    if ( optind < argc ) {
//		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

