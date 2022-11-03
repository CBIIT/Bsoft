/**
@file	bmultislice.cpp
@brief	Multi-slice simulation of the electron imaging process.
@author Bernard Heymann
@date	Created: 20030118
@date	Modified: 20221024
**/

#include "rwmolecule.h"
#include "rwimg.h"
#include "mg_multislice.h"
#include "mg_ctf.h"
#include "model_map.h"
#include "model_util.h"
#include "molecule_to_map.h"
#include "mol_util.h"
#include "Matrix.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmultislice [options] input.pdb/input_pot.img output.img",
"---------------------------------------------------------------",
"Simulates electron imaging using a multi-slice approach.",
"The output file contains the projection image.",
"The calculation of the atomic potential can be done in real space (fast but less accurate).",
"or in reciprocal space (slow but more accurate), both using atomic electron scattering curves.",
" ",
"Actions:",
"-Defocus 1.2,1.0,47      Defocus average & deviation, and astigmatism angle (required to apply CTF).",
"-poisson                 Add poisson noise.",
"-gaussian 0.2            Add gaussian noise with the given signal-to-noise ratio.",
"-MTF 17.5                Set the MTF decay constant (required to impose MTF).",
"-exit                    Calculate exit wave and final image.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype b              Data type of the projection image (default floating point).",
"-size 100,80,70          Size (default automatic voxels).",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-sampling 1,1.3,1.1      Sampling in angstrom/voxel (default 1,1,1).",
"-thickness 20            Slice thickness in angstrom (default 10).",
"-resolution 8.5          Resolution limits (default Nyquest frequency).",
"-electrondose 35.5       Set the electron dose (default 20 angstrom/pixel).",
"-Bfactor 24.5            Set the overall temperature factor (default 0).",
"-type real               Type of potential calculation (default reciprocal, other gauss).",
" ",
#include "use_ctf.inc"
" ",
"Input:",
"-parameters param.star   Parameter file name for electron scattering curves(default atom_prop.star).",
" ",
"Output:",
"-potential pot.img       Output atomic potential image(s).",
"-wave wav.img            Output exit wave image.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Float);			// Default data type for projection
	Vector3<long>	size = {1024,1024,1};		// Simulated projection size
	Vector3<double>	ori;						// Coordinate origin placement
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam = {1,1,1};				// Sampling in angstrom
	int				set_sampling(0);			// Flag to set sampling
	int				exit_wave(0);				// Flag to calculate exit wave
	double			thickness(10);				// Slice thickness in angstrom
	double			resolution(0);				// High resolution limit
	double			def_avg(0);					// Defocus average (angstrom)
	double			def_dev(0);					// Defocus deviation (angstrom)
	double			ast_angle(-1);	 			// Astigmatism angle
	int				poisson(0);					// No poisson noise added
	double			gauss(1e10);				// Gaussian signal-to-noise ratio, > 1000 means no noise
	double			dose(20);					// Electron dose
	double			kmtf(0);					// MTF decay constant
	double			Bfactor(0);					// Overall temperature factor
	int				pottype(0);					// Reciprocal space potential calculation
	Bstring    		atom_select("all");
	Bstring			potentialfile;				// Atomic potential filename
	Bstring			exitwavefile;				// Exit wave filename
	Bstring			paramfile;					// Use default parameter file for atomic properties

	double			v;
	JSvalue			jsctf(JSobject);

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "exit" ) exit_wave = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				ori = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" ) {
			sam = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "thickness" )
			if ( ( thickness = curropt->value.real() ) < 0.001 )
				cerr << "-thickness: A slice thickness must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->values(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: At least the defocus average must be specified!" << endl;
			else {
				if ( def_avg < 100 ) def_avg *= 1e4;// Assume um
				if ( def_dev < 100 ) def_dev *= 1e4;// Assume um
				ast_angle *= M_PI/180;				// Assume degrees
			}
		}
#include "ctf.inc"
		if ( curropt->tag == "electrondose" )
			if ( ( dose = curropt->value.real() ) < 0.001 )
				cerr << "-electrondose: The electron dose must be specified!" << endl;
		if ( curropt->tag == "poisson" )
			poisson = 1;
		if ( curropt->tag == "gaussian" )
			if ( ( gauss = curropt->value.real() ) < 0.0000001 )
				cerr << "-gaussian: The signal-to-noise ratio must be specified!" << endl;
		if ( curropt->tag == "MTF" )
			if ( ( kmtf = curropt->value.real() ) < 0.001 )
				cerr << "-MTF: The MTF decay constant must be specified!" << endl;
		if ( curropt->tag == "Bfactor" )
			if ( ( Bfactor = curropt->value.real() ) < 0.001 )
				cerr << "-Bfactor: A temperature factor must be specified!" << endl;
		if ( curropt->tag == "type" ) {
			pottype = 0;
			if ( curropt->value.contains("rea") ) pottype = 1;
			if ( curropt->value.contains("gau") ) pottype = 2;
		}
		if ( curropt->tag == "Parameter" )
			paramfile = curropt->filename();
		if ( curropt->tag == "potential" )
			potentialfile = curropt->filename();
		if ( curropt->tag == "wave" )
			exitwavefile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	CTFparam		cp = ctf_from_json(jsctf);
	cp.defocus_average(def_avg);
	cp.astigmatism(def_dev, ast_angle);

	// Read the atomic structure file
	Bstring			filename(argv[optind]);
	Bmolgroup*		molgroup = NULL;
//	Bmolgroup*		molgroup = read_molecule(filename, atom_select, paramfile);
	Bmodel*			model = read_model(filename, paramfile);
	
	Bimage*			p = NULL; 
	Bimage*			ppot = NULL; 
	Bimage*			pgrate = NULL; 
	Bimage*			pwav = NULL;
	
	if ( set_origin == 2 ) ori = size/2;
	ori[2] = 0;

	if ( molgroup ) {
		ppot = img_calc_potential(molgroup, size, ori, sam, thickness,
			resolution, Bfactor, paramfile, pottype);
		if ( potentialfile.length() ) {
			pgrate = ppot->copy();
			ppot->complex_to_real();
			if ( verbose )
				cout << "Writing the atomic potential image(s)" << endl;
			ppot->change_type(nudatatype);
			write_img(potentialfile, ppot, 0);
			delete ppot;
		} else {
			pgrate = ppot;
		}
		optind++;
	} else if ( model ) {
		long					nslices(size[2]);
		sam[2] = thickness;
		if ( size[2] < 2 ) nslices = size[1]*sam[2]/thickness;
		ppot = new Bimage(Float, TComplex, size[0], size[1], 1, nslices);
		ppot->origin(ori);
		ppot->sampling(sam);
		ppot->fourier_type(Standard);
		img_electron_scattering(model, ppot, cp, 0, 0, paramfile, 0);
		ppot->phase_shift(ori);
		ppot->fft(FFTW_BACKWARD, 2);	// Full normalization 
		ppot->fourier_type(NoTransform);
		ppot->multiply(POTPREFAC);
		if ( potentialfile.length() ) {
			pgrate = ppot->copy();
			ppot->complex_to_real();
//			ppot->complex_to_intensities();
			if ( verbose )
				cout << "Writing the atomic potential image(s)" << endl;
			ppot->change_type(nudatatype);
			write_img(potentialfile, ppot, 0);
			delete ppot;
		} else {
			pgrate = ppot;
		}
		optind++;
	} else {
		pgrate = read_img(argv[optind++], 1, -1);
		if ( set_sampling ) pgrate->sampling(sam);
		pgrate->simple_to_complex();
	}
	
	if ( pgrate == NULL ) {
		cerr << "Error: No coordinate or image file given!" << endl;
		bexit(-1);
	}
	
	p = pgrate;

	double			dose_per_pixel = dose*p->sampling(0)[0]*p->sampling(0)[1];
	if ( exit_wave ) {
		img_calc_phase_grating(pgrate, cp.volt());

//			write_img("pg.mrc", pgrate, 0);

		p = img_calc_multi_slice(pgrate, thickness, cp.volt());
//		p = img_calc_multi_slice_correct(pgrate, thickness, cp.volt());
	
		
		if ( exitwavefile.length() ) {
			pwav = p->copy();
			pwav->fft(FFTW_BACKWARD, 2);
			pwav->complex_to_intensities();
			pwav->change_type(nudatatype);
			if ( verbose )
				cout << "Writing the exit wave image" << endl;
			write_img(exitwavefile, pwav, 0);
			delete pwav;
		}
		
//		if ( def_avg )
//			img_apply_complex_CTF(p, def_avg, def_dev, ast_angle, volt, Cs, Cc, sin(amp_fac), alpha, energy_spread);
		if ( def_avg )
			img_apply_complex_CTF(p, cp);
		
		p->fft(FFTW_BACKWARD, 2);
		if ( p->compound_type() == TComplex ) p->complex_to_intensities();
		p->fourier_type(NoTransform);
//		p->rescale(dose_per_pixel, 0);
		
		// The inherent SNR is calculated as the ratio of the image variance to the dose
		// where the dose gives the expected variance of the correlated poisson noise
		if ( verbose )
			cout << "Inherent signal-to-noise ratio: " << p->standard_deviation()*p->standard_deviation()/dose_per_pixel << endl;
		
		if ( poisson )
			p->noise_poisson(p->average());
		
		if ( gauss < 1000 )
			p->noise_gaussian(0, p->standard_deviation()/sqrt(gauss));
		
		if ( kmtf > 0 )
			p->fspace_weigh_B_factor(4*kmtf);
	}
	
	if ( optind < argc ) {
		if ( p->fourier_type() > NoTransform ) p->fft(FFTW_BACKWARD, 1, Intensity);
//		if ( p->compound_type() == TComplex ) p->complex_to_intensities();
		p->change_type(nudatatype);
		if ( verbose )
			cout << "Writing the simulated projection image" << endl;
		write_img(argv[optind], p, 0);
	}
	
	if ( p != pgrate ) delete pgrate;
	delete p;
	
	molgroup_kill(molgroup);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

