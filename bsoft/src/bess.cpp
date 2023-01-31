/**
@file	bess.cpp
@brief	Electron scattering simulation
@author	Bernard Heymann
@date	20190724 - 20220414

clang++ -o bin/bess src/bess.cpp -I. -I/usr/local/include -I/Users/bernard/b20/bsoft/include -L/Users/bernard/b20/bsoft/lib -lbsoft -std=c++11 -I/Users/bernard/b20/fftw-3.3.6-pl2/include 
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "model_map.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_util.h"
#include "ctf.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			img_ssnr(Bimage* p, long window, double res_hi, double sampling_ratio);

// Usage assistance
const char* use[] = {
" ",
"Usage: bess [options] input.pdb output.mrc",
"------------------------------------------",
"Simulates electron microscope imaging given a molecular structure.",
" ",
"Actions:",
"-center                  Center coordinates before calculations.",
"-View 0.3,-0.5,0.8,33    View to rotate the molecule to.",
"-dose 20,2.8,0.1         Simulate dose and radiation damage effects:",
"                         dose/fluence, critical dose parameters (all e/Å2).",
"-snr 5                   Estimate SSNR over a summation window.",
"-noabberrations          Do not apply aberrations.",
"-ewald upper             Apply Ewald sphere shift (upper, lower, combine).",
"-Bfactor 44,25.3         B-factor application: B-factor (A^2), high resolution limit (A,optional)",
"                         Multiplied in reciprocal space by exp(-B-factor/4 * s^2).",
"-back                    Transform the output back to real space.",
"-convert real            Convert the complex transform: real, imag, Amp, Int (default not).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-size 250,230            Output image size (default from coordinates).",
"-images 120              Output number of images (default 1).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-origin 136,123          Origin placement within image (default center).",
"-Defocus 1.2,1.0,47      Defocus average & deviation, and astigmatism angle (default 0,0,0).",
"-ratio 2.5               Radial sampling to Cartesian sampling ratio (default 1).",
" ",
#include "use_ctf.inc"
" ",
"Input:",
"-record ctf.json         Image acquisition parameter file.",
"-atoms atomprop.star     Input atom properties file.",
"-residues resprop.star   Input residue properties file.",
" ",
"Output:",
"-output particle.json    Output parameter file.",
" ",
NULL
};


int		main(int argc, char** argv)
{
	// Initialize variables
	int				center(0);				// Flag to center coordinates
	int				ab_flag(1);				// Flag to apply aberration weights
	int				ewald_flag(0);			// Flag to apply Ewald sphere shift: 1=one, 2=combine
	View2<double>	view;					// View to generate
	bool 			set_backtransform(0);	// Flag for back transformation
    ComplexConversion	conv(NoConversion);		// Conversion from complex transform
	Vector3<long>	size;					// Image size
	long			nimg(1);				// Number of images
	Vector3<double>	sam(1,1,1);    			// Sampling in angstrom/voxel side
	Vector3<double>	origin;					// Coordinate origin placement
	int				set_origin(0);			// Flag to set origin
	double			def_avg(0);				// In angstrom
	double			def_dev(0);				// In angstrom
	double			ast_angle(0);	 		// Used to limit astigmatism
	vector<double>	dose = {1,0};		// Accumulated dose and critical dose parameters
	double			sampling_ratio(1);		// For SNR estimation
	double			snr_window(0);			// Summing window for SNR estimation
	double			Bfactor(0);				// B-factor decay
	Bstring			elements;				// Selected elements
	Bstring			atom_select("all");		// Selection
	Bstring 		atompropfile;			// Atom properties file
	Bstring 		respropfile;			// Residue properties file
	Bstring			paramfile;
	string			recordfile;				// Image acquisition parameter file
	string			outfile;				// JSON parameter output file

	double			v;
	JSvalue			jsctf(JSobject);

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	
	// Find the record file first
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "record" ) {
			recordfile = curropt->filename().c_str();
			cout << "Reading " << recordfile << endl;
			JSparser	parser(recordfile);
			JSvalue		js = parser.parse();
			JSvalue		jsctf = js;
			if ( !jsctf.exists("Volt") ) {
				string		str("*['ctf']");
				jsctf = *(js(str)[0]);
			}
			cout << jsctf << endl;
			if ( !jsctf.exists("Volt") ) {
				cerr << "CTF parameters not found in file " << recordfile << endl;
				bexit(-1);
			}
		}
	}
	
	// Parse other options
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "noabberrations" ) ab_flag = 0;
		if ( curropt->tag == "ewald" ) {
			ewald_flag = 1;
			if ( curropt->value[0] == 'l' ) ewald_flag = -1;
			if ( curropt->value[0] == 'c' ) ewald_flag = 2;
		}
		if ( curropt->tag == "View" )
			view = curropt->view2();
		if ( curropt->tag == "snr" )
			if ( ( snr_window = curropt->value.integer() ) < 1 )
 				cerr << "-snr: A window larger than 0 must be specified." << endl;
		if ( curropt->tag == "back" )
			set_backtransform = 1;
		if ( curropt->tag == "convert" )
			conv = curropt->complex_conversion();
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "images" )
			if ( ( nimg = curropt->value.integer() ) < 1 )
				cerr << "-images: A number of images must be specified!" << endl;
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
#include "ctf.inc"
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->values(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: At least the defocus average must be specified!" << endl;
			else {
//				if ( fabs(def_avg) < 20 ) def_avg *= 1e4;	// Assume um
//				if ( fabs(def_dev) < 10 ) def_dev *= 1e4;	// Assume um
				ast_angle *= M_PI/180;				// Assume degrees
			}
		}
		if ( curropt->tag == "dose" ) {
			dose = curropt->value.split_into_doubles(",");
			if ( dose.size() < 1 )
				cerr << "-dose: At least a dose must be specified!" << endl;
		}
		if ( curropt->tag == "ratio" )
			if ( ( sampling_ratio = curropt->value.real() ) < 0.1 )
				cerr << "-ratio: A ratio must be specified!" << endl;
		if ( curropt->tag == "Bfactor" )
			if ( ( Bfactor = curropt->value.real() ) < 1 )
				cerr << "-Bfactor: A B factor value must be specified!" << endl;
		if ( curropt->tag == "atoms" )
			atompropfile = curropt->filename();
		if ( curropt->tag == "residues" )
			respropfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().c_str();
    }
	option_kill(option);
	
	if ( dose[0] < 1 && nimg > 1 ) {
		cerr << "Warning: Without a dose specified, only one image will be generated." << endl;
		nimg = 1;
	}
	
	double			ti = timer_start();

	CTFparam		cp = ctf_from_json(jsctf);
	cp.defocus_average(def_avg);
	cp.astigmatism(def_dev, ast_angle);
	
	Bstring			coorfile(argv[optind++]);
	Bmodel*			model = NULL;
	Bstring*		file_list = NULL;
	Vector3<double> box;
	
	if ( coorfile.length() ) {
		file_list = coorfile.split(",");
		model = read_model(file_list, paramfile);
		
		if ( !model ) {
			cerr << "Error: Problem with coordinate file " << coorfile << ", exiting!" << endl;
			bexit(-1);
		}
	}
	
	if ( !model ) {
		cerr << "No coordinates read!" << endl;
		bexit(-1);
	}
	
	string			imgfile;
	if ( optind < argc ) imgfile = argv[optind++];
	
	if ( center )
		models_shift(model, -models_center_of_coordinates(model));

	if ( view[2] < 1 )
		model_rotate(model, view);

	vector<Vector3<double>>	bounds = models_calculate_bounds(model);

	if ( size.volume() < 100 ) {
		model->calculate_bounds();
		size = sam * (bounds[1] - bounds[0]);
		size[2] = 1;
	}
	
	Bimage*			p = new Bimage(Float, TComplex, size, nimg);
	p->file_name(imgfile);
	p->sampling(sam);
	if ( set_origin != 1 ) origin = p->size()/2;
	p->fourier_type(Standard);
	p->image->view(view);

	Vector3<double>		start = -origin*sam;
	Vector3<double>		end = start + sam*size;
	if ( size[2] < 2 ) {	// Projection
		start[2] = bounds[0][2] - 1;
		end[2] = bounds[1][2] + 1;
	}
	long				nsel = models_select_within_bounds(model, start, end);
	if ( verbose ) {
		cout << "Model bounds:" << tab << bounds[0] << tab << bounds[1] << endl;
		cout << "Volume bounds:" << tab << start << tab << end << endl;
		cout << "Components selected:" << tab << nsel << endl;
	}

	model_show_selection(model);

		img_electron_scattering(model, p, cp, dose[0], dose[1], atompropfile, ab_flag | 2*(ewald_flag>0));
	
	if ( ewald_flag == 2  ) p->combine_ewald();

//	if ( dose.size() )
//		p->fspace_weigh_accumulated_dose(dose);

	if ( snr_window )
		img_ssnr(p, snr_window, p->sampling(0)[0], sampling_ratio);

	if ( Bfactor ) p->fspace_weigh_B_factor(Bfactor, 0);
	
	if ( set_backtransform ) {
		p->phase_shift(origin);
		p->fft(FFTW_BACKWARD, 0, conv);
	}

	jsctf = ctf_to_json(cp);
	(*p)["ctf"] = jsctf;
	
	// Write output file
	if ( imgfile.length() )
		write_img(imgfile, p, 0);
	
	if ( outfile.size() )
		p->meta_data().write(outfile);

 	delete model;
	if ( p ) delete p;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}



int			img_ssnr(Bimage* p, long window, double res_hi, double sampling_ratio)
{
	Bstring			filename(p->file_name());

	Bimage*			psum = p->fspace_subset_sums(window, 0);

	Bplot*			plot = psum->fspace_subset_ssnr(window, res_hi, sampling_ratio, 2);

	Bstring			psfile = filename.base() + "_subssnr.ps";
	
	ps_plot(psfile, plot);
	
	delete plot;

	psum->progressive_sum();
	psum->next->progressive_sum();

	plot = psum->fspace_subset_ssnr(window, res_hi, sampling_ratio, 3);

	psfile = filename.base() + "_progssnr.ps";
	
	ps_plot(psfile, plot);

	return 0;
}
