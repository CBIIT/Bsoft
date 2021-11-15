/**
@file	bsim.cpp
@brief	A tool to generate a project with random orientations for a molecule and to simulate TEM images from given orientations.
@author Bernard Heymann
@date	Created: 20030805
@date 	Modified: 20041021
**/

#include "mg_processing.h"
#include "mg_multislice.h"
#include "rwmg.h"
#include "rwmolecule.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bsim [options] in.star [in2.star]",
"----------------------------------------",
"Simulates TEM images from atomic coordinates.",
"The program has different entry and exit points:",
"A new STAR file can be generated (use -generate or -symmetry option and no input STAR file)",
"   or parameters can be read in from existing STAR files.",
"The provision of a coordinate file leads to calculation of potential slices.",
"Input of STAR file(s) without a coordinate file implies that the potentials",
"   have already been calculated and will be used for multislice simulation.",
" ",
"Actions:",
"-generate 10,3,50        Generate a project, specifying number of fields, number of",
"                         micrographs per field and number particles per micrograph.",
"-symmetry C5             Generates images based on a grid in the asymmetric unit.",
"-coordinates file.pdb    Input coordinate file, generate simulated images in project.",
"-image                   Generate a final image and parameter file (default not).",
"-poisson                 Add poisson noise (default not).",
"-gaussian 0.2            Add gaussian noise with the given signal-to-noise ratio.",
"-MTF 17.5                Set the MTF decay constant (required to impose MTF).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Parameters for generation:",
"-translate 2.4           Standard deviation for random translation (default 0).",
"-origin 120              Set the particle image origin (default 50 pixels).",
"-angles 8,6              Step sizes for theta and phi in the asymmetric unit, one value sets both.",
"                         Use with the -symmetry option (default 10,10).",
"-Pixelsize 5.5           Set the particle pixel size (default 1 angstrom/pixel).",
"-electrondose 35.5       Set the electron dose (default 20 electrons/angstrom).",
"-Cs 2.0                  Set the spherical aberration, Cs (default 2 mm).",
"-Cc 2.0                  Set the chromatic aberration, Cc (default 2.0 mm).",
"-alpha 0.15              Set the illumination half-angle (default 0.1 milliradians).",
"-energyspread 2e-5       Set the effective energy spread (default 1e-5).",
"-Volt 100                Set the acceleration voltage (default 120 kV).",
"-Amplitude 0.07          Set the amplitude contrast (fraction, default 0.07).",
"-Defocus 0.5,3.4         Defocus range (required to apply CTF).",
"-fieldbasename field     Field-of-view base name (default field).",
"-fieldnumber 5           Field starting number (default 1).",
"-mgbasename mg           Micrograph base name (default mg).",
"-mgnumber 23             Micrograph starting number (default 1).",
"-partbasename part       Particle base name (default part).",
"-partnumber 158          Particle starting number (default 1).",
" ",
"Parameters for simulation:",
"-fieldname field         Field-of-view name (for selecting a specific field).",
"-mgname mg               Micrograph base name (for selecting a specific micrograph).",
"-partselect 132          Particle selection for potential generation (default all).",
"-type real               Type of potential calculation (default reciprocal, other gauss).",
"-size 100,80,120         Size of enclosing box (default based on input coordinates).",
"-thickness 20            Slice thickness in angstrom (default 10 A).",
"-resolution 8.5          Resolution limits (default Nyquest frequency).",
"-Bfactor 24.5            Set the overall temperature factor (default 0 A^2).",
" ",
"Input:",
"-parameters param.star   Parameter file name for electron scattering curves (default atom_prop.star).",
"-Water water.pdb         Block of water to embed the molecule into (default none).",
" ",
"Output:",
"-output mg.star          Parameter file (use with generate option).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int				genfield(0), genmg(0), genpart(0);	// Number of randomly oriented molecules to generate
	Bstring			symmetry_string;				// Default: none
	double			theta_step(M_PI*10.0/180.0);	// Angular step size for theta
	double			phi_step(M_PI*10.0/180.0);		// Angular step size for phi
	double			tsigma(0);				// No random shift in origins
	double			img_origin(0);			// Particle image origin
	Vector3<double>	pixel_size;				// Take pixel size from input files
	double			dose(20);				// Electron dose
	double			Cs(2e7);				// Spherical aberration
	double			Cc(2e7);				// Chromatic aberration (angstrom)
	double			Volts(120000);			// Acceleration voltage
	double			Amp(0.07);				// Amplitude contrast contribution
	double			def_min(0), def_max(0);	// Defocus range (angstrom) default none
//	DataType 		datatype = Float;		// Default data type for projection
	int				pottype(0);				// Reciprocal space potential calculation
	Vector3<long>	size;					// Simulated box size
	double			thickness(10);			// Slice thickness in angstrom
	double			resolution(0);			// High resolution limit
	double			Bfactor(0);				// Overall temperature factor
	double			alpha(0.0001);			// Illumination half-angle (radians)
	double			energy_spread(1e-5);	// Effective energy spread
	int				image_flag(0);			// Flag to generate the final image and parameter file
	int				poisson(0);				// No poisson noise added
	double			gauss(1e10);			// Gaussian signal-to-noise ratio, > 1000 means no noise
	double			kmtf(0);				// MTF decay constant
    Bstring    		atom_select("ALL");
	Bstring			fieldbase("field");		// Field-of-view base name
	Bstring			fieldname;				// Selected field-of-view name
	Bstring			mgbase("mg");			// Micrograph base name
	Bstring			mgname;					// Selected micrograph name
	Bstring			partbase("part");		// Particle base name
	int				fieldnumber(1);			// Starting field number
	int				mgnumber(1);			// Starting micrograph number
	int				partnumber(1);			// Starting particle ID
	int				partselect(0);			// Selected particle ID
	Bstring			coorfile;				// Input coordinate file
	Bstring			waterfile;				// Input solvent coordinate file
	Bstring			paramfile;				// Use default parameter file for atomic properties
	Bstring			outfile;				// Output parameter file
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "generate" )
        	if ( curropt->values(genfield, genmg, genpart) < 3 )
				cerr << "-generate: Three numbers must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				if ( i < 2 ) phi_step = theta_step;
				else phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "translate" )
			if ( ( tsigma = curropt->value.real() ) < 1e-30 )
				cerr << "-translate: A translation standard deviation must be specified!" << endl;
		if ( curropt->tag == "origin" )
			if ( ( img_origin = curropt->value.real() ) < 1e-30 )
				cerr << "-origin: The image origin must be specified!" << endl;
		if ( curropt->tag == "Pixelsize" )
			pixel_size = curropt->scale();
//			if ( ( pixel_size = curropt->value.real() ) < 1e-30 )
//				cerr << "-Pixelsize: The particle pixel size must be specified!" << endl;
		if ( curropt->tag == "electrondose" )
			if ( ( dose = curropt->value.real() ) < 1e-30 )
				cerr << "-electrondose: The electron dose must be specified!" << endl;
		if ( curropt->tag == "Cs" ) {
			if ( ( Cs = curropt->value.real() ) < 1e-30 )
				cerr << "-Cs: A Cs value must be specified!" << endl;
			else if ( Cs < 1e3 ) Cs *= 1e7;						// Assume mm
		}
		if ( curropt->tag == "Cc" ) {
			if ( ( Cc = curropt->value.real() ) < 1e-30 )
				cerr << "-Cc: A Cc value must be specified!" << endl;
			else if ( Cc < 1e3 ) Cc *= 1e7;			// Assume mm
		}
		if ( curropt->tag == "alpha" ) {
			if ( ( alpha = curropt->value.real() ) < 1e-30 )
				cerr << "-alpha: The illumination half-angle must be specified!" << endl;
			else {
				if ( alpha > 0.01 ) alpha /= 1000;	// Assume milliradians
				if ( alpha < 0 ) alpha = 0.0001;
			}
		}
		if ( curropt->tag == "energyspread" ) {
			if ( ( energy_spread = curropt->value.real() ) < 1e-30 )
				cerr << "-energyspread: A energy spread must be specified!" << endl;
			else {
				if ( energy_spread > 0.001 ) energy_spread /= 1e5;	// Assume eV
				if ( energy_spread < 0 ) energy_spread = 0;
			}
		}
		if ( curropt->tag == "Volt" ) {
			if ( ( Volts = curropt->value.real() ) < 1e-30 )
				cerr << "-Volt: A voltage must be specified!" << endl;
			else if ( Volts < 1e3 ) Volts *= 1e3;				// Assume kilovolts
		}
		if ( curropt->tag == "Amplitude" ) {
			if ( ( Amp = curropt->value.real() ) < 1e-30 )
				cerr << "-Amplitude: A fraction must be specified!" << endl;
			else if ( Amp > 1 ) Amp /= 100;						// Assume percentage
		}
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->values(def_min, def_max) < 1 )
				cerr << "-Defocus: Both defocus minimum and maximum must be specified!" << endl;
			if ( def_max < def_min ) def_max = def_min;
			if ( def_min < 100 && def_max < 100 ) {	// Assume um
				def_min *= 1e4;
				def_max *= 1e4;
			}
		}
//		if ( curropt->tag == "datatype" )
//			datatype = curropt->datatype();
		if ( curropt->tag == "type" ) {
			pottype = 0;
			if ( curropt->value.contains("rea") ) pottype = 1;
			if ( curropt->value.contains("gau") ) pottype = 2;
		}
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "thickness" )
			if ( ( thickness = curropt->value.real() ) < 1e-30 )
				cerr << "-thickness: A slice thickness must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "Bfactor" )
			if ( ( Bfactor = curropt->value.real() ) < 0.001 )
				cerr << "-Bfactor: A temperature factor must be specified!" << endl;
		if ( curropt->tag == "image" )
			image_flag = 1;
		if ( curropt->tag == "poisson" )
			poisson = 1;
		if ( curropt->tag == "MTF" )
			if ( ( kmtf = curropt->value.real() ) < 1e-30 )
				cerr << "-MTF: The MTF decay constant must be specified!" << endl;
		if ( curropt->tag == "gaussian" )
			if ( ( gauss = curropt->value.real() ) < 1e-30 )
				cerr << "-gaussian: The signal-to-noise ratio must be specified!" << endl;
		if ( curropt->tag == "fieldbasename" ) {
			if ( curropt->value.length() > 0 ) fieldbase = curropt->value;
			else cerr << "-fieldbasename: A field-of-view base name must be specified!" << endl;
		}
		if ( curropt->tag == "fieldname" ) {
			if ( curropt->value.length() > 0 ) fieldname = curropt->value;
			else cerr << "-fieldname: A field-of-view name must be specified!" << endl;
		}
		if ( curropt->tag == "mgbasename" ) {
			if ( curropt->value.length() > 0 ) mgbase = curropt->value;
			else cerr << "-mgbasename: A micrograph base name must be specified!" << endl;
		}
		if ( curropt->tag == "mgname" ) {
			if ( curropt->value.length() > 0 ) mgname = curropt->value;
			else cerr << "-mgname: A micrograph name must be specified!" << endl;
		}
		if ( curropt->tag == "partbasename" ) {
			if ( curropt->value.length() > 0 ) partbase = curropt->value;
			else cerr << "-partbasename: A particle base name must be specified!" << endl;
		}
		if ( curropt->tag == "fieldnumber" )
			if ( ( fieldnumber = curropt->value.integer() ) < 1 )
				cerr << "-fieldnumber: A field-of-view starting number must be specified!" << endl;
		if ( curropt->tag == "mgnumber" )
			if ( ( mgnumber = curropt->value.integer() ) < 1 )
				cerr << "-mgnumber: A micrograph starting number must be specified!" << endl;
		if ( curropt->tag == "partnumber" )
			if ( ( partnumber = curropt->value.integer() ) < 1 )
				cerr << "-partnumber: A particle starting number must be specified!" << endl;
		if ( curropt->tag == "partselect" )
			if ( ( partselect = curropt->value.integer() ) < 1 )
				cerr << "-partselect: A particle number must be specified!" << endl;
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->filename();
		if ( curropt->tag == "Water" )
			waterfile = curropt->filename();
		if ( curropt->tag == "Parameter" )
			paramfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
    
	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	Bproject*		project = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( file_list ) {
		project = read_project(file_list);
		string_kill(file_list);
	}

	if ( img_origin <= 0 ) img_origin = size[0]/2;
	
	if ( !project ) {
		if ( symmetry_string.length() > 0 ) {
			project = project_generate_asu(symmetry_string, pixel_size, img_origin,
				theta_step, phi_step, Volts, Cs, Cc, alpha, energy_spread, Amp, 
				def_min, dose, mgbase, partbase);
		} else if ( genfield*genmg*genpart > 0 ) {
			project = project_generate(genfield, genmg, genpart, pixel_size, img_origin,
				Volts, Cs, Cc, alpha, energy_spread, Amp, def_min, def_max, dose, 
				tsigma, fieldbase, mgbase, partbase, fieldnumber, mgnumber, partnumber);
		}
	}
	
	if ( !project ) {
		cerr << "Error: No STAR file generated or read!" << endl;
		bexit(-1);
	}
	
	if ( outfile.c_str() ) {
		write_project(outfile, project, 0, 0);
	}
	
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		water = NULL;
	if ( coorfile.length() ) {
		molgroup = read_molecule(coorfile, atom_select, paramfile);
		if ( !molgroup ) {
			error_show(coorfile.c_str(), __FILE__, __LINE__);
			bexit(-1);
		}
		molecule_update_comment(molgroup, argc, argv);
		if ( waterfile.c_str() ) {
			water = read_molecule(waterfile, atom_select, paramfile);
			if ( !water ) {
				error_show(waterfile.c_str(), __FILE__, __LINE__);
				bexit(-1);
			}
		}
		project_generate_potential(molgroup, water, project, fieldname, mgname, partselect,
				size, thickness, resolution, Bfactor, pottype, paramfile);
	}
	
	Bstring		outfile2;
	Bstring		insert("_tfn.");
	if ( image_flag ) {
		// Generates an image with an average of one and the CTF applied
		// Thickness and resolution must correspond to the values used for potential calculation
		project_generate_image(project, thickness, resolution);
	
		if ( outfile.length() ) {
			// Adjust for dose and apply noise and envelope functions
			project_apply_distortions(project, poisson, gauss, kmtf);
//			outfile2 = insert_in_filename(outfile, insert, '_');
			outfile2 = outfile.pre_rev('.') + insert + outfile.post_rev('.');
			write_project(outfile2, project, 0, 0);
		}
	}
	
	molgroup_kill(molgroup);
	molgroup_kill(water);
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

