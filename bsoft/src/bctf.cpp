/**
@file	bctf.cpp
@brief	A program to determine and correct for the CTF (contrast transfer function) in electron micrographs.
@author Bernard Heymann
@date	Created: 19970715
@date	Modified: 20230127
**/

#include "rwimg.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "mg_ctf_sim.h"
#include "ps_ctf_plot.h"
#include "ps_micrograph.h" 
#include "mg_processing.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bctf [options] input.star/input.img [input2.star/output.img]",
"-------------------------------------------------------------------",
"Applies, determines or corrects for the CTF (contrast transfer function) in electron micrographs.",
" ",
"Actions:",
"-action correct          Action: 1=flip, 2=apply, 3=correct, 4=wienerfilter,",
"                         5=baseline, 6=baseline2, 7=baseflip, 8=basecorrect,",
"                         11=prepare, 12=fit, 13=prepfit (default 0=none).",
"-contrast                Inverts contrast by negating the CTF.",
"-invertaxis              Inverts the tilt axis to correct for improper setup.",
"-filter                  Filter extremes before doing anything with the image.",
"-background              Correct background after applying CTF (default not).",
"-fitastigmatism 2        Fit astigmatism for individual micrographs (1) or combined (2) (default not).",
"-isotropy                Assess the power spectrum isotropy.",
"-toparticle              Transfer micrograph CTF parameters to particle records (default not).",
//"-noaberration odd        Delete aberration weights (all, odd or even).",
"-statistics              Display CTF parameter statistics.",
" ",
"Simulations:",
"-ctf 512,512,10          Generate a CTF image of this size (z>1 for focal series).",
"-center                  Center the simulated image.",
"-power                   Convert to intensities.",
//"-series 1.2,1.4,0.05     Calculate a focal series function: start, end, increment.",
"-simulate                Simulate a tilted CTF for each x or y-line.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-mgpath dir/subdir       Set the micrograph file paths.",
"-pspath dir/subdir       Set the power spectrum file paths.",
"-partpath dir/subdir     Set the particle file paths.",
"-tile 1024,1024,1        Size of power spectrum generated during preparation (default 512,512,1).",
"-micrograph              Use micrograph file given in parameter file instead of particle file.",
"-frames                  Use micrograph frames file given in parameter file instead of particle file.",
"-sampling 1.5,1.5,1.5    Sampling (angstrom/pixel, a single value sets all three).",
"-resolution 27.5,125.3   Resolution limits for CTF determination and application (default 0.1,1e6).",
"-wiener 0.15             Wiener factor for CTF correction (default 0).",
" ",
#include "use_ctf.inc"
" ",
"Imaging parameters:",
"-Defocus 1.2,1.0,47      Defocus average & deviation, and astigmatism angle (default 2 um, 0, 0).",
"-Astigmatism 0.3,-34     Set defocus deviation and astigmatism angle.",
"-axis 74.7               Tilt axis angle relative to x-axis (default 0 or from parameter file).",
"-tilt 15                 Tilt angle (default from parameter file).",
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
"Fitting parameters:",
"-Range 0.6,3.5,0.2       Defocus minimum, maximum and increment for fitting (default 0.1,20,0.1 um).",
" ",
"Input:",
"-json file.json          Input JSON file with CTF parameters.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-jsonout file.json       Output CTF parameters to a JSON file.",
"-PSaverage file.mrc      Output average power spectrum.",
"-Postscript file.ps      Postscript output file name for CTF fit.",
"-Histogram file.ps       Postscript output file name for defocus histogram.",
"-Plotastigmatism file.ps Postscript output file name for astigmatism.",
"-Average file.ps         Postscript output file name for average of all CTF curves.",
"-Zeroes file.ps          Postscript output file name for defocus-vs-zeroes plot.",
"-Envelope file.ps        Envelope curve output file name.",
"-Pointspread file.ps     Point spread curve output file name.",
"-Function ctf.mrc        Contrast transfer function output file name.",
" ",
"Examples:",
"To apply the CTF:",
"bctf -verbose 7 -action apply -Defocus 2.3 input.img distorted.img",
" ",
"To prepare a power spectrum from a micrograph and to fit the CTF:",
"bctf -v 7 -datatype float -action prepfit -sam 1.842 -out new.star mg0001.tif mg0001_ps.pif",
"	This is a good way to generate the initial STAR file for each micrograph",
"	of a single particle analysis effort, as well as a power spectrum file.",
"	Subsequently, the automatic fit should be checked and corrected in bshow.",
" ",
"To correct for the CTF by simply flipping the phases in every second Thon ring:",
"bctf -verbose 7 -action flip -Defocus 1.4 part.pif part_flip.pif",
"	Applying the same flipping operation twice will recover the original image.",
" ",
"To correct for the CTF:",
"bctf -v 7 -action correct -Defocus 1.9 -wiener 0.1 part.pif part_corr.pif",
"	The equation used for the correction is:",
"		U = I*C/(C^2 + w)",
"	where:",
"		U: corrected image",
"		I: distorted image",
"		C: CTF",
"		w: Wiener factor, default 0.2",
" ",
"To correct for the CTF with baseline compensation:",
"bctf -v 7 -act baseline -Def 1.2 -wien 0.2 -base 1,1.5,-2.6,12.9,30,118 part.pif part_corr.pif",
"	The equation used for the correction is:",
"		U = I*C/(B*C^2 + w)",
"	where:",
"		U: corrected image",
"		I: distorted image",
"		C: CTF",
"		B: baseline",
"		w: Wiener factor, default 0.2",
"	Baseline equations:",
"		1.	a0 + a1*s + a2*s^2 + a3*s^3 + a4*s^4",
"		2.	a0 + a1*e^(a2*s^2) + a3*e^(a4*s^2)",
"		3.	a0 + a1*e^(a2*s^0.5 + a3*s^2)",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int				invert_axis(0);			// Flag to invert micrograph matrices
	int				toparticle(0);			// Flag to transfer parameters to particles
	Bstring			mgpath;					// Micrograph file path
	Bstring			pspath;					// Power spectrum file path
	Bstring			partpath;				// Particle file path
	Vector3<double>	sam;    				// Pixel size
	Vector3<long>	tile_size(512,512,1);	// Size of power spectrum
	double			resolution_lo(0);		// Low resolution limit
	double			resolution_hi(0);		// High resolution limit
	int 			action(0); 				// Default no CTF operation
	int				ctf_flag(0);			// Flag for CTF image: bit 1 = center, bit 2 = power
	double			wiener(0);				// Wiener for CTF correction
	double			def_avg(2e4);			// In angstrom
	double			def_dev(0);				// In angstrom
//	double			def_ser_start(0), def_ser_end(0), def_ser_inc(0);	// In angstrom
	double			ast_angle(-1);	 		// Used to limit astigmatism
	double			tilt_axis(0);			// Tilt axis angle in radians
	double			tilt_angle(0);			// Tilt angle start in radians
	int				isotropy(0);			// Flag to assess isotropy
//	int				noab(0);				// Flag to delete aberration weights
	int				stats(0);				// Flag to display CTF parameter statistics
	int				basetype(1);			// Baseline type: 1=poly, 2=double_gauss, 3=EMAN
	vector<double>	base = {1,0,0,0,0};		// Baseline coefficients
	int				envtype(1);				// Envelope type: 1=gauss, 2=gauss+, 3=double_gauss, 4=double_gauss+
	vector<double>	env = {1,-1,0,0,0};		// Envelope coefficients
	double			def_start(1e3), def_end(2e5), def_inc(1e3);	// Defocus fitting range and increment
	int				setDefocus(0);			// Flags
	int				setbase(0);
	int				setenv(0);
	int				flags(0);				// Flags for using mg (1), filter extremes (2), background (4), astigmatism (8), use frames (16), invert contrast (32)
	int				fitastig(0);			// Flag to fit astigmatism individually (1) or global (2)
	int				simulate(0);			// Flag to simulate a tilted micrograph
	Vector3<long>	size(1,1,1);			// Size of CTF output file
	Bstring			outfile;				// Output parameter file name
	Bstring			jsin, jsout;			// JSON files
	Bstring			psaverage;				// Output power spectrum average file name
	Bstring			psctffile;				// Output postscript file for CTF fit
	Bstring			pshisfile;				// Output postscript file for defocus histogram
	Bstring			psastfile;				// Astigmatism plot file
	Bstring			psavgfile;				// Output postscript file for average CTF
	Bstring			pszeroesfile;			// Output postscript file for zeroes
	Bstring			psenvfile;				// Output postscript file for envelope
	Bstring			pspointfile;			// Output postscript file for point spread function
	Bstring			funcfile;				// Output image file for CTF
	int				wave_aberration(0);		// Flag to output wave aberration function
	
	double			v;
	JSvalue			jsctf(JSobject);
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "action" )
			action = curropt->ctf_action();
		if ( curropt->tag == "invertaxis" )
			invert_axis = 1;
		if ( curropt->tag == "toparticle" ) toparticle = 1;
		if ( curropt->tag == "ctf" )
			size = curropt->size();
		if ( curropt->tag == "center" ) ctf_flag |= 1;
		if ( curropt->tag == "power" ) ctf_flag |= 2;
		if ( curropt->tag == "tile" ) {
			tile_size = curropt->size();
			if ( tile_size.volume() < 1 )
				cerr << "-tile: All three dimensions must be specified." << endl;
		}
		if ( curropt->tag == "micrograph" )
			flags |= 1;
		if ( curropt->tag == "frames" )
			flags |= 17;
		if ( curropt->tag == "filter" )
	        flags |= 2;
//		if ( curropt->tag == "background" )
//			flags |= 4;
		if ( curropt->tag == "contrast" )
			flags |= INVERT;
		if ( curropt->tag == "fitastigmatism" ) {
			fitastig = curropt->value.integer();
			if ( fitastig == 1 ) flags |= 8;
			else if ( fitastig == 3 ) flags |= 32;
		}
		if ( curropt->tag == "isotropy" )
			isotropy = 1;
		if ( curropt->tag == "simulate" )
			simulate = 1;
		if ( curropt->tag == "statistics" )
			stats = 1;
		if ( curropt->tag == "mgpath" ) {
			mgpath = curropt->value;
			if ( mgpath.length() < 1 )
				cerr << "-mgpath: The micrograph file path must be specified!" << endl;
			else
				if ( mgpath[-1] != '/' ) mgpath += "/";
		}
		if ( curropt->tag == "pspath" ) {
			pspath = curropt->value;
			if ( pspath.length() < 1 )
				cerr << "-pspath: The power spectrum file path must be specified!" << endl;
			else
				if ( pspath[-1] != '/' ) pspath += "/";
		}
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "axis" ) {
			tilt_axis = curropt->value.real();
			tilt_axis = angle_set_negPI_to_PI(tilt_axis*M_PI/180.0);
		}
		if ( curropt->tag == "tilt" ) {
    	    tilt_angle = curropt->value.real();
			tilt_angle = angle_set_negPI_to_PI(tilt_angle*M_PI/180.0);
        }
		if ( curropt->tag == "resolution" )
			if ( curropt->values(resolution_hi, resolution_lo) < 1 )
				cerr << "-resolution: At least one resolution limit must be specified!" << endl;
		if ( curropt->tag == "wiener" ) {
			if ( ( wiener = curropt->value.real() ) < 0.000001 )
				cerr << "-wiener: A Wiener factor must be specified!" << endl;
			else {
				if ( wiener < 0.01 ) wiener = 0.01;
//				if ( wiener > 1 ) wiener = 1;
			}
		}
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->real_units(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: At least the defocus average must be specified!" << endl;
			else {
				ast_angle *= M_PI/180;				// Assume degrees
				setDefocus = 1;
			}
		}
		if ( curropt->tag == "Astigmatism" ) {
			if ( curropt->real_units(def_dev, ast_angle) < 1 )
				cerr << "-Astigmatism: A defocus value must be specified!" << endl;
			else {
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
//		if ( curropt->tag == "series" )
//			if ( curropt->real_units(def_ser_start, def_ser_end, def_ser_inc) < 3 )
//				cerr << "-series: All three values must be specified!" << endl;
#include "ctf.inc"
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
		if ( curropt->tag == "Range" )
			if ( curropt->real_units(def_start, def_end, def_inc) < 2 )
				cerr << "-Range: At least two values must be specified!" << endl;
		if ( curropt->tag == "json" )
			jsin = curropt->filename();
		if ( curropt->tag == "jsonout" )
			jsout = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "PSaverage" ) {
			psaverage = curropt->filename();
			fitastig = 2;
		}
		if ( curropt->tag == "Postscript" )
			psctffile = curropt->filename();
		if ( curropt->tag == "Histogram" )
			pshisfile = curropt->filename();
		if ( curropt->tag == "Plotastigmatism" )
			psastfile = curropt->filename();
		if ( curropt->tag == "Average" )
			psavgfile = curropt->filename();
		if ( curropt->tag == "Zeroes" )
			pszeroesfile = curropt->filename();
		if ( curropt->tag == "Envelope" )
			psenvfile = curropt->filename();
		if ( curropt->tag == "Pointspread" )
			pspointfile = curropt->filename();
		if ( curropt->tag == "Function" ) {
			funcfile = curropt->filename();
//			wave_aberration = 1;
		}
    }
	option_kill(option);
	
	double			ti = timer_start();
	
	if ( jsin.length() ) jsctf = JSparser(jsin.c_str()).parse();
	CTFparam		cp = ctf_from_json(jsctf);
	cp.defocus_average(def_avg);
	cp.astigmatism(def_dev, ast_angle);
	if ( setbase ) cp.baseline(basetype, base);
	if ( setenv ) cp.envelope(envtype, env);

//			ctf_to_json(cp).write("test.json");

	if ( pszeroesfile.length() ) if ( pszeroesfile.contains(".ps") )
		ps_ctf_defocus_zeroes(pszeroesfile, cp.volt(), cp.Cs(), cp.amp_shift());

	Bimage*			p = NULL;

	if ( psenvfile.length() ) {
		cp.show();
		double		freq_step = (resolution_hi)? 0.001/resolution_hi: 0.0005;
		ps_ctf_plot(psenvfile, cp, 1000, freq_step);
	}

	if ( pspointfile.length() ) {
		cp.show();
		double		freq_step = (resolution_hi)? 0.001/resolution_hi: 0.0005;
		ps_point_spread(pspointfile, cp, 1000, freq_step);
	}

	if ( funcfile.length() ) {
		if ( size.volume() < 2  && optind < argc ) {
			if ( file_type(argv[optind]) == Image ) {
				p = read_img(argv[optind], 0, 0);
				if ( p ) {
					size = p->size();
					delete p;
				}
			}
		}
		if ( size.volume() > 1 ) {
			if ( action < 1 ) action = 2;	// Only CTF function - no baseline or envelope
/*			if ( def_ser_inc ) {
//				p = img_ctf_gradient(cp, def_ser_start, def_ser_end, def_ser_inc,
//						size, sam, resolution_lo, resolution_hi);
				p = img_ctf_focal_series(cp, def_ser_start, def_ser_end, def_ser_inc,
						size, sam, resolution_lo, resolution_hi);
			} else if ( wave_aberration ) {*/
			if ( wave_aberration ) {
				p = img_wave_aberration(cp, size, sam);
			} else if ( action > 2 ) {
				p = img_ctf_calculate(cp, action, wiener, size,
						sam, resolution_lo, resolution_hi);
			} else {
				p = img_ctf_calculate(cp, (action==1), wiener, size,
						sam, resolution_lo, resolution_hi);
			}
//			p->complex_to_real();
			if ( ctf_flag & 1 ) p->center_wrap();
//			if ( ctf_flag & 2 ) p->square();
			if ( ctf_flag & 2 ) p->complex_to_intensities();
			p->change_type(nudatatype);
			write_img(funcfile, p, 0);
			delete p;
		} else {
			cerr << "Error: A size needs to be given for the CTF function image!" << endl;
			bexit(-1);
		}
	}
	
	Bproject*			project = NULL;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	Bstring*			file_list = NULL;
	Bstring				basename, filename, newfilename;
	
	if ( optind < argc ) {
		if ( file_type(argv[optind]) == Micrograph ) {
			while ( optind < argc ) string_add(&file_list, argv[optind++]);
			project = read_project(file_list);
			string_kill(file_list);
			if ( project_count_mg_particles(project) < 1 ) flags |= 1;
			if ( project->field->mg->ctf ) {
				cp.update(project->field->mg->ctf);
				cp.show();
				cp.show_baseline();
				cp.show_envelope();
				jsctf = ctf_to_json(cp);
			}
		} else {
			filename =  argv[optind++];
			if ( simulate && fabs(tilt_angle) > 0.001 ) {
				p = read_img(filename, 1, -1);
				Bimage*		psim = img_ttf_simulate(p, cp, action, wiener, tilt_angle, 0, tilt_axis);
				if ( optind < argc )
					write_img(argv[optind++], psim, 0);
				delete p;
				delete psim;
			}
			p = read_img(filename, 0, -1);
			if ( !p ) bexit(-1);
			if ( nudatatype == Unknown_Type )
				nudatatype = p->data_type();
//			img_ctf_fit_prepare(p, 10);
			project = new Bproject;
//			basename = (Bstring(p->file_name())).base();
			basename = filename.base();
			if ( p->sizeZ() < 2 ) {	// Create a micrograph record
				field = field_add(&project->field, basename);
				mg = micrograph_add(&field->mg, basename);
				mg->ctf = new CTFparam;
				mg->ctf->update(cp);
				mg->pixel_size = p->sampling(0);
				mg->tilt_angle = tilt_angle;
				mg->tilt_axis = tilt_axis;
				mg->origin = p->image->origin();
				if ( action == 12 ) {
//					mg->fps = p->file_name();
					mg->fps = filename;
				} else if ( p->images() > 1 ) {
//					mg->fpart = p->file_name();
					mg->fpart = filename;
					for ( i=1; i<=p->images(); i++ ) {
						part = particle_add(&part, i);
						if ( !mg->part ) mg->part = part;
						part->mg = mg;
					}
				} else {
//					mg->fmg = p->file_name();
					mg->fmg = filename;
					flags |= 1;
				}
			} else {	// Create a reconstruction record
				project->select = 1;
				rec = reconstruction_add(&project->rec, basename);
				rec->ctf = new CTFparam;
				rec->ctf->update(cp);
				rec->voxel_size = p->sampling(0);
				rec->origin = p->image->origin();
				if ( p->images() > 1 ) {
//					rec->fpart = p->file_name();
					rec->fpart = filename;
					for ( i=1; i<=p->images(); i++ ) {
						part = particle_add(&part, i);
						if ( !rec->part ) rec->part = part;
					}
				} else {
//					rec->frec = p->file_name();
					rec->frec = filename;
					flags |= 1;
				}
			}
			if ( optind < argc ) newfilename = argv[optind];
			delete p;
			if ( access(newfilename.c_str(), F_OK) == 0 ) {
				cerr << "Error: The output file, " << newfilename << ", already exists!" << endl;
				bexit(-1);
			}
		}
	}
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		if ( jsout.length() )
			ctf_to_json(cp).write(jsout.c_str());
		bexit(-1);
	}

	if ( mgpath.length() )
		project_set_micrograph_path(project, mgpath);
	
	// Set the sampling for all the micrographs
	if ( sam[0] )
		project_set_mg_pixel_size(project, sam);

	project_update_ctf(project, jsctf);

	// Set the defocus values for all the micrographs
	if ( setDefocus )
		project_set_defocus(project, def_avg, def_dev, ast_angle);
	else if ( def_dev )
		project_set_astigmatism(project, def_dev, ast_angle);

	// Set the tile parameters for all the micrographs
	if ( tilt_angle )
		project_set_tilt(project, tilt_axis, tilt_angle);

	// Invert the geometry in case it is incorrectly set up
	if ( invert_axis )
		project_invert_tilt_axis(project);

	// Set other parameters
	if ( setbase == 1 )
		project_set_baseline_type(project, basetype);
	if ( setenv == 1 )
		project_set_envelope_type(project, envtype);
	if ( setbase > 1 )
		project_set_baseline(project, basetype, base);
	if ( setenv > 1 )
		project_set_envelope(project, envtype, env);

	if ( toparticle )
		project_CTF_to_part(project);
	
//	cout << " action = " << action << endl;
	
	// Determine the CTF in all micrographs
	if ( action > 0 && action < 11 )
		project_ctf(project, action, resolution_lo, resolution_hi, tile_size, wiener, 
				nudatatype, partpath, newfilename, flags);
	else if ( action > 10 )
		project_ctf_prepare(project, action, resolution_lo, resolution_hi, tile_size,
				def_start, def_end, def_inc, pspath, newfilename, flags);

	if ( isotropy )
		project_powerspectrum_isotropy(project, resolution_lo, resolution_hi);

	if ( jsout.length() )
		jsctf.write(jsout.c_str());

	if ( psctffile.length() )
		project_plot_ctf(project, psctffile);

	if ( pshisfile.length() )
		ps_defocus_histogram(pshisfile, project);
	
	if ( psastfile.length() )
		ps_astigmatism_plot(psastfile, project);
	
	if ( psavgfile.length() )
		project_ctf_average(project, psavgfile);
	
   // Write an average power spectrum if a file name is given
	Bimage*		ps = NULL;
 	if ( fitastig == 2 ) {
		cp.update(project->field->mg->ctf);
		JSvalue		jd = project_defocus_range(project);
		if ( jd["defocus_average"].real() > 10 ) cp.defocus_average(jd["defocus_average"].real());
		else cp.defocus_average(2e4);
		ps = project_powerspectrum_average(project, cp.defocus_average());
		if ( ps ) {
			img_ctf_fit_astigmatism(ps, 0, cp, resolution_lo, resolution_hi);
			if ( verbose )
				cout << "Astigmatism deviation and angle:     " << cp.defocus_deviation() << 
					tab << cp.astigmatism_angle()*180.0/M_PI << endl;
			project_set_astigmatism(project, cp.defocus_deviation(), cp.astigmatism_angle());
		} else {
			cerr << "Error: Failed to calculate an average power spectrum!" << endl;
		}
	}
	
	if ( stats ) project_ctf_statistics(project);
	
 	if ( ps && psaverage.length() )
    	write_img(psaverage, ps, 0);

	delete ps;

	// Write an output STAR format file if a name is given
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

