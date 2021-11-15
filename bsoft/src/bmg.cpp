/**
@file	bmg.cpp
@brief	Manipulating project and micrograph structures
@author Bernard Heymann
@date	Created: 20020826
@date	Modified: 20210706
**/

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_particle_select.h"
#include "ps_micrograph.h"
#include "rwmg.h"
#include "mg_ctf.h"
#include "file_util.h"
#include "utilities.h"
#include "options.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmg [options] input.star/input.img [input.star/input.img]",
"----------------------------------------------------------------",
"Manipulates micrograph parameter files.",
"If the input files are images, a new parameter file will be created.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-getfield field1         Get a single field from a parameter file.",
"-getmicrograph mg1       Get a single micrograph from a parameter file.",
"-selected                Show selected particles.",
"-histogram 80,0.01       Histogram of magnification values of selected particles.",
"-add                     Add micrograph coordinates and particle origins.",
"-Euler                   Convert orientations from views to Euler angles.",
"-View                    Convert orientations from Euler angles to views.",
"-flip y                  Reverse the particle coordinates in the micrographs along the given axis.",
"-oriflip y               Reverse the particle origins along the given axis.",
"-toview -2.5,5.5,9,45    Transform particle orientations relative to this view.",
"-toeuler 45,30,10        Transform particle orientations relative to these Euler angles.",
"-remove                  Do not write non-selected micrographs into the parameter file.",
" ",
"Actions dealing with file references:",
"-check                   Check if files can be found.",
"-filenames               Delete file names of files not found.",
" ",
"Actions for generating parameters from images:",
"-extract part            Extract views from images: mg, frame, rec, part, fil.",
"-intensities             Calculate micrograph intensities.",
" ",
"Actions modifiying images:",
"-write                   Write particle views into images.",
"-bin 3                   Bin micrographs by given kernel size, new file names with a \"_b<n>\" insert.",
"-changepixelsize 3.68    Change the pixel size and all associated size parameters.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-numberperfield 3        Number of micrographs per field-of-view (default 1).",
"-fieldID dbgd            Set the field ID (multiple micrographs will get numbered ID's).",
"-mgpath dir/subdir       Set micrograph file paths.",
"-framepath dir/subdir    Set micrograph frames file paths.",
"-pspath dir/subdir       Set power spectrum file paths.",
"-partpath dir/subdir     Set particle file paths.",
"-filpath dir/subdir      Set filament file paths.",
"-samemgpart              Set micrograph and particle file names to be the same.",
"-Scansampling 3.6        Set scan sampling or pixel size (um/pixel).",
"-mgpixelsize 1.2,1.15    Set micrograph pixel size (angstrom/pixel, one value sets all).",
"-framepixelsize 1.2,1.1  Set frame pixel size (angstrom/pixel, one value sets all).",
"-partpixelsize 2.3,2.2,1 Set particle pixel size (angstrom/pixel, one value sets all).",
"-boxsize 55              Set particle box size (default from parameter file).",
"-mgorigin 1256,1256,1    Set micrograph origins (voxels).",
"-partorigin 215,215,1    Set particle origins (voxels).",
" ",
#include "use_dose.inc"
" ",
"Actions setting microscope and CTF parameters:",
"-Magnification 38500     Set magnification.",
"-Defocus 1.8,0.3,-34     Set defocus average, deviation, and astigmatism angle (um, degrees).",
"-Astigmatism 0.3,-34     Set defocus deviation and astigmatism angle (um, degrees).",
"-Cs 2.0                  Set spherical aberration, Cs (mm).",
"-Volt 100                Set acceleration voltage (kV).",
"-Amplitude 0.07          Set amplitude contrast (fraction)",
"-alpha 0.15              Set illumination half-angle (milliradians).",
"-basetype 2              Baseline type: (default 1)",
"                         Type 1: Polynomial with 5 coefficients.",
"                         Type 2: Double Gaussian with 5 coefficients.",
"                         Type 3: EMAN style with 4 coefficients.",
"-baseline 1,1.5,-2.6,12.9,30,118 Baseline type and 4 or 5 coefficients: (default 1,0,0,0,0)",
"-envtype 3               Envelope type: (default 4)",
"                         Type 1: Single gaussian (2 coefficients).",
"                         Type 2: Single gaussian with constant (3 coefficients).",
"                         Type 3: Double gaussian (4 coefficients).",
"                         Type 4: Double gaussian with constant (5 coefficients).",
"-reset defocus           Reset particle defocus to micrograph defocus.",
" ",
"Input:",
"-input input.star        Input master parameter file to merge in data blocks from other input files.",
"-revert part.star,15     Replace file names with those in this file,",
"                         with a flag indicating: 1=mg, 2=ps, 4=rec, 8=part.",
"-replacectf ctf.star     Replace CTF parameters with those in this file.",
"-replacepart part.star   Replace particle parameters with those in this file.",
"-validate mg.xsd         Validate XML input file(s) with this schema.",
"-ppx update              Process files in directory \"ppx\":",
"                         update: Update particle records in the parameter file.",
"                         list: List existing ppx files.",
"                         absent: List absent ppx files.",
" ",
"Output:",
"-output output.star      Output parameter file (required any time \"-split\" is used to indicate the type).",
"-split 3                 Split micrograph or reconstruction data blocks into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": micrograph ID's are used as file names (no -output required).",
"                         Argument: \"field\": field ID's are used as file names (no -output required).",
"-dump file.txt           Dump all the particle origins and views.",
"-Postscript ori.ps       Plot the micrograph origins.",
" ",
NULL
};

int main (int argc, char **argv)
{
	// Initialize optional variables
	int				use_rec(0);				// Flag to process reconstructions
	int				fom_index(0);			// Index of FOM to act on
	Bstring			mgpath;					// Micrograph file path
	Bstring			framepath;				// Micrograph frames file path
	Bstring			pspath;					// Power spectrum file path
	Bstring			partpath;				// Particle file path
	Bstring			filpath;				// Filament file path
	Bstring			field_id;				// User-specified field ID
	Bstring			get_field;				// Field ID to get
	Bstring			get_mg;					// Micrograph ID to get
//	Bstring			magfile;				// STAR format file with map magnification parameters
	int				show_selected(0);		// Don't show selected by default
	int				bins(0);				// Bins must be >0 to get histogram
	double			histo_inc(0.01);		// Histogram increment
	int 			nseries(0);				// Number of micrographs in a series from a field
	int 			add_origins(0);			// Add origins to micrograph coordinates
	int				bin(0);					// Binning for micrographs
	Bstring			extract;				// Create records from images
	int 			write_views(0);			// Use particle views from parameter file
	int 			mg_part_flag(0);		// Different mg & part file names
	double			magnification(0);		// Magnification
	double			sampling(0);			// Take scan sampling from input files
	Vector3<double>	mg_pixel;				// Take micrograph pixel size from input files
	Vector3<double>	frame_pixel;			// Take micrograph pixel size from input files
	Vector3<double>	part_pixel;				// Take particle pixel size from input files
	Vector3<double>	new_pixel_size;			// New pixel size to rescale all length parameters
	Vector3<double>	mgorigin;				// Micrograph origin
	Vector3<double>	partorigin;				// Particle origin
	int				box_size(0);			// Particle box size
	double			def_avg(0);				// No defocus specified
	double			def_dev(0);				// Default no astigmatism
	double			ast_angle(0);			// Default astigmatism angle
	JSvalue			dose_frac(JSobject);	// Container for dose fractionation parameters
	double			Cs(0);					// Spherical aberration
	double			Volts(0);				// Acceleration voltage
	double			Amp(0);					// Amplitude contrast contribution
	double			alpha(0);				// Illumination half-angle (radians)
	int				basetype(1);			// Baseline type: 1=poly, 2=double_gauss, 3=EMAN
	double			base[5] = {1,0,0,0,0};	// Baseline coefficients
	int				envtype(4);				// Envelope type: 1=gauss, 2=gauss+, 3=double_gauss, 4=double_gauss+
	double			env[5] = {1,-10,0,0,0};	// Envelope coefficients
	int				setbase(0);
	int				setenv(0);
	int				setcohenv(0);			// Flag to set partial coherence envelope
	Bstring			reset;                  // Reset a particle parameter from its micrograph
	int 			to_euler(0);			// No conversion
	int 			to_views(0);			// No conversion
	int				loc_flip(0);			// No location flipping
	int				ori_flip(0);			// No origin flipping
	View			view;					// View for new reference
	Euler			euler;					// Euler angles for new reference
	Bstring			masterfile;				// Input master parameter file name
	Bstring			revertfilenames;		// Parameter file name with old file names
	int				name_flag(0);			// File names to reveert to old versions
	Bstring			replacectffile;			// Parameter file name with desired CTF parameters
	Bstring			replacepartfile;		// Parameter file name with desired particle parameters
	Bstring			xsdfile;				// XML schema file
	Bstring			ppx("none");			// Processing ppx files
	Bstring			outfile;				// Output parameter file name
	Bstring			dumpfile;				// File to dump particle info to
	Bstring			psfile;					// File to plot micrograph origins
	int				split(0);				// Output one big STAR file
	int				read_flags(0);			// Flags to pass to the parameter file reading function
	int				write_flags(0);			// Flags to pass to the parameter file writing function
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" ) use_rec = 1;
		if ( curropt->tag == "getfield" ) {
			get_field = curropt->value;
			if ( get_field.length() < 1 )
				cerr << "-getfield: A field ID must be specified!" << endl;
		}
		if ( curropt->tag == "getmicrograph" ) {
			get_mg = curropt->value;
			if ( get_mg.length() < 1 )
				cerr << "-getmicrograph: A micrograph ID must be specified!" << endl;
		}
		if ( curropt->tag == "selected" )
			show_selected = 1;
		if ( curropt->tag == "histogram" )
			if ( curropt->values(bins, histo_inc) < 1 ) {
				cerr << "-histogram: The number of bins must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "numberperfield" )
			if ( ( nseries = curropt->value.integer() ) < 1 ) {
				cerr << "-numberperfield: The number of micrographs in a field-of-view must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "add" )
			add_origins = 1;
		if ( curropt->tag == "toview" ) {
			view = curropt->view();
		}
		if ( curropt->tag == "toeuler" ) {
			euler = curropt->euler();
			view = euler.view();
        }
		if ( curropt->tag == "remove" ) write_flags |= 2;
		if ( curropt->tag == "check" ) read_flags |= 8;
		if ( curropt->tag == "filenames" ) read_flags |= 24;
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: A binning value must be specified!" << endl;
		if ( curropt->tag == "extract" ) extract = curropt->value;
		if ( curropt->tag == "intensities" ) read_flags |= 32;
		if ( curropt->tag == "write" )
			write_views = 1;
		if ( curropt->tag == "changepixelsize" )
			new_pixel_size = curropt->scale();
		if ( curropt->tag == "fieldID" ) {
			field_id = curropt->value;
			if ( field_id.length() < 1 )
				cerr << "-fieldID: A field ID must be specified!" << endl;
		}
		if ( curropt->tag == "mgpath" ) {
			mgpath = curropt->value;
			if ( mgpath.length() < 1 )
				cerr << "-mgpath: The micrograph file path must be specified!" << endl;
			else
				if ( mgpath[-1] != '/' ) mgpath += "/";
		}
		if ( curropt->tag == "framepath" ) {
			framepath = curropt->value;
			if ( framepath.length() < 1 )
				cerr << "-framepath: The micrograph frames file path must be specified!" << endl;
			else
				if ( framepath[-1] != '/' ) framepath += "/";
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
		if ( curropt->tag == "filpath" ) {
			filpath = curropt->value;
			if ( filpath.length() < 1 )
				cerr << "-filpath: The filament file path must be specified!" << endl;
			else
				if ( filpath[-1] != '/' ) filpath += "/";
		}
		if ( curropt->tag == "samemgpart" )
			mg_part_flag = 1;
		if ( curropt->tag == "Magnification" )
			if ( ( magnification = curropt->value.real() ) < 0.00001 )
				cerr << "-Magnification: The magnification must be specified!" << endl;
		if ( curropt->tag == "Scansampling" ) {
			if ( ( sampling = curropt->value.real() ) < 0.0001 )
				cerr << "-Scansampling: The scan sampling must be specified!" << endl;
			else if ( sampling < 1e3 ) sampling *= 1e4;			// Assume um
		}
		if ( curropt->tag == "mgpixelsize" )
			mg_pixel = curropt->scale();
		if ( curropt->tag == "framepixelsize" )
			frame_pixel = curropt->scale();
		if ( curropt->tag == "partpixelsize" )
			part_pixel = curropt->scale();
		if ( curropt->tag == "boxsize" )
        	if ( ( box_size = curropt->value.integer() ) < 1 )
				cerr << "-boxsize: A particle box size in pixels must be specified" << endl;
		if ( curropt->tag == "mgorigin" )
			mgorigin = curropt->origin();
		if ( curropt->tag == "partorigin" )
			partorigin = curropt->origin();
#include "dose.inc"
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->values(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: A defocus value must be specified!" << endl;
			else {
				if ( def_avg < 1e3 ) def_avg *= 1e4;			// Assume um
				if ( def_dev < 1e3 ) def_dev *= 1e4;			// Assume um
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
		if ( curropt->tag == "Astigmatism" ) {
			if ( curropt->values(def_dev, ast_angle) < 1 )
				cerr << "-Defocus: A defocus value must be specified!" << endl;
			else {
				if ( def_dev < 1e3 ) def_dev *= 1e4;			// Assume um
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
		if ( curropt->tag == "Cs" ) {
			if ( ( Cs = curropt->value.real() ) < 0.001 )
				cerr << "-Cs: A Cs value must be specified!" << endl;
			else if ( Cs < 1e3 ) Cs *= 1e7;						// Assume mm
		}
		if ( curropt->tag == "Volt" ) {
			if ( ( Volts = curropt->value.real() ) < 1 )
				cerr << "-Volt: A voltage must be specified!" << endl;
			else if ( Volts < 1e3 ) Volts *= 1e3;				// Assume kilovolts
		}
		if ( curropt->tag == "Amplitude" ) {
			if ( ( Amp = curropt->value.real() ) < 0.00001 )
				cerr << "-Amplitude: A fraction must be specified!" << endl;
			else if ( Amp > 1 ) Amp /= 100;						// Assume percentage
		}
		if ( curropt->tag == "alpha" ) {
			if ( ( alpha = curropt->value.real() ) < 0.00000001 )
				cerr << "-alpha: The illumination half-angle must be specified!" << endl;
			else
				if ( alpha > 0.01 ) alpha /= 1000;	// Assume milliradians
		}
		if ( curropt->tag == "basetype" ) {
			basetype = curropt->value.integer();
			if ( basetype < 1 || basetype > 3 ) {
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
		if ( curropt->tag == "reset" )
			reset = curropt->value;
		if ( curropt->tag == "Euler" )
			to_euler = 1;				// Convert orientations to Euler angles
		if ( curropt->tag == "View" )
			to_views = 1;				// Convert orientations to views
		if ( curropt->tag == "flip" ) {
			if ( curropt->value.contains("x") ) loc_flip += 1;
			if ( curropt->value.contains("y") ) loc_flip += 2;
			if ( curropt->value.contains("z") ) loc_flip += 4;
		}
		if ( curropt->tag == "oriflip" ) {
			if ( curropt->value.contains("x") ) ori_flip += 1;
			if ( curropt->value.contains("y") ) ori_flip += 2;
			if ( curropt->value.contains("z") ) ori_flip += 4;
		}
		if ( curropt->tag == "input" )
			masterfile = curropt->filename();
		if ( curropt->tag == "revert" ) {
			revertfilenames = curropt->value.pre(',');
			name_flag = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "replacectf" )
			replacectffile = curropt->filename();
		if ( curropt->tag == "replacepart" )
			replacepartfile = curropt->filename();
		if ( curropt->tag == "validate" )
			xsdfile = curropt->filename();
		if ( curropt->tag == "ppx" ) ppx = curropt->value;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( curropt->value.contains("field") || curropt->value.contains("FIELD") ) split = -9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
		}
		if ( curropt->tag == "dump" )
			dumpfile = curropt->filename();
    }
	option_kill(option);
	
	if ( pspath.length() < 1 && mgpath.length() ) pspath = mgpath;

	Bproject*		project = NULL;
	Bproject*		project2 = NULL;
	
	if ( masterfile.length() )
		project = read_project(masterfile, read_flags);
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	if ( file_type(*file_list) == Micrograph )
		project2 = read_project(file_list, xsdfile, read_flags);
	else
		project2 = project_create_from_images(file_list, extract);
		
	string_kill(file_list);

//	cout << "---" << project2->filename << endl;
	
	if ( project ) {
		project_update(project, project2, fom_index);
		project_kill(project2);
	} else
		project = project2;
	
	if ( use_rec ) project->select = 1;
	
	if ( revertfilenames.length() ) {
		project2 = read_project(revertfilenames, xsdfile, read_flags);
		project_revert_filenames(project, project2, name_flag);
		project_kill(project2);
	}
	
	if ( replacectffile.length() ) {
		project2 = read_project(replacectffile, xsdfile, read_flags);
		project_merge_CTF_parameters(project, project2);
		project_kill(project2);
	}
	
	if ( replacepartfile.length() ) {
		project2 = read_project(replacepartfile, xsdfile, read_flags);
		project_merge_part_parameters(project, project2);
		project_kill(project2);
	}
	
	if ( ppx[0] == 'u' ) project_update_from_ppx(project);
	else if ( ppx[0] == 'l' ) {
		part_set_sequential(project);
		project_list_ppx(project, 2);
	} else if ( ppx[0] == 'a' ) {
		part_set_sequential(project);
		project_list_ppx(project, 1);
	}
	
	if ( get_field.length() )
		project_select_field(project, get_field);
	
	if ( get_mg.length() )
		project_select_micrograph(project, get_mg);
	
	if ( reset.length() ) project_reset(project, reset);
	
	if ( bin > 1 ) {
		project_bin_micrographs(project, bin, mgpath, partpath);
	} else if ( new_pixel_size[0] ) {
		project_change_pixel_size(project, new_pixel_size, mgpath, partpath);
	} else {
		if ( mgpath.length() )
			project_set_micrograph_path(project, mgpath);
		if ( framepath.length() )
			project_set_frame_path(project, framepath);
		if ( pspath.length() )
			project_set_powerspectrum_path(project, pspath);
		if ( partpath.length() )
			project_set_particle_path(project, partpath);
		if ( filpath.length() )
			project_set_filament_path(project, filpath);
	}
	
	if ( mg_part_flag )
		project_equal_mg_part_files(project);

	if ( extract.contains("part") )
		project_set_views_from_images(project);

	if ( write_views )
		project_set_views_in_images(project);

	if ( nseries )
		project_set_field_id(project, nseries, field_id);

	if ( magnification )
		project_set_magnification(project, magnification);
	
	if ( sampling )
		project_set_scan_sampling(project, sampling);
	
	if ( mgorigin[0] )
		project_set_micrograph_origins(project, mgorigin);
		
	if ( partorigin[0] )
		project_set_particle_origins(project, partorigin);
		
	if ( mg_pixel[0] )
		project_set_mg_pixel_size(project, mg_pixel);
		
	if ( frame_pixel[0] )
		project_set_frame_pixel_size(project, frame_pixel);
		
	if ( part_pixel[0] )
		project_set_part_pixel_size(project, part_pixel);
		
	if ( box_size )
		project_set_particle_box_size(project, box_size);
	
	if ( def_avg )
		project_set_defocus(project, def_avg, def_dev, ast_angle);
	else if ( def_dev )
		project_set_astigmatism(project, def_dev, ast_angle);
		
	if ( dose_frac.size() )
		project_set_dose(project, dose_frac);
	
	if ( Volts )
		project_set_volts(project, Volts);
	
	if ( Cs )
		project_set_Cs(project, Cs);
		
	if ( Amp )
		project_set_amp_shift(project, sin(Amp));
	
	if ( alpha )
		project_set_alpha(project, alpha);

	if ( setbase )
		project_set_baseline(project, basetype, base);
		
	if ( setenv )
		project_set_envelope(project, envtype, env);
	
	if ( setcohenv )
		project_set_coherence_envelope(project);
	
	if ( add_origins )
		project_add_origins_to_coords(project);
	
	if ( loc_flip )
		project_flip_particle_coordinates(project, loc_flip);

	if ( ori_flip )
		project_flip_origins(project, ori_flip);
	
	if ( view.vector_size() )
		project_rotate_particle_views(project, view);
		
	if ( verbose & (VERB_RESULT|VERB_LABEL) )
		project_show_hierarchy(project);
	
	if ( show_selected )
		project_show_selected(project);
		
	if ( bins )
		project_show_mag_histogram(project, bins, histo_inc);
		
	if ( to_views && !to_euler )
		project->euler_flag = 0;
	
	if ( to_euler && !to_views )
		project->euler_flag = 1;
	
	if ( psfile.length() ) {
		Bstring		title("Micrograph origins");
		ps_mg_origins(psfile, title, project);
	}
	
	if ( dumpfile.length() ) project_dump(project, dumpfile);

//	
	if ( project ) {
		if ( split < 0 )
			project_split_field_write(project);
		else if ( outfile.length() || split == 9 ) {
			project->split = split;
			write_project(outfile, project, write_flags);
		}
	}
	
	project_kill(project);

	bexit(0);
}
