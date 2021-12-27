/**
@file	bmgconvert.cpp
@brief	Convert between micrograph parameter formats
@author Bernard Heymann
@author David Belnap
@author James Conway
@author Juha Huiskonen
@date	Created: 20061101
@date	Modified: 20211223

**/

#include "mg_processing.h"
#include "rwmg.h"
#include "rwmgSTAR.h"
#include "rwmgRELION.h"
#include "rwmgXML.h"
#include "mg_ctf.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

#define DATLINELENGTH 120

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bproject*	read_project_conv(Bstring* file_list);
int			write_project_conv(Bstring& filename, Bproject* project, int flags);
int			read_project_dat(Bstring& filename, Bproject* project);
int			write_project_dat(Bstring& filename, Bproject* project, int flags);
int			read_project_crd(Bstring& filename, Bproject* project);
int			write_project_crd(Bstring& filename, Bproject* project, int flags);
int			project_modify_parameters(Bproject* project, char flags);

// Usage assistence
const char* use[] = {
" ",
"Usage: bmgconvert [options] input.star [input2.star]",
"----------------------------------------------------",
"Converts micrograph parameter files.",
"Supported formats:",
"	PIC CRD files",
"	PFT3DR DAT files",
"	Relion STAR files",
"DAT and CRD files only hold one micrograph per file:",
"	multiple micrographs are written as numbered DAT or CRD files.",
"	use -mgselect to select only one micrograph.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1      Sampling (A/pixel, one value sets all).",
"-euler                   Convert views to Euler angles.",
"-omega                   Add or subtract 90 degrees to or from the OMEGA angle.",
"-magnification           Invert the magnification factor:  1/mag.",
"-mgpath dir/subdir       Set micrograph file paths.",
"-framepath dir/subdir    Set micrograph frames file paths.",
"-pspath dir/subdir       Set power spectrum file paths.",
"-partpath dir/subdir     Set particle file paths.",
"-filpath dir/subdir      Set filament file paths.",
"-extension pif           Set the particle image file format.",
"-remove                  Do not write non-selected micrographs into the parameter file.",
" ",
"Input:",
"-particles part.mrcs     Input file with a stack of particles (Relion).",
" ",
"Output:",
"-output file.star        Output parameter file.",
" ",
NULL
};

int			main(int argc, char** argv)
{
	// Initializing variables
	Vector3<double>	sam;				// A/pixel
	int				euler(0);			// Flag to specify Euler angles
	char			flags(0);			// Flags to modify some input parameters
	Bstring			partfile;			// Relion stack of particles file
	Bstring			mgpath;				// Micrograph file path
	Bstring			framepath;			// Micrograph frames file path
	Bstring			pspath;				// Power spectrum file path
	Bstring			partpath;			// Particle file path
	Bstring			filpath;			// Filament file path
	Bstring			partext;			// Particle image file format
	Bstring			outfile;			// Output parameter file
	int				write_flags(0);		// Flags to pass to the parameter file writing function

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
//			if ( ( sampling = curropt->value.real() ) < 0.1 )
//				cerr << "-sampling: A pixel size must be specified" << endl;
		if ( curropt->tag == "euler" )
			euler = 1;
		if ( curropt->tag == "magnification" )
			flags += 1;
		if ( curropt->tag == "omega" )
			flags += 2;
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
		if ( curropt->tag == "extension" ) {
			partext = curropt->value;
			if ( partext.length() < 1 )
				cerr << "-extension: The particle file extension must be specified!" << endl;
		}
		if ( curropt->tag == "particles" )
			partfile = curropt->filename();
		if ( curropt->tag == "remove" ) write_flags |= 2;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project_conv(file_list);
	
	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No project generated!" << endl;
		bexit(-1);
	}

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

	if ( sam[0] > 0.1 )
		project_set_mg_pixel_size(project, sam);

	if ( flags )
		project_modify_parameters(project, flags);
	
	if ( euler ) project->euler_flag = 1;
	
	if ( partfile.length() || partext.length() )
		project_split_particles(project, partfile, partpath, partext);
	
	// Write an output parameter file if a name is given
	if ( outfile.length() )
		write_project_conv(outfile, project, write_flags);
	
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

Bproject*	read_project_conv(Bstring* file_list)
{
	if ( !file_list ) {
		cerr << "Error: No micrograph parameter filename!" << endl;
		return NULL;
	}
	
	Bstring			ext;
	Bstring*		thisfile;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Parameter filenames:";
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	int					err(0);
	FileType			type;
	Bproject*			project = new Bproject;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	project->filename = *file_list;
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		for ( field = project->field; field && field->next; field = field->next ) ;
		if ( field ) for ( mg = field->mg; mg && mg->next; mg = mg->next ) ;
		for ( rec = project->rec; rec && rec->next; rec = rec->next) ;
		ext = thisfile->extension();
		err = -1;
		if ( ext.contains("star") ) {
			type = file_type(*thisfile);
			if ( type == Micrograph )
				err = read_project_star2(*thisfile, project);
			else if ( type == MgRelion )
				err = read_project_relion(*thisfile, project);
			else
				cerr << "Error: The STAR file is not a valid form Bsoft can read!" << endl;
		} else if ( ext.contains("xml") ) {
			err = read_project_xml(*thisfile, project);
		} else if ( ext.contains("dat") ) {
			err = read_project_dat(*thisfile, project);
		} else if ( ext.contains("crd") ) {
			err = read_project_crd(*thisfile, project);
		} else {
			cerr << "Error: Extension " << ext << " not valid for parameter files!" << endl;
		}
		if ( !err ) {
			if ( !field ) field = project->field;
			else field = field->next;
			if ( field ) {
				if ( !mg ) mg = field->mg;
				else mg = mg->next;
//				if ( mg ) field_resolve_file_access(field, mg, *thisfile);
			}
			if ( !rec ) rec = project->rec;
			else rec = rec->next;
//			if ( rec ) reconstruction_resolve_file_access(rec, *thisfile);
		} else {
			cerr << "	Parameter file " << *thisfile << " not read." << endl;
		}
	}

	if ( err < 0 ) {
		error_show(file_list->c_str(), __FILE__, __LINE__);
		project_kill(project);
		return NULL;
	}
	
	project_update_first_zero(project);
	
//	project_check(project);
	
//	if ( verbose & VERB_PROCESS )
//		project_display_counts(project);
	
	return project;
}

int			write_project_conv(Bstring& filename, Bproject* project, int flags)
{
	int				err(0);

	if ( filename.empty() ) {
		cerr << "Error: No micrograph parameter filename!" << endl;
		return -1;
	}
	
	project_update_first_zero(project);
	
//	project_check(project);
	
//	if ( verbose & VERB_PROCESS )
//		project_display_counts(project);
	
	Bstring			ext = filename.extension();

	if ( ext.contains("star") && project->split == 9 ) {
		err = write_project_star2(filename, project, flags & 2, flags & 4);
		return err;
	}
	
	if ( verbose & VERB_LABEL )
		cout << "Parameter filename: " << filename << endl;

	if ( ext.contains("dat") ) {
		err = write_project_dat(filename, project, flags);
    } else if ( ext.contains("crd") ) {
		err = write_project_crd(filename, project, flags);
	} else {
		err = write_project(filename, project, flags);
	}
	
	if ( err < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	return err;
}


/**
@brief 	Reading micrograph parameters from a DAT file.
@param 	&filename		file name (or comma-delimited list).
@param 	*project		initialized project structure.
@return int				error code (<0 means failure).
**/
int			read_project_dat(Bstring& filename, Bproject* project)
{
	int					err(0), nx(0), ny(0), n(0);
	char				inputLine[1024];
	int					unitsInHeaderI;
	Bstring				unitsInHeaderC;   // "angstrom", "nm" or "pixel"
	double				px, volt, amp, def, dev, ast, cs;

	Bstring				field_id = filename.base();
	Bstring				mg_id = field_id;
	Bfield*				field = field_add(&project->field, field_id);
	Bmicrograph*		mg = micrograph_add(&field->mg, mg_id);
	Bparticle*			part = NULL;
	CTFparam*			ctf = mg->ctf = new CTFparam;

	if ( verbose ) cout << "Reading a DAT file:             " << filename << endl;
	
	project->comment = "# DAT file: " + filename + "\n";

	/*========================
       Read DAT-format file
    ========================*/
	ifstream		fdat(filename.c_str());
	if ( fdat.fail() ) {
		cerr << "File " << filename << " does not exist!.  Exit program." << endl;
		bexit(-1);
	}
	
	fdat.getline(inputLine, 1023);
	mg->fpart = inputLine;
	mg->fpart = mg->fpart.no_lead_space();
	mg->fpart = mg->fpart.pre(' ');
	mg->fpart = mg->fpart.remove('\n');
    if ( mg->fpart.empty() )  {
		cerr << "ERROR reading line 1 of PFT-format DAT file.  Exit program" << endl;
		exit(-1);
    }
	
	// General micrograph parameters
	sscanf(strchr(inputLine+mg->fpart.length(), ' '), " %d %d %d", &nx, &ny, &n);
//	cout << "\"" << inputLine << "\"" << nx << " " << ny << " " << n << endl;
	
	fdat.getline(inputLine, 1023);
	if ( !fdat.good() )  {
        cout << "ERROR reading line 2 of PFT-format DAT file.  Exit program" << endl;
        exit(-1);
	}
	if ( sscanf(inputLine, "%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf\n", &px, &unitsInHeaderI,
				&volt, &amp, &def, &dev, &ast, &cs) == 8)  {
		switch(unitsInHeaderI)  {
			case 0:
				unitsInHeaderC = "pixels";
				break;
			case 1:
				unitsInHeaderC = "angstroms";
				break;
			case 2:
				unitsInHeaderC = "nm";
			default:
			break;
		}
	} else if ( sscanf(inputLine, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &px, &volt,
					&amp, &def, &dev, &ast, &cs) == 7)  {
	} else {
			cout << "ERROR interpreting line 2 of PFT-format DAT file.  Exit program" << endl;
			exit(-1);
    }

    if ( volt < 1000.0 ) volt *= 1000;	// Volt units for STAR-format file
    cs *= 10000000;						// Angstrom units for STAR-format file
	if ( def < 100 ) def *= 1e4;		// Angstrom units for STAR-format file
	if ( dev < 100 ) dev *= 1e4;		// Angstrom units for STAR-format file
	def = (def + dev)/2;
	dev -= def;
	if ( dev < 0 ) dev = -dev;
	
	if ( nx && ny ) mg->box_size = Vector3<int>(nx, ny, 1);
	mg->pixel_size = Vector3<double>(px, px, 1);
	ctf->volt(volt);
	ctf->amp_shift(asin(amp));
	ctf->defocus_average(def);
	ctf->defocus_deviation(dev);
	ctf->astigmatism_angle(ast*M_PI/180.0);
	ctf->Cs(cs);
	
    char				line[DATLINELENGTH];
	int					numpart(0), pid;
	float				theta, phi, omega, ox, oy, mag, fom1, fom2, fom3;
    
	while ( fdat.getline(line, DATLINELENGTH) && fdat.good() ) {
		if ( sscanf(line, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f", &pid, &theta, &phi, &omega, &ox, &oy, &mag, &fom1, &fom2, &fom3) > 6 )  {
			part = particle_add(&part, pid);
			if ( !mg->part ) mg->part = part;
			part->ori[0] = ox;
			part->ori[1] = oy;
			part->mag = mag;
			part->fom[0] = fom1;
			part->fom[1] = fom2;
			part->fom[2] = fom3;
			part->view = Euler(-omega*M_PI/180.0, 
					theta*M_PI/180.0, phi*M_PI/180.0).view();
			numpart++;
		}
    }
	
	fdat.close();         // close input DAT file

	if ( n && numpart != n )
		cerr << "Warning: The number of particles (" << numpart <<
			") does not agree with the number in the header (" << n << ")" << endl;
	
	if ( verbose )
		cout << "Number of particles read:       " << numpart << endl << endl;
	
	return err;
}

/**
@brief 	Writing micrograph parameters to a DAT file.
@param 	*filename		file name.
@param 	*project		project structure.
@param 	flags			write flags to select micrographs.
@return int				error code (<0 means failure).
**/
int			write_project_dat(Bstring& filename, Bproject* project, int flags)
{
	int					err(0);
	int 				mg_select(flags & 2);
	int					i, f, numpart(0);
	double				def_min, def_max, omega;
	Euler				euler;
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	CTFparam*			ctf;
	ofstream			fdat;
	Bstring				outname = filename;

	if ( verbose ) cout << "Writing a DAT file:             " << filename << endl;

	if ( project_count_micrographs(project) < 2 ) mg_select = 1;

	for ( i=1, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, i++ )  if ( mg_select < 1 || mg->select > 0 ) {
			if ( mg_select < 1 ) outname = filename.pre_rev('.') + Bstring(i, "_%06d.") + filename.post_rev('.');
			for ( part = mg->part; part; part = part->next ) numpart++;
			ctf = mg->ctf;
			if ( !ctf ) ctf = mg->ctf = new CTFparam;
			def_min = ctf->defocus_average() - ctf->defocus_deviation();
			def_max = ctf->defocus_average() + ctf->defocus_deviation();
			fdat.open(outname.c_str());
			if ( fdat.fail() )  {                       // make sure input file exists
				cerr << "Unable to create file " << filename << "!.  Exit program." << endl;
				bexit(-1);
			}
			fdat << mg->fpart << " " << mg->box_size[0] << " " << mg->box_size[1] << " " << numpart << endl;
			fdat << fixed << right << setprecision(4) << setw(9) << mg->pixel_size[0] << ", "
				<< 1 << ", "
				<< setprecision(1) << setw(8) << ctf->volt() << ", "
				<< setprecision(4) << setw(10) << sin(ctf->amp_shift()) << ", "
				<< setprecision(4) << setw(9) << def_max/1e4 << ", "
				<< setprecision(4) << setw(9) << def_min/1e4 << ", "
				<< setprecision(3) << setw(9) << ctf->astigmatism_angle()*180.0/M_PI << ", "
				<< setprecision(2) << setw(7) << ctf->Cs()/1e7 << endl;
			for ( part = mg->part; part; part = part->next ) {
				euler = Euler(part->view);
				omega = -euler.psi()*180.0/M_PI;
				if ( omega < 0 ) omega += 360;
				fdat << setw(7) << part->id << ", "
					<< setprecision(3) << setw(7) << euler.theta()*180.0/M_PI << ", "
					<< setprecision(3) << setw(7) << euler.phi()*180.0/M_PI << ", "
					<< setprecision(3) << setw(7) << omega << ", "
					<< setprecision(3) << setw(7) << part->ori[0] << ", "
					<< setprecision(3) << setw(7) << part->ori[1] << ", "
					<< setprecision(3) << setw(5) << part->mag;
				for ( f=0; f<3; f++ ) fdat << ", " << setprecision(3) << setw(5) << part->fom[f];
				fdat << endl;
			}
			fdat.close();
		}
	}

	if ( verbose )
		cout << "Number of particles written:    " << numpart << endl << endl;
	
	return err;
}

/**
@brief 	Reading micrograph parameters from a CRD file.
@param 	&filename		file name (or comma-delimited list).
@param 	*project		initialized project structure.
@return int				error code (<0 means failure).
**/
int			read_project_crd(Bstring& filename, Bproject* project)
{
	int					err(0);
    char				line[DATLINELENGTH];
	char*				value;
	double				d;
	
	Bstring				field_id = filename.base();
	Bstring				mg_id = field_id;
	Bfield*				field = field_add(&project->field, field_id);
	Bmicrograph*		mg = micrograph_add(&field->mg, mg_id);
	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;
	Bstring				part_prefix;
	Bstring				part_ext;
	
	if ( verbose ) cout << "Reading a CRD file:             " << filename << endl;

	project->comment = "# CRD file: " + filename + "\n";

    /*========================
	 Read CRD-format file
	 ========================*/
	ifstream		fcrd(filename.c_str());
	if ( fcrd.fail() ) {
		cerr << "File " << filename << " does not exist!.  Exit program." << endl;
		bexit(-1);
	}
	
	fcrd.getline(line, DATLINELENGTH);
	
	// General micrograph parameters
	while ( fcrd.getline(line, DATLINELENGTH) && !strstr(line, "$END") ) {
		value = strchr(line, '=') + 1;
		while ( isspace(value[0]) ) value++; 
		if ( strstr(line, "PIC_FILENAME") ) {
			mg->fmg = value;
			mg->fmg = mg->fmg.remove('\n');
			mg->fmg = mg->fmg.remove(',');
			mg->fmg = mg->fmg.remove('\'');
		} else if ( strstr(line, "OUTFILE_PREFIX") ) {
			part_prefix = value;
			part_prefix = part_prefix.remove('\n');
			part_prefix = part_prefix.remove(',');
			part_prefix = part_prefix.remove('\'');
		} else if ( strstr(line, "OUTFILE_TYPE") ) {
			part_ext = value;
			part_ext = part_ext.remove('\n');
			part_ext = part_ext.remove(',');
			part_ext = part_ext.remove('\'');
		} else if ( strstr(line, "ANGSTROMS") ) {
			sscanf(value, "%lf", &d);
			if ( d > 0 ) mg->pixel_size[0] = mg->pixel_size[1] = d;
		} else if ( strstr(line, "RADIUS_PICK") ) {
			sscanf(value, "%lf", &d);
			mg->box_size[0] = mg->box_size[1] = 2*d;
			mg->box_size[2] = 1;
		} else if ( strstr(line, "RADIUS_BAD") ) {
			sscanf(value, "%lf", &mg->bad_radius);
		}
	}
	
	if ( part_prefix.length() )
		mg->fpart = part_prefix + "." + part_ext;
	
	int					numpart(0), numbad(0), pid;
	float				x, y;
    
	while ( fcrd.getline(line, DATLINELENGTH) ) {
		if ( sscanf(line, "%d %f %f", &pid, &x, &y) > 2 )  {
			if ( pid > 0 ) {
				part = particle_add(&part, pid);
				if ( !mg->part ) mg->part = part;
				part->loc[0] = x;
				part->loc[1] = y;
				numpart++;
			} else {
				bad = (Bbadarea *) add_item((char **) &bad, sizeof(Bbadarea));
				if ( !mg->bad ) mg->bad = bad;
				bad->loc[0] = x;
				bad->loc[1] = y;
				numbad++;
			}
		}
    }
	
	fcrd.close();         // close input CRD file
	
	if ( verbose )
		cout << "Number of particles written:    " << numpart << " (" << numbad << ")" << endl << endl;
	
	return err;
}

/**
@brief 	Writing micrograph parameters to a CRD file.
@param 	*filename		file name.
@param 	*project		project structure.
@param 	flags			write flags to select micrographs.
@return int				error code (<0 means failure).
**/
int			write_project_crd(Bstring& filename, Bproject* project, int flags)
{
	int					err(0);
	int 				mg_select(flags & 2);
	int					i, numpart(0), numbad(0);
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;
	Bbadarea*			bad = NULL;
	ofstream			fcrd;
	Bstring				outname = filename;
	Bstring				prefix, type;
	
	if ( verbose ) cout << "Writing a CRD file:             " << filename << endl;
	
	if ( project_count_micrographs(project) < 2 ) mg_select = 1;
	
	for ( i=1, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, i++ )  if ( mg_select < 1 || mg->select > 0 ) {
			if ( mg_select < 1 ) outname = filename.pre_rev('.') + Bstring(i, "_%06d.") + filename.post_rev('.');
			if ( mg->fpart.length() ) {
				prefix = mg->fpart.pre_rev('.');
				type = mg->fpart.post_rev('.');
			}
			fcrd.open(outname.c_str());
			if ( fcrd.fail() )  {                       // make sure input file exists
				cerr << "Unable to create file " << filename << "!.  Exit program." << endl;
				bexit(-1);
			}
			fcrd << " $TRIMNEWPARAMETERS" << endl;
			fcrd << " CRD_VERSION     = 3," << endl;
			fcrd << " PIC_FILENAME    = '" << mg->fmg << "'," << endl;
			fcrd << " PIC_FLIPFLAG    = F," << endl;
			fcrd << " FLIPALLROWSFLAG = F," << endl;
			fcrd << " FLIPALLCOLSFLAG = F," << endl;
			fcrd << " INVERT_INTENSTY = F," << endl;
			fcrd << " ANGSTROMS       = " << mg->pixel_size[0] << "," << endl;
			fcrd << " CONTRAST        = 1.380000" << endl;
			fcrd << " BRIGHTNESS      = 0.120000" << endl;
			fcrd << " OUTFILE_PREFIX  = '" << prefix << "'," << endl;
			fcrd << " OUTFILE_START   = 1," << endl;
			fcrd << " OUTFILE_WIDTH   = 3," << endl;
			fcrd << " OUTFILE_TYPE    = '" << type << "'," << endl;
			fcrd << " OD_CONVERT      = F," << endl;
			fcrd << " RADIUS          = 50," << endl;
			fcrd << " RADIUS_BAD      = " << mg->bad_radius << "," << endl;
			fcrd << " RADIUS_PICK     = " << mg->box_size[0]/2.0 << "," << endl;
			fcrd << " FADE            = 9," << endl;
			fcrd << " EXTRACT_DX      = 501," << endl;
			fcrd << " EXTRACT_DY      = 501," << endl;
			fcrd << " FINAL_DX        = 401," << endl;
			fcrd << " FINAL_DY        = 401," << endl;
			fcrd << " FINAL_MEAN      = 127," << endl;
			fcrd << " FINAL_STDDEV    = 40," << endl;
			fcrd << " SKIPMASKFLAG    = F," << endl;
			fcrd << " SKIPGRADIENTFLAG= T" << endl;
			fcrd << " SELECTEXTENDED  = T" << endl;
			fcrd << " SELECTEXTLENMODE= FIXED" << endl;
			fcrd << " SELECTEXTFLENGTH= 40" << endl;
			fcrd << " $END" << endl;
			for ( part = mg->part; part; part = part->next, numpart++ )
				fcrd << right << setw(6) << part->id << setw(8) << part->loc[0] << setw(8) << part->loc[1] << endl;
			for ( bad = mg->bad; bad; bad = bad->next, numbad++ )
				fcrd << right << setw(6) << -numbad-1 << setw(8) << bad->loc[0] << setw(8) << bad->loc[1] << endl;
			fcrd.close();
		}
	}
	
	if ( verbose )
		cout << "Number of particles written:    " << numpart << " (" << numbad << ")" << endl << endl;
	
	return err;
}



/**
@brief 	Modifies some micrograph parameters.
@param 	*project		project structure.
@param 	flags			bit 0 = omega 90 flag, bit 1 = mag flag.
@return int				error code (<0 means failure).
**/
int			project_modify_parameters(Bproject* project, char flags)
{
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;

	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			for ( part = mg->part; part; part = part->next ) {
				if ( flags & 1 ) part->mag = 1.0/part->mag;
				if ( flags & 2 ) part->view[3] -= M_PI_2;
			}
		}
	}
	
	return 0;
}

