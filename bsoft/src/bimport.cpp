/**
@file	bimport.cpp
@brief	Generate parameter files from a cisTEM hierarchy
@author Bernard Heymann
@date	Created: 20210412
@date	Modified: 20210413

**/

#include "mg_processing.h"
#include "rwmg.h"
#include "mg_ctf.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bproject*	generate_project_from_files(Bstring& framepath, Bstring& mgpath, Bstring& partpath, Bstring& pspath, long mg_per_field, int flags);

// Usage assistence
const char* use[] = {
" ",
"Usage: bimport [options] Assets",
"-------------------------------",
"Generates a micrograph parameter file from directories of files.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1      Sampling (A/pixel, one value sets all).",
"-assets Assets           Generate from cisTEM assets (modified by specifying other paths).",
"-framepath movies,tif    Path for movie frames and extension.",
"-mgpath images,mrc       Path for micrographs and extension.",
"-partpath images,pif     Path for particles and extension.",
"-pspath powerspec,map    Path for power spectra and extension.",
"-numberperfield 3        Number of micrographs per field-of-view (default 1).",
//"-separate                Individual particle files instead of stacks.",
" ",
//"Input:",
//"-particles part.mrcs     Input file with a stack of particles (Relion).",
//" ",
"Output:",
"-output file.star        Output parameter file.",
" ",
NULL
};

int			main(int argc, char** argv)
{
	// Initializing variables
	Vector3<double>	sam;				// A/pixel
	Bstring			assets;				// Path for cisTEM assets
	Bstring			framepath;			// Path for movies
	Bstring			mgpath;				// Path for micrographs
	Bstring			partpath;			// Path for particles
	Bstring			pspath;				// Path for power spectra
	long			mg_per_field(1);	// Micrographs per field
	int				flags(0);			// 1=separate particles
	Bstring			outfile;			// Output parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "assets" ) {
			assets = curropt->value;
			framepath = assets + "/Movies,tif";
			mgpath = assets + "/Images,mrc";
			partpath = assets + "/ParticleStacks,mrc";
			pspath = assets + "/Images/Spectra,mrc";
		}
		if ( curropt->tag == "framepath" )
			framepath = curropt->value;
		if ( curropt->tag == "mgpath" )
			mgpath = curropt->value;
		if ( curropt->tag == "partpath" )
			partpath = curropt->value;
		if ( curropt->tag == "pspath" )
			pspath = curropt->value;
		if ( curropt->tag == "numberperfield" )
			if ( ( mg_per_field = curropt->value.integer() ) < 1 ) {
				cerr << "-numberperfield: The number of micrographs in a field-of-view must be specified!" << endl;
				bexit(-1);
			}
		if ( curropt->tag == "separate" ) flags = 1;
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);

	double			ti = timer_start();
	
	Bproject*		project = generate_project_from_files(framepath, mgpath, partpath, pspath, mg_per_field, flags);
	
	if ( project == NULL )  {
		cerr << "Error: No project generated!" << endl;
		bexit(-1);
	}
	
	if ( sam[0] > 0.1 )
		project_set_mg_pixel_size(project, sam);

	// Write an output parameter file if a name is given
	if ( outfile.length() )
		write_project(outfile, project);
	
	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

Bstring		path_and_extension(Bstring& fn)
{
	Bstring		ext;
	if ( fn.contains(",") ) {
		ext = fn.post(',');
		fn = fn.pre(',');
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG path_and_extension: " << fn << tab << ext << endl;
	
	return ext;
}

long		project_add_files(Bproject* project, vector<string>& files, int ftype, long number)
{
	if ( files.size() < 1 ) return 0;
	
	if ( verbose )
		cout << "Adding files of type " << ftype << endl;
	
	long			i, j, nmg(0), mg_per_field(1);
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bstring			mgid;
	
	if ( number > 1 && ( ftype == 1 || ftype == 2 ) )
		mg_per_field = number;
	
	if ( !project->field ) {
		for ( auto fn: files ) {
			if ( verbose )
				cout << fn << endl;
			mgid = Bstring(fn).base();
			if ( nmg%mg_per_field == 0 ) {
				mg = NULL;
				field = field_add(&field, mgid);
				if ( !project->field ) project->field = field;	// First field
			}
			mg = micrograph_add(&mg, mgid);
			if ( !field->mg ) field->mg = mg;
		}
	}
	
	for ( i=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, ++i ) {
			if ( i >= files.size() ) break;
			if ( verbose )
				cout << i << tab << files[i] << endl;
			Bimage*		p = read_img(files[i], 0, -1);
			if ( !p ) {
				cerr << "File " << files[i] << " not found!" << endl;
			} else {
				if ( ftype == 1 ) {
					mg->fframe = files[i];
					mg->pixel_size = p->image->sampling();
					mg->origin = p->image->origin();
				} else if ( ftype == 2 ) {
					mg->fmg = files[i];
					mg->pixel_size = p->image->sampling();
					mg->origin = p->image->origin();
				} else if ( ftype == 3 ) {
					mg->fpart = files[i];
					mg->box_size = p->size();
					for ( j=0; j<p->images(); ++j ) {
						part = particle_add(&part, j+1);
						if ( !mg->part ) mg->part = part;
						part->pixel_size = p->image[j].sampling();
						part->ori = p->image[j].origin();
						part->view = p->image[j].view();
					}
				} else if ( ftype == 4 ) {
					mg->fps = files[i];
				}
				delete p;
			}
			nmg++;
		}
	}
	
	if ( verbose )
		cout << nmg << " micrographs processed" << endl << endl;

	return nmg;
}

Bproject*	generate_project_from_files(Bstring& framepath, Bstring& mgpath, Bstring& partpath, Bstring& pspath, long mg_per_field, int flags)
{
	Bstring			frameext = path_and_extension(framepath);
	Bstring			mgext = path_and_extension(mgpath);
	Bstring			partext = path_and_extension(partpath);
	Bstring			psext = path_and_extension(pspath);
//cout << "-" << framepath << "-" << frameext << "-" << endl;
	vector<string>	frame_list = file_list(framepath.str(), frameext.str());
//cout << "-" << mgpath << "-" << mgext << "-" << endl;
	vector<string>	mg_list = file_list(mgpath.str(), mgext.str());
	vector<string>	part_list = file_list(partpath.str(), partext.str());
	vector<string>	ps_list = file_list(pspath.str(), psext.str());
	
	if ( verbose ) {
		cout << "Creating file references:" << endl;
		if ( frame_list.size() )
			cout << "Movie references:                " << framepath
				<< " (" << frameext << ", " << frame_list.size() << ")" << endl;
		if ( mg_list.size() )
			cout << "Micrograph references:           " << mgpath
				<< " (" << mgext << ", " << mg_list.size() << ")" << endl;
		if ( part_list.size() )
			cout << "Particle references:             " << partpath
				<< " (" << partext << ", " << part_list.size() << ")" << endl;
		if ( ps_list.size() )
			cout << "Power spectra references:        " << pspath
				<< " (" << psext << ", " << ps_list.size() << ")" << endl;
	}
	
	// Ckecks
	long			err(0);
	if ( frame_list.size() && mg_list.size() && frame_list.size() != mg_list.size() ) {
		cerr << "The number of movie frame files and micrograph files must be the same!" << endl;
		err--;
	}
	if ( mg_list.size() && ps_list.size() && ps_list.size() != mg_list.size() ) {
		cerr << "The number ofpower spectrum files and micrograph files must be the same!" << endl;
		err--;
	}
	
	if ( err ) {
		cerr << "Too many errors! Abort" << endl;
		bexit(-1);
	}
	
	Bproject*		project = new Bproject;
	
	if ( frame_list.size() ) project_add_files(project, frame_list, 1, mg_per_field);
	
	if ( mg_list.size() ) project_add_files(project, mg_list, 2, mg_per_field);
	
	if ( part_list.size() ) project_add_files(project, part_list, 3, mg_per_field);
	
	if ( ps_list.size() ) project_add_files(project, ps_list, 4, mg_per_field);
	
	
	return project;
}

