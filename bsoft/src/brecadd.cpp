/**
@file	brecadd.cpp
@brief	Program to average selected reconstructions
@author Bernard Heymann
@date	Created: 20070419 
@date	Modified: 20180427
**/

#include "rwimg.h"
#include "rwmg.h"
#include "mg_reconstruct.h"
#include "mg_processing.h"
#include "mg_particle_select.h"
#include "symmetry.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: brecadd [options] input.star [input.star]",
"------------------------------------------------",
"Averages reconstructions or subvolumes of reconstructions based on",
"the selection in parameter files.",
" ",
"Actions:",
"-all                     Reset selection to all particles before other selections.",
"-group 3022              Select particles belonging to this group.",
"-setnumber 3             Select given selection number and set others to zero.",
"-multiple 2              Generate a number of maps from alternating particle images (default 1).",
"-size 50,50,70           Specify a size to extract new particles from a map.",
"-resolution 15           Specify a resolution to do a reciprocal space summation.",
"-std                     Output standard deviation image rather than variance image.",
"-symmetry C5             Apply point group symmetry after reconstruction.",
"-rescale -0.1,5.2        Rescale output images to average and standard deviation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-average output.img      Output image (default avg.map).",
"-fom fom.img             Output variance/standard deviation image (default fom.map).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int 			reset(0);					// Flag to reset selection
	int				group(0);					// Group to select
	int				set_number(-1);				// Selection number to keep
	int 			nmaps(1); 					// Number of maps to reconstruct
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			nuavg(0), nustd(0);			// Average and standard deviation for rescaling
	int				calcfom(0);					// Flag to calculate variance=1, or std=2
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Vector3<long>	size;						// Size to extract new particles
	double			resolution(0);				// Resolution for frequency space reconstruction
	Bsymmetry		sym;						// Point group
	Bstring			outfile;					// Output parameter file name
	Bstring			avgfile("avg.map");			// Output image file name
	Bstring			fomfile("fom.map");			// Output FOM file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) reset = 1;
		if ( curropt->tag == "group" )
			if ( ( group = curropt->value.integer() ) < 1 )
				cerr << "-group: A group must be specified!" << endl;
		if ( curropt->tag == "setnumber" )
			if ( ( set_number = curropt->value.integer() ) < 1 )
				cerr << "-setnumber: A selection number must be specified!" << endl;
		if ( curropt->tag == "multiple" )
			if ( ( nmaps = curropt->value.integer() ) < 1 )
				cerr << "-multiple: The number of maps must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.1 )
				cerr << "-resolution: The resolution must be specified!" << endl;
		if ( curropt->tag == "std" ) calcfom += 1;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "average" )
        	avgfile = curropt->filename();
		if ( curropt->tag == "output" )
        	outfile = curropt->filename();
		if ( curropt->tag == "fom" ) {
        	fomfile = curropt->filename();
			calcfom += 1;
		}
    }
	option_kill(option);
	
	if ( avgfile.length() < 3 ) {
		cerr << "Error: No output map file name specified!" << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();

	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*			project = read_project(file_list);
	string_kill(file_list);

	if ( !project )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	project->select = 1;
	
	if ( reset ) part_reset_selection(project, 3);
	
	if ( group > 0 ) part_select_group(project, group);

	if ( set_number > 0 ) part_set_selection(project, set_number);

	if ( nmaps > 1 ) part_set_multi_maps(project, set_number, nmaps);

	int					n;
	Bimage*				psum = NULL;
	Bstring				id;
	Bstring				filename;
	Breconstruction*	rec;
	
	for ( n=1; n<=nmaps; n++ ) {
		if ( nmaps > 1 ) set_number = n;
		if ( resolution > 0 )
			psum = project_reconstruct_3D(project, set_number, size, resolution);
		else
			psum = project_reconstruct_3D(project, set_number, calcfom, size);
		if ( !psum ) {
			cerr << "Error: No average calculated!" << endl;
			bexit(-1);
		}
		if ( sam.volume() > 0 ) psum->sampling(sam);
		if ( sym.point() > 101 ) {
			if ( verbose )
				cout << "Applying symmetry " << sym.label() << endl;
			psum->symmetrize(sym, 1);
		}
		if ( avgfile.length() ) {
			if ( nustd > 0 ) psum->rescale_to_avg_std(nuavg, nustd);
			psum->change_type(nudatatype);
			id = avgfile.pre_rev('.');
			if ( nmaps > 1 ) id += Bstring(n, "_%03d");
			filename = id + "." + avgfile.post_rev('.');
//			cout << "writing " << filename << endl;
			write_img(filename, psum, 0);
			if ( verbose )
				cout << "Adding " << id << " to the reconstruction list" << endl;
			rec = reconstruction_add(&project->rec, id);
			rec->frec = filename;
			rec->voxel_size = psum->sampling(0);
			rec->origin = psum->image->origin();
			if ( sym.point() > 101 ) rec->symmetry = sym.label();
			rec->select = 1;
		}
		if ( calcfom && psum->next ) {
			filename = fomfile;
			if ( nmaps > 1 ) filename = fomfile.pre_rev('.') + Bstring(n, "_%03d.") + fomfile.post_rev('.');
			write_img(filename, psum->next, 0);
		}
		delete psum;
	}
	
	if ( outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

