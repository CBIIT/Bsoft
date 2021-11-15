/**
@file	bbackproj.cpp
@brief	3D reconstruction by backprojection
@author Bernard Heymann
@date	Created: 20010414
@date	Modified: 20150814
**/

#include "rwimg.h"
#include "mg_processing.h"
#include "mg_reconstruct.h"
#include "rwmg.h"
#include "symmetry.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bbackproj [options] input.star [input.star]",
"--------------------------------------------------",
"Reconstructs 3D maps from single particles by backprojection.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-edge                    Smooth map with an oval edge.",
" ",
"Parameters for configuration:",
"-classes 1,3-5           Class selection: can be all, selected, or any combination of numbers (default selected).",
"-fullmap                 Output one full map per class (default if -halfmaps not used).",
"-halfmaps                Output 2 maps from half sets.",
"-threads 2               Threads per class (default 1, must be even for -halfmaps).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-size 120,120,80         Output 3D volume size (default from input file).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-scale 1.2               Scale or magnification of reconstruction compared to original images (default 1).",
"-resolution 3.4          Resolution (angstrom).",
"-symmetry D6             Symmetry: Point group identifier.",
" ",
"Output:",
"-reconstruction file.map Output 3D reconstruction.",
"-output file.star        Output parameter file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<long>	map_size;				// New 3D volume size
	Vector3<double>	sam;    				// Sampling
	double			scale(1);				// Scale or magnification of reconstruction compared to original images
	double			resolution(0); 			// Undefined resolution
	Bsymmetry		sym;					// No symmetry specified
	double			nuavg, nustd(0); 		// Values for rescaling
	int 			smooth_oval(0);			// Flag for oval smoothing
	Bstring			classes("select");		// Class selection string
	long 			nmaps(0); 				// Number of maps per class
	long 			nthreads(1); 			// Number of threads per map
	long			imap(0);				// Select one of multiple maps to reconstruct
	Bstring			outfile;				// Output parameter file
	Bstring			reconsfile;				// Output reconstruction file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "size" )
			map_size = curropt->size();
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.0001 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: Resolution must be specified" << endl;
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "edge" )
			smooth_oval = 1;
		if ( curropt->tag == "classes" )
			classes = curropt->value;
		if ( curropt->tag == "fullmap" ) nmaps += 1;
		if ( curropt->tag == "halfmaps" ) nmaps += 2;
		if ( curropt->tag == "threads" )
			if ( ( nthreads = curropt->value.integer() ) < 1 )
				cerr << "-threads: A number of threads must be specified!" << endl;
 		if ( curropt->tag == "reconstruction" )
			reconsfile = curropt->filename();
		if ( curropt->tag == "output" )
				outfile = curropt->filename();
   }
	option_kill(option);
	
	double		ti = timer_start();
	
	if ( nmaps < 1 ) nmaps = 1;
	if ( !nmaps%2 ) nmaps = 2;
	if ( nmaps > 3 ) nmaps = 3;
	if ( nmaps & 2 && nthreads%2 ) nthreads++;
	
	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*			project = read_project(file_list);
	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( sam[0] > 0 ) {
		if ( sam[0] < 0.1 ) sam = Vector3<double>(1,1,1);
		project_set_mg_pixel_size(project, sam);
	}

	if ( verbose ) {
		cout << "3D backprojection:" << endl;
	}
	
	long				nclasses(1);
	
	nclasses = project_configure_for_reconstruction(project, classes, nmaps, nthreads);
	
	long				ntotal = nclasses*nthreads;
	
	long				i;
	Bimage**			pacc = new Bimage*[ntotal];
	Bimage**			prec = new Bimage*[nclasses*nmaps];
//	View				ref_view;
	Vector3<double>		start;
	Bstring				filename = reconsfile;
	Bstring				id;
	Breconstruction*	rec = NULL;

	Vector3<long>		ms;
	if ( map_size.volume() <= 0 ) {
		ms = project_set_reconstruction_size(project, scale, 0);
		map_size = {ms[0], ms[1], ms[2]};
	}

	if ( map_size.volume() <= 0 ) {
		error_show("Error in bbackproj", __FILE__, __LINE__);
		cerr << "No particles to determine the size from!" << endl << endl;
		return -1;
	}

	if ( verbose ) {
		cout << "Map size:                       " << map_size << endl;
		cout << "Resolution limit:               " << resolution << " A" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << endl;
	}

	fft_plan			planf = fft_setup_plan(map_size[0], map_size[1], 1, FFTW_FORWARD, 1);
	fft_plan			planb = fft_setup_plan(map_size[0], map_size[1], 1, FFTW_BACKWARD, 1);

	// Do the backprojection, weighing each image
#ifdef HAVE_GCD
	dispatch_apply(ntotal, dispatch_get_global_queue(0, 0), ^(size_t i){
		pacc[i] = project_back_projection(project, i+1, map_size, sam,
					scale, resolution, planf, planb);
	});
#else
#pragma omp parallel for
	for ( i=0; i<ntotal; i++ )
		pacc[i] = project_back_projection(project, i+1, map_size, sam,
					scale, resolution, planf, planb);
#endif

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);

	if ( verbose )
		cout << "Combining " << nclasses*nmaps << " reconstructions" << endl;

	// Accumulate partial reconstructions
#ifdef HAVE_GCD
	dispatch_apply(nclasses*nmaps, dispatch_get_global_queue(0, 0), ^(size_t i){
		prec[i] = img_backprojection_accumulate(pacc, i, nmaps, nthreads);
		if ( sym.point() > 101 ) prec[i]->symmetrize(sym, 1);
		if ( smooth_oval ) prec[i]->edge(1, map_size, start, 5, FILL_USER, 0);
		if ( nustd > 0 ) prec[i]->rescale_to_avg_std(nuavg, nustd);
		prec[i]->change_type(nudatatype);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nclasses*nmaps; i++ ) {
		prec[i] = img_backprojection_accumulate(pacc, i, nmaps, nthreads);
		if ( sym.point() > 101 ) prec[i]->symmetrize(sym, 1);
		if ( smooth_oval ) prec[i]->edge(1, map_size, start, 5, FILL_USER, 0);
		if ( nustd > 0 ) prec[i]->rescale_to_avg_std(nuavg, nustd);
		prec[i]->change_type(nudatatype);
	}
#endif

	for ( i=0; i<ntotal; i++ ) delete pacc[i];
	delete[] pacc;

	// Write an output reconstruction if an output filename is given
	for ( i=0; i<nclasses*nmaps; i++ ) {
		imap = 0;
		if ( nmaps == 2 ) imap = i%2 + 1;
		else if ( nmaps == 3 ) imap = i%3;
		if ( prec[i] && reconsfile.length() ) {
			if ( sym.point() > 101 ) prec[i]->symmetry(sym.label().str());
			id = reconsfile.pre_rev('.');
			if ( imap ) id += Bstring(imap, "_%02d");
			filename = reconsfile;
			if (imap ) filename = reconsfile.pre_rev('.') + Bstring(imap, "_%02d.") + reconsfile.post_rev('.');
			write_img(filename, prec[i], 0);
			rec = reconstruction_add(&project->rec, id);
			rec->frec = filename;
			rec->voxel_size = project->field->mg->pixel_size/scale;
			rec->origin = prec[i]->image->origin();
			if ( sym.point() > 101 ) rec->symmetry = sym.label();
		}
		delete prec[i];
	}	
	
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	delete[] prec;
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

