/**
@file	breconstruct.cpp
@brief	3D reconstruction from single particle images
@author	Bernard Heymann
@date	Created: 20010403
@date	Modified: 20220722
**/

#include "mg_processing.h"
#include "mg_reconstruct.h"
#include "mg_ctf.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "symmetry.h"
#include "Matrix.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern int	thread_limit;	// Thread limit

// Usage assistance
const char* use[] = {
" ",
"Usage: breconstruct [options] input.star [input.star]",
"-----------------------------------------------------",
"Reconstructs 2D and 3D maps from single particle images.",
"Note: -classes and -bootstrap are mutually exclusive.",
" ",
"Actions:",
"-mode 1                  Symmetry mode: 0=during reconstruction (default), 1=after reconstruction,",
"                         2=pick a random symmetry view for each particle.",
"-ewald                   Ewald sphere integration (default central section).",
"-TwoD                    Generate a 2D reconstruction (default 3D).",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-rescale -0.1,5.2        Rescale reconstruction to an average and standard deviation.",
" ",
"Parameters for configuration:",
"-bootstrap               Use particle selection numbers to weigh their contributions.",
"-filaments               Generate one map per filament.",
"-classes 1,3-5           Class selection: can be all, selected, or any combination of numbers (default selected).",
"-fullmap                 Output one full map per class (default if -halfmaps not used).",
"-halfmaps                Output 2 maps from half sets.",
"-threads 2               Total number of threads (default 1, must be even for -halfmaps).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-size 100,120,90         Size of reconstruction (default from images).",
"-sampling 1.5            Reconstruction sampling (A/pixel; default from first particle image).",
"-scale 1.2               Scale or magnification of reconstruction compared to original images (default 1).",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
"-interpolation weighted  Interpolation type: nearest (default), weighted, trilinear.",
"-pad 3                   Image padding factor (default 2).",
"-symmetry D6             Symmetry: Point group identifier",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
" ",
"Output:",
"-reconstruction file.ext Reconstruction file name.",
"-fom file.ext            Figure-of-merit file name.",
"-output file.star        Output parameter file name.",
"-Postscript file.ps      Postscript output filename.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int 			sym_mode(0);				// 0=during, 1=after, 2=random
	int 			flags(0);					// Flags: 1=rescale, 2=2D, 4=bootstrap, 8=ewald
	double			nuavg, nustd(0); 			// Values for rescaling
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	Vector3<long>	map_size;					// Size of reconstruction
	int				interp_type(0);				// Interpolation type
	Vector3<double>	scale(1,1,1);				// Scale or magnification of reconstruction compared to original images
	int				ctf_action(0);				// Default no CTF operation
	double			wiener(0);					// Wiener for CTF correction, flip if 0
	int				pad_factor(2);				// Image padding factor
	double 			resolution(0); 				// Must be set > 0 to limit resolution
	Bsymmetry		sym;						// No symmetry specified
	Bstring			classes("select");			// Class selection string
	long 			nmaps(0); 					// Number of maps per class
	long			imap(0);					// Select one of multiple maps to reconstruct
	int				filaments(0);				// Flag to generate separate filament reconstructions
	Bstring			outfile;					// Output parameter file
	Bstring			reconsfile;					// Output reconstruction file
	Bstring			fomfile;					// Figure-of-merit file
	Bstring			ps_file;					// Output postscript file
	
	thread_limit = 1;		// Default number of threads

	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "mode" ) {
			sym_mode = curropt->value.integer();
			if ( sym_mode < 0 || sym_mode > 2 )
				cerr << "-mode: A symmetry mode of 0, 1 or 2 must be specified!" << endl;
		}
		if ( curropt->tag == "ewald" )
			flags |= 8;
		if ( curropt->tag == "TwoD" )
			flags |= 2;
		if ( curropt->tag == "CTF" )
			ctf_action = curropt->ctf_action();
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			else if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
			else
				flags |= 1;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "size" )
			map_size = curropt->size();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "scale" )
			scale = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "interpolation" ) {
			if ( curropt->value[0] == 'n' ) interp_type = 0;
			if ( curropt->value[0] == 'w' ) interp_type = 1;
			if ( curropt->value[0] == 't' ) interp_type = 2;
		}
		// Configuration
		if ( curropt->tag == "bootstrap" )
			flags |= 4;
		if ( curropt->tag == "classes" )
			classes = curropt->value;
		if ( curropt->tag == "fullmap" ) nmaps += 1;
		if ( curropt->tag == "halfmaps" ) nmaps += 2;
		if ( curropt->tag == "filaments" ) filaments = 1;
		if ( curropt->tag == "wiener" ) {
			if ( ( wiener = curropt->value.real() ) < 0.000001 )
				cerr << "-wiener: A Wiener factor must be specified!" << endl;
			else {
				if ( wiener < 0.01 ) wiener = 0.01;
//				if ( wiener > 1 ) wiener = 1;
			}
		}
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: A positive factor must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "reconstruction" )
			reconsfile = curropt->filename();
		if ( curropt->tag == "fom" )
			fomfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
 		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
   }
	option_kill(option);

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
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( ctf_action > 2 && wiener < 0.01 ) wiener = 0.2;
	
//	project->field->mg->ctf->show();

	if ( verbose ) {
		if ( flags & 2 )
			cout << "2D reciprocal space reconstruction:" << endl;
		else
			cout << "3D reciprocal space reconstruction:" << endl;
	}
	
	long 				npart = project_count_mg_part_selected(project);
	if ( npart < 1 ) {
		cerr << "Error: No particles selected!" << endl;
		bexit(-1);
	}

	if ( thread_limit < 2 ) thread_limit = 2;
	if ( thread_limit > npart/2 ) thread_limit = npart/2;
	
	if ( nmaps < 1 ) nmaps = 1;
	if ( !nmaps%2 ) nmaps = 2;
	if ( nmaps > 3 ) nmaps = 3;
	
	int					nclasses(1);

	if ( filaments ) {
		nclasses = part_set_filament_maps(project);
		nmaps = 1;
	} else {
		nclasses = project_configure_for_reconstruction(project, classes, nmaps, thread_limit);
	}

	int					ntotal = nclasses*thread_limit;

	if ( nmaps & 2 && thread_limit%2 ) thread_limit++;
//	cout << " thread_limit = " << thread_limit << endl;

	Bimage**			pacc = new Bimage*[ntotal];
	Bimage**			prec = new Bimage*[nclasses*nmaps];
	View				ref_view;
	Bstring				filename = reconsfile;
	Bstring				id;
	Breconstruction*	rec = NULL;
	Bparticle*			part = part_find_first(project);

	if ( map_size.volume() <= 0 )
		map_size = project_set_reconstruction_size(project, sam, scale[0], flags & 2);
		
	if ( resolution < 0.1 ) resolution = sam[0];

	if ( map_size.volume() <= 0 ) {
		error_show("Error in breconstruct", __FILE__, __LINE__);
		cerr << "No particles to determine the size from!" << endl << endl;
		return -1;
	}

	long				ft_size = part->mg->box_size[0];
	if ( map_size[0] > ft_size ) ft_size = map_size[0];
	
//	ft_size = part_ft_size(ft_size, scale[0], pad_factor);
	ft_size = part_ft_size(ft_size, 1, pad_factor);
	
	if ( verbose ) {
		cout << "Map size:                       " << map_size << endl;
		cout << "Map sampling:                   " << sam << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Symmetry mode:                  " << sym_mode << endl;
		cout << "Resolution limit:               " << resolution << " A" << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Interpolation type:             " << interp_type << endl;
		cout << "CTF application type:           " << ctf_action << endl;
		if ( flags & 4 ) cout << "Bootstrapping on" << endl;
		if ( flags & 8 ) cout << "Ewald sphere correction on" << endl;
		cout << "Padding factor:                 " << pad_factor << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl << endl;
	}
	
	long				memreq = (ntotal+nmaps)*5*map_size.volume()*sizeof(float);
	memory_check(memreq);

	fft_plan			plan = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	
#ifdef HAVE_GCD
	dispatch_apply(ntotal, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bparticle*	partlist = project_selected_partlist(project, i+1, flags & 4);
		pacc[i] = particle_reconstruct(partlist, sym, sym_mode,
				resolution, scale, sam, map_size, ft_size, plan,
				interp_type, ctf_action, wiener, flags, (i==0));
		particle_kill(partlist);
		if ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) )
			cout << "List " << i+1 << " done: " << timer_report(ti) << endl;
	});
#else
#pragma omp parallel for
	for ( i=0; i<ntotal; i++ ) {
		Bparticle*	partlist = project_selected_partlist(project, i+1, flags & 4);
		pacc[i] = particle_reconstruct(partlist, sym, sym_mode,
				resolution, scale, sam, map_size, ft_size, plan, 
				interp_type, ctf_action, wiener, flags, (i==0));
		particle_kill(partlist);
		if ( verbose & ( VERB_TIME | VERB_PROCESS | VERB_RESULT ) )
			cout << "List " << i+1 << " done: " << timer_report(ti) << endl;
	}
#endif
//	cout << "F0=" << pacc[0]->complex(0).real() << endl;
	
	fft_destroy_plan(plan);

	if ( verbose )
		cout << "Weighing " << nclasses*nmaps << " reconstructions" << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(nclasses*nmaps, dispatch_get_global_queue(0, 0), ^(size_t i){
		prec[i] = img_reconstruction_sum_weigh(pacc, i, nmaps, thread_limit, resolution);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nclasses*nmaps; i++ )
		prec[i] = img_reconstruction_sum_weigh(pacc, i, nmaps, thread_limit, resolution);
#endif
//	cout << "F0=" << prec[0]->complex(0).real() << endl;
	
	for ( i=0; i<ntotal; i++ ) if ( pacc[i] ) delete pacc[i];
	delete[] pacc;

	vector<double> 	fsccut{0.143, 0.3, 0.5, 0.8};	// FSC cutoff values
	vector<double>	dprcut;							// DPR cutoff values
	if ( nmaps > 1 ) {
		for ( i=nmaps-2; i<nclasses*nmaps; i+=nmaps ) {
//		for ( i=nmaps-2; i<nmaps; i+=nmaps ) {
			Bplot*	plot = prec[i]->fsc(prec[i+1], resolution);
			plot->resolution_display(fsccut, dprcut);
			if ( ps_file.length() ) ps_plot(ps_file, plot);
			delete plot;
		}
	}
	
	if ( verbose )
		cout << "Transforming " << nclasses*nmaps << " reconstructions to real space" << endl;

	plan = fft_setup_plan(prec[0]->size(), FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(nclasses*nmaps, dispatch_get_global_queue(0, 0), ^(size_t i){
		prec[i]->fft_back(plan);
		prec[i]->statistics();
		prec[i]->calculate_background();
		if ( sym_mode == 1 ) prec[i]->symmetrize(sym, 1);
		if ( nustd > 0 ) prec[i]->rescale_to_avg_std(nuavg, nustd);
		prec[i]->change_type(nudatatype);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nclasses*nmaps; i++ ) {
		prec[i]->fft_back(plan);
		prec[i]->statistics();
		prec[i]->calculate_background();
		if ( sym_mode == 1 ) prec[i]->symmetrize(sym, 1);
		if ( nustd > 0 ) prec[i]->rescale_to_avg_std(nuavg, nustd);
		prec[i]->change_type(nudatatype);
	}
#endif
	
	fft_destroy_plan(plan);

	// Write an output reconstruction if an output filename is given
	for ( i=0; i<nclasses*nmaps; i++ ) {
		imap = 0;
		if ( nmaps == 2 ) imap = i%2 + 1;
		else if ( nmaps == 3 ) imap = i%3;
		if ( prec[i] && reconsfile.length() ) {
			if ( interp_type == 2 ) prec[i]->filter_sinc();
			if ( fomfile.length() ) {
				filename = fomfile.pre_rev('.');
				if ( nclasses > 1 ) filename += Bstring(i%nmaps + 1, "_%02d");
				if ( imap ) filename += Bstring(imap, "_%02d");
				filename += "." + fomfile.post_rev('.');
				write_img(filename, prec[i]->next, 0);
			}
			prec[i]->symmetry(sym.label().str());
			id = reconsfile.pre_rev('.');
			if ( nclasses > 1 ) id += Bstring(i/nmaps + 1, "_%02d");
			if ( imap ) id += Bstring(imap, "_%02d");
			filename = id + "." + reconsfile.post_rev('.');
			if ( verbose )
				cout << "Writing " << filename << endl;
			write_img(filename, prec[i], 0);
			rec = reconstruction_add(&project->rec, id);
			rec->frec = filename;
			rec->voxel_size = prec[i]->image->sampling();
			rec->origin = prec[i]->image->origin();
			rec->symmetry = sym.label();
		}
		delete prec[i];
	}	
	
	if ( project && outfile.length() ) {
		if ( verbose )
			cout << "Writing " << outfile << endl;
		write_project(outfile, project, 0, 0);
	}
	
	delete[] prec;
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}


