/**
@file	btomrec.cpp
@brief	Disk-based 3D reconstruction for a tomography series
@author	Bernard Heymann
@date	Created: 20031205
@date	Modified: 20200505
**/

#include "mg_tomo_rec.h"
#include "mg_align.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: btomrec [options] input.star [input2.star]",
"-------------------------------------------------",
"Reconstructs a tomogram from an aligned tilt series in Fourier space.",
"Requires large amounts of disk space.",
" ",
"Actions:",
"-removemarkers 14        Mask out markers with this radius (pixels) before reconstruction.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
"-slab 25,40              Reconstruct a slab of these slices (default all).",
"-transform 2D            Back transform reconstruction: none (default), full, slices.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling of input images (A/pixel; a single value can be given).",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
"-scale 1.2               Scale or magnification of reconstruction compared to original images (default 1).",
"-size 100,120,90         Size of reconstruction (default from images).",
"-interpolation weighted  Interpolation type: nearest (default), weighted, trilinear.",
"-pad 3                   Image padding factor (default 2).",
"-edge 15                 An edge smoothing width (default 0).",
"-fill 127                Fill value for erasing/painting markers (default average).",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
" ",
"Output:",
"-reconstruction file.ext Reconstruction file name.",
"-output file.star        Output STAR file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	double			marker_radius(0);			// Radius to mask out markers
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double			nuavg(0), nustd(0); 		// Values for rescaling
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	Vector3<long>	size;						// Size of reconstruction
	int				slab_start(0), slab_end(-1);// Slab definition
	int				interp_type(0);				// Interpolation type
	int				pad_factor(2);				// Image padding factor
	double			edge_width(0);				// An edge smoothing width
	double			scale(1);					// Scale or magnification of reconstruction compared to original images
	double 			resolution(0); 				// Must be set > 0 to limit resolution
	int				transform(0);				// No back transform
	int 			fill_type(FILL_AVERAGE);	// Type of fill value
	double			fill(0);					// Fill value for new areas
	int				ctf_action(0);				// Default no CTF operation
	double			wiener(0.2);				// Wiener for CTF correction
	Bstring			reconsfile;					// Output reconstruction file
	Bstring			outfile;					// Output STAR file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "removemarkers" )
			if ( ( marker_radius = curropt->value.real() ) < 1 )
				cerr << "-removemarkers: A marker radius must be specified!" << endl;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
		if ( curropt->tag == "slab" )
			if ( curropt->values(slab_start, slab_end) < 1 )
				cerr << "-slab: A slab size must be specified!" << endl;
		if ( curropt->tag == "transform" ) {
			if ( curropt->value[0] == 'n' ) transform = 0;
			if ( curropt->value[0] == 'f' ) transform = 1;
			if ( curropt->value[0] == 's' ) transform = 2;
		}
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "scale" )
			if ( ( scale = curropt->value.real() ) < 0.001 )
				cerr << "-scale: A scale must be specified!" << endl;
		if ( curropt->tag == "interpolation" ) {
			interp_type = curropt->value.integer();
			if ( curropt->value[0] == 'n' ) interp_type = 0;
			if ( curropt->value[0] == 'w' ) interp_type = 1;
			if ( curropt->value[0] == 't' ) interp_type = 2;
		}
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: An integer factor must be specified!" << endl;
		if ( curropt->tag == "edge" )
			if ( ( edge_width = curropt->value.real() ) < 0 )
				cerr << "-edge: An edge smoothing width must be specified!" << endl;
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "CTF" )
			ctf_action = curropt->ctf_action();
		if ( curropt->tag == "wiener" ) {
			if ( ( wiener = curropt->value.real() ) < 0.000001 )
				cerr << "-wiener: A Wiener factor must be specified!" << endl;
			else {
				if ( wiener < 0.01 ) wiener = 0.01;
//				if ( wiener > 1 ) wiener = 1;
			}
		}
		if ( curropt->tag == "reconstruction" )
			reconsfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

#ifdef HAVE_GCD
	fftwf_init_threads();
	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	

	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( sam[0] > 0 ) {
		if ( sam[0] < 0.1 ) sam[0] = 1;
		project_set_mg_pixel_size(project, sam);
	}

	Bimage*			prec = NULL;

	if ( size[2] < 2 ) {
		prec = mg_tomo_reconstruct2D(project, 10, size, scale, resolution);
	} else if ( slab_end > 0 ) {
		prec = project_fourier_reconstruction_slab(project, resolution,
					scale, size, slab_start, slab_end, marker_radius,
					fill_type, fill, ctf_action, wiener);
	} else {
		prec = project_tomo_reconstruct(project, resolution,
					scale, size, interp_type, pad_factor, 
					edge_width, marker_radius,
					fill_type, fill, ctf_action, wiener);
		prec->change_transform_size(size);
	}

	Bstring			id("1");
	Bmicrograph*	mg = field_find_zero_tilt_mg(project->field);
	Vector3<double>	ori(mg->origin);
	Vector3<double>	trans(prec->image->origin() - ori*scale);
	if ( verbose )
		cout << "Translating markers:      " << trans << endl;
	
	if ( transform == 1 ) {
		if ( verbose )
			cout << "Full Fourier transform" << endl;
		if ( prec->compound_type() == TComplex ) prec->fft_back();
	} else if ( transform == 2 ) {
		if ( verbose )
			cout << "Fourier transforming slices" << endl;
		img_backtransform_slices(prec);
	}
	
	if ( prec && reconsfile.length() ) {
		if ( nustd > 0 ) prec->rescale_to_avg_std(nuavg, nustd);
		prec->change_type(nudatatype);
		if ( verbose )
			cout << "Writing " << reconsfile << endl;
		write_img(reconsfile, prec, 0);
	}

	if ( project && outfile.length() ) {
		if ( !project->rec )
			reconstruction_add(&project->rec, id);
		project->rec->frec = reconsfile;
		project->rec->voxel_size = project->field->mg->pixel_size/scale;
		project->rec->voxel_size[2] = project->rec->voxel_size[0];
		project->rec->origin = prec->image->origin();
		project->rec->scale = Vector3<double>(scale, scale, scale);
		project->rec->box_size = project->field->mg->box_size * scale;
		project->rec->filament_width = project->field->mg->filament_width * scale;
		project->rec->fil_node_radius = project->field->mg->fil_node_radius * scale;
		project->rec->bad_radius = project->field->mg->bad_radius * scale;
		project->rec->sf_radius = project->field->mg->sf_radius * scale;
		project->rec->mark_radius = project->field->mg->mark_radius * scale;
//		marker_transform(project->rec->mark, tf);
		marker_scale(project->rec->mark, project->rec->scale);
		marker_shift(project->rec->mark, trans);
		if ( verbose )
			cout << "Writing " << outfile << endl;
		write_project(outfile, project, 0, 0);
	}

	delete prec;
	project_kill(project);

#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

