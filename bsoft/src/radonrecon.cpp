/**
@file	radonrecon.cpp
@brief	Program to do radon transforms.
@author Salvatore Lanzavecchia, Francesca Cantele and Pier Luigi Bellon
         Dip. Chimica Strutturale e Stereochimica Inorganica, Via Venezian 21, 20133 Milano, Italy
@author	Bernard Heymann
         Rm 1515, 50 South Dr., NIH, Bethesda, MD, 20892, USA

@date	Created: 2003 07 04
@date	Modified: 20160728 (BH)
**/

#include "mg_processing.h"
#include "rwmg.h"
#include "rwimg.h"
#include "img_radon.h"
#include "symmetry.h"
#include "utilities.h"
#include "options.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

// Usage assistence
const char *use[] = {
  " ",
  "Usage: radonrecon [options] input.star [input.star]",
  "---------------------------------------------------",
  "Program to reconstruct a quarter (PI*PI) radon transform from projections.",
  "Projections must be square.",
  " ",
  "Parameters:",
  "-verbose 7               Verbosity of output",
  "-datatype u              Force writing of a new data type (default output floating point)",
  "-origin 32,32            Projections origin (default center of projections)",
  "-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
  "-kernel 15,4             Kernel size and power (default 11,2)",
  "-threshold 19.4          Threshold for accepting a particle image (default 0)",
  "-symmetry D5             Number of symmetry equivalents (default 1)",
  " ",
  "Output:",
  "-output file.star        Output parameter file.",
  "-reconstruction rec.map  Reconstruction file name.",
  "-mask mask.map           Binary mask for POCS filter.",
  " ",
  NULL
};

int
main (int argc, char **argv)
{
	// Initialize variables
	DataType 	newdatatype = Unknown_Type;	// Conversion to new data type for the reconstruction
	int 		rec_size(0);				// Reconstruction size - default from particle image size
	Vector3<double> origin;					// Origin of reconstruction
	Vector3<double> sampling;				// Units for the three axes (A/pixel)
	int 		n_theta(0);					// Number of angular samples per 2*PI
	int 		nkernel(11);				// Kernel width
	int 		kernel_power(2);			// Kernel power
	double 		threshold(0);				// Threshold to accept particle images
	Bsymmetry	sym;						// Point group
	Bstring		outfile;					// File name for parameter output
	Bstring		reconsfile;					// File name for reconstruction
	Bstring		maskfile;					// File name for mask file output

	int sizeTable = 8;

	int optind;
	Boption *option = get_option_list (use, argc, argv, optind);
	Boption *curropt;
	for (curropt = option; curropt; curropt = curropt->next)  {
		if ( curropt->tag == "datatype" )
			newdatatype = curropt->datatype();
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "kernel" )
			if ( curropt->values(nkernel, kernel_power) < 1)
				cerr << "-kernel: At least the kernel size must be specified!" << endl;
		if ( curropt->tag == "threshold" )
			if ( ( threshold = curropt->value.real() ) < 1e-30)
				cerr << "-threshold: The threshold must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "reconstruction" )
			reconsfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
    }
	option_kill (option);

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


	if ( sampling[0] > 0 ) {
		if ( sampling[0] < 0.1 ) sampling = Vector3<double>(1,1,1);
		project_set_mg_pixel_size (project, sampling);
    }

	Bimage *prec = img_radon_reconstruction (project, sym, maskfile, rec_size, n_theta, sizeTable, threshold, origin, nkernel, kernel_power);

	if (prec && reconsfile.length()) {
		prec->statistics();
		prec->change_type(newdatatype);
		write_img (reconsfile, prec, 0);
    }

	delete prec;

	if (project && outfile.length())
		write_project (outfile, project, 0, 0);

	project_kill (project);

	return 0;
}
