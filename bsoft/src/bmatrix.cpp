/**
@file	bmatrix.cpp
@brief	Program to process matrices.
@author Bernard Heymann
@date	Created: 20010723
@date	Modified: 20200916
**/

#include "cluster.h"
#include "Matrix.h"
#include "matrix_linear.h"
#include "matrix_util.h"
#include "math_util.h"
#include "Bimage.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmatrix [options] matrix_in.dat matrix_out.dat",
"-----------------------------------------------------",
"Processes matrices.",
" ",
"Actions:",
"-random 4,7              Create a random matrix of the given size (rows & columns).",
"-negate                  Change the sign of the matrix.",
"-invert                  Invert by LU decomposition.",
"-svd                     Invert by singular value decomposition.",
"-window 5                Window over which to do permutations for linear sorting.",
"-preference -10          Homogenous preference for clustering.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-ln_1_R                  Use the -ln(1-R) conversion of the input matrix.",
"-lambda 0.6              Damping factor for clustering (default 0.5).",
" ",
"Parameters for image output:",
"-datatype b              Force writing of a new data type.",
"-scale 3                 Scale the image (must be an integer).",
"-color                   Convert to a heat map.",
" ",
"Input:",
"-vector file.dat         Input a vector to multiply with the matrix.",
" ",
"Output:",
"-outvector file.dat      Output the vector.",
"-outimage file.img       Output the matrix as an image.",
" ",
NULL
};

int 		main(int argc, char** argv)
{
	// Initialize variables
	long			m(0), n(0);				// Dimensions for random matrix
	int 			negate(0);				// Flag to negate
	int 			invert(0);				// Flag for LU decomp inversion
	int				svd(0);					// Flag for singular value decomposition
	int 			window(0);				// Window for linear sorting
	double			pref(-1e37);			// Preference value for clustering
	double			lambda(0.5);			// Damping factor for clustering
	int 			log_1_R(0);
	Bstring			vector_file;			// Vector input file name
	Bstring			vecout_file;			// Vector output file name
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	long			img_scale(0);			// Image enlargement
	int				img_color(0);			// Flag to color the image
	Bstring			img_file;				// Output image file name
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "random" ) {
			if ( ( i = curropt->values(m, n) ) < 1 )
				cerr << "-random: A matrix size must be specified!" << endl;
			else
				if ( i < 2 ) n = m;
		}
		if ( curropt->tag == "negate" )
			negate = 1;
		if ( curropt->tag == "invert" )
			invert = 1;
		if ( curropt->tag == "svd" )
			svd = 1;
		if ( curropt->tag == "window" ) {
			if ( ( window = curropt->value.integer() ) < 1 )
				cerr << "-window: A window size must be specified!" << endl;
			else
				if ( window < 2 ) window = 2;
		}
		if ( curropt->tag == "preference" )
			if ( ( pref = curropt->value.real() ) < -1e36 )
				cerr << "-preference: A preference value must be specified!" << endl;
		if ( curropt->tag == "ln_1_R" )
			log_1_R = 1;
		if ( curropt->tag == "lambda" )
			if ( ( lambda = curropt->value.real() ) < 0.001 )
				cerr << "-lambda: A damping factor must be specified!" << endl;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "scale" )
			if ( ( img_scale = curropt->value.integer() ) < 2 )
				cerr << "-scale: An integer must be specified!" << endl;
		if ( curropt->tag == "color" ) img_color = 1;
		if ( curropt->tag == "vector" )
			vector_file = curropt->filename();
		if ( curropt->tag == "outvector" )
			vecout_file = curropt->filename();
		if ( curropt->tag == "outimage" )
			img_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	Matrix		matrix, vector;
	Bstring		filename;

	if ( n > 0 ) {
		matrix = Matrix(m, n);
		matrix.randomize();
	} else {
		filename = argv[optind++];
		matrix = Matrix(filename);
	}
	
	if ( verbose )
		cout << "Matrix size:" << tab << matrix.rows() << " x " << matrix.columns() << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Input matrix:" << endl << matrix << endl;
	
	if ( vector_file.length() )
		vector = Matrix(vector_file);
	
	if ( negate )
		matrix = -matrix;
	
	double		det;
	if ( invert ) {
		det = matrix.LU_decomposition();
		if ( verbose & VERB_PROCESS ) {
			cout << "Determinant:                        " << det << endl;
			cout << "Inverse matrix:" << endl << matrix << endl;
		}
	}
	
	if ( svd ) {
		matrix.singular_value_decomposition();
		if ( verbose & VERB_PROCESS ) {
			cout << "Inverse matrix:" << endl << matrix << endl;
		}
	}

	if ( log_1_R ) matrix_log_1_R(matrix);
	
	int*		order = NULL;
	if ( window ) {
		order = matrix_find_linear_sequence(matrix, window);
//		matrix_rearrange(matrix, order);
		delete[] order;
	}
	
	long		ncluster(0);
	if ( pref > -1e36 ) {
		for ( i=0; i<matrix.rows(); i++ )
			matrix[i][i] = pref;
		affin_prop_clustering(matrix, 500, 50, lambda, ncluster);
	}
	
	if ( vector.rows() )
		vector = matrix * vector;

	if ( verbose && vector.rows() ) {
		cout << "Output vector:" << endl << vector << endl;
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Output matrix:" << endl << matrix << endl;
	}
	
	if ( vecout_file.length() )
		vector.write(vecout_file);
	
	if ( optind < argc ) {
		filename = argv[optind];
		matrix.write(filename);
	}
	
//	Bimage*		pmat = NULL;
//	Bimage*		pt = NULL;
	if ( img_file.length() ) {
//		pmat = img_from_matrix(matrix, img_scale);
		Bimage*		pmat = new Bimage(matrix, img_scale);
		if ( img_color ) {
			pmat->change_type(UCharacter);
			Bimage*		pt = pmat->color_spectrum(pmat->minimum(), pmat->maximum());
			delete pt;
		}
		if ( nudatatype > Unknown_Type ) pmat->change_type(nudatatype);
		write_img(img_file, pmat, 0);
		delete pmat;
	}
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}



