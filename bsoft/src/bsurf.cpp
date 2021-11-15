/**
@file	bsurf.cpp 
@brief	AFM image analysis and reconstruction
@author Bernard Heymann
@date	Created: 19990124
@date 	Modified: 20160604
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bsurf [options] input.tif output.map",
"--------------------------------------------",
"Calculates a 3D surface from an AFM image.",
" ",
"Actions:",
"-expand 10               Surface expansion: value = new z dimension.",
"-threshold 1.6,1         Calculate a 2D height image using this threshold",
"                         and direction: 0=bottom up, 1=top down.",
"-invert                  Invert inside & outside (default not).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-positions 2.5,6.3       Z positions of surface minimum and maximum (angstrom).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-sampling 1.5,1.5,1.5    Dimension scaling (default 1,1,1 A/pixel).",
"-resolution 0.5,2.3      Upper & lower sigma bounds (default 3.236,3.236).",
"-Density 0.74            Density within the surface (default 0.81 Da/A3).",
" ",
"Input:",
"-Stdev sd.map            Standard deviation map (default no map).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam;    			// Pixel size
    long 			nz(1);						// New z-dimension set to preclude expansion
	double			threshold(-1e37);			// Threshold for height image
	int				dir(0);						// Direction for height image
    double 			density(0.74);				// Density within the surface
    double			resmin(1+sqrt(5.0));		// Resolution minimum = 2 phi = 3.236
    double			resmax(1+sqrt(5.0));		// Resolution maximum = 2 phi = 3.236
    double			posmin(0), posmax(0); 		// No min and max for surface expansion
    int 			setinvert(0); 				// Don't invert inside-outside
    Bstring			sdfile;						// Standard deviation file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "expand" )
			if ( curropt->values(nz, dir) < 1 )
				cerr << "-expand: A z dimension must be specified!" << endl;
		if ( curropt->tag == "threshold" )
			if ( ( threshold = curropt->value.real() ) < -1e36 )
				cerr << "-threshold: A threshold must be specified!" << endl;
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "positions" )
			if ( curropt->values(posmin, posmax) < 2 )
				cerr << "-positions: Minimum and maximum positions must be specified!" << endl;
		if ( curropt->tag == "resolution" ) {
			if ( curropt->values(resmin, resmax) < 1 )
				cerr << "-resolution: A resolution must be specified!" << endl;
			if ( resmax < resmin ) resmax = resmin;
		}
		if ( curropt->tag == "Density" ) {
			if ( ( density = curropt->value.real() ) < 0.001 )
				cerr << "-Density: A density must be specified!" << endl;
			if ( density < 1e-10 ) density = 0.74;
		}
		if ( curropt->tag == "Stdev" )
			sdfile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
    // Read the input file if a file name is given
	Bimage*		p = read_img(argv[optind++], 1, -1);
	Bimage*		psd = NULL;
	Bimage*		pnu = NULL;
	
	p->change_type(Float);
	if ( sam.volume() ) p->sampling(sam);
	
	// Expand the image to a 3D surface density
	if ( nz > 1 ) {
		p->show_minimum(posmin);		// Positions in angstrom
		p->show_maximum(posmax);
		
		// Read the standard deviation file if given
		if ( sdfile.length() ) {
			psd = read_img(sdfile, 1, -1);
			
			psd->show_minimum(0.5*resmin);
			psd->show_maximum(0.5*resmax);
		
    		// Convert all data to floating point and check the size
			psd->change_type(Float);
			if ( p->size() != psd->size() ) {
				cerr << "Error: Dimensions of the standard deviation file " << sdfile << " incorrect!" << endl;
				bexit(-1);
			}
		}
		
		// Calculate the surface
		pnu = p->topograph_to_surface(psd, nz, density, resmin);
		delete p;
		p = pnu;
	}
	
	if ( threshold > -1e36 && p->sizeZ() > 1 ) {
		pnu = p->surface_to_topograph(threshold, dir);
		delete p;
		p = pnu;
	}
	
	// Invert if requested
	if ( setinvert ) p->invert();

	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
    // Write an output file if a file name is given
    if ( argc > optind && strspn(argv[optind],"-") != 1 ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete psd;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

