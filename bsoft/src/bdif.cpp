/**
@file	bdif.cpp
@brief	Analyze diffraction patterns
@author Bernard Heymann
@date	Created: 20050217
@date	Modified: 20160403
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
"Usage: bdif [options] [input.img] output.img",
"--------------------------------------------",
"Analyzes diffraction patterns.",
" ",
"Actions:",
"-findorigin              Find the origin of the diffraction pattern.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 60.4,-32.2,44    Origin of image (voxels).",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-resolution 4.5,130      Resolution range for correlation (default 0 - 1e6 angstrom).",
"-angle 36.8              Angle for directional averages (default 0).",
" ",
NULL
};

int 	main(int argc, char* argv[])
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> 	sam;				// Units for the three axes (A/pixel)
	int				find_origin(0);
	int				set_origin(0);
	Vector3<double>	origin;			// Origin in file
	double			angle(0);					// Angle for directional averages
	double			hires(0), lores(0);			// Limiting resolution range (hires must be > 0 to be set)

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "findorigin" )
			find_origin = 1;
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
		if ( curropt->tag == "angle" )
        	angle = curropt->value.real() * M_PI/180.0;
    }
	option_kill(option);
        
	double		ti = timer_start();
	
    Bimage* 	p = read_img(argv[optind++], 1, -1);
	Bimage*		pmask = NULL;
	
	if ( !p ) {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	if ( nudatatype == Unknown_Type ) nudatatype = p->data_type(); // Preserve the old type
	
	if ( nudatatype >= Float ) p->change_type(nudatatype);		// Convert for efficiency
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
		
	if ( sam.volume() ) p->sampling(sam);

	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0]*10;
	if ( lores < hires ) lores = hires*100;
	
	double			radius(p->sizeX()/4.0);
	if ( find_origin ) p->find_center(pmask, hires, lores, radius, 0, 1);

	long			i, len(p->sizeX()/2);
	Vector3<double>	end(cos(angle)*len, sin(angle)*len, 0);
	Bimage*			prad = p->radial(0, len, 1);
	Bimage*			pdir1 = p->extract_line(0, p->image->origin(), end, 5);;
	end = Vector3<double>(sin(angle)*len, cos(angle)*len, 0);
	Bimage*			pdir2 = p->extract_line(0, p->image->origin(), end, 5);;
	
	cout << "Radius\ts\tAvg\tDir1\tDir2" << endl;
	for ( i=1; i<prad->sizeX(); i++ )
		cout << i << tab << p->sizeX()*p->sampling(0)[0]/i << tab << (*prad)[i] << tab << (*pdir1)[i] << tab << (*pdir2)[i] << endl;
	cout << endl;
	
	delete prad;
	delete pdir1;
	delete pdir2;
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	delete pmask;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

