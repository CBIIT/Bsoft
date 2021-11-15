/**
@file	bscale.cpp 
@brief	Program to find the scaling of map with respect to a reference
@author Bernard Heymann
@date	Created: 20160317
@date 	Modified: 20170106
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
"Usage: bscale [options] input.pif output.map",
"--------------------------------------------",
"Scales a map with respect to a reference.",
" ",
"Actions:",
"-bin 3                   Bin the images before comparison.",
"-average 7               Averaging/smoothing filter: kernel size.",
"-symmetry C5             Find point group symmetry equivalent.",
"-fit                     Fit the map to the reference.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-sampling 1.5,1.5,1.5    Dimension scaling (default 1,1,1 A/pixel).",
"-range 0.95,1.08,0.01    Scaling range and step size (default 0.99,1.01,0.01).",
" ",
"Input:",
"-reference temp.map      Reference map to scale to.",
"-mask mask.mrc           Mask (must be the same size as the reference).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	int				bin(0);						// Binning
	Bsymmetry		sym;					// Point group
	int				fit(0);
	double			nuavg(0), nustd(0); 		// Values for rescaling
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam;    			// Pixel size
	long			average_kernel(0);			// Average filter kernel size
    double 			smin(0.99), smax(1.01);		// Scaling search range
	double			step(0);					// Scaling search step
    Bstring			reffile;					// Reference file name
    Bstring			maskfile;					// Mask file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "fit" ) fit = 1;
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
			if ( nustd <= 0 )
				cerr << "-rescale: A positive standard deviation must be specified!" << endl;
		}
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
		if ( curropt->tag == "average" )
			if ( ( average_kernel = curropt->value.integer() ) < 1 )
				cerr << "-average: The kernel edge size must be specified!" << endl;
		if ( curropt->tag == "range" )
			if ( curropt->values(smin, smax, step) < 2 )
				cerr << "-range: Minimum and maximum scales must be specified!" << endl;
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
    // Read the input file if a file name is given
	Bimage*			p = read_img(argv[optind++], 1, -1);
	Bimage*			pref = NULL;
	Bimage*			pmask = NULL;
 	Bimage*			pscale = NULL;
	Matrix3			mat(1);
       
	p->change_type(Float);
	if ( sam.volume() ) p->sampling(sam);

	if ( set_origin > 1 ) p->origin(p->size()/2);
	else if ( set_origin ) p->origin(origin);
	
	if ( bin > 1 ) {
		pscale = p->bin_around_origin(bin);
		delete p;
		p = pscale;
		pscale = NULL;
	}
	
	if ( average_kernel ) p->filter_average(average_kernel);
	
	if ( reffile.length() ) {
		pref = read_img(reffile, 1, -1);

		if ( average_kernel ) pref->filter_average(average_kernel);
	
		if ( maskfile.length() )
			pmask = read_img(maskfile, 1, -1);

		if ( sym.point() > 101 && sym.point() != 432 ) {
			mat = p->symmetry_equivalent(pref, pmask, sym);
			if ( sym.point() < 200 ) p->rotate(-p->image->origin()+p->size()/2, mat);
			else p->rotate(mat);
//		} else {
//			mat = img_align(p, pref, pmask);
//			p->rotate(-p->image->origin()+p->size()/2, mat);
		}
	
		if ( step ) pscale = p->scale_to_reference(pref, pmask, smin, smax, step);
		else if ( fit ) pscale = p->scale_to_reference(pref, pmask);
		else pscale = p->scale_to_same_size(pref);
	}
	
	if ( pscale ) {
		delete p;
		p = pscale;
	}

	delete pref;
	delete pmask;

	if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
	
    // Write an output file if a file name is given
    if ( argc > optind ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}
	
	delete p;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}
/*
Matrix3		img_align(Bimage* p, Bimage* pref, Bimage* pmask)
{
	Matrix3			mat(1);
	double			accuracy(2/pref->sizeX());
	
	while ( dang >= accuracy ) {
		ccbest = -1;
		cout << "∆ang = " << dang*180.0/M_PI << tab << "res = " << hires << endl;
		cout << "Shift\t\t\tAngle\tCC" << endl;
		for ( angle = angmin; angle <= angmax; angle += dang ) {
			mat = Matrix3(axis, angle);
			pt = transform(pref->size(), vscale, ori, translate, mat, FILL_BACKGROUND);
//			cc = pt->correlate(pref, 0, pref->sizeX()/2, pmask);
			shift = pt->find_shift(pref, NULL, hires, 0, 10, 0, 1, cc);
			if ( ccbest < cc ) {
				ccbest = cc;
				best_angle = angle;
				best_shift = shift;
//				if ( verbose & VERB_PROCESS )
//					cout << "Better CC:                      " << cc << " (" << angle*180.0/M_PI << ")" << endl;
			}
			delete pt;
			if ( verbose )
				cout << shift << tab << angle*180.0/M_PI << tab << cc << endl;
		}
		dang *= 0.75;
		angmin = best_angle - 2*dang;
		angmax = best_angle + 2.1*dang;
		hires = dang*pref->real_size()[0];
	}	
	
	return mat;
}
*/
