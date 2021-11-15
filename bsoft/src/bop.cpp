/**
@file	bop.cpp
@brief	A program to operate on image pairs
@author Bernard Heymann
@date	Created: 19990219
@date	Modified: 20200108
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bop [options] file1.img file2.img outfile.img",
"----------------------------------------------------",
"Binary operations for 2D and 3D images.",
"Only one operation can be done at a time.",
" ",
"Actions:",
"-slices                  Interpret multiple 2D images as z-slices of a single 3D image",
"-images                  Interpret slices of a single 3D image as 2D images",
"-add 1.4,3.6             Add two images, scaling and shifting the second.",
"-multiply -5.1,0.5       Multiply two images, scaling and shifting the second.",
"-divide -1.8,-2.4        Divide the first image by the second, scaling and shifting the second.",
"-largest 3.1,0.5         Maximum of two images, scaling and shifting the second.",
"-smallest 3.1,0.5        Minimum of two images, scaling and shifting the second.",
"-truncate -0.3,10        Truncate result to min and max.",
"-variance 0.15,-0.5      Scale based on variance ratio, optionally shifting the second.",
"-fit 6                   Linear fit, excluding a percentage of voxels as outliers, output difference.",
"-coefficient             Calculates correlation coefficients between two sets of images.",
"-correlate corr.img      Correlate two images and generate a correlation image.",
"-Histomatch 20           Match histogram of image 2 to image 1, output modified image 2.",
"-blend 7                 Blend from image 1 to image 2, output a new multi-image file.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-radius minr,maxr        Correlate or fit only between min. and max. radius (pixel).",
" ",
"Input:",
"-Mask mask.tif           Mask file to use for limiting comparisons to a defined region.",
" ",
NULL
};

int			main(int argc, char* argv[])
{
	// Initialize all settings
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double> sampling;					// Units for the three axes (A/pixel)
	int				znswitch(0);				// 0=not, 1=n2z, 2=z2n
	int 			setadd(0);
	int 			setmultiply(0);
	int 			setdivide(0);
	int 			setlargest(0);
	int 			setsmallest(0);
	int				setvar(0);
	int 			setlinear(0);
	double 			excl_voxels(0);			// Percentage voxels to excl_voxels from linear fit
	int				setcc(0);
	int				radiusflag(0);
	int 			setminmax(0);
	int 			sethisto(0), bins(0);
	double			scale(1);
	double			shift(0);
//	double			threshold(0);
//	double			fill(0);
	double			min(0), max(0);
	double			minr(0), maxr(0);
	int 			blend_number(0);
	Bstring			corrfile;
	Bstring			maskfile;			// Mask file name
	    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sampling = curropt->scale();
		if ( curropt->tag == "slices" )
			znswitch = 1;
		if ( curropt->tag == "images" )
			znswitch = 2;
		if ( curropt->tag == "add" ) {
        	if ( curropt->values(scale, shift) < 2 )
				cerr << "-add: Both scale and shift must be specified!" << endl;
			else setadd = 1;
		}
		if ( curropt->tag == "multiply" ) {
        	if ( curropt->values(scale, shift) < 2 )
				cerr << "-multiply: Both scale and shift must be specified!" << endl;
			else setmultiply = 1;
		}
		if ( curropt->tag == "divide" ) {
        	if ( curropt->values(scale, shift) < 2 )
				cerr << "-divide: Both scale and shift must be specified!" << endl;
			else setdivide = 1;
		}
		if ( curropt->tag == "largest" ) {
        	if ( curropt->values(scale, shift) < 2 )
				cerr << "-largest: Both scale and shift must be specified!" << endl;
			else setlargest = 1;
		}
		if ( curropt->tag == "smallest" ) {
        	if ( curropt->values(scale, shift) < 2 )
				cerr << "-smallest: Both scale and shift must be specified!" << endl;
			else setsmallest = 1;
		}
		if ( curropt->tag == "truncate" ) {
        	if ( curropt->values(min, max) < 2 )
				cerr << "-truncate: Both min and max must be specified!" << endl;
			else setminmax = 1;
		}
		if ( curropt->tag == "variance" ) {
        	if ( curropt->values(scale, shift) < 1 )
				cerr << "-variance: A variance ratio be specified!" << endl;
			else setvar = setadd = 1;
		}
 		if ( curropt->tag == "fit" ) {
       		if ( ( excl_voxels = curropt->value.real() ) < 0 )
				cerr << "-fit: A percentage of voxels to exclude must be specified!" << endl;
			else
				setlinear = 1;
		}
		if ( curropt->tag == "coefficient" )
        	setcc = 1;
		if ( curropt->tag == "radius" ) {
        	if ( curropt->values(minr, maxr) < 2 )
				cerr << "-radius: Both radii must be specified!" << endl;
			else radiusflag = 1;
		}
		if ( curropt->tag == "Histomatch" ) {
        	if ( ( bins = curropt->value.integer() ) < 1 )
				cerr << "-Histomatch: The number of bins must be specified!" << endl;
			else
				sethisto = 1;
		}
		if ( curropt->tag == "blend" )
        	if ( ( blend_number = curropt->value.integer() ) < 1 )
				cerr << "-blend: A number of images must be specified!" << endl;
 		if ( curropt->tag == "correlate" )
			corrfile = curropt->filename();
 		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
    // Read the input files
	int 		dataflag = 0;
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p1 = read_img(argv[optind++], dataflag, -1);
	if ( p1 == NULL )  {
		cerr << "Error: No first input file read!" << endl;
		bexit(-1);
	}

	if ( znswitch == 1 )
		if ( p1->images_to_slices() < 0 ) bexit(-1);
	
	if ( znswitch == 2 ) 
		if ( p1->slices_to_images() < 0 ) bexit(-1);
	
    p1->change_type(Float);
	if ( p1->standard_deviation() <= 0 ) p1->statistics();
	if ( sampling.volume() > 0 ) p1->sampling(sampling);
	
    Bimage* 	p2 = read_img(argv[optind++], dataflag, -1);
	if ( p2 == NULL )  {
		cerr << "Error: No second input file read!" << endl;
		bexit(-1);
	}

	if ( znswitch == 1 )
		if ( p2->images_to_slices() < 0 ) bexit(-1);
	
	if ( znswitch == 2 ) 
		if ( p2->slices_to_images() < 0 ) bexit(-1);
	
    p2->change_type(Float);
	if ( p2->standard_deviation() <= 0 ) p2->statistics();
	if ( sampling.volume() > 0 ) p2->sampling(sampling);
	
	if ( !p1->data_pointer() ) {
		cerr << "Error: Data missing from image " << p1->file_name() << endl;
		bexit(-1);
	}
	if ( !p2->data_pointer() ) {
		cerr << "Error: Data missing from image " << p2->file_name() << endl;
		bexit(-1);
	}
	
	Bimage*		pmask = NULL;
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, -1);
	
	Bimage* 	pnu = NULL;
	if ( blend_number ) {
		pnu = p1->blend(p2, blend_number);
		delete p1;
		p1 = pnu;
	}

	if ( setvar )
		scale = p1->standard_deviation()/(p2->standard_deviation() * sqrt(scale));
	
//	cout << "scale=" << scale << " shift=" << shift << endl;
	if ( fabs(scale - 1) > 1e-10 ) p2->multiply(scale);
	
	if ( shift ) p2->add(shift);

	if ( setadd )	    		// Add two images
    	p1->add(p2);
    	
	if ( setmultiply )	    	// Multiply two images
    	p1->multiply(p2);
	
	if ( setdivide ) {	    	// Divide the first image by the second image
		if ( p2->images() == 1 ) p1->divide_one(p2);
    	else p1->divide(p2);
	}
		
	if ( setlargest )	    	// Largest of two images
    	p1->largest(p2);
		
	if ( setsmallest )	    	// Smallest of two images
    	p1->smallest(p2);
	
	int				n;
	double			cc(0);
	Vector3<double>	origin(p1->image->origin());
	if ( setlinear ) {	    	// Do a linear fit
		if ( !pmask && radiusflag ) {
			pmask = new Bimage(UCharacter, p1->compound_type(), p1->size(), p1->images());
			pmask->origin(origin);
			pmask->mask_shell(origin, minr, maxr);
//			write_img("mask.pif", pmask);
		}
    	cc = p1->linear_fit(p2, pmask, excl_voxels);
		if ( verbose & VERB_RESULT ) {
			cout << "Image\tR-factor" << endl;
			for ( n=0; n<p1->images(); n++ )
				cout << n+1 << tab << p1->image[n].FOM() << endl;
			cout << endl;
		}
	}
	
	if ( corrfile.length() ) {	    // Calculate a correlation image
		cc = p1->correlate(p2, minr, maxr, pmask, 1);
		p2->change_type(nudatatype);
    	write_img(corrfile, p2, 0);
	} else if ( setcc ) {
		cc = p1->correlate(p2, minr, maxr, pmask, 0);
	}
	
	if ( sethisto ) {
		if ( bins > 10 )
			p1->histomatch(p2, bins);
		else
			cerr << "Too few bins!" << endl;
	}
		
	if ( setminmax )
		p1->truncate_to_min_max(min, max);
	
    // Write an output file if a file name is given
    if ( argc > optind && strspn(argv[optind],"-") != 1 ) {
		p1->change_type(nudatatype);
    	write_img(argv[optind], p1, 0);
	}
	
	delete p1;
	delete p2;
	delete pmask;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit((int)(1000*cc));
}

