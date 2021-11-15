/**
@file	bradial.cpp
@brief	Polar images and radial averages
@author Bernard Heymann
@date	Created: 20000620
@date	Modified: 20210611
**/

#include "rwimg.h"
#include "ps_plot.h"
#include "Matrix.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bradial [options] input.img output.img",
"---------------------------------------------",
"Calculates radial profiles of 2D and 3D images.",
"The origin within the image will be used by default to calculate a radial profile.",
"To use the radial fit, the sampling (voxel sizes) for the two images",
"	must be within 20 % correct. Successive calibrations with the same input",
"   will not necessarily give the exact same factors on output because ",
"   the simplex algorithm with random initial parameters is used.",
" ",
"Actions:",
"-subtractbackground      Calculate and subtract the background before calculating radial profiles.",
"-radial                  Radial average profile.",
"-elliptical 1.2,20       Elliptical average profile with given ellipticity and angle.",
//"-volume 0.66             Radial volume profile above a threshold.",
"-Radialps 15.5,1,3       Radial power spectrum to a given resolution (in angstrom) ",
"                         with flags and sampling ratio:",
"                         Flag bits: 1=zero origin, 2=average multiple power spectra, 8=logarithm.",
"-polar 156,24,35         Polar image: annuli, phi steps, theta steps.",
"-cylindrical 124,36      Cylindrical image: annuli, phi steps.",
"-shells                  Calculate radial shells projected along the z-axis.",
"-expand 150,150,1        Size to generate a 2D or 3D image from a radial profile (use with -origin).",
"-average                 Average multiple images before writing an output image.",
"-logarithm               Calculate logarithm of radial power spectrum.",
"-reslice -z+xy           Reslice (switch axes) after all other operations.",
"-coverage 1.4            Shell coverage for the given threshold.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (default from image; a single value can be given).",
"-origin 0,-10,30         Set the origin for radial profile (default image origin).",
"-minmax 4,123            Minimum and maximum radii for radial and volume profiles (angstrom).",
"-step 2.68               Step size for radial profiles (default pixel size in angstrom).",
"-wrap                    Flag to wrap around image boundaries (default not).",
" ",
"Input:",
"-fit file.img            Map to fit radial profile of input image to.",
"-Mask mask.tif           Mask to limit calculations.",
" ",
"Output:",
"-Postscript rad.ps       Postscript output file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;    
	Vector3<double>	origin;						// Origin
	int				set_origin(0);				// Flag to set origin
	int 			subtract_background(0);		// Flag to subtract backgrounds
	int 			setradial(0);
	int				avg_flag(0);				// Flag for averaging radial profiles
	int				ps_flags(0);				// Radial power spectrum flags, bits: 1=zero origin, 2=average, 8=log
	Vector3<long>	expand_size;				// Size for expansion
	double			resolution(0); 				// Resolution limit for power spectrum
	double			sampling_ratio(1);			// Sampling ratio for power spectrum
	double			minrad(0), maxrad(0); 		// Radii for radial profiles (angstrom)
	double			ellipticity(1);				// Elliptical averaging parameters
	double			angle(0);
//	double			threshold(-1e37);			// Threshold for volume profiles
	int 			nannuli(-1), ntheta(-1), nphi(-1); // Polar dimensions
	int				set_shells(0);				// Radial shells flag
	int 			setreslice(0);				// Reslice after all other operations
	Bstring			order("xyz");
	int				coverage(0);				// Flag for shell coverage
	double			threshold(-1e37);			// Threshold for shell coverage
	double			step(1);					// Radial step size
	int				wrap(0);					// Flag to wrap
	Bstring			mapfile;					// Map to fit to
	Bstring			maskfile;					// Mask file
	Bstring			outfile;					// Output image file
	Bstring			psfile;						// Postscript output file
	
	int				i, j, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "subtractbackground" )
        	subtract_background = 1;
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
		if ( curropt->tag == "polar" )
			if ( curropt->values(nannuli, nphi, ntheta) < 3 )
				cerr << "-polar: All three dimensions must be specified!" << endl;
		if ( curropt->tag == "cylindrical" )
			if ( curropt->values(nannuli, nphi) < 2 )
				cerr << "-cylindrical: Both dimensions must be specified!" << endl;
		if ( curropt->tag == "shells" )
			set_shells = 1;
		if ( curropt->tag == "radial" ) setradial = 1;
		if ( curropt->tag == "elliptical" ) {
			if ( curropt->values(ellipticity, angle) < 2 )
				cerr << "-elliptical: Both ellipticity and angle must be specified!" << endl;
			else {
				setradial = 1;
				angle *= M_PI/180;
			}
		}
		if ( curropt->tag == "volume" ) {
			if ( ( threshold = curropt->value.real() ) < -1e30 )
				cerr << "-volume: A threshold must be specified!" << endl;
			else
				setradial = 1;
		}
		if ( curropt->tag == "Radialps" )
			if ( curropt->values(resolution, ps_flags, sampling_ratio) < 1 )
				cerr << "-Radialps: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "expand" )
			expand_size = curropt->size();
		if ( curropt->tag == "average" ) { avg_flag = 1; ps_flags |= 2; }
		if ( curropt->tag == "logarithm" ) ps_flags |= 8;
		if ( curropt->tag == "reslice" ) {
			if ( curropt->value.length() < 3 )
				cerr << "-reslice: The order of x, y and z must be specified!" << endl;
			else {
				order = curropt->value;
				setreslice = 1;
			}
		}
		if ( curropt->tag == "coverage" ) {
			coverage = 1;
			threshold = curropt->value.real();
		}
		if ( curropt->tag == "fit" ) {
			mapfile = curropt->filename();
			if ( mapfile.length() ) setradial = 1;
		}
		if ( curropt->tag == "minmax" )
			if ( curropt->values(minrad, maxrad) < 1 )
				cerr << "-minmax: At least a minimum must be specified!" << endl;
		if ( curropt->tag == "step" )
			if ( ( step = curropt->value.real() ) < 0.1 )
				cerr << "-step: A step size must be specified!" << endl;
		if ( curropt->tag == "wrap" )
			wrap = 1;
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	int 		dataflag(0);
	if ( setradial || resolution || optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	if ( !p->data_pointer() ) {
		delete p;
		bexit(0);
	}
	
	// Get the output image file name
	if ( optind < argc ) outfile = argv[optind];
	
//	if ( nudatatype == Unknown_Type )
//		nudatatype = p->data_type();		// Preserve the old type

	if ( subtract_background ) {
		if ( p->background(long(0)) == 0 ) p->calculate_background();
		p->subtract_background();
	}
	
	long			nx, n, k;
	double			f;
	double*			fit_result = NULL;
	Vector3<double> scale(1,1,1);
	Vector3<double>	translate;
	Matrix3			mat(1);
	Bimage*			pmap = NULL;
	Bimage*			prad = NULL;
	Bimage*			pmaprad = NULL;
	Bimage*			pmask = NULL;
	Bstring			title;
	Bplot*			plot = NULL;
	RGB<float>		color;
	
	if ( mapfile.length() )
		pmap = read_img(mapfile, 1, -1);
		
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, -1);
	
	// Sampling is always positive
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	} else origin = p->image->origin();
	
	long			minpix(minrad/(step*p->sampling(0)[0]));
	long			maxpix(maxrad/(step*p->sampling(0)[0]));
	
	if ( setradial ) {
		prad = p;
		if ( p->sizeY() > 1 ) {
			prad = p->radial(minpix, maxpix, step, ellipticity, angle, pmask, wrap);
		} else
			p->change_type(Float);

		nx = prad->sizeX();
		if ( pmap ) {
			n = 3;
			title = "Radial fit";
		} else {
			n = prad->images() + 1;
			title = "Radial plot";
		}
		
//		cout << "Radial size = " << nx << endl;

		plot = new Bplot(1, nx, n);
		plot->title(title);
		plot->page(0).title(title);
		plot->page(0).columns(n);
		for ( i=0; i<n; i++ ) plot->page(0).column(i).number(i);
		plot->page(0).column(0).label("Distance(A)");
		plot->page(0).column(0).axis(1);
		for ( i=1; i<n; i++ ) {
			plot->page(0).column(i).type(2);
			plot->page(0).column(i).label("Average");
			plot->page(0).column(i).axis(3);
		}
		if ( pmap ) {
			plot->page(0).column(1).color(1,0,0);
			plot->page(0).column(2).color(0,0,1);
		}
		plot->page(0).axis(1).min(minrad);
		plot->page(0).axis(1).max(maxrad);
		plot->page(0).axis(3).min(prad->minimum());
		plot->page(0).axis(3).max(prad->maximum());

		if ( pmap ) {
			pmaprad = pmap;
			if ( pmap->sizeY() > 1 ) {
				pmaprad = pmap->radial(minpix, maxpix, step, ellipticity, angle, pmask, wrap);
			} else
				p->change_type(Float);
			fit_result = prad->radial_fit(pmaprad);
			if ( plot ) {
				plot->page(0).axis(3).min(pmaprad->minimum());
				plot->page(0).axis(3).max(pmaprad->maximum());
				for ( i=0; i<nx; i++ ) (*plot)[i] = prad->sampling(0)[0]*i;
				for ( i=0, k=nx; i<nx; i++, k++ ) (*plot)[k] = (*pmaprad)[i];
				for ( i=0; i<nx; i++, k++ ) {
					f = i/fit_result[0];
					j = (int) f;
					f -= j;
					if ( j >= 0 && j < nx-1 )
						(*plot)[k] = fit_result[1]*((1-f)*(*prad)[j] + f*(*prad)[j+1]) + fit_result[2];
				}
			}
			scale[0] = scale[1] = scale[2] = fit_result[0];
			delete[] fit_result;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG bradial: Scale = " << scale[0] << endl;
			if ( outfile.length() > 0 ) {
				translate = pmap->image->origin() - origin;
				p->transform(scale, origin, translate, mat, FILL_BACKGROUND, 0);
				p->sampling(pmap->sampling(0));
			}
		} else {
			if ( verbose & VERB_PROCESS ) {
				cout << "Radial averages:" << endl;
				cout << "Radius";
				for ( j=0; j<prad->images(); j++ ) cout << tab << j+1;
				cout << endl;
			}
//		cout << "Radial size = " << nx << endl;
			for ( i=0; i<nx; i++ ) {
				if ( prad->sampling(0)[0]*i >= minrad ) {
					cout << prad->sampling(0)[0]*i;
					for ( j=0; j<prad->images(); j++ ) cout << tab << (*prad)[nx*j+i];
					cout << endl;
				}
			}
			if ( plot ) {
				for ( i=0; i<nx; i++ ) (*plot)[i] = prad->sampling(0)[0]*i;
				for ( i=0, k=nx; i<nx*prad->images(); i++, k++ ) (*plot)[k] = (*prad)[i];
			}
			delete p;
			p = prad;		// Assign for output
			prad = NULL;
		}
	} else if ( resolution > 0 ) {
		if ( p->compound_type() != TComplex ) p->fft();
		if ( ps_flags & 1 ) p->zero_fourier_origin();
		prad = p->fspace_radial_power(resolution, sampling_ratio);
		delete p;
		p = prad;		// Assign for output
		prad = NULL;
		if ( ps_flags & 2 ) p->average_images();
		if ( ps_flags & 8 ) p->logarithm();
		if ( psfile.length() ) {
			title = "Radial power spectrum";
			n = p->images() + 1;
			nx = p->sizeX();
			plot = new Bplot(1, nx, n);
			plot->title(title);
			plot->page(0).title(title);
			plot->page(0).columns(n);
			for ( i=0; i<n; i++ ) plot->page(0).column(i).number(i);
			plot->page(0).column(0).label("s(1/A)");
			plot->page(0).column(0).axis(1);
			plot->page(0).axis(1).min(0);
			plot->page(0).axis(1).max(1/resolution);
			for ( i=1; i<n; i++ ) {
				plot->page(0).column(i).type(2);
				plot->page(0).column(i).label(Bstring(i, "%d"));
				plot->page(0).column(i).axis(3);
				if ( n > 2 ) color.spectrum(i,1,n-1);
				plot->page(0).column(i).color(color.r(),color.g(),color.b());
			}
			if ( ps_flags & 8 ) {
				plot->page(0).axis(3).flags(2);
				plot->page(0).axis(3).label("log Power");
			} else {
				plot->page(0).axis(3).min(0);
//				plot->page(0).axis(3).max(p->maximum());
//				plot->page(0).axis(3).max((*p)[1]);
				f = (*p)[1];
				for ( i=1; i<nx; ++i )
					if ( f < (*p)[i] ) f = (*p)[i];
				plot->page(0).axis(3).max(f);
				plot->page(0).axis(3).label("Power");
			}
			for ( i=0; i<nx; i++ ) (*plot)[i] = p->sampling(0)[0]*i;
			for ( i=0, k=nx; i<nx*p->images(); i++, k++ ) (*plot)[k] = (*p)[i];
		}
	} else if ( expand_size.volume() > 0 ) {
		pmap = p->radial_to_full(expand_size, origin);
		delete p;
		p = pmap;
		pmap = NULL;
	} else if ( nannuli > -1 ) {
		if ( ntheta > 0 && p->sizeZ() > 1 )
			pmap = p->cartesian_to_spherical(nannuli, nphi, ntheta);
		else
			pmap = p->cartesian_to_cylindrical(nannuli, nphi);
		delete p;
		p = pmap;
		pmap = NULL;
	} else if ( set_shells ) {
		p->radial_shells();
	} else if ( coverage ) {
		p->radial_coverage(threshold, step);
	}

	if ( psfile.length() && plot )
		ps_plot(psfile, plot);

	if ( avg_flag ) p->average_images();


	if ( outfile.length() ) {
		if ( setreslice ) p->reslice(order);
		p->change_type(nudatatype);
		write_img(outfile, p, 0);
	}
	
	delete p;
	order = 0;
	if ( pmap ) delete pmap;
	if ( prad ) delete prad;
	if ( pmaprad ) delete pmaprad;
	if ( pmask ) delete pmask;
	if ( plot ) delete plot;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

