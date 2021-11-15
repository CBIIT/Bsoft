/**
@file	bproject.cpp
@brief	Projecting a 3D map and calculating comparison statistics of the projections.
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20200318
**/

#include "rwimg.h"
#include "matrix_util.h"
#include "symmetry.h"
#include "ps_views.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Matrix		img_compare_projections(Bimage* p, double SNR, double hires, double lores, double shift_limit);
int			img_compare_projections_lowpass(Bimage* p, double SNR, double shift_limit);
JSvalue		js_views(View* views);

// Usage assistance
const char* use[] = {
" ",
"Usage: bproject [options] input.img output.img",
"----------------------------------------------",
"Projects a 3D map and calculating comparison statistics of the projections.",
"If the input is a 2D image, it is assumed to be projections.",
" ",
"Actions:",
"-axis y,2                Project down a major axis (x, y, z).",
"                         and a flag for minimum (1) or maximum (2).",
"-kernel 6,2              Reciprocal space projection: kernel size and power.",
"-average 7,5,3           Averaging/smoothing filter before projection: kernel size.",
"-rescale -0.1,5.2        Rescale data to average and standard deviation (after -truncate).",
"-View 0.3,-0.5,0.8,33    View to generate symmetry-related projections.",
"-Tilt -65,72,2.5,55      Tilt series of projections with angle min, max, step and axis angle.",
"-side 15                 Generate side view projections within the given angle from the equator.",
"-compare 0.02,15         Compare projections, with the given SNR imposed and shift limit (default off).",
"-random 24               Generate a number of random projections.",
"-asu                     Change to views to fit into the asymmetric unit.",
"-edgewidth 80            Gaussian width of edge to background (default 0).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-origin 10,-10,20        Origin for rotation (voxels, default 0,0,0).",
"-sampling 2,3.5,1        Sampling (angstrom/voxel, a single value sets all three).",
"-resolution 15.6,200     High and low resolution limits (default 0,1e30).",
"-fill 127                Fill value for resizing: average (default), background, or value.",
"-symmetry D6             Symmetry: Point group identifier.",
"-angles 8,6              Step sizes for theta and phi in the asymmetric unit, one value sets both.",
"-noinplanerotation       No in-plane rotations applied to projections (default applied).",
"-nonormalization         No normalization applied to projections (default applied).",
"-unitcell 10,23,77,90,90,90 Unit cell parameters.",
" ",
"Output:",
"-Plotviews plot.ps       Output postscript file with a plot of projection vectors.",
"-json views.json         Output JSON file with views.",
"-Matrix cc.dat           Projection comparison matrix.",
"-Map cc.map              Projection comparison matrix in an image.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType		nudatatype(Unknown_Type);	// Conversion to new type
	char			axis('n');					// Projection axis
	int				flags(1);					// Flags to modify axis-projection
	int				kernel_width(0);			// Reciprocal space kernel width
	int				kernel_power(2);			// Reciprocal space kernel power
	Vector3<long>	average_kernel;				// Average filter kernel size
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	Vector3<double>	origin;						// Origin
	int				set_origin(0);				// Flag to set origin
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	double			hires(0), lores(1e30);		// High resolution limit
	Bstring			symmetry_string("C1");		// Default: asymmetric or C1 point group
	double			theta_step(0);				// Angular step size for theta
	double			phi_step(0);				// Angular step size for phi
	double			ang_min(0), ang_max(0), ang_step(0), ang_axis(0); // Tilt series
	double			side_ang(-1);				// Side view variation angle
    double   		edge_width(0);    	    	// Gaussian width of edge
	int				setfill(0);
	double			fill(0);
	int 			fill_type(FILL_BACKGROUND);	// Fill type for resizing
	double			compare_snr(0);				// SNR to compare projections
	double			shift_limit(-1);			// Shift limit for projection comparison
	int				lowpass(0);					// Flag to do low-pass range comparison
	int				asu(0);						// Flag to change views to ASU
	int				nviews(0);					// Number of views
	int				spacegroup(0);  	    	// Spacegroup
	UnitCell		uc(0,0,0,M_PI_2,M_PI_2,M_PI_2);
	int				view_flag(1);				// Flag for projection generation
	View			theview;					// View to generate symmetry-related projections
	int				symviews(0);				// Flag to generate symmetry-related projections
	int				ps_flag(0);					// Flag for postscript output
	Bstring			ps_file;
	Bstring			js_file;
	Bstring			matrix_file;
	Bstring			map_file;
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( curropt->tag == "axis" ) {
			if ( ( axis = curropt->value[0] ) < 'a' )
				cerr << "-axis: An axis for projection must be specified!" << endl;
			if ( curropt->value.post(',').length() ) {
				i = 2 * curropt->value.post(',').integer();
				flags |= i;
			}
		}
		if ( curropt->tag == "kernel" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-kernel: At least the kernel size must be specified!" << endl;
		if ( curropt->tag == "average" ) {
			if ( ( i = curropt->values(average_kernel[0],
					average_kernel[1], average_kernel[2]) ) < 1 )
				cerr << "-average: The kernel edge size must be specified!" << endl;
			else
				if ( i < 2 ) average_kernel[2] = average_kernel[1] = average_kernel[0];
		}
		if ( curropt->tag == "rescale" )
        	if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
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
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "fill" ) {
			if ( ! ( fill = curropt->value.real() ) )
				cerr << "-fill: A fill value must be specified!" << endl;
			else
				setfill = 1;
		}
 		if ( curropt->tag == "compare" )
			if ( curropt->values(compare_snr, shift_limit, lowpass) < 1 )
				cerr << "-compare: A SNR value must be specified!" << endl;
 		if ( curropt->tag == "asu" )
       	    asu = 1;
 		if ( curropt->tag == "edgewidth" )
			if ( ( edge_width = curropt->value.real() ) < 0.001 )
				cerr << "-edgewidth: An edge width must be specified!" << endl;
		if ( curropt->tag == "random" )
			if ( ( nviews = curropt->value.integer() ) < 1 )
				cerr << "-random: A number of views must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified!" << endl;
			else {
				theta_step *= M_PI/180.0;
				if ( i < 2 ) phi_step = theta_step;
				else phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "unitcell" )
        	uc = curropt->unit_cell();
		if ( curropt->tag == "noinplanerotation" )
        	view_flag |= 2;
		if ( curropt->tag == "nonormalization" )
        	view_flag |= 4;
		if ( curropt->tag == "View" ) {
			theview = curropt->view();
			symviews = ps_flag = 1;
		}
		if ( curropt->tag == "Tilt" ) {
			if ( curropt->values(ang_min, ang_max, ang_step, ang_axis) < 3 )
				cerr << "-Tilt: Three angles must be specified!" << endl;
			else {
				ang_min *= M_PI/180.0;
				ang_max *= M_PI/180.0;
				ang_step *= M_PI/180.0;
				ang_axis *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "side" ) {
			if ( ( side_ang = curropt->value.real() ) < 0.001 )
				cerr << "-side: One angle must be specified!" << endl;
			else
				side_ang *= M_PI/180.0;
		}
		if ( curropt->tag == "Plotviews" )
			ps_file = curropt->filename();
		if ( curropt->tag == "json" )
			js_file = curropt->filename();
		if ( curropt->tag == "Matrix" )
			matrix_file = curropt->filename();
		if ( curropt->tag == "Map" )
			map_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

#ifdef HAVE_GCD
	fftwf_init_threads();
//	fftwf_plan_with_nthreads(system_processors());
#endif
	if ( verbose )
		cout << "Number of threads:              " << system_processors() << endl;
	
	Bsymmetry 	sym(symmetry_string);
	View*		views = NULL;
	
	if ( axis == 'n' ) {
		if ( nviews ) {
			views = random_views(nviews);
		} else if ( side_ang > -1 ) {
			views = side_views(sym, side_ang, theta_step, phi_step);
		} else if ( ang_step ) {
			views = tilt_views(ang_min, ang_max, ang_step, ang_axis);
		} else if ( symviews ) {
			views = symmetry_get_all_views(sym, theview);
		} else if ( theta_step && phi_step ) {
			views = asymmetric_unit_views(sym, theta_step, phi_step, view_flag);
			ps_flag = 1;
		}
		if ( views ) {
			nviews = count_list((char *) views);
			if ( verbose )
				cout << "Generating " << nviews << " projections" <<endl << endl;
			if ( asu ) change_views_to_asymmetric_unit(sym, views);
			if ( ps_file.length() )
				ps_views(ps_file, symmetry_string, views, ps_flag);
			if ( js_file.length() ) {
				string		fn(js_file.c_str());
				JSvalue		js = js_views(views);
				js.write(fn);
			}
		}
	} else if ( axis == 'z' ) {
		views = new View(0,0,1,0);
	} else if ( axis == 'y' ) {
		views = new View(0,1,0,0);
	} else if ( axis == 'x' ) {
		views = new View(1,0,0,0);
	}
	
	if ( optind >= argc ) {
		if ( views ) kill_list((char *) views, sizeof(View));
		bexit(0);
	}

	FSI_Kernel*		kernel = NULL;
	if ( kernel_width > 1 )
		kernel = new FSI_Kernel(kernel_width, kernel_power);

	// Read image file
	int 		dataflag(0);
	if ( optind < argc - 1 ) dataflag = 1;
	if ( compare_snr ) dataflag = 1;
	
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();		// Preserve the old type
	
	if ( spacegroup ) p->space_group(spacegroup);
	
	if ( sam.volume() ) p->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	p->calculate_background();
	if ( setfill )  p->background(fill);
	
	if ( uc.check() ) p->unit_cell(uc);

	if ( average_kernel[0] > 0 ) p->filter_average(average_kernel);
	
	Bimage* 	proj = p;
	if ( axis != 'n' ) {
		proj = p->project(axis, flags);
	} else if ( p->sizeZ() > 1 && views ) {
		if ( kernel )
			proj = p->project(views, hires, kernel);
		else
			proj = p->project(views, !(view_flag&4));
	}
	
	if ( views ) kill_list((char *) views, sizeof(View));

	Vector3<long>	size((int)(proj->sizeX() - 2*edge_width), 
		(int)(proj->sizeY() - 2*edge_width), (int)(proj->sizeZ() - 2*edge_width));
	origin[0] = origin[1] = origin[2] = edge_width;
    if ( edge_width ) proj->edge(1, size, origin, edge_width, fill_type, fill);

	if ( nustd > 0 ) proj->rescale_to_avg_std(nuavg, nustd);

//	Matrix			cc;
//	Bimage*			pmat;
	if ( compare_snr > 0 ) {
		if ( lowpass ) {
			img_compare_projections_lowpass(proj, compare_snr, shift_limit);
		} else {
			Matrix 	cc = img_compare_projections(proj, compare_snr, hires, lores, shift_limit);
			if ( matrix_file.length() ) cc.write(matrix_file);
			if ( map_file.length() ) {
				Bimage*		pmat = new Bimage(cc, 1);
				write_img(map_file, pmat, 0);
				delete pmat;
			}
		}
	}
	
	if ( optind < argc ) {
		proj->change_type(nudatatype);
		write_img(argv[optind], proj, 0);
	}
	
	if ( p != proj ) delete p;	
	delete proj;
	if ( kernel ) delete kernel;

#ifdef HAVE_GCD
	fftwf_cleanup_threads();
#endif
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/*
	Compares all test images against one reference projection (i)
*/
double	img_compare_one_projection(Bimage* p, Bimage* pref, long i,
				double hires, double lores, double radius, double sigma,
				Matrix* cc, fft_plan planf, fft_plan planb)
{
	long		j;
	double		cc1(0);
	Bimage* 	px;
	
	Bimage*		prefx = pref->extract(i);
	
	for ( j=0; j<p->images(); j++ ) {
		px = p->extract(j);
		px->find_shift(prefx, NULL, hires, lores, radius, sigma, 0, planf, planb, cc1);
		(*cc)[i][j] = cc1;
//		cc[j] = (px->image->origin() - prefx->image->origin()).length();
//		cout << "Shift:    " << (px->image->origin() - prefx->image->origin()).length() << endl;
		delete px;
	}
	
	delete prefx;
	
	return 0;
}

Matrix*		img_compare_projections_driver(Bimage* p, double SNR, double hires, double lores, double shift_limit)
{
	if ( lores <= 0 ) lores = 1e30;
	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	if ( hires > lores ) swap(hires, lores);
	if ( shift_limit < 0 ) shift_limit = 1e30;
	
	p->change_type(Float);
	
	Bimage*			pt = p->copy();
	
	if ( SNR < 101 ) pt->noise_gaussian(0, pt->standard_deviation()/sqrt(SNR));
	
//	write_img("pt.pif", pt);

	long			n(p->images());
	double			sigma(0);
	
	Matrix*			cc = new Matrix(n, n);
	
	fft_plan		planf = p->fft_setup(FFTW_FORWARD, 1);
	fft_plan		planb = p->fft_setup(FFTW_BACKWARD, 1);

#ifdef HAVE_GCD
	dispatch_apply(n, dispatch_get_global_queue(0, 0), ^(size_t i){
		img_compare_one_projection(pt, p, i, hires, lores, shift_limit, sigma, cc, planf, planb);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<n; i++ )
		img_compare_one_projection(pt, p, i, hires, lores, shift_limit, sigma, cc, planf, planb);
#endif

	fft_destroy_plan(planf);
	fft_destroy_plan(planb);

	delete pt;
	
	return cc;
}


/**
@brief 	Compares a set of projections from a 3D density map.
@param 	*p			the 2D projections.
@param 	SNR			SNR to impose on test projections.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	shift_limit	limit on cross-correlation shift.
@return Matrix 		matrix with projection comparisons.
**/
Matrix		img_compare_projections(Bimage* p, double SNR, double hires, double lores, double shift_limit)
{
	Matrix*			ccp = img_compare_projections_driver(p, SNR, hires, lores, shift_limit);
	Matrix			cc = *ccp;
	
	long  			i, j, m, n(p->images()), nw(0), im;
	double			ccmax(0), ccmax1, ccavg(0), ccstd(0);

	if ( verbose & VERB_FULL )
		cout << "Image1\tImage2\tCC" << endl;
	for ( j=0; j<n; ++j ) {
		im = 0;
		ccmax1 = -1;
		for ( i=0; i<n; ++i ) {
			if ( ccmax1 < cc[i][j] ) {
				ccmax1 = cc[i][j];
				im = i;
			}
			if ( i == j ) ccmax += cc[i][j];
			else {
				ccavg += cc[i][j];
				ccstd += cc[i][j]*cc[i][j];
			}
			if ( verbose & VERB_FULL )
				cout << setprecision(5) << i << tab << j << tab << cc[i][j] << endl;
		}
		if ( j != im ) nw++;
		if ( verbose )
			cout << j << tab << im << tab << ccmax1 << endl;
	}
	
	m = n*(n-1);
	ccmax /= n;
	ccavg /= m;
	ccstd = ccstd/m - ccavg * ccavg;
	if ( ccstd > 0 ) ccstd = sqrt(ccstd);
	else ccstd = 0;
	
	if ( verbose ) {
		cout << "Correlating all projections:" << endl;
		cout << "Imposed SNR:                    " << SNR << endl;
		cout << "Resolution limits:              " << hires << " " << lores << endl;
		cout << "Shift limit:                    " << shift_limit << endl;
		cout << "CC maximum:                     " << ccmax << endl;
		cout << "CC average:                     " << ccavg << endl;
		cout << "CC standard deviation:          " << ccstd << endl;
		cout << "Fraction incorrect:             " << nw*1.0/n << endl;
		cout << "Relative maximum sigma units:   " << (ccmax - ccavg)/ccstd << endl;
		cout << "Fisher's z-transform:           " << fishers_z_transform(ccmax) << endl;
		if ( verbose & VERB_PROCESS )
			cout << cc << endl;
	}
	
	return cc;
}

/**
@brief 	Compares a set of projections from a 3D density map.
@param 	*p			the 2D projections.
@param 	SNR			SNR to impose on test projections.
@param 	shift_limit	limit on cross-correlation shift.
@return int 				0.
**/
int			img_compare_projections_lowpass(Bimage* p, double SNR, double shift_limit)
{

	unsigned long  	i, j, k, m, n(p->images());
	double			hires, lores(p->real_size()[0]);
	double			ccmax, ccavg, ccstd, cct;
	Matrix*			ccp;
	Matrix			cc;
	
	if ( verbose )
		cout << "ResLim\tk\tCCmax\tCCavg\tCCstd\tCCth" << endl;
	for ( hires = p->sampling(0)[0]; hires < lores; hires *= 2 ) {
		ccp = img_compare_projections_driver(p, SNR, hires, lores, shift_limit);
		cc = *ccp;
		
		ccmax = ccavg = ccstd = 0;
		for ( m=i=0; i<n; i++ ) {
			for ( j=0; j<n; j++, m++ ) {
				if ( i == j ) ccmax += cc[i][j];
				else {
					ccavg += cc[i][j];
					ccstd += cc[i][j]*cc[i][j];
				}
				if ( verbose & VERB_FULL )
					cout << setprecision(5) << i << tab << j << tab << cc[i][j] << endl;
			}
		}
	
		m -= n;
		ccmax /= n;
		ccavg /= m;
		ccstd = ccstd/m - ccavg * ccavg;
		if ( ccstd > 0 ) ccstd = sqrt(ccstd);
		else ccstd = 0;
		
		k = (unsigned long) (p->real_size()[0]/hires);
		cct = p->sampling(0)[0]/hires;
		cct *= M_PI*cct;
		if ( cct > 1 ) cct = 1;
		cct = sqrt(SNR/(SNR + cct));
		
		delete ccp;
		
		if ( verbose )
			cout << hires << tab << k << tab << ccmax << tab << ccavg << tab << ccstd << tab << cct << endl;
	}
			
	return 0;
}

JSvalue		js_views(View* views)
{
	View*			v;
	vector<double>	vv;
	JSvalue			js(JSarray);
	
	for ( v = views; v; v = v->next ) {
		vv = v->array();
		vv[3] *= 180.0/M_PI;
		js.push_back(JSvalue(vv));
	}
		
	return js;
}

