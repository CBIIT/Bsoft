/**
@file	bhelix.cpp
@brief	Dealing with images with helical symmetry.
@author Bernard Heymann
@date	Created: 20021127
@date	Modified: 20210322
**/

#include "rwimg.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bhelix [options] input.img output.img",
"--------------------------------------------",
"Symmetrizes and projects images with helical symmetry.",
" ",
"Actions:",
"-normalize               Normalize symmetrized output, use for maps close to symmetry (default not).",
"-background              Calculate new background values.",
"-project 20              Generate projections with this angular step size (over 360 degrees).",
"-kernel 6,2              Reciprocal space projection: kernel size and power.",
"-elliptical 37,5.5       Elliptical distortion: major axis angle and shift (pixels).",
"-rescale -0.1,5.2        Rescale data to average and standard deviation.",
" ",
"Actions for symmetrizing:",
"-helix 25.3,67.1,2,5     Symmetrize: Helical rise and angle, dyad axis (1/2) and cyclic symmetry.",
"                         (default 10 angstrom, 45 degrees, no dyad axis, C1).",
"-rule 237.5,-2,7         Symmetrize by selection rule: Helical repeat (angstrom), turns and subunits per repeat.",
"-super 25.3,67.1,22,15.5 Superhelix: Helical rise and angle, and offset in xy plane.",
"-seam 1.5                Symmetrize with the given subunit shift along the seam.",
" ",
"Actions for determining helical parameters:",
"-findrs 5,8,0.1,10,20,1  Find helical parameters in real space:",
"                         rise start, end and increment, rotation angle start, end and increment.",
"-findcc 10,22.4,1.5      Find helical parameters by cross-correlation:",
"                         rotation angle start, end and increment.",
"-segment 10,22.4,1.5,18  Correlate helical segments: rotation angle start,",
"                         end, increment and segment thickness.",
" ",
"Actions for constructing a filament from projections:",
"-views seq               Sets the view sequential or random.",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 110,50,44        Origin for helix (default 0,0,0 or from image).",
"-zlimits 8,47            Range of slices along the helical axis to use (default 0,inf).",
"-radius 22               Radial limit for symmetrization and finding parameters (default from image).",
"-resolution 15.6,200     Resolution limits for projection and cross-correlation (default 0).",
"-bin 2                   Binning for finding parameters (default 1).",
" ",
"Output:",
"-Postscript file.ps      Postscript output file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			set_sampling(0);
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	double			helix_rise(0);				// Rise per asymmetric unit
	double			helix_angle(0);				// Rotation angle per asymmetric unit
	int				dyad_axis(1);				// No dyad axis
	int				cyclic(1);					// Cyclic symmetry
	double			helix_repeat(0);			// Helical repeat
	int				turns(0), subunits(0);		// Turns and subunits per repeat
	double			seam_shift(0);				// Shift along seam (0=no seam)
	Vector3<double>	offset;						// Superhelix offset
	double			project_angle(0);			// Projection angular step size
	int				kernel_width(0);			// Reciprocal space kernel width
	int				kernel_power(2);			// Reciprocal space kernel power
	double			rise_start(0), rise_end(0);	// Rise range for finding parameters
	double			rise_step(0);				// Rise increment for finding parameters
	double			angle_start(0), angle_end(0);	// Angle range for finding parameters
	double			angle_step(0);				// Angle increment for finding parameters
	int				thickness(0);				// Segment thickness for correlation
	int				find_bin(1);				// Binning for searching
	double			hires(0), lores(0);			// Resolution limits for cross-correlation
	double			nuavg(0), nustd(0);			// For rescaling
	Vector3<double>	origin;						// Helix origin
	int				set_origin(0);				// Flag to set origin
	int 			calc_background(0);			// Flag to calculate background values
	int 			norm_flag(0);				// Normalization flag
	int				zmin(0), zmax(1000000);		// Limits along the helical axis
	double			radius(0);					// Radial limit for symmetrization
	double			ell_angle(0), ell_shift(0);	// Elliptical distortion parameters
	int				set_views(0);				// Flag to reconstruct from projections
	Bstring			ps_file;					// Postscript output file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" ) {
			sam = curropt->scale();
			set_sampling = 1;
		}
		if ( curropt->tag == "helix" ) {
			if ( curropt->values(helix_rise, helix_angle, dyad_axis, cyclic) < 2 )
				cerr << "-helix: Both rise and angle must be specified" << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "rule" )
			if ( curropt->values(helix_repeat, turns, subunits) < 3 )
				cerr << "-rule: Repeat, turns and subunits must be specified" << endl;
		if ( curropt->tag == "super" ) {
			if ( curropt->values(helix_rise, helix_angle, offset[0], offset[1]) < 3 )
				cerr << "-super: Rise, angle and offset must be specified" << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "seam" )
			if ( ( seam_shift = curropt->value.real() ) < 0.1 )
				cerr << "-seam: A unit shift must be specified" << endl;
		if ( curropt->tag == "project" ) {
			if ( ( project_angle = curropt->value.real() ) < 0.1 )
				cerr << "-project: A projection angle step must be specified" << endl;
			else
				project_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "kernel" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-kernel: At least the kernel size must be specified!" << endl;
		if ( curropt->tag == "findcc" )
			if ( curropt->values(angle_start, angle_end, angle_step) < 3 )
				cerr << "-findcc: All three angles must be specified" << endl;
		if ( curropt->tag == "findrs" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 6 )
				cerr << "-findrs: All six rise and angle parameters must be specified" << endl;
			else {
				rise_start = d[0];
				rise_end = d[1];
				rise_step = d[2];
				angle_start = d[3];
				angle_end = d[4];
				angle_step = d[5];
			}
		}
		if ( curropt->tag == "segment" )
			if ( curropt->values(angle_start, angle_end, angle_step, thickness) < 4 )
				cerr << "-segment: All three angles and a thickness must be specified!" << endl;
		if ( curropt->tag == "elliptical" ) {
			if ( curropt->values(ell_angle, ell_shift) < 2 )
				cerr << "-elliptical: Both major axis angle and shift must be specified!" << endl;
			else
				ell_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "background" )
			calc_background = 1;
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "zlimits" )
			if ( curropt->values(zmin, zmax) < 2 )
				cerr << "-zlimits: Both limits on the helical axis must be specified" << endl;
		if ( curropt->tag == "normalize" )
			norm_flag = 1;
		if ( curropt->tag == "radius" )
			if ( ( radius = curropt->value.real() ) < 1 )
				cerr << "-radius: A radius must be specified" << endl;
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: Resolution limits must be specified!" << endl;
		if ( curropt->tag == "bin" )
			if ( ( find_bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: A bin value must be specified!" << endl;
 		if ( curropt->tag == "views" ) {
			if ( curropt->value[0] == 's' ) set_views = 1;
			if ( curropt->value[0] == 'r' ) set_views = 2;
 		}
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
   }
	option_kill(option);

	angle_start *= M_PI/180.0;
	angle_end *= M_PI/180.0;
	angle_step *= M_PI/180.0;
	
	double			ti = timer_start();
	
	Bstring			symmetry_string("H");
	Bsymmetry 		sym(symmetry_string);
	View*			views = NULL;
	Bplot*			plot = NULL;
	
	FSI_Kernel*		kernel = NULL;
	if ( kernel_width > 1 )
		kernel = new FSI_Kernel(kernel_width, kernel_power);

	int 			dataflag = 1;
	Bimage*			p = read_img(argv[optind++], dataflag, -1);
	Bimage* 		proj = NULL;
	
	if ( !p ) {
		cerr << "Error: No input file!" << endl;
		bexit(-1);
	}
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
	
	if ( set_sampling ) p->sampling(sam);
	
	if ( calc_background ) p->calculate_background();
	
	if ( subunits > 0 ) {
		helix_rise = helix_repeat*1.0L/subunits;
		helix_angle = TWOPI*turns*1.0L/subunits;
	}
	
	if ( angle_step > 0 ) {
		if ( thickness > 0 ) {
			plot = p->helix_segment_correlation(thickness,
				angle_start, angle_end, angle_step,
				find_bin, hires, lores, radius);
		} else if ( rise_step > 0 ) {
			plot = p->find_helix_parameters(rise_start, rise_end, rise_step,
				angle_start, angle_end, angle_step, find_bin, radius);
		} else {
			plot = p->find_helix_parameters(angle_start, angle_end, angle_step,
				find_bin, hires, lores, radius);
		}
	}
	
	if ( helix_rise && helix_angle ) {
		if ( offset[0] ) p->convert_to_helix(helix_rise, helix_angle, offset);
		else if ( seam_shift ) plot = p->seamed_helix_symmetrize(helix_rise, helix_angle,
				seam_shift, dyad_axis, zmin, zmax, radius, norm_flag);
		else p->helix_symmetrize(helix_rise, helix_angle,
				dyad_axis, zmin, zmax, radius, norm_flag);
	} else if ( p->sizeZ() > 1 && project_angle ) {
		views = asymmetric_unit_views(sym, TWOPI, project_angle, 0);	
		if ( kernel )
			proj = p->project(views, hires, kernel);
		else
			proj = p->project(views, 1);
		kill_list((char *) views, sizeof(View));
		delete p;
		p = proj;
	}

	if ( cyclic > 1 )
		p->symmetrize_cyclic(cyclic, norm_flag);

	if ( ell_shift )
		p->distort_elliptically(ell_angle, ell_shift);
		
	if ( ps_file.length() && plot ) ps_plot(ps_file, plot);
	
	delete plot;
	
	Bimage* 	pf = NULL;
	if ( set_views ) {
		pf = p->filament_from_projections(hires, set_views-1);
		delete p;
		p = pf;
	}

	if ( optind < argc ) {
		if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
		p->change_type(nudatatype);
		write_img(argv[optind], p, 0);
	}
	
	delete p;
	if ( kernel ) delete kernel;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

