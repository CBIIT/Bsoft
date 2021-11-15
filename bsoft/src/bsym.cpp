/**
@file	bsym.cpp
@brief	Program to generate symmetry axes for point group symmetries
@author Bernard Heymann
@date	Created: 20001119
@date	Modified: 20191104
**/

#include "rwimg.h"
#include "rwsymop.h"
#include "symmetry.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			img_write_symmetry_views(Bimage* p, Bsymmetry& sym, View ref_view,
				Bstring& filename, DataType datatype, double avg, double std);

// Usage assistance
const char* use[] = {
" ",
"Usage: bsym [options] in.img out.img",
"------------------------------------",
"Analyzes for and applies point group symmetry to images.",
"Generates symmetry axes for point group symmetries.",
" ",
"Actions:",
"-invert                  Invert density in the image.",
"-rescale -0.1,5.2        Rescale symmetrized data to average and standard deviation.",
"-axis 3,1                Rotate to desired axis, with flag for alternate (requires -symmetry).",
"-handedness 1.4          Check handedness above threshold (density is positive).",
"-check C,D,I             Check symmetry in correctly oriented map.",
"-find 1.5,2,1            Find orientation: angle step, binning, and",
"                         flag to skip major axis (use with -symmetry).",
"-change C5,C6            Change point group symmetry.",
"-Views                   Output all symmetry-related views, (filename: out_??.img).",
"-asu 3                   Output a single or multi-level mask of asymmetric units (use with -symmetry).",
"                            > 0, generate the corresponding asymetric unit.",
"                            == 0, generate all asymetric units.",
"-edge oval               Apply an edge smoothing of the indicated type.",
"-replicate               Replicate the asymmetric unit (use with -symmetry).",
"-nonorm                  Do not normalize after symmetrization (use with -symmetry).",
"-show                    Show operational symmetry matrices (use with -symmetry).",
"-pdb                     Show PDB symmetry matrices (use with -symmetry).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype f              Force writing of a new image data type.",
"-origin 0,22.5,30        Set the origin for rotation.",
"-reference 0,1.5,-0.2,35 Reference symmetry axis and rotation angle (default 0,0,1,0).",
"-resolution 15.6,200     Resolution limits for cross-correlation (default 0).",
"-symmetry C5             Point group symmetry.",
"-fill 0.02               Fill value (default image background).",
"-radius 34.5,1.2         Radius in voxels of repeats for changing symmetry.",
"                         and slope along z to adjust the radius.",
" ",
"Input:",
"-Template temp.pif       Template to find symmetry equivalent.",
"-mask mask.mrc           Mask (must be the same size as the template).",
" ",
"Output:",
"-Parameters out.star     Output STAR file with symmetry axes and orders.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	int 			setinvert(0);				// Flag to invert density before operations
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	long			axis(0), axis_flag(0);		// Axis to rotate to
	Bstring			check_sym;					// String of symmetries to check for
	double			find_angle(0);				// Angle step size for symmetry search
	int				find_bin(1);				// Binning for searching
	int				flags(0);					// Only search for minor axes
	double			radius(0);					// Radius to derive shift for symmetry change
	double			z_slope(0);					// Slope along z to adjust the radius
	Vector3<double>	origin;						// Origin
	int				set_origin(0); 				// Flag to set origin
	View			ref_view;					// Reference view
	double			hires(0), lores(0);			// Resolution limits for cross-correlation
	Bsymmetry		sym;						// Point group
	Bsymmetry		symnu;						// New point group
	int				norm_flag(1);				// Normalization flag
	int 			fill_type(FILL_BACKGROUND);	// Use background
	double			fill(0);					// Default fill is zero
	int				write_all_views(0);			// Flag to output al symmetry-related views
	int 			hand_flag(0);				// Flag for checking of handedness
	double			hand_threshold(0); 			// Threshold for checking handedness
	int				asu_mask(-1);				// Flag to generate a mult-level mask
	int				setedge(-1);				// Edge type; default none
	int				replicate(0);				// Flag to replicate the asymmetric unit
	int				show(0);					// Flag to show operational matrices
	Bstring			tempfile;					// Input template for symmetry equivalents
    Bstring			maskfile;					// Mask file name
	Bstring			pgfile;						// Output point group file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "invert" )
			setinvert = 1;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "axis" ) 
			if ( curropt->values(axis, axis_flag) < 1 )
				cerr << "-axis: An axis order must be specified!" << endl;
		if ( curropt->tag == "check" ) check_sym = curropt->value;
		if ( curropt->tag == "find" ) {
			if ( curropt->values(find_angle, find_bin, flags) < 1 )
				cerr << "-find: An angle step size must be specified!" << endl;
			else
				find_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "change" ) {
			sym = Bsymmetry(curropt->value.pre(','));
			symnu = Bsymmetry(curropt->value.post(','));
		}
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "reference" )
			ref_view = curropt->view();
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: Resolution limits must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "Views" )
			write_all_views = 1;
		if ( curropt->tag == "asu" )
			if ( ( asu_mask = curropt->value.integer() ) < 0 )
				cerr << "-asu: An asymmetric unit index must be specified!" << endl;
		if ( curropt->tag == "edge" ) {
			if ( curropt->value.length() < 1 ) {
				cerr << "-edge: An edge shape must be specified." << endl;
			} else {
				if ( tolower(curropt->value[0]) == 'r' ) setedge = 0;
				if ( tolower(curropt->value[0]) == 'o' ) setedge = 1;
				if ( tolower(curropt->value[0]) == 'c' ) setedge = 2;
			}
		}
		if ( curropt->tag == "replicate" )
			replicate = 1;
		if ( curropt->tag == "nonorm" )
			norm_flag = 0;
		if ( curropt->tag == "show" ) show |= 1;
		if ( curropt->tag == "pdb" ) show |= 2;
		if ( curropt->tag == "handedness" ) {
			if ( ( hand_threshold = curropt->value.real() ) )
				cerr << "-handedness: A density threshold must be specified!" << endl;
			else
				hand_flag = 1;
		}
		if ( curropt->tag == "radius" )
			if ( curropt->values(radius, z_slope) < 1 )
				cerr << "-radius: A radius in voxels must be specified!" << endl;
		if ( curropt->tag == "Template" )
			tempfile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "Parameters" )
			pgfile = curropt->filename();
    }
	option_kill(option);	
	
	double		ti = timer_start();

	if ( show & 1 ) sym_show_operational_matrices(sym);
	if ( show & 2) sym_show_pdb_matrices(sym);
	
	if ( pgfile.length() )
		write_pointgroup(pgfile, sym, ref_view);
	
	// Read image file
	int 		dataflag(0);
	Bimage*		p = NULL;
	if ( hand_flag || optind < argc - 1 ) dataflag = 1;
	if ( optind < argc ) {
		p = read_img(argv[optind++], dataflag, -1);
		if ( !p ) {
			cerr << "Error: File \"" << argv[optind-1] << "\" was not read!" << endl;
			bexit(-1);
		}
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();
	} else
		bexit(0);
	
	if ( setinvert ) p->invert();
	
	if ( fill_type == FILL_AVERAGE ) fill = p->average();
	if ( fill_type == FILL_BACKGROUND ) {
		if ( !p->background(long(0)) ) p->calculate_background();
		fill = p->background(long(0));
	}
	
	p->background(fill);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->size()/2);
		else p->origin(origin);
	}
	
	Bimage*			ptemp = NULL;
	if ( tempfile.length() ) {
		ptemp = read_img(tempfile, 1, 0);
		if ( !ptemp ) {
			cerr << "Error:File \"" << tempfile << "\" was not read!" << endl;
			bexit(-2);
		}
	}
	
	Bstring			filename;
	Bimage*			pmask = NULL;
	Vector3<double>	start;
	Matrix3			mat;
	
	if ( optind < argc ) {
		filename = argv[optind];
		if ( p ) {
			if ( check_sym.length() ) {
				p->check_point_group(check_sym);
			}
			if ( sym.point() > 101 ) {
				if ( write_all_views ) {
					img_write_symmetry_views(p, sym, ref_view,
						filename, nudatatype, nuavg, nustd);
				}
				if ( asu_mask >= 0 ) {
					pmask = p->levelmask_asymmetric_units(sym, asu_mask);
					delete p;
					p = pmask;
				} else if ( find_angle > 0 ) {
					p->find_point_group(sym, find_angle, find_bin, hires, lores, flags);
				} else if ( replicate > 0 ) {
					p->replicate_asymmetric_unit(sym);
				} else if ( symnu.point() > 101 ) {
					p->change_symmetry(sym, symnu, radius, z_slope);
				} else if ( axis > 0 ) {
					p->rotate_to_axis(sym, axis, axis_flag);
				} else if ( ptemp ) {
					if ( maskfile.length() )
						pmask = read_img(maskfile, 1, -1);
					mat = p->symmetry_equivalent(ptemp, pmask, sym);
					p->rotate(mat);
				} else {
					p->symmetrize(sym, ref_view, norm_flag);
				}
			}
		}
		if ( setedge >= 0 ) p->edge(setedge, p->size(), start, 2, FILL_BACKGROUND, 0);
		if ( nustd ) p->rescale_to_avg_std(nuavg, nustd);
		p->change_type(nudatatype);
		write_img(filename, p, 0);
	}
	
	delete p;
	delete ptemp;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Rotates an image to all symmetry-related views and writes out the resultant images.
@param 	*p				image to be rotated and saved.
@param 	*sym			point group symmetry.
@param 	ref_view		reference view (default 0,0,1,0).
@param 	&filename		output filename (converted to name_??.img)
@param 	datatype		new data type.
@param 	avg				new average.
@param 	std				new standard deviation.
@return int				0.

	The resultant image is rescaled to the desired average and standard
	deviation, and converted to the desired data type. If the standard
	deviation is zero, no rescaling is applied. The output filenames
	are numbered with an underscore and two digits, starting at 1.

**/
int			img_write_symmetry_views(Bimage* p, Bsymmetry& sym, View ref_view,
				Bstring& filename, DataType datatype, double avg, double std)
{
	int				i;
	Bstring			outname;
	Bimage*			prot = NULL;
	View*			v;
	View*			view = symmetry_get_all_views(sym, ref_view);

	for ( i=1, v=view; v; v=v->next, i++ ) {
		prot = p->rotate(p->size(), *v);
		outname = filename.pre_rev('.') + Bstring(i, "_%02d.") + filename.post_rev('.');
		if ( std ) prot->rescale_to_avg_std(avg, std);
		prot->change_type(datatype);
		write_img(outname, prot, 0);
		delete prot;
	}
	
	kill_list((char *) view, sizeof(View));

	return 0;
}


