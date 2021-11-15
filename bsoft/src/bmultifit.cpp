/**
@file	bmultifit.cpp
@brief	Searching for a template in a map and returning multiple hits in a model.
@author Bernard Heymann
@date	Created: 20021027
@date	Modified: 20160604
**/

#include "model_multifit.h"
#include "ps_views.h"
#include "linked_list.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bmultifit [options] input.img view.grd",
"---------------------------------------------",
"Searches orientation space for a set of solutions for the fit of a template within a 3D map.",
"Modes: Global: Search an asymmetric unit with a fixed angular step size.",
"       Symmetry: Search all symmetry-related views using the view from the -View option.",
"The output image contains 4-valued views, only supported by the GRD format.",
"Only one image in a multi-image file will be used (see -image option).",
" ",
"Actions:",
"-mode sym                Mode: global (default) or symmetric.",
"-subset 25,5             View subset selection from a start offset and with a given size.",
"-fom 0.35                FOM cutoff to use (default set to FOMmax/2).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 110,50,44        Origin for rotations of the template (default taken from template image).",
"-angles 8.8,2.5,5.2      Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-resolution 15.6,200     Resolution limits for cross-correlation (default 0).",
"-View 0.3,-0.5,0.8,33    View to initiate search: vector {xyz} and angle (use with -symmetry option).",
"-symmetry D6             Symmetry: Point group identifier (use with -View option).",
"-image 2                 Image number in a multi-image file (default 0=first image).",
"-bin 3                   Binning of input, template and mask by given kernel size,",
"                         output image remains same size as input.",
" ",
"Input:",
"-Template image.map      Template to search for.",
"-Mask mask.tif           Mask in reciprocal space to apply during cross-correlation.",
" ",
"Output:",
"-output file.star        Model output file with best orientations.",
"-Postscript plot.ps      Output postscript file with a plot of projection vectors.",
"-FOM cc.map              Output correlation coefficient map.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	Bstring			mode("global");					// Mode, default global
	int				view_start(0);					// Subset start
	int				view_subset(0);					// Subset size
	double			cutoff(0);						// Cutoff set to 0 to calculate automatically
	int				img_num(0);						// Image number in file
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Vector3<double>	origin;							// Template origin
	int				set_origin(0);					// Flag to set origin
	double			alpha(-1);						// Rotation around view vector
	double			alpha_step = M_PI_4;			// Angular step size for alpha
	double			theta_step = M_PI_4;			// Angular step size for theta
	double			phi_step = M_PI_4;				// Angular step size for phi
	Bstring			symmetry_string("C1");			// Default: asymmetric or C1 point group
	double			hires(0), lores(0);				// Limiting resolution for cross-correlation
	View			currview(0,0,-10,0);			// View to initiate search from and use as current view
	View			ref_view;						// Reference view
	Vector3<long>	bin = {1,1,1};					// No binning
	Bstring			template_file;
	Bstring			mask_file;
	Bstring			model_file;
	Bstring			ps_file;
	Bstring			fom_file;
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "mode" ) {
			mode = curropt->value.lower();
			if ( mode.length() < 1 )
				cerr << "-mode: A mode must be specified" << endl;
		}
		if ( curropt->tag == "subset" )
			if ( curropt->values(view_start, view_subset) < 2 )
				cerr << "-subset: A view subset start and size must be specified" << endl;
		if ( curropt->tag == "fom" )
			if ( ( cutoff = curropt->value.real() ) < 0.0001 )
				cerr << "-fom: A cutoff value must be specified!" << endl;
		if ( curropt->tag == "image" )
			if ( ( img_num = curropt->value.integer() ) < 1 )
				cerr << "-image: An image number must be specified" << endl;
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
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(alpha_step, theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified" << endl;
			else {
				if ( i < 2 ) phi_step = theta_step = alpha_step;
				else if ( i < 3 ) phi_step = theta_step;
				alpha_step *= M_PI/180.0;
				theta_step *= M_PI/180.0;
				phi_step *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution must be specified!" << endl;
		if ( curropt->tag == "View" )
			currview = curropt->view();
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "bin" ) {
			if ( ( bin[0] = curropt->value.integer() ) < 1 )
				cerr << "-bin: A bin size must be specified!" << endl;
			else
				bin[2] = bin[1] = bin[0];
		}
		if ( curropt->tag == "Template" )
			template_file = curropt->filename();
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "output" )
			model_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
		if ( curropt->tag == "FOM" )
			fom_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	// Do the symmetry-related stuff first
	View*			view = NULL;
	
	Matrix3			mat = ref_view.matrix();
	Bsymmetry 		sym(symmetry_string);
	sym.transform(mat);
//	for ( i=0; i<sym.operations(); i++ )
//		sym.op[i].axis = mat * sym.op[i].axis;
	
	if ( mode.contains("glo") ) {
		view = asymmetric_unit_views(sym, theta_step, phi_step, 1);
	} else if ( mode.contains("sym") ) {
		if ( currview[2] < -1 ) currview[2] = 1;
		view = symmetry_get_all_views(sym, currview);
		alpha_step = TWOPI;
	}

	if ( view_subset > 0 ) if ( view_list_subset(&view, view_start, view_subset) < 1 ) {
		cerr << "Error: At least one view must be selected!" << endl;
		bexit(-1);
	}

	if ( ps_file.length() && view ) ps_views(ps_file, view);

	if ( verbose )
		cout << "Number of views:                " << count_list((char *) view) << endl;
	
	// Read image file
	Bimage*			p = NULL;
	Bimage*			ptemp = NULL;
	Bimage*			pmask = NULL;
	Vector3<double>	scale(1,1,1);
	Vector3<long>	size;
	
	if ( optind < argc ) {
		p = read_img(argv[optind++], 1, img_num);
		size = p->size();
		if ( bin[0]*bin[1]*bin[2] > 1 ) p->bin(bin);
	}
	
	if ( !p ) {
		cerr << "Error: No input file given!" << endl;
		bexit(-1);
	}

	if ( template_file.length() ) {
		ptemp = read_img(template_file, 1, 0);
		if ( sam.volume() > 0 ) ptemp->sampling(sam);
		if ( set_origin ) {
			if ( set_origin == 2 ) ptemp->origin(ptemp->size()/2);
			else ptemp->origin(origin);
		} else origin = ptemp->image->origin();
		if ( bin[0]*bin[1]*bin[2] > 1 ) ptemp->bin(bin);
	}
	
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, -1);
		if ( bin[0]*bin[1]*bin[2] > 1 ) pmask->bin(bin);
	}
	
	Bmodel*			model = NULL;
	Bimage*			pfit = NULL;
	
	if ( ptemp ) {
		if ( model_file.length() ) {
			model = model_from_densities(p, ptemp, view, alpha, alpha_step, hires, lores, pmask, cutoff);
				write_model(model_file, model);
			model_kill(model);
		} else {
//			pfit = img_multifit(p, ptemp, view, alpha, alpha_step, hires, lores, pmask, cutoff);
			pfit = p->search_volume(ptemp, view, alpha, alpha_step, hires, lores, pmask, cutoff);
			if ( optind < argc )
				write_img(argv[optind], pfit, 0);
			if ( fom_file.length() )
				write_img(fom_file, pfit->next, 0);
//			delete pfit->next;
			delete pfit;
		}
		delete ptemp;
	}

	kill_list((char *) view, sizeof(View));
	delete p;
	delete pmask;

	if ( !pfit ) {
		cerr << "Error: No output file generated!" << endl;
		bexit(-1);
	}

	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

