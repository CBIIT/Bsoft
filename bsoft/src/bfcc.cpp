/**
@file	bfcc.cpp
@brief	Fast cross-correlation search for the best fit of a 3D map to a template.
@author Bernard Heymann
@date	Created: 20130523
@date	Modified: 20160604
**/

#include "rwimg.h"
#include "mg_processing.h"
#include "mg_particle_select.h"
#include "rwmg.h"
#include "linked_list.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
double		img_fcc_search(Bimage* p, Bimage* ptemp, View* view,
				FSI_Kernel* kernel, double hires, double lores, Bimage* pmask);
double		img_fcc_for_view(Bimage* p, Bimage* ptemp, View view,
				FSI_Kernel* kernel, double hires, double lores, Bimage* pmask);
double		img_cc_for_view(Bimage* p, Bimage* ptemp, View view,
				FSI_Kernel* kernel, double hires, double lores, double searchrad,
				Bimage* pmask);

// Usage assistance
const char* use[] = {
" ",
"Usage: bfcc [options] input.img/input.star output.img",
"-----------------------------------------------------",
"Fast cross-correlation for the best fit of a 3D map to a template.",
"Modes: Global: Search an asymmetric unit with a fixed angular step size.",
"       Directional: Do a search around a view within an angular distance (-limit option).",
"       Refine: Refine iteratively around a view to some accuracy (-accuracy option).",
"       Symmetry: Search all symmetry-related views using a given view.",
"For modes other than global, the view is taken from the command line (-View option),",
"	or from the parameter file, or as a last resort from the image header.",
"If a parameter file is used as input, the input image is selected based on ",
"	the reconstruction ID and the particle image number.",
"Only one image in a multi-image file will be used (see -image option).",
"The output file contains the transformed map only if the input is an image (i.e., not a parameter file).",
"Output of the transformed map or template can be specified with the",
"	-particle and -newtemplate options.",
" ",
"Actions:",
"-global                  Global search before other modes.",
"-mode refine             Local modes: directional, refine or symmetric.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 110,50,44        Template image origin (default taken from template file).",
"-partorigin 50,80,47     Particle image origin (default taken from particle file).",
"-kernel 6,2              Kernel for reciprocal space interpolation: size and power.",
"-resolution 15.6,200     Resolution limits for cross-correlation (default 0).",
"-shiftlimit 5.2          Radial limit to search for shift in cross-correlation map (default 1/4 of map).",
" ",
"Selections:",
"-setsequential           Set particle selection numbers sequentially.",
"-recid A435              The 3D reconstruction ID (should correspond to that in the parameter file).",
"-image 2                 Image number in a multi-image file or particle number (default 1=first image).",
"-select 5                Particle selection number in parameter file (default not used).",
" ",
"View parameters:",
"-angles 8.8,2.5,5.2      Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-side 15                 Generate side view projections within the given angle from the equator.",
" ",
"Parameters for the global and symmetric modes:",
"-symmetry D6             Symmetry: Point group identifier.",
" ",
"Parameters for the directional and refinement modes:",
"-View 0.3,-0.5,0.8,33    View to initiate search: vector {xyz} and angle.",
" ",
"Parameters for the directional mode:",
"-limit 17.5              Radial limit around given view (default 20 degrees).",
" ",
"Parameters for the refinement mode:",
"-accuracy 0.3            Accuracy (default 0.5 degrees).",
" ",
"Input:",
"-Template image.map      Template to search for.",
"-Mask mask.tif           Mask in reciprocal space to apply during cross-correlation.",
" ",
"Output:",
//"-Postscript plot.ps      Output postscript file with a plot of global view vectors.",
"-output file.star        Parameter output file with best orientation.",
"-particle newpart.pif    New transformed particle file.",
"-newtemplate newtemp.pif New transformed template file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	int				global(0);					// Flag for global search
	Bstring			mode;						// Mode: directional, refine or symmetric
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Bstring			rec_id;						// Reconstruction ID
	int				img_num(1);					// Image number in file
	int				sel_num(-1);				// Selection number in parameter file
	int				set_seq_sel(0);				// Flag to set selection number sequentially
	Vector3<double>	sam;				// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// Template origin
	Vector3<double>	part_ori;					// Particle origin
	int				set_origin(0);				// Flag to set origin
	int				kernel_width(6);			// Reciprocal space kernel width
	int				kernel_power(2);			// Reciprocal space kernel power
	double			alpha_step(M_PI_4);			// Angular step size for alpha
	double			theta_step(M_PI_4);			// Angular step size for theta
	double			phi_step(M_PI_4);			// Angular step size for phi
	double			side_ang(-1);				// Side view variation angle
	double			angle_limit(M_PI/9.0);		// Angular limit for directional search
	Bsymmetry		sym;					// Default: asymmetric or C1 point group
	double			hires(0), lores(0);			// Limiting resolution for cross-correlation
	double			search_radius(-1);			// Use default search radius
	int				currview_set(0);			// Flag to indicate current view is set
	View			currview;					// View to initiate search from and use as current view
	View			bestview;					// Eventually the best view
	View			ref_view;					// Reference view
	double			accuracy(M_PI/360.0);		// Accuracy for refinement default = 0.5 degrees
	Bstring			template_file;				// Reference template file name
	Bstring			mask_file;					// Reciprocal space mask file name
	Bstring			ps_file;					// Postscript output file name
	Bstring			param_file;					// Parameter file name
	Bstring			newpart_file;				// Transformed particle output file name
	Bstring			newtemp_file;				// Transformed template output file name
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "global" ) global = 1;
		if ( curropt->tag == "mode" ) {
			if ( curropt->value.length() > 0 ) {
				if ( curropt->value[0] == 'd' || curropt->value[0] == 'D' ) mode = "directional";
				if ( curropt->value[0] == 'r' || curropt->value[0] == 'R' ) mode = "refine";
				if ( curropt->value[0] == 's' || curropt->value[0] == 'S' ) mode = "symmetric";
			} else
				cerr << "-mode: A mode must be specified" << endl;
		}
		if ( curropt->tag == "setsequential" ) set_seq_sel = 1;
		if ( curropt->tag == "image" )
			if ( ( img_num = curropt->value.integer() ) < 1 )
				cerr << "-image: An image number must be specified" << endl;
		if ( curropt->tag == "select" )
			if ( ( sel_num = curropt->value.integer() ) < 0 )
				cerr << "-select: A selection number must be specified" << endl;
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
		if ( curropt->tag == "partorigin" )
			part_ori = curropt->origin();
		if ( curropt->tag == "kernel" )
			if ( curropt->values(kernel_width, kernel_power) < 1 )
				cerr << "-kernel: At least the kernel size must be specified!" << endl;
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
		if ( curropt->tag == "side" ) {
			if ( ( side_ang = curropt->value.real() ) < 0.1 )
				cerr << "-side: One angle must be specified!" << endl;
			else
				side_ang *= M_PI/180.0;
		}
		if ( curropt->tag == "limit" ) {
			if ( ( angle_limit = curropt->value.real() ) < 0.1 )
				cerr << "-limit: An angular limit must be specified!" << endl;
			else
				angle_limit *= M_PI/180.0;
		}
		if ( curropt->tag == "View" ) {
			currview = curropt->view();
			currview_set = 1;
		}
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution must be specified!" << endl;
		if ( curropt->tag == "shiftlimit" )
			if ( ( search_radius = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A shift search radius must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "accuracy" ) {
			if ( ( accuracy = curropt->value.real() ) < 0.1 )
				cerr << "-accuracy: An accuracy in degrees must be specified!" << endl;
			else
				accuracy *= M_PI/180.0;
		}
		if ( curropt->tag == "recid" )
			rec_id = curropt->value;
		if ( curropt->tag == "Template" )
			template_file = curropt->filename();
		if ( curropt->tag == "Mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
		if ( curropt->tag == "output" )
			param_file = curropt->filename();
		if ( curropt->tag == "particle" )
			newpart_file = curropt->filename();
		if ( curropt->tag == "newtemplate" )
			newtemp_file = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	int					read_img_num = img_num - 1;
	int					cnt(1), found;
	Bimage*				p = NULL;
	Bimage*				ptemp = NULL;
	Bimage*				pmask = NULL;
	Bstring				filename;
	Vector3<double>		scale(1,1,1);
	Vector3<long>		size;
	Bproject*			project = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
//	if ( hires > lores ) swap(hires, lores);

	if ( kernel_width < 4 ) kernel_width = 4;
	if ( kernel_power < 2 ) kernel_power = 2;
	FSI_Kernel*			kernel = new FSI_Kernel(kernel_width, kernel_power);
	
	// Read image file
	if ( optind < argc ) {
		filename = argv[optind++];
		if ( file_type(filename) == Micrograph ) {
			project = read_project(filename);
			part_set_FOM(project, 0, 0);
			if ( set_seq_sel )
				for ( rec = project->rec; rec; rec = rec->next )
					for ( part = rec->part; part; part = part->next ) part->sel = cnt++;
			rec = project->rec;
			if ( sel_num > 0 ) {		// Finds the first particle with this selection number
				for ( found=0, rec = project->rec; rec && !found; ) {
					for ( part = rec->part; part && !found; ) {
						if ( part->sel == sel_num ) found = 1;
						else part = part->next;
					}
					if ( !found ) rec = rec->next;
				}
			} else {			// Finds the specified particle and image number
				if ( rec_id.length() )
					for ( rec = project->rec; rec && rec_id != rec->id; rec = rec->next ) ;
				if ( rec ) {
					for ( part = rec->part; part && part->id != img_num; part = part->next ) ;
				} else {
					cerr << "Error: No reconstruction with ID = " << rec_id << endl;
					bexit(-1);
				}
			}
			if ( rec ) {
				filename = rec->frec;
				if ( rec->fpart.length() ) filename = rec->fpart;
				if ( part ) {
					img_num = part->id;
					read_img_num = part->id - 1;
					if ( !currview_set ) {
						currview = part->view;
						currview_set = 1;
					}
//					cout << "part->ori=" << part->ori << endl;
					if ( part_ori.length() < 1 && part->ori.length() ) part_ori = part->ori;
					if ( part->fpart.length() ) {
						filename = part->fpart;
						read_img_num = 0;
					}
				} else {
					cerr << "Error: Particle not defined!" << endl;
					bexit(-1);
				}
				if ( !currview_set ) {
					currview = rec->view;
					currview_set = 1;
				}
			} else {
				cerr << "Error: Reconstruction not defined!" << endl;
				bexit(-1);
			}
		}
		if ( filename.length() < 1 ) {
			cerr << "Error: No image file name found!" << endl;
			bexit (-1);
		}
	}
	
//	cout << "part_ori=" << part_ori << endl;

	if ( filename.length() < 1 ) {
		cerr << "Error: No image file name specified!" << endl;
	} else {
		p = read_img(filename, 1, read_img_num);
		if ( p ) {
			if ( sam.volume() > 0 ) p->sampling(sam);
			if ( part_ori.length() ) p->origin(part_ori);
			else part_ori = p->image->origin();
		} else {
			cerr << "Error: %s not read!" << filename << endl;
		}
	}

	if ( template_file.length() < 1 ) {
		cerr << "Error: No input template file specified!" << endl;
	} else {
		ptemp = read_img(template_file, 1, 0);
		if ( ptemp ) {
			if ( sam.volume() > 0 ) ptemp->sampling(sam);
			if ( set_origin ) {
				if ( set_origin == 2 ) ptemp->origin(ptemp->size()/2);
				else ptemp->origin(origin);
			} else origin = ptemp->image->origin();
		} else {
			cerr << "Error: " << template_file << " not read!" << endl;
		}
	}

	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, -1);
		if ( !pmask ) {
			cerr << "Error: " << mask_file << " not read!" << endl;
		}
	}

	if ( search_radius < 1 ) search_radius = ptemp->sizeX()/4;
	
	p->image->origin(origin);
	
	if ( verbose ) {
		cout << "Finding a template in a map:" << endl;
		cout << "Template:                       " << ptemp->file_name() << endl;
		cout << "Map:                            " << p->file_name() << " (" << read_img_num << ")" << endl;
		if ( pmask ) cout << "Mask:                           " << pmask->file_name() << endl;
		cout << "Template origin:                " << ptemp->image->origin() << endl;
		cout << "Map origin:                     " << p->image->origin() << endl;
		cout << "Resolution range:               ";
		if ( lores <= 0 ) cout << "inf";
		else cout << lores;
		cout << " - " << hires << " A" << endl;
		cout << "Shift search radius:            " << search_radius << " pixels" << endl;
		cout << "Kernel width and power:         " << kernel_width << " " << kernel_power << endl << endl;
		cout << "Angle\tShift\t\t\tView\t\t\t\tCC" << endl;
	}
	

	Bimage*			pc = p->pack_two_in_complex(ptemp);
	if ( !pc ) bexit(-1);
	
	delete p;
	delete ptemp;

//	cout << "Doing forward transform" << endl;
	pc->fft(FFTW_FORWARD, 0);

//	cout << "Unpacking transform" << endl;
	ptemp = pc->unpack_combined_transform();
	p = pc;
	
//	cout << "Phase shifting" << endl;
	p->origin(origin);
	ptemp->origin(origin);
	p->phase_shift_to_origin();
	ptemp->phase_shift_to_origin();
	
	View*			view = NULL, *view2 = NULL;
	Vector3<double>	currshift, bestshift;
	double			cc, best_cc = -1e37;
	
	Matrix3			mat = ref_view.matrix();
	sym.transform(mat);
//	for ( i=0; i<sym.operations(); i++ )
//		sym.op[i].axis = mat * sym.op[i].axis;
	
	bestview = currview;

	if ( global || !currview_set ) {
//		cout << "Doing global" << endl;
		
		if ( side_ang < 0 )
			view = asymmetric_unit_views(sym, theta_step, phi_step, 1);
		else
			view = side_views(sym, side_ang, theta_step, phi_step);
		view2 = view_list_expand_angles(view, -M_PI, M_PI - alpha_step/2, alpha_step);
		kill_list((char *) view, sizeof(View));
		view = view2;
		
		if ( verbose & VERB_FULL )
			show_views(view);

//		if ( ps_file.length() && view ) ps_views(ps_file, symmetry_string, view, 0);
	
//		if ( verbose )
//			cout << "Number of global views:         " << count_list((char *) view) << endl << endl;
	
		best_cc = img_fcc_search(p, ptemp, view, kernel, hires, lores, pmask);
		bestview = p->image->view();

		best_cc = img_cc_for_view(p, ptemp, bestview, kernel, hires, lores,
			search_radius, pmask);
		bestshift = ptemp->image->origin() - p->image->origin();
		p->phase_shift(bestshift);
		
		if ( verbose )
			cout << alpha_step*180.0/M_PI << tab << bestshift << tab << bestview << tab << best_cc << endl;
		
		kill_list((char *) view, sizeof(View));
		view = NULL;
	}

	if ( mode == "directional" ) {
		view = views_within_limits(bestview, theta_step, phi_step, alpha_step, angle_limit, M_PI);
		if ( verbose ) {
			cout << "Central view:                   " << bestview << endl;
			cout << "Number of directional views:    " << count_list((char *) view) << endl << endl;
		}
	} else if ( mode == "symmetric" ) {
		view = symmetry_get_all_views(sym, bestview);
		if ( verbose )
			cout << "Number of symmetry views:       " << count_list((char *) view) << endl << endl;
	}

	if ( mode.length() ) {
		if ( mode == "directional" || mode == "symmetric" ) {
			best_cc = img_fcc_search(p, ptemp, view, kernel, hires, lores, pmask);
			bestview = p->image->view();

			kill_list((char *) view, sizeof(View));
		} else if ( mode == "refine" ) {
			currview = bestview;
			best_cc = -1;
			while ( alpha_step >= accuracy ) {
				view = views_for_refinement(bestview, alpha_step);
				cc = img_fcc_search(p, ptemp, view, kernel, hires, lores, pmask);
				currview = p->image->view();
				if ( currview.residual(bestview) < 1e-6 ) alpha_step /= 2;	// Contract around best point
				if ( best_cc < cc ) {
					best_cc = cc;
					bestview = currview;
//					bestshift = currshift;
				}

				if ( verbose )
					cout << alpha_step*180.0/M_PI << tab << bestshift << tab << bestview << tab << best_cc << endl;

				kill_list((char *) view, sizeof(View));
			}
		}
		
		best_cc = img_cc_for_view(p, ptemp, bestview, kernel, hires, lores,
			search_radius, pmask);
		bestshift += ptemp->image->origin() - p->image->origin();
	}

	if ( verbose )
		cout << "Best view:\t" << setprecision(4) <<
			bestshift[0] << tab << bestshift[1] << tab << bestshift[2] << tab <<
			bestview[0] << tab << bestview[1] << tab << bestview[2] << tab << 
			bestview.angle()*180/M_PI << tab << best_cc << endl;

	delete p;
	delete pmask;
	delete ptemp;
	delete kernel;
	
	if ( newpart_file.length() < 1 && optind < argc ) newpart_file = argv[optind];

	if ( newtemp_file.length() ) {  // Output transformed template
		ptemp = read_img(template_file, 1, 0);
//		size = ptemp->size();
		mat = bestview.matrix();
		mat = mat.transpose();
		ptemp->calculate_background();
		ptemp->transform(scale, origin, -bestshift, mat, FILL_BACKGROUND, 0);
		ptemp->view(bestview);
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();	// Preserve the old type
		ptemp->change_type(nudatatype);
		write_img(newtemp_file, ptemp, 0);
		delete ptemp;
	}
	
	if ( newpart_file.length() ) {	// Output transformed image
		p = read_img(filename, 1, read_img_num);
//		size = p->size();
		mat = bestview.matrix();
//		origin += bestshift;	// Origin for rotation
		p->calculate_background();
		p->transform(scale, origin-bestshift, bestshift, mat, FILL_BACKGROUND, 0);
		p->image->view(0,0,1,0);
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();	// Preserve the old type
		p->change_type(nudatatype);
		write_img(newpart_file, p, 0);
		delete p;
	}
	
	if ( param_file.length() ) {
		p = read_img(filename, 0, 0);
		if ( !project ) {
			project = new Bproject;
			if ( rec_id.length() < 1 ) {
				rec_id = p->file_name();
				rec_id = rec_id.pre_rev('.');
			}
			rec = reconstruction_add(&project->rec, rec_id);
			part = particle_add(&rec->part, img_num);
		}
		if ( sam[0] > 0 ) rec->voxel_size = sam;
		else rec->voxel_size = p->sampling(0);
		rec->symmetry = sym.label();
		if ( part ) {
			part->fpart = filename;
			if ( newpart_file.length() ) {
				part->fpart = newpart_file;
				part->ori = origin;
			} else {
				part->view = bestview;
				part->ori = origin - bestshift;
			}
			for ( part = rec->part; part; part = part->next )
				if ( part->id == img_num ) part->fom[0] = best_cc;
				else part->fom[0] = 0;
		} else {
			rec->frec = filename;
			if ( newpart_file.length() ) {
				rec->frec = newpart_file;
				rec->origin = origin;
			} else {
				rec->view = bestview;
				rec->origin = origin - bestshift;
			}
			for ( rec = project->rec; rec; rec = rec->next )
				if ( rec->id == rec_id ) rec->fom = best_cc;
				else rec->fom = 0;
		}
		write_project(param_file, project, 0, 0);
		delete p;
	}

	project_kill(project);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Fast cross-correlation searches a 2D/3D density map for a template.
@param 	*p			the Fourier transformed image.
@param 	*ptemp		the Fourier transformed template to be searched for.
@param 	view		view.
@param 	*kernel		interpolation kernel.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		mask for cross-correlation (ignored if NULL).
@return double		correlation coefficient.

	The template is rotated to the view and cross-correlated to find
	a set of high-scoring fits.
	The views must be calculated externally to allow for custom sets.
	The best view is returned in the image view record.

**/
double		img_fcc_search(Bimage* p, Bimage* ptemp, View* view,
				FSI_Kernel* kernel, double hires, double lores, Bimage* pmask)
{
	if ( !p ) return -1e37;

	long 			n, nviews, ibest;
	double			best(-1e37);

	View*			view_arr = view_array(view, nviews);

	double*			cc = new double[nviews];
	
	if ( verbose & VERB_PROCESS  )
		cout << "Number of views:                " << nviews << endl;
	
#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t k){
		cc[k] = img_fcc_for_view(p, ptemp, view_arr[k], kernel, hires, lores, pmask);
	});
#else
#pragma omp parallel for
	for ( n=0; n<nviews; n++ )
		cc[n] = img_fcc_for_view(p, ptemp, view_arr[n], kernel, hires, lores, pmask);
#endif

	double			ccmin(1), ccmax(-1), ccavg(0), ccstd(0);
	
	for ( n=ibest=0; n<nviews; n++ ) {
		if ( best < cc[n] ) {
			best = cc[n];
			ibest = n;
		}
		if ( ccmin > cc[n] ) ccmin = cc[n];
		if ( ccmax < cc[n] ) ccmax = cc[n];
		ccavg += cc[n];
		ccstd += cc[n]*cc[n];
	}

	ccavg /= n;
	ccstd = ccstd/n - ccavg*ccavg;
	if ( ccstd > 0 ) ccstd = sqrt(ccstd);
	else ccstd = 0;

	if ( verbose & VERB_PROCESS )
		cout << "Min, max, avg, std:             " << ccmin << " " << ccmax << " " << ccavg << " " << ccstd << endl;

	p->view(view_arr[ibest]);
	p->image->FOM(best);
	
	delete[] view_arr;
	delete[] cc;
	
	return best;
}

/*
@brief Searches a 2D/3D density map for a template using a specific view.
@param 	*p				the Fourier transformed image.
@param 	*ptemp			the Fourier transformed template to be searched for.
@param 	view			view.
@param 	*cc				interpolation kernel.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	*pmask			mask for cross-correlation (ignored if NULL).
@return double			correlation coefficient.

	The template is rotated to the view and cross-correlated to find
	a set of high-scoring fits.
	The views must be calculated externally to allow for custom sets.
**/
double		img_fcc_for_view(Bimage* p, Bimage* ptemp, View view,
				FSI_Kernel* kernel, double hires, double lores, Bimage* pmask)
{
	Matrix3			mat = view.matrix();
//	mat = mat.transpose();
	
	long			i, n(0);
	double			minrad = (lores>0)? ptemp->real_size()[0]/lores: 0;
	double			maxrad = (hires>0)? ptemp->real_size()[0]/hires: ptemp->sizeX()/2;
	if ( maxrad > ptemp->sizeX()/2 ) maxrad = ptemp->sizeX()/2;

//	cout << "minrad=" << minrad << " maxrad=" << maxrad << endl;
	
	Complex<double>	v1, v2, cv;
	Vector3<double>	c, cr;
	Vector3<long>	min = Vector3<int>((int) (-maxrad - 0.5), (int) (-maxrad - 0.5), (int) (-maxrad - 0.5));
	Vector3<long>	max = Vector3<int>((int) (maxrad + 0.5), (int) (maxrad + 0.5), (int) (maxrad + 0.5));
	Vector3<long>	ci;
	double			d, pwr1(0), pwr2(0), cc(0);
	
//	cout << "min=" << min << " max=" << max << endl;
	
	for ( c[2]=min[2]; c[2]<=max[2]; c[2]+=1 ) {
		for ( c[1]=min[1]; c[1]<=max[1]; c[1]+=1 ) {
			for ( c[0]=min[0]; c[0]<=max[0]; c[0]+=1 ) {
				d = c.length();
				if ( d >= minrad && d <= maxrad ) {
					ci = c;
					if ( ci[0] < 0 ) ci[0] += p->sizeX();
					if ( ci[1] < 0 ) ci[1] += p->sizeY();
					if ( ci[2] < 0 ) ci[2] += p->sizeZ();
					i = (ci[2]*p->sizeY() + ci[1])*p->sizeX() + ci[0];
					if ( !pmask || ( pmask && (*pmask)[i] ) ) {
						v1 = p->complex(i);
						cr = mat * c;
						v2 = ptemp->fspace_interpolate(0, cr, kernel);
						pwr1 += v1.power();
						pwr2 += v2.power();
						cv += v1*v2.conj();
						n++;
					}
				}
			}
		}
	}
	
	d = pwr1*pwr2;
	if ( d > 0 ) cc = cv.amp()/sqrt(d);

	if ( verbose & VERB_FULL )
		cout << view << tab << n << tab << cc << endl;
	
	return cc;
}

/*
@brief Searches a 2D/3D density map for a template using a specific view.
@param 	*p				the Fourier transformed image.
@param 	*ptemp			the Fourier transformed template to be searched for.
@param 	view			view.
@param 	*cc				interpolation kernel.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param 	searchrad		correlation map search radius.
@param 	*pmask			mask for cross-correlation (ignored if NULL).
@return double			correlation coeficient.

	The template is rotated to the view and cross-correlated to find
	a set of high-scoring fits.
	The views must be calculated externally to allow for custom sets.
	The shift is encoded in the image origin.
**/
double		img_cc_for_view(Bimage* p, Bimage* ptemp, View view,
				FSI_Kernel* kernel, double hires, double lores, double searchrad,
				Bimage* pmask)
{
	Matrix3			mat = view.matrix();
//	mat = mat.transpose();
	
	long			i;
	double			minrad = (lores>0)? ptemp->real_size()[0]/lores: 0;
	double			maxrad = ptemp->real_size()[0]/hires;
	if ( maxrad > ptemp->sizeX()/2 ) maxrad = ptemp->sizeX()/2;
	
	Complex<double>	v1, v2;
	Vector3<double>	c, cr;
	Vector3<long>	min = Vector3<int>((int) (-maxrad - 0.5), (int) (-maxrad - 0.5), (int) (-maxrad - 0.5));
	Vector3<long>	max = Vector3<int>((int) (maxrad + 0.5), (int) (maxrad + 0.5), (int) (maxrad + 0.5));
	Vector3<long>	ci;
	double			d, pwr1(0), pwr2(0);
	
	Bimage*			pcc = new Bimage(Float, TComplex, ptemp->sizeX(),
						ptemp->sizeY(), ptemp->sizeZ(), 1);
	
	for ( c[2]=min[2]; c[2]<=max[2]; c[2]+=1 ) {
		for ( c[1]=min[1]; c[1]<=max[1]; c[1]+=1 ) {
			for ( c[0]=min[0]; c[0]<=max[0]; c[0]+=1 ) {
				d = c.length();
				if ( d >= minrad && d <= maxrad ) {
					ci = c;
					if ( ci[0] < 0 ) ci[0] += p->sizeX();
					if ( ci[1] < 0 ) ci[1] += p->sizeY();
					if ( ci[2] < 0 ) ci[2] += p->sizeZ();
					i = (ci[2]*p->sizeY() + ci[1])*p->sizeX() + ci[0];
					if ( !pmask || ( pmask && (*pmask)[i] ) ) {
						v1 = p->complex(i);
						cr = mat * c;
						v2 = ptemp->fspace_interpolate(0, cr, kernel);
						pwr1 += v1.power();
						pwr2 += v2.power();
						pcc->set(i, v1*v2.conj());
					}
				}
			}
		}
	}
	
	double			scale = 1/sqrt(pwr1*pwr2);
	
	for ( i=0; i<pcc->data_size(); i++ ) pcc->set(i, (*pcc)[i] * scale);

	pcc->fft(FFTW_BACKWARD, 0);
	
	pcc->complex_to_real();

	pcc->find_peak(searchrad, 0);
	double			cc = pcc->image->FOM();

	pcc->refine_peak();
	
//	shift = mat * shift;
	
	if ( verbose & VERB_FULL )
		cout << view << tab << pcc->image->origin() << tab << cc << endl;

	p->origin(ptemp->image->origin() + pcc->image->origin());
	
	delete pcc;
	
	return cc;
}



