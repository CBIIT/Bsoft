/**
@file	bfind.cpp
@brief	Searches orientation space for the best fit of a 3D map to a template.
@author Bernard Heymann
@date	Created: 20021027
@date	Modified: 20200401
**/

#include "rwimg.h"
#include "mg_processing.h"
#include "mg_particle_select.h"
#include "rwmg.h"
#include "symmetry.h"
#include "ps_views.h"
#include "linked_list.h"
#include "file_util.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"

#include <sys/stat.h>
#include <fcntl.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
vector<Bimage*>		setup_images(Bstring& map_file, Bstring& template_file, 
				Bstring& rs_mask, Bstring& fs_mask,
				int read_img_num, Vector3<long>& bin, Vector3<double> sam,
				int set_origin, Vector3<double>& origin, Vector3<double>& part_ori);

// Usage assistance
const char* use[] = {
" ",
"Usage: bfind [options] input.img/input.star output.img",
"------------------------------------------------------",
"Searches orientation space for the best fit of a 3D map to a template.",
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
"Specify the output of the transformed map or template with the",
"	-particle and -newtemplate options.",
" ",
"Actions:",
"-global 2                Global search before other modes with binning.",
"-mode refine             Local modes: directional, refine or symmetric.",
"-setsequential           Set particle selection numbers sequentially.",
"-subset 25,5             View subset selection from a start offset and with a given size.",
"-rescale -0.1,5.2        Rescale output particle to average and standard deviation.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; default from input file; a single value can be given).",
"-origin 110,50,44        Template image origin (default taken from template file).",
"-partorigin 50,80,47     Particle image origin (default taken from particle file).",
"-angles 8.8,2.5,5.2      Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-side 15                 Generate side view projections within the given angle from the equator.",
"-resolution 15.6,200     Resolution limits for cross-correlation (default 0).",
"-shiftlimit 5.2          Radial limit to search for shift in cross-correlation map (default 1/4 of map).",
"-recid A435              The 3D reconstruction ID (should correspond to that in the parameter file).",
"-image 2                 Image number in a multi-image file or particle number (default 1=first image).",
"-select 5                Particle selection number in parameter file (default not used).",
"-bin 3                   Binning of input, template and mask by given kernel size,",
"                         output image remains same size as input (only for local modes).",
" ",
"Parameters for the global and symmetric modes:",
"-symmetry D6             Symmetry: Point group identifier.",
" ",
"Paramweters for the directional and refinement modes:",
"-View 0.3,-0.5,0.8,33    View to initiate search: vector {xyz} and angle.",
" ",
"Parameters for the directional mode:",
"-limit 17.5              Radial limit around given view (default 20 degrees).",
" ",
"Parameters for the refinement mode:",
"-accuracy 0.3            Smallest angle step size (default 0.5 degrees).",
" ",
"Input:",
"-Template image.map      Template to search for.",
"-mask rs_mask.mrc        Mask in real space to apply before cross-correlation.",
"-Mask fs_mask.mrc        Mask in reciprocal space to apply during cross-correlation.",
" ",
"Output:",
"-Postscript plot.ps      Output postscript file with a plot of global view vectors.",
"-output file.star        Parameter output file with best orientation.",
"-ppx                     Write temporary particle parameter files to directory \"ppx\".",
"-json file.json          JSON single particle output file with best orientation.",
"-particle newpart.pif    New transformed particle file.",
"-newtemplate newtemp.pif New transformed template file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	Vector3<long>	global_bin;				// Flag and binning for global search
	Vector3<long>	bin(1,1,1);				// No binning for local modes
	Bstring			mode;					// Mode: directional, refine or symmetric
	double			nuavg(0), nustd(0);		// Rescaling to average and stdev
	int				view_start(0);			// Subset start
	int				view_subset(0);			// Subset size
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Bstring			rec_id;					// Reconstruction ID
	int				img_num(1);				// Image number in file
	int				sel_num(-1);			// Selection number in parameter file
	int				set_seq_sel(0);			// Flag to set selection number sequentially
	Vector3<double>	sam;    				// Units for the three axes (A/pixel)
	Vector3<double>	temp_ori;				// Template origin
	Vector3<double>	part_ori;				// Particle origin
	int				set_origin(0);			// Flag to set the template origin
	double			alpha_step(M_PI_4);		// Angular step size for alpha
	double			theta_step(M_PI_4);		// Angular step size for theta
	double			phi_step(M_PI_4);		// Angular step size for phi
	double			side_ang(-1);			// Side view variation angle
	double			angle_limit(M_PI/9.0);	// Angular limit for directional search
	Bsymmetry		sym;				// Default: asymmetric or C1 point group
	double			hires(0), lores(0);		// Limiting resolution for cross-correlation
	double			search_radius(-1);		// Use default search radius
	int				currview_set(0);		// Flag to indicate current view is set
	View			currview;				// View to initiate search from and use as current view
	View			bestview;				// Eventually the best view
	View			ref_view;				// Reference view
	double			accuracy(0);			// Accuracy for refinement: default at Nyquest
	int				flags(0);				// Flags for processing options
	Bstring			template_file;			// Reference template file name
	Bstring			rs_mask;				// Real space mask file name
	Bstring			fs_mask;				// Reciprocal space mask file name
	Bstring			ps_file;				// Postscript output file name
	Bstring			param_file;				// Parameter file name
	Bstring			jsfile;					// JSON file for one particle
	Bstring			newpart_file;			// Transformed particle output file name
	Bstring			newtemp_file;			// Transformed template output file name
	
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "global" ) {
			if ( ( global_bin[0] = curropt->value.integer() ) < 1 )
				cerr << "-global: A bin size must be specified!" << endl;
			else {
				if ( global_bin[0] < 1 ) global_bin[0] = 1;
				global_bin[2] = global_bin[1] = global_bin[0];
			}
		}
		if ( curropt->tag == "mode" ) {
			if ( curropt->value.length() > 0 ) {
				if ( curropt->value[0] == 'd' || curropt->value[0] == 'D' ) mode = "directional";
				if ( curropt->value[0] == 'r' || curropt->value[0] == 'R' ) mode = "refine";
				if ( curropt->value[0] == 's' || curropt->value[0] == 'S' ) mode = "symmetric";
			} else
				cerr << "-mode: A mode must be specified" << endl;
		}
		if ( curropt->tag == "setsequential" ) set_seq_sel = 1;
		if ( curropt->tag == "subset" )
			if ( curropt->values(view_start, view_subset) < 2 )
				cerr << "-subset: A view subset start and size must be specified" << endl;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
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
				set_origin |= 2;
			} else {
				temp_ori = curropt->origin();
				set_origin |= 1;
			}
		}
		if ( curropt->tag == "partorigin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin |= 8;
			} else {
				part_ori = curropt->origin();
				set_origin |= 4;
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
		if ( curropt->tag == "side" ) {
			if ( ( side_ang = curropt->value.real() ) < 0.001 )
				cerr << "-side: One angle must be specified!" << endl;
			else
				side_ang *= M_PI/180.0;
		}
		if ( curropt->tag == "limit" ) {
			if ( ( angle_limit = curropt->value.real() ) < 0.001 )
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
			if ( ( search_radius = curropt->value.real() ) < 0.001 )
				cerr << "-shiftlimit: A shift search radius must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			sym = curropt->symmetry();
		if ( curropt->tag == "bin" ) {
			if ( ( bin[0] = curropt->value.integer() ) < 1 )
				cerr << "-bin: A bin size must be specified!" << endl;
			else {
				if ( bin[0] < 1 ) bin[0] = 1;
				bin[2] = bin[1] = bin[0];
			}
		}
		if ( curropt->tag == "accuracy" ) {
			if ( ( accuracy = curropt->value.real() ) < 0.001 )
				cerr << "-accuracy: An accuracy in degrees must be specified!" << endl;
			else
				accuracy *= M_PI/180.0;
		}
		if ( curropt->tag == "ppx" ) flags |= WRITE_PPX | CHECK_PPX;
		if ( curropt->tag == "recid" )
			rec_id = curropt->value;
		if ( curropt->tag == "Template" )
			template_file = curropt->filename();
		if ( curropt->tag == "mask" )
			rs_mask = curropt->filename();
		if ( curropt->tag == "Mask" )
			fs_mask = curropt->filename();
		if ( curropt->tag == "Postscript" )
			ps_file = curropt->filename();
		if ( curropt->tag == "output" )
			param_file = curropt->filename();
		if ( curropt->tag == "json" )
			jsfile = curropt->filename();
		if ( curropt->tag == "particle" )
			newpart_file = curropt->filename();
		if ( curropt->tag == "newtemplate" )
			newtemp_file = curropt->filename();
    }
	option_kill(option);
	
	double				ti = timer_start();

	int					read_img_num = img_num - 1;
	int					cnt(1), found(0);
	Bimage*				p = NULL;
	Bimage*				ptemp = NULL;
	Bstring				filename;
	Vector3<double>		scale(1,1,1);
	Vector3<long>		size;
	Bproject*			project = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	// Sets up the input image filename
	// If the input file is a parameter file, find the reconstruction or particle to process
	if ( optind < argc ) {
		filename = argv[optind++];
		if ( file_type(filename) == Micrograph ) {
			project = read_project(filename);
			project->select = 1;
			part_set_FOM(project, 0, 0);
			if ( set_seq_sel ) part_set_sequential(project);
//				for ( rec = project->rec; rec; rec = rec->next )
//					for ( part = rec->part; part; part = part->next ) part->sel = cnt++;
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
					if ( verbose ) {
						cout << "Processing particle " << part->id << ":" << endl;
						cout << "Origin:                         " << part->ori << endl;
						cout << "View:                           " << part->view << endl;
//						cout << "Current view:                   " << currview << endl;
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

	FOMType 	fom_tag[NFOM] = {FOM};
	if ( part && ( flags & CHECK_PPX ) ) {
		// Directory for individual particle parameter files
		if ( flags & WRITE_PPX ) {
			mkdir("ppx", O_CREAT );
			chmod("ppx", 0755);
		}
		if ( ppx_check(part, fom_tag) ) {
			cerr << filename << " already done! Aborting" << endl;
			bexit(0);
		}
	}

//	cout << "part_ori=" << part_ori << endl;
	
	View*			view = NULL;
	Vector3<double>	origin, currshift, bestshift;
	double			cc, best_cc(-1e37);
	
	Matrix3			mat = ref_view.matrix();
	sym.transform(mat);
//	for ( i=0; i<sym.operations(); i++ )
//		sym.op[i].axis = mat * sym.op[i].axis;
	
	bestview = currview;

	if ( global_bin[0] > 0 || !currview_set ) {
		
		if ( side_ang < 0 )
			view = asymmetric_unit_views(sym, theta_step, phi_step, alpha_step, 1);
		else
			view = side_views(sym, side_ang, theta_step, phi_step, alpha_step);
		
		if ( verbose & VERB_FULL )
			show_views(view);

		if ( view_subset > 0 ) if ( view_list_subset(&view, view_start, view_subset) < 1 ) {
			cerr << "Error: At least one view must be selected!" << endl;
			bexit(-1);
		}

		if ( ps_file.length() && view ) ps_views(ps_file, sym.label(), view, 0);
	
		if ( verbose )
			cout << "Number of global views:         " << count_list((char *) view) << endl << endl;
	
		origin = temp_ori;

		vector<Bimage*>		imgvec = setup_images(filename, template_file, rs_mask, fs_mask, 
					read_img_num, global_bin, sam, set_origin, origin, part_ori);
		
		p = imgvec[0];

		best_cc = p->search_views(imgvec[1], view, hires, lores, 
			search_radius/global_bin[0], imgvec[2], bestview, bestshift);
		
		kill_list((char *) view, sizeof(View));
		view = NULL;
		
		if ( global_bin[0] != bin[0] )
			bestshift *= Vector3<double>(global_bin[0]/bin[0], global_bin[1]/bin[1], global_bin[2]/bin[2]);

		p->meta_data_update();
		
		if ( jsfile.length() ) p->meta_data().write(jsfile.c_str());

		for ( auto it = imgvec.begin(); it != imgvec.end(); ++it ) delete *it;
	}

	if ( mode.length() ) {
		if ( bin[0] > 1 ) search_radius /= bin[0];
		origin = temp_ori;

		vector<Bimage*>		imgvec = setup_images(filename, template_file, rs_mask, fs_mask, 
					read_img_num, bin, sam, set_origin, origin, part_ori);

		p = imgvec[0];

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

		if ( mode == "directional" || mode == "symmetric" ) {
			best_cc = p->search_views(imgvec[1], view, hires, lores, 
				search_radius, imgvec[2], bestview, bestshift);
			kill_list((char *) view, sizeof(View));
		} else if ( mode == "refine" ) {
			currview = bestview;
			best_cc = -1;
			cnt = 0;
			if ( accuracy < 2.0/p->sizeX() ) accuracy = 2.0/p->sizeX();
			while ( alpha_step >= accuracy ) {
				if ( hires > alpha_step * p->real_size()[0] )
					 hires = alpha_step * p->real_size()[0];
				view = views_for_refinement(bestview, alpha_step);
				cc = p->search_views(imgvec[1], view, hires, lores, 
					search_radius, imgvec[2], bestview, currshift);
				if ( currview.residual(bestview) < 1e-6 ||
					fabs(cc - best_cc) < 1e-6 ||
					cnt > 3 ) {
						alpha_step /= 2;	// Contract around best point
						cnt = 0;
				}
				if ( best_cc < cc ) {
					best_cc = cc;
					bestview = currview;
					bestshift = currshift;
				}
				kill_list((char *) view, sizeof(View));
				cnt++;
			}
		}
		
		p->meta_data_update();
		
		if ( jsfile.length() ) p->meta_data().write(jsfile.c_str());

		for ( auto it = imgvec.begin(); it != imgvec.end(); ++it ) delete *it;
	}

	if ( bin[0] ) bestshift *= Vector3<double>(bin[0], bin[1], bin[2]);
		
	if ( newpart_file.length() < 1 && optind < argc ) newpart_file = argv[optind];

	currshift = -bestshift;

	if ( newtemp_file.length() ) {  // Output transformed template
		mat = bestview.matrix();
		ptemp = read_img(template_file, 1, 0);
		if ( set_origin & 1 ) ptemp->origin(temp_ori);
		else if ( set_origin & 2) ptemp->origin(ptemp->size()/2);
		mat = mat.transpose();
		ptemp->calculate_background();
		ptemp->transform(scale, ptemp->image->origin(), currshift, mat, FILL_BACKGROUND, 0);
		ptemp->view(bestview);
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();	// Preserve the old type
		ptemp->change_type(nudatatype);
		write_img(newtemp_file, ptemp, 0);
		delete ptemp;
	}
	
	Bimage*		ps = NULL;
	
	if ( newpart_file.length() ) {	// Output transformed image
//		cout << "origin=" << origin << " shift=" << bestshift << " view=" << bestview << endl;
		mat = bestview.matrix();
		p = read_img(filename, 1, read_img_num);
		ptemp = read_img(template_file, 0, 0);
		if ( p->size().distance(ptemp->size()) > 0.5 || p->real_size().distance(ptemp->real_size()) > 0.1 ) {
			ps = p->scale_to_same_size(ptemp);
			delete p;
			p = ps;
		}
		delete ptemp;
		origin += currshift;	// Origin for rotation
		p->calculate_background();
//		cout << "origin=" << origin << " shift=" << bestshift << " view=" << bestview << endl;
		p->transform(scale, origin, bestshift, mat, FILL_BACKGROUND, 0);
		p->image->view(0,0,1,0);
		if ( nudatatype == Unknown_Type ) nudatatype = p->data_type();	// Preserve the old type
		p->change_type(nudatatype);
		if ( nustd > 0 ) p->rescale_to_avg_std(nuavg, nustd);
		write_img(newpart_file, p, 0);
		delete p;
	}
	
	p = read_img(filename, 0, 0);
	if ( !project ) project = new Bproject;
	if ( !rec ) {
		if ( rec_id.length() < 1 ) {
			rec_id = p->file_name();
			rec_id = rec_id.pre_rev('.');
		}
		rec = reconstruction_add(&project->rec, rec_id);
		part = particle_add(&rec->part, img_num);
	}
	if ( sam[0] > 0 ) rec->voxel_size = sam;
	else if ( p ) rec->voxel_size = p->sampling(0);
	delete p;
	rec->symmetry = sym.label();
	if ( part ) {
		if ( rec->fpart.length() < 1 ) part->fpart = filename;
		if ( newpart_file.length() ) {
			part->fpart = newpart_file;
		} else {
			part->view = bestview;
			part->ori = origin + currshift;
		}
		part->fom[0] = best_cc;
		if ( flags & WRITE_PPX ) {
			Bstring		log_name = ppx_filename(rec->id, part->id);
			write_particle(log_name, part, 0, 0, fom_tag);
		}
//		for ( part = rec->part; part; part = part->next )
//			if ( part->id == img_num ) part->fom[0] = best_cc;
//			else part->fom[0] = 0;
	} else {
		rec->frec = filename;
		if ( newpart_file.length() ) {
			rec->frec = newpart_file;
		} else {
			rec->view = bestview;
			rec->origin = origin + currshift;
		}
		for ( rec = project->rec; rec; rec = rec->next )
			if ( rec->id == rec_id ) rec->fom = best_cc;
			else rec->fom = 0;
	}
	
	if ( param_file.length() ) {
		write_project(param_file, project, 0, 0);
	}

	project_kill(project);
		
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}


vector<Bimage*>		setup_images(Bstring& map_file, Bstring& template_file, 
				Bstring& rs_mask, Bstring& fs_mask,
				int read_img_num, Vector3<long>& bin, Vector3<double> sam,
				int set_origin, Vector3<double>& origin, Vector3<double>& part_ori)
{
	int			err(0);

	Bimage*		p = NULL;
	Bimage*		ptemp = NULL;
	Bimage*		prs_mask = NULL;
	Bimage*		pfs_mask = NULL;
	
	if ( map_file.length() < 1 ) {
		cerr << "Error: No image file name specified!" << endl;
		err--;
	} else {
		p = read_img(map_file, 1, read_img_num);
		if ( p ) {
			if ( sam.volume() > 0 ) p->sampling(sam);
			if ( set_origin & 4 ) p->origin(part_ori);
			else if ( set_origin & 8 ) p->origin(p->size()/2);
			part_ori = p->image->origin();
		} else {
			cerr << "Error: " << map_file << " not read!" << endl;
			err--;
		}
	}

	if ( template_file.length() < 1 ) {
		cerr << "Error: No input template file specified!" << endl;
		err--;
	} else {
		ptemp = read_img(template_file, 1, 0);
		if ( ptemp ) {
			if ( sam.volume() > 0 ) ptemp->sampling(sam);
			if ( set_origin & 1 ) ptemp->origin(origin);
			else if ( set_origin & 2 ) ptemp->origin(ptemp->size()/2);
			origin = ptemp->image->origin();
		} else {
			cerr << "Error: " << template_file << " not read!" << endl;
			err--;
		}
	}

	if ( rs_mask.length() ) {
		prs_mask = read_img(rs_mask, 1, 0);
		if ( !prs_mask ) {
			cerr << "Error: " << rs_mask << " not read!" << endl;
			err--;
		}
		ptemp->multiply(prs_mask);
		delete prs_mask;
	}
	
	if ( fs_mask.length() ) {
		pfs_mask = read_img(fs_mask, 1, 0);
		if ( !pfs_mask ) {
			cerr << "Error: " << fs_mask << " not read!" << endl;
			err--;
		}
	}
	
	if ( err ) {
		cerr << "Too many errors! (" << -err << ")" << endl;
		bexit(-1);
	}
	
	Bimage*		ps = NULL;
	
	if ( p->size().distance(ptemp->size()) > 0.5 && p->real_size().distance(ptemp->real_size()) > 0.1 ) {
		ps = p->scale_to_same_size(ptemp);
		delete p;
		p = ps;
	}

	if ( bin[0] > 1 ) {
		if ( p ) p->bin(bin);
		if ( ptemp ) ptemp->bin(bin);
		if ( pfs_mask ) pfs_mask->bin(bin);
	} else {
		bin = Vector3<long>(1,1,1);
	}
	
	if ( verbose ) {
		if ( p )
			cout << "Map file:                       " << p->file_name() << " (" << read_img_num << ")" << endl;
		if ( ptemp )
			cout << "Template file:                  " << ptemp->file_name() << endl;
		if ( pfs_mask )
			cout << "Mask file:                      " << pfs_mask->file_name() << endl;
		cout << "Size:                           " << p->size() << endl;
		cout << endl;
	}
	
	vector<Bimage*>		imgvec;
	imgvec.push_back(p);
	imgvec.push_back(ptemp);
	imgvec.push_back(pfs_mask);
	
	return imgvec;
}

