/**
@file	bmaxpow.cpp
@brief	Determining orientations by maximum power of 3D reconstruction from single particle images.
@author	Bernard Heymann
@date	Created: 20080424
@date	Modified: 20190516
**/

#include "rwmg.h"
#include "mg_processing.h"
#include "mg_particles.h"
#include "mg_particle_select.h"
#include "mg_ctf.h"
#include "rwimg.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

#include <sys/stat.h>
#include <fcntl.h>


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			project_maximum_power(Bproject* project, Bsymmetry& sym,
				double theta_step, double phi_step, double alpha_step,
				int part_select, double hi_res, double lo_res, double scale, 
				Vector3<long> size, Bimage* pmask, int pad_factor, 
                int ctf_action, double wiener, int flags);

// Usage assistance
const char* use[] = {
" ",
"Usage: bmaxpow [options] input.star [input.star]",
"------------------------------------------------",
"Finds orientations of particle images based on maximizing reciprocal space reconstruction power.",
" ",
"Actions:",
"-all                     Reset selection to all particles before other selections (default not).",
"-center                  Center particles before determining orientation (default not).",
"-symmetry C5             Point group symmetry to find orientations.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-select 14               Selection number of particles to process (default all selected).",
"-size 100,120,90         Size of reconstruction (default from images).",
"-angles 8.8,2.5,5.2      Step size for alpha, theta and phi, one value sets all (default 45 degrees).",
"-resolution 10,500       Resolution limits for finding centers (angstrom, default 20,1000).",
"-pad 3                   Image padding factor (default 0, no padding).",
"-CTF flip                Apply CTF correction to images before reconstruction (default not).",
"-wiener 0.15             Wiener factor for CTF correction (default 0.2).",
" ",
"Input:",
"-reference file.pif      Centering reference file (turns on centering).",
"-mask mask.mrc           Real space 2D mask to be applied to particles.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-ppx                     Write temporary particle parameter files to directory \"ppx\".",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int 			reset(0);				// Keep selection as read from file
	int				center(0);				// Flag to turn centering on
	int				part_select(-1);		// Process all selected particles
	Vector3<long>	map_size;				// Size of reconstruction
	int				pad_factor(0);			// Image padding factor
	int				ctf_action(0);			// Default no CTF operation
	double			wiener(0.2);			// Wiener for CTF correction
	double			hires(20);				// High resolution limit
	double			lores(1000);			// Low resolution limit
	double			theta_step(M_PI_4);		// Theta angular step size.
	double			phi_step(M_PI_4);		// Phi angular step size.
	double			alpha_step(M_PI_4);		// Alpha angular step size (in-plane).
	Bstring			symmetry_string;		// No symmetry specified
	Bstring			reffile;				// Reference image for centering
	Bstring			maskfile;				// Mask to be applied to particles
	Bstring			outfile;				// Output parameter file
	int				flags(FULL_ASU);		// Flags for processing options

	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" )
			reset = 1;
		if ( curropt->tag == "center" )
			center = 1;
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 1 )
				cerr << "-select: A selection number must be specified!" << endl;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "size" )
			map_size = curropt->size();
		if ( curropt->tag == "angles" ) {
			if ( ( i = curropt->values(alpha_step, theta_step, phi_step) ) < 1 )
				cerr << "-angles: An angle step size must be specified!" << endl;
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
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "pad" )
			if ( ( pad_factor = curropt->value.integer() ) < 0 )
				cerr << "-pad: An integer factor must be specified!" << endl;
		if ( curropt->tag == "CTF" )
			ctf_action = curropt->ctf_action();
		if ( curropt->tag == "wiener" ) {
			if ( ( wiener = curropt->value.real() ) < 0.000001 )
				cerr << "-wiener: A Wiener factor must be specified!" << endl;
			else {
				if ( wiener < 0.01 ) wiener = 0.01;
//				if ( wiener > 1 ) wiener = 1;
			}
		}
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "ppx" ) flags |= WRITE_PPX | CHECK_PPX;
    }
	option_kill(option);

	double		ti = timer_start();

	if ( hires > lores ) swap(hires, lores);

	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);
	
	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}

	if ( reset ) part_reset_selection(project, 3);
	
	Bimage*			pref = NULL;
	int				max_iter(10);
	
	if ( reffile.length() ) {
		pref = read_img(reffile, 1, 0);
		project_set_particle_centers(project, pref, part_select, hires, lores);
	} else if ( center ) {
		pref = project_find_particle_centers(project, max_iter, part_select, hires, lores);
	}
	
	if ( pref ) delete pref;
	
	Bimage*			pmask = NULL;
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, 0);

	Bsymmetry		sym(symmetry_string);
	if ( sym.point() > 101 )
		project_maximum_power(project, sym, theta_step, phi_step, alpha_step,
				part_select, hires, lores, 1, map_size, pmask, pad_factor, 
				ctf_action, wiener, flags);
	
	if ( pmask ) delete pmask;

	if ( outfile.length() )
		write_project(outfile, project, 0, 0);
	
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

double		img_reconstruct_one(Bimage* p, Vector3<long> size, Bsymmetry& sym,
				View view, double hi_res, double lo_res, double scale)
{
	long 			i, n;
	double 			fom(0);
	Vector3<double>	vscale(scale, scale, scale);
	
	Bimage*			prec = new Bimage(Float, TComplex, size, 1);
	
	prec->fspace_pack_2D(p, view, sym, hi_res, lo_res, vscale, 1, 0);

	prec->fspace_reconstruction_weigh();

	float*			power = (float *) prec->next->data_pointer();
	float* 			weight = (float *) prec->next->next->data_pointer();

	for ( i=n=0; i<prec->image_size(); i++ ) {
		if ( weight[i] > 2.5 ) {
			fom += prec->complex(i).power()/power[i];
			n++;
		}
	}
	
	if ( n ) fom /= n;
	else cout << "No overlapping voxels!" << endl;
	
	delete prec;

	return fom;
}

double		img_reconstruct_one_v3(Bimage* p, Vector3<long> size, Bsymmetry& sym,
				View view, double hi_res, double lo_res, double scale)
{
	long 			i, n;
	double 			pwr(0), avg(0), w(0), fom(0);
	Vector3<double>	vscale(scale, scale, scale);
	
	Bimage*			prec = new Bimage(Float, TComplex, size, 1);
	
	prec->fspace_pack_2D(p, view, sym, hi_res, lo_res, vscale, 1, 0);

	float*			power = (float *) prec->next->data_pointer();
	float* 			weight = (float *) prec->next->next->data_pointer();

	for ( i=n=0; i<prec->image_size(); i++ ) {
		if ( weight[i] > 1.5 ) {
			avg = prec->complex(i).power();
			pwr = power[i];
			w = pwr - avg/weight[i];
			if ( w > SMALLFLOAT ) {
				fom += (avg - pwr)/w;
				n++;
			}
		}
	}
	if ( n ) fom /= n;
	else cout << "No overlapping voxels!" << endl;
	
	delete prec;

	return fom;
}

double		img_reconstruct_one_v2(Bimage* p, Vector3<long> size, Bsymmetry& sym,
				View view, double hi_res, double lo_res, double scale)
{
	long 			i, n;
	double 			pwr(0), avg(0), w(0), fom(0);
	Vector3<double>	vscale(scale, scale, scale);
	
	Bimage*			prec = new Bimage(Float, TComplex, size, 1);
	
	prec->fspace_pack_2D(p, view, sym, hi_res, lo_res, vscale, 1, 0);

//	prec->fspace_reconstruction_weigh();

	float*			power = (float *) prec->next->data_pointer();
	float* 			weight = (float *) prec->next->next->data_pointer();

	for ( i=n=0; i<prec->image_size(); i++ ) {
		if ( weight[i] > 1.5 ) {
			avg += prec->complex(i).power();
			pwr += power[i];
			w += weight[i];
//			fom += power[i] - prec->complex(i).power();
//			n++;
		}
	}
//	if ( n ) fom /= n;
//	else cout << "No overlapping voxels!" << endl;
	if ( w ) fom = (avg - pwr)/(pwr - avg/w);
	else cout << "No overlapping voxels!" << endl;
	
//	fom = 1/(1+fom);

	delete prec;

	return fom;
}

double		img_reconstruct_one_v1(Bimage* p, Vector3<long> size, Bsymmetry& sym,
				View view, double hi_res, double lo_res, double scale)
{
	long 			i, n;
	double 			pwr(0), avg(0), w(0), fom(0);
	Vector3<double>	vscale(scale, scale, scale);
	
	Bimage*			prec = new Bimage(Float, TComplex, size, 1);
	
	prec->fspace_pack_2D(p, view, sym, hi_res, lo_res, vscale, 1, 0);

	prec->fspace_reconstruction_weigh();

	float*			power = (float *) prec->next->data_pointer();
	float* 			weight = (float *) prec->next->next->data_pointer();

	for ( i=n=0; i<prec->image_size(); i++ ) {
		if ( weight[i] > 1.5 ) {
			avg += prec->complex(i).power();
			pwr += power[i];
			w += weight[i];
//			fom += power[i] - prec->complex(i).power();
			n++;
		}
	}
//	if ( n ) fom /= n;
//	else cout << "No overlapping voxels!" << endl;
	if ( n ) fom = (n*avg - pwr)/(pwr - avg);
	else cout << "No overlapping voxels!" << endl;
	
	fom = 1/(1+fom);

	delete prec;

	return fom;
}

/**
@brief 	Determines the orientation of a particle image by maximum power.  
@param 	*p			particle image.
@param 	size		reconstruction size.
@param 	*sym		point group symmetry.
@param 	theta_step	theta angular step size.
@param 	phi_step	phi angular step size.
@param 	alpha_step	alpha angular step size (in-plane).
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit (infinite if 0).
@param 	scale		scale of reconstruction and particle magnification.
@param 	pad_factor	factor that determines image padding.
@return	double 	 	best FOM.

	For each view, a reconstruction is generated using the given symmetry.
	The reconstruction with the maximum power is selected and the
	corresponding view assigned to the particle.
		P = sum(|F|^2)

**/
double		img_find_maximum_power(Bimage* p, Vector3<long> size, Bsymmetry& sym,
				double theta_step, double phi_step, double alpha_step,
				double hi_res, double lo_res, double scale, int pad_factor)
{
	if ( !p->data_pointer() ) return 0;
	
	if ( hi_res < 2*p->sampling(0)[0]/scale ) hi_res = 2*p->sampling(0)[0]/scale;
	
	if ( verbose & VERB_FULL )
		cout << "Finding the orientation based on maximum power up to " << hi_res << " A resolution" << endl;
	
	long			iv, nviews;
	double			fomavg, bestfom;
	View*			v;
	View*			view = asymmetric_unit_views(sym, theta_step, phi_step, 1);
	View*			view2 = view_list_expand_angles(view, -M_PI, M_PI - alpha_step/2, alpha_step);
	kill_list((char *) view, sizeof(View));
	view = view2;
	View			bestview;
	
	for ( nviews=0, v=view; v; v=v->next ) nviews++;
	
//	cout << "number of views = " << nviews << endl;

	View*			view_arr = new View[nviews];
	double*			fom_arr = new double[nviews];

	for ( iv=0, v=view; v; v=v->next, iv++ ) view_arr[iv] = *v;
	
	if ( verbose & VERB_FULL )
		cout << "x\ty\tz\ta\tFOM" << endl;

#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t iv){
		fom_arr[iv] = img_reconstruct_one(p, size, sym, view_arr[iv], hi_res, lo_res, scale);
	});
#else
#pragma omp parallel for
	for ( iv = 0; iv < nviews; iv++ )
		fom_arr[iv] = img_reconstruct_one(p, size, sym, view_arr[iv], hi_res, lo_res, scale);
#endif

	for ( iv = 0, bestfom = -1e30, fomavg = 0; iv < nviews; iv++ ) {
		fomavg += fom_arr[iv];
		if ( bestfom < fom_arr[iv] ) {
			bestfom = fom_arr[iv];
			bestview = view_arr[iv];
		}
	}
	
	fomavg /= nviews;
	
	if ( verbose & VERB_FULL )
		cout << bestview << tab << fomavg << endl;
	
	delete[] view_arr;
	delete[] fom_arr;
	kill_list((char *) view, sizeof(View));
	
	p->view(bestview);
	p->image->FOM(fomavg);

	return bestfom;
}

double		part_find_maximum_power(Bparticle* part,
				Vector3<long> size, Bimage* pmask, long ft_size, Bsymmetry& sym,
				double theta_step, double phi_step, double alpha_step,
				double hi_res, double lo_res, double scale, int pad_factor, 
                int ctf_action, double wiener, int flags)
{
	FOMType 		fom_tag[NFOM] = {FOM_CC, FOM_CC_AVG};

	if ( flags & CHECK_PPX )
		if ( ppx_check(part, fom_tag) )
			return part->fom[0];
	
	double			part_scale(1);
	Bimage*			p = NULL;
	Bmicrograph*	mg = part->mg;
	
	if ( ( p = read_img(mg->fpart, 1, part->id - 1) ) == NULL ) {
		error_show("part_find_maximum_power", __FILE__, __LINE__);
		return -1;
	}

    CTFparam        em_ctf;
    if ( mg->ctf ) em_ctf.update(mg->ctf);
    if ( part->def > 0 ) em_ctf.defocus_average(part->def);

	if ( pmask ) p->multiply(pmask);

	if ( part->ori[0] <= 0 ) {
		if ( p->image->origin()[0] <= 0 ) part->ori = p->size()/2;
		else part->ori = p->image->origin();
	}
	p->image->origin(part->ori);
					
	p->view(part->view);
				
	if ( size[2] < 2 ) {
		if ( p->image->view()[2] >= 0 ) p->image->view(0,0,1,0);
		else p->image->view(0,0,-1,0);
	}
					
	if ( part->mag > 0 ) part_scale = scale/part->mag;
	else  part_scale = scale;

	p->rescale_to_avg_std(0, 1);
	p->calculate_background();
	
	p->pad(ft_size, FILL_BACKGROUND, p->background(long(0)));
	
	p->fft();
				
	p->phase_shift_to_origin();

	p->sampling(mg->pixel_size);

	if ( mg->ctf && ctf_action )
		img_ctf_apply(p, em_ctf, ctf_action, wiener, 0, hi_res);

	part->fom[0] = img_find_maximum_power(p, size, sym, theta_step, 
			phi_step, alpha_step, hi_res, lo_res, part_scale, pad_factor);

	part->view = p->image->view();
	part->fom[1] = p->image->FOM();
				
	delete p;
//	cout << "Particle done" << endl;										

	if ( flags & WRITE_PPX ) {
		Bstring		ppx_name = ppx_filename(part->mg->id, part->id);
		write_particle(ppx_name, part, 0, 0, fom_tag);
		ppx_name = 0;
	}
	
	return part->fom[0];
}

/**
@brief 	Determines the orientation of each particle from the maximum power of its reconstruction.  
@param 	*project 	image processing parameter structure.
@param 	*sym		point group symmetry.
@param 	theta_step	theta angular step size.
@param 	phi_step	phi angular step size.
@param 	alpha_step	alpha angular step size (in-plane).
@param 	part_select	selection number from the selection column.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit (infinite if 0).
@param 	scale		scale of reconstruction.
@param 	size		size of reconstruction.
@param 	*pmask		mask to eliminate unwanted parts (can be NULL).
@param 	pad_factor	factor that determines image padding.
@param 	ctf_action	type of CTF calculated (1-8).
@param 	wiener		Wiener factor (fraction).
@param 	flags		option flags.
@return	int			0.

	Each particle is transformed to a view on a grid and a reconstruction
	generated. The view associated with the reconstruction with the highest
	power is accepted.

**/
int			project_maximum_power(Bproject* project, Bsymmetry& sym,
				double theta_step, double phi_step, double alpha_step,
				int part_select, double hi_res, double lo_res, double scale, 
				Vector3<long> size, Bimage* pmask, int pad_factor, 
                int ctf_action, double wiener, int flags)
{
	if ( pad_factor < 1 ) pad_factor = 1;
	if ( pad_factor > 8 ) pad_factor = 8;
	
	Bfield*			field = project->field;
	Bparticle*		part = part_find_first(project);

	if ( !part ) {
		error_show("Error in project_maximum_power", __FILE__, __LINE__);
		cerr << "No particles found! " << endl << endl;
		return 0;
	}

	Bmicrograph*	mg = part->mg;

	if ( !mg->fpart.length() ) {
		error_show("Error in project_maximum_power", __FILE__, __LINE__);
		cerr << "No particle image file name given for micrograph " << mg->id << endl << endl;
		return 0;
	}

	Bimage* 		p = read_img(mg->fpart, 0, 0);
	p->sampling(mg->pixel_size);

	if ( size.volume() < 1 ) {
		size[0] = p->sizeX();
		if ( size[0] < p->sizeY() ) size[0] = p->sizeY();
		if ( size[0] < p->sizeZ() ) size[0] = p->sizeZ();
		size[2] = size[1] = size[0];
	}
	
	if ( size[0] > 2*p->real_size()[0]/hi_res ) 
		size[2] = size[1] = size[0] = 2*p->real_size()[0]/hi_res;
	
	long 	ft_size = findNextPowerOf((long)(pad_factor*size[0]*scale), 2);
//	if ( ft_size < p->sizeX() ) ft_size = p->sizeX();
	if ( ft_size < p->sizeX() ) ft_size = findNextPowerOf(p->sizeX(), 2);

	delete p;
	
	if ( size[0] > 500 && hi_res > 4*mg->pixel_size[0] ) size /= 2;
	
//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "Finding particle orientations by maximum power:" << endl;
		cout << "Map size:                       " << size << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Angular steps:                  " << alpha_step*180.0/M_PI
			<< " " << theta_step*180.0/M_PI << " " << phi_step*180.0/M_PI << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Resolution limits:              " << hi_res << " - " << lo_res << " A" << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "CTF application type:           " << ctf_action << endl;
		if ( pmask )
		cout << "Real space mask:                " << pmask->file_name() << endl << endl;
	}

	// Directory for individual particle parameter files
	if ( flags & WRITE_PPX ) {
		mkdir("ppx", O_CREAT );
		chmod("ppx", 0755);
	}
	
	long 			nrec(0), psel;
	long			nsel = project_count_mg_part_selected(project);

	if ( verbose )
		cout << "Part\tx\ty\tz\ta\tFOMavg\tFOMmax\tRatio" << endl;
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			for ( part=mg->part; part; part=part->next ) {
				psel = 0;
				if ( part->sel ) {
					if ( part_select < 0 ) {
						psel = 1;
					} else if ( part->sel == part_select ) {
						psel = 1;
					}
				}
				if ( psel ) {	// Use only selected images
					part_find_maximum_power(part, size, pmask, ft_size, sym,
						theta_step, phi_step, alpha_step, hi_res, lo_res, scale, 
                        pad_factor, ctf_action, wiener, flags);
					if ( verbose )
						cout << part->id << tab << setprecision(4) <<
							part->view[0] << tab <<
							part->view[1] << tab << part->view[2] << tab <<
							part->view.angle()*180.0/M_PI << tab <<
							part->fom[1] << tab << part->fom[0] << tab <<
							part->fom[0]/part->fom[1] << endl;

					nrec++;
					if ( verbose & ( VERB_TIME | VERB_PROCESS ) )
						cout << "Complete:                       " << nrec*100.0/nsel << "%\r" << flush;
				}
			}
		}
	}
	
	return 0;
}				

