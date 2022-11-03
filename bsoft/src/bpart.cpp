/**
@file	bpart.cpp
@brief	Process single particle images in various ways.
@author Bernard Heymann
@date	Created: 20080424
@date	Modified: 20220531
**/

#include "mg_processing.h"
#include "mg_particles.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "rwmg.h"
#include "mg_particle_select.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpart [options] input.star [input.star]",
"----------------------------------------------",
"Manipulates single particle images.",
" ",
"Actions:",
"-reconstructions         Operate on reconstruction parameters rather than micrographs.",
"-all                     Reset selection to all particles before other selections.",
"-noaberration odd        Delete aberration weights (all, odd or even).",
"-find 20                 Finds the centers of the particle images (iterations).",
"-shift                   Center the actual particle images (new images with \"_cen\" in names).",
"-align                   Output aligned particle images (new images with \"_aln\" in names).",
"-rescale -0.1,5.2        Rescale particle images to average and standard deviation.",
"-settilt -30,45          Set micrograph tilt, axis and particle defocus.",
"-calctilt                Calculate micrograph tilt parameters from particle defocus.",
"-phasedifference p.grd   Phase difference from reference: Output \"even\" and \"odd\" files.",
"-fit 2.3,2,10000         Fit phase image to resolution limit; flag: 0=full, 1=odd, 2=even; fitting iterations.",
//"-ewald ew.grd            Ewald phase map from particles.",
"-ewald ew.grd            Ewald correlation map from particles.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-select 14               Selection number of particles to process (default all selected).",
"-resolution 10,500       Resolution limits for cross-correlation (angstrom, default 20,1000).",
"-partpath dir/subdir     Set the particle file paths.",
" ",
"Input:",
"-reference map.pif       Reference for cross-correlation.",
"-mask mask.pif           Input 3D mask file to mask particles.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-composite file.pif      Composite image file.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
	int				use_rec(0);					// Flag to process reconstructions
	int 			reset(0);					// Keep selection as read from file
	int				noab(0);					// Flag to delete aberration weights
	int 			find(0);					// Iterations to find the centers
	int				shift(0);					// Flag to shift the actual images
	int				align(0);					// Flag to align images
	double			nuavg(0), nustd(0); 		// Rescaling to average and stdev
	int				calc_tilt(0);				// Flag to calculate micrograph tilt parameters
	double			axis(0), tilt(0);			// Tilt axis and angle to set
	int				part_select(-1);			// Process all selected particles
	double			hires(0);					// High resolution limit
	double			lores(1000);				// Low resolution limit
	int				fit_flag(0);				// Fit type: 0=full, 1=odd, 2=even
	double			fit_res(0);					// Fit resolution limit
	long			fit_iter(0);				// Fit iterations
	Bstring			partpath;					// Particle file path
	Bstring			maskfile;					// Particle mask file
	Bstring			outfile;					// Output parameter file
	Bstring			reffile;					// Reference image file
	Bstring			compfile;					// Composite image file
	Bstring			phifile;					// Phase difference image file base name
	Bstring			ewaldfile;					// Ewald phase ouput file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "reconstructions" )
			use_rec = 1;
		if ( curropt->tag == "all" )
			reset = 1;
		if ( curropt->tag == "noaberration" ) {
			if ( curropt->value[0] == 'o' ) noab = 1;
			if ( curropt->value[0] == 'e' ) noab = 2;
			if ( curropt->value[0] == 'a' ) noab = 3;
		}
		if ( curropt->tag == "find" )
			if ( ( find = curropt->value.integer() ) < 1 )
				cerr << "-find: A number of iterations must be specified!" << endl;
		if ( curropt->tag == "shift" )
			shift = 1;
		if ( curropt->tag == "align" )
			align = 1;
		if ( curropt->tag == "rescale" )
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both average and standard deviation must be specified!" << endl;
		if ( curropt->tag == "calctilt" ) calc_tilt = 1;
		if ( curropt->tag == "settilt" ) {
			if ( curropt->values(tilt, axis) < 1 )
				cerr << "-settilt: A tilt angle must be specified!" << endl;
			else {
				tilt *= M_PI/180.0;
				axis *= M_PI/180.0;
			}
		}
		if ( curropt->tag == "phasedifference" )
			phifile = curropt->filename();
		if ( curropt->tag == "fit" )
    	    if ( curropt->values(fit_res, fit_flag, fit_iter) < 1 )
				cerr << "-fit: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "ewald" )
			ewaldfile = curropt->filename();
		if ( curropt->tag == "select" )
			if ( ( part_select = curropt->value.integer() ) < 0 )
				cerr << "-select: A selection number must be specified!" << endl;
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "composite" )
			compfile = curropt->filename();
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

	if ( use_rec ) project->select = 1;
	
	if ( reset ) part_reset_selection(project, 3);

	if ( noab ) project_delete_aberration(project, noab);

	Bimage*			pcomp = NULL;
	Bimage*			pmask = NULL;
	
	if ( find ) 
		pcomp = project_find_particle_centers(project, find, part_select, hires, lores);
	
	if ( align )
		project_align_particles(project, part_select, nuavg, nustd);
	else if ( shift )
		project_center_particles(project, part_select, nuavg, nustd);

	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		if ( pmask ) {
			project_mask_particles(project, pmask, partpath);
			delete pmask;
		}
	}

	if ( calc_tilt )
		project_tilt_from_particle_defocus(project);

	if ( fabs(tilt) > 0.01 )
		project_set_particle_defocus_from_tilt(project, axis, tilt);

	Bimage*			pref = NULL;
	if ( reffile.length() )
		pref = read_img(reffile, 1, 0);

	FSI_Kernel*		kernel = NULL;
	Bimage*			pphi = NULL;
	if ( pref && phifile.length() ) {
		kernel = new FSI_Kernel(8, 2);
		pphi = project_aberration_phase_difference(project, pref, hires, kernel);
		Bimage*		podd = pphi->complex_split();
		Bstring		filename = phifile.pre_rev('.') + "_even." + phifile.post_rev('.');
		write_img(filename, pphi, 0);
		filename = phifile.pre_rev('.') + "_odd." + phifile.post_rev('.');
		write_img(filename, podd, 0);
		if ( fit_res && pphi ) {
			map<string,CTFparam>	cpa = project_ctf_optics_groups(project);
			vector<map<pair<long,long>,double>> we = img_aberration_fit(pphi, 0, fit_res, 2, fit_iter);
			vector<map<pair<long,long>,double>> wo = img_aberration_fit(podd, 0, fit_res, 1, fit_iter);
			long			i(0);
			for ( auto& cp: cpa ) {
				cp.second.aberration_weights(we[i]);
				cp.second.aberration_weights(wo[i++]);
				if ( verbose )
					cout << i << tab << cp.second.aberration_weight_string() << endl;
			}
			project_update_ctf_aberration(project, cpa, 1);
			Bimage*		pf = pphi->copy();
			img_create_aberration(pf, we, 2);
			filename = phifile.pre_rev('.') + "_even_fit." + phifile.post_rev('.');
			write_img(filename, pf, 0);
			pphi->subtract(pf);
			filename = phifile.pre_rev('.') + "_even_fitdif." + phifile.post_rev('.');
			write_img(filename, pphi, 0);
			img_create_aberration(pf, wo, 1);
			filename = phifile.pre_rev('.') + "_odd_fit." + phifile.post_rev('.');
			write_img(filename, pf, 0);
			podd->subtract(pf);
			filename = phifile.pre_rev('.') + "_odd_fitdif." + phifile.post_rev('.');
			write_img(filename, podd, 0);
			delete pf;
		}
		delete pref;
		delete pphi;
		delete podd;
		delete kernel;
	}
	
	if ( pref && ewaldfile.length() ) {
		if ( !kernel ) kernel = new FSI_Kernel(8, 2);
//		Bimage*		pew = project_ewald_phase(project, pref, hires, kernel);
		Bimage*		pew = project_ewald_correlation(project, pref, hires, kernel);
		write_img(ewaldfile, pew, 0);
		delete pew;
	}

	if ( outfile.length() )
		write_project(outfile, project, 0, 0);
	
	if ( pcomp && compfile.length() )
		write_img(compfile, pcomp, 0);
	
	delete pcomp;
	project_kill(project);

	if ( verbose & VERB_TIME )
		timer_report(ti);

	bexit(0);
}

