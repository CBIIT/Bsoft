/**
@file	bico.cpp
@brief	Program to analyze icosahedral subtomograms
@author Bernard Heymann
@date	Created: 20180427
@date	Modified: 20201124

**/

#include "mg_processing.h"
#include "mg_reconstruct.h"
#include "rwmg.h"
#include "rwimg.h"
#include "rwsymop.h"
#include "symmetry.h"
#include "linked_list.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

fft_plan	planf, planb;

Bimage*		img_vertex_reconstruct(Bproject* project, Vector3<long> size, double maxres);
long		project_ico_worst_vertex(Bproject* project, Bimage* pref, 
				Bimage* pmask, Bimage* pfsmask, double hires, double lores);
long		project_ico_opposite_vertex(Bproject* project, Bimage* pmask, 
				Bimage* pfsmask, double hires, double lores);
long		project_ico_best_c5(Bproject* project, Bimage* pref, 
				Bimage* pmask, Bimage* pfsmask, double hires, double lores);
long		project_assign_special_vertex_orientations(Bproject* project, Bproject* project_specvert);
long		project_analyze_vertices(Bproject* project, double threshold);
long		project_analyze_vertices2(Bproject* project, double significance);
long		project_analyze_vertices3(Bproject* project, long maxcut);

// Usage assistance
const char* use[] = {
" ",
"Usage: bico [options] tomo.star",
"-------------------------------",
"Analyzes icosahedral subtomograms.",
" ",
"Actions:",
"-find opp                Type of correlation: worst, opposite, c5.",
"-analyze 0.1             Calculate statistics for vertices.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-resolution 15.6,200     Resolution limits for cross-correlation (default 0,inf).",
"-maxresolution 9.5       Maximum reconstruction resolution limit (default Nyquist).",
" ",
"Input:",
"-reference image.map     Reference for cross-correlation.",
"-mask image.map          Real space mask.",
"-Mask fsmask.mrc         Frequency space mask (missing data).",
"-special vert.star       Parameter file with special vertex orientations.",
" ",
"Output:",
"-output file.star        Output parameter file.",
"-reconstruction file.ext Reconstruction file name.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	// Initialize variables
	Bstring			type;				// Type of analysis: worst or opposite
	int				analyze(0);			// Flag to calculate statistics
	double			threshold(-1);		// Threshold for finding odd vertices
	double			hires(0), lores(0);	// Limiting resolution for cross-correlation
	double			maxres(0);			// Maximum reconstruction resolution
	Bstring			ref_file;			// Reference template file name
	Bstring			mask_file;			// Real space mask
	Bstring			fsmask_file;		// Frequency space missing mask
	Bstring			specvert_file;		// Special vertex parameter file
	Bstring			outfile;			// Output parameter file name
	Bstring			reconsfile;			// Output reconstruction file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "find" )
			type = curropt->value;
		if ( curropt->tag == "analyze" ) {
			analyze = 1;
			if ( ( threshold = curropt->value.real() ) < -0.99 )
				cerr << "-analyze: A threshold must be specified!" << endl;
		}
		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified!" << endl;
		if ( curropt->tag == "maxresolution" )
			if ( ( maxres = curropt->value.real() ) < 0.1 )
				cerr << "-maxresolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "reference" )
			ref_file = curropt->filename();
		if ( curropt->tag == "mask" )
			mask_file = curropt->filename();
		if ( curropt->tag == "Mask" )
			fsmask_file = curropt->filename();
		if ( curropt->tag == "special" )
			specvert_file = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
 		if ( curropt->tag == "reconstruction" )
			reconsfile = curropt->filename();
   }
	option_kill(option);	
	
	double		ti = timer_start();

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

	Vector3<long> 	size;
	Bimage*			pref = NULL;
	if ( ref_file.length() ) {
		pref = read_img(ref_file, 1, 0);
	} else if ( !analyze ) {
		pref = project_reconstruct_3D(project, -1, size, maxres);
		write_img("reference.mrc", pref, 0);
	}
	
	if ( !analyze && !pref ) {
		cerr << "Error: No template map!" << endl;
		bexit(-1);
	}
	
	Bimage*			pmask = NULL;
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, 0);
		pref->multiply(pmask);
	}

	Bimage*			pfsmask = NULL;
	if ( fsmask_file.length() ) {
		pfsmask = read_img(fsmask_file, 1, 0);
		if ( pfsmask->image->origin().length() < 1 )
			pfsmask->center_wrap();		// Make sure the origin is in the center
	}

	if ( maxres > hires ) maxres = hires;
 
 	int				done(0);
	Bimage*			pavg = NULL;
	
	if ( type.length() ) {
		if ( type[0] == 'w' ) {
			done = project_ico_worst_vertex(project, pref, pmask, pfsmask, hires, lores);
		} else if ( type[0] == 'o' ) {
			done = project_ico_opposite_vertex(project, pmask, pfsmask, hires, lores);
		} else if ( type == "c5" ) {
			done = project_ico_best_c5(project, pref, pmask, pfsmask, hires, lores);
		}
		if ( done ) pavg = img_vertex_reconstruct(project, size, maxres);
	}

	if ( specvert_file.length() ) {
		Bproject*		project_specvert = read_project(specvert_file);
		project_assign_special_vertex_orientations(project, project_specvert);
		project_kill(project_specvert);
	}
	
	if ( analyze ) project_analyze_vertices3(project, threshold);

	
	if ( project && outfile.length() )
		write_project(outfile, project, 0, 0);
	
	long			i(0);
	Bstring			filename(reconsfile);
	if ( pavg && reconsfile.length() ) {
		for ( Bimage* p = pavg; p; p = p->next ) {
				if ( pavg->next ) filename = reconsfile.pre_rev('.') + Bstring(++i, "_%02d.") + reconsfile.post_rev('.');
			if ( verbose )
				cout << "Writing " << filename << endl;
			write_img(filename, p, 0);
		}
	}

	project_kill(project);
	delete pref;
	delete pavg;
	delete pmask;
	delete pfsmask;

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

Bimage*		img_read_orient_mask(Bstring& filename, Bimage* pmask, View view)
{
	Bimage*		p = read_img(filename, 1, 0);
	
	if ( !p ) {
		error_show("img_read_mask_orient", __FILE__, __LINE__);
		return p;
	}
	
	p->rotate(view);
	
	if ( pmask ) p->multiply(pmask);
	
	p->image->view(view);
	
	return p;
}

double		img_cc(Bimage* p1, Bimage* p2, Bimage* pfsmask, double hires, double lores)
{
//	double			cc(1);
	Bimage*			pfsmrot = NULL;
	if ( pfsmask ) {
		pfsmrot = pfsmask->copy();
		if ( pfsmrot->image->origin().length() < 1 )
			pfsmrot->center_wrap();
		pfsmrot->rotate(p1->image->view());
		pfsmrot->zero_origin();
	}
	
	p1->find_shift(0, p2, pfsmrot, hires, lores, 0.5, planf, planb);
/*	
	if ( cc <= 1e-10 ) {
		cerr << "Correlation coefficient too small! " << cc << endl;
		bexit(-1);
	}*/
	
	delete pfsmrot;
	
	return p1->image->FOM();
}

Bimage*		img_vertex_reconstruct(Bproject* project, Vector3<long> size, double maxres)
{
	Bimage*		prec = project_reconstruct_3D(project, 1, size, maxres);

    delete prec->next;

	prec->next = project_reconstruct_3D(project, 2, size, maxres);
	
    delete prec->next->next;
    prec->next->next = NULL;
    
    return prec;
}

long		particle_worst_vertex(Bparticle* vert, Bimage* pref, Bimage* pmask, Bimage* pfsmask, double hires, double lores)
{
	long				i, psel(0);
	double				cc_worst(1);
	Bparticle*			worst = NULL;
	
	vector<Bparticle*>	varr(12);

	for ( i=0, worst = vert; worst && i<12; worst = worst->next, ++i )
		varr[i] = worst;

#ifdef HAVE_GCD
	dispatch_apply(11, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bparticle*		v = varr[i];
		Bimage*			p = img_read_orient_mask(v->fpart, pmask, v->view);
		v->fom[0] = img_cc(p, pref, pfsmask, hires, lores);
		delete p;
	});
#else
#pragma omp parallel for
	for ( i=0; i<12; ++i ) {
		Bparticle*		v = varr[i];
		Bimage*			p = img_read_orient_mask(v->fpart, pmask, v->view);
		v->fom[0] = img_cc(p, pref, pfsmask, hires, lores);
		delete p;
	}
#endif

	for ( i=0; vert && i<12; vert = vert->next, ++i ) {
		if ( verbose & VERB_FULL )
			cout << vert->sel << tab << vert->fom[0] << endl;
		if ( cc_worst > vert->fom[0] ) {
			psel = vert->sel;
			cc_worst = vert->fom[0];
			worst = vert;
		}
		vert->sel = 1;
	}

	if ( worst )
		worst->sel = 2;
	
	return psel;
}

long		particle_cc_opposite_vertices(Bparticle* vert, Bimage* pmask, 
				Bimage* pfsmask, double hires, double lores)
{
	long			i, j, psel(0);
	double			at(M_PI - 0.02), cc, cc_worst(1);
	Bparticle*		vert1 = vert;
	Bparticle*		vert2;
	Bparticle*		worst = NULL;
	Bparticle*		worst2 = NULL;
	Bimage			*p1, *p2;

	for ( i=0, vert2 = vert1; vert2 && i<12; vert2 = vert2->next, ++i )
		vert2->fom[0] = -2;
	
	for ( i=0; vert && i<11; vert = vert->next, ++i ) if ( vert->fom[0] < -1 ) {
		p1 = img_read_orient_mask(vert->fpart, pmask, vert->view);
		for ( j=i+1, vert2 = vert->next; vert2 && j<12; vert2 = vert2->next, ++j )
			if ( vert->view.backward().angle(vert2->view.backward()) > at ) break;
		if ( j<12 && vert2 ) {
			p2 = img_read_orient_mask(vert2->fpart, pmask, vert2->view);
			cc = img_cc(p1, p2, pfsmask, hires, lores);
			delete p2;
			vert->fom[0] = vert2->fom[0] = cc;
			if ( verbose & VERB_FULL )
				cout << vert->sel << " - " << vert2->sel << tab << cc << endl;
			if ( cc_worst > cc ) {
				psel = vert->sel;
				cc_worst = cc;
				worst = vert;
				worst2 = vert2;
			}
		}
		delete p1;
	}

	for ( i=0, vert2 = vert1; vert2 && i<12; vert2 = vert2->next, ++i )
		vert2->sel = 1;
	
	if ( worst )
		worst->sel = worst2->sel = 2;
	
	return psel;
}

long		particle_worst_opposite_vertex(Bparticle* vert, Bimage* pref, 
				Bimage* pmask, Bimage* pfsmask, double hires, double lores)
{
	long		i, psel(0);
	double		cc, cc1(0);
	Bparticle*	vert1 = NULL;
	Bparticle*	vert2;
	Bimage*		p;

	for ( i=1, vert2 = vert; vert2 && i<=12; vert2 = vert2->next, ++i ) {
		if ( vert2->sel == 2 ) {
			p = img_read_orient_mask(vert2->fpart, pmask, vert2->view);
			cc = img_cc(p, pref, pfsmask, hires, lores);
			delete p;
			if ( vert1 ) {
				if ( cc1 > cc ) {
					vert1->sel = 1;
				} else {
					vert2->sel = 1;
					psel = i;
				}
				break;
			} else {
				vert1 = vert2;
				cc1 = cc;
				psel = i;
			}
		}
	}
		
	return psel;
}

long		particle_best_c5(Bparticle* part, Bimage* pref, Bimage* pmask, 
				Bimage* pfsmask, double hires, double lores)
{
	long				i, psel(0);	
	double*				cc = new double[5];
	
	Bsymmetry 			sym("C5");
	View*				views = symmetry_get_all_views(sym, part->view);
	View*				v;
	vector<View>		view_arr(5);
	for ( i=0, v=views; v; v=v->next, ++i ) view_arr[i] = *v;
	kill_list((char *) views, sizeof(View));

#ifdef HAVE_GCD
	dispatch_apply(4, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bimage*			p = img_read_orient_mask(part->fpart, pmask, view_arr[i]);
		cc[i] = img_cc(p, pref, pfsmask, hires, lores);
		delete p;
	});
#else
#pragma omp parallel for
	for ( i=0; i<5; ++i ) {
		Bimage*			p = img_read_orient_mask(part->fpart, pmask, view_arr[i]);
		cc[i] = img_cc(p, pref, pfsmask, hires, lores);
		delete p;
	}
#endif

	part->fom[0] = -1;
	for ( i=0; i<5; ++i ) {
		if ( verbose & VERB_FULL )
			cout << i << tab << cc[i] << endl;
		if ( part->fom[0] < cc[i] ) {
			psel = i+1;
			part->fom[0] = cc[i];
			part->view = view_arr[i];
		}
	}
	
	delete[] cc;
	
	return psel;
}

/**
@brief 	Determines the special vertex for icosahedral symmetry and returns two maps.
@param 	*project	parameter structure with all parameters.
@param	*pref		reference map.
@param	*pmask		real space mask.
@param	*pfsmask	frequency space mask.
@param 	hires		maximum resolution for correlation.
@param 	lores		minimum resolution for correlation.
@return long		1 if done, 0 on error.

	The function correlates each vertex with the reference, and selects the one
	out of the twelve for each particle with the worst coefficient as the 
	likely special vertex.
	Two averages are calculated, the first for the non-special vertices,
	and the second for the special vertices.
**/
long		project_ico_worst_vertex(Bproject* project, Bimage* pref, 
				Bimage* pmask, Bimage* pfsmask, double hires, double lores)
{
	if ( hires < 2*project->rec->voxel_size[0] )
		hires = 4*project->rec->voxel_size[0];
	if ( lores < hires ) lores = 1e6;
		
	if ( !pref ) {
		cerr << "Error: No template map!" << endl;
		return 0;
	}
	
	planf = pref->fft_setup(FFTW_FORWARD, 1);
	planb = pref->fft_setup(FFTW_BACKWARD, 1);
	
	long				i, psel;
	Vector3<long> 		size;
	Breconstruction*	rec;
	Bparticle*			part;
	Bparticle*			vert;
	
	if ( verbose ) {
		cout << "Determining the special vertices:" << endl;
		cout << "Resolution limits:               " << hires << " - " << lores << " A" << endl;
		if ( pref )
			cout << "Reference:                       " << pref->file_name() << endl;
		if ( pmask )
			cout << "Real space mask:                 " << pmask->file_name() << endl;
		if ( pfsmask )
			cout << "Frequency space mask:            " << pfsmask->file_name() << endl;
	}

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose ) {
			cout << "Reconstruction: " << rec->frec << endl;
			cout << "Particle\tSpecial vertex" << endl;
		}
		for ( i=1, vert = part = rec->part; part; part = part->next, ++i ) {
			if ( i == 12 ) {
				psel = particle_worst_vertex(vert, pref, pmask, pfsmask, hires, lores);
				if ( verbose )
					cout << vert->group << tab << psel << endl;
				vert = part->next;
				i = 0;
			}
		}
	}
	
//	delete pref;

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	return 1;
}

/**
@brief 	Determines the special vertex for icosahedral symmetry and returns two maps.
@param 	*project	parameter structure with all parameters.
@param	*pmask		real space mask.
@param	*pfsmask	frequency space mask.
@param 	hires		maximum resolution for correlation.
@param 	lores		minimum resolution for correlation.
@return long		1 if done, 0 on error.

	The function compares opposite vertices of the twelve for every particle
	to identify the pair that is most dissimilar. 
	One of these is a candidate for the special vertex.
	A reference map is then calculated from all the rest of the vertices.
	The previously selected pair is then correlated agiants the reference
	and the one with the worst coefficient selected as final candidate.
	Two averages are calculated, the first for the non-special vertices,
	and the second for the special vertices.
**/
long		project_ico_opposite_vertex(Bproject* project, Bimage* pmask, 
				Bimage* pfsmask, double hires, double lores)
{
	if ( hires < 2*project->rec->voxel_size[0] )
		hires = 4*project->rec->voxel_size[0];
	if ( lores < hires ) lores = 1e6;
	
	long				i, psel;
	Breconstruction*	rec;
	Bparticle*			part;
	Bparticle*			vert;

	if ( verbose ) {
		cout << "Determining the special vertices:" << endl;
		cout << "Resolution limits:               " << hires << " - " << lores << " A" << endl;
		if ( pmask )
			cout << "Real space mask:                 " << pmask->file_name() << endl;
		if ( pfsmask )
			cout << "Frequency space mask:            " << pfsmask->file_name() << endl;
	}

	planf = fft_setup_plan(project->rec->box_size, FFTW_FORWARD, 1);
	planb = fft_setup_plan(project->rec->box_size, FFTW_BACKWARD, 1);

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose ) {
			cout << "Reconstruction: " << rec->frec << endl;
			cout << "Particle\tVertex pair" << endl;
		}
		for ( i=1, vert = part = rec->part; part; part = part->next, ++i ) {
			if ( i == 12 ) {
				psel = particle_cc_opposite_vertices(vert, pmask, pfsmask, hires, lores);
				if ( verbose )
					cout << vert->group << tab << psel << endl;
				vert = part->next;
				i = 0;
			}
		}
	}

	Vector3<long> 		size;
	Bimage*				pavg = project_reconstruct_3D(project, 1, size, hires);
/*
	if ( mask_file.length() ) {
		pmask = read_img(mask_file, 1, 0);
		pavg->multiply(pmask);
	}
*/	
	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose ) {
			cout << "Reconstruction: " << rec->frec << endl;
			cout << "Particle\tSpecial vertex" << endl;
		}
		for ( i=1, vert = part = rec->part; part; part = part->next, ++i ) {
			if ( i == 12 ) {
				psel = particle_worst_opposite_vertex(vert, pavg, pmask, pfsmask, hires, lores);
				if ( verbose )
					cout << vert->group << tab << psel << endl;
				vert = part->next;
				i = 0;
			}
		}
	}

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);
    
    delete pavg;

	return 1;
}

/**
@brief 	Determines the best 5-fold orientation of the vertex wrt a reference.
@param 	*project	parameter structure with all parameters.
@param	*pref		reference map.
@param	*pmask		real space mask.
@param	*pfsmask	frequency space mask.
@param 	hires		maximum resolution for correlation.
@param 	lores		minimum resolution for correlation.
@return long		1 if done, 0 on error.

	The function compares opposite vertices of the twelve for every particle
	to identify the pair that is most dissimilar. 
	One of these is a candidate for the special vertex.
	A reference map is then calculated from all the rest of the vertices.
	The previously selected pair is then correlated agiants the reference
	and the one with the worst coefficient selected as final candidate.
	Two averages are calculated, the first for the non-special vertices,
	and the second for the special vertices.
**/
long		project_ico_best_c5(Bproject* project, Bimage* pref, 
				Bimage* pmask, Bimage* pfsmask, double hires, double lores)
{
	if ( hires < 2*project->rec->voxel_size[0] )
		hires = 4*project->rec->voxel_size[0];
	if ( lores < hires ) lores = 1e6;
	
	long				psel, n(0), nc(0);
	Breconstruction*	rec;
	Bparticle*			part;

	if ( !pref ) {
		cerr << "Error: No template map!" << endl;
		return 0;
	}

	if ( verbose ) {
		cout << "Determining the best 5-fold for each vertex:" << endl;
		cout << "Resolution limits:               " << hires << " - " << lores << " A" << endl;
		cout << "Reference:                       " << pref->file_name() << endl;
		if ( pmask )
			cout << "Real space mask:                 " << pmask->file_name() << endl;
		if ( pfsmask )
			cout << "Frequency space mask:            " << pfsmask->file_name() << endl;
	}

	planf = fft_setup_plan(project->rec->box_size, FFTW_FORWARD, 1);
	planb = fft_setup_plan(project->rec->box_size, FFTW_BACKWARD, 1);

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose ) {
			cout << "Reconstruction: " << rec->frec << endl;
			cout << "Particle\tBestC5" << endl;
		}
		for ( part = rec->part; part; part = part->next ) {
			if ( part->sel == 2 ) {
				psel = particle_best_c5(part, pref, pmask, pfsmask, hires, lores);
				if ( verbose )
					cout << part->group << tab << psel << endl;
				n++;
				if ( psel > 1 ) nc++;
			}
		}
	}
	
	if ( verbose ) 
		cout << "Number of view changes:          " << nc << " (" << nc*100.0/n << "%)" << endl;

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	return 1;
}

long		project_assign_special_vertex_orientations(Bproject* project, Bproject* project_specvert)
{
	long				i, npart(0);

	Breconstruction*	rec, *rec_sv;
	Bparticle*			part, *vert;

	if ( verbose )
		cout << "Transferring the special vertex orientations" << endl;

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		for ( rec_sv = project_specvert->rec; rec_sv; rec_sv = rec_sv->next )
			if ( rec_sv->id == rec->id && rec->frec == rec_sv->frec ) break;
		if ( rec_sv && rec_sv->part ) {
			if ( verbose )
				cout << "Reconstruction: " << rec->frec << endl;
			vert = rec_sv->part;
			for ( part = rec->part; part; part = part->next ) {
				if ( part->id == vert->group ) {
					for ( i=0; vert && i<12; vert = vert->next, ++i ) {
						cout << part->id << tab << vert->id << tab << vert->sel << endl;
						if ( vert->sel == 2 ) {
							part->view = vert->view;
							npart++;
						}
					}
				} else {
					cerr << "Error: reconstruction " << rec->id << " particle " << part->id << endl;
				}
			}
		}
	}
	
	if ( verbose )
		cout << "Number of particles updated: " << npart << endl;
	
	return npart;
}


long		project_analyze_vertices(Bproject* project, double threshold)
{
	long				i, i2, imin, npart(0), nvert(0), nodd, noddtot(0);
	long				h, nh(40), bins(nh+1);
	double				avg, max, min, min2, fommax(2), scale(nh/fommax);
	vector<double>		fom(12,0);
	vector<long>		fomhis(bins,0);

	Breconstruction*	rec;
	Bparticle*			part, *vert;

	if ( verbose )
		cout << "Analyzing vertices:" << endl;

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose )
			cout << "Reconstruction: " << rec->id << endl;
		for ( part = vert = rec->part; part; part = vert, npart++ ) {
			avg = 0;
			min = min2 = 1;
			max = -1;
			for ( vert = part, i = i2 = 0; vert && vert->group == part->group;
						vert = vert->next, ++i ) {
//				cout << vert->group << tab << vert->id << tab << vert->fom[0] << endl;
				fom[i] = vert->fom[0];
				if ( min > vert->fom[0] ) {
					min = vert->fom[0];
					imin = i;
				}
				if ( max < vert->fom[0] ) max = vert->fom[0];
				if ( vert->sel == 1 ) {
					avg += vert->fom[0];
				} else {
					i2 = i;
					min2 = vert->fom[0];
				}
			}
			avg /= 11;
			for ( vert = part, i = nodd = 0; vert && vert->group == part->group;
						vert = vert->next, ++i ) {
//				fom[i] = (fom[i] - min)/(avg - min);
//				fom[i] = fom[i]/max;
				fom[i] = (fom[i] - min)/avg;
				if ( fom[i] <= threshold ) nodd++;
				h = long(scale*fom[i] + 0.5);
				if ( h < 0 ) h = 0;
				if ( h > nh ) h = nh;
				fomhis[h]++;
			}
			if ( verbose )
				cout << part->group << tab << avg << tab << i2 << tab << imin << tab << nodd << endl;
			nvert += i;
			noddtot += nodd;
			if ( i < 12 )
				cerr << "Error: Too few vertices! (" << i << ")" << endl;
		}
	}
	
	if ( verbose ) {
		cout << "Number of particles:     " << npart << endl;
		cout << "Number of odd vertices:  " << noddtot << endl;
	}
	
	cout << "Bin\tFOM\tCount" << endl;
	for ( h=0; h<bins; ++h )
		cout << h << tab << h/scale << tab << fomhis[h] << endl;
	
	return npart;
}


long		project_analyze_vertices2(Bproject* project, double significance)
{
	long				i, npart(0), nvert(0), non, nop, nontot(0), noptot(0);
	double				avg, sig, t;
	double				avgall(0), sigall(0);
//	double				alpha(1.796);
	double				alpha(significance);

	Breconstruction*	rec;
	Bparticle*			part, *vert;

	if ( verbose )
		cout << "Analyzing vertices:" << endl;

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose )
			cout << "Reconstruction: " << rec->id << endl;
		for ( part = vert = rec->part; part; part = vert, npart++ ) {
			// Calculate statistics for majority vertices (11)
			avg = 0;
			sig = 0;
			for ( vert = part; vert && vert->group == part->group;
						vert = vert->next ) {
				if ( vert->sel == 1 ) {
					avg += vert->fom[0];
					sig += vert->fom[0] * vert->fom[0];
				}
			}
			sig = (sig - avg/11)/10;
			sig /= 11;
			if ( sig > 0 ) sig = sqrt(sig);	// Standard deviation estimate
			else sig = 0;
			avg /= 11;						// Average
			avgall += avg;
			sigall += sig;
			
			for ( vert = part, i = non = nop = 0; vert && vert->group == part->group;
						vert = vert->next, ++i ) {
				t = (vert->fom[0] - avg)/sig;	// T statistic
				if ( t < -alpha ) non++;		// Negative outlier
				if ( t > alpha ) nop++;			// Positive outlier
			}
			if ( verbose )
				cout << part->group << tab << avg << tab << sig << tab << avg/sig << tab << non << tab << nop << endl;
			nvert += i;
			nontot += non;
			noptot += nop;
			if ( i < 12 )
				cerr << "Error: Too few vertices! (" << i << ")" << endl;
		}
	}
	
	avgall /= npart;
	sigall /= npart;
	
	if ( verbose ) {
		cout << "Number of particles:     " << npart << endl;
		cout << "Average:                 " << avgall << endl;
		cout << "Standard deviation average:   " << sigall << endl;
		cout << "Number of negative outliers:  " << nontot << endl;
		cout << "Number of positive outliers:  " << noptot << endl;
	}
	
	return npart;
}

long		project_analyze_vertices3(Bproject* project, long maxcut)
{
	long				i, nv(12), npart(0), cut;
	double				d, avg(0);
	vector<double>		fom(nv), h(nv, 0);
	
	if ( maxcut > nv ) maxcut = nv;

	Breconstruction*	rec;
	Bparticle*			part, *vert;

	if ( verbose )
		cout << "Analyzing vertices:" << endl;

	for ( rec = project->rec; rec; rec = rec->next ) if ( rec->part ) {
		if ( verbose )
			cout << "Reconstruction: " << rec->id << endl;
		for ( part = vert = rec->part; part; part = vert, npart++ ) {
			
			for ( vert = part, i = 0; vert && vert->group == part->group;
						vert = vert->next, ++i )
				fom[i] = vert->fom[0];
			sort(fom.begin(), fom.end());

			// Attempting to detect no special vertex
			for ( i=2, d=1e30; i<5; ++i )
				if ( d > fom[i] - fom[i-1] ) d = fom[i] - fom[i-1];
			if ( d > fom[1] - fom[0] ) {
				cut = 0;
			} else {
				for ( cut = 1; cut<maxcut; ++cut )
					if ( fom[cut+1] - fom[cut] < fom[cut] - fom[cut-1] ) break;
			}
			
			h[cut]++;
			
			cout << "Particle: " << part->group << tab << cut << endl;
			
			for ( i = 0; i<maxcut; ++i )
				cout << i+1 << tab << fom[i] << tab << fom[i+1] - fom[i] << endl;
		}
	}

	for ( cut = 0; cut<maxcut; ++cut ) avg += cut*h[cut];
	avg /= npart;

	cout << "Histogram of special vertex counts:" << endl;
	cout << "Vertices\tCount\t%\tBinomial" << endl;
	for ( cut = 0; cut<maxcut; ++cut )
		cout << cut << tab << h[cut] << tab << h[cut]*100.0/npart << tab << 
			100 * number_of_combinations(nv, cut) * pow(avg/nv, cut) * pow(1-avg/nv, nv-cut) << endl;
	cout << "Number of particles:     " << npart << endl;
	cout << "Average per particle:    " << avg << endl << endl;
	
	return npart;
}

