/**
@file	mg_class.cpp
@brief	Classifies raw single particle images with respect to multiple models
@author Bernard Heymann
@date	Created: 20010222
@date	Modified: 20150424
**/

#include "Bimage.h"
#include "rwimg.h"
#include "mg_class.h"
#include "mg_processing.h"
#include "mg_select.h"
#include "mg_ctf.h"
#include "Matrix3.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Classifies particle images with respect to a series of reference maps.
@param 	*project		project parameter structure.
@param 	resolution_hi 	upper resolution limit.
@param 	resolution_lo 	lower resolution limit.
@param 	fom_type		FOM type: 0=CC, 1=R, 2=PD.
@param 	fom_cut			input FOM threshold to classify a particle.
@param 	*kernel			frequency space interpolation kernel.
@param 	ctf_apply		apply CTF to projections.
@param 	img_out			image output options: 0=none, 1=projections, 2=differences, 3=both
@return long				particles selected.

	For every particle image, a projection is made from every reference map
	according to the input orientation parameters and compared to the
	particle image. The FOM calculated is a real space correlation coefficient.

**/
long   mg_classify(Bproject* project, double resolution_hi, double resolution_lo,
					int fom_type, double fom_cut, FSI_Kernel* kernel, int ctf_apply, int img_out)
{
	// Main loop over all reference maps
	bool			invert(0);
	long			i, j, m, nsel(0);
	double			maxfom, CC, PD(0), R;
	Vector3<double>	translate, origin;
	Matrix3 		mat(1);
	Bstring			filename;
	Bstring			name_base;
	Bstring			insert("_diff.");
	Bimage*			p = NULL;
	Bimage*			pmap = NULL;
	Bimage*			proj = NULL;
	Bimage*			pdiff = NULL;
	Bimage*			mapproj = NULL;

	long			nmap = count_list((char *)project->reference);
	Bstring*		map = project->reference;
	
	Bfield*			field;
	Bmicrograph*	mg;
	Bparticle*		part;
	long			npart = project_count_mg_particles(project);

	int*			num = new int[nmap+1];
	float*			fom = new float[npart*nmap];
	float*			fom_avg = new float[nmap+1];
	float*			fom_std = new float[nmap+1];
	
	for ( i=0; i<nmap+1; i++ ) fom_avg[i] = fom_std[i] = num[i] = 0;
	for ( i=0; i<npart*nmap; i++ ) fom[i] = 0;

	if ( verbose ) {
		cout << "Classifying particles with respect to maps:" << endl;
		cout << "Number of maps:                 " << nmap << endl;
		cout << "Resolution range:               " << resolution_hi << " - " << resolution_lo << " Å" << endl;
		cout << "FOM type:                       " << fom_type << endl;
		cout << "FOM cutoff:                     " << fom_cut << endl;
		if ( ctf_apply )
			cout << "Applying CTF parameters" << endl;
		if ( img_out%2 == 1 )
			cout << "Output projections" << endl;
		if ( img_out > 1 )
			cout << "Output differences" << endl;
	}

	pmap = read_img(*(project->reference), 0, 0);

	fft_plan		planf_2D = fft_setup_plan(pmap->sizeX(), pmap->sizeY(), 1, FFTW_FORWARD, 1);
	fft_plan		planb_2D = fft_setup_plan(pmap->sizeX(), pmap->sizeY(), 1, FFTW_BACKWARD, 1);
	
	delete pmap;
	
	if ( verbose & VERB_RESULT ) {
		cout << endl << "Map\tFile\tPart\t";
		if ( fom_type == 1 ) cout << "R" << endl;
		else if ( fom_type == 2 ) cout << "dPhi" << endl;
		else cout << "CC" << endl;
	}
	for ( i=0, map = project->reference; map; map = map->next, i++ ) {
		pmap = read_img(map->str(), 1, 0);
		if ( !pmap )
			return error_show("mg_classify", __FILE__, __LINE__);
		pmap->change_type(Float);
		pmap->statistics();
		pmap->calculate_background();
		origin = pmap->image->origin();
		if ( origin[0] <= 0 || origin[0] >= pmap->sizeX() ) origin[0] = pmap->sizeX()/2;
		if ( origin[1] <= 0 || origin[1] >= pmap->sizeY() ) origin[1] = pmap->sizeY()/2;
		if ( origin[2] <= 0 || origin[2] >= pmap->sizeZ() ) origin[2] = pmap->sizeZ()/2;
		pmap->origin(origin);
		if ( kernel ) {
			pmap->fft();
			pmap->phase_shift_to_origin();
		}
		name_base = Bstring(i+1, "map%02d");
		m = 0;
		for ( field=project->field; field; field=field->next ) {
			for ( mg=field->mg; mg; mg=mg->next ) {
				p = read_img(mg->fpart, 0, -1);
				if ( !p )
					return error_show("mg_classify", __FILE__, __LINE__);
				proj = p->copy_header(p->images());
				delete p;
//				filename = p->file_name();
				filename = mg->fpart;
				filename = filename.pre_rev('.') + "_" + name_base + "." + filename.post_rev('.');
				proj->file_name(filename.str());
				proj->data_type(Float);
				proj->data_alloc();
				pdiff = proj->copy_header(proj->images());
				filename = filename.pre_rev('.') + insert + filename.post_rev('.');
				pdiff->file_name(filename.str());
				pdiff->data_alloc();
				for ( j=0, part=mg->part; part; part=part->next, j++ ) {
					p = read_img(mg->fpart, 1, part->id-1);
					if ( !p )
						return error_show("Error in mg_classify", __FILE__, __LINE__);
					p->change_type(Float);
					p->statistics();
					if ( !p->rescale_to_avg_std(0, 1) ) {
						// Calculate the projected reference image
						if ( verbose & VERB_FULL )
							cout << "View and origin: " << part->view << tab << part->ori << endl;
						mat = part->view.matrix();
						if ( part->mag ) mat *= part->mag;
						translate[0] = part->ori[0] - pmap->sizeX()/2;
						translate[1] = part->ori[1] - pmap->sizeY()/2;
						if ( kernel ) {
							mapproj = pmap->central_section(mat, resolution_hi, kernel); 
							mapproj->phase_shift_to_center();
							if ( ctf_apply )
								img_ctf_apply_to_proj(mapproj, *(mg->ctf), part->def, 1e6, resolution_hi, invert, planf_2D, planb_2D);
							mapproj->fft_back(planb_2D);
							mapproj->shift(translate);
							mapproj->correct_background();
						} else {
							mapproj = pmap->rotate_project(mat, translate, pmap->sizeX()/2.0);
							if ( ctf_apply )
								img_ctf_apply_to_proj(mapproj, *(mg->ctf), part->def, 1e6, resolution_hi, invert, planf_2D, planb_2D);
						}
						mapproj->statistics();
						mapproj->rescale_to_avg_std(0, 1);
						// Compare two 2D images ===> FOM
						p->sampling(part->pixel_size);
						p->origin(part->ori);
						proj->image[j].origin(part->ori);
						pdiff->image[j].origin(part->ori);
						proj->image[j].view(part->view);
						pdiff->image[j].view(part->view);
						proj->replace(j, mapproj);
						if ( fom_type == 0 ) {
							CC = p->correlate(mapproj);
							fom[m*nmap+i] = CC;
						} else if ( fom_type == 2 ) {
							PD = mapproj->average_phase_difference(p,
								resolution_hi, resolution_lo, 1);
							fom[m*nmap+i] = cos(PD);
						}
						R = mapproj->linear_fit(p, NULL, 0);
						mapproj->invert();
						if ( fom_type == 1 ) fom[m*nmap+i] = 1 - R;
						pdiff->replace(j, mapproj);
						delete p;
						delete mapproj;
						if ( verbose & VERB_RESULT ) {
							cout << *map << tab << mg->fpart << tab << part->id << tab;
							if ( fom_type == 1 ) cout << R << endl;
							else if ( fom_type == 2 ) cout << PD*180.0/M_PI << endl;
							else cout << CC << endl;
						}
					}
					m++;
				}
				if ( img_out%2 == 1 ) {
					proj->statistics();
					cout << "Writing " << proj->file_name() << endl;
					write_img(proj->file_name(), proj, 0);
					mg->fpart = proj->file_name();
				}
				if ( img_out > 1 ) {
					pdiff->statistics();
					cout << "Writing " << pdiff->file_name() << endl;
					write_img(pdiff->file_name(), pdiff, 0);
					mg->fpart = pdiff->file_name();
				}
				delete proj;
				delete pdiff;
			}
		}
		delete pmap;
	}

    fft_destroy_plan(planf_2D);
    fft_destroy_plan(planb_2D);

	// Get the best FOM and assign the particle to that map
	if ( verbose & VERB_RESULT ) {
		cout << endl << "File\tPart\tMap\tFOM=";
		if ( fom_type == 1 ) cout << "1-R" << endl;
		else if ( fom_type == 2 ) cout << "cos(dPhi)" << endl;
		else cout << "CC" << endl;
	}
	m = nsel = 0;
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			for ( part=mg->part; part; part=part->next ) if ( part->sel ) {
				maxfom = -1e37;
				part->sel = 0;
				for ( i=0; i<nmap; i++ ) {
					if ( verbose & VERB_RESULT )
						cout << mg->fpart << tab << part->id << tab << i+1 << tab << fom[m*nmap+i] << endl;
					if ( maxfom < fom[m*nmap+i] ) {
						maxfom = fom[m*nmap+i];
						if ( maxfom >= fom_cut ) part->sel = i+1;
					}
				}
				num[part->sel]++;
				fom_avg[part->sel] += maxfom;
				fom_std[part->sel] += maxfom*maxfom;
				part->fom[0] = maxfom;
				if ( part->sel ) nsel++;
				m++;
			}
		}
	}
	
	if ( verbose )
		cout << "Map\tNumber\t%\tAvg\tStd" << endl;
	for ( i=0; i<=nmap; i++ ) {
		if ( num[i] ) {
			fom_avg[i] /= num[i];
			fom_std[i] = fom_std[i]/num[i] - fom_avg[i]*fom_avg[i];
			if ( fom_std[i] > 0 )
				fom_std[i] = sqrt(fom_std[i]);
			else
				fom_std[i] = 0;
		}
		if ( verbose )
			cout << i << tab << num[i] << tab << num[i]*100.0/npart << tab << fom_avg[i] << tab << fom_std[i] << endl;
	}
	
	delete[] num;
	delete[] fom;
	delete[] fom_avg;
	delete[] fom_std;

	return nsel;
}

