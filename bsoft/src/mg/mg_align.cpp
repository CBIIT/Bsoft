/**
@file	mg_align.cpp
@brief	Functions to align micrographs or coordinates from micrographs and apply the resultant transformation.
@author Bernard Heymann and Samuel Payne
@date	Created: 20000505
@date	Modified: 20211004
**/

#include "Bimage.h"
#include "mg_align.h"
#include "rwimg.h"
#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_ctf_fit.h"
#include "mg_select.h"
#include "mg_tomography.h"
#include "marker.h"
#include "matrix_linear.h"
#include "simplex.h"
#include "ps_marker.h"
#include "linked_list.h"
#include "utilities.h"

#include <sys/stat.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Aligns the particle coordinates of a series of micrographs.
@param 	*project	project parameter structure.
@param 	refset			reference image or set of coordinates (-1 means undefined).
@return int 				0.

	The sets of micrograph coordinates of a focal series in an image 
	processing parameter structure are fitted to each other, giving a 
	rotation angle, and shifting and scaling in x and y for each pair 
	of micrographs.
	There must be the same number of particles in each micrograph
	of a focal series.

**/
int			mg_align_coordinates(Bproject* project, int refset)
{
	if ( refset < 0 ) refset = 0;

	int 				i;
	long				nmg, npart, nmg_max = 1000;
	Transform			t;

	Bfield* 			field;
	Bmicrograph*		mg;
	Bparticle*			part;

	Bmarker**			set = new Bmarker*[nmg_max];
	Bmarker*			m = NULL;
	
	Vector3<double>		origin;

	if ( verbose & VERB_LABEL )
		cout << "Aligning micrograph particle coordinates for series of micrographs." << endl;
	if ( verbose & VERB_PROCESS ) 
		cout << "Particles\tX-shift\tY-shift\tX-scale\tY-scale\tAngle" << endl; 

	for ( field=project->field; field; field=field->next ) {
		for ( i=1, mg=field->mg; mg->next && i < refset; mg=mg->next ) i++;
		origin = micrograph_get_nominal_origin(mg);
		nmg = field_count_micrographs(field);
		npart = micrograph_count_particles(field->mg);
		for ( i=0, mg=field->mg; mg; mg=mg->next, i++ ) {
			for ( part=mg->part; part; part=part->next ) {
				m = (Bmarker *) add_item((char **) &m, sizeof(Bmarker));
				if ( !set[i] ) set[i] = m;
				m->loc = part->loc;
				m->sel = 1;
				m->fom = 1;
			}
		}
		for ( i=0, mg=field->mg->next; mg; mg=mg->next, i++ ) if ( i != refset ) {
			t = markers_find_rottrans(set[refset], set[i], 0.01);
			mg->rot_angle = t.angle;
			mg->origin = origin + t.trans;
			mg->scale = t.scale;
			if ( verbose & VERB_PROCESS ) 
				cout << npart << tab << t.trans[0] << tab << t.trans[1] << tab << t.trans[2] << tab <<
					t.scale[0] << tab << t.scale[1] << tab << t.scale[2] << tab << t.angle*180.0/M_PI << endl;
		}
		for ( i=0; i<nmg; i++ )
			kill_list((char *) set[i], sizeof(Bmarker));
	}
			
	delete[] set;
	
	return 0;
}

/*
@brief 	Aligns the second image to the first.
@param 	*p1				first image.
@param 	*p2				second image (transformed)
@param 	tile_size		3-valued vector for the size of sub-images.
@param 	res_lo			low resolution limit for cross-correlation.
@param 	res_hi			high resolution limit for cross-correlation.
@param 	max_shift		maximum shift allowed (default 1/4 of tile).
@param 	filter_flag		flag to filter micrograph extremes.
@param 	refine_flag 	flag to turn on refinement of shift.
@return Transform		structure with shift, scale, rotation angle, and R factor.

	The second image is transformed to fit on the first image.

*/
Transform	img_align(Bimage* p1, Bimage* p2, Vector3<long> tile_size, 
				double res_lo, double res_hi, double max_shift, int filter_flag, int refine_flag)
{
	long				i, iter, maxiter(5);
	Vector3<double>		shift;
	Transform			t, best_t;
	best_t.fom = 1e37;
	
	Vector3<long> 		start1, common_size, step_size;
	
	// Find the largest area covering both images
	common_size = p1->size().min(p2->size());
	common_size = common_size.max(1);
	
	// Ensure the tile fits into the images
	tile_size = tile_size.max(1);
	tile_size = tile_size.min(common_size);
	if ( tile_size[0] >= common_size[0] && common_size[0] > 1 ) {
		tile_size[0] = (int) (0.4*common_size[0]);
//		start1[0] = (common_size[0] - tile_size[0])/2;
	}
	if ( tile_size[1] >= common_size[1] && common_size[1] > 1 ) {
		tile_size[1] = (int) (0.4*common_size[1]);
//		start1[1] = (common_size[1] - tile_size[1])/2;
	}
	if ( tile_size[2] >= common_size[2] && common_size[2] > 1 ) {
		tile_size[2] = (int) (0.4*common_size[2]);
//		start1[2] = (common_size[2] - tile_size[2])/2;
	}

	if ( max_shift < 1 ) max_shift = tile_size[0]/4;
		
	Vector3<double>		half_tile;
	half_tile[0] = (tile_size[0] - 1)/2;
	half_tile[1] = (tile_size[1] - 1)/2;
	half_tile[2] = (tile_size[2] - 1)/2;
	
	Bimage*				pex1 = NULL;
	Bimage*				pex2 = NULL;
	Bmarker*			set1 = NULL;
	Bmarker*			set2 = NULL;
	Bmarker*			m1 = NULL;
	Bmarker*			m2 = NULL;
	Bstring				fn1(p1->file_name()), fn2(p2->file_name());
	Bstring				filename(fn1.base() + "_" + fn2.base() + ".ps");
	
	if ( verbose & VERB_LABEL ) {
		cout << "Aligning images ";
		if ( refine_flag ) cout << "with ";
		else cout << "without ";
		cout << "shift refinement." << endl << endl;
	}
	
	if ( common_size[0] > 2*tile_size[0] )
		start1[0] = (common_size[0]
				- (int) (tile_size[0]*(common_size[0]/tile_size[0])))/2;
	if ( common_size[1] > 2*tile_size[1] )
		start1[1] = (common_size[1]
				- (int) (tile_size[1]*(common_size[1]/tile_size[1])))/2;
	if ( common_size[2] > 2*tile_size[2] )
		start1[2] = (common_size[2]
				- (int) (tile_size[2]*(common_size[2]/tile_size[2])))/2;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Common size:                    " << common_size << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Maximum shift:                  " << max_shift << endl;
		cout << "Reference start:                " << start1 << endl << endl;
	}
	
	if ( res_hi < 0.1 || res_hi > tile_size[0]*p1->sampling(0)[0] )
		res_hi = 4.0*p1->sampling(0)[0];
	
	if ( filter_flag ) {
		p1->filter_extremes();
		p2->filter_extremes();
	}
	
	pex1 = p1->extract_tiles(0, start1, common_size, tile_size, step_size, 0);
	
//	Vector3<long>*	tile_start = new Vector3<long>[pex1->images()];
	vector<Vector3<long>>	tile_start(pex1->images());
	
	for ( i=0; i<pex1->images(); i++ ) {
		m1 = (Bmarker *) add_item((char **) &m1, sizeof(Bmarker));
		m2 = (Bmarker *) add_item((char **) &m2, sizeof(Bmarker));
		if ( !set1 ) set1 = m1;
		if ( !set2 ) set2 = m2;
		m1->sel = m2->sel = 1;
		m1->fom = m2->fom = 1;
		m1->loc = pex1->image[i].origin() + half_tile;
		tile_start[i] = pex1->image[i].origin();
	}
	
	for ( iter=0; iter<maxiter && best_t.fom>0.1; iter++ ) {
		pex2 = p2->extract_tiles(0, tile_start, tile_size);
		// Determines the shift of the second image onto the first
		pex1->find_shift(pex2, NULL, res_hi, res_lo, max_shift, 0, refine_flag);
		for ( i=0, m2=set2; i<pex1->images(); i++, m2=m2->next ) {
			shift = pex1->image[i].origin() - pex2->image[i].origin();
			if ( shift[0] > half_tile[0] ) shift[0] -= tile_size[0];
			if ( shift[1] > half_tile[1] ) shift[1] -= tile_size[1];
			if ( shift[2] > half_tile[2] ) shift[2] -= tile_size[2];
			// To get the location in the second image, the shift must be subtracted
			m2->loc = tile_start[i] + half_tile - shift;
		}
		t = markers_find_rottrans(set1, set2, 0.01);
		if ( verbose & VERB_PROCESS )
			cout << "Fit " << iter+1 << ":\t" << p1->file_name() << tab << p2->file_name() << tab 
				<< t.trans << tab
				<< t.scale << tab 
				<< t.angle*180.0/M_PI << tab << t.fom << endl;
		if ( best_t.fom > t.fom ) best_t = t;
		else iter = maxiter;
		for ( i=0, m1=set1, m2=set2; i<pex1->images(); i++, m1=m1->next, m2=m2->next ) {
			m2->loc[0] = best_t.scale[0]*(m1->loc[0]*cos(best_t.angle) -
					m1->loc[1]*sin(best_t.angle)) + t.trans[0] + 0.5;
			m2->loc[1] = best_t.scale[1]*(m1->loc[0]*sin(best_t.angle) +
					m1->loc[1]*cos(best_t.angle)) + t.trans[1] + 0.5;
			m2->loc[2] = best_t.scale[2]*m1->loc[2] + t.trans[2];
			tile_start[i] = m2->loc - half_tile;
		}
		delete pex2;
	}
	
	if ( verbose & VERB_RESULT ) 
		cout << "Best fit:\t" <<
				p1->file_name() << tab << p2->file_name() << tab <<
				best_t.trans << tab <<
				best_t.scale << tab <<
				best_t.angle*180.0/M_PI << tab << best_t.fom << endl;
	
	delete pex1;

	ps_marker_errors(filename, set1, set2, t, p1->size(), 20);

	kill_list((char *) set1, sizeof(Bmarker));
	kill_list((char *) set2, sizeof(Bmarker));
	
	return best_t;
}

Bimage*		img_sum_subset(Bimage* p, Bstring& subset, double dose_per_frame, double time_per_frame)
{
	long			n, ns(0);
	double			dose(0), exposure(0);
	vector<int>		numsel = select_numbers(subset, p->images());

	Bimage*			p1;
	Bimage*			pa = p->copy_header(1);
	pa->data_type(Float);
	pa->sampling(p->sampling(0));
	pa->data_alloc_and_clear();
	
	for ( n=0; n<p->images(); n++ ) {
		dose += dose_per_frame;
		exposure += time_per_frame;
		if ( numsel[n] ) {
			if ( verbose & VERB_PROCESS )
				cout << "Adding frame " << n << endl;
			p1 = p->extract(n);
			p1->change_type(Float);
			pa->add(p1);
			(*pa)["dose"] = dose;
			(*pa)["exposure"] = exposure;
			delete p1;
			ns++;
		}
	}
	
	double			var1(p->standard_deviation()*p->standard_deviation());
	double			var(pa->standard_deviation()*pa->standard_deviation());
	double			snr(0);
	if ( var < ns*ns*var1 ) snr = (var - ns*var1)/(ns*var1 - var/ns);
	
	pa->image->FOM(snr);
	
	return pa;
}


/**
@brief 	Aligns a focal series of micrographs.
@param 	*project		parameter structure.
@param 	refset			reference image or set of coordinates (< 1 means undefined).
@param 	tile_size		3-valued vector for the size of sub-images.
@param 	res_lo			low resolution limit for cross-correlation.
@param 	res_hi			high resolution limit for cross-correlation.
@param 	max_shift		maximum shift allowed (default 1/4 of tile).
@param 	filter_flag		flag to filter micrograph extremes.
@param 	refine_flag 		flag to turn on refinement of shift.
@return int				0, <0 on error.

	A series of micrograph images specified in an image processing 
	parameter structure are aligned by segmented cross-correlation. The 
	micrograph data blocks are assumed to be arranged with a series in
	consequent data blocks. The micrographs are segmented into tiles 
	and the tile shifts with respect to each other determined by 
	cross-correlation. The shifts are assumed to most accurately 
	represent the displacement of the center of one tile with respect 
	to the center of the corresponding tile in the other micrograph.
	The resultant sets of coordinates are fitted to each
	other, giving a 3-value shift vector, a 3-value scale vector,
	and a rotation angle for each pair of micrographs.
	A reference micrograph is chosen as:
		1.	the first micrograph with particle coordinates
		2.	otherwise, the first micrograph
	If coordinates are supplied for particles in the reference micrograph,
	the determined transformation parameters are applied and written into
	the other micrograph structures.

**/
int 		mg_align_micrographs(Bproject* project, int refset,
				Vector3<long> tile_size, double res_lo, double res_hi,
				double max_shift, int filter_flag, int refine_flag)
{
	long 				i;
	
	Bfield* 			field;
	Bmicrograph*		mg, *mg_ref;
	
	Bimage*				p = NULL;
	Bimage*				p2 = NULL;
	
	Transform			t;
	
	if ( verbose ) {
		cout << "Aligning micrographs for micrograph series ";
		if ( refine_flag ) cout << "with ";
		else cout << "without ";
		cout << "shift refinement." << endl << endl;
	}
	
	for ( field=project->field; field; field=field->next ) if ( field->mg ) {
		for ( i=1, mg = mg_ref = field->mg; mg->next && i < refset; mg=mg->next, i++ ) ;
		if ( mg ) mg_ref = mg;
		p = read_img(mg_ref->fmg, 1, mg->img_num);
		if ( !p ) {
			cerr << "Error in mg_align_micrographs: File " << mg_ref->fmg << " not read!" << endl;
			return -1;
		}
		if ( mg_ref->pixel_size[0] > 0.1 ) p->sampling(mg_ref->pixel_size);
		mg_ref->origin = p->size()/2;
		t = Transform();
		mg_ref->rot_angle = 0;
		mg_ref->scale = t.scale;
		for ( mg=field->mg; mg; mg=mg->next ) {
			if ( mg_ref != mg ) {
				p2 = read_img(mg->fmg, 1, mg->img_num);
				if ( !p ) {
					cerr << "Error in mg_align_micrographs: File " << mg->fmg << " not read!" << endl;
					return -1;
				}
				if ( mg->pixel_size[0] > 0.1 ) p2->sampling(mg->pixel_size);
				mg->origin = mg_ref->origin;
				t = Transform();
				// The alignment transforms the second image
//				cout << "aligning:" << tab << p->file_name() << tab << p2->file_name() << endl;
				t = img_align(p, p2, tile_size, res_lo, res_hi, max_shift, filter_flag, refine_flag);
				mg->rot_angle = t.angle;
				mg->origin += t.trans;
				mg->scale = t.scale;
				mg->matrix = Matrix3(t.axis, -t.angle);
				mg_apply_transform(mg_ref, mg);
				delete p2;
			}
		}
		delete p;
	}
	
	return 0;
}

/*
@author  Samuel Payne and Bernard Heymann
@brief 	Finds local minimums in a grey scale image and uses those as features
@param 	*p				image (modified to floating point).
@param 	extract_method	The method to determine the center of the particle
@param 	target_num 		Max number of desired features (-1 will return all)
@param 	threshold 		Threshold value used to determine whether
						a pixel is considered a feature.
@return Bmarker*		a set containing xyz coordinates of features

	The input image is expected to have been cleaned in the following manner.
	With each suggestion is listed the corresponding bsoft function call. 
	1.  Any extremes such as white lettering, or black edges have been removed.
		bfilter -e filname filename
	2. Edge filter to smooth out the edges.
		bedge -x<xsize>,<ysize>,<zsize> -o150,150,0 -g75 filename filename
	3. Bandpass filtered to reduce the noise.  
		bfilter -b1000,2000 -df filename filename
	
	The algorithm finds blocks of dark grey. The features are black dots 
	in the image, which could be particles, but don't have to be.  As long 
	as the picking is consistent, it does not matter what is picked.  
	(it could be a pickle).  The method used to find the center of the feature 
	depends on the input. (0=Center of Mass, 1=Center of Area, 2=Darkest pixel)

**/
Bmarker*	img_extract_features(Bimage* p, int extract_method,
					int target_num, double threshold)
{
	// The extreme considered as a feature depends on the relative threshold
	int				sign = 1;
	if ( threshold < p->average() ) sign = -1;

	long	   		i, j, x, y, z, n;
	long	  		count(0);
	Vector3<long>	coor;
	double			temp;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Using the ";
		switch ( extract_method ) {
			case 1: cout << "center of mass "; break;
			case 2: cout << "center of volume "; break;
			case 3: cout << "highest pixel density "; break;
			default: cout << "center of mass "; break;
		}
		cout << "to determine the feature coordinates in " << p->file_name() << endl;
	}
	
	// Definition of features based on regions beyond a threshold
	Bimage*			pmask = p->regions(threshold, 0);
	count = (long) pmask->maximum();
	
	if ( verbose & VERB_LABEL )
		cout << "Features extracted:             " << count << endl << endl;
	
	if ( count < 1 ) {
		delete pmask;
		return NULL;
	}	
			
	Bmarker*	set = NULL;
	Bmarker*	m = NULL;
	Bmarker**	lut = new Bmarker*[count];
	
	for ( i=0; i<count; i++ ) {
		m = (Bmarker *) add_item((char **) &m, sizeof(Bmarker));
		if ( !set ) set = m;
		m->id = i+1;
		m->sel = i+1;
		lut[i] = m;
	}
		
    for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			coor[2] = z;
			for ( y=0; y<p->sizeY(); y++ ) {
				coor[1] = y;
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					coor[0] = x;
					if ( (*pmask)[i] > 0 ) {	//meaning that it belongs to a group
						j = (long) ((*pmask)[i]-1);
						m = lut[j];
//						cout << "Region " << i << ": " << (*pmask)[i] << endl;
						temp = fabs((*p)[i]);
						if ( extract_method == 1 ) {
							m->loc += coor * temp;
							m->fom += temp;
						} else if ( extract_method == 2 ) {
							m->loc += coor;
							m->fom += 1;
						}
						if ( sign*(*p)[i] > sign*m->fom ) {
							m->fom = (*p)[i];
							if ( extract_method == 3 ) {
								m->loc = coor;
							}
						}
						m->img_num = n;
					}
				}				
			}
		}
	}

	delete pmask;
	delete[] lut;
	
	cout << "Feature marker list generated" << endl;
	
	markers_limit(target_num, sign, &set);
	
//	if ( verbose & VERB_FULL )
		cout << "Feature\tx\ty\tz\tex\tey\tez\tFOM" << endl;
	
	for ( m=set; m; m=m->next ) {
		if ( extract_method == 1 || extract_method == 2 )
			m->loc /= m->fom;
		m->err = m->loc - p->size() /2;
//		if ( verbose & VERB_FULL )
			cout << m->id << tab << m->loc << tab <<
				m->err << tab << m->fom << endl;
	}	
	if ( verbose & VERB_FULL )
		cout << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "\nDEBUG img_extract_features: The Feature set" << endl;
		for ( m=set; m; m=m->next )
			cout << "Feature " << m->id << ": " << m->loc << tab << m->fom << endl;
		pmask = p->mask_by_threshold(threshold);
		Bstring		filename(p->file_name());
		filename = filename.base() + "_features.pif";
		write_img(filename, pmask, 0);
		delete pmask;
	}

	
	return set;
}


/**
@author  Samuel Payne
@brief	Aligns images by feature extraction
@param 	*project		the set of image parameters
@param 	max_features	the maximum number of features to extract
@param 	res_low			the low resolution limit
@param 	res_high		the high resolution limit
@param 	thresh			threshold used for feature extraction
@param 	extract_method	method used for finding the center of particles
@return int 				error code.

Calculates the transformation parameters for each set images, by 
	picking features in the images and finding the best fit for the 
	matching of the features.

	The features are picked and returned using function 
	img_extract_features().  The set of features is sent to
	find_transform_params() that finds the transformation parameters
	that relate the images to each other.  The shift in xy directions, 
	the scale in xy directions, and rotation angle is returned.  The 
	program works best if the images have at least 40 distinct features.

**/
int         mg_align_feature_extraction(Bproject* project, int max_features,
				double res_low, double res_high, double thresh, int extract_method)
{
	long 				i, refset(-1), nmg;
	Vector3<long>		size;
	Vector3<double>  	origin;
	double 				width;
	Transform*			tlist = NULL;
	Transform*			t = NULL;

	Bfield* 			field;
	Bmicrograph*		mg, *mg_ref = NULL;
	
	long				nmg_max = field_count_micrographs(project->field);
	if ( nmg_max < 5 ) nmg_max = 5;

	Bmarker**			set = new Bmarker*[nmg_max];
		
	Bimage*				p = NULL;
	
	if ( verbose & VERB_LABEL )
		cout << "Aligning micrographs for micrograph series using features." << endl;

	for ( field=project->field; field; field=field->next ) {
		nmg = field_count_micrographs(field);
		if ( refset < 0 ) refset = 0;
		for ( i=0, mg=field->mg; mg; mg=mg->next, i++ ) {
			if ( mg->part ) refset = i;
			if ( i == refset ) mg_ref = mg;
			p = read_img(mg->fmg, 1, 0);
			mg->origin = p->size()/2;
			// Before features are extracted, the image must be filtered for extreme
			// values, the edge smoothed and finally bandpass-filtered.
			p->filter_extremes();
			width = 0.05*p->sizeX();
			if ( p->sizeY() > p->sizeX() ) width = 0.05*p->sizeY();
			origin = p->size() * 0.1;
			if ( origin[0] > p->sizeX() - 1 ) origin[0] = p->sizeX() - 1;
			if ( origin[1] > p->sizeY() - 1 ) origin[1] = p->sizeY() - 1;
			if ( origin[2] > p->sizeZ() - 1 ) origin[2] = p->sizeZ() - 1;
			size[0] = (int) (p->sizeX() - 2*origin[0]);
			size[1] = (int) (p->sizeY() - 2*origin[1]);
			size[2] = (int) (p->sizeZ() - 2*origin[2]);
			if ( size[0] > p->sizeX() ) size[0] = p->sizeX();
			if ( size[1] > p->sizeY() ) size[1] = p->sizeY();
			if ( size[2] > p->sizeZ() ) size[2] = p->sizeZ();
			p->shape(0, size, origin, width, FILL_AVERAGE, p->average());
			p->fspace_bandpass(res_high, res_low, 0);
			set[i] = img_extract_features(p, extract_method, max_features, thresh);
			if ( set[i] ==  NULL) {
				error_show("Error in mg_align_feature_extraction", __FILE__, __LINE__);
				cerr << "Unable to extract any features from " << p->file_name() << endl;
				return -1;
			}
			size = p->size();
			delete p;
		}
		tlist = markers_map_and_find_transform(set, nmg, refset, size);
		if ( !t->fom ) {
			error_show("Error in mg_align_feature_extraction", __FILE__, __LINE__);
			cerr << "Transformation parameters not found!" << endl;
			return -1;
		}
		for ( i=0, mg=field->mg, t=tlist; mg; mg=mg->next, t=t->next, i++ ) if ( i != refset ) {
			mg->rot_angle = t->angle;
			mg->origin = mg->origin + t->trans;
			mg->scale = t->scale;
			mg_apply_transform(mg_ref, mg);
		}
		for ( i=0, mg=field->mg; mg; mg=mg->next, i++ )
			kill_list((char *)set[i], sizeof(Bmarker));
	}
	
	delete[] set;
	kill_list((char *) tlist, sizeof(Transform));
	
	return 0;
}

/**
@brief 	Applies a transformation to particle coordinates in a reference micrograph and writes them into an application micrograph structure.
@param 	*mg_ref		micrograph used as reference.
@param 	*mg_apply	micrograph to apply transformation to.
@return int 			0.

	The transformation parameters specified in the second micrograph
	are used.

**/
int			mg_apply_transform(Bmicrograph* mg_ref, Bmicrograph* mg_apply)
{
	if ( mg_ref == mg_apply ) return 0;
	
	long				npref = micrograph_count_particles(mg_ref);
	long				npapp = micrograph_count_particles(mg_apply);
	long				nfref = filament_node_count(mg_ref->fil);
	long				nfapp = filament_node_count(mg_apply->fil);
	Bparticle			*part_ref, *part_app;
	Bbadarea			*bad_ref, *bad_app;
	Bfilament			*fil_ref, *fil_app;
	Bfilnode			*fn_ref, *fn_app;
	Bmarker				*mark_ref, *mark_app;
	
//	double				angle = mg_ref->rot_angle - mg_apply->rot_angle;
//	double				cos_ang = cos(angle), sin_ang = sin(angle);
//	Vector3<float>		shift = mg_ref->origin - mg_apply->origin;
//	Vector3<float>		scale = mg_ref->scale / mg_apply->scale;

	double				angle = mg_apply->rot_angle - mg_ref->rot_angle;
	double				cos_ang = cos(angle), sin_ang = sin(angle);
	Vector3<double>		shift = mg_apply->origin - mg_ref->origin;
	Vector3<double>		scale = mg_apply->scale / mg_ref->scale;
	Vector3<double>		loc;

//	if ( npref < 1 ) {
//		cerr << "Warning: No particle coordinates for reference micrograph " << mg_ref->id << endl;
//		return 0;
//	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying transform parameters to micrographs:" << endl;
		cout << "Reference micrograph:           " << mg_ref->id << endl;
		cout << "Transform micrograph:           " << mg_apply->id << endl;
		cout << "Shift:                          " << shift << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Rotation angle:                 " << angle*180.0/M_PI << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Applying transform parameters to micrographs" << endl << endl;

	mg_apply->box_size = mg_ref->box_size;
	mg_apply->bad_radius = mg_ref->bad_radius;
	mg_apply->filament_width = mg_ref->filament_width;
	mg_apply->fil_node_radius = mg_ref->fil_node_radius;
	mg_apply->mark_radius = mg_ref->mark_radius;
	
	if ( npapp != npref ) {
		particle_kill(mg_apply->part);
		mg_apply->part = part_app = NULL;
		for ( part_ref=mg_ref->part; part_ref; part_ref=part_ref->next ) {
			part_app = particle_add(&part_app, part_ref->id);
			if ( !mg_apply->part ) mg_apply->part = part_app;
		}
	}
	
	if ( mg_ref->part ) {
		for ( part_ref=mg_ref->part, part_app=mg_apply->part; part_ref; 
				part_ref=part_ref->next, part_app=part_app->next ) {
			part_app->id = part_ref->id;
			loc = scale*(part_ref->loc + shift);
			part_app->loc[0] = loc[0] * cos_ang - loc[1] * sin_ang;
			part_app->loc[1] = loc[0] * sin_ang + loc[1] * cos_ang;
		}
	}
	
	if ( mg_apply->bad ) {
		kill_list((char *) mg_apply->bad, sizeof(Bbadarea));
		mg_apply->bad = NULL;
	}
	
	if ( mg_ref->bad ) {
		for ( bad_ref=mg_ref->bad; bad_ref; bad_ref=bad_ref->next ) {
			bad_app = (Bbadarea *) add_item((char **) &mg_apply->bad, sizeof(Bbadarea));
			loc = scale*(bad_ref->loc + shift);
			bad_app->loc[0] = loc[0] * cos_ang - loc[1] * sin_ang;
			bad_app->loc[1] = loc[0] * sin_ang + loc[1] * cos_ang;
		}
	}
	
	if ( nfapp != nfref ) {
		filament_kill(mg_apply->fil);
		mg_apply->fil = fil_app = NULL;
		for ( fil_ref=mg_ref->fil; fil_ref; fil_ref=fil_ref->next ) {
			fil_app = filament_add(&fil_app, fil_ref->id);
			if ( !mg_apply->fil ) mg_apply->fil = fil_app;
			fil_app->node = fn_app = NULL;
			for ( fn_ref = fil_ref->node; fn_ref; fn_ref = fn_ref->next ) {
				fn_app = filament_node_add(&fn_app, fn_ref->id);
				if ( !fil_app->node ) fil_app->node = fn_app;
			}
		}
	}
	
	if ( mg_ref->fil ) {
		for ( fil_ref=mg_ref->fil, fil_app=mg_apply->fil; fil_ref; 
				fil_ref=fil_ref->next, fil_app=fil_app->next ) {
			fil_app->id = fil_ref->id;
			for ( fn_ref=fil_ref->node, fn_app=fil_app->node; fn_ref; 
					fn_ref=fn_ref->next, fn_app=fn_app->next ) {
				fn_app->id = fn_ref->id;
				loc = scale*(fn_ref->loc + shift);
				fn_app->loc[0] = loc[0] * cos_ang - loc[1] * sin_ang;
				fn_app->loc[1] = loc[0] * sin_ang + loc[1] * cos_ang;
			}
		}
	}
	
	if ( mg_apply->mark ) {
		kill_list((char *) mg_apply->mark, sizeof(Bmarker));
		mg_apply->mark = NULL;
	}
	
	if ( mg_ref->mark ) {
		for ( mark_ref=mg_ref->mark; mark_ref; mark_ref=mark_ref->next ) {
			mark_app = (Bmarker *) add_item((char **) &mg_apply->mark, sizeof(Bmarker));
			loc = scale*(mark_ref->loc + shift);
			mark_app->loc[0] = loc[0] * cos_ang - loc[1] * sin_ang;
			mark_app->loc[1] = loc[0] * sin_ang + loc[1] * cos_ang;
			mark_app->sel = mark_ref->sel;
			mark_app->fom = mark_ref->fom;
		}
	}
	
	return 0;
}

/**
@brief 	Merges corresponding particle images in each focal series.
@param 	*project		project parameter structure.
@param 	use_old_origins	flag to use old origins rather than cross-correlation.
@return int 				error code.
**/
int			mg_merge_focal_series(Bproject* project, int use_old_origins)
{
	int 			i, ifield(0);
	double			cc;
	
	Bimage*			p = NULL;				// Particle image
	Bimage*			pt = NULL;				// Temporary image
	Bimage*			pref = NULL;			// Particle reference image
	Bimage*			ptr = NULL;				// Temporary reference image
	Bimage*			psum = NULL;			// Particle sum image

//	Vector3<float>*	shift = NULL;			// Shift vectors
//	Vector3<float>*	origin = NULL;			// Origin vectors
	
	Bfield*			field;
	Bmicrograph		*mg, *mg2, *mg_ref;
	Bparticle		*part, *part_ref;
	Bstring			filename;
	Bstring			insert("_merged.");
	Vector3<double>	origin, scale(1,1,1), axis(0,0,1), shift;
	Matrix3			mat;
	
	for ( ifield = 0, field=project->field; field; field=field->next, ifield++ ) {
		// Find the closest-to-focus particle image
		for ( mg_ref = mg = field->mg; mg; mg = mg->next )
			if ( mg_ref->ctf && mg->ctf )
				if ( mg_ref->ctf->defocus_average() < mg->ctf->defocus_average() ) mg_ref = mg;
		pref = read_img(mg_ref->fpart, 1, -1);
		if ( !pref )
			return error_show("Error in mg_merge_focal_series", __FILE__, __LINE__);
		pref->origin(pref->size()/2);
		pref->change_type(Float);
		psum = pref->copy();
		for ( mg = field->mg; mg; mg = mg->next ) if ( mg != mg_ref ) {
			p = read_img(mg->fpart, 1, -1);
			if ( !p )
				return error_show("Error in mg_merge_focal_series", __FILE__, __LINE__);
			p->origin(p->size()/2);
			p->change_type(Float);
/*			origin = new Vector3<float>[p->images()];
			if ( use_old_origins ) {
				shift = new Vector3<float>[p->images()];
				for ( i=0, part = mg->part, part_ref = mg_ref->part; part;
						i++, part = part->next, part_ref = part_ref->next ) {
					shift[i] = part_ref->ori - part->ori;
					origin[i] = part->ori;
				}
			} else {
//				shift = img_find_shift(pref, p, NULL, p->sampling()[0]*2, 1e10, p->sizeX()/4.0, 0, 1, NULL);
//				for ( i=0, part = mg->part; part; i++, part = part->next )
//					origin[i] = mg->origin - shift[i];
				p->find_shift(pref, NULL, p->sampling()[0]*2, 1e10, p->sizeX()/4.0, 0, 1);
				for ( i=0, part = mg->part; part; i++, part = part->next )
					origin[i] = mg->origin + pref->image[i].origin() - p->image[i].origin();
			}
			img_unique_shift_global_rotate(p, shift, origin, mg->rot_angle - mg_ref->rot_angle);
			delete[] shift;
			delete[] origin;
//			img_add(psum, p, 1, 0);
*/
			mat = Matrix3(axis, mg->rot_angle - mg_ref->rot_angle);
			for ( i=0, part = mg->part, part_ref = mg_ref->part; part;
						i++, part = part->next, part_ref = part_ref->next ) {
				pt = p->extract(i);
				if ( use_old_origins ) {
					shift = part_ref->ori - part->ori;
					origin = part->ori;
				} else {
					ptr = pref->extract(i);
					shift = pt->find_shift(ptr, NULL, ptr->sampling(0)[0]*2, 1e10, ptr->sizeX()/4.0, 0, 1, cc);
					origin = ptr->image->origin() - shift;
					delete ptr;
				}
				pt->transform(scale, origin, shift, mat, FILL_BACKGROUND, 0);
				p->replace(i, pt);
			}
			psum->add(p);
			delete p;
		}
		delete pref;
		filename = mg_ref->fpart;
		filename = filename.base() + insert + filename.post_rev('.');
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG mg_merge_focal_series: mg=" << mg_ref->id << " filename=" << filename << endl;
		write_img(filename, psum, 0);
		delete psum;
		for ( mg = field->mg; mg; ) {
			mg2 = mg->next;
			if ( mg != mg_ref ) micrograph_kill(mg);
			mg = mg2;
		}
		field->mg = mg_ref;
		mg_ref->next = NULL;
		mg_ref->block = ifield;
		mg_ref->fpart = filename;
	}

	return 0;
}

double		field_write_aligned_average(Bfield* field, Bimage* pgr, 
				Bstring& imgfile, DataType datatype, Bstring& subset)
{
	long	i, n;
	Bmicrograph*	mg;
	Vector3<double>	ori;
	Bimage* 		p;
	Bimage* 		pavg = NULL;
	Bimage* 		pstd = NULL;
	Vector3<double> scale(1,1,1);
	Matrix3 		mat(1);
	Vector3<double>	nori;

	for ( n=0, mg = field->mg; mg; mg = mg->next, n++ ) ;
//	int*			numsel = new int[n];
//	select_numbers(subset, n, numsel);
	vector<int>		numsel = select_numbers(subset, n);
	
	for ( i=n=0, mg = field->mg; mg; mg = mg->next, i++ ) if ( numsel[i] ){
		p = read_img(mg->fmg, 1, mg->img_num);
		if ( datatype == Unknown_Type ) datatype = p->data_type();
		if ( datatype > p->data_type() ) p->change_type(datatype);
		if ( pgr ) p->multiply(pgr);
		p->calculate_background();
		ori = p->size()/2;
		p->transform(scale, nori, (ori - mg->origin), mat, FILL_BACKGROUND, 0);
		p->change_type(datatype);
		if ( !pavg ) pavg = p->copy();
		else pavg->add(p);
		p->square();
		if ( !pstd ) pstd = p->copy();
		else pstd->add(p);
		delete p;
		n++;
	}
	
//	delete[] numsel;
	
	if ( n ) {
		pavg->multiply(1.0L/n);
		pstd->multiply(1.0L/n);
	} else {
		field->select = 0;
		return -1;
	}
	
	double			snr = pavg->standard_deviation()*pavg->standard_deviation();
	
	mg = field->mg;
	micrograph_kill(mg->next);
	mg->next = NULL;
	
	View			view = View(mg->matrix);
	pavg->image->origin(ori);
	pavg->image->view(view);
	if ( mg->pixel_size[0] ) pavg->sampling(mg->pixel_size);
	
	if ( imgfile.length() && pavg ) {
		write_img(imgfile, pavg, 0);
		mg->fmg = imgfile;
		mg->img_num = 0;
	}
	
	pavg->square();
	pstd->subtract(pavg);
	snr /= pstd->average();

	delete pavg;
	delete pstd;
	
	field->fom = snr;
		
	return snr;
}

/**
@brief 	Calculates aligned micrograph images and write them to a file.
@param 	*project   		micrograph project.
@param 	*pgr			gain reference.
@param 	&imgfile		output image file name.
@param 	datatype   		output data type.
@return double			0.

	Only the origin is adjusted.

**/
double		project_write_aligned_images(Bproject* project, Bimage* pgr,
				Bstring& imgfile, DataType datatype)
{
	long			i, nmg;
	Bfield*			field = project->field;
	Bmicrograph*	mg = field->mg;
	Bimage* 		p, *pone, *paln;
	Vector3<double>	origin;
	Vector3<double> scale(1,1,1);
	Matrix3 		mat(1);
	Vector3<double>	nori;
	View			view;
	
	for ( field = project->field; field; field = field->next ) {
		for ( nmg=0, mg = field->mg; mg; mg = mg->next ) nmg++;
		mg = field->mg;
		p = read_img(mg->fmg, 0, -1);
		if ( datatype == Unknown_Type ) datatype = p->data_type();
		paln = new Bimage(datatype, p->compound_type(), p->sizeX(), p->sizeY(), 1, nmg);
		if ( p->sizeZ() > 1 || p->images() > 1 ) {
			delete p;
			p = read_img(mg->fmg, 1, -1);
			if ( p->sizeZ() > 1 ) p->slices_to_images();
		} else {
			delete p;
			p = NULL;
		}
		if ( pgr ) p->multiply(pgr);
		if ( project->field->next )
			imgfile = mg->fmg.base() + "_aln." + mg->fmg.post_rev('.');
		for ( i=0, mg = field->mg; mg; mg = mg->next, i++ ) {
			view = View(mg->matrix);
			if ( p ) pone = p->extract(i);
			else pone = read_img(mg->fmg, 1, mg->img_num);
			if ( datatype > pone->data_type() ) pone->change_type(datatype);
			origin = pone->size()/2;
			if ( origin != mg->origin ) {
				pone->calculate_background();
				pone->transform(scale, nori, (origin - mg->origin), mat, FILL_BACKGROUND, 0);
			}
			pone->change_type(datatype);
			paln->replace(i, pone);
			paln->image[i].origin(origin);
			paln->image[i].view(view);
			if ( mg->pixel_size[0] ) paln->sampling(mg->pixel_size);
			if ( imgfile.length() ) {
				mg->fmg = imgfile;
				mg->img_num = i;
				mg->origin = origin;
			}
			delete pone;
		}
		if ( imgfile.length() && paln )
			write_img(imgfile, paln, 0);
		delete paln;
	}
	
	return 0;
}

/**
@brief 	Calculates aligned micrograph images and write them to a file.
@param 	*project		micrograph project.
@param 	*pgr			gain reference.
@param 	&imgfile		output image file name.
@param 	datatype		output data type.
@param 	&subset			subset to average (all if empty)
@return double			average SNR.

	Only the origin is adjusted.

**/
double		project_write_aligned_averages(Bproject* project, Bimage* pgr,
				Bstring& imgfile, DataType datatype, Bstring& subset)
{
	Bfield*			field = project->field;
	long			n, nfield(0);
	double			snr_avg(0);
	
	if ( verbose & VERB_RESULT ) {
		cout << "Writing micrograph average to " << imgfile << endl;
		if ( subset.length() )
			cout << "Subset:                 " << subset << endl;
	}

	Bstring			path;
	if ( imgfile.contains("/") ) path = imgfile.pre_rev('/') + "/";

	for ( field = project->field; field; field = field->next ) nfield++;
	
	Bfield**		farr = new Bfield*[nfield];
	for ( n=0, field = project->field; field; field = field->next, n++ )
		farr[n] = field;
	
#ifdef HAVE_GCD
	dispatch_apply(nfield, dispatch_get_global_queue(0, 0), ^(size_t i){
		Bmicrograph*	mg = farr[i]->mg;
		Bstring			newname(imgfile);
		if ( project->field->next ) {
			newname = mg->fmg.pre_rev('.') + "_avg." + mg->fmg.post_rev('.');
			if ( path.length() ) newname = path + newname.post_rev('/');
		}
		field_write_aligned_average(farr[i], pgr, newname, datatype, subset);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nfield; i++ ) {
		Bmicrograph*	mg = farr[i]->mg;
		Bstring			newname(imgfile);
		if ( project->field->next ) {
			newname = mg->fmg.pre_rev('.') + "_avg." + mg->fmg.post_rev('.');
			if ( path.length() ) newname = path + newname.post_rev('/');
		}
		field_write_aligned_average(farr[i], pgr, newname, datatype, subset);
	}
#endif

	delete[] farr;

	if ( verbose & VERB_RESULT )
		cout << "Field\tSNR" << endl;

	for ( n=0, field = project->field; field; field = field->next ) {
		snr_avg += field->fom;
		n++;
		if ( verbose & VERB_RESULT )
			cout << field->id << tab << field->fom << endl;
	}
	
	if ( n ) snr_avg /= n;
		
	return snr_avg;
}

double		mg_write_frame_sum(Bmicrograph* mg, Bimage* pgr,
				DataType datatype, Bstring& subset, double sampling_ratio, int flag)
{
	if ( !mg->frame ) return 0;
	if ( sampling_ratio < 1 ) sampling_ratio = 1;
	
	long			n;
	Bframe*			frame;

	Bimage*			p = read_img(mg->fframe, 1, -1);

	if ( p->sizeZ() > 1 ) p->slices_to_images();
	p->sampling(mg->frame_pixel_size);

	if ( pgr ) p->multiply(pgr);

	if ( flag & 1 ) p->histogram_counts(2);

	long			nimg = p->set_subset_selection(subset);
	
	if ( verbose )
		cout << "Micrograph: " << mg->id << " with " << nimg << " frames at " << mg->dose/p->images() << " e/Å2/frame" << endl;

	for ( n=0, frame = mg->frame; frame; frame = frame->next, ++n )
		p->image[n].origin(frame->shift + p->size()/2);

	Bimage*			psum = p->fspace_shift_sum();

	double			res_hi(p->image->sampling()[0]*2);
	Bplot*			plot = psum->fspace_ssnr(nimg, res_hi, sampling_ratio);

	psum->fft_back();
	
	double			dose(0), exposure(0), dose_per_frame(mg->dose/p->images()), time_per_frame(mg->exposure/p->images());
	for ( long i=0; i<p->images(); ++i ) {
		dose += dose_per_frame;
		exposure += time_per_frame;
		if ( p->image[i].select() ) {
			mg->dose = dose;
			mg->exposure = exposure;
		}
	}
	(*psum)["dose"] = dose;
	(*psum)["exposure"] = exposure;

	Bstring			ext(mg->fframe.extension());
	if ( ext.contains("dm") ) ext = "mrc";
	if ( ext.contains("tif") ) ext = "mrc";
	if ( ext.contains("eer") ) ext = "mrc";

	mg->fmg = mg->fframe.base() + "_sum." + ext;
	mg->pixel_size = psum->sampling(0);
	mg->origin = psum->size()/2;
	psum->origin(mg->origin);
	mg->intensity = psum->average();
	
	if ( verbose )
		cout << "Writing aligned frame sum: " << mg->fmg << endl;
	
	write_img(mg->fmg, psum, 0);
	
	delete psum;
	
	Bstring			psfile = mg->fframe.base() + "_ssnr.ps";
	ps_plot(psfile, plot);
	delete plot;

	return 0;
}

/**
@brief 	Sums aligned micrograph frames and write them to files.
@param 	*project		micrograph project.
@param 	*pgr			gain reference.
@param 	datatype		output data type.
@param 	&subset			subset to average (all if empty).
@param 	sampling_ratio	radial sampling ratio (1 or larger).
@param 	flag			flag to calculate counts from histogram.
@return double			average SNR.

	Only the origin is adjusted.

**/
double		project_write_frame_sums(Bproject* project, Bimage* pgr,
				DataType datatype, Bstring& subset, double sampling_ratio, int flag)
{
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Writing micrograph frame sums" << endl;
		if ( subset.length() )
			cout << "Subset:                 " << subset << endl;
	}

	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_write_frame_sum(mg, pgr, datatype, subset, sampling_ratio, flag);
			
	return 0;
}

double		mg_align_frames(Bmicrograph* mg, long ref_num, long window, long step,
				Bimage* pgr, Bimage* pmask, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width, 
				long bin, Bstring& subset, int flag)
{
	if ( bin < 1 ) bin = 1;
	
	int				mode(flag&16);
	
	Bimage*			p = read_img(mg->fframe, bin, -1);
	Bstring			ext(mg->fframe.extension());
	if ( ext.contains("eer") ) bin = 1;

	if ( mg->frame_pixel_size.volume() ) p->sampling(mg->frame_pixel_size);
//	cout << mg->frame_pixel_size << endl;
//	cout << p->sampling(0) << endl;
	
	if ( mg->origin.length() > 0 ) p->origin(mg->origin);
	
	if ( p->sizeZ() > 1 ) p->slices_to_images();
		
	if ( p->images() < 2 ) {
		delete p;
		return 0;
	}
	
//	p->information();

	if ( pgr ) {
		if ( verbose )
			cout << "Multiplying with " << pgr->file_name() << endl;
		p->multiply(pgr);
	}
	
	if ( bin > 1 ) {
		if ( verbose )
			cout << "Binning by " << bin << endl;
		p->bin(bin);
//		mg->pixel_size = p->image->sampling();
	}
	
//	write_img("pb.mrc", p, 0);
	
	if ( verbose )
		cout << "Aligning frames from micrograph " << mg->id << endl;

	if ( flag & 1 ) p->histogram_counts(2);

	Vector3<long>	aln_bin(1,1,1);
	if ( hi_res > 3*p->image->sampling()[0] ) aln_bin = Vector3<long>(2,2,1);
	
	vector<Vector3<double>>	sh = p->align(ref_num, window, step, pmask, hi_res, lo_res, shift_limit,
						edge_width, gauss_width, aln_bin, mode);

	long			i;
	double			d, cc_avg(0), shift_avg(0), shift_var(0);
	Vector3<double>	shift, pshift;
	Bframe*			frame = NULL;
	
	for ( i=0, frame = mg->frame; i<p->images(); ++i, frame = frame->next ) {
		if ( !frame ) frame = frame_add(&mg->frame, i+1);
		frame->shift[0] = sh[i][0];
		frame->shift[1] = sh[i][1];
		frame->fom = sh[i][2];
		if ( i ) {
			shift = (frame->shift - pshift)*mg->pixel_size;
			d = shift.length();
			shift_avg += d;
			shift_var += d*d;
		}
		cc_avg += sh[i][2];
		pshift = frame->shift;
	}
	shift_avg /= p->images()-1;
	shift_var /= p->images()-1;
	shift_var -= shift_avg*shift_avg;
	mg->fom = cc_avg /= p->images();
	
	if ( verbose ) {
		cout << "Shift average per frame:       " << shift_avg << " A" << endl;
		cout << "Shift variance per frame:      " << shift_var << " A2" << endl;
		if ( mg->exposure ) {
			cout << "Exposure:                      " << mg->exposure << " s" << endl;
			cout << "Average movement:              " << 
				(shift_avg/mg->exposure)*p->images() << " A/s" << endl << endl;
		}
	}

	delete p;
	
	return mg->fom;
}

/**
@brief 	Aligns a series of micrographs by cross-correlation.
@param 	*project		project parameter structure.
@param 	ref_img			reference frame number (starts from 0).
@param 	window			moving sum window (default 1, no moving sum).
@param 	step			moving sum interval (default 1).
@param 	*pgr			gain reference.
@param 	*pmask			reciprocal space mask, 0's and 1's.
@param 	origin			tilt origin.
@param 	hi_res			high resolution limit.
@param 	lo_res			low resolution limit.
@param 	shift_limit		maximum shift from nominal origin of image.
@param 	edge_width		edge smoothing width (not done if 0).
@param 	gauss_width		edge decay width.
@param 	bin				integer bin factor.
@param	subset			a subset to sum.
@param 	flag			options flag.
@return double			root-mean-square of offsets.

	Each micrograph frame is cross-correlated with the reference
	frame and the shift determined.
	Options encoded in the flag:
	1	rescale image based on histogram.
	2	weigh by accumulated dose.
	4	write aligned frames with insert "_aln".
	8	write aligned frame sum with insert "_sum".
	16	initial alignment: local rather than progressive.

**/
double		project_align_frames(Bproject* project, int ref_img, long window, long step,
				Bimage* pgr, Bimage* pmask, Vector3<double> origin, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width,
				long bin, Bstring& subset, int flag)
{
	if ( bin < 1 ) bin = 1;
	if ( window < 1 ) window = 1;
	
	Bfield*			field = project->field;
	Bmicrograph*	mg;
	long			nmg = project_count_micrographs(project);
	double			d(0);

//	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
	if ( verbose ) {
		cout << "Aligning micrograph frames by cross-correlation:" << endl;
		cout << "Number of micrographs:          " << nmg << endl;
		cout << "Reference frame:                " << ref_img << endl;
		cout << "Moving sum window and step:     " << window << tab << step << endl;
		cout << "Resolution limits:              " << hi_res << " - " << lo_res << " A" << endl;
		cout << "Shift limit:                    " << shift_limit << endl;
		cout << "Edge masking width & smoothing: " << edge_width << " " << gauss_width << endl;
		cout << "Binning:                        " << bin << endl;
		if ( pgr )
			cout << "Gain reference file:            " << pgr->file_name() << endl;
		if ( pmask )
			cout << "Using mask file:                " << pmask->file_name() << endl;
		if ( flag&16 )
			cout << "Initial alignment mode:         local" << endl;
		else
			cout << "Initial alignment mode:         progressive" << endl;
		cout << endl;
	}

	for ( nmg=0, field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next, nmg++ ) {
			mg_align_frames(mg, ref_img, window, step, pgr, pmask, hi_res, lo_res, shift_limit,
				edge_width, gauss_width, bin, subset, flag);
			d += mg->fom;
		}
	}
	
	d /= nmg;

	if ( verbose & VERB_RESULT )
		cout << "Overall average shift per frame: " << d << endl << endl;
	
	return d;
}


double		mg_align_series(Bmicrograph* mg, Bmicrograph* mg_ref, Bimage* pref, 
				Bimage* pgr, Bimage* pmask, double hi_res, double lo_res, double shift_limit, 
				double edge_width, double gauss_width,
				long bin, fft_plan planf, fft_plan planb)
{
	if ( mg == mg_ref ) return 0;
	
	double			cc(0);
	Vector3<long>	obin(bin, bin, 1);
	Vector3<double>	edge_origin(edge_width, edge_width, 0), shift;
	Vector3<long>	edge_size = pref->size() - (long) (2*edge_width);
	edge_size = edge_size.max(1);
	
	Bimage*			p = read_img(mg->fmg, 1, mg->img_num);

	if ( pgr ) p->multiply(pgr);
	
	p->calculate_background();
	
	if ( edge_width > 0 )
		p->shape(0, edge_size, edge_origin, gauss_width, FILL_BACKGROUND, 0);
	
	p->origin(mg_ref->origin);
	mg->pixel_size = mg_ref->pixel_size;
	p->sampling(pref->sampling(0));
	
	if ( bin > 1 ) p->bin(bin);
	
	shift = p->find_shift(pref, pmask, hi_res, lo_res, shift_limit, 0, 1, planf, planb, cc);
	if ( bin > 1 ) shift *= obin;
	
	mg->origin = mg_ref->origin + shift;
	mg->fom = cc;
	
	delete p;
	
	mg->scale = Vector3<float>(1,1,1);
	mg->matrix = Matrix3(1);
	mg_apply_transform(mg_ref, mg);

	return mg->fom;
}

double		field_align_series(Bfield* field, int ref_img, Bimage* pgr, Bimage* pmask,
				Vector3<double> origin, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width,
				long bin, Bstring& subset, int flag)
{
	long			i, j, nmg, nsel(0);
	Bmicrograph*	mg, *mg_ref;
	
	for ( i=0, mg_ref = field->mg; mg_ref && ( i<ref_img ); mg_ref = mg_ref->next, i++ ) ;
	if ( !mg_ref ) {
		mg_ref = field->mg;
		ref_img = 0;
	}

	for ( nmg=0, mg = field->mg; mg; mg = mg->next ) nmg++;
	vector<int>		numsel = select_numbers(subset, nmg);
	
	for ( i=nsel=0; i<nmg; i++ ) if ( numsel[i] ) nsel++;
	
	if ( nsel < 2 ) {
		cerr << "Error: There are less than two micrographs in the field-of-view!" << endl;
		return 0;
	}
	
	Bmicrograph**	mg_arr = new Bmicrograph*[nsel];

	for ( i=j=0, mg = field->mg; mg; mg = mg->next, i++ )
		if ( numsel[i] ) mg_arr[j++] = mg;
	
	Bimage*			pref = read_img(mg_ref->fmg, 1, mg_ref->img_num);
	if ( origin.length() < 1 ) origin = pref->size()/2;
	pref->origin(origin);
	mg_ref->origin = origin;
	if ( !mg_ref->pixel_size[0] ) mg_ref->pixel_size = pref->sampling(0);
	pref->sampling(mg_ref->pixel_size);
	
	Vector3<double>	edge_origin(edge_width, edge_width, 0);
	Vector3<long>	edge_size = pref->size() - (long) (2*edge_width);
	edge_size = edge_size.max(1);
	
	pref->calculate_background();
	
	if ( edge_width > 0 )
		pref->edge(0, edge_size, edge_origin, gauss_width, FILL_BACKGROUND);

	if ( bin > 1 ) pref->bin(bin);
	
	if ( shift_limit < 0 ) shift_limit = pref->sizeX()/10;

	if ( verbose & VERB_RESULT )
		cout << "Field " << field->id << endl << "Micrograph\tdx\tdy\tdz\tCC" << endl;

	fft_plan		planf = pref->fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = pref->fft_setup(FFTW_BACKWARD, 0);
	
#ifdef HAVE_GCD
	dispatch_apply(nsel, dispatch_get_global_queue(0, 0), ^(size_t i){
		mg_align_series(mg_arr[i], mg_ref, pref, pgr, pmask, hi_res, lo_res,
			shift_limit, edge_width, gauss_width, bin, planf, planb);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nsel; i++ )
		mg_align_series(mg_arr[i], mg_ref, pref, pgr, pmask, hi_res, lo_res,
			shift_limit, edge_width, gauss_width, bin, planf, planb);
#endif

    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	delete[] mg_arr;
	delete pref;

	double			d, da(0);

	for ( nsel=i=0, mg = field->mg; mg; mg = mg->next, i++ ) if ( numsel[i] ) {
		if ( verbose & VERB_RESULT )
			cout << mg->id << tab <<
					mg->origin - mg_ref->origin << tab <<
					mg->fom << endl;
		if ( mg != mg_ref ) {
			d = mg->origin.distance(mg_ref->origin)/(fabs(i-ref_img));
			da += d;
			nsel++;
		}
	}
	
	if ( nsel ) da /= nsel;
	
	if ( verbose & VERB_RESULT )
		cout << "Average shift per micrograph:    " << da << endl;

	field->fom = da;

	if ( flag & 4 ) {
		Bstring			ext = mg->fmg.extension();
		if ( ext.contains("dm") ) ext = "mrc";
		if ( ext.contains("tif") ) ext = "mrc";
		DataType		dt(Float);
		Bstring			label("Written by bseries");
		Bstring			favg = field->mg->fmg.base() + "_avg." + ext;
		field_write_aligned_average(field, pgr, mg->fmg, dt, subset);
	}
	
	return da;
}

/**
@brief 	Aligns a series of micrographs by cross-correlation.
@param 	*project		project parameter structure.
@param 	ref_img			reference micrograph number (starts from 0).
@param 	*pgr			gain reference.
@param 	*pmask			reciprocal space mask, 0's and 1's.
@param 	origin			tilt origin.
@param 	hi_res			high resolution limit.
@param 	lo_res			low resolution limit.
@param 	shift_limit		maximum shift from nominal origin of image.
@param 	edge_width		edge smoothing width (not done if 0).
@param 	gauss_width		edge decay width.
@param 	bin				3-value vector of integer bin factors.
@param 	&subset			subset to average (all if empty)
@param 	flag			options flag.
@return double			root-mean-square of offsets.

	Each micrograph in the series is cross-correlated with the reference
	micrograph and the shift determined.
	Options encoded in the flag:
	1	rescale image based on histogram.
	2	weigh by accumulated dose.
	4	write aligned frames with insert "_aln".
	8	write aligned frame sum with insert "_avg".

**/
double		project_align_series(Bproject* project, int ref_img, Bimage* pgr,
				Bimage* pmask, Vector3<double> origin, double hi_res, double lo_res,
				double shift_limit, double edge_width, double gauss_width,
				long bin, Bstring& subset, int flag)
{
	Bfield*			field = project->field;
	long			nfield(0);
	long			nmg = project_count_micrographs(project);

	for ( field = project->field; field; field = field->next ) nfield++;
	

//	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
	if ( verbose ) {
		cout << "Aligning micrographs in series by cross-correlation:" << endl;
		cout << "Number of images:               " << nmg << endl;
		cout << "Resolution limits:              " << hi_res << " - " << lo_res << " A" << endl;
		cout << "Shift limit:                    " << shift_limit << endl;
		cout << "Edge masking width & smoothing: " << edge_width << " " << gauss_width << endl;
		cout << "Binning:                        " << bin << endl;
		if ( pmask )
			cout << "Using mask file:                " << pmask->file_name() << endl;
		if ( subset.length() )
			cout << "Subset:                         " << subset << endl;
		cout << endl;
	}

	double			d(0);
	for ( field = project->field; field; field = field->next ) {
		field_align_series(field, ref_img, pgr, pmask, origin, hi_res, lo_res,
			shift_limit, edge_width, gauss_width, bin, subset, flag);
		d += field->fom;
	}
	
	d /= nfield;

	if ( verbose & VERB_RESULT )
		cout << "Overall average shift per frame: " << d << endl << endl;
	
	return d;
}

Bimage*		mg_tomo_reconstruct(Bmicrograph** mgarr, long mg_id, long mg_min, long mg_max,
				double hi_res, double scale, Vector3<long> size,
				long ft_size, fft_plan planf, fft_plan planb)
{
	long   			i;
	Bimage* 		p = NULL;

	Bmicrograph*	mg_res = mgarr[mg_id];
	Bmicrograph*	mg;
	
	if ( !mg_res ) {
		error_show("Error in mg_tomo_reconstruct", __FILE__, __LINE__);
		cerr << "Micrograph number " << mg_id << " not found!" << endl << endl;
		return NULL;
	}

	if ( !mg_res->select ) {
		if ( verbose & VERB_FULL )
			cout << "Skipping micrograph number " << mg_id << endl << endl;
		return NULL;
	}
	
	if ( mg_max < mg_min || ( mg_max == mg_id && mg_min == mg_id ) ) {
		error_show("Error in mg_tomo_reconstruct", __FILE__, __LINE__);
		cerr << "Micrograph range is incorrect: " << mg_min << " - " << mg_max << endl << endl;
		return NULL;
	}
	
	if ( mg_res->fmg.length() < 0 ) {
		error_show("Error in mg_tomo_reconstruct", __FILE__, __LINE__);
		cerr << "No micrograph image file name given for micrograph " << mg_res->id << endl;
		return NULL;
	}
	
	if ( hi_res < 2*mg_res->pixel_size[0] ) hi_res = 2*mg_res->pixel_size[0];
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mg_tomo_reconstruct: size=" << size << " ft_size=" << ft_size << endl;

	// Header structure for the final reconstruction
	Bimage* 		prec = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	prec->sampling(mg_res->pixel_size);
	prec->origin(prec->size()/2);
	prec->fourier_type(Standard);
	
	prec->next = new Bimage(Float, TSimple, prec->size(), prec->images());

	long   			data_size = (long) prec->size().volume();
	
	long 			nrec(0), nmg(mg_max-mg_min), nover(0);
//	long			interp_type(0), pad_factor(0);
//	double			pad_ratio(ft_size*1.0/prec->sizeX());
	Vector3<long>	translate;
	
	for ( i=mg_min; i<=mg_max; i++ ) if ( mgarr[i]->select && i != mg_id ) {
		mg = mgarr[i];
		if ( verbose & VERB_FULL )
			cout << "Reading image " << mg->img_num << " (micrograph " << mg->id << ")" << endl;
		if ( ( p = read_img(mg->fmg, 1, mg->img_num) ) == NULL ) {
			error_show("mg_tomo_reconstruct", __FILE__, __LINE__);
			return NULL;
		}
		
		if ( mg->origin.length() < 1 )
			mg->origin = p->size()/2;
	
		translate = (prec->size() - p->size())/2;
//		cout << translate << endl;
		
		p->change_type(Float);
		p->image->origin(mg->origin[0], mg->origin[1], 0);
		p->sampling(mg->pixel_size);
		p->statistics();
//		p->rescale_to_avg_std(0,1);
		p->resize(prec->size(), translate, FILL_AVERAGE, 0);
	

		p->fft(planf);
		p->phase_shift_to_origin();
		
//		cout << "preparation for reconstruction done" << endl;
		
		nover += prec->fspace_pack_2D_into_central_section(p, ft_size, scale, hi_res, 0, mg->matrix, mg_res->matrix);
		
		delete p;
					
		nrec++;
					
		if ( verbose & ( VERB_TIME | VERB_PROCESS ) )
			cout << "Complete:                       " << nrec*100.0/nmg << "%\r" << flush;
		
	}
	
	if ( nrec < 1 ) {
		delete prec;
		error_show("Error in mg_tomo_reconstruct: No micrographs used in this reconstruction!", __FILE__, __LINE__);
		return NULL;
	}

	float*			fom = (float *) prec->next->data_pointer();

	if ( verbose & VERB_FULL ) {
		cout << "Contributions:                    " << nover << " (" << nover*100.0/data_size << "%)" << endl;
		cout << "Weighing reconstruction.          " << endl;
	}
	
	for ( i=0; i<data_size; i++ ) {
		if ( fom[i] ) prec->set(i, prec->complex(i) / fom[i]);
		fom[i] = fom[i]/(fom[i] + 1);
	}

	prec->set(0, 0);
	prec->origin(0, 0, 0);
	prec->phase_shift_to_center();
	prec->fft_back(planb);
	prec->statistics();
	
	return prec;
}

/**
@brief 	Reconstructs the zero-tilt image of a tilt-series.
@param 	*project		project parameter structure.
@param 	dimg			number of adjacent images in reconstructions.
@param 	size			reconstruction size.
@param 	scale			reconstruction scale.
@param 	hi_res			high resolution limit (angstrom).
@return Bimage*			new 2D image.

**/
Bimage*		mg_tomo_reconstruct2D(Bproject* project, long dimg,
				Vector3<long> size, double scale, double hi_res)
{
	if ( scale < 0.1 ) scale = 1;
	
	long			i, nmg(0);
	
	Bmicrograph*	mg0 = field_find_zero_tilt_mg(project->field);
	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);

	if ( hi_res < 2*mg0->pixel_size[0] )
		hi_res = 4*mg0->pixel_size[0];
	
	if ( dimg < 1 ) {
		for ( i=0; i<nmg; ++i )
			if ( fabs(mgarr[i]->tilt_angle) <= 15*M_PI/180.0 )
				break;
		dimg = abs(mg0->img_num - mgarr[i]->img_num);
	}
	
	Bimage*			p = read_img(mg0->fmg, 0, mg0->img_num);
	if ( size.volume() < 1 ) size = p->size()/scale;
	size[2] = 1;
	
	long 	ft_size = p->sizeX();
	if ( ft_size < p->sizeY() ) ft_size = p->sizeY();
	ft_size = findNextPowerOf(ft_size, 2);
	
	delete p;

//	size = Vector3<long>(ft_size,ft_size,1);

	if ( verbose ) {
		cout << "Reconstructing a 2D image:" << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Reconstruction size:            " << size << endl;
		cout << "Reconstruction scale:           " << scale << endl;
		cout << "Number of adjacent images:      " << dimg << endl;
		cout << "Resolution limit:               " << hi_res << " A" << endl;
		cout << endl;
	}

//	fft_plan		planf = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
//	fft_plan		planb = fft_setup_plan(ft_size, ft_size, 1, FFTW_BACKWARD, 1);
	fft_plan		planf = fft_setup_plan(size, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(size, FFTW_BACKWARD, 1);
	
	Bimage*			prec = mg_tomo_reconstruct(mgarr, mg0->img_num,
						mg0->img_num - dimg, mg0->img_num + dimg,
						hi_res, scale, size, ft_size, planf, planb);
	
	delete[] mgarr;
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	return prec;
}

double		mg_reconstruct_align(Bmicrograph** mgarr, long nmg,
				long mg_id, long dimg, double resolution, double shift_limit,
				Vector3<long> size, double edge_width, double gauss_width,
				long ft_size, fft_plan planf, fft_plan planb)
{
	long			imin, imax, thickness(size[2]);
	double			cc;
	Bmicrograph*	mg_res = mgarr[mg_id];
	
	if ( shift_limit < 0 ) shift_limit = size[0]/4;
	
	for ( imin=(mg_id>dimg)? mg_id-dimg: 0; mgarr[imin]->fom < -0.5 && imin<=mg_id; imin++ ) ;
	
	for ( imax=(mg_id+dimg<nmg)? mg_id+dimg: nmg-1; mgarr[imax]->fom < -0.5 && imax>=mg_id; imax-- ) ;
	
//	cout << mg_id << ": " << imin << " - " << imax << endl;

	Bimage*			prec = mg_tomo_reconstruct(mgarr, mg_id, imin, imax,
						resolution, 1, size, ft_size, planf, planb);

	if ( !prec ) {
		cerr << "Error: Reconstruction for " << mg_id << " failed!" << endl;
		bexit(-1);
	}
	
//	cout << "reconstruction origin = " << prec->image->origin() << endl;
	
//	cout << mg_id << ": " << imin << " - " << imax << endl;


	Bimage*			p = read_img(mg_res->fmg, 1, mg_res->img_num);
	if ( mg_res->origin.length() < 1 ) mg_res->origin = p->size()/2;
	p->image->origin(mg_res->origin);
	p->image->sampling(mg_res->pixel_size);
	
	// Resize to a smaller sub-images
	Vector3<long>	translate = (prec->size() - p->size())/2;
	p->resize(prec->size(), translate, FILL_AVERAGE, 0);

	Vector3<double>	edge_origin(edge_width, edge_width, 0);
	Vector3<double>	edge_size(p->sizeX() - (long)edge_width, p->sizeY() - (long)edge_width, 1);
	if ( edge_width > 0 ) {
		if ( thickness > 0 ) {
			img_clear_extraneous_areas(p, mg_res->tilt_axis, mg_res->tilt_angle,
				thickness, edge_width);
			img_clear_extraneous_areas(prec, mg_res->tilt_axis, mg_res->tilt_angle,
				thickness, edge_width);
		} else {
			p->edge(0, edge_size, edge_origin, gauss_width, FILL_AVERAGE, 0);
			prec->edge(0, edge_size, edge_origin, gauss_width, FILL_AVERAGE, 0);
		}
	}

	p->image->origin(prec->image->origin());
	p->sampling(prec->sampling(0));
	
	Vector3<double>	shift = p->find_shift(prec, NULL, resolution, 0, shift_limit, 0, 1, planf, planb, cc);
//	cout << "origin = " << p->image->origin() << endl;
//	cout << "mg origin = " << mg_res->origin << endl;
	
	// Change in origin from previous one
	Vector3<double>		oldori = mg_res->origin;
	shift -= translate;
	mg_res->origin = prec->image->origin() + shift;
	mg_res->fom = cc;
	shift = mg_res->origin - oldori;
	
	if ( verbose )
		cout << mg_id << ": " << imin << " - " << imax << tab <<
			mg_res->origin << tab << shift << tab << mg_res->fom << endl;
/*
	long			k(19);
	if ( mg_id == k ) {
		Bstring			name(k, "r%d.pif");
		write_img(name, prec, 0);
		name = Bstring(k, "t%d.pif");
		write_img(name, p, 0);
		cout << "writing " << name << endl;
		bexit(0);
	}
*/
	
	delete p;
	delete prec;

	return shift.length();
}


/**
@brief 	Aligns a series of micrographs by cross-correlation sequentially.
@param 	*project		project parameter structure.
@param 	thickness		reconstruction thickness.
@param 	iter			number of alignment iterations.
@param 	dchange			threshold change in origin.
@param 	dimg			number of adjacent images in reconstructions.
@param 	resolution		high resolution limit (angstrom).
@param 	shift_limit		limit on shift search (pixels).
@param 	edge_width		smoothing edge width.
@param 	gauss_width		smoothing edge decay.
@return long				0, <0 on failure.

	Each pair of adjacent micrographs in the series is cross-correlated
	and the relative shift determined. The relative shifts are adjusted
	relative to the reference micrograph, defined as the one closest
	to a zero degree tilt.
	The images are stretched to compensate for tilt difference.
	The relationship between an euler representation of the view and
	the tilt axis and tilt angle is:
		tilt_axis = phi - 90 = - psi - 90
		tilt_angle = theta

**/
long			project_tomo_align(Bproject* project, long thickness, long iter, double dchange,
				long dimg, double resolution, double shift_limit, double edge_width, double gauss_width)
{
	Bmicrograph*	mg0 = field_find_zero_tilt_mg(project->field);

	if ( !mg0->fmg.length() ) {
		error_show("Error in project_tomo_align ", __FILE__, __LINE__);
		cerr << "No micrograph image or transform file name given for micrograph "
			<< mg0->id << endl << endl;
		return -1;
	}

	Bimage*			p = read_img(mg0->fmg, 0, mg0->img_num);

//	long 	ft_size = p->sizeX();
//	if ( ft_size < p->sizeY() ) ft_size = p->sizeY();
//	ft_size = findNextPowerOf(ft_size, 2);
	
	// Select a sub-image fitting into the frame with a power of 2
	long 	ft_size = p->sizeX();
	if ( ft_size > p->sizeY() ) ft_size = p->sizeY();
	ft_size = findNextPowerOf(ft_size/2, 2);
	
	if ( thickness < 10 ) thickness = ft_size/5;
	
	delete p;

	if ( resolution < 2*mg0->pixel_size[0] ) resolution = 4*mg0->pixel_size[0];
	
	long			it, nmg(0), i, j, n, mg_id0(0);
	double			dori(10);
	Vector3<long>	size(ft_size,ft_size,thickness);

	if ( shift_limit < 0 ) shift_limit = ft_size/4;
	
	Bmicrograph**	mgarr = project_micrograph_array(project, nmg);
	for ( i=0; i<nmg; i++ ) {
		mgarr[i]->fom = -1;
		if ( mg0 == mgarr[i] ) mg_id0 = i;
	}
	mg0->fom = 1;

//	ft_size = part_ft_size(size[0], 1, 2);
	
	if ( verbose ) {
		cout << "Aligning images in a tilt series by cross-correlation:" << endl;
		cout << "Maximum number of iterations:   " << iter << endl;
		cout << "Number of images:               " << nmg << endl;
		cout << "Fourier transform size:         " << ft_size << " x " << ft_size << endl;
		cout << "Reconstruction size:            " << size << endl;
		cout << "Number of adjacent images:      " << dimg << endl;
		cout << "Resolution limit:               " << resolution << " A" << endl;
		cout << "Shift limit:                    " << shift_limit << " pixels" << endl;
		cout << "Smoothing edge and width:       " << edge_width << " " << gauss_width << " pixels" << endl;
		cout << endl;
	}

	fft_plan		planf = fft_setup_plan(ft_size, ft_size, 1, FFTW_FORWARD, 1);
	fft_plan		planb = fft_setup_plan(ft_size, ft_size, 1, FFTW_BACKWARD, 1);
	
	for ( it = 0; it < iter && dori > dchange; it++ ) {
		if ( verbose ) {
			cout << "Iteration " << it+1 << ":" << endl;
			cout << "Micrograph: range" << tab << "Origin" << tab << "Shift" << tab << "FOM" << endl;
		}
		n = 0;
		dori = 0;
		for ( i=mg_id0-1, j=mg_id0+1; i>=0 || j<nmg; i--, j++ ) {
			if ( i >= 0 ) {
				dori += mg_reconstruct_align(mgarr, nmg, i, dimg, resolution,
					shift_limit, size, edge_width, gauss_width, ft_size, planf, planb);
				n++;
			}
			if ( j < nmg ) {
				dori += mg_reconstruct_align(mgarr, nmg, j, dimg, resolution,
					shift_limit, size, edge_width, gauss_width, ft_size, planf, planb);
				n++;
			}
		}
		dori /= n;
		if ( verbose )
			cout << "Iteration:\t" << it+1 << tab << dori << endl;
	}
	
	delete[] mgarr;
    fft_destroy_plan(planf);
    fft_destroy_plan(planb);

	Bmicrograph*	mg;
	if ( verbose & VERB_RESULT ) {
		cout << "\nImage\tdx\tdy\tdz" << endl;
		for ( mg = project->field->mg; mg; mg = mg->next )
			cout << mg->id << tab << mg->origin[0] << tab << mg->origin[1]
				<< tab << mg->origin[2] << endl;
		cout << endl;
	}
	
	return 0;
}



double		field_ssnr(Bfield* field, fft_plan plan)
{
	long			i;
	Bmicrograph*	mg;
	Vector3<double>	origin;
	Bimage* 		p, *ps;
	Bimage* 		pavg = NULL;
	Bimage* 		pstd = NULL;
	Vector3<double> scale(1,1,1), nori;
	Matrix3 		mat(1);
	Vector3<long>	size;
	
	for ( i=0, mg = field->mg; mg; mg = mg->next, i++ ) {
		p = read_img(mg->fmg, 1, mg->img_num);
		p->calculate_background();
		origin = p->size()/2;
		if ( size.volume() < 1 ) {	// Make sure it is square
			size = p->size();
			if ( size[0] < size[1] ) size[1] = size[0];
			else size[0] = size[1];
		}
		ps = p->transform(size, scale, nori, (origin - mg->origin), mat, FILL_BACKGROUND, 0);
		ps->fft(plan);
		if ( !pavg ) pavg = ps->copy();
		else pavg->add(ps);
		ps->complex_to_intensities();
		if ( !pstd ) pstd = ps->copy();
		else pstd->add(ps);
		delete p;
		delete ps;
	}
	
	if ( i ) {
		pavg->multiply(1.0L/i);
		pstd->multiply(1.0L/i);
	}
	
	pavg->complex_to_intensities();
	pstd->subtract(pavg);
	pstd->truncate_to_min_max(1e-10, pstd->maximum());
	pavg->divide(pstd, 1, 0);
	
	field->fom = pavg->average();

	Bstring		imgfile = field->id + "_snr.map";
	write_img(imgfile, pavg, 0);

	double 		rad_step(10), maxrad(size[0]/2);
	Bimage* 	prad = pavg->radial(0, maxrad, rad_step, 1);
	imgfile = field->id + "_snr.txt";
	write_img(imgfile, prad, 0);

	delete pavg;
	delete pstd;
	delete prad;
	
	return 0;
}

double		project_ssnr(Bproject* project)
{
	Bfield*			field = project->field;
	long			n, nfield(0);

	for ( field = project->field; field; field = field->next ) nfield++;
	
	Bfield**		farr = new Bfield*[nfield];
	for ( n=0, field = project->field; field; field = field->next, n++ )
		farr[n] = field;

	if ( verbose )
		cout << "Calculating the SSNR for each field:" << endl;
		
	Bimage*			p = read_img(project->field->mg->fmg, 0, 0);

	fft_plan		plan = p->fft_setup(FFTW_FORWARD, 0);
	
	delete p;
	
//#ifdef HAVE_GCD
//	dispatch_apply(nfield, dispatch_get_global_queue(0, 0), ^(size_t i){
//		field_ssnr(farr[i], plan);
//	});
//#else
//#pragma omp parallel for
	for ( long i=0; i<nfield; i++ ) {
		field_ssnr(farr[i], plan);
	}
//#endif

    fft_destroy_plan(plan);
	delete[] farr;

	if ( verbose ) {
		cout << "Field\tSNR" << endl;
		for ( field = project->field; field; field = field->next )
			cout << field->id << tab << field->fom << endl;
	}
	
	return 0;
}

double		progressive_snr_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); ++i ) {
		df = f[i] - simp.parameter(0)*(1 - exp(-simp.parameter(1)*x[i]));
		R += df*df;
	}
	
	R /= simp.points();
	R /= simp.dependent_variance();
	
	return R;
}

// y = (snr*cd)*(1-exp(-x/cd)
double 		fit_progressive_snr(vector<double>& x, vector<double>& y, double& snr, double& cd )
{
	Bsimplex			simp(1, 2, 0, x.size(), x, y);
	
	snr = fabs(y.back());
//	cd = fabs(snr*x[0]/y[0]);
	cd = 0.1;
	
//	cout << snr << tab << cd << endl;

	simp.parameter(0, snr);
	simp.parameter(1, cd);
	simp.limits(0, 1e-6, 10*snr);
	simp.limits(1, 1e-6, 10*cd);
	
	double			R = simp.run(10000, 1e-5, progressive_snr_R);
	
	snr = simp.parameter(0)*simp.parameter(1);
	cd = 1/simp.parameter(1);
	
	return R;
}
/*
double		fit_dose_tolerance(vector<double>& s, vector<double>& c, double f);
{
	Bsimplex			simp(1, 1, 0, c.size(), s, c);
	
	snr = fabs(y.back());
	cd = fabs(snr*x[0]/y[0]);
	
//	cout << snr << tab << cd << endl;

	simp.parameter(0, snr);
	simp.parameter(1, 1/cd);
	simp.limits(0, 1e-6, 5*snr);
	simp.limits(1, 0.2/cd, 5/cd);
	
	double			R = simp.run(10000, 1e-5, dose_tolerance_R);
	
	snr = simp.parameter(0)*simp.parameter(1);
	cd = 1/simp.parameter(1);
	
	return R;
}
*/

#include <fstream>
Bstring		progsnr_file, pssnr_curves;

double		fit_individual_progressive_snr(Bplot* plot, vector<double>& minima, long nimg, double dose_per_subset)
{
	long			i, j, jj, ii, nc(0);
	double			nv, mx, snr, cd, c, ca(0), cv(0), R;
	vector<double>	x(nimg,0);
	vector<double>	v(nimg,0);
//	vector<double>	s, c;

    ofstream		ftxt(progsnr_file.str());
    if ( ftxt.fail() ) return -1;
    ofstream		fcurves(pssnr_curves.str());
    if ( fcurves.fail() ) return -1;
    
    vector<vector<double>>	curves;

	for ( i=0; i<nimg; ++i )
		x[i] = dose_per_subset*(i+1);
		
	curves.push_back(x);

	ftxt << "s\tResolution\tSNR0\tCd\t1/Cd\tR\t#" << setprecision(3) << endl;
	if ( verbose )
		cout << "s\tResolution\tSNR0\tCd\tR\t#\tc" << setprecision(3) << endl;
	i = jj = 1;
	mx = 0;
	for ( auto m: minima ) {
		mx = (mx + m)/2;
		for ( ; jj<plot->rows(); ++jj ) if ( (*plot)[jj] > m ) break;
		for ( nv=0; i<jj; ++i ) {
			for ( j=0, ii=i+plot->rows(); j<v.size(); ++j, ii+=plot->rows() )
				v[j] += (*plot)[ii];
			nv += 1;
		}
		if ( nv ) {
			for ( j=0; j<v.size(); ++j ) v[j] /= nv;
			R = fit_progressive_snr(x, v, snr, cd);
			curves.push_back(v);
//			for ( j=0; j<v.size(); ++j ) {
//				snrc = snr*cd*(1-exp(-x[j]/cd));
//				cout << x[j] << tab << v[j] << tab << snrc << endl;
//			}
//			if ( v.back() < 0.01 || R > 0.2 ) continue;
			c = 1/(mx*cd);
			if ( snr > 0.01 && R < 0.1 ) {
				ca += c;
				cv += c*c;
				nc++;
//				s.push_back(mx);
//				c.push_back(cd);
			}
			ftxt << mx << tab << 1/mx << tab << snr << tab << cd << tab << 1/cd << tab << R << tab << nv << endl;
			if ( verbose )
				cout << mx << tab << 1/mx << tab << snr << tab << cd << tab << R << tab << nv << tab << c << endl;
		}
		mx = m;
	}
	
	ca /= nc;
	cv /= nc;
	cv -= ca*ca;
	if ( cv > 0 ) cv = sqrt(cv);
	else cv = 0;
	if ( verbose )
		cout << "Coefficient: " << ca << " (" << cv << ", " << nc << ")" << endl << endl;
	
/*
	if ( c.size() ) {
		R = fit_dose_tolerance(s, c, f);
		if ( verbose )
			cout << "Dose tolrance coefficient: " << f << endl;
	}
*/

	fcurves << "Dose (e/A2)";
	for ( j=1; j<curves.size(); ++j ) fcurves << tab << j;
	fcurves << endl;
	for ( i=0; i<x.size(); ++i ) {
		fcurves << curves[0][i];
		for ( j=1; j<curves.size(); ++j )
			fcurves << tab << curves[j][i];
		fcurves << endl;
	}

	ftxt.close();
	fcurves.close();
	
	return R;
}


double		progressive_snr_R2(Bsimplex& simp)
{
	long			i, j;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0, j=f.size(); i<simp.points(); ++i, ++j ) {
		df = f[i] - (simp.parameter(0)/x[i])*(1 - exp(-simp.parameter(1)*x[i]*x[j]));
		R += df*df;
	}
	
	R /= simp.points();
	R /= simp.dependent_variance();
	
	return R;
}

// y = (snr/(cs))*(1-exp(-csd))
double 		fit_progressive_snr2(vector<double>& x, vector<double>& y, double& snr, double& c )
{
	Bsimplex			simp(2, 2, 0, y.size(), x, y);
	
	snr = 0.1;	// SNR0/c
	c = 2;
	
//	cout << snr << tab << cd << endl;

	simp.parameter(0, snr);
	simp.parameter(1, c);
	simp.limits(0, 1e-6, 10*snr);
	simp.limits(1, 0.1*c, 10*c);
	
	double			R = simp.run(10000, 1e-5, progressive_snr_R2);
	
	snr = simp.parameter(0)*simp.parameter(1);
	c = simp.parameter(1);
	
	return R;
}

double		fit_full_progressive_snr(Bplot* plot, vector<double>& minima, long nimg, double dose_per_subset)
{
	long			i, j, k, k1, ii, jj, kk, ip;
	double			snr, snrc, cd, R;
	vector<double>	x(2*nimg*minima.size(),0);
	vector<double>	v(nimg*minima.size(),0);

	for ( i=k=0, k1=v.size(), jj=1; i<minima.size(); ++i ) {
		for ( ii = jj; jj<plot->rows(); ++jj ) if ( (*plot)[jj] > minima[i] ) break;
//		cout << ii << tab << jj << tab << jj-ii << endl;
		for ( j=0; j<nimg; ++j, ++k, ++k1 ) {
			if ( i ) x[k] = (minima[i] + minima[i-1])/2;
			else x[k] = minima[0]/2;
			x[k1] = dose_per_subset*(j+1);
//			ip=ii+(j+1)*plot->rows();
			for ( kk=ii, ip=ii+(j+1)*plot->rows(); kk<jj; ++kk, ++ip )
				v[k] += (*plot)[ip];
			v[k] /= jj-ii;
			if ( v[k] < 0 ) v[k] = 0;
//			if ( verbose )
//				cout << x[k] << tab << 1/x[k] << tab << x[k1] << tab << v[k] << endl;
		}
	}

	R = fit_progressive_snr2(x, v, snr, cd);
	
	if ( verbose ) {
		cout << "SNR0 = " << snr << endl;
		cout << "c = " << cd << endl;
		cout << "R = " << R << endl;
		cout << "s\tResolution\tDose\tSNR]\tSNRcalc" << setprecision(3) << endl;
		for ( i=k=0, k1=v.size(); i<minima.size(); ++i ) {
			for ( j=0; j<nimg; ++j, ++k, ++k1 ) {
				snrc = snr/(cd*x[k])*(1-exp(-cd*x[k]*x[k1]));
				cout << x[k] << tab << 1/x[k] << tab << x[k1] << tab << v[k] << tab << snrc << endl;
			}
		}
	}

	return R;
}

int		img_filter_spikes(Bimage* p, double ratio)
{
	long			i, ns(0), ds(p->image_size()*p->images());
	double			v, va;
	Complex<double>	cv;
	
    for ( i=0; i<ds; i++ ) {
		va = p->kernel_neighbor_average(i, 1);
		if ( va < 0 ) cerr << "Error: " << i << tab << va << endl;
		if ( va <= 0 ) va = 1e-6;
		cv = p->complex(i);
		v = cv.power();
		if ( v > ratio*va ) {
			p->set(i, cv*va/v);
			ns++;
		}
	}
	
	if ( verbose )
		cout << "Spikes removed: " << ns << " (" << ns*100.0L/ds << " %)" << endl;
	
	return ns;
}


int			mg_frames_snr(Bmicrograph* mg, double res_hi, long window, Bstring& subset, double sampling_ratio, int flag)
{
	if ( !mg->frame ) return 0;
	
	long			n;
	Bframe*			frame;

	Bimage*			p = read_img(mg->fframe, 0, -1);

//	long			memreq = p->images()*(1.0+1.0/window)*p->size().volume()*sizeof(Complex<float>);
//	memory_check(memreq);
	
	delete p;

	p = read_img(mg->fframe, 1, -1);
	p->sampling(mg->frame_pixel_size);
	if ( p->sizeZ() > 1 ) p->slices_to_images();
	
	if ( res_hi < p->image->sampling()[0] * 2 ) res_hi = p->image->sampling()[0] * 2;

	if ( flag & 1 ) p->histogram_counts(2);

	long			nimg = p->set_subset_selection(subset);
	
	if ( verbose )
		cout << "Micrograph: " << mg->id << " with " << nimg << " frames at " << mg->dose/p->images() << " e/Å2/frame" << endl;

	for ( n=0, frame = mg->frame; frame; frame = frame->next, ++n )
		p->image[n].origin(frame->shift + p->size()/2);

	Bimage*			psum = p->fspace_subset_sums(window, 1|(flag&2));
	
	delete p;
	
	(*psum)["dose"] = mg->dose;
	(*psum)["exposure"] = mg->exposure;
	
//	psum->sum_images();
	
	if ( flag & 4 ) img_filter_spikes(psum, 100);
	
//	write_img("tt.mrc", psum, 0);
	
//	bexit(0);
	
	Bplot*			plot = psum->fspace_subset_ssnr(window, res_hi, sampling_ratio, (flag&2));

	Bstring			psfile = mg->fframe.base() + "_subssnr.ps";
	
	ps_plot(psfile, plot);
	
	delete plot;

	psum->progressive_sum();
	psum->next->progressive_sum();

	plot = psum->fspace_subset_ssnr(window, res_hi, sampling_ratio, 1);

	psfile = mg->fframe.base() + "_progssnr.ps";
	
	ps_plot(psfile, plot);
	
	long			i, j;
	vector<double>	v(plot->rows());
	for ( i=0, j=plot->rows()*(plot->columns()-1); i<plot->rows(); ++i, ++j ) v[i] = (*plot)[j];
	
	if ( !mg->ctf ) mg->ctf = new CTFparam;
	if ( mg->ctf->defocus_average() < 10 ) {
		double			step_size = sampling_ratio/psum->real_size()[0];
		long			rmin = (long) (1/(30*step_size));
		long			rmax = (long) (1/(2*res_hi*step_size));
		ctf_find_defocus(v, *(mg->ctf), rmin, rmax, step_size, 1e3, 2e5, 1e3);
	}

	if ( verbose )
		cout << "Defocus:                        " << mg->ctf->defocus_average()*1e-4 << " um" << endl;
	
	vector<double>	minima = mg->ctf->zeroes(0.8/res_hi);
	
	double			dose_per_subset(mg->dose/psum->images());

	if ( verbose )
		cout << "Dose per subset = " << dose_per_subset << endl;
	if ( dose_per_subset <= 0 ) {
		cerr << "Error: The dose is not specified!" << endl;
		bexit(-1);
	}

	progsnr_file = mg->fframe.base() + "_progssnr.txt";
	pssnr_curves = mg->fframe.base() + "_pssnr_curves.txt";
	fit_individual_progressive_snr(plot, minima, psum->images(), dose_per_subset);

//	fit_full_progressive_snr(plot, minima, psum->images(), dose_per_subset);

	delete psum;
	delete plot;

	return 0;
}

/**
@brief 	Calculates the estimated SNR from aligned movie frames.
@param 	*project		project parameter structure.
@param 	res_hi			high resolution limit.
@param 	window			number of frames to sum for each curve.
@param 	&subset			subset to average (all if empty).
@param 	sampling_ratio	radial sampling ratio (1 or larger).
@param 	flag			flag to convert to counts (1) and calculate a progressive sum (2), filter extremes (4).
@return int				0, <0 on failure.

	The SNR is calculated as for reconstructions.

**/
int			project_frames_snr(Bproject* project, double res_hi,
				long window, Bstring& subset, double sampling_ratio, int flag)
{
	if ( sampling_ratio < 1 ) sampling_ratio = 1;
	
	Bfield*			field = NULL;
	Bmicrograph*	mg = NULL;
	
	if ( verbose ) {
		cout << "Calculating the estimated SSNR:" << endl;
		cout << "High resolution limit:          " << res_hi << " A" << endl;
		cout << "Summing window:                 " << window << endl;
		cout << "Sampling ratio:                 " << sampling_ratio << endl;
		cout << endl;
	}
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_frames_snr(mg, res_hi, window, subset, sampling_ratio, flag);
	
	return 0;
}

int			mg_frame_shift_analysis(Bmicrograph* mg, long window, double resolution)
{
	if ( !mg->frame ) return 0;
	if ( resolution < 1 ) resolution = 10;

	long			i, j, k, nf(0);
	double			s(M_PI*mg->frame_pixel_size[0]/resolution), E;
	vector<double>	csum;
	Vector3<double>	pshift;
	Bframe*			frame;
	
	for ( nf=j=0, frame = mg->frame; frame; frame = frame->next, ++nf ) {
		if ( nf ) {
			if ( nf%window == 0 ) {
				csum.push_back(0);
				j++;
			}
			csum[j] += frame->shift.distance(pshift);
		} else {
			csum.push_back(frame->shift.distance(frame->next->shift));
		}
		pshift = frame->shift;
	}
	
	cout << "Micrograph: " << mg->id << endl;
	cout << "Frames\tShift\tE(s=" << 1/resolution << ")" << endl;
	for ( j=0, i=1, k=window; j<csum.size(); ++j, i+=window, k+=window ) {
		if ( k > nf ) k = nf;
		E = sin(s*csum[j])/(s*csum[j]);
		E = 0.5*(1+E*E);
		cout << i << "-" << k << tab << csum[j] << tab << E << endl;
	}

	return 0;
}

/**
@brief 	Calculates the estimated SNR from aligned movie frames.
@param 	*project		project parameter structure.
@param 	window			number of frames to sum for each curve.
@param 	resolution		resolution for calculating envelope.
@return int				0, <0 on failure.

	The SNR is calculated as for reconstructions.

**/
int			project_frame_shift_analysis(Bproject* project, long window, double resolution)
{
	if ( window < 1 ) window = 1;
	if ( resolution < 1 ) resolution = 10;
	
	Bfield*			field = NULL;
	Bmicrograph*	mg = NULL;
	
	if ( verbose ) {
		cout << "Frame shift analysis:" << endl;
		cout << "Summing window:                 " << window << endl;
		cout << "Resolution:                     " << resolution << " A" << endl;
		cout << endl;
	}
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			mg_frame_shift_analysis(mg, window, resolution);
	
	return 0;
}
