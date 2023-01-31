/**
@file	Bimage_segment.cpp
@brief	Methods for segmentation.
@author	Bernard Heymann and Samuel Payne
@date	Created: 20010516
@date	Modified: 20210118 (BH)
**/

#include "rwimg.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Segments an image into contiguous regions.
@param 	threshold 	the level at which things are ignored.
@param 	sign		sign controlling direction of thresholding.
@return Bimage* 	level mask with indexed regions.

	The image is segmented into contiguous regions where a region is 
	defined as all those voxels exceeding the given threshold and
	adjacent to each other.
	This method uses a list to keep track of voxels assigned to a region.
	In subsequent iterations the non-assigned voxels next to those in
	the list are assigned if they exceed the threshold, and their indices
	are kept in a new list. The new list is transferred to the old list
	and the process iteratively continued until the new list contains
	no more voxels. The indices of all the regions are packed into
	an integer data block within a new image.
	The new image maps the indices of the regions starting from one to as
	many regions as were found (the maximum gives the number of regions).
	The new image has zeroes for all voxels outside the regions.

**/
Bimage* 	Bimage::regions(double threshold, int sign)
{
	// The extreme considered as a region depends on the threshold and sign
	if ( sign < 0 ) sign = -1;
	else if ( sign > 0 ) sign = 1;
	else {
		sign = 1;
		if ( threshold < avg ) sign = -1;
	}
	
	long	   		i, count(0), rsize;
	double			abs_threshold(sign*threshold);

	if ( verbose & VERB_PROCESS ) {
		cout << "Finding regions:" << endl;
		cout << "Threshold:                      " << threshold << " (" << sign << ")" << endl;
	}
	
	Bimage*			pmask = new Bimage(Integer, TSimple, size(), n);
	pmask->file_name(file_name());
	pmask->sampling(sampling(0));
	for ( i=0; i<n; ++i ) pmask->image[i] = image[i];

	//Using a mask, all the regions are grouped and numbered.
	for ( i=0; i<datasize; ++i ) {
		if ( (*pmask)[i] == 0 ) {
			if ( sign*(*this)[i] > abs_threshold ) {
				count++;
				rsize = region_assign(pmask, i, count, threshold, sign);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG Bimage::regions: i=" << i << " count=" << 
						count << " rsize=" << rsize << endl;
			}
		}
	}
	
	pmask->statistics();
	
	if ( verbose & VERB_PROCESS )
		cout << "Regions:                        " << count << endl << endl;
	
	return pmask;
}


/**
@author	Bernard Heymann and Samuel Payne
@brief 	Finds all the pixels that are part of the same region.
@param 	*pmask			a mask holding region assignments.
@param 	idx				the first voxel of a region.
@param 	region_number	the region number that voxels are assigned.
@param 	threshold 		the level to define background.
@param 	sign			sign controlling direction of thresholding.
@return int				0.

	This method uses a list to keep track of voxels assigned to a region.
	In subsequent iterations the non-assigned voxels next to those in
	the list are assigned if they exceed the threshold, and their indices
	are kept in a new list. The new list is transferred to the old list
	and the process iteratively continued until the new list contains
	no more voxels.

**/
long		Bimage::region_assign(Bimage* pmask, long idx,
					long region_number, double threshold, int sign)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::region_assign: Starting with region " <<
			region_number << " idx=" << idx << endl;
	
	if ( sign < 0 ) sign = -1;
	else if ( sign > 0 ) sign = 1;
	else {
		sign = 1;
		if ( threshold < avg ) sign = -1;
	}
	
	double			abs_threshold(sign*threshold);
	long	   		i, j, nn, kx, ky, kz;
	Vector3<long>	lo, hi;

	// Set up two lists:
	//		First list is for voxels to loop through for checking
	//		Second list is for new voxels to be checked in a subsequent iteration
	long	   		number, new_number, count(1);
	long	   		list_size = 2*(x + y)*z;
	if ( z > 1 ) list_size += 2*x*y;
	int*			list = new int[list_size];
	int*			new_list = new int[list_size];
	
	// Initialize the first item in the list
	pmask->set(idx, region_number);
	number = 1;
	list[0] = idx;
	nn = idx / (x*y*z);
	
	// Find neighbors in shells around already found neighbors
	while ( number > 0 ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG Bimage::region_assign: region=" << region_number << " list size=" << number << endl;
		new_number = 0;
		for ( i=0; i<number; ++i ) {
			lo = kernel_low(list[i]);
			hi = kernel_high(list[i]);
	
			// Loop through the kernel and construct the new list
			for ( kz=lo[2]; kz<=hi[2]; kz++ ) {
				for ( ky=lo[1]; ky<=hi[1]; ky++ ) {
					for ( kx=lo[0]; kx<=hi[0]; kx++ ) {
						j = index(kx, ky, kz, nn);
						if ( ((*pmask)[j] == 0) && (sign*(*this)[j] > abs_threshold) ) {
							pmask->set(j, region_number);
							count++;
							new_list[new_number] = j;
							new_number++;
							if ( new_number >= list_size ) {
								cerr << "Error: The list size exceeded its limit (" 
									<< list_size << ")" << endl;
							}
						}
					}
				}
			}
		}
		number = new_number;
		for ( i=0; i<number; ++i ) list[i] = new_list[i];
	}
	
	delete[] list;
	delete[] new_list;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::region_assign: Ending with region " << region_number << endl;
	
	return count;
}

/**
@brief 	Segments a map through a series of thresholds and reports results.
@param 	threshold_first	the level to pick initial regions.
@param 	threshold_last	the lowest level to include voxels in regions.
@param 	threshold_step 	the incremental change in threshold.
@return int				0.
**/
int			Bimage::region_threshold_series(double threshold_first,
				double threshold_last, double threshold_step)
{
	int			sign(1);
	double		threshold, a(0);
	Bimage* 	pmask;
	
	if ( threshold_step == 0 ) {
		if ( threshold_last == threshold_first )
			threshold_step = 1;
		else
			threshold_step = threshold_last - threshold_first;
	} else {
		if ( threshold_step > 0 && threshold_first > threshold_last )
			threshold_step = -threshold_step;
		if ( threshold_step < 0 && threshold_first < threshold_last )
			threshold_step = -threshold_step;
	}
	
	if ( threshold_step < 0 ) sign = -1;
	threshold_last += threshold_step/2;
	
	cout << "Finding the number of regions at each threshold:" << endl;
	cout << "Thresh\tRegions\tAvgSize" << endl;
	for ( threshold = threshold_first; sign*threshold <= sign*threshold_last; 
			threshold+=threshold_step ) {
		pmask = regions(threshold, 0);
		a = pmask->levelmask_average_region_size();
		cout << threshold << tab << pmask->maximum() << tab << a << endl;
		delete pmask;
	}
	cout << endl;
	
	return 0;
}

/**
@brief 	Creates and expands a region map from a starting to ending threshold value.
@param 	*pmask			region map (if NULL, generate from high threshold).
@param 	threshold_hi 	the level to pick initial regions.
@param 	threshold_lo	the lowest level to include voxels in regions.
@param 	threshold_step 	the incremental change in threshold.
@param 	fill_borders	flag to assign borders between regions (default not).
@return long	 			number of regions.

	A region map is calculated at the threshold furthest from the average.
	The region map is then expanded by lowering the threshold gradually 
	and assigning newly included voxels to neighboring regions (i.e., flooding).
	Voxels with neighbours assigned to two or more different regions are tagged
	as indeterminate to indicate borders between regions. These border voxels
	are counted and can be used to estimate the extent of interfaces.
	The indices of all the final regions are packed into an integer data block
	within a new image.
	The image is assumed to have high values for objects (i.e., density is white).

**/
long 		Bimage::region_flood(Bimage* pmask, double threshold_hi,
				double threshold_lo, double threshold_step, int fill_borders)
{
	if ( !pmask ) {
		error_show("No region mask in Bimage::region_flood", __FILE__, __LINE__);
		return 0;
	}
	
	long	   		i;
	double			threshold, dt;

	if ( threshold_hi < threshold_lo )
		swap(threshold_hi, threshold_lo);
	if ( threshold_step == 0 )
		threshold_step = (threshold_lo - threshold_hi)/100;
	else if ( (threshold_lo - threshold_hi)/threshold_step < 0 )
		threshold_step *= -1;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Segmenting an image by flooding from a threshold or region set:" << endl;
		if ( pmask )
			cout << "Region set:                     " << pmask->file_name() << endl;
		cout << "Starting threshold:             " << threshold_hi << endl;
		cout << "Finishing threshold:            " << threshold_lo << endl;
		cout << "Threshold step:                 " << threshold_step << endl;
	}
    
	if ( !check_if_same_size(pmask) ) {
		error_show("Error in Bimage::region_flood", __FILE__, __LINE__);
		cerr << "The region mask must have the same dimensions as the image!" << endl;
		return 0;
	}

	if ( verbose & VERB_LABEL )
		cout << "Starting regions:                " << pmask->maximum() << endl;
	
	Bimage*			pt = pmask->copy();
	
	// Progressively decrease threshold and find new adjacent voxels belonging to
	// each region - disputes over voxel assignment are settled later
	for ( threshold=threshold_hi; threshold>threshold_lo; threshold+=threshold_step ) {
		dt = threshold_hi - threshold_lo;
		if ( dt > 0 ) dt = 100*(threshold_hi - threshold)/dt;
		else dt = 100;
		if ( verbose )
			cout << "Processing threshold:           " << threshold << " (" << dt << " %)     \r" << flush;
	
#ifdef HAVE_GCD
		dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t it){
			if ( (*this)[it] > threshold )
				pt->set(it, pmask->check_neighbors(it));
			else
				pt->set(it, 0);
		});
#else
#pragma omp parallel for
		for ( long it=0; it<datasize; it++ ) {
			if ( (*this)[it] > threshold )
				pt->set(it, pmask->check_neighbors(it));
			else
				pt->set(it, 0);
		}
#endif
		for ( i=0; i<datasize; ++i ) pmask->set(i, (*pt)[i]);
	}
	
	if ( verbose )
		cout << endl;
	
	long	   		nborder(0), dnb(1), maxtries(20);
	
	for ( i=0; i<datasize; ++i ) if ( (*pmask)[i] < 0 ) nborder++;
	if ( verbose & VERB_LABEL )
		cout << "Border voxels:                  " << nborder << endl;
	
	// Assign the borders between regions to one of the regions
	if ( fill_borders ) {
		if ( verbose & VERB_PROCESS )
			cout << "Assigning borders:" << endl;
			
		while ( maxtries-- > 0 && nborder > 0 && dnb > 0 ) {
			if ( verbose )
				cout << "Processing border voxels:       " << nborder << "     \r" << flush;
			dnb = nborder;
			nborder = 0;

#ifdef HAVE_GCD
			dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t it){
				if ( (*pmask)[it] < 0 ) {
					size_t	j = kernel_max_neigbor(it, 1);
					pmask->set(it, (*pmask)[j]);
					if ( (*pmask)[it] < 0 ) {
						j = kernel_max_neigbor(it, 2);
						pmask->set(it, (*pmask)[j]);
					}
				}
			});
#else
#pragma omp parallel for
			for ( long it=0; it<datasize; it++ ) {
				if ( (*pmask)[it] < 0 ) {
					long	j = kernel_max_neigbor(it, 1);
					pmask->set(it, (*pmask)[j]);
					if ( (*pmask)[it] < 0 ) {
						j = kernel_max_neigbor(it, 2);
						pmask->set(it, (*pmask)[j]);
					}
				}
			}
#endif
			for ( long it=0; it<datasize; it++ ) if ( (*pmask)[it] < 0 ) nborder++;
			
			dnb -= nborder;
		}
		
		if ( verbose )
			cout << "Processing border voxels:       " << nborder << endl;
		if ( verbose && nborder > 0 )
			cerr << "Warning: " << nborder << " border voxels were not assigned!" << endl;
	}
	
	long			nseg = (long)(pmask->maximum() + 1.9);
	int*			vol;
	
	if ( verbose & VERB_FULL ) {
		vol = new int[nseg];
		for ( i=0; i<nseg; ++i ) vol[i] = 0;
		for ( i=0; i<datasize; ++i )
			if ( (*pmask)[i] > -1 ) vol[(long)(*pmask)[i]]++;
		cout << "Region\tVoxels\tVolume\tMass" << endl;
		for ( i=1; i<nseg; ++i )
			cout << i << tab << vol[i] << tab << vol[i]*sampling(0).volume()
				<< tab << vol[i]*sampling(0).volume()*RHO << endl;
		delete[] vol;
	}

	if ( verbose )
		cout << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::region_flood: Finished with segmentation" << endl;
	
	return nseg;
}

/**
@brief 	Check neighbors in region map and assign a value if a neighbor is assigned.
@param 	idx		index in multi-image.
@return long		new value assigned.
**/
long		Bimage::check_neighbors(long idx)
{
	if ( (*this)[idx] != 0 ) return (long)(*this)[idx];
	
	long			i, xx, yy, zz, nn(0), nv(0);

	//calculate the bounds on the kernel
	Vector3<long>	lo = kernel_low(idx);
	Vector3<long>	hi = kernel_high(idx);
	
	// Loop through the kernel and check neighbors
	for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
		for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
			for ( xx=lo[0]; xx<=hi[0]; ++xx ) {
				i = index(0,xx,yy,zz,nn);
				if ( i != idx ) {
					if ( (*this)[i] > 0 ) {
						if ( nv == 0 )
							nv = (long)(*this)[i];
						else if ( nv != (*this)[i] )
							nv = -100;	// Border voxel
					}
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::check_neighbors: Ending with region " << (int)(*this)[idx] << endl;
	
	return nv;
}

/**
@brief 	Generates a segmented image from peaks above a threshold value.
@param 	threshold 		the lowest level to include voxels in regions.
@param 	flag	 		Track towards maxima (0) or minima (1).
@return Bimage* 			integer image containing tracking indices.

	Every voxel with a value above (below) the threshold is tagged by a pointer
	pointing to the highest (lowest) value within a kernel of 3x3x3.
	A peak is defined as pointing to itself, while every other voxel
	points towards a peak.

**/
Bimage*		Bimage::track_gradient(double threshold, int flag)
{
	long 			nn;

	if ( verbose & VERB_LABEL ) {
		cout << "Tracking the gradient:" << endl;
		cout << "Threshold:                      " << threshold << endl;
		cout << "Direction:                      ";
		if ( flag ) cout << "Towards maxima" << endl;
		else cout << "Towards minima" << endl;
	}
    
	Bimage*			ptrack = new Bimage(Integer, TSimple, size(), n);
	for ( nn=0; nn<n; nn++ ) ptrack->image[nn] = image[nn];
	ptrack->fill(-1);

#ifdef HAVE_GCD
	dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t i){
		if ( flag ) {
			if ( (*this)[i] >= threshold ) ptrack->set(i, kernel_max(i, 1));
		} else {
			if ( (*this)[i] <= threshold ) ptrack->set(i, kernel_min(i, 1));
		}
	});
#else
#pragma omp parallel for
	for ( long i=0; i<datasize; ++i ) {
		if ( flag ) {
			if ( (*this)[i] >= threshold ) ptrack->set(i, kernel_max(i, 1));
		} else {
			if ( (*this)[i] <= threshold ) ptrack->set(i, kernel_min(i, 1));
		}
	}
#endif
	
	return ptrack;
}

/**
@brief 	Generates a segmented image from peaks above a threshold value.
@param 	kernel_size		size of kernel to determine peaks.
@param 	threshold 		the lowest level to include voxels in regions.
@param 	flood	 		flag to flood to the threshold.
@param 	wrap	 		flag to wrap the kernel around image boundaries.
@return Bimage* 			integer image containing region index numbers.

	Every voxel with a value above the threshold is tagged by a pointer 
	pointing to the highest value within a kernel. 
	A peak is defined as pointing to itself, while every other voxel
	points towards a peak.
	The voxels in the mask are then iteratively assigned to the peaks
	they point to.
	The indices of all the final regions are packed into an integer data block
	within a new image.

**/
Bimage*		Bimage::region_peaks(long kernel_size, double threshold, int flood, int wrap)
{
	if ( kernel_size < 3 ) kernel_size = 3;
	
	long 			i, nn, dn;
	long			hk(kernel_size/2), npeak(0);

	if ( verbose & VERB_LABEL ) {
		cout << "Segmenting an image by flooding from peaks:" << endl;
		cout << "Kernel size:                    " << 2*hk+1 << endl;
		cout << "Threshold:                      " << threshold << endl;
		if ( wrap ) cout << "With wrapping" << endl;
	}
    
	Bimage*			pmask = new Bimage(Integer, TSimple, size(), n);
	for ( nn=0; nn<n; nn++ ) pmask->image[nn] = image[nn];

	int*			pnt = new int[datasize];
	for ( i=0; i<datasize; ++i ) pnt[i] = -1;
	
#ifdef HAVE_GCD
	dispatch_apply(datasize, dispatch_get_global_queue(0, 0), ^(size_t i){
		if ( (*this)[i] >= threshold ) {
			if ( wrap ) pnt[i] = kernel_max_wrap(i, hk);
			else pnt[i] = kernel_max(i, hk);
		}
	});
#else
#pragma omp parallel for
	for ( i=0; i<datasize; ++i ) {
		if ( (*this)[i] >= threshold ) {
			if ( wrap ) pnt[i] = kernel_max_wrap(i, hk);
			else pnt[i] = kernel_max(i, hk);
		}
	}
#endif

	for ( i=npeak=0; i<datasize; ++i )
		if ( pnt[i] == i ) pmask->set(i, ++npeak);

	if ( flood ) {	
		nn = 0;
		dn = 1;
		while ( dn ) {
			dn = nn;
			for ( i=0; i<datasize; ++i )
				if ( (*this)[i] >= threshold && (*pmask)[i] < 1 )
					pmask->set(i, (*pmask)[pnt[i]]);
			for ( i=nn=0; i<datasize; ++i ) if ( (*pmask)[i] > 0 ) nn++;
			dn = nn - dn;
		}
	}

	delete[] pnt;
	
	pmask->maximum(npeak);
	
	if ( verbose & VERB_LABEL )
		cout << "Number of peaks:                " << npeak << endl << endl;
	
	return pmask;
}

/**
@brief 	Identifies contiguous regions above a given threshold (blobs) and 
	eliminate those less than a minimum size.
@param 	threshold 	threshold to define contiguous regions.
@param 	min_size	minimum number of voxels to keep a region.
@param 	max_size	maximum number of voxels to keep a region.
@param 	setvalue 	value to assign to deleted regions (typically set to threshold).
@param 	sign		sign controlling direction of thresholding.
@return long 			number of voxels set to threshold.

	All contiguous regions in an image are indexed and their sizes 
	calculated. Whether a region is contiguous is defined by a
	threshold. Regions smaller than a minimum number of voxels are
	set to the threshold value. The input image is modified and
	left as a floating point image.

**/
long 		Bimage::blobs(double threshold, double min_size, double max_size,
				double setvalue, int sign)
{
	if ( min_size > max_size ) swap(min_size, max_size);
	
	if ( sign < 0 ) sign = -1;
	else if ( sign > 0 ) sign = 1;
	else {
		sign = 1;
		if ( threshold < avg ) sign = -1;
	}
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		if ( sign > 0 )
			cout << "Finding blobs above " << threshold << endl;
		else
			cout << "Finding blobs below " << threshold << endl;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Size range to keep:             " << min_size << " - " << max_size << endl;
		cout << "Set value:                      " << setvalue << endl << endl;
    }
    
	Bimage*			pmask = regions(threshold, sign);
	long			count = (long) pmask->maximum();
	
	long   		i, f, nvox(0), nblobs(0);

	int*			region_size = new int[count+1];

	for ( i=0; i<=count; ++i ) region_size[i] = 0;
	
	for ( i=0; i<datasize; ++i ) region_size[(long)(*pmask)[i]]++;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Region\tSize";
		for ( i=0; i<=count; ++i ) {
			if ( i%5 == 0 ) cout << endl << i;
			cout << tab << region_size[i];
		}
		cout << endl << endl;
	}
	
	if ( ( min_size > 0 ) || ( max_size < datasize ) ) {
		for ( i=1; i<=count; ++i )
			if ( region_size[i] < min_size || region_size[i] > max_size ) nblobs++;
	
		for ( i=0; i<datasize; ++i ) {
			f = (long)(*pmask)[i];
			if ( f > 0 ) {
				if ( ( region_size[f] < min_size ) || ( region_size[f] > max_size ) ) {
					set(i, setvalue);
					nvox++;
				}
			}
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of blobs deleted:        " << nblobs << endl;
		cout << "Number of voxels set:           " << nvox << endl << endl;
	}
	
	delete[] region_size;
	
//	write_img("pf.pif", pmask);
	
	delete pmask;
	
	statistics();
	
	return nvox;
}

/**
@brief 	Filters the extremes out of a micrograph image.
@return int 				0.

	Segmentation is used to identify large contiguous regions of high
	and low values, and these regions are set to the average.

**/
int			Bimage::filter_extremes()
{
	statistics();
	
	double		threshold = avg + 3*std;
	if ( threshold >= max ) threshold = (max - avg)*0.75 + avg;
	
	blobs(threshold, 100, 1e30, avg, 0);

	threshold = avg - 3*std;
	if ( threshold <= min ) threshold = (min - avg)*0.75 + avg;
	
	blobs(threshold, 1000, 1e30, avg, 0);
	
	return 0;
}

/**
@brief 	Filters the extremes out of an image.
@param 	mod_flag		modification flag (0=set to avg, 1=set to min_max).
@return int 			0.

	A histogram of an image is calculated.  The first minimum in the first
	quarter of the histogram and the last minimum in the last quarter
	are taken to define the small and large outliers.  Pixel outside
	these minima are then either set to average or to minimum or maximum
	(depending on the mod_flag argument) to remove the outliers.

**/
int			Bimage::filter_extremes(int mod_flag)
{
	double		tmin, tmax;
	
	histogram_minmax(tmin, tmax);
	
    if ( mod_flag )
		truncate_to_min_max(tmin, tmax);
	else
		truncate_to_avg(tmin, tmax);
	
	return 0;
}

/**
@brief 	Filters the extremes out of an image by replacing with adjacent averages.
@param 	tmin		minimum.
@param 	tmax		maximum.
@param 	kernel		kernel edge size.
@return int 			0.

	Pixels smaller than the minimum or larger than the maximum are set to
	the average within a defined kernel.

**/
int			Bimage::filter_extremes(double tmin, double tmax, int kernel)
{
	if ( !data_pointer() ) return -1;
	
    if ( tmin < min ) tmin = min;
    if ( tmax > max ) tmax = max;
	if ( tmin > tmax ) swap(tmin, tmax);
	if ( kernel < 3 ) kernel = 3;
   
    long   			i, kh(kernel/2);
    double  	    v1;
    
	if ( verbose & VERB_PROCESS ) {
	    cout << "Truncating to:                  " << tmin << " " << tmax << endl;
	    cout << "Kernel size:                    " << kernel << endl << endl;
	}
	
    for ( i=0; i<datasize; i++ ) {
 		v1 = (*this)[i];
		if ( v1 < tmin || v1 > tmax )
			set(i, kernel_average(i, kh, tmin, tmax));
	}

	statistics();

	return 0;
}


/**
@brief Replaces maxima above a threshold using local averages.
@param	threshold 	threshold.
@return long 		number of voxels replaced.

	The image is first segemented into regions, each defined as a 
	contiguous cluster of voxels with values above the threshold. 
	Each region is encoded in an integer image, with all the values
	in this image set to the indices of the regions, or zero elsewhere.
	The regions are iteratively filled in, where in every iteration, 
	the border values of the regions are replaced by the average 
	of the neighbouring voxels outside the region. After each iteration, 
	every region is shrunk by excluding the newly replaced voxels.
**/
long		Bimage::replace_maxima(double threshold)
{
	Bimage*			pmask = regions(threshold, 0);

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Replacing maxima above:         " << threshold << endl;
		cout << "Number of maxima:               " << pmask->maximum() << endl << endl;
	}
	
	if ( pmask->maximum() < 1 ) {  // No regions found = no modification
		delete pmask;
		return 0;
	}
	
//	long   			i, j, k, m, nn, xx, yy, zz, iter, nvr(0);
//	long			kx, ky, kz, nval, tx, ty, tz;
//	long   			nr((long) (pmask->maximum() + 1.9)), rs;
	long   			i, j, nn, xx, yy, zz, nvr(0);
	long   			nr((long) (pmask->maximum() + 1.9));
	double			val;
	Vector3<long>	lo, hi, region_size;
	
	vector<int>				num(nr,0);
	vector<int>				nimg(nr,0);
	vector<Vector3<float>>	region(nr);
	vector<Vector3<int>>	rmin(nr);
	vector<Vector3<int>>	rmax(nr);
	vector<float>			a(nr,0);
	vector<int>				na(nr,0);
	
	for ( j=1; j<nr; j++ ) rmin[j] = size();
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; ++zz ) {
			for ( yy=0; yy<y; ++yy ) {
				for ( xx=0; xx<x; ++xx, ++i ) {
					j = (long) (*pmask)[i];
					if ( j ) {
						num[j]++;
						nimg[j] = nn;
						region[j][0] += xx;
						region[j][1] += yy;
						region[j][2] += zz;
						if ( rmin[j][0] > xx ) rmin[j][0] = xx;
						if ( rmin[j][1] > yy ) rmin[j][1] = yy;
						if ( rmin[j][2] > zz ) rmin[j][2] = zz;
						if ( rmax[j][0] < xx ) rmax[j][0] = xx;
						if ( rmax[j][1] < yy ) rmax[j][1] = yy;
						if ( rmax[j][2] < zz ) rmax[j][2] = zz;
						val = kernel_average(i, 3, -1e30, threshold);
						if ( val > min ) {
							a[j] += val;
							na[j]++;
						}
					}
				}
			}
		}
	}

	if ( verbose & VERB_STATS )
		cout << "Regions:\n#\tImage\tx\ty\tz\tSize\tAvg" << endl;
	for ( j=1; j<nr; j++ ) {
		nn = nimg[j];
		if ( num[j] ) region[j] /= num[j];
		if ( na[j] ) a[j] /= na[j];
		if ( verbose & VERB_STATS )
			cout << j << tab << nn << tab << region[j] << tab << num[j] << tab << a[j] << endl;
	}
	
	for ( i=0; i<datasize; ++i ) {
		j = (long) (*pmask)[i];
		if ( j ) set(i, a[j]);
	}
	
/*
	if ( verbose & VERB_STATS )
		cout << "Regions:\n#\tImage\tx\ty\tz\tSize" << endl;
	for ( j=1; j<nr; j++ ) {
		nn = nimg[j];
		if ( num[j] ) region[j] /= num[j];
		if ( navg[j] ) avg[j] /= navg[j];
		if ( verbose & VERB_STATS )
			cout << j << tab << nn << tab << region[j] << tab << num[j] << endl;
		rmin[j] = rmin[j].max(0);
		rmax[j] = rmax[j].min(size()-1);
		region_size = rmax[j] - rmin[j] + 1;
		rs = (long) region_size.volume();
		vector<float>	temp(rs,0);
		nval = iter = 0;
		val = 0;
		while ( nval < num[j] ) {
			for ( k=0; k<rs; k++ )
				temp[k] = min;
			for ( zz=rmin[j][2], tz=0; zz<=rmax[j][2]; ++zz, tz++ ) {
				for ( yy=rmin[j][1], ty=0; yy<=rmax[j][1]; ++yy, ty++ ) {
					for ( xx=rmin[j][0], tx=0; xx<=rmax[j][0]; ++xx, tx++ ) {
						i = index(0,xx,yy,zz,nn);
						set(i, avg[j]);
						lo = kernel_low(i);
						hi = kernel_high(i);
						if ( (*pmask)[i] == j ) {
							val = 0;
							m = 0;
							for ( kz=lo[2]; kz<=hi[2]; kz++ ) {
								for ( ky=lo[1]; ky<=hi[1]; ky++ ) {
									for ( kx=lo[0]; kx<=hi[0]; kx++ ) {
										k = index(0,kx,ky,kz,nn);
										if ( (*pmask)[k] < 1 ) {
											val += (*this)[k];
											m++;
										}
									}
								}
							}
							if ( m ) {
								k = (tz*region_size[1] + ty)*region_size[0] + tx;
								temp[k] = val/m;
							}
						}
					}
				}
			}
			for ( zz=rmin[j][2], k=tz=0; zz<=rmax[j][2]; ++zz, tz++ ) {
				for ( yy=rmin[j][1], ty=0; yy<=rmax[j][1]; ++yy, ty++ ) {
					for ( xx=rmin[j][0], tx=0; xx<=rmax[j][0]; ++xx, tx++, k++ ) {
						if ( temp[k] > min ) {
							i = index(0,xx,yy,zz,nn);
							set(i, temp[k]);
							if ( (*pmask)[i] == j ) pmask->set(i, 0);
							nval++;
						}
					}
				}
			}
			iter++;
			if ( verbose & VERB_FULL )
				cout << iter << ":\t" << nval << endl;
		}
		nvr += nval;
	}*/
	if ( verbose & ( VERB_STATS | VERB_FULL ) )
		cout << endl;

	delete pmask;
	
	return nvr;
}

/**
@brief 	Finds the density threshold associated with a particular molecular weight.
@param 	img_num		sub-image number (first = 0,).
@param 	mol_weight 	molecular weight.
@param 	rho			protein density in Da/A3.
@return float 		threshold.

	A threshold value is determined which contour a density such that
	it corresponds to 100% of the given molecular weight.
	An image is assumed to have density represented by higher (white) values.

**/
double 		Bimage::mass_threshold(long img_num, double mol_weight, double rho)
{
	if ( rho <= 0 ) rho = RHO;
	
	long   			i, j, nn, nvox(0);
	long			dnvox(-1);
	long   			imgsize(x*y*z);
	
	double			mass, dmass, vs(voxel_size());
	double			volume(mol_weight/rho);
	double			threshold(0);
	double			thresh_step;
	double			vol_diff;
	double			min_vol(vs);
	double			min_step(1e-6*std);
	
	if ( verbose & VERB_LABEL ) {
		cout << "Finding the threshold for image " << img_num+1 << ":" << endl;
		cout << "Protein density:                " << rho << " Da/A3" << endl;
		cout << "Molecular weight:               " << mol_weight << " Da" << endl;
		cout << "Molecular volume:               " << volume << " A3" << endl;
		cout << "Voxel size:                     " << vs << " A3" << endl << endl;
    }
    
    if ( volume > real_size().volume() ) {
		cerr << "Warning: The image volume is smaller than the molecular volume!" << endl << endl;
		return 0;
    }
	
	if ( verbose )
		cout << "Image\tThresh\tVolRes\tIntMass\tMW/mass" << endl;
	for ( nn=0; nn<n; nn++ ) {
		if ( verbose & VERB_FULL )
			cout << " Threshold\tThreshStep\tVolume(A3)\tDifference\tMass" << endl;
		threshold = avg;
		thresh_step = std/2;
		nvox = 0;
		mass = 0;
		vol_diff = 1e30;
		while ( fabs(vol_diff) >= min_vol && fabs(thresh_step) > min_step ) {
			dnvox = nvox;
			dmass = mass;
			if ( vol_diff > 0 ) {
				if ( thresh_step < 0 ) thresh_step *= -0.5;
			} else {
				if ( thresh_step > 0 ) thresh_step *= -0.5;
			}
			threshold += thresh_step;
			for ( i=nn*imgsize, j=i+imgsize, nvox = 0, mass = 0; i<j; ++i ) {
				if ( (*this)[i] >= threshold ) {
					nvox++;
					mass += (*this)[i];
				}
			}
			vol_diff = nvox*vs - volume;
			if ( verbose & VERB_FULL )
				cout << threshold << tab << thresh_step << tab 
					<< nvox*vs << tab << vol_diff << tab << mass << endl;
			dnvox = nvox - dnvox;
			dmass = mass - dmass;
		}
		if ( verbose )
			cout << nn+1 << tab << threshold << tab <<
				fabs(vol_diff)/volume << tab << mass << tab << mol_weight/mass << endl;
	}
	if ( verbose )
		cout << endl;

/*	
	if ( verbose & VERB_LABEL ) {
		cout << "Threshold:                      " << threshold << endl;
		cout << "Relative volume residual:       " << fabs(vol_diff)/volume << endl;
		cout << "Integrated mass:                " << mass << endl;
		cout << "Molecular weight/mass:          " << mol_weight/mass << endl << endl;
	}
*/	
	return threshold;
}

/**
@brief 	Calculates the mass from the density threshold.
@param 	img_num		sub-image number (first = 0).
@param 	threshold 	density threshold.
@param 	rho			protein density in Da/A3.
@return double 		threshold.

	An image is assumed to have density represented by higher (white) values.

**/
double 		Bimage::mass_at_threshold(long img_num, double threshold, double rho)
{
	if ( rho <= 0 ) rho = RHO;
	
	long			i, j, nvox;
	long   		imgsize(x*y*z);
	
	for ( i=img_num*imgsize, j=i+imgsize, nvox=0; i<j; ++i )
		if ( (*this)[i] >= threshold ) nvox++;

	double			vs(voxel_size());
	double			volume(nvox*vs);
	double			mass(volume*rho);

	if ( verbose & VERB_LABEL ) {
		cout << "Protein density:                " << rho << " Da/A3" << endl;
		cout << "Threshold:                      " << threshold << endl;
		cout << "Volume:                         " << volume << " A3" << endl;
		cout << "Molecular weight:               " << mass << " Da" << endl << endl;
	}
	
	return mass;
}

/**
@brief 	Calculates the internal volume of a shell.
@param 	threshold 		threshold to define density.
@return Bimage*			internal volume mask.

**/
Bimage*		Bimage::internal_volume(double threshold)
{
	long   			i, j, jn, jz, jy, nn;
	long   			xx, yy, zz;
	long			slicesize(x*y), imgsize(x*y*z);
	long			mode, kmax;
	double			d0, dk, dd;
	Vector3<long>	lo, hi, nk(3,3,3), k;

	int				sign(1);
	if ( threshold < avg ) sign = -1;
	
	Bimage*			pv = copy();
	pv->fill(0);

	for ( i=nn=0; nn<n; nn++ ) {
		jn = nn*imgsize;
		for ( zz=0; zz<z; ++zz ) {
			lo[2] = (zz>nk[2])? zz-nk[2]: 0;
			hi[2] = (zz<z-nk[2])? zz+nk[2]: zz;
			cout << "z = " << zz << endl;
			for ( yy=0; yy<y; ++yy ) {
				lo[1] = (yy>nk[1])? yy-nk[1]: 0;
				hi[1] = (yy<y-nk[1])? yy+nk[1]: yy;
				for ( xx=0; xx<x; ++xx, ++i ) {
					lo[0] = (xx>nk[0])? xx-nk[0]: 0;
					hi[0] = (xx<x-nk[0])? xx+nk[0]: xx;
					mode = kmax = 0;
					if ( sign*(*this)[i] >= threshold ) {
						mode = 1;
						d0 = (Vector3<double>(xx, yy, zz) -
							image[nn].origin()).length();
					}
					jz = jn + lo[2]*slicesize;
					dd = 0;
					for ( k[2]=lo[2]; k[2]<=hi[2]; ++k[2], jz += slicesize ) {
						jy = jz + lo[1]*x;
						for ( k[1]=lo[1]; k[1]<=hi[1]; ++k[1], jy += x ) {
							j = jy + lo[0];
							for ( k[0]=lo[0]; k[0]<=hi[0]; ++k[0], ++j ) {
								if ( mode ) {
									dk = (image[nn].origin() - k).length();
									if ( sign*(*this)[j] < threshold )
										dd += d0 - dk;
									if ( dk < d0 && sign*(*this)[j] < threshold )
										kmax++;
//										pv->set(j, 1);
//								} else {
//									if ( kmax < (*pv)[j] ) kmax = (*pv)[j];
								}
							}
						}
					}
//					if ( kmax > 6 ) pv->set(i, 1);
					if ( dd > 12 ) pv->set(i, 1);
				}
			}
		}
	}
	
	return pv;
}

/**
@brief 	Calculates the internal volume of a shell.
@param 	threshold 		threshold to define density.
@param 	mask_out_freq	mask output frequency.
@return Bimage*			internal volume mask.

**/
Bimage*		Bimage::internal_volume(double threshold, int mask_out_freq)
{
	long   			i, j, jn, jz, jy, nn;
	long   			xx, yy, zz, kx, ky, kz;
	long   			iter, dvol, vol(0);
	long			slicesize(x*y), imgsize(x*y*z);
	double			vs(voxel_size());
	Vector3<long>	lo, hi;
	Vector3<long>	lo2, hi2;

	Bstring			fn;
	
	int				nb, nbcut(6), sign(1);
	if ( z < 2 ) nbcut = 2;
	if ( threshold < avg ) sign = -1;
	
	Bimage*			pv = copy();
	pv->fill(0);

	for ( nn=0; nn<n; nn++ ) {
		jn = nn*imgsize;
		lo2 = image[nn].origin();
		i = index(lo2,nn);
		lo = image[nn].origin() - 1;
		lo = lo.max(0);
		hi = image[nn].origin() + 1;
		hi = hi.min(size()-1);
		pv->set(i, 2);
		dvol = vol = 1;
		if ( verbose ) {
			cout << "Image " << nn+1 << ":" << endl;
			cout << "Origin: " << lo2 << endl;
			cout << "Iter\tVol\tdVol\tVol(A3)\tdVol(A3)" << endl;
		}
		iter = 0;
		while ( dvol > 0 ) {
			dvol = 0;
			iter++;
			for ( zz=lo[2]; zz<=hi[2]; ++zz ) {
				lo2[2] = (zz)? zz-1: 0;
				hi2[2] = (zz<z-1)? zz+1: zz;
				for ( yy=lo[1]; yy<=hi[1]; ++yy ) {
					lo2[1] = (yy)? yy-1: 0;
					hi2[1] = (yy<y-1)? yy+1: yy;
					for ( xx=lo[0]; xx<=hi[0]; ++xx ) {
						lo2[0] = (xx)? xx-1: 0;
						hi2[0] = (xx<x-1)? xx+1: xx;
						i = index(0,xx,yy,zz,nn);
						nb = 0;
						if ( (*pv)[i] < 1 && sign*(*this)[i] <= threshold ) {
							jz = jn + lo2[2]*slicesize;
							for ( kz=lo2[2]; kz<=hi2[2]; kz++, jz += slicesize ) {
								jy = jz + lo2[1]*x;
								for ( ky=lo2[1]; ky<=hi2[1]; ky++, jy += x ) {
									j = jy + lo2[0];
									for ( kx=lo2[0]; kx<=hi2[0]; kx++, j++ ) {
										if ( i != j && (*pv)[j] > 1 ) nb++;
									}
								}
							}
							if ( ( nb >= vol ) || ( nb > nbcut ) ) {
								pv->set(i, 1);
								dvol++;
							} else pv->set(i, 0);
						}
					}
				}
			}
			lo2 = lo;
			hi2 = hi;
			for ( zz=lo2[2]; zz<=hi2[2]; ++zz ) {
				for ( yy=lo2[1]; yy<=hi2[1]; ++yy ) {
					for ( xx=lo2[0]; xx<=hi2[0]; ++xx ) {
						i = index(0,xx,yy,zz,nn);
						if ( (*pv)[i] == 1 ) {
							pv->set(i, 2);
							if ( lo[0] == xx && xx > 0 ) lo[0]--;
							if ( lo[1] == yy && yy > 0 ) lo[1]--;
							if ( lo[2] == zz && zz > 0 ) lo[2]--;
							if ( hi[0] == xx && xx < x - 1 ) hi[0]++;
							if ( hi[1] == yy && yy < y - 1 ) hi[1]++;
							if ( hi[2] == zz && zz < z - 1 ) hi[2]++;
						}
					}
				}
			}
			vol += dvol;
			if ( verbose )
				cout << iter << tab << vol << tab << dvol << tab 
					<< vol*vs << tab << dvol*vs << endl;
			if ( ( mask_out_freq > 0 && iter%mask_out_freq == 0 ) ||
					( dvol <= 0 ) ) {
				if ( n == 1 ) fn = Bstring(iter, "v%04ld.mrc");
				else fn = Bstring(nn, "v%02ld") + Bstring(iter, "_%04ld.mrc");
//				write(fn);
				write_img(fn, pv, 0);
			}
		}
	}
	
	return pv;
}

/**
@brief 	Segments an image based on K-means.
@param 	nregion 		number of regions.
@param 	max_iter 	maximum number of iterations.
@param 	ratio	 	balance between density and distance.
@return Bimage* 		segmentation mask.

	The metric to choose region membership is based on the minimum of:
		d = |<c>-c|/s + r|<v>-v|
	for:
		c: 		coordinates
		<c>:	region average coordinates.
		s:		image size.
		v:		density.
		<v>:	region average density.
		r:		ratio.

**/
Bimage*		Bimage::kmeans_segment(long nregion, long max_iter, double ratio)
{
	random_seed();
	
	long			i, j, j4, h, k, m;
	double			v, dv, mdv, da(1e30);
	double			r(ratio/std);
	Vector3<long> 	coor;
	Vector3<double> vec, start, end(size()-1);

	long*			nm = new long[nregion];
	double*			a = new double[4*nregion];	// 3 values for coordinates and one for density
	double*			an = new double[4*nregion];

	Bimage*			pseg = copy_header();
	pseg->data_type(Integer);
	pseg->data_alloc_and_clear();
	
	if ( verbose ) {
		cout << "Segmentation using K-means:" << endl;
		cout << "Number of regions:             " << nregion << endl;
		cout << "Density:distance ratio:        " << ratio << endl << endl;
	}
	
	// Initial averages
	for ( j=j4=0; j<nregion; j++, j4+=4 ) {
		vec = vector3_random(start, end);
		for ( h=0; h<3; h++ ) a[j4+h] = vec[h];
		a[j4+3] = get(0, vec);
		if ( verbose & VERB_FULL )
			cout << j+1 << tab << vec << tab << a[j4+3] << endl;
	}
	
	// Iterations
	if ( verbose )
		cout << "#\tChange" << endl;
	for ( i=0; i<max_iter && da > 0.001; ++i ) {
		for ( j=0; j<nregion; j++ ) nm[j] = 0;
		for ( j=0; j<4*nregion; j++ ) an[j] = 0;
	
		for ( k=0; k<datasize; k++ ) {
			v = (*this)[k];
			coor = coordinates(k);
			for ( j=j4=m=0, mdv=1e30; j<nregion; j++, j4+=4 ) {
				if ( i ) dv = r * fabs(a[j4+3] - v);
				else dv = 0;	// First iteration only distance
				for ( h=0; h<3; h++ ) dv += fabs(a[j4+h] - coor[h])/size()[h];
				if ( mdv > dv ) {
					mdv = dv;
					m = j;
				}
			}
			an[4*m+3] += v;
			for ( h=0; h<3; h++ ) an[4*m+h] += coor[h];
			nm[m]++;
			pseg->set(k, m);
		}

		for ( j=j4=0, da=0; j<nregion; j++, j4+=4 ) {
			if ( nm[j] < 1 ) {
				vec = vector3_random(start, end);
				for ( h=0; h<3; h++ ) an[j4+h] = vec[h];
				an[j4+3] = get(0, vec);
			}
			if ( verbose & VERB_FULL )
				cout << j+1 << tab << nm[j];
			for ( h=0; h<4; h++ ) {
				if ( nm[j] ) an[j4+h] /= nm[j];
				da += fabs(a[j4+h] - an[j4+h]);
				a[j4+h] = an[j4+h];
				if ( verbose & VERB_FULL )
					cout << tab << a[j4+h];
			}
			if ( verbose & VERB_FULL )
				cout << endl;
		}
		da /= 4*nregion;
		if ( verbose )
			cout << i+1 << tab << da << endl;
	}
	
	delete[] nm;
	delete[] a;
	delete[] an;
	
	return pseg;
}

/**
@brief 	Initializing voxels and edges for graph-based segmentation.
@param 	connect_type	connection type: 0=direct neighbors, 1=all neighbors.
@return GSgraph			graph with segment designations.
@remarks	R. Nock, F. Nielsen: Statistical Region Merging.
	IEEE Trans. Pattern Anal. Mach. Intell. 26(11): 1452-1458 (2004)

	Edges are set up with 6 or 26 neighbors.

**/
GSgraph		Bimage::graph_setup(int connect_type)
{
	long			i, xx, yy, zz, x1, y1, z1;
	double			a;
	GSgraph			g;

	// Sets up voxels
	for ( i=0; i<image_size(); ++i ) {
		a = (*this)[i];
		if ( compound_type() == TRGB )
			a = rgb(i).average();
		g.add_voxel(i, 0, 1, a);
	}

	// Sets up edges
	for ( i=zz=0; zz<z; ++zz ) {
		for ( yy=0; yy<y; ++yy ) {
			for ( xx=0; xx<x; ++xx, ++i ) {
				if ( xx < x - 1 ) {
					x1 = xx + 1;
					if ( connect_type ) {
						for ( z1=(zz>0)?zz-1:zz; z1<zz+1 && z1<z; z1++ ) {
							for ( y1=(yy>0)?yy-1:yy; y1<yy+1 && y1<y; y1++ ) {
								g.add_edge(i, (z1*y + y1)*x + x1, 0);
							}
						}
					} else {
						g.add_edge(i, (zz*y + yy)*x + x1, 0);
					}
				}
				if ( yy < y - 1 ) {
					y1 = yy + 1;
					if ( connect_type ) {
						for ( z1=(zz>0)?zz-1:zz; z1<zz+1 && z1<z; z1++ ) {
							g.add_edge(i, (z1*y + y1)*x + xx, 0);
						}
					} else {
						g.add_edge(i, (zz*y + y1)*x + xx, 0);

					}
				}
				if ( zz < z - 1 ) {
					z1 = zz + 1;
					g.add_edge(i, (z1*y + yy)*x + xx, 0);
				}
			}
		}
	}
	
	if ( compound_type() == TSimple ) {
		for ( i=0; i<g.edge_count(); ++i )
			g.edge(i).weight(fabs((*this)[g.edge(i)[0]] - (*this)[g.edge(i)[1]]));
	} else if ( compound_type() == TRGB ) {
		for ( i=0; i<g.edge_count(); ++i )
			g.edge(i).weight(fabs((rgb(g.edge(i)[0]) - rgb(g.edge(i)[1])).rms()));
	}
	
	g.edge_sort();
	
	return g;
}

/**
@brief 	Graph-based segmentation of an image.
@param 	type			segmentation type: 1=threshold, 2=statistical region merging.
@param 	connect_type	connection type: 0=direct neighbors, 1=all neighbors.
@param 	complexity		a value determining the number of segments.
@param 	min_size		minimum segment size.
@return GSgraph			graph with segment designations.

	An array of edges between neighboring voxels is set up and sorted in
	non-decreasing order. The voxels are then aggregated into regions based
	on two selectable criteria:
	simple		edge difference with adjustable threshold.
	srm			statistical region merging.
	Regions below a given cutoff size are merged.
	Only the first sub-image is segmented.

@remarks	R. Nock, F. Nielsen: Statistical Region Merging.
	IEEE Trans. Pattern Anal. Mach. Intell. 26(11): 1452-1458 (2004)

**/
GSgraph		Bimage::graph_segment(int type, int connect_type,
				double complexity, long min_size)
{
	long			dim = (z > 1)? 3: 2;
	double			threshold(std/complexity);
	
	if ( verbose ) {
    	cout << "Segmenting image:" << endl;
		cout << "Algorithm:                      ";
		if ( type == 1 ) {
			cout << "Simple" << endl;
			cout << "Threshold:                      " << threshold << endl;
		} else {
			cout << "Statistical region merging" << endl;
			cout << "Complexity:                     " << complexity << endl;
		}
		cout << "Connectivity:                   ";
		if ( connect_type < 1 ) cout << 2*dim << endl;
		else if ( dim < 3 ) cout << 8 << endl;
		else cout << 26 << endl;
		cout << "Minimum size:                   " << min_size << endl;
	}
	
	GSgraph			g = graph_setup(connect_type);
	
	if ( verbose )
		cout << "Number of edges:                " << g.edge_count() << endl;
	
	long			i, j, nrc(image_size());
	
	// Find segment memberships
	if ( type == 1 ) {
		nrc = g.region_merging(threshold);
	} else {
		nrc = g.statistical_region_merging(threshold);
	}
	
	// Merge small regions
	if ( min_size > 0 ) {
		if ( verbose )
			cout << "Regions before merging:         " << nrc << endl;
		nrc = g.region_merge_small(nrc, min_size);
	}
	
	g.region_number();

	for ( i=0; i<image_size(); ++i ) g.voxel(i).clear();

	for ( i=0; i<image_size(); ++i ) {
		j = g.region_find(i);
		g.voxel(j).add_to_average((*this)[i]);
		g.voxel(j).increment_voxels();
	}
	
	for ( i=0; i<image_size(); ++i )
		if ( g.voxel(j).voxels() ) g.voxel(i).divide_average_by_voxels();

	if ( verbose )
		cout << "Number of regions:              " << nrc << endl << endl;

	g.show_regions();
	
	return g;
}

/**
@brief 	Converting a graph-based segmentation to an image.
@param 	g				graph segmentation.
@return Bimage*			segmented image.

**/
Bimage*		Bimage::graph_segments_to_image(GSgraph& g)
{
	long		i, j;
	
	if ( verbose & VERB_PROCESS )
		cout << "Converting a graph segmentation to an image" << endl << endl;
	
	Bimage*		pseg = new Bimage(Float, TSimple, size(), 1);
	
	for ( i=0; i<image_size(); ++i ) {
		j = g.region_find(i);
		pseg->set(i, g.voxel(j).average());
	}
	
	pseg->sampling(sampling(0));
	pseg->origin(image->origin());
	pseg->statistics();

	return pseg;
}

/**
@brief 	Converting a graph-based segmentation to a multi-level mask.
@param 	g				graph segmentation.
@return Bimage*			segmented image.

**/
Bimage*		Bimage::graph_segments_to_mask(GSgraph& g)
{
	if ( verbose & VERB_PROCESS )
		cout << "Converting a graph segmentation to a multi-level mask" << endl << endl;
	
	Bimage*		pseg = new Bimage(Integer, TSimple, size(), 1);
	(*pseg)["type"] = "mask";

	for ( long i=0; i<image_size(); ++i )
		pseg->set(i, g.voxel(i).joins());

	pseg->sampling(sampling(0));
	pseg->origin(image->origin());
	pseg->statistics();

	return pseg;
}

/*
	For a pixel i:
		calculate kernel extents, size=k
		calculate superpixel centers, O(k): (4+c)*add
		calculate superpixel weights, O(k): (c-1)*add + (c+1)*mul
*/
int			img_update_segment(long i, Bimage* p, Bimage* pmask, vector<long> vstep,
				double colorweight, vector<Bsuperpixel>& seg)
{
	long			j, xx, yy, zz;
	double			d, icw(1/colorweight);
	double			m(p->standard_deviation()*p->standard_deviation()/colorweight);
	vector<double>	dvec = seg[i].coordinates();
	vector<long>	vec(dvec.begin(), dvec.end());

	vector<long>	kmin = vec - vstep;
	kmin = vecmax(kmin, 0);
//	kmin = vecmin(kmin, p->size());
	if ( kmin[0] > p->sizeX() ) kmin[0] = p->sizeX();
	if ( kmin[1] > p->sizeY() ) kmin[1] = p->sizeY();
	if ( kmin[2] > p->sizeZ() ) kmin[2] = p->sizeZ();
	
	vector<long>	kmax = vec + vstep;
//	kmax = vecmin(kmax, p->size());
	if ( kmax[0] > p->sizeX() ) kmax[0] = p->sizeX();
	if ( kmax[1] > p->sizeY() ) kmax[1] = p->sizeY();
	if ( kmax[2] > p->sizeZ() ) kmax[2] = p->sizeZ();

	seg[i].clear();
	seg[i].weight(m);

	for ( zz=kmin[2]; zz<kmax[2]; ++zz ) {
		for ( yy=kmin[1]; yy<kmax[1]; ++yy ) {
			for ( xx=kmin[0]; xx<kmax[0]; ++xx ) {
				j = p->index(xx, yy, zz);
				if ( i == (*pmask)[j]) {
					seg[i].add(xx, yy, zz, p->values(j));
				}
			}
		}
	}
	
	seg[i].average();

	for ( zz=kmin[2]; zz<kmax[2]; ++zz ) {
		for ( yy=kmin[1]; yy<kmax[1]; ++yy ) {
			for ( xx=kmin[0]; xx<kmax[0]; ++xx ) {
				j = p->index(xx, yy, zz);
				if ( i == (*pmask)[j]) {
					d = seg[i].difference2(p->values(j));
					d *= icw;
					if ( seg[i].weight() < d ) seg[i].weight(d);
				}
			}
		}
	}
	
	return 0;
}

int			Bimage::superpixels_update(Bimage* pmask, vector<long> vstep,
				double colorweight, vector<Bsuperpixel>& seg)
{
	long				nseg(seg.size());
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG superpixels_update" << endl;

#ifdef HAVE_GCD
	dispatch_apply(nseg, dispatch_get_global_queue(0, 0), ^(size_t i){
		img_update_segment(i, this, pmask, vstep, colorweight, seg);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nseg; ++i )
		img_update_segment(i, this, pmask, vstep, colorweight, seg);
#endif

	return 0;
}

/*
	For a pixel i:
		For each neighboring superpixel center m:
			fom = wd * distance^2 + color_difference^2 / superpixel_weight
			set the mask to the minimum fom
	O(NNEIGHBOR)
			(6+c)*add + (7+c)*mul + 1*div
*/
int			img_assign_pixel(long i, Bimage* p, Bimage* pmask, double wd, 
				vector<Bsuperpixel>& seg)
{
	long			j((*pmask)[i]), k, m;
	double			dd, dc, d, fom(1e30);
	Vector3<long>	coor3 = p->coordinates(i*p->channels());
//	vector<double>	coor = {double(coor3[0]), double(coor3[1]), double(coor3[2])};
	vector<double>	scoor;

	vector<double>	val = p->values(i);
	
//	cout << i << tab << j << tab << coor3 << endl;

	for ( k=-1, m=j; k<NNEIGHBOR;  ) {
		scoor = seg[m].coordinates() - coor3;
		dd = length2(scoor);
		dc = seg[m].difference2(val);
		d = dd*wd + dc/seg[m].weight();
		if ( fom >= d ) {
			fom = d;
			pmask->set(i, m);
		}
		if ( ++k == NNEIGHBOR || ( m = seg[j].neighbor(k) ) < 0 ) break;
	}
	
	return 0;
}

int			img_assign_segments(Bimage* p, Bimage* pmask, vector<long> vstep,
				vector<Bsuperpixel>& seg)
{
	long			i, nc(0);	
	double			dr(volume(vstep)), wd(1.0L/dr);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG: img_assign_segments" << endl;

	Bimage*			pmc = pmask->copy();
	
#ifdef HAVE_GCD
	dispatch_apply(p->image_size(), dispatch_get_global_queue(0, 0), ^(size_t i){
		img_assign_pixel(i, p, pmask, wd, seg);
	});
#else
#pragma omp parallel for
	for ( i=0; i<p->image_size(); ++i )
		img_assign_pixel(i, p, pmask, wd, seg);
#endif
	
	for ( i=0; i<pmask->image_size(); ++i ) if ( (*pmask)[i] - (*pmc)[i] ) nc++;
	
	delete pmc;

	return nc;
}

int			img_impose_segments(Bimage* p, Bimage* pmask, vector<Bsuperpixel>& seg)
{
	long			i, j, cc;
	
	for ( i=j=0; i<p->image_size(); ++i )
		for ( cc=0; cc<p->channels(); cc++, j++ )
			p->set(j, seg[(long)(*pmask)[i]].channel(cc));
	
	return 0;
}

long		segment_lowest_neighbor(vector<Bsuperpixel>& seg, long i, long cc)
{
	long		j, k, ln(i);
	double		min(seg[i].channel(cc));
	for ( j=0; j<NNEIGHBOR && ( k = seg[i].neighbor(j) ) >= 0; ++j ) {
		if ( seg[k].channel(cc) < min ) {
			min = seg[k].channel(cc);
			ln = k;
		}
	}
	return ln;
}

int			img_impose_lowest_neighbor(Bimage* p, Bimage* pmask, vector<Bsuperpixel>& seg)
{
	long			i, j, is, cc;
	
	for ( i=j=0; i<p->image_size(); ++i )
		for ( cc=0; cc<p->channels(); cc++, j++ ) {
			is = segment_lowest_neighbor(seg, (*pmask)[i], cc);
			p->set(j, seg[is].channel(cc));
		}
	
	return 0;
}

long		segment_highest_neighbor(vector<Bsuperpixel>& seg, long i, long cc)
{
	long		j, k, ln(i);
	double		max(seg[i].channel(cc));
	for ( j=0; j<NNEIGHBOR && ( k = seg[i].neighbor(j) ) >= 0; ++j ) {
		if ( seg[k].channel(cc) > max ) {
			max = seg[k].channel(cc);
			ln = k;
		}
	}
	return ln;
}

int			img_impose_highest_neighbor(Bimage* p, Bimage* pmask, vector<Bsuperpixel>& seg)
{
	long			i, j, is, cc;
	
	for ( i=j=0; i<p->image_size(); ++i )
		for ( cc=0; cc<p->channels(); cc++, j++ ) {
			is = segment_highest_neighbor(seg, (*pmask)[i], cc);
			p->set(j, seg[is].channel(cc));
		}
	
	return 0;
}

int			img_impose_difference_from_highest_neighbor(Bimage* p, Bimage* pmask, vector<Bsuperpixel>& seg)
{
	long			i, j, is, cc;
	double			v;
	
	for ( i=j=0; i<p->image_size(); ++i )
		for ( cc=0; cc<p->channels(); cc++, j++ ) {
			is = segment_highest_neighbor(seg, (*pmask)[i], cc);
			v = seg[(long)(*pmask)[i]].channel(cc);
			p->set(j, seg[is].channel(cc) - v);
		}
	
	return 0;
}

long		segment_setup_neighbors(long i, vector<Bsuperpixel>& seg, double dmax)
{
	long			nseg(seg.size());
	long			j, nn(0);
	double			d;
	
	for ( j=0; j<nseg; j++ ) if ( j != i ) {
		d = distance(seg[i].coordinates(), seg[j].coordinates());
		if ( d < dmax ) {
			seg[i].add_neighbor(j);
			nn++;
		}
	}
	
	return nn;
}

long		segments_setup_neighbors(vector<Bsuperpixel>& seg, long step)
{
	long			nseg(seg.size());
	long			i, j, nn(0);
	double			dmax(1.9*step);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG segments_setup_neighbors: step=" << step << endl;

	for ( i=0; i<nseg; ++i ) seg[i].clear_neighbors();

#ifdef HAVE_GCD
	dispatch_apply(nseg, dispatch_get_global_queue(0, 0), ^(size_t i){
		segment_setup_neighbors(i, seg, dmax);
	});
#else
#pragma omp parallel for
	for ( i=0; i<nseg; ++i )
		segment_setup_neighbors(i, seg, dmax);
#endif

	for ( i=0; i<nseg; ++i )
		for ( j=0; j<NNEIGHBOR && seg[i].neighbor(j) >= 0; ++j ) nn++;

	if ( verbose & VERB_FULL ) {
		for ( i=0; i<nseg; ++i ) {
			cout << seg[i].index() << tab << seg[i].coordinates();
			for ( j=0; j<NNEIGHBOR && seg[i].neighbor(j) >= 0; ++j ) cout << tab << seg[i].neighbor(j);
			cout << endl;
		}
	}
	
	return nn;
}
			
/**
@brief 	Create superpixels from a multilevel mask.
@param	cc			number of channels in the original image.
@param 	step		distance limit to determine neigbors.
@return vector<Bsuperpixel>	array of superpixels.
**/
vector<Bsuperpixel>	Bimage::superpixels_from_mask(long cc, long step)
{
	long			nseg = (long) (max + 1.9);
	
	if ( verbose )
		cout << "Setting up " << nseg << " segments" << endl;
	
	vector<Bsuperpixel>	seg(nseg);

	long			i, j, nn(0);
	Vector3<long>	coor3;
	vector<double>	v(cc, 0);
	
	for ( i=0; i<nseg; ++i ) {
		seg[i].index(i);
		seg[i].channels(cc);
	}
	
	for ( i=0; i<image_size(); ++i ) {
		j = (*this)[i];
		coor3 = coordinates(i);
		seg[j].add(coor3[0], coor3[1], coor3[2], v);
	}
		
	for ( i=0; i<nseg; ++i )
		seg[i].average();

	nn = segments_setup_neighbors(seg, step);
	
	if ( verbose )
		cout << "Total neigbors: " << nn << endl << endl;

	return seg;
}

/**
@brief 	Segment the image into superpixels.
@param 	step		initial superpixel intervals.
@param	colorweight	weight of color differences compared to spatial distances.
@param	iterations	maximum number of iterations.
@param	stop		stopping condition as a percent of voxel changes.
@return vector<Bsuperpixel>	array of superpixels.

	The segment array contents are:
	0		count
	1-3		coordinates
	4+		channels
	last	color weight for the segment - maximum squared distance

	The mask is linked to the original image for return.
**/
vector<Bsuperpixel>	Bimage::superpixels(long step, double colorweight, long iterations, double stop)
{
	Bimage*			pmask = tile_mask(step);

	long			i, nc(image_size()), ncp(0), nstop(image_size()*stop/100);
	
	if ( verbose ) {
		cout << "SLIC superpixel segmentation:" << endl;
		cout << "Step size:                           " << step << endl;
		cout << "Color weight:                        " << colorweight << endl;
		cout << "Maximum number of iterations:        " << iterations << endl;
		cout << "Stopping condition:                  " << stop << " % (" << nstop << ")" << endl << endl;
	}
	
	vector<Bsuperpixel>	seg = pmask->superpixels_from_mask(c, step);

	vector<long>		vstep = {step,step,step};
	for ( i=0; i<3; ++i ) vstep[i] *= 2;

	if ( vstep[2] > z ) vstep[2] = z;

	superpixels_update(pmask, vstep, colorweight, seg);

//	write_img("t.tif", pmask, 0);

	time_t			t = time(NULL);
	
	if ( verbose )
		cout << "Iter\tChange\t∆\tTime" << endl;
	
	for ( i=0; i<iterations && nc > nstop; ++i ) {
		nc = img_assign_segments(this, pmask, vstep, seg);	// Slow step
		if ( verbose )
			cout << i+1 << tab << nc << tab << nc-ncp << tab << (long)(time(NULL) - t) << endl;
//		write_img("t.tif", pmask);
		ncp = nc;
		superpixels_update(pmask, vstep, colorweight, seg);
		segments_setup_neighbors(seg, step);
	}

	next = pmask;

	return seg;
}

Bimage*		img_unbin_update(Bimage* p)
{
	if ( !p->next ) return 0;
	
	long			i, j;
	Vector3<long>	c;
	Bimage*			pb = NULL;
	
	for ( pb = p; pb->next; p = pb, pb = pb->next ) ;
	
	for ( i=0; i<p->image_size(); ++i ) {
		c = p->coordinates(i)/2;
		j = pb->index(c, 0);
		p->set(i, (*pb)[j]);
	}
	
	delete pb;
	p->next = NULL;
	
	return p;
}

vector<Bsuperpixel>	Bimage::superpixels(long step, double colorweight, long iterations, long bin_level, double stop)
{
	Bimage*			pmask = tile_mask(step);
		
	if ( verbose ) {
		cout << "SLIC superpixel segmentation:" << endl;
		cout << "Step size:                           " << step << endl;
		cout << "Bin levels:                          " << bin_level << endl;
		cout << "Color weight:                        " << colorweight << endl;
		cout << "Maximum number of iterations:        " << iterations << endl;
		cout << "Stopping condition:                  " << stop << " %" << endl << endl;
	}
	
	vector<Bsuperpixel>	seg = pmask->superpixels_from_mask(c, step);

	long			i, bin(1);
	double			w = ( z > 1 )? 1.0/8.0: 1.0/4.0;
	vector<long>	vstep = {step,step,step}, vbstep(3);
	for ( i=0; i<3; ++i ) vstep[i] *= 2;
	if ( vstep[2] > z ) vstep[2] = z;

	superpixels_update(pmask, vstep, colorweight, seg);

//	write_img("t.tif", pmask, 0);

	time_t			t = time(NULL);

	Bimage*			pb = this;
	Bimage*			pm = pmask;
	for ( i=1; i<bin_level; ++i ) {
		bin *= 2;
//		cout << "creating bin level " << i << endl;
		pb->next = pb->bin_copy(2);
		pm->next = pm->bin_copy(2);
		pb = pb->next;
		pm = pm->next;
		pm->multiply(w);
	}

	write_img("t.tif", pm, 0);

	long			nc, ncp(0), nstop(pb->image_size()*stop/100);
	long			bstep(step/bin);
	for ( i=0; i<3; ++i ) vbstep[i] = vstep[i]/bin;
	if ( vbstep[2] < 1 ) vbstep[2] = 1;
	for ( auto it = seg.begin(); it != seg.end(); ++it )
		it->scale(1.0/bin);

	while ( bin ) {
		nc = pb->image_size();
		
//		cout << "size = " << pb->size() << endl;
//		cout << "step = " << bstep << tab << vbstep[1] << endl;
		
		if ( verbose ) {
			cout << "Binning:   " << bin << endl;
			cout << "Iter\tChange\t∆\tTime" << endl;
		}

		for ( i=0; i<iterations && nc > nstop; ++i ) {
			nc = img_assign_segments(pb, pm, vbstep, seg);	// Slow step
			if ( verbose )
				cout << i+1 << tab << nc << tab << nc-ncp << tab << (long)(time(NULL) - t) << endl;
			ncp = nc;
			pb->superpixels_update(pm, vbstep, colorweight, seg);
			segments_setup_neighbors(seg, bstep);
		}
		
		bin /= 2;
		
		if ( bin ) {
			bstep *= 2;
			for ( i=0; i<3; ++i ) {
				vbstep[i] *= 2;
				if ( size()[i] > 1 ) ncp *= 2;
			}
			if ( z == 1 ) vbstep[2] = 1;
			for ( auto it = seg.begin(); it != seg.end(); ++it )
				it->scale(2.0);
			for ( pb = this; pb->next && pb->next->next; pb = pb->next ) ;
			delete pb->next;
			pb->next = NULL;
			pm = img_unbin_update(pmask);
//			pb->superpixels_update(pm, vbstep, colorweight, seg);
//			segments_setup_neighbors(seg, bstep);
		}
		
		iterations /= 2;
	}

	next = pmask;
	
	if ( verbose ) cout << endl;

	return seg;
}

/**
@brief 	Impose superpixel features onto an image.
@param 	*pmask		mask defining superpixels.
@param	seg			array of superpixels.
@param	impose		flag to select feature
@return 0			.

	Features to select:
	0	none
	1	segment average
	2	lowest neighboring segment
	3	difference from lowest neighbor
	4	higest neigboring segment
	5	difference from highest neighbor
**/
int			Bimage::impose_superpixels(Bimage* pmask, vector<Bsuperpixel>& seg, int impose)
{
	if ( impose < 1 ) return 0;
	if ( !pmask ) {
		cerr << "Error in impose_superpixels: No multi-level mask!" << endl;
		return -1;
	}
	
	long			i, j, cc, is, isn;
	double			v, vv;

	change_type(Float);
//	pmask->information();
	
	if ( next ) delete next;
	next = copy();
	
	if ( verbose )
		cout << "Imposing segments using feature type " << impose << endl;
	
	for ( i=j=0; i<image_size(); ++i ) {
		for ( cc=0; cc<channels(); cc++, j++ ) {
			is = isn = (long)(*pmask)[i];
			if ( impose == 2 || impose == 3 ) {
				isn = segment_lowest_neighbor(seg, is, cc);
			} else if ( impose == 4 || impose == 5 ) {
				isn = segment_highest_neighbor(seg, is, cc);
			}
			if ( impose == 3 || impose == 5 ) {
				v = seg[isn].channel(cc) - seg[is].channel(cc);
				vv = seg[isn].variance(cc) - seg[is].variance(cc);
			} else {
				v = seg[isn].channel(cc);
				vv = seg[isn].variance(cc);
			}
			set(j, v);
			next->set(j, vv);
		}
	}
	
	return 0;
}


/**
@brief 	Color an image based on the volumes of regions in a multi-level mask.
@param 	*pmask		mask to limit areas for segmentation.
@return Bimage*		new colored image.

**/
/*Bimage*		Bimage::region_color_by_size(Bimage* pmask)
{	
	Bimage*		pcol = new Bimage(UCharacter, TRGB, size(), n);
	pcol->sampling(sampling(0));
	pcol->origin(image->origin());
	RGB<unsigned char>*	color = (RGB<unsigned char> *) pcol->data_pointer();

	if ( verbose )
		cout << "Coloring regions by size" << endl;
		
	long				i, j, nseg((long)(pmask->maximum() + 1.9));
	long				vmin(datasize), vmax(0);
	double				vavg(0), vstd(0);
	RGB<unsigned char>	black(0,0,0);
	int*				vol = new int[nseg];
	
	for ( i=0; i<nseg; ++i ) vol[i] = 0;
	
	for ( i=0; i<datasize; ++i )
		if ( (*pmask)[i] > -1 ) vol[(long)(*pmask)[i]]++;
	
	for ( i=1; i<nseg; ++i ) {
		if ( vmin > vol[i] ) vmin = vol[i];
		if ( vmax < vol[i] ) vmax = vol[i];
		vavg += vol[i];
		vstd += vol[i] * vol[i];
	}
	vavg /= nseg - 1;
	vstd /= nseg - 1;
	vstd -= vavg * vavg;
	vstd = sqrt(vstd);
	if ( vmin < vavg - 3*vstd ) vmin = vavg - 3*vstd;
	if ( vmax > vavg + 3*vstd ) vmax = vavg + 3*vstd;
	
	if ( verbose )
		cout << "Coloring range:                 " << vmin << " - " << vmax << endl;
	
	for ( i=0; i<datasize; ++i ) {
		j = (*pmask)[i];
		if ( j > 0 )
			color[i].spectrum(vol[j], vmin, vmax);
		else
			color[i] = black;
	}
	
	delete[] vol;
	
	return pcol;
}
*/

