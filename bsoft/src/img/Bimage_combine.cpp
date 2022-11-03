/**
@file	Bimage_combine.cpp
@brief	Functions to combine two images in various ways
@author Bernard Heymann
@date	Created: 19990219
@date	Modified: 20150902
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Sums an array of images with their FOM blocks.
@param 	m			number of images in the array.
@param 	**p			array of images.

	The images must all have the same dimensions.

**/
void		Bimage::sum(long m, Bimage** p)
{
	long			i, j, k, cc, ds(x*y*z*n);
	
	clear();
	
	vector<double>			b(n);
	for ( i=0; i<n; i++ ) b[i] = 0;
	
	for ( j=0; j<m; j++ ) {
		for ( i=0; i<ds; i++ )
			for ( cc=0, k=i*c; cc<c; cc++, k++ )
				add(k, (*p[j])[k]);
		for ( i=0; i<n; i++ ) b[i] += p[j]->image[i].background();
	}

	for ( i=0; i<n; i++ ) image[i].background(b[i]);
	
	statistics();
}

/**
@brief 	Catenates an array of images of the same size into a multi-image structure.
@param 	m			number of images in the array.
@param 	**p			array of images.

	The images must all have the same dimensions.

**/
void		Bimage::catenate(long m, Bimage** p)
{
	long			i, j, k, l, ds;
	long			imgsize(c*x*y*z);

	for ( i=n=0; i<m; i++ ) if ( p[i] ) n += p[i]->n;
	
	if ( n < 1 ) {
		error_show("Error: No images to concatenate!", __FILE__, __LINE__);
		return;
	}

	data_type(p[0]->datatype);
	compoundtype = p[0]->compoundtype;
	c = p[0]->c;
	x = p[0]->x;
	y = p[0]->y;
	z = p[0]->z;

	images(n);
	
	data_alloc_and_clear();

	if ( verbose & VERB_PROCESS )
		cout << "Catenating " << n << " images" << endl;

	for ( i=j=0; i<m; i++ ) if ( p[i] ) {
		ds = imgsize*p[i]->n;
		for ( k=j*imgsize, l=0; l<ds; k++, l++ ) set(k, (*p[i])[l]);
		for ( k=0; k<p[i]->n; k++, j++ )
			image[j] = p[i]->image[k];
	}
	
	statistics();
}


/**
@brief 	Blends the two images, creating a new set of sub-images.
@param 	*p			second image.
@param 	number		number of images in the series.
@return Bimage* 	new image structure, NULL if error.

	A number of images are created by blending the two input images in
	different ratios. The input images become the first and last sub-images of 
	the new image structure, with the intermediate images changing over from
	the first to the last:
			new_data = (1-fraction)*data1 + fraction*data2
	where fraction = index/(number - 1)
	At least 3 new sub-images are packed into the new image. 
	All of the header information in the first image is copied into the new image.
	Both images are converted to floating point.

**/
Bimage* 	Bimage::blend(Bimage* p, long number)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::blend", __FILE__, __LINE__);
		return NULL;
	}
	
    if ( number < 3 ) number = 3;
	
    if ( verbose & VERB_LABEL )
	    cout << "Blending two images to create a series of " << number << " sub-images" << endl;
	
    long			i, j, nn;
	double			fraction1, fraction2;
	
	calculate_background();
	p->calculate_background();
    
	// Copy the first image into the new image structure with a new number of sub-images
	Bimage* 		pnu = copy_header(number);
	pnu->data_type(Float);
	pnu->data_alloc();
	
//	cout << pnu->meta_data() << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_morph_blend: datasize=" << datasize << endl;
	
    for ( i=nn=0; nn<number; nn++ ) {
		fraction2 = nn*1.0/(number - 1.0);
		fraction1 = 1 - fraction2;
		if ( verbose & VERB_PROCESS )
			cout << "Calculating image " << nn+1 << ": dnew = " << fraction1
				<< " x image1 + " << fraction2 << " x image2" << endl;
		for ( j=0; j<datasize; j++, i++ )
    		pnu->set(i, fraction1*(*this)[j] + fraction2*(*p)[j]);
//		pnu->image[nn].background(fraction1*image->background() +
//				fraction2*p->image->background());
		pnu->image[nn].background(fraction1*background(long(0)) +
								  fraction2*p->background(long(0)));
	}
	if ( verbose & VERB_PROCESS ) cout << endl;
	
	pnu->statistics();

    return pnu;
}

/**
@brief 	Places a small image into a large image.
@param	nn			sub-image.
@param 	*p			image to place.
@param 	loc			location in large image of small image origin.
@param 	radius		radial mask to transfer small image.
@param 	scale		density scale to apply to second image.
@param 	shift		density shift to apply to second image.
@param 	operation	operation to apply.
@return int			0, <0 if error.

	The small image is placed with its origin at given origin in the large image.
	The second image is scaled and shifted before placing into the first:
		image1 = image1 + image2*scale + shift
	Both images are converted to floating point.
	The operation can be selected:
		0	simple addition.
		1	replace if smaller.
		2	replace if larger.
	Requirement: Both images must have the same pixel size.

**/
int 		Bimage::place(long nn, Bimage* p, Vector3<double> loc,
				double radius, double scale, double shift, int operation)
{
	if ( radius <= 0 ) radius = (x > y)? x: y;

    long			i, si, xx, yy, zz;
	long			sz, sy, sx;
	double			sx2, sy2, sz2, r2(radius*radius), v;
	Vector3<double>	start = loc - p->image[nn].origin() + 0.5;
	Vector3<long>	lo(start);
	Vector3<long>	hi(p->size());
	hi += lo;
	lo = lo.max(0);
	hi = hi.max(1);
	hi = hi.min(size());
	
	if ( verbose & VERB_FULL ) {
		cout << "Placing a small image into a larger image:" << endl;
		cout << "Small image size:               " << p->size() << endl;
		cout << "Location:                       " << loc << endl;
		cout << "Low limits:                     " << lo << endl;
		cout << "High limits:                    " << hi << endl << endl;
	}
	
	for ( zz=lo[2]; zz<hi[2]; zz++ ) {
		sz = (long) zz - start[2];
		sz2 = p->image[nn].origin()[2] - sz;
		sz2 *= sz2;
		for ( yy=lo[1]; yy<hi[1]; yy++ ) {
			sy = (long) yy - start[1];
			sy2 = p->image[nn].origin()[1] - sy;
			sy2 *= sy2;
			for ( xx=lo[0]; xx<hi[0]; xx++ ) {
				sx = (long) xx - start[0];
				sx2 = p->image[nn].origin()[0] - sx;
				sx2 *= sx2;
				if ( sx2 + sy2 + sz2 <= r2 ) {
					i = index(xx, yy, zz, nn);
					si = p->index(sx, sy, sz);
					v = (*p)[si] * scale + shift;
					switch ( operation ) {
						case 1: if ( v < (*this)[i] ) set(i, v); break;
						case 2: if ( v > (*this)[i] ) set(i, v); break;
						default: add(i, v);
					}
				}
			}
		}
	}
	
	return 0; 
}

/**
@brief 	Packs a tile into a new composite image with addition within overlap.
@param 	*p			image = tiles.
@param 	nn			tile image (can be a sub-image).
@return int			0.

	The tiles are added to the image.
	The contributions at each voxel are counted in a linked image.
	The tile placement is in the sub-image origins.

**/
int			Bimage::place_with_addition(Bimage* p, long nn)
{
	if ( !next ) next = new Bimage(Float, compoundtype, size(), n);
	
	long			i, j, cc, tx, ty, tz, xx, yy, zz;
	float*			count = (float *) next->data_pointer();

	if ( verbose & VERB_FULL ) {
		cout << "Adding tile:" << endl;
		cout << "Location:                       " << p->image[nn].origin() << endl;
	}
	
	zz = p->image[nn].origin()[2];
	for ( tz=0; tz<p->z && zz<z; tz++, zz++ ) {
		yy = p->image[nn].origin()[1];
		for ( ty=0; ty<p->y && yy<y; ty++, yy++ ) {
			xx = p->image[nn].origin()[0];
			for ( tx=0; tx<p->x && xx<x; tx++, xx++ ) {
				i = p->index(0, tx, ty, tz, nn);
				j = index(0, xx, yy, zz, 0);
				for ( cc=0; cc<c; cc++, i++, j++ ) {
					add(j, (*p)[i]);
					count[j] += 1;
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Packs a tile into a new composite image with weighted overlap.
@param 	*p			image = tiles.
@param 	nn			tile image (can be a sub-image).
@return int			0.

	The overlaps between tiles are filled in with a weighted average based
	on a linear transition from one to the other.
	The tile placement is in the sub-image origins.

**/
int			Bimage::place_with_overlap(Bimage* p, long nn)
{
	long			i, j, cc, tx, ty, tz, xx, yy, zz;
	double			wvol;
	Vector3<double>	d, w, imin, imax(p->size());

	// Determine overlap in all dimensions
	for ( i=0; i<p->n; i++ ) if ( i != nn ) {
		d = p->image[nn].origin() - p->image[i].origin();
		if ( d.abs() < p->size() ) {
			if ( d[0] > p->x/2 ) imin[0] = p->x - d[0];
			if ( d[0] < -p->x/2 ) imax[0] = -d[0];
			if ( d[1] > p->y/2 ) imin[1] = p->y - d[1];
			if ( d[1] < -p->y/2 ) imax[1] = -d[1];
			if ( d[2] > p->z/2 ) imin[2] = p->z - d[2];
			if ( d[2] < -p->z/2 ) imax[2] = -d[2];
		}
	}
	imin = imin.max(0);
	imax = imax.min(p->size());
	
	Vector3<double>		wmin(1,1,1);
	Vector3<double>		wmax(1,1,1);
	
	for ( i=0; i<3; ++i ) {
		if ( imin[i] ) wmin[i] = 1.0L/imin[i];
		if ( imax[i] ) wmax[i] = 1.0L/(p->size()[i] - imax[i]);
	}

	if ( verbose & VERB_FULL ) {
		cout << "Placing tile:" << endl;
		cout << "Location:                       " << p->image[nn].origin() << endl;
		cout << "Tile size:                      " << p->size() << endl;
		cout << "Overlap minima:                 " << imin << endl;
		cout << "Overlap maxima:                 " << imax << endl << endl;
	}
	
	zz = p->image[nn].origin()[2];
	for ( tz=0; tz<p->z; tz++, zz++ ) {
		w[2] = 1;
		if ( tz < imin[2] ) w[2] = tz*wmin[2];
		else if ( tz > imax[2] ) w[2] = (p->z - tz)*wmax[2];
		yy = p->image[nn].origin()[1];
		for ( ty=0; ty<p->y; ty++, yy++ ) {
			w[1] = 1;
			if ( ty < imin[1] ) w[1] = ty*wmin[1];
			else if ( ty > imax[1] ) w[1] = (p->y - ty)*wmax[1];
			xx = p->image[nn].origin()[0];
			for ( tx=0; tx<p->x; tx++, xx++ ) {
				w[0] = 1;
				if ( tx < imin[0] ) w[0] = tx*wmin[0];
				else if ( tx > imax[0] ) w[0] = (p->x - tx)*wmax[0];
				i = p->index(0, tx, ty, tz, nn);
				j = index(0, xx, yy, zz, 0);
				wvol = w[0]*w[1]*w[2];
//				wvol = 1;
				for ( cc=0; cc<c; cc++, i++, j++ ) add(j, (*p)[i] * wvol);
			}
		}
	}

	return 0;
}

/**
@brief 	Packs a tile into a new composite image, retaining only the central part.
@param 	*p			image = tiles.
@param 	nn			tile image (can be a sub-image).
@return int			0.

	The overlaps between tiles are divided between the neighboring tiles.
	The tile placement is in the sub-image origins.

**/
int			Bimage::place_central_part(Bimage* p, long nn)
{
	long			i, j, cc, tx, ty, tz, xx, yy, zz, w;
	Vector3<long>	imin, imax(p->size());
	Vector3<double>	d;
	
	// Determine overlap in all dimensions
	for ( i=0; i<p->n; i++ ) if ( i != nn ) {
		d = p->image[nn].origin() - p->image[i].origin();
		for ( j=0; j<3; ++j ) {
			if ( d[j] && fabs(d[j]) <= p->size()[j] ) {
				if ( d[j] > 0 ) imin[j] = (long)(p->size()[j] - d[j])/2;
				else imax[j] = (long)(p->size()[j] - d[j])/2;
			}
		}
	}
	
	// Special cases: The edges of the image
	for ( i=0; i<3; ++i )
		if ( p->image[nn].origin()[i] + p->size()[i] >= size()[i] )
			imax[i] = size()[i] - p->image[nn].origin()[i];		
	
	if ( verbose & VERB_FULL ) {
		cout << "Placing tile:" << endl;
		cout << "Location:                       " << p->image[nn].origin() << endl;
		cout << "Tile size:                      " << p->size() << endl;
		cout << "Overlap minima:                 " << imin << endl;
		cout << "Overlap maxima:                 " << imax << endl << endl;
	}

//	if ( p->image[nn].origin()[0] < 1 ) cout << nn << tab << imin << endl;
	
	zz = long(p->image[nn].origin()[2] + 0.5);
	for ( tz=0; tz<p->z; tz++, zz++ ) {
		w = 1;
		if ( tz < imin[2] ) w = 0;
		else if ( tz > imax[2] ) w = 0;
		if ( w ) {
			yy = long(p->image[nn].origin()[1] + 0.5);
			for ( ty=0; ty<p->y; ty++, yy++ ) {
				w = 1;
				if ( ty < imin[1] ) w = 0;
				else if ( ty > imax[1] ) w = 0;
				if ( w ) {
					xx = long(p->image[nn].origin()[0] + 0.5);
					for ( tx=0; tx<p->x; tx++, xx++ ) {
						w = 1;
						if ( tx < imin[0] ) w = 0;
						else if ( tx > imax[0] ) w = 0;
						if ( w ) {
							i = p->index(tx, ty, tz, nn);
							j = index(xx, yy, zz);
							for ( cc=0; cc<c; cc++, i++, j++ ) set(j, (*p)[i]);
						}
					}
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Assembles overlapping tiles into this image.
@param 	*pt				multi-image containing tiles.
@param 	flag			flag to set overlap handling.
@return int				0.

	The origin in each tile specifies the starting location.
	Where tiles overlap, the integration is determined by the flag:
	0 = simple addition with averaging afterwards.
	1 = a gradual transition from one tile to the other.
	2 = place just the central part of each tile.
	The original data in the image is overwritten.

**/
int			Bimage::assemble_tiles(Bimage* pt, int flag)
{
	clear();
	
	for ( long nn=0; nn<pt->n; nn++ ) {
		if ( flag == 1 ) {
			place_with_overlap(pt, nn);
		} else if ( flag == 2 ) {
			place_central_part(pt, nn);
		} else {
			place_with_addition(pt, nn);
		}
	}
	
	if ( flag == 0 && next ) {
		for ( long i=0; i<data_size(); ++i )
			if ( (*next)[i] > 1 ) set(i, (*this)[i]/(*next)[i]);
		delete next;
		next = NULL;
	}

	return 0;
}

/**
@brief 	Linear least squares fit of two images.
@param 	*p			second image.
@param 	*pmask		mask to limit calculation to a certain region.
@param 	max_exclude	maximum percentage of outlying points to exclude.
@return double 		first R factor.

	The data blocks from two images are fit by a simple linear least squares
	regression algorithm with exclusion of a percentage of outliers:
		image2 = intercept + slope * image1
	The first image is modified to return the difference:
		difference = intercept + slope * image1 - image2
	The residual returned is:
		R = sqrt(sum(difference^2) / sum((image2-avg2)^2))
	Note: A linear fit is not symmetric with respect to the two input
	data sets - the order of the input images determine the output.
	The two data blocks must have the same size and are converted to
	floating point.

**/
double		Bimage::linear_fit(Bimage* p, Bimage* pmask, double max_exclude)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::linear_fit", __FILE__, __LINE__);
		return -1.0;
	}
	
	if ( max_exclude >= 1 ) max_exclude /= 100;	// Assume it is a percentage
	if ( max_exclude > 0.5 ) max_exclude = 0.5; // At least half of the data must be fit
	
    long   			i, j, nn, nm, nfit, best_excl(0);
    long   			imgsize(x*y*z*c), num(imgsize);
    int     	    bin;
    double	  		sx, sx2, sy, sxy, sd, dy, denominator, a, b, R, bestR(1e30), best_a(0), best_b(0);
	double			excl_voxels, diff, maxdiff, diff_cutoff;
    
	int				hist[1024];
	for ( i=0; i<1024; i++ ) hist[i] = 0;
	int*			inc_mask = new int[imgsize];
	for ( i=0; i<imgsize; i++) inc_mask[i] = 1;
	
	if ( pmask ) pmask->change_type(UCharacter);
    
    if ( verbose & VERB_LABEL ) {
	    cout << "Linear least squares fit:       image2 = intercept + slope x image1" << endl;
	    cout << "Maximum voxels to exclude:      " << 100.0*max_exclude << " %" << endl;
	}

	if ( verbose & VERB_FULL )
		cout << " %Set\t Voxels  \t%Excluded\tIntercept\tSlope\t\tR" << endl;
	else if ( verbose & VERB_PROCESS )
		cout << "Intercept\t  Slope\t  R-factor\tExcluded\t%" << endl;
	
	for ( nn=0; nn<n; nn++) {
		best_excl = 0;
		best_a = best_b = 0;
		bestR = 1e30;
		if ( pmask ) {
			nm = nn;
			if ( pmask->images() < n ) nm = 0;
			for ( num=j=0, i=nm*imgsize; j<imgsize; j++,i++ ) {
				if ( (*pmask)[i] ) {
					inc_mask[i] = 1;
					num++;
				} else inc_mask[j] = -1;
			}
		}
		for ( excl_voxels = 0; excl_voxels < max_exclude+0.01; excl_voxels += 0.01 ) {
			nfit = 0;
			a = b = sx = sx2 = sy = sxy = sd = dy = maxdiff = 0;
			for ( j=0, i=nn*imgsize; j<imgsize; j++, i++ ) if ( inc_mask[j] > 0 ) {
				sx += (*this)[i];
				sy += (*p)[i];
				sx2 += (*this)[i]*(*this)[i];
				sxy += (*this)[i]*(*p)[i];
				nfit++;
			}
			denominator = nfit*sx2 - sx*sx;
			if ( denominator ) {
				a = (sx2*sy   - sx*sxy)/denominator;
				b = (nfit*sxy - sx*sy) /denominator;
			}
			sy /= nfit;
    
			for ( j=0, i=nn*imgsize; j<imgsize; j++, i++ ) if ( inc_mask[j] >= 0 ) {
				diff = a + b*(*this)[i] - (*p)[i];
				if ( inc_mask[j] ) {
					dy += diff*diff;
					sd += ((*p)[i]-sy)*((*p)[i]-sy);
				}
				diff = fabs(diff);
				if ( maxdiff < diff ) maxdiff = diff;
			}
			R = 1e30;
			if ( dy < 1e-30 ) R = 0;
			else if ( sd ) R = sqrt(dy/sd);
			if ( bestR > R ) {
				bestR = R;
				best_a = a;
				best_b = b;
				best_excl = num - nfit;
			}
		
			if ( verbose & VERB_FULL )
				cout<< 100*excl_voxels << tab << nfit << tab 
					<< 100.0-nfit*100.0/num << tab << a << tab << b << tab << R << endl;
	
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG img_linear_fit: maxdiff=" << maxdiff << endl;
			if ( maxdiff < 1e-30 ) maxdiff = 1;
			for ( i=0; i<1024; i++ ) hist[i] = 0;
			for ( j=0, i=nn*imgsize; j<imgsize; j++, i++ ) if ( inc_mask[j] >= 0 ) {
				diff = fabs(a + b*(*this)[i] - (*p)[i]);
				bin = (int) (1000*diff/maxdiff);
				if ( bin > 999 ) bin = 999;
				hist[bin]++;
			}
		
			i = 999;
			j = hist[i];
			while ( i > 1 && j < num*excl_voxels ) {
				i--;
				j += hist[i];
			}
			diff_cutoff = i*maxdiff/1000.0;
		
			for ( j=0, i=nn*imgsize; j<imgsize; j++, i++ ) if ( inc_mask[j] >= 0 ) {
				diff = fabs(a + b*(*this)[i] - (*p)[i]);
				inc_mask[j] = 1;
				if ( diff > diff_cutoff ) inc_mask[j] = 0;
			}
		}
	
		for ( j=0, i=nn*imgsize; j<imgsize; j++, i++ )
			set(i, best_a + best_b*(*this)[i] - (*p)[i]);
		
		image[nn].FOM(bestR);
		
		if ( verbose & VERB_PROCESS )
			cout << best_a << tab << best_b << tab << bestR << tab << best_excl << tab << best_excl*100.0/num << " %" << endl;
	}
	
	delete[] inc_mask;
	
	if ( verbose >= VERB_PROCESS ) cout << endl;
    
	return bestR;
}

/**
@brief 	Fits two images by matching the histogram of the second to the first.
@param 	*p			second image.
@param 	bins		number of bins in the histograms.
@return int 		0, <0 if error.

	Both images are converted to floating point.

**/
int 		Bimage::histomatch(Bimage* p, long bins)
{
	if ( !check_if_same_size(p) ) {
		error_show("Bimage::histomatch", __FILE__, __LINE__);
		return -1;
	}
	
    long			i;
	long 			h1, h2;
	double 			denom, fraction;
	int*			hist1 = new int[bins];
	int*			hist2 = new int[bins];
	int*			map = new int[bins];
	double*			dens = new double[bins];
	
	for ( i=0; i<bins; i++ )
		hist1[i] = hist2[i] = 0;
    
    if ( verbose & VERB_LABEL )
	    cout << "Matching histograms in " << bins << " bins:" << endl << endl;
	
    for ( i=0; i<datasize; i++ ) {
		h1 = (long) (bins*((*this)[i] - min)/(max - min));
		h2 = (long) (bins*((*p)[i] - p->minimum())/(p->maximum() - p->minimum()));
		if ( h1 >= bins ) h1 = bins - 1;
		if ( h2 >= bins ) h2 = bins - 1;
		hist1[h1]++;
		hist2[h2]++;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "Histograms done" << endl;

	for ( i=bins-1; i>0; i-- ) {
		hist1[i-1] += hist1[i];
		hist2[i-1] += hist2[i];
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "Integration done" << endl;	
	
	for ( h2=0; h2<bins; h2++ ) {
		h1 = 0;
		while ( ( h1 < bins ) && ( hist1[h1] > hist2[h2] ) ) h1++;
		if ( h1 > 0 ) h1--;
		if ( h2 < bins - 1 )
			denom = hist1[h1+1] - hist1[h1] - hist2[h2+1] + hist2[h2];
		else
			denom = hist1[h1] - hist1[h1-1] - hist2[h2] + hist2[h2-1];
		if ( denom != 0 )
			fraction = (hist2[h2] - hist1[h1])/denom;
		else
			fraction = 0.5;
		map[h2] = h1;
		dens[h2] = (h2 + fraction)*(max - min)/bins + min;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Mapping calculated" << endl;
		cout << "Bin\tHist1\tHist2\tMapping\tDensity" << endl;
		for ( i=0; i<bins; i++ )
			cout << i << tab << hist1[i] << tab << hist2[i] << tab 
				<< map[i] << tab << dens[i] << endl;
	}	
	
    for ( i=0; i<datasize; i++ ) {
		h2 = (int) (bins*((*p)[i] - p->minimum())/(p->maximum() - p->minimum()));
		if ( h2 >= bins ) h2 = bins - 1;
		set(i, dens[map[h2]]);
	}
	
	statistics();
	
	if ( verbose & VERB_DEBUG )
		cout << "Histogram matching done" << endl;
	
	delete[] hist1;
	delete[] hist2;
	delete[] map;
	delete[] dens;
	
	return 0;
}

/**
@brief 	Replaces values for >x/2 with the given image.
@param 	*p			second image.
@return int 			0, <0 if error.

**/
int 		Bimage::replace_half(Bimage* p)
{
	long 		nn, xx, cc, i, j, yz(y*z);
	
	for ( nn=0; nn<n; ++nn )
		for ( j=0; j<yz; ++j )
			for ( xx=x/2, i=((nn*yz + j)*x + xx)*c; xx<x; ++xx )
				for ( cc=0; cc<c; ++cc, ++i )
					set(i, (*p)[i]);
	
	statistics();
	
	return 0;
}

