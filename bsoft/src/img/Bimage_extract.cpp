/**
@file	Bimage_extract.cpp
@brief	Library routines to extract parts of images
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20210609
**/

#include "Bimage.h"
#include "spline.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Extracts one sub-image into new image structure.
@param 	nn			image number to extract.
@return Bimage*		the new image structure, NULL if copy failed.
**/
Bimage*		Bimage::extract(long nn)
{
	long			j, k, imgsize(x*y*z*c);
	
	Bimage*			img = copy_header(1);
	img->image[0] = image[nn];
	img->data_alloc();
	
	for ( j=nn*imgsize, k=0; k<imgsize; j++, k++ ) img->set(k, (*this)[j]);
	
	return img;
}

/**
@brief 	Extracts a set of sub-images into new image structure.
@param 	n1			first sub-image to extract.
@param 	n2			last sub-image to extract.
@return Bimage*		the new image structure, NULL if copy failed.
**/
Bimage*		Bimage::extract(long n1, long n2)
{
	if ( n1 > n2 ) swap(n1, n2);
	if ( n1 >= n ) n1 = n-1;
	if ( n2 >= n ) n2 = n-1;
	
	long			i, j, k, nn, ne, imgsize(x*y*z*c);
	
	Bimage*			img = copy_header(n2-n1+1);
	img->data_alloc();
	
	for ( ne=0, nn=n1; nn<=n2; nn++, ne++ ) {
		img->image[ne] = image[nn];
		for ( i=0, j=nn*imgsize, k=ne*imgsize; i<imgsize; i++, j++, k++ )
			img->set(k, (*this)[j]);
	}
	
	return img;
}

/**
@brief Extracts a region of one sub-image into new image structure.
@param 	nn			image number to extract.
@param 	coords		extraction start.
@param 	ext_size	extraction size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return Bimage*		the new image structure, NULL if copy failed.
**/
Bimage*		Bimage::extract(long nn, Vector3<long> coords, Vector3<long> ext_size,
						int fill_type, double fill)
{
	long			i, j, xx, yy, zz, cc;
	long			oldx, oldy, oldz;
	
	Bimage*			img = copy_header(1);
	img->size(ext_size);
	img->data_alloc();

	if ( fill_type == FILL_AVERAGE ) fill = avg;
	if ( fill_type == FILL_BACKGROUND ) {
		if ( fabs(background(nn)) < 1e-20  ) calculate_background(nn);
		fill = background(nn);
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::extract: image size = " << size() << endl;
		cout << "DEBUG Bimage::extract: extraction size = " << img->size() << endl;
		cout << "DEBUG Bimage::extract: fill = " << fill << endl;
	}

	for ( zz=0, oldz=coords[2]; zz<img->z; zz++, oldz++ ) {
		for ( yy=0, oldy=coords[1]; yy<img->y; yy++, oldy++ ) {
			for ( xx=0, oldx=coords[0]; xx<img->x; xx++, oldx++ ) {
				i = img->index(xx, yy, zz);
				if ( within_boundaries(oldx, oldy, oldz) ) {
					j = index(oldx, oldy, oldz, nn);
					for ( cc=0; cc<c; cc++, i++, j++ ) img->set(i, (*this)[j]);
				} else {
					for ( cc=0; cc<c; cc++, i++ ) img->set(i, fill);
				}
			}
		}
	}
	
	img->origin(coords);
	
	return img;
}

/**
@brief 	Extracts a region of one sub-image into new image structure.
@param 	nn			image number to extract.
@param 	loc			extraction location.
@param 	size		extraction size.
@param 	ori			extraction origin.
@param	mat			orientation matrix.
@return Bimage*		the new image structure, NULL if copy failed.
**/
Bimage*		Bimage::extract(long nn, Vector3<double> loc, Vector3<long> size,
				Vector3<double> ori, Matrix3 mat)
{
	if ( ori[0] == 0 && ori[1] == 0 && ori[2] == 0 )
		ori = size/2;
	
	long			i, xx, yy, zz, cc;
	Vector3<double>	v, vold;
	
	Bimage*			img = copy_header(1);
	img->size(size);
	img->image[0] = image[nn];
	img->origin(ori);
	img->data_alloc();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::extract: size=" << size << endl;
	
	for ( i=zz=0; zz<img->z; zz++ ) {
		v[2] = double(zz) - ori[2];
		for ( yy=0; yy<img->y; yy++ ) {
			v[1] = double(yy) - ori[1];
			for ( xx=0; xx<img->x; xx++ ) {
				v[0] = double(xx) - ori[0];
				vold = (mat * v) + loc;
				for ( cc=0; cc<c; cc++, i++ )
					img->set(i, interpolate(cc, vold, nn, background(nn)));
			}
		}
	}
	
	return img;
}

/**
@brief 	Extracts a region of one sub-image into new image structure with wrapping.
@param 	nn			image number to extract.
@param 	loc			extraction location.
@param 	size		extraction size.
@param 	ori			extraction origin.
@param	mat			orientation matrix.
@return Bimage*		the new image structure, NULL if copy failed.
**/
Bimage*		Bimage::extract_wrap(long nn, Vector3<double> loc, Vector3<long> size,
				Vector3<double> ori, Matrix3 mat)
{
	long			i, xx, yy, zz, cc;
	Vector3<double>	v, vold;
	
	Bimage*			img = copy_header(1);
	img->size(size);
	img->image[0] = image[nn];
	img->origin(ori);
	img->data_alloc();
	
	for ( i=zz=0; zz<img->z; zz++ ) {
		v[2] = zz - ori[2];
		for ( yy=0; yy<img->y; yy++ ) {
			v[1] = yy - ori[1];
			for ( xx=0; xx<img->x; xx++ ) {
				v[0] = xx - ori[0];
				vold = (mat * v) + loc;
				for ( cc=0; cc<c; cc++, i++ )
					img->set(i, interpolate_wrap(cc, vold, nn));
			}
		}
	}
	
	return img;
}

/**
@brief 	Extracts a shell from an image into a new image.
@param 	nn			image number from which to extract, -1 indicates all images.
@param 	minrad		minimum shell radius.
@param 	maxrad		maximum shell radius.
@return Bimage*		extracted image.

	A single image (only one sub-image) is extracted from a given sub-image
	in an image structure, starting at a specified point and with a 
	specified size.
	The old data is not affected.
	Statistics for the extracted image are calculated.

**/
Bimage*		Bimage::extract_shell(long nn, double minrad, double maxrad)
{
	long   		i, j, xx, yy, zz, cc;
	double			rx2, ry2, rz2, r2;
	
	if ( minrad > maxrad ) swap(minrad, maxrad);
	if ( minrad < 0 ) minrad = 0;
	if ( maxrad > x ) maxrad = x - 1;
	
	double			minrad2 = minrad*minrad;
	double			maxrad2 = maxrad*maxrad;
	
	Bimage*			pex = copy_header(1);
	pex->image[0] = image[nn];
	pex->data_alloc();
	pex->fill(avg);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Extracting subimage:            " << nn << endl;
		cout << "Shell:                          " << minrad << " " << maxrad << endl;
		cout << endl;
	}
	
	for ( zz=0; zz<z; zz++ ) {
		if ( z > 1 ) {
			rz2 = zz - image[nn].origin()[2];
			rz2 *= rz2;
		} else rz2 = 0;
		for ( yy=0; yy<y; yy++ ) {
			ry2 = yy - image[nn].origin()[1];
			ry2 *= ry2;
			for ( xx=0; xx<x; xx++ ) {
				rx2 = xx - image[nn].origin()[0];
				rx2 *= rx2;
				r2 = rx2 + ry2 + rz2;
				if ( r2 >= minrad2 && r2 <= maxrad2 ) {
					i = pex->index(xx, yy, zz);
					j = index(xx, yy, zz, nn);
					for ( cc=0; cc<c; cc++, i++, j++ ) pex->set(i, (*this)[j]);
				}
			}
		}
	}
	
	pex->statistics();
	
	return pex;
}

/**
@brief 	Generates a set of tile coordinates for an image.
@param 	&start			3-vector start for first tile to be extracted.
@param 	&region			3-vector size of part of image to be extracted (0 = whole image).
@param 	&tile_size		3-vector size of extracted image.
@param 	&step_size		3-vector size of intervals between tiles.
@param 	exceed			flag to allow tiles to exceed the input image size.
@return vector<Vector3<long>>	set of coordinates.

	Calculating the coordinates of tiles of a specified size and within
	a specified region within the image. The overlap is specified by
	the step size and a flag indicates the option to exceed the
	bounds of the image.
	The old data is not affected.

**/
vector<Vector3<long>>	Bimage::tile_coordinates(Vector3<long>& start,
				Vector3<long>& region, Vector3<long>& tile_size,
				Vector3<long>& step_size, int exceed)
{
	vector<Vector3<long>>	coords;

	if ( tile_size.length() < 1 ) {
		error_show("Bimage::tile_coordinates", __FILE__, __LINE__);
		return coords;
	}
	
	long			i, xx, yy, zz;
	
	if ( region.length() < 2 ) region = size();
	
	if ( start[0] > x ) start[0] = 0;
	if ( start[1] > y ) start[1] = 0;
	if ( start[2] > z ) start[2] = 0;
	
	if ( region[0] > x - start[0] ) region[0] = x - start[0];
	if ( region[1] > y - start[1] ) region[1] = y - start[1];
	if ( region[2] > z - start[2] ) region[2] = z - start[2];
	
	if ( !exceed ) {
		region -= tile_size - 1;
		region = region.max(tile_size);
	}

//	tile_size = tile_size.min(region);
	
	if ( step_size[0] <= 0 ) step_size[0] = tile_size[0];
	if ( step_size[1] <= 0 ) step_size[1] = tile_size[1];
	if ( step_size[2] <= 0 ) step_size[2] = tile_size[2];
	
//	calculate_background();
	
	// Count the number of tiles to be extracted
	long			ntiles(0);
	Vector3<long>	tile_pattern((region+step_size-1)/step_size);
	for ( zz=0; zz<region[2]; zz+=step_size[2] )
		for ( yy=0; yy<region[1]; yy+=step_size[1] )
			for ( xx=0; xx<region[0]; xx+=step_size[0] )
				ntiles++;

	coords.resize(ntiles);
	
	for ( i=zz=0; zz<region[2]; zz+=step_size[2] )
		for ( yy=0; yy<region[1]; yy+=step_size[1] )
			for ( xx=0; xx<region[0]; xx+=step_size[0], i++ )
				coords[i] = Vector3<long>(xx + start[0], yy + start[1], zz + start[2]);

	if ( verbose & VERB_FULL ) {
		cout << "Generating " << ntiles << " tile coordinates:" << endl;
		cout << "Start:                          " << start << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Step size:                      " << step_size << endl;
		cout << "Tile pattern:                   " << tile_pattern << endl;
		cout << endl;
	}
	
	return coords;
}

/**
@brief 	Generates a set of tile coordinates to fit in the image dimensions.
@param 	tile_size			3-vector size of extracted image.
@param 	&step_size			3-vector size of intervals between tiles.
@return vector<Vector3<long>>	set of coordinates.

	Calculating the coordinates of tiles of a specified size to fit within
	the bounds of the image. The overlap is specified by
	the given step size, which is adjusted to fit the image bounds.
	The old data is not affected.

**/
vector<Vector3<long>>	Bimage::tile_coordinates(Vector3<long> tile_size, Vector3<long>& step_size)
{
	vector<Vector3<long>>	coords;

	if ( tile_size.volume() < 1 ) {
		error_show("Bimage::tile_coordinates", __FILE__, __LINE__);
		return coords;
	}
	
	long			i, xx, yy, zz;
	
	if ( step_size[0] <= 0 ) step_size[0] = tile_size[0];
	if ( step_size[1] <= 0 ) step_size[1] = tile_size[1];
	if ( step_size[2] <= 0 ) step_size[2] = tile_size[2];
	
	Vector3<long>	nt((size() - 1)/step_size + 1);
	nt = nt.min(size());
	
	Vector3<long>	of(size() - (tile_size + step_size*(nt - 1))), pof;
	if ( verbose & VERB_FULL )
		cout << step_size << tab << of << endl;
	
	while ( of.length() > 1 && of != pof ) {
		pof = of;
		step_size += of/nt;
		of = size() - (tile_size + step_size*(nt - 1));
		if ( verbose & VERB_FULL )
			cout << step_size << tab << of << endl;
	}
	
	Vector3<long>	t, fin(size()-tile_size+step_size);
	for ( i=zz=0; zz<fin[2]; zz+=step_size[2] ) {
		t[2] = zz;
		if ( zz + tile_size[2] > z ) t[2] = z - tile_size[2];
		for ( yy=0; yy<fin[1]; yy+=step_size[1] ) {
			t[1] = yy;
			if ( yy + tile_size[1] > y ) t[1] = y - tile_size[1];
			for ( xx=0; xx<fin[0]; xx+=step_size[0], i++ ) {
				t[0] = xx;
				if ( xx + tile_size[0] > x ) t[0] = x - tile_size[0];
				coords.push_back(t);
			}
		}
	}

	if ( verbose & VERB_FULL ) {
		cout << "Generating " << coords.size() << " tile coordinates:" << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Step size:                      " << step_size << endl;
		cout << "Tiles:                          " << nt << endl;
		cout << "Offset:                         " << of << endl;
		cout << endl;
		for ( auto it: coords )
			cout << it << endl;
	}
	
	return coords;
}

/**
@brief 	Extracts a set of tiles at specified positions from an image into a new image.
@param 	nn			image number from which to extract, -1 indicates all images.
@param 	coords 		coordinates for the tile origins.
@param 	tile_size	3-vector size of extracted image.
@return Bimage* p	extracted set of sub-images.

	A set of tiles of specified size are extracted from a given sub-image 
	in an image structure, using given coordinates for the tiles, which have
	to fit into the image size. The origins of the tiles are inserted into 
	the sub-image origin fields.
	The old data is not affected.

**/
Bimage*		Bimage::extract_tiles(long nn, vector<Vector3<long>>& coords, Vector3<long> tile_size)
{
	if ( tile_size.length() < 1 ) {
		error_show("Bimage::extract_tiles", __FILE__, __LINE__);
		cerr << "Error: Size of subimage is less than one!" << endl;
		return NULL;
	}
	
//	tile_size = tile_size.min(size());
	
	if ( nn >= n ) nn = n - 1;
	
	long			ntiles(coords.size());
	Bimage*			pex = copy_header(ntiles);
	
	pex->size(tile_size);
	pex->page_size(pex->size());
	
	pex->data_alloc();
//	pex->fill(avg);
	
	if ( verbose & VERB_FULL ) {
		cout << "Extracting " << pex->images() << " tiles using coordinates:" << endl;
		cout << "Tile size:                      " << pex->size() << endl;
		cout << endl;
	}
	
#ifdef HAVE_GCD
	dispatch_apply(ntiles, dispatch_get_global_queue(0, 0), ^(size_t nex){
		Bimage*		pt = extract(nn, coords[nex], tile_size, FILL_BACKGROUND);
		pex->replace(nex, pt);
		delete pt;
	});
#else
#pragma omp parallel for
	for ( long nex=0; nex<ntiles; nex++ ) {
		Bimage*		pt = extract(nn, coords[nex], tile_size, FILL_BACKGROUND);
		pex->replace(nex, pt);
		delete pt;
	}
#endif

	pex->statistics();
		
	return pex;
}

/**
@brief 	Extracts a set of tiles from an image into a new image.
@param 	nn			image number from which to extract.
@param 	start		3-vector start for first tile to be extracted.
@param 	region		3-vector size of part of image to be extracted (0 = whole image).
@param 	tile_size	3-vector size of extracted image.
@param	step_size	3-vector size of tile intervals.
@param 	exceed		flag to allow tiles to exceed the input image size.
@return Bimage*		extracted set of sub-images.

	A set of tiles of specified size are extracted from a given sub-image
	in an image structure, starting from a point and generating
	as many tiles as would fit into the image size given. If the image size
	given is zero, then the whole image is used. The origins of the tiles
	are inserted into the sub-image origin fields.
	The old data is not affected.

**/
Bimage*		Bimage::extract_tiles(long nn, Vector3<long> start,
				Vector3<long> region, Vector3<long> tile_size,
				Vector3<long> step_size, int exceed)
{
	vector<Vector3<long>>	coords = tile_coordinates(start, region, tile_size, step_size, exceed);

	if ( verbose & VERB_FULL )
		cout << "Extracting " << coords.size() << " tiles:" << endl << endl;
	
	Bimage*		pex = extract_tiles(nn, coords, tile_size);
	
	return pex;
}

long		find_optimal_tiling(long target, long tile_size, double fraction=0.2)
{
	if ( fraction > 0.5 ) fraction = 0.5;
	
	long		overlap(fraction*tile_size), i, imin(overlap/2), imax(overlap*2);
	long		step(tile_size - overlap);
	long		num(target/step);
	long		d, dmin(target);
	
	for ( i = imin; i < imax; i++ ) {
		step = tile_size - i;
		num = target/step;
		d = num*step + i - target;
		if ( d >= 0 && dmin > d ) {
			overlap = i;
			dmin = d;
			if ( verbose & VERB_FULL )
				cout << overlap << tab << num << tab << step << tab << d << endl;
		}
	}

	return overlap;
}

/**
@brief 	Extracts a set of tiles from an image into a new image.
@param 	nn			image number from which to extract.
@param 	tile_size	3-vector size of extracted image.
@param 	fraction	overlap fraction to aim for.
@return Bimage*		extracted set of sub-images.

	A set of tiles of specified size are extracted from a given sub-image
	in an image structure. The tiles are located with overlap (~20% of tile width)
	with the aim of covering the whole image but not exceeding it.
	The origins of the tiles are inserted into the sub-image origin fields.
	The old data is not affected.

**/
Bimage*		Bimage::extract_tiles(long nn, Vector3<long> tile_size, double fraction)
{
	tile_size = tile_size.min(size());
	
	Vector3<long>		num, overlap, step, icoor;
	Vector3<long>		maxori(size() - tile_size);
	
	for ( long j=0; j<3; j++ )
		if ( size()[j] > 1 ) overlap[j] = find_optimal_tiling(size()[j], tile_size[j], fraction);

	step = tile_size - overlap;
	
	num = size()/step;
	num = num.max(1);

	if ( verbose & VERB_FULL ) {
		cout << size() << endl;
		cout << num * step + overlap << endl;
	}
	
	long				i, ntiles(num.volume());
	
	vector<Vector3<long>>	coords(ntiles);
	
	for ( i=icoor[2]=0; icoor[2]<num[2]; icoor[2]++ )
		for ( icoor[1]=0; icoor[1]<num[1]; icoor[1]++ )
			for ( icoor[0]=0; icoor[0]<num[0]; icoor[0]++, i++ ) {
				coords[i] = icoor * step;
				coords[i] = coords[i].min(maxori);
			}
				

	if ( verbose & VERB_FULL ) {
		cout << "Extracting " << ntiles << " tiles:" << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << "Number of tiles:                " << num << endl;
		cout << "Overlap:                        " << overlap << endl << endl;
	}
	
	Bimage*		pex = extract_tiles(nn, coords, tile_size);
	
	return pex;
}

/**
@brief 	Extracts a stack of tiles at a specified position from an image into a new image.
@param 	coords 		coordinates for the tile origins.
@param 	tile_size	3-vector size of extracted image.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return Bimage* p		extracted stack of tiles.

	A stack of tiles of specified size are extracted from all the sub-images
	in an image structure, using given coordinates for the tiles.
	The origins of the tiles are inserted into the sub-image origin fields.
	The old data is not affected.

**/
Bimage*		Bimage::extract_tile_stack(Vector3<long> coords, Vector3<long> tile_size, int fill_type, double fill)
{
	Bimage*		p = copy_header();
	p->size(tile_size);
	p->page_size(size());
	p->data_alloc();
//	p->fill(bkg);

	for ( long nn=0; nn<n; ++nn ) {
		Bimage*		pt = extract(nn, coords, tile_size, fill_type, fill);
//		pt->origin(pt->size()/2);
		p->replace(nn, pt);
		delete pt;
	}

	return p;
}


/**
@brief 	Extracts stacks of tiles at a specified positions from an image into an array of new images.
@brief 	Extracts a set of tiles at specified positions from an image into a new image.
@param 	coords 			coordinates for the tile origins.
@param 	tile_size		3-vector size of extracted image.
@return Bimage**			extracted stacks of tiles.

	Stacks of tiles of specified size are extracted from all the sub-images
	in an image structure, using given coordinates for the tiles.
	The origins of the tiles are inserted into the sub-image origin fields.
	The old data is not affected.

**/
Bimage**	Bimage::extract_tile_stacks(vector<Vector3<long>>& coords, Vector3<long> tile_size)
{
	if ( tile_size.length() < 1 ) {
		error_show("Bimage::extract_tile_stacks", __FILE__, __LINE__);
		cerr << "Error: Size of tile is less than one!" << endl;
//		return NULL;
	}
	
	calculate_background();
	
	tile_size = tile_size.min(size());
	
	long			ntiles(coords.size());
	
	Bimage**		parr = new Bimage*[ntiles];

	if ( verbose & VERB_FULL ) {
		cout << "Extracting " << ntiles << " tiles using coordinates:" << endl;
		cout << "Tile size:                      " << tile_size << endl;
		cout << endl;
	}
	
#ifdef HAVE_GCD
	dispatch_apply(ntiles, dispatch_get_global_queue(0, 0), ^(size_t t){
		parr[t] = extract_tile_stack(coords[t], tile_size, FILL_BACKGROUND);
	});
#else
#pragma omp parallel for
	for ( long t=0; t<ntiles; ++t ) {
		parr[t] = extract_tile_stack(coords[t], tile_size, FILL_BACKGROUND);
	}
#endif

	return parr;
}


/**
@brief 	Extracts a line from an image into a new image.
@param 	nn			image number from which to extract.
@param 	start		3-vector start of line.
@param 	end			3-vector end of line.
@param 	width		width of integration perpendicular to line.
@return Bimage*		extracted line image (1D).

	The line values extracted are the interpolated values along the vector
	defined by start and end coordinates..

**/
Bimage*		Bimage::extract_line(long nn, Vector3<double> start, Vector3<double> end, long width)
{
	if ( width < 1 ) width = 1;
	
	long			i, j, wh(width/2);
	Vector3<double>	vec(end - start), vecp, coor;
	long			len(vec.length() + 1);
	
	vec.normalize();
	vecp = Vector3<double>(vec[1], vec[0], vec[2]);
	
	Bimage*			pline = new Bimage(Float, TSimple, len, 1, 1, 1);
	
	for ( i=0; i<len; i++ ) {
		coor = start + vec*i;
		for ( j=-wh; j<=wh; j++ )
			pline->add(i, interpolate(coor + vecp*j, nn));
		pline->set(i, (*pline)[i]/width);
	}
	
	pline->statistics();
	
	return pline;
}

/**
@brief 	Extracts a tetrahedral part of the image.
@param 	*tet		four 3-value vectors defining the tetrahedron.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return int			0.

	All voxels outside a tetrahedron defined by four points or vectors 
	are set to a given fill value.
	The new data replaces the old data.

**/
Bimage*		Bimage::extract_tetrahedron(Vector3<double>* tet, int fill_type, double fill)
{
	if ( fill_type == FILL_AVERAGE ) fill = avg;
	
    long		     	i, j, k, nn, xx, yy, zz, cc;
	
	Vector3<double>		center;
	Vector3<double> 	normal[4];
	double				origin[4];
	
	for ( i=0; i<4; i++ ) center = center + tet[i];
	center /= 4.0;
	
	normal[0] = tet[0].normal(tet[1], tet[2]);
	normal[1] = tet[1].normal(tet[2], tet[3]);
	normal[2] = tet[2].normal(tet[3], tet[0]);
	normal[3] = tet[3].normal(tet[0], tet[1]);
	
	// The distance of a point from a plane is given by the inner product sum of
	// the normal to the plane and the vector from that point to any point on the plane
	for ( i=0; i<4; i++ ) {
		origin[i] = normal[i].scalar(tet[i]);
		if ( normal[i].scalar(center) < origin[i] ) {
			normal[i] = -normal[i];
			origin[i] = -origin[i];
		}
	}
	
	Bimage*				ptet = copy();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Extracting a tetrahedron:" << endl;
		cout << "Vertices:                       " << tet[0] << endl;
		cout << "                                " << tet[1] << endl;
		cout << "                                " << tet[2] << endl;
		cout << "                                " << tet[3] << endl;
		cout << "Center:                         " << center << endl;
		cout << "Normals:                        " << normal[0] << " (" << origin[0] << ")" << endl;
		cout << "                                " << normal[1] << " (" << origin[1] << ")" << endl;
		cout << "                                " << normal[2] << " (" << origin[2] << ")" << endl;
		cout << "                                " << normal[3] << " (" << origin[3] << ")" << endl;
		cout << endl; 
	} else if ( verbose & VERB_LABEL )
	    cout << "Extracting a tetrahedron" << endl << endl;
	
	for ( i=nn=0; nn<n; nn++ ) {
		if ( fill_type == FILL_BACKGROUND ) fill = background(nn);
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				for ( xx=0; xx<x; xx++, i++ ) {
					for ( j=0; j<4; j++ ) {
						if ( normal[j][0]*xx + normal[j][1]*yy + normal[j][2]*zz < origin[j] )
							for ( cc=0, k=i*c; cc<c; cc++, k++ )
								ptet->set(k, fill);
					}
				}
			}
		}
	}	
	
	return 0;
}

/**
@brief 	Extracts orthogonal views around a voxel.
@param 	nn			image number to extract.
@param 	voxel		voxel of intersection.
@param 	ext_size	size of slices to extract.
@return Bimage*		extracted slice for 2D and 3 slices for 3D.

	Only the desired region is extracted from the original image.
	The fill value is taken from the image background.

**/
Bimage*		Bimage::orthogonal_slices(long nn, Vector3<long> voxel,
				Vector3<long> ext_size)
{
	if ( voxel[0] < 0 || voxel[1] < 0 || voxel[2] < 0 ) return NULL;
	if ( voxel[0] >= x || voxel[1] >= y || voxel[2] >= z ) return NULL;
	
	if ( ext_size.volume() < 1 ) ext_size = size();
	
	int 				xf, yf;
	long				i, j, iy, xx, yy, zz, cc;
	long				xo, yo, zo;
	long				h, nh(1);
	if ( z > 1 ) nh = 3;
	long				nx = (long) (ext_size[0]);
	long				ny = (long) (ext_size[1]);
	
	Vector3<long>		start(voxel - (ext_size * 0.5));
	Vector3<long>		end(start + ext_size);
	Vector3<double>		ori;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Extracting orthogonal views:" << endl;
		cout << "Start:                          " << start << endl;
		cout << "End:                            " << end << endl;
		cout << "Fill value:                     " << image[nn].background() << endl;
		cout << endl;
	}

    Bimage*     		porth = new Bimage(datatype, compoundtype, nx, ny, 1, nh);
	porth->sampling(sampling(0));
	porth->fill(image[nn].background());
	
	for ( h=0; h<nh; h++ ) {
		porth->background(h, background(nn));
		if ( h%2 == 0 ) ori[0] = image[nn].origin()[0] - start[0];
		else ori[0] = image[nn].origin()[2] - start[2];
		if ( h < 2 ) ori[1] = image[nn].origin()[1] - start[1];
		else ori[1] = image[nn].origin()[2] - start[2];
		porth->origin(h, ori);
		if ( h == 0 ) porth->image[h].view(0,0,1,0);
		else if ( h == 1 ) porth->image[h].view(1,0,0,0);
		else porth->image[h].view(0,1,0,0);
		xo = voxel[0];
		yo = voxel[1];
		zo = voxel[2];
		for ( yy=0; yy<ny; yy++ ) {
			yf = 1;
			if ( h < 2 ) {
				yo = (long) (yy + start[1]);
				if ( yo < 0 ) yf = 0;
				else if ( yo >= y ) yf = 0;
			} else {
				zz = yy;
				zo = (long) (zz + start[2]);
				if ( zo < 0 ) yf = 0;
				else if ( zo >= z ) yf = 0;
			}
			iy = (h*ny + yy)*nx;
			for ( xx=0; xx<nx; xx++ ) {
				xf = 1;
				if ( h%2 == 0 ) {
					xo = (long) (xx + start[0]);
					if ( xo < 0 ) xf = 0;
					else if ( xo >= x ) xf = 0;
				} else {
					zz = xx;
					zo = (long) (zz + start[2]);
					if ( zo < 0 ) xf = 0;
					else if ( zo >= z ) xf = 0;
				}
				i = (iy + xx)*porth->c;
				if ( xf && yf ) {
					j = index(xo, yo, zo, nn);
					for ( cc=0; cc<porth->c; cc++, i++, j++ ) porth->set(i, (*this)[j]);
				}
			}
		}
    }

	porth->statistics();
	
    return porth;
}

/**
@brief 	Extracts orthogonal views around a voxel and create a montage.
@param 	voxel		voxel of intersection.
@param 	ext_size	size of slices to extract.
@param 	pad			padding size.
@param 	fill_type	FILL_AVERAGE, FILL_BACKGROUND, FILL_USER.
@param 	fill		fill value.
@return Bimage*		extracted slice for 2D and 3 slices for 3D.

	Only the desired region is extracted from the original image.
	The fill value is taken from the image background.

**/
Bimage*		Bimage::orthogonal_montage(Vector3<long> voxel, Vector3<long> ext_size, int pad, int fill_type, double fill)
{
	if ( voxel[0] < 0 || voxel[1] < 0 || voxel[2] < 0 ) return NULL;
	if ( voxel[0] >= x || voxel[1] >= y || voxel[2] >= z ) return NULL;
	
	if ( ext_size.volume() < 1 ) ext_size = size();
	
	bool 				xf, yf;
	long				i, j, iy, iz, xx, yy, zz, cc, nn;
	long				xo, yo, zo;
	
	Vector3<long>		start(voxel - (ext_size * 0.5));
	Vector3<long>		end(start + ext_size);
	Vector3<double>		ori;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Extracting orthogonal views:" << endl;
		cout << "Voxel:                          " << voxel << endl;
		cout << "Start:                          " << start << endl;
		cout << "End:                            " << end << endl;
		cout << "Padding:                        " << pad << endl;
		cout << endl;
	}

    Bimage*     		porth = new Bimage(datatype, compoundtype, ext_size[0] + ext_size[2] + pad,
    						ext_size[1] + ext_size[2] + pad, 1, n);
    						
	for ( nn=0; nn<n; ++nn ) {
		if ( fill_type == FILL_BACKGROUND ) fill = image[nn].background();
		else if ( fill_type == FILL_AVERAGE ) fill = image[nn].average();
//		cout << "Fill: " << fill << endl;
		porth->background(nn, background(nn));
		porth->sampling(sampling(nn));
		ori = image[nn].origin();
		ori[1] += ext_size[2] + pad;
		porth->image[nn].origin(ori);
		xo = voxel[0];
		yo = voxel[1];
		zo = voxel[2];
		iz = nn*porth->sizeY();
		// xy slice
		for ( yy=0; yy<porth->y; ++yy ) {
			iy = (iz + yy)*porth->x;
			yf = 0;
			if ( yy < ext_size[2] ) {
				yo = voxel[1];
				zz = yy;
				zo = (long) (zz + start[2]);
				if ( zo >= 0 && zo < z ) yf = 1;
			} else if ( yy >= ext_size[2] + pad ) {
				zo = voxel[2];
				yo = (long) (yy + start[1] - ext_size[2] - pad);
				if ( yo >= 0 && yo < y ) yf = 1;
			}
			for ( xx=0; xx<porth->x; ++xx ) {
				i = (iy + xx)*porth->c;
				xf = 0;
				if ( xx < ext_size[0] ) {
					xo = (long) (xx + start[0]);
					if ( xo >= 0 && xo < x ) xf = 1;
				} else if ( xx >= ext_size[0] + pad ) {
					xo = voxel[0];
					zz = xx - ext_size[0] - pad;
					zo = (long) (zz + start[2]);
					if ( zo >= 0 && zo < z && yy >= ext_size[2] ) xf = 1;
				}
				j = index(xo, yo, zo, nn);
				for ( cc=0; cc<porth->c; ++cc, ++i, ++j ) {
					if ( xf && yf ) porth->set(i, (*this)[j]);
					else porth->set(i, fill);
				}
			}
		}
    }

	porth->statistics();
	
    return porth;
}


int			Bimage::extract_show_chunk(Bimage* pshow, int aflag, long i, long len)
{
	long				nn(sn);
	long				zz(sz);
	double				scale(ss);

	if ( compoundtype == TComplex ) {
		zz += z/2;
		if ( zz >= z ) zz -= z;
	}
	
	if ( scale >= 1 ) aflag = 0;	// No averaging for enlargements
	
	int 				threshold(0);
	if ( fabs(smax - smin) < 1e-37 ) threshold = 1;
	
	long				ic, j, jz, xx, yy, cc, m, xm(x/2), ym(y/2);
	long				xo, yo;
	double				xf, yf;
	double				iscale(1.0/scale), shift = (scale < 1)? (iscale - 1)/2: 0;
	
	cc = (pshow->c < c)? c: pshow->c;
 	vector<double>		value(cc,0);
	double				dispval[4] = {0,0,0,0};
    double				dscale = 255.0/(smax - smin);
	long				ds(pshow->x*pshow->y);
	RGB<double>			rgb;
	CMYK<double>		cmyk;
	
	// Note that the y-axis is flipped
	xx = i%pshow->x;
	yy = i/pshow->x;
	jz = (nn*z + zz)*y;

	m = i+len;
	if ( m > ds ) m = ds;
	for ( ; i<m; i++, xx++ ) {
		if ( xx >= pshow->x ) {
			xx = 0;
			yy++;
		}
		if ( compoundtype == TComplex ) {
			xf = xx*iscale - xm;
			yf = (pshow->y - yy - 1)*iscale - ym;
			if ( xf < 0 ) xf += x;
			if ( yf < 0 ) yf += y;
		} else {
			xf = xx*iscale;
			yf = (pshow->y - yy - 1)*iscale;
		}
		if ( xf >= 0 || xf <= x-1 || yf >=0 || yf <= y-1 ) {
			yo = (long) (yf + shift);
			xo = (long) (xf + shift);
			if ( datatype == Bit ) {
				j = (jz + yo)*page_size()[0] + xo;
				pshow->set(i, 255 * (*this)[j]);
				if ( dscale < 0 ) pshow->set(i, 255 - (*pshow)[i]);
			} else {
				if ( aflag ) {
					for ( cc=0; cc<c; cc++ ) value[cc] = average2D(cc, xf, yf, zz, nn, iscale);
				} else {
					for ( j = ((jz + yo)*x + xo)*c, cc=0; cc<c; cc++, j++ ) 
						value[cc] = (*this)[j];
				}
				if ( compoundtype == TComplex ) {
					dispval[0] = value[0]*value[0] + value[1]*value[1];
				} else if ( compoundtype == TCMYK ) {
					if ( datatype == UCharacter ) for ( cc=0; cc<4; cc++ ) value[cc] /= 255;
					cmyk = CMYK<double>(value);
					rgb = RGB<double>(cmyk);
					for ( cc=0; cc<3; cc++ ) dispval[cc] = 255*rgb[cc];
				} else {
					for ( cc=0; cc<pshow->c; cc++ ) dispval[cc] = value[cc];
				}
				if ( threshold ) {
					for ( cc=0; cc<pshow->c && cc<3; cc++ ) {
						if ( dispval[cc] > smin ) dispval[cc] = 255;
						else dispval[cc] = 0;
					}
				} else {
					for ( cc=0; cc<pshow->c; cc++ ) {
						if ( cc < 3 ) dispval[cc] = floor(dscale*(dispval[cc]-smin)+0.5);
						if ( dispval[cc] > 255 ) dispval[cc] = 255;
						else if ( dispval[cc] < 0 ) dispval[cc] = 0;
					}
				}
//				cout << dispval[3] << tab;
				if ( pshow->c > 3 ) dispval[3] = 255;
//				if ( dispval[3] > 100 ) dispval[3] = 255;
				for ( cc=0, ic = i*pshow->c; cc<pshow->c; cc++, ic++ )
					pshow->set(ic, dispval[cc]);
			}
		}
    }
//	cout << endl;

//	delete [] value;
	
    return 0;
}

/**
@brief 	Converts a slice from an image to a 2D plane for display.
@param 	aflag		averaging flag.
@return Bimage*		extracted slice.

	Only the desired image slice is extracted from the original image and
	resized based on the given scale argument.
	The image, slice and scale values are encoded in the image structure.
	The dynamic range is rescaled using the image display minimum and maximum:
		new_data = data*255/(max-min)
	with truncation of the data below 0 and above 255.
	Bit data are converted to 0 (black) and 255 (white).
	RGB data are rescaled as for gray scale images.
	Complex data types are converted to intensities.

**/
Bimage*		Bimage::extract_show(int aflag)
{
	CompoundType		ctype;
	
	switch ( compoundtype ) {
		case TSimple: ctype = TSimple; break;
		case TComplex: ctype = TSimple; break;
		case TVector2:
		case TVector3:
		case TMulti:
		case TCMYK:
		case TRGB: ctype = TRGB; break;
		case TView:
		case TRGBA: ctype = TRGBA; break;
		default: ctype = compoundtype;
	}
	
    Bimage*     	pshow = new Bimage(UCharacter, ctype, 
							(long) (ss*x), (long) (ss*y), 1, 1);
	pshow->sampling(sampling(0)/ss);
	
	long			ds = pshow->x*pshow->y;
//	long			chunk_size = get_chunk_size(ds, c);
	long			chunk_size = get_chunk_size(ds);

#ifdef HAVE_GCD
	dispatch_apply((ds - 1)/chunk_size + 1, dispatch_get_global_queue(0, 0), ^(size_t i){
		extract_show_chunk(pshow, aflag, i*chunk_size, chunk_size);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<ds; i+=chunk_size )
		extract_show_chunk(pshow, aflag, i, chunk_size);
#endif

    return pshow;
}


/**
@brief 	Extracts a region from an image to magnify.
@param 	nn			image number to extract.
@param 	center		center of region to magnify.
@param 	ext_size	size of region to magnify.
@param 	scale		dimensional scaling.
@return Bimage*		extracted slice for 2D and 3 slices for 3D.

	Only the desired region is extracted from the original image and
	resized based on the given scale argument.
	The dynamic range is rescaled using the image display minimum and 
	maximum:
		new_data = data*255/(max-min)
	with truncation of the data below 0 and above 255.
	Bit data are converted to 0 (black) and 255 (white).
	RGB data are rescaled as for gray scale images.
	Complex data types are converted to intensities.

**/
Bimage*		Bimage::extract_magnify(long nn, Vector3<long> center,
				Vector3<long> ext_size, double scale)
{
	if ( compoundtype == TComplex ) {
		center = vector3_set_PBC(center - size()/2, size());
	} else {
		if ( center[0] < 0 || center[1] < 0 || center[2] < 0 ) return NULL;
		if ( center[0] >= x || center[1] >= y || center[2] >= z ) return NULL;
	}
	
	int 				threshold(0), xf, yf;
	if ( fabs(smax - smin) < 1e-37 ) threshold = 1;
	
	long				i, j, iy, jz, xx, yy, zz, cc;
	long				xo, yo, zo;
	long				h, nh(1);
	if ( z > 1 ) nh = 3;
	long				nx = (long) (scale*ext_size[0]);
	long				ny = (long) (scale*ext_size[1]);
	double				iscalex = 1.0/scale;
	double				iscaley = 1.0/scale;
	double				iscalez = 1.0/scale;
	RGB<double>			rgb;
	CMYK<double>		cmyk;
	
	Vector3<long>		start = center - (ext_size * 0.5);
	Vector3<long>		end = start + ext_size;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::extract_magnify: start=" << start << endl;
		cout << "DEBUG Bimage::extract_magnify: end=" << end << endl;
		cout << "DEBUG Bimage::extract_magnify: size: nx=" << nx << " ny=" << ny << " c=" << c << endl;
	}
	
	CompoundType		ctype;
	
	switch ( compoundtype ) {
		case TSimple: ctype = TSimple; break;
		case TComplex: ctype = TSimple; break;
		case TVector3:
		case TCMYK:
		case TRGB: ctype = TRGB; break;
		case TView:
		case TRGBA: ctype = TRGBA; break;
		default: ctype = compoundtype;
	}
	
    Bimage*     		pmag = new Bimage(UCharacter, ctype, nx, ny, 1, nh);
	pmag->sampling(sampling(nn)/scale);
	
 	double*				value = new double[c];
	double				dispval[4];    
    
    double				dscale = 255.0/(smax - smin);
	
	// Note that the y-axis is flipped
	for ( h=0; h<nh; h++ ) {
		xo = center[0];
		yo = center[1];
		zo = center[2];
		for ( yy=0; yy<ny; yy++ ) {
			yf = 1;
			if ( h < 2 ) {
				yo = (long) (yy*iscaley + start[1]);
				if ( compoundtype == TComplex ) {
					if ( yo < 0 ) yo += y;
					else if ( yo >= y ) yo -= y;
				} else {
					if ( yo < 0 ) yf = 0;
					else if ( yo >= y ) yf = 0;
				}
			} else {
				zz = yy;
				zo = (long) (zz*iscalez + start[2]);
				if ( compoundtype == TComplex ) {
					if ( zo < 0 ) zo += z;
					else if ( zo >= z ) yo -= z;
				} else {
					if ( zo < 0 ) yf = 0;
					else if ( zo >= z ) yf = 0;
				}
			}
			iy = (h*ny + ny - 1 - yy)*nx;
			for ( xx=0; xx<nx; xx++ ) {
				xf = 1;
				if ( h%2 == 0 ) {
					xo = (long) (xx*iscalex + start[0]);
					if ( compoundtype == TComplex ) {
						if ( xo < 0 ) xo += x;
						else if ( xo >= x ) xo -= x;
					} else {
						if ( xo < 0 ) xf = 0;
						else if ( xo >= x ) xf = 0;
					}
				} else {
					zz = xx;
					zo = (long) (zz*iscalez + start[2]);
					if ( compoundtype == TComplex ) {
						if ( zo < 0 ) zo += z;
						else if ( zo >= z ) zo -= z;
					} else {
						if ( zo < 0 ) xf = 0;
						else if ( zo >= z ) xf = 0;
					}
				}
				if ( xf && yf ) {
					i = (iy + xx)*pmag->c;
					jz = (nn*z + zo)*y;
					if ( datatype == Bit ) {
						j = (jz + yo)*page_size()[0] + xo;
						pmag->set(i, 255 * (*this)[j]);
						if ( dscale < 0 ) pmag->set(i, 255 - (*pmag)[i]);
					} else {
						j = ((jz + yo)*x + xo)*c;
						for ( cc=0; cc<c; cc++, j++ ) value[cc] = (*this)[j];
						if ( compoundtype == TComplex ) {
							dispval[0] = value[0]*value[0] + value[1]*value[1];
						} else if ( compoundtype == TCMYK ) {
							if ( datatype == UCharacter ) for ( cc=0; cc<4; cc++ ) value[cc] /= 255;
							cmyk = CMYK<double>(value);
							rgb = RGB<double>(cmyk);
							for ( cc=0; cc<3; cc++ ) dispval[cc] = 255*rgb[cc];
						} else {
							for ( cc=0; cc<pmag->c; cc++ ) dispval[cc] = value[cc];
						}
						if ( threshold ) {
							for ( cc=0; cc<pmag->c; cc++ ) {
								if ( dispval[cc] > smin ) dispval[cc] = 255;
								else dispval[cc] = 0;
							}
						} else {
							for ( cc=0; cc<pmag->c; cc++ ) {
								dispval[cc] = floor(dscale*(dispval[cc]-smin)+0.5);
								if ( dispval[cc] > 255 ) dispval[cc] = 255;
								else if ( dispval[cc] < 0 ) dispval[cc] = 0;
							}
						}
						for ( cc=0; cc<pmag->c; cc++, i++ ) pmag->set(i, dispval[cc]);
					}
				}
			}
		}
    }
	
	delete [] value;
	
    return pmag;
}

/**
@brief 	Extracts a given slice or slices from an image.
@param 	nz			slice number to extract.
@return Bimage*		extracted slice(s).

	Only the desired slices are extracted from the original image.
	The dynamic range is rescaled using the image display minimum and 
	maximum:
		new_data = data*255/(max-min)
	with truncation of the data below 0 and above 255.
	Bit data are converted to 0 (black) and 255 (white).
	RGB data are rescaled as for gray scale images.
	Complex data types are converted to intensities.

**/
Bimage*		Bimage::extract_slice(long nz)
{
	int 				threshold(0);
	if ( fabs(smax - smin) < 1e-37 ) threshold = 1;
	
	long				i, j, xx, yy, cc, nn;
	RGB<double>			rgb;
	CMYK<double>		cmyk;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::extract_slice: nz=" << nz << endl;
	
	CompoundType		ctype;
	
	switch ( compoundtype ) {
		case TSimple: ctype = TSimple; break;
		case TComplex: ctype = TSimple; break;
		case TVector3:
		case TCMYK:
		case TRGB: ctype = TRGB; break;
		case TView:
		case TRGBA: ctype = TRGBA; break;
		default: ctype = compoundtype;
	}
	
    Bimage*     		pslice = new Bimage(UCharacter, ctype, x, y, 1, n);
	pslice->sampling(sampling(0));
	
 	double*				value = new double[c];
	double				dispval[4];    
    
    double				dscale = 255.0/(smax - smin);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::extract_slice: dscale=" << dscale << endl;

	// Note that the y-axis is flipped
	for ( i=nn=0; nn<n; nn++ ) {
		j = (nn*z + nz)*x*y*c;
		for ( yy=0; yy<y; yy++ ) {
			for ( xx=0; xx<x; xx++ ) {
					if ( datatype == Bit ) {
						j *= page_size()[0];
						pslice->set(i, 255 * (*this)[j]);
						if ( dscale < 0 ) pslice->set(i, 255 - (*pslice)[i]);
					} else {
						for ( cc=0; cc<c; cc++, j++ ) value[cc] = (*this)[j];
						if ( compoundtype == TComplex ) {
							dispval[0] = value[0]*value[0] + value[1]*value[1];
						} else if ( compoundtype == TCMYK ) {
							if ( datatype == UCharacter ) for ( cc=0; cc<4; cc++ ) value[cc] /= 255;
							cmyk = CMYK<double>(value);
							rgb = RGB<double>(cmyk);
							for ( cc=0; cc<3; cc++ ) dispval[cc] = 255*rgb[cc];
						} else {
							for ( cc=0; cc<pslice->c; cc++ ) dispval[cc] = value[cc];
						}
						if ( threshold ) {
							for ( cc=0; cc<pslice->c; cc++ ) {
								if ( dispval[cc] > smin ) dispval[cc] = 255;
								else dispval[cc] = 0;
							}
						} else {
							for ( cc=0; cc<pslice->c; cc++ ) {
								dispval[cc] = floor(dscale*(dispval[cc]-smin)+0.5);
								if ( dispval[cc] > 255 ) dispval[cc] = 255;
								else if ( dispval[cc] < 0 ) dispval[cc] = 0;
							}
						}
						for ( cc=0; cc<pslice->c; cc++, i++ ) pslice->set(i, dispval[cc]);
					}
			}
		}
    }
	
	delete [] value;
	
    return pslice;
}

/**
@brief 	Extracts a filament defined by a series of coordinates.
@param 	img_num		image number from which to extract filament.
@param 	width		width of filament image to extract.
@param 	axis		helical axis alignment: x=1, y=2, z=3.
@param 	nspline		number of coordinates in spline curve.
@param 	*spline		spline curve.
@return Bimage*		one filament image.

	A single filament is extracted and returned.

**/
Bimage*		Bimage::extract_filament(long img_num, double width,
				int axis, long nspline, Vector3<double>* spline)
{
	if ( !d.uc ) return NULL;
	if ( nspline < 1 ) return NULL;
	
	if ( axis < 1 || axis > 3 ) {	// Default axis selection
		if ( z > 1 ) axis = 3;
		else axis = 1;
	}
		
	Vector3<long>	filsize((long)width, (long)width, (long)width);
	filsize[axis-1] = nspline;
	filsize = filsize.min(size());
	
	Bimage*			pfil = new Bimage(Float, compoundtype, filsize, 1);
	pfil->sampling(sampling(0));
	pfil->label(label());

	long			i, xx, yy, zz, s(0);
	double			w, halfwidth(width/2);
	Vector3<double>	vref, v, vf;
	Matrix3			mat(1);
	
	vref[axis-1] = 1;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage:extract_filament: img_num=" << img_num << " width=" << width << " nspline=" << nspline << endl;
		cout << "DEBUG Bimage:extract_filament: axis=" << axis << " size=" << filsize << endl;
	}
	
	for ( i=zz=0; zz<pfil->z; zz++ ) {
		if ( axis == 3 ) s = zz;
		else if ( pfil->z > 1 ) vf[2] = zz - halfwidth;
		for ( yy=0; yy<pfil->y; yy++ ) {
			if ( axis == 2 ) s = yy;
			else if ( pfil->y > 1 ) vf[1] = yy - halfwidth;
			for ( xx=0; xx<pfil->x; xx++, i++ ) {
				if ( axis == 1 ) s = xx;
				else if ( pfil->x > 1 ) vf[0] = xx - halfwidth;
				if ( s < nspline ) {
					if ( s > 0 ) v = spline[s] - spline[s-1];
					else v = spline[1] - spline[0];
					v.normalize();
					mat = Matrix3(vref, v);
//					cout << mat << endl;
					v = mat * vf;
					v += spline[s];
					w = interpolate(v, img_num, background(img_num));
					pfil->add(i, w);
				}
			}
		}
	}
	
	return pfil;
}

