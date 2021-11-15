/**
@file	Bimage_bin.cpp
@brief	Library routines for binning images
@author Bernard Heymann
@date	Created: 19990904
@date	Modified: 20161017
**/

#include "Bimage.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


int 		Bimage::integer_interpolation(int integer_factor)
{
	return integer_interpolation(integer_factor, 21);
}

/**
@brief 	Interpolates by an integer scale with a density-preserving overlapping kernel.
@param 	integer_factor	integer interpolation factor.
@param 	odd				flag to ensure the dimensions are odd.
@return int 			0.

	An image is interpolated by integer scaling (i.e., 2, 3, 4-fold or 
	more) sometimes referred to as a form of binning.  A kernel is used
	such that it overlaps with its neighbouring positions.  Voxels where
	neighbouring kernel positions overlap contribute to 2, 4 or 8 new
	voxels based on the number of overlapping kernel positions.  Only
	the central voxel is unique to a kernel position.  The kernel is 
	calculated as:
		w(i,j,k) = (1/s^2n)*(s-|s-1-i|)*(s-|s-1-j|)*(s-|s-1-k|)
	where	s is the integer interpolation factor.
			n is the number of dimensions (1D, 2D or 3D).
	The flag determines whether the dimensions are forced to be even or odd:
		0		No forcing
		1		x odd
		2		x even
		4		y odd
		8		y even
		16		z odd
		32		z even
		21		all odd
		42		all even
	and any other combination.

**/
int 		Bimage::integer_interpolation(int integer_factor, int odd)
{
	if ( integer_factor < 2 ) return 0;
	
	// Calculate the new size and kernel weights
	long   			i, j, ix, iy, iz, cc;
	long			nn, wx, wy, wz, xx, yy, zz;
	int				dimensions(0);
	Vector3<long>	nusize(1,1,1);
	Vector3<long>	oldnomori(x/2, y/2, z/2);
	Vector3<long>	nunomori;
	Vector3<long>	kernelsize(1,1,1);
	Vector3<double>	nusam(image->sampling());

	if ( x > 1 ) {
		kernelsize[0] = 2*integer_factor - 1;
		nusize[0] = (x + 1)/integer_factor;
		if ( (odd & 1) && nusize[0]%2 == 0 ) nusize[0]--; // Size must be odd
		if ( (odd & 2) && nusize[0]%2 == 1 ) nusize[0]--; // Size must be even
		nunomori[0] = nusize[0]/2;
		nusam[0] *= integer_factor;
		dimensions++;
	}
	if ( y > 1 ) {
		kernelsize[1] = 2*integer_factor - 1;
		nusize[1] = (y + 1)/integer_factor;
		if ( (odd & 4) && nusize[1]%2 == 0 ) nusize[1]--; // Size must be odd
		if ( (odd & 8) && nusize[1]%2 == 1 ) nusize[1]--; // Size must be even
		nunomori[1] = nusize[1]/2;
		nusam[1] *= integer_factor;
		dimensions++;
	}
	if ( z > 1 ) {
		kernelsize[2] = 2*integer_factor - 1;
		nusize[2] = (z + 1)/integer_factor;
		if ( (odd & 16) && nusize[2]%2 == 0 ) nusize[2]--; // Size must be odd
		if ( (odd & 32) && nusize[2]%2 == 1 ) nusize[2]--; // Size must be even
		nunomori[2] = nusize[2]/2;
		nusam[2] *= integer_factor;
		dimensions++;
	}
	
	double			kernel_scale = 1.0/pow(1.0*integer_factor, 2.0*dimensions);
	vector<double>	kernel(kernelsize[0]*kernelsize[1]*kernelsize[2]);
	
	for ( i=zz=0; zz<kernelsize[2]; zz++ ) {
		wz = (int) (integer_factor - fabs(integer_factor - 1.0 - zz));
		for ( yy=0; yy<kernelsize[1]; yy++ ) {
			wy = (int) (integer_factor - fabs(integer_factor - 1.0 - yy));
			for ( xx=0; xx<kernelsize[0]; xx++, i++ ) {
				wx = (int) (integer_factor - fabs(integer_factor - 1.0 - xx));
				kernel[i] = wx*wy*wz*kernel_scale;
			}
		}
	}
	
	long			datasize = n*nusize[0]*nusize[1]*nusize[2]*c;
	double*			nudata = new double[datasize];
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Special integer interpolation:" << endl;
		cout << "Integer factor:                 " << integer_factor << endl;
		cout << "Interpolation kernel size:      " << kernelsize << endl;
		cout << "New image size:                 " << nusize << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Special integer interpolation" << endl << endl;

	if ( verbose & VERB_FULL ) {
		cout << "Kernel interpolation: Kernel" << endl;
		for ( i=zz=0; zz<kernelsize[2]; zz++ ) {
			for ( yy=0; yy<kernelsize[1]; yy++ ) {
				for ( xx=0; xx<kernelsize[0]; xx++, i++ ) {
					cout << " " << kernel[i];
				}
				cout << endl;
			}
		}
		cout << endl;
	}
	
	// Do the kernel interpolation
	long 			oldx, oldy, oldz, xkernel, ykernel, zkernel;
	long 			zlo, zhi, ylo, yhi, xlo, xhi;
	vector<double>	sum(c);
	
	for ( nn=0; nn<n; nn++ ) {
		for ( iz=0; iz<nusize[2]; iz++ ) {
			oldz = (iz - nunomori[2])*integer_factor + oldnomori[2] - kernelsize[2]/2;
			zlo = oldz;
			if ( zlo < 0 ) zlo = 0;
			zhi = oldz + kernelsize[2];
			if ( zhi > z ) zhi = z;
			for ( iy=0; iy<nusize[1]; iy++ ) {
				oldy = (iy - nunomori[1])*integer_factor + oldnomori[1] - kernelsize[1]/2;
				ylo = oldy;
				if ( ylo < 0 ) ylo = 0;
				yhi = oldy + kernelsize[1];
				if ( yhi > y ) yhi = y;
				for ( ix=0; ix<nusize[0]; ix++ ) {
					oldx = (ix - nunomori[0])*integer_factor + oldnomori[0] - kernelsize[0]/2;
					xlo = oldx;
					if ( xlo < 0 ) xlo = 0;
					xhi = oldx + kernelsize[0];
					if ( xhi > x ) xhi = x;
					for ( cc=0; cc<c; cc++ ) sum[cc] = 0;
					for ( zz=zlo; zz<zhi; zz++ ) {
						zkernel = zz - oldz;
						for ( yy=ylo; yy<yhi; yy++ ) {
							ykernel = yy - oldy;
							for ( xx=xlo; xx<xhi; xx++ ) {
								xkernel = xx - oldx;
								i = index(0, xx, yy, zz, nn);
								j = (zkernel*kernelsize[1] + ykernel)*kernelsize[0] + xkernel;
								for ( cc=0; cc<c; cc++ ) sum[cc] += (*this)[i++] * kernel[j];
								
							}
						}
					}
					i = (((nn*nusize[2] + iz)*nusize[1] + iy)*nusize[0] + ix)*c;
					for ( cc=0; cc<c; cc++ ) nudata[i++] = sum[cc];
				}
			}
		}
		image[nn].origin((image[nn].origin() - oldnomori)/integer_factor + nunomori);
		if ( verbose & VERB_STATS )
			cout << "Image " << nn+1 << " origin:                " << image[nn].origin() << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "Kernel interpolation: Sums done" << endl;
	
	// Set new image parameters
	size(nusize);
	page_size(nusize);
	sampling(nusam);
	
	// Transfer the data to a new block and free the old
	data_alloc();
	for ( i=0; i<datasize; i++ ) set(i, nudata[i]);
	
	delete[] nudata;
	
	statistics();
	
	return 0;
}

int			Bimage::bin(long i, Vector3<long> bk, Bimage* pb)
{
	long   			j, cc, xx, yy, zz;
	long			nn(i/(pb->image_size()*c));
	Vector3<long>	coor = pb->coordinates(i);
	Vector3<long>	bs(bk * coor);
	Vector3<long>	bf(bs + bk);
	long			binsize(0);
	vector<double>	binsum(c);
	for ( cc=0; cc<c; cc++ ) binsum[cc] = 0;

	for ( zz=bs[2]; zz<bf[2] && zz<z; zz++ ) {
		for ( yy=bs[1]; yy<bf[1] && yy<y; yy++ ) {
			for ( xx=bs[0]; xx<bf[0] && xx<x; xx++ ) {
				binsize++;
				j = index(0, xx, yy, zz, nn);
				for ( cc=0; cc<c; cc++ ) binsum[cc] += (*this)[j++];
			}
		}
	}
	
//	for ( cc=0; cc<c; cc++ ) pb->set(i++, binsum[cc]/binsize);
	for ( cc=0; cc<c; cc++ ) pb->set(i++, binsum[cc]);

	return 0;
}

/**
@brief 	Bins by an integer size in place.
@param 	bk		3-value vector of integer bin factors.
@return int 	0.

	An image is binned by an integer size.

**/
int 		Bimage::bin(Vector3<long> bk)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::bin: bk=" << bk << endl;
	
	Bimage*		pb = bin_copy(bk);
	
	// Set new image parameters
	size(pb->size());
	page_size(size());

	for ( long i=0; i<n; i++ ) {
		image[i].origin(pb->image[i].origin());
		image[i].sampling(pb->image[i].sampling());
	}
	
//	cout << "alloc_size=" << alloc_size() << endl;
	// Transfer the data to a new block and free the old
	data_alloc();	
//	cout << "datasize=" << datasize << endl;
	for ( long i=0; i<datasize; i++ ) set(i, (*pb)[i]);	
	delete pb;
	
	statistics();

	return 0;
}

/**
@brief 	Bins by an integer size and returns a new image.
@param 	bk			3-value vector of integer bin factors.
@return Bimage* 	new image.

	An image is binned by an integer size.

**/
Bimage*		Bimage::bin_copy(Vector3<long> bk)
{
	if ( bk[0] < 1 ) bk[0] = 1;
	if ( bk[1] < 1 ) bk[1] = 1;
	if ( bk[2] < 1 ) bk[2] = 1;
	if ( x < bk[0] ) bk[0] = x;
	if ( y < bk[1] ) bk[1] = y;
	if ( z < bk[2] ) bk[2] = z;
	
	if ( bk[0]*bk[1]*bk[2] < 2 ) return copy();
	
	Vector3<long>	nusize;
	nusize[0] = (x + bk[0] - 1)/bk[0];
	nusize[1] = (y + bk[1] - 1)/bk[1];
	nusize[2] = (z + bk[2] - 1)/bk[2];
	
	Bimage*				pb = copy_header();
	
	pb->data_type(Float);
	pb->size(nusize);

	for ( long i=0; i<pb->images(); i++ ) {
		pb->image[i].origin(pb->image[i].origin()/bk);
		pb->image[i].sampling(image[i].sampling()*bk);
	}
	
	pb->data_alloc();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Binning:" << endl;
		cout << "Bin size:                       " << bk << endl;
		cout << "New image size:                 " << nusize << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Binning" << endl << endl;

#ifdef HAVE_GCD
	dispatch_apply(pb->data_size()/c, dispatch_get_global_queue(0, 0), ^(size_t i){
		bin(i*c, bk, pb);
	});
#else
#pragma omp parallel for
	for ( long i=0; i<pb->data_size(); i+=c )
		bin(i, bk, pb);
#endif

	if ( verbose & VERB_DEBUG )
		cout << "Binning: Sums done" << endl;
	
	pb->statistics();
	
//	pb->information();
//	delete pb;
	
	return pb;
}

/**
@brief 	Bins by an integer size making sure the origin falls on a binned voxel, and returns a new image.
@param 	bin			integer bin factor.
@return Bimage* 	new image.

	An image is binned by an integer size, square in 2D and cubic in 3D.

**/
Bimage*		Bimage::bin_around_origin(int bin)
{
	if ( bin%2 == 0 ) bin++;	// bin must be odd
	
	Bimage*			pt = copy();
	
	pt->filter_average(bin);
	
	long			i, j, xx, yy, zz, xo, yo, zo, xj, yj, zj, nn;
	Vector3<long>	nusize(x/bin + 1, y/bin + 1, z/bin + 1);
	Vector3<double>	shift;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Binning:                        " << bin << endl;
		cout << "New binned size:                " << nusize << endl << endl;
	}
	
	Bimage*			pb = new Bimage(Float, compoundtype, nusize, n);
	pb->sampling(image->sampling()*bin);

	for ( i=nn=0; nn<n; ++nn ) {
		pb->image[nn].origin(image[nn].origin()/bin);
		shift = image[nn].origin() - pb->image[nn].origin() * bin;
		if ( verbose )
			cout << nn << tab << shift << endl;
		zo = shift[2];
		for ( zz=0; zz<pb->sizeZ(); ++zz, zo+=bin ) {
			zj = zo;
			if ( zj < 0 ) zj = 0;
			if ( zj >= z ) zj = z - 1;
			yo = shift[1];
			for ( yy=0; yy<pb->sizeY(); ++yy, yo+=bin ) {
				yj = yo;
				if ( yj < 0 ) yj = 0;
				if ( yj >= y ) yj = y - 1;
				xo = shift[0];
				for ( xx=0; xx<pb->sizeX(); ++xx, xo+=bin, ++i ) {
					xj = xo;
					if ( xj < 0 ) xj = 0;
					if ( xj >= x ) xj = x - 1;
					j = index(0, xj, yj, zj, nn);
					pb->set(i, (*pt)[j]);
				}
			}
		}
	}
	
	delete pt;
	
	return pb;
}


/**
@brief 	Bins by an integer size, selecting the kernel median.
@param 	binning 	integer bin factor.
@return int 		0.

	An image is binned by an integer size, square in 2D and cubic in 3D.

**/
int 		Bimage::median_bin(int binning)
{
	if ( binning < 2 ) return -1;
	
	if ( compound_type() > TSimple ) {
		error_show("Bimage::median_bin", __FILE__, __LINE__);
		cerr << "Error: Interpolation of compound data sets not supported!" << endl << endl;
		return -1;
	}
	
	// Calculate the new size and set up
	long				i, j, k;
	long				nn, ix, iy, iz, xx, yy, zz, dimensions(0);
	Vector3<long>		nusize(1,1,1);
	Vector3<long>		shift;
	Vector3<long>		blocksize(1,1,1);
	Vector3<double>		nusam(image->sampling());

	if ( x > 1 ) {
		blocksize[0] = 2*binning - 1;
		nusize[0] = (x + 1)/binning;
		if ( 2*(nusize[0]/2) == nusize[0] ) nusize[0]--; // Size must be odd
		shift[0] = nusize[0] - x/binning;
		nusam[0] *= binning;
		dimensions++;
	}
	if ( y > 1 ) {
		blocksize[1] = 2*binning - 1;
		nusize[1] = (y + 1)/binning;
		if ( 2*(nusize[1]/2) == nusize[1] ) nusize[1]--; // Size must be odd
		shift[1] = nusize[1] - y/binning;
		nusam[1] *= binning;
		dimensions++;
	}
	if ( z > 1 ) {
		blocksize[2] = 2*binning - 1;
		nusize[2] = (z + 1)/binning;
		if ( 2*(nusize[2]/2) == nusize[2] ) nusize[2]--; // Size must be odd
		shift[2] = nusize[2] - z/binning;
		nusam[2] *= binning;
		dimensions++;
	}
	
	long				blocktotalsize = blocksize[0]*blocksize[1]*blocksize[2];
	vector<double>		block(blocktotalsize);

	long				datasize = n*nusize[0]*nusize[1]*nusize[2]*c;
	double*				nudata = new double[datasize];
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Binning using a median filter:" << endl;
		cout << "Block size:                     " << blocksize << endl;
		cout << "New image size:                 " << nusize << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Binning using a median filter" << endl << endl;

	// Do the median binning interpolation
	long				oldx, oldy, oldz, xblock, yblock, zblock, nblock, imedian;
	long				zlo, zhi, ylo, yhi, xlo, xhi;
	for ( i=nn=0; nn<n; nn++ ) {
		for ( iz=0; iz<nusize[2]; iz++ ) {
			oldz = iz*binning - shift[2]; 		// Origin of the block
			zlo = oldz;
			if ( zlo < 0 ) zlo = 0;
			zhi = oldz + blocksize[2];
			if ( zhi > z ) zhi = z;
			for ( iy=0; iy<nusize[1]; iy++ ) {
				oldy = iy*binning - shift[1];
				ylo = oldy;
				if ( ylo < 0 ) ylo = 0;
				yhi = oldy + blocksize[1];
				if ( yhi > y ) yhi = y;
				for ( ix=0; ix<nusize[0]; ix++, i++ ) {
					oldx = ix*binning - shift[0];
					xlo = oldx;
					if ( xlo < 0 ) xlo = 0;
					xhi = oldx + blocksize[0];
					if ( xhi > x ) xhi = x;
					nblock = 0;
					for ( zz=zlo; zz<zhi; zz++ ) {
						zblock = zz - oldz;
						for ( yy=ylo; yy<yhi; yy++ ) {
							yblock = yy - oldy;
							for ( xx=xlo; xx<xhi; xx++ ) {
								xblock = xx - oldx;
								k = ((nn*z + zz)*y + yy)*x + xx;
								j = (zblock*blocksize[1] + yblock)*blocksize[0] + xblock;
								block[j] = (*this)[k];
								nblock++;
							}
						}
					}
					imedian = nblock/2;
					partition(block, nblock, imedian);
					nudata[i] = block[imedian];
				}
			}
		}
		image[nn].origin(image[nn].origin() * (1.0/binning));
		if ( verbose & VERB_STATS )
			cout << "Image " << nn+1 << " origin:                " << image[nn].origin() << endl;
	}

	if ( verbose & VERB_DEBUG )
		cout << "Median binning interpolation done" << endl;
	
	// Set new image parameters
	size(nusize);
	page_size(nusize);
	sampling(nusam);
	
	// Transfer the data to a new block and free the old
	data_alloc();
	for ( i=0; i<datasize; i++ ) set(i, nudata[i]);
	
	delete[] nudata;
	
	statistics();
	
	return 0;
}

