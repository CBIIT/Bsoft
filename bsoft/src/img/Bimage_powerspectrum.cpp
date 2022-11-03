/**
@file	Bimage_powerspectrum.cpp
@brief	Functions for calculating and using power spectra in electron micrographs
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20220113
**/

#include "Bimage.h"
#include "simplex.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			set_tile_size(Vector3<long>& tile_size, Vector3<long> img_size)
{
	// Minimum and maximum tile sizes allowed
	long			min_size(64);
	long			max_size(4096);
	
	tile_size = tile_size.max(min_size);
	tile_size = tile_size.min(img_size);
	tile_size = tile_size.min(max_size);

	if ( tile_size[0] != tile_size[1] ) {
		if ( tile_size[0] > tile_size[1] ) tile_size[0] = tile_size[1];
		else tile_size[1] = tile_size[0];
		if ( img_size[2] > 1 ) tile_size[2] = tile_size[0];
	}
	
	return 0;
}
	

/**
@brief 	Calculates a power spectrum.
@param 	flags			1=norm, 2=avg, 4=shift, 8=log.
@return int				error code.

	All the sub-images are Fourier transformed.
	The flags variable controls options base don which bits are set:
		1	normalize image before transformation
		2	average all power spectra
		4	shift the origin to the center
		8	calculate the logarithm of the power spectrum
		16	do edge smoothing to get rid of the cross artifact

**/
int			Bimage::power_spectrum(int flags)
{
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::power_spectrum: flags=" << flags << endl;
		cout << "DEBUG Bimage::power_spectrum: tiles=" << n << endl;
	}
	
	color_to_simple();
	change_type(Float);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating a power spectrum:" << endl;
		cout << "Flags:                          " << flags << endl;
		cout << "Tiles:                          " << n << endl << endl;
	}
	
	Bimage*		pe = NULL;
	
	if ( flags & 16 )
		pe = extract_edge_difference();
	else {
		Vector3<double>		start(1,1,0);
		Vector3<long> 		rect(size()-start*2);
		edge(0, rect, start, 1, FILL_AVERAGE);
	}

	if ( fft() )
		return error_show("Bimage::power_spectrum", __FILE__, __LINE__);

	if ( flags & 1 ) zero_fourier_origin();

	if ( flags & 16 ) {		// Very slow! 3D does not seem to work
		if ( verbose & VERB_PROCESS )
			cout << "Decomposing for edge smoothing" << endl;
			
		pe->fft();
		pe->zero_fourier_origin();

		long		i, nn, xx, yy, zz;
		double		cx = TWOPI/(double)x;
		double		cy = TWOPI/(double)y;
		double		cz = TWOPI/(double)z;
		
		vector<double>	dx(x);
		vector<double>	dy(y);
		vector<double>	dz(z,0);
		
		for ( xx=0; xx<x; ++xx ) dx[xx] = cos(cx*xx);
		for ( yy=0; yy<y; ++yy ) dy[yy] = cos(cy*yy);
		if ( z > 1 ) for ( zz=0; zz<z; ++zz ) dz[zz] = cos(cz*zz);
	
		for ( nn=0; nn<n; ++nn ) {
//			cout << n << endl;
			for ( i=nn*x*y*z, zz=0; zz<y; ++zz ) {
 				for ( yy=0; yy<y; ++yy ) {
					for ( xx=0; xx<x; ++xx, ++i ) {
   						if ( z == 1 )
							set(i, complex(i) - pe->complex(i) * (0.5/(2.0-dx[xx]-dy[yy])));
						else
							set(i, complex(i) - pe->complex(i) * (0.5/(3.0-dx[xx]-dy[yy]-dz[zz])));
					}
				}
			}
		}
	
		delete pe;
	}
	
	complex_to_intensities();
	
	fourier_type(NoTransform);

	if ( flags & 2 ) average_images();
//	cout << "number of images = " << images() << endl;
	
	statistics();
	
	origin(0,0,0);
	
	show_maximum((complex(1)).power());
	
	if ( flags & 4 ) center_wrap();
	
	if ( flags & 8 ) logarithm();			// Get logarithm of intensities
	
	return 0;
}

/**
@brief 	Prepares a tiled power spectrum from an image for determining CTF parameters.
@param	img_num			sub-image to transform.
@param 	tile_size		tile size (if (0,0,0) don't tile).
@param 	flags			1=norm, 2=avg, 4=shift, 8=log.
@return Bimage*			power spectrum.

	A large single image (a micrograph) is converted to a number of tiles
	packed into a multi-image structure.
	All the sub-images are Fourier transformed and the power spectra calculated.
	The flag indicates if the images are normalized, averaged, shifted and
	the logarithm calculated

**/
Bimage*		Bimage::powerspectrum_tiled(long img_num, Vector3<long> tile_size, int flags)
{
	// If the image is square it is assumed to be a CCD image or extracted to be square
	// and no padding around the edges is necessary
	int 			pad(0);
	long			max_size(4096);
	
	// If the image is not square, it is assumed to be a micrograph with severe edge
	// features that will affect the power spectrum
	if ( x != y && x > max_size && y > max_size) {
		pad = x/10;
		if ( pad < y/10 ) pad = y/10;
	}
	
	Vector3<long> 	start, ext_size, step_size;

	set_tile_size(tile_size, size());

	ext_size[0] = tile_size[0]*((x - 2*pad)/tile_size[0]);
	ext_size[1] = tile_size[1]*((y - 2*pad)/tile_size[1]);
	ext_size[2] = tile_size[2]*((z - 2*pad)/tile_size[2]);
	ext_size = ext_size.max(1);

	start[0] = (x - ext_size[0])/2;
	start[1] = (y - ext_size[1])/2;
	start[2] = (z - ext_size[2])/2;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::powerspectrum_tiled: img_num=" << img_num << endl;
		cout << "DEBUG Bimage::powerspectrum_tiled: ext_size=" << ext_size << endl;
		cout << "DEBUG Bimage::powerspectrum_tiled: start=" << start << endl;
		cout << "DEBUG Bimage::powerspectrum_tiled: tile_size=" << tile_size << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Tile size:                      " << tile_size << endl << endl;
	
	// If a single image assume a micrograph that needs to be tiled
	Bimage* 		pex = extract_tiles(img_num, start, ext_size, tile_size, step_size, 0);
	
	pex->power_spectrum(flags);
	
	return pex;
}

Bimage*		Bimage::powerspectrum_tiled_exact(long img_num, Vector3<long> tile_size, int flags)
{
	Vector3<long> 	start, ext_size(size()), step_size;

	if ( verbose & VERB_PROCESS )
		cout << "Tile size:                      " << tile_size << endl << endl;
	
	Bimage* 		pex = extract_tiles(img_num, start, ext_size, tile_size, step_size, 0);
	
	pex->power_spectrum(flags);
	
	return pex;
}


/**
@brief 	Prepares a tiled powerspectrum from a tilted image for determining CTF parameters.
@param	img_num			sub-image to transform.
@param 	tile_size		tile size (if (0,0,0) don't tile).
@param 	tilt_axis		tilt axis angle (in radians).
@param 	tilt_offset		offset perpendicular to tilt axis (in pixels).
@param 	flags			1=norm, 2=avg, 4=shift, 8=log.
@return Bimage*			power spectrum.

	A large single image (a micrograph) is converted to a number of tiles
	along the tilt axis packed into a multi-image structure.
	All the sub-images are Fourier transformed and the power spectra calculated.
	The flag indicates if the images are normalized, averaged, shifted and
	the logarithm calculated

**/
Bimage*		Bimage::powerspectrum_tilt_axis(long img_num, Vector3<long> tile_size,
					double tilt_axis, double tilt_offset, int flags)
{
	long			i, xx, yy, ntiles(0);
	double			tan_a(tan(tilt_axis));
	double			inv_tan_a(1.0/tan_a);
	double			ratang(atan2(y, x));
	Vector3<long> 	origin(size()/2);
	vector<Vector3<long>>	tile;

	set_tile_size(tile_size, size());
	
	if ( fabs(tilt_axis) < ratang || fabs(tilt_axis) > M_PI - ratang ) {
		origin[1] += (long)tilt_offset;
		ntiles = x/tile_size[0];
		tile.resize(ntiles);
		for ( i=0, xx = origin[0] - (ntiles * tile_size[0])/2; i<ntiles; i++, xx += tile_size[0] ) {
			tile[i][0] = xx;
			tile[i][1] = (long) (origin[1] + (xx - origin[0])*tan_a);
		}
	} else {
		origin[0] += (long)tilt_offset;
		ntiles = y/tile_size[1];
		tile.resize(ntiles);
		for ( i=0, yy = origin[1] - (ntiles * tile_size[1])/2; i<ntiles; i++, yy += tile_size[1] ) {
			tile[i][0] = (long) (origin[0] + (yy - origin[1])*inv_tan_a);
			tile[i][1] = yy;
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::powerspectrum_tilted: img_num=" << img_num << endl;
		cout << "DEBUG Bimage::powerspectrum_tilted: origin=" << origin << endl;
		cout << "DEBUG Bimage::powerspectrum_tilted: tile_size=" << tile_size << endl;
		cout << "DEBUG Bimage::powerspectrum_tilted: ntiles=" << ntiles << endl;
//		for ( i=0; i<ntiles; i++ )
//			cout << "DEBUG Bimage::powerspectrum_tilted: " << tile[i][0] << " " << tile[i][1] << endl;
	}
	
	// If a single image assume a micrograph that needs to be tiled
	Bimage* 			pex = extract_tiles(img_num, tile, tile_size);

	Vector3<long> 		translate;
	Vector3<long> 		nusize(1,1,1);
	if ( pex->sizeX() != pex->sizeY() ) {
		nusize[0] = pex->sizeX();
		if ( pex->sizeY() > 1 && nusize[0] > pex->sizeY() ) nusize[0] = pex->sizeY();
		if ( pex->sizeZ() > 1 && nusize[0] > pex->sizeZ() ) nusize[0] = pex->sizeZ();
		if ( pex->sizeY() > 1 ) nusize[1] = nusize[0];
		if ( pex->sizeZ() > 1 ) nusize[2] = nusize[0];
		pex->resize(nusize, translate, FILL_USER, 0);
	}

//	write_img("pex.pif", pex);

	pex->power_spectrum(flags);

//	write_img("pex.pif", pex);
	
	return pex;
}

Bimage*		Bimage::defocus_scale(long nn, double df, double df2, double iCL2, int fill_type)
{
	Bimage*			ps = new Bimage(Float, TSimple, size(), 1);
	
	long			i, j, xx, yy;
	double			interval(1.0/(image[nn].sampling()[0]*x));
	double			s1, s2, f1(iCL2*df), f2(2*iCL2*df2), len, fill(image[nn].background());
	vector<double>	vs(x,0);
	
	for ( size_t i=0; i<vs.size(); ++i ) {
		s2 = interval*i;
		s2 *= s2;
		s1 = f1 - sqrt(f1*f1 + s2*s2 - f2*s2);
		vs[i] = sqrt(s1)/interval;
	}
	
	Vector3<double>	dis, coor, ori(size()/2);
	
	for ( i=0, yy=0; yy<y; ++yy ) {
		dis[1] = double(yy) - ori[1];
		for ( xx=0; xx<x; ++xx, ++i ) {
			dis[0] = double(xx) - ori[0];
			len = dis.length();
			coor = ori;
			if ( len ) {
				j = long(len);
				f2 = len - j;
				f1 = 1 - f2;
				coor += dis * (f1*vs[j] + f2*vs[j+1])/len;
			}
			ps->set(i, interpolate(coor, nn, fill));
		}
	}
	
	return	ps;
}

/**
@brief 	Prepares a tiled powerspectrum from a tilted image for determining CTF parameters.
@param	img_num			sub-image to transform.
@param 	tile_size		tile size (if (0,0,0) don't tile).
@param 	tilt_axis		tilt axis angle (in radians).
@param 	tilt_angle		tilt angle (in radians).
@param	defocus			average defocus to adjust for change in focus.
@param	iCL2			inverse of product of spherical aberration and wavelenght squared.
@param 	flags			1=norm, 2=avg, 4=shift, 8=log.
@return Bimage*			power spectrum.

	A large single image (a micrograph) is converted to a number of tiles
	packed into a multi-image structure.
	All the sub-images are Fourier transformed and the power spectra calculated.
	The power spectra are scaled based on the tilt and average defocus of the image.
	The flag indicates if the images are normalized, averaged, shifted and
	the logarithm calculated

**/
Bimage*		Bimage::powerspectrum_tilted(long img_num, Vector3<long> tile_size,
				double tilt_axis, double tilt_angle, double defocus, double iCL2, int flags)
{
	if ( image->origin()[0] <= 0 || image->origin()[1] <= 0 )
		image->origin(size()/2);

	set_tile_size(tile_size, size());
	
	long				nn;
	Vector3<long> 		start, region(size()), step_size(tile_size/2);
	Vector3<long> 		edge_rect(tile_size/1.5), edge_start((tile_size-edge_rect)/2);
	edge_rect = edge_rect.max(1);
	
//	cout << tile_size << tab << edge_rect << tab << edge_start << endl;

	// The coordinates are the lower-left corners of the tiles
	vector<Vector3<long>>	coor = tile_coordinates(start, region, tile_size, step_size, 0);

	Bimage*				pex = extract_tiles(img_num, coor, tile_size);

	pex->edge(1, edge_rect, edge_start, edge_start[0]/2.0, FILL_AVERAGE, 0);
	
	pex->power_spectrum(5);
	pex->calculate_background();

	double				ddef, ca(cos(tilt_axis)), sa(sin(tilt_axis)), tt(tan(tilt_angle));
	Bimage*				pt;
	Vector3<double>		vec, origin(pex->size()/2), scale(1,1,1), translate;
	Matrix3				mat(1);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::powerspectrum_tilted: origin=" << image->origin() << endl;
		cout << "DEBUG Bimage::powerspectrum_tilted: ca=" << ca << " sa=" << sa << " tt=" << tt << endl;
	}
	
	// The scale is calculated so that it compensates for the tile defocus
	// The defocus difference is subtracted from the reference defocus
	// This is equivalent to negating the pixel vector or negating the tilt axis
//	if ( verbose )
//		cout << "Tile\tx\ty\t∆f" << endl;
	for ( nn=0; nn<pex->images(); nn++ ) {
		vec = (origin - image->origin() + coor[nn])*image->sampling();
//		ddef = (vec[0]*sa - vec[1]*ca)*tt;
		ddef = (vec[1]*ca - vec[0]*sa)*tt;
//		if ( iCL2 > 0.002 ) {
//			scale[0] = scale[1] = sqrt(defocus/(defocus + ddef));
			scale[0] = scale[1] = sqrt(1 + ddef/defocus);
			pt = pex->transform(nn, pex->size(), scale, origin, translate, mat, FILL_BACKGROUND);
//		} else {
//			pt = pex->defocus_scale(nn, defocus, defocus + ddef, iCL2, FILL_BACKGROUND);
//		}
		pex->replace(nn, pt);
//		if ( verbose )
//			cout << nn << tab << vec[0] << tab << vec[1] << tab << ddef << endl;
		delete pt;
	}
	
	if ( flags & 2 )
		pex->average_images();
	
	pex->origin(pex->size()/2);
	
	return pex;
}

/**
@brief 	Prepares a tiled powerspectrum from a tilted image for determining CTF parameters.
@param 	tile_size		tile size (if (0,0,0) don't tile).
@param 	tilt_axis		tilt axis angle (in radians).
@param 	tilt_angle		tilt angle (in radians).
@param 	tilt_offset		offset perpendicular to tilt axis (in pixels).
@param	defocus			average defocus to adjust for change in focus.
@param	iCL2			inverse of product of spherical aberration and wavelenght squared.
@param 	flags			1=norm, 2=avg, 4=shift, 8=log, 16=add.
@return Bimage*			power spectrum.

	A large single image (a micrograph) is converted to a number of tiles
	packed into a multi-image structure.
	All the sub-images are Fourier transformed and the power spectra calculated.
	The power spectra are scaled based on the tilt and average defocus of the image.
	The flag indicates if the images are normalized, averaged, shifted and
	the logarithm calculated

**/
Bimage*		Bimage::powerspectrum_tiled_and_tilted(Vector3<long> tile_size,
				double tilt_axis, double tilt_angle, double tilt_offset, 
				double defocus, double iCL2, int flags)
{
	long			nn, nimg(n);
	Bimage*			ps1 = NULL;
	if ( flags & 16 ) nimg = 1;

	set_tile_size(tile_size, size());
	
	Bimage*			ps = new Bimage(Float, TSimple, tile_size, nimg);	
	ps->sampling(sampling(0));
	
	for ( nn = 0; nn<n; nn++ ) {
		if ( fabs(tilt_angle) > 1e-6 ) {
			if ( defocus > 0 )
				ps1 = powerspectrum_tilted(nn, tile_size, tilt_axis, tilt_angle, defocus, iCL2, flags);
			else
				ps1 = powerspectrum_tilt_axis(nn, tile_size, tilt_axis, tilt_offset, flags);
		} else
			ps1 = powerspectrum_tiled(nn, tile_size, flags);
		if ( flags & 16 ) {
			ps->add(ps1);
		} else {
			ps->replace(nn, ps1);
		}
		delete ps1;
	}
	
	if ( flags & 16 ) ps->multiply(1.0/nn);
	
	ps->sampling(sampling(0));
	ps->statistics();
	ps->file_name(file_name());
	ps->origin(ps->size()/2);

	return ps;
}

double		isotropy_R(Bsimplex& simp)
{
	long			i;
	double			R(0), df;
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=0; i<simp.points(); i++ ) {
		df = f[i] - (simp.parameter(0) + simp.parameter(1)*cos(2*(simp.parameter(2)+x[i])));
		R += df*df;
	}
	
	R = sqrt(R/i);
			
	return R;
}

/**
@brief 	Calculates a measure of anisotropy in a poer spectrum.
@param	n				sub-image number.
@param 	&lores			low resolution limit.
@param 	&hires			high resolution limit
@return vector<double>		3-vlaue vector: power average and deviation and maximum power angle.
	The power between the indicated resolution shells are averaged for
	each angle and fitted to an equation for anisotropy:
		P = Pavg + Pdev*cos(2(a-phi))
	where phi is the direction of maximum power.
**/
vector<double>	Bimage::powerspectrum_isotropy(long n, double& lores, double& hires)
{
	if ( lores <= 0 || lores > real_size()[0] ) lores = real_size()[0];
	if ( lores < hires ) swap(lores, hires);
	if ( hires < 2*sampling(0)[0] ) hires = 2*sampling(0)[0];
	
	long			i, j, xx, yy;
	long			kmin(real_size()[0]/lores);
	long			kmax(real_size()[0]/hires);
	long			na(TWOPI*kmax), h(x/2);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Power spectrum anisotropy:" << endl;
		cout << "Origin:                         " << image[n].origin() << endl;
		cout << "Resolution:                     " << hires << " - " << lores << " A" << endl;
		cout << "Pixel radii:                    " << kmin << " - " << kmax << endl;
		cout << "Angles:                         " << na << endl;
	}
	
	double			dx, dy, r, a(0), f, v, y_avg(0), y_var(0), xy_scale(x/y);
	vector<double>	vx(na,0), vy(na,0);
		
	for ( i=n*size().volume(), yy=0; yy<y; ++yy ) {
		dy = (yy - image[n].origin()[1])*xy_scale;
		if ( dy > h ) dy -= x;
		for ( xx=0; xx<x; ++xx, ++i ) {
			dx = xx - image[n].origin()[0];
			if ( dx > h ) dx -= x;
			r = sqrt(dx*dx + dy*dy);
			if ( r >= kmin && r <= kmax ) {
				v = (*this)[i];
				a = (na/TWOPI)*atan2(dy, dx);
				if ( a < 0 ) a += na;
				j = long(a);
				if ( j > na-1 ) cerr << "j too large!" << endl;
				f = a - j;
				vx[j] += 1-f;
				vy[j] += (1-f)*v;
				j++;
				vx[j] += f;
				vy[j] += f*v;
			}
		}
	}

	for ( i=0; i<na; ++i) {
//		cout << i << tab << x[i] << tab << y[i] << endl;
		vy[i] /= vx[i];
		vx[i] = TWOPI*i*1.0/na;
		y_avg += vy[i];
		y_var += vy[i]*vy[i];
	}
	y_avg /= na;
	y_var /= na;
	y_var -= y_avg*y_avg;
	
//	cout << "avg=" << y_avg << " var=" << y_var << endl;

	Bsimplex		simp(1, 3, 0, na, vx, vy);
	
	simp.parameter(0, y_avg);
	simp.parameter(1, 2*sqrt(y_var));
	simp.parameter(2, 0);
	simp.limits(0, 0.7*y_avg, 1.5*y_avg);
	simp.limits(1, y_var, 5*sqrt(y_var));
	simp.limits(2, -M_PI, M_PI);
	
	double			R = simp.run(10000, 1e-5, isotropy_R);
	
	double			aniso = simp.parameter(1)/simp.parameter(0);
	vector<double>	fit(3,0);
	fit[0] = simp.parameter(0);
	fit[1] = simp.parameter(1);
	fit[2] = simp.parameter(2);
	if ( fit[2] < -M_PI_2 ) fit[2] += M_PI;
	if ( fit[2] > M_PI_2 ) fit[2] -= M_PI;
	
//	cout << "limit_lo=" << simp.limit_low(1) << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Fitted function:                " << fit[0] << " + " << fit[1] <<
			" * cos(2*(a + " << fit[2]*180/M_PI << "))" << endl;
		cout << "Anisotropy:                     " << aniso << endl;
		cout << "R:                              " << R << endl << endl;
	}

	if ( verbose & VERB_FULL ) {
		cout << "Angle\tPower\tFit" << endl;
		for ( i=0; i<na; ++i)
			cout << vx[i]*180.0/M_PI << tab << vy[i] << tab << fit[0] + fit[1]*cos(2*(vx[i]+fit[2])) << endl;
		cout << endl;
	}
	
	return fit;
}

/*
int			Bimage::powerspectrum_edge_smoothed(int flags)
{
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::power_spectrum: flags=" << flags << endl;
		cout << "DEBUG Bimage::power_spectrum: tiles=" << n << endl;
	}
	
	color_to_simple();
	change_type(Float);
	
	if ( verbose & VERB_FULL ) {
		cout << "Calculating a power spectrum:" << endl;
		cout << "Flags:                           " << flags << endl;
		cout << "Tiles:                           " << n << endl;
	}

	Bimage*		pe = extract_edge_difference();

	if ( fft() )
		return error_show("Bimage::power_spectrum", __FILE__, __LINE__);

	if ( flags & 1 ) zero_fourier_origin();

	pe->fft();
	pe->zero_fourier_origin();

	long		i, nn, xx, yy;
	double		dx, dy, v;
	double		cx = TWOPI/(double)x;
	double		cy = TWOPI/(double)y;
	
	for ( nn=0; nn<n; ++nn ) {
		for ( i=nn*x*y*z, yy=0; yy<y; ++yy ) {
   			dy = cos(cy*yy);
			for ( xx=0; xx<x; ++xx, ++i ) {
   				dx = cos(cx*xx);
				set(i, complex(i) - pe->complex(i) * (0.5/(2.0-dx-dy)));
			}
		}
	}
	
	delete pe;

	complex_to_intensities();
	
	fourier_type(NoTransform);

	if ( flags & 2 ) average_images();
//	cout << "number of images = " << images() << endl;
	
	statistics();
	
	origin(0,0,0);
	
	show_maximum((complex(1)).power());
	
	if ( flags & 4 ) center_wrap();
	
	if ( flags & 8 ) logarithm();			// Get logarithm of intensities
	
	return 0;
}
*/
