/**
@file	Bimage_reconstruct.cpp
@brief	2D and 3D reconstruction from single particle images
@author	Bernard Heymann
@date	Created: 20010403
@date	Modified: 20200127
**/

#include "Bimage.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Interpolates a 2D image for packing ito a 3D reciprocal space volume.  
@param 	cv			complex value from 2D transform.
@param 	m			location in 3D relative to origin.
@param 	part_weight	weight to assign to value (usually 1).
@param 	interp_type	interpolation type (0=nearest neighbor, 1=weighted nearest neigbor, 2=trilinear).
@return int			0, <0 on error.

**/
int			Bimage::fspace_2D_interpolate(Complex<float> cv, Vector3<double> m,
				double part_weight, int interp_type)
{
	if ( !next ) next = new Bimage(Float, TSimple, size(), n);
	if ( !next->next ) next->next = new Bimage(Float, TSimple, size(), n);
	if ( !next->next->next ) next->next->next = new Bimage(Float, TSimple, size(), n);

	long 				xx, yy, zz, j, jx, jy, jz;
	double				w(part_weight);
	Vector3<long>		coor;
	Vector3<double>		dist;
	
	float*				power = (float *) next->data_pointer();
	float* 				weight = (float *) next->next->data_pointer();
	float* 				weight2 = (float *) next->next->next->data_pointer();
	
	if ( interp_type < 2 ) {			// Nearest neighbour
		coor = Vector3<long>((long) floor(m[0] + 0.5),
			(long) floor(m[1] + 0.5),
			(long) floor(m[2] + 0.5));
		if ( interp_type == 1 ) {		// Weighted nearest neighbour
			dist = m - coor;
			w *= 1 - dist.length();
		}
		if ( w > 0 ) {
			j = index_wrap(coor);
			add(j, cv * w);
			power[j] += w*cv.power();
			weight[j] += w;
			weight2[j] += w*w;
		}
	} else if ( interp_type == 2 ) {	// Trilinear interpolation
		coor = Vector3<long>((long) floor(m[0]),
			(long) floor(m[1]),
			(long) floor(m[2]));
		dist = m - coor;
		for ( zz=0; zz<2; zz++ ) {
			jz = coor[2] + zz;
			if ( jz < 0 ) jz += z;
			else if ( jz >= z ) jz -= z;
			dist[2] = 1.0 - dist[2];
			for ( yy=0; yy<2; yy++ ) {
				jy = coor[1] + yy;
				if ( jy < 0 ) jy += y;
				else if ( jy >= y ) jy -= y;
				dist[1] = 1.0 - dist[1];
				for ( xx=0; xx<2; xx++ ) {
					jx = coor[0] + xx;
					if ( jx < 0 ) jx += x;
					else if ( jx >= x ) jx -= x;
					dist[0] = 1.0 - dist[0];
					w = dist.volume() * part_weight;
					j = index(jx, jy, jz);
					add(j, cv * w);
					power[j] += w*cv.power();
					weight[j] += w;
					weight2[j] += w*w;
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Packs a 2D Fourier transform into a 3D reciprocal space volume.  
@param 	*p			2D particle image transform.
@param 	mat			rotation matrix.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit (infinite if 0).
@param 	scale		scale of reconstruction and particle magnification.
@param 	part_weight	weight of particle (usually 1).
@param 	interp_type	interpolation type (0=nearest neighbor, 1=weighted nearest neigbor, 2=trilinear).
@return int			0, <0 on error.

	The rotation matrix is used to determine the plane in reciprocal space
	to which the 2D transform data is added. The map is assumed to be cubic
	and the 2D transform square. The orientation parameters must be written
	into the image structure. 

**/
int			Bimage::fspace_pack_2D(Bimage* p, Matrix3 mat, double hi_res, 
				double lo_res, Vector3<double> scale, double part_weight, int interp_type)
{
	double			min_rad_sq = (lo_res)? sampling(0)[0]/lo_res: 0;
	double			max_rad_sq = (hi_res)? sampling(0)[0]/hi_res: 0.5;
	if ( max_rad_sq > 0.5 ) max_rad_sq = 0.5;
	min_rad_sq *= x;
	max_rad_sq *= x;
	max_rad_sq += 1;		// Add one pixel width additional interpolation to be removed later
	min_rad_sq *= min_rad_sq;
	max_rad_sq *= max_rad_sq;
	
	long 			i, xx, yy;
	long			hx = (p->x - 1)/2, hy = (p->y - 1)/2;
	double			d, w;
	Vector3<double>	m, iv;
	Vector3<double> vscale(x/(scale[0]*p->x), y/(scale[1]*p->y), 1);
	
	mat[2] *= z*1.0/x;	 
	
	if ( verbose & VERB_FULL )
		cout << "Packing an image into reciprocal space up to " << hi_res << " A resolution" << endl;

	if ( verbose & VERB_DEBUG )
		cout << mat << endl;
		
	for ( yy=i=0; yy<p->y; yy++ ) {
		iv[1] = yy;
		if ( yy > hy ) iv[1] -= (double)p->y;
		iv[1] *= vscale[1];
		for ( xx=0; xx<p->x; xx++, i++ ) {
			iv[0] = xx;
			if ( xx > hx ) iv[0] -= (double)p->x;
			iv[0] *= vscale[0];
			d = iv.length2();
			if ( d >= min_rad_sq && d <= max_rad_sq ) {
				d = max_rad_sq - d;
				m = mat * iv;
				w = part_weight;
				if ( d < 1 ) w *= sqrt(d);
				fspace_2D_interpolate(p->complex(i), m,
					w, interp_type);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Packs a 2D Fourier transform into a 3D reciprocal space volume.  
@param 	*p			2D particle image transform.
@param 	asu_view	view of asymmetric unit.
@param 	*sym		point group symmetry.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit (infinite if 0).
@param 	scale		scale of reconstruction and particle magnification.
@param 	part_weight	weight of particle (usually 1).
@param 	interp_type	interpolation type (0=nearest neighbor, 1=weighted nearest neigbor, 2=trilinear).
@return int			0, <0 on error.

	The rotation matrix is used to determine the plane in reciprocal space
	to which the 2D transform data is added. The map is assumed to be cubic
	and the 2D transform square. The orientation parameters must be written
	into the image structure.

**/
int			Bimage::fspace_pack_2D(Bimage* p, View asu_view, Bsymmetry& sym, double hi_res, 
				double lo_res, Vector3<double> scale, double part_weight, int interp_type)
{
	View*			v;
	View*			view = symmetry_get_all_views(sym, asu_view);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG Bimage::fspace_pack_2D: sym=" << sym.label() << endl;
		cout << "DEBUG Bimage::fspace_pack_2D: view=" << asu_view << endl;
		cout << "DEBUG Bimage::fspace_pack_2D: part_weight=" << part_weight << endl;
	}
	
	for ( v=view; v; v=v->next )
		fspace_pack_2D(p, v->matrix(), hi_res, lo_res, scale, part_weight, interp_type);
	
	kill_list((char *) view, sizeof(View));
	
	return 0;
}

/**
@brief 	Packs a 2D Fourier transform into a central section of a 3D reciprocal space volume.  
@param 	*p			2D particle image transform.
@param 	ft_size		Fourier transform size.
@param 	scale		scale of reconstruction and particle magnification.
@param 	hi_res		high resolution limit.
@param 	lo_res		low resolution limit (infinite if 0).
@param 	matr		central section orientation matrix.
@param 	mat			image orientation matrix.
@return int			0, <0 on error.

	The rotation matrix is used to determine the plane in reciprocal space
	to which the 2D transform data is added in reference to the rotation
	matrix of the central section. The map is assumed to be cubic
	and the 2D transform square. The orientation parameters must be written
	into the image structure. 

**/
long		Bimage::fspace_pack_2D_into_central_section(Bimage* p,
				long ft_size, double scale, double hi_res, double lo_res, 
				Matrix3 matr, Matrix3 mat)
{
	if ( !p->data_pointer() ) return 0;
	
	check_resolution(hi_res);
	
	if ( verbose & VERB_FULL )
		cout << "Packing an image into reciprocal space up to " << hi_res << " A resolution" << endl << endl;
	
	long 			i, j, xx, yy, kx, ky;
	long 			ix, iy, ixx, iyy, nover(0);
	double 			w, dx, dy, d2;
	Vector3<double>	m, d(0,0,1), iv;
	
//	Vector3<double> vscale(x/(scale*ft_size), y/(scale*ft_size), 1/scale);

	float*			fom = (float *) next->data_pointer();
	
	if ( verbose & VERB_FULL )
		cout << mat << endl;
	
	double			maxrad = p->image->sampling()[0]/hi_res;
	double			maxrad_sq = maxrad*maxrad;	// Maximum radius squared on map scale
	double			ymin = floor(-maxrad*p->y);
	double			ymax = -ymin;
	double			xmin = floor(-maxrad*p->x);
	double			xmax = -xmin;
	
//	cout << mat << endl << matr << endl;
	
	// Vector to calculate z-shift
///	Vector3<double>	vn(mat[2][0]/mat[2][2], mat[2][1]/mat[2][2], 0);
//	Vector3<double>	vn(mat[2][2]/mat[2][0], mat[2][2]/mat[2][1], 0);
//	cout << vn << tab << vscale << endl;
	
///	matr = matr.transpose();
//	matr = vscale * matr;
	
	matr = mat * matr.transpose();
//	matr = matr.transpose() * mat;
	Vector3<double>	vn(1/matr[0][0], 1/matr[1][1], 0);
//	cout << "vector = " << vn << endl;

	for ( iv[1]=ymin; iv[1]<=ymax; iv[1]+=1 ) {
		if ( iv[1] >= 0 ) yy = (long)iv[1];
		else yy = (long)(iv[1] + p->y);
		dy = iv[1]/p->y;
		for ( iv[0]=xmin; iv[0]<=xmax; iv[0]+=1 ) {
			if ( iv[0] >= 0 ) xx = (long)iv[0];
			else xx = (long)(iv[0] + p->x);
			dx = iv[0]/p->x;
			d2 = dx*dx + dy*dy;
			if ( d2 <= maxrad_sq ) {
				m = matr * iv;
///				m = mat * iv;
///				m[2] = m.scalar(vn);
///				m = matr * m;
//				d[2] = fabs(m[2]);
				if ( m[2] <= 0.5 ) {
					m *= vn;
					ix = (long) floor(m[0]);
					iy = (long) floor(m[1]);
					d[0] = m[0] - ix;
					d[1] = m[1] - iy;
					i = yy*p->x + xx;
					for ( ky=0; ky<2; ky++ ) {
						iyy = iy + ky;
						if ( iyy < 0 ) iyy += y;
						d[1] = 1 - d[1];
						for ( kx=0; kx<2; kx++ ) {
							ixx = ix + kx;
							if ( ixx < 0 ) ixx += x;
							d[0] = 1 - d[0];
							j = iyy*x + ixx;
							w = d[0]*d[1];
							fom[j] += w;
							add(j, p->complex(i) * w);
						}
					}
					nover++;
				}
			}
		}
	}

	return nover;
}


/**
@brief 	Packs a 3D Fourier transform into a 3D reciprocal space volume.
@param 	*p			3D particle image transform.
@param 	hi_res		high resolution limit.
@param 	threshold	threshold to exclude low intensities.
@return int			0, <0 on error.

	The image is added up to the high resolution limit and excluding
	low intensities as defined by the threshold.

**/
int			Bimage::fspace_pack_3D(Bimage* p, double hi_res, double threshold)
{
	if ( !next ) next = new Bimage(Float, TSimple, size(), n);
	if ( !next->next ) next->next = new Bimage(Float, TSimple, size(), n);
	if ( !next->next->next ) next->next->next = new Bimage(Float, TSimple, size(), n);

	check_resolution(hi_res);
	
	long				i, xx, yy, zz;
	double				rx, ry, rz, r2, maxrad, maxrad2, pwr, w;
	Vector3<long>		h((size()-1)/2);
	Vector3<double>		freq_scale(1.0/real_size());
	Complex<float>		cv;

	float*				power = (float *) next->data_pointer();
	float* 				weight = (float *) next->next->data_pointer();
	float* 				weight2 = (float *) next->next->next->data_pointer();

	double				rad_scale(real_size()[0]), rad_scale2(rad_scale*rad_scale);
	
	if ( hi_res > 0.1 )
		maxrad = rad_scale/hi_res;
	else
		maxrad = size().max()/2.0;
	maxrad2 = maxrad*maxrad;
	
	for ( i=zz=0; zz<z; zz++ ) {
		rz = zz;
		if ( zz > h[2] ) rz -= z;
		rz *= freq_scale[2];
		rz *= rz;
		for ( yy=0; yy<y; yy++ ) {
			ry = yy;
			if ( yy > h[1] ) ry -= y;
			ry *= freq_scale[1];
			ry *= ry;
			for ( xx=0; xx<x; xx++, i++ ) {
				rx = xx;
				if ( xx > h[0] ) rx -= x;
				rx *= freq_scale[0];
				rx *= rx;
				cv = p->complex(i);
				pwr = cv.power();
				if ( pwr >= threshold ) {
					r2 = rad_scale2*(rx + ry + rz);
					if ( r2 <= maxrad2 ) {
						w = 1;
						add(i, cv);
						power[i] += pwr;
						weight[i] += w;
						weight2[i] += w*w;
					}
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Adds all components to a reconstruction.
@param	*p			frequency space reconstruction to add.
@return long		0.

	The FOM block contains the sum of powers, the next image conatins the
	weight sum, and the next image FOM block contains the sum of the
	weight squared.

**/
long		Bimage::fspace_reconstruction_add(Bimage* p)
{
//	cout << "adding" << endl;
	
	add(p);
	
	if ( next && p->next ) {
		next->add(p->next);
		
		if ( next->next && p->next->next ) {
			next->next->add(p->next->next);
			
			if ( next->next->next && p->next->next->next ) next->next->next->add(p->next->next->next);
		}
	}
	
	return 0;
}

/**
@brief 	Weighs a reconstruction.
@return long		coverage.

	The FOM block contains the sum of powers, the next image conatins the
	weight sum, and the next image FOM block contains the sum of the
	weight squared.

**/
long		Bimage::fspace_reconstruction_weigh()
{
	long		 		i, cov(0), ds(x*y*z);
	Complex<double>		cv;

	float*				power = (float *) next->data_pointer();
	float* 				weight = (float *) next->next->data_pointer();
	float* 				weight2 = (float *) next->next->next->data_pointer();
	
//	weight[0] = 0;
	
	for ( i=0; i<ds; i++ ) {
		if ( weight[i] > SMALLFLOAT ) {
			weight2[i] = weight[i] - weight2[i]/weight[i];
			set(i, complex(i) / weight[i]);
			if ( weight2[i] > 1 ) power[i] /= weight2[i];
			else if ( weight[i] > 1 ) power[i] /= weight[i];
			if ( power[i] < 0 ) power[i] = 0;
			cov++;
		} else {
			set(i, cv); 
			power[i] = 0;
			weight[i] = weight2[i] = SMALLFLOAT;
		}
	}
	
	power[0] = 1; 			// The zero frequency contains no information
		
	return cov;
}

/**
@brief 	Calculates Fourier shell statistics.
@param 	resolution		high resolution limit.
@param	sampling_ratio	frequency space sampling (default 1 pixel/sample).
@return int				0.

	The average FOM and number of FOM values in each resolution shell is determined.

**/
int 		Bimage::fspace_reconstruction_stats(double resolution, double sampling_ratio)
{
//	cout << "Compound type: " << compoundtype << endl;
//	cout << "Fourier type:  " << fouriertype << endl;
	
	long			maxrad = fspace_maximum_radius(resolution, sampling_ratio);
	
	Bimage*			prad = fspace_radial_power(resolution, sampling_ratio);
	Bimage*			pfom = prad->next = next->radial(0, maxrad, sampling_ratio, 1);
	Bimage*			pw = prad->next->next = next->next->radial(0, maxrad, sampling_ratio, 1);
	Bimage*			pw2 = prad->next->next->next = next->next->next->radial(0, maxrad, sampling_ratio, 1);

	long			i;
	double			rad_scale(real_size()[0]/sampling_ratio);
	double			s, noise, signal, w, w2;
	double			res, fsc, fscp(1), snr, snrp(10), sc(0), res_fsc(0), res_snr(0);
	double			fsc_cut(0.3), snr_cut(0.5);

	if ( verbose ) {
		cout << "Radius\ts(1/A)\tRes(A)\t"
			<< "Weight\tPower\tSignal\tNoise\tSNR\tFSC" << endl;
		for ( i=1; i<maxrad; i++ ) {
			s = i/rad_scale;
			res = rad_scale/i;
			fsc = snr = signal = noise = 0;
			w = (*pw)[i];
			if ( w > 1e-3 ) {
				w2 = (*pw2)[i];
//				w2 = w - (*pw2)[i]/w;
//				if ( w < 1 ) w = 1;
//				if ( w2 < 1 ) w2 = 1;
				noise = ((*pfom)[i] - (*prad)[i])/w2;
				signal = (w*(*prad)[i] - (*pfom)[i])/w2;
//				if ( noise < 0 )
//					noise = (*prad)[i] - signal;
//				if ( noise < 0 )
//					cout << i << tab << w2 << tab << (*prad)[i] << tab << (*pfom)[i] << tab << (*prad)[i]/(*pfom)[i] << endl;
//				if ( signal > 0 )signal -= sqrt(noise*signal);
//				else signal = 0;
//				signal = ((w-2-2/w)*(*prad)[i] - (1+2/w)*(*pfom)[i])/w2;
//				if ( signal < 0 ) signal = 0;
				if ( noise > 0 ) snr = signal/noise;
//				fsc = (w-(*pfom)[i]/(*prad)[i])/(w-1);
				fsc = snr/(snr+1);
				if ( fsc < 0 ) fsc = 0;
			}
			if ( fscp >= fsc_cut && fsc < fsc_cut ) {
				sc = s - (fsc_cut - fsc)/((fscp - fsc)*rad_scale);
				res_fsc = 1/sc;
			}
			if ( snrp >= snr_cut && snr < snr_cut ) {
				sc = s - (snr_cut - snr)/((snrp - snr)*rad_scale);
				res_snr = 1/sc;
			}
			fscp = fsc;
			snrp = snr;
			cout << i << tab << setprecision(4) << s << tab
				<< setprecision(2) << res << tab
				<< setprecision(2) << w << tab
				<< setprecision(6) << (*prad)[i] << tab
				<< setprecision(6) << signal << tab
				<< setprecision(6) << noise << tab
				<< setprecision(4) << snr << tab << fsc << endl;
		}
		cout << endl;
		cout << "Resolution estimates:" << endl;
		cout << "SNR(0.5):                       " << res_snr << " A" << endl;
		cout << "FSC(0.3):                       " << res_fsc << " A" << endl << endl;
	}
	
	delete prad;
	
	return 0;
}

/**
@brief 	Calculates the SNR map.
@return long		coverage.

	The FOM block contains the sum of powers, the next image conatins the
	weight sum, and the next image FOM block contains the sum of the
	weight squared.

**/
long		Bimage::fspace_reconstruction_snr()
{
	long		 		i, cov(0), ds(x*y*z);
	double				signal, noise;

	float*				power = (float *) next->data_pointer();
	float* 				weight = (float *) next->next->data_pointer();
	
	for ( i=0; i<ds; i++ ) {
		signal = complex(i).power();
		noise = power[i] - signal;
		if ( noise < signal ) noise = signal;
//		if ( signal < SMALLFLOAT && noise < SMALLFLOAT ) 
	//		cout << i << tab << signal << tab << noise << endl;
		if ( weight[i] > 1 && signal > SMALLFLOAT && noise > SMALLFLOAT ) {
			power[i] = (weight[i]*signal - power[i])/noise;
			cov++;
		} else
			power[i] = 0;
	}
	
	power[0] = 1; 			// The zero frequency contains no information
		
	return cov;
}
