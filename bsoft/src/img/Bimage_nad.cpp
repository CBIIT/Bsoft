/**
@file	Bimage_nad.cpp
@brief	Denoising by nonlinear anisotropic diffusion: Coherence and edge enhancing diffusion.
@author Achilleas Frangakis
@author Bernard Heymann
@date	Created: 20020803
@date	Modified: 20161210 (BH)
**/

#include "Bimage.h"
#include "Matrix.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
 Calculates a structure tensor for a 2D image. 
*/
Bimage*		Bimage::structure_tensor_2D(long nn)
{
	if ( z > 1 ) {
		cerr << "Error in Bimage::structure_tensor_2D: The image must be 2D!" << endl;
		return NULL;
	}
	
	long			i, k, xx, yy;
	double			dx, dy;
	
	Bimage*			pd = new Bimage(Float, 3, x, y, 1, 1);
	pd->image->origin(image->origin());
	pd->sampling(sampling(0));
	
	for ( i=nn*image_size(), k=yy=0; yy<y; yy++ ) {
		for ( xx=0; xx<x; xx++, i++ ) {
			dx = dy = -avg;
			if ( xx > 0 ) dx = -(*this)[i-1];
			if ( xx < x-1 ) dx += (*this)[i+1];
			else dx += avg;
			if ( yy > 0 ) dy = -(*this)[i-x];
			if ( yy < y-1 ) dy += (*this)[i+x];
			else dy += avg;
			dx /= 2;
			dy /= 2;
			pd->set(k++, dx*dx);
			pd->set(k++, dx*dy);
			pd->set(k++, dy*dy);
//			cout << xx << tab << yy << tab << dx << tab << dy << endl;
		}
	}
	
	return pd;
}

/*
 Calculates a structure tensor for a 3D volume. 
*/
Bimage*		Bimage::structure_tensor(long zs, long zw)
{
	long			i, k, xx, yy, zz, zf;
	long			nxy(x*y);
	double			dx, dy, dz;
	
	// Padding for the top and bottom slices
	if ( zs > 0 ) zw++;
	if ( zs+zw < z ) zw++;
	
	Bimage*			pd = new Bimage(Float, 6, x, y, zw, 1);
	pd->image->origin(image->origin());
	pd->sampling(sampling(0));
	
	// Start point for the first slice
	if ( zs > 0 ) {
		pd->image->origin(0,0,1);
		zs--;
	}
	
	for ( zz=zs, zf=zs+zw, k=0, i=nxy*zs; zz<zf; zz++ ) {
		for ( yy=0; yy<y; yy++ ) {
			for ( xx=0; xx<x; xx++, i++ ) {
				dx = dy = dz = -avg;
				if ( xx > 0 ) dx = -(*this)[i-1];
				if ( xx < x-1 ) dx += (*this)[i+1];
				else dx += avg;
				if ( yy > 0 ) dy = -(*this)[i-x];
				if ( yy < y-1 ) dy += (*this)[i+x];
				else dy += avg;
				if ( zz > 0 ) dz = -(*this)[i-nxy];
				if ( zz < z-1 ) dz += (*this)[i+nxy];
				else dz += avg;
				dx /= 2;
				dy /= 2;
				dz /= 2;
				pd->set(k++, dx*dx);
				pd->set(k++, dx*dy);
				pd->set(k++, dx*dz);
				pd->set(k++, dy*dy);
				pd->set(k++, dy*dz);
				pd->set(k++, dz*dz);
//				cout << xx << tab << yy << tab << zz << tab << dx << tab << dy << tab << dz << endl;
			}
		}
	}
	
	return pd;
}

/*
Calculates the diffusion tensor of EED or CED using the structure tensor.
*/
int			Bimage::diffusion_tensor_2D(double lambda, double C, double alpha)
{
	if ( z > 1 ) {
		cerr << "Error in Bimage::diffusion_tensor_2D: The image must be 2D!" << endl;
		return -1;
	}
	
	long			i, j, imgsize(x*y);
	double			lam[2] = {1,1};		/* eigenvalues of diffusion tensor */
//	double*			d;					/* vector with unordered eigenvalues */
	Matrix			a(2,2);				/* real tensor matrix */
	double			beta(1 - alpha);	/* time saver */
	double			thexp(16), grd, gol;
	
	for ( i=0; i<imgsize; i++ ) {
		j = i*c;
		a[0][0] = (*this)[j++];
		a[0][1] = a[1][0] = (*this)[j++];
		a[1][1] = (*this)[j++];
//		cout << i << tab << a[0] << tab << a[1] << tab << a[2] << endl;
		
		grd = a[0][0] + a[1][1];
		vector<double>	d = a.jacobi_rotation();
		a.eigen_sort(d);
		
		if ( C ) {
			// for plane-like structures, calculate the second eigenvalue
			lam[0] = lam[1] = alpha;
			if ( d[1] )
				lam[1] += beta * exp (- C / (d[1]*d[1]));
		} else {
			if ( grd > 0.0) {
				gol = sqrt(grd)/lambda;
				lam[0] = 1.0 - exp (-3.31488 / pow(gol, thexp));
			} else {
				lam[0] = lam[1] = 1;
			}
		}
		
//		delete[] d;
		
		/* principal axis backtransformation */
		j = i*c;
		set(j++, lam[0]*a[0][0]*a[0][0] + lam[1]*a[0][1]*a[0][1]);
		set(j++, lam[0]*a[0][0]*a[1][0] + lam[1]*a[0][1]*a[1][1]);
		set(j++, lam[0]*a[1][0]*a[1][0] + lam[1]*a[1][1]*a[1][1]);		
		//		cout << i;
		//		for ( int k=0; k<3; k++ ) cout << tab << (*this)[j++];
		//		cout << endl;
	}
	
	return 0;
}

/*
Calculates the diffusion tensor of EED or CED using the structure tensor.
Principal axis backtransformation of a symmetric 3x3 matrix. 
A = U * diag(lam[0], lam[1], lam[2]) * U_transpose with U = (v1 | v2 | v3)     
*/
int			Bimage::diffusion_tensor(double lambda, double C, double alpha)
{
	long			i, j, imgsize(x*y*z);
	double			lam[3] = {1,1,1};	/* eigenvalues of diffusion tensor */
//	double*			d;					/* vector with unordered eigenvalues */
	Matrix			a(3,3);				/* real tensor matrix */
	double			beta(1 - alpha);	/* time saver */
	double			thexp(16), grd, gol, denom;
	
	for ( i=0; i<imgsize; i++ ) {
		j = i*c;
		a[0][0] = (*this)[j++];
		a[0][1] = a[1][0] = (*this)[j++];
		a[0][2] = a[2][0] = (*this)[j++];
		a[1][1] = (*this)[j++];
		a[1][2] = a[2][1] = (*this)[j++];
		a[2][2] = (*this)[j++];
		//		cout << i << tab << a[0] << tab << a[1] << tab << a[2] << tab << a[4] << tab << a[5] << tab << a[8] << endl;
		
		grd = a[0][0] + a[1][1] + a[2][2];
//		d = a.jacobi_rotation();
		vector<double>	d = dsyevq3(a);
		a.eigen_sort(d);
//		cout << d[0] << tab << d[1] << tab << d[2] << endl;
//		cout << a << endl;
		
		if ( C ) {
			// for plane-like structures, calculate the second eigenvalue
			lam[0] = lam[1] = lam[2] = alpha;
			denom = d[0] - d[1];
			if ( denom && denom > d[2] && ( denom > d[1] - d[2] ) )
				lam[1] += beta * exp (- C / (denom*denom));
			/* calculate third eigenvalue */
			denom = d[0] - d[2];
			if ( denom ) lam[2] += beta * exp (- C / (denom*denom));
		} else {
			if ( grd > 0.0) {
				gol = sqrt(grd)/lambda;
				lam[0] = lam[1] = 1.0 - exp (-3.31488 / pow(gol, thexp));
			} else {
				lam[0] = lam[1] = 1;
			}
		}
		
//		delete[] d;
		
		/* principal axis backtransformation */
		j = i*c;
		set(j++, lam[0]*a[0][0]*a[0][0] + lam[1]*a[0][1]*a[0][1] + lam[2]*a[0][2]*a[0][2]);
		set(j++, lam[0]*a[0][0]*a[1][0] + lam[1]*a[0][1]*a[1][1] + lam[2]*a[0][2]*a[1][2]);
		set(j++, lam[0]*a[0][0]*a[2][0] + lam[1]*a[0][1]*a[2][1] + lam[2]*a[0][2]*a[2][2]);
		set(j++, lam[0]*a[1][0]*a[1][0] + lam[1]*a[1][1]*a[1][1] + lam[2]*a[1][2]*a[1][2]);
		set(j++, lam[0]*a[1][0]*a[2][0] + lam[1]*a[1][1]*a[2][1] + lam[2]*a[1][2]*a[2][2]);
		set(j++, lam[0]*a[2][0]*a[2][0] + lam[1]*a[2][1]*a[2][1] + lam[2]*a[2][2]*a[2][2]);
		
		//		cout << i;
		//		for ( int k=0; k<6; k++ ) cout << tab << (*this)[j++];
		//		cout << endl;
	}
	
	return 0;
}

int			Bimage::diffuse_2D(Bimage* pdn, Bimage* pd, double ht, long nn)
{
	if ( z > 1 ) {
		cerr << "Error in Bimage::diffuse_2D: The image must be 2D!" << endl;
		return -1;
	}
	
	long			i, j, k, xx, yy, ii;
	long			iw, xw, yw, ix, iy, xf1, yf1;
	double			w[9], v;
	double			r(ht/4);
	
	if ( verbose & VERB_FULL )
		cout << "Calculating the new image chunk" << endl;
	
	xf1 = pd->x - 1;
	yf1 = pd->y - 1;
	
	for ( i=yy=0, j=nn*image_size(); yy<y; yy++ ) {
		for ( xx=0; xx<x; xx++, i++, j++ ) {
			for ( iw=0; iw<9; iw++ ) w[iw] = 0;
				
			k = pd->c*i;
			
//				w[12] = w[14] = 2*(*pd)[k]   - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
				w[3] = w[5] = 2*(*pd)[k]   - fabs((*pd)[k+1]);
//				w[10] = w[16] = 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
				w[1] = w[7] = 2*(*pd)[k+2] - fabs((*pd)[k+1]);
//				w[4]  = w[22] = 2*(*pd)[k+5] - fabs((*pd)[k+4]) - fabs((*pd)[k+2]);
				
//				w[9]  = w[17] =  (*pd)[k+1] + fabs((*pd)[k+1]);
				w[0]  = w[8] =  (*pd)[k+1] + fabs((*pd)[k+1]);
//				w[15] = w[11] = -(*pd)[k+1] + fabs((*pd)[k+1]);
				w[6] = w[2] = -(*pd)[k+1] + fabs((*pd)[k+1]);
				
				if ( xx > 0 ) {
					k = pd->c*(i - 1);
//					w[12] += 2*(*pd)[k] - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
					w[3] += 2*(*pd)[k] - fabs((*pd)[k+1]);
				}
				if ( xx < xf1 ) {
					k = pd->c*(i + 1);
//					w[14] += 2*(*pd)[k] - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
					w[5] += 2*(*pd)[k] - fabs((*pd)[k+1]);
				}
				
				if ( yy > 0 ) {
					k = pd->c*(i - pd->x);
//					w[10] += 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
					w[1] += 2*(*pd)[k+2] - fabs((*pd)[k+1]);
				}
				if ( yy < yf1 ) {
					k = pd->c*(i + pd->x);
//					w[16] += 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
					w[7] += 2*(*pd)[k+2] - fabs((*pd)[k+1]);
				}
				
				if ( xx > 0 && yy > 0 ) {
					k = pd->c*(i - pd->x - 1) + 1;
//					w[9] += (*pd)[k] + fabs((*pd)[k]);
					w[0] += (*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && yy > 0 ) {
					k = pd->c*(i - pd->x + 1) + 1;
//					w[11] += -(*pd)[k] + fabs((*pd)[k]);
					w[2] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx > 0 && yy < yf1 ) {
					k = pd->c*(i + pd->x - 1) + 1;
//					w[15] += -(*pd)[k] + fabs((*pd)[k]);
					w[6] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && yy < yf1 ) {
					k = pd->c*(i + pd->x + 1) + 1;
//					w[17] += (*pd)[k] + fabs((*pd)[k]);
					w[8] += (*pd)[k] + fabs((*pd)[k]);
				}
				
			/* modify weights to prevent flux across boundaries */
			if ( xx==0 ) w[6] = w[3] = w[0] = 0.0; 
			if ( xx==x-1 ) w[2] = w[5] = w[8] = 0.0; 
			if ( yy==0 ) w[0] = w[1] = w[2] = 0.0; 
			if ( yy==y-1 ) w[8] = w[7] = w[6] = 0.0; 
				
			//				cout << xx << tab << yy <<;
			//				for ( iw=0; iw<9; iw++ ) cout << tab << w[iw]*r;
			//				cout << endl;
				
			/* evolution */
			v = (*this)[i];
			for ( yw=0; yw<2; yw++ ) {
				iy = yy + yw - 1;
				if ( iy >= 0 && iy < y ) for ( xw=0; xw<2; xw++ ) {
					ix = xx + xw - 1;
					if ( ix >= 0 && ix < x ) {
						iw = 2*yw + xw;
						if ( w[iw] ) {
							ii = iy*x + ix;
							v += w[iw]*r*((*this)[ii] - (*this)[i]);
						}
					}
				}
			}
			pdn->set(j, v);
		}
	}
	
	return 0;
}

int			Bimage::diffuse(Bimage* pdn, Bimage* pd, double ht, long zs)
{
	long			i, k, xx, yy, zz, zf, ii, ik;
	long			iw, xw, yw, zw, ix, iy, iz, zs1, xf1, yf1, zf1;
	long			nxy(x*y);
	double			w[27], v;
	double			r(ht/4);
	
	if ( verbose & VERB_FULL )
		cout << "Calculating the new image chunk" << endl;
	
	xf1 = pd->x - 1;
	yf1 = pd->y - 1;
	zs1 = zs;
	ik = 0;
	if ( zs > 0 ) {	// The diffusion tensor reference starts one slice in
		ik = nxy;
		zs1 = zs - 1;
	}
	zf = zs1 + pd->z;
	zf1 = zf - 1;
	if ( zf < z ) zf -= 1;
	if ( zf > z ) zf = z;
	
	if ( verbose & VERB_DEBUG )
		cout << "pd->sizeZ()=" << pd->z << " zs=" << zs << " zs1=" << zs1 << " zf=" << zf << " zf1=" << zf1 << endl;
	
	for ( zz=zs, i=nxy*zs; zz<zf; zz++ ) {
		for ( yy=0; yy<y; yy++ ) {
			for ( xx=0; xx<x; xx++, i++, ik++ ) {
				for ( iw=0; iw<27; iw++ ) w[iw] = 0;
				
				k = pd->c*ik;
				
				w[12] = w[14] = 2*(*pd)[k]   - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
				w[10] = w[16] = 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
				w[4]  = w[22] = 2*(*pd)[k+5] - fabs((*pd)[k+4]) - fabs((*pd)[k+2]);
				
				w[9]  = w[17] =  (*pd)[k+1] + fabs((*pd)[k+1]);
				w[15] = w[11] = -(*pd)[k+1] + fabs((*pd)[k+1]);
				
				w[3]  = w[23] =  (*pd)[k+2] + fabs((*pd)[k+2]);
				w[5]  = w[21] = -(*pd)[k+2] + fabs((*pd)[k+2]);
				
				w[1]  = w[25] =  (*pd)[k+4] + fabs((*pd)[k+4]);
				w[7]  = w[19] = -(*pd)[k+4] + fabs((*pd)[k+4]);
				
				if ( xx > 0 ) {
					k = pd->c*(ik - 1);
					w[12] += 2*(*pd)[k] - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
				}
				if ( xx < xf1 ) {
					k = pd->c*(ik + 1);
					w[14] += 2*(*pd)[k] - fabs((*pd)[k+1]) - fabs((*pd)[k+2]);
				}
				
				if ( yy > 0 ) {
					k = pd->c*(ik - pd->x);
					w[10] += 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
				}
				if ( yy < yf1 ) {
					k = pd->c*(ik + pd->x);
					w[16] += 2*(*pd)[k+3] - fabs((*pd)[k+1]) - fabs((*pd)[k+4]);
				}
				
				if ( zz > zs1 ) {
					k = pd->c*(ik - nxy);
					w[4] += 2*(*pd)[k+5] - fabs((*pd)[k+4]) - fabs((*pd)[k+2]);
				}
				if ( zz < zf1 ) {
					k = pd->c*(ik + nxy);
					w[22] += 2*(*pd)[k+5] - fabs((*pd)[k+4]) - fabs((*pd)[k+2]);
				}
				
				if ( xx > 0 && yy > 0 ) {
					k = pd->c*(ik - pd->x - 1) + 1;
					w[9] += (*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && yy > 0 ) {
					k = pd->c*(ik - pd->x + 1) + 1;
					w[11] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx > 0 && yy < yf1 ) {
					k = pd->c*(ik + pd->x - 1) + 1;
					w[15] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && yy < yf1 ) {
					k = pd->c*(ik + pd->x + 1) + 1;
					w[17] += (*pd)[k] + fabs((*pd)[k]);
				}
				
				if ( xx > 0 && zz > zs1 ) {
					k = pd->c*(ik - nxy - 1) + 2;
					w[3] += (*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && zz > zs1 ) {
					k = pd->c*(ik - nxy + 1) + 2;
					w[5] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx > 0 && zz < zf1 ) {
					k = pd->c*(ik + nxy - 1) + 2;
					w[21] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( xx < xf1 && zz < zf1 ) {
					k = pd->c*(ik + nxy + 1) + 2;
					w[23] += (*pd)[k] + fabs((*pd)[k]);
				}
				
				if ( yy > 0 && zz > zs1 ) {
					k = pd->c*(ik - pd->x - nxy) + 4;
					w[1] += (*pd)[k] + fabs((*pd)[k]);
				}
				if ( yy < yf1 && zz > zs1 ) {
					k = pd->c*(ik + pd->x - nxy) + 4;
					w[7] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( yy > 0 && zz < zf1 ) {
					k = pd->c*(ik - pd->x + nxy) + 4;
					w[19] += -(*pd)[k] + fabs((*pd)[k]);
				}
				if ( yy < yf1 && zz < zf1 ) {
					k = pd->c*(ik + pd->x + nxy) + 4;
					w[25] += (*pd)[k] + fabs((*pd)[k]);
				}
				
				/* modify weights to prevent flux across boundaries */
				if ( xx==0 ) w[15] = w[12] = w[9] = w[3] = w[21] = 0.0; 
				if ( xx==x-1 ) w[11] = w[14] = w[17] = w[5] =w[23] = 0.0; 
				if ( yy==0 ) w[9] = w[10] = w[11] = w[1] = w[19] = 0.0; 
				if ( yy==y-1 ) w[17] = w[16] = w[15] = w[7] = w[25] = 0.0; 
				if ( zz==0 ) w[7] = w[4] = w[1] = w[5] = w[3] = 0.0; 
				if ( zz==z-1 )  w[23] = w[22] = w[21] = w[1] =w[25] = 0.0;
				
				//				cout << xx << tab << yy << tab << zz;
				//				for ( iw=0; iw<27; iw++ ) cout << tab << w[iw]*r;
				//				cout << endl;
				
				/* evolution */
				v = (*this)[i];
				for ( zw=0; zw<3; zw++ ) {
					iz = zz + zw - 1;
					if ( iz >= 0 && iz < z ) for ( yw=0; yw<3; yw++ ) {
						iy = yy + yw - 1;
						if ( iy >= 0 && iy < y ) for ( xw=0; xw<3; xw++ ) {
							ix = xx + xw - 1;
							if ( ix >= 0 && ix < x ) {
								iw = 9*zw + 3*yw + xw;
								if ( w[iw] ) {
									ii = (iz*y + iy)*x + ix;
									v += w[iw]*r*((*this)[ii] - (*this)[i]);
								}
							}
						}
					}
				}
				pdn->set(i, v);
			}
		}
	}
	
	return 0;
}

int			Bimage::nad_chunk_2D(Bimage* pdn, double lambda, double C,
				double alpha, double ht, long nn)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::nad_chunk_2D: nn = " << nn << endl;
	
	Bimage*			pd = structure_tensor_2D(nn);
	
	pd->diffusion_tensor_2D(lambda, C, alpha);
	
	diffuse_2D(pdn, pd, ht, nn);
	
	delete pd;
	
	return 0;
}

int			Bimage::nad_chunk(Bimage* pdn, double lambda, double C,
				double alpha, double ht, long zs, long zw)
{
	if ( zs+zw > z ) zw = z - zs;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bimage::nad_chunk: Chunk start = " << zs << " and size = " << zw << endl;
	
	Bimage*			pd = structure_tensor(zs, zw);
	
	pd->diffusion_tensor(lambda, C, alpha);
	
	diffuse(pdn, pd, ht, zs);
	
	delete pd;
	
	return 0;
}

/**
@brief 	Denoises a 2D image by non-linear anisotropic diffusion.
@param 	ht			time step size, 0 < ht <= 0.25.
@param 	lambda		lamda parameter for EED.
@param 	C			coherence parameter for CED.
@param 	alpha		.
@return Bimage*		new image.

	The diffusion tensor is calculated with the aim of enhancing edges (EED)
	or planes (CED).

**/
Bimage*		Bimage::nad_2D(double ht, double lambda, double C, double alpha)
{
	change_type(Float);
	
	Bimage*			pdn = copy();
	
	for ( long nn=0; nn<n; ++nn )
		nad_chunk_2D(pdn, lambda, C, alpha, ht, nn);

	return pdn;
}


/**
@brief 	Denoises a 3D density map by non-linear anisotropic diffusion.
@param 	ht			time step size, 0 < ht <= 0.25.
@param 	zw			slab size for piece-wise denoising.
@param 	lambda		lamda parameter for EED.
@param 	C			coherence parameter for CED.
@param 	alpha		.
@return Bimage*		new image.

	The diffusion tensor is calculated with the aim of enhancing edges (EED)
	or planes (CED).

**/
Bimage*		Bimage::nad(double ht, long zw, double lambda, double C, double alpha)
{
	change_type(Float);
	
	Bimage*			pdn = copy();
	
#ifdef HAVE_GCD
	dispatch_apply((z - 1)/zw + 1, dispatch_get_global_queue(0, 0), ^(size_t zz){
		nad_chunk(pdn, lambda, C, alpha, ht, zz*zw, zw);
	});
#else
#pragma omp parallel for
	for ( long zz=0; zz<z; zz+=zw )
		nad_chunk(pdn, lambda, C, alpha, ht, zz, zw);
#endif

	return pdn;
}


