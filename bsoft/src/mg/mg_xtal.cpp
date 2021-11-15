/**
@file	mg_xtal.cpp
@brief	Functions to process crystallographic data
@author	Bernard Heymann
@date	Created: 20061110
@date	Modified: 20181030
**/

#include "mg_processing.h"
#include "mg_xtal.h"
#include "matrix_linear.h"
#include "Matrix.h"
#include "linked_list.h"
#include "Complex.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the unit cell vectors for a 2D crystal.
@param 	*mg			micrograph.
@return int			0, <0 on error.

	Finds the unit cell vectors u and v by solving the equation:
		x = uh + vk
	where x is the location of the reflection or structure factor,
	and h and k are the associated Miller indices.

**/
int			mg_unitcell_vectors(Bmicrograph* mg)
{
	long			i;
	vector<double>	bx(2,0), by(2,0);
	Matrix			a(2,2);
	Bstrucfac*		sf;
	
	for ( i = 0, sf = mg->sf; sf; sf = sf->next, i++ ) {
		a[0][0] += sf->index[0]*sf->index[0];
		a[0][1] += sf->index[0]*sf->index[1];
		a[1][1] += sf->index[1]*sf->index[1];
		bx[0] += sf->index[0]*sf->loc[0];
		bx[1] += sf->index[1]*sf->loc[0];
		by[0] += sf->index[0]*sf->loc[1];
		by[1] += sf->index[1]*sf->loc[1];
	}
	a[1][0] = a[0][1];
	
	if ( i < 2 ) {
		cerr << "Error: Too few structure factors to calculate unit cell vectors!" << endl;
		return -1;
	}
	
	a.LU_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	
	mg->hvec[0] = bx[0];
	mg->hvec[1] = by[0];
	mg->kvec[0] = bx[1];
	mg->kvec[1] = by[1];
		
	return 0;
}

/**
@brief 	Generates reflections given the unit cell vectors.
@param 	*mg				micrograph.
@param	real_size		physical image size.
@param 	resolution		resolution limit.
@return long			number of reflections generated, <0 on error.

	The structure factor location is given by:
		x = uh + vk
	where u and v are the unit cell vectors,
	and h and k are the associated Miller indices.

**/
long		mg_generate_reflections(Bmicrograph* mg, Vector3<double> real_size, double resolution)
{
	if ( resolution < 2*mg->pixel_size[0] )
		resolution = 2*mg->pixel_size[0];
	
	double				hkl_lim = real_size[0]/resolution;
	if ( real_size[1] > real_size[0] ) hkl_lim = real_size[1]/resolution;
	double				hkl_lim2 = hkl_lim*hkl_lim + 0.001;
	
	if ( verbose )
		cout << "Generating reflections to " << resolution << " Å" << endl;
		
	long				n(0), h, k, l, hmin, kmin, lmin, hmax, kmax, lmax;
	double				h2, k2, l2, hlen, klen;
	Vector3<double>		loc;
	Bstrucfac*			sf = NULL;
	
	hlen = mg->hvec.length();
	hmax = (long) (hkl_lim/hlen);
	hmin = -hmax;

	klen = mg->kvec.length();
	kmax = (long) (hkl_lim/klen);					
	kmin = -kmax;
	
	lmin = lmax = 0;
	
	kill_list((char *) mg->sf, sizeof(Bstrucfac));
	mg->sf = NULL;
	
	if ( verbose ) {
		cout << "real_size=" << real_size << endl;
		cout << "hvec=" << mg->hvec << endl;
		cout << "kvec=" << mg->kvec << endl;
		cout << "hmax=" << hmax << " kmax=" << kmax << endl;
		cout << "mg->origin=" << mg->origin << endl;
	}
	
	for ( l=lmin; l<=lmax; l++ ) {
		l2 = l*l;
		for ( k=kmin; k<=kmax; k++ ) {
			k2 = k*klen;
			k2 *= k2;
			for ( h=hmin; h<=hmax; h++ ) {
				h2 = h*hlen;
				h2 *= h2;
				if ( h2 + k2 + l2 <= hkl_lim2 ) {
					loc[0] = h*mg->hvec[0] + k*mg->kvec[0] + l*mg->lvec[0];
					loc[1] = h*mg->hvec[1] + k*mg->kvec[1] + l*mg->lvec[1];
					loc[2] = h*mg->hvec[2] + k*mg->kvec[2] + l*mg->lvec[2];
					sf = (Bstrucfac *) add_item((char **) &sf, sizeof(Bstrucfac));
					if ( !mg->sf ) mg->sf = sf;
					sf->index[0] = h;
					sf->index[1] = k;
					sf->index[2] = l;
					sf->loc = loc;
					sf->fom = sf->sel = 1;
					n++;
				}
			}
		}
	}
	
	return n;
}

/**
@brief 	Masks the image using the list of reflections.
@param 	*p			complex image.
@param 	*sflist		reflection list.
@param 	radius		radius around reflection to mask.
@return int			error code.
**/
int			img_mask_reflections(Bimage* p, Bstrucfac* sflist, double radius)
{
	long			i, datasize = (long) p->size().volume();
//	Vector3<long>	mid(p->size()/2);
	Vector3<double>	loc;
	Bstrucfac*		sf = NULL;

	Bimage*			pmask = new Bimage(Float, p->compound_type(), p->size(), 1);
	Complex<double>	cv;

	if ( p->compound_type() == TSimple ) {	// Power spectrum
		for ( sf = sflist; sf; sf = sf->next )
			pmask->sphere(sf->loc, radius, 0, FILL_USER, 1);
		for ( i=0; i<datasize; i++ )
			if ( (*pmask)[i] < 1 ) p->set(i, 0);
	} else if ( p->compound_type() == TComplex ) {	// Standard Fourier transform
		for ( sf = sflist; sf; sf = sf->next ) {
			loc = vector3_set_PBC(sf->loc, p->size());
			pmask->sphere(loc, radius, 0, FILL_USER, 1);
		}
		for ( i=0; i<datasize; i++ )
			if ( (*pmask)[i] < 1 ) p->set(i, cv);;
	}
	
	delete pmask;
	
    return 0;
}


