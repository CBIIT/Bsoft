/**
@file	model_multifit.cpp
@brief	Searching for a template in a map and returning multiple hits in a model.
@author Bernard Heymann
@date	Created: 20021027
@date	Modified: 20151127
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "model_multifit.h"
#include "symmetry.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Searches a 2D/3D density map for a template.
@param 	*p			the image.
@param 	*ptemp		the template to be searched for.
@param 	*view		views.
@param 	alpha		rotation around view vector, <0 = use 2*PI (radians).
@param 	alpha_step	angular step size around view vector (radians).
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		mask for cross-correlation (ignored if NULL).
@param 	threshold	threshold value, if 0, threshold = FOMmax/2.
@return Bmodel*		list of solutions.

	The template is rotated and cross-correlated to find the best fit.
	The views must be calculated externally to allow for custom sets.

**/
Bmodel*		model_from_densities(Bimage* p, Bimage* ptemp, View* view,
				double alpha, double alpha_step, double hires, double lores, Bimage* pmask, double threshold)
{
	int				mode = 0;		// Mode: 0=global, 1=refine, 2=symmetry
	long			i, n = 0, nviews;
	double			alpha_min = 0;
	double			alpha_max = TWOPI - alpha_step/2;
	double			a;
	View*			v;
	
	if ( alpha >= 0 ) {			// Case of refinement
		mode = 1;
		alpha_min = alpha - alpha_step;
		alpha_max = alpha + alpha_step;
	} else if ( alpha_step >= 1.999*M_PI ) mode = 2;	// Symmetry
	
	if ( !p ) return NULL;
	
	p->change_type(Float);
	ptemp->change_type(Float);
	p->rescale_to_avg_std(0, 1);
	ptemp->rescale_to_avg_std(0, 1);
	p->calculate_background();
	ptemp->calculate_background();

	Vector3<double>	origin(ptemp->image->origin());
	
	if ( p->sizeZ() == 1 ) origin[2] = 0;
	
	if ( origin.sum() == 0 ) {
		origin = p->size()/2;
		p->origin(origin);
	}
	
	for ( nviews=0, v=view; v; v=v->next ) 
		for ( a=alpha_min; a<=alpha_max; a+=alpha_step ) nviews++;
		
	View*			view_arr = new View[nviews];
	Bmodel**		mod_arr = new Bmodel*[nviews];
	
	for ( n=0, v=view; v; v=v->next ) 
		for ( a=alpha_min; a<=alpha_max; a+=alpha_step, n++ ) {
			view_arr[n] = *v;
			if ( mode < 2 ) view_arr[n][3] = a;
		}
		
	if ( verbose ) {
		cout << "Finding a template in a map:" << endl;
		cout << "Number of views:                " << nviews << endl;
		cout << "Alpha step size:                " << alpha_step*180.0/M_PI << endl;
		cout << "Origin:                         " << origin << endl;
	}

	Bstring			id("1");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;

	model->mapfile(p->file_name());
	model->identifier(base(model->mapfile()));

	Bcomptype*		ct = model->add_type(id);
	ct->file_name(ptemp->file_name());
	
	if ( verbose & VERB_RESULT )
		cout << "#\tx\ty\tz\tvx\tvy\tvz\tva\tCC" << endl;

//#pragma omp parallel for private(comp)
	for ( n=0; n<nviews; n++ ) {
		mod_arr[n] = model_from_densities_for_view(p, ptemp, view_arr[n], hires, lores, pmask, threshold);
		for ( comp = mod_arr[n]->comp; comp; comp = comp->next ) {
			if ( verbose & VERB_RESULT )
				cout << n << tab << comp->location()[0] << tab << comp->location()[1] << tab << 
					comp->location()[2] << tab << comp->view() << tab << comp->FOM() << endl;
		}
	}
	
	double		fommin = 1, fommax = 0, fomavg = 0, fomstd = 0;
	for ( n=i=0; n<nviews; n++ ) {
		for ( comp = mod_arr[n]->comp; comp; comp = comp->next ) {
			i++;
			comp->identifier() = to_string(i);
			comp->type(ct);
			comp->select(1);
			if ( fommin > comp->FOM() ) fommin = comp->FOM();
			if ( fommax < comp->FOM() ) fommax = comp->FOM();
			fomavg += comp->FOM();
			fomstd += comp->FOM()*comp->FOM();
		}
		if ( model->comp ) {
			for ( comp = model->comp; comp->next; comp = comp->next ) ;
			comp->next = mod_arr[n]->comp;
		} else {
			model->comp = mod_arr[n]->comp;
		}
		mod_arr[n]->comp = NULL;
		model_kill(mod_arr[n]);
	}

	delete[] view_arr;
	delete[] mod_arr;
	
	if ( verbose ) {
		if ( i ) {
			fomavg /= i;
			fomstd = fomstd/i - fomavg*fomavg;
			if ( fomstd > 0 ) fomstd = sqrt(fomstd);
			else fomstd = 0;
		}
		cout << "FOM minimum and maximum:        " << fommin << " " << fommax << endl;
		cout << "FOM average:                    " << fomavg << endl;
		cout << "FOM standard deviation:         " << fomstd << endl << endl;
	}
	
	return model;
}

/**
@brief 	Searches a 2D/3D density map for a template using a specific view.
@param 	*p			the image.
@param 	*ptemp		the template to be searched for.
@param 	view		view.
@param 	hires		high resolution limit.
@param 	lores		low resolution limit.
@param 	*pmask		mask for cross-correlation (ignored if NULL).
@param 	threshold	threshold value, if 0, threshold = FOMmax/2.
@return Bmodel*		list of solutions.

	The template is rotated to the view and cross-correlated to find
	a set of high-scoring fits.
	The views must be calculated externally to allow for custom sets.

**/
Bmodel*		model_from_densities_for_view(Bimage* p, Bimage* ptemp, View view,
				double hires, double lores, Bimage* pmask, double threshold)
{
	Matrix3			mat = view.matrix();
//	mat = mat.transpose();
	
	Bimage*			prot = ptemp->rotate(p->size(), mat);
	
	prot->view(view);
	if ( p->sizeZ() == 1 ) prot->project('z');

	Bimage* 		pcc = p->cross_correlate(prot, hires, lores, pmask);

	Bmodel*			model = model_from_peaks(pcc, threshold, 1);
	Bcomponent*		comp = NULL;

	for ( comp = model->comp; comp; comp = comp->next ) {
		comp->view(View2<float>(view[0],view[1],view[2],view[3]));
		comp->select(1);
	}

	delete prot;
	delete pcc;
	
	return model;
}

/**
@brief 	Generates a model from peaks in a map.
@param 	*p			cross-correlation map.
@param 	threshold	threshold value, if 0, threshold = FOMmax/2.
@param 	wrap		flag to wrap around image boundaries.
@return Bmodel*		list of peaks.

	Peaks found in the map are returned in the FOM block.
	A component is generated from each peak with a FOM above the threshold.

**/
Bmodel*		model_from_peaks(Bimage* p, double threshold, int wrap)
{
	Bimage*			pfom = p->find_peaks(3);
	
	if ( pfom->maximum() <= 0 ) return NULL;

	if ( threshold <= 0 ) threshold = pfom->maximum()/2;
	
    long			i, j, x, y, z, n;
	Vector3<long>	h(p->size()/2);
	Vector3<double>	origin;
		
//	Bstring			id(p->file_name().post_rev('/').pre_rev('.'));
	Bstring			id("Peaks");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	Bstring			nutype("PEAK");
	
	model->mapfile() = p->file_name();
	model->FOM(1);

	for ( i=j=n=0; n<p->images(); n++ ) {
		origin = p->image[n].origin();
		for ( z=0; z<p->sizeZ(); z++ ) {
			for ( y=0; y<p->sizeY(); y++ ) {
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					if ( (*pfom)[i] >= threshold ) {
//						id = Bstring(++j, "%d");
//						comp = component_add(&comp, id);
//						if ( !model->comp ) model->comp = comp;
						if ( comp ) comp = comp->add(++j);
						else model->comp = comp = new Bcomponent(++j);
						comp->type(model->add_type(nutype));
						comp->location(Vector3<double>(x, y, z) - origin);
						if ( wrap ) {
							if ( comp->location()[0] > h[0] ) comp->location()[0] -= p->sizeX();
							if ( comp->location()[1] > h[1] ) comp->location()[1] -= p->sizeY();
							if ( comp->location()[2] > h[2] ) comp->location()[2] -= p->sizeZ();
						}
						comp->scale(p->sampling(0));
						comp->select(n+1);
						comp->FOM((*pfom)[i]);
					}
				}
			}
		}
	}
	
	delete pfom;
	
	if ( verbose )
		cout << "Components generated:           " << j << endl << endl;
	
	return model;
}

