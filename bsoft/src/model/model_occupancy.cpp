/**
@file	model_occupancy.cpp
@brief	Library routines to count components in maps
@author Daniel Nemecek and Bernard Heymann
@date	Created: 20091202
@date	Modified: 20151110
**/

#include "model_occupancy.h"
#include "binomial.h"
#include "ps_plot.h"
#include "img_combine.h"
#include "model_extract_build.h"
#include "model_compare.h"
#include "model_util.h"
#include "math_util.h"
#include "utilities.h"
#include "linked_list.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

//Internal function prototypes
double		component_coverage(Bimage* p, int img_num, Bcomponent* comp, double threshold);

/**
@author Daniel Nemecek and Bernard Heymann
@brief     Calculates the occupancy of components in a map.
@param 	*model		model structure to be modified.
@param	*pmask
@param 	mol_weight	molecular weight to determine threshold.
@param 	rho			protein density in Da/A3.
@param 	cutoff		coverage cutoff to use for assigning occupancy.
@param 	invert_flag		flag to invert map density.
@return int					0.

	The map must have positive density (higher values are considered density).
	The coverage for a component is defined as the fraction of voxels within
	a sphere around the component location that is above the threshold.
	The threshold can be derived from the molecular weight or is taken as
	the sum of the map average and standard deviation.
	The radius of the sphere is taken from the component radius.

**/
int			model_occupancy(Bmodel* model, Bimage* pmask, double mol_weight, double rho, double cutoff, int invert_flag)
{
	double			threshold;
	double			f, sigma;
	Vector3<double>	shift, loc;
	Bmodel*			mp;
	Bcomponent*		comp;
	Bimage*			p = NULL;

	if ( verbose & VERB_PROCESS ) {
		cout << "Model component occupancy:" << endl;
		cout << "Molecular weight:               " << mol_weight << " Da (rho = " << rho << " Da/Å)" << endl;
		cout << "Component coverage cutoff:      " << cutoff << endl << endl;
	}
	
	for ( mp = model; mp; mp = mp->next ) {
		p = read_img(mp->mapfile(), 1, mp->image_number());
		
		if ( invert_flag ) p->invert();

		//Determine an intensity treshold for the subimage
		if ( mol_weight > 0 )
			threshold = p->mass_threshold(0, mol_weight, rho);
		else
			threshold = p->average() + p->standard_deviation();
			
        if ( verbose & VERB_PROCESS )
			cout << "Intensity Threshold for MW = " << mol_weight << " is " << threshold << endl;

        //Mask the shell or other region in the subimage
//        if ( pmask ) img_multiply(p, pmask, 1, 0);
		if ( pmask ) p->multiply(pmask);
//		write_img("masked_img.pif", p);
	
		//Evaluate intensity of voxels for each component and decide its occupancy
        if ( verbose & VERB_PROCESS )
			cout << "Model\tComponent\tdx\tdy\tdz\tDensity\tCoverage\tOccupancy" << endl;
		for ( comp = mp->comp; comp; comp = comp->next ) {  
			loc = comp->location()/p->sampling(0) + p->image->origin();  // location of the component in the map (Å->vox)
//			comp->density(img_density_in_marker(p, loc, comp->radius()/p->sampling(0)[0], sigma)); // weighted average of intensity
			comp->density(p->density(0, loc, comp->radius()/p->sampling(0)[0], sigma)); // weighted average of intensity
			f = component_coverage(p, 0, comp, threshold);
			if ( f >= cutoff ) comp->select(1);
			else comp->select(0);
			comp->FOM(f);
			shift = comp->location() - comp->force();
			if ( verbose & VERB_PROCESS )
				cout << mp->identifier() << tab << comp->identifier() << tab << fixed << 
					setprecision(2) << shift[0] << tab << shift[1] << tab << shift[2] << tab <<
					setprecision(4) << comp->density() << tab << comp->FOM() << tab << comp->select() << endl;
		}
        if ( verbose & VERB_PROCESS )
			cout << endl;
	
		delete p;
	}
	
	return 0;
}

/**
@author Daniel Nemecek and Bernard Heymann
@brief	Calculates the occupancy distribution of models.
@param 	*model			model structure to be modified.
@param 	cutoff			coverage cutoff to determine occupancy.
@param 	nfit			number of binomial curves to fit.
@param 	&ncomp			maximum number of components in a model.
@param 	*prob			weight and probability array (2*nfit).
@param 	&R				pointer to fit residual.
@return vector<double>&	occupancy distribution histogram and fits.

	The component coverage must already be calculated and stored in the
	FOM property of each component. A component is considered occupied
	if its coverage exceeds the given cutoff value.
	The distribution, error and fit array is set up with 3+nfit columns:
		Column1:	Distribution histogram
		Column2:	Standard deviations
		Column3:	Overall binomial fit curve (sum of remaining columns)
		Column4+:	nfit individual binomial curves

**/
vector<double>	model_occupancy_distribution(Bmodel* model, double cutoff, int nfit,
				long& ncomp, vector<double>& prob, double& R)
{
	long			i, j, k, h, occ, bins, nsets(3+nfit);
	double			c, s;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	ncomp = model_maxnum_components(model);
	bins = ncomp + 1;
			
	vector<double>	dist(bins);
	vector<double>	disterrfit(bins*nsets);
	vector<double>	sum(bins);
	vector<double>	sum2(bins);
	
	for ( i=k=0; i<bins; i++ ) {
		dist[i] = sum[i] = sum2[i] = 0;
		for ( j=0; j<nsets; j++, k++ ) disterrfit[k] = 0;
	}
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( occ=0, comp = mp->comp; comp; comp = comp->next )
			if ( comp->FOM() >= cutoff ) occ++;
		disterrfit[occ]++;
	}
	for ( s=0, i=0; i<bins; i++ ) s += disterrfit[i];   
	for ( i=0; i<bins; i++ ) disterrfit[i] /= s;
	
	
	for ( c = cutoff - 0.2; c <= cutoff + 0.2; c += 0.1 ) {
		for ( i=0; i<bins; i++ ) dist[i] = 0;
		for ( mp = model; mp; mp = mp->next ) {
			for ( occ=0, comp = mp->comp; comp; comp = comp->next )
				if ( comp->FOM() >= c ) occ++;
			dist[occ]++;
		}
		for ( s=0, i=0; i<bins; i++ ) s += dist[i];
		for ( i=0; i<bins; i++ ) {
			dist[i] /= s;
			sum[i] += dist[i];
			sum2[i] += dist[i]*dist[i];
		}
	}

	for ( i=0, j=bins; i<bins; i++, j++ ) {
		s = sum2[i]/5 - sum[i]*sum[i]/25;
		if ( s > 0 ) disterrfit[j] = sqrt(s);
	}
	
	R = binomial_fit(disterrfit, ncomp, nfit, prob);
	
	// Order distribution fits in ascending probability order
	for ( i=0, k=nfit; i<nfit-1; i++, k++ ) {
		for ( j=i+1, h=j+nfit; j<nfit; j++, h++ ) {
			if ( prob[k] > prob[h] ) {
				swap(prob[i], prob[j]);
				swap(prob[k], prob[h]);
			}
		}
	}

	for ( i=0, j=3; i<nfit; i++, j++ ) {
		for ( k=0, h=j*bins; k<bins; k++, h++ )
			disterrfit[h] = Bnpk(ncomp, prob[nfit+i], k, prob[i]);
		for ( k=0, h=2*bins; k<bins; k++, h++ )
			disterrfit[h] += disterrfit[j*bins+k]; 
	}
			
	// Print results
	if ( verbose ) {
		cout << endl << "Binomial distributions:" << endl;
		cout << "#\tWeight\tProbability" << endl;
		for ( i=0; i<nfit; i++)
			cout << i+1 << tab << prob[i] << tab << prob[nfit+i] << endl;
		cout << "Residual:                       " << R << endl;
		cout << endl << "Occup\t Models\t Std\t Fit";
		for ( j=0; j<nfit; j++) cout << "\tB" << j+1;
		cout << endl;
 
		for ( i=0; i<bins; i++ ) {
			cout << i << tab << disterrfit[i] << tab << disterrfit[bins+i] << tab << disterrfit[2*bins+i];
			for ( j=0; j<nfit; j++)
				cout << "\t " << disterrfit[(3+j)*bins+i];
			cout << endl;
		}
		cout << endl;
	}
		
	return disterrfit;
}

/**
@author Daniel Nemecek and Bernard Heymann
@brief 	Refines component views and positions by cross-correlation.
@param 	*model		model.
@param 	*pmask2    	mask for the input image.
@param 	*ptemp		density template.
@param 	*pmask		cross-correlation mask.
@param 	hires		high resolution limit for cross-correlation.
@param 	lores		low resolution limit for cross-correlation.
@param 	max_shift	maximum shift in coordinates (angstrom).
@return int			0.

	The density origin is positioned on the component.
	The component views must already be set.
	The density and search radii are derived from the radius of each component.
	When a shift is out-of-range, returns to the initial position.
	The size of the template determines the search area.

**/
int			model_refine_comp_for_occupancy(Bmodel* model, Bimage* pmask2, Bimage* ptemp, Bimage* pmask,
									 double hires, double lores, double max_shift)
{
	if ( max_shift <= 0 ) max_shift = 1e37;
	
	Bimage*			p = read_img(model->mapfile(), 1, model->image_number()); // reads only the appropriate subimage from the image file for the given model
	if ( p == NULL ) {
		cerr << "Error: No model map file read! (" << model->mapfile() << ")" << endl;
		bexit(-1);
	}
	else if ( pmask2 ) 
//	    img_multiply(p, pmask2, 1, 0);
		p->multiply(pmask2);
	
	
	p->change_type(Float);
	ptemp->change_type(Float);

	if ( hires > lores ) swap(hires, lores);
	
	double	searchradius = model->comp->radius()*2/ptemp->sampling(0)[0];	// radius for cross-correlation (max. shift about 2r)
	double	radius = model->comp->radius()/ptemp->sampling(0)[0];			// radius for determination of intensity [avg,std] within a given component

	
	if ( verbose ) {
		cout << "Refinement of component locations for model " << model->identifier() << ":" << endl;
		cout << "Radius for determination of average intensity: " << radius << " voxels" << endl;
		cout << "Radius for cross-correlation search:  " << searchradius << " voxels" << endl << endl;
		cout << "Refining component positions:" << endl;
		cout << "Resolution limits:              " << hires << " " << lores << " A" << endl;
		cout << "Maximum shift:                  " << max_shift << " A" << endl;
		cout << "Search radius:                  " << searchradius << endl;
		cout << "Density radius:                 " << radius << endl << endl;
	}
	
	int				ncomp(0);
	long			n, imax(0);
	double			cc, dist12, dist10, dist20;
	double			sigma1, sigma2;
	double			ccmax, s, d, ccavg(0), d1avg(0), d2avg(0);
	Vector3<double>	zero_origin;
	Vector3<int>	size(ptemp->size());
	Vector3<double>	origin(ptemp->image->origin());
	Vector3<double>	map_origin(p->image->origin());
	Vector3<double>	scale(ptemp->sampling(0));
	Vector3<double>	r, loc, comp_shift, shift_avg;
	Vector3<double>	shift;
	Bimage*			pcomp = NULL;
	Bimage*			pone = NULL;
	Bimage*			pmaskrot = NULL;
	
	Bcomponent*		comp;
	Bcomponent*		comp2;
	Matrix3			mat;
	int*			num = new int[ptemp->images()];
	for ( n=0; n<ptemp->images(); n++ ) num[n] = 0;
	
	double			diameter; // in angstroms
	
	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		ccmax = -1;
		loc = comp->location()/scale + map_origin;  // location of the component in the map (Å->vox)
		comp->density(p->density(0, loc, radius, sigma1)); // weighted average of intensity
		d1avg += comp->density();
		sigma2 = sigma1;
		comp->force(comp->location());
		comp_shift = 0;
		mat = comp->view().matrix();  // transformation matrix (local->global?)
		//	Extracts a density associated with a set of coordinates in a model. (size from the density template, only one component)
		pcomp = p->extract(0, loc, size, origin, mat);
		// Transform of the image mask
		if ( pmask ) pmaskrot = pmask->extract_wrap(0, size, mat);
		for ( n=0; n<ptemp->images(); n++ ) { // cycle over subimages in the template
			pone = ptemp->extract(n); //extracts image of the template
			pone->origin(origin);
			// find a shift to match the extracted image density and the template
			shift = pcomp->find_shift(pone, pmaskrot, hires, lores, searchradius, 0, 0, cc);
			if ( verbose & VERB_FULL )
				cout << "Shift = " << shift << endl;
			if ( ccmax < cc ) { // if the cross-correlation coef is higher than before
				ccmax = cc;
				imax = n;
				comp_shift = mat * shift; // component shift in global coordinates
			}
			delete pone;
		}
		delete pcomp;
		delete pmaskrot;
 		
		loc += comp_shift;			 // new location of the component (in voxels)
		comp->location(scale * (loc - map_origin));
		comp_shift = comp->location() - comp->force();  // shift of the component in Å
		num[imax]++;			// stores number of copies of each subimage of the template
		p->density(0, loc, radius, sigma2); // weighted average of intensity
		if ( verbose & VERB_FULL )
			cout << "Sigma1=" << sigma1 << "\t Sigma2=" << sigma2 << endl;
		s = comp_shift.length();

		// Test for overlapping components
		for ( comp2 = model->comp; comp2 && comp2 != comp; comp2 = comp2->next ) {
			dist12 = comp->location().distance(comp2->location());
			diameter = comp->radius() + comp2->radius();
			if ( dist12 < diameter ) {
				dist10 = comp->location().distance(comp->force()); 
				dist20 = comp2->location().distance(comp2->force());
				if ( dist10 < dist20 )	{ // zero the longer shift
					comp2->location(comp2->force());
				} else {				
					comp->location(comp->force());
					s = 0;
				}
			}
		}
		
		// shift only when within the maximal allowed radius and sigma of intensity is smaller
		if ( s > 0 && s < max_shift && sigma1 > sigma2 ) {		
			comp->select(imax);			// select alignment to the best template subimage
			comp->FOM(ccmax);			// update fom to the highest cc
		} else {
			comp->location(comp->force());
			if ( verbose & VERB_FULL )
				cout << "Component " << ncomp+1 << " was not shifted!" << endl;
		}
		
	}

	if ( verbose )
		cout << "Comp\tdx\tdy\tdz\tShift\tCC\tDens1\tDens2\tdDens" << endl;

	for ( ncomp=0, comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		d = comp->density();
//		comp->density(img_density_in_marker(p, comp->location()/scale + map_origin, radius, sigma2));
		comp->density(p->density(0, comp->location()/scale + map_origin, radius, sigma2)); // weighted average of intensity
		d2avg += comp->density();
		comp_shift = comp->location() - comp->force();
		shift_avg += comp_shift;
		ccavg += comp->FOM();
		if ( verbose )
			cout << ncomp+1 << tab << comp_shift[0] << tab << comp_shift[1] 
				<< tab << comp_shift[2] << tab << comp_shift.length() << tab 
				<< comp->FOM() << tab << d << tab << comp->density() << tab << comp->density() - d << endl;
	}
	
	shift_avg /= ncomp;		// average shift (vector)
	s = shift_avg.length();		// average length of the shift
	ccavg /= ncomp;			// average cc	
	d1avg /= ncomp;			// average density in the original components
	d2avg /= ncomp;			// average density in the shifted components
	
	if ( verbose ) {
		cout << "Average\t" << shift_avg[0] << tab << shift_avg[1] << tab 
			<< shift_avg[2] << tab << s << tab << ccavg << tab << d1avg 
			<< tab << d2avg << tab << d2avg - d1avg << endl << endl;
		cout << "Component type distribution:\n#\tCount" << endl;
		for ( n=0; n<ptemp->images(); n++ ) cout << n+1 << tab << num[n] << endl;
		cout << endl;
	}
	
	delete[] num;
	delete p;
	
	return 0;
}

/*
@author Daniel Nemecek and Bernard Heymann
@brief	Calculates the coverage of the component in a map.
@param 	*p			image.
@param 	img_num		sub-image number.
@param	comp		component.
@param 	threshold	threshold to count covered voxels.
@return	double		fraction of examined sphere covered.

	The map must have positive density (higher values are considered density).
	The coverage for a component is defined as the fraction of voxels within
	a sphere around the component location that is above the threshold.
	The radius of the sphere is taken from the component radius.
**/
double		component_coverage(Bimage* p, int img_num, Bcomponent* comp, double threshold)
{
	int					n(0), np(0);
	long		i, x, y, z;
	double				a(1), b(1), c(1), fx, fy, fz;
	double				value, fraction(0);
	Vector3<int>		start;
	Vector3<int>		end;

	start[0] = (int) (p->image->origin()[0] + (comp->location()[0] - comp->radius())/p->sampling(0)[0]);
	start[1] = (int) (p->image->origin()[1] + (comp->location()[1] - comp->radius())/p->sampling(0)[1]);
	start[2] = (int) (p->image->origin()[2] + (comp->location()[2] - comp->radius())/p->sampling(0)[2]);
	end[0] = (int) (p->image->origin()[0] + (comp->location()[0] + comp->radius())/p->sampling(0)[0]);
	end[1] = (int) (p->image->origin()[1] + (comp->location()[1] + comp->radius())/p->sampling(0)[1]);
	end[2] = (int) (p->image->origin()[2] + (comp->location()[2] + comp->radius())/p->sampling(0)[2]);
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG component_coverage: start=" << start << endl;
		cout << "DEBUG component_coverage: end=" << end << endl;
	}
	
	if ( ( start[0] < 0 ) && ( end[0] < 0 ) ) return 0;
	if ( ( start[1] < 0 ) && ( end[1] < 0 ) ) return 0;
	if ( ( start[2] < 0 ) && ( end[2] < 0 ) ) return 0;
	if ( ( start[0] >= p->sizeX() ) && ( end[0] >= p->sizeX() ) ) return 0;
	if ( ( start[1] >= p->sizeY() ) && ( end[1] >= p->sizeY() ) ) return 0;
	if ( ( start[2] >= p->sizeZ() ) && ( end[2] >= p->sizeZ() ) ) return 0;
	start = start.max(0);
	end = end.max(0);
	if ( start[0] >= p->sizeX() ) start[0] = p->sizeX() - 1;
	if ( start[1] >= p->sizeY() ) start[1] = p->sizeY() - 1;
	if ( start[2] >= p->sizeZ() ) start[2] = p->sizeZ() - 1;
	if ( end[0] >= p->sizeX() ) end[0] = p->sizeX() - 1;
	if ( end[1] >= p->sizeY() ) end[1] = p->sizeY() - 1;
	if ( end[2] >= p->sizeZ() ) end[2] = p->sizeZ() - 1;	
	if ( start[0] > end[0] ) swap(start[0], end[0]);
	if ( start[1] > end[1] ) swap(start[1], end[1]);
	if ( start[2] > end[2] ) swap(start[2], end[2]);
	
	Vector3<int>		size = end - start;
	Vector3<double>		center;
	center[0] = start[0] + size[0]/2.0;
	center[1] = start[1] + size[1]/2.0;
	center[2] = start[2] + size[2]/2.0;
	if ( size[0] ) a = 2.0/size[0];
	if ( size[1] ) b = 2.0/size[1];
	if ( size[2] ) c = 2.0/size[2];
	a *= a;
	b *= b;
	c *= c;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG component_coverage: start=" << start << endl;
		cout << "DEBUG component_coverage: end=" << end << endl;
		cout << "DEBUG component_coverage: center=" << center << endl;
	}

	for ( n=0, z=start[2]; z<=end[2]; z++ ) {
		fz = z - center[2];
		fz = c*fz*fz;
		for ( y=start[1]; y<=end[1]; y++ ) {
			fy = y - center[1];
			fy = b*fy*fy;
			for ( x=start[0]; x<=end[0]; x++ ) {
				fx = x - center[0];
				fx = a*fx*fx;
				if ( fx + fy + fz <= 1 ) {
					i = p->index(0, x, y, z, img_num);
					value = (*p)[i];
					if ( value > threshold ) np++;
					n++;
 				}
			}
		}
	}
	
	if ( n ) fraction = np*1.0L/n;
	
	return fraction;
}



