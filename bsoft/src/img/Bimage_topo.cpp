/**
@file	Bimage_topo.cpp
@brief	Functions to convert between 2D topology images to 3D surfaces.
@author Bernard Heymann
@date	Created: 19990124
@date	Modified: 20170613
**/

#include "Bimage.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Converts a 2D AFM image to a 3D density map.
@param 	*psd		2D standard deviation image.
@param 	nz			z dimension of the new 3D map.
@param 	density		density inside the surface (Da/A3).
@param 	resolution	sigma bounds.
@return int			0.

	The 2D image is expanded in the z direction according to:
	                          voxel_density
	dens[x,y,z] = --------------------------------------
	              1+exp(factor*(z-zdis[x,y])/sigma[x,y])
	where	voxel_density is the total density attributed to a voxel
			factor = 1.618
			z - zdis[x,y] is the z distance from the surface
			sigma is a standard deviation term to describe the surface variability.
	The sigma term can be set automatically or given by an input standard deviation map.
	The position of the surface is calculated as:
	zdis[x,y] = loz + (hiz - loz)*(image[x,y] - min)/(max - min)
	where	loz and hiz are the lowest and highest points of the surface in the 3D map
			min and max are the image minimum and maximum values

**/
Bimage*		Bimage::topograph_to_surface(Bimage* psd, long nz, double density, double resolution)
{
	// Density within the surface
	if ( density < 1e-10 ) density = 1;
	double			voxel_density(density*sampling(0).volume());
	int				local_psd(0);
	long			i, j, xx, yy, zz, slice_size(x*y);
    double			factor(GOLDEN);
	double			mn, mx, zdis, scale;
	
	if ( psd == NULL ) {
		local_psd = 1;
		psd = copy_header();
		mn = mx = resolution*0.5;
		psd->show_minimum(resolution*0.5);
		psd->show_maximum(resolution*0.5);
		scale = 1;
		psd->data_alloc();
		for ( i=0; i<slice_size; i++ )
			psd->set(i, psd->show_minimum());
	} else {
		mn = psd->minimum();
		mx = psd->maximum();
		scale = (psd->show_maximum() - psd->show_minimum())/(mx - mn);
		for ( i=0; i<slice_size; i++ )
			psd->set(i, scale*(*psd)[i] + psd->show_minimum());
	}
	
	Bimage*			psurf = new Bimage(Float, TSimple, x, y, nz, 1);
	psurf->sampling(sampling(0));
	
	if ( show_minimum() == 0 && show_maximum() == 0 ) {
		show_minimum(resolution);				
		show_maximum((nz - 1)*image->sampling()[2] - resolution);
	}
	
    double			loz(show_minimum()/image->sampling()[2]);
	scale = (show_maximum()-show_minimum())/(image->sampling()[2]*(max-min));
    if ( show_minimum() > show_maximum() ) factor *= -1;
	
    if ( verbose & VERB_PROCESS ) {
		cout << "Converting a 2D image into a 3D surface density" << endl;
		cout << "New z dimension:                " << nz << endl;
		cout << "Dimension scale:                " << sampling(0) << " A/pixel" << endl;
		cout << "Density scale:                  " << scale << " A/(unit density)" << endl;
		cout << "Z min & max:                    " << show_minimum() << " A  " << show_maximum() << " A" << endl;
		cout << "Z thickness:                    " << show_maximum() - show_minimum() << " A" << endl;
		cout << "Resolution:                     " << resolution << " A" << endl;
		if ( psd )
    		cout << "Sigma min & max:                " << psd->show_minimum() << " " << psd->show_maximum() << " A" << endl;
		cout << "Density within surface:         " << density << " Da/A3" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Converting a 2D image into a 3D surface density" << endl << endl;

	for ( i=yy=0; yy<y; yy++ ) {
		for ( xx=0; xx<x; xx++, i++ ) {
       		zdis = loz + scale*((*this)[i]-min);
    	    for ( j=i, zz=0; zz<psurf->z; zz++, j+=slice_size )
    			psurf->set(j, voxel_density/(1+exp(factor*(zz-zdis)/(*psd)[i])));
   		}
	}
    
	if ( local_psd ) delete psd;
    
	psurf->statistics();
	
    return psurf;
}

/**
@brief 	Converts a 3D image to a 2D height image.
@param 	threshold	threshold to define the surface (assuming positive density).
@param	dir			direction: 0=bottom up, 1=top down.
@return Bimage*		2D height image.

	The threshold defines the surface of the positive density in the 3D image.
	Each line of voxels in the z-direction is scanned for the voxel
	exceeding the threshold with the highest z-index.

**/
Bimage*		Bimage::surface_to_topograph(double threshold, int dir)
{
	long   			i, j, xx, yy, zz, nn;
	
	Bimage*			p = new Bimage(Float, TSimple, x, y, 1, n);
	
    if ( verbose & VERB_PROCESS ) {
		cout << "Converting a 3D image into a 2D height image" << endl;
		cout << "Threshold:                      " << threshold << endl << endl;
	}
	
	for ( i=nn=0; nn<n; nn++ ) {
		for ( zz=0; zz<z; zz++ ) {
			for ( yy=0; yy<y; yy++ ) {
				for ( xx=0; xx<x; xx++, i++ ) {
					j = (nn*y+yy)*x+xx;
					if ( !(*p)[j] ) {
						if ( dir ) {
							if ( (*this)[i] >= threshold ) p->set(j, zz);
						} else {
							if ( (*this)[i] <= threshold ) p->set(j, zz);
						}
					}
				}
			}
		}
	}
		
	return p;
}

/**
@brief 	Rotates a 3D map and calculates the height along the z-axis.
@param 	mat				3x3 rotation or skewing matrix.
@param 	translate		3-value vector for translation after transformation.
@param 	threshold		density threshold to consider as object.
@return Bimage*			new 2D height image.

	The 3D map is rotated around its center as origin and the height 
	calculated along the z-direction.
	The resultant 2D image is translated if the translation vector is 
	non-zero.
	The rotation origin is obtained from the map origin.

**/
Bimage*		Bimage::rotate_height(Matrix3 mat, Vector3<double> translate, double threshold)
{
	long			i, xx, yy, zz;
	Vector3<double>	old, nuvec;
	Vector3<double>	u(image->sampling());
	
	double			radial_cutoff = x;
	if ( radial_cutoff < 1 ) radial_cutoff = x;
//	double			rad_sq = radial_cutoff*radial_cutoff;
	
	Vector3<double>	oldorigin(image->origin());

	Vector3<double>	nuorigin = oldorigin + translate;
		
	Matrix3			matinv = mat.transpose();
	
	if ( verbose & VERB_FULL ) {
		cout << "Rotating around the map center and calculating height:" << endl;
		cout << "Rotation origin:                " << oldorigin << endl;
		cout << matinv << endl;
		if ( translate[0] != 0 || translate[1] != 0 || translate[2] != 0 )
			cout << "Translation:                    " << translate << endl;
//		cout << "Radial cutoff:                  " << radial_cutoff << " pixels" << endl;
		cout << "Threshold:                      " << threshold << endl;
		cout << endl;
	}

	Bimage*			pheight = new Bimage(Float, TSimple, x, y, 1, 1);
	pheight->sampling(u[0], u[1], 1);
	pheight->origin(nuorigin[0], nuorigin[1], 0);

	double			d, dz, bottom(u[2]*z);
	
	for ( i=yy=0; yy<y; yy++ ) {
		nuvec[1] = yy - nuorigin[1];
		for ( xx=0; xx<x; xx++, i++ ) {
			nuvec[0] = xx - nuorigin[0];
			for ( zz=z-1; zz; zz-- ) {
				nuvec[2] = zz - nuorigin[2];
				old = mat * nuvec + oldorigin;
				d = interpolate(old);
				if ( d > threshold ) {
					dz = u[2]*zz;
					if ( (*pheight)[i] <= 0 ) pheight->set(i, dz);
					if ( bottom > dz ) bottom = dz;
				}
			}
		}
	}
	
	for ( i=0; i<x*y; ++i ) if ( (*pheight)[i] ) pheight->add(i, -bottom);
	
	pheight->statistics();
	
	return pheight;
}

/**
@brief 	Calculates a set of height images from a 3D density map.
@param 	*views			linked list of views.
@param 	threshold		density threshold to consider as object.
@return Bimage* 		2D height images as sub-images.

	A set of height images is calculated according to a list of views.

**/
Bimage* 	Bimage::height(View* views, double threshold)
{
	if ( !views ) {
		error_show("Error in Bimage::height: No views defined!", __FILE__, __LINE__);
		return NULL;
	}
	
	calculate_background();
	
	long			i;
	Vector3<double>	translate;
	Vector3<double>	u(image->sampling());
	u[2] = 1;

	long			nviews = count_list((char *)views);
	View*			vlist = new View[nviews];
//	vector<View>	vlist(nviews);
	View*			v;
	
	for ( v=views, i=0; v; v = v->next, i++ ) vlist[i] = *v;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Calculating height images:" << endl;
		cout << "Number of images:               " << nviews << endl;
		cout << "Threshold:                      " << threshold << endl << endl;
	}

	Bimage*			pheight = new Bimage(Float, TSimple, x, y, 1, nviews);
	pheight->file_name(file_name());
	pheight->sampling(u);

#ifdef HAVE_GCD
	dispatch_apply(nviews, dispatch_get_global_queue(0, 0), ^(size_t i){
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	one_height = rotate_height(mat, translate, threshold);
		one_height->shift_background(0);
		pheight->replace(i, one_height);
//		pheight->image[i].origin(one_height->image->origin());
		pheight->origin(i, one_height->image->origin());
		pheight->image[i].view(vlist[i]);
		delete one_height;
	});
#else
#pragma omp parallel for
	for ( long i=0; i<nviews; i++ ) {
//		if ( verbose & VERB_LABEL )
//			cout << "Projection " << i+1 << ":" << endl;
		Matrix3	mat = vlist[i].matrix();
		Bimage*	one_height = rotate_height(mat, translate, threshold);
		one_height->shift_background(0);
		pheight->replace(i, one_height);
//		pheight->image[i].origin(one_height->image->origin());
		pheight->origin(i, one_height->image->origin());
		pheight->image[i].view(vlist[i]);
		delete one_height;
	}
#endif

	delete[] vlist;
	
	pheight->statistics();
	
	return pheight;
}
