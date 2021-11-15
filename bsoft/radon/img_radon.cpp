/**
@file	img_radon.cpp
@brief	Radon transform functions.
@author Salvatore Lanzavecchia, Francesca Cantele and Pier Luigi Bellon
         Dip. Chimica Strutturale e Stereochimica Inorganica, Via Venezian 21, 20133 Milano, Italy
@author	Bernard Heymann
         Rm 1515, 50 South Dr., NIH, Bethesda, MD, 20892, USA
@date	Created: 2001 05 19
@date	Modified: 20160728 (BH)
**/

#include "kernlut.h"
#include "gal_ricoI.h"
#include "sin_gal.h"
#include "pocs_tool.h"
#include "proj_tool.h"
#include "radon_util.h"
#include "img_radon.h"

#include "mg_processing.h"
#include "Bimage.h"
#include "symmetry.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

/**
@brief 	Calculates the Radon transform of an image.
@param 	*p				image to be transformed.
@param	type			type of transformation.
@param	nkernel			kernel size.
@param	kernel_power	kernel exponent (usually 2).
@param	padd			padding flag (0=no padding, 1=padding twofold).
@param	ntheta			number of angles in the 2*PI range.
@return	Bimage*			transformed image.

	Computes the Radon transform of a 3D image with all sides equal
	(a cube) and a power of 2. The transform is calculated in spherical
	coordinates where ntheta is the number of sampling points in the 
	range 0 - 2*PI along the two angular axes. Because half of this 
	range is redundant (0-PI is enough), a smaller output can be 
	generated.
        Type:
	1 = From 3D structure to single axis proj., output is a gallery of ntheta 2D image with size (x,x)
	2 = From single axis proj. to quarter radon transform, output is a 3D image with size (x,ntheta/2,ntheta/2)
	3 = From 3D structure to full radon transform, output is a 3D image with size (x,ntheta,ntheta)
	4 = From 3D structure to quarter radon transform, output is a 3D image with size (x,ntheta/2,ntheta/2)

**/
Bimage *
img_radon_transform (Bimage * p, int type, int nkernel, int kernel_power, int padd, int ntheta)
{
  char order[12] = "xyz";
  Vector3<double> shift;

  kernel *ker = kernel_prep (nkernel, kernel_power);

  p->change_type(Float);

  if (type != 2)
    sphere (p);
  else
    mean_to_0(p);

  int nang = 180;

  if (type == 1 || type == 3)
    nang = 360;

  if (verbose & (VERB_LABEL | VERB_PROCESS))
    cout << "Forward radon transform:" << endl;

  if (verbose & VERB_PROCESS)
    {
      cout << "Kernel size and power:          " << nkernel << " " << kernel_power << endl;
      cout << "Type:                           " << type << endl;
      cout << "Number of angles:               " << ntheta << endl;
      cout << "Padding:                        " << padd << endl;
      if (type != 3)
	cout << "Transform size:                 PI x PI\n" << endl;
      else
	cout << "Transform size:                 2 PI x 2 PI\n" << endl;
    }

  Bimage *prad = p;
  Bimage *prad2 = NULL;

  if (type == 2)
    {
      ntheta = prad->sizeZ();
//      prad->z /= 2;
    }
  else
    {
      strcpy (order, "zxy");
		p->reslice(order);
//write_img("p.pif", p);

      prad = sin_gal (p, nang, shift, ker, padd, ntheta);
//write_img("prad.pif", prad);
      strcpy (order, "xzy");
		prad->reslice(order);
    }

  if (type != 1)
    {
      prad2 = sin_gal (prad, nang, shift, ker, padd, ntheta);
      if (prad != p) delete prad;
      prad = prad2;

      strcpy (order, "yzx");
		prad->reslice(order);
    }

  free (ker->kern);
  free (ker);

  return (prad);
}

/**
@brief	Calculates a 3D image from a Radon transform.
@param 	*p				radon transform to be back-transformed.
@param	type			type of transformation.
@param	nkernel			kernel size.
@param	kernel_power	kernel exponent (usually 2).
@param	padd			padding flag (0=no padding, 1=padding twofold).
@return	Bimage*			transformed image.

	Type:
	1 = From single axis proj. to 3D structure, output is a 3D image with size (z,z,z)
	2 = From quarter radon transform to single axis proj., output is a gallery of ntheta 2D image with size (z,z)
	3 = From full radon transform to 3D structure, output is a 3D image with size (z,z,z)
	4 = From quarter radon transform to 3D structure, output is a 3D image with size (z,z,z)

**/
Bimage *
img_radon_inverse_transform (Bimage * p, int type, int nkernel, int kernel_power, int padd)
{
  char order[12] = "xyz";

  kernel *ker = kernel_prep (nkernel, kernel_power);

  p->change_type(Float);

  mean_to_0(p);

  if (verbose & (VERB_LABEL | VERB_PROCESS))
    cout << "Backward radon transform:" << endl;

  if (verbose & VERB_PROCESS)
    {
      cout << "Kernel size and power:          " << nkernel << " " << kernel_power << endl;
      cout << "Type:                           " << type << endl;
      cout << "Padding:                        " << padd << endl << endl;
    }

  Bimage *pinv = p;
  Bimage *pinv2;

  if (type == 2)
    {
      pinv = copy_I_part_rad (p);
      copy_II_part_rad (pinv);
//      p = pinv;
    }

  if (type != 1)
    {
      strcpy (order, "zxy");
		p->reslice(order);

      pinv2 = gal_ricoI (pinv, type, ker, padd);
	  if ( p != pinv ) delete pinv;
	  pinv = pinv2;
    }

  if (type != 2)
    {

      strcpy (order, "xzy");
		pinv->reslice(order);

      pinv2 = gal_ricoI (pinv, type, ker, padd);
      if (p != pinv) delete pinv;
      pinv = pinv2;

      strcpy (order, "yzx");
		pinv->reslice(order);
    }

  free (ker->kern);
  free (ker);

  return (pinv);
}

/**
@brief	Filters a radon transform using the POCS method.
@param 	*p				radon transform.
@param 	n_cyc_out		outer cycles, swapping r,phi and r,theta planes
@param 	n_cyc_in		inner cycles, within r,phi and r,theta planes
@param 	rad_3D			limiting radius in 3D.
@param 	rad_plane		limiting radius in plane.
@param 	support			flag to impose finiteness in real space.
@param 	*pmask			mask of dimension (ntheta/2,ntheta/2).
@return	Bimage*			transformed image.

	Filters a Radon transform to impose consistency and/or to fill holes.
	It can be used in two ways:
	a) When the Radon transform is not completely filled, it fills the holes
		described in the mask file produced by proj_to_radon.
	b) When the Radon transform is filled and noisy, it impose consistency.
	The mask must have the same dimensions as a plane in the radon transform.

**/
int
img_radon_pocs_filter (Bimage * p, int n_cyc_out, int n_cyc_in, double rad_3D, double rad_plane, int support, Bimage * pmask)
{
  char order[12] = "xyz";
  int i;

  p->change_type(Float);

  mean_to_0(p);

  if (verbose & (VERB_LABEL | VERB_PROCESS))
    cout << "POCS filtering of a radon transform:" << endl;

  if (verbose & VERB_PROCESS)
    {
      cout << "Outer and inner cycles:         " << n_cyc_out << " " << n_cyc_in << endl;
      cout << "Limiting radii, 3D and plane:   " << rad_3D << " " << rad_plane << endl;
      cout << "Finiteness support:             " << support << endl << endl;
    }

  // load Radon transform like theta,phi,rho //
  for (i = 0; i < n_cyc_out; i++)
    {

      if (support)
        support_f (p, rad_3D);

      //  execute T2 --> rho,theta,phi //
      strcpy (order, "zxy");
		p->reslice(order);

      filter_X_mask (p, pmask, n_cyc_in, rad_3D, 0);

      // go to rho,phi,theta //
      strcpy (order, "xzy");
		p->reslice(order);

      filter_X_mask (p, pmask, n_cyc_in, rad_plane, 1);

      // come back at theta,phi,rho //
      strcpy (order, "zyx");
		p->reslice(order);
    }

  if (support)
    support_f (p, rad_3D);

  return 0;
}

/**
@brief	Reconstructs a Radon transform from a set of projections.
@param 	*project		image processing parameter structure.
@param 	&sym			point group symmetry.
@param 	&file_mask		output file name of mask of dimension (ntheta,ntheta).
@param 	rec_size		reconstruction size (x,y,z)
@param 	ntheta			number of angles in the 2*PI range.
@param 	table_size		lookup table size.
@param 	threshold		threshold for rejecting images.
@param 	origin			origin reference for shifts.
@param 	nkernel			kernel size.
@param 	kernel_power	kernel exponent (usually 2).
@return	Bimage*			transformed image.

	The parameters are defined in the hierarchical project structure.
	The radon transform of each image is calculated and all its 
	symmetry-related views are written into the reconstruction volume.
	A mask image is calculated for the angular coverage of orientation
	space and used to weigh the reconstruction.

**/
Bimage *
img_radon_reconstruction (Bproject * project, Bsymmetry& sym, Bstring& file_mask, int rec_size, int ntheta, int table_size, double threshold, Vector3<double> origin, int nkernel, int kernel_power)
{
  int i, padd;
  double x, y, z, f, t, t1;

  Bfield *field = project->field;
  Bmicrograph *mg = field->mg;
  Bparticle *part = mg->part;

  if ( mg->fpart.length() < 1 )
    {
      cout << "No particle image file name given!" << endl << endl;
      exit (-1);
    }

  kernel *ker = kernel_prep (nkernel, kernel_power);

  Bimage *p = NULL;
  if (rec_size < 1)
	{
		p = read_img (mg->fpart, 0, -1);
		rec_size = findNextPowerOf(p->sizeX(), 2);
		delete p;
	}

  ntheta = rec_size * 2;

  if ( origin[0] < 1 ) origin[0] = rec_size / 2;
  if ( origin[1] < 1 ) origin[1] = rec_size / 2;
  if ( origin[2] < 1 ) origin[2] = rec_size / 2;

  if (verbose & (VERB_LABEL | VERB_PROCESS))
    cout << "Radon transform reconstruction from particle images:" << endl;

  if (verbose & VERB_PROCESS)
    {
      cout << "Reconstruction size:            " << ntheta / 2 << " " << ntheta / 2 << " " << rec_size << endl;
      cout << "Symmetry:                       " << sym.label() << endl;
      cout << "Kernel size and power:          " << nkernel << " " << kernel_power << endl;
      cout << "Number of angles:               " << ntheta << endl;
      cout << "Threshold:                      " << threshold << endl;
      if (file_mask.length())
	cout << "Mask file:                      " << file_mask << endl << endl;
      else
	cout << endl;
    }

  LUTable *lut = create_table (rec_size, table_size);

  Bimage *prec = new Bimage (Float, TSimple, ntheta / 2, ntheta / 2, rec_size, 1);
  prec->image->origin(origin);

  Bimage *pmask = new Bimage (Float, TSimple, ntheta / 2, ntheta / 2, 1, 1);

  if (verbose & VERB_PROCESS)
    cout << "Reconstruction dimensions:     " << prec->size() << endl;

  float *priga;
  if ((priga = (float *) calloc (prec->sizeZ(), sizeof (float))) == 0)
    {
      perror ("img_radon_reconstruction: priga error calloc");
      exit (0);
    }

  char order[12] = "yxz";
  int psel = 0;
  Vector3<double> d;
  Bimage *prad;
  Euler euler;
  View *view, *v;
//  int nview = sym.order();

  for (field = project->field; field; field = field->next)
    {
      for (mg = field->mg; mg; mg = mg->next)
	{
	  for (part = mg->part; part; part = part->next)
	    {
	      psel = 0;
	      if (threshold)
		{
		  if (part->fom[0] >= threshold)
		    psel = 1;
		}
	      else
		{
		  if (part->sel)
		    psel = 1;
		}
	      if (psel)		// Use only selected images
		{
		  p = read_img (mg->fpart, 1, part->id - 1);
		  if (!p)
		    {
		      cerr << "Error: Image file " << mg->fpart << " not read!" << endl;
		      exit (-1);
		    }
			p->change_type(Float);
		  p->background(0.0);
		  p->image->origin(part->ori);
		  img_resize_to_next_power2(p, FILL_BACKGROUND, p->background(long(0)));

		  padd = 1;
		  if ( p->image[0].origin()[0] > 0 && p->image[0].origin()[1] )
		      d = origin - p->image->origin();
		  prad = sin_gal (p, 180, d, ker, padd, ntheta);
		  delete p;

			prad->reslice(order);
		  view = symmetry_get_all_views (sym, part->view);

		for ( v=view; v; v=v->next )
		    {
		      euler = Euler(*v);
		      for (i = 0; i < ntheta / 2; i++)
			{
			  t1 = (float) i *2. * M_PI / ntheta;	// sinogram step //
			  x = 1.;
			  y = 0.;
			  z = 0.;
			  rotZ_cart (x, y, z, t1, &x, &y, &z);
			  rotZ_cart (x, y, z, euler.psi(), &x, &y, &z);
			  rotY_cart (x, y, z, euler.theta(), &x, &y, &z);
			  rotZ_cart (x, y, z, euler.phi(), &x, &y, &z);
			  f = asin (y);
			  t = atan2 (x, z);
			  write_line (prad, prec, i, f, t, pmask, lut, priga);
			}
		    }
		  delete prad;
		  kill_list((char *) view, sizeof(View));
		}
	    }
	}
    }

  free (priga);

  kill_table (lut);

  weigh_radon_transf (prec, pmask);

  if (file_mask.length())
    write_img(file_mask, pmask, 0);

  delete pmask;

  free (ker->kern);
  free (ker);

  return (prec);
}

/**
@brief 	Resizes without interpolation or rescaling to the next power of 2.
@param 	p				image (modified).
@param 	fill_type		FILL_AVERAGE, FILL_BACKGROUND, FILL_USER
@param 	fill			value to fill in new regions.
@return int				0.

	An image is resized to the next power of two in each dimension
	greater than 1 with translation and filling of new regions with 
	a given value.
	The new data replaces the old data.
**/
int			img_resize_to_next_power2(Bimage* p, int fill_type, double fill)
{
	Vector3<int>		newsize(1,1,1);
	
	if ( p->sizeX() > 1 ) newsize[0] = findNextPowerOf(p->sizeX(), 2);
	if ( p->sizeY() > 1 ) newsize[1] = findNextPowerOf(p->sizeY(), 2);
	if ( p->sizeZ() > 1 ) newsize[2] = findNextPowerOf(p->sizeZ(), 2);
	
	if ( newsize[0] == p->sizeX() && newsize[1] == p->sizeY() && newsize[2] == p->sizeZ() )
		return 0;
	
//	cout << "resizing" << endl;
	
	Vector3<int>	translate(newsize[0]/2 - p->sizeX()/2, newsize[1]/2 - p->sizeY()/2, newsize[2]/2 - p->sizeZ()/2);
	
	p->resize(newsize, translate, fill_type, fill);
	
	return 0;
}
