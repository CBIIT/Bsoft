/**
@file	radon_util.cpp
@author P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "radon_util.h"
#include "utilities.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

/**
@brief 	Rotates a vector of cartesian coordinates along the Z axis
@param 	x			x coordinate of starting vector
@param 	y			y coordinate of starting vector
@param 	z			z coordinate of starting vector
@param 	ang			rotation angle (in radiants)
@param 	*xx			x coordinate of rotated vector
@param 	*yy			y coordinate of rotated vector
@param 	*zz			z coordinate of rotated vector
@return int			0
**/
int
rotZ_cart (double x, double y, double z, double ang, double *xx, double *yy, double *zz)
{
  *xx = cos (ang) * x - sin (ang) * y;
  *yy = sin (ang) * x + cos (ang) * y;
  *zz = z;

  return 0;
}

/**
@brief 	Rotates a vector of cartesian coordinates along the Y axis
@param 	x			x coordinate of starting vector
@param 	y			y coordinate of starting vector
@param 	z			z coordinate of starting vector
@param 	ang			rotation angle (in radiants)
@param 	*xx			x coordinate of rotated vector
@param 	*yy			y coordinate of rotated vector
@param 	*zz			z coordinate of rotated vector
@return int			0
**/
int
rotY_cart (double x, double y, double z, double ang, double *xx, double *yy, double *zz)
{
  *xx = cos (ang) * x + sin (ang) * z;
  *yy = y;
  *zz = -sin (ang) * x + cos (ang) * z;

  return 0;
}

/**
@brief 	Rotates a vector of cartesian coordinates along the X axis
@param 	x			x coordinate of starting vector
@param 	y			y coordinate of starting vector
@param 	z			z coordinate of starting vector
@param 	ang			rotation angle (in radiants)
@param 	*xx			x coordinate of rotated vector
@param 	*yy			y coordinate of rotated vector
@param 	*zz			z coordinate of rotated vector
@return int			0
**/
int
rotX_cart (double x, double y, double z, double ang, double *xx, double *yy, double *zz)
{
  *xx = x;
  *yy = cos (ang) * y - sin (ang) * z;
  *zz = sin (ang) * y + cos (ang) * z;

  return 0;
}

/**
@brief 	Calculates the angle (in radians) from 2 vectors
@param 	a1			first coordinate of first vector
@param 	b1			second coordinate of first vector
@param 	a2			first coordinate of second vector
@param 	b2			second coordinate of second vector
@return double		angle
**/
double
ang_one_two (double a1, double b1, double a2, double b2)
{
  double coseno;
  double sb1, sb2;

  sb1 = sin (b1);
  sb2 = sin (b2);
  coseno = sb1 * cos (a1) * sb2 * cos (a2) + sb1 * sin (a1) * sb2 * sin (a2) + cos (b1) * cos (b2);

  if (coseno > 1)
    coseno = 1.;
  if (coseno < -1)
    coseno = -1.;

  double ang = acos (coseno);

  return (ang);
}

/**

@param	*p	
@return int			0
**/
int
sphere (Bimage * p)
{
  int i, j, k, d, dmax, cont, ll, hx, hy, hz;
  int rmax, r_int;
  float mean, div, rl, rk, fct;
  float wpix, wmean, i2, j2, k2;
  float *pr, *pr1, *pr2;

  pr = (float *) p->data_pointer();

	hx = p->sizeX() / 2;
	hy = p->sizeY() / 2;
	hz = p->sizeZ() / 2;
  ll = 2;
  dmax = hz - 1;
  pr1 = pr;
  cont = 0;
  mean = 0;
  for (k = -hz; k < hz; k++)
    for (j = -hy; j < hy; j++)
      for (i = -hx; i < hx; i++)
	{
	  d = k * k + j * j + i * i;
	  d = (int) sqrt (d);
	  if (d >= dmax - ll && d < dmax)
	    {
	      cont++;
	      mean += *pr1;
	    }
	  pr1++;
	}
  mean /= cont;
  rmax = (hz - 1);
  r_int = (hz - ll * 2 - 1);
  div = (float) (2 * ll);
  div = M_PI / div;
  pr1 = pr;
  for (k = -hz; k < hz; k++)
    {
      k2 = k * k;
      for (i = -hy; i < hy; i++)
	{
	  i2 = i * i;
	  for (j = -hx; j < hx; j++)
	    {
	      j2 = k2 + i2 + j * j;
	      fct = sqrt (j2);
	      if (fct >= r_int && fct <= rmax)
		{
		  rl = fct - r_int;
		  rk = rmax - fct;
		  wpix = .5 + .5 * cos (rl * div);
		  wmean = .5 + .5 * cos (rk * div);
		  pr2 = pr1;
		  *pr1 = *pr2 * wpix + mean * wmean - mean;
		}
	      else
		{
		  if (fct > rmax)
		    *pr1 = 0.;
		  else
		    *pr1 -= mean;
		}
	      pr1++;
	    }
	}
    }

  if (verbose & VERB_FULL)
    cout << "count=" << cont << " mean=" << mean << endl << endl;

  return 0;
}

/**
@brief	Set each value equal to old_value-(first image pixel)
@param	*p     		image.
@return	int			0.
**/
int mean_to_0(Bimage *p)
{
  unsigned long i, ntot;
  float mean,*pr;


  pr = (float *) p->data_pointer();
  mean = *pr;
  ntot = p->sizeX() * p->sizeY() * p->sizeZ() * p->images();

  if ( verbose & VERB_PROCESS )
    printf("Shifting the values range of all images by %g\n\n", -mean);

  for (i = 0; i < ntot; i++,pr++)
     *pr -= mean;

  return(0);
}
