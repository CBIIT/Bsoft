/**
@file	kernlut.cpp
@author P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy
@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "kernlut.h"

#include "utilities.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

/**
@brief 	Prepare coefficients table for interpolation
@param	*ker	table of coefficients with its parameters
@return int     	0
**/
int
kernlut (kernel * ker)
{
  float *kern1;
  double s1, s2, c, arg1, arg2, ink, pink;
  int i, j, rnkern = ker->nkern, rnk = ker->nk, rnkm1;
  double fnkm1, arg, incr;

  if (ker->nk % 2 == 0)
    rnkm1 = ker->nk / 2 - 1;
  else
    rnkm1 = (ker->nk - 1) / 2;
  i = 1;
  ink = (float) i / (float) ker->nk;
  pink = M_PI * ink;
  incr = (float) 1 / (float) rnkern;
  kern1 = ker->kern + rnkm1;
  fnkm1 = (float) rnkm1;
  *kern1 = 1.;
  kern1 = ker->kern + rnk;
  arg = incr;
  for (i = 1; i < rnkern; i++)
    {
      arg1 = M_PI * (arg + fnkm1);
      s1 = ink * sin (arg1);
      arg2 = arg1 * ink;
      for (j = 0; j < rnk; j++)
	{
	  s2 = sin (arg2);
	  c = cos (arg2);
	  *kern1 = (s1 / s2) * pow (c, ker->power);
	  kern1++;
	  arg2 -= pink;
	  s1 = -s1;
	}
      arg += incr;
    }
  kern1 += rnkm1 + 1;
  *kern1 = 1.;
  return 0;
}


kernel *
kernel_prep (int nk, int power)
{
  kernel *ker = NULL;

  if ((ker = (kernel *) calloc (1, sizeof (kernel))) == NULL)
    {
      perror ("kernel_prep: ker error calloc");
      exit (0);
    }

  if (verbose & VERB_FULL)
    cout << "Preparing the kernel: nk = " << nk << ",  power = " << power << endl;

  ker->nk = nk;
  ker->power = power;

  ker->nkern = 512;
  ker->nker = ker->nkern * ker->nk + ker->nk;
  if ((ker->kern = (float *) calloc (ker->nker, sizeof (float))) == 0)
    {
      perror ("kernel_prep: ker->kern error calloc");
      exit (0);
    }

  if (ker->power == -1)
    {
      ker->bordo = 1;
      ker->nk = 1;
    }				// linear interp. //
  else if ((ker->nk % 2) == 0)
    ker->bordo = ker->nk / 2 - 1;
  else
    ker->bordo = (ker->nk - 1) / 2;

  if (ker->power >= 0)
    kernlut (ker);

  return (ker);
}
