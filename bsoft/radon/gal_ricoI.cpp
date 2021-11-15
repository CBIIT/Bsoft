/**
@file	gal_ricoI.cpp
@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy
@date	Created: 7 04 2003
@date	Modified: 11 08 2011 (BH)
**/

#include "gal_ricoI.h"
#include "fft_tool.h"

#include "Bimage.h"
#include "utilities.h"

int polcar_f (Bimage * p, float *pfr, float *pfi, int bordo, int n_radii, int n_theta, kernel * ker, int padd);
int polcart_fourier (float *pfri, float r, float phi, float *val, int bordo, int n_radii, int n_theta, kernel * ker);
int padding (Bimage * p, Bimage * ppiano, int z);
int copy180 (Bimage * p, float *pfr, float *pfi, int bordo);
int shift_phase (Bimage * p, float npix);


/**
@brief 	Reconstruct an image from a sinogram
@param	*p			parameters of the image
@param	tipo		type of reconstruction (can be 1/2/3/4)
@param	*ker		table of coefficients for interpolation
@param	padd		use of padding (0/1)
@return int			0

	float *psource	= sinogram
	float *pdest	= reconstructed image

**/
Bimage *
gal_ricoI (Bimage * p, int tipo, kernel * ker, int padd)
{
  int k;
  int n_radii, n_theta;
  float	t;
  float *pr, *pfr, *pfi, *pr1, *pk1;

  int bordo = ker->bordo;

  Bimage *prad = p->copy_header (0);
  prad->sizeY(prad->sizeX());

  float *pdest = (float *) malloc (prad->sizeX() * prad->sizeY() * prad->sizeZ() * sizeof (float));
  prad->data_assign((unsigned char *) pdest);

  float *psource = (float *) p->data_pointer();

  int pf_x = p->sizeX() / 2 + 2 * ker->bordo;
  int pf_y = p->sizeY() + 2 * ker->bordo;
  if (padd)
    pf_x = p->sizeX() * 2 * padd + 2 * ker->bordo;

  if ((pfr = (float *) malloc (2 * pf_x * pf_y * sizeof (float))) == 0)
    {
      perror ("sin_gal: pfr error calloc");
      exit (0);
    }
  pfi = pfr + pf_x * pf_y;

  Bimage *ppiano = new Bimage (Float, TSimple, p->sizeX() * (padd + 1), p->sizeY(), 1, 1);
  float *piano;

  Bimage *ppol = new Bimage (Float, TSimple, prad->sizeX(), prad->sizeY(), 1, 1);

  char order[12] = "yxz";
  for (k = 0; k < prad->sizeZ(); k++)
    {
      pr = psource + p->sizeX() * p->sizeY() * k;
      piano = (float *) ppiano->data_pointer();

      if (padd)
	padding (p, ppiano, k);
      else
	for (pr1 = pr, pk1 = piano; pr1 < pr + p->sizeX() * p->sizeY(); pr1++, pk1++)
	  *pk1 = *pr1;

      if (tipo != 4)
	ppiano->sizeY(ppiano->sizeY() / 2);
 		ppiano->reslice(order);

      fft_1D_forward (ppiano);

      shift_phase (ppiano, (float) (ppiano->sizeY() / 2));

      memset ((char *) pfr, 0, 2 * pf_x * pf_y * sizeof (float));

      copy180 (ppiano, pfr, pfi, bordo);

      n_theta = ppiano->sizeX();
      n_radii = ppiano->sizeY() / 2;

 //     swap(ppiano->x, ppiano->y);
	  t = ppiano->sizeX();
	  ppiano->sizeX(ppiano->sizeY());
	  ppiano->sizeY(t);

      memset (ppol->data_pointer(), 0, prad->sizeX() * prad->sizeY() * sizeof (float));
      polcar_f (ppol, pfr, pfi, bordo, n_radii, n_theta, ker, padd);

      rephase_orig (ppol);

      ritoab (ppol);

      zeroes (ppol);

      fft_2D_backward (ppol);

      pr = (float *) ppol->data_pointer();

      if (padd)
	for (pr1 = pr; pr1 < pr + prad->sizeX() * prad->sizeY(); pr1++)
	  *pr1 *= 2. * padd;

      memcpy (pdest + prad->sizeX() * prad->sizeY() * k, pr, prad->sizeX() * prad->sizeY() * sizeof (float));

      if (tipo != 4)
	{
	  ppiano->sizeY(ppiano->sizeY() * 2);
	  ppiano->data_assign((unsigned char *)realloc (ppiano->data_pointer(), ppiano->sizeX() * ppiano->sizeY() * sizeof (float)));
	}
    }

  delete ppol;
  delete ppiano;
  free ((float *) pfr);

  return (prad);
}

/**

@param	*p	
@param	*pfr
@param	*pfi
@param	bordo
@param	n_radii
@param	n_theta
@param	*ker
@param	padd
@return int     0
**/
int
polcar_f (Bimage * p, float *pfr, float *pfi, int bordo, int n_radii, int n_theta, kernel * ker, int padd)
{
  int i, j, hx, hy;
  float r, r2, r2max, phi, vr, vi, memo, round;
  float vr1, vi1;
  float *pr1;

  float *pr = (float *) p->data_pointer();

	hx = p->sizeX() / 2;
	hy = p->sizeY() / 2;
  round = 1. / (2 * ker->nkern);
  r2max = (hx - 1) * (hx - 1);
  for (i = -hy + 1; i < hy; i++)
    for (j = 0; j < hx; j++)
      {
	r2 = i * i + j * j;
	if (r2 <= r2max)
	  {
	    if (j != 0)
	      phi = atan2 ((double) i, (double) j);
	    else if (i >= 0)
	      phi = +M_PI / 2.0;
	    else
	      phi = -M_PI / 2.0;
	    r = sqrt (r2);
	    if (padd)
	      r *= 2 * padd;
	    r += round;
	    if (phi < 0)
	      phi += 2 * M_PI;
	    if (phi >= 2 * M_PI)
	      phi -= 2 * M_PI;
	    polcart_fourier (pfr, r, phi, &vr1, bordo, n_radii, n_theta, ker);
	    polcart_fourier (pfi, r, phi, &vi1, bordo, n_radii, n_theta, ker);
	    memo = 1;
	    if (phi >= M_PI)
	      memo = -1;
	    vr = vr1;
	    vi = vi1 * memo;
	    if (i >= 0)
	      {
		pr1 = pr + 2 * j * p->sizeX() + 2 * i;
		*pr1 = vr;
		pr1 += p->sizeX();
		*pr1 = vi;
		if (i == 0)
		  {
		    *++pr1 = vr;
		    pr1 -= p->sizeX();
		    *pr1 = vi;
		  }
	      }
	    else
	      {
		pr1 = pr + 2 * j * p->sizeX() - 2 * i + 1;
		*pr1 = vi;
		pr1 += p->sizeX();
		*pr1 = vr;
	      }
	  }
      }

  return 0;
}

/**

@param	*pfri
@param	r
@param	phi
@param	*val
@param	bordo
@param	n_radii
@param	n_theta
@param	*ker
@return int     0
**/
int
polcart_fourier (float *pfri, float r, float phi, float *val, int bordo, int n_radii, int n_theta, kernel * ker)
{
  int ncol = n_theta + 2 * bordo;
  float *itbx, *itby;
  double x, y, toty, totx, xyt, f512 = ker->nkern;
  int i, j, k, rnk = ker->nk, nkm1;
  float *pf1, *pf2;

  if (phi < 0)
    phi += 2 * M_PI;
  if (phi > 2 * M_PI)
    phi -= 2 * M_PI;
  if (phi >= M_PI)
    phi -= M_PI;
  x = phi * (float) n_theta / M_PI + bordo;
  y = r + bordo;
  x += 1. / (2. * ker->nkern);

  if (bordo == 1)
    {
      i = (int) floor (x);
      j = (int) floor (y);
      totx = x - i;
      toty = y - j;
      pf1 = pfri + i + j * ncol;
      *val = pf1[0] * (1 - totx) * (1 - toty) + pf1[1] * totx * (1 - toty) + pf1[ncol] * (1 - totx) * toty + pf1[ncol + 1] * totx * toty;
    }
  else
    {
      if (ker->nk % 2 == 0)
	nkm1 = ker->nk / 2 - 1;
      else
	nkm1 = (ker->nk - 1) / 2;
      xyt = floor (y);
      j = (int) floor (f512 * (y - xyt));
      j *= rnk;
      itby = ker->kern + j;
      i = (int) xyt;
      i -= nkm1;

      pf1 = pfri + i * ncol;

      xyt = floor (x);
      j = (int) floor (f512 * (x - xyt));
      j *= rnk;
      itbx = ker->kern + j;
      i = (int) xyt;
      i -= nkm1;
      pf1 = pf1 + i;

      pf2 = pf1;
      toty = 0;
      for (j = 0; j < rnk; j++)
	{
	  totx = 0;
	  for (k = 0; k < rnk; k++)
	    {
	      totx += *itbx * *pf2++;
	      itbx++;
	    }
	  toty += *itby * totx;
	  itby++;
	  itbx -= rnk;
	  pf2 += ncol;
	  pf2 -= rnk;
	}
      *val = toty;
    }
  return 0;
}


/**
@brief 	Padds an image
@param	*p		starting image parameters
@param	*ppiano	padded image
@param	z
@return int	0
**/
int
padding (Bimage * p, Bimage * ppiano, int z)
{
  float *pk1, *pr1, *pk2;
  int i, j, offset;

  float *pr = (float *) p->data_pointer();
  float *piano = (float *) ppiano->data_pointer();

  memset (piano, 0, ppiano->sizeX() * ppiano->sizeY() * sizeof (float));
  offset = (ppiano->sizeX() - p->sizeX()) / 2;
  pk1 = piano + offset;
  pr1 = pr + z * p->sizeX() * p->sizeY();
  for (i = 0; i < p->sizeY(); i++)
    {
      pk2 = pk1;
      for (j = 0; j < p->sizeX(); j++)
	*pk2++ = *pr1++;
      pk1 += ppiano->sizeX();
    }

  return 0;
}

/**

@param	*p	
@param	*pfr
@param	*pfi
@param	bordo
@return int	0
**/
int
copy180 (Bimage * p, float *pfr, float *pfi, int bordo)
{
  float *pr1, *pf1, *pr2, *pf2;
  int i, j, ncol, nrow;

  float *pr = (float *) p->data_pointer();

  ncol = p->sizeX() + 2 * bordo;
  nrow = p->sizeY() / 2 + 2 * bordo;

  pf1 = pfr + bordo + bordo * ncol;
  for (pr1 = pr; pr1 < pr + p->sizeX(); pr1 += 2 * p->sizeX())
    {
      pf2 = pf1;
      for (pr2 = pr1; pr2 < pr1 + p->sizeX(); pr2++)
	*pf2++ = *pr2;
      pf1 += ncol;
    }

  for (pr1 = pr + 2 * p->sizeX(); pr1 < pr + p->sizeX() * p->sizeY(); pr1 += 2 * p->sizeX())
    {
      pf2 = pf1;
      for (pr2 = pr1; pr2 < pr1 + p->sizeX(); pr2++)
	*pf2++ = *pr2;
      pf1 += ncol;
    }

  pf1 = pfi + bordo + bordo * ncol;
  for (pr1 = pr + p->sizeX(); pr1 < pr + 2 * p->sizeX(); pr1 += 2 * p->sizeX())
    {
      pf2 = pf1;
      for (pr2 = pr1; pr2 < pr1 + p->sizeX(); pr2++)
	*pf2++ = 0.;
      pf1 += ncol;
    }

  for (pr1 = pr + 3 * p->sizeX(); pr1 < pr + p->sizeX() * p->sizeY(); pr1 += 2 * p->sizeX())
    {
      pf2 = pf1;
      for (pr2 = pr1; pr2 < pr1 + p->sizeX(); pr2++)
	*pf2++ = *pr2;
      pf1 += ncol;
    }

  for (i = 1; i <= bordo; i++)
    {
      pf1 = pfr + bordo + (bordo - i) * ncol;
      pf2 = pfr + bordo + (bordo + i) * ncol;
      for (j = 0; j < p->sizeX(); j++)
	*pf1++ = *pf2++;
    }
  for (i = 0; i < bordo; i++)
    {
      pf1 = pfr + bordo + (bordo + p->sizeY() / 2 - 1 - i) * ncol;
      pf2 = pfr + bordo + (bordo + p->sizeY() / 2 + i) * ncol;
      for (j = 0; j < p->sizeX(); j++)
	*pf2++ = 0.;
    }
  for (j = 0; j < nrow; j++)
    {
      pf1 = pfr + j * ncol;
      pf2 = pfr + j * ncol + p->sizeX();
      for (i = 0; i < bordo; i++)
	*pf1++ = *pf2++;
    }
  for (j = 0; j < nrow; j++)
    {
      pf1 = pfr + j * ncol + p->sizeX() + bordo;
      pf2 = pfr + j * ncol + bordo;
      for (i = 0; i < bordo; i++)
	*pf1++ = *pf2++;
    }

  for (i = 1; i <= bordo; i++)
    {
      pf1 = pfi + bordo + (bordo - i) * ncol;
      pf2 = pfi + bordo + (bordo + i) * ncol;
      for (j = 0; j < p->sizeX(); j++)
	*pf1++ = -*pf2++;
    }
  for (i = 0; i < bordo; i++)
    {
      pf1 = pfi + bordo + (bordo + p->sizeY() / 2 - 1 - i) * ncol;
      pf2 = pfi + bordo + (bordo + p->sizeY() / 2 + i) * ncol;
      for (j = 0; j < p->sizeX(); j++)
	*pf2++ = 0.;
    }
  for (j = 0; j < nrow; j++)
    {
      pf1 = pfi + j * ncol;
      pf2 = pfi + j * ncol + p->sizeX();
      for (i = 0; i < bordo; i++)
	*pf1++ = -*pf2++;
    }
  for (j = 0; j < nrow; j++)
    {
      pf1 = pfi + j * ncol + p->sizeX() + bordo;
      pf2 = pfi + j * ncol + bordo;
      for (i = 0; i < bordo; i++)
	*pf1++ = -*pf2++;
    }
  return 0;
}

/**

@param	*p	
@param	npix
@return int		0
**/
int
shift_phase (Bimage * p, float npix)
{
  float *pr0, *pr1, *pi1, *ptot, *pcol;
  double a, cs, sn, cs1, sn1;
  double ft, ang;

  float *pr = (float *) p->data_pointer();

  ft = npix;
  ft *= M_PI / (p->sizeY() / 2);
  ang = ft;

  cs1 = cos (ang);
  sn1 = sin (ang);
  cs = cs1;
  sn = sn1;
  ptot = pr + p->sizeX() * p->sizeY();
  for (pr0 = pr + 2 * p->sizeX(); pr0 < ptot; pr0 += 2 * p->sizeX())
    {
      pcol = pr0 + p->sizeX();
      for (pr1 = pr0; pr1 < pcol; pr1++)
	{
	  pi1 = pr1 + p->sizeX();
	  a = *pr1 * cs - *pi1 * sn;
	  *pi1 = *pr1 * sn + *pi1 * cs;
	  *pr1 = a;
	}
      a = cs * cs1 - sn * sn1;
      sn = cs * sn1 + sn * cs1;
      cs = a;
    }
  return 0;
}
