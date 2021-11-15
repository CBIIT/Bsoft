/**
@file	sin_gal.cpp
@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 23 01 2006 (BH)
**/

#include "sin_gal.h"
#include "fft_tool.h"
#include "Bimage.h"
#include "utilities.h"

int copy_sin (Bimage * p);
int expand (Bimage * p);
int merge (Bimage * p, int padd);
int cart_pol_Fourier_rad (Bimage * p, int n_theta, int padd, int range, int pf_x, int pf_y, float *pfr, float *pfi, kernel * ker);
int separate_ri (Bimage * p, float *pfr, float *pfi, int bordo);
int rephase_1D (Bimage * p);
int map_shift (Bimage * p, float dx, float dy);

/**

@param	*p	
@param	tipo
@param	shift
@param	*ker
@param	padd
@param	n_theta
@return Bimage*
**/
Bimage *
sin_gal (Bimage * p, int tipo, Vector3<float> shift, kernel * ker, int padd, int n_theta)
{
  int i, j;
//  float	t;
  float *pr, *pfr, *pfi;
	
  p->change_type(Float);

  Bimage *prad = p->copy_header ();

  prad->sizeY(n_theta / 2);
  if (tipo != 180)
    prad->sizeY(n_theta);

//  float *pdest = (float *) malloc (prad->sizeX() * prad->sizeY() * prad->sizeZ() * sizeof (float));
//  prad->data_assign((unsigned char *) pdest);

  float *psource = (float *) p->data_pointer();
	
	float *pdest = (float *) prad->data_alloc();

//	cout << "prad size = " << prad->size() << endl;
	
  int pf_x = p->sizeX() / 2 + ker->bordo;
  int pf_y = p->sizeY();
  if (padd)
    {
      pf_x = padd * p->sizeX() + ker->bordo;
      pf_y = 2 * p->sizeY() * padd;
    }

  if ((pfr = (float *) calloc (2 * pf_x * pf_y, sizeof (float))) == 0)
    {
      perror ("sin_gal: pfr error calloc");
      exit (0);
    }
  pfi = pfr + pf_x * pf_y;

  Bimage *ppiano = new Bimage (Float, TSimple, p->sizeX() * (padd + 1), p->sizeY() * (padd + 1), 1, 1);
  float *piano = (float *) ppiano->data_pointer();

  Bimage *ppol = new Bimage (Float, TSimple, prad->sizeX(), prad->sizeY(), 1, 1);

//	cout << "ppiano size" << ppiano->size() << endl;
//	cout << "ppol size" << ppol->size() << endl;
//	cout << "starting sin_gal" << endl;
	
  char order[12] = "yxz";
  long		p_plane(p->sizeX() * p->sizeY());
  long		prad_plane(prad->sizeX() * prad->sizeY());
  long		piano_plane(ppiano->sizeX() * ppiano->sizeY());
	
  for (i = 0; i < p->sizeZ(); i++)
    {
		
 //     swap(ppol->x, ppol->y);
	  j = ppol->sizeX();
	  ppol->sizeX(ppol->sizeY());
	  ppol->sizeY(j);
      if (tipo == 360)
	ppol->sizeX(ppol->sizeX() / 2);

      pr = psource + p_plane * i;
//		cout << "slice " << i << " pointer = " << pr - psource << endl;

		// need to reset because the pointer changes
		piano = (float *) ppiano->data_pointer();

      for (j = 0; j < p_plane; j++)
			piano[j] = pr[j];
      for (j = p_plane; j < piano_plane; j++)
			piano[j] = 0.;

      if (padd)
			merge (ppiano, padd);

     fft_2D_forward (ppiano);

      zeroes (ppiano);

      abtori (ppiano);

      // shift 2D in FT //
      map_shift (ppiano, shift[0], shift[1]);

      // Shift origin
      rephase_orig (ppiano);

      // Transpose 2D
 		ppiano->reslice(order);

      separate_ri (ppiano, pfr, pfi, ker->bordo);

      cart_pol_Fourier_rad (ppol, n_theta / 2, padd, 180, pf_x, pf_y, pfr, pfi, ker);

      rephase_1D (ppol);

       fft_1D_backward (ppol);

 		ppol->reslice(order);
	  
      if (tipo == 360)
	{
	  ppol->sizeY(ppol->sizeY() * 2);
	  ppol->data_assign((unsigned char *)realloc (ppol->data_pointer(), ppol->sizeX() * ppol->sizeY() * sizeof (float)));
		copy_sin (ppol);
	}
	
      pr = (float *) ppol->data_pointer();
      memcpy (pdest + prad_plane * i, pr, prad_plane * sizeof (float));
    }

//	cout << "ending sin_gal" << endl;
	
  delete ppol;
  delete ppiano;
  free (pfr);

  return (prad);
}


/**

@param	*p	
@return Bimage *		
**/
Bimage *
copy_I_part_rad (Bimage * p)
{
  float *pr1, *pr2;
  int i, j, k, k2, piano, piano2;

  Bimage *p2 = new Bimage (Float, TSimple, p->sizeX() * 2, p->sizeY(), p->sizeZ(), 1);
  piano = p->sizeX() * p->sizeY();
  piano2 = p2->sizeX() * p2->sizeY();
  for (k = 0; k < p->sizeZ(); k++)
    {
      pr1 = (float *) p->data_pointer() + k * piano;
      pr2 = (float *) p2->data_pointer() + k * piano2;
      for (j = 0; j < p->sizeY(); j++, pr2 += p->sizeX())
	for (i = 0; i < p->sizeX(); i++)
	  *pr2++ = *pr1++;

      k2 = (p->sizeZ() - k) % p->sizeZ();
      pr1 = (float *) p->data_pointer() + k * piano;
      pr2 = (float *) p2->data_pointer() + k2 * piano2 + p->sizeX();
      for (j = 0; j < p->sizeY(); j++, pr2 += p->sizeX())
	for (i = 0; i < p->sizeX(); i++)
	  *pr2++ = *pr1++;
    }
  return (p2);
}

/**

@param	*p	
@return int			0
**/
int
copy_II_part_rad (Bimage * p)
{
  float *pr1, *pr2, *pr3, *pr4;
  int i, j, k, j1, j2;

  expand (p);

  float *pr = (float *) p->data_pointer();

  pr1 = pr;
  pr2 = pr + p->sizeX() * p->sizeY() / 2;
  for (k = 0; k < p->sizeY() / 2; k++)
    {
      j1 = p->sizeX() / 4;
      for (j = 0; j < p->sizeX(); j++)
	{
	  j1 = (p->sizeX() + p->sizeX() / 4 - j) % p->sizeX();
	  j2 = (p->sizeX() + p->sizeX() / 4 + j) % p->sizeX();
	  pr3 = pr1 + j1;
	  pr4 = pr2 + j2;
	  for (i = 0; i < p->sizeZ(); i++)
	    {
	      *pr4 = *pr3;
	      pr4 += p->sizeX() * p->sizeY();
	      pr3 += p->sizeX() * p->sizeY();
	    }
	}
      pr1 += p->sizeX();
      pr2 += p->sizeX();
    }
  return 0;
}


/**

@param	*p	
@return int			0
**/
int
expand (Bimage * p)
{
  float *pr1, *pr2;
  unsigned long j, k;

  unsigned long datasize = p->sizeX() * p->sizeY() * p->sizeZ() * sizeof (float);
  unsigned long planesize = p->sizeX() * p->sizeY();
  float *pv = (float *) p->data_pointer();
  float *pr = (float *) malloc (datasize * 2);

  memcpy (pr, pv, datasize);
//  free (pv);

  for (k = p->sizeZ() - 1; k > 0; k--)
    {
      pr1 = pr + k * planesize;
      pr2 = pr + k * planesize * 2;
      for (j = 0; j < planesize; j++)
	*pr2++ = *pr1++;
    }

  p->data_assign((unsigned char *) pr);
  p->sizeY(p->sizeY()*2);
//  free(pv);

  return 0;
}

/**
@brief 	Pads an image and shifts it to the middle.
@param	*p	
@param	padd		0=no padding, 1=padding
@return int			0
**/
int
merge (Bimage * p, int padd)
{
  int i, j, offset, nr_old;
  float *pr1, *pr2;

  float *pr = (float *) p->data_pointer();

//		cout << "merge: pad=" << padd << endl;
	
   nr_old = p->sizeY() / (padd + 1);
  offset = (int) (p->sizeY() * 3.0 / 8);
  if (padd == 1)
    offset = p->sizeY() / 4;

  pr2 = pr;
  pr1 = pr + p->sizeX() * offset + offset;
  for (i = 0; i < nr_old; i++)
    {
      for (j = 0; j < nr_old; j++)
	*pr1++ = *pr2++;
      pr1 += 2 * offset;
    }
  j = p->sizeX() * p->sizeY() / 4;
  pr1 = pr;
  for (i = 0; i < j; i++)
    *pr1++ = 0.;

  return 0;
}

/**

@param	*p	
@param	n_theta
@param	padd
@param	range
@param	pf_x
@param	pf_y
@param	*pfr
@param	*pfi
@param	*ker
@return int			0
**/
int
cart_pol_Fourier_rad (Bimage * p, int n_theta, int padd, int range, int pf_x, int pf_y, float *pfr, float *pfi, kernel * ker)
{
  int i, j, k, jmin, jmax, kmin, kmax, row, col, nkm1, rnk, rowmax;
  float *pr1, *pf1, *pr2, *pf2, *itbx, *itby;
  float fat, x, y, r, theta, xyt, totx, toty, memo, weight, round;

  float *pr = (float *) p->data_pointer();

  round = 1. / (2 * ker->nkern);
  weight = 1.;
  if (padd == 1)
    weight = 4.;
  if (padd == 2)
    weight = 16.;
  rnk = ker->nk;		// Kernel size
  if (range == 360)
    fat = 2 * M_PI / n_theta;	// fat is simply radians/degrees
  else
    fat = M_PI / n_theta;
	nkm1 = ker->bordo;
//  if (ker->nk % 2 == 0)
//    nkm1 = ker->nk / 2 - 1;
//  else
//    nkm1 = (ker->nk - 1) / 2;
  rowmax = pf_y;
  if (padd == 1)
    rowmax /= 2;
  if (padd == 2)
    rowmax /= 4;
  for (row = 0; row < rowmax; row += 2)
    {
      pr1 = pr + row * n_theta;
      pr2 = pr1 + n_theta;
      r = (float) row *.5;
      if (padd)
			r *= 2*padd;
      for (col = 0; col < n_theta; col++)
	{
	  theta = col * fat;
	  memo = 1;
	  while (theta > M_PI)
	    theta -= 2 * M_PI;	// Range -PI to PI
	  while (theta < -M_PI)
	    theta += 2 * M_PI;
	  if (theta > M_PI / 2)	// Range -PI/2 to PI/2
	    {
	      theta -= M_PI;
	      memo = -1;
	    }
	  if (theta < -M_PI / 2)
	    {
	      theta += M_PI;
	      memo = -1;
	    }
	  x = r * cos (theta) + ker->bordo + round;
	  y = r * sin (theta) + pf_y / 2 + round;

	  xyt = floor (y);
	  j = (int) (ker->nkern * (y - xyt));	// Fraction * size of kernel ( = 512 )
	  j *= rnk;								// Starting pount in LUT
	  itby = ker->kern + j;					// Pointer to starting point in LUT
	  i = (int) xyt;
	  i -= nkm1;
	  j = i + rnk;
	  if (i < 0 || j > pf_y)
	    {
	      if (i < 0)
		{
		  jmin = -i;
		  jmax = rnk;
		  itby -= i;
		  i = 0;
		}
	      else
		{
		  jmin = 0;
		  jmax = pf_y - i;
		}
	    }
	  else
	    {
	      jmin = 0;
	      jmax = rnk;
	    }

	  pf1 = pfr + i * pf_x;
	  pf2 = pfi + i * pf_x;

	  xyt = floor (x);
	  j = (int) (ker->nkern * (x - xyt));
	  j *= rnk;
	  itbx = ker->kern + j;
	  i = (int) xyt;
	  i -= nkm1;
	  j = i + rnk;
	  if (i < 0 || j > pf_x)
	    {
	      if (i < 0)
		{
		  kmin = -i;
		  kmax = rnk;
		  itbx -= i;
		  i = 0;
		}
	      else
		{
		  kmin = 0;
		  kmax = pf_x - i;
		}
	    }
	  else
	    {
	      kmin = 0;
	      kmax = rnk;
	    }

	  pf1 = pf1 + i;
	  pf2 = pf2 + i;

	  toty = 0;
	  for (j = jmin; j < jmax; j++)
	    {
	      totx = 0.;
	      for (k = kmin; k < kmax; k++)
				totx += *itbx++ * *pf1++;
	      toty += *itby++ * totx;
	      itbx -= kmax - kmin;
	      pf1 += pf_x - (kmax - kmin);
	    }
	  *pr1++ = toty * weight;

	  itby -= jmax - jmin;
	  toty = 0;
	  for (j = jmin; j < jmax; j++)
	    {
	      totx = 0.;
	      for (k = kmin; k < kmax; k++)
				totx += *itbx++ * *pf2++;
	      toty += *itby++ * totx;
	      itbx -= kmax - kmin;
	      pf2 += pf_x - (kmax - kmin);
	    }
	  *pr2++ = toty * memo * weight;
	}
    }
  return 0;
}

/**

@param	*p	
@param	*pfr
@param	*pfi
@param	bordo
@return int			0
**/
int
separate_ri (Bimage * p, float *pfr, float *pfi, int bordo)
{
  int i, j, row, ncol, hx, hy;
  float *pr1, *pfr1, *pfr2, *pfi1, *pfi2;

  float *pr = (float *) p->data_pointer();

  hx = p->sizeX()/2;
  hy = p->sizeY()/2;
  ncol = hx + bordo;
  pfr1 = pfr + ncol + bordo;
  pfi1 = pfi + ncol + bordo;
  for (i = -hy + 1; i < hy; i++)
    {
      pfr2 = pfr1;
      row = abs (i);
      pr1 = pr + row * p->sizeX() * 2;
      if (i < 0)
	pr1 += p->sizeX() + 1;
      for (j = 0; j < hx; j++)
	{
	  *pfr2++ = *pr1;
	  pr1 += 2;
	}
      pfi2 = pfi1;
      pr1 = pr + row * p->sizeX() * 2 + 1;
      if (i < 0)
	pr1 += p->sizeX() - 1;
      for (j = 0; j < hx; j++)
	{
	  *pfi2++ = *pr1;
	  pr1 += 2;
	}
      pfr1 += ncol;
      pfi1 += ncol;
    }
  pfr1 = pfr + ncol + bordo + 1;
  pfr2 = pfr + ncol * (p->sizeY() - 1) + bordo - 1;
  for (i = -hy + 1; i < hy; i++)
    {
      for (j = 0; j < bordo; j++)
	*pfr2-- = *pfr1++;
      pfr1 += ncol - bordo;
      pfr2 -= ncol - bordo;
    }
  pfi1 = pfi + ncol + bordo + 1;
  pfi2 = pfi + ncol * (p->sizeY() - 1) + bordo - 1;
  for (i = -hy + 1; i < hy; i++)
    {
      for (j = 0; j < bordo; j++)
	*pfi2-- = -*pfi1++;
      pfi1 += ncol - bordo;
      pfi2 -= ncol - bordo;
    }
  return 0;
}

/**

@param	*p	
@return int			0
**/
int
rephase_1D (Bimage * p)
{
  float *pr1, *pr2;
  float meno;
  int ncol2;

  float *pr = (float *) p->data_pointer();

  ncol2 = 2 * p->sizeX();
  meno = 1.;
  for (pr1 = pr + p->sizeX(); pr1 < pr + ncol2; pr1++)
    *pr1 = 0.;
  for (pr1 = pr + ncol2; pr1 < pr + p->sizeX() * p->sizeY(); pr1 += ncol2)
    {
      meno *= -1.;
      for (pr2 = pr1; pr2 < pr1 + ncol2; pr2++)
	*pr2 *= meno;
    }

  return 0;
}

/**

@param	*p	
@param	dx
@param	dy
@return int			0
**/
int
map_shift (Bimage * p, float dx, float dy)
{
  if (dx == 0 && dy == 0)
    return 0;

  int j, k;
  float *Rp, *Rm, *Ip, *Im;
  float a1, a2, cs1, sn1, cs2, sn2, csa, sna;
  float vr, vi, val;

  float *pr = (float *) p->data_pointer();

  a1 = -dx * M_PI / (p->sizeY() / 2);
  a2 = -dy * M_PI / (p->sizeX() / 2);
  cs1 = cos (a1);
  sn1 = sin (a1);
  cs2 = cos (a2);
  sn2 = sin (a2);

  // shift on x //
  csa = 1.;
  sna = 0.;
  for (j = 0; j < p->sizeY(); j += 2)
    {
      Rp = pr + j * p->sizeX();
      Im = Rp + 1;
      Ip = Rp + p->sizeX();
      Rm = Ip + 1;
      for (k = 0; k < p->sizeX(); k += 2)
	{
	  vr = *Rp * csa - *Ip * sna;
	  vi = *Rp * sna + *Ip * csa;
	  *Rp = vr;
	  *Ip = vi;
	  vr = *Rm * csa - *Im * sna;
	  vi = *Rm * sna + *Im * csa;
	  *Rm = vr;
	  *Im = vi;
	  Rp += 2;
	  Ip += 2;
	  Rm += 2;
	  Im += 2;
	}
      val = csa * cs1 - sna * sn1;
      sna = csa * sn1 + sna * cs1;
      csa = val;
    }

  // shift on y //
  csa = 1.;
  sna = 0.;
  for (j = 0; j < p->sizeX(); j += 2)
    {
      Rp = pr + j;
      Im = Rp + 1;
      Ip = Rp + p->sizeX();
      Rm = Ip + 1;
      for (k = 0; k < p->sizeY(); k += 2)
	{
	  vr = *Rp * csa - *Ip * sna;
	  vi = *Rp * sna + *Ip * csa;
	  *Rp = vr;
	  *Ip = vi;
	  vr = *Rm * csa + *Im * sna;
	  vi = -*Rm * sna + *Im * csa;
	  *Rm = vr;
	  *Im = vi;
	  Rp += 2 * p->sizeX();
	  Ip += 2 * p->sizeX();
	  Rm += 2 * p->sizeX();
	  Im += 2 * p->sizeX();
	}
      val = csa * cs2 - sna * sn2;
      sna = csa * sn2 + sna * cs2;
      csa = val;
    }
  return 0;
}
