/**
@file	fft_tool.cpp
@author	Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "fft_tool.h"
#include "Bimage.h"

#include "utilities.h"

int
fft_1D_forward (Bimage * p)
{
//	cout << "forward 1D fft" << endl;
	
  scramrig (p);
  farfarig (p, 1);
  mixrig (p, 1);

  return 0;
}

int
fft_1D_backward (Bimage * p)
{
//	cout << "backward 1D fft" << endl;
	
  mixrig (p, -1);
  scramrig (p);
  farfarig (p, -1);

  return 0;
}

int
fft_2D_forward (Bimage * p)
{
  char order[12] = "yxz";

//	cout << "forward 2D fft" << endl;
	
   fft_1D_forward (p);
	p->reslice(order);
  fft_1D_forward (p);

  return 0;
}

int
fft_2D_backward (Bimage * p)
{
  char order[12] = "yxz";

//	cout << "backward 2D fft" << endl;
	
  fft_1D_backward (p);
	p->reslice(order);
  fft_1D_backward (p);

  return 0;
}

/**
@brief 	With MIXRIG and SCRAMRIG fast Fourier transform a 2D image
@param	*p		image
@param	ifl		direction of transformation (-1 = direct / 1 = inverse)
@return int		0
**/
int
farfarig (Bimage * p, int ifl)
{
  double a, b, c, e, dp = 1., q, uno = 1.;
  int j, i, r, u = 2, rifl;
  float *ppr, *ppi, *pqr, *pqi, *plim;
  int rncol, rnrow, ripw, ipwy;

  float *pr = (float *) p->data_pointer();

  for (i = 2, ipwy = 0; i < p->sizeY(); ipwy++, i *= 2);

  rncol = p->sizeX();
  rnrow = p->sizeY();
  ripw = ipwy;
  rifl = ifl;

  for (i = 0; i < ripw; i++)
    {
      c = 1;
      e = 0;
      q = rifl * sqrt ((uno - dp) * 0.5);
      dp = sqrt ((uno + dp) * 0.5);
      if (i == 0)
	dp = (-dp);
      for (r = 0; r < u; r += 2)
	{
	  for (j = r; j < rnrow; j += u * 2)
	    {
	      ppr = pr + j * rncol;
	      ppi = ppr + rncol;
	      plim = ppr + u * rncol;
	      pqi = plim + rncol;
	      for (pqr = plim; pqr < plim + rncol; pqr++)
		{
		  a = c ** pqr + e ** pqi;
		  b = e ** pqr - c ** pqi;
		  *pqr = *ppr - a;
		  *ppr += a;
		  *pqi = *ppi + b;
		  *ppi -= b;
		  ppr++;
		  ppi++;
		  pqi++;
		}
	    }
	  a = e * dp + c * q;
	  c = c * dp - e * q;
	  e = a;
	}
      u += u;
    }
  return 0;
}

/**
@brief 	With FARFARIG and SCRAMRIG fast Fourier transform a 2D image
@param	*p		image
@param	ifl		direction of transformation (0 = direct / 1 = inverse)
@return int		0
**/
int
mixrig (Bimage * p, int ifl)
{
  double a, b, u, v, fat;
  int j;
  float *prj, *pij, *prk, *pik, *plim;
  double dp, q, c, e, ad1, ad, ad2;
  int rncol, rnr2, rnrow;

  float *pr = (float *) p->data_pointer();

  rncol = p->sizeX();
  rnrow = p->sizeY();
  rnr2 = p->sizeY() / 2;
  ad2 = rnrow;
  ad = rnrow;
  if (ifl == -1)
    {
      ad = 0.5;
      ad2 = 1;
    }
  ad1 = 1 / (ad * 2);
  ad = 1 / ad;
  ad2 = 1 / ad2;
  dp = cos (M_PI / rnr2);
  q = sin (M_PI / rnr2) * ifl;

  fat = ad2;
  prj = pr;
  pij = pr + rncol;
  plim = pij;
  while (prj < plim)
    {
      a = *prj;
      b = *pij;
      *prj++ = (a + b) * fat;
      *pij++ = (a - b) * fat;
    }

  c = ifl;
  e = 0;
  fat = ad1;
  for (j = 2; j < rnr2; j += 2)
    {
      a = e * dp + c * q;
      c = c * dp - e * q;
      e = a;
      prk = pr + (rnrow - j) * rncol;
      pik = prk + rncol;
      prj = pr + j * rncol;
      pij = prj + rncol;
      plim = pik;
      while (prk < plim)
	{
	  a = *prj + *prk;
	  b = (*pij + *pik) * c - (*prj - *prk) * e;

	  u = *pij - *pik;
	  v = (*pij + *pik) * e + (*prj - *prk) * c;

	  *prj++ = (a + b) * fat;
	  *pij++ = (u - v) * fat;
	  *prk++ = (a - b) * fat;
	  *pik++ = -(u + v) * fat;
	}
    }

  fat = ad;
  prk = pr + rnr2 * rncol;
  pik = prk + rncol;
  plim = pik;
  while (prk < plim)
    {
      *pik *= -fat;
      *prk *= fat;
      pik++;
      prk++;
    }
  return 0;
}

/**
@brief 	With FARFARIG and MIXRIG fast Fourier transform a 2D image
@param	*p		image
@return int		0

	Rearrange each row in the complex image.

**/
int
scramrig (Bimage * p)
{
  int j, i, k = 0, a, ncol2 = 2 * p->sizeX();
  int *prj, *prk, *pij, *pik, *plim;
  int rncol = p->sizeX(), nr2 = p->sizeY() / 2;

  float *pr = (float *) p->data_pointer();

  for (j = 1; j < nr2; j++)
    {
      i = 2;
      while (k >= nr2 / i)
	{
	  k -= nr2 / i;
	  i += i;
	}
      k += nr2 / i;
      if (k > j)
	{
	  prj = (int *) pr + j * ncol2;
	  pij = prj + rncol;
	  prk = (int *) pr + k * ncol2;
	  pik = prk + rncol;
	  plim = pij;
	  while (prj < plim)
	    {
	      a = *prj;
	      *prj++ = *prk;
	      *prk++ = a;
	      a = *pij;
	      *pij++ = *pik;
	      *pik++ = a;
	    }
	}
    }
  return 0;
}

/**
@brief 	Changes direct Fourier transform representation
@param	*p		image
@return int		0

	This subroutine changes the way in which a direct transform 
	is represented; from the conventional form Txryr, Txiyr,
	Txryi, Txiyi (here called a,b,c,d ) it produces a new 
	representation in which real and imaginary parts
	are calculated, for the frequencies of the type f(Kx,Ky)
	and f(Kx,-Ky) .      
	Coefficients are addressed by pointers in this way:
	pp1 --> real f(Kx,Ky)       pp2 --> imaginary f(Kx,Ky)
	pp3 --> imaginary f(Kx,-Ky) pp4 --> real f(Kx,-Ky)

**/
int
abtori (Bimage * p)
{
  float *p1, *pp1, *pp2, *pp3, *pp4, *ptot;
  float a, c;
  int ncol2;

  float *pr = (float *) p->data_pointer();

  ptot = pr + p->sizeX() * p->sizeY();
  ncol2 = 2 * p->sizeX();
  pp2 = pr + p->sizeX();
  pp3 = pr + 1;
  pp4 = pr + 1 + p->sizeX();
  for (p1 = pr; p1 < ptot; p1 += ncol2)
    {
      for (pp1 = p1; pp1 < p1 + p->sizeX(); pp1 += 2)
	{
	  a = *pp1;
	  *pp1 = *pp1 - *pp4;
	  *pp4 += a;
	  c = *pp2 - *pp3;
	  *pp2 += *pp3;
	  *pp3 = c;
	  pp2 += 2;
	  pp3 += 2;
	  pp4 += 2;
	}
      pp2 += p->sizeX();
      pp3 += p->sizeX();
      pp4 += p->sizeX();
    }
  return 0;
}

/**
@brief 	Changes direct Fourier transform representation
@param	p		image
@return int		0

	This subroutine changes the representation of the
	direct transform: from real & imaginary parts to the 
	conventional one; here a,b,c,d ,pointed by pp1,pp2,
	pp3 and pp4, have the following meaning:
	a=TXrYr    b=TXiYi    c=TXrYi    d=TXiYr
	where TX(Y)r(i) means real (imaginary) part of the coefficients 
	obtained by a direct transform computed along the X (Y) direction.
	A quartet a,b,c,d is defined for each point of the reciprocal
	space in the first quadrant (i.e. points with positive Kx & Ky
	coordinates.     

**/
int
ritoab (Bimage * p)
{
  float *p1, *pp1, *pp2, *pp3, *pp4;
  float rp, im;

  float *pr = (float *) p->data_pointer();

  pp2 = pr + p->sizeX();
  pp3 = pr + 1;
  pp4 = pr + 1 + p->sizeX();
  for (p1 = pr; p1 < pr + p->sizeX() * p->sizeY(); p1 += 2 * p->sizeX())
    {
      for (pp1 = p1; pp1 < p1 + p->sizeX(); pp1 += 2)
	{
	  rp = *pp1;
	  *pp1 += *pp4;
	  *pp1 /= 2;
	  *pp4 -= rp;
	  *pp4 /= 2;
	  im = *pp3;
	  *pp3 = *pp2 - im;
	  *pp3 /= 2;
	  *pp2 += im;
	  *pp2 /= 2;
	  pp2 += 2;
	  pp3 += 2;
	  pp4 += 2;
	}
      pp2 += p->sizeX();
      pp3 += p->sizeX();
      pp4 += p->sizeX();
    }

  return 0;
}

/**
@param	*p		image
@return int		0
**/
int
zeroes (Bimage * p)
{
  float *pb1, *pc1, *pd1;

  float *pr = (float *) p->data_pointer();

  pd1 = pr + p->sizeX();
  *pd1 = 0.;			// origin  Tri=Tir=Tii=0 //
  for (pd1 = pr + p->sizeX(); pd1 < pr + 2 * p->sizeX(); pd1 += 2)	// axis Ky  Tri=0 Tii=0   //
    {
      pb1 = pd1 + 1;
      *pd1 = 0.;
      *pb1 = 0.;
    }
  for (pc1 = pr + 1; pc1 < pr + p->sizeX() * p->sizeY(); pc1 += 2 * p->sizeX())	// axis Kx  Tir=Tii=0 //
    {
      pb1 = pc1 + p->sizeX();
      *pc1 = 0.;
      *pb1 = 0.;
    }
  return 0;
}

/**
@brief 	Shifting the origin by half the complex image size
@param	*p		image
@return int		0

	Coefficients are addressed by pointers in this way:
	pp1 --> real f(Kx,Ky)      pp2 --> imaginary f(Kx,-Ky)
	pp3 --> imaginary f(Kx,Ky) pp4 --> real f(Kx,-Ky)

**/
int
rephase_orig (Bimage * p)
{
  float *p1, *pp1, *pp2, *pp3, *pp4;
  float cs, csx, csy, meno;

  float *pr = (float *) p->data_pointer();

  meno = -1.;
  csx = -1.;
  pp2 = pr + p->sizeX();
  pp3 = pr + 1;
  pp4 = pr + 1 + p->sizeX();
  for (p1 = pr; p1 < pr + p->sizeX() * p->sizeY(); p1 += 2 * p->sizeX())
    {
      csy = -1.;
      csx *= meno;
      for (pp1 = p1; pp1 < p1 + p->sizeX(); pp1 += 2)
	{
	  csy *= meno;
	  cs = csy * csx;
	  *pp1 *= cs;
	  *pp3 *= cs;
	  *pp4 *= cs;
	  *pp2 *= cs;
	  pp2 += 2;
	  pp3 += 2;
	  pp4 += 2;
	}
      pp2 += p->sizeX();
      pp3 += p->sizeX();
      pp4 += p->sizeX();
    }
  return 0;
}
