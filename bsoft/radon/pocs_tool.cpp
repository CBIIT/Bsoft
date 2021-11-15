/**
@file	pocs_tool.cpp
@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "pocs_tool.h"
#include "fft_tool.h"

int filter_sin_mask (Bimage * p, int n_cicli, float r_max, Bimage * pmask, int piano, int tipo, int first);
int restore_with_mask (Bimage * p, Bimage * pcopia, Bimage * pmask, int linea, int tipo);
int copy_sin (Bimage * p);
int copy_phi (Bimage * pelab, Bimage * p, int i);
int filter_bessel (Bimage * p, int r, int first);


/**

@param	*p			Radon transform
@param	n_cicli		number of cicles to do
@param	r_max		maximun ray of the object
@param	*pmask		mask
@param	piano		number of plane processed
@param	tipo		can be 0 or 1
@param	first
@return int     	0
**/
int
filter_sin_mask (Bimage * p, int n_cicli, float r_max, Bimage * pmask, int piano, int tipo, int first)
{
  int j;

  Bimage *pcopia = p->copy();

  for (j = 0; j < n_cicli; j++)
    {
      fft_2D_forward (p);

      filter_bessel (p, (int) r_max, first);
      zeroes (p);

      fft_2D_backward (p);

      if (pmask)
	restore_with_mask (p, pcopia, pmask, piano, tipo);
    }

  delete pcopia;

  return 0;
}

/**

@param	*p			Radon transform
@param	*pmask		mask
@param	n_cicli		number of cicles to do
@param	r_maxx		maximun ray of the object
@param	passata		can be 0 or 1
@return int     	0
**/
int
filter_X_mask (Bimage * p, Bimage * pmask, int n_cicli, float r_maxx, int passata)
{
  int i, first = 2;
  int planesize = p->sizeX() * p->sizeY();
  float *pp;
  float fat, r_max, theta;

  float *pv = (float *) p->data_pointer();

  Bimage *pelab = new Bimage (Float, TSimple, p->sizeX(), 2 * p->sizeY(), 1, 1);
	
  if (passata == 0)
    for (i = 0; i < p->sizeZ(); i++)
      {
	r_max = r_maxx;
	pp = pv + i * planesize;
	memcpy (pelab->data_pointer(), pp, planesize * sizeof (float));
	copy_sin (pelab);
	filter_sin_mask (pelab, n_cicli, r_max, pmask, i, 0, first);
	memcpy (pp, pelab->data_pointer(), planesize * sizeof (float));
      }
  else
    for (i = 0; i < p->sizeZ(); i++)
      {
	theta = M_PI * (i - (int)p->sizeZ() / 2) / p->sizeZ();
	fat = fabs (sin (theta));
	r_max = r_maxx * fat;
	pp = pv + i * planesize;
	memcpy (pelab->data_pointer(), pp, planesize * sizeof (float));
	copy_phi (pelab, p, i);
	filter_sin_mask (pelab, n_cicli, r_max, pmask, i, 1, first);
	memcpy (pp, pelab->data_pointer(), planesize * sizeof (float));
      }

  delete pelab;

  return 0;
}

/**

@param	*p			Radon transform
@param	*pcopia		temporary array
@param 	*pmask		mask
@param	linea
@param	tipo
@return int     	0
**/
int
restore_with_mask (Bimage * p, Bimage * pcopia, Bimage * pmask, int linea, int tipo)
{
  float *pm1;

  float *pr1, *pr2;
  int i, j, k, icol, val;

  float *copia = (float *) pcopia->data_pointer();
  float *elab = (float *) p->data_pointer();
  float *mask = NULL;
  if (pmask)
    mask = (float *) pmask->data_pointer();

  if (tipo == 0)
    {
      pm1 = mask + linea * p->sizeX();
      pr1 = copia;
      pr2 = elab;
      for (i = 0; i < p->sizeY(); i++)
	{
	  k = i % (p->sizeY() / 2);
	  val = (int) (pm1[k]);
	  for (j = 0; j < p->sizeX(); j++)
	    {
	      if (val != 0)
		*pr2 = *pr1;
	      pr1++;
	      pr2++;
	    }
	}
      pr1 = copia;
      pr2 = elab;
    }
  else
    {
      pm1 = mask + linea;
      pr1 = copia;
      pr2 = elab;
      for (i = 0; i < p->sizeY(); i++)
	{
	  k = i % (p->sizeY() / 2);
	  icol = k * p->sizeX();
	  val = (int) (pm1[icol]);
	  for (j = 0; j < p->sizeX(); j++)
	    {
	      if (val != 0)
		*pr2 = *pr1;
	      pr1++;
	      pr2++;
	    }
	}
      pr1 = copia;
      pr2 = elab;
    }
  return 0;
}

/**
@brief 	Copy the second half of a sinogram from the first one
@param	*p			sinogram
@return int     	0
**/
int
copy_sin (Bimage * p)
{
  float *pr1, *pr2, *pr3, *pr4;
  int j, k;

  float *pr = (float *) p->data_pointer();

  pr1 = pr;
  pr2 = pr + p->sizeX() * p->sizeY() / 2;

  for (k = 0; k < p->sizeY() / 2; k++)
    {
      for (j = 0; j < p->sizeX(); j++)
	{
	  pr3 = pr1 + (p->sizeX() - j) % p->sizeX();
	  pr4 = pr2 + j;
	  *pr4 = *pr3;
	}
      pr1 += p->sizeX();
      pr2 += p->sizeX();
    }
  return 0;
}

/**

@param	*pelab
@param	*p	
@param	i
@return int     		0
**/
int
copy_phi (Bimage * pelab, Bimage * p, int i)
{
  int j, nt2;
  float *pr1, *pp, *pp1;

  float *pr = (float *) pelab->data_pointer();
  float *pv = (float *) p->data_pointer();

  nt2 = pelab->sizeX() * pelab->sizeY() / 2;
  j = (p->sizeZ() - i) % p->sizeZ();
  pp = pv + j * nt2;
  pr1 = pr + nt2;
  for (pp1 = pp; pp1 < pp + nt2; pp1++)
    *pr1++ = *pp1;
  return 0;
}

/**

@param	*p	
@param	r
@param	first
@return int     		0
**/
int
filter_bessel (Bimage * p, int r, int first)
{
  float *pr1, *pr2;
  int i, j;
  float m, rmax, rmaxx;
  int x, ymax, imax;

  float *pr = (float *) p->data_pointer();

  rmax = (float) r / (float) p->sizeY();
  imax = first * 2;

  rmaxx = (float) (p->sizeX() / 2 - 1) / (2. * M_PI * rmax);
  m = rmaxx / (float) (p->sizeX() / 2 - first);
  pr[1] = 0.;
  pr[p->sizeX()] = 0.;
  pr[p->sizeX() + 1] = 0.;
  for (i = 2; i < imax; i += 2)
    {
      pr1 = pr + i;
      pr2 = pr1 + p->sizeX();
      *pr1++ = 0.;
      *pr1 = 0.;
      *pr2++ = 0.;
      *pr2 = 0.;
    }
  for (i = imax; i < p->sizeX(); i += 2)
    {
      x = (i / 2 - first);
      ymax = (int) ceil (m * (float) x);
      if (ymax == 0)
	ymax = 1;
      if (ymax > (p->sizeY() / 2))
	ymax = p->sizeY() / 2;
      for (j = 0; j < ymax; j++)
	{
	  pr1 = pr + j * 2 * p->sizeX() + i;
	  pr2 = pr1 + p->sizeX();
	  *pr1++ = 0.;
	  *pr1 = 0.;
	  *pr2++ = 0.;
	  *pr2 = 0.;
	}
    }
  return 0;
}

/**

@param	*p	
@param	rmax
@return int     		0
**/
int
support_f (Bimage * p, float rmax)
{
  int i, j, j1, j2;
  float weight;
  float *pr1;

  float *pr = (float *) p->data_pointer();

//  cout << "in support_f" << endl;
	
  j1 = p->sizeZ() / 2 - (int) rmax;
  j2 = p->sizeZ() / 2 + (int) rmax;

  double ang_fac = M_PI / (2.0 * j1);

  for (j = 0; j < j1; j++)
    {
      pr1 = pr + j * p->sizeX() * p->sizeY();
      weight = sin (ang_fac * j);
      for (i = 0; i < p->sizeX() * p->sizeY(); i++)
	pr1[i] *= weight;
    }
  for (j = j2 + 1; j < p->sizeZ(); j++)
    {
      pr1 = pr + j * p->sizeX() * p->sizeY();
      weight = sin (ang_fac * (p->sizeZ() - j));
      for (i = 0; i < p->sizeX() * p->sizeY(); i++)
	pr1[i] *= weight;
    }
  return 0;
}
