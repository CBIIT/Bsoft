/**
@file	proj_tool.cpp
@author P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "proj_tool.h"
#include "radon_util.h"
#include "utilities.h"

// Declaration of global variables
extern int verbose;		// Level of output to the screen

/**

@param	ncol
@param	sizeT
@return int			0
**/
LUTable *
create_table (int ncol, int sizeT)
{
  int i;

  LUTable *lut = (LUTable *) malloc (sizeof (LUTable));
  lut->ncol = ncol;
  lut->sizeT = sizeT;

  int size_tab_1 = (ncol * sizeT + 1) * (sizeT + 1);
  int size_tab_2 = size_tab_1 * ncol / 2;

  int **tab1;
  float **tab2;

  if ((tab1 = (int **) calloc (size_tab_1, sizeof (int *))) == 0)
    {
      perror ("create_table: tab1 error calloc");
      exit (0);
    }

  for (i = 0; i < size_tab_1; i++)
    if ((tab1[i] = (int *) calloc (2, sizeof (int))) == 0)
      {
	perror ("create_table: tab1 error calloc");
	exit (0);
      }

  if ((tab2 = (float **) calloc (size_tab_2, sizeof (float *))) == 0)
    {
      perror ("create_table: tab2 error calloc");
      exit (0);
    }

  for (i = 0; i < size_tab_2; i++)
    if ((tab2[i] = (float *) calloc (3, sizeof (float))) == 0)
      {
	perror ("create_table: tab1 error calloc");
	exit (0);
      }

  lut->tab1 = tab1;
  lut->tab2 = tab2;

  double fat, a, b, a1, b1, a2, b2, d;
  double alpha, beta, ang, arg;
  float w2, weight, stepb;
  int n, n1, n2, n_theta, n_theta_2, n_theta_4, step;

  w2 = 9.2023;

  n_theta = 2 * ncol;
  n_theta_2 = n_theta / 2;
  n_theta_4 = n_theta / 4;
  step = sizeT;
  fat = 2. * M_PI / n_theta;
//  delta = 1. / fat;
//  soglia = fat * 180. / M_PI;

  n1 = n2 = 0;
  stepb = 1. / (float) step;
  for (beta = -(float) n_theta_4; beta <= (float) n_theta_4 + .001; beta += stepb)
    {
      for (alpha = 0; alpha <= 1.001; alpha += stepb)
	{
	  b1 = (int) floor (beta);
	  b2 = b1 + 1;
	  n = 0;
	  for (b = b1; b <= b2; b += 1)
	    {
	      if (b != 0 && b != n_theta_2)
		{
		  d = 1.0 / sin (b * fat);
		  d = fabs (d);
		  a1 = (int) floor (alpha - d);
		  a2 = (int) floor (alpha + d) + 1;
		}
	      else
		{
		  a1 = 0;
		  a2 = n_theta - 1;
		}
	      for (a = a1; a <= a2; a += 1.0)
		{
		  ang = ang_one_two (alpha * fat, beta * fat, a * fat, b * fat);
		  if (ang < fat - 0.001)
		    {
		      arg = ang * n_theta / 2.;
		      if (arg >= 0.001)
			weight = (sin (arg) / arg) * exp (-(arg * arg) / w2);
		      else
			weight = 1;
		      ang *= 180. / M_PI;
		      tab2[n2][0] = a;
		      tab2[n2][1] = b + n_theta_4;
		      tab2[n2][2] = weight;
		      n2++;
		      n++;
		    }
		}
	    }
	  tab1[n1][0] = n;
	  tab1[n1][1] = n2;
	  n1++;
	}
    }

  if (verbose & VERB_FULL)
    cout << "LUTable sizes:                  " << n1 << " (" << size_tab_1 <<
		") " << n2 << " (" << size_tab_2 << ")" << endl;

  return (lut);
}

int
kill_table (LUTable * lut)
{
  int i;
  int size_tab_1 = (lut->ncol * lut->sizeT + 1) * (lut->sizeT + 1);
  int size_tab_2 = size_tab_1 * lut->ncol / 2;

  int **tab1 = lut->tab1;
  float **tab2 = lut->tab2;

  for (i = 0; i < size_tab_1; i++)
    free (tab1[i]);

  free (tab1);

  for (i = 0; i < size_tab_2; i++)
    free (tab2[i]);

  free (tab2);

  free (lut);

  return 0;
}


/**

@param	*p	
@param	*prec
@param	angolo		angle index
@param	csi
@param	eta
@param	*pmask
@param 	*lut
@param	*priga
@return int			0
**/
int
write_line (Bimage * p, Bimage * prec, int angolo, float csi, float eta, Bimage * pmask, LUTable * lut, float *priga)
{
  int i, j, k, l, memo, memo0;
  float *pg, *pg1, *pg2, *pv1, *pr1;
  int eta_in_pix, csi_in_pix;
  int place, nlines, ind_eta, ind_csi, ind;
  float eta0, eta_in_fl, csi_in_fl, weight;

  int sizeT = lut->sizeT;
  int **ptab1 = lut->tab1;
  float **ptab2 = lut->tab2;

  float *pr = (float *) p->data_pointer();
  float *pv = (float *) prec->data_pointer();
  float *pw1 = (float *) pmask->data_pointer();

  pg = pr;
  pr1 = priga;

  pg1 = pg + angolo;
  for (l = 0; l < p->sizeY(); l++)
    {
      pg2 = pg1 + l * p->sizeX();
      pr1[l] = *pg2;
    }

  memo0 = 0;
  if (csi < 0)
    {
      csi = -csi;
      eta -= M_PI;
      memo0 = 1;
    }
  while (eta < 0)
    eta += 2. * M_PI;

  eta_in_fl = eta * (float) p->sizeX() / M_PI;
  csi_in_fl = csi * (float) p->sizeX() / M_PI;
  eta0 = floor (eta_in_fl);

  ind_eta = (int) rint ((eta_in_fl - eta0) * (float) sizeT);
  ind_csi = (int) rint (csi_in_fl * (float) sizeT);
  ind = ind_eta + ind_csi * (sizeT + 1);

  nlines = ptab1[ind][0];
  place = ptab1[ind][1];

  for (j = 0; j < nlines; j++)
    {
      k = place + j;
      eta_in_pix = (int) ptab2[k][0] + (int) eta0;
      csi_in_pix = (int) ptab2[k][1];
      weight = ptab2[k][2];

      memo = memo0;
      while (eta_in_pix < 0)
	eta_in_pix += 2 * p->sizeX();
      if (csi_in_pix < 0)
	{
	  csi_in_pix *= -1;
	  eta_in_pix += p->sizeX();
	  memo += 1;
	}
      if (eta_in_pix >= 2 * p->sizeX())
	eta_in_pix -= 2 * p->sizeX();
      if (eta_in_pix >= p->sizeX() && csi_in_pix != 0)
	{
	  eta_in_pix -= p->sizeX();
	  csi_in_pix = (p->sizeX() - csi_in_pix) % p->sizeX();
	}
      if (eta_in_pix >= p->sizeX() && csi_in_pix == 0)
	{
	  eta_in_pix -= p->sizeX();
	  csi_in_pix = (p->sizeX() - csi_in_pix) % p->sizeX();
	  memo += 1;
	}

      memo = memo % 2;
      pw1[eta_in_pix * p->sizeX() + csi_in_pix] += weight;

      pv1 = pv + eta_in_pix * p->sizeX() + csi_in_pix;
      pv1[0] += pr1[0] * weight;
      pv1 += p->sizeX() * p->sizeY();
      if (memo == 0)
	for (i = 1; i < p->sizeY(); i++)
	  {
	    *pv1 += pr1[i] * weight;
	    pv1 += p->sizeX() * p->sizeY();
	  }
      else
	for (i = p->sizeY() - 1; i > 0; i--)
	  {
	    *pv1 += pr1[i] * weight;
	    pv1 += p->sizeX() * p->sizeY();
	  }
    }
  return 0;
}

/**
@brief 	Weighs a radon transform with a mask.
@param	*p			radon transform.
@param	*pmask		mask calculated for the radon transform.
@return int     	0

	The mask is searched for the minimum value to assess completeness.
	Each value in the mask above a threshold is inverted (the weight to be applied).
	Values below the threshold is set to zero.
	The threshold is set to 0.9.

**/
int
weigh_radon_transf (Bimage * p, Bimage * pmask)
{
  float *pw1, *pv1, *pv2, *pv3;
  int j, k;
  float soglia, val, min;

  float *pr = (float *) p->data_pointer();
  float *pw = (float *) pmask->data_pointer();

  min = 1000;
  soglia = 0.90;		// Threshold to accept a voxel
  pw1 = pw;
  for (j = 0; j < p->sizeX() * p->sizeY(); j++)
    {				// Finds minimum and sets values below threshold to zero
      val = *pw1;
      if (val < min)
	min = val;
      if (val > soglia)
	*pw1 = 1. / val;
      else
	*pw1 = 0.;
      pw1++;
    }

  if (verbose & VERB_PROCESS)
    {
      cout << "Minimum weight:                 " << min << endl;
      if (min > .01)
	cout << "Radon transform filled" << endl;
      else
	cout << "Radon transform unfilled" << endl;
    }

  for (k = 0; k < p->sizeZ(); k++)	// Weighs each z-plane with mask
    {
      pv1 = pr + p->sizeX() * p->sizeY() * k;
      pv2 = pv1 + p->sizeX() * p->sizeY();
      for (pv3 = pv1, pw1 = pw; pv3 < pv2; pv3++, pw1++)
	*pv3 *= *pw1;
    }

  for (j = 0, pw1 = pw; j < p->sizeX() * p->sizeY(); j++, pw1++)	// Sets the mask values to zero and one
    if (*pw1)
      *pw1 = 1;

  return 0;
}
