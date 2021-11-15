/**
@file	model_occupancy.h
@brief	Library routines to count components in maps.
@author Daniel Nemecek and Bernard Heymann
@date	Created: 20091202
@date	Modified: 20190201
**/


#include "rwimg.h"
#include "rwmodel.h"

// Function prototypes
int			model_occupancy(Bmodel* model, Bimage* pmask, double mol_weight, double rho, double cutoff, int invert_flag);
vector<double>	model_occupancy_distribution(Bmodel* model, double cutoff, int nfit,
				long& ncomp, vector<double>& prob, double& R);
int         model_refine_comp_for_occupancy(Bmodel* model, Bimage* pmask2, Bimage* ptemp, Bimage* pmask,
                           double hires, double lores, double max_shift);
