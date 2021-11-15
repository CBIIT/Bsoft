/**
@file	ps_model.h
@brief	Header file for postscript tools dealing with models
@author Bernard Heymann
@date	Created: 20090203
@date	Modified: 20190201
**/

#include "rwmodel.h"

//Function prototypes
int			ps_model_views(Bstring& filename, Bmodel* model, int combined);
int			ps_model_symmetry_views(Bstring& filename, Bmodel* model, string& symmetry_string, int combined);
int			ps_model_fom_histogram(Bstring& filename, Bmodel* model);
int			ps_model_occupancy(Bmodel* model, double cutoff, int bins, int nfit, 
				vector<double>& distrib, vector<double>& prob, double R, Bstring& filename);

