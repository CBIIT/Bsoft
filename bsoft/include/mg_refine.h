/**
@file	mg_refine.h
@brief	Reciprocal space refinement of orientation parameters of particle images.
@author  Bernard Heymann
@date	Created: 20070115
@date	Modified: 20190207
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "symmetry.h"

// Function prototypes
long		mg_refine_orientations(Bproject* project, Bstring& reffile, Bstring& maskfile,
				Bstring& sym_string, int part_select, int max_iter, 
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, int fom_type, vector<double> weight,
				double hi_res, double lo_res, int kernel_width, int kernel_power, 
				double edge_radius, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags);
long		project_refine_orientations(Bproject* project, Bstring& reffile, Bstring& maskfile,
				Bstring& sym_string, int part_select, int max_iter, 
				double alpha_step, double accuracy, double shift_step, 
				double shift_accuracy, int fom_type, vector<double> weight,
				double hi_res, double lo_res, int kernel_width, int kernel_power,
				double edge_radius, double def_std, double shift_std, 
				double view_std, double max_angle, double max_mag, int flags);

