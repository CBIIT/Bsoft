/**
@file	model_extract_build.h
@brief	Functions to extract subvolumes and build new maps.
@author Bernard Heymann
@date	Created: 20060411
@date	Modified: 20190821
**/

#include "rwmodel.h"
#include "rwimg.h"

// Function prototypes 
int			model_refine_components(Bmodel* model, Bstring* ct_names, Bimage* ptemp,
				Bimage* pmask, Bimage* pfsmask, int max_iter, double viewstep, double rotstep,
				double hires, double lores, double accuracy, double max_shift,
				double max_view_angle, double max_rot_angle, int shift_flag);
int			model_refine_link_positions(Bmodel* model, Bimage* ptemp, Bimage* pmask,
				Bimage* pfsmask, double hires, double lores, double max_shift, int shift_flag, double bias);
Bimage*		model_average_component_density(Bmodel* model, Vector3<long> size, Vector3<double> origin, int npt);
Bimage*		model_extract_component_densities(Bmodel* model, Vector3<long> size, Vector3<double> origin);
Bimage*		model_average_link_density(Bmodel* model, Vector3<long> size, Vector3<double> origin);
Bimage*		model_build_from_component_density(Bmodel* model, Vector3<long> size, 
				Vector3<double> origin, int flags);
Bimage*		model_build_from_link_density(Bmodel* model, Bstring& linkmap, 
				Vector3<long> size, Vector3<double> origin, int link_select, int flags);

