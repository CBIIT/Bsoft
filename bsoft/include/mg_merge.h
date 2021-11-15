/**
@file	mg_merge.h
@brief	Header file for functions to merge images.
@author David Belnap & Bernard Heymann
@date	Created: 20030410
@date	Modified: 20030901 (BH)
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "Bstring.h"


// Function prototypes
int			mg_particle_merge_series(Bproject* project, int mg_ref_select, int mg_index, 
				float mg_rot_ang, int mg_ori_select, Bstring outimg);
int			mg_particle_unmerge(Bproject* project, Bproject* orientations, float fom_diff);

