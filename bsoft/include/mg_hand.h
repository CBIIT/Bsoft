/**
@file	mg_hand.h
@brief	Author: David Belnap
@date	Created: 20011003
@date	Modified: 20050906
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
int  project_get_handedness(Bimage* penantiomer, Bproject* project, double* mg_ang, 
			int* mg_index, int* mg_select, double rad_min, double rad_max,
			double res_min, double res_max, double AmB_min, double AB_min, 
			int diff_out, int origins2, Bstring outimg);
int   hand_select_consist(Bproject* project, double* mg_ang, int* mg_index, int* mg_select, int sel_consist);

