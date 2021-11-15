/**
@file	pocs_tool.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "rwimg.h"

int filter_X_mask(Bimage *p, Bimage *pmask, int n_cicli, float r_maxx, int passata);
int support_f(Bimage *p, float rmax);
