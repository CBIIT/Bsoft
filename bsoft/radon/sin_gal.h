/**
@file	sin_gal.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 23 01 2006
**/

#include "kernlut.h"
#include "Bimage.h"
#include "Vector3.h"

Bimage *sin_gal(Bimage *p, int tipo, Vector3<float> shift, kernel *ker, int padd, int n_theta);

Bimage*	copy_I_part_rad(Bimage * p);
int copy_II_part_rad(Bimage *p);
