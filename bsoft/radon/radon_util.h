/**
@file	radon_util.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "rwimg.h"

int rotZ_cart(double x, double y, double z, double ang, double *xx, double *yy, double *zz);
int rotY_cart(double x, double y, double z, double ang, double *xx, double *yy, double *zz);
int rotX_cart(double x, double y, double z, double ang, double *xx,double *yy,double *zz);
double ang_one_two(double a1, double b1, double a2, double b2);
int sphere(Bimage *p);
int mean_to_0(Bimage *p);
