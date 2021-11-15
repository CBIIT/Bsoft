/**
@file	proj_tool.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "rwimg.h"

struct	LUTable {
	int ncol, sizeT;
	int** tab1;
	float** tab2;
} ;

LUTable *create_table(int ncol, int sizeT);
int kill_table(LUTable *lut);
int write_line(Bimage *p, Bimage *prec, int angolo, float csi, float eta, Bimage *pmask, LUTable *lut, float *priga);
int weigh_radon_transf(Bimage *p, Bimage *pmask);
