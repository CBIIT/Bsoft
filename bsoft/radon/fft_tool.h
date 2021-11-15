/**
@file	fft_tool.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#include "rwimg.h"

int fft_1D_forward(Bimage *p);
int fft_1D_backward(Bimage *p);
int fft_2D_forward(Bimage *p);
int fft_2D_backward(Bimage *p);

int farfarig(Bimage *p, int ifl);
int mixrig(Bimage *p, int ifl);
int scramrig(Bimage *p);
int abtori(Bimage *p);
int ritoab(Bimage *p);
int zeroes(Bimage *p);
int rephase_orig(Bimage *p);
