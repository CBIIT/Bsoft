/**
@file	kernlut.h

@author	P.L. Bellon, F. Cantele and S. Lanzavecchia
	Dip. Chimica Strutturale e Stereochimica Inorganica
	Via Venezian 21, 20133 Milano, Italy

@date	Created: 7 04 2003
@date	Modified: 07 07 2005
**/

#ifndef _kernel_
typedef struct
   {
   int nk;			// Width of kernel
   int power;
   int nkern;		// Number of divisions in sampling
   int nker;		// Total number of elements in kernel
   int bordo;		// Half of the kernel width
   float *kern;		// Kernel values
   } kernel;
#define _kernel_
#endif

int kernlut(kernel *ker);
kernel *kernel_prep(int nk, int power);
