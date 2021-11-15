/**
@file	rwSPE.h
@brief	Header file for reading and writing SPE files
@author Bernard Heymann
@date	Created: 20081105
@date	Modified: 20111217

	Format: SPE Princeton Instruments CCD image file format
**/

#include "rwimg.h"

#define SPESIZE   4100	// Size of the SPE header (constant)


struct SPEhead {             // file header for SPE data
		short pad1[10];			//	0	0
        char date[10];			//  10	20	Date
		short pad2[6];			//	15	30
        short xdim;				//  21	42	X size
		short pad3[32];			//	22	44
		short datatype;			//	54	108	0=float, 1=int, 2=short, 3=unsigned short
		short pad4[45];			//	55	110
		char Comments[400];		//	100	200
		short pad5[28];			//	300	600
        short ydim;				//  328	656	Y size
		char pad6[788];			//	329	658
        short NumFrames;		//  723	1446	Number of frames
		float MaxIntensity;		//	725 1450	max intensity of data
		float MinIntensity;		//	727	1454	min intensity of data 
        short pad7[1321];		//  729	1458	Padding
} ;

// I/O prototypes
int 		readSPE(Bimage* p, int readdata, int img_select);
int 		writeSPE(Bimage* p);
