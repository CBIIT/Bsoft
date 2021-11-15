/**
@file	rwBIORAD.h
@brief	Header file for reading and writing BioRad files
@date	Created: 19990427
@date	Modified: 20111217

	Format: BioRad 3D confocal microscopy image file format@author Bernard Heymann
**/

#include "rwimg.h"

#define BIORADSIZE   76	// Size of the BioRad header (constant)


struct BIORADhead {             // file header for BioRad data
        short nx;               //  0	0	X size
        short ny;               //  1	2	Y size
        short nz;               //  2	4	Z size or number of pictures
        short black;            //  3	6	Contrast ramp 0, pixel value black
        short white;            //  4	8	Contrast ramp 255, pixel value white
        short pad1[2];          //  5	10	Padding
        short mode;             //  7	14	Byte = 1, Word = 0
        short pad2[17];         //  8	16	Padding
        short mark;             //  25	50	Tag for merging or marking, usually 0
        short color;            //  26	52	Color selection
        short magic;            //  27	54	File identification = 12345
        short black2;           //  28	56	Contrast ramp 0, 2nd picture merged
        short white2;           //  29	58	Contrast ramp 255, 2nd picture merged
        short color2;           //  30	60	Color selection for 2nd picture merged
        short pad3;             //  31	62	Padding
        short lens;             //  32	64	Lens power (integer) or flag indic. fp
        float mag;              //  33	66	Magnification factor
        short pad4[3];          //  35	70	Padding
} ;

// I/O prototypes
int 		readBIORAD(Bimage* p, int readdata);
int 		writeBIORAD(Bimage* p);
