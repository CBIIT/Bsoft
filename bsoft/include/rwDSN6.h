/**
@file	rwDSN6.h
@brief	Header file for reading and writing DSN6 files
@author Bernard Heymann
@date	Created: 20011226
@date	Modified: 200111217

	Format: 3D crystallographic image file format for the DSN6 format
**/

#include "rwimg.h"

#define DSN6SIZE   512	// Constant size of the DSN6 header

struct DSN6head {               	// file header for DSN6 data
	short x_start;  				//  0	 0	unit cell offset
	short y_start;  				//  1	 2
	short z_start;  				//  2	 4
	short x_extent; 				//  3    6	image size
	short y_extent; 				//  4    8
	short z_extent; 				//  5   10
	short x_sampling;				//  6	12	unit cell size (voxels)
	short y_sampling;				//  7	14
	short z_sampling;				//  8	16
	short a;                		//	9   18	unit cell dimensions in A
	short b;                		// 10	20	(*cell_scale)
	short c;                		// 11	22
	short alpha;            		// 12	24	unit cell angles in degrees
	short beta;             		// 13   26	(*cell_scale)
	short gamma;            		// 14	28
	short product;            		// 15	30	scale*(253-3)/(max-min)
	short plus;            			// 16	32	(3*max-253*min)/(max-min)
	short cell_scale;            	// 17	34	cell constant scaling factor
	short scale;            		// 18	36	product scaling factor
	short extra[237]; 				// 19	38 - 511
} ;


// I/O prototypes
int 		readDSN6(Bimage* p, int readdata);
int 		writeDSN6(Bimage* p);
