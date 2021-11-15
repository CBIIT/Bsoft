/**
@file	rwGOODFORD.h
@brief	Header file for reading and writing Peter Goodford's GRID files
@author Bernard Heymann
@date	Created: 20000924
@date	Modified: 20111217

	Format: 3D electrostatic potential format for UHBD
**/

#include "rwimg.h"

#define GOODFORDSIZE     168	// Size of the GOODFORD header (constant)
#define GFSLICESIZE      20 	// Size of the slice header (constant)

struct GOODFORDhead {       // file header for GOODFORD data
	int 	pad1;			//	0		FORTRAN start integer = 160
	char	title[72]; 		//	4		File title
	float	scale;			//	76		multiplier for output
	int 	pad2;			//	80		= 0
	int 	grdflg; 		//	84		= 1, type of grid
	int 	pad3;			//	88		= 0
	int 	pad4;			//	92		= km
	int 	pad5;			//	96		= 1
	int 	pad6;			//	100 	= km
	int 	im; 			//	104 	x-dimension
	int 	jm; 			//	108 	y-dimension
	int 	km; 			//	112 	z-dimension
	float 	h;				//	116 	Sampling/spacing (angstrom/voxel)
	float 	ox;				//	120 	Origin (in angstrom)
	float 	oy;				//	124
	float 	oz;				//	128
	int 	pad7[8];		//	132 	= 0
	int 	pad8;			//	164 	FORTRAN end integer = 160
} ;

struct GFslice_head {
	int 	pad1;			//	FORTRAN start integer = 12
	int 	k;				//	Slice number
	int 	im; 			//	x-dimension
	int 	jm; 			//	y-dimension
	int 	pad2;			//	FORTRAN end integer = 12
} ;

// I/O prototypes
int 		readGOODFORD(Bimage* p, int readdata);
int 		writeGOODFORD(Bimage* p);


