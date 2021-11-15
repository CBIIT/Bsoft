/**
@file	rwPIF.h
@brief	Header file for reading and writing PIF files
@author Bernard Heymann
@date	Created: 19991112
@date	Modified: 20111218

	Format: 3D crystallographic image file format for the Purdue package
**/

#include "rwimg.h"

#define PIFSIZE   1024	// Minimum size of the PIF header

struct PIFhead {               			// file header for PIF data
        char file_id[8];                //  0   0	magic bytes
		char realscalefactor[16];		//	1	8	floatint -> real
		int numimages;					//	2  24	number of images
		int endianness; 				//	3  28	0=little, 1=big
		char genprogram[32];			//	4  32	last program
		int htype;						//	5  64	1=all images same dimensions
		int nx; 						//	6  68	image size
        int ny;                 		//  7  72
        int nz;                 		//  8  76
        int mode;               		//  9  80	data type (large number)
//		char label[428];				// 10  84	extra to 512
		int even;						// 11  84
		int mrcX;						// 12  88
		int mrcY;						// 13  92
		int mrcZ;						// 14  96
		char floatScaleStr[16];			// 15 100
		char a0[13];					// 16 116
		char a1[13];					// 17 129
		char a2[13];					// 18 142
		char a3[13];					// 19 155
		char a4[13];					// 20 168
		char a5[13];					// 21 181
		char dummy[62];					// 22 194	extra to 256
		char label[256];				// 23 256	extra to 512
};

struct PIFcolor {						//	colormap
	short r[256];
	short g[256];
	short b[256];
} ;

struct PIFimagehead {					// header for each image
		int nx; 						//	0	0	size
		int ny; 						//	1	4
		int nz; 						//	2	8
		int mode;						//	3  12	data type
		int bkgnd;						//	4  16	background value
		int packradius; 				//	5  20	radius of boxed image
        int nxstart;            		//  6		unit cell offset
        int nystart;            		//  7
        int nzstart;            		//  8
        int mx;                 		//  9		unit cell size in voxels
        int my;                 		// 10  40
        int mz;                 		// 11
        int xlength;                	// 12		cell dimensions in A
        int ylength;                	// 13		(converted to floating point)
        int zlength;                	// 14
        int alpha;            			// 15		cell angles in degrees
        int beta;             			// 16		(converted to floating point)
        int gamma;            			// 17
        int mapc;               		// 18		column axis
        int mapr;               		// 19		row axis
        int maps;               		// 20  80	section axis
        int min;	             		// 21		minimum density value (->fp)
        int max;	             		// 22   	maximum density value (->fp)
        int mean;           	  		// 23		average density value (->fp)
		int stddev;						// 24		standard deviation (->fp)
        int ispg;               		// 25 100	space group number
        int nsymbt;             		// 26		bytes used for sym. ops. table
		int xorigin;					// 27		x origin
		int yorigin;					// 28		y origin
		char title[80]; 				// 29 116	title/description (628)
		char timestamp[32]; 			// 30 196	date & time last modified (708)
		char micrograph[16];			// 31 228	unique micrograph number (740)
		char scannumber[8]; 			// 32 244	scan number of micrograph (756)
		int aoverb;						// 33 252
		int map_abang;					// 34 256
		int dela; 						// 35 260
		int delb;						// 36
		int delc; 						// 37
		int t_matrix[6];				// 38
		int dthe; 						// 39
		int dphi_90;					// 40 300
		int symmetry; 					// 41
		int binfactor;					// 42		image compression factor
		int a_star;						// 43		reciprocal space unit cell
		int b_star;						// 44
		int c_star;						// 45
		int alp_star; 					// 46
		int bet_star; 					// 47
		int gam_star; 					// 48
		int zorigin;					// 49		z origin
		int pixelsize_x;				// 50 340   (converted to floating point)
		int pixelsize_y;				// 51		
		int pixelsize_z;				// 52		
		int view_x;						// 53		View vector
		int view_y; 					// 54		(converted to floating point)
		int view_z; 					// 55
		int view_angle; 				// 56		View rotation angle (->fp)
		char extra[144];				// 57 368
} ;

struct PIFimagehead_new {					// header for each image - new PIF specification
		int nx; 						//	0	0	size
		int ny; 						//	1	4
		int nz; 						//	2	8
		int mode;						//	3  12	data type
		int bkgnd;						//	4  16	background value
		int packradius; 				//	5  20	radius of boxed image
        int nxstart;            		//  6		unit cell offset
        int nystart;            		//  7
        int nzstart;            		//  8
        int mx;                 		//  9		unit cell size in voxels
        int my;                 		// 10  40
        int mz;                 		// 11
        int xlength;                	// 12		cell dimensions in A
        int ylength;                	// 13		(converted to floating point)
        int zlength;                	// 14
        int alpha;            			// 15		cell angles in degrees
        int beta;             			// 16		(converted to floating point)
        int gamma;            			// 17
        int mapc;               		// 18		column axis
        int mapr;               		// 19		row axis
        int maps;               		// 20  80	section axis
        int min;	             		// 21		minimum density value (->fp)
        int max;	             		// 22   	maximum density value (->fp)
        int mean;           	  		// 23		average density value (->fp)
		int stddev;						// 24		standard deviation (->fp)
        int ispg;               		// 25 100	space group number
        int nsymbt;             		// 26		bytes used for sym. ops. table
		int xorigin;					// 27		x origin
		int yorigin;					// 28		y origin
		char title[80]; 				// 29 116	title/description (628)
		char timestamp[32]; 			// 30 196	date & time last modified (708)
		char micrograph[16];			// 31 228	unique micrograph number (740)
		char scannumber[8]; 			// 32 244	scan number of micrograph (756)
		int aoverb;						// 33 252
		int map_abang;					// 34 256
		int dela; 						// 35 260
		int delb;						// 36
		int delc; 						// 37
		int t_matrix[6];				// 38
		int dthe; 						// 39
		int dphi_90;					// 40 300
		int symmetry; 					// 41
		int binfactor;					// 42		image compression factor
		int a_star;						// 43		reciprocal space unit cell
		int b_star;						// 44
		int c_star;						// 45
		int alp_star; 					// 46
		int bet_star; 					// 47
		int gam_star; 					// 48
		int pixelSize;					// 49
		int pixSizUnits;				// 50
		int res_out;					// 51		
		int ctfMode;					// 52		
		int tempFactor;					// 53
		int wiener; 					// 54
		int borderType; 				// 55
		int borderWidth; 				// 56
		int datar4_min; 				// 57
		int datar4_max; 				// 58
		int rad_lo;						// 59
		int rad_hi;						// 60
		int transX;						// 61
		int transY;						// 62
		int transZ;						// 63
		int isTrunc;					// 64
		int origX;						// 65 400
		int origY;						// 66
		int origZ;						// 67
		char extra[72];					// 68 412
		int pixelsize_x;				// 69 484   (converted to floating point)
		int pixelsize_y;				// 70		
		int pixelsize_z;				// 71		
		int view_x;						// 72		View vector
		int view_y; 					// 73		(converted to floating point)
		int view_z; 					// 74
		int view_angle; 				// 75 508	View rotation angle (->fp)
} ;

// I/O prototypes
int 		readPIF(Bimage* p, int readdata, int img_select);
int 		writePIF(Bimage* p);

