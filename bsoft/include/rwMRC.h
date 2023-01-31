/**
@file	rwMRC.h
@brief	Header file for reading and writing MRC files
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20111217

	Format: 3D crystallographic image file format for the MRC package
**/

#include "rwimg.h"

#define MRCSIZE    1024	// Minimum size of the MRC header (when nsymbt = 0)

struct MRCheadold {                		// file header for MRC data
        int nx;                 		//  0   0	image size
        int ny;                 		//  1   4
        int nz;                 		//  2   8
        int mode;               		//  3		0=uchar,1=short,2=float
        int nxStart;            		//  4		unit cell offset
        int nyStart;            		//  5
        int nzStart;            		//  6
        int mx;                 		//  7		unit cell size in voxels
        int my;                 		//  8
        int mz;                 		//  9
        float a;                		// 10   40	cell dimensions in A
        float b;                		// 11
        float c;                		// 12
        float alpha;            		// 13		cell angles in degrees
        float beta;             		// 14   
        float gamma;            		// 15
        int mapc;               		// 16		column axis
        int mapr;               		// 17		row axis
        int maps;               		// 18		section axis
        float amin;             		// 19		minimum density value
        float amax;             		// 20   80	maximum density value
        float amean;            		// 21		average density value
        int ispg;               		// 22		space group number
        int nsymbt;             		// 23		bytes used for sym. ops. table
        float extra[29];        		// 24		user-defined info
        float xOrigin;          		// 53		phase origin in pixels
        float yOrigin;          		// 54
        int nlabl;              		// 55		number of labels used
        char labels[10][80];    		// 56-255	10 80-character labels
} ;

struct MRChead2000 {                		// file header for MRC data
        int nx;                 		//  0   0	image size
        int ny;                 		//  1   4
        int nz;                 		//  2   8
        int mode;               		//  3		0=char,1=short,2=float
        int nxStart;            		//  4		unit cell offset
        int nyStart;            		//  5
        int nzStart;            		//  6
        int mx;                 		//  7		unit cell size in voxels
        int my;                 		//  8
        int mz;                 		//  9
        float a;                		// 10   40	cell dimensions in A
        float b;                		// 11
        float c;                		// 12
        float alpha;            		// 13		cell angles in degrees
        float beta;             		// 14   
        float gamma;            		// 15
        int mapc;               		// 16		column axis
        int mapr;               		// 17		row axis
        int maps;               		// 18		section axis
        float amin;             		// 19		minimum density value
        float amax;             		// 20   80	maximum density value
        float amean;            		// 21		average density value
        int ispg;               		// 22		space group number
        int nsymbt;             		// 23		bytes used for sym. ops. table
        float extra[25];        		// 24		user-defined info
        float xOrigin;          		// 49		phase origin in pixels
        float yOrigin;          		// 50	200
        float zOrigin;          		// 51
        char map[4];  	        		// 52	    identifier for map file ("MAP ")
        char machst[4];         		// 53		machine stamp
        float arms;             		// 54	    RMS deviation
        int nlabl;              		// 55		number of labels used
        char labels[800];				// 56-255	10 80-character labels
} ;

struct MRChead {                		// file header for MRC data
        int nx;                 		//  0   0	image size
        int ny;                 		//  1   4
        int nz;                 		//  2   8	mz*#vol
        int mode;               		//  3		0=schar,1=short,2=float,3=complex_ushort,4=complex_float,6=ushort
        int nxStart;            		//  4		unit cell offset
        int nyStart;            		//  5
        int nzStart;            		//  6
        int mx;                 		//  7		unit cell size in voxels
        int my;                 		//  8
        int mz;                 		//  9		#vol = nz/mz
        float a;                		// 10   40	cell dimensions in A
        float b;                		// 11
        float c;                		// 12
        float alpha;            		// 13		cell angles in degrees
        float beta;             		// 14   
        float gamma;            		// 15
        int mapc;               		// 16		column axis
        int mapr;               		// 17		row axis
        int maps;               		// 18		section axis
        float amin;             		// 19		minimum density value (amin > amax => not detn)
        float amax;             		// 20   80	maximum density value
        float amean;            		// 21		average density value (< min(amin,amax) => not detn)
        int ispg;               		// 22		space group number/ 2/3/4D config: 0=2D,sg=3D,sg+400=4D
        int nsymbt;             		// 23	92	bytes used for extended header
        float extra1[2];        		// 24		user-defined info
		char exttype[4];				// 26		extended type (CCP4,MRCO,SERI,AGAR,FEI1,JSON)
		int nversion;        			// 27	108	version: year*10 + sub-version [0,9]
		float extra2[10];        		// 28		user-defined info
		int imodStamp;					// 38	152	IMOD stamp: 1146047817
		int imodFlags;					// 39		1=signed bytes, 16=4-bit values
		float extra3[9];        		// 40		user-defined info
		float xOrigin;          		// 49		phase origin in angstroms
        float yOrigin;          		// 50	200
        float zOrigin;          		// 51
        char map[4];  	        		// 52	    identifier for map file ("MAP ")
        char machst[4];         		// 53		machine stamp
        float arms;             		// 54	    RMS deviation (<0 => not detn)
        int nlabl;              		// 55		number of labels used
        char labels[800];				// 56-255	10 80-character labels
} ;

// I/O prototypes
int 		readMRC(Bimage* p, int readdata, int img_select);
int 		writeMRC(Bimage* p);
