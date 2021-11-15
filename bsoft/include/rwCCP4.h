/**
@file	rwCCP4.h
@brief	Header file for reading and writing CCP4 files
@author Bernard Heymann
@date	Created: 19990410
@date	Modified: 20010401

	Format: 3D crystallographic image file format for the CCP4 package
**/

#include "rwimg.h"

#define CCP4SIZE   1024	// Minimum size of the CCP4 header (when nsymbt = 0)

struct CCP4head {               		// file header for CCP4 data
        int nx;                 		//  0   0	image size
        int ny;                 		//  1   4
        int nz;                 		//  2   8
        int mode;               		//  3		0=signed char,1=short,2=float
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
        int lskflg; 	        		// 24	    skew transformation flag, none=0, foll=1
        float skwmat[9];        		// 25-33    Skew matrix S (order S11,S12,S13,S21...) if lskflg>0.
        float skwtrn[3];        		// 34-36    Skew translation if lskflg>0.
        float extra[15];        		// 37-51	user-defined info
        char map[4];  	        		// 52	    identifier for map file ("MAP ")
        char machst[4];         		// 53		machine stamp
        float arms;             		// 54	    RMS deviation
        int nlabl;              		// 55		number of labels used
        char labels[800];   			// 56-255	10 80-character labels
} ;


// I/O prototypes
int			set_CCP4_machine_stamp(char* machine_stamp);
int 		readCCP4(Bimage* p, int readdata);
int 		writeCCP4(Bimage* p);
