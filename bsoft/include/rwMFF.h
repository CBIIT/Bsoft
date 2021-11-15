/**
@file	rwMFF.h
@brief	Header file for reading and writing the What If MFF files
@author Bernard Heymann
@date	Created: 19990424
@date	Modified: 20111217

	Format: 3D crystallographic image file format for the program WHAT IF
**/

#include "rwimg.h"

#define MFFSIZE	    268 		// Size of the MFF header (constant)

struct MFFhead {                // file header for What If MFF data
        int pad1;  	    		//  0	    	padding
        char title[80];     	//  1-20	80-character title
        int pad2[2];  	    	// 21-22	padding
        float a;	    		// 23    	Unit cell dimensions
        float b;	    		// 24
        float c;	    		// 25
        float alpha; 	    	// 26    	Unit cell angles
        float beta; 	    	// 27
        float gamma; 	    	// 28
        int pad3[2];  	    	// 29-30	padding
        float stepu;	    	// 31	    	Step size in direction U
        float stepv;	    	// 32	    	Step size in direction V
        float stepw;	    	// 33	    	Step size in direction W
        int pad4[2];  	    	// 34-35	padding
        int idiru;  	    	// 36	    	Order for direction U
        int idirv;  	    	// 37	    	Order for direction V
        int idirw;  	    	// 38	    	Order for direction W
        int pad5[2];  	    	// 39-40	padding
        int ndiva;  	    	// 41	    	Unit cell size in pixels
        int ndivb;  	    	// 42
        int ndivc;  	    	// 43
        int pad6[2];  	    	// 44-45	padding
        int nx;     	    	// 46	    	Size of data
        int ny;     	    	// 47
        int nz;     	    	// 48
        int pad7[2];  	    	// 49-50	padding
        float uvworg[3];    	// 51-53    	Fractional coordinates of first point
        int pad8[2];  	    	// 54-55	padding
        float uvwmat[3][3]; 	// 56-64    	Fractionalization matrix
        int pad9[2];  	    	// 65-66	padding
} ;


// I/O prototypes
int 		readMFF(Bimage* p, int readdata);
int 		writeMFF(Bimage* p);
