/**
@file	rwGRD.h
@brief	Header file for reading and writing Basel GRD files
@author Bernard Heymann
@date	Created: 19990410
@date	Modified: 20210304

	Format: 3D crystallographic image file format (Basel)
**/

#include "rwimg.h"

#define GRDSIZE     512	// Size of the GRD header (constant)

struct GRDhead {                // file header for GRD data
        int magic;  	        //  0	0   Magic number (42 for new spec)
        int version;  	        //  1	4   Old: 1/2; New: 1??
        int mode;   	        //  2	8   Bsoft data and compound types
        int doffset;	        //  3  12   Data offset (0 => 512 bytes)
        int nx;                 //  4  16	Image size
        int ny;                 //  5  20
        int nz;                 //  6  24
        int nn;          	  	//  7  28	Number of images
        int channels;           //  8  32	Number of channels
        int slowest;            //  9	Dimension order (deprecated)
        float a;            	// 10  40	Cell dimensions in A
        float b;            	// 11
        float c;            	// 12
        float alpha;            // 13		Cell angles
        float beta;             // 14   
        float gamma;            // 15
        float amin;             // 16		Minimum density value
        float amax;             // 17  	    Maximum density value
        float amean;            // 18		Average density value
        float arms;          	// 19	    RMS deviation
		float ux;				// 20  80	Sampling/voxel size
		float uy;				// 21
		float uz;				// 22
		float bkg;				// 23		Background
		float ox;				// 24		Origin
		float oy;				// 25
		float oz;				// 26
		float vx;				// 27		View
		float vy;				// 28
		float vz;				// 29
		float va;				// 30
		long rle;				// 31		Run-length encoded size (0=no compression)
        float extra[31];		// 33-63	user-defined info
		char datetime[16];		// 64 256	Date and time as YYYYMMDDhhmmss
		char symmetry[32];		// 68 272	Symmetry string
		char label[208];		// 76 304	Text label
} ;


// I/O prototypes
int 		readGRD(Bimage* p, int readdata, int img_select);
int 		writeGRD(Bimage* p, int flags);
