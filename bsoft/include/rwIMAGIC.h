/**
@file	rwIMAGIC.h
@brief	Header file for reading and writing Image Science's Imagic files
@author Bernard Heymann
@date	Created: 19990424
@date	Modified: 20111217

	Format: 2D image file format for the program Imagic (Image Science)
**/

#include "rwimg.h"

#define IMAGICSIZE 1024	// Size of the IMAGIC header for each image

struct IMAGIChead {             // file header for IMAGIC data
    	int imn;    	    	//  0	    	image location number (1,2,...)
    	int ifn;    	    	//  1	    	# images following
    	int ierror; 	    	//  2	    	error code: error if >0
    	int nhfr;   	    	//  3	    	# header records per image
    	int ndate;  	    	//  4	    	creation day
    	int nmonth; 	    	//  5	    	creation month
    	int nyear;  	    	//  6	    	creation year
    	int nhour;  	    	//  7	    	creation hour
    	int nminut; 	    	//  8	    	creation minute
    	int nsec;   	    	//  9	    	creation second
    	int npix2;  	    	// 10	    	# 4-byte reals in image
    	int npixel; 	    	// 11	    	# image elements
    	int ixlp;	    		// 12	    	lines per image (Y)
    	int iylp; 	    		// 13	    	pixels per line (X)
    	char type[4];	    	// 14	    	image type
    	int ixold;	    		// 15	    	top-left X coordinate
    	int iyold;	    		// 16	    	top-left Y coordinate
    	float avdens; 	    	// 17	    	average
    	float sigma; 	    	// 18	    	standard deviation
    	float varian; 	    	// 19	    	variance
    	float oldavd;	    	// 20	    	old average
    	float densmax; 	    	// 21	    	maximum
    	float densmin; 	    	// 22	    	minimum
		int complex;			// 23			label indicating that data is complex
		float cxlength;			// 24			cell dimension in Angstr. (x-direction)
		float cylength;			// 25			cell dimension in Angstr. (y-direction)
		float czlength;			// 26			cell dimension in Angstr. (z-direction)
		float calpha;			// 27			cell angle alpha
		float cbeta;			// 28			cell angle beta
    	char name[80]; 	    	// 29-48    	image name
		float cgamma;			// 49			cell angle gamma
		int mapc;				// 50			axis corresponding to columns
		int mapr;				// 51			axis corresponding to rows
		int maps;				// 52			axis corresponding to sections
		int ispg;				// 53			space group
		int nxstart;			// 54			number of 1st column in map
		int nystart;			// 55			number of 1st row in map
		int nzstart;			// 56			number of 1st section in map
		int nxintv;				// 57			number of intervals along X
		int nyintv;				// 58			number of intervals along Y
		int nzintv;				// 59			number of intervals along Z
		int izlp;				// 60			number of 2D planes in 3D data
		int i4lp;				// 61			number of images
		int i5lp;				// 62
		int	i6lp;				// 63
		float alpha;			// 64			Euler angle alpha
		float beta;				// 65			Euler angle beta
		float gamma;			// 66			Euler angle gamma
		int imavers;			// 67			IMAGIC-5 version (yyyymmdd)
		int realtype;			// 68			floating point type, machine stamp
    	int buffer[29]; 		// 69-97    	Variables that control the buffering
		int ronly;				// 98			flag in calling program to open file readonly
		float angle;			// 99			last rotation angle
		float rcp;				// 100			rotational correlation peak
		int ixpeak;				// 101			shift correlation peak in X direction
		int iypeak;				// 102			shift correlation peak in Y direction
		float ccc;				// 103			cross correlation peak hight
		float errar;			// 104			error in angular reconstitution
		float err3d;			// 105			error in 3D reconstruction
		int ref;				// 106			(multi-) reference number
		float classno;			// 107			class number in classification
		float locold;			// 108			location number before CUT-IMAGE (boxing)
		float oldavd2;			// 109			old average density
		float oldsigma;			// 110			old sigma
		float xshift;			// 111			last shift in X direction
		float yshift;			// 112			last shift in Y direction
		float numcls;			// 113			number of class members 
		float ovqual;			// 114			overall class quality 
		float eangle;			// 115			equivalent angle
		float exshift;			// 116			equivalent shift in X direction
		float eyshift;			// 117			equivalent shift in Y direction
		float cmtotvar;			// 118			used in MSA/IMAGECOOS
		float informat;			// 119			Gauss norm / real*FT Space information of the data set
		int numeigen;			// 120			number of eigen values in MSA
		int niactive;			// 121			number of active images in MSA
		float resolx;			// 122			Angstrom per pixel/voxel in X direction
		float resoly;			// 123			Angstrom per pixel/voxel in Y direction
		float resolz;			// 124			Angstrom per pixel/voxel in Z direction
		float alpha2;			// 125			Euler angle alpha (from projection matching)
		float beta2;			// 126			Euler angle beta (from projection matching)
		float gamma2;			// 127			Euler angle gamma (from projection matching)
		float nmetric;			// 128			Metric used in MSA calculations
		float actmsa;			// 129			a flag indicating whether the "image" is
		float coosmsa[69];		// 130-198		co-ordinates of "image" along factorial axis 
    	char history[228];     	// 199-255  	history
} ;

/*
Machine stamp:
16777216 for VAX/VMS 
33686018 for DEC/OSF,ULTRIX, LINUX, MS Windows 
67372036 for SiliconGraphics, SUN, HP, IBM
*/

struct IMAGIChead_old {             // file header for IMAGIC data
    	int imn;    	    	//  0	    	image location number (1,2,...)
    	int ifn;    	    	//  1	    	# images following
    	int ierror; 	    	//  2	    	error code: error if >0
    	int nhfr;   	    	//  3	    	# header records per image
    	int ndate;  	    	//  4	    	creation day
    	int nmonth; 	    	//  5	    	creation month
    	int nyear;  	    	//  6	    	creation year
    	int nhour;  	    	//  7	    	creation hour
    	int nminut; 	    	//  8	    	creation minute
    	int nsec;   	    	//  9	    	creation second
    	int npix2;  	    	// 10	    	# 4-byte reals in image
    	int npixel; 	    	// 11	    	# image elements
    	int ixlp;	    		// 12	    	lines per image (Y)
    	int iylp; 	    		// 13	    	pixels per line (X)
    	char type[4];	    	// 14	    	image type
    	int ixold;	    		// 15	    	top-left X coordinate
    	int iyold;	    		// 16	    	top-left Y coordinate
    	float avdens; 	    	// 17	    	average
    	float sigma; 	    	// 18	    	standard deviation
    	float varian; 	    	// 19	    	variance
    	float oldavd;	    	// 20	    	old average
    	float densmax; 	    	// 21	    	maximum
    	float densmin; 	    	// 22	    	minimum
//    	double sum; 	    	// 23+24		sum of densities
//    	double squares;  		// 25+26		sum of squares
    	float dummy[4]; 		// 23-26		dummy place holder
    	char lastpr[8];	    	// 27+28    	last program writing file
    	char name[80]; 	    	// 29-48    	image name
    	float extra_1[8]; 		// 49-56    	additional parameters
    	float eman_alt; 		// 57    		EMAN: equiv to psi & PFT omega
    	float eman_az; 			// 58    		EMAN: equiv to theta
    	float eman_phi; 		// 59    		EMAN: equiv to phi
    	float extra_2[69]; 		// 60-128    	additional parameters
		float euler_alpha;		// 129			Euler angles:	psi
		float euler_beta;		// 130							theta
		float euler_gamma;		// 131							phi
		float proj_weight;		// 132			weight of each projection
    	float extra_3[66];  	// 133-198    	additional parameters
    	char history[228];     	// 199-255  	history
} ;


// I/O prototypes
int 		readIMAGIC(Bimage* p, int readdata, int img_select);
int 		writeIMAGIC(Bimage* p);
