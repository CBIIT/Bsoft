/**
@file	rwEM.h
@brief	Header file for reading and writing EM files
@author Bernard Heymann
@date	Created: 19990418
@date	Modified: 20111217

	Format: 3D image file format for the EM package
**/

#include "rwimg.h"

#define EMSIZE      512	// Size of the EM header (constant)

struct EMhead {				// file header for EM data
	char stamp[4];  	    //  0   0	machine & data stamp
	int nx;                 //  1   4	image size
	int ny;                 //  2   8
	int nz;                 //  3  12
	char label[80];	        //  4  16	80-character comment
	int volt;				//	   96	Acceleration voltage (kV)
	int Cs; 				//	  100	Spherical aberration (micron)
	int aperture;			//	  104	Aperture (mrad)
	int mag;				//	  108	Nominal magnification wrt film
	int post_mag;			//	  112	Additional magnification at CCD (x1000)
	int exp_time;			//	  116	Exposure time (ms)
	int nsampling;			//	  120	Number of sampling points
	int pixel_size; 		//	  124	CCD or scanner pixel size (nm)
	int pad1;				//	  128
	int field_length;		//	  132	= sampling points x pixel size
	int defocus;			//	  136	Defocus determined from power spectrum (angstrom)
	int astigmatism;		//	  140	Magnitude of astigmatism (angstrom)
	int ast_angle;			//	  144	Angle of astigmatism (deg x 1000)
	int drift;				//	  148	Drift during exposure time (angstrom)
	int drift_angle;		//	  152	Angle of drift (deg x 1000)
	int defocus_inc;		//	  156	Defocus increment for focal series (angstrom)
	int pad2;				//	  160
	int pad3;				//	  164
	int tilt_angle; 		//	  168	Tilt angle (deg x 1000)
	int tilt_dir;			//	  172	Tilt direction (deg x 1000)
	char extra[336];		//    176	user-defined info
} ;

// I/O prototypes
int 		readEM(Bimage* p, int readdata);
int 		writeEM(Bimage* p);
