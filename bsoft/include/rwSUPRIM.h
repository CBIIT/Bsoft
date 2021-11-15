/**
@file	rwSUPRIM.h
@brief	Header file for reading and writing SUPRIM files
@author Bernard Heymann
@date	Created: 19990930
@date	Modified: 20050705

	Format: 3D image file format for the SUPRIM package
**/

#include "rwimg.h"

/* REGISTER ASSIGNMENT */

#define F_AVG        0 
#define BACKGROUND   0
#define ANG_SAMP     1
#define RAD_SAMP     2
#define O_ROW_SZ     3
#define O_COL_SZ     4
#define ROW_ORG      5			// x origin
#define COL_ORG      6			// y origin
#define MASK_RAD     7
#define F_NOTE_R     8
#define F_NOTE_I     9
#define SAMP_DIST   10			// Sampling
#define NSLICES     11			// z dimension
#define ENL_BP      12
#define X_FMAX      12
#define THE_BP      13			// Theta Euler angle
#define PHI_BP      14			// Phi Euler angle
#define FOUR_TR_CEN 15
#define GSBAR       15 
#define X_NIMA      16
#define SYMP2       17
#define X_NPIX      17
#define X_NCLU      18
#define X_NCTR      19
#define X_NAME      20
#define X_LNAM      21
#define X_NTRM      21
#define X_NFAC      22
#define X_NTYP      23
#define FRTHDIM     24
#define X_NIMI      26
#define X_UFAC      27
#define X_WALL      28
#define SLICE_ORG   29			// z origin
#define NUM_POLYS   31
#define X_FBEG      32
#define PSI_BP      40			// Psi Euler angle (added JBH 20011102)
#define X_TRP       70
#define R_ANGX      81
#define R_ANGY      82
#define R_ANGZ      83
#define R_NPROJ     84
#define XGLOB	    90
#define YGLOB	    91
#define ROTGLOB     92
#define XX_CELL     93
#define XY_CELL     94
#define YX_CELL     95
#define YY_CELL     96
#define IR_WINID    97
#define O_SLICE_SZ  98
#define CALIB       99
#define NAV        100

#define SUPRIMSIZE 548			// Minimum size of the SUPRIM header (excluding trace)

#define MAXREG     128

typedef union {
    int l;
    float f;
} REG;

struct SUPRIMhead {             // file header for SUPRIM data
    	int nrow;				//  0	0	y dimension
    	int ncol;				//	1	4	x dimension
    	int format; 			//	2	8	word size (1/2/4/8)
    	int intern; 			//	3  12	data type
    	int type;				//	4  16	file type
    	float min;				//	5  20	minimum
    	float max;				//	6  24	maximum
    	float av;				//	7  28	average
    	float sd;				//	8  32	standard deviation
    	REG   reg[MAXREG];		//  9  36	additional registers (see above for index definitions)
		char trace[1024];		// 10 548	trace (history and date)
} ;


// I/O prototypes
int 		readSUPRIM(Bimage* p, int readdata);
int 		writeSUPRIM(Bimage* p);
