/**
@file	rwimg.h
@brief	Header file for 2D and 3D image I/O
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20180419
**/

#include "Bimage.h"
#include "Vector3.h"
#include "View.h"
#include "UnitCell.h"

#include <fstream>

#define SWAPTRIG	65535	// Threshold file z size above which bytes are swapped

// Function prototypes
Bimage* 	read_img(char* filename, int readdata, int img_select);
Bimage* 	read_img(Bstring filename, int readdata, int img_select);
Bimage* 	read_img(string filename, int readdata, int img_select);
int 		write_img(const char* filename, Bimage* p, int compression);
int 		write_img(Bstring filename, Bimage* p, int compression);
int 		write_img(string filename, Bimage* p, int compression);
int			img_convert_fourier(Bimage* p, FourierType newtransform);

