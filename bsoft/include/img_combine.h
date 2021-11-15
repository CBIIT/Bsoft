/**
@file	img_combine.h
@brief	Header file for combining two images in various ways
@author Bernard Heymann
@date	Created: 20000430
@date	Modified: 20190208
**/

#include "rwimg.h"

// Function prototypes
Bimage* 	img_add(Bstring* file_list, int flags);
Bimage*	 	img_setup_combined(Bstring* file_list, long& nimg, int cat=0);
Bimage*		img_catenate(Bstring* file_list, Bstring& rawstring, DataType newdatatype, 
				Vector3<long> nusize, int setZslices=0, int fill_type=0, double fill=0, 
				double newavg=0, double newstd=0);
Bimage* 	img_add_weighed(Bstring* file_list, vector<double> weight,
					double newavg=0, double newstd=0, int flags=0);
