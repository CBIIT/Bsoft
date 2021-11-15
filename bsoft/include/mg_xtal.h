/**
@file	mg_xtal.h
@brief	Header file for functions to process crystallographic data.
@author Bernard Heymann and Samuel Payne
@date	Created: 20061110
@date	Modified: 20180410
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
int			mg_unitcell_vectors(Bmicrograph* mg);
long		mg_generate_reflections(Bmicrograph* mg, Vector3<double> real_size, double resolution);
int			img_mask_reflections(Bimage* p, Bstrucfac* sflist, double radius);

