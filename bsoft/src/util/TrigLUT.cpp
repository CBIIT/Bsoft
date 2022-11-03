/**
@file	TrigLUT.cpp
@brief	Trigonometric look-up tables
@author Bernard Heymann
@date	Created: 20220713
@date	Modified: 20220715
**/

#include "TrigLUT.h"

#define TRIGACC	100000

// Definition of global variables 
extern int 	verbose;

TrigLUT		cos_lut(TRIGACC, cos);
TrigLUT		sin_lut(TRIGACC, sin);
TrigLUT		atan2_lut(TRIGACC, atan);


