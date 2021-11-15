/**
@file	qsort_functions.h
@brief	Header file for utility functions for use in qsort calls.
@author Bernard Heymann
@date	Created: 20010516
@date	Modified: 20100724
**/

// Structure for ranking correlation coefficients
#ifndef _int_float_
#define _int_float_
struct int_float { int i; float f; } ;
#endif


// Function prototypes
int			QsortSmallToLargeFloat(const void *x, const void *y);
int			QsortLargeToSmallFloat(const void *x, const void *y);
int			QsortSmallToLargeDouble(const void *x, const void *y);
int			QsortLargeToSmallDouble(const void *x, const void *y);
int			QsortSmallToLargeIntFloat(const void *x,const void *y);
int			QsortLargeToSmallIntFloat(const void *x,const void *y);

