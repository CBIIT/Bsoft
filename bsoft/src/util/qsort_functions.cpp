/**
@file	qsort_functions.cpp
@brief	Library utility functions for use in qsort calls.
@author Bernard Heymann
@date	Created: 20010516
@date	Modified: 20100724
**/

#include "qsort_functions.h"

/**
@brief 	Utility function for sorting floating point values in qsort.
@param 	*x		first floating point value.
@param 	*y		second floating point value.
@return 			-1 if x < y and 1 otherwise.
**/
int			QsortSmallToLargeFloat(const void *x, const void *y)
{
	float* 		r1 = (float*) x;
	float* 		r2 = (float*) y;
	if( *r1 < *r2 ) return -1;
	else return 1;
}

/**
@brief 	Utility function for sorting floating point values in qsort.
@param 	*x		first floating point value.
@param 	*y		second floating point value.
@return 			-1 if x > y and 1 otherwise.
**/
int			QsortLargeToSmallFloat(const void *x, const void *y)
{
	float* 		r1 = (float*) x;
	float* 		r2 = (float*) y;
	if( *r1 > *r2 ) return -1;
	else return 1;
}

/**
@brief 	Utility function for sorting double precision floating point values in qsort.
@param 	*x		first floating point value.
@param 	*y		second floating point value.
@return 			-1 if x < y and 1 otherwise.
**/
int			QsortSmallToLargeDouble(const void *x, const void *y)
{
	double* 		r1 = (double*) x;
	double* 		r2 = (double*) y;
	if( *r1 < *r2 ) return -1;
	else return 1;
}

/**
@brief 	Utility function for sorting double precision floating point values in qsort.
@param 	*x		first floating point value.
@param 	*y		second floating point value.
@return 			-1 if x > y and 1 otherwise.
**/
int			QsortLargeToSmallDouble(const void *x, const void *y)
{
	double* 		r1 = (double*) x;
	double* 		r2 = (double*) y;
	if( *r1 > *r2 ) return -1;
	else return 1;
}

/**
@brief 	Utility function for sorting int float pairs in qsort.
@param 	*x		first pair.
@param 	*y		second pair.
@return 			-1 if x < y, 1 if x > y, and 0 otherwise.
**/
int			QsortSmallToLargeIntFloat(const void *x,const void *y)
{
	int_float*		r1 = (int_float *) x;
	int_float*		r2 = (int_float *) y;

	if ( r1->f < r2->f ) return -1;
	else if ( r1->f > r2->f ) return 1;
	else return 0;
}

/**
@brief 	Utility function for sorting int float pairs in qsort.
@param 	*x		first pair.
@param 	*y		second pair.
@return 			-1 if x > y, 1 if x < y, and 0 otherwise.
**/
int			QsortLargeToSmallIntFloat(const void *x,const void *y)
{
	int_float*		r1 = (int_float *) x;
	int_float*		r2 = (int_float *) y;

	if ( r1->f > r2->f ) return -1;
	else if ( r1->f < r2->f ) return 1;
	else return 0;
}

