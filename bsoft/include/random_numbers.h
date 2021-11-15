/**
@file	random_numbers.h
@brief	Header file for functions for creating random images
@author Bernard Heymann
@date	Created: 19990703
@date	Modified: 20151113
**/

#include <Vector3.h>



// Function prototypes
long 		get_rand_max();
long 		random_seed();
float* 		random_array_uniform(long n, double min, double max);
double 		random_gaussian(double avg, double std);
float* 		random_array_gaussian(long n, double avg, double std);
double		random_poisson(double avg);
float*		random_array_poisson(int n, double avg);
double 		random_logistical(double avg, double std);
float* 		random_array_logistical(long n, double avg, double std);
Vector3<double> 	vector3_random_unit_sphere();
Vector3<double> 	vector3_random(const double min, const double max);
Vector3<double> 	vector3_random(Vector3<double> min, Vector3<double> max);
Vector3<double> 	vector3_random(const double length);
Vector3<double> 	vector3_random_gaussian(double avg, double std);
Vector3<double> 	vector3_xy_random_gaussian(double avg, double std);

