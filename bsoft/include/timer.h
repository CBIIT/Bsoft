/**
@file	timer.h 
@brief	Utilities for timing functions 
@author Bernard Heymann 
@date	Created: 20010316
@date	Modified: 20100726 
**/
 
//#include <time.h>
//#include <sys/time.h>


// Function prototypes 
double		getwalltime();
double		getcputime();
double		timer_start();
double		timer_report(double lasttime);

