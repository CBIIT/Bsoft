/**
@file	timer.cpp
@brief	Utilities for timing functions
@author Bernard Heymann 
@date	Created: 20010316
@date	Modified: 20221103
**/

#include <ctime>
#include <sys/time.h>
//#include <chrono>

#include "timer.h" 
#include "utilities.h" 
 
// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Returns the current local time.
@return tm* 			the current local time struct.

 The time format is in the tm struct.

*/
tm*			get_localtime()
{
	time_t 			time_sec = time(NULL);
	return localtime(&time_sec);
}

/**
@brief 	Returns the current time.
@return double 			the current time in seconds.

 The time format is in seconds since Jan. 1, 1970, including microseconds.

*/
double		getwalltime()
{
	timeval 		thetime;
	
	gettimeofday(&thetime, NULL);
	
	return (double) (thetime.tv_sec + 1e-6L*thetime.tv_usec);
}
/*
double		getwalltime()
{
	return chrono::system_clock::now();
}
*/

/**
@brief 	Returns the clock time.
@return double 			the clock time in seconds.


	The time format is in clock seconds since the start of the process.

**/
double		getcputime()
{
	clock_t			thetime = clock();
	
	return (double) (thetime/CLOCKS_PER_SEC);
}

/**
@brief 	Starts timer and prints the time.
@return double 			the current time.
**/
double		timer_start()
{
	double 		thetime = getwalltime();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG timer_start: time = " << thetime << endl;
	
	return thetime;
}

/**
@brief 	Reports the time elapsed since a given time.
@param	lasttime		the given initial time.
@return double 			the current time.
**/
double		timer_report(double lasttime)
{
	double 		difftime = getwalltime() - lasttime;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG timer_report: difftime = " << difftime << endl;
	
	if ( verbose & VERB_TIME )
		cout << "Time elapsed: " << difftime << " s" << endl;
	
	return difftime;
}

