/**
@file	Quaternion.cpp
@brief	Quaternion functions
@author Bernard Heymann 
@date	Created: 20051017
@date	Modified: 20111124
**/

#include "Quaternion.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, Quaternion& q) {
	output.setf(ios::fixed, ios::floatfield);
	output.precision(4);
	output << q.scalar() << tab << q.axis();
	return output;
}



