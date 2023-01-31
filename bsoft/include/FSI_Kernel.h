/**
@file	FSI_Kernel.h
@brief	Header file for frequency space interpolation kernel.
@author	Bernard Heymann
@date	Created: 20051102
@date	Modified: 20190212
**/

#include "Vector3.h"
#include "Complex.h"
#include "utilities.h"

#ifndef _FSI_Kernel_
class FSI_Kernel {
private:
	long 			wd;			// Width of kernel
	long 			power;		// Interpolation power
	long 			ndiv;		// Number of divisions in sampling
	long 			size;		// Total number of elements in kernel
	vector<double>	data;		// Kernel values
public:
	FSI_Kernel() { wd = power = ndiv = size = 0; }
/**
@brief 	Calculates a reciprocal space interpolation lookup table.
@param 	kernel_width			table width.
@param 	kernel_power			interpolation power.
**/
	FSI_Kernel(long kernel_width, long kernel_power) {
		long			i, j, ik, wd_2;
		double			s1, s2, c, arg, arg1, arg2, iwd, piwd, incr;
	
		power = kernel_power;

		if ( power < 0 ) wd = 1;
    	else wd = kernel_width;
		
		wd_2 = (wd - 1) / 2;
		ndiv = 512;
		size = (ndiv + 1) * wd;
		data.resize(size);

		for ( i = 0; i < size; i++ ) data[i] = 0;
	
		if ( power < 0 ) {
			for ( i = 0; i < size; i++ ) data[i] = 1;
		} else {
			iwd = 1.0L / wd;
			piwd = M_PI * iwd;
			incr = 1.0L / ndiv;
			arg = incr;
			data[wd_2] = 1;
			for ( i = 1, ik = wd; i < ndiv; i++ ) {
				arg1 = M_PI * (arg + wd_2);
				s1 = iwd * sin (arg1);
				arg2 = arg1 * iwd;
				for ( j = 0; j < wd; j++, ik++ ) {
					s2 = sin (arg2);
					c = cos (arg2);
					data[ik] = (s1 / s2) * pow (c, power);
					arg2 -= piwd;
					s1 = -s1;
				}
				arg += incr;
			}
			ik += wd_2 + 1;
			data[ik] = 1;
		}
	}
	long		width() { return wd; }
	long		half_width() { return (wd - 1)/2; }
	long		divisions() { return ndiv; }
	double		operator[](long i) const { return data[i]; }
	void		show() {
		long				i, j, ik;
		
		cout << "Kernel parameters:" << endl;
		cout << "Width:                          " << wd << endl;
		cout << "Half width:                     " << half_width() << endl;
		cout << "Power:                          " << power << endl;
		cout << "Divisions:                      " << ndiv << endl;
		cout << "Data size:                      " << size << endl << endl;

		cout << "Kernel content:" << endl;
		for ( i=ik=0; i<=ndiv; i++ ) {
			cout << i;
			for ( j=0; j<wd; j++, ik++ ) {
				cout << tab << data[ik];
			}
			cout << endl;
		}
		cout << endl;
	}
} ;

#define _FSI_Kernel_
#endif


