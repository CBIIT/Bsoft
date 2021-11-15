/**
@file	versions.h
@brief	Library functions to encode versions of all included packages
@author Bernard Heymann
@date	Created: 20200610
@date	Modified: 20200610
**/

#ifndef _VERSIONS_

#include "fft.h"
#include "tiffvers.h"
#include "jversion.h"
#include "png.h"
#include <tk.h>
#include "utilities.h"

/**
@brief 	Displaying library versions.
@return int				0
**/
int			show_library_versions()
{
	cout << "Library versions:" << endl;
	cout << "----------------" << endl << endl;
	cout << "Bsoft library: " << BVERSION << endl << endl;
	cout << "FFTW3: " << FFTW3_VERSION << endl << endl;
	cout << TIFFLIB_VERSION_STR << endl << endl;
	cout << "JPEG: " << JVERSION << endl << endl;
	cout << "PNG: " << PNG_HEADER_VERSION_STRING << endl << endl;
	cout << "Tcl/Tk: " << TK_PATCH_LEVEL << endl << endl;
	return 0;
}

#define _VERSIONS_
#endif

