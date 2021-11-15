/**
@file	tcltk_bimage.h
@brief	A shared object to load Bsoft image files in TCL/Tk
@author	Bernard Heymann
@date	Created: 20010210
@date	Modified: 20130924
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "rwimg.h"

// Function prototypes
int			image_processing(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[]);

