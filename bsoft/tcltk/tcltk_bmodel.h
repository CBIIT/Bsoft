/**
@file	tcltk_bmodel.h
@brief	A shared object to manage model parameter files in TCL/Tk
@author	Bernard Heymann
@date	Created: 20071002
@date	Modified: 20130924
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "rwmodel.h"
#include "rwimg.h"
#include "mg_processing.h"

// Function prototypes
int			model_processing(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[]);


