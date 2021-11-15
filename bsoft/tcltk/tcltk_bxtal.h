/**
@file	tcltk_bmg.h
@brief	A shared object to manage micrograph parameter files in TCL/Tk
@author	Bernard Heymann
@date	Created: 20030813
@date	Modified: 20130726
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "mg_processing.h"

// Function prototypes
int			do_unitcell_vectors(Bproject* project, int objc, Tcl_Obj *CONST objv[]);
Tcl_Obj*	do_spot(Bproject* project, int objc, Tcl_Obj *CONST objv[]);


