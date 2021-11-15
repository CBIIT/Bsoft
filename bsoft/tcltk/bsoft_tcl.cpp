/**
@file	bsoft_tcl.cpp
@brief	A stub function to interface the Bsoft library within TCL/Tk
@author Bernard Heymann
@date	Created: 20010210
@date	Modified: 20010302
**/

#include <tk.h>


/* this function is dynamically loaded by the "load" command */
extern "C" int Bsoft_tcl_Init(Tcl_Interp *interp)
{
  return TCL_OK;
}

extern "C" int Bsoft_tcl_SafeInit(Tcl_Interp *interp)
{
  return TCL_OK;
}
