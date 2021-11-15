/**
@file	bshow.cpp
@brief	A shared object to load Bsoft image files in TCL/Tk
@author Bernard Heymann
@date	Created: 20010210
@date	Modified: 20200514
**/

// Tk must be included before anything else to remedy symbol conflicts
#include <tk.h>

#include "tcltk_bimage.h"
#include "tcltk_bmg.h"
#include "tcltk_bmodel.h"

#include "file_util.h"
#include "timer.h"
#include "utilities.h"

// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string	command;		// Command line

// Globals necessary to make the Tcl/Tk interface work
static Tcl_ObjCmdProc	systemMemoryCmd;
static Tcl_ObjCmdProc	bitSizeCmd;
static Tcl_ObjCmdProc	fileTypeCmd;
static Tcl_ObjCmdProc	spectrumCmd;
static Tcl_ObjCmdProc 	BimageCmd;
static Tcl_ObjCmdProc 	BmgCmd;
static Tcl_ObjCmdProc 	BmodelCmd;
Bimage* 				imglist = NULL;
Bproject* 				project = NULL;
Bmodel* 				model = NULL;

Bimage*					imgtemp = NULL;
		
/* These functions are dynamically loaded by the "load" command */
extern "C" int Bshow_Init(Tcl_Interp *interp)
{
	Tcl_CreateObjCommand(interp, "systemMemory", systemMemoryCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "bitSize", bitSizeCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "fileType", fileTypeCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "spectrum", spectrumCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bimage", BimageCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bmg", BmgCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bmodel", BmodelCmd, NULL, NULL);
	return TCL_OK;
}

extern "C" int Bshow_SafeInit(Tcl_Interp *interp)
{
	Tcl_CreateObjCommand(interp, "systemMemory", systemMemoryCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "bitSize", bitSizeCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "fileType", fileTypeCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "spectrum", spectrumCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bimage", BimageCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bmg", BmgCmd, NULL, NULL);
	Tcl_CreateObjCommand(interp, "Bmodel", BmodelCmd, NULL, NULL);
	return TCL_OK;
}

static int systemMemoryCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_ResetResult(interp);

	Tcl_Obj*			returnObj = Tcl_NewLongObj(system_memory());
	
	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}

static int bitSizeCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_ResetResult(interp);

	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	int					bits = 8*sizeof(long);
	
	Tcl_SetIntObj(returnObj, bits);
	
	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}

static int fileTypeCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_ResetResult(interp);

	char*				filename = Tcl_GetStringFromObj(objv[1], NULL);

	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	int					ft = file_type(filename);
	
//	cout << "File type = " << ft << endl;

	Tcl_SetIntObj(returnObj, ft);
	
	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}

static int spectrumCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	Tcl_ResetResult(interp);

	Tcl_Obj*			returnObj = Tcl_NewObj();
	
	double				value(0), min(0), max(1);
	RGB<unsigned char>	color;
	
	if ( objc > 1 ) Tcl_GetDoubleFromObj(NULL, objv[1], &value);
	if ( objc > 2 ) Tcl_GetDoubleFromObj(NULL, objv[2], &min);
	if ( objc > 3 ) Tcl_GetDoubleFromObj(NULL, objv[3], &max);
	
	color.spectrum(value, min, max);
	
	Bstring				hexcolor(color.hex());
	Tcl_SetStringObj(returnObj, hexcolor.c_str(), hexcolor.length());

	Tcl_SetObjResult(interp, returnObj);
	
	return TCL_OK;
}


/**
@brief 	Implements the "Bimage" command in Tcl/Tk to access image files through Bsoft.
@param 	clientdata	some pointer to some data structure (not used).
@param 	*interp		a Tcl interpreter within Tcl.
@param 	objc		number of arguments passed (+1).
@param 	*objv[]		arguments passed as Tcl objects.
@return static int		Tcl result.

	Calls the "image_processing" function.

**/
static int BimageCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	return image_processing(interp, objc, objv);
}

/**
@brief 	Implements the "Bmg" command in Tcl/Tk to access micrograph parameter files through Bsoft.
@param 	clientdata	some pointer to some data structure (not used).
@param 	*interp		a Tcl interpreter within Tcl.
@param 	objc		number of arguments passed (+1).
@param 	*objv[]		arguments passed as Tcl objects.
@return static int		Tcl result.

	Calls the "project_processing" function.

**/
static int BmgCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	return project_processing(interp, objc, objv);
}

/**
@brief 	Implements the "Bmodel" command in Tcl/Tk to access model parameter files through Bsoft.
@param 	clientdata	some pointer to some data structure (not used).
@param 	*interp		a Tcl interpreter within Tcl.
@param 	objc		number of arguments passed (+1).
@param 	*objv[]		arguments passed as Tcl objects.
@return static int		Tcl result.

	Calls the "model_processing" function.

**/
static int BmodelCmd(ClientData clientdata, Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
	return model_processing(interp, objc, objv);
}
