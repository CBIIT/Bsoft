/**
@file	rwPostScript.h
@brief	Header file for (reading and) writing postscript image files
@author Bernard Heymann
@date	Created: 20010614
@date	Modified: 20220124

	Format: Postscript image file format
**/

#include "rwimg.h"

// I/O prototypes
int 		writePostScriptImage(Bimage* p);
int 		writePostScriptImage(Bimage* p, vector<Vector3<double>>& t);
