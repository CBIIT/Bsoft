/**
@file	ps_marker.h
@brief	Header file for postscript tools dealing with Bmarkers
@author Samuel Payne and Bernard Heymann
@date	Created: 20010516
@date	Modified: 20200310 (BH)
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "marker.h"
#include "Transform.h"

//Function prototypes
int			ps_marker_plots(ofstream* fps, Bstring& title, Bproject* project);
//int			ps_marker_errors(Bstring& filename, Bstring& title, Bproject* project);
int			ps_marker_errors(Bstring& filename, Bmarker* set1, Bmarker* set2, 
				Transform t, Vector3<long> size, double err_scale);
int			ps_marker_errors(ofstream* fps, Bmarker* set1, Bmarker* set2, 
				Transform t, Vector3<long> size, double err_scale);
int 		ps_marker_match(Bmarker* set1, Bmarker* set2, Bstring& filename); 
int			ps_dual_zcompare(Bstring& filename, Bproject* project);
