/**
@file	ps_views.h
@brief	Header file for postscript tools dealing with views
@author Bernard Heymann
@date	Created: 20011127
@date	Modified: 20201125
**/

#include <fstream>
#include "Bstring.h"
#include "View.h"
#include "View2.h"
#include "Euler.h"
#include "symmetry.h"

//Function prototypes
int 		ps_views(Bstring& filename, Bstring& symmetry_string, View* view, int flags);
int 		ps_views(ofstream* fps, Bstring& symmetry_string, View* view, int flags);
int 		ps_views2(ofstream* fps, string& symmetry_string, list<View2<float>>& view, int flags);
int 		ps_views(Bstring& filename, View* view);
int 		ps_views(ofstream* fps, View* view);
int			ps_sets_of_views(Bstring& filename, Bstring& title, int nv, View* views, int ns, double* fom);
int			ps_phi_theta_plot(ofstream* fps, int left, int bottom, int width, int height, int ncol);


