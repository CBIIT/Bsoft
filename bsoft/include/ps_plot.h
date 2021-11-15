/**
@file	ps_plot.h 
@brief	Header file for postscript output functions.
@author Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20200310
**/

#ifndef _ps_plot_

#include "Bplot.h"
#include "utilities.h"
#include <fstream>


// Function prototypes
ofstream*	ps_open_and_init(Bstring filename, Bstring title, int npages, 
					int width, int height);
ofstream*	ps_open_and_init(Bstring filename, Bplot* plot);
int			ps_close(ofstream* fps);
int			ps_plot(Bstring filename, Bplot* plot);
Bplot*		ps_read(Bstring& filename);
int 		ps_graph(ofstream* fps, Bplot* plot, int page_number);
int 		ps_scale(ofstream* fps, double x1, double y1, double x2, double y2,
					double min, double max, double increment,
					double tick_length, double tick_angle,
					int digits, int fontsize, int inverse);
Bplot*		plot_curve(unsigned long nrow, double* c0, double* c1, double* c2, double* c3);
int			ps_define_arrowline(ofstream* fps);

#define _ps_plot_
#endif
