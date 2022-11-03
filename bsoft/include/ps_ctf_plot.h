/**
@file	ps_ctf_plot.h 
@brief	Header file for postscript output functions.
@author Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20160229
**/
 
#include "ctf.h"
#include "Bplot.h"

// Function prototypes
int			ps_show_ctf_param(Bplot* plot, CTFparam& em_ctf);
int 		ps_ctf_plot(Bstring& filename, CTFparam& em_ctf, size_t n, double freq_step);
int 		ps_ctf_plot(long n, double* rps, double interval,
				CTFparam* em_ctf, Bstring& filename);
int 		ps_ctf_defocus_zeroes(Bstring& filename, double volts, double Cs, double amp_shift);
int 		ps_point_spread(Bstring& filename, CTFparam& em_ctf, size_t n, double freq_step);


