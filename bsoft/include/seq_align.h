/**
@file	seq_align.h
@brief	Header file for functions to generate and analyze dot plots
@author Bernard Heymann
@date	Created: 20001029
@date	Modified: 20110411
**/

#include "rwmolecule.h"
#include "rwresprop.h"
#include "Matrix.h"

// Function prototypes 
long		seq_pair_align(Bmolecule* mol1, Bmolecule* mol2, double gapopen, 
				double gapextend, Bresidue_matrix* simat);
int			seq_find_best_offset(Bmolecule* mol1, Bmolecule* mol2, 
				long& nres, Bresidue_matrix* simat);
Matrix		seq_dot_plot(Bmolecule* mol1, Bmolecule* mol2, Bresidue_matrix* simat);
Matrix		seq_dot_plot_mov_avg(Matrix dot_plot, int window);
int			seq_dot_plot_interpret(Matrix dot_plot);
double		seq_dot_plot_best_segments(Matrix dot_plot);
int			seq_show_segments(Bmolecule* mol1, Bmolecule* mol2, double threshold, Matrix dot_plot);

