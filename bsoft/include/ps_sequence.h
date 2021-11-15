/**
@file	ps_sequence.h 
@brief	Header file for postscript output for sequence analysis functions.
@author Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20210426
**/
 
#include "Complex.h"
#include "utilities.h"

// Function prototypes
int 		ps_seq_info(Bstring& filename, Bstring& title, int nres, int length,
				double* info, double* nseq, Complex<float>* per, double* freq, char* pattern);
int 		ps_seq_info(Bstring& filename, Bstring& title, int nres,
				vector<double>& info, vector<double>& nseq,
				vector<Complex<float>>& per, vector<double>& freq, string& pattern);
int 		ps_seq_hydrophob(Bstring& filename, Bstring& title, int length,
				double* Hphob, int* HPseg, double* nseq, Complex<float>* per);
int 		ps_seq_hydrophob(Bstring& filename, Bstring& title,
				vector<double>& Hphob, vector<int>& HPseg,
				vector<double>& nseq, vector<Complex<float>>& per);

		
