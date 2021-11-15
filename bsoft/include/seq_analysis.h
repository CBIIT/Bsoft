/**
@file	seq_analysis.h 
@brief	Header file for sequence analysis functions
@author Bernard Heymann 
@date	Created: 19990123
@date	Modified: 20210426
**/
 
#include "rwmolecule.h"
#include "rwresprop.h"
#include "Complex.h"
#include "Matrix.h"
#include "utilities.h"
  
// Function prototypes 
long		seq_limit(Bmolgroup* molgroup, Bstring& refseq);
Matrix	 	seq_aligned_identity(Bmolgroup* molgroup);
Matrix	 	seq_aligned_similarity(Bmolgroup* molgroup, double threshold, Bresidue_matrix* simat);
long		seq_select(Bmolgroup* molgroup, long minlen, long maxlen);
long		seq_select(Bmolgroup* molgroup, Matrix mat, long ref, double cutoff);
long		seq_delete(Bmolgroup* molgroup, Matrix mat);
string		seq_aligned_profile(Bmolgroup* molgroup);
int 	 	seq_aligned_information(Bmolgroup* molgroup, int window, Bstring& psfile);
int 	 	seq_aligned_hydrophobicity(Bmolgroup* molgroup,
				int window, double threshold, Bstring& hphobfile, Bstring& psfile);
vector<Complex<float>>	seq_frequency_analysis(long win, long start, long end, vector<double>& data);
Matrix		seq_correlated_mutation(Bmolgroup* molgroup,
					Bstring& refseqid, double cutoff, Bstring& simfile);
