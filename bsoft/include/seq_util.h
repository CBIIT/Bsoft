/**
@file	seq_util.h 
@brief	Header file for sequence utilities 
@author Bernard Heymann 
@date	Created: 20001029
@date	Modified: 20220713
**/
 
#include "Bstring.h"
#include "rwmolecule.h"

// Function prototypes 
int			seq_from_residues(Bmolgroup* molgroup);
int			seq_show(Bmolgroup* molgroup);
int			seq_mass(Bmolgroup* molgroup);
vector<double>	seq_elements(Bmolgroup* molgroup, Bstring& paramfile);
int 		seq_complement_all(Bmolgroup* molgroup);
int 		seq_translate_all(Bmolgroup* molgroup, int frame, Bstring& gcname);
long 		seq_find_dna(Bmolgroup* molgroup, Bstring& seq);
long 		seq_find_protein(Bmolgroup* molgroup, Bstring& seq);
Bstring		seq_find_protein_in_dna(Bmolgroup* molgroup, Bstring& seq, int seqlenmin, int seqlenmax, 
				int side1, int side2, double threshold, Bstring& gcfile);
int			getcode3(char c, char* cod);
char		getcode1(char* acode); 
int 		complement_sequence(Bstring& nucseq);
char		get_complement(char nuc);
Bstring 	sequence_translate(Bstring& nucseq, long frame, Bstring& gencode);

