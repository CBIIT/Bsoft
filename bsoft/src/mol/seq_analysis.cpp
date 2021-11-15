/**
@file	seq_analysis.cpp
@brief	Analysis of protein sequences 
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20210426
**/
 
#include "seq_analysis.h"
#include "seq_util.h" 
#include "ps_sequence.h"
#include "moving_average.h"
#include "linked_list.h"
#include "utilities.h" 
 
// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


double			log_2(double a)
{
	if ( a <= 0 ) return 0;
	else return log(a)/log(2.0L);
}

/**
@brief 	Limits the selection to the reference sequence in an aligned set.
@param 	*molgroup 		the set of sequences.
@param 	&refseq 			reference sequence identifier.
@return long 				number of selected residues.

**/
long		seq_limit(Bmolgroup* molgroup, Bstring& refseq)
{
	long			i, n(0);
	Bmolecule*		mol;
	
	for ( mol = molgroup->mol; mol; mol = mol->next )
		if ( mol->id.contains(refseq) ) break;
		
	if ( !mol ) {
		cerr << "Sequence " << refseq << " not found!" << endl;
		return molgroup->maxlen;
	}
	
	for ( i=0; i<mol->seq.length(); ++i ) {
		if ( mol->seq[i] == '-' ) {
			molgroup->seqflag[i] = 0;
		} else {
			molgroup->seqflag[i] = 1;
			n++;
		}
	}
	
	return n;
}

/**
@brief 	Calculates the pairwise identities between aligned sequences.
@param 	*molgroup 		the set of sequences.
@return Matrix 			the matrix of identities.

The identity between two sequences is defined as:
			   number of identical residues
	identity = ----------------------------
					   overlap
where the overlap is the number of positions with residues in both sequences.

**/
Matrix	 	seq_aligned_identity(Bmolgroup* molgroup)
{
	long   			i, j, k, n(0), nid, overlap, nmol(0), maxlen;
	long			idsum(0), idssum(0), overlapsum(0), overlapssum(0);
	double			idavg, overlapavg;
	Bmolecule*		mol1;
	Bmolecule*		mol2;
	
	if ( verbose & VERB_LABEL )
	    cout << "Aligned identity analysis:" << endl;

	if ( verbose & VERB_PROCESS )
		cout << "Seq1\tSeq2\tIdentity\tnID\tOverlap\tName1\tName2" << endl;
	
	for ( mol1 = molgroup->mol; mol1; mol1 = mol1->next ) nmol++;
	
	// Initialize an image structure to hold the results
	Matrix		mat(nmol, nmol);
	
	for ( i=1, mol1 = molgroup->mol->next; mol1; i++, mol1 = mol1->next ) {
		for ( j=0, mol2 = molgroup->mol; mol2 != mol1 && j<i; j++, mol2 = mol2->next ) {
			n++;
			nid = overlap = 0;
			maxlen = mol1->seq.length();
			if ( maxlen > mol2->seq.length() ) maxlen = mol2->seq.length();
			for ( k=0; k<maxlen; k++ ) {
				if ( molgroup->seqflag[k] > 0 &&
						mol1->seq[k] != '-' && 
						mol2->seq[k] != '-' ) {
					overlap++;
					if ( mol1->seq[k] == mol2->seq[k] ) nid++;
				}
			}
			idsum += nid;
			idssum += nid*nid;
			overlapsum += overlap;
			overlapssum += overlap*overlap;
//			k = i*nmol + j;
//			if ( overlap ) pimg->set(k, nid*1.0/overlap);
			if ( overlap ) mat[i][j] = mat[j][i] = nid*1.0/overlap;
			if ( verbose & VERB_PROCESS )
				cout << i+1 << tab << j+1 << tab << fixed << setprecision(3) << setw(8) <<
					mat[i][j] << tab << nid << tab << overlap << tab <<
					mol1->id << tab << mol2->id << endl;
		}
	}
	
	idavg = idsum*1.0/n;
	overlapavg = overlapsum*1.0/n;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Average identical residues:  " << idavg << " (" << 
			sqrt(idssum*1.0/n - idavg*idavg) << ")" << endl;
		cout << "Average overlap:             " << overlapavg << " (" << 
			sqrt(overlapssum*1.0/n - overlapavg*overlapavg) << ")" << endl << endl;
	}
	
    return mat;
}

/**
@brief 	Calculates the pairwise similarities between aligned sequences.
@param 	*molgroup 		the set of sequences.
@param 	threshold 		threshold to accept residues as similar.
@param 	*simat			residue similarity matrix.
@return Matrix			the matrix of similarities.

The similarity between two sequences is defined as:
				   sum(residue similarity)
	similarity   = -----------------------
						  overlap
						  number of residues with similarity > threshold
	fraction similarity = ----------------------------------------------
											overlap
where the overlap is the number of positions with residues in both sequences.
The residue similarity is taken from a residue substitution matrix.
The default substitution matrix is BLOSUM62.

**/
Matrix	 	seq_aligned_similarity(Bmolgroup* molgroup, double threshold, Bresidue_matrix* simat)
{
    long 				i, j, k, thelen, overlap, nmol(0);
    long				ii, jj;
	double				simsum(0), percentsim(0);
	Bmolecule			*mol1, *mol2;
	
	long				nres = simat->n;
	Bstring				code1 = simat->c;
	float*				sim_mat = simat->m;
    
	for ( mol1 = molgroup->mol; mol1; mol1 = mol1->next ) nmol++;
	
	if ( verbose & VERB_LABEL ) {
	    cout << "Aligned similarity analysis:" << endl;
	    cout << "Similar residue threshold:      " << threshold << endl;
	    cout << "Number of sequences:            " << nmol << endl << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Seq1\tSeq2\tSimilarity\tFraction\tOverlap\tName1\tName2" << endl;
	
	// Initialize an image structure to hold the results
	Matrix		mat(nmol, nmol);
	
	for ( i=1, mol1 = molgroup->mol->next; mol1; i++, mol1 = mol1->next ) {
//		seq1 = mol1->seq.c_str();
		for ( j=0, mol2 = molgroup->mol; mol2 != mol1; j++, mol2 = mol2->next ) {
//			seq2 = mol2->seq.c_str();
//			thelen = strlen(seq1);
//          if ( thelen > strlen(seq2) ) thelen = strlen(seq2);
            thelen = mol1->seq.length();
            if ( thelen > mol2->seq.length() ) thelen = mol2->seq.length();
            simsum = percentsim = 0.0;
			overlap = 0;
			if ( verbose & VERB_DEBUG )
				cout << "i=" << i << " j=" << j << " thelen=" << thelen << endl;
			for ( k=0; k<thelen; k++ ) {
				if ( molgroup->seqflag[k] > 0 &&
						( mol1->seq[k] != '-' ) && ( mol2->seq[k] != '-' ) ) {
//					ii = find_char_index(code1, mol1->seq[k]);
//					jj = find_char_index(code1, mol2->seq[k]);
					ii = code1.index(mol1->seq[k]);
					jj = code1.index(mol2->seq[k]);
					if ( ii > -1 && jj > -1 ) {
	                	simsum += sim_mat[ii*nres+jj];
						if ( sim_mat[ii*nres+jj] >= threshold ) 
							percentsim += 1;
						overlap++;
					}
				}
				if ( verbose & VERB_DEBUG )
					cout << k << tab << simsum << tab << overlap << endl;
			}
//			if ( overlap ) {
//				pimg->set(i*nmol+j, simsum*1.0/overlap);
//				percentsim /= overlap;
//			} else pimg->set(i*nmol+j, 0);
			if ( overlap ) {
				mat[i][j] = simsum*1.0/overlap;
				percentsim /= overlap;
			} else mat[i][j] = 0;
			mat[j][i] = mat[i][j];
			if ( verbose & VERB_PROCESS )
				cout << i+1 << tab << j+1 << tab << fixed << setprecision(3) << setw(8) <<
					mat[i][j] << tab << setw(8) << percentsim << tab << overlap << tab <<
					mol1->id << tab << mol2->id << endl;
        }
	}
	cout << endl;
    
    return mat;
}

/**
@brief 	Selects sequences within a range of lengths.
@param 	*molgroup 		the set of sequences.
@param 	minlen 			minimum length.
@param 	maxlen 			maximum length.
@return long				number of sequences retained.
**/
long		seq_select(Bmolgroup* molgroup, long minlen, long maxlen)
{
	long			i, n(0), len;
	Bmolecule*		mol;

	if ( verbose )
		cout << "Selecting sequences of length " << minlen << " - " << maxlen << endl;

	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( i=len=0; i<mol->seq.length(); ++i )
			if ( mol->seq[i] != '-' ) len++;
		if ( len >= minlen && len <=maxlen ) {
			mol->sel = 1;
			n++;
		} else {
			mol->sel = 0;
		}
	}

	if ( verbose )
		cout << "Selected:                  " << n << endl << endl;

	return n;
}

/**
@brief 	Selects sequences based on a comparison matrix of aligned sequences.
@param 	*molgroup 		the set of sequences.
@param 	mat 				comparison matrix.
@param 	ref		 		reference sequence number (starting at 1).
@param 	cutoff	 		threshold for selecting sequences.
@return long				number of sequences retained.
**/
long		seq_select(Bmolgroup* molgroup, Matrix mat, long ref, double cutoff)
{
	if ( ref > mat.rows() ) return mat.rows();
	
	long			i(ref-1), j, n(0);
	Bmolecule*		mol;

	for ( j = 1, mol = molgroup->mol; j < ref && mol; ++j ) ;
	
	if ( verbose ) {
		cout << "Selecting based on sequence     " << ref << ":" << endl;
		cout << "Sequence ID:                    " << mol->id << endl;
		cout << "Cutoff:                         " << cutoff << endl << endl;
	}
	
	for ( j = 0, mol = molgroup->mol; j < mat.rows() && mol; mol = mol->next, ++j )
		if ( i != j && mat[i][j] < cutoff ) {
			mol->sel = 0;
		} else {
			mol->sel = 1;
			n++;
		}

	if ( verbose )
		cout << "Selected:                  " << n << endl << endl;
	
	return n;
}
/*
long		seq_select(Bmolgroup* molgroup, Matrix mat, long ref, double cutoff)
	{
		if ( ref > mat.rows() ) return -1;
		
		long			i(ref-1), j;
		Bmolecule*		mol, *pmol;

		for ( j = 1, mol = molgroup->mol; j < ref && mol; ++j ) ;
		
		if ( verbose ) {
			cout << "Selecting based on sequence     " << ref << ":" << endl;
			cout << "Sequence ID:                    " << mol->id << endl;
			cout << "Cutoff:                         " << cutoff << endl << endl;
		}
		
	//	cout << mat << endl;

		for ( j = 0, mol = molgroup->mol; j < mat.rows() && mol; mol = mol->next, ++j )
			if ( i != j && mat[i][j] < cutoff ) mol->sel = 0;
			else mol->sel = 1;

		for ( j = 0, mol = molgroup->mol; mol; mol = mol->next ) if ( mol->sel ) ++j;
		
	//	cout << "reference = " << i << endl;
	//	cout << "selected = " << j << endl;

		for ( j = 0; j < mat.rows();  ) {
			if ( i != j && mat[i][j] < cutoff ) {
				mat = mat.delete_row_column(j);
				if ( i > j ) i--;
			} else j++;
		}

		for ( mol = pmol = molgroup->mol; mol; ) {
			if ( mol->sel < 1 ) {
				if ( mol == molgroup->mol ) {
					molgroup->mol = pmol = mol->next;
					molecule_kill(mol);
					mol = pmol;
				} else {
					pmol->next = mol->next;
					molecule_kill(mol);
					mol = pmol->next;
				}
			} else {
				pmol = mol;
				mol = mol->next;
			}
		}

		for ( j = 0, mol = molgroup->mol; mol; mol = mol->next ) ++j;

		if ( verbose ) {
			cout << "Number of sequences retained:   " << j << endl << endl;
			cout << "Number of sequences retained:   " << mat.rows() << endl << endl;
		}
		
		return mat.rows();
	}
*/

/**
@brief 	Deletes non-selected sequences and corresponding elelments of a comparison matrix.
@param 	*molgroup 		the set of sequences.
@param 	mat 				comparison matrix.
@return long				number of sequences retained.
**/
long		seq_delete(Bmolgroup* molgroup, Matrix mat)
{
	long			i;
	Bmolecule*		mol, *pmol;

	for ( i = 0, mol = molgroup->mol; mol; mol = mol->next )
		if ( mol->sel == 0 ) ++i;

	if ( verbose )
		cout << "Deleting " << i << " sequences" << endl << endl;
	
//	cout << mat << endl;

	if ( mat.rows() ) {
		vector<long>	del(mat.rows(), 0);
		for ( i = 0, mol = molgroup->mol; i < mat.rows() && mol; mol = mol->next, ++i )
			if ( mol->sel == 0 ) del[i] = 1;
		for ( i = mat.rows() - 1; i >= 0; --i )
			if ( del[i] ) mat = mat.delete_row_column(i);
	}

	for ( mol = pmol = molgroup->mol; mol; ) {
		if ( mol->sel < 1 ) {
			if ( mol == molgroup->mol ) {
				molgroup->mol = pmol = mol->next;
				molecule_kill(mol);
				mol = pmol;
			} else {
				pmol->next = mol->next;
				molecule_kill(mol);
				mol = pmol->next;
			}
		} else {
			pmol = mol;
			mol = mol->next;
		}
	}

	for ( i = 0, mol = molgroup->mol; mol; mol = mol->next ) ++i;

	if ( verbose ) {
		cout << "Number of sequences retained:   " << i << endl << endl;
		cout << "Number of sequences retained:   " << mat.rows() << endl << endl;
	}
	
	return mat.rows();
}

/**
@brief 	Generates a PROSITE format profile from an aligned set of sequences.
@param 	*molgroup 		the set of sequences.
@return string			profile in PROSITE format.

At each position in the alignment, the number of distinct residue types
are counted. If there are more than 3 residue types represented at a
position, or there is a gap, it is designated as variable by an "x".
The profile finally contains 1-3 residue type possibilities for highly
conserved positions interspersed by variable length gaps.

**/
string		seq_aligned_profile(Bmolgroup* molgroup)
{
	bool			notdone;
    long 			i, j, k, x, g;
	Bmolecule		*mol;
	char			reslist[3];
	
	long			len(3*molgroup->maxlen);
	string			profseq(len,' ');
	string			profile;
	
	if ( verbose & VERB_LABEL )
	    cout << "Profile:" << endl;
	
	for ( i=k=0; i<molgroup->maxlen; i++, k+=3 ) {
		notdone = 1;
		for ( j=0; j<3; j++ ) reslist[j] = ' '; 
		for ( mol = molgroup->mol; mol; mol = mol->next ) {
			if ( mol->seq[i] == '-' ) {
				for ( j=k; j<k+3; j++ ) profseq[j] = '-';
				notdone = 0;
			} else {
				for ( j=0; j<3 && reslist[j] != ' ' && mol->seq[i] != reslist[j]; j++ ) ;
				if ( j<3 ) reslist[j] = mol->seq[i];
				else {
					for ( j=k; j<k+3; j++ ) profseq[j] = 'x';
					notdone = 0;
				}
			}
			if ( notdone ) for ( j=0; j<3; j++ ) profseq[k+j] = reslist[j];
		}
	}
	
	for ( i=g=x=k=0; i<molgroup->maxlen; i++, k+=3 ) {
		if ( profseq[k] == '-' ) g++;
		else if ( profseq[k] == 'x' ) x++;
		else {
			if ( profile.length() ) profile += '-';
			if ( g + x > 0 ) {
				if ( g == 0 ) {
					profile += "x(" + to_string(x) + ")-";
				} else {
					profile += "x(" + to_string(x) + "," + to_string(g) + ")-";
				}
			}
			if ( profseq[k+1] == ' ' ) {
				profile += profseq[k];
			} else {
				profile += '[';
				for ( j=0; j<3 && profseq[k+j] != ' '; j++ ) profile += profseq[k+j];
				profile += ']';
			}
			g = x = 0;
		}
	}
	
	if ( verbose )
		cout << profile << endl;
	
    return profile;
}

/**
@brief 	Calculates the sequence logo representation for an alignment.
@param 	*molgroup 		the set of sequences.
@param 	window			window for calculating the moving average.
@param 	&psfile			the postscript file name.
@return int 				0.

The information content of each position in an alignment is calculated
as:
	information = log_2(n) - sum(pi * log_2(pi) )
					 fi
	pi          =  -------
				   sum(fi)
	fi          =  frequency of residue type i at this position
	n           =  sum(fi) if sum(fi) < 20, otherwise n = 20
A moving average of the information is calculated over a given window
to smooth the resultant data.
The sequence logo representation for the occurrence of every residue type
at every position is generated and written into a postscript file.

**/
int 	 	seq_aligned_information(Bmolgroup* molgroup, int window, Bstring& psfile)
{
	long			i, j, k, ii, n, nres;
    double       	fgap;
	
	Bstring			temp;
	Bresidue_type*	rt = get_residue_properties(temp);
	Bresidue_type*	rt_curr = NULL;
	for ( nres = 0, rt_curr = rt; rt_curr; rt_curr = rt_curr->next ) nres++;
	string			code1(nres, ' ');
	for ( i = 0, rt_curr = rt; rt_curr; rt_curr = rt_curr->next, i++ )
		code1[i] = rt_curr->c;
	
	kill_list((char *) rt, sizeof(Bresidue_type));
    
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_aligned_information: code=" << code1 << endl;
    
	long			nmol(0);
	Bmolecule*		mol = molgroup->mol;
	for ( nmol = 0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	long			npos(0);
	for ( i=0; i<molgroup->maxlen; i++ ) npos += molgroup->seqflag[i];
	
	long   			igap = code1.find_first_of('-');
	vector<double>	fsum(npos,0);
	vector<double>	info(npos,0);
    vector<double>	freq(nres*npos);
	string			pattern(nres*npos,' ');
	string			seq(npos,' ');
	
    if ( verbose & VERB_LABEL ) {
		cout << "Information content analysis:" << endl;
		cout << "Sequences:                      " << nmol << endl;
		cout << "Alignment positions:            " << npos << endl;
		cout << "Moving average window:          " << window << endl;
	}
	
    // Calculating the frequency of each residue type occurring at each
    // position in the alignment
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
    	j = k = 0;
    	while ( j < mol->seq.length() ) {
			if ( molgroup->seqflag[j] ) {
				ii = code1.find_first_of(mol->seq[j]);
    	        if ( ii > -1 ) freq[nres*k+ii] += 1;
				k++;
			}
    	    j++;
    	}
    }
	
    // Calculating information content for a position in the alignment
	multimap<double,char> 	ri;
	double					v;
	k = 0;
	mol = molgroup->mol;
    for ( j=0; j<molgroup->maxlen; j++ ) {
        if ( molgroup->seqflag[j] ) {
	    	fgap = 0;
	    	ri.clear();
    		for ( i=0; i<nres; i++ ) {
				if ( i == igap )
					fgap += freq[nres*k+i];
				else
    		    	fsum[k] += freq[nres*k+i];
			}
			info[k] = 0;
			if ( fsum[k] > 0 ) {
    			for ( i=0; i<nres; i++ ) if ( i != igap ) {
					v = freq[nres*k+i]*log_2(freq[nres*k+i])/fsum[k];
    	    		ri.insert(pair<double,char>(v,code1[i]));
    	    		info[k] += v;
				}
				if ( fsum[k] > 20 )
					info[k] += log(20.0/fsum[k])/log(2.0);
			}
			i = 0;
			for ( auto it: ri ) {
				pattern[nres*k+i] = it.second;
				freq[nres*k+i] = it.first;
				i++;
			}
			seq[k] = mol->seq[j];
			k++;
		}
    }
	
	// Calculate a moving average
	vector<double>		movavg = moving_average(info, window);
	
	// Print out the list of results
    if ( verbose & VERB_LABEL )
		cout << "  Pos\tBits\tMovAvg\tnSeq\tPattern" << endl;
    for ( k=0; k<npos; k++ ) {
    	cout << k+1 << tab << info[k] << tab << movavg[k] << tab << fsum[k] << tab;
		for ( i=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 )
			cout << pattern[nres*k+i];
		cout << endl;
    }
	
    for ( k=0; k<npos; k++ ) {
		if ( k > 0 ) cout << "-";
		for ( i=n=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 ) n++;
		if ( n < 1 ) {
			cout << "x";
		} else {
			if ( n > 1 ) cout << "[";
			for ( i=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 )
				cout << pattern[nres*k+i];
			if ( n > 1 ) cout << "]";
		}
    }
	cout << endl;
	
	vector<Complex<float>>	per = seq_frequency_analysis(16, 4, 4, info);
	
	vector<Complex<float>>	per_avg = moving_average_complex(per, window);
	
	Bstring				title("Information");
	if ( psfile.length() )
		ps_seq_info(psfile, title, nres, movavg, fsum, per_avg,
			freq, pattern);
	
	return 0;
}

/**
@brief 	Calculates the average hydrophobicity at every position in an alignment.
@param 	*molgroup 		the set of sequences.
@param 	window			moving average window.
@param 	threshold		fraction of sequences with a residue in a position.
@param 	&hphobfile		parameter file.
@param 	&psfile			postscript output file.
@return int 	 			0.

The default hydrophobicity scale is the GES scale.

**/
int 	 	seq_aligned_hydrophobicity(Bmolgroup* molgroup,
				int window, double threshold, Bstring& hphobfile, Bstring& psfile)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_aligned_hydrophobicity: " << hphobfile << endl;
	
    long   			i, j, nres;
	Bresidue_type*	rt = get_residue_properties(hphobfile);
	Bresidue_type*	rt_curr = NULL;
	for ( nres = 0, rt_curr = rt; rt_curr; rt_curr = rt_curr->next ) nres++;
	string			code1(nres,'*');
	for ( i = 0, rt_curr = rt; rt_curr; rt_curr = rt_curr->next, i++ )
		code1[i] = rt_curr->c;
	
    if ( verbose & VERB_LABEL )
		cout << "Average aligned hydrophobicity:" << endl;
    
	vector<double>	nseq(molgroup->maxlen,0);
	vector<double>	avgHP(molgroup->maxlen,0);
	vector<int>		HPseg(molgroup->maxlen,0);
	Bmolecule*		mol;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
        for ( j=0; j<mol->seq.length(); j++ ) {
			if ( mol->seq[j] != '-' ) {
				for ( rt_curr = rt; rt_curr && mol->seq[j] != rt_curr->c; rt_curr = rt_curr->next ) ;
				if ( mol->seq[j] == rt_curr->c ) {
					avgHP[j] += -rt_curr->hphob; 	// Note: Negative for scales of hydrophilicity
					nseq[j]++;
				}
			}
        }
    }

	kill_list((char *) rt, sizeof(Bresidue_type));
	
    for ( j=0; j<molgroup->maxlen; j++ )
        if ( nseq[j] > 0 ) avgHP[j] /= nseq[j];
    
    if ( verbose & VERB_PROCESS )
		cout << "nRes\tavgHP\tintHP\tnSeq\tRes" << endl;
	
	double		sum(0);
	i = 0;
	mol = molgroup->mol;
    for ( j=0; j<molgroup->maxlen; j++ ) {
        if ( molgroup->seqflag[j] ) {
			i++;
			sum += avgHP[j];
			cout << i << tab << avgHP[j] << tab << sum << tab << nseq[j] << tab << mol->seq[j] << endl;
		}
	}
	
    if ( verbose & VERB_PROCESS )
		cout << endl << "Segment\tStart\tEnd\tLength\tsumHP" << endl;
	
	int 		k, m, sum_flag(0), start[1000], end[1000], newstart, newend, nmol(0);
	double		max, sumHP[1000];
	
	for ( nmol = 0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( threshold < 0.5 ) threshold = 0.5;
	if ( threshold > 1 ) threshold = 1;
	threshold *= nmol;
	i = 0;
	sumHP[0] = 0;
    for ( j=1; j<molgroup->maxlen; j++ ) {
        if ( molgroup->seqflag[j] ) {
			if ( nseq[j] >= threshold && nseq[j-1] < threshold ) {
				if ( sumHP[i] > 20 ) i++;
				start[i] = j;
				sum_flag = 1;
				sumHP[i] = 0;
			}
			if ( sum_flag && nseq[j] < threshold ) {
				end[i] = j - 1;
				sum_flag = 0;
				newstart = start[i];
				newend = end[i];
				max = sumHP[i];
				for ( k=start[i]; k<end[i]; k++ ) {
					sum = 0;
					for ( m=k; m<=end[i]; m++ ) {
						sum += avgHP[m];
						if ( max < sum ) {
							max = sum;
							newstart = k;
							newend = m;
						}
					}
				}
				start[i] = newstart;
				end[i] = newend;
				sumHP[i] = max;
				if ( max > 20 ) {
					for ( k=start[i]; k<=end[i]; k++ ) HPseg[k] = 1;
					cout << i+1 << tab << start[i]+1 << tab << end[i]+1 << tab << 
						end[i] - start[i] + 1 << tab << sumHP[i] << endl;
				}
			}
			if ( sum_flag ) sumHP[i] += avgHP[j];
		}
	}
	cout << endl;
	
	// Calculate a moving average
	vector<double>			movavg = moving_average(avgHP, window);
	
	vector<Complex<float>>	per = seq_frequency_analysis(16, 4, 4, avgHP);
	
	vector<Complex<float>>	per_avg = moving_average_complex(per, window);
	
	Bstring				title("Hydrophobicity");
	if ( psfile.length() )
		ps_seq_hydrophob(psfile, title, movavg, HPseg, nseq, per_avg);
	
	return 0;
}

/**
@brief 	Fourier transforms a vector for frequency analysis.
@param 	win				window size.
@param 	start			start within window.
@param 	end				end within window.
@param 	*data			sequence.
@return int 	 			0.

A brute force Fourier transform is done.

**/
vector<Complex<float>>	seq_frequency_analysis(long win, long start, long end, vector<double>& data)
{
	if ( start > end ) swap(start, end);
	if ( start < 0 ) start = 0;
	if ( end > win - 1 ) end = win - 1;
	
    if ( verbose )
		cout << "Frequency analysis:" << endl;
       
	if ( verbose & VERB_RESULT ) {
		cout << "Window:                         " << win << endl;
		cout << "Start - end:                    " << start << " - " << end << endl;
	}
	
    long         	j, k, l, length(data.size());
	long			half_win = win/2;
	double			real, imag, amp, phi, dphi;
	
    vector<Complex<float>>	transform_data(length*(end-start+1));
	
    for ( j=0; j<length-win; j++ ) {
		cout << j+half_win;
		for ( k=start; k<=end; k++ ) {
			real = 0;
			imag = 0;
            for ( l=0; l<win; l++ ) {
				phi = MIN2PI*k*l/win;
                real += data[l+j]*cos(phi);
                imag += data[l+j]*sin(phi);
			}
			amp = sqrt(real*real + imag*imag);
            phi = atan2(imag, real);
			if ( amp < 0.001 ) phi = 0;
            dphi = phi - TWOPI*j*k/win;
			dphi = angle_set_negPI_to_PI(dphi);
			transform_data[length*(k-start)+j+half_win] = Complex<float>(amp*cos(dphi), amp*sin(dphi));
			if ( verbose & VERB_RESULT )
				cout << tab << amp << tab << dphi*180/M_PI;
		}
		if ( verbose & VERB_RESULT )
			cout << endl;
    }
	cout << endl;
	
    return transform_data;
}

vector<double>		seq_aligned_weight(Bmolgroup* molgroup)
{
	char			m1, m2, igap = '-';
	long   			i, k, l, nseq, overlap, identity;
	
	Bmolecule		*mol1, *mol2;
	
	for ( nseq=0, mol1 = molgroup->mol; mol1; mol1 = mol1->next ) nseq++;
	
	vector<double>	weight(nseq*nseq, 0);
	
	// Weighting is based on the difference between identity and overlap
	for ( k=1, mol1 = molgroup->mol->next; k<nseq; k++, mol1 = mol1->next ) {
		for ( l=0, mol2 = molgroup->mol; l<k; l++, mol2 = mol2->next ) {
			overlap = identity = 0;
			for ( i=0; i<molgroup->maxlen; i++ ) {
				m1 = mol1->seq[i];
				m2 = mol2->seq[i];
				if ( m1 != igap && m2 != igap ) {
					overlap++;
					if ( m1 == m2 ) identity++;
				}
			}
			weight[k*nseq+l] = (overlap - identity)*1.0/molgroup->maxlen;
			if  ( verbose & VERB_DEBUG )
				cout << "k=" << k << " l=" << l << " weight=" << weight[k*nseq+l] << endl;
		}
	}
	
	return weight;
}

/**
@brief 	Correlated mutation analysis of an alignment.
@param 	*molgroup 		the set of aligned sequences.
@param 	refseqid		reference sequence to report on.
@param 	cutoff			cutoff for reporting correlated mutations.
@param 	&simfile		similarity matrix file.
@return Matrix 			the analysis result matrix.

Reference: Gobel, Sander & Schneider (1994) Proteins 18, 309-317.
Mutation (residue variation) correlation is defined as:
					1
	r(i,j) =  ------------- sum(w(k,l)*(s(i,k,l) - <s(i)>)*(s(j,k,l) - <s(j)>))
			  m^2*o(i)*o(j)
	where:
		m:         number of sequences
		o(i):      standard deviation of similarities at alignment position i
		w(k,l):    weight for sequences k and l
				   (1 - fractional identity: see function seq_aligned_identity)
		s(i,k,l):  similarity for alignment position i between sequences k and l
		<s(i)>:    average similarity at alignment position i
Individual high-scoring correlations (using the given cutoff value) are reported
as follows:
	Res1	Num1	Res2	Num2	Total	Corr
	T	9	I	17	210	 0.631
	TAIIIVVVIVVVIVIIIIIII
	IILLLLLLLLLLLLLLLLLLL
The first 4 values gives the type and alignment position of the correlating residues.
The total is the number of comparisons made: maximally m*(m-1)/2
The last number is the correlation coefficient.
The following two lines gives the corresponding residues at the two alignment positions
for all the sequences, allowing the user to see on what basis this is a high correlation.

**/
Matrix		seq_correlated_mutation(Bmolgroup* molgroup,
					Bstring& refseqid, double cutoff, Bstring& simfile)
{
	long 				h, i, j, k, l;
	long				refseqnum, length(molgroup->maxlen), number(0);
	
	Bmolecule			*molref, *mol1, *mol2;
	for ( number=0, mol1 = molgroup->mol; mol1; mol1 = mol1->next ) number++;
	for ( refseqnum=0, molref = molgroup->mol; molref; molref = molref->next, ++refseqnum )
		if ( molref->id.contains(refseqid) ) break;
	
	if ( !molref ) {
		cerr << "Sequence " << refseqid << " not found!" << endl;
		return Matrix();
	}
	
	Bstring				refseq(molref->seq);
	vector<char>		selseq(length*number, 0);

	if  ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_correlated_mutation: length=" << length << " number=" << number << endl;
		
	int 				m1, m2, ir, ntot;
	int 				total = number*(number-1)/2;		// Only lower triangle
	double				sum, ssum, sigma2, maxCC(-1e30);
	vector<double>		avg(length, 0);
	vector<double>		std(length, 0);
	vector<double>		num(length, 0);
	vector<double>		score(length*number*number, 0);
	
	Bresidue_matrix*	simat = get_residue_matrix(simfile);

	int					nres = simat->n;
	Bstring				code1 = simat->c;
	float*				sim_mat = simat->m;
    
	int 				igap = code1.index('-');
	
	if  ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_correlated_mutation: total=" << total << endl;
	
	vector<double>		weight = seq_aligned_weight(molgroup);
	
	for ( i=ir=0; i<length; i++ ) {
		if ( molgroup->seqflag[i] ) {
			sum = ssum = 0;
			m1 = code1.index(molref->seq[i]);
			for ( k=1, mol1 = molgroup->mol->next; k<number; k++, mol1 = mol1->next ) {
				m1 = code1.index(mol1->seq[i]);
				for ( l=0, mol2 = molgroup->mol; l<k; l++, mol2 = mol2->next ) {
					m2 = code1.index(mol2->seq[i]);
					h = (ir*number + k)*number + l;
						// Exclude all gaps
					if ( m1 != igap && m2 != igap ) {
						score[h] = sim_mat[m2*nres + m1];
						sum += score[h];
						ssum += score[h]*score[h];
						num[ir]++;
					} else score[h] = -999;
//					if ( m1 == m2 ) score[h] = -999;
				}
			}
			avg[ir] = std[ir] = 0;
			if ( total ) {
				avg[ir] = sum/total;
				std[ir] = sqrt(ssum/total - avg[ir]*avg[ir]);
			}
			refseq[ir] = refseq[i];
			for ( j=0, mol1 = molgroup->mol; j<number; j++, mol1 = mol1->next )
				selseq[ir*number+j] = mol1->seq[i];
			if  ( verbose & VERB_DEBUG )
				cout << "i=" << ir << " avg=" << avg[ir] << " std=" << std[ir] << " num=" << num[ir] << endl;
			ir++;
		}
	}
	length = ir;

	if ( verbose & VERB_PROCESS ) {
		cout << "Correlated mutation analysis:" << endl;
		cout << "Matrix size:                    " << length << " x " << length << endl;
		cout << "Reference sequence:             " << molref->id << " (" << refseqnum + 1 << ")" << endl << endl;
		cout << endl << "Res1\tNum1\tRes2\tNum2\tTotal\tCorr" << endl;
	}

	Matrix				mat(length, length);
	
	ir = 0;
	sigma2 = 0;
	for ( i=0; i<length; i++ ) {
		for ( j=0; j<=i; j++ ) {
			ntot = 0;
			if ( num[i] > 2 && num[j] > 2 ) {
				if ( std[i] && std[j] ) {
					for ( k=1; k<number; k++ ) {
						for ( l=0; l<k; l++ ) {
							m1 = (i*number + k)*number + l;
							m2 = (j*number + k)*number + l;
							if ( score[m1] > -100 && score[m2] > -100 ) {
								mat[i][j] += weight[k*number+l]*
									(score[m1] - avg[i])*(score[m2] - avg[j]);
								ntot++;
							}
						}
					}

					mat[i][j] /= total*std[i]*std[j];
				}
			}
			if ( i == j ) mat[i][j] = 1;
			else if ( maxCC < mat[i][j] ) maxCC = mat[i][j];
			sigma2 += mat[i][j] * mat[i][j]; 		// Sum of squared correlations
			if ( mat[i][j] > cutoff && i != j ) {
				cout << refseq[j] << tab << j+1 << tab << 
					refseq[i] << tab << i+1 << tab << ntot << tab << mat[i][j] << endl;
				if ( verbose & VERB_PROCESS ) {
					for ( k=0; k<number; k++ )
						cout << selseq[j*number+k];
					cout << endl;
					for ( k=0; k<number; k++ )
						cout << selseq[i*number+k];
					cout << endl;
				}
				ir++;
			}
			mat[j][i] = mat[i][j];
		}
	}
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Correlations reported:           " << ir << endl;
		cout << "Maximum off-diagonal coefficient: " << maxCC << endl << endl;
	}
	
	residue_matrix_kill(simat);
	
	return mat;
}

