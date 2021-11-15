/**
@file	rwFASTA.cpp
@brief	Library routines to read and write FASTA sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20210422
**/

#include "rwFASTA.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads FASTA format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup molecule group.
@return int 				number of molecules read (<0 if reading failed).

	FASTA format:
	>1hiwa TRIMERIC HIV-1 MATRIX PROTEIN   MOLECULE: HIV-1 MATRIX PROTEIN;   CHAIN:
	VLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTIAVLY
	CVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAAD

**/
int			readFASTA(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readFASTA: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    long			i, sl, nres, n(0);
	long			m(0);
	string			s;
	Bstring			bs;
	Bmolecule*		mol = NULL;
	string			seq;
	
	nres = 0;
	while ( !fseq.eof() ) {
		getline(fseq, s);
		sl = s.length();
		if ( sl ) {
			if ( s.back() == '\n' ) {	// Remove trailing newline
				s.pop_back();
				sl = s.length();
			}
			if ( s[0] == '>' ) {
				if ( mol ) {
					seq = seq.substr(0,m);
					mol->seq = seq.c_str();
					mol->nres = nres;
 				}
				bs = s.substr(1);
				mol = molecule_add(&molgroup->mol, bs);
				cout << mol->id << endl;
				n++;
				nres = m = 0;
				seq.clear();
			} else {
				seq.append(s);
	            for ( i=0; i<sl; i++ ) {
    	            if ( isalpha(s[i]) ) {
        	            seq[m] = toupper(s[i]);
						nres++;
                	    m++;
	                }
					if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
						seq[m] = '-';
						m++;
					}
                }
            }
		}
    }
    
	if ( mol->nres < 1 ) {		// Clean up last sequence
		seq = seq.substr(0,m);
		mol->seq = seq.c_str();
		mol->nres = nres;
	}
    
    fseq.close();
	
    return n;
}


/**
@brief 	Writes FASTA format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules written (<0 if writing failed).
**/
int			writeFASTA(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeFASTA: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			i(0), j, n;
	Bmolecule*		mol;
	string			seq, s80;
    
    for ( mol = molgroup->mol; mol; mol = mol->next, ++i ) {
		if ( ( n = mol->seq.length() ) ) seq = mol->seq.str();
		else if ( ( n = mol->naseq.length() ) ) seq = mol->naseq.str();
		else n = 0;
		fseq << ">" << mol->id << endl;
        for ( j=0; j<n; j+=80 ) {
			s80 = seq.substr(j,80);
			fseq << s80 << endl;
        }
    }
    
    fseq.close();
    
    return i;
}

