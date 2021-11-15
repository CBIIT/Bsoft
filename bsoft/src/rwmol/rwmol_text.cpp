/**
@file	rwmol_text.cpp
@brief	Library routines to read and write molecule files in plain text
@author Bernard Heymann
@date	Created: 20050419
@date	Modified: 20210423
**/

#include "rwmolecule.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a text file.
@param 	&filename		the file name.
@param 	*molgroup		the molecule group.
@return int 				number of molecules read (<0 if reading failed).
**/
int 		read_mol_text(Bstring& filename, Bmolgroup *molgroup)
{
  	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_mol_text: filename=" << filename << endl;
	
	ifstream    	fseq(filename.c_str());
    if ( fseq.fail() ) {
		cout << "File not opened!" << endl;
		return -1;
	}
    
	long			i, nres(0), nna(0), m(0), nseq(0);
	Bstring			molname;
	char			na[8] = "ACGTU";
    string			s, seq;
	Bmolecule*		mol = NULL;
	
 	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			if ( isalnum(s[0]) ) {
				molname = s;
				mol = molecule_add(&molgroup->mol, molname);
				seq.clear();
				m = nna = 0;
				nseq++;
			} else {
				s = s.substr(10);
				seq.append(s);
	            for ( i=0; i<s.length(); i++ ) {
    	            if ( isalpha(s[i]) ) {
        	            seq[m] = toupper(s[i]);
						if ( strchr(na, seq[m]) ) nna++;
						nres++;
                	    m++;
	                }
					if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
						seq[m] = '-';
						m++;
					}
                }
			}
		} else if ( m && mol ) {
			seq = seq.substr(0,m);
			if ( 1.5*nna > nres ) {
				mol->naseq = seq;
				mol->nbase = nres;
			} else {
				mol->seq = seq;
				mol->nres = nres;
			}
		}
	}

    fseq.close();	

	return nseq;
}

/**
@brief 	Writes a molecule group to a text file.
@param 	&filename		the file name.
@param 	*molgroup		the molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int 		write_mol_text(Bstring& filename, Bmolgroup *molgroup)
{
  	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_mol_text: filename=" << filename << endl;
	
    ofstream    	fseq(filename.c_str());

    if ( fseq.fail() ) return -1;

    int     	i(0), j, k, n;
	string		seq;
	Bmolecule*	mol;

    for ( mol = molgroup->mol; mol; mol = mol->next, ++i ) {
		fseq << mol->id << endl;
		n = mol->naseq.length();
        seq = mol->naseq.str();
        for ( j=0; j<n; j+=60 ) {
			fseq << setw(9) << right << j+1;
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seq.substr(k,10);
            fseq << endl;
        }
		n = mol->seq.length();
		seq = mol->seq.str();
        for ( j=0; j<n; j+=60 ) {
			fseq << setw(9) << right << j+1;
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seq.substr(k,10);
            fseq << endl;
        }
        fseq << endl;
    }

    fseq.close();
    
	return i;
}
