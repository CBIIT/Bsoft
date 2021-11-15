/**
@file	rwEMBL.cpp
@brief	Library routines to read and write EMBL sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20210422
**/

#include "rwEMBL.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads EMBL format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules read (<0 if reading failed).

	EMBL format:
	ID   AQP1_HUMAN
	DE   AQP1_HUMAN, 527 bases, EB1B8D0 checksum.
	SQ             527 BP
	    ---------- ---------- ---------- ---------- ---------- ----MASEFK

**/
int 	readEMBL(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readEMBL: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    string			s, tag, seq;
	Bstring			bs;
    long			i;
	long			nid(0), nseq(0), nres(0);
	int				n(0), m(0), seq_flag(0);
	Bmolecule*		mol = NULL;

	// Pass to read the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		tag = s.substr(0,2);
		if ( s.length() > 5 ) s = s.substr(5);
		else s.clear();
		if ( tag == "ID" ) {
			bs = s;
			nid++;
		} else if ( tag == "//" ) {
			seq = seq.substr(0,m);
			mol = molecule_add(&molgroup->mol, bs);
 			mol->seq = seq.c_str();
			mol->nres = nres;
			nres = 0;
            seq_flag = 0;
			seq.clear();
 			m = 0;
            n++;
		} else if ( tag == "SQ" ) {
			seq_flag = 1;
 			nseq++;
        } else if ( seq_flag ) {
			seq.append(s);
			for ( i=0; i<s.length(); ++i ) {
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
    
	fseq.close();

	if ( nid != nseq ) {
		cerr << "Error: Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")!" << endl;
		return -1;
	}
	
    return nseq;
}

/**
@brief 	Writes EMBL format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules written (<0 if writing failed).
**/
int 	writeEMBL(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeEMBL: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			i(0), j, k, n;
	Bmolecule*		mol;
	string			seq;
    
    for ( mol = molgroup->mol; mol; mol = mol->next, ++i ) {
		if ( ( n = mol->seq.length() ) ) seq = mol->seq.str();
		else if ( ( n = mol->naseq.length() ) ) seq = mol->naseq.str();
		else n = 0;
        fseq << "ID   " << mol->id << endl;
        fseq << "DE   " << mol->id << ", " << n << " bases, 0 checksum" << endl;
        fseq << "SQ   " << n << " BP" << endl;
        for ( j=0; j<n; j+=60 ) {
			fseq << "    ";
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seq.substr(k,10);
            fseq << endl;
        }
        fseq << "//" << endl;
    }
    
	fseq.close();
    
    return i;
}

