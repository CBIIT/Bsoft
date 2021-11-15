/**
@file	rwPIR.cpp
@brief	Library routines to read and write PIR sequence files
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20210423
**/

#include "rwPIR.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads PIR format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules read (<0 if reading failed).

	PIR format:
	>P1;BAC1_HALS1
	BAC1_HALS1, 308 bases, DC74A5E6 checksum.
	 ---------- ---MDPIALT AAVGADLLGD GRPETLWLGI GTLLMLIGTF
	 ASAAD---*

**/
int			readPIR(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPIR: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    string			s, seq;
    Bstring			seqname;
	long			i, nid(0), nseq(0), nres, m(0);
	int				seq_flag(0);
	Bmolecule*		mol;

//	int 		gaps = 1;								// Retain gaps
//	if ( strstr(molgroup->select, "NOGAP") ) gaps = 0;	// Omit gaps
	
	// Pass to read the data
	nres = nid = nseq = 0;
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.substr(0,3) == ">P1" ) {
			seqname = s.substr(4);
            seq_flag = 1;
			nid++;
		} else if ( seq_flag == 1 ) {
			seq_flag++;
		} else if ( seq_flag > 1 ) {
			seq.append(s);
            for ( i=0; i<s.length(); i++ ) {
                if ( isalpha(s[i]) ) {
                    seq[m] = toupper(s[i]);
					nres++;
                    m++;
                }
                if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
//					if ( gaps ) {
	                    seq[m] = '-';
    	                m++;
//					}
                }
				if ( s[i] == '*' ) {
					seq = seq.substr(0,m);
					mol = molecule_add(&molgroup->mol, seqname);
					mol->seq = seq.c_str();
					mol->nres = nres;
					nres = 0;
            		seq_flag = 0;
					seq.clear();
					m = 0;
            		nseq++;
				}
            }
		}
    }
    
	fseq.close();
	
	if ( nseq != nid ) {
		cout << "Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")" << endl;
		return -1;
	}
	
    return nid;
}


/**
@brief 	Writes PIR format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules written (<0 if writing failed).
**/
int			writePIR(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePIR: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    int				i(0);
    long			j, k, n;
	Bmolecule*		mol;
	string			seq;
    
    for ( mol = molgroup->mol; mol; mol = mol->next, ++i ) {
		if ( ( n = mol->seq.length() ) ) seq = mol->seq.str();
		else if ( ( n = mol->naseq.length() ) ) seq = mol->naseq.str();
		else n = 0;
		fseq << ">P1;" << mol->id << endl;
		fseq << mol->id << ", " <<  mol->nres << " bases, 0 checksum" << endl;
        for ( j=0; j<n; j+=50 ) {
			for ( k=j; k<j+50 && k<n; k+=10 )
				fseq << " " << seq.substr(k,10);
            fseq << endl;
        }
        fseq << "*" << endl;
    }
    
    fseq.close();
    
    return i;
}
