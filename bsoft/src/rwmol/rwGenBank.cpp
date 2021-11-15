/**
@file	rwGenBank.cpp
@brief	Library routines to read and write GenBank sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20210423
**/

#include "rwGenBank.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads GenBank format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup molecule group.
@return int 				number of molecules read (<0 if reading failed).

	GenBank format:
	LOCUS       NP_039940                284 aa            linear   VRL 13-FEB-2003
	DEFINITION  UL6 [Human herpesvirus 5].
	ACCESSION   NP_039940
	ORIGIN      
	        1 mhakmngwag vrlvthclnt rsrtyvalnm lafartprgv psclfnkvwv sryalvlilm
	       61 vcasesstsw avtsnrlpnc stitttagqd aelhgpapls cnvtqwgrye ngstpvlwct
		  121 lwgsrtrvsl ghrvafgcsw ktffiynvse ssggty
	//

**/
int 	readGenBank(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readGenBank: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    string			s, seq;
    Bstring			seqname;
    size_t			i;
	int				type(0);		// Type of sequence: DNA=1, protein=2
	long			nid(0), nseq(0), nres(0), m(0);
	int				seq_flag(0);
	Bmolecule*		mol;
	
//	int 		gaps = 1;								// Retain gaps
//	if ( strstr(molgroup->select, "NOGAP") ) gaps = 0;	// Omit gaps
	
	// Pass to read the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			vector<string>	vs = split(s);
			if ( vs[0] == "LOCUS" ) {
				if ( s.find("bp") != string::npos || s.find("DNA") != string::npos ) type = 1;
				else type = 2;
				seqname = vs[1];
				seq_flag = 1;
				nid++;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readGenBank: type=" << type << endl;
			} else if ( vs[0] == "ORIGIN" ) {
				seq_flag++;
			} else if ( vs[0] == "//" ) {
				seq = seq.substr(0,m);
				mol = molecule_add(&molgroup->mol, seqname);
				if ( type == 1 ) {
					mol->naseq = seq.c_str();
					mol->nbase = nres;
				} else {
					mol->seq = seq.c_str();
					mol->nres = nres;
				}
				nres = 0;
				seq_flag = 0;
				seq.clear();
				m = 0;
				nseq++;
			} else if ( seq_flag ) {
				s = s.substr(10);
				seq.append(s);
	            for ( i=0; i<s.length(); i++ ) {
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
    
	fseq.close();
	
	if ( nseq != nid ) {
		cout << "Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")" << endl;
		return -1;
	}
	
    return nid;
}


/**
@brief 	Writes GenBank format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int 	writeGenBank(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGenBank: filename=" << filename << endl;
	
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
        fseq << "LOCUS       " << mol->id << endl;
//		fseq << mol->id << "," << mol->nres << " bases, 0 checksum" << endl;
        fseq << "ORIGIN" << endl;
        for ( j=0; j<n; j+=60 ) {
			fseq << setw(9) << right << j+1;
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seq.substr(k,10);
            fseq << endl;
        }
        fseq << "//" << endl << endl;
    }
    
	fseq.close();
    
    return i;
}

