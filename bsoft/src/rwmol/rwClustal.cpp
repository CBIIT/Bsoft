/**
@file	rwClustal.cpp
@brief	Library routines to read and write Clustal sequence files
@author Bernard Heymann
@date	Created: 20030309
@date	Modified: 20210422
**/

#include "rwClustal.h"
#include "utilities.h"
#include <fstream>
#include <sstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Clustal format sequence files.
@param 	&filename	sequence file name.
@param	*molgroup 	molecule group.
@return int 			number of molecules read (<0 if reading failed).

	Clustal format:
	
	NP_039940       -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY
	P16720          -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY
	DAA00163        -----MHAKMNGWAGVRL--VTHCLNTRSRTYVALNMLAFARTPRGVPSCLFNKVWVSRY

**/
int 	readClustal(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readClustal: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

	string			s, seqname, aseq;
	Bstring			bs;
    long     		m(0), nseq(0);
	Bmolecule*		mol = NULL;
	
	// Read the header line
	getline(fseq, s);
	
	// Read the rest of the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			istringstream	ss(s);
			ss >> seqname >> aseq;
			if ( isalnum(seqname[0]) ) {
				bs = seqname;
				if ( m == 0 ) {
					mol = molecule_add(&molgroup->mol, bs);
					nseq++;
				}
				mol->seq = mol->seq + aseq;
				mol = mol->next;
			}
			seqname.clear();
		} else {
			if ( !mol && nseq > 0 ) m++;		// Test for end of first section
			mol = molgroup->mol;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readClustal: sequences read = " << nseq << " lines per sequence = " << m << endl;
		
	fseq.close();
	
    return nseq;
}

/**
@brief 	Writes Clustal format sequence files.
@param 	&filename		sequence file name.
@param	*molgroup 		molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int 	writeClustal(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeClustal: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			m, n, nseq;
    long			max_name_len(10), max_seq_len(60);
	Bmolecule*		mol;
	Bstring			molname;
	string			s, seq;
    
	for ( nseq = 0, mol = molgroup->mol; mol; mol = mol->next ) {
		nseq++;
		molname = mol->id.pre(' ');
		if ( max_name_len < molname.length() ) max_name_len = molname.length();
	}
	if ( max_name_len > 10 ) max_seq_len = 50;
	max_name_len += 6;
	
	time_t		ti = time(NULL);

	fseq << "CLUSTAL alignment file, Written by Bsoft on " << asctime(localtime(&ti)) << endl << endl;

	for ( n = 0; n < molgroup->maxlen; n += max_seq_len ) {
		fseq << endl;
		for ( mol = molgroup->mol; mol; mol = mol->next ) {
			if ( mol->seq.length() ) seq = mol->seq.str();
			else if ( mol->naseq.length() ) seq = mol->naseq.str();
			m = max_seq_len;
			if ( molgroup->maxlen - n < max_seq_len ) m = molgroup->maxlen - n;
			s.clear();
			if ( n < seq.length() ) s = seq.substr(n, max_seq_len);
			if ( s.length() < m ) s.append(m-s.length(), '-');
			molname = mol->id.pre(' ');
			fseq << setw(max_name_len) << molname << s << endl;
		}
		fseq << endl;
	}
    
	fseq.close();
    
    return nseq;
}

