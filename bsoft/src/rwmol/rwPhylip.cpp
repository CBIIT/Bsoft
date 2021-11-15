/**
@file	rwPhylip.cpp
@brief	Library routines to read and write Phylip sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20210423
**/

#include "rwPhylip.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Phylip format sequence files.
@param 	&filename		sequence file name.
@param	*molgroup 		molecule group.
@return int 				number of molecules read (<0 if reading failed).

	Phylip format:
	21   1419
	WMBEX6    MT------------------------------------------------
	P89429    MA------------------------------------------------
	P09302    MAEITSLFNNSS--------------------------------------

**/
int			readPhylip(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPhylip: filename=" << filename << endl;
	
    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
	Bstring			seqname;
	string			s;
    long			i(0), j, k(0), n(0), nseq(0), maxlen;
	Bmolecule*		mol;
	
//	int 		gaps = 1;								// Retain gaps
//	if ( strstr(molgroup->select, "NOGAP") ) gaps = 0;	// Omit gaps
	
	fseq >> nseq >> maxlen;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPhylip: nseq=" << nseq << " maxlen=" << maxlen << endl;
	
	// Set up molecule structures and read the first parts of the sequences
	vector<string>   	seq(nseq);

	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			if ( i < nseq ) {
				seqname = s;
				seqname = seqname.pre(' ');
				mol = molecule_add(&molgroup->mol, seqname);
			} else {
				for ( j=0, mol = molgroup->mol; mol && j<n; mol = mol->next ) j++;
			}
			if ( mol ) {
				seq[n].append(s.substr(10));
			}
			i++; n++;
		}
		if ( n >= nseq )	// Empty line and reset sequence counter
			n = 0;
	}
	
	for ( n=0, mol = molgroup->mol; n<nseq && mol; ++n, mol = mol->next ) {
		s = seq[n];
		for ( i=0, k=0; i<s.length(); ++i ) {
			if ( isalpha(s[i]) ) s[k++] = s[i];
			if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) s[k++] = '-';
		}
		s = s.substr(0,k);
		mol->seq = s.c_str();
	}
		
	fseq.close();
	
    return nseq;
}

/**
@brief 	Writes Phylip format sequence files.
@param 	&filename		sequence file name.
@param	*molgroup 		molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int			writePhylip(Bstring& filename, Bmolgroup* molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePhylip: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			m, n;
	int				nseq;
	Bmolecule*		mol;
	string			molname, seq, s;
    
	for ( nseq = 0, mol = molgroup->mol; mol; mol = mol->next ) nseq++;
	
	fseq << setw(6) << nseq << setw(7) << molgroup->maxlen << endl;

	for ( n = 0; n < molgroup->maxlen; n += 50 ) {
		for ( mol = molgroup->mol; mol; mol = mol->next ) {
			molname = mol->id.substr(0,10).str();
			if ( mol->seq.length() ) seq = mol->seq.str();
			else if ( mol->naseq.length() ) seq = mol->naseq.str();
			m = 50;
			if ( molgroup->maxlen - n < 50 ) m = molgroup->maxlen - n;
			s.clear();
			if ( n < seq.length() ) s = seq.substr(n, 50);
			if ( s.length() < m ) s.append(m-s.length(), '-');
			if ( n == 0 ) fseq << setw(10) << molname << " " << s << endl;
			else fseq << setw(11) << " " << s << endl;
		}
		fseq << endl;
	}
    
	fseq.close();
    
    return nseq;
}

