/**
@file	seq_util.cpp
@brief	Sequence utility functions 
@author Bernard Heymann
@date	Created: 20001029
@date	Modified: 20220713
**/
 
#include "rwgencode.h"
#include "seq_analysis.h"
#include "seq_util.h" 
#include "linked_list.h"
#include "utilities.h" 

#include <map>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
  
/* One-letter to 3-letter mapping */
/*const char*   res_code[] = { 
        "-", "GAP", 
        "A", "ALA", 
        "B", "ASX", 
        "C", "CYS", 
        "D", "ASP", 
        "E", "GLU", 
        "F", "PHE", 
        "G", "GLY", 
        "H", "HIS", 
        "I", "ILE", 
        "K", "LYS", 
        "L", "LEU", 
        "M", "MET", 
        "N", "ASN", 
        "P", "PRO", 
        "Q", "GLN", 
        "R", "ARG", 
        "S", "SER", 
        "T", "THR", 
        "V", "VAL", 
        "W", "TRP", 
        "Y", "TYR", 
        "Z", "GLX", 
        "X", "UNK", 
} ; 
*/
const map<char, string> res_code = {
        {'-', "GAP"}, 
        {'*', "UNK"}, 
        {'A', "ALA"}, 
        {'B', "ASX"}, 
        {'C', "CYS"}, 
        {'D', "ASP"}, 
        {'E', "GLU"}, 
        {'F', "PHE"}, 
        {'G', "GLY"}, 
        {'H', "HIS"}, 
        {'I', "ILE"}, 
        {'K', "LYS"}, 
        {'L', "LEU"}, 
        {'M', "MET"}, 
        {'N', "ASN"}, 
        {'P', "PRO"}, 
        {'Q', "GLN"}, 
        {'R', "ARG"}, 
        {'S', "SER"}, 
        {'T', "THR"}, 
        {'V', "VAL"}, 
        {'W', "TRP"}, 
        {'Y', "TYR"}, 
        {'Z', "GLX"}, 
        {'X', "UNK"}
} ;

/**
@brief 	Obtain sequences from residues.
@param 	*molgroup		set of molecules.
@return int				0
**/
int			seq_from_residues(Bmolgroup* molgroup)
{
	long			maxnum;
	Bstring			seq;
	Bmolecule*		mol;
	Bresidue*		res;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		if ( mol->res ) {
			for ( maxnum=mol->nres, res = mol->res; res; res = res->next )
				if ( maxnum < res->num ) maxnum = res->num;
			mol->seq = Bstring('-', maxnum);
			mol->nres = 0;
			for ( res = mol->res; res; res = res->next, mol->nres++ )
				mol->seq[res->num-1] = getcode1(res->type);
		}
	}
	
	return 0;
}

/**
@brief 	Shows all molecular sequences.
@param 	*molgroup		set of sequences.
@return int				0
**/
int			seq_show(Bmolgroup* molgroup)
{
	long			i, maxnum;
	Bstring			seq;
	Bmolecule*		mol;
	Bresidue*		res;
	Bsecondary*		sec;
	
	if ( verbose & VERB_LABEL )
		cout << endl << "Sequences:" << endl << endl;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		cout << "Molecule " << mol->id << ": " << mol->nbase << " nucleotides, "
			<< mol->nres << " residues" << endl;
		if ( mol->naseq.length() ) cout << mol->naseq << endl;
		if ( mol->res ) {
			for ( maxnum=mol->nres, res = mol->res; res; res = res->next )
				if ( maxnum < res->num ) maxnum = res->num;
			seq = Bstring('-', maxnum);
			for ( res = mol->res; res; res = res->next )
				seq[res->num-1] = getcode1(res->type);
			cout << seq << endl;
			seq = 0;
		} else if ( mol->seq.length() ) cout << mol->seq << endl;
		if ( mol->sec ) {
			for ( maxnum=0, sec=mol->sec; sec; sec=sec->next )
				if ( maxnum < sec->last->num ) maxnum = sec->last->num;
			seq = Bstring('-', maxnum);
			for ( sec=mol->sec; sec; sec=sec->next ) {
				if ( sec->type < Strand )
					for ( i=sec->first->num-1; i<sec->last->num; i++ ) seq[i] = 'H';
				if ( sec->type == Strand )
					for ( i=sec->first->num-1; i<sec->last->num; i++ ) seq[i] = 'E';
			}
			cout << seq << endl;
			seq = 0;
		}
		cout << endl;
	}
	
	return 0;
}
 
/**
@brief 	Shows the masses of all molecular sequences.
@param 	*molgroup		set of sequences.
@return int				0
**/
int			seq_mass(Bmolgroup* molgroup)
{
	long			i, n;
	double			mass;
	Bstring			temp;
	Bmolecule*		mol;
	Bresidue*		res;
	Bresidue_type*	rt_list = get_residue_properties(temp);
	Bresidue_type*	rt = NULL;
		
	if ( verbose & VERB_LABEL )
		cout << "Calculating molecular weights:" << endl << endl;
	
	cout << "Molecule\tMass" << endl;
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mass = 0;
		n = mol->seq.length();
		if ( n ) {
			for ( i=0; i<n; i++ ) {
				for ( rt = rt_list; rt && rt->c != mol->seq[i]; rt = rt->next ) ;
				if ( rt ) mass += rt->mass;
			}
		} else if ( mol->res ) {
			for ( res = mol->res; res; res = res->next ) {
				for ( rt = rt_list; rt && rt->cod != res->type; rt = rt->next ) ;
				if ( rt ) mass += rt->mass;
			}
		}
		mass -= 18*(n-1);
		cout << mol->id << tab << mass << endl;
	}

	cout << endl;

	kill_list((char *) rt_list, sizeof(Bresidue_type));
	
	return 0;
}

/**
@brief 	Shows the elemental composition of all molecular sequences.
@param 	*molgroup		set of sequences.
@param 	&paramfile		file of residue parameters.
@return vector<double>	array of element numbers: HCNOS
**/
vector<double>	seq_elements(Bmolgroup* molgroup, Bstring& paramfile)
{
	long			i, j, n, nat, natr;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Bresidue_type*	rt_list = get_residue_properties(paramfile);
	Bresidue_type*	rt = NULL;
	vector<double>	el(5,0), elr(5,0), elall(5,0), el1(5,0);
		
	if ( verbose & VERB_LABEL )
		cout << "Calculating the elemental composition:" << endl << endl;
	
	cout << "Molecule\tH\tC\tN\tO\tS\tTotal" << endl;
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( j=0; j<5; ++j ) el[j] = elr[j] = 0;
		n = mol->seq.length();
		if ( n ) {
			for ( i=0; i<n; i++ ) {
				for ( rt = rt_list; rt && rt->c != mol->seq[i]; rt = rt->next ) ;
				if ( rt ) for ( j=0; j<5; ++j ) el[j] += rt->comp[j];
			}
		}
//		} else if ( mol->res ) {
/*		if ( mol->res ) {
			for ( res = mol->res; res; res = res->next ) {
				char		c1 = getcode1(res->type);
				for ( rt = rt_list; rt && rt->c != c1; rt = rt->next ) ;
				if ( rt ) for ( j=0; j<5; ++j ) elr[j] += rt->comp[j];
			}
		}*/
		if ( mol->res ) {
			for ( res = mol->res; res; res = res->next ) {
				char		c1 = getcode1(res->type);
				for ( rt = rt_list; rt && rt->c != c1; rt = rt->next ) ;
				for ( j=0; j<5; ++j ) el1[j] = 0;
				for ( atom = res->atom; atom; atom = atom->next ) if ( atom->el[1] == 0 ) {
					if ( atom->el[0] == 'H' ) elr[0]++;
					if ( atom->el[0] == 'C' ) elr[1]++;
					if ( atom->el[0] == 'N' ) elr[2]++;
					if ( atom->el[0] == 'O' ) elr[3]++;
					if ( atom->el[0] == 'S' ) elr[4]++;
					if ( atom->el[0] == 'H' ) el1[0]++;
					if ( atom->el[0] == 'C' ) el1[1]++;
					if ( atom->el[0] == 'N' ) el1[2]++;
					if ( atom->el[0] == 'O' ) el1[3]++;
					if ( atom->el[0] == 'S' ) el1[4]++;
				}
/*				cout << c1;
				for ( j=0; j<5; ++j ) cout << tab << el1[j];
				cout << endl;
				cout << c1;
				for ( j=0; j<5; ++j ) cout << tab << rt->comp[j];
				cout << endl;*/
			}
		}
		cout << mol->id;
		for ( j=nat=natr=0; j<5; ++j ) {
			nat += el[j];
			natr += elr[j];
			elall[j] += el[j];
			cout << tab << el[j];
		}
		cout << tab << nat << endl;
//		cout << mol->id;
//		for ( j=0; j<5; ++j ) cout << tab << elr[j];
//		cout << tab << natr << endl;
	}

	cout << "Total";
	for ( j=n=0; j<5; n+=elall[j], ++j ) cout << tab << elall[j];
	cout << tab << n << endl << "Percentage:";
	for ( j=0; j<5; ++j ) cout << tab << elall[j]*100.0/n;
	cout << endl << endl;

	kill_list((char *) rt_list, sizeof(Bresidue_type));
	
	return elall;
}
 

/**
@brief 	Complements all nucleotide sequences.
@param 	*molgroup 		the molecule group.
@return int 			0.

	Search through a list of 1-3 mappings for the desired letter.

**/
int 		seq_complement_all(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	
	if ( verbose )
		cout << endl << "Complementing all nucleotide sequences" << endl << endl;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) 	// Complement all molecules
		if ( mol->naseq.length() ) complement_sequence(mol->naseq);
	
	return 0;
}

/**
@brief 	Translates all nucleotide sequences to protein sequences.
@param 	*molgroup the molecule group.
@param 	frame			the frame for translation.
@param 	&gcfile			file with genetic code.
@return int 			0.

	Each nucleic acid sequence in the molecule group is translated to the
	protein sequence.

**/
int 		seq_translate_all(Bmolgroup* molgroup, int frame, Bstring& gcfile)
{
	Bmolecule*		mol;
	
	Bstring			gencode = get_genetic_code(gcfile);
	
	if ( verbose & VERB_LABEL )
		cout << endl << "Translating all nucleotide sequences to amino acid sequences in frame " << frame << endl << endl;
	
	molgroup->maxlen = 0;
	for ( mol = molgroup->mol; mol; mol = mol->next ) if ( mol->naseq.length() ) {
		mol->seq = sequence_translate(mol->naseq, frame, gencode);
		mol->nres = mol->seq.length();
		if ( molgroup->maxlen < mol->nres ) molgroup->maxlen = mol->nres;
	}
	
	return 0;
}

/**
@brief 	Finds a nucleotide sequence.
@param 	*molgroup the molecule group.
@param 	&seq			sequence to find.
@return long 			position.
**/
long 		seq_find_dna(Bmolgroup* molgroup, Bstring& seq)
{
	Bmolecule*		mol;
	Bstring			sseq = seq.upper();
	Bstring			mol_naseq;
	long			i, offset(0);
	
	if ( verbose )
		cout << "Searching for: " << sseq << endl << endl;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mol_naseq = mol->naseq.upper();
		if ( mol_naseq.length() ) {
			offset = mol_naseq.find(sseq, 0);
			if ( offset > -1 ) {
				cout << mol->id << ": " << offset << endl;
				for ( i=0; i<offset; i++ ) cout << " ";
				cout << sseq << endl << mol_naseq << endl;
			} else {
				complement_sequence(mol_naseq);
				offset = mol_naseq.find(sseq, 0);
				if ( offset > -1 ) {
					cout << mol->id << ": " << offset << " (complement)" << endl;
					for ( i=0; i<offset; i++ ) cout << " ";
					cout << sseq << endl << mol_naseq << endl;
				} else {
					cout << mol->id << ": not found" << endl;
				}
			}
		} else
			cout << "Molecule " << mol->id << ": Nucleic acid sequence not found" << endl;
	}
	
	return offset;
}

/**
@brief 	Finds an amino acid sequence.
@param 	*molgroup 		the molecule group.
@param 	&seq			sequence to find.
@return long 			position.
**/
long 		seq_find_protein(Bmolgroup* molgroup, Bstring& seq)
{
	Bmolecule*		mol;
	Bstring			sseq = seq.upper();
	Bstring			mol_seq;
	long			offset(0);
	
	if ( verbose )
		cout << "Searching for: " << sseq << endl << endl;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mol_seq = mol->seq.upper();
		if ( mol_seq.length() ) {
			offset = mol_seq.find(sseq, 0);
			if ( offset > -1 ) cout << mol->id << ": " << offset << endl;
			else cout << mol->id << ": not found" << endl;
		} else
			cout << "Molecule " << mol->id << ": Protein sequence not found" << endl;
	}
	
	return offset;
}

/**
@brief 	Finds the coding region for an amino acid sequence.
@param 	*molgroup 		the molecule group.
@param 	&seq			sequence to find.
@param 	seqlenmin		sequence length minimum.
@param 	seqlenmax		sequence length maximum.
@param 	side1			preceding sequence length to include.
@param 	side2			succeeding sequence length to include.
@param 	threshold		threshold for reporting possible hits.
@param 	&gcfile			file with genetic code.
@return Bstring 		coding sequence.

	All molecules in the group are searched in all 6 possible frames.

**/
Bstring		seq_find_protein_in_dna(Bmolgroup* molgroup, Bstring& seq, int seqlenmin, int seqlenmax, 
				int side1, int side2, double threshold, Bstring& gcfile)
{
	long			n = seq.length();
	if ( n < 1 ) return "";
	
	Bmolecule*		mol;
	Bmolecule*		bestmol = NULL;
	if ( seqlenmin < n/2 ) seqlenmin = n/2;
	if ( seqlenmax <= n ) seqlenmin = 2*n;
	
	int				linelength = 80;
	int 			f, i, j, nres = 0, start, length, frame, term, seqlen = 0;
	int 			bestscore = 0, bestframe = 0;
	int 			bestindex = 0, bestseqlen = 0;
	Bstring			bestseq;
	char			format[100];
	snprintf(format, 100, "%%.%ds%c", linelength, '\n');
	
	Bstring			gencode = get_genetic_code(gcfile);
	
	if ( verbose ) {
		cout << endl << "Finding an amino acid sequence in translated nucleotide sequences:" << endl;
		cout << "Sequence length range:          " << seqlenmin << " - " << seqlenmax << endl;
		cout << "Score threshold (residues):     " << threshold << endl;
	}
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) { 	// Search all molecules
		for ( f=0; f<6; f++ ) { 			// Search all translational frames
			Bstring		nucseq = mol->naseq;
			frame = f;
			if ( f > 2 ) {
				frame = f - 3;
				complement_sequence(nucseq);
			}
			Bstring		aaseq = sequence_translate(nucseq, frame, gencode);
			nres = aaseq.length();
			
			vector<int>	score(nres,0);
			
			for ( i=0; i<=nres - n; i++ ) {	// Loop over starting positions
				if ( i == 0 || aaseq[i] == '*' ) {
					seqlen = 0;
					while ( ( i+seqlen+1 < nres ) && ( aaseq[i+seqlen+1] != '*' ) )
						seqlen++;
				}
				if ( seqlen >= seqlenmin && seqlen <= seqlenmax ) {
					term = 0;
					for ( j=0; j<n && j<nres && !term; j++ ) {
						if ( aaseq[i+j] == seq[j] ) score[i]++;
						if ( aaseq[i+j] == '*' ) term = 1;
					}
					if ( !term && threshold > 0 && score[i] >= threshold ) {
						if ( verbose )
							cout << "Frame: " << f << "  Position: " << 3*i + frame << 
								"   Length: " << seqlen << "   Score: " << score[i] << 
								" (" << score[i]*100.0/seq.length() << " %)" << endl;
						for ( j=0; j<n; j+= linelength )
							printf(format, &aaseq[i + j]);
					}
					if ( !term && bestscore < score[i] ) {
						bestscore = score[i];
						bestseqlen = seqlen;
						bestmol = mol;
						bestframe = f;
						bestindex = 3*i + frame;
						start = (i - side1)*3 + frame;
						while ( start < 0 ) start += 3;
						length = (n + side1 + side2)*3;
						if ( length + start > nucseq.length() )
							length = nucseq.length() - start;
						bestseq = nucseq.substr(start, length);
					}
				}
			}
		}
	}
	
	cout << "Best sequence =" << bestseq << endl;
	
	if ( bestmol && verbose ) {
		cout << endl << "Best sequence score:            " << bestscore << 
			" (" << bestscore*100.0/seq.length() << "%)" << endl;
		cout << "Molecule:                       " << bestmol->id << endl;
		cout << "Length:                         " << bestseqlen << endl;
		if ( bestframe < 3 ) cout << "Frame:                          " << bestframe << endl;
		else cout << "Frame:                          " << bestframe - 3 << " (complement)" << endl;
		cout << "Index:                          " << bestindex << endl;
		cout << "Nucleotide sequence:" << endl;
		for ( i=0; i<bestseq.length(); i+= linelength )
			printf(format, bestseq.c_str() + i);
		Bstring		aaseq = sequence_translate(bestseq, 0, gencode);
		cout << "Amino acid sequence:" << endl;
		for ( i=0; i<aaseq.length(); i+= linelength )
			printf(format, aaseq.c_str() + i);
		cout << "Requested sequence:" << endl;
		for ( i=0; i<seq.length(); i+= linelength )
			printf(format, seq.c_str() + i);
	}
	cout << endl;
	
	return bestseq;
}

/**
@brief 	Converts a one-letter amino acid designation to the three-letter equivalent.
@param 	c			the desired amino acid code letter
@param 	*cod		the corresponding three-letter code
@return int			0.

	Search through a list of 1-3 mappings for the desired letter.

**/
int			getcode3(char c, char* cod) 
{ 
	c = toupper(c); 
	
	string		rc(res_code.at(c));
	for ( int i=0; i<3; ++i ) cod[i] = rc[i];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG getcode3: c=" << c << " cod=" << cod << endl;
		 
	return 0; 
} 
 
/**
@brief 	Converts a three-letter amino acid designation to the one-letter equivalent.
@param 	*acode 		the desired amino acid three-letter code
@return char		the corresponding one-letter code

	Search through a list of 3-1 mappings for the desired three-letter code.

**/
char		getcode1(char* acode) 
{ 
	size_t			i;
	char			c(' ');
	 
	for ( i=0; i<3 && i<strlen(acode); i++ ) acode[i] = toupper(acode[i]);
	acode[i] = 0;
	
	for ( auto it=res_code.begin(); it!=res_code.end(); ++it )
		if ( it->second == acode ) c = it->first;
	 
	return c;
} 
 
/**
@brief 	Complements a nucleotide sequence in place.
@param 	&nucseq		nucleotide sequence to be translated.
@return int 		0.
**/
int 		complement_sequence(Bstring& nucseq)
{
	int 		i, j;
	char		nuc;
	
	j = nucseq.length();
	for ( i=0; i<(j+1)/2; i++ ) {
		nuc = get_complement(nucseq[j-i-1]);
		nucseq[j-i-1] = get_complement(nucseq[i]);
		nucseq[i] = nuc;
	}
	
	return 0;
}

/**
@brief 	Get the Watson-Crick complement of a nucleotide base.
@param 	nuc			nucleotide.
@return char 		complementing nucleotide.
**/
char		get_complement(char nuc)
{
	char		nunuc = 'X';
	
	switch ( nuc ) {
		case 'A': nunuc = 'T'; break;
		case 'G': nunuc = 'C'; break;
		case 'C': nunuc = 'G'; break;
		case 'T': nunuc = 'A'; break;
		case 'U': nunuc = 'A'; break;
		default: break;
	}

	return nunuc;
}

/**
@brief 	Translates a nucleotide sequence to a protein sequence.
@param 	&nucseq		nucleotide sequence to be translated.
@param 	frame		coding frame.
@param 	&gencode	genetic code: array of amino acids.
@return Bstring 	translated protein sequence.
**/
Bstring 	sequence_translate(Bstring& nucseq, long frame, Bstring& gencode)
{
	long 			i, j;
	long 			nnuc = nucseq.length();
	long 			nres = (nnuc - frame)/3;
	Bstring			aaseq(nres, ' ');
	
	for ( i=0, j=frame; i<nres && j<nnuc; i++, j+=3 )
		aaseq[i] = gencode[index_from_codon(&nucseq[j])];
	
	if ( i < nres ) aaseq[i] = 0;
	
	return aaseq;
}


