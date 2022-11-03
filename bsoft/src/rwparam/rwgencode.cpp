/**
@file	rwgencode.cpp
@brief	Library routines to read and write genetic codes
@author Bernard Heymann
@date	Created: 20030316
@date	Modified: 20210328
**/

#include "rwgencode.h"
#include "star.h"
#include "mol_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
Bstring		read_gencode_star(Bstring& propfile);


/**
@brief 	Gets the genetic code from a parameter file.
@param 	&filename	file name (if empty, use a default file).
@return Bstring				64 character array of residues (+1 0 termination character).
**/
Bstring		get_genetic_code(Bstring& filename)
{
	if ( verbose & VERB_DEBUG )	
		cout << "DEBUG get_genetic_code: Initializing atomic parameters" << endl;
	
	Bstring			code;
		
	// Atom parameter file
	Bstring			gcfile = "gencode.star";
	if ( filename.c_str() ) gcfile = filename;
	
	Bstring			propfile = parameter_file_path(gcfile);
	Bstring			ext = gcfile.extension();
	if ( ext.length() ) {
		if ( ext.contains("star") )
			code = read_gencode_star(propfile);
	}
	
	if ( code.length() != 64 ) {
		if ( verbose )
			cout << "Genetic code file " << gcfile << " not opened! Using default code" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_genetic_code: " << code << endl;
		
	return code;
}

/**
@brief 	Writing genetic code.
@param 	&filename	file name.
@param 	&gc			64-byte string with the genetic code.
@return int 				0.
**/
int 		write_genetic_code(Bstring& filename, Bstring& gc)
{
	
	return 0;
}

Bstring		read_gencode_star(Bstring& propfile)
{
 	Bstar			star;
	
 	if ( star.read(propfile.str()) < 0 )
		error_show(propfile.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return NULL;
	}

	int				i, j, index;
	Bstring			gc(' ', 64L);

	for ( auto ib: star.blocks() ) {
		for ( auto il: ib.loops() ) {
			if ( ( i = il.find(RESPROP_CODE1) ) >= 0 &&
					( j = il.find(RESPROP_CODON) ) >= 0 ) {
				if ( il.data().size() != 64 ) {
					cerr <<  "Error: File " << propfile << " does contain only " << il.data().size() << " codons, 64 required!" << endl;
					return "";
				}
				for ( auto ir: il.data() ) {
					index = index_from_codon(ir[j].c_str());
					gc[index] = ir[i][0];
				}
			}
		}
	}
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_gencode_star: " << gc << endl;
		
	return gc;
}

/**
@brief 	Gets the index associated with a specific codon.
@param	*codon	the codon.
@return int					index.
**/
//int			index_from_codon(const char* codon)
int			index_from_codon(const char codon[3])
{
	int				i, j, index(0), e;
	char			acgt[8] = "ACGT";
	char			upcodon[4];
	
	for ( i=0; i<3; i++ ) {
		upcodon[i] = toupper(codon[i]);
		e = (int) pow(4.0, 2.0 - i);
		for ( j=0; j<4; j++ ) {
			if ( upcodon[i] == acgt[j] ) index += e*j;
		}
	}
	
	return index;
}
