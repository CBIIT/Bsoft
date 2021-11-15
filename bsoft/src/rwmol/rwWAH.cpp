/**
@file	rwWAH.cpp
@brief	Library routines to read and write Wayne Hendrickson coordinate files
@author Bernard Heymann
@date	Created: 20050217
@date	Modified: 20120212
**/

#include "rwWAH.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a Wayne Hendrickson coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules read (<0 if reading failed).

	WAH format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	   SER 1071CA     8.38771  32.21584 115.00745   0.00000   0.00000
	   ARG 1072C     11.26685  34.89560 115.61476   0.00000   0.00000

**/
int 	readWAH(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readWAH: filename=" << filename << endl;
	
    ifstream		fwah(filename.c_str());
    if ( fwah.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    int				atomnum=0, resnum, prev_resnum = -1, n = 0, nmol = 0;
    char			aline[MAXLINELEN] = "A";
    char			restype[20], atomtype[20];
	
	Bmolecule*		mol = molecule_add(&molgroup->mol, aline);
	Bresidue*		res = NULL;
	Batom*  		atom = NULL;
	
    if ( verbose )
		cout << "Title:                          " << mol->id << endl;
		
	while ( !fwah.eof() ) {
    	fwah.getline(aline, MAXLINELEN);
		strncpy(restype, aline + 3, 4);
		resnum = get_integer(aline + 7, 4);
		strncpy(atomtype, aline + 11, 4);
		restype[4] = 0; 	// Residue names only 4 characters
		atomtype[4] = 0; 	// Atom names only 4 characters
		if ( n < 1 || resnum != prev_resnum ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readWAH: resnum=" << resnum << " restype=" << restype << endl;
			if ( mol->res ) res = residue_add(&res, restype);
			else res = residue_add(&mol->res, restype);
			res->num = resnum;
			prev_resnum = resnum;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readWAH: atomnum=" << atomnum << " atomtype=" << atomtype << endl;
		if ( !res->atom ) atom = atom_add(&res->atom, atomtype);
		else atom = atom_add(&atom, atomtype);
		atom->num = ++atomnum;
		atom->coord[0] = get_float(aline + 15, 10);
		atom->coord[1] = get_float(aline + 25, 10);
		atom->coord[2] = get_float(aline + 35, 10);
		atom->q = get_float(aline + 45, 10);
		atom->b = get_float(aline + 55, 10);
		atom->sel = 1;	    /* All atoms read selected */
		n++;
    }
    
 	if ( n ) nmol = 1;
 	else nmol = -1;

	fwah.close();
    
    return nmol;
}

/**
@brief 	Writes a Wayne Hendrickson coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules written (<0 if writing failed).

	WAH format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	   SER 1071CA     8.38771  32.21584 115.00745   0.00000   0.00000
	   ARG 1072C     11.26685  34.89560 115.61476   0.00000   0.00000

**/
int 	writeWAH(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeWAH: filename=" << filename << endl;
	
	if ( !molgroup ) return -1;
	
    ofstream		fwah(filename.c_str());
    if ( fwah.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			nw = 0;
	int				nmol = 0;
	char			string[MAXLINELEN];
    
	if ( verbose & VERB_DEBUG )
		cout << "New file name:                  " << filename << endl;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		nw = 0;
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( atom->sel ) {
					snprintf(string, MAXLINELEN, "   %-4s%4d%-4s%10.5f%10.5f%10.5f%10.5f%10.5f\n",
						res->type, res->num, atom->type,
						atom->coord[0], atom->coord[1], atom->coord[2], atom->q, atom->b);
					fwah << string;
					nw++;
				}
			}
		}
		if ( nw > 0 ) nmol++;
    }
    
    fwah.close();
    
    return nmol;
}

