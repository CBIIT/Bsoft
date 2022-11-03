/**
@file	rwmol_star.cpp
@brief	Library routines to read and write molecule files in STAR format
@author Bernard Heymann
@date	Created: 19980822 
@date	Modified: 20210329
**/

#include "star.h"
#include "rwmolecule.h"
#include "rwmol_star.h"
#include "mol_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a STAR format file.
@param 	&filename		the file name.
@param 	*molgroup		the molecule group.
@return long				number of molecules read (<0 if reading failed).
**/
long		read_mol_star(Bstring& filename, Bmolgroup *molgroup)
{
 	Bstar			star;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_mol_star: filename=" << filename << endl;

 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return 0;
	}

	long			i, nmol(0), nres(0), natom(0), resnum(0);
	Bstring			molname, resname;
	Bmolecule*		mol = NULL;
	Bresidue*		res = NULL;
	Batom*			atom = NULL;

	for ( auto ib: star.blocks() ) {
		molname = ib.tag();
		if ( molname.length() < 1 ) molname = "A";
		mol = molecule_add(&molgroup->mol, molname);
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_mol_star: Molecule " << mol->id << endl;
		mol->seq = ib.at(MOLECULE_SEQUENCE).c_str();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_mol_star: sequence " << mol->seq << endl;
		mol->nres = ib.real(MOLECULE_LENGTH);
		if ( mol->seq.length() > 0 )
			if ( mol->nres < 1 ) mol->nres = mol->seq.length();
		for ( auto il: ib.loops() ) {
			if ( il.find(ATOM_RESNUMBER) >= 0 ) {
//				il.show_tags();
				mol->nres = 0;
				for ( auto ir: il.data() ) {
//					for ( auto it: ir ) cout << it << endl;
					if ( ( i = il.find(ATOM_RESNUMBER) ) >= 0 ) {
						resnum = to_integer(ir[i]);
						if ( !resnum ) break;
						if ( !res || res->num != resnum ) {
							res = (Bresidue *) add_item((char **) &res, sizeof(Bresidue));
							if ( !mol->res ) mol->res = res;
							res->num = resnum;
							atom = NULL;
							nres++;
							mol->nres++;
//							cout << resnum << endl;
						}
						atom = (Batom *) add_item((char **) &atom, sizeof(Batom));
						if ( !res->atom ) res->atom = atom;
						atom->sel = 1;
						natom++;
					}
					if ( ( i = il.find(ATOM_RESINSERT) ) >= 0 ) {
						res->insert[0] = ir[i][0];
						if ( !res->insert[0] || res->insert[0] == '?' ) res->insert[0] = ' ';
					}
					if ( ( i = il.find(ATOM_RESIDUE) ) >= 0 ) {
						resname = ir[i].c_str();
						strncpy(res->type, ir[i].c_str(), 6);
					}
					if ( ( i = il.find(ATOM_NUMBER) ) >= 0 )
						atom->num = to_integer(ir[i]);
					if ( ( i = il.find(ATOM_SYMBOL) ) >= 0 )
						strncpy(atom->el, ir[i].c_str(), 2);
					if ( ( i = il.find(ATOM_TOPTYPE) ) >= 0 )
						atom_clean_type(atom, ir[i].c_str());
					if ( ( i = il.find(ATOM_X) ) >= 0 )
						atom->coord[0] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Y) ) >= 0 )
						atom->coord[1] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Z) ) >= 0 )
						atom->coord[2] = to_real(ir[i]);
					atom->q = 1;
					if ( ( i = il.find(ATOM_OCCUPANCY) ) >= 0 )
						atom->q = to_real(ir[i]);
					if ( ( i = il.find(ATOM_BFACTOR) ) >= 0 )
						atom->b = to_real(ir[i]);
					if ( ( i = il.find(ATOM_CHARGE) ) >= 0 )
						atom->chrg = to_real(ir[i]);
				}
			}
		}
		nmol++;
		molname[0]++;
		if ( molname[0] > 'Z' ) molname[0] = 'A';
	}
		
	if ( nmol < 1 ) nmol--;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG read_mol_star: " << nmol << " molecules" << endl;
		cout << "DEBUG read_mol_star: " << nres << " residues" << endl;
		cout << "DEBUG read_mol_star: " << natom << " atoms" << endl;
	}
	
	return nmol;
}

/**
@brief 	Writes a molecule group to a STAR format file.
@param 	&filename		the file name.
@param 	*molgroup		the molecule group.
@return long 				number of molecules written (<0 if writing failed).
**/
long 		write_mol_star(Bstring& filename, Bmolgroup *molgroup)
{
 	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_mol_star: filename=" << filename << endl;

 	Bstar			star;

	star.comment(molgroup->comment.str());
	star.line_length(120);
	
	long			natom(0), nmol(0);
	string			chain("A");
	
	// Convert the atomic data
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		BstarBlock&		block = star.add_block(mol->id.str());
		block[MOLECULE_NAME] = mol->id.str();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG write_mol_star: sequence " << mol->seq << endl;
		if ( ( mol->nbase > 0 ) && mol->naseq.length() ) {
			block[MOLECULE_LENGTH] = to_string(mol->nbase);
			block[MOLECULE_SEQUENCE] = mol->naseq.str();
		}
		if ( ( mol->nres > 0 ) && mol->seq.length() ) {
			block[MOLECULE_LENGTH] = to_string(mol->nres);
			block[MOLECULE_SEQUENCE] = mol->seq.str();
		}
//		chain = mol->id.str();
		natom = 0;
		for ( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) natom++;
		if ( natom ) {
			BstarLoop&		loop = block.add_loop();
			loop.tags()[ATOM_NUMBER] = 0;
			loop.tags()[ATOM_SYMBOL] = 1;
			loop.tags()[ATOM_TOPTYPE] = 2;
			loop.tags()[ATOM_RESIDUE] = 3;
			loop.tags()[ATOM_RESNUMBER] = 4;
			loop.tags()[ATOM_RESINSERT] = 5;
			loop.tags()[ATOM_CHAIN] = 6;
			loop.tags()[ATOM_X] = 7;
			loop.tags()[ATOM_Y] = 8;
			loop.tags()[ATOM_Z] = 9;
			loop.tags()[ATOM_OCCUPANCY] = 10;
			loop.tags()[ATOM_BFACTOR] = 11;
			loop.tags()[ATOM_CHARGE] = 12;
			for ( res = mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					vector<string>&	vs = loop.add_row(13);
					vs[0] = to_string(atom->num);
					vs[1] = atom->el;
					vs[2] = atom->type;
					vs[3] = res->type;
					vs[4] = to_string(res->num);
					if ( !res->insert[0] || res->insert[0] == ' ' ) vs[5] = "?";
					else vs[5] = res->insert;
					vs[6] = chain;
					vs[7] = to_string(atom->coord[0]);
					vs[8] = to_string(atom->coord[1]);
					vs[9] = to_string(atom->coord[2]);
					vs[10] = to_string(atom->q);
					vs[11] = to_string(atom->b);
					vs[12] = to_string(atom->chrg);				}
			}
		}
		nmol++;
		chain[0]++;
	}
			
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_mol_star: " << filename << endl;

	star.write(filename.str());
	
	return nmol;
}
