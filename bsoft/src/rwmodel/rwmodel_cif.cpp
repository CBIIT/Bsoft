/**
@file	rwmodel_cif.cpp
@brief	Library routines to read and write molecule files in CIF format
@author Bernard Heymann
@date	Created: 19980822 
@date	Modified: 20220208
**/

#include "star.h"
#include "rwmolecule.h"
#include "rwmodel_cif.h"
#include "mol_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a CIF format file.
@param 	&filename		the file name.
@param 	&paramfile		parameter file.
@return long				number of molecules read (<0 if reading failed).
**/
Bmodel*		read_model_cif(Bstring& filename, Bstring& paramfile)
{
 	Bstar			star;
	
	if ( verbose )
		cout << "Reading file:                   " << filename << endl;

 	if ( star.read(filename.str()) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return 0;
	}

	long			i, j, nmol(0), nres(0), natom(0);
	string			pdb, typestr, sid, gid, chain, el, atomtype, restype, restype2, resnum, resnum2, seq;
	map<string,int>	sheet_strands;
	map<string,int>	sheet_order;

	Bmodel*			model_list = NULL;
	Bmodel*			model = NULL;
	Bcomponent*		comp = NULL;
//	Bcomponent*		comp2 = NULL;
//	Blink*			link = NULL;
	Bgroup*			group = NULL;

	for ( auto ib: star.blocks() ) {
//		model->identifier(ib.tag());
//		if ( verbose & VERB_DEBUG )
//			cout << "DEBUG read_model_cif: Molecule " << model->identifier() << endl;
/*		molname = ib.tag();
		if ( molname.length() < 1 ) molname = "A";
		mol = molecule_add(&molgroup->mol, molname);
		mol->seq = ib.at(MOLECULE_SEQUENCE).c_str();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_model_cif: sequence " << mol->seq << endl;
		mol->nres = ib.real(MOLECULE_LENGTH);
		if ( mol->seq.length() > 0 )
			if ( mol->nres < 1 ) mol->nres = mol->seq.length();*/
		for ( auto il: ib.loops() ) {
			if ( il.find(MOL_ID) >= 0 ) {
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(MOL_ID) ) >= 0 )
						chain = ir[i];
					if ( model ) model = model->add(chain);
					else model_list = model = new Bmodel(chain);
//					if ( verbose & VERB_DEBUG )
						cout << "DEBUG read_model_cif: Molecule " << model->identifier() << endl;
				}
			}
			if ( il.find(HELIX_ID) >= 0 ) {
//				il.show_tags();
				j = 0;
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(HELIX_ID) ) >= 0 )
						gid = to_string(++j) + " " + ir[i];
					if ( ( i = il.find(HELIX_CHAIN) ) >= 0 )
						chain = ir[i];
					if ( ( i = il.find(HELIX_RESIDUE1) ) >= 0 )
						restype = ir[i];
					if ( ( i = il.find(HELIX_RESIDUE2) ) >= 0 )
						restype2 = ir[i];
					if ( ( i = il.find(HELIX_RESNUM1) ) >= 0 )
						resnum = ir[i];
					if ( ( i = il.find(HELIX_RESNUM2) ) >= 0 )
						resnum2 = ir[i];
					model = model_list->find(chain);
					group = model->add_group(gid);
					group->group_type("HELIX ");
					group->type("1");								// Type of helix
					// Chain residue1 residue2
					group->description(chain + " " + restype + " " + restype2);
					seq = resnum + "-" + resnum2;					// Residue numbers
					group->sequence(seq);
				}
			}
			if ( il.find(SHEET_ORDER_ID) >= 0 ) {
//				il.show_tags();
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(SHEET_ORDER_ID) ) >= 0 )
						sid = ir[i];
					if ( ( i = il.find(SHEET_ORDER_STRAND_ID) ) >= 0 )
						gid = ir[i] + " " + sid;
					if ( ( i = il.find(SHEET_ORDER_SENSE) ) >= 0 ) {
						if ( ir[i].find("anti-parallel") != string::npos ) sheet_order[gid] = -1;
						else sheet_order[gid] = 1;
					}
//					cout << gid << tab << ir[i] << tab << sheet_order[gid] << endl;
				}
			}
			if ( il.find(SHEET_ID) >= 0 ) {
//				il.show_tags();
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(SHEET_ID) ) >= 0 )
						sid = ir[i];
					if ( ( i = il.find(SHEET_STRAND_ID) ) >= 0 )
						gid = ir[i] + " " + sid;
					if ( ( i = il.find(SHEET_CHAIN) ) >= 0 )
						chain = ir[i];
					if ( ( i = il.find(SHEET_RESIDUE1) ) >= 0 )
						restype = ir[i];
					if ( ( i = il.find(SHEET_RESIDUE2) ) >= 0 )
						restype2 = ir[i];
					if ( ( i = il.find(SHEET_RESNUM1) ) >= 0 )
						resnum = ir[i];
					if ( ( i = il.find(SHEET_RESNUM2) ) >= 0 )
						resnum2 = ir[i];
					model = model_list->find(chain);
					group = model->add_group(gid);
					group->group_type("SHEET ");
					group->type(to_string(sheet_order[gid]));		// Direction of strand
					// Chain residue1 residue2 strands_in_sheet
					group->description(chain + " " + restype + " " + restype + " ");
					seq = resnum + "-" + resnum2;					// Residue numbers
					group->sequence(seq);
					if ( sheet_strands.find(sid) != sheet_strands.end() )
						sheet_strands[sid] += 1;
					else
						sheet_strands[sid] = 1;
				}
				for ( model = model_list; model; model = model->next ) {
					for ( group = model->group; group; group = group->next ) {
						if ( group->group_type() == "SHEET " ) {
							sid = group->identifier();
							sid = sid.substr(sid.find(" ")+1);
//							cout << sid << tab << sheet_strands[sid] << endl;
							group->description() += to_string(sheet_strands[sid]);
						}
					}
				}
			}
			if ( il.find(ATOM_NUMBER) >= 0 ) {
//				il.show_tags();
				for ( auto ir: il.data() ) {
//					for ( auto it: ir ) cout << it << endl;
					if ( ( i = il.find(ATOM_CHAIN) ) >= 0 )
						chain = ir[i];
					if ( !model || model->identifier() != chain ) {
						model = model_list->find(chain);
						comp = model->comp;
					}
					if ( ( i = il.find(ATOM_NUMBER) ) >= 0 ) {
						if ( comp ) comp = comp->add(ir[i]);
						else comp = model->comp = new Bcomponent(ir[i]);
						comp->density(1);
						comp->select(1);
						natom++;
					}
					if ( ( i = il.find(ATOM_PDB) ) >= 0 )
						if ( ir[i] == "HETATM" ) comp->select(2);
//					if ( ( i = il.find(ATOM_PDB) ) >= 0 )
//						cout << "-" << ir[i] << "-" << endl;
					if ( ( i = il.find(ATOM_SYMBOL) ) >= 0 )
						el = ir[i];
					if ( ( i = il.find(ATOM_TOPTYPE) ) >= 0 )
						atomtype = ir[i];
					if ( ( i = il.find(ATOM_RESNUMBER) ) >= 0 )
						resnum = ir[i];
//					if ( ( i = il.find(ATOM_RESINSERT) ) >= 0 ) {
//						res->insert[0] = ir[i][0];
//						if ( !res->insert[0] || res->insert[0] == '?' ) res->insert[0] = ' ';
//					}
					if ( ( i = il.find(ATOM_RESIDUE) ) >= 0 )
						restype = ir[i];
					typestr = el + " " + atomtype + " " + restype + " " + chain + " " + resnum;
					comp->type(model->add_type(atomtype));
					comp->description(typestr);
					if ( ( i = il.find(ATOM_X) ) >= 0 )
						comp->location()[0] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Y) ) >= 0 )
						comp->location()[1] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Z) ) >= 0 )
						comp->location()[2] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_OCCUPANCY) ) >= 0 )
						comp->density(to_real(ir[i]));
					if ( ( i = il.find(ATOM_BFACTOR) ) >= 0 )
						comp->FOM(to_real(ir[i]));
					if ( ( i = il.find(ATOM_CHARGE) ) >= 0 )
						comp->charge(to_real(ir[i]));
				}
			}
		}
//		nmol++;
//		molname[0]++;
//		if ( molname[0] > 'Z' ) molname[0] = 'A';
	}
		
	if ( nmol < 1 ) nmol--;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG read_model_cif: " << nmol << " molecules" << endl;
		cout << "DEBUG read_model_cif: " << nres << " residues" << endl;
		cout << "DEBUG read_model_cif: " << natom << " atoms" << endl;
	}
	
	return model_list;
}

/**
@brief 	Writes a molecule group to a CIF format file.
@param 	&filename		the file name.
@param 	*model			the model.
@return long 				number of molecules written (<0 if writing failed).
**/
long 		write_model_cif(Bstring& filename, Bmodel* model)
{
 	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_cif: filename=" << filename << endl;

 	Bstar			star;

//	star.comment(molgroup->comment.str());
	star.line_length(120);
	
//	long			natom(0), nmol(0);
//	string			chain("A");
	
	Bcomponent*		comp = NULL;
//	Bcomponent*		comp2 = NULL;
//	Blink*			link = NULL;
//	Bgroup*			group = NULL;

	BstarBlock&		block = star.add_block(model->identifier());
//	block[MOLECULE_NAME] = mol->id.str();
//	if ( verbose & VERB_DEBUG )
//		cout << "DEBUG write_model_cif: sequence " << mol->seq << endl;

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

	for ( comp = model->comp; comp; comp = comp->next ) {
		vector<string>	vd = split(comp->description());
		vector<string>&	vs = loop.add_row(13);
		vs[0] = comp->identifier();
		vs[1] = vd[0];
		vs[2] = vd[1];
		vs[3] = vd[2];
		vs[4] = vd[4];
//		if ( !res->insert[0] || res->insert[0] == ' ' ) vs[5] = "?";
//		else vs[5] = res->insert;
		vs[5] = "?";
		vs[6] = vd[3];
		vs[7] = to_string(comp->location()[0]);
		vs[8] = to_string(comp->location()[1]);
		vs[9] = to_string(comp->location()[2]);
		vs[10] = to_string(comp->density());
		vs[11] = to_string(comp->FOM());
		vs[12] = to_string(comp->charge());
	}
			
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_cif: " << filename << endl;

	star.write(filename.str());
	
	return 1;
}
