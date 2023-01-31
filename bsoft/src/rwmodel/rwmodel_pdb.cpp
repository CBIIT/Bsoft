/**
@file	rwmodel_pdb.cpp
@brief	Library routines to read and write PDB coordinate files
@author Bernard Heymann
@date	Created: 20211231
@date	Modified: 20220324
**/

#include "rwmodel.h"
#include "model_links.h"
#include "model_util.h"
#include "linked_list.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bmodel* 	readPDB(Bstring& filename, long n);
int			writePDB(Bstring& filename, Bmodel* model);

/**
@brief 	Reads molecular model parameters.
@param 	*file_list		list of model parameter file names.
@param 	&paramfile		parameter file.
@return Bmodel*			model parameters.

	Each file is considered a separate model.
	
	The component description field contains the following detail as a space-separated list:
		element
		atom type
		residue type
		chain id
		residue number

**/
Bmodel*		read_model_pdb(Bstring* file_list, Bstring& paramfile)
{
	Bstring*		filename;
    Bstring    		atom_select("all");
	Bstring			id;
	
	int				i;
	RGBA<float>		rgba(1,1,1,1);
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_molecule: " << *file_list << endl;
	
	for ( i=1, filename = file_list; filename; filename = filename->next, i++ ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << *filename << endl;
		mp = readPDB(*filename, i);
		if ( mp ) {
			if ( model ) model->add(mp);
			else model = mp;
		}
//		mp->model_type(mp->identifier());
	}
	
	for ( i=0, mp = model; mp; mp = mp->next ) i++;
	
	if ( verbose )
		cout << "Models read: " << i << endl;

	models_process(model, model_setup_links);
	
	return model;
}

Bmodel* 	readPDB(Bstring& filename, long n)
{
	if ( verbose )
		cout << "Reading file:                   " << filename << endl;


	// Open pdb file read only
	ifstream		fmod(filename.c_str());
	if ( fmod.fail() ) {
		cerr << "Error: File " << filename << " not opened!" << endl;
		return NULL;
	}

	string			s, recordname, chain("-"), grpch, resnumstr;
	int				i, nmol(0), natom(0), ngrp(0);
	int				resnum, atomnum, ba;
//	char			insert, prev_insert = ' ';
	string			gid, seq, el, restype, atomtype, typestr;
//	Bstring			molname(" ");
	double			d;
	vector<Bcomponent*>	al(100000);

//	UnitCell		unitcell = molgroup->unitcell;

	Bmodel*			model = new Bmodel(n);
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	Bgroup*			group = NULL;
	
	while ( !fmod.eof() ) {
		getline(fmod, s);
//		cout << s << tab << s.length() << endl;
		recordname = s.substr(0,6);
//		s = s.substr(6);
//		if ( recordname == "HEADER" ) {
//			model->identifier(s.substr(60));
//		}
		if ( recordname == "HELIX " ) {
			grpch = s[19];
//			gid = s.substr(11, 3);
			gid = s.substr(6, 4) + " " + s.substr(11, 3);
			if ( verbose & VERB_FULL )
				cout << "Adding helix: " << gid << endl;
			group = model->add_group(gid);
			group->group_type(recordname);
			group->type(s.substr(38, 2));					// Type of helix
			// Chain residue1 residue2
			group->description(grpch + " " + s.substr(15, 3) + " " + s.substr(27, 3));
			seq = s.substr(21, 4) + "-" + s.substr(33, 4);	// Residue numbers
			group->sequence(seq);
//				sec->num = to_integer(s.substr(6, 4));
			ngrp++;
		}
		if ( recordname == "SHEET " ) {
			grpch = s[21];
			gid = s.substr(6, 4) + " " + s.substr(11, 3);
			if ( verbose & VERB_FULL )
				cout << "Adding strand: " << gid << endl;
			group = model->add_group(gid);
			group->group_type(recordname);
			group->type(s.substr(38, 2));					// Direction of strand
			// Chain residue1 residue2 strands_in_sheet
			group->description(grpch + " " + s.substr(17, 3) + " " + s.substr(28, 3) + " " + s.substr(14, 2));
			seq = s.substr(22, 4) + "-" + s.substr(33, 4);	// Residue numbers
			group->sequence(seq);
			ngrp++;
		}
		if ( recordname == "TURN  " ) {
			grpch = s[19];
			gid = s.substr(7, 4) + s.substr(11, 3);
			if ( verbose & VERB_FULL )
				cout << "Adding turn: " << gid << endl;
			group = model->add_group(gid);
			group->group_type(recordname);
//			group->type(s.substr(38, 2));					// ??
			group->description(grpch);
			seq = s.substr(20, 4) + "-" + s.substr(31, 4);	// Residue numbers
			group->sequence(seq);
			ngrp++;
		}
		if ( ( recordname == "ATOM  " || recordname == "HETATM" ) && s[16] != 'B' ) {
			if ( s[21] != chain[0] ) {
				if ( natom ) {
					seq += "-" + resnumstr;
					group->sequence(seq);
					if ( verbose & VERB_FULL )
						cout << seq << endl;
				}
				chain = s[21];
				if ( verbose & VERB_FULL )
					cout << "Adding molecule: " << chain << endl;
				group = model->add_group(chain);
				group->group_type("CHAIN");
//				group->type(s.substr(38, 2));	// Type of molecule?
				seq = s.substr(22, 4);			// First residue
				nmol++;
				ngrp++;
			}
			atomnum = to_integer(s.substr(6, 5));
			if ( atomnum < 1 ) atomnum = natom + 1;
			resnumstr = s.substr(22, 4);
			resnum = to_integer(resnumstr);
			if ( resnum < 1 ) resnum = natom + 1;
			if ( comp ) comp = comp->add(atomnum);
			else comp = model->comp = new Bcomponent(atomnum);
			comp->location()[0] = to_real(s.substr(30, 8));
			comp->location()[1] = to_real(s.substr(38, 8));
			comp->location()[2] = to_real(s.substr(46, 8));
			comp->density(to_real(s.substr(54, 6)));		// occupancy
			comp->charge(to_real(s.substr(78, 6)));			// charge
			comp->FOM(to_real(s.substr(60, 6)));			// b
			comp->select(1);
			if ( recordname == "HETATM" ) comp->select(2);
			atomtype = s.substr(12, 4);
			atomtype = remove_spaces2(atomtype);
			el = s.substr(76, 2);
			if ( el[0] == ' ' ) el = el.substr(1);
			restype = s.substr(17, 3);
//			typestr = s.substr(12, 14);
			typestr = el + " " + atomtype + " " + restype + " " + chain + " " + resnumstr;
			comp->type(model->add_type(atomtype));
			comp->description(typestr);

//			insert = s[26];
//			if ( s.length() > 78 ) atom->chrg = to_integer(s.substr(78, 6));
//			if ( verbose & VERB_DEBUG )
//				cout << "DEBUG readPDB: " << atom->num << " " << res->num << " " << atom->coord << " " << atom->chrg << endl;
			al[atomnum] = comp;		/* Atom list */
			if ( natom <= atomnum ) natom = atomnum + 1;
		} else if ( recordname == "CONECT" ) {  // Assumed to always be after all ATOM records
			atomnum = ba = to_integer(s.substr(6, 5));
//			cout << ba << endl;
			if ( atomnum < natom ) for ( i=11; ba && s.length() > i + 4; i+=5 ) {
				ba = to_integer(s.substr(i, 5));
				if ( ba <= natom && ba > atomnum && al[atomnum] && al[ba] ) {
					if ( verbose & VERB_DEBUG )
						printf ("DEBUG readPDB: Atom %d bound to %d\n", atomnum, ba);
					comp = al[atomnum];
					comp2 = al[ba];
					d = comp->location().distance(comp2->location());
					if ( link ) link = link->add(comp, comp2, d, 1);
					else link = model->link = new Blink(comp, comp2, d, 1);
				}
			}
/*		} else if ( recordname == "CRYST1" ) {
			unitcell[0] = to_real(s.substr(6, 9));
			unitcell[1] = to_real(s.substr(15, 9));
			unitcell[2] = to_real(s.substr(24, 9));
			unitcell[3] = to_real(s.substr(33, 7));
			unitcell[4] = to_real(s.substr(40, 7));
			unitcell[5] = to_real(s.substr(47, 7));
			molgroup->sgstring = s.substr(54, 11).c_str();
			molgroup->spacegroup = to_integer(s.substr(65, 4));*/
		}
	}

	// Final residue
	seq += "-" + resnumstr;
	group->sequence(seq);
	if ( verbose & VERB_FULL )
		cout << seq << endl;

	fmod.close();

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG readPDB: Number of molecules  = " << nmol << endl;
		cout << "DEBUG readPDB: Number of groups = " << ngrp << endl;
	}

	return model;
}


/**
@brief 	Writes molecular model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	split		flag to split into individual files and number of digits for the insert.
@return int			files written.

	if ( split=0 ) only the first model is written
**/
int			write_model_pdb(Bstring& filename, Bmodel* model, int split)
{
	int				n(1);
	Bstring			name(filename);
	Bmodel*			mp = NULL;
	char			format[32];
	
	if ( split ) {
		snprintf(format, 32, "_%%0%dd.", split);
 		for ( n=1, mp = model; mp; mp = mp->next, n++ ) {
			if ( model->next )
				name = filename.pre_rev('.') + Bstring(n, format) + filename.post_rev('.');
			writePDB(name, mp);
		}
	} else {
		writePDB(name, model);
	}
	
	return  n;
}

int			writePDB(Bstring& filename, Bmodel* model)
{
	if ( verbose )
		cout << "Writing file:                   " << filename << endl;

	int				i, j, nw(0);
	string			atomtag, atomtype("?");
	string			chain("A"), restype("UNK"), resinsert(" "), alt(" "), el, sid("    ");
	long			atomnum(0), resnum(1), a1, a2;
	
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	Bgroup*			group = NULL;

	ofstream		fmod(filename.c_str());
	if ( fmod.fail() ) return 0;

	time_t			ti = time(NULL);
	tm 				tm = *localtime(&ti);
    
	Bstring			title = command_line();
//	fmod << left << setw(50) << "HEADER    Written by Bsoft" << asctime(localtime(&ti));
	fmod << left << setw(50) << "HEADER    Written by Bsoft" << put_time(&tm, "%d-%b-%y ") << model->identifier() << endl;
	fmod << "TITLE     " << title << endl;

	for ( group = model->group; group; group = group->next ) {
//		cout << group->identifier() << tab << group->description() << tab <<
//			group->group_type() << tab << group->type() << tab << group->sequence() << endl;
		vector<string>	vi = split(group->identifier());
		vector<string>	vs = split(group->description());
		vector<string>	vq = split(group->sequence(), '-');
//		cout << vq[0] << tab << vq[1] << endl;
		if ( group->group_type().find("HELIX") != string::npos ) {
			fmod << "HELIX " << right << setw(4) << vi[0] << " " <<
				left << setw(4) << vi[1] <<
				left << setw(4) << vs[1] << vs[0] <<
				right << setw(5) << vq[0] <<
				left << "  " <<
				setw(4) << vs[2] << vs[0] <<
				right << setw(5) << vq[1] <<
				left << " " <<
				right << setw(2) << group->type() << endl;
		}
		if ( group->group_type().find("SHEET") != string::npos ) {
			fmod << "SHEET " << right << setw(4) << vi[0] << " " <<
				right << setw(3) << vi[1] <<
				right << setw(2) << vs[3] << " " <<
				left << setw(4) << vs[1] << vs[0] <<
				right << setw(4) << vq[0] <<
				left << "  " <<
				setw(4) << vs[2] << vs[0] <<
				right << setw(4) << vq[1] <<
				left << " " <<
				right << setw(2) << group->type() << endl;
		}
		if ( group->group_type().find("TURN") != string::npos ) {
			fmod << "TURN  " << right << setw(4) << vi[0] << " " <<
				left << setw(4) << vi[1] <<
				setw(4) << vs[1] << vs[0] <<
				right << setw(4) << vq[0] <<
				left << "  " <<
				setw(4) << vs[2] << vs[0] <<
				right << setw(4) << vq[1] <<
				left << " " << endl;
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: groups done" << endl;

	for ( comp = model->comp; comp; comp = comp->next ) {
//		chain = mol->id[0];
		if ( comp->select() ) {
			atomnum = stoi(comp->identifier());
			if ( comp->type() )
				atomtype = comp->type()->identifier().substr(0,4);
			vector<string>	vs = split(comp->description());
			if ( vs.size() > 0 ) el = vs[0];
			if ( vs.size() > 2 ) restype = vs[2];
			if ( vs.size() > 3 ) chain = vs[3];
			if ( vs.size() > 4 ) resnum = to_integer(vs[4]);
			if ( comp->select() < 2 ) {
				atomtag = "ATOM  ";
				if ( el.length() < 1 && comp->type() )
//					el = comp->type()->identifier().substr(0,1) + " ";
					el = " " + comp->type()->identifier().substr(0,1);
			} else {
				atomtag = "HETATM";
				if ( el.length() < 1 && comp->type() )
					el = comp->type()->identifier().substr(0,2);
			}
			fmod << setw(6) <<
				atomtag << right << setw(5) << atomnum << " " << left << setw(4) << atomtype << alt <<
				setw(3) << restype << " " << chain << right << setw(4) << resnum << resinsert << "   " <<
				fixed << setprecision(3) << setw(8) << comp->location()[0] << setw(8) <<
				comp->location()[1] << setw(8) << comp->location()[2] <<
				setprecision(2) << setw(6) << comp->density() << setw(6) << comp->FOM() << "      " <<
				setw(4) << sid << right << setw(2) << el << setprecision(1) << setw(6) << comp->charge() << endl;
			nw++;
		}
	}

	vector<vector<int>> cnct(100000);
	
	for ( link = model->link; link; link = link->next ) {
		a1 = stoi((link->comp[0]->identifier()));
		a2 = stoi((link->comp[1]->identifier()));
		cnct[a1].push_back(a2);
		cnct[a2].push_back(a1);
	}
	
	for ( i=1; i<100000; ++i ) if ( cnct[i].size() ) {
		fmod << "CONECT" << setw(5) << i;
		for ( j=0; j<10 && j<cnct[i].size(); ++j ) {
			if ( j == 4 ) fmod << endl << "CONECT" << setw(5) << i;
			fmod << setw(5) << cnct[i][j];
		}
		fmod << endl;
	}

	fmod << "TER" << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_pdb: atoms written: " << nw << endl;

	fmod.close();
	
	return (nw>0);
}

