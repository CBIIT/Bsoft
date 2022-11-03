/**
@file	rwPDB.cpp
@brief	Library routines to read and write PDB coordinate files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20110623
**/

#include "rwPDB.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a PDB format coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules read (<0 if reading failed).

	PDB fixed format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	ATOM      1  N   ASP     1     -45.145   1.941 -34.322  1.00  0.00
	ATOM      2  OD2 ASP     1     -44.830   3.862 -30.017  1.00  0.00
	ATOM      3  OD1 ASP     1     -45.547   1.817 -29.798  1.00  0.00

**/
int 	readPDB(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPDB: filename=" << filename << endl;
	
    // Open pdb file read only
    ifstream		fpdb(filename.c_str());
    if ( fpdb.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}
    
    int				readflag, linenum(0);
    int				i, nmol(0), natom(0), nsec(0), nr(0);
    int				resnum, resnum2, prev_resnum = -10000;
    int				atomnum, ba, helixclass;
	char			chain, insert, prev_insert = ' ';
	char			sgstring[32];
	char			label[12], restype[12], atomtype[12], satomnum[8];
    char			aline[MAXLINELEN], *aptr;
	Bstring			molname(" ");
	double			d;
	
	Bmolecule*		mol = NULL;
	Bresidue*		res = NULL;
	Bsecondary*		sec = NULL;
	Batom*			atom = NULL;
	Batom*			atom2 = NULL;
	Bbond*			bond = NULL;
	UnitCell		unitcell = molgroup->unitcell;
	
	Batom**			al = new Batom*[100000];
	memset(al, 0, 100000*sizeof(Batom*));
	
	nmol = natom = 0;
	chain = '%';

    while ( !fpdb.eof() ) {
    	fpdb.getline(aline, MAXLINELEN);
		linenum++;
		readflag = 0;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG readPDB: " << aline << endl;
		if ( strncmp(aline, "TER", 3) == 0 ||  strncmp(aline, "END", 3) == 0 ) chain = '%';
		if ( strncmp(aline, "CRYST1", 6) == 0 ) {
			i = sscanf(aline, "%6c%9lf%9lf%9lf%7lf%7lf%7lf%11c%4d", 
					label, &unitcell[0], &unitcell[1], &unitcell[2], 
					&unitcell[3], &unitcell[4], &unitcell[5], 
					sgstring, &molgroup->spacegroup);
			molgroup->sgstring = sgstring;
//			unitcell.alf += M_PI/180.0;
//			unitcell.bet += M_PI/180.0;
//			unitcell.gam += M_PI/180.0;
			if ( i < 9 ) molgroup->sgstring = "P 1";
		}
		if ( ( strncmp(aline, "ATOM", 4) == 0 ) || ( strncmp(aline, "HETATM", 6) == 0 ) ) {
			if ( molgroup->select.contains("ATOM") ) {
				if ( strstr(aline,"ATOM") ) readflag = 1;
			} else if ( molgroup->select.contains("HETATM") ) {
				if ( strstr(aline,"HETATM") ) readflag = 1;
			} else if ( molgroup->select.contains("CA") ) {
				if ( strncmp(&aline[13],"CA", 2) == 0 ) readflag = 1;
			} else readflag = 1;
			if ( strstr(aline, "nan") || strstr(aline, "NAN") ) {
				cerr << "Warning: NAN at line number " << linenum << endl;
				readflag = 0;
			}
			if ( aline[16] == 'B' ) readflag = 0;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG readPDB: readflag=" << readflag << endl;
			if ( readflag ) {
//				if ( !isspace(aline[21]) && aline[21] != chain ) {
				if ( aline[21] != chain ) {
					chain = aline[21];
					molname[0] = chain;
					mol = molecule_add(&mol, molname);
					if ( !molgroup->mol ) molgroup->mol = mol;
					nmol++;
					res = NULL;
					prev_resnum = -10000;
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG readPDB: chain=" << chain << " nmol=" << nmol << endl;
				}
				memset(satomnum, 0, 8);
				atomnum = get_integer(aline + 6, 5);
				if ( atomnum < 1 ) atomnum = natom + 1;
				strncpy(atomtype, aline + 12, 4);
				strncpy(restype, aline + 17, 3);
				resnum = get_integer(aline + 22, 4);
//				if ( resnum < 1 ) resnum = natom + 1;
				insert = aline[26];
				restype[3] = 0; 	// Residue names only 3 characters
				atomtype[4] = 0; 	// Atom names only 4 characters
				if ( natom < 1 || resnum != prev_resnum || insert != prev_insert ) {
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG readPDB: resnum=" << resnum << " restype=" << restype << endl;
					if ( !mol ) {
						error_show("No molecule structure allocated!", __FILE__, __LINE__);
						return -1;
					}
					res = residue_add(&res, restype);
					if ( !mol->res ) mol->res = res;
					res->num = resnum;
					res->insert[0] = insert;
					prev_resnum = resnum;
					prev_insert = insert;
					atom = NULL;
				}
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readPDB: atomnum=" << atomnum << " atomtype=" << atomtype << endl;
	    		atom = atom_add(&atom, atomtype);
	    		if ( !res->atom ) res->atom = atom;
				atom->num = atomnum;
				atom->coord[0] = get_float(aline + 30, 8);
				atom->coord[1] = get_float(aline + 38, 8);
				atom->coord[2] = get_float(aline + 46, 8);
				atom->q = get_float(aline + 54, 6);
				atom->b = get_float(aline + 60, 6);
				atom->chrg = get_float(aline + 78, 6);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readPDB: " << atom->num << " " << res->num << " " << atom->coord << " " << atom->chrg << endl; 
		    	atom->sel = 1;			/* All atoms read selected */
				if ( strncmp(aline, "HETATM", 6) == 0 ) atom->sel = 2;
				al[atomnum] = atom;		/* Atom list */
				if ( natom <= atomnum ) natom = atomnum + 1;
				nr++;
			}
		}
		if ( strncmp(aline, "CONECT", 6) == 0 ) {   // Assumed to always be after all ATOM records
			atomnum = ba = get_integer(aline + 6, 5);
			if ( atomnum < natom ) for ( aptr=aline+11; ba && strlen(aptr) > 4; aptr+=5 ) {
				ba = get_integer(aptr, 5);
				if ( ba <= natom && ba > atomnum && al[atomnum] && al[ba] ) {
					if ( verbose & VERB_DEBUG )
						printf ("DEBUG readPDB: Atom %d bound to %d\n", atomnum, ba);
					atom = al[atomnum];
					atom2 = al[ba];
					d = atom->coord.distance(atom2->coord);
					bond = bond_add(&bond, atom, atom2, d, 1);
					if ( !molgroup->bond ) molgroup->bond = bond;
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readPDB: atoms read = " << nr << endl;

	delete[] al;
	
	if ( natom ) nmol++;
	else nmol = -1;

//	fseek(fpdb, 0, SEEK_SET);
	fpdb.clear();
	fpdb.seekg(0, ios::beg);
	
	mol = molgroup->mol;

	chain = 0;
	sec = NULL;
    while ( !fpdb.eof() ) {
    	fpdb.getline(aline, MAXLINELEN);
		if ( strncmp(aline, "HELIX", 5) == 0 ) {
			chain = aline[19];
//			if ( sec ) {
//				while ( mol && chain != mol->id[0] ) mol = mol->next;
//				if ( !mol )
					for ( mol = molgroup->mol; mol && chain != mol->id[0]; mol = mol->next ) ;
//			}
			if ( mol ) {
				sec = (Bsecondary *) add_item((char **) &mol->sec, sizeof(Bsecondary));
				sec->num = get_integer(aline + 6, 4);
				strncpy(sec->id, aline + 11, 3);
				helixclass = get_integer(aline + 38, 2);
				switch ( helixclass ) {
					case 1: sec->type = RightHandedAlpha; break;
					case 2: sec->type = RightHandedOmega; break;
					case 3: sec->type = RightHandedPi; break;
					case 4: sec->type = RightHandedGamma; break;
					case 5: sec->type = RightHanded310; break;
					case 6: sec->type = LeftHandedAlpha; break;
					case 7: sec->type = LeftHandedOmega; break;
					case 8: sec->type = LeftHandedGamma; break;
					case 9: sec->type = TwoSevenRibbon; break;
					case 10: sec->type = Polyproline; break;
				}
				resnum = get_integer(aline + 21, 4);
				resnum2 = get_integer(aline + 33, 4);
				for ( res=mol->res; res; res=res->next ) {
					if ( res->num == resnum )  sec->first = res;
					if ( res->num == resnum2 ) sec->last = res;
				}
				nsec++;
				if ( !sec->first )
					cerr << "Warning: First residue " << resnum << " not found!" << endl;
				if ( !sec->last )
					cerr << "Warning: Last residue " << resnum2 << " not found!" << endl;
			}
		}
		if ( ( strncmp(aline, "ATOM", 4) == 0 ) || ( strncmp(aline, "HETATM", 6) == 0 ) )
			break;
    }

//	fseek(fpdb, 0, SEEK_SET);
	fpdb.clear();
	fpdb.seekg(0, ios::beg);
	
	mol = molgroup->mol;

	chain = 0;
	sec = NULL;
    while ( !fpdb.eof() ) {
    	fpdb.getline(aline, MAXLINELEN);
		if ( strncmp(aline, "SHEET", 5) == 0 ) {
			chain = aline[21];
//			if ( sec ) {
//				while ( mol && chain != mol->id[0] ) mol = mol->next;
//				if ( !mol )
					for ( mol = molgroup->mol; mol && chain != mol->id[0]; mol = mol->next ) ;
//			}
			if ( mol ) {
				sec = (Bsecondary *) add_item((char **) &mol->sec, sizeof(Bsecondary));
				sec->num = get_integer(aline + 6, 4);
				strncpy(sec->id, aline + 11, 3);
				sec->type = Strand;
				sec->sense = get_integer(aline + 38, 2);
				sec->nstrands = get_integer(aline + 14, 2);
				resnum = get_integer(aline + 22, 4);
				resnum2 = get_integer(aline + 33, 4);
				for ( res=mol->res; res; res=res->next ) {
					if ( res->num == resnum )  sec->first = res;
					if ( res->num == resnum2 ) sec->last = res;
				}
				nsec++;
				if ( !sec->first )
					cerr << "Warning: First residue " << resnum << " not found!" << endl;
				if ( !sec->last )
					cerr << "Warning: Last residue " << resnum2 << " not found!" << endl;
			}
		}
		if ( ( strncmp(aline, "ATOM", 4) == 0 ) || ( strncmp(aline, "HETATM", 6) == 0 ) )
			break;
    }
	
// 	fseek(fpdb, 0, SEEK_SET);
	fpdb.clear();
	fpdb.seekg(0, ios::beg);
	
	mol = molgroup->mol;

	chain = 0;
	sec = NULL;
    while ( !fpdb.eof() ) {
    	fpdb.getline(aline, MAXLINELEN);
		if ( strncmp(aline, "TURN", 4) == 0 ) {
			chain = aline[19];
//			if ( sec ) {
//				while ( mol && chain != mol->id[0] ) mol = mol->next;
//				if ( !mol )
					for ( mol = molgroup->mol; mol && chain != mol->id[0]; mol = mol->next ) ;
//			}
			if ( mol ) {
				sec = (Bsecondary *) add_item((char **) &mol->sec, sizeof(Bsecondary));
				sec->num = get_integer(aline + 7, 4);
				strncpy(sec->id, aline + 11, 3);
				sec->type = Turn;
				resnum = get_integer(aline + 20, 4);
				resnum2 = get_integer(aline + 31, 4);
				for ( res=mol->res; res; res=res->next ) {
					if ( res->num == resnum )  sec->first = res;
					if ( res->num == resnum2 ) sec->last = res;
				}
				nsec++;
			}
		}
		if ( ( strncmp(aline, "ATOM", 4) == 0 ) || ( strncmp(aline, "HETATM", 6) == 0 ) )
			break;
    }
	
	fpdb.close();
	
	molgroup_stats(molgroup);

	if ( verbose & VERB_DEBUG ) {
    	cout << "DEBUG readPDB: Number of molecules  = " << nmol << endl; 
		cout << "DEBUG readPDB: Secondary structures = " << nsec << endl;
	}
	
    return nmol;
}

/**
@brief 	Writes a PDB format coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules written (<0 if writing failed).

	PDB fixed format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	ATOM      1  N   ASP     1     -45.145   1.941 -34.322  1.00  0.00
	ATOM      2  OD2 ASP     1     -44.830   3.862 -30.017  1.00  0.00
	ATOM      3  OD1 ASP     1     -45.547   1.817 -29.798  1.00  0.00

**/
int 	writePDB(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: filename=" << filename << endl;
	
	if ( !molgroup ) return -1;
	if ( !molgroup->mol ) return -1;
	if ( !molgroup->mol->res ) return -1;
	
    ofstream		fpdb(filename.c_str());
    if ( fpdb.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    int     	i, j, nw(0), nmol(0), a1, a2, pnum;
	char		atomtag[8], atomtype[12];
	char		alt[4] = " ", sid[8] = "    ", esm[4] = "  ", chain;
    
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: New file name:                  " << filename << endl;
	
//    if ( ( fpdb = fopen(filename.c_str(), "w") ) == NULL ) return -1;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: File " << filename << " opened" << endl;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Bsecondary*		sec;
	Batom*			atom;
	Bbond*			bond;
	UnitCell		unitcell = molgroup->unitcell;
	
	time_t			t = time(NULL);
	tm*				ltm = localtime(&t);
	char			stm[16];
	strftime(stm, sizeof(stm), "%d-%b-%y ", ltm);
	Bstring			title = command_line();
//	fpdb << left << setw(50) << "HEADER    Written by Bsoft" << asctime(localtime(&t));
//	fpdb << left << setw(50) << "HEADER    Written by Bsoft" <<
//		ltm->tm_mday << "-" << ltm->tm_mon << "-" << ltm->tm_year << endl;
	fpdb << left << setw(50) << "HEADER    Written by Bsoft" << stm << molgroup->id << endl;
	fpdb << "TITLE     " << title << endl;
//	if ( molgroup->comment.length() )
//		fpdb << "REMARK " << molgroup->comment << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: remarks done" << endl;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		chain = mol->id[0];
		for ( sec = mol->sec; sec; sec = sec->next ) {
			if ( sec->type < Strand && sec->first && sec->last ) {
//				fpdb << "HELIX %4d %3s %3s %c %4d%1s %3s %c %4d%1s%2d\n",
				fpdb << "HELIX " << right << setw(4) << sec->num << " " <<
					left << setw(4) << sec->id <<
					setw(4) << sec->first->type << chain <<
					right << setw(5) << sec->first->num <<
					left << sec->first->insert << " " <<
					setw(4) << sec->last->type << chain <<
					right << setw(5) << sec->last->num <<
					left << sec->last->insert <<
					right << setw(2) << sec->type << endl;
			}
		}
	}	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: helices done" << endl;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		chain = mol->id[0];
		for ( sec = mol->sec; sec; sec = sec->next ) {
			if ( sec->type == Strand && sec->first && sec->last ) {
//				fpdb << "SHEET %4d %3s%2d %3s %c%4d%1s %3s %c%4d%1s%2d\n",
				fpdb << "SHEET " << right << setw(4) << sec->num << " " <<
					left << setw(3) << sec->id <<
					right << setw(2) << sec->nstrands << " " <<
					left << setw(4) << sec->first->type << chain <<
					right << setw(4) << sec->first->num <<
					left << sec->first->insert << " " <<
					setw(4) << sec->last->type << chain <<
					right << setw(4) << sec->last->num <<
					left << sec->last->insert <<
					right << setw(2) << sec->sense << endl;
			}
		}
	}	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: sheets done" << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		chain = mol->id[0];
		for ( sec = mol->sec; sec; sec = sec->next ) {
			if ( sec->type == Turn ) {
//				fpdb << "TURN  %4d %3s %3s %c%4d%1s %3s %c%4d%1s\n",
				fpdb << "TURN  " << right << setw(4) << sec->num << " " <<
					left << setw(4) << sec->id <<
					setw(4) << sec->first->type << chain <<
					right << setw(4) << sec->first->num <<
					left << sec->first->insert << " " <<
					setw(4) << sec->last->type << chain <<
					right << setw(4) << sec->last->num <<
					left << sec->last->insert << endl;
			}
		}
	}	

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: turns done" << endl;
	
//    fpdb << "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%11s%4d\n",
    fpdb << "CRYST1" << fixed << setprecision(3) << setw(9) << unitcell.a() << setw(9) << 
		unitcell.b() << setw(9) << unitcell.c() << setprecision(2) << setw(7) << 
		unitcell.alpha()*180.0/M_PI << setw(7) << unitcell.beta()*180.0/M_PI << 
		setw(7) << unitcell.gamma()*180.0/M_PI << setw(11) << 
		molgroup->sgstring << setw(4) << molgroup->spacegroup << endl;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: crystal parameters done" << endl;

	pnum = molgroup->mol->res->atom->num - 1;
//	cout << "First atom number = " << pnum << endl;
	strcpy(sid, "    ");
	sid[4] = 0;
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		chain = mol->id[0];
		nw = 0;
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( atom->num <= pnum ) atom->num = pnum + 1;
				if ( atom->num > 99999 ) atom->num = 1;
				pnum = atom->num;
				if ( atom->sel ) {
					if ( atom->sel < 2 ) strncpy(atomtag, "ATOM  ", 6);
					else strncpy(atomtag, "HETATM", 6);
					atomtag[6] = 0;
					if ( strlen(atom->type) < 4 )
						snprintf(atomtype, 12, " %s", atom->type);
					else
						snprintf(atomtype, 12, "%s", atom->type);
					atomtype[4] = 0;
					strncpy(esm, atom->el, 2);
					esm[2] = 0;
					if ( strlen(esm) < 2 ) esm[1] = ' ';
//					fpdb << "%6s%5d %-4s%1s%3s %c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%6.3f\n",
					fpdb << setw(6) <<
						atomtag << right << setw(5) << atom->num << " " << left << setw(4) << atomtype << alt <<
						setw(3) << res->type << " " << chain << right << setw(4) << res->num << res->insert << "   " <<
						setprecision(3) << setw(8) << atom->coord[0] << setw(8) << 
						atom->coord[1] << setw(8) << atom->coord[2] << 
						setprecision(2) << setw(6) << atom->q << setw(6) << atom->b << "      " <<
						setw(4) << sid << right << setw(2) << esm << setprecision(3) << setw(6) << atom->chrg << endl;
					nw++;
				}
			}
		}
		if ( nw > 0 ) {
			nmol++;
    		fpdb << "TER" << endl;
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: atoms done: " << nw << endl;
/*
	int*			al = new int[1000000];
	for ( i=0; i<1000000; i++ ) al[i] = 0;
	
	for ( bond = molgroup->bond; bond; bond = bond->next ) {
		a1 = 10*bond->atom1->num;
		a2 = 10*bond->atom2->num;
		for ( i=0; i<10 && al[a1+i]; i++ ) ;
		if ( i<10 ) al[a1+i] = bond->atom2->num;
		for ( i=0; i<10 && al[a2+i]; i++ ) ;
		if ( i<10 ) al[a2+i] = bond->atom1->num;
	}
	
	for ( i=1; i<100000; i++ ) if ( al[i*10] ) {
		fpdb << "CONECT" << setw(5) << i;
		for ( j=0; j<10 && al[i*10+j]; j++ ) {
			if ( j == 4 ) fpdb << endl << "CONECT" << setw(5) << i;
			fpdb << setw(5) << al[i*10+j];
		}
		fpdb << endl;
	}
	
	delete[] al;
*/
//	vector<vector<int>> cnct(100000, vector<int>(10, 0));
	vector<vector<int>> cnct(100000);
	
	for ( bond = molgroup->bond; bond; bond = bond->next ) {
		a1 = bond->atom1->num;
		a2 = bond->atom2->num;
		cnct[a1].push_back(a2);
		cnct[a2].push_back(a1);
	}
	
	for ( i=1; i<100000; ++i ) if ( cnct[i].size() ) {
		fpdb << "CONECT" << setw(5) << i;
		for ( j=0; j<10 && j<cnct[i].size(); ++j ) {
			if ( j == 4 ) fpdb << endl << "CONECT" << setw(5) << i;
			fpdb << setw(5) << cnct[i][j];
		}
		fpdb << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: bonds done" << endl;
  
	fpdb << "END" << endl;
	    
    fpdb.close();

//	molgroup_stats(molgroup);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writePDB: File " << filename << " closed" << endl;
	
    return nmol;
}

