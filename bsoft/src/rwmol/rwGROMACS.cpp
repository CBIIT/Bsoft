/**
@file	rwGROMACS.cpp
@brief	Library routines to read and write GROMACS coordinate files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20060121
**/

#include "rwGROMACS.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a Gromacs coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules read (<0 if reading failed).

	Gromacs format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	MD of 2 waters, t= 0.0
	    6
	    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
	    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
	    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180

**/
int 	readGROMACS(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readGROMACS: filename=" << filename << endl;
	
    ifstream		fgro(filename.c_str());
    if ( fgro.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

	size_t			i;
	long			j;
    int				readflag(0), atomnum, resnum, prev_resnum = -1;
    int				n(0), natom, nmol(0);
    char			aline[MAXLINELEN];
    char			restype[20], atomtype[20];
	Vector3<double>	coord, vel;
	
    fgro.getline(aline, MAXLINELEN);
	aline[MAXLINELEN-1] = 0;
	
	Bmolecule*		mol = molecule_add(&molgroup->mol, aline);
	
    if ( verbose )
		cout << "Title:                          " << mol->id << endl;
		
    fgro.getline(aline, MAXLINELEN);
    sscanf(aline, "%d", &natom);
	
    if ( verbose )
		cout << "Number of atoms:                " << natom << endl;
    
	Bresidue*	res = NULL;
	Batom*  	atom = NULL;
    int     	linenum = 2;
	
    while ( !fgro.eof() && n < natom ) {
    	fgro.getline(aline, MAXLINELEN);
		linenum++;
		i = strcspn(aline,"\n");
		if ( i > 0 ) aline[i] = ' ';
    	if ( i > 40 ) readflag = 1;
		if ( strstr(aline, "nan") || strstr(aline, "NAN") ) {
			cerr << "Warning: NAN at line number " << linenum << endl;
			readflag = 0;
		}
		if ( readflag ) {
			coord = 0;
			vel = 0;
			sscanf(aline, "%5d%5c%5c%5d%8lf%8lf%8lf%8lf%8lf%8lf", 
				&resnum, restype, atomtype, &atomnum,
				&coord[0], &coord[1], &coord[2], &vel[0], &vel[1], &vel[2]);
			if ( resnum < 1 ) resnum = n + 1;
			if ( atomnum < 1 ) atomnum = n + 1;
	    	restype[5] = 0; 	// Residue names only 5 characters
			if ( resnum != prev_resnum ) {
				if ( mol->res ) res = residue_add(&res, restype);
				else res = residue_add(&mol->res, restype);
				res->num = prev_resnum = resnum;
			}
	    	atomtype[5] = 0; 	// Atom names only 5 characters
			for ( i=0, j=0; i<strlen(atomtype); i++ )
				if ( !isspace(atomtype[i]) ) atomtype[j++] = atomtype[i];
			atomtype[j] = 0;	// Atom type string shifted to left
			if ( !res->atom ) atom = atom_add(&res->atom, atomtype);
			else atom = atom_add(&atom, atomtype);
			atom->num = atomnum;
			atom->coord = coord * 10;
			atom->vel = vel * 10;
    	    atom->q = 1;
    	    atom->b = 10;
			if ( verbose & VERB_DEBUG )
    	    	cout << atom->num << " " << res->num << " " << atom->coord << endl; 

	    	atom->sel = 1;	    /* All atoms read selected */
    	    n++;
		}
    }
	
	if ( n ) nmol = 1;
 	else nmol = -1;
   
	fgro.getline(aline, MAXLINELEN);

	Vector3<double>		box;
    sscanf( aline, "%lf %lf %lf", &box[0], &box[1], &box[2]);
    if ( verbose & VERB_PROCESS )
		cout << "Box size:                       " << box << " nm" << endl;
    
	molgroup->box = box * 10;
    
    fgro.close();
    
    return nmol;
}

/**
@brief 	Writes a Gromacs coordinate file.
@param 	&filename	coordinate file name.
@param	*molgroup molecule group.
@return int 				number of molecules written (<0 if writing failed).

	Gromacs format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	MD of 2 waters, t= 0.0
	    6
	    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
	    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
	    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180

**/
int 	writeGROMACS(Bstring& filename, Bmolgroup *molgroup)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGROMACS: filename=" << filename << endl;
	
	if ( !molgroup ) return -1;
	
    ofstream		fgro(filename.c_str());
    if ( fgro.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    long			i, j, nw(0), natom(0);
	char			atomtype[8];
	memset(atomtype, 0, 8);
    
	if ( verbose & VERB_DEBUG )
		cout << "New file name:                  " << filename << endl;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for ( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				natom++;

//  time_t  	t = time(NULL);
//	fgro << "Gromacs file written by Bsoft: " << asctime(localtime(&t)) << endl;
	char		comment[MAXLINELEN];
	for ( i=0; i<molgroup->comment.length() && molgroup->comment[i] != '\n' && i<80; i++ ) 
		comment[i] = molgroup->comment[i];
	comment[i++] = '\n';
	comment[i] = 0;
	fgro << comment;
    fgro << " " << natom << endl;
    
	int				nmol = 0;
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		nw = 0;
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( atom->sel ) {
					memset(atomtype, ' ', 5);
					for ( i=(long)strlen(atom->type)-1, j=4; i>=0 && j>=0; i-- )
						if ( !isspace(atom->type[i]) ) atomtype[j--] = atom->type[i];
					fgro << setw(5) << res->num << setw(5) << res->type << 
						setw(5) << atomtype << setw(5) << atom->num <<
						fixed << setprecision(3) << setw(8) << 0.1*atom->coord[0] << 
						setw(8) << 0.1*atom->coord[1] << setw(8) << 0.1*atom->coord[2] <<
						setprecision(4) << setw(8) << 0.1*atom->vel[0] << 
						setw(8) << 0.1*atom->vel[1] << setw(8) << 0.1*atom->vel[2] << endl;
					nw++;
				}
			}
		}
		if ( nw > 0 ) nmol++;
    }
    
    fgro << fixed << setprecision(2) << setw(9) << 0.1*molgroup->box[0] << 
		setw(9) << 0.1*molgroup->box[1] << setw(9) << 0.1*molgroup->box[2] << endl;
    
    fgro.close();
    
    return nmol;
}

