/**
@file	rwatomprop.cpp
@brief	Library routines to read and write atom properties
@author Bernard Heymann
@date	Created: 19991114
@date	Modified: 20210214
**/

#include "rwatomprop.h"
#include "star.h"
#include "mol_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
Batomtype*		read_atom_prop_star(Bstring& filename);
int 			write_atom_prop_star(Bstring& filename, Batomtype* at_first);

/* Atomic parameters */
static struct atom_type  {   
	char	name[8];		// Atom name/symbol
	float	mass;			// Atom mass
	float	oxid;			// Oxidation state
	float	bond;			// Bond length
	float	vdw;			// VdW radius
	float	red,grn,blu;	// RGB colour
	float	scat[9];		// Atom scattering factor parameters (a1,a2,a3,a4,b1,b2,b3,b4,c)
} atom_data[] = {
    {"H",  1.007,  0.0, 1.0, 1.100, 1.0,1.0,1.0, {0.307,0.202,0.0,0.0,22.97,3.98,0.0,0.0,0.0}},
    {"C", 12.0110, 0.0, 1.5, 1.548, 0.7,0.7,0.7, {1.494,0.937,0.0,0.0,23.22,3.79,0.0,0.0,0.0}},
    {"N", 14.0067, 0.0, 1.5, 1.400, 0.0,0.0,1.0, {1.263,0.872,0.0,0.0,17.65,2.94,0.0,0.0,0.0}},
    {"O", 15.9994, 0.0, 1.5, 1.348, 1.0,0.0,0.0, {1.085,0.822,0.0,0.0,14.17,2.40,0.0,0.0,0.0}},
    {"P", 30.9738, 0.0, 1.5, 1.880, 1.0,0.0,1.0, {3.590,1.696,0.0,0.0,28.40,3.65,0.0,0.0,0.0}},
    {"S", 32.0600, 0.0, 1.5, 1.808, 1.0,1.0,0.0, {3.319,1.649,0.0,0.0,23.06,3.22,0.0,0.0,0.0}},
    {" ",       0,   0,   0,     0,   0,  0,  0, {0.0  ,0.0  ,0.0,0.0, 0.0, 0.0, 0.0,0.0,0.0}}
} ;

/**
@brief 	Reading atomic properties from parameter files.
@param 	&filename	file name (if empty, use a default file).
@return Batomtype*	atom property structure, NULL on failure.
**/
Batomtype*		get_atom_properties(Bstring& filename)
{
	if ( verbose & VERB_DEBUG )	
		cout << "DEBUG get_atom_properties: Initializing atomic parameters from " << filename << endl;
	
	Batomtype*		at_curr = NULL;
	Batomtype*		at_first = NULL;
		
	int 			i, j;
		
	// Atom parameter file
	Bstring			atfile;
	Bstring			propfile;
	if ( filename.length() ) atfile = filename;
	else atfile = "atom_prop.star";
	
	if ( access(atfile.c_str(), R_OK) == 0 )
		propfile = atfile;
	else
		propfile = parameter_file_path(atfile);

	filename = propfile;
	
	Bstring			ext = atfile.extension();
	if ( ext.length() ) {
		if ( ext.contains("star") )
			at_first = read_atom_prop_star(propfile);
//		else
//			at_first = read_atom_properties(propfile);
	}
	
	if ( !at_first ) {
		if ( verbose )
			cout << "Atom property file " << atfile << " not opened! Using default properties" << endl;
		for ( i=0; atom_data[i].mass && atom_data[i].vdw; i++ ) {
			at_curr = atom_type_add(&at_curr, atom_data[i].name);
			if ( i == 0 ) at_first = at_curr;
			at_curr->mass = atom_data[i].mass;
			at_curr->oxid = atom_data[i].oxid;
			at_curr->bond = atom_data[i].bond;
			at_curr->vdw = atom_data[i].vdw;
			at_curr->color[0] = (char) (255*atom_data[i].red);
			at_curr->color[1] = (char) (255*atom_data[i].grn);
			at_curr->color[2] = (char) (255*atom_data[i].blu);
			for ( j=0; j<4; j++ ) {
				at_curr->sfa[j] = atom_data[i].scat[j];
				at_curr->sfb[j] = atom_data[i].scat[j+4];
			}
			at_curr->sfc = atom_data[i].scat[8];
		}
	}
	
	for ( i=0, at_curr = at_first; at_curr; at_curr = at_curr->next, i++ ) ;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_atom_properties: Number of atom types = " << i << endl;
		
	return at_first;
}

/**
@brief 	Writing atomic properties from parameter files.
@param 	&filename	file name.
@param	*at			the atom property structure.
@return int			error code (<0 means failure).
**/
int 			write_atom_properties(Bstring& filename, Batomtype* at)
{
	int			err(0);
	Bstring		ext = filename.extension();

	if ( ext.contains("star") )
		err = write_atom_prop_star(filename, at);
	else {
		cerr <<  "Error: Extension \"" << ext << "\" not valid for parameter files!" << endl;
		cerr <<  "	Parameter file \"" << filename << "\" not written." << endl;
		err = -1;
	}
	
	if ( err < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	return err;
}

/**
@brief 	Adds an atom type structure to a linked list.
@param	**atype		pointer to any atom type in the list.
@param	*aname		atom type name.
@return Batomtype* 	new atom type.

	The function allocates memory for a new atom type structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Batomtype*		atom_type_add(Batomtype** atype, const char* aname)
{
	Batomtype* 		this_atype = *atype;
	Batomtype* 		new_atype = new Batomtype;
	memset(new_atype, 0, sizeof(Batomtype));
	
	if ( aname ) strncpy(new_atype->name, aname, 8);
	if ( aname ) {
		while ( isspace(aname[0]) ) aname++;
		strncpy(new_atype->name, aname, 6);
		if ( isalpha(aname[0]) ) {
			new_atype->el[0] = toupper(aname[0]);
			if ( isalpha(aname[1]) ) {
				if ( islower(aname[1]) ) new_atype->el[1] = aname[1];
			}
		}
	}

	if ( !this_atype )
		*atype = new_atype;
	else {
		while ( this_atype->next ) this_atype = this_atype->next;
		this_atype->next = new_atype;
	}
	
	return new_atype;
}

/**
@brief 	Destroys an atom type.
@param 	atype			the atom type.
@return int 			0.
**/
int			atom_type_kill(Batomtype* atype)
{
	if ( !atype ) return 0;
	
	Batomtype*	atype2;
	
	for ( ; atype; ) {
		atype2 = atype->next;
		delete atype;
		atype = atype2;
	}
	
	return 0;
}


Batomtype* 	read_atom_prop_star(Bstring& filename)
{
 	Bstar			star;
	
 	if ( star.read(filename.c_str()) < 0 ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		return NULL;
	}

	BstarBlock		block = star.block(0);
	BstarLoop		loop = block.loop(0);
	
	// Check if it contains atom parameters
	if ( loop.find(ATOM_TYPE_SYMBOL) < 0 ) {
		cerr <<  "Error: Star to atom properties conversion unsuccessful!" << endl;
		return NULL;
	}
	
	Batomtype*		at = NULL;
	Batomtype*		at_first = NULL;
	Bstring			symbol;
	long			j, k;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_atom_prop_star: Reading " << filename << endl;
	
	// Convert the atomic data
	for ( auto il = block.loops().begin(); il != block.loops().end(); ++il ) {
		if ( ( j = il->find(ATOM_TYPE_SYMBOL) ) >= 0 ) {
			for ( auto ir = il->data().begin(); ir != il->data().end(); ++ir ) {
				symbol = (*ir)[j];
				if ( at ) at = atom_type_add(&at, symbol.c_str());
				else at = atom_type_add(&at_first, symbol.c_str());
				if ( ( k = il->find(ATOM_TYPE_NUMBER) ) >= 0 ) at->z = to_integer((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_MASS) ) >= 0 ) at->mass = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_OXIDATION) ) >= 0 ) at->oxid = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_RADIUS_BOND) ) >= 0 ) at->bond = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_RADIUS_VDW) ) >= 0 ) at->vdw = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_A1) ) >= 0 ) at->sfa[0] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_A2) ) >= 0 ) at->sfa[1] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_A3) ) >= 0 ) at->sfa[2] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_A4) ) >= 0 ) at->sfa[3] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_A5) ) >= 0 ) at->sfa[4] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_B1) ) >= 0 ) at->sfb[0] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_B2) ) >= 0 ) at->sfb[1] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_B3) ) >= 0 ) at->sfb[2] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_B4) ) >= 0 ) at->sfb[3] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_B5) ) >= 0 ) at->sfb[4] = to_real((*ir)[k]);
				if ( ( k = il->find(ATOM_TYPE_SCAT_C) ) >= 0 ) at->sfc = to_real((*ir)[k]);
			}
		}
	}
		
	return at_first;
}

map<string, int>	atom_type_tags()
{
	map<string, int>	tag;
	tag[ATOM_TYPE_SYMBOL] = 0;
	tag[ATOM_TYPE_NUMBER] = 1;
	tag[ATOM_TYPE_MASS] = 2;
	tag[ATOM_TYPE_OXIDATION] = 3;
	tag[ATOM_TYPE_RADIUS_BOND] = 4;
	tag[ATOM_TYPE_RADIUS_VDW] = 5;
	tag[ATOM_TYPE_SCAT_A1] = 6;
	tag[ATOM_TYPE_SCAT_A2] = 7;
	tag[ATOM_TYPE_SCAT_A3] = 8;
	tag[ATOM_TYPE_SCAT_A4] = 9;
	tag[ATOM_TYPE_SCAT_A5] = 10;
	tag[ATOM_TYPE_SCAT_B1] = 11;
	tag[ATOM_TYPE_SCAT_B2] = 12;
	tag[ATOM_TYPE_SCAT_B3] = 13;
	tag[ATOM_TYPE_SCAT_B4] = 14;
	tag[ATOM_TYPE_SCAT_B5] = 15;
	tag[ATOM_TYPE_SCAT_C] = 16;

	return tag;
}

int 			write_atom_prop_star(Bstring& filename, Batomtype* at_first)
{
	if ( !at_first ) {
		cerr <<  "Error: No atom parameters to write!" << endl;
		return -1;
	}

 	Bstar			star;
	BstarBlock&		block = star.add_block("Atomic_parameters");
	BstarLoop&		loop = block.add_loop();

	loop.tags() = atom_type_tags();

	Batomtype*		at;
	
	star.line_length(120);
	
	// Count the number of atom types
	int				n, err(0);
	for ( n=0, at = at_first; at; at = at->next, n++ ) {
		vector<string>&	vs = loop.add_row(6);
		vs[0] = at->name;
		vs[1] = to_string(at->z);
		vs[2] = to_string(at->mass);
		vs[3] = to_string(at->oxid);
		vs[4] = to_string(at->bond);
		vs[5] = to_string(at->vdw);
		vs[6] = to_string(at->sfa[0]);
		vs[7] = to_string(at->sfa[1]);
		vs[8] = to_string(at->sfa[2]);
		vs[9] = to_string(at->sfa[3]);
		vs[10] = to_string(at->sfa[4]);
		vs[11] = to_string(at->sfb[0]);
		vs[12] = to_string(at->sfb[1]);
		vs[13] = to_string(at->sfb[2]);
		vs[14] = to_string(at->sfb[3]);
		vs[15] = to_string(at->sfb[4]);
		vs[16] = to_string(at->sfc);
	}

	err = star.write(filename.c_str());
	
	if ( err < 0 ) return err;
	
	return n;
}


