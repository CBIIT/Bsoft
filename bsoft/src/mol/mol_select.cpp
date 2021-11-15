/**
@file	mol_select.cpp
@brief	Library routines to select atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20211029
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "rwresprop.h"
#include "mol_select.h"
#include "mol_util.h"
#include "seq_util.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "Matrix.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Returns the number of atoms selected.
@param 	*molgroup 	molecule group structure.
@return long 		number of atoms selected.
**/
long		molgroup_atoms_selected(Bmolgroup* molgroup)
{
	long			n(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				if ( atom->sel ) ++n;

	if ( verbose )
		cout << "Atoms selected:                 " << n << endl << endl;
	
	return n;
}

/**
@brief 	Sets selection based on a specification.
@param 	*molgroup 	molecule group structure.
@param 	selstr 		selection specification.
@return long 		number of atoms selected.
**/
long		molgroup_select(Bmolgroup* molgroup, Bstring selstr)
{
	if ( selstr == "all" )
		return molgroup_select_all(molgroup);
	
	long			i, j, len(selstr.length());
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	Bstring			mol_id;
	Bstring			res_id;
	Bstring			atom_id;

	molgroup_deselect_all(molgroup);

	for ( i=0; i<len; ++i ) {
		if ( selstr[i] == '#' ) {
			j = ++i;
			while ( i<len && selstr[i] != ':' && selstr[i] != '@' ) ++i;
			mol_id = selstr.substr(j, i-1);
		}
		if ( selstr[i] == ':' ) {
			j = ++i;
			while ( i<len && selstr[i] != '@' ) ++i;
			res_id = selstr.substr(j, i-1);
		}
		if ( selstr[i] == '@' ) {
			atom_id = selstr.substr(++i, len-1);
		}
	}
		
	if ( verbose ) {
		cout << "Selection identifiers:" << endl;
		if ( mol_id.length() )
			cout << "Model:                          " << mol_id << endl;
		if ( res_id.length() )
			cout << "Residue:                        " << res_id << endl;
		if ( atom_id.length() )
			cout << "Atom:                           " << atom_id << endl;
//		cout << endl;
	}

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		if ( mol_id.length() ) {
			if ( mol->id == mol_id ) mol->sel = 1;
			else mol->sel = 0;
		}
		for( res = mol->res; res; res = res->next ) {
			if ( res_id.length() ) {
//				cout << "-" << res->type << "-" << endl;
				if ( strstr(res_id.c_str(), res->type) ) i = 1;
				else i = 0;
			}
			if ( i ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					atom->sel = 1;
					if ( atom_id.length() ) {
						if ( strstr(atom->type, atom_id.c_str()) ) atom->sel = 1;
						else atom->sel = 0;
					}
				}
			}
		}
	}

	return molgroup_atoms_selected(molgroup);
}

/**
@brief 	Sets selection for all atoms.
@param 	*molgroup 	molecule group structure.
@return long 		number of atoms.
**/
long		molgroup_select_all(Bmolgroup* molgroup)
{
	long			n(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mol->sel = 1;
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next, n++ )
				atom->sel = 1;
	}
	
	return n;
}

/**
@brief 	Unsets selection for all atoms.
@param 	*molgroup 	molecule group structure.
@return long 		number of atoms.
**/
long		molgroup_deselect_all(Bmolgroup* molgroup)
{
	long			n(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next, n++ )
				atom->sel = 0;
	
	return n;
}

/**
@brief	Selects molecules.
@param 	*molgroup 	molecule group structure.
@param 	chains		comma-separated list of molecule ids.
@return int 		chains selected.
**/
int			molgroup_select_chains(Bmolgroup* molgroup, Bstring chains)
{
	int				nsel(0);
	Bmolecule*		mol;
	
	if ( verbose )
		cout << "Selecting chains: " << chains << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		if ( chains.contains(mol->id) ) {
			mol->sel = 1;
			nsel++;
		} else {
			mol->sel = 0;
		}
	}
	
	return nsel;
}

/**
@brief	Selects atoms in a ring.
@param 	*molgroup 	molecule group structure.
@param 	rmin		minimum radius.
@param 	rmax		maximum radius.
@return int 		0.
**/
int 		molgroup_coor_select_ring(Bmolgroup* molgroup, double rmin, double rmax)
{
    if ( ( rmin < 0 ) && ( rmax > 1000 ) ) return 0;
	
	int 		sel;
	double		d;
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting ring:                 " << rmin << " " << rmax << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			sel = 0;
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->sel = 0;
				d = sqrt(atom->coord[0]*atom->coord[0]+atom->coord[1]*atom->coord[1]);
				if ( ( d > rmin ) && ( d < rmax ) ) sel = 1;
			}
			if ( sel )
				for ( atom = res->atom; atom; atom = atom->next )
					atom->sel = 1;
			
		}
    }
	
	return 0;
}

/**
@brief	Selects atoms within a box.
@param 	*molgroup 	molecule group structure.
@param 	min			three-valued vector of minima.
@param 	max			three-valued vector of maxima.
@return int 		0.
**/
int 		molgroup_coor_select(Bmolgroup* molgroup, Vector3<double> min, Vector3<double> max)
{
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Selecting atoms:                " << min << " - " << max << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->sel = 1;
				if ( 	atom->coord[0] < min[0] || atom->coord[0] > max[0] ||
						atom->coord[1] < min[1] || atom->coord[1] > max[1] ||
						atom->coord[2] < min[2] || atom->coord[2] > max[2]	)
					atom->sel = 0;
			}
    	}
	}
	
	return 0;
}

/**
@brief	Selects a number of atoms.
@param 	*molgroup 	molecule group structure.
@param 	number		number of atoms to select.
@return int 			0.
**/
int 		molgroup_select(Bmolgroup* molgroup, long number)
{
	molgroup_deselect_all(molgroup);
	
	random_seed();

	vector<Batom*>	atarr = atom_get_array(molgroup);
	long		natom(atarr.size()), nsel(0);
	double		irm = natom*1.0L/get_rand_max();

	if ( verbose & VERB_PROCESS )
		cout << "Randomly selecting " << number << " atoms" << endl;
	
	long		j;
	while ( nsel < number ) {
		j = irm*random();
		if ( j < natom && atarr[j]->sel == 0 ) {
			atarr[j]->sel = 1;
			nsel++;
		}
	}
		
	return 0;
}


