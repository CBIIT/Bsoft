/**
@file	rwmolecule.cpp
@brief	Library routines to read and write molecule file, including sequences and coordinates
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20220427
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "seq_util.h"
#include "linked_list.h"
#include "utilities.h"

#include "rwmol_star.h"
#include "rwmol_text.h"
#include "rwClustal.h"
#include "rwEMBL.h"
#include "rwFASTA.h"
#include "rwGenBank.h"
#include "rwGROMACS.h"
#include "rwPDB.h"
#include "rwPhylip.h"
#include "rwPIR.h"
#include "rwWAH.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
int 		molgroup_sequence_check(Bmolgroup* molgroup);

/**
@brief 	Initializes and allocates a new molecule group.
@return Bmolgroup*		the new molecule group, NULL if initialization failed.

	The selection string is set to "all".
	The spacegroup is set to 1, the space group string to "P 1".
	The point group is set to "C1".
	The unit cell is set to 1,1,1,90,90,90.

**/
Bmolgroup*  molgroup_init()
{
    Bmolgroup*	molgroup = new Bmolgroup;
	memset(molgroup, 0, sizeof(Bmolgroup));
	
	if ( !molgroup ) return NULL;

	molgroup->select = "all";
	molgroup->sel = 1;
	
	molgroup->spacegroup = 1;
	molgroup->unitcell = UnitCell(1, 1, 1, M_PI_2, M_PI_2, M_PI_2);
	
	molgroup->sgstring = "P 1";
	molgroup->pointgroup = "C1";
	
	return molgroup;
}

/**
@brief 	Adds a molecule to a linked list.
@param 	**mol		pointer to any molecule in the list.
@param 	*name			molecule name.
@return Bmolecule* 			new molecule.

	The function allocates memory for a new molecule structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bmolecule*	molecule_add(Bmolecule** mol, char* name)
{
	Bstring		thename(name);
	return molecule_add(mol, thename);
}

Bmolecule*	molecule_add(Bmolecule** mol, Bstring& name)
{
	Bmolecule* 		this_mol = *mol;
	Bmolecule* 		nu_mol = new Bmolecule;
	memset(nu_mol, 0, sizeof(Bmolecule));
	
	if ( name.length() ) nu_mol->id = name;
	nu_mol->sel = 1;
	
	if ( !this_mol )
		*mol = nu_mol;
	else {
		while ( this_mol->next ) this_mol = this_mol->next;
		this_mol->next = nu_mol;
	}
	
	return nu_mol;
}

/**
@brief 	Adds a residue to a linked list.
@param 	**res		pointer to any residue in the list.
@param	*type	residue type.
@return Bresidue* 			new residue.

	The function allocates memory for a new residue structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bresidue*	residue_add(Bresidue** res, const char* type)
{
	Bstring		thetype(type);
	return residue_add(res, thetype);
}

Bresidue*	residue_add(Bresidue** res, Bstring& type)
{
	Bresidue* 		this_res = *res;
	Bresidue* 		nu_res = new Bresidue;
	memset(nu_res, 0, sizeof(Bresidue));
	
	if ( type.length() ) {
		if ( type.contains("***") ) strcpy(nu_res->type, "UNK");
		else strncpy(nu_res->type, type.c_str(), 6);
	}
	
	if ( !this_res )
		*res = nu_res;
	else {
		while ( this_res->next ) this_res = this_res->next;
		this_res->next = nu_res;
	}
	
	return nu_res;
}

/**
@brief 	Adds an atom to a linked list.
@param 	**atom		pointer to any atom in the list.
@param	*type	atom type.
@return Batom* 				new atom.

	The function allocates memory for a new atom structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Batom*		atom_add(Batom** atom, const char* type)
{
	Bstring		thetype(type);
	return atom_add(atom, thetype);
}

Batom*		atom_add(Batom** atom, Bstring& type)
{
	Batom* 		this_atom = *atom;
	Batom* 		nu_atom = new Batom;
	memset(nu_atom, 0, sizeof(Batom));
	
	atom_clean_type(nu_atom, type.c_str());

	nu_atom->q = 1;
	nu_atom->mass = 1;
	nu_atom->sel = 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG atom_add: type=" << nu_atom->type << " el=" << nu_atom->el << endl;
	
	if ( !this_atom )
		*atom = nu_atom;
	else {
		while ( this_atom->next ) this_atom = this_atom->next;
		this_atom->next = nu_atom;
	}
	
	return nu_atom;
}

Batom*		atom_copy(Batom* atom)
{
	Batom* 		nu_atom = new Batom;
	
	nu_atom->next = NULL;	 		// Next atom in linked list
    nu_atom->num = atom->num;    	// Atom number
    memcpy(nu_atom->el, atom->el, 4);		// Element
    memcpy(nu_atom->type, atom->type, 8);	// Atom type
	nu_atom->tnum = atom->tnum;		// Atom type number (for reference purposes)
    nu_atom->coord = atom->coord;	// Coordinates
    nu_atom->q = atom->q;			// Occupancy
    nu_atom->b = atom->b;    		// Temperature factor
    nu_atom->mass = atom->mass;		// Atom mass
    nu_atom->chrg = atom->chrg;		// Charge on atom
    nu_atom->sel = atom->sel;    	// Select flag - extra
	nu_atom->vel = atom->vel;		// Velocities (for molecular dynamics)
	nu_atom->F = atom->F;			// Force vector (for molecular dynamics)

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG atom_copy: type=" << nu_atom->type << " el=" << nu_atom->el << endl;
	
	return nu_atom;
}

/**
@brief 	Counts the number of residues in a molecule group.
@param 	*molgroup			the molecule group.
@return long				number of residues.
**/
long		residue_count(Bmolgroup* molgroup)
{
	long			nres(0);
	Bmolecule*		mol;
	Bresidue*		res;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next ) nres++;
	
	return nres;
}

/**
@brief 	Counts the number of atoms in a molecule group.
@param 	*molgroup			the molecule group.
@return long				number of atoms.
**/
long		atom_count(Bmolgroup* molgroup)
{
	long			natom(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) natom++;
	
	return natom;
}

/**
@brief 	Cleans up the type string and assigns an element code to an atom.
@param 	*atom			atom.
@param	*type	atom type.
@return int					0.

	The first two alphanumeric characters of the type string is used to
	determine the element.

**/
int			atom_clean_type(Batom* atom, const char* type)
{
	if ( !type ) return 0;
	
	while ( isspace(type[0]) ) type++;
	if ( type[0] == '"' ) type++;
	
	strncpy(atom->type, type, 8);
	
	if ( atom->el[0] ) return 0;
	
	if ( isdigit(type[0]) ) {
		if ( type[1] == 'H' ) atom->el[0] = 'H';
		else atom->el[0] = 'C';
	} else if ( isalpha(type[0]) ) {
		atom->el[0] = toupper(type[0]);
		if ( isalpha(type[1]) ) {
			if ( islower(type[1]) ) atom->el[1] = type[1];
			// Special cases
			if ( atom->el[0] == 'M' )
				atom->el[1] = tolower(type[1]);
		}
	} else {
			cout << "Warning: No element assigned for atom " << atom->type << endl;
	}
	
	return 0;
}

/**
@brief 	Adds a bond to a linked list.
@param 	**bond		pointer to any bond in the list.
@param 	*atom1		atom1 of bond.
@param 	*atom2		atom2 of bond.
@param 	l				reference bond length.
@param 	k				bond strength.
@return Bbond* 				new bond.

	The function allocates memory for a new bond structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bbond*		bond_add(Bbond** bond, Batom* atom1, Batom* atom2, double l, double k)
{
	if ( l < 0.01 ) l = atom1->coord.distance(atom2->coord);
	if ( k < 1e-10 ) k = 1;
	
	Bbond* 			this_bond = *bond;
	Bbond* 			nu_bond = new Bbond;
	memset(nu_bond, 0, sizeof(Bbond));
	
	nu_bond->atom1 = atom1;
	nu_bond->atom2 = atom2;
	nu_bond->l = l;
	nu_bond->k = k;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG bond_add: strength=" << k << endl;
	
	if ( !this_bond )
		*bond = nu_bond;
	else {
		while ( this_bond->next ) this_bond = this_bond->next;
		this_bond->next = nu_bond;
	}
	
	return nu_bond;
}

/**
@brief 	Adds an angle to a linked list.
@param 	**angle		pointer to any angle in the list.
@param 	*atom1		atom1 of angle.
@param 	*atom2		atom2 of angle (central atom).
@param 	*atom3		atom3 of angle.
@param 	a				reference angle.
@param 	k				angle strength.
@return Bangle* 			new angle.

	The function allocates memory for a new angle structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bangle*		angle_add(Bangle** angle, Batom* atom1, Batom* atom2, Batom* atom3, double a, double k)
{
	Bangle* 		this_angle = *angle;
	Bangle* 		nu_angle = new Bangle;
	memset(nu_angle, 0, sizeof(Bangle));
	
	Vector3<double>	d1 = atom1->coord - atom2->coord;
	Vector3<double>	d2 = atom3->coord - atom2->coord;
	
	if ( a < 0.1 ) a = d1.angle(d2);
				
	nu_angle->atom1 = atom1;
	nu_angle->atom2 = atom2;
	nu_angle->atom3 = atom3;
	nu_angle->a = a;
	nu_angle->k = k;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG angle_add: size=" << a << " strength=" << k << endl;
	
	if ( !this_angle )
		*angle = nu_angle;
	else {
		while ( this_angle->next ) this_angle = this_angle->next;
		this_angle->next = nu_angle;
	}
	
	return nu_angle;
}

/**
@brief 	Destroys a molecule group linked list.
@param 	*molgroup		the molecule group linked list.
@return int 			0.
**/
int 		molgroup_list_kill(Bmolgroup* molgroup)
{
	if ( molgroup == NULL ) return 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG molgroup_list_kill: comment length = " << 
			molgroup->comment.length() << endl;
		cout << "DEBUG molgroup_list_kill: maxlen = " << molgroup->maxlen << 
			" (" << molgroup->maxlen*sizeof(int) << ")" << endl;
	}
	
	Bmolgroup*	molgroup2;
	
	for ( ; molgroup; ) {
		molgroup2 = molgroup->next;
		molgroup_kill(molgroup);
		molgroup = molgroup2;
	}

	return 0;
}

/**
@brief 	Destroys a molecule group.
@param 	*molgroup		the molecule group.
@return int 			0.
**/
int 		molgroup_kill(Bmolgroup* molgroup)
{
	if ( molgroup == NULL ) return 0;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG molgroup_kill: name = " << molgroup->id << endl;
		cout << "DEBUG molgroup_kill: comment length = " << 
			molgroup->comment.length() << endl;
		cout << "DEBUG molgroup_kill: maxlen = " << molgroup->maxlen << 
			" (" << molgroup->maxlen*sizeof(int) << ")" << endl;
	}
	
	Bmolecule	*mol, *mol2;
	
	for ( mol = molgroup->mol; mol; ) {
		mol2 = mol->next;//cout << mol << " - " << mol2 << endl;
		molecule_kill(mol);
		mol = mol2;
	}
	
	if ( molgroup->seqflag )
		delete[] molgroup->seqflag;
	
	bond_kill(molgroup->bond);
	angle_kill(molgroup->angle);
	
	molgroup->id = 0;
	molgroup->filename = 0;
	molgroup->comment = 0;
	molgroup->pointgroup = 0;
	molgroup->sgstring = 0;
	molgroup->select = 0;
	
	delete molgroup;

	return 0;
}

/**
@brief 	Destroys a molecule.
@param 	*mol		the molecule.
@return int 			0.
**/
int 		molecule_kill(Bmolecule* mol)
{
	if ( mol == NULL ) return 0;
	
	Bresidue	*res, *res2;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG molecule_kill: name = " << mol->id << endl;
		if ( mol->seq.length() ) cout << "DEBUG molecule_kill: seq length = " << 
			mol->seq.length() << endl;
		if ( mol->naseq.length() ) cout << "DEBUG molecule_kill: naseq length = " << 
			mol->naseq.length() << endl;
	}
	
	for ( res = mol->res; res; ) {
		res2 = res->next;
		residue_kill(res);
		res = res2;
	}
	
	if ( mol->sec ) kill_list((char *) mol->sec, sizeof(Bsecondary));
	
	mol->id = 0;
	mol->seq = 0;
	mol->naseq = 0;
	
	delete mol;
	
	return 0;
}
	
/**
@brief 	Destroys a residue.
@param 	*res		the residue.
@return int 			0.
**/
int 		residue_kill(Bresidue* res)
{
	if ( res == NULL ) return 0;
	
//	kill_list((char *) res->atom, sizeof(Batom));
	
	Batom		*atom, *atom2;
	
	for ( atom = res->atom; atom; ) {
		atom2 = atom->next;
		delete atom;
		atom = atom2;
	}
	
	delete res;
	
	return 0;
}
	
/**
@brief 	Dealocates a list of bonds.
@param 	*bond			first bond in the list.
@return int					0.

	All bonds downstream are deallocated.

**/
int			bond_kill(Bbond* bond)
{
	Bbond*		bond2;
	
	while ( bond ) {
		bond2 = bond->next;
		delete bond;
		bond = bond2;
	}
	
	return 0;
}

/**
@brief 	Dealocates a list of angles.
@param 	*angle		first angle in the list.
@return int					0.

	All angles downstream are deallocated.

**/
int			angle_kill(Bangle* angle)
{
	Bangle*		angle2;
	
	while ( angle ) {
		angle2 = angle->next;
		delete angle;
		angle = angle2;
	}
	
	return 0;
}

/**
@brief 	Copies a molecule group list.
@param 	*molgroup			the molecule group list.
@return Bmolgroup*			new molecule group list.

	All molecule groups are copied to a completely new list.

**/
Bmolgroup*	molgroup_list_copy(Bmolgroup* molgroup)
{
	Bmolgroup* 		molgroupcopy = NULL;
	Bmolgroup* 		mgc = NULL;
	Bmolgroup* 		mgc2 = NULL;
	
	for ( ; molgroup; molgroup=molgroup->next ) {
		mgc = molgroup_copy(molgroup);
		if ( !molgroupcopy ) molgroupcopy = mgc;
		else mgc2->next = mgc;
		mgc2 = mgc;
	}
	
	return molgroupcopy;
}

/**
@brief 	Copies a molecule group.
@param 	*molgroup			the molecule group.
@return Bmolgroup*			new molecule group.

	All parts of a molecule group are copied to a completely new structure
	hierarchy, except sequence flag array.

**/
Bmolgroup*	molgroup_copy(Bmolgroup* molgroup)
{
	Bmolgroup* 		molgroupcopy = molgroup_init();
	
	Bmolecule*		mol;
	Bmolecule*		molcopy = NULL, *newmol;

	molgroupcopy->next = NULL;
	molgroupcopy->num = molgroup->num;
	molgroupcopy->id = molgroup->id;
	molgroupcopy->filename = molgroup->filename;
	molgroupcopy->comment = molgroup->comment;
	molgroupcopy->spacegroup = molgroup->spacegroup;
	molgroupcopy->pointgroup = molgroup->pointgroup;
	molgroupcopy->sgstring = molgroup->sgstring;
	molgroupcopy->unitcell = molgroup->unitcell;
	molgroupcopy->maxlen = molgroup->maxlen;
	molgroupcopy->maxnalen = molgroup->maxnalen;
	molgroupcopy->min = molgroup->min;
	molgroupcopy->max = molgroup->max;
	molgroupcopy->box = molgroup->box;
	molgroupcopy->select = molgroup->select;
	molgroupcopy->sel = molgroup->sel;
	molgroupcopy->fom = molgroup->fom;

	if ( molgroup->seqflag ) {
		molgroupcopy->seqflag = new int[molgroup->maxlen];
		memcpy(molgroupcopy->seqflag, molgroup->seqflag, molgroup->maxlen*sizeof(int));
	}
	molgroupcopy->mol = NULL;
	molgroupcopy->bond = NULL;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		newmol = molecule_copy(mol);
		if ( molcopy ) molcopy->next = newmol;
		else molgroupcopy->mol = newmol;
		molcopy = newmol;
	}

	molgroupcopy->bond = molgroup_bond_list_copy(molgroup, molgroupcopy);

	return molgroupcopy;
}

/**
@brief 	Copies a molecule.
@param 	*mol		the molecule to be copied.
@return Bmolecule* 		the new molecule, NULL if copy failed.

	Generates a new molecule with the same structure as the given molecule.

**/
Bmolecule* 	molecule_copy(Bmolecule* mol)
{
	if ( mol == NULL ) return NULL;
	
	Bmolecule* 		nu_mol = new Bmolecule;
	memset(nu_mol, 0, sizeof(Bmolecule));

	nu_mol->num = mol->num;				// Molecule number
	nu_mol->id = mol->id;					// Molecule identifier
	nu_mol->nres = mol->nres;				// Number of residues
	nu_mol->nbase = mol->nbase;			// Number of bases
	nu_mol->first_codon = mol->first_codon;	// First amino acid codon in NA sequence
	nu_mol->seq = mol->seq;				// Protein sequence
	nu_mol->naseq = mol->naseq;			// Nucleic acid sequence
	nu_mol->min = mol->min;				// Coordinate minima
	nu_mol->max = mol->max;				// Coordinate maxima
	nu_mol->Bmin = mol->Bmin;				// B-factor statistics
	nu_mol->Bmax = mol->Bmax;
	nu_mol->Bave = mol->Bave;
	nu_mol->Bstd = mol->Bstd;
	nu_mol->sel = mol->sel;				// Selection flag
	nu_mol->fom = mol->fom;				// Figure-of-merit
	
	Bsecondary	*sec, *nu_sec;
	Bresidue	*res, *nu_res;
	Batom		*atom, *nu_atom;
	
	for ( res = mol->res; res; res = res->next ) {
		nu_res = residue_add(&nu_mol->res, res->type);
		nu_res->num = res->num;					// Residue number
		strncpy(nu_res->insert, res->insert, 2);	// Insertion code
		strncpy(nu_res->type, res->type, 6);		// Residue type
		nu_atom = NULL;
		for ( atom = res->atom; atom; atom = atom->next ) {
//			nu_atom = atom_add(&nu_res->atom, atom->type);
//			memcpy(nu_atom, atom, sizeof(Batom));
//			nu_atom->next = NULL;
			if ( nu_res->atom ) nu_atom = nu_atom->next = atom_copy(atom);
			else nu_res->atom = nu_atom = atom_copy(atom);
		}
	}
	
	for ( sec = mol->sec; sec; sec = sec->next ) {
		nu_sec = (Bsecondary *) add_item((char **) &nu_mol->sec, sizeof(Bsecondary));
		memcpy(nu_sec, sec, sizeof(Bsecondary));
		nu_sec->next = NULL;
		for ( res=mol->res, nu_res=nu_mol->res; res; res=res->next, nu_res=nu_res->next ) {
			if ( sec->first == res ) nu_sec->first = nu_res;
			if ( sec->last == res ) nu_sec->last = nu_res;
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molecule_copy: nu_mol->id = " << nu_mol->id << endl;
	
	return nu_mol;
}

/**
@brief 	Copies a molecule and assign to a new pointer in the molgroup.
@param 	*molgroup		the molecule group.
@param 	*mol		the molecule to be copied.
@return Bmolecule* 		the new molecule, NULL if copy failed.

	Adds a new molecule to the molecule group identical to the given
	molecule and returns a pointer to the new molecule.

**/
Bmolecule* 	mol_copy_and_add_to_molgroup(Bmolgroup* molgroup, Bmolecule* mol)
{
	if ( molgroup == NULL ) return NULL;
	if ( mol == NULL ) return NULL;
	
	Bmolecule*	nu_mol = molecule_copy(mol);

//	char		molid = 'A';
//	char		notset(1);
	
	if ( molgroup->mol ) {
		for ( mol=molgroup->mol; mol->next; mol=mol->next ) ;
		mol->next = nu_mol;
/*		if ( nu_mol->id[0] != ' ' ) {
			molid--;
			while ( notset ) {
				notset = 0;
				molid++;
				for ( mol=molgroup->mol; mol; mol=mol->next )
					if ( mol->id[0] == molid ) notset = 1;
			}
			if ( molid < 'A' || molid > 'Z' ) molid = 'A';
			nu_mol->id[0] = molid;
		}*/
	} else
		molgroup->mol = nu_mol;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mol_copy_and_add_to_molgroup: nu_mol->id = " << nu_mol->id << endl;
	
	return nu_mol;
}

/**
@brief 	Copies a bond list.
@param 	*molgroup 			molecule group structure.
@param 	*molgroupcopy		molecule group structure to copy bonds to.
@return Bbond* 				new bond list.

	A copy of the molecule group bond list is generated and returned. 

**/
Bbond*		molgroup_bond_list_copy(Bmolgroup* molgroup, Bmolgroup* molgroupcopy)
{
	long			maxid;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	Bbond*			bond;
	Bbond*			bond_copy = NULL;
	Bangle*			angle;
	Bangle*			angle_copy = NULL;
	
	molgroupcopy->bond = NULL;
	molgroupcopy->angle = NULL;
	
	for ( maxid=0, mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				if ( maxid < atom->num ) maxid = atom->num;
	
	maxid++;
	
	Batom**			acl = new Batom *[maxid];

    for ( mol=molgroupcopy->mol; mol; mol=mol->next )
		for( res=mol->res; res; res=res->next )
			for ( atom=res->atom; atom; atom=atom->next )
				acl[atom->num] = atom;
	
	for ( bond = molgroup->bond; bond; bond = bond->next ) {
		bond_copy = bond_add(&bond_copy, acl[bond->atom1->num], 
				acl[bond->atom2->num], bond->l, bond->k);
		if ( !molgroupcopy->bond ) molgroupcopy->bond = bond_copy;
	}

	for ( angle = molgroup->angle; angle; angle = angle->next ) {
		angle_copy = angle_add(&angle_copy, acl[angle->atom1->num], 
				acl[angle->atom2->num], acl[angle->atom3->num], angle->a, angle->k);
		if ( !molgroupcopy->angle ) molgroupcopy->angle = angle_copy;
	}

	delete[] acl;
	
	return molgroupcopy->bond;
}

/**
@brief 	Converts a molecule group list to a single molecule group.
@param 	*molgroup 	molecule group list.
@return int					0.

	The input molecule group list is replace by a single molecule group. 

**/
int			molgroup_from_molgroup_list(Bmolgroup* molgroup)
{
	Bmolgroup*		mg = molgroup;
	Bmolecule*		mol = molgroup->mol;
	Bbond*			bond = molgroup->bond;
	Bangle*			angle = molgroup->angle;
	
	for ( mg = molgroup->next; mg; mg = mg->next ) {
		if ( mol ) {
			for ( ; mol->next; mol = mol->next ) ;
			mol->next = mg->mol;
		} else {
			mol = molgroup->mol = mg->mol;
		}
		mg->mol = NULL;
		if ( bond ) {
			for ( ; bond->next; bond = bond->next ) ;
			bond->next = mg->bond;
		} else {
			bond = molgroup->bond = mg->bond;
		}
		mg->bond = NULL;
		if ( angle ) {
			for ( ; angle->next; angle = angle->next ) ;
			angle->next = mg->angle;
		} else {
			angle = molgroup->angle = mg->angle;
		}
		mg->angle = NULL;
	}
	molgroup_list_kill(molgroup->next);
	molgroup->next = NULL;
	molgroup_atom_renumber(molgroup, 1);
	
	return 0;
}

/**
@brief 	The generalized function for reading molecular files.
@param 	*filename	the file name.
@param 	*atom_select	a selection string.
@param	*paramfile	parameter file name.
@return Bmolgroup*				new molecule group, NULL if reading failed.

	All sequence and atomic coordinate information is read from a file into
	an internal hierarchy of structures in linked lists:
		Bmolgroup   molecule group or collection of molecules
		Bmolecule   linked list of molecules in the group
		Bresidue	linked list of residues in a molecule
		Batom		linked list of atoms in a residue
		Bbond		linked list of bonds in the molecule group
	The selection string is used to select for specific atom types:
		CA			C-alpha atoms only
	The parameter file is used to load atomic properties, such as mass
		and charge. The default file is bsoft/parameters/atom_prop.star.
	The input format is based on the file name extension.

**/
Bmolgroup*  read_molecule(const char *filename, const char *atom_select, const char* paramfile)
{
	Bstring		thefile(filename);
	Bstring		atomsel(atom_select);
	Bstring		theparam(paramfile);
	return read_molecule(thefile, atomsel, theparam);
}

Bmolgroup*  read_molecule(Bstring& filename, Bstring& atom_select, Bstring& paramfile)
{
	if ( filename.empty() ) {
		error_show("No molecule filename!", __FILE__, __LINE__);
		return NULL;
	}
	
	Bstring			ext = filename.extension();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_molecule: extension=" << ext << endl;
	
	Bstring			thefile = filename.pre(':');
	
	if ( thefile.length() < 1 ) {
		error_show("No molecule file given!", __FILE__, __LINE__);
		return NULL;
	}
	
	long			i, n(-1);
	
	Bmolgroup*		molgroup = molgroup_init();
	molgroup->filename = thefile;
	molgroup->select = atom_select.upper();
	
	if ( verbose & VERB_LABEL ) {
	    cout << "Reading file:                   " << molgroup->filename << endl;
    	cout << "Selecting:                      " << molgroup->select << endl << endl;
	}
	
//	cout << thefile << "\t" << sizeof(Bmolgroup) << endl;
	
	if ( ext.contains("star") || ext.contains("cif") )
    	n = read_mol_star(thefile, molgroup);
	else if ( ext.contains("txt") )
		n = read_mol_text(thefile, molgroup);
	else if ( ext.contains("aln") )
		n = readClustal(thefile, molgroup);
	else if ( ext.contains("embl") )
		n = readEMBL(thefile, molgroup);
	else if ( ext.contains("fasta") )
		n = readFASTA(thefile, molgroup);
	else if ( ext.contains("gb") || ext.contains("gp") || ext.contains("gen") )
		n = readGenBank(thefile, molgroup);
    else if ( ext.contains("gro") )
    	n = readGROMACS(thefile, molgroup);
    else if ( ext.contains("pdb") || ext.contains("ent") )
    	n = readPDB(thefile, molgroup);
	else if ( ext.contains("phy") )
		n = readPhylip(thefile, molgroup);
	else if ( ext.contains("pir") )
		n = readPIR(thefile, molgroup);
	else if ( ext.contains("wh") || ext.contains("wah") )
		n = readWAH(thefile, molgroup);
	else {
		cerr << "Error: File type not supported!" << endl;
		n = -1;
	}
	
	if ( n < 0 ) {
		error_show(thefile.c_str(), __FILE__, __LINE__);
		molgroup_kill(molgroup);
		return NULL;
	}
	
	if ( molgroup->id.length() < 1 ) molgroup->id = molgroup->mol->id;
	
	molgroup_sequence_check(molgroup);
	
	molecule_get_masses(molgroup, paramfile);
	
	molgroup_stats(molgroup);
    
	if ( molgroup->maxlen > 0 ) {
		molgroup->seqflag = new int[molgroup->maxlen];
		for ( i=0; i<molgroup->maxlen; i++ ) molgroup->seqflag[i] = 1;
	}
	
    return molgroup;
}

/**
@brief	Reads and catenates multiple molecule files.
@param 	file_list 			list of file names.
@param 	set_pbc 			flag to fit within periodic boundaries.
@param 	box 				periodic boundary box.
@param 	atom_select 		atomic selection.
@param 	paramfile 			atomic parameters.
@return Bmolgroup*			new molecule group.
**/
Bmolgroup*	read_molecule(Bstring* file_list, int set_pbc, Vector3<double> box, 
				Bstring atom_select, Bstring paramfile)
{
	Bstring*		filename;
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		molgroup1 = NULL;
	Bmolecule*		mol = NULL;
	Bbond*			bond;
	Bangle*			angle;
	
	for ( filename = file_list; filename; filename = filename->next ) {
		molgroup1 = read_molecule(*filename, atom_select, paramfile);
		if ( set_pbc ) {
			if ( box[0] > 0 && box[1] > 0 && box[2] > 0 ) molgroup1->box = box;
			molgroup_resolve_pbc(molgroup1);
		}
		if ( molgroup ) {
			if ( mol ) {
				for ( ; mol->next; mol = mol->next ) ;
				mol->next = molgroup1->mol;
			} else {
				mol = molgroup->mol = molgroup1->mol;
			}
			molgroup1->mol = NULL;
			if ( bond ) {
				for ( ; bond->next; bond = bond->next ) ;
				bond->next = molgroup1->bond;
			} else {
				bond = molgroup->bond = molgroup1->bond;
			}
			molgroup1->bond = NULL;
			if ( angle ) {
				for ( ; angle->next; angle = angle->next ) ;
				angle->next = molgroup1->angle;
			} else {
				angle = molgroup->angle = molgroup1->angle;
			}
			molgroup1->angle = NULL;
			molgroup_kill(molgroup1);
		} else {
			molgroup = molgroup1;
			mol = molgroup->mol;
			bond = molgroup->bond;
			angle = molgroup->angle;
		}
	}

	return molgroup;
}

/**
@brief 	Writes a molecule group.

	The output format is based on the file name extension.

@param 	*filename		the file name.
@param	*molgroup the molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int 		write_molecule(char *filename, Bmolgroup *molgroup)
{
	Bstring		thefile(filename);
	return write_molecule(thefile, molgroup);
}

int 		write_molecule(Bstring& filename, Bmolgroup *molgroup)
{
	if ( !filename.length() ) {
		error_show("No output file name given!", __FILE__, __LINE__);
		return -1;
	}

	if ( !molgroup->mol ) {
		cerr << "Error: No molecules to write!" << endl;
		return -1;
	}
	
	int				n;
	
	Bstring			ext = filename.extension();
	
	Bstring			thefile(filename.pre(':'));
	
    if ( verbose & VERB_LABEL )
		cout << "Writing file:                   " << thefile << endl;
    
    if ( ext.contains("star") || ext.contains("cif") )
    	n = write_mol_star(thefile, molgroup);
	else if ( ext.contains("txt") )
		n = write_mol_text(thefile, molgroup);
	else if ( ext.contains("aln") )
		n = writeClustal(thefile, molgroup);
	else if ( ext.contains("embl") )
		n = writeEMBL(thefile, molgroup);
	else if ( ext.contains("fasta") )
		n = writeFASTA(thefile, molgroup);
	else if ( ext.contains("gb") || ext.contains("gp") || ext.contains("gen") )
		n = writeGenBank(thefile, molgroup);
    else if ( ext.contains("gro") )
    	n = writeGROMACS(thefile, molgroup);
    else if ( ext.contains("pdb") || ext.contains("ent") )
    	n = writePDB(thefile, molgroup);
	else if ( ext.contains("phy") )
		n = writePhylip(thefile, molgroup);
	else if ( ext.contains("pir") )
		n = writePIR(thefile, molgroup);
	else if ( ext.contains("wh") || ext.contains("wah") )
		n = writeWAH(thefile, molgroup);
	else {
		cerr << "Error: File type not supported!" << endl;
		n = -1;
	}
    
	if ( n < 0 ) {
		error_show(thefile.c_str(), __FILE__, __LINE__);
		return n;
	}
	
	molgroup->filename = filename;
	
    if ( verbose & VERB_PROCESS )
		cout << "Molecules written:              " << n << endl;
	
	molgroup_stats(molgroup);
    
    return n;
}

/**
@brief 	Writes a molecule group list.

	The output files are numbered if the list constains more than one molecule group.

@param 	*filename		the file name.
@param	*molgroup the molecule group.
@return int 				number of molecules written (<0 if writing failed).
**/
int			molgroup_list_write(Bstring& filename, Bmolgroup* molgroup)
{
	int			i, n(0);
	Bstring		onename;
	
	if ( molgroup->next ) {
		for ( i=1; molgroup; molgroup = molgroup->next, i++ ) {
			onename = filename.pre_rev('.') + Bstring(i, "_%04d.") + filename.post_rev('.');
			n += write_molecule(onename, molgroup);
		}
	} else {
		n = write_molecule(filename, molgroup);
	}
	
	return n;
}

/**
@brief 	Counts the total number of molecules in a molecule group.
@param	*molgroup the molecule group.
@return long 				number of molecules.
**/
long 		molgroup_count_molecules(Bmolgroup* molgroup)
{
	long		n(0);
	Bmolecule*	mol;

	for ( mol=molgroup->mol; mol; mol=mol->next ) n++;

	return n;
}

/**
@brief 	Counts the total number of residues in a molecule group.
@param	*molgroup the molecule group.
@return long 				number of residues.
**/
long 		molgroup_count_residues(Bmolgroup* molgroup)
{
	long		n(0);
	Bmolecule*	mol;

	for ( mol=molgroup->mol; mol; mol=mol->next )
		n += mol_count_residues(mol);

	return n;
}

/**
@brief 	Counts the total number of residues in a molecule.
@param 	*mol		the molecule.
@return long 				number of residues.
**/
long 		mol_count_residues(Bmolecule *mol)
{
	long		n(0);
	Bresidue*	res;
	
	if ( mol->res )
		for ( res=mol->res; res; res=res->next ) n++;
	else if ( mol->seq.length() )
		n = mol->seq.length();

	return n;
}

/**
@brief 	Counts the total number of atoms in a molecule group.
@param	*molgroup the molecule group.
@return long 				number of atoms.
**/
long 		molgroup_count_atoms(Bmolgroup* molgroup)
{
	long		n(0);
	Bmolecule*	mol;

	for ( mol=molgroup->mol; mol; mol=mol->next )
		n += mol_count_atoms(mol);

	return n;
}

/**
@brief 	Counts the total number of atoms in a molecule.
@param 	*mol		the molecule.
@return long 				number of atoms.
**/
long 		mol_count_atoms(Bmolecule *mol)
{
	long		n(0);
	Bresidue*	res;
	Batom*  	atom;

	for ( res=mol->res; res; res=res->next )
		for ( atom=res->atom; atom; atom=atom->next ) n++;

	return n;
}

/*
@brief 	Checks the sequences in a molecule group, and generate them if needed.

	The sequence is first converted to upper case and then checked to see
	if it is a nucleic acid or protein sequence. The nucleic acid and
	protein sequences are referenced through different pointers.

@param	*molgroup the molecule group.
@return int 				number of molecules (<0 if writing failed).
**/
int 		molgroup_sequence_check(Bmolgroup* molgroup)
{
	if ( !molgroup ) return -1;
	
	Bmolecule*	mol;
	Bresidue*  	res;

	char*		temp;
	long 		i, j, n, nmol(0), gaps(1);
	
	if ( molgroup->select.contains("NOGAP") ) gaps = 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_sequence_check: select=" << molgroup->select << " gaps=" << gaps << endl;
	
    if ( verbose & VERB_PROCESS )
		cout << "Checking sequence" << endl << endl;
	
	molgroup->maxlen = 0;
	for ( mol = molgroup->mol; mol; mol = mol->next, nmol++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG molgroup_sequence_check: sequence " << mol->seq << endl;
		res = mol->res;
		if ( ( n = mol->naseq.length() ) ) {
			temp = new char[n+1];
			if ( gaps ) {
				for ( i=0, j=0; i<n; i++ ) {
					if ( isalpha(mol->naseq[i]) ) temp[j++] = toupper(mol->naseq[i]);
					if ( mol->naseq[i] == '-' || mol->naseq[i] == '.' ) temp[j++] = '-';
				}
			} else {
				for ( i=0, j=0; i<n; i++ ) {
					if ( isalpha(mol->naseq[i]) ) temp[j++] = toupper(mol->naseq[i]);
				}
			}
			temp[j] = 0;
			mol->naseq = temp;
			delete[] temp;
			mol->nbase = j;
		}
		if ( ( n = mol->seq.length() ) ) {
			temp = new char[n+1];
			if ( gaps ) {
				for ( i=0, j=0; i<n; i++ ) {
					if ( isalpha(mol->seq[i]) ) temp[j++] = toupper(mol->seq[i]);
					else if ( mol->seq[i] == '-' || mol->seq[i] == '.' ) temp[j++] = '-';
				}
			} else {
				for ( i=0, j=0; i<n; i++ ) {
					if ( isalpha(mol->seq[i]) ) temp[j++] = toupper(mol->seq[i]);
				}
			}
			temp[j] = 0;
			mol->seq = temp;
			delete[] temp;
			mol->nres = j;
		} else if ( res ) {
			for ( n = 0, res = mol->res; res; res = res->next, n++);
			mol->nres = n;
			temp = new char[n+1];
			for ( n = 0, res = mol->res; res; res = res->next, n++ )
				temp[n] = getcode1(res->type);
			temp[n] = 0;
			mol->seq = temp;
			delete[] temp;
			mol->nres = n;
		}
		if ( molgroup->maxlen < mol->seq.length() )
			molgroup->maxlen = mol->seq.length();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG molgroup_sequence_check: length=" << mol->nres << " sequence " << mol->seq << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_sequence_check: maxlen=" << molgroup->maxlen << endl;
		
	return nmol;
}

/**
@brief 	Removes redundant gaps from an alignment.

	All positions in an alignment with only gaps are removed.

@param	*molgroup the molecule group.
@return int 				0.
**/
int			molgroup_consolidate_gaps(Bmolgroup* molgroup)
{
	long		nmol = molgroup_count_molecules(molgroup);
	long		i, n, m;
	Bstring		newseq;
	Bmolecule*	mol;
	
	for ( m=i=0; i<molgroup->maxlen; i++ ) {
		for ( n=0, mol=molgroup->mol; mol; mol=mol->next ) {
			if ( i >= mol->seq.length() ) n++;
			else if ( mol->seq[i] == '-' ) n++;
		}
		if ( n >= nmol ) {
			molgroup->seqflag[i] = 0;
			m++;
		}
	}
	
	m = molgroup->maxlen - m;
	
	for ( mol=molgroup->mol; mol; mol=mol->next ) {
//		newseq = Bstring(m, '\0');
		newseq = Bstring('\0', m);
		for ( n=i=0; i<molgroup->maxlen; i++ ) {
			if ( molgroup->seqflag[i] )
				newseq[n++] = mol->seq[i];
		}
		mol->seq = newseq;
		mol->nres = m;
	}
	
	delete[] molgroup->seqflag;
	
	molgroup->maxlen = m;
	
	molgroup->seqflag = new int[m];
	for ( i=0; i<m; i++ ) molgroup->seqflag[i] = 1;
	
	return 0;
}

/**
@brief 	Calculates the statistics of a molecule group.
@param	*molgroup molecule group.
@param 	show			flag to show statistics.
@return long 				number of molecules (<0 if writing failed).
**/
long 		molgroup_stats(Bmolgroup* molgroup, int show)
{
	if ( !molgroup ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_stats: Calculating molecule group statistics" << endl;

	long 		nmol(0), natom_all(0), nres_all(0), nbase_all(0);
	Bmolecule*	mol;
	double		Bmin_all(1000), Bmax_all(-1000);
	
	molgroup->min[0] = molgroup->min[1] = molgroup->min[2] = 1000;
	molgroup->max[0] = molgroup->max[1] = molgroup->max[2] = -1000;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		natom_all += mol_stats(mol);
		nres_all += mol->nres;
		nbase_all += mol->nbase;
		if ( Bmin_all > mol->Bmin ) Bmin_all = mol->Bmin;
		if ( Bmax_all < mol->Bmax ) Bmax_all = mol->Bmax;
		molgroup->min = molgroup->min.min(mol->min);
		molgroup->max = molgroup->max.max(mol->max);
		nmol++;
	}
	
	long		nbond = count_list((char *) molgroup->bond);
	
	if ( molgroup->box.volume() < 0.1 )
		molgroup->box = molgroup->max - molgroup->min;
	
	if ( show ) {
		cout << "Statistics for molecule group:" << endl;
		cout << "Molecules:                      " << nmol << endl;
		cout << "Maximum sequence length:        " << molgroup->maxlen << endl;
		if ( nbase_all )
			cout << "Bases:                          " << nbase_all << endl;
		if ( nres_all )
			cout << "Residues:                       " << nres_all << endl;
		if ( natom_all )
			cout << "Atoms:                          " << natom_all << endl;
		if ( nbond )
			cout << "Bonds specfied:                 " << nbond << endl;
		if ( natom_all > 0 ) {
			cout << "Coordinate minima and maxima:   " << molgroup->min << " " << molgroup->max << endl;
			cout << "Bounding box size:              " << molgroup->box << endl;
			cout << "Bfactor min & max:              " << Bmin_all << " " << Bmax_all << endl;
		}
		cout << endl;
	}
	
	return nmol;
}

long 		molgroup_stats(Bmolgroup* molgroup)
{
	int			show(0);
	if ( verbose & VERB_PROCESS ) show = 1;
	return molgroup_stats(molgroup, show);
}

/**
@brief 	Calculates the statistics of a molecule.
@param 	*mol		the molecule.
@param 	show			flag to show statistics.
@return long				number of atoms (<0 if writing failed).
**/
long 		mol_stats(Bmolecule* mol, int show)
{
	if ( !mol ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mol_stats: Calculating molecule statistics: " << mol->id << endl;

	long			natom(0), nres(0), nhelix(0), nstrand(0);
	Bsecondary*		sec = NULL;
	Bresidue*		res;
	Batom*			atom;
	double			Bsum(0), Bssum(0);
	
	mol->Bave = mol->Bstd = 0;
	if ( mol->res ) {
		mol->min = Vector3<double>(1000, 1000, 1000);
		mol->max = Vector3<double>(-1000, -1000, -1000);
		mol->Bmin = 1000;
		mol->Bmax = -1000;

		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( atom->sel ) {
					mol->min = mol->min.min(atom->coord);
					mol->max = mol->max.max(atom->coord);
					if ( mol->Bmin > atom->b ) mol->Bmin = atom->b;
					if ( mol->Bmax < atom->b ) mol->Bmax = atom->b;
					Bsum += atom->b;
					Bssum += atom->b*atom->b;
					natom++;
				}
			}
			nres++;
		}
		if ( natom ) {
			mol->Bave = Bsum/natom;
			mol->Bstd = sqrt(Bssum/natom - mol->Bave*mol->Bave);
		}
	}
    
	if ( mol->max[0] - mol->min[0] < 1 ) {
		mol->max[0] = mol->min[0] = (mol->max[0] + mol->min[0])/2;
		mol->max[0] += 0.5;
		mol->min[0] -= 0.5;
	}
	if ( mol->max[1] - mol->min[1] < 1 ) {
		mol->max[1] = mol->min[1] = (mol->max[1] + mol->min[1])/2;
		mol->max[1] += 0.5;
		mol->min[1] -= 0.5;
	}
	if ( mol->max[2] - mol->min[2] < 1 ) {
		mol->max[2] = mol->min[2] = (mol->max[2] + mol->min[2])/2;
		mol->max[2] += 0.5;
		mol->min[2] -= 0.5;
	}
	
	if ( mol->sec ) {
		for ( sec=mol->sec; sec; sec=sec->next )
			if ( sec->type == Strand ) nstrand++;
			else nhelix++;
	}
		
	if ( show ) {
		cout << "Statistics for molecule:        " << mol->id << endl;
		if ( mol->nbase )
			cout << "Bases:                          " << mol->nbase << endl;
		if ( mol->nres )
			cout << "Sequence length:                " << mol->nres << endl;
		if ( nres )
			cout << "Residues:            	         " << nres << endl;
		if ( natom )
			cout << "Atoms:                          " << natom << endl;
		if ( nhelix )
			cout << "Helices:            	         " << nhelix << endl;
		if ( nstrand )
			cout << "Strands:            	         " << nstrand << endl;
		if ( natom > 0 ) {
			cout << "Coordinate minima and maxima:   " << mol->min << " " << mol->max << endl;
			cout << "Bfactor min, max, avg, & std:   " << 
				mol->Bmin << " " << mol->Bmax << " " << mol->Bave << " " << mol->Bstd << endl;
		}
		cout << endl;
	}
	
	return natom;
}

long 		mol_stats(Bmolecule* mol)
{
	int			show(0);
	if ( verbose & VERB_FULL ) show = 1;
	return mol_stats(mol, show);
}

/**
@brief 	Puts a set of strings and time in the main comment of a molecule group.

	This is designed to pack the command line into a string followed by
	a second string for the time.

@param 	*molgroup 	the molecule group.
@param 	n				the number of strings.
@param	**strings		an array of strings.
@return int 				string length of the new comment.
**/
int 		molecule_update_comment(Bmolgroup* molgroup, int n, char** strings)
{
	int			i; 
	Bstring		nu_comment("");
	Bmolgroup*	mg;
	
	for ( i=0; i<n; i++ )
		nu_comment = nu_comment + strings[i] + " ";
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molecule_update_comment: " << nu_comment << endl;
	 
	for ( mg = molgroup; mg; mg = mg->next ) 
		mg->comment = nu_comment;

	return molgroup->comment.length();
}

/**
@brief 	Gets atomic masses from a parameter file.
@param 	*molgroup 	the molecule group.
@param 	&paramfile	parameter file name.
@return int 				0.
**/
int 		molecule_get_masses(Bmolgroup* molgroup, Bstring& paramfile)
{
	Batomtype*		atompar = get_atom_properties(paramfile);
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Batomtype*  	at;
	
	int				t;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->tnum = -1;
				for ( t=0, at = atompar; at && atom->tnum < 0; at = at->next, t++ ) {
					if ( strncmp(atom->el, at->el, 2) == 0 ) {
//						cout << "atom=" << atom->el << ", type=" << at->el << "." << endl;
						atom->tnum = t;
						atom->mass = at->mass;
					}
				}
				if ( atom->tnum < 0 ) {		// Sometimes people use incorrect case
					for ( t=0, at = atompar; at && atom->tnum < 0; at = at->next, t++ ) {
						if ( strncasecmp(atom->el, at->el, 2) == 0 ) {
//							cout << "atom=" << atom->el << ", type=" << at->el << "." << endl;
							atom->tnum = t;
							atom->mass = at->mass;
						}
					}
				}
				if ( verbose ) {
					if ( atom->tnum < 0 )
						cerr << "Warning: No type asigned for atom " << atom->num 
							<< " (" << atom->type << ", " << atom->el << ")!" << endl;
					if ( atom->mass <= 0 )
						cerr << "Warning: No mass asigned for atom " << atom->num 
							<< " (" << atom->type << ", " << atom->el << ")!" << endl;
				}
			}
		}
	}
	
	atom_type_kill(atompar);
	
	return 0;
}

int			bond_exists(Bbond* bondlist, Batom* atom1, Batom* atom2)
{
	if ( !bondlist ) return 0;

	Bbond*			bond = NULL;
	
	for ( bond = bondlist; bond; bond = bond->next ) {
		if ( bond->atom1 == atom1 && bond->atom2 == atom2 ) return 1;
		if ( bond->atom1 == atom2 && bond->atom2 == atom1 ) return 1;
	}
	
	return 0;
}

/**
@brief 	Generates a bond list based on atom separation.
@param 	*molgroup 	the molecule group.
@param 	maxlength		maximum bond length.
@param 	wrap			wrap around periodic boundaries if !=0.
@return Bbond* 				new bond list.
**/
Bbond*		molgroup_bond_list_generate(Bmolgroup* molgroup, double maxlength, int wrap)
{	
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Bond list already defined!" << endl;
		return molgroup->bond;
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating bonds:" << endl;
		cout << "Bond length cutoff:             " << maxlength << " A" << endl;
	}
	
	long			natom(0);
	long			i, ii, x, y, z, xx, yy, zz, ix, iy, iz, xs, ys, xe, ye, ze;
	long			nbond(0);
	Vector3<double>	box = molgroup->box;
	Vector3<double>	d;
	double			dist, reflength;
	Vector3<double>	sampling(maxlength, maxlength, maxlength);
	Vector3<int>	size((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/size[i] + 0.001;
	long			boxsize = (long) size.volume();
	Bbond*			bond = NULL;
	Bbond*			bondlist = NULL;
	Batom			*atom, *atom2;
	Latom			*latom, *latom2;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG molgroup_bond_list_generate: box=" << box << endl;
		cout << "DEBUG molgroup_bond_list_generate: size=" << size << endl;
		cout << "DEBUG molgroup_bond_list_generate: sampling=" << sampling << endl;
	}
	
	Latom**			alist = molgroup_atom_mesh_lists(molgroup, size, sampling);
	
	for ( i=z=0; z<size[2]; z++ ) {
		for ( y=0; y<size[1]; y++ ) {
			for ( x=0; x<size[0]; x++,i++ ) {
				for ( latom = alist[i]; latom; latom = latom->next ) {
					natom++;
					ze = ( size[2] < 2 )? z: z+1;
					for ( zz=z; zz<=ze; zz++ ) {
						iz = zz;
						if ( iz < 0 ) iz += (long) size[2];
						if ( iz >= size[2]) iz -= (long) size[2];
						ys = ( z == zz )? y: y-1;
						ye = ( size[1] < 2 )? y: y+1;
//						for ( yy=y-1; yy<=y+1; yy++ ) {
						for ( yy=ys; yy<=ye; yy++ ) {
							iy = yy;
							if ( iy < 0 ) iy += (long) size[1];
							if ( iy >= size[1] ) iy -= (long) size[1];
							xs = ( z==zz && y==yy )? x: x-1;
							xe = ( size[0] < 2 )? x: x+1;
//							for ( xx=x-1; xx<=x+1; xx++ ) {
							for ( xx=xs; xx<=xe; xx++ ) {
								ix = xx;
								if ( ix < 0 ) ix += (long) size[0];
								if ( ix >= size[0] ) ix -= (long) size[0];
								ii = (iz*(long)size[1] + iy)*(long)size[0] + ix;
								if ( i == ii ) latom2 = latom->next;
								else latom2 = alist[ii];
								for ( ; latom2; latom2 = latom2->next ) {
									atom = latom->atom;
									atom2 = latom2->atom;
									if ( wrap )
										d = vector3_difference_PBC(atom->coord, atom2->coord, box);
									else
										d = atom->coord - atom2->coord;
									dist = d.length();
									if ( strncmp(atom->el, "H\0", 2) && strncmp(atom2->el, "H\0", 2) )
										reflength = maxlength;
									else
										reflength = 1.1;
									if ( dist > 0.1*reflength && dist <= reflength ) {
										if ( bond_exists(bondlist, atom, atom2) == 0 ) {
											bond = bond_add(&bond, atom, atom2, dist, 1);
											if ( !bondlist ) bondlist = bond;
											nbond++;
											if ( verbose & VERB_FULL )
												cout << "Bond " << nbond << ":\t" 
													<< atom->num << " - " << atom2->num << " = "
													<< dist << endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for ( i=0; i<boxsize; i++ ) kill_list((char *) alist[i], sizeof(Latom));
	delete[] alist;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of atoms:                " << natom << endl;
		cout << "Number of bonds generated:      " << nbond << endl << endl;
	}
	
	molgroup->bond = bondlist;

	return bondlist;
}

/**
@brief 	Generates an intramolecular distance-based bond list.

	This function defines bonds on distance and within molecules.
	If the molecule group already has a bond list, no new bonds are generated. 

@param 	*molgroup molecule group structure.
@param 	bondlength	maximum bond length.
@param 	wrap			wrap around periodic boundaries if !=0.
@return Bbond* 				new bond list.
**/
Bbond*		mol_bond_list_generate(Bmolgroup* molgroup, double bondlength, int wrap)
{
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Bond list already defined!" << endl;
		return molgroup->bond;
	}
	
	long 		natom(0), nbond(0);
	double		dist;
	Bbond*		bond = NULL;
	Bbond*		bondlist = NULL;
	
	Bmolecule*	mol;
	Bresidue*	res, *res2;
	Batom*  	atom, *atom2;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Generating bonds:" << endl;
		cout << "Bond length cutoff:             " << bondlength << " A" << endl;
	}
	
    for ( natom=0, mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				natom++;
				for( res2 = res; res2; res2 = res2->next ) {
					for ( atom2 = res2->atom; atom2; atom2 = atom2->next ) {
						if ( atom < atom2 ) {
							if ( wrap )
								dist = (vector3_difference_PBC(atom->coord, atom2->coord, molgroup->box)).length();
							else
								dist = (atom->coord - atom2->coord).length();
							if ( dist > 0.1*bondlength && dist < 1.1*bondlength ) {
								if ( bond ) bond = bond_add(&bond, atom, atom2, bondlength, 1);
								else bond = bond_add(&bondlist, atom, atom2, bondlength, 1);
								if ( verbose & VERB_FULL )
									cout << "Bond " << nbond << ": " << dist << " " 
										<< atom->coord[2] << " " << atom2->coord[2] << endl;
								nbond++;
							}
						}
					}
				}
			}
		}
	}
			
	if ( verbose & VERB_PROCESS )
		cout << "Number of bonds generated:      " << nbond << endl << endl;
	
	molgroup->bond = bondlist;
	
	return bondlist;
}

/**
@brief 	Converts molecules in a molecule group to individual molecule groups.

	A new linked list of molecule groups is created and the links to the
	individual molecules set. 

@param 	*molgroup molecule group structure (modified).
@return int					0.
**/
int			molecules_to_molgroups(Bmolgroup* molgroup)
{
	Bmolgroup*	mg = molgroup;
	Bmolecule*	mol;
	
	for ( mol = molgroup->mol->next; mol; mol = mol->next ) {
		mg->next = molgroup_init();
		mg = mg->next;
		mg->comment = molgroup->comment;
		mg->id = mol->id;
		mg->pointgroup = molgroup->pointgroup;
		mg->mol = mol;
	}
	
	for ( mg = molgroup; mg; mg = mg->next ) if ( mg->mol ) mg->mol->next = NULL;
	
	return 0;
}

