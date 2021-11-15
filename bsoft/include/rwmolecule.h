/**
@file	rwmolecule.h
@brief	Header file for reading reflection files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20170509
**/

#include "UnitCell.h"
#include "Vector3.h"
#include "Bstring.h"

/* Constants */
#define MAXSEQLEN	1000000	// Maximum sequence length

#ifndef _Bmolecule_
/************************************************************************
@Object: struct Batom
@Description:
	Atom property structure.
@Features:
	PDB line data structure: items labeled extra not in the PDB specification.
*************************************************************************/
struct Batom {
    Batom*			next;	 	// Next atom in linked list
    int				num;    	// Atom number
    char			el[4];		// Atom type
    char			type[8];	// Atom type
	int				tnum;		// Atom type number (for reference purposes)
    Vector3<double>	coord;		// Coordinates
    double			q,b;    	// Occupancy and temperature factor
    double			mass;		// Atom mass
    double			chrg;		// Charge on atom
    int				sel;    	// Select flag - extra
	Vector3<double>	vel;		// Velocities (for molecular dynamics)
	Vector3<double>	F;			// Force vector (for molecular dynamics)
} ;

/************************************************************************
@Object: struct Bresidue
@Description:
	Residue definition and property structure.
@Features:
	A residue is viewed as equivalent to the monomeric units incorportated
	into a polymer (amino acids for proteins and nucleotides for nucleic
	acids). For small molecules, the residue and molecule definitions
	are identical.
*************************************************************************/
struct Bresidue {
	Bresidue*	next;		// Next residue in linked list
	int			num;		// Residue number
	char		insert[2];	// Insertion code
	char		type[6];	// Residue type
	Batom*		atom;		// First atom in residue
} ;

/************************************************************************
@Object: enum SecondaryType
@Description:
	Protein secondary structure type specifier.
@Features:
	This is derived from the PDB definitions of secondary structure.
*************************************************************************/
enum SecondaryType {
	RightHandedAlpha = 1,
	RightHandedOmega = 2,
	RightHandedPi = 3,
	RightHandedGamma = 4,
	RightHanded310 = 5,
	LeftHandedAlpha = 6,
	LeftHandedOmega = 7,
	LeftHandedGamma = 8,
	TwoSevenRibbon = 9,
	Polyproline = 10,
	Strand = 11,
	Turn = 12
} ;

/************************************************************************
@Object: struct Bsecondary
@Description:
	Protein secondary structure definition.
@Features:
	Every definition of a secondary structure element includes a type
	and the sequence range.
*************************************************************************/
struct Bsecondary {
	Bsecondary*	next;			// Next secondary structure item in linked list
	SecondaryType type;			// Class
	int			num;			// Serial number
	char		id[4];			// Identifier: Different series for helix, sheet and turn
	int			sense;			// Strand direction
	int			nstrands;		// Number of strands in the sheet
	Bresidue*	first;			// First residue
	Bresidue*	last;			// Last residue
} ;

/************************************************************************
@Object: struct Bmolecule
@Description:
	Molecule definition and property structure.
@Features:
	A molecule is viewed as a single covalently bound set of atoms with 
	properties such as gene and protein sequences, and composed
	of a set of residues containing atoms.
*************************************************************************/
struct Bmolecule {
	Bmolecule*		next;			// Next molecule in linked list
	int				num;			// Molecule number
	Bstring			id;				// Molecule identifier
	long 			nres;			// Number of residues
	long 			nbase;			// Number of bases
	long			first_codon;	// First amino acid codon in NA sequence
	Bstring			seq;			// Protein sequence
	Bstring			naseq;			// Nucleic acid sequence
	Vector3<double>	min; 			// Coordinate minima
	Vector3<double>	max;	 		// Coordinate maxima
	double			Bmin,Bmax;		// B-factor statistics
	double			Bave,Bstd;
	int				sel;			// Selection flag
	double			fom;			// Figure-of-merit for molecule
	Bresidue*		res;			// First residue in molecule
	Bsecondary*		sec;			// Secondary structure
} ;

/************************************************************************
@Object: struct Bbond
@Description:
	All bonds defined for a molecule group.
@Features:
	Intra- and intermolecular bonds are defined by pointers to two atoms each.
	The structure also holds bond length and strength properties.
*************************************************************************/
struct Bbond {
	Bbond*		next;	 	// Next bond in linked list
	Batom*		atom1;		// First atom
	Batom*		atom2;		// Second atom
	double		l;			// Reference bond length
	double		k;			// Bond strength
} ;

/************************************************************************
@Object: struct Bangle
@Description:
	All angles defined for a molecule group.
@Features:
	The angle between two bonds are defined by the three connected atoms.
	The structure also holds reference angle length and strength properties.
*************************************************************************/
struct Bangle {
	Bangle*		next;	 	// Next angle in linked list
	Batom*		atom1;		// First atom
	Batom*		atom2;		// Second atom
	Batom*		atom3;		// Third atom
	double		a;			// Reference angle
	double		k;			// Angle strength
} ;


/************************************************************************
@Object: struct Bmolgroup
@Description:
	A collection of molecules.
@Features:
	This holds any set of molecules such as:
	1.	Members of a complex.
	2.	Homologues used for sequence or structural similarity analysis.
	3.	All molecules in a molecular dynamics setup.
*************************************************************************/
struct Bmolgroup {
	Bmolgroup*		next;			// Next molecule group in linked list
	Bstring			filename;		// File name
	int				num;			// Molecule group number
	Bstring			id;				// Molecule group identifier
	Bstring			comment;		// Comment string
	int				spacegroup;		// Crystallographic space group number
	Bstring			pointgroup;		// Point group symmetry for molecule group
	Bstring			sgstring;		// Crystallographic space group notation
	UnitCell		unitcell;		// Crystallographic unit cell
	long			maxlen; 		// Maximum protein sequence length
	long			maxnalen; 		// Maximum nucleic acid sequence length
	Vector3<double>	min; 			// Coordinate minima
	Vector3<double>	max;	 		// Coordinate maxima
	Vector3<double>	box;	 		// Bounding box details
	Bstring			select;			// Selection options
	int				sel;			// Select flag
	double			fom;			// Figure-of-merit for molecule group
	int*			seqflag; 		// Position selection flag
	Bmolecule*		mol;			// First molecule in group
	Bbond*			bond;			// First bond for group
	Bangle*			angle;			// First angle for group
} ;
#define _Bmolecule_
#endif

/* Function prototypes */
Bmolgroup*  molgroup_init();
Bmolecule*	molecule_add(Bmolecule** mol, char* name);
Bmolecule*	molecule_add(Bmolecule** mol, Bstring& name);
Bresidue*	residue_add(Bresidue** res, const char* type);
Bresidue*	residue_add(Bresidue** res, Bstring& type);
Batom*		atom_add(Batom** atom, const char* type);
Batom*		atom_add(Batom** atom, Bstring& type);
long		residue_count(Bmolgroup* molgroup);
long		atom_count(Bmolgroup* molgroup);
int			atom_clean_type(Batom* atom, const char* type);
Bbond*		bond_add(Bbond** bond, Batom* atom1, Batom* atom2, double l, double k);
Bangle*		angle_add(Bangle** angle, Batom* atom1, Batom* atom2, Batom* atom3, double a, double k);
int 		molgroup_list_kill(Bmolgroup* molgroup);
int 		molgroup_kill(Bmolgroup* molgroup);
int 		molecule_kill(Bmolecule* mol);
int 		residue_kill(Bresidue* res);
int			bond_kill(Bbond* bond);
int			angle_kill(Bangle* angle);
Bmolgroup*	molgroup_list_copy(Bmolgroup* molgroup);
Bmolgroup*	molgroup_copy(Bmolgroup* molgroup);
Bmolecule* 	molecule_copy(Bmolecule* mol);
Bmolecule* 	mol_copy_and_add_to_molgroup(Bmolgroup* molgroup, Bmolecule* mol);
Bbond*		molgroup_bond_list_copy(Bmolgroup* molgroup, Bmolgroup* molgroupcopy);
int			molgroup_from_molgroup_list(Bmolgroup* molgroup);
Bmolgroup*  read_molecule(const char *filename, const char *select, const char* paramfile);
Bmolgroup*  read_molecule(Bstring& filename, Bstring& atom_select, Bstring& paramfile);
Bmolgroup*	read_molecule(Bstring* file_list, int set_pbc, Vector3<double> box, 
				Bstring atom_select, Bstring paramfile);
int 		write_molecule(char *filename, Bmolgroup* molgroup);
int 		write_molecule(Bstring& filename, Bmolgroup *molgroup);
int			molgroup_list_write(Bstring& filename, Bmolgroup* molgroup);
long 		molgroup_count_molecules(Bmolgroup* molgroup);
long 		molgroup_count_residues(Bmolgroup* molgroup);
long 		mol_count_residues(Bmolecule *mol);
long 		molgroup_count_atoms(Bmolgroup* molgroup);
long 		mol_count_atoms(Bmolecule *mol);
int			molgroup_consolidate_gaps(Bmolgroup* molgroup);
long 		molgroup_stats(Bmolgroup* molgroup, int show);
long 		molgroup_stats(Bmolgroup* molgroup);
long 		mol_stats(Bmolecule* mol, int show);
long 		mol_stats(Bmolecule* mol);
int 		molecule_update_comment(Bmolgroup* molgroup, int n, char** strings);
int 		molecule_get_masses(Bmolgroup* molgroup, Bstring& paramfile);
Bbond*		molgroup_bond_list_generate(Bmolgroup* molgroup, double maxlength, int wrap);
Bbond*		mol_bond_list_generate(Bmolgroup* molgroup, double bondlength, int wrap);
int			molecules_to_molgroups(Bmolgroup* molgroup);

