/**
@file	rwatomprop.h
@brief	Header file for reading atom property files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20110810
**/

#include <Bstring.h>

#ifndef _Batomtype_
/************************************************************************
@Object: struct Batomtype
@Description:
	A structure for an atom type.
@Features:
	This defines all the properties of an atom type.
*************************************************************************/
struct Batomtype {
	Batomtype*		next;		// Link to next atom type if not NULL
	int				z;			// Atomic number
	char			el[4];		// Element
	char			name[8];	// Atom name
	float			mass;		// Atomic mass
	float			oxid;		// Oxidation state
	float			bond;		// Bond length
	float			vdw;		// VdW radius
	unsigned char	color[4];	// RGBA colour
	float			sfa[5];		// Cromer-Mann coefficients
	float			sfb[5];
	float			sfc;
	long			number;		// Number of this type
} ;

#define _Batomtype_
#endif

/* Function prototypes */
Batomtype*	get_atom_properties(Bstring& filename);
int 		write_atom_properties(Bstring& filename, Batomtype* at);
Batomtype*	atom_type_add(Batomtype** atype, const char* aname);
int			atom_type_kill(Batomtype* atype);


