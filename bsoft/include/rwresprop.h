/**
@file	rwresprop.h
@brief	Header file for reading residue property files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20100916
**/

#include "Bstring.h"

/* Constants */
#define MAXRES	25 		// Maximum number of different residues

#ifndef _Brestype_
/************************************************************************
@Object: struct Bresidue_type
@Description:
	Residue type property structure.
@Features:
	A residue type may be identified by either a one- or three-letter code.
*************************************************************************/
struct Bresidue_type {
	Bresidue_type*	next;	// Link to next residue type if not NULL
	char		c;			// Single letter representation
	char		cod[3];		// Three-letter representation
	float		mass;		// Residue mass
	float		vol;		// Residue volume
	float		ext;		// Residue extension
	float		extsd;		// Residue extension standard deviation
	float		charge; 	// Residue charge
	float		hphob;		// Residue hydrophobicity
	float		comp[5];	// Composition: HCNOS
} ;

/************************************************************************
@Object: struct Bresidue_matrix
@Description:
	Residue relationship matrix.
@Features:
	A pairwise residue relationship matrix to encode a property such 
	as similarity.
*************************************************************************/
/* Residue relationship matrix (similarity) */
struct Bresidue_matrix {
	int 		n;			// Number of symbols for matrix
    Bstring		c;			// Symbol list for matrix (n characters)
	float*		m; 			// nxn Matrix (similarity)
} ;

#define _Brestype_
#endif

/* Function prototypes */
Bresidue_type*	get_residue_properties(Bstring& filename);
int 			write_residue_properties(Bstring& filename, Bresidue_type* rt);
Bresidue_matrix*	get_residue_matrix(Bstring& filename);
int				residue_matrix_kill(Bresidue_matrix* resmat);


