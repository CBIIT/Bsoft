/**
@file	rwmd.h
@brief	Header to read and write molecular dynamics parameters in STAR format
@author Bernard Heymann
@date	Created: 20030919
@date	Modified: 20060424
**/

#include "rwatomprop.h"
#include "Vector3.h"

#ifndef _Bmd_
#define _Bmd_
/************************************************************************
@Object: struct Bbondtype
@Description:
	A structure for a bond type.
@Features:
	This defines all the properties of a bond type.
*************************************************************************/
struct Bbondtype {
	Bbondtype*	next;
	char		type1[8];		// Atom type 1
	char		type2[8];		// Atom type 2
	float		covlength;		// Covalent bond length
	float		vdwdist;		// Van der Waals interaction distance
	long		number;			// Number of this type
	float		std;			// Standard deviation of covalent bond length
} ;

/************************************************************************
@Object: struct Bangletype
@Description:
	A structure for an angle type.
@Features:
	This defines all the properties of an angle type.
*************************************************************************/
/*struct Bangletype {
	Bangletype*	next;
	char		type1[8];		// Atom type 1
	char		type2[8];		// Atom type 2
	char		type3[8];		// Atom type 3
	float		angle;			// Angle
	long		number;			// Number of this type
	float		std;			// Standard deviation of angle
} ;
*/
/************************************************************************
@Object: struct Bmd
@Description:
	A structure used for molecular mechanics.
@Features:
	This defines variables and constants for doing molecular dynamics or
	other molecular mechanics.
*************************************************************************/
struct Bmd {
	double		timestep;		// Integration time step for MD
	double		Kfriction;		// Friction coefficient
	double		Kbond;			// Bond energy constant
	double		Kangle;			// Angle energy constant
	double		Kelec;			// Electrostatic energy constant
	double		Kvdw;			// Van der Waals energy constant
	double		Ksep;			// Separation distance energy constant
	double		Kpoint;			// Point force energy constant
	double		Kmap;			// Density map energy constant
	double		VdWcoeff1;		// First Van der Waals coefficient (default 1/12)
	double		VdWcoeff2;		// Second Van der Waals coefficient (default 1/6)
	double		sepdist;		// Separation distance grid sampling
	double		cutoff;			// Distance cutoff for non-bonded calculations
	double		pointdecay;		// Decay constant for the point force
	int			bondsteps;		// Number of sampling intervals along a bond
	int			wrap;			// Flag to turn periodic boundaries on
	double		Ebond;			// Bond energy
	double		Eangle;			// Angle energy
	double		Eelec;			// Electrostatic energy
	double		Evdw;			// Van der Waals energy
	double		Esep;			// Separation distance energy (overlap penalty)
	double		Epoint;			// Point force energy
	double		Emap;			// Density map associated energy
	double		Ekin;			// Kinetic energy
	double		Epot;			// Potential energy
	Vector3<float>	point;		// Center of point force
	Batomtype*	atom;			// Linked list of atom types
	Bbondtype*	bond;			// Linked list of bond types
//	Bangletype*	angle;			// Linked list of angle types
} ;
#endif

// Function prototypes
Bmd*		md_init();
Bmd*		md_init_with_types();
Bmd*		read_md_parameters(Bstring& filename);
int			write_md_parameters(Bstring& filename, Bmd* md);
int			md_kill(Bmd* md);

