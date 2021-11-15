/**
@file	mol_compare.cpp
@brief	Library routines used to compare sets of atomic coordinates
@author Bernard Heymann
@date	Created: 20021020
@date	Modified: 20150209
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "mol_symmetry.h"
#include "mol_compare.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "seq_align.h"
#include "Matrix3.h"
#include "matrix_linear.h"
#include "Matrix.h"
#include "symmetry.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Copies and rotates the molecule group and compares it with the original.
@param 	*molgroup		molecule group.
@param 	t				rotation operation.
@return double			RMSD.
**/
double		molgroup_rotate_and_compare(Bmolgroup* molgroup, Transform t)
{
	Bmolgroup*		mol_rot;
	double			R;

	t.axis.normalize();
	mol_rot = molgroup_copy(molgroup);
	molgroup_coor_rotate(mol_rot, t);
	R = molgroup_calc_brute_rmsd(molgroup, mol_rot);
	molgroup_kill(mol_rot);
	
	return R;
}

/**
@brief	Determines the transformation between two groups of identical molecules.
@param 	*molgroup1		first molecule group.
@param 	*molgroup2		second molecule group.
@return Transform		transform structure.
**/
Transform	molgroup_find_transformation(Bmolgroup* molgroup1, Bmolgroup* molgroup2)
{
	long			i, j, n(0);
//	double			bx[4], by[4], bz[4], v[4] = {0,0,0,1};
	vector<double>	bx(4), by(4), bz(4), v(4,0);
	Matrix			a(4,4);
	
	v[3] = 1;
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	Bmolecule		*mol1, *mol2;
	Bresidue		*res1, *res2;
	Batom			*atom1, *atom2;

	Vector3<double> 	com = molgroup_center_of_mass(molgroup1);
	Vector3<double> 	c1, c2;
	
	for ( n=0, mol1=molgroup1->mol, mol2=molgroup2->mol; mol1 && mol2; mol1=mol1->next, mol2=mol2->next ) {
		for ( res1=mol1->res, res2=mol2->res; res1 && res2; res1=res1->next, res2=res2->next ) {
			for ( atom1=res1->atom, atom2=res2->atom; atom1 && atom2; atom1=atom1->next, atom2=atom2->next ) {
				if ( strncmp(atom1->el, atom2->el, 2) == 0 ) {
					c1 = atom1->coord - com;
					v[0] = c1[0];
					v[1] = c1[1];
					v[2] = c1[2];
					c2 = atom2->coord - com;
					for ( i=0; i<4; i++ ) {
						bx[i] += v[i]*c2[0];
						by[i] += v[i]*c2[1];
						bz[i] += v[i]*c2[2];
						for ( j=0; j<=i; j++ ) a[i][j] += v[i]*v[j];
					}
					n++;
				}
			}
		}
	}
	
	for ( i=0; i<3; i++ )
		for ( j=i+1; j<4; j++ ) a[i][j] = a[j][i];
	
//	a.LU_decomposition();
	a.singular_value_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	a.multiply_in_place(bz);

	Matrix3			mat(bx[0],bx[1],bx[2],by[0],by[1],by[2],bz[0],bz[1],bz[2]);
//	mat = mat.transpose();
	mat.normalize();
	Transform		t(mat);
	t.origin = com;
	t.trans = Vector3<double>(bx[3], by[3], bz[3]);

	if ( verbose & VERB_PROCESS ) {
		cout << "Atoms compared:                 " << n << endl;
		matrix3_show_hp(mat);
		cout << endl;
	}

	return t;
}

/**
@brief	Determines the transformation between two identical molecules.
@param 	*mol1			first molecule.
@param 	*mol2			second molecule.
@param 	offset			offset of second sequence with respect to first.
@return Transform		transform structure.

	The transformation:
		coord2 = rot_mat * (coord1 - origin) + origin + shift
	is solved.
	The algorithm is set up in parts, solving first for:
		coord2 = rot_mat * coord1 + shift_temp
	The last term is given by:
		shift_temp = shift + origin - rot_mat * origin
	The shift vector must be parallel to the rotation axis, and is
	determined as:
		shift = axis * |shift_temp| * cos(alpha)
	where alpha is the angle between shift_temp and the rotation axis.
	The origin is then calculated from:
		origin = (shift_temp - shift) * inverse(id_mat - rot_mat)
	Note that the origin still has one degree of freedom: It can be
	anywhere along the rotation axis.

**/
Transform	mol_find_transformation(Bmolecule* mol1, Bmolecule* mol2, int offset)
{
	long			i, j, n(0);
//	double			bx[4], by[4], bz[4], v[4] = {0,0,0,1};
	vector<double>	bx(4), by(4), bz(4), v(4,0);
	Matrix			a(4,4);
	
	v[3] = 1;
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	Batom			*atom1, *atom2;
	Bresidue*		res1=mol1->res;
	Bresidue*		res2=mol2->res;
	if ( offset < 0 ) for ( i=0; i<-offset; i++ ) res2=res2->next;
	else for ( i=0; i<offset; i++ ) res1=res1->next;
	
	res1=mol1->res;
	res2=mol2->res;
	if ( offset < 0 ) for ( i=0; i<-offset; i++ ) res2=res2->next;
	else for ( i=0; i<offset; i++ ) res1=res1->next;

	Vector3<double> 	com = mol_center_of_mass(mol1);
	Vector3<double> 	c1, c2;
	
	for ( n=0; res1 && res2; res1=res1->next, res2=res2->next ) {
		for ( atom1=res1->atom, atom2=res2->atom; atom1 && atom2; atom1=atom1->next, atom2=atom2->next ) {
			if ( strncmp(atom1->el, atom2->el, 2) == 0 ) {
				c1 = atom1->coord - com;
				v[0] = c1[0];
				v[1] = c1[1];
				v[2] = c1[2];
				c2 = atom2->coord - com;
				for ( i=0; i<4; i++ ) {
					bx[i] += v[i]*c2[0];
					by[i] += v[i]*c2[1];
					bz[i] += v[i]*c2[2];
					for ( j=0; j<=i; j++ ) a[i][j] += v[i]*v[j];
				}
				n++;
			}
		}
	}
	
	for ( i=0; i<3; i++ )
		for ( j=i+1; j<4; j++ ) a[i][j] = a[j][i];
	
//	a.LU_decomposition();
	a.singular_value_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	a.multiply_in_place(bz);

	Matrix3			mat(bx[0],bx[1],bx[2],by[0],by[1],by[2],bz[0],bz[1],bz[2]);
//	Transform		t = transform_from_matrix3(mat);
	mat.normalize();
	Transform		t(mat);
	t.origin = com;
	t.trans = Vector3<double>(bx[3], by[3], bz[3]);

//	matrix3_show_hp(mat);
	
	return t;
}

/**
@brief	Calculates the root-mean-square-deviation between two molecule groups.  
@param 	*molgroup1		first molecule group.
@param 	*molgroup2		second molecule group.
@return	double			root-mean-square-deviation.

	The two molecule groups must have identical structures.

**/
double		molgroup_calculate_rmsd(Bmolgroup* molgroup1, Bmolgroup* molgroup2)
{
	long			n(0);
	double			d, ad(0), sd(0), R(0);
	Bmolecule		*mol1, *mol2;
	Bresidue		*res1, *res2;
	Batom			*atom1, *atom2;
	
//	if ( verbose & VERB_FULL )
//		cout << "Atom\tDistance" << endl;
	for ( mol1=molgroup1->mol, mol2=molgroup2->mol; mol1 && mol2; mol1=mol1->next, mol2=mol2->next ) {
		for ( res1=mol1->res, res2=mol2->res; res1 && res2; res1=res1->next, res2=res2->next ) {
			for ( atom1=res1->atom, atom2=res2->atom; atom1 && atom2; atom1=atom1->next, atom2=atom2->next ) {
				if ( strncmp(atom1->el, atom2->el, 2) == 0 ) {
					d = atom1->coord.distance2(atom2->coord);
					sd += d;
					d = sqrt(d);
					ad += d;
					n++;
//					if ( verbose & VERB_FULL )
//						cout << atom1->num << tab << d << endl;
				}
			}
		}
	}
	
	if ( n ) {
		ad /= n;
		sd /= n;
		if ( sd > 0 ) R = sqrt(sd);
		sd -= ad*ad;
		if ( sd > 0 ) sd = sqrt(sd);
		else sd = 0;
	}
	
	Vector3<double>	comd = molgroup_center_of_mass(molgroup2) - molgroup_center_of_mass(molgroup1);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of comparisons:          " << n << endl;
		cout << "Center-of-mass difference:      " << comd << " " << comd.length() << endl;
		cout << "Distance avg & std:             " << ad << " " << sd << " A" << endl;
		cout << "RMSD:                           " << R << " A" << endl << endl;
	}
	
	return R;
}

/**
@brief	Calculates the RMSD between two molecules.
@param 	*mol1			first molecule.
@param 	*mol2			second molecule.
@return double 			RMSD.

	The root-mean-square-deviation between two sets of corrdinates is given by:
		R = sqrt(sum(length(coord1-coord2))/number).

**/
double		mol_calculate_rmsd(Bmolecule* mol1, Bmolecule* mol2)
{
	long			n(0);
	double			d, ad(0), sd(0), R(0);
	Bresidue		*res1, *res2;
	Batom			*atom1, *atom2;
	
	for ( res1=mol1->res, res2=mol2->res; res1 && res2; res1=res1->next, res2=res2->next ) {
		for ( atom1=res1->atom, atom2=res2->atom; atom1 && atom2; atom1=atom1->next, atom2=atom2->next ) {
			if ( strncmp(atom1->el, atom2->el, 2) == 0 ) {
				d = atom1->coord.distance2(atom2->coord);
				sd += d;
				ad += sqrt(d);
				n++;
			}
		}
	}
	
	if ( n ) {
		ad /= n;
		sd /= n;
		if ( sd > 0 ) R = sqrt(sd);
		sd -= ad*ad;
		if ( sd > 0 ) sd = sqrt(sd);
		else sd = 0;
	}
	
	Vector3<double>	comd = mol_center_of_mass(mol2) - mol_center_of_mass(mol1);
	
	if ( verbose & VERB_FULL ) {
		cout << "Center-of-mass difference:      " << comd << " " << comd.length() << endl;
		cout << "Distance avg & std:             " << ad << " " << sd << " A" << endl;
		cout << "RMSD:                           " << R << " A" << endl << endl;
	}
	
	return R;
}

/**
@brief	Calculates the RMSD based on nearest atoms.
@param 	*molgroup1		first molecule group.
@param 	*molgroup2		second molecule group.
@return double			RMSD.
**/
double		molgroup_calc_brute_rmsd(Bmolgroup* molgroup1, Bmolgroup* molgroup2)
{
	Bmolecule		*mol1, *mol2;
	Bresidue		*res1, *res2;
	Batom			*atom1, *atom2;
	long			n(0);
	double			d, dd, R(0);

    for ( mol1 = molgroup1->mol; mol1; mol1 = mol1->next ) {
		for( res1 = mol1->res; res1; res1 = res1->next ) {
			for ( atom1 = res1->atom; atom1; atom1 = atom1->next ) {
				d = 1e30;
				for ( mol2 = molgroup2->mol; mol2; mol2 = mol2->next ) {
					for( res2 = mol2->res; res2; res2 = res2->next ) {
						for ( atom2 = res2->atom; atom2; atom2 = atom2->next ) {
							dd = atom1->coord.distance(atom2->coord);
							if ( d > dd ) d = dd;
						}
					}
				}
				R += d*d;
				n++;
			}
		}
	}

	R = sqrt(R/n);
	
	return R;
}

/**
@brief 	Calculates the distance matrix between the residues in two molecules.
@param 	*m1				first molecules structure.
@param 	*m2				second molecules structure.
@return Matrix			distance matrix.

	The matrix is calculated from the pairwise distances between residues. 

**/
Matrix		mol_distance_matrix(Bmolecule* m1, Bmolecule* m2)
{
	long			i, j;
	Bresidue		*res1, *res2;
	Batom			*atom1, *atom2;
	
	Matrix			mat(m1->nres, m2->nres);

	for ( i=0, res1 = m1->res; res1; res1 = res1->next, ++i ) {
		for ( atom1 = res1->atom; atom1 && strncmp(atom1->type, "CA", 2); atom1 = atom1->next ) ;
		if ( atom1 ) {
			for ( j=0, res2 = m2->res; res2; res2 = res2->next, ++j ) {
				for ( atom2 = res2->atom; atom2 && strncmp(atom2->type, "CA", 2); atom2 = atom2->next ) ;
				if ( atom2 )
					mat[i][j] = atom1->coord.distance(atom2->coord);
			}
		}
	}
	
	return mat;
}
