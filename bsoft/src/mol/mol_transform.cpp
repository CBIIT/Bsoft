/**
@file	mol_transform.cpp
@brief	Library routines used for atomic coordinate transformations
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20070614
**/

#include "rwmolecule.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
Bmolgroup*	molgroup_rotate(Bmolgroup* molgroup, Matrix3 mat, Vector3<double> origin, Vector3<double> trans);
Bmolecule*	mol_rotate(Bmolecule* mol, Matrix3 mat, Vector3<double> origin, Vector3<double> trans);

/**
@brief 	Translates a molecule group.
@param 	*molgroup		molecule group structure.
@param 	shift	three-valued translation vector.
@return int						0.
**/
int 		molgroup_coor_shift(Bmolgroup* molgroup, Vector3<double> shift)
{
	if ( shift.length() < 1e-30 ) return 0;

	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL )
		cout << "Shifting coordinates:           " << shift << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord + shift;
			}
		}
    }
	
	return 0;
}

/**
@brief 	Translates a molecule.
@param 	*mol		molecule structure.
@param 	shift		three-valued translation vector.
@return int 				0.
**/
int 		mol_coor_shift(Bmolecule* mol, Vector3<double> shift)
{
	if ( shift.length() < 1e-30 ) return 0;

	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL )
		cout << "Shifting coordinates:           " << shift << endl;
	
    for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			atom->coord = atom->coord + shift;
		}
    }
	
//	mol_stats(mol);
	
	return 0;
}

/**
@brief 	Rotates and translates a molecule group.

	The new coordinates are:
		coord_new = (coord - origin)*rot_mat + origin + shift.

@param 	*molgroup molecule group structure.
@param 	t			transform structure.
@return int 				0.
**/
int 		molgroup_coor_rotate(Bmolgroup* molgroup, Transform t)
{
	if ( ( t.angle == 0 || t.axis.length() < 1e-10 ) && ( t.trans.length() < 1e-10 ) )return 0;

	Matrix3		mat = Matrix3(t.axis, t.angle);
	
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL ) {
		cout << endl << "Rotating and translating coordinates:" << endl;
		cout << "Rotation axis and angle:        " << t.axis << " " << t.angle*180.0/M_PI << endl;
		cout << "Rotation origin:                " << t.origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << t.trans << endl << endl;
	}
		
	Vector3<double>		shift = t.origin + t.trans;
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord - t.origin;
				atom->coord = mat * atom->coord;
				atom->coord = atom->coord + shift;
		    }
		}
    }
	
	return 0;
}

/**
@brief 	Rotates and translates a molecule.

	The new coordinates are:
		coord_new = (coord - origin)*rot_mat + origin + shift.

@param 	*mol		molecule structure.
@param 	t			transform structure.
@return int 				0.
**/
int 		mol_coor_rotate(Bmolecule* mol, Transform t)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mol_coor_rotate: t.angle=" << t.angle << endl;
		
	if ( t.angle < 1e-20 ) return 0;

	if ( t.axis.length() < 1e-10 ) return 0;

	Matrix3		mat = Matrix3(t.axis, t.angle);
	
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL ) {
		cout << endl << "Rotating and translating coordinates:" << endl;
		cout << "Rotation axis and angle:        " << t.axis << " " << t.angle*180.0/M_PI << endl;
		cout << "Rotation origin:                " << t.origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << t.trans << endl << endl;
	}
		
	Vector3<double>		shift = t.origin + t.trans;
    for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			atom->coord = atom->coord - t.origin;
			atom->coord = mat * atom->coord;
			atom->coord = atom->coord + shift;
		}
    }
	
//	mol_stats(mol);
	
	return 0;
}

/**
@brief 	Rotates a molecule group to a specified view and translates it afterwards.

	A new rotated and translated molecule group is generated.

@param 	*molgroup		molecule group structure (unmodified).
@param 	view				view to rotate to.
@param 	origin	origin of rotation.
@param 	trans	3-valued translation vector.
@return Bmolgroup*				new molecule group.
**/
Bmolgroup*	molgroup_rotate_to_view(Bmolgroup* molgroup, View view, Vector3<double> origin, Vector3<double> trans)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating and translating coordinates:" << endl;
		cout << "Rotation to view:               " << view << endl;
	}
	
	Matrix3			mat = view.matrix();

	return molgroup_rotate(molgroup, mat, origin, trans);
}

/**
@brief 	Rotates a molecule group from a specified view and translates it afterwards.

	A new rotated and translated molecule group is generated.

@param 	*molgroup		molecule group structure (unmodified).
@param 	view				view to rotate from.
@param 	origin	origin of rotation.
@param 	trans	3-valued translation vector.
@return Bmolgroup*				new molecule group.
**/
Bmolgroup*	molgroup_rotate_from_view(Bmolgroup* molgroup, View view, Vector3<double> origin, Vector3<double> trans)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating and translating coordinates:" << endl;
		cout << "Rotation from view:             " << view << endl;
	}
	
	Matrix3			mat = view.matrix();
	mat = mat.transpose();
	
	return molgroup_rotate(molgroup, mat, origin, trans);
}

Bmolgroup*	molgroup_rotate(Bmolgroup* molgroup, Matrix3 mat, Vector3<double> origin, Vector3<double> trans)
{
	Bmolgroup*		mol_rot = molgroup_copy(molgroup);

	Vector3<double>	shift = origin + trans;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	if ( verbose & VERB_FULL ) {
		cout << "Origin:                         " << origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << trans << endl << endl;
	}
	
    for ( mol = mol_rot->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord - origin;
				atom->coord = mat * atom->coord;
				atom->coord = atom->coord + shift;
		    }
		}
    }
	
	return mol_rot;
}

/**
@brief 	Rotates a molecule to a specified view and translates it afterwards.

	A new rotated and translated molecule is generated.

@param 	*mol			molecule structure (unmodified).
@param 	view				view to rotate to.
@param 	origin	origin of rotation.
@param 	trans	3-valued translation vector.
@return Bmolecule*				new molecule group.
**/
Bmolecule*	mol_rotate_to_view(Bmolecule* mol, View view, Vector3<double> origin, Vector3<double> trans)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating and translating coordinates:" << endl;
		cout << "Rotation to view:               " << view << endl;
	}

	Matrix3			mat = view.matrix();

	return mol_rotate(mol, mat, origin, trans);
}

/**
@brief 	Rotates a molecule from a specified view and translates it afterwards.

	A new rotated and translated molecule is generated.

@param 	*mol			molecule structure (unmodified).
@param 	view				view to rotate from.
@param 	origin	origin of rotation.
@param 	trans	3-valued translation vector.
@return Bmolecule*				new molecule group.
**/
Bmolecule*	mol_rotate_from_view(Bmolecule* mol, View view, Vector3<double> origin, Vector3<double> trans)
{
	if ( verbose & VERB_FULL ) {
		cout << "Rotating and translating coordinates:" << endl;
		cout << "Rotation from view:             " << view << endl;
	}
	
	Matrix3			mat = view.matrix();
	mat = mat.transpose();
	
	return mol_rotate(mol, mat, origin, trans);
}

Bmolecule*	mol_rotate(Bmolecule* mol, Matrix3 mat, Vector3<double> origin, Vector3<double> trans)
{
	Bmolecule*		mol_rot = molecule_copy(mol);
	
	Vector3<double>	shift = origin + trans;
	
	Bresidue*		res;
	Batom*  		atom;
	
	if ( verbose & VERB_FULL ) {
		cout << "Origin:                         " << origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << trans << endl << endl;
	}
	
    for( res = mol_rot->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			atom->coord = atom->coord + origin;
			atom->coord = mat * atom->coord;
			atom->coord = atom->coord + shift;
		}
    }
	
	return mol_rot;
}

/**
@brief 	Rotates and translates a molecule group.

	The new coordinates are:
		coord_new = (coord - origin)*rot_mat + origin + shift.

@param 	*molgroup molecule group structure.
@param 	t			transform structure.
@return int 				0.
**/
int 		molgroup_coor_transform(Bmolgroup* molgroup, Transform t)
{
	Matrix3		mat = Matrix3(t.axis, t.angle);
	
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL ) {
		cout << "Transforming coordinates:" << endl;
		cout << "Rotation axis and angle:        " << t.axis << " " << t.angle*180.0/M_PI << endl;
		cout << "Rotation origin:                " << t.origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << t.trans << endl;
		cout << "Scale:                          " << t.scale << endl << endl;
	}
	
	mat = t.scale * mat;
	
	Vector3<double>		shift = t.origin + t.trans;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord - t.origin;
				atom->coord = mat * atom->coord;
				atom->coord = atom->coord + shift;
		    }
		}
    }
	
	return 0;
}

/**
@brief 	Rotates and translates a molecule.

	The new coordinates are:
		coord_new = (coord - origin)*rot_mat + origin + shift.

@param 	*mol		molecule structure.
@param 	t			transform structure.
@return int 				0.
**/
int 		mol_coor_transform(Bmolecule* mol, Transform t)
{
	Matrix3		mat = Matrix3(t.axis, t.angle);
	
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_FULL ) {
		cout << "Transforming and translating coordinates:" << endl;
		cout << "Rotation axis and angle:        " << t.axis << " " << t.angle*180.0/M_PI << endl;
		cout << "Rotation origin:                " << t.origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation:                    " << t.trans << endl;
		cout << "Scale:                          " << t.scale << endl << endl;
	}
		
	mat = t.scale * mat;
	
	Vector3<double>		shift = t.origin + t.trans;
	
    for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			atom->coord = atom->coord - t.origin;
			atom->coord = mat * atom->coord;
			atom->coord = atom->coord + shift;
		}
    }
	
//	mol_stats(mol);
	
	return 0;
}

/**
@brief 	Inverts the coordinates of a molecule through a given point.
@param 	*molgroup molecule group structure.
@param 	point		3-valued inversion point.
@return int 				0.
**/
int 		molgroup_coor_invert(Bmolgroup* molgroup, Vector3<double> point)
{
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_PROCESS )
		cout << "Inverting coordinates:          " << point << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord[0] = 2*point[0] - atom->coord[0];
				atom->coord[1] = 2*point[1] - atom->coord[1];
				atom->coord[2] = 2*point[2] - atom->coord[2];
			}
		}
    }
	
	return 0;
}

/**
@brief 	Recombines coordinates split across periodic box boundaries.
@param 	*molgroup 		molecule group structure.
@return int 			0.

	The center of mass of each molecule is calculated to determine the
	shift needed to be applied to get the center of mass within the box.

**/
int  		molgroup_resolve_pbc(Bmolgroup* molgroup)
{
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	double		summass;
	Vector3<double>		csum, shift;
	Vector3<double>		box = molgroup->box;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Resolving periodic box coordinates:" << endl;
		cout << "Box:                            " << molgroup->box << endl << endl;
	}
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		csum = 0;
		summass = 0;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( summass == 0 ) {
					csum = atom->coord * atom->mass;
				} else {
					shift = csum / summass;
					shift += vector3_difference_PBC(atom->coord, shift, box);
					csum += shift * atom->mass;
				}
				summass += atom->mass;
			}
		}
		if ( summass ) csum /= summass;
		if ( csum[0] < 0 ) csum[0] += box[0];
		if ( csum[1] < 0 ) csum[1] += box[1];
		if ( csum[2] < 0 ) csum[2] += box[2];
		if ( csum[0] > molgroup->box[0] ) csum[0] -= box[0];
		if ( csum[1] > molgroup->box[1] ) csum[1] -= box[1];
		if ( csum[2] > molgroup->box[2] ) csum[2] -= box[2];
		if ( verbose & VERB_FULL )
			cout << "COM for molecule " << mol->id << ", mass " << summass << ": " << csum << endl;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				shift = atom->coord - csum;
				if ( 2*shift[0] > box[0] ) shift[0] = -box[0];
				else if ( 2*shift[0] < -box[0] ) shift[0] = box[0];
				else shift[0] = 0;
				if ( 2*shift[1] > box[1] ) shift[1] = -box[1];
				else if ( 2*shift[1] < -box[1] ) shift[1] = box[1];
				else shift[1] = 0;
				if ( 2*shift[2] > box[2] ) shift[2] = -box[2];
				else if ( 2*shift[2] < -box[2] ) shift[2] = box[2];
				else shift[2] = 0;
				atom->coord = atom->coord + shift;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Distributes coordinates across periodic box boundaries.
@param 	*molgroup molecule group structure.
@return int 				0.
**/
int  		molgroup_pack_in_periodic_box(Bmolgroup* molgroup)
{
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Wrapping coordinates within a periodic box:" << endl;
		cout << "Box:                            " << molgroup->box << endl << endl;
	}
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = vector3_set_PBC(atom->coord, molgroup->box);
			}
		}
	}
	
	return 0;
}

/**
@brief 	Translates a molecule within a periodic box.

	The periodic box is defined in the molecule group structure.

@param 	*molgroup molecule group structure.
@param 	shift		three-valued translation vector.
@return int 				0.
**/
int 		molgroup_coor_shift_PBC(Bmolgroup* molgroup, Vector3<double> shift)
{
	if ( shift.length2() < 1e-30 ) return 0;

	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Shifting coordinates:           " << shift << endl;
		cout << "Wrapping within box:            " << molgroup->box << endl;
	}
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord + shift;
				atom->coord = vector3_set_PBC(atom->coord, molgroup->box);
			}
		}
    }
	
	return 0;
}

/**
@brief 	Translates and rotates a molecule within a periodic box.

	The molecule group is first rotated around the origin and then shifted.
	The periodic box is defined in the molecule group structure.

@param 	*molgroup molecule group structure.
@param 	origin		3-valued origin for rotation.
@param 	mat			3x3 rotation matrix.
@param 	shift		3-valued translation vector.
@return int 				0.
**/
int 		molgroup_coor_shift_rotate_PBC(Bmolgroup* molgroup, Vector3<double> origin, 
				Matrix3 mat, Vector3<double> shift)
{
	if ( shift.length2() < 1e-30 ) return 0;

	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Rotating and translating coordinates:" << endl;
		cout << "Rotation origin:                " << origin << endl;
		cout << "Rotation matrix:" << endl;
		cout << mat << endl;
		cout << "Determinant:                    " << mat.determinant() << endl;
		cout << "Translation after rotation:     " << shift << endl;
		cout << "Wrapping within box:            " << molgroup->box << endl << endl;
	}
	
	shift += origin;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord = atom->coord - origin;
				atom->coord = mat * atom->coord;
				atom->coord = atom->coord + shift;
				atom->coord = vector3_set_PBC(atom->coord, molgroup->box);
			}
		}
    }
	
	return 0;
}

/**
@brief     Translates a set of coordinates to the center of mass.

	The function molgroup_center_of_mass is used to calculate the center of mass.
	The function molgroup_coor_shift is used to shift the coordinates.

@param 	*molgroup 	molecule group structure.
@return int 					0.
**/
int 		molgroup_shift_to_center_of_mass(Bmolgroup* molgroup)
{
	Vector3<double>		com = molgroup_center_of_mass(molgroup);
	
	com = -com;
	
	molgroup_coor_shift(molgroup, com);
	
	return 0;
}

/**
@brief     Translates a set of coordinates to the center of mass.

	The function mol_center_of_mass is used to calculate the center of mass.
	The function mol_coor_shift is used to shift the coordinates.

@param 	*mol			molecule structure.
@return int 					0.
**/
int 		mol_shift_to_center_of_mass(Bmolecule* mol)
{
	Vector3<double>		com = mol_center_of_mass(mol);
	
	com = -com;
	
	mol_coor_shift(mol, com);
	
	return 0;
}

/**
@brief     Translates a set of coordinates to a defined center-of-mass.

	The function molgroup_center_of_mass is used to calculate the center of mass.
	The function molgroup_coor_shift is used to shift the coordinates.

@param 	*molgroup 	molecule group structure.
@param 	location		desired center-of-mass coordinates.
@return int 					0.
**/
int			molgroup_place_at_coordinates(Bmolgroup* molgroup, Vector3<double> location)
{
	Vector3<double>		com = molgroup_center_of_mass(molgroup);
	
	molgroup_coor_shift(molgroup, (location - com));
	
	return 0;
}


