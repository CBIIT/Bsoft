/**
@file	mol_param.cpp
@brief	Functions to extract parameters from coordinate files
@author Bernard Heymann
@date	Created: 20050304
@date	Modified: 20090226
**/

#include "rwmolecule.h"
#include "rwmd.h"
#include "mol_bonds.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Obtains atomic types from molecule groups.
@param 	*molgroup	molecule group.
@param 	elements		flag to take atom types from element records.
@return Batomtype*			atom property structure, NULL on failure.
**/
Batomtype*	molgroup_get_atom_types(Bmolgroup* molgroup, int elements)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Batomtype*		at = NULL;
	Batomtype*		atlist = NULL;
	char			string[32];
	
	if ( verbose )
		cout << "Atom types:" << endl;

	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( elements ) strcpy(string, atom->el);
				else strcpy(string, atom->type);
				if ( atlist )
					for ( at = atlist; at && strcmp(at->name, string); at = at->next ) ;
				if ( !at ) {
					at = atom_type_add(&atlist, string);
					at->z = 1;
					at->el[0] = string[0];
					at->mass = 1;
				}
				at->number++;
			}
		}
	}
	
	if ( verbose )
		cout << "Type\tNumber" << endl;
	for ( at = atlist; at; at = at->next )
		cout << at->name << tab << at->number << endl;
	cout << endl;
	
	return atlist;
}

/**
@brief 	Calculates atomic parameters from molecule groups.
@param 	*molgroup	molecule group.
@param 	elements		flag to take atom types from element records.
@param 	show_bond		flag to show bonds.
@param 	show_angle		flag to show angles.
@return Bmd*				molecular dynamics object, NULL on failure.
**/
Bmd*		md_calculate_parameters(Bmolgroup* molgroup, int elements, int show_bond, int show_angle)
{
	Bmd*			md = md_init();
	
	Bbond*			bond = NULL;
//	Bangle*			angle = NULL;
	Bbondtype*		bt = NULL;
//	Bangletype*		at = NULL;
	
	double			d;

	if ( !md->atom ) md->atom = molgroup_get_atom_types(molgroup, elements);
	
	if ( !molgroup->bond ) molgroup->bond = mol_bond_list_generate(molgroup, 1.8, 0);
	
//	if ( !molgroup->angle ) md_generate_angle_list(molgroup, md);
	
	if ( verbose )
		cout << "Bond types:" << endl;

	for ( bond = molgroup->bond; bond; bond = bond->next ) {
		bt = md_find_bond_type(bond->atom1, bond->atom2, md->bond);
		if ( !bt ) bt = (Bbondtype *) add_item((char **) &md->bond, sizeof(Bbondtype));
		if ( elements ) {
			strncpy(bt->type1, bond->atom1->el, 2);
			strncpy(bt->type2, bond->atom2->el, 2);
		} else {
			strncpy(bt->type1, bond->atom1->type, 8);
			strncpy(bt->type2, bond->atom2->type, 8);
		}
		d = bond->atom1->coord.distance(bond->atom2->coord);
		if ( show_bond ) cout << bond->atom1->type << tab << bond->atom2->type << tab << d << endl;
		bt->covlength += d;
		bt->std += d*d;
		bt->number++;
	}
	
	if ( verbose )
		cout << "Type1\tType2\tNumber\tLength\tSTD" << endl;
	for ( bt = md->bond; bt; bt = bt->next ) if ( bt->number ) {
		bt->covlength /= bt->number;
		bt->std = bt->std/bt->number - bt->covlength*bt->covlength;
		if ( bt->std > 0 ) bt->std = sqrt(bt->std);
		else bt->std = 0;
		if ( verbose )
			cout << bt->type1 << tab << bt->type2 << tab << bt->number << tab << bt->covlength << tab << bt->std << endl;
	}

	if ( verbose )
		cout << endl << "Angle types:" << endl;
/*
	for ( angle = molgroup->angle; angle; angle = angle->next ) {
		at = md_find_angle_type(angle->atom1, angle->atom2, angle->atom3, md->angle);
		if ( show_angle ) cout << angle->atom1->type << tab << angle->atom2->type << tab << angle->atom3->type << tab << angle->a*180.0/M_PI << endl;
		if ( !at ) {
			at = (Bangletype *) add_item((char **) &md->angle, sizeof(Bangletype));
			if ( elements ) {
				strncpy(at->type1, angle->atom1->el, 2);
				strncpy(at->type2, angle->atom2->el, 2);
				strncpy(at->type3, angle->atom3->el, 2);
			} else {
				strncpy(at->type1, angle->atom1->type, 8);
				strncpy(at->type2, angle->atom2->type, 8);
				strncpy(at->type3, angle->atom3->type, 8);
			}
		}
		at->angle += angle->a;
		at->std += angle->a*angle->a;
		at->number++;
	}
	
	for ( at = md->angle; at; at = at->next ) {
		at->angle /= at->number;
		at->std = at->std/at->number - at->angle*at->angle;
		if ( at->std > 0 ) at->std = sqrt(at->std);
		else at->std = 0;
	}

	if ( verbose ) {
		cout << "Type1\tType2\tType3\tNumber\tAngle\tSTD" << endl;
		for ( at = md->angle; at; at = at->next )
			cout << at->type1 << tab << at->type2 << tab << at->type3 << tab << 
				at->number << tab << at->angle*180.0/M_PI << tab << at->std*180.0/M_PI << endl;
		cout << endl;
	}
*/
	return md;
}
