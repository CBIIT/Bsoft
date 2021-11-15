/**
@file	rwmd.cpp
@brief	Library routines to read and write molecular dynamics parameters in STAR format
@author Bernard Heymann
@date	Created: 20030919
@date	Modified: 20080924
**/

#include "rwmd.h"
#include "rwstar.h"
#include "mol_tags.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Initializes a molecular dynamics structure.
@return Bmd*				molecular dynamics structure.

	The coefficients for the Lennard-Jones Van der Waals interaction are
	defined as:
		VdWcoeff1 = 1/12
		VdWcoeff2 = 1/6

**/
Bmd*		md_init()
{
	Bmd*			md = new Bmd;
	memset(md, 0, sizeof(Bmd));
	
	md->cutoff = 5;
	md->sepdist = 4;
	md->VdWcoeff1 = 1.0/12.0;
	md->VdWcoeff2 = 1.0/6.0;
	
	return md;
}

/**
@brief 	Initializes a molecular dynamics structure.
@return Bmd*				molecular dynamics structure.

	The coefficients for the Lennard-Jones Van der Waals interaction are
	defined as:
		VdWcoeff1 = 1/12
		VdWcoeff2 = 1/6

**/
Bmd*		md_init_with_types()
{
	Bmd*			md = md_init();
	
	char			type[8] = "A";
	
	add_item((char **) &md->atom, sizeof(Batomtype));
	strncpy(md->atom->name, type, 8);
	md->atom->el[0] = type[0];
	
	add_item((char **) &md->bond, sizeof(Bbondtype));
	strncpy(md->bond->type1, type, 8);
	strncpy(md->bond->type2, type, 8);
	md->bond->covlength = 1;
/*
	add_item((char **) &md->angle, sizeof(Bangletype));
	strncpy(md->angle->type1, type, 8);
	strncpy(md->angle->type2, type, 8);
	strncpy(md->angle->type3, type, 8);
	md->angle->angle = M_PI*120.0/180.0;

	md->angle = new Bangletype(type, type, type);
	md->angle->angle(M_PI*120.0/180.0);
*/
	return md;
}

/**
@brief 	Reading molecular dynamics parameters from STAR files.
@param 	&filename		file name (or comma-delimited list).
@return Bmd*			molecular dynamics structure.
**/
Bmd*		read_md_parameters(Bstring& filename)
{
	int				i;
	Bstar*          star = init_star();
	star->line_length = 160;                // Set the output line length

	if ( read_star(filename.c_str(), star) != 0 ) return NULL;

	Bstar_block*	block = star->block;
	Bstar_item*		item = NULL;
	Bstring*		data;
	
	if ( verbose & VERB_FULL )
		cout << "Starting to convert molecular dynamics parameters from the STAR database" << endl;
	
	Bmd*			md = md_init();
	Batomtype*		atom = NULL;
	Bbondtype*		bond = NULL;
//	Bangletype*		angle = NULL;
	
	int				natom = item_get_number(star, ATOM_TYPE_SYMBOL);
	int				nbond = item_get_number(star, BOND_TYPE_SYMBOL1);
	int				nangle = item_get_number(star, ANGLE_TYPE_SYMBOL1);
	
	if ( verbose & VERB_FULL )
		cout << nbond << " bond types and " << nangle << " angle types found" << endl << endl;
	
	if ( natom ) {
		item = item_find(block, ATOM_TYPE_SYMBOL);
		if ( item )
			for ( data=item->data; data; data=data->next )
				atom = atom_type_add(&md->atom, data->c_str());
		item = item_find(block, ATOM_TYPE_NUMBER);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->z = data->integer();
		item = item_find(block, ATOM_TYPE_MASS);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->mass = data->real();
		item = item_find(block, ATOM_TYPE_OXIDATION);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->oxid = data->real();
		item = item_find(block, ATOM_TYPE_RADIUS_BOND);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->bond = data->real();
		item = item_find(block, ATOM_TYPE_RADIUS_VDW);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->vdw = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_A1);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfa[0] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_A2);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfa[1] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_A3);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfa[2] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_A4);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfa[3] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_A5);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfa[4] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_B1);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfb[0] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_B2);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfb[1] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_B3);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfb[2] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_B4);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfb[3] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_B5);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfb[4] = data->real();
		item = item_find(block, ATOM_TYPE_SCAT_C);
		if ( item )
			for ( data=item->data, atom=md->atom; data && atom; data=data->next, atom=atom->next )
				atom->sfc = data->real();
	}

	if ( nbond ) {
		block = block_find_with_tag(star, BOND_TYPE_SYMBOL1);
		for ( i=0; i<nbond; i++ ) {
			if ( bond ) bond = (Bbondtype *) add_item((char **) &bond, sizeof(Bbondtype));
			else bond = (Bbondtype *) add_item((char **) &md->bond, sizeof(Bbondtype));
		}
		item = item_find(block, BOND_TYPE_SYMBOL1);
		if ( item )
			for ( data=item->data, bond=md->bond; data && bond; data=data->next, bond=bond->next )
				strncpy(bond->type1, data->c_str(), 8);
		item = item_find(block, BOND_TYPE_SYMBOL2);
		if ( item )
			for ( data=item->data, bond=md->bond; data && bond; data=data->next, bond=bond->next )
				strncpy(bond->type2, data->c_str(), 8);
		item = item_find(block, BOND_TYPE_LENGTH);
		if ( item )
			for ( data=item->data, bond=md->bond; data && bond; data=data->next, bond=bond->next )
				bond->covlength = data->real();
		item = item_find(block, BOND_TYPE_VDWDIST);
		if ( item )
			for ( data=item->data, bond=md->bond; data && bond; data=data->next, bond=bond->next )
				bond->vdwdist = data->real();
	}
/*
	if ( nangle ) {
		block = block_find_with_tag(star, ANGLE_TYPE_SYMBOL1);
		for ( i=0; i<nangle; i++ ) {
			if ( angle ) angle = (Bangletype *) add_item((char **) &angle, sizeof(Bangletype));
			else angle = (Bangletype *) add_item((char **) &md->angle, sizeof(Bangletype));
		}
		item = item_find(block, ANGLE_TYPE_SYMBOL1);
		if ( item )
			for ( data=item->data, angle=md->angle; data && angle; data=data->next, angle=angle->next )
				strncpy(angle->type1, data->c_str(), 8);
		item = item_find(block, ANGLE_TYPE_SYMBOL2);
		if ( item )
			for ( data=item->data, angle=md->angle; data && angle; data=data->next, angle=angle->next )
				strncpy(angle->type2, data->c_str(), 8);
		item = item_find(block, ANGLE_TYPE_SYMBOL3);
		if ( item )
			for ( data=item->data, angle=md->angle; data && angle; data=data->next, angle=angle->next )
				strncpy(angle->type3, data->c_str(), 8);
		item = item_find(block, ANGLE_TYPE_ANGLE);
		if ( item )
			for ( data=item->data, angle=md->angle; data && angle; data=data->next, angle=angle->next )
				angle->angle = data->real()*M_PI/180.0;
	}
*/
	kill_star(star);

	return md;
}

/**
@brief 	Writing molecular dynamics parameters to a STAR file.
@param 	&filename		file name (or comma-delimited list).
@param 	*md				molecular dynamics parameters structure.
@return int				error code (<0 means failure).
**/
int			write_md_parameters(Bstring& filename, Bmd* md)
{
	int				natom, nbond, nangle, loop, err(0);
	Bstring			comment("Bsoft molecular dynamics parameter file\n");
	Bstar*          star = init_star();
	Bstar_block*	block = NULL;
	star->line_length = 160;                // Set the output line length
	star->comment += comment;

	Batomtype*		atom;
	Bbondtype*		bond;
//	Bangletype*		angle;
	
	for ( natom=0, atom=md->atom; atom; atom=atom->next, natom++ ) ;
	for ( nbond=0, bond=md->bond; bond; bond=bond->next ) nbond++;
//	for ( nangle=0, angle=md->angle; angle; angle=angle->next ) nangle++;
	
	if ( natom ) {
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		block->tag = "atom_types";
		item_put_list(block, ATOM_TYPE_SYMBOL, (char *) md->atom, 
			(char *)&md->atom->name - (char *)md->atom, "%s");
		item_put_list(block, ATOM_TYPE_NUMBER, (char *) md->atom, 
			(char *)&md->atom->z - (char *)md->atom, "%4d");
		item_put_list(block, ATOM_TYPE_MASS, (char *) md->atom, 
			(char *)&md->atom->mass - (char *)md->atom, "%7.3f");
		item_put_list(block, ATOM_TYPE_OXIDATION, (char *) md->atom, 
			(char *)&md->atom->oxid - (char *)md->atom, "%7.3f");
		item_put_list(block, ATOM_TYPE_RADIUS_BOND, (char *) md->atom, 
			(char *)&md->atom->bond - (char *)md->atom, "%7.3f");
		item_put_list(block, ATOM_TYPE_RADIUS_VDW, (char *) md->atom, 
			(char *)&md->atom->vdw - (char *)md->atom, "%7.3f");
		item_put_list(block, ATOM_TYPE_SCAT_A1, (char *) md->atom, 
			(char *)&md->atom->sfa[0] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_A2, (char *) md->atom, 
			(char *)&md->atom->sfa[1] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_A3, (char *) md->atom, 
			(char *)&md->atom->sfa[2] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_A4, (char *) md->atom, 
			(char *)&md->atom->sfa[3] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_A5, (char *) md->atom, 
			(char *)&md->atom->sfa[4] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_B1, (char *) md->atom, 
			(char *)&md->atom->sfb[0] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_B2, (char *) md->atom, 
			(char *)&md->atom->sfb[1] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_B3, (char *) md->atom, 
			(char *)&md->atom->sfb[2] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_B4, (char *) md->atom, 
			(char *)&md->atom->sfb[3] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_B5, (char *) md->atom, 
			(char *)&md->atom->sfb[4] - (char *)md->atom, "%7.4f");
		item_put_list(block, ATOM_TYPE_SCAT_C, (char *) md->atom, 
			(char *)&md->atom->sfc - (char *)md->atom, "%7.4f");
		loop = item_index(block, ATOM_TYPE_SYMBOL);
		loop_set_identifier(block, loop, 1, "atom_type");
	}
	
	if ( nbond ) {
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		block->tag = "bond_types";
		item_put_list(block, BOND_TYPE_SYMBOL1, (char *) md->bond, 
			(char *)&md->bond->type1 - (char *)md->bond, "%s");
		item_put_list(block, BOND_TYPE_SYMBOL2, (char *) md->bond, 
			(char *)&md->bond->type2 - (char *)md->bond, "%s");
		item_put_list(block, BOND_TYPE_LENGTH, (char *) md->bond, 
			(char *)&md->bond->covlength - (char *)md->bond, "%7.3f");
		item_put_list(block, BOND_TYPE_VDWDIST, (char *) md->bond, 
			(char *)&md->bond->vdwdist - (char *)md->bond, "%7.3f");
		loop = item_index(block, BOND_TYPE_SYMBOL1);
		loop_set_identifier(block, loop, 1, "bond_type");
	}
/*
	if ( nangle ) {
		block = (Bstar_block *) add_item((char **)&star->block, sizeof(Bstar_block));
		block->tag = "angle_types";
		item_put_list(block, ANGLE_TYPE_SYMBOL1, (char *) md->angle, 
			(char *)&md->angle->type1 - (char *)md->angle, "%s");
		item_put_list(block, ANGLE_TYPE_SYMBOL2, (char *) md->angle, 
			(char *)&md->angle->type2 - (char *)md->angle, "%s");
		item_put_list(block, ANGLE_TYPE_SYMBOL3, (char *) md->angle, 
			(char *)&md->angle->type3 - (char *)md->angle, "%s");
		item_put_angle_list(block, ANGLE_TYPE_ANGLE, (char *) md->angle, 
			(char *)&md->angle->angle - (char *)md->angle, "%7.2f");
		loop = item_index(block, ANGLE_TYPE_SYMBOL1);
		loop_set_identifier(block, loop, 1, "angle_type");
	}
*/
	err = write_star(filename.c_str(), star);
	
	kill_star(star);

	return err;
}

/*
@brief 	Adds a bond type to a linked list.
@param 	**bond			pointer to any bond type in the list.
@param 	*symbol1		atom1 of bond.
@param 	*symbol2		atom2 of bond.
@param 	covlength		reference covalent bond length.
@param 	vdwdist			reference Van der Waals distance.
@return Bbondtype* 		new bond type.

	The function allocates memory for a new bond type structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.

**/
Bbondtype*	bond_type_add(Bbondtype** bond, char* symbol1, char* symbol2, float covlength, float vdwdist)
{
	Bbondtype* 		this_bond = *bond;
	Bbondtype* 		new_bond = new Bbondtype;
	memset(new_bond, 0, sizeof(Bbondtype));
	
	strncpy(new_bond->type1, symbol1, 8);
	strncpy(new_bond->type2, symbol2, 8);
	new_bond->covlength = covlength;
	new_bond->vdwdist = vdwdist;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG bond_type_add: covlength=" << covlength << " vdwdist=" << vdwdist << endl;
	
	if ( !this_bond )
		*bond = new_bond;
	else {
		while ( this_bond->next ) this_bond = this_bond->next;
		this_bond->next = new_bond;
	}
	
	return new_bond;
}

/*
@brief	Adds an angle type to a linked list.
@param	**angle	pointer to any angle type in the list.
@param	*symbol1		atom1 of angle.
@param	*symbol2		atom2 of angle (central atom).
@param	*symbol3		atom3 of angle.
@param	float a			reference angle.
@return	Bangletype* 	new angle type.

	The function allocates memory for a new angle type structure.
	If the content of the pointer is null, the new structure is
	the first in the list. Otherwise, the end of the list is found
	and the new structure added to it.
**/
/*Bangletype*	angle_type_add(Bangletype** angle, char* symbol1, char* symbol2, char* symbol3, float a)
{
	Bangletype* 	this_angle = *angle;
	Bangletype* 	new_angle = new Bangletype;
	memset(new_angle, 0, sizeof(Bangletype));
	
	strncpy(new_angle->type1, symbol1, 8);
	strncpy(new_angle->type2, symbol2, 8);
	strncpy(new_angle->type3, symbol3, 8);
	new_angle->angle = a;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG angle_type_add: angle = " << a << endl;
	
	if ( !this_angle )
		*angle = new_angle;
	else {
		while ( this_angle->next ) this_angle = this_angle->next;
		this_angle->next = new_angle;
	}
	
	return new_angle;
}
*/
/**
@brief 	Deallocates a molecular dynamics structure.
@param 	*md				molecular dynamics structure.
@return int				0.

	All linked lists within the structure is also deallocated.

**/
int			md_kill(Bmd* md)
{
	if ( !md ) return 0;
	
	kill_list((char *) md->bond, sizeof(Bbondtype));
	
//	kill_list((char *) md->angle, sizeof(Bangletype));

	delete md;
	
	return 0;
}
