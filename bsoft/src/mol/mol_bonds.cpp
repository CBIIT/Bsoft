/**
@file	mol_bonds.cpp
@brief	Functions for molecular dynamics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20180226
**/

#include "mol_bonds.h"
#include "rwmd.h"
#include "rwmolecule.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates a distance-based bond list.
@param 	*molgroup 	molecule group structure.
@param 	*md			global molecular dynamics structure.
@return Bbond* 		new bond list.

	This function assumes very little and defines bonds purely on distance.
	This means that the bond distances must already have been defined well.
	The bond length is defined in the molecular dynamics structure.
	If the molecule group already has a bond list, no new bonds are generated. 
**/
Bbond*		md_generate_bond_list(Bmolgroup* molgroup, Bmd* md)
{
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Bond list already defined!" << endl;
		return molgroup->bond;
	}
	
	int 		i, j, natom, nbond;
	double		dist, bondlength;
	Bbond*		bond = NULL;
	Bbond*		bondlist = NULL;
	
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom, *atom2;
	
    for ( natom=0, mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) natom++;
			
	Batom**		atomlist = new Batom*[natom];

    for ( i=0, mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) {
				atomlist[i++] = atom;
//				atom->sel = 0;
			}
			
	if ( verbose ) {
		cout << "Generating a bond list for" << natom << " atoms" << endl;
		if ( md->wrap ) cout << "With wrapping" << endl;
	}
	
	for ( i=0, nbond=0; i<natom-1; i++ ) {
    	atom = atomlist[i];
		for ( j=i+1; j<natom; j++ ) {
			atom2 = atomlist[j];
			bondlength = md_find_bond_length(atom, atom2, md->bond);
			if ( md->wrap ) 
				dist = (vector3_difference_PBC(atom->coord, atom2->coord, molgroup->box)).length();
			else
				dist = atom->coord.distance(atom2->coord);
			if ( dist > 0.5*bondlength && dist < 1.1*bondlength ) {
				if ( bond ) bond = bond_add(&bond, atom, atom2, bondlength, 1);
				else bond = bond_add(&bondlist, atom, atom2, bondlength, 1);
				nbond++;
		    }
		}
    }
	
	if ( verbose & VERB_FULL )
		md_show_bonds(molgroup);
	
	if ( verbose )
		cout << "Number of bonds generated:      " << nbond << endl << endl;
	
	delete[] atomlist;
	
	molgroup->bond = bondlist;
	
	return bondlist;
}

/**
@brief 	Generates an intramolecular distance-based bond list.
@param 	*molgroup 	molecule group structure.
@param 	*md			global molecular dynamics structure.
@return Bbond* 		new bond list.

	This function defines bonds on distance and within molecules.
	This means that the bond distances must already have been defined well.
	The bond length is defined in the molecular dynamics structure.
	If the molecule group already has a bond list, no new bonds are generated. 
**/
Bbond*		md_generate_molecular_bond_list(Bmolgroup* molgroup, Bmd* md)
{
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Bond list already defined!" << endl;
		return molgroup->bond;
	}
	
	long 		natom = 0, nbond = 0;
	double		dist, bondlength;
	Bbond*		bond = NULL;
	Bbond*		bondlist = NULL;
	
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom, *atom2;
	
    for ( natom=0, mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom->next; atom = atom->next ) {
				natom++;
				for ( atom2 = atom->next; atom2; atom2 = atom2->next ) {
					bondlength = md_find_bond_length(atom, atom2, md->bond);
					if ( md->wrap )
						dist = (vector3_difference_PBC(atom->coord, atom2->coord, molgroup->box)).length();
					else
						dist = atom->coord.distance(atom2->coord);
					if ( dist > 0.1*bondlength && dist < 1.1*bondlength ) {
						if ( bond ) bond = bond_add(&bond, atom, atom2, bondlength, 1);
						else bond = bond_add(&bondlist, atom, atom2, bondlength, 1);
						nbond++;
					}
				}
			}
		}
	}
			
	if ( verbose & VERB_FULL )
		md_show_bonds(molgroup);
	
	if ( verbose )
		cout << "Number of bonds generated:      " << nbond << endl << endl;
	
	molgroup->bond = bondlist;
	
	return bondlist;
}

/*
@brief 	Finds the top shortest bonds.
@param 	*bondlist	list of bonds.
@param 	number		number of shortest bonds to find.
@return Bbond* 		short bond list.

	This function finds the top number of shortests bonds with the aim
	of assigning bonds based valence. The current bond length must be
	given by the bond strength variable, k, in the bond structure.
**/
Bbond*		md_bondlist_get_top(Bbond* bondlist, int number)
{
	if ( !bondlist ) return NULL;
	
	int			i, imax = 0;
	double		max = 0;
	double*		val = new double[number];
	Bbond**		bondarray = new Bbond*[number];
	Bbond*		newbond = NULL, *bond;
	
	for ( i=0; i<number; i++ ) {
		val[i] = 1e37;
		bondarray[i] = NULL;
	}
	
	for ( bond = bondlist; bond; bond = bond->next ) {
		for ( i=imax=0, max=0; i<number; i++ ) if ( max < val[i] ) {
			max = val[i];
			imax = i;
		}
		if ( val[imax] > bond->k ) {
			val[imax] = bond->k;
			bondarray[imax] = bond;
		}
	}
	
	for ( i=0; i<number; i++ ) if ( bondarray[i] )
		bond_add(&newbond, bondarray[i]->atom1, bondarray[i]->atom2, bondarray[i]->l, 1);
	
	delete[] val;
	delete[] bondarray;
	
	return newbond;
}

/**
@brief 	Generates a bond list based on a single valence and the shortest current bond lengths.
@param 	*molgroup 	molecule group structure.
@param 	*md			global molecular dynamics structure.
@param 	valence		atom valence.
@return Bbond* 		new bond list.

	This function attempts to assign bonds based on atomic valence by
	finding the top number of shortests bonds. The current bond length is
	encoded in the bond strength variable, k, in the bond structure.
**/
Bbond*		md_generate_bond_list_with_valence(Bmolgroup* molgroup, Bmd* md, int valence)
{
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Bond list already defined!" << endl;
		return molgroup->bond;
	}
	
	double			cutoff = md->cutoff;
	
	int				i, ii, x, y, z, xx, yy, zz, ix, iy, iz;
	Vector3<double>	box = molgroup->box;
	Vector3<double>	sampling(cutoff, cutoff, cutoff);
	Vector3<int>	size((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/size[i] + 0.001;
	long			boxsize = (long) size.volume();
	double			dist, bondlength;
	Batom			*atom, *atom2;
	Latom			*latom, *latom2;
	Bbond*			tbondlist = NULL;
	Bbond*			tbond = NULL;
	Bbond*			bond = NULL;
	
	if ( verbose )
		cout << "Sampling intervals for cutoff " << cutoff << ": " << sampling << endl;
	
	Latom**			alist = molgroup_atom_mesh_lists(molgroup, size, sampling);
	
	for ( z=0; z<size[2]; z++ ) {
		for ( y=0; y<size[1]; y++ ) {
			for ( x=0; x<size[0]; x++ ) {
				i = (z*size[1] + y)*size[0] + x;
				tbond = tbondlist = NULL;
				for ( latom = alist[i]; latom; latom = latom->next ) {
					atom = latom->atom;
					for ( zz=z-1; zz<=z+1; zz++ ) {
						iz = zz;
						if ( iz < 0 ) iz += size[2];
						if ( iz >= size[2]) iz -= size[2];
						for ( yy=y-1; yy<=y+1; yy++ ) {
							iy = yy;
							if ( iy < 0 ) iy += size[1];
							if ( iy >= size[1] ) iy -= size[1];
							for ( xx=x-1; xx<=x+1; xx++ ) {
								ix = xx;
								if ( ix < 0 ) ix += size[0];
								if ( ix >= size[0] ) ix -= size[0];
								ii = (iz*size[1] + iy)*size[0] + ix;
								for ( latom2 = alist[ii]; latom2; latom2 = latom2->next ) {
									atom2 = latom2->atom;
									bondlength = md_find_bond_length(atom, atom2, md->bond);
									if ( md->wrap )
										dist = (vector3_difference_PBC(atom->coord, atom2->coord, box)).length();
									else
										dist = atom->coord.distance(atom2->coord);
									if ( dist > 0.1*bondlength && dist < 1.3*bondlength ) {
										if ( tbond ) tbond = bond_add(&tbond, atom, atom2, bondlength, dist);
										else tbond = bond_add(&tbondlist, atom, atom2, bondlength, dist);
									}
								}
							}
						}
					}
				}
				tbond = md_bondlist_get_top(tbondlist, valence);
				if ( bond ) bond->next = tbond;
				else bond = molgroup->bond = tbond;
				if ( bond ) for ( ; bond->next; bond = bond->next ) ;
				bond_kill(tbondlist);
			}
		}
	}
	
	for ( i=0; i<boxsize; i++ ) kill_list((char *) alist[i], sizeof(Latom));
	delete[] alist;
	
	if ( verbose & VERB_FULL )
		md_show_bonds(molgroup);
	
	return molgroup->bond;
}

/**
@brief 	Assigns covalent bond lengths to bonds.
@param 	*bondlist	list of bonds.
@param 	*bondtype	list of bond types.
@return double		fraction of bond types found.

	The bond lengths are taken from a bond type list.
**/
double		md_bond_list_set_parameters(Bbond* bondlist, Bbondtype* bondtype)
{
	long			n(0), nf(0);
	double			f(0);
	Bbond*			bond;
	Bbondtype*		bt;
	
	for ( bond = bondlist; bond; bond = bond->next, ++n ) {
		bt = md_find_bond_type(bond->atom1, bond->atom2, bondtype);
		if ( bt ) {
			bond->l = bt->covlength;
			nf++;
		}
	}
	
	f = nf*1.0L/n;
	
	if ( verbose )
		cout << "Fraction of bond types found:   " << f << endl;
	
	return f;
}

/**
@brief 	Shows the bonds and bond lengths.
@param 	*molgroup	molecule group.
@return int			number of bonds.

	Uses the bond list defined for the molecule group.
**/
int			md_show_bonds(Bmolgroup* molgroup)
{
	int				nbond;
	Bbond*			bond;

	cout << "Bond#\tAtom1\tAtom2\tLength\tStrength" << endl;
	for ( nbond=0, bond = molgroup->bond; bond; bond = bond->next, nbond++ )
		cout << nbond+1 << tab << bond->atom1->num << tab << 
				bond->atom2->num << tab << bond->l << tab << bond->k << endl;
	cout << endl;
	
	return nbond;
}

/**
@brief 	Shows the angles and angular sizes.
@param 	*molgroup	molecule group.
@return int			number of angles.

	Uses the angle list defined for the molecule group.
**/
int			md_show_angles(Bmolgroup* molgroup)
{
	int				nang;
	Bangle*			angle = NULL;

	cout << "Angle#\tAtom1\tAtom2\tAtom3\tAngle\tStrength" << endl;
	for ( nang=0, angle = molgroup->angle; angle; angle = angle->next, nang++ )
		cout << nang+1 << tab << angle->atom1->num << tab << angle->atom2->num << tab << 
			angle->atom3->num << tab << angle->a*180.0/M_PI << tab << angle->k << endl;
	cout << endl;
	
	return nang;
}

/**
@brief 	Shows the number of bonds and the valency per atom.
@param 	*molgroup	molecule group.
@return int			number of bonds.

	Uses the bond list defined for the molecule group.
**/
int			md_show_bond_stats(Bmolgroup* molgroup)
{
	int			i, nbond;
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	Bbond*		bond;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for ( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				atom->sel = 0;

	for ( nbond=0, bond = molgroup->bond; bond; bond = bond->next, nbond++ ) {
		bond->atom1->sel++;
		bond->atom2->sel++;
	}
	
	cout << "Number of bonds:  " << nbond << endl;
	
	int			hist[100];
	for ( i=0; i<100; i++ ) hist[i] = 0;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for ( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				hist[atom->sel]++;

	cout << "Val\tAtoms\tBonds" << endl;
	for ( i=0; i<10; i++ ) cout << i << tab << hist[i] << tab << i*hist[i] << endl;
	
	return nbond;
}

/*
@brief 	Generates an angle list.
@param 	*molgroup 	molecule group structure.
@param 	*md			global molecular dynamics structure (can be NULL).
@return Bangle* 	new angle list.

	This function uses the bond list for the molecule group to find bonds
	to the same atom and set an angle structure with the corresponding angle.
	The angle is defined in the molecular dynamics structure.
**/
/*Bangle*		md_generate_angle_list(Bmolgroup* molgroup, Bmd* md)
{
	Bbond*			bondlist = molgroup->bond;
	if ( !bondlist ) return NULL;
	
	Batom			*atom1, *atom2, *atom3;
	Bbond			*bond1, *bond2;
	Bangle*			anglelist = NULL;
	Bangle*			angle = NULL;
	double			set_angle = 0;
	int				nangle = 0, found;
	Vector3<double>	v1, v2;
	
	for ( bond1=bondlist; bond1->next; bond1=bond1->next ) {
		for ( bond2=bond1->next; bond2; bond2=bond2->next ) {
			found = 0;
			if ( bond1->atom1 == bond2->atom1 ) {
				atom1 = bond1->atom2;
				atom2 = bond1->atom1;
				atom3 = bond2->atom2;
				found = 1;
			} else if ( bond1->atom1 == bond2->atom2 ) {
				atom1 = bond1->atom2;
				atom2 = bond1->atom1;
				atom3 = bond2->atom1;
				found = 1;
			} else if ( bond1->atom2 == bond2->atom1 ) {
				atom1 = bond1->atom1;
				atom2 = bond1->atom2;
				atom3 = bond2->atom2;
				found = 1;
			} else if ( bond1->atom2 == bond2->atom2 ) {
				atom1 = bond1->atom1;
				atom2 = bond1->atom2;
				atom3 = bond2->atom1;
				found = 1;
			}
			if ( found ) {
				if ( md ) set_angle = md_find_angle(atom1, atom2, atom3, md->angle);
				else {
					v1 = atom1->coord - atom2->coord;
					v2 = atom3->coord - atom2->coord;
					set_angle = v1.angle(v2);
				}
				angle = angle_add(&angle, atom1, atom2, atom3, set_angle, 1);
				if ( !anglelist ) anglelist = angle;
				nangle ++;
			}
		}
	}
	
	molgroup->angle = anglelist;
	
	if ( verbose & VERB_FULL )
		md_show_angles(molgroup);
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of angles generated:     " << nangle << endl << endl;
	
	return anglelist;
}
*/
/*
@brief 	Assigns angles from angle types.
@param 	*anglelist	list of angles.
@param 	*angletype	list of angle types.
@return int			0.

	The angles are taken from an angle type list.
**/
/*int			md_angle_list_set_parameters(Bangle* anglelist, Bangletype* angletype)
{
	Bangle*			angle;
	Bangletype*		at;
	
	for ( angle = anglelist; angle; angle = angle->next ) {
		at = md_find_angle_type(angle->atom1, angle->atom2, angle->atom3, angletype);
		if ( at ) angle->a = at->angle;
		else if ( verbose & VERB_FULL )
			cout << "Angle type not found: " << angle->atom1->type << tab << angle->atom2->type << tab << angle->atom3->type << endl;
	}
	
	return 0;
}
*/
/**
@brief 	Finds the reference covalent bond length.
@param 	*atom1		first atom.
@param 	*atom2		second atom.
@param 	*bondtype	bond type list.
@return double		bond length.

	The bond type is found in the bond type list, and the bond length is
	set from the bond type structure.
**/
double		md_find_bond_length(Batom* atom1, Batom* atom2, Bbondtype* bondtype)
{
	double			bondlength = 0;
	Bbondtype*		bt;
	
	for ( bt = bondtype; bt && bondlength < 0.01; bt = bt->next ) {
		if ( strstr(atom1->type, bt->type1) && strstr(atom2->type, bt->type2) )
			bondlength = bt->covlength;
		else if ( strstr(atom1->type, bt->type2) && strstr(atom2->type, bt->type1) )
			bondlength = bt->covlength;
		else if ( strcmp(atom1->el, bt->type1) == 0 && strcmp(atom2->el, bt->type2) == 0 )
			bondlength = bt->covlength;
		else if ( strcmp(atom1->el, bt->type2) == 0 && strcmp(atom2->el, bt->type1) == 0 )
			bondlength = bt->covlength;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG md_find_bond_length: atom1=" << atom1->type << 
			" atom2=" << atom2->type << " bondlength=" << bondlength << endl;
	
	return bondlength;
}

/**
@brief 	Finds the bond type.
@param 	*atom1		first atom.
@param 	*atom2		second atom.
@param 	*bondtype	bond type list.
@return Bbondtype*	bond type structure.

	The bond type is found in the bond type list.
**/
Bbondtype*	md_find_bond_type(Batom* atom1, Batom* atom2, Bbondtype* bondtype)
{
	int				l1, l2;
	Bbondtype*		bt;
	Bbondtype*		bt_found = NULL;
	
	for ( bt = bondtype; bt && !bt_found; bt = bt->next ) {
		l1 = strlen(bt->type1);
		l2 = strlen(bt->type2);
		if ( strncmp(atom1->type, bt->type1, l1) == 0 && strncmp(atom2->type, bt->type2, l2) == 0 )
			bt_found = bt;
		else if ( strncmp(atom1->type, bt->type2, l2) == 0 && strncmp(atom2->type, bt->type1, l1) == 0 )
			bt_found = bt;
		else if ( strcmp(atom1->el, bt->type1) == 0 && strcmp(atom2->el, bt->type2) == 0 )
			bt_found = bt;
		else if ( strcmp(atom1->el, bt->type2) == 0 && strcmp(atom2->el, bt->type1) == 0 )
			bt_found = bt;
	}
	
	if ( !bt_found && ( verbose & VERB_FULL ) )
		cout << "Bond type not found: atom1=" << atom1->type << " atom2=" << atom2->type << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG md_find_bond_type: atom1=" << atom1->type << 
			" atom2=" << atom2->type << " length=" << bt_found->covlength << endl;
	
	return bt_found;
}

/**
@brief 	Shows the bond type composition.
@param 	*molgroup		molecule group.
@param 	*bondtype		bond type list.
@return int				0.

	For each bond, the bond type is found in the bond type list.
	The number variable for each bond type is used for counting.
**/
int			md_show_bond_types(Bmolgroup* molgroup, Bbondtype* bondtype)
{
	long			nbt, i;
	double			d;
	Bbond*			bond;
	Bbondtype*		bt, *bt2;
	
	cout << "Bond types:" << endl;
	
	for ( nbt=0, bt = bondtype; bt; bt = bt->next, nbt++ ) bt->number = 0;

	double*			dist = new double[nbt];
	double*			dev = new double[nbt];
	
	for ( i=0; i<nbt; i++ ) dist[i] = dev[i] = 0;
	
	for ( bond = molgroup->bond; bond; bond = bond->next) {
		bt = md_find_bond_type(bond->atom1, bond->atom2, bondtype);
		for ( i=0, bt2 = bondtype; bt2 && bt2!=bt; bt2 = bt2->next, i++ ) ;
		d = bond->atom1->coord.distance(bond->atom2->coord);
		if ( i < nbt ) {
			dist[i] += d;
			dev[i] += d*d;
			bt->number++;
		}
	}
	
	cout << "Type1\tType2\tNumber\tRefLen\tLength\tDeviation" << endl;
	for ( i=0, bt = bondtype; bt; bt = bt->next, i++ ) {
		if ( bt->number ) {
			dist[i] /= bt->number;
			dev[i] = dev[i]/bt->number - dist[i]*dist[i];
			if ( dev[i] > 0 ) dev[i] = sqrt(dev[i]);
			else dev[i] = 0;
			cout << bt->type1 << tab << bt->type2 << tab << bt->number << tab << 
				bt->covlength << tab << dist[i] << tab << dev[i] << endl;
		}
	}
	cout << endl;
	
	delete[] dist;
	delete[] dev;

	return 0;
}

/**
@brief 	Calculates the angle between three atoms.
@param 	*atom1			first atom.
@param 	*atom2			second atom.
@param 	*atom3			third atom.
@return double			angle.

	The angle associated with the second atom and between the bonds linking
	the first and third atoms is calculated.

**/
double		md_angle(Batom* atom1, Batom* atom2, Batom* atom3)
{
	Vector3<double>	d1 = atom1->coord - atom2->coord;
	Vector3<double>	d2 = atom3->coord - atom2->coord;
	
	double			angle = d1.angle(d2);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG md_angle: atom1=" << atom1->type << " atom2=" << 
			atom2->type << " atom3=" << atom3->type << " angle=" << angle << endl;
	
	return angle;
}

/*
@brief 	Finds the reference angle type between three atoms.

	The angle type is found in the angle type list, and the angle type is
	returned.

@param 	*atom1			first atom.
@param 	*atom2			second atom.
@param 	*atom3			third atom.
@param 	*angletype		angle type list.
@return Bangletype*		angle type.
**/
/*Bangletype*	md_find_angle_type(Batom* atom1, Batom* atom2, Batom* atom3, Bangletype* angletype)
{
	int				l1, l2, l3;
	Bangletype*		at;
	Bangletype*		at_found = NULL;
	Bangletype*		at_best = NULL;

	double			angle = md_angle(atom1, atom2, atom3);
	
	for ( at = angletype; at; at = at->next ) {
		l1 = strlen(at->type1);
		l2 = strlen(at->type2);
		l3 = strlen(at->type3);
		if ( strncmp(atom2->type, at->type2, l2) == 0 ) {
			if ( strncmp(atom1->type, at->type1, l1) == 0 && strncmp(atom3->type, at->type3, l3) == 0 )
				at_found = at;
			else if ( strncmp(atom1->type, at->type3, l3) == 0 && strncmp(atom3->type, at->type1, l1) == 0 )
				at_found = at;
		} else if ( strncmp(atom2->el, at->type2, l2) == 0 ) {
			if ( strncmp(atom1->el, at->type1, l1) == 0 && strncmp(atom3->el, at->type3, l3) == 0 )
				at_found = at;
			else if ( strncmp(atom1->el, at->type3, l3) == 0 && strncmp(atom3->el, at->type1, l1) == 0 )
				at_found = at;
		}
		if ( at_found ) {
			if ( at_best ) {
				if ( fabs(at_found->angle - angle) < fabs(at_best->angle - angle) )
					at_best = at_found;
			} else {
				at_best = at_found;
			}
		}
	}
	
	if ( !at_best && ( verbose & VERB_FULL ) )
		cout << "Angle type not found: atom1=" << atom1->type << " atom2=" << 
			atom2->type << " atom3=" << atom3->type << endl;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG md_find_angle_type: atom1=" << atom1->type << " atom2=" << 
			atom2->type << " atom3=" << atom3->type << " angle=" << at_best->angle << endl;
	
	return at_best;
}
*/
/*
@brief 	Finds the reference angle between three atoms.
@param 	*atom1			first atom.
@param 	*atom2			second atom.
@param 	*atom3			third atom.
@param 	*angletype		angle type list.
@return double			angle.

	The angle type is found in the angle type list, and the angle is
	set from the angle type structure.
**/
/*double		md_find_angle(Batom* atom1, Batom* atom2, Batom* atom3, Bangletype* angletype)
{
	double			angle = 0;
	Bangletype*		at = md_find_angle_type(atom1, atom2, atom3, angletype);

	if ( at ) angle = at->angle;
	
	if ( angle < 0.1 )
		angle = md_angle(atom1, atom2, atom3);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG md_find_angle: atom1=" << atom1->type << " atom2=" << 
			atom2->type << " atom3=" << atom3->type << " angle=" << angle << endl;
	
	return angle;
}
*/
/**
@brief 	Calculates deviations from the reference bond lengths and angles.
@param 	*molgroup	molecule group.
@param 	wrap		flag to wrap around periodic boundaries.
@return double		angle.
**/
double		md_calculate_deviations(Bmolgroup* molgroup, int wrap)
{
	int				n = 0;
	double			dist, bond_dev = 0, da, angle_dev = 0;
	Vector3<double>	d1, d2;
	Bbond*			bond;
	Bangle*			angle;
	Bbond*			bondlist = molgroup->bond;
	Bangle*			anglelist = molgroup->angle;
	
	if ( verbose & VERB_FULL )
		cout << "Bond#\tAtom1\tAtom2\tdLength\tStrength" << endl;
	for ( n=0, bond = bondlist; bond; bond = bond->next, n++ ) {
		if ( wrap )
			dist = (vector3_difference_PBC(bond->atom1->coord, bond->atom2->coord, molgroup->box)).length();
		else
			dist = bond->atom1->coord.distance(bond->atom2->coord);
		dist -= bond->l;
		bond_dev += dist*dist;
		if ( verbose & VERB_FULL )
			cout << n+1 << tab << bond->atom1->num << tab << 
				bond->atom2->num << tab << dist << tab << bond->k << endl;
	}
	
	if ( n ) {
		if ( bond_dev > 0 )
			bond_dev = sqrt(bond_dev/n);
		else
			bond_dev = 0;
	}
	
	cout << "Bond deviation:                 " << bond_dev << " A" << endl;
	
	if ( verbose & VERB_FULL )
		cout << "Angle#\tAtom1\tAtom2\tAtom3\tdAngle\tStrength" << endl;
	for ( n=0, angle = anglelist; angle; angle = angle->next, n++ ) {
		if ( wrap ) {
			d1 = vector3_difference_PBC(angle->atom2->coord, angle->atom1->coord, molgroup->box);
			d2 = vector3_difference_PBC(angle->atom2->coord, angle->atom3->coord, molgroup->box);
		} else {
			d1 = angle->atom2->coord - angle->atom1->coord;
			d2 = angle->atom2->coord - angle->atom3->coord;
		}
		da = angle_set_negPI_to_PI(fabs(d1.angle(d2) - angle->a));
		angle_dev += da*da;
		if ( verbose & VERB_FULL )
			cout << n+1 << tab << angle->atom1->num << tab << angle->atom2->num << 
				tab << angle->atom3->num << tab << da*180.0/M_PI << tab << angle->k << endl;
	}
		
	if ( n ) {
		if ( angle_dev > 0 )
			angle_dev = sqrt(angle_dev/n);
		else
			angle_dev = 0;
	}
	
	cout << "Angle deviation:                " << angle_dev*180.0/M_PI << " degrees" << endl << endl;
	
	return bond_dev;
}

/**
@brief 	Calculates the radial deviations from the center of the box.
@param 	*molgroup	molecule group.
@return double		angle.
**/
int			md_calculate_radial_deviation(Bmolgroup* molgroup)
{
	int				natom;
	double			radius, radavg=0, radstd=0;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Vector3<double>	center = (molgroup->min + molgroup->max) * 0.5;
	
    for ( natom=0, mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) {
				radius = atom->coord.distance(center);
				radavg += radius;
				radstd += radius*radius;
				natom++;
			}
			
	radavg /= natom;
	radstd = radstd/natom - radavg*radavg;
	if ( radstd > 0 ) radstd = sqrt(radstd);
	else radstd = 0;
	
	cout << "Radial avg & std:               " << radavg << " " << radstd << endl << endl;
	
	return 0;
}

