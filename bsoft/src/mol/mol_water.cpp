/**
@file	mol_water.cpp
@brief	Generating and managing water
@author Bernard Heymann
@date	Created: 20001014
@date	Modified: 20110811
**/

#include "rwmolecule.h"
#include "mol_bonds.h"
#include "Matrix3.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Generates one water molecule at a given location.
@param 	**mollist		molecule list.
@param 	*watername		molecule name.
@param 	Ocoord			oxygen coordinates.
@return Bmolecule*			new water molecule.
**/
Bmolecule*	mol_generate_one_water(Bmolecule** mollist, char* watername, Vector3<double> Ocoord)
{
	double			rand_max = 1.0*get_rand_max();	
	double			OHbond = 0.9572;
	double			HOHangle = 109.47*M_PI/180.0;
	char			resname[12] = "HOH";
	char			atomname[12] = "O";
	Vector3<double>	H1vec, H2vec;
	Matrix3			mat(1);
	
	Bmolecule*		mol = molecule_add(mollist, watername);
	Bresidue*		res = residue_add(&mol->res, resname);
	res->num = 1;
	
	strcpy(atomname, "O");
	Batom*			atom = atom_add(&res->atom, atomname);
	atom->el[0] = 'O';
	atom->num = 1;
	atom->mass = 16;
	atom->chrg = -0.834;
	atom->coord = Ocoord;
	
	strcpy(atomname, "HO1");
	atom = atom_add(&atom, atomname);
	atom->el[0] = 'H';
	atom->num = 2;
	atom->mass = 1;
	atom->chrg = 0.417;
	H1vec[0] = random()/rand_max - 0.5;
	H1vec[1] = random()/rand_max - 0.5;
	H1vec[2] = random()/rand_max - 0.5;
	H1vec.normalize();
	atom->coord = (H1vec * OHbond) + Ocoord;
		
	strcpy(atomname, "HO2");
	atom = atom_add(&atom, atomname);
	atom->el[0] = 'H';
	atom->num = 3;
	atom->mass = 1;
	atom->chrg = 0.417;
	H2vec[0] = random()/rand_max - 0.5;
	H2vec[1] = random()/rand_max - 0.5;
	H2vec[2] = random()/rand_max - 0.5;
	H2vec.normalize();
	
	H2vec = H1vec.cross(H2vec);
	mat = Matrix3(H2vec, HOHangle);
	H2vec = mat * H1vec;
	atom->coord = (H2vec * OHbond) + Ocoord;
	
	return mol;
}

/**
@brief 	Generates a block of water based on a regular lattice.

	The number of water molecules generated is calculated as:
		n = volume * 0.03346.

@param 	size	size of block.
@param 	type			type of lattice, 2=rectangular, 3=tetrahedral.
@return Bmolgroup*			new molecule group.
**/
Bmolgroup*	molgroup_generate_regular_water(Vector3<double> size, int type)
{
	if ( size.volume() < 1 ) {
		cerr << "Error: A box size must be specified!" << endl;
		return NULL;
	}
	
	Bmolgroup*		molgroup = molgroup_init();
	
	char			c = 65;		// Character for molecule ID
//	double			rand_max = 1.0*get_rand_max();
	double			x, y, z;
	
	Vector3<double>	Ocoord, start;
	Vector3<double>	interval(3.1, 3.1, 3.1);			// Interval between water molecules in each direction
	if ( type == 3 ) {
		interval[0] = 3.48;
		interval[1] = interval[0]*sqrt(3.0)/2.0;
		interval[2] = interval[0]*sqrt(2.0/3.0);
	}
	
	size[0] = interval[0]*((int) (size[0]/interval[0])) + 0.001;
	size[1] = interval[1]*((int) (size[1]/interval[1])) + 0.001;
	size[2] = interval[2]*((int) (size[2]/interval[2])) + 0.001;

	double			volume = size.volume();
	long			n = (long) (volume*0.03346);
	
	Bmolecule*		mol = NULL;
	
	char			watername[20] = "A";
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating regular water" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Interval between molecules:     " << interval << " A" << endl;
		cout << "Estimated number of waters:     " << n << endl << endl;
	}
	
	if ( molgroup->mol ) for ( mol = molgroup->mol; mol->next; mol = mol->next ) ;
	
	n = 0;
	for ( z=interval[2]/2; z<size[2]; z+=interval[2] ) {
		if ( type == 3 ) start[1] = 1 - start[1];
		start[0] = start[1];
		for ( y=interval[1]/2+start[1]*interval[0]/(2*sqrt(3.0)); y<size[1]; y+=interval[1] ) {
			if ( type == 3 ) start[0] = 1 - start[0];
			for ( x=(2-start[0])*interval[0]/2; x<size[0]; x+=interval[0] ) {
				watername[0] = c;
				Ocoord[0] = x;
				Ocoord[1] = y;
				Ocoord[2] = z;
				mol = mol_generate_one_water(&mol, watername, Ocoord);
				if ( !molgroup->mol ) molgroup->mol = mol;
				n++;
				c++;
				if ( c > 90 ) c = 65;
			}
		}
	}
	
	molgroup_stats(molgroup);
	
	molgroup->box[0] = molgroup->max[0] + interval[0]/2;
	molgroup->box[1] = molgroup->max[1] + interval[1]/2;
	molgroup->box[2] = molgroup->max[2] + interval[2]/2;
	
	molgroup->box = size;
	
	volume = molgroup->box.volume();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of waters:               " << n << endl;
		cout << "Volume:                         " << molgroup->box << " = " << volume << " A3" << endl;
		cout << "Density:                        " << n*1.0/volume << 
			" molecules/A3 (" << n/(volume*0.03346) << ")" << endl << endl;
	}
	
	return molgroup;
}

/**
@brief 	Generates a block of water with random placement.

	The number of water molecules generated is calculated as:
		n = volume * 0.03346.

@param 	size	size of block.
@return Bmolgroup*			new molecule group.
**/
Bmolgroup*	molgroup_generate_random_water(Vector3<double> size)
{
	if ( size.volume() < 1 ) {
		cerr << "Error: A box size must be specified!" << endl;
		return NULL;
	}
	
	random_seed();
	
	int				i;
	char			c = 65;		// Character for molecule ID
	double			rand_max = 1.0*get_rand_max();
	
	double			volume = size.volume();
	long			n = (long) (volume*0.03346);
	Vector3<double>	Ocoord;
	
	Bmolecule*		mol = NULL;
	Batom*			atom;
	
	char			watername[20] = "A";
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Generating random water" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Volume:                         " << size << " = " << volume << " A3" << endl;
		cout << "Number of waters:               " << n << endl << endl;
	}
	
	Bmolgroup*		molgroup = molgroup_init();
	
	if ( molgroup->mol ) for ( mol = molgroup->mol; mol->next; mol = mol->next ) ;
	
	for ( i=0; i<n; i++ ) {
		watername[0] = c;
		Ocoord[0] = size[0]*random()/rand_max;
		Ocoord[1] = size[1]*random()/rand_max;
		Ocoord[2] = size[2]*random()/rand_max;
		mol = mol_generate_one_water(&mol, watername, Ocoord);
		if ( !molgroup->mol ) molgroup->mol = mol;
		for ( atom = mol->res->atom; atom; atom = atom->next )
			atom->coord = vector3_set_PBC(atom->coord, size);
		c++;
		if ( c > 90 ) c = 65;
	}
	
	molgroup_stats(molgroup);
	
	molgroup->box = size;
	
	volume = molgroup->box.volume();
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of waters:               " << n << endl;
		cout << "Volume:                         " << molgroup->box << " = " << volume << " A3" << endl;
		cout << "Density:                        " << n*1.0/volume << 
			" molecules/A3 (" << n/(volume*0.03346) << ")" << endl << endl;
	}
	
	return molgroup;
}

/**
@brief 	Generates a bond list for a block of waters.
@param 	*molgroup	molecule group.
@return Bbond*				new bond list.
**/
Bbond*		water_bond_list(Bmolgroup* molgroup)
{
//	if ( molgroup->bond ) {
//		cerr << "Warning: Bond list already defined!" << endl;
//		return molgroup->bond;
//	}
	
	if ( molgroup->bond ) {
		if ( verbose ) cerr << "Warning: Deleting original bond list!" << endl;
		bond_kill(molgroup->bond);
	}
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Bbond*			bondlist = NULL;
	Bbond*			bond = NULL;
	double			OHbond = 0.9572;
	int				nbond = 0;
	
    for ( nbond = 0, mol = molgroup->mol; mol; mol = mol->next ) {
		res = mol->res;
		atom = res->atom;
		if ( bond ) bond = bond_add(&bond, atom, atom->next, OHbond, 1);
		else bond = bond_add(&bondlist, atom, atom->next, OHbond, 1);
		bond = bond_add(&bond, atom, atom->next->next, OHbond, 1);
		nbond += 2;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of bonds generated:      " << nbond << endl;
	
	molgroup->bond = bondlist;
	
	return bondlist;
}

/**
@brief 	Generates a bond angle list for a block of waters.
@param 	*molgroup	molecule group.
@return Bangle*				new bond angle list.
**/
Bangle*		water_angle_list(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Bangle*			angle = NULL;
	double			OHangle = 109.47*M_PI/180.0;
	int				nangle = 0;

	molgroup->angle = NULL;
	
    for ( nangle = 0, mol = molgroup->mol; mol; mol = mol->next ) {
		res = mol->res;
		atom = res->atom;
		angle = angle_add(&angle, atom->next, atom, atom->next->next, OHangle, 1);
		if ( !molgroup->angle ) molgroup->angle = angle;
		nangle ++;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of angles generated:     " << nangle << endl;
	
	return molgroup->angle;
}

/**
@brief 	Calculates a radial distribution function for water molecules.
@param 	*molgroup		molecule group.
@param 	interval		interval between bins.
@param 	cutoff			distance cutoff.
@return int				0.
**/
int			molgroup_calc_water_rdf(Bmolgroup* molgroup, double interval, double cutoff)
{
	int				i, ii, j, t, x, y, z, xx, yy, zz, ix, iy, iz;
	Vector3<double>	box = molgroup->box;
	double			dist;
	Vector3<double>	d, sampling(cutoff, cutoff, cutoff);
	Vector3<int>	size((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/size[i] + 0.001;

	Bmolecule*		mol;
	Bresidue*		res;
	
	Latom**			alist = molgroup_atom_mesh_lists(molgroup, size, sampling);
	Latom			*latom, *latom2;
	
	double			mult = 1.0/interval;
	int				n = (int) (mult*cutoff);
	int*			rdf = new int[n*3];
	for ( i=0; i<n*3; i++ ) rdf[i] = 0;
	
	for ( i=z=0; z<size[2]; z++ ) {
		for ( y=0; y<size[1]; y++ ) {
			for ( x=0; x<size[0]; x++, i++ ) {
				for ( latom = alist[i]; latom; latom = latom->next ) {
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
								for ( latom2 = alist[ii]; latom2; latom2 = latom2->next ) 
										if ( latom2->atom != latom->atom ) {
									d = vector3_difference_PBC(latom2->atom->coord, latom->atom->coord, box);
									dist = d.length();
									if ( dist < cutoff ) {
										j = (int) (mult*dist + 0.5);
										if ( j < n ) {
											if ( latom->atom->el[0] == 'H' && latom2->atom->el[0] == 'H' ) t = 0;
											else if ( latom->atom->el[0] == 'H' || latom2->atom->el[0] == 'H' ) t = 1;
											else t = 2;
											rdf[3*j+t] += 1;
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

	for ( i=0; i<size.volume(); i++ ) kill_list((char *) alist[i], sizeof(Latom));
	delete[] alist;
	
	int			nwater = 0;
	for ( nwater = 0, mol = molgroup->mol; mol; mol = mol->next )
		for ( res = mol->res; res; res = res->next ) nwater++;
			
	cout << "Calculating the radial distribution function:" << endl;
	
	double		vol_per_water = 180/6.022;	// Volume occupied by water = 29.9 A^3
//	double		vol_within_cutoff = (4.0/3.0)*M_PI*cutoff*cutoff*cutoff;
//	double		water_within_cutoff = vol_within_cutoff/vol_per_water;
	long		ndist[3] = {0,0,0};
	double		rdftot[3] = {0,0,0};
	double		norm[3], vol_surface;
	double		vol_surf_fac = (4.0/3.0)*M_PI*interval*interval*interval;
	cout << endl << "r(A)\tH-H\tRDF(HH)\tH-O\tRDF(HO)\tO-O\tRDF(OO)" << endl;
	for ( i=1; i<n; i++ ) {
		cout << i*interval;
		vol_surface = vol_surf_fac*(3*i*i + 0.25);
		norm[2] = vol_per_water/(vol_surface*2*nwater);		// Normalizer for RDF(OO)
		norm[0] = norm[1] = 0.25*norm[2];					// Volume per H is half and 2 times more H's
		for ( t=0; t<3; t++ ) {
			cout << tab << rdf[3*i+t] << tab << rdf[3*i+t]*norm[t];
			ndist[t] += rdf[3*i+t];
			rdftot[t] += rdf[3*i+t]*norm[t];
		}
		cout << endl;
	}
	cout << "Number of distances included = " << ndist[0] << " " << ndist[1] << " " << ndist[2] << endl;
	cout << "RDF sum = " << rdftot[0] << " " << rdftot[1] << " " << rdftot[2] << endl;
	cout << "RDF avg = " << rdftot[0]/n << " " << rdftot[1]/n << " " << rdftot[2]/n << endl << endl;
	
	delete[] rdf;
	
	return 0;
}

