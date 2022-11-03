/**
@file	mol_edit.cpp
@brief	Library routines used for atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20211029
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "rwresprop.h"
#include "mol_transform.h"
#include "mol_select.h"
#include "mol_edit.h"
#include "mol_util.h"
#include "seq_util.h"
#include "linked_list.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief     Sets atom types to element names.
@param 	*molgroup 	molecule group structure.
@return long					number of atoms.
**/
long		molgroup_set_atom_types_to_elements(Bmolgroup* molgroup)
{
	long			n(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;

    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				strcpy(atom->type, atom->el);

	if ( verbose & VERB_PROCESS )
		cout << "Atoms renamed    :              " << n << endl << endl;
	
	return n;
}

/**
@brief	Removes all hydrogens from a molecular system.
@param 	*molgroup 	molecule group structure.
@return int			number of hydrogens removed.
**/
int			molgroup_remove_hydrogens(Bmolgroup* molgroup)
{
	long			n(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom			*atom, *atom2;

    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = atom2 = res->atom; atom; ) {
				if ( atom->tnum == 0 ) {
					if ( atom == atom2 ) {
						atom2 = res->atom = atom->next;
						delete atom;
						atom = atom2;
					} else {
						atom2->next = atom->next;
						delete atom;
						atom = atom2->next;
					}
					n++;
				} else {
					atom2 = atom;
					atom = atom2->next;
				}
			}

	if ( verbose & VERB_PROCESS )
		cout << "Hydrogens removed:              " << n << endl << endl;
	
	return n;
}

bool		bond_find(Bmolgroup* molgroup, Batom* atom1, Batom* atom2)
{
	Bbond*		bond;
	
	for ( bond = molgroup->bond; bond; bond = bond->next )
		if ( bond->atom1 == atom1 && bond->atom2 == atom2 ) return 1;
	
	return 0;
}

/**
@brief	Adds disulfide bonds when sulfur atoms are close to each other.
@param 	*molgroup 	molecule group structure.
@param	distance	maximum separation between sulfur atoms.
@return int			number of disulfide bonds added.
**/
int			molgroup_add_disulfides(Bmolgroup* molgroup, double distance)
{
	long			n(0);
	double			d;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Latom*			alist = NULL;
	Latom*			la = NULL;
	Latom*			la2 = NULL;

	if ( verbose )
		cout << "Adding disulfide bonds:" << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) if ( strstr(res->type, "CYS") ) {
			for ( atom = res->atom; atom; atom = atom->next ) if ( strstr(atom->type, "SG" ) ) {
				la = (Latom *) add_item((char **) &la, sizeof(Latom));
				if ( !alist ) alist = la;
				la->atom = atom;
				n++;
			}
		}
	}
	
	if ( verbose )
		cout << "Cysteine sulfurs found:         " << n << endl;
	
	for ( n=0, la = alist; la->next; la = la->next ) {
		for ( la2 = la->next; la2; la2 = la2->next ) {
			d = la->atom->coord.distance(la2->atom->coord);
			if ( d <= distance ) {
				if ( !bond_find(molgroup, la->atom, la2->atom) ) {
					bond_add(&molgroup->bond, la->atom, la2->atom, d, 1);
					if ( verbose )
						cout << la->atom->num << tab << la2->atom->num << tab << d << endl;
					n++;
				}
			}
		}
	}

	kill_list((char *) alist, sizeof(Latom));

	if ( verbose & VERB_PROCESS )
		cout << "Disulfide bonds created:        " << n << endl << endl;
	
	return n;
}

/**
@brief	Splits a set of coordinates into slices with a given thickness.
@param 	*molgroup 	molecule group structure.
@param 	slice_thickness	slice thickness (in angstrom).
@param 	&nslices			pointer to the number of slices generated
@return Bmolgroup**				set of molecule groups.
**/
Bmolgroup**	molgroup_split_into_slices(Bmolgroup* molgroup, double slice_thickness, int& nslices)
{
	molgroup_stats(molgroup);
	
	int				i;
	double			min = molgroup->min[2];
	char			molname[20] = "";
	char			restype[20] = "UNK";
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

	nslices = 1;
	if ( slice_thickness < molgroup->max[2] - molgroup->min[2] ) {
		min = slice_thickness*floor(molgroup->min[2]/slice_thickness);
		nslices = (int) ((molgroup->max[2] - min)/slice_thickness + 1);
	}
	
	Bmolgroup**		slice_molgroup = new Bmolgroup*[nslices];
	Bmolecule*		slice_mol;
	Bresidue**		slice_res = new Bresidue*[nslices];
	Batom**  		slice_atom = new Batom*[nslices];
	
	for ( i=0; i<nslices; i++ )	{
		slice_molgroup[i] = NULL;
		slice_res[i] = NULL;
		slice_atom[i] = NULL;
	}
	
	if ( verbose ) {
		cout << "Splitting molecule group:" << endl;
		cout << "Number of slices:               " << nslices << endl;
		cout << "Slice thickness:                " << slice_thickness << " A" << endl << endl;
	}
	
	for ( i=0; i<nslices; i++ ) {
		snprintf(molname, 20, "%d", i+1);
		slice_molgroup[i] = molgroup_init();
		slice_mol = molecule_add(&slice_molgroup[i]->mol, molname);
		slice_res[i] = residue_add(&slice_mol->res, restype);
//		slice_molgroup[i]->comment += molgroup->comment;
	}
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		if ( verbose & VERB_PROCESS )
			cout << "Splitting molecule " << mol->id << endl;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				i = (int) ((atom->coord[2] - min)/slice_thickness);
				if ( i < 0 || i >= nslices ) {
					error_show("Error in molgroup_split_into_slices", __FILE__, __LINE__);
					cerr << "Index out of range: " << i << endl;
					return NULL;
				}
//				if ( slice_atom[i] )
//					slice_atom[i] = atom_add(&slice_atom[i], atom->type);
//				else
//					slice_atom[i] = atom_add(&slice_res[i]->atom, atom->type);
//				memcpy(slice_atom[i], atom, sizeof(Batom));
				if ( slice_res[i]->atom )
					slice_atom[i] = slice_atom[i]->next = atom_copy(atom);
				else
					slice_res[i]->atom = slice_atom[i] = atom_copy(atom);
//				slice_atom[i]->next = NULL;
			}
		}
	}
	
	delete[] slice_res;
	delete[] slice_atom;

	for ( i=0; i<nslices; i++ ) {
		if ( verbose & VERB_PROCESS )
			cout << "Slice " << i+1 << ":" << endl;
		molgroup_stats(slice_molgroup[i]);
	}
	
	return slice_molgroup;
}

/**
@brief     Inserts one set of molecules into another.

	Molecules overlapping in the receiving molecule group are deleted.
	The footprint of the molecules being inserted is calculated on a grid
	and all atoms within this footprint is tested for deletion.
	Note: The molecule list is transferred from the insertion group to 
		the main group and the insertion group is deallocated.

@param 	*molgroup 	molecule group structure to be modified.
@param 	*molinsert 	molecule group structure to insert. (deallocated)
@param 	distance			cutoff distance to remove atoms.
@return int						0.
**/
int			molgroup_insert(Bmolgroup* molgroup, Bmolgroup* molinsert, double distance) 
{
	if ( distance <= 0 ) distance = 2;  // Default
	if ( distance < 1 ) distance = 1;   // Limits on cutoff distance in angstrom
	if ( distance > 5 ) distance = 5;
	
	long	i;
	int				delete_res, nresdel(0), nmoldel(0);
	long			ii, x, y, z, xx, yy, zz, ix, iy, iz;
	Vector3<double>	sampling(distance, distance, distance);
	Bmolecule*		mol, *prevmol;
	Bresidue*		res, *prevres;
	Batom*  		atom;
	Vector3<double>	box = molgroup->box;
	Vector3<int>	gridsize((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	gridsize = gridsize.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/gridsize[i] + 0.001;
	long	gridvol = (long) gridsize.volume();
	
	if ( verbose )
		cout << "Inserting a molecule group and deleting overlapping residues" << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Distance cutoff:                " << distance << " A" << endl;

	Latom*			latom;
	Latom**			alist = molgroup_atom_mesh_lists(molinsert, gridsize, sampling);
	
	// Find the atoms under the footprint to be deleted
	for ( mol = prevmol = molgroup->mol; mol; ) {
		for( res = prevres = mol->res; res; ) {
			delete_res = 0;
			for ( atom = res->atom; atom; atom = atom->next ) {
				x = (long) ((atom->coord[0] - molgroup->min[0])/sampling[0]);
				y = (long) ((atom->coord[1] - molgroup->min[1])/sampling[1]);
				z = (long) ((atom->coord[2] - molgroup->min[2])/sampling[2]);
				if ( x >=0 && x < gridsize[0] && y >= 0 && y < gridsize[1] && z >= 0 && z < gridsize[2] ) {
					i = (z*gridsize[1] + y)*gridsize[0] + x;
					for ( zz=z-1; zz<=z+1; zz++ ) {
						iz = zz;
						if ( iz < 0 ) iz += gridsize[2];
						if ( iz >= gridsize[2]) iz -= gridsize[2];
						for ( yy=y-1; yy<=y+1; yy++ ) {
							iy = yy;
							if ( iy < 0 ) iy += gridsize[1];
							if ( iy >= gridsize[1] ) iy -= gridsize[1];
							for ( xx=x-1; xx<=x+1; xx++ ) {
								ix = xx;
								if ( ix < 0 ) ix += gridsize[0];
								if ( ix >= gridsize[0] ) ix -= gridsize[0];
								ii = (iz*gridsize[1] + iy)*gridsize[0] + ix;
								for ( latom = alist[ii]; latom; latom = latom->next ) {
									if ( atom->coord.distance(latom->atom->coord) < distance )
										delete_res = 1;
								}
							}
						}
					}
				}
			}
			if ( delete_res ) {
				if ( verbose & VERB_FULL )
					cout << "Removing a residue" << endl;
				if ( res == mol->res ) {
					mol->res = prevres = res->next;
					residue_kill(res);
					res = mol->res;
				} else {
					prevres->next = res->next;
					residue_kill(res);
					res = prevres->next;
				}
				nresdel++;
			} else {
				prevres = res;
				if ( res ) res = res->next;
			}
		}
		if ( !mol->res ) {
			if ( verbose & VERB_FULL )
				cout << "Removing a molecule" << endl;
			if ( molgroup->mol == mol ) {
				molgroup->mol = prevmol = mol->next;
				molecule_kill(mol);
				mol = molgroup->mol;
			} else {
				prevmol->next = mol->next;
				molecule_kill(mol);
				mol = prevmol->next;
			}
			nmoldel++;
		} else {
			prevmol = mol;
			if ( mol ) mol = mol->next;
		}
	}

	for ( i=0; i<gridvol; i++ ) kill_list((char *) alist[i], sizeof(Latom));
	delete[] alist;

	if ( verbose & VERB_PROCESS )
		cout << "Molecules and residues deleted: " << nmoldel << " " << nresdel << endl;
	
	// Add the new molecules
	if ( molgroup->mol ) {
		for ( mol = molgroup->mol; mol->next; mol = mol->next ) ;
		mol->next = molinsert->mol;
	} else  molgroup->mol = molinsert->mol;
	molinsert->mol = NULL;
	
	if ( verbose & VERB_PROCESS )
		molgroup_density(molgroup) ;
	
	molgroup_stats(molgroup);
	
	return 0;
}

/**
@brief	Applies random displacements to coordinates.
@param 	*molgroup 	molecule group structure to be modified.
@param 	random_max	maximum displacement.
@return int			0.
**/
int			molgroup_randomize(Bmolgroup* molgroup, double random_max)
{
	if ( random_max < 1e-30 ) return 0;

	random_seed();
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	double			irm = random_max*2.0/get_rand_max();
	
	if ( verbose & VERB_PROCESS )
		cout << "Randomizing coordinates to maximum: " << random_max << " angstrom" << endl << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord[0] += irm*random() - random_max;
				atom->coord[1] += irm*random() - random_max;
				atom->coord[2] += irm*random() - random_max;
			}
		}
    }
	
	return 0;
}

/**
@brief	Applies random displacements to coordinates.
@param 	*molgroup 	molecule group structure to be modified.
@param 	B			B factor.
@return int			0.
**/
int			molgroup_randomize_B(Bmolgroup* molgroup, double B)
{
	if ( B < 1e-30 ) return 0;
	
	random_seed();
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	double			stdev(sqrt(B));
	
	if ( verbose & VERB_PROCESS )
		cout << "Randomizing coordinates with Bfactor: " << B << " Å^2" << endl << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord += vector3_random_gaussian(0, stdev);
			}
		}
    }
	
	return 0;
}

/**
@brief	Applies random displacements to a selected number of coordinates.
@param 	*molgroup 	molecule group structure to be modified.
@param 	number		number of coordinates to displace.
@param 	stdev		standard devaition of displacement.
@return int			0.
**/
int			molgroup_random_displace_number(Bmolgroup* molgroup, long number, double stdev)
{
	if ( stdev < 1e-30 ) return 0;
	
	molgroup_select(molgroup, number);

	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	if ( verbose & VERB_PROCESS )
		cout << "Randomizing " << number << " coordinates to standard deviation " << stdev << " angstrom" << endl << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) if ( atom->sel )
				atom->coord += vector3_random_gaussian(0, stdev);
    
    molgroup_select_all(molgroup);
	
	return 0;
}


/**
@brief 	Removes overlapping atoms.
@param 	*molgroup	molecule group.
@param 	mindist		minimum distance allowed between atoms.
@return int 			number of atoms removed, <0 on error.

	The input molecule group is checked for any atom pairs closer than
	a minimum allowed distance. The second atom of an overlapping pair
	is removed. This is intended to clean up after symmetry operations
	that generate overlapping pseudo-atoms lying on symmetry axes.

**/
int			molgroup_remove_overlapping_atoms(Bmolgroup* molgroup, double mindist)
{
	if ( mindist < 0.001 ) return 0;
	
	int				nrem(0);
	double			dist;
	Bmolecule		*mol, *mol2, *mol3;
	Bresidue		*res, *res2, *res3;
	Batom			*atom, *atom2, *atom3;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_remove_overlapping_atoms: mindist=" << mindist << endl;
	
    for ( mol = molgroup->mol; mol && mol->next; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				for ( mol2 = mol3 = mol; mol2; mol2 = mol2 ? mol2->next : NULL ) {
					for( res2 = res3 = mol2->res; res2; res2 = res2 ? res2->next : NULL ) {
						for ( atom2 = atom3 = res2->atom; atom2; atom2 = atom2 ? atom2->next : NULL ) {
							if ( atom != atom2 ) {
								dist = atom->coord.distance(atom2->coord);
								if ( dist < mindist ) {
									if ( atom2 == res2->atom ) {
										res2->atom = atom2->next;
										delete atom2;
										atom2 = res2->atom;
									} else {
										atom3->next = atom2->next;
										delete atom2;
										atom2 = atom3;
									}
									nrem++;
								}
							}
							atom3 = atom2;
						}
						if ( !res2->atom ) {
							if ( res2 == mol2->res ) {
								mol2->res = res2->next;
								residue_kill(res2);
								res2 = mol2->res;
							} else {
								res3->next = res2->next;
								residue_kill(res2);
								res2 = res3;
							}
						}
						res3 = res2;
					}
					if ( !mol2->res ) {
						if ( mol2 == molgroup->mol ) {
							molgroup->mol = mol2->next;
							molecule_kill(mol2);
							mol2 = molgroup->mol;
						} else {
							mol3->next = mol2->next;
							molecule_kill(mol2);
							mol2 = mol3;
						}
					}
					mol3 = mol2;
				}
			}
		}
	}
	
	if ( molgroup->mol == NULL ) {
		error_show("Error in molgroup_remove_overlapping_atoms: All molecules deleted!", __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Overlapping atoms removed:      " << nrem << endl << endl;
	
	molgroup_stats(molgroup);
	
	return nrem;
}
	
/**
@brief 	Places pseudo-atoms on bonds.
@param 	*molgroup 		the molecule group.
@param 	atoms_per_bond	number of pseudo-atoms per bond.
@param 	wrap			wrap around periodic boundaries if !=0.
@return Bbond* 			new bond list.

	The requested number of new atoms are placed to coincide with each bond.

**/
int			molgroup_bond_pseudo_atoms(Bmolgroup* molgroup, int atoms_per_bond, int wrap)
{
	int				i;
	Vector3<double>	v;
	Bresidue*		res = residue_add(&molgroup->mol->res, "PSD");
	Batom*			atom;
	Bbond*			bond;
	
	for ( bond=molgroup->bond; bond; bond=bond->next ) {
		if ( wrap ) v = vector3_difference_PBC(bond->atom2->coord, bond->atom1->coord, molgroup->box);
		else v = bond->atom2->coord - bond->atom1->coord;
		for ( i=1; i<=atoms_per_bond; i++ ) {
			atom = atom_add(&res->atom, "A");
			atom->coord = bond->atom1->coord + (v * (i*1.0L/(atoms_per_bond+1)));
			if ( wrap ) atom->coord = vector3_set_PBC(atom->coord, molgroup->box);
		}
	}
	
	return 0;
}

char*		grid_from_molgroup(Bmolgroup* molgroup, Vector3<double> min, double sampling, Vector3<int> size)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	long			x, y, z;
	long	i, vol = (long)size.volume();
	
	char*			grid = new char[vol];
	for ( i=0; i<vol; i++ ) grid[i] = 0;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				x = (long) ((atom->coord[0] - min[0])/sampling);
				y = (long) ((atom->coord[1] - min[1])/sampling);
				z = (long) ((atom->coord[2] - min[2])/sampling);
				grid[(z*size[1] + y)*size[0] + x] = 1;
			}
		}
	}
	
	return grid;
}

char*		grid_from_molecule(Bmolecule* mol, Vector3<double> min, double sampling, Vector3<int> size)
{
	Bresidue*		res;
	Batom*			atom;
	long			x, y, z;
	long	i, vol = (long)size.volume();
	
	char*			grid = new char[vol];
	for ( i=0; i<vol; i++ ) grid[i] = 0;
	
	for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			x = (long) ((atom->coord[0] - min[0])/sampling);
			y = (long) ((atom->coord[1] - min[1])/sampling);
			z = (long) ((atom->coord[2] - min[2])/sampling);
			grid[(z*size[1] + y)*size[0] + x] = 1;
		}
	}
	
	return grid;
}

/**
@brief 	Deletes overlapping molecules.
@param 	*molgroup 	the molecule group.
@return int					0.
**/
int			molgroup_prune_molecules(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	
	long	i;
	double			sampling = 4;
	Vector3<int>	size;	
	size[0] = (int) (molgroup->box[0]/sampling + 1);
	size[1] = (int) (molgroup->box[1]/sampling + 1);
	size[2] = (int) (molgroup->box[2]/sampling + 1);	
	long	vol = (long)size.volume();

	char*			grid1 = NULL;
	char*			grid = grid_from_molecule(molgroup->mol, molgroup->min, sampling, size);
	molgroup->mol->sel = 1;
	
	for ( mol = molgroup->mol->next; mol; mol = mol->next ) {
		mol->sel = 0;
		grid1 = grid_from_molecule(mol, molgroup->min, sampling, size);
		for ( i=0; i<vol; i++ ) {
			if ( grid1[i] ) {
				if ( grid[i] ) mol->sel++;
				else grid[i] = 1;
			}
		}
		delete[] grid1;
		if ( mol->sel ) mol->sel = 0;
		else mol->sel = 1;
	}
	
	delete[] grid;
	
	molgroup_delete_deselected_molecules(molgroup);

	return 0;
}

/**
@brief     Deletes deselected molecules.
@param 	*molgroup 	molecule group structure to be modified.
@return long					number of molecules deleted.
**/
long		molgroup_delete_deselected_molecules(Bmolgroup* molgroup)
{
	long			nmol(0), ndel(0);
	Bmolecule*		mol, *mol2;
	
    for ( mol = mol2 = molgroup->mol; mol;  ) {
		if ( mol->sel < 1 ) {
			if ( verbose & VERB_FULL )
				cout << "Deleting molecule " << mol->id << endl;
			if ( mol == molgroup->mol ) {
				molgroup->mol = mol->next;
				molecule_kill(mol);
				mol = mol2 = molgroup->mol;
			} else {
				mol2->next = mol->next;
				molecule_kill(mol);
				mol = mol2->next;
			}
			ndel++;
		} else {
			mol2 = mol;
			mol = mol->next;
		}
		nmol++;
	}
	
	if ( verbose )
		cout << "Molecules deleted:              " << ndel << " (" << ndel*100.0/nmol << " %)" << endl << endl;
	
	return ndel;
}

/**
@brief     Prunes overlapping atoms based on a distance criterion.

	The first atom in any pair of overlapping atoms is kept.

@param 	*molgroup 	molecule group structure to be modified.
@param 	mindist			distance criterion.
@return long					number of remaining atoms.
**/
long		molgroup_prune_overlapping_atoms(Bmolgroup* molgroup, double mindist)
{

	vector<Batom*>	atomarray = atom_get_array(molgroup);

	long			i, j, natom(atomarray.size()), ndel(0);

	for ( i=0; i<natom-1; i++ ) atomarray[i]->sel = 1;
	
	for ( i=0; i<natom-1; i++ ) if ( atomarray[i]->sel ) {
		for ( j=i+1; j<natom; j++ ) if ( atomarray[j]->sel ) {
			if ( atomarray[i]->coord.distance(atomarray[j]->coord) < mindist ) {
				atomarray[j]->sel = 0;
				ndel++;
			}
		}
	}

	if ( verbose )
		cout << "Atoms marked for deletion:      " << ndel << " (" << ndel*100.0/natom << ")" << endl << endl;
	
	ndel = molgroup_delete_deselected_atoms(molgroup);
	
	return natom - ndel;
}

/**
@brief     Deletes deselected atoms.
@param 	*molgroup 	molecule group structure to be modified.
@return long					number of atoms deleted.
**/
long		molgroup_delete_deselected_atoms(Bmolgroup* molgroup)
{
	long			natom(0), ndel(0);
	Bmolecule*		mol, *mol2;
	Bresidue*		res, *res2;
	Batom*  		atom, *atom2;
	
    for ( mol = mol2 = molgroup->mol; mol;  ) {
		for( res = res2 = mol->res; res;  ) {
			for ( atom = atom2 = res->atom; atom;  ) {
				if ( atom->sel < 1 ) {
					if ( verbose & VERB_FULL )
						cout << "Deleting atom " << atom->num << endl;
					if ( atom == res->atom ) {
						res->atom = atom->next;
						delete atom;
						atom = atom2 = res->atom;
					} else {
						atom2->next = atom->next;
						delete atom;
						atom = atom2->next;
					}
					ndel++;
				} else {
					atom2 = atom;
					atom = atom->next;
				}
				natom++;
			}
			if ( !res->atom ) {
				if ( verbose & VERB_FULL )
					cout << "Deleting residue " << res->num << endl;
				if ( res == mol->res ) {
					mol->res = res->next;
					residue_kill(res);
					res = res2 = mol->res;
				} else {
					res2->next = res->next;
					residue_kill(res);
					res = res2->next;
				}
			} else {
				res2 = res;
				res = res->next;
			}
		}
		if ( !mol->res ) {
			if ( verbose & VERB_FULL )
				cout << "Deleting molecule " << mol->id << endl;
			if ( mol == molgroup->mol ) {
				molgroup->mol = mol->next;
				molecule_kill(mol);
				mol = mol2 = molgroup->mol;
			} else {
				mol2->next = mol->next;
				molecule_kill(mol);
				mol = mol2->next;
			}
		} else {
			mol2 = mol;
			mol = mol->next;
		}
	}
	
	if ( verbose )
		cout << "Atoms deleted:                  " << ndel << " (" << ndel*100.0/natom << ")" << endl << endl;
	
	return ndel;
}

/**
@brief 	Moves overlapping molecules away from each other.

	The overlap of molecules are assessed by projecting atom positions
	onto a grid. Overlapping molecules are then moved away from each other
	along a vector through their centers-of-mass.

@param 	*molgroup 	the molecule group.
@param 	sampling		grid sampling (angstrom).
@param 	lambda		damping factor.
@return int					0.
**/
int			molgroup_untangle_molecules(Bmolgroup* molgroup, double sampling, double lambda)
{
	Bmolecule		*mol, *mol2;

	long			iter, i, x, y, z, m, m2, nmol;
	double			d, dmax, all_shift = 1;
	Vector3<int>	size;	
	Vector3<double>	com, com2, u, c, v;
	char			*grid, *grid2;

	for ( nmol=0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( verbose )
		cout << "Untangling " << nmol << " molecules" << endl;
		
	Vector3<double>*	F = new Vector3<double>[nmol];

	for ( iter=1; all_shift > 0 && iter < 100; iter++ ) {
		all_shift = 0;
		size[0] = (int) (molgroup->box[0]/sampling + 1);
		size[1] = (int) (molgroup->box[1]/sampling + 1);
		size[2] = (int) (molgroup->box[2]/sampling + 1);	
		for ( m=0, mol = molgroup->mol; mol->next; mol = mol->next, m++ ) {
			com = mol_center_of_mass(mol);
			grid = grid_from_molecule(mol, molgroup->min, sampling, size);
			for ( m2=m+1, mol2 = mol->next; mol2; mol2 = mol2->next, m2++ ) {
				com2 = mol_center_of_mass(mol2);
				grid2 = grid_from_molecule(mol2, molgroup->min, sampling, size);
				u = com2 - com;
				u.normalize();
				c = (com + com2)/2;
				dmax = 0;
				for ( i=z=0; z<size[2]; z++ ) {
					for ( y=0; y<size[1]; y++ ) {
						for ( x=0; x<size[0]; x++, i++ ) {
							if ( grid[i] & grid2[i] ) {
								v = Vector3<double>(x, y, z);
								v = (v * sampling + molgroup->min - c);
								d = v.scalar(u);
//								d = fabs(v.scalar(u));
								if ( dmax < d ) dmax = d;
							}
						}
					}
				}
				delete[] grid2;
				if ( dmax ) {
					dmax *= lambda;
					v = u * dmax;
					F[m] -= v;
					F[m2] += v;
					all_shift += dmax;
				}
			}
			delete[] grid;
		}
		for ( m=0, mol = molgroup->mol; mol; mol = mol->next, m++ ) {
			mol_coor_shift(mol, F[m]);
			F[m] = 0;
		}
		molgroup_stats(molgroup);
		molgroup->box = molgroup->max - molgroup->min;
		if ( verbose )
			cout << iter << tab << all_shift << endl;
	}

	delete[] F;

	return 0;
}

/**
@brief 	Moves overlapping molecules away from each other.

	The overlap of molecules are assessed by projecting atom positions
	onto a grid. Overlapping molecules are then moved away from each other
	along a vector through their centers-of-mass.

@param 	*molgroup 	the molecule group.
@param 	sampling		grid sampling (angstrom).
@param 	lambda		damping factor.
@return int					0.
**/
int			molgroup_untangle_groups(Bmolgroup* molgroup, double sampling, double lambda)
{
	Bmolgroup		*mg, *mg2;

	long			iter, i, x, y, z, m, m2, nmg;
	double			d, dmax, all_shift = 1;
	Vector3<int>	size;	
	Vector3<double>	min, max, box, com, com2, u, c, v;
	char			*grid, *grid2;

	for ( nmg=0, mg = molgroup; mg; mg = mg->next ) nmg++;

	if ( verbose ) {
		cout << "Untangling " << nmg << " molecule groups" << endl;
		cout << "Sampling:                       " << sampling << " A" << endl;
		cout << "Damping factor:                 " << lambda << endl;
	}
		
	Vector3<double>*	F = new Vector3<double>[nmg];

	if ( verbose )
		cout << "Iter\tShift" << endl;
	for ( iter=1; all_shift > 0 && iter < 100; iter++ ) {
		all_shift = 0;
		min = molgroup->min;
		max = molgroup->max;
		for ( mg = molgroup; mg; mg = mg->next ) {
			min = min.min(mg->min);
			max = max.max(mg->max);
		}
		box = max - min;
		size[0] = (int) (box[0]/sampling + 1);
		size[1] = (int) (box[1]/sampling + 1);
		size[2] = (int) (box[2]/sampling + 1);	
		for ( m=0, mg = molgroup; mg->next; mg = mg->next, m++ ) {
			com = molgroup_center_of_mass(mg);
			grid = grid_from_molgroup(mg, min, sampling, size);
			for ( m2=m+1, mg2 = mg->next; mg2; mg2 = mg2->next, m2++ ) {
				com2 = molgroup_center_of_mass(mg2);
				grid2 = grid_from_molgroup(mg2, min, sampling, size);
				u = com2 - com;
				u.normalize();
				c = (com + com2)/2;
				dmax = 0;
				for ( i=z=0; z<size[2]; z++ ) {
					for ( y=0; y<size[1]; y++ ) {
						for ( x=0; x<size[0]; x++, i++ ) {
							if ( grid[i] & grid2[i] ) {
								v = Vector3<double>(x, y, z);
								v = (v * sampling + min - c);
								d = v.scalar(u);
//								d = fabs(v.scalar(u));
								if ( dmax < d ) dmax = d;
							}
						}
					}
				}
				delete[] grid2;
				if ( dmax ) {
					dmax *= lambda;
					v = u * dmax;
					F[m] -= v;
					F[m2] += v;
					all_shift += dmax;
				}
			}
			delete[] grid;
		}
		for ( m=0, mg = molgroup; mg; mg = mg->next, m++ ) {
			if ( F[m].length() ) molgroup_coor_shift(mg, F[m]);
			F[m] = 0;
			molgroup_stats(mg);
			mg->box = mg->max - mg->min;
		}
		if ( verbose )
			cout << iter << tab << all_shift << endl;
	}

	delete[] F;

	return 0;
}
