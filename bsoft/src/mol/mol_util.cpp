/**
@file	mol_util.cpp
@brief	Library routines used for atomic coordinates
@author Bernard Heymann
@date	Created: 19980214
@date	Modified: 20220713
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "rwresprop.h"
#include "mol_util.h"
#include "seq_util.h"
#include "linked_list.h"
#include "Matrix.h"
#include "utilities.h"

#include <map>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the average Bfactors for residues in molecules.
@param 	*molgroup 	molecule group structure.
@return int 		0.
**/
int  		molgroup_Bfactors(Bmolgroup* molgroup)
{
	int				nat = 0;
	double			avgB = 0, sdB = 0, atomB = 0, dist, maxdist = 0;
	char			atomtype[10];
	Vector3<double>	ca;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Residue\tNatoms\tAvgB\tStDevB" << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			nat = 0;
			avgB = sdB = 0;
			for ( atom = res->atom; atom; atom = atom->next ) {
				nat++;
				avgB += atom->b;
				sdB += atom->b*atom->b;
			}
			avgB /= nat;
			sdB = sqrt(sdB/nat - avgB*avgB);
			cout << res->num << tab << nat << tab << avgB << tab << sdB << endl;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Residue\tMaxDist\tBfactor\tAtom" << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				if ( strstr(atom->type, "CA" ) ) {
					ca = atom->coord;
					strncpy(atomtype, atom->type, 8);
					atomB = atom->b;
				} else if ( strncmp(atom->type, "H", 1 ) &&
						strncmp(atom->type, " H", 2 ) &&
						!strstr(atom->type, "N " ) &&
						!strstr(atom->type, "C " ) &&
						!strstr(atom->type, "O " ) ) {
					dist = (atom->coord - ca).length();
					if ( maxdist < dist ) {
						maxdist = dist;
						strncpy(atomtype, atom->type, 8);
						atomB = atom->b;
					}
				}
			}
			cout << res->num << tab << maxdist << tab << atomB << tab << atomtype << endl;
			maxdist = 0;
		}
	}
	
	return 0;
}

/**
@brief 	Prints the sequence of a molecule as a single letter code string.
@param 	*molgroup 	molecule group structure.
@return int 		0.
**/
int 		molgroup_print_sequence(Bmolgroup* molgroup)
{
	Bmolecule*	mol = molgroup->mol;

	if ( mol->seq.length() < 1 )
		seq_from_residues(molgroup);

	if ( verbose & VERB_LABEL )
		cout << "Sequence:" << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		cout << "Molecule: " << mol->id << endl << mol->seq << endl;
		if ( verbose & VERB_PROCESS )
			cout << endl << "Length = " << mol->nres << endl << endl;
	}
	
	return 0;
}

/**
@brief	Shifts all coordinates to a given radius.
@param 	*molgroup 	molecule group structure.
@param 	radius		spherical radius.
@return int 		0.
**/
int 		molgroup_set_radius(Bmolgroup* molgroup, double radius)
{
	if ( radius < 0.1 ) return 0;
	
	long		natom;
	double		dist, sum, sumsq, avg, std;
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Shifting to radius:             " << radius << endl << endl;
	
    for ( sum=sumsq=natom=0, mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next, natom++ ) {
				dist = atom->coord.length();
				sum += dist;
				sumsq += dist*dist;
				atom->coord *= radius/dist;
			}
    	}
	}
	
	avg = sum/natom;
	std = sumsq/natom - avg*avg;
	if ( std > 0 ) std = sqrt(std);
	else std = 0;
	
	if ( verbose & VERB_PROCESS )
		cout << "Radius avg and std:             " << avg << " " << std << endl << endl;
	
	molgroup_stats(molgroup);
	
	return 0;
}

/**
@brief	Renames molecules.
@param 	*molgroup 	molecule group structure.
@param 	first_name	letter of first molecule.
@return int 		0.

	Letters are assigned to the molecules in sequence starting from the
	given first letter. If it extends beyond 'Z', it restarts at 'A'.

**/
int			molgroup_rename(Bmolgroup* molgroup, char first_name)
{
	if ( first_name == 0 ) return 0;
	
	char		letter = first_name;
	Bmolecule*	mol;
	
	for ( mol=molgroup->mol; mol; mol=mol->next ) {
		mol->id[0] = letter++;
		if ( letter > 'Z' ) letter = 'A';
	}
	
	return 0;
}

/**
@brief	Renumbers residues.
@param 	*molgroup 	molecule group structure.
@param 	first		number of first residue.
@return int 		0.
**/
int 		molgroup_residue_renumber(Bmolgroup* molgroup, int first)
{
    if ( first < 1 ) return 0;
	
    int     	resnum;
	Bmolecule*	mol;
	Bresidue*	res;
	
	if ( verbose & VERB_LABEL )
		cout << "Renumbering residues from:      " << first << endl;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		resnum = first;
		for( res = mol->res; res; res = res->next ) {
			res->num = resnum;
			res->insert[0] = ' ';
			resnum++;
		}
	}
	
	return 0;
}
    
/**
@brief	Renumbers all atoms, disregarding molecule distinction.
@param 	*molgroup 	molecule group structure.
@param 	first		number of first atom.
@return int 		0.
**/
int 		molgroup_atom_renumber(Bmolgroup* molgroup, int first)
{
    if ( first < 1 ) return 0;
	
    int     	atomnum = first;
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*		atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Renumbering atoms from:         " << first << endl;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->num = atomnum;
				atomnum++;
			}
		}
	}
	
	return 0;
}
    
/**
@brief	Sets occupancies to new values..
@param 	*molgroup 	molecule group structure.
@param 	range_first first residue in range to change.
@param 	range_last	last residue in range to change.
@param 	occupancy	value to set occupancy to.
@return int 		0.
**/
int 		molgroup_coor_reset_occupancy(Bmolgroup* molgroup, int range_first, 
				int range_last, double occupancy)
{
	Bmolecule*	mol;
	Bresidue*	res;
	Batom*  	atom;
	
	if ( verbose & VERB_LABEL )
		cout << "Reseting occupancy to range:    " << range_first << " " <<range_last << endl;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			if ( res->num >= range_first && res->num <= range_last )
			for ( atom = res->atom; atom; atom = atom->next ) {
    			atom->q = occupancy;
			}
		}
	}
	
	return 0;
}

/**
@brief	Calculates the weight of each molecule.
@param 	*molgroup 	molecule group structure.
@return double 		total mass of the molecule group.
**/
double 		molgroup_weight_from_atoms(Bmolgroup* molgroup)
{
	if ( !molgroup->mol ) return 0;
	if ( !molgroup->mol->res ) return 0;
	
	int 			n;
	double			summass(0), totmass(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

	if ( verbose & VERB_PROCESS )
		cout << "Molecule\tMW(Da)\tAtoms" << endl;
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		summass = 0;
		n = 0;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				n++;
				summass += atom->mass;
			}
	    }
		if ( verbose & VERB_PROCESS )
			cout << mol->id << tab << summass << tab << n << endl;
		totmass += summass;
	}
	
	if ( verbose ) {
		cout << "Total molecular weight:         " << totmass << " Da" << endl;
		cout << "Total volume:                   " << totmass/RHO << " A3" << endl << endl;
	}
	
	return totmass;
}

/**
@brief	Calculates the weight of each molecule.
@param 	*molgroup 	molecule group structure.
@return double 		total mass of the molecule group.
**/
double 		molgroup_weight_from_sequence(Bmolgroup* molgroup)
{
	if ( !molgroup->mol ) return 0;
	if ( !molgroup->mol->seq.length() ) return 0;

	Bstring			paramfile;
	Bresidue_type*	respar = get_residue_properties(paramfile);
	
	long 			n;
	double			summass(0), totmass(0);
	Bmolecule*		mol;
	Bresidue_type*	rt;

	if ( verbose & VERB_PROCESS )
		cout << "Molecule\tMW(Da)\tResidues" << endl;
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		summass = 0;
	    if ( mol->seq.length() ) {
			for ( n=0; n<mol->nres; ++n ) {
				for ( rt = respar; rt && rt->c != mol->seq[n]; rt = rt->next ) ;
				if ( rt ) summass += rt->mass - 18;
			}
	    }
		if ( verbose & VERB_PROCESS )
			cout << mol->id << tab << summass << tab << n << endl;
		totmass += summass;
	}
	
	if ( verbose ) {
		cout << "Total molecular weight:         " << totmass << " Da" << endl;
		cout << "Total volume:                   " << totmass/RHO << " A3" << endl << endl;
	}
	
	return totmass;
}

/**
@brief	Calculates the center of mass of a molecule group.
@param 	*molgroup 		molecule group structure.
@return Vector3<double> 	3-valued center of mass vector.
**/
Vector3<double> 	molgroup_center_of_mass(Bmolgroup* molgroup)
{
	double			summass(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	acoord, com, vec;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				acoord = atom->coord * atom->mass;
				com += acoord;
				vec += acoord*acoord;
				summass += atom->mass;
			}
	    }
	}
	
	if ( summass ) {
		com /= summass;
		vec = vec/summass - com*com;
		vec = vec.square_root();
		vec.normalize();
	} else cout << "Error: No atom masses found!" << endl << endl;
	
	if ( verbose & VERB_FULL ) {
		cout << "Total molecular weight:         " << summass << " Da" << endl;
		cout << "Center of mass:                 " << com << " A" << endl;
		cout << "Major axis:                     " << vec << endl;
	}
	
	return com;
}
	
/**
@brief	Calculates the center of mass of a molecule group.
@param 	*molgroup 		molecule group structure.
@return Vector3<double> 	3-valued center of mass vector.
**/
Vector3<double> 	molgroup_selected_center_of_mass(Bmolgroup* molgroup)
{
	double			summass(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	acoord, com, vec;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) if ( mol->sel ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) if ( atom->sel ) {
				acoord = atom->coord * atom->mass;
				com += acoord;
				vec += acoord*acoord;
				summass += atom->mass;
			}
	    }
	}
	
	if ( summass ) {
		com /= summass;
		vec = vec/summass - com*com;
		vec = vec.square_root();
		vec.normalize();
	} else cout << "Error: No atom masses found!" << endl << endl;
	
	if ( verbose & VERB_FULL ) {
		cout << "Total molecular weight:         " << summass << " Da" << endl;
		cout << "Center of mass:                 " << com << " A" << endl;
		cout << "Major axis:                     " << vec << endl;
	}
	
	return com;
}
	
/**
@brief	Calculates the center of mass of a molecule group.
@param 	*molgroup molecule group structure.
@return Vector3<double> 			3-valued center of mass vector.
**/
Vector3<double> 	molgroup_show_center_of_mass(Bmolgroup* molgroup)
{
	Vector3<double> 	com = molgroup_center_of_mass(molgroup);
	
	cout << "Center of mass:                 " << com << " A" << endl;
	
	return com;
}

/**
@brief	Calculates the center of mass of a molecule.
@param 	*mol		molecule structure.
@return Vector3<double> 			3-valued center of mass vector.
**/
Vector3<double> 	mol_center_of_mass(Bmolecule* mol)
{
	double			summass(0);
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	acoord, com, vec;
	
	for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			acoord = atom->coord * atom->mass;
			com += acoord;
			vec += acoord*acoord;
			summass += atom->mass;
	    }
	}
	
	if ( summass ) {
		com /= summass;
		vec = vec/summass - com*com;
		vec = vec.square_root();
		vec.normalize();
	} else cerr << "Error: No atom masses found!" << endl << endl;
	
	if ( verbose & VERB_FULL ) {
		cout << "Molecule:                       " << mol->id << endl;
		cout << "Molecular weight:               " << summass << " Da" << endl;
		cout << "Center of mass:                 " << com << " A" << endl;
		cout << "Major axis:                     " << vec << endl << endl;
	}
	
	return com;
}

/**
@brief	Calculates the center of mass of a molecule.
@param 	*mol		molecule structure.
@return Vector3<double> 			3-valued center of mass vector.
**/
Vector3<double> 	mol_show_center_of_mass(Bmolecule* mol)
{
	Vector3<double> 	com = mol_center_of_mass(mol);
	
	cout << "Molecule:                       " << mol->id << endl;
	cout << "Center of mass:                 " << com << " A" << endl;
	
	return com;
}

/**
@brief	Calculates the principal axes of a molecule group.
@param 	*molgroup			molecule group structure.
@param 	*eigenvec	eigen vectors (can be NULL).
@return Vector3<double>				3-valued vector of principal axes.
**/
Vector3<double> 	molgroup_principal_axes(Bmolgroup* molgroup, Vector3<double>* eigenvec)
{
	double			summass(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	eigenval, vec, vec2, vecx;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mol_principal_axes(mol, eigenvec);
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				vec += atom->coord * atom->mass;					// Sums
				vec2 += (atom->coord*atom->coord) * atom->mass;		// Square sums
				vecx[0] += atom->mass*atom->coord[0]*atom->coord[1];	// Cross-term sums
				vecx[1] += atom->mass*atom->coord[0]*atom->coord[2];
				vecx[2] += atom->mass*atom->coord[1]*atom->coord[2];
				summass += atom->mass;
			}
	    }
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Molecule group:" << endl;
		cout << "Molecular weight:               " << summass << " Da" << endl;
	}
	
	if ( summass < 1 ) {
		cout << "Error: No atom masses found!" << endl << endl;
		return vec;
	}
	
	vec /= summass;
	vec2 /= summass;
	vecx /= summass;
	
	return principal_axes(vec, vec2, vecx, eigenvec);
}

/**
@brief	Calculates the principal axes of a molecule.
@param 	*mol				molecule structure.
@param 	*eigenvec	eigen vectors (can be NULL).
@return Vector3<double>				3-valued vector of principal axes.
**/
Vector3<double> 	mol_principal_axes(Bmolecule* mol, Vector3<double>* eigenvec)
{
	double			summass(0);
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	eigenval, vec, vec2, vecx;
	
	for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			vec += atom->coord * atom->mass;					// Sums
			vec2 += (atom->coord*atom->coord) * atom->mass;		// Square sums
			vecx[0] += atom->mass*atom->coord[0]*atom->coord[1];	// Cross-term sums
			vecx[1] += atom->mass*atom->coord[0]*atom->coord[2];
			vecx[2] += atom->mass*atom->coord[1]*atom->coord[2];
			summass += atom->mass;
	    }
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Molecule:                       " << mol->id << endl;
		cout << "Molecular weight:               " << summass << " Da" << endl;
	}

	if ( summass < 1 ) {
		cout << "Error: No atom masses found!" << endl << endl;
		bexit(-1);
	}
	
	vec /= summass;
	vec2 /= summass;
	vecx /= summass;
	
	return principal_axes(vec, vec2, vecx, eigenvec);
}

/**
@brief	Calculates how close the coordinate set is to a spherical shape.
@param 	*molgroup 	molecule group structure.
@param	da			angular step size. (in radians)
@return double		sphericity.

	For each solid angle, the algorithm finds the coordinates with the
	maximum distance from the center-of-mass.

**/
double		molgroup_sphericity(Bmolgroup* molgroup, double da)
{
	if ( da > 1 ) da *= M_PI/180.0;
	
	long			i, j, k;
	double			t, p, dp, pm, st, dist, r, dm(0), avg(0), std(0), sph(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Vector3<double>	c;

	Vector3<double>	com = molgroup_center_of_mass(molgroup);
	
	if ( verbose )
		cout << "Calculating sphericity with a step size of "
			<< da*180.0/M_PI << " degrees" << endl << endl;

	vector<long>	n, na;
	
	for ( i = 0, t = 0; t <= M_PI; t += da ) {
		st = sin(t);
		if ( st > 1e-3 ) {
			dp = da/sin(t);
			pm = TWOPI - dp/2;
		} else {
			dp = 10;
			pm = TWOPI;
		}
		for ( j=0, p = 0; p<pm; p += dp, i++, j++ )
			if ( verbose & VERB_FULL )
				cout << i << tab << t << tab << p << endl;
		na.push_back(i);
		n.push_back(j);
		if ( verbose & VERB_FULL )
			cout << t << tab << n.back() << tab << na.back() << endl;
	}

	vector<long>	d(i, 0);
	
    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) {
				c = atom->coord - com;
				dist = c.normalize();
				if ( dm < dist ) dm = dist;
				t = acos(c[2]);
				j = long(t/da + 0.5);
				if ( j >= n.size() ) j = n.size() - 1;
				i = j;
				if ( i < 0 || i >= n.size() )
					cerr << "i = " << i << endl;
				k = 0;
				if ( j ) {
					i = na[j-1];
					r = sqrt(c[0]*c[0]+c[1]*c[1]);
					if ( r > da/2 ) {
						dp = da/r;
						p = atan2(c[1],c[0]);
//						if ( p < -dp ) p += TWOPI;
//						if ( p >= TWOPI - dp ) p -= TWOPI;
						if ( p < 0 ) p += TWOPI;
						if ( p >= TWOPI ) p -= TWOPI;
						k = long(p/dp);
						if ( k >= n[j] ) k -= n[j];
//						if ( k < 0 )
//							cerr << "k = " << k << endl;
						i += k;
					}
				}
				if ( i >= d.size() ) {
					cerr << "Error: array out of bounds!" << endl;
					cerr << i << tab << j << tab << k << endl;
				}
				if ( d[i] < dist ) d[i] = dist;
			}
			
	for ( i = 0; i < d.size(); ++i ) {
		for ( j = 0; j < na.size() && na[j] <= i; ++j ) ;
		if ( j ) k = i - na[j-1];
		else k = 0;
		t = da*j;
		st = sin(t);
		if ( st > 1e-3 )
			p = (da/st)*k;
		else
			p = 0;
		if ( verbose & VERB_PROCESS )
			cout << i << tab << j << tab << t << tab << p << tab << d[i] << endl;
		avg += d[i];
		std += d[i] * d[i];
	}
	
	avg /= d.size();
	std /= d.size();
	std = sqrt(std - avg*avg);
	sph = 1/(1+std);
	
	if ( verbose ) {
		cout << "Number of solid angle bins:     " << d.size() << endl;
		cout << "Maximal distance:               " << dm << endl;
		cout << "Average bin maximal distance:   " << avg << endl;
		cout << "Standard deviation:             " << std << endl;
		cout << "Sphericity:                     " << sph << endl << endl;
	}
	
	return sph;
}

/**
@brief	Calculates the number and mass density of a molecular system.
@param 	*molgroup 	molecule group structure.
@return double 		density: Da/A^3.

	Requirement: The box size must be correctly specified in the
	molecular group structure.
	The number of molecules, residues and atoms are calculated and the
	density for each number reported.
	The total mass is calculated and the mass density reported as
	Da/A^3 and g/ml.

**/
double		molgroup_density(Bmolgroup* molgroup)
{
	long			nmol, nres, natom;
	double			mass = 0;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
    for ( nmol=nres=natom=0, mol = molgroup->mol; mol; mol = mol->next, nmol++ )
		for( res = mol->res; res; res = res->next, nres++ )
			for ( atom = res->atom; atom; atom = atom->next, natom++ ) mass += atom->mass;
			
	double			volume = molgroup->box.volume();
	
	double			density = mass/volume;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Number of molecules and atoms:  " << nmol << " " << nres << " " << natom << endl;
		cout << "Volume:                         " << molgroup->box << " = " << volume << " A3" << endl;
		cout << "Density:                        " << density << " Da/A3" << density*1.66058 << " g/ml" << endl;
		cout << "Density:                        " << nmol*1.0/volume << " molecules/A3  "
			<< nres*1.0/volume << " residues/A3  " << natom*1.0/volume << " atoms/A3" << endl << endl;
	}
	
	return density;
}

/**
@brief	Calculates the volume of a molecular system.
@param 	*molgroup 	molecule group structure.
@param 	&paramfile		atomic parameter file.
@param 	wrap				flag to turn on wrapping.
@return double					volume in angstrom^3.

	The volume occupied by all the atoms is estimated from the footprint
	under the Van der Waals volume of each atom.

**/
double		molgroup_volume(Bmolgroup* molgroup, Bstring& paramfile, int wrap)
{
	Batomtype*		atompar = get_atom_properties(paramfile);
	Bresidue_type*	respar = get_residue_properties(paramfile);
	
//	write_atom_properties("at.star", atompar);
	
	if ( molgroup->box.volume() < 1 )
		molgroup->box = molgroup->max - molgroup->min;
	
	long			i, x, y, z, ix, iy, iz, rad;
	double			volume = 0, VdWvol = 0, mass = 0, resmass, sampling = 1, edge = 3, dist, dist2;
	double			voxel_volume = sampling*sampling*sampling;
	Vector3<int>	gridsize, lo, hi;
	Vector3<double>	c, d;
	if ( wrap ) edge = 0;
	gridsize[0] = (int) ((molgroup->max[0] - molgroup->min[0] + 2*edge)/sampling + 1);
	gridsize[1] = (int) ((molgroup->max[1] - molgroup->min[1] + 2*edge)/sampling + 1);
	gridsize[2] = (int) ((molgroup->max[2] - molgroup->min[2] + 2*edge)/sampling + 1);
//	Vector3<int>	h = gridsize * 0.5;
	long			size = (long) gridsize.volume();
	
	char*			grid = new char[size];
	for ( i=0; i<size; i++ ) grid[i] = 0;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Batomtype*  	at;
	Bresidue_type*	rt;
	
	cout << "Calculating the Van der Waals volume:" << endl;
	cout << "    in volume:                  " << gridsize << " = " << size << endl;
	if ( wrap ) cout << "with wrapping." << endl;
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom, resmass=0; atom; atom = atom->next ) {
				resmass += atom->mass;
				for ( at = atompar; at && strncmp(atom->el, at->el, 2)!=0; at = at->next ) ;
				dist = 2/sampling;	// Atom types not found tend to be larger
				if ( at ) {
					if ( at->vdw < 0.8 )
						cout << "VdW distance too short! atom=" << atom->el << " dist=" << at->vdw << endl;
					else
						dist = at->vdw/sampling;
				}
				dist2 = dist*dist;
				rad = (long) (dist + 1);
				c = atom->coord - molgroup->min;
				c += edge;
				c /= sampling;
				lo = Vector3<int>((int)c[0], (int)c[1], (int)c[2]) - rad;
				hi = Vector3<int>((int)c[0], (int)c[1], (int)c[2]) + rad;
				if ( !wrap ) {
					lo = lo.max(0);
					hi = hi.min(gridsize - 1);
				}
				for ( z=lo[2]; z<=hi[2]; z++ ) {
					iz = z;
					if ( iz < 0 ) iz += gridsize[2];
					if ( iz >= gridsize[2] ) iz -= gridsize[2];
					d[2] = c[2] - z;
					for ( y=lo[1]; y<=hi[1]; y++ ) {
						iy = y;
						if ( iy < 0 ) iy += gridsize[1];
						if ( iy >= gridsize[1] ) iy -= gridsize[1];
						d[1] = c[1] - iy;
						for ( x=lo[0]; x<=hi[0]; x++ ) {
							ix = x;
							if ( ix < 0 ) ix += gridsize[0];
							if ( ix >= gridsize[0] ) ix -= gridsize[0];
							d[0] = c[0] - ix;
							if ( d.length2() <= dist2 )
								grid[(iz*gridsize[1] + iy)*gridsize[0] + ix] = 1;
						}
					}
				}
			}
			mass += resmass;
			for ( rt = respar; rt && strncmp(rt->cod, res->type, 3)!=0; rt = rt->next ) ;
			if ( rt ) volume += rt->vol;
			else volume += resmass/RHO;
		}
	}
	
	cout << "Mass:                           " << mass << " Da" << endl;
	cout << "Volume from residue volumes:    " << volume << " A3 (" << 
		volume*100.0/molgroup->box.volume() << " %)" << endl;
	cout << "Density:                        " << mass/volume << " Da/A3 " 
		<< 1.66058*mass/volume << " g/ml" << endl << endl;

	for ( i=0, VdWvol=0; i<size; i++ ) VdWvol += grid[i];
	VdWvol *= voxel_volume;
	
	cout << "Van der Waals volume:           " << VdWvol << " A3 (" 
		<< VdWvol*100.0/(voxel_volume*size) << " %)" << endl;
	cout << "Packing density:                " << VdWvol/volume << endl << endl;
	
	delete[] grid;
	
	kill_list((char *) atompar, sizeof(Batomtype));
	kill_list((char *) respar, sizeof(Bresidue_type));
	
	return volume;
}

/**
@brief	Calculates the elemental composition of a molecular system.
@param 	*molgroup 	molecule group structure.
@param 	&paramfile		atomic parameter file.
@return int						0.
**/
int			molgroup_composition(Bmolgroup* molgroup, Bstring& paramfile)
{
	Batomtype*		atompar = get_atom_properties(paramfile);
	
	long			i, t, nel(0), nres(0), el[1024], an(0);
	double			mass[1024], mw(0), extinc(0);
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Batomtype*  	at;

	map<char, long>	num;
	
	for ( t=0, at = atompar; at && t < 1022; at = at->next, t++ )
		mass[t] = el[t] = 0;

    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		if ( mol->seq.length() ) {
			for ( i=0; i<mol->seq.length(); i++, nres++ ) {
 				if ( num.find(mol->seq[i]) == num.end() )
  					num[mol->seq[i]] = 1;
				else
					num[mol->seq[i]]++;
			}
		}
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next ) {
				el[atom->tnum]++;
				mass[atom->tnum] += atom->mass;
				mw += atom->mass;
				nel++;
			}
	}

	for ( t=0, at = atompar; at && t < 1022; at = at->next, t++ )
		an += el[t]*at->z;
	
	// From Mach, Middaugh and Lewis (1992) Analytical Biochemistry 200:74-80
	// From Pace et al. (1995) Protein Science 4:2411-2423.
	// Cystines: 125 /M/cm
	extinc += num['W'] * 5500;
	extinc += num['Y'] * 1490;
	// From different source
	extinc += num['F'] * 200;
	
	if ( nel ) {
		cout << "Elemental composition:" << endl;
		cout << "Element\tNumber\tMass" << endl;
		for ( t=0, at = atompar; at && t < 1022; at = at->next, t++ )
			if ( el[t] ) cout << at->el << tab << el[t] << tab << mass[t] << endl;
		cout << endl;
		cout << endl;
		cout << "Total molecular weight:         " << mw << " Da" << endl;
		cout << "Total atomic number:            " << an << endl << endl;
	}
	
	if ( nres ) {
		cout << "Residue composition:" << endl;
		cout << "Residue\tNumber" << endl;
		for ( auto it=num.begin(); it!=num.end(); ++it )
			cout << it->first << tab << it->second << endl;
		cout << endl;
		cout << "Extinction coefficient (280 nm): " << extinc << " /M/cm" << endl << endl;
	}
	
	kill_list((char *) atompar, sizeof(Batomtype));
	
	return 0;
}

/**
@brief 	Calculates the elemental composition.
@param 	*molgroup		molecule group.
@return JSvalue			composition.
**/
JSvalue		molgroup_elements(Bmolgroup* molgroup)
{
	JSvalue			comp(JSobject);

	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

    for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next )
				if ( comp.exists(atom->el) )
					comp[atom->el] += 1;
				else
					comp[atom->el] = 1;

//	cout << comp << endl;
	
	return comp;
}

/**
@brief 	Calculates the elemental composition from residue compositions.
@param 	*molgroup		molecule group.
@param 	&paramfile		residue parameter file.
@return JSvalue			composition.
**/
JSvalue		molgroup_elements(Bmolgroup* molgroup, Bstring& paramfile)
{
	long			i, j, n(0);
	Bmolecule*		mol;
	Bresidue_type*	rt_list = get_residue_properties(paramfile);
	Bresidue_type*	rt = NULL;
	vector<double>	el(5,0), eltot(5,0);

	if ( verbose & VERB_LABEL )
		cout << "Calculating the elemental composition:" << endl << endl;
	
	if ( verbose )
		cout << "Molecule\tH\tC\tN\tO\tS\tTotal" << endl;
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
//		cout << mol->seq << endl;
		for ( i=0; i<mol->seq.length(); i++ ) {
			for ( rt = rt_list; rt && rt->c != mol->seq[i]; rt = rt->next ) ;
			if ( rt )
				for ( j=0; j<5; ++j )
					el[j] += rt->comp[j];
		}
		if ( verbose )
			cout << mol->id;
		for ( j=n=0; j<5; ++j ) {
			n += el[j];
			eltot[j] += el[j];
			if ( verbose )
				cout << tab << el[j];
			el[j] = 0;
		}
		if ( verbose )
			cout << tab << n << endl;
	}
	if ( verbose )
		cout << "Total";
	for ( j=n=0; j<5; ++j ) {
		n += eltot[j];
		if ( verbose )
			cout << tab << eltot[j];
		el[j] = 0;
	}
	if ( verbose )
		cout << tab << n << endl;

	kill_list((char *) rt_list, sizeof(Bresidue_type));
	
	if ( verbose ) {
		cout << "%";
		for ( j=0; j<5; ++j )
			cout << tab << eltot[j]*100.0/n;
		cout << endl;
	}
	
	JSvalue			comp(JSobject);

	comp["H"] = eltot[0];
	comp["C"] = eltot[1];
	comp["N"] = eltot[2];
	comp["O"] = eltot[3];
	comp["S"] = eltot[4];
	
	if ( verbose & VERB_PROCESS )
		cout << comp << endl;

	return comp;
}

/**
@brief 	Default protein composition adjusted by mass.
@param 	mass			molecular weight.
@return JSvalue			composition.
**/
JSvalue		protein_composition_default(double mass)
{
	JSvalue			comp(JSobject);
	
	comp["H"] = 0.498;
	comp["C"] = 0.320;
	comp["N"] = 0.085;
	comp["O"] = 0.095;
	comp["S"] = 0.002;
	
	double			mc = comp["H"].real() +
						12*comp["C"].real() +
						14*comp["N"].real() +
						16*comp["O"].real() +
						32*comp["S"].real();

	double			sc = mass/mc;
//	cout << mass << tab << mc << tab << sc << endl;
	
	for ( auto it = comp.object_begin(); it != comp.object_end(); ++it )
		it->second *= sc;

	if ( verbose & VERB_PROCESS )
		cout << comp << endl;
	
	return comp;
}

/**
@brief 	Calculates a radial distribution function.
@param 	*molgroup		molecule group.
@param 	interval		interval between bins.
@param 	cutoff			distance cutoff.
@param 	wrap			flag to use periodic boundaries.
@return int				0.

	The atoms are not distinguished by type or properties.

**/
int			molgroup_radial_density(Bmolgroup* molgroup, double interval, double cutoff, int wrap)
{
	int				i, ii, j, x, y, z, xx, yy, zz, ix, iy, iz;
	int				xlo, ylo, zlo, xhi, yhi, zhi;
	Vector3<double>	box = molgroup->box;
	double			dist;
	Vector3<double>	d, sampling(cutoff, cutoff, cutoff);
	Vector3<int>	size((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/size[i] + 0.001;

	Latom**			alist = molgroup_atom_mesh_lists(molgroup, size, sampling);
	Latom			*latom, *latom2;
	
	double			mult = 1.0/interval;
	int				n = (int) (mult*cutoff+1);
	int*			rdf = new int[n];
	for ( i=0; i<n; i++ ) rdf[i] = 0;
	
	for ( i=z=0; z<size[2]; z++ ) {
		zlo = z - 1;
		zhi = z + 1;
		if ( !wrap ) {
			if ( zlo < 0 ) zlo = 0;
			if ( zhi >= size[2] ) zhi = size[2] - 1;
		}
		for ( y=0; y<size[1]; y++ ) {
			ylo = y - 1;
			yhi = y + 1;
			if ( !wrap ) {
				if ( ylo < 0 ) ylo = 0;
				if ( yhi >= size[1] ) yhi = size[1] - 1;
			}
			for ( x=0; x<size[0]; x++, i++ ) {
				xlo = x - 1;
				xhi = x + 1;
				if ( !wrap ) {
					if ( xlo < 0 ) xlo = 0;
					if ( xhi >= size[0] ) xhi = size[0] - 1;
				}
				for ( latom = alist[i]; latom; latom = latom->next ) {
					for ( zz=zlo; zz<=zhi; zz++ ) {
						iz = zz;
						if ( wrap ) {
							if ( iz < 0 ) iz += size[2];
							if ( iz >= size[2]) iz -= size[2];
						}
						for ( yy=ylo; yy<=yhi; yy++ ) {
							iy = yy;
							if ( wrap ) {
								if ( iy < 0 ) iy += size[1];
								if ( iy >= size[1] ) iy -= size[1];
							}
							for ( xx=xlo; xx<=xhi; xx++ ) {
								ix = xx;
								if ( wrap ) {
									if ( ix < 0 ) ix += size[0];
									if ( ix >= size[0] ) ix -= size[0];
								}
								ii = (iz*size[1] + iy)*size[0] + ix;
								for ( latom2 = alist[ii]; latom2; latom2 = latom2->next ) 
										if ( latom2->atom != latom->atom ) {
									d = vector3_difference_PBC(latom2->atom->coord, latom->atom->coord, box);
									dist = d.length();
									if ( dist <= cutoff ) {
										j = (int) (mult*dist);
										if ( j < n ) rdf[j]++;
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
	
	cout << "Calculating the radial distribution function within " << cutoff << " A" << endl;
	
	cout << "r(A)\tRDF" << endl;
	for ( i=0; i<n; i++ ) cout << i*interval << tab << rdf[i] << endl;
	cout << endl;
	
	delete[] rdf;
	
	return 0;
}

/**
@brief	Generates lists of atoms based on a grid.

	The coordinates must be positive to fit into a grid starting at {0,0,0}.

@param 	*molgroup 	molecule group structure to be modified.
@param 	size		size of grid.
@param 	sampling	spacing in each dimension.
@return Latom**					array of linked lists of atoms.
**/
Latom**		molgroup_atom_mesh_lists(Bmolgroup* molgroup, Vector3<int> size, Vector3<double> sampling)
{
	long			i, x, y, z;
	
	Vector3<double>	origin = molgroup->min;

	long			gridvol = (long) size.volume();
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Latom**			alist = new Latom*[gridvol];
	Latom**			alistp = new Latom*[gridvol];
	
	for ( i=0; i<gridvol; i++ ) alist[i] = alistp[i] = NULL;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				x = (long) ((atom->coord[0] - origin[0])/sampling[0]);
				y = (long) ((atom->coord[1] - origin[1])/sampling[1]);
				z = (long) ((atom->coord[2] - origin[2])/sampling[2]);
				if ( x < 0 ) x = 0;
				if ( y < 0 ) y = 0;
				if ( z < 0 ) z = 0;
				if ( x >= size[0] ) x = size[0] - 1;
				if ( y >= size[1] ) y = size[1] - 1;
				if ( z >= size[2] ) z = size[2] - 1;
				i = (z*size[1] + y)*size[0] + x;
				alistp[i] = (Latom *) add_item((char **) &alistp[i], sizeof(Latom));
				alistp[i]->atom = atom;
				if ( !alist[i] ) alist[i] = alistp[i];
			}
		}
	}
	
	delete[] alistp;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_atom_mesh_lists: Done!" << endl;
	
	return alist;
}

/**
@brief	Generates an array of pointers to atom structures.
@param 	*molgroup 		molecule group structure to be modified.
@return vector<Batom*>	array of pointers to atoms.
**/
vector<Batom*>	atom_get_array(Bmolgroup* molgroup)
{
	long			n;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	vector<Batom*>	atomarray;

    for ( n=0, mol = molgroup->mol; mol; mol = mol->next ) if ( mol->sel )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom; atom = atom->next, n++ )
				atomarray.push_back(atom);


	return atomarray;
}

