/**
@file	mol_map_energy.cpp
@brief	Functions to calculate estimates of the fitting of molecules to maps.
@author Bernard Heymann
@date	Created: 20041230
@date	Modified: 20111030
**/

#include "rwmolecule.h"
#include "rwimg.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Energy calculation of a molecule fit into a map.  
@param 	*molgroup		molecule group.
@param 	*map			map.
@param 	Kmap			map energy constant.
@return double			energy.
**/
double		molgroup_map_energy(Bmolgroup* molgroup, Bimage* map, double Kmap)
{
	map->change_type(Float) ;
	
	Bmolgroup*		mg;
	Bmolecule*		mol;	
	Bresidue*		res;
	Batom*			atom;
	long			i, n(0), nmg, nm;
	long			x, y, z;
	long			slicesize = map->sizeX()*map->sizeY();
	double			E(0), Emg, Em;
	double			invunit = 0.5*Kmap/map->sampling(0)[0];
	Vector3<double>	F;
	
	for ( Emg=0, nmg=0, mg = molgroup; mg; mg = mg->next ) {
		for ( Em=0, nm=0, mol = mg->mol; mol; mol = mol->next ) {
			for( res = mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					x = (long) (atom->coord[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
					y = (long) (atom->coord[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
					z = (long) (atom->coord[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
					if ( x>0 && x<map->sizeX()-1 && y>0 && y<map->sizeY()-1 && z>0 && z<map->sizeZ()-1 ) {
						i = map->index(0, x, y, z, 0);
						F = Vector3<double>((*map)[i+1] - (*map)[i-1], (*map)[i+map->sizeX()] - (*map)[i-map->sizeX()], 
							(*map)[i+slicesize] - (*map)[i-slicesize]);
						Em -= (*map)[i];
						atom->F += F * invunit;
					}
					nm++;
				}
			}
			if ( nm ) mol->fom = (Em/nm + map->average())/map->standard_deviation() + 10;
			else mol->fom = 100;
			Emg += Em;
			nmg += nm;
		}
		if ( nmg ) mg->fom = (Emg/nmg + map->average())/map->standard_deviation() + 10;
		else mg->fom = 100;
		E += Emg;
		n += nmg;
	}
	
	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Molecule group outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap*((E/n + map->average())/map->standard_deviation() + 10);
	
	return E;
}

/**
@brief 	Energy calculation of a molecule fit into a map.  
@param 	*mol			molecule.
@param 	*map			map.
@param 	Kmap			map energy constant.
@return double			energy.
**/
double		mol_map_energy(Bmolecule* mol, Bimage* map, double Kmap)
{
	map->change_type(Float) ;
	
	Bresidue*		res;
	Batom*			atom;
	long			i, n(0);
	long			x, y, z;
	long			slicesize = map->sizeX()*map->sizeY();
	double			E(0);
	double			invunit = 0.5*Kmap/map->sampling(0)[0];
	Vector3<double>	F;
	
    for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			x = (long) (atom->coord[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
			y = (long) (atom->coord[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
			z = (long) (atom->coord[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
			if ( x>0 && x<map->sizeX()-1 && y>0 && y<map->sizeY()-1 && z>0 && z<map->sizeZ()-1 ) {
				i = map->index(0, x, y, z, 0);
				F = Vector3<double>((*map)[i+1] - (*map)[i-1], (*map)[i+map->sizeX()] - (*map)[i-map->sizeX()], 
					(*map)[i+slicesize] - (*map)[i-slicesize]);
				E -= (*map)[i];
				atom->F += F * invunit;
			}
			n++;
		}
    }
	
	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Molecule outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap*((E/n + map->average())/map->standard_deviation() + 10);
	
	return E;
}

/**
@brief 	Calculation of the fit of a molecule group into a map as a correlation coefficient.  
@param 	*molgroup		molecule group.
@param 	*map			map.
@return double			correlation coefficient.
**/
double		molgroup_map_correlation(Bmolgroup* molgroup, Bimage* map)
{
	map->change_type(Float) ;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	long			i, n(0);	
	long			datasize = map->data_size();
	long			x, y, z;
	double			sumc2(0), sumd2(0), sumcd(0);
	map->next = new Bimage(Float, TSimple, map->size(), map->images());
	float*			calc = (float *) map->next->data_pointer();
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				x = (long) (atom->coord[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
				y = (long) (atom->coord[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
				z = (long) (atom->coord[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
				if ( x>=0 && x<map->sizeX() && y>=0 && y<map->sizeY() && z>=0 && z<map->sizeZ() ) {
					i = map->index(0, x, y, z, 0);
					calc[i] += atom->mass;
				}
			}
		}
	}
	
	for ( i=0; i<datasize; i++ ) {
		if ( calc[i] ) {
			sumc2 += calc[i]*calc[i];
			sumd2 += (*map)[i]*(*map)[i];
			calc[i] *= (*map)[i];
			sumcd += calc[i];
			n++;
		}
	}
	
	double		den = sqrt(sumc2*sumd2);
	double		CC = sumcd/den;
	double		norm = n/den;

	for ( i=0; i<datasize; i++ )
		if ( calc[i] ) calc[i] *= norm;
	
	if ( verbose & VERB_FULL )
		cout << "Correlation:                    " << CC << " (" << n << ")" << endl << endl;
	
	return CC;
}

/**
@brief 	Calculation of the fit of a molecule into a map as a correlation coefficient.  
@param 	*mol			molecule.
@param 	*map			map.
@return double			correlation coefficient.
**/
double		mol_map_correlation(Bmolecule* mol, Bimage* map)
{
	map->change_type(Float) ;
	
	Bresidue*		res;
	Batom*			atom;
	long			i, n(0);	
	long			datasize = map->data_size();
	long			x, y, z;
	double			sumc2(0), sumd2(0), sumcd(0);
	map->next = new Bimage(Float, TSimple, map->size(), map->images());
	float*			calc = (float *) map->next->data_pointer();
	
    for( res = mol->res; res; res = res->next ) {
		for ( atom = res->atom; atom; atom = atom->next ) {
			x = (long) (atom->coord[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
			y = (long) (atom->coord[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
			z = (long) (atom->coord[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
			if ( x>=0 && x<map->sizeX() && y>=0 && y<map->sizeY() && z>=0 && z<map->sizeZ() ) {
				i = map->index(0, x, y, z, 0);
				calc[i] += atom->mass;
			}
		}
	}
	
	for ( i=0; i<datasize; i++ ) {
		if ( calc[i] ) {
			sumc2 += calc[i]*calc[i];
			sumd2 += (*map)[i]*(*map)[i];
			calc[i] *= (*map)[i];
			sumcd += calc[i];
			n++;
		}
	}
	
	double		den = sqrt(sumc2*sumd2);
	double		CC = sumcd/den;
	
	if ( verbose & VERB_FULL )
		cout << "Correlation for " << mol->id << ":      " << CC << " (" << n << ")" << endl << endl;
	
	return CC;
}

/**
@brief 	Energy calculation of a molecule fit into a map.  
@param 	*molgroup		molecule group.
@param 	*map			map.
@param 	Kmap			map energy constant.
@param 	steps			number of steps along a bond.
@return double			energy.


	The energy is the negative of the density at each step along a bond, rho, 
	plus a fudge factor to make it positive:
		E = Kmap * ((-sum(rho)/n + avg)/std + 10)
	where n is the number of voxels sampled, avg is the average of the map,
	and std is the standard deviation of the map.
	The force associated with this energy is the gradient at each step:
		Fx = Kmap * (rho(x+1) - rho(x-1))/(2*u)
	where u is the voxel size. The contribution of a force to an atom at 
	a step is weighted by the fractional distance the step is away from
	the atom along the bond. 

**/
double		molgroup_bond_fit_map_energy(Bmolgroup* molgroup, Bimage* map, double Kmap, int steps)
{
	map->change_type(Float) ;
	
	Bbond*			bond;
	Batom			*atom1, *atom2;
	long			i, n(0);
	long			slicesize = map->sizeX()*map->sizeY();
	long			x, y, z;
	double			frac, df, frac1, s1, s2, E(0);
	double			invunit = 0.5*Kmap/map->sampling(0)[0];
	Vector3<double>	F, F1, F2, loc;
	
	for ( bond = molgroup->bond; bond; bond = bond->next ) {
		atom1 = bond->atom1;
		atom2 = bond->atom2;
		df = 1.0L/steps;
		s1 = s2 = 0;
		F1 = F2 = 0;
		for ( frac=0; frac<=1; frac+=df ) {
			frac1 = 1 - frac;
			loc = atom1->coord * frac1 + atom2->coord * frac;
			x = (long) (loc[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
			y = (long) (loc[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
			z = (long) (loc[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
			if ( x>0 && x<map->sizeX()-1 && y>0 && y<map->sizeY()-1 && z>0 && z<map->sizeZ()-1 ) {
				i = map->index(0, x, y, z, 0);
				F = Vector3<double>((*map)[i+1] - (*map)[i-1], (*map)[i+map->sizeX()] - (*map)[i-map->sizeX()], 
					(*map)[i+slicesize] - (*map)[i-slicesize]);
				F1 += F * frac1;
				F2 += F * frac;
				if ( !isfinite((*map)[i]) ) cout << "Infinite value at " << x << " " << y << " " << z << endl;
				E -= (*map)[i];
				s1 += frac1;
				s2 += frac;
			}
			n++;
		}
		atom1->F += F1 * (invunit/s1);
		atom2->F += F2 * (invunit/s2);
	}

	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Molecule outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap * ((E/n + map->average())/map->standard_deviation() + 10);
	
	return E;
}



