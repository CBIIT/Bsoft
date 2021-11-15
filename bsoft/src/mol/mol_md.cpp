/**
@file	mol_md.cpp
@brief	Functions for molecular dynamics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20110816
**/

#include "mol_md.h"
#include "mol_bonds.h"
#include "rwmolecule.h"
#include "random_numbers.h"
#include "Matrix.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Molecular dynamics using the leapfrog integrator.
@param 	*molgroup		molecular structure.
@param 	*md				molecular dynamics parameters.
@param 	max_iter		maximum number of iterations to run.
@param 	velocitylimit	limit on velocity per time step.
@return double 			energy.

	Leapfrog integration for any coordinate x, velocity vx and force Fx:
		x(t+1) = x(t) + vx(t+1) * dt
		vx(t+1) = (Fx(t) * dt/m + vx(t)) * kf
		where
			kf:	friction constant (1=no friction)
			dt:	time step
			m: atomic mass
	The velocity is limited each time step to damp chaotic oscillations.

**/
double		md_leapfrog(Bmolgroup* molgroup, Bmd* md, int max_iter, double velocitylimit)
{
	if ( max_iter < 1 ) return 1e10;
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Bbond*			bondlist = molgroup->bond;
	Bangle*			anglelist = molgroup->angle;
		
	int 			cycle = 0;
	double			E(0), dE(0), vel;
	double			timestep = md->timestep;
	double			Kfriction = md->Kfriction;
	double			minvel = -velocitylimit, maxvel = velocitylimit;	// Maximum velocity per time step
	
	random_seed();
	
	if ( !bondlist ) {
		error_show("Error in md_leapfrog: No bond list!", __FILE__, __LINE__);
		return 1e10;
	}
	
	if ( !anglelist && md->Kangle ) {
		error_show("Error in md_leapfrog: No angle list!", __FILE__, __LINE__);
		return 1e10;
	}
	
	double			volume = molgroup->box.volume();
	if ( volume < 1 ) {
		molgroup->box = molgroup->max;
		molgroup->box[0] += 1.7;
		molgroup->box[1] += 1.7;
		molgroup->box[2] += 1.7;
		volume = molgroup->box.volume();
	}
	
	if ( verbose ) {
		cout << "Molecular dynamics using the leapfrog integrator:" << endl;
		cout << "Periodic box:                   " << molgroup->box << " " << volume << " A3" << endl;
		cout << "Maximum iterations:             " << max_iter << endl;
		cout << "Time step:                      " << timestep << endl;
		cout << "Kfriction:                      " << Kfriction << endl;
		cout << "Kbond:                          " << md->Kbond << endl;
		cout << "Kangle:                         " << md->Kangle << endl;
		cout << "Kelectrostatic:                 " << md->Kelec << endl;
		cout << "KVanderWaals:                   " << md->Kvdw << endl;
		cout << "Non-bonded cutoff distance:     " << md->cutoff << " A" << endl;
	
		cout << endl << "Cycle\tEbond\tEangle\tEelec\tEvdw\tEpoint\tEnergy\tdE\tEkin" << endl;
	}
	
//	while ( cycle < max_iter && fabs(dE) > 1e-37 ) {
	while ( cycle < max_iter && isfinite(E) ) {
		
		dE = E; 			// Save the old energy value
			// Initialize for the new energy values
		E = md->Ebond = md->Eangle = md->Eelec = md->Evdw = md->Ekin = md->Epoint = 0;
				
		// Zero the force on each atom
		md_zero_forces(molgroup);
		
		// Bond forces
		md->Ebond = md_bond_forces(molgroup, md->Kbond, md->wrap);
		
		// Angle forces
		md->Eangle = md_angular_forces(molgroup, md->Kangle, md->wrap);
		
		// Nonbonded forces
		md_nonbonded_forces(molgroup, md);

		md->Epoint = md_point_force(molgroup, md->point, md->Kpoint, md->pointdecay);
		
		// Calculate the total potential energy
		E = md->Ebond + md->Eangle + md->Eelec + md->Evdw + md->Epoint;
		dE = E - dE;		// Get the change in potential energy for this cycle
		
		if ( dE <= 0 ) Kfriction = md->Kfriction;
		else Kfriction = md->Kfriction*exp(-dE);
		
		// Calculate new atomic coordinates using the leapfrog method and the kinetic energy
		for ( mol = molgroup->mol; mol; mol = mol->next ) {
			for( res = mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					atom->vel += atom->F * (timestep/atom->mass);
					atom->vel *= Kfriction;
					atom->vel = vector3_scalar_range(atom->vel, minvel, maxvel);
					atom->coord += atom->vel * timestep;
					if ( md->wrap ) atom->coord = vector3_set_PBC(atom->coord, molgroup->box);
					vel = atom->vel.length();
					md->Ekin += 0.5*atom->mass*vel*vel;
				}
			}
		}
		
		// Show the current energies
		cycle++;
		
		if ( verbose )
			cout << cycle << tab << md->Ebond << tab << md->Eangle << tab << 
				md->Eelec << tab << md->Evdw << tab << md->Epoint << tab << 
				E << tab << dE << tab << md->Ekin << endl;
		else
			cout << cycle << tab << E << endl;
			
		if ( !isfinite(E) ) {
			cerr << endl << "***** Warning: The system exploded!!! *****" << endl << endl;
			return E;
		}

	}
	
	if ( verbose )
		cout << endl;

	return E;
}

/**
@brief 	Zero all atomic forces.
@param 	*molgroup		molecular structure.
@return int				0.
**/
int			md_zero_forces(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	for ( mol = molgroup->mol; mol; mol = mol->next )
		for( res = mol->res; res; res = res->next )
			for ( atom = res->atom; atom->next; atom = atom->next )
				atom->F = 0;
	
	return 0;
}

/**
@brief 	Calculates the covalent bond length forces and energy.
@param 	*molgroup		molecular structure.
@param 	Kbond			bond energy constant.
@param 	wrap			flag to wrap around periodic boundaries.
@return double			total bond length energy.

	The energy is defined as a harmonic function around the reference 
	bond length, |ro|:
		Eb = Kb*(|r|-|ro|)^2
	The force is the derivative of the energy:
		Fb = -2*Kb*(|r|-|ro|)*r/|r|
	where r is the distance vector and Kb is the bond energy constant.

**/
double		md_bond_forces(Bmolgroup* molgroup, double Kbond, int wrap)
{
	if ( Kbond <= 0 ) return 0;
	
	double			dist, dev, fac, energy(0);
	Vector3<double>	d, Force;
	Bbond*			bond;
	Bbond*			bondlist = molgroup->bond;
	
	for ( bond = bondlist; bond; bond = bond->next ) {
		if ( wrap )
			d = vector3_difference_PBC(bond->atom1->coord, bond->atom2->coord, molgroup->box);
		else
			d = bond->atom1->coord - bond->atom2->coord;
		dist = d.length();
		dev = dist - bond->l;
		energy += Kbond*dev*dev;
		fac = -2*Kbond*dev/dist;
		Force = d * fac;
		bond->atom1->F += Force;
		bond->atom2->F -= Force;
	}

	return energy;
}

/**
@brief 	Calculates the covalent bond angular forces and energy.
@param 	*molgroup		molecular structure.
@param 	Kangle			bond angle energy constant.
@param 	wrap			flag to wrap around periodic boundaries.
@return double			total bond angle energy.

	The energy is defined as a harmonic function around the reference 
	bond angle, a0:
		Ea = Ka*(cos(a0)-r1*r2/(|r1|*|r2|))^2
	The force is the derivative of the energy on the first and last atoms:
		Fa1 = 2*Ka*(cos(a0)-r1*r2/(|r1|*|r2|))/(|r1|*|r2|) * ((r1*r2/|r1|)*r1-r2)
		Fa3 = 2*Ka*(cos(a0)-r1*r2/(|r1|*|r2|))/(|r1|*|r2|) * ((r1*r2/|r2|)*r2-r1)
	where r1 is the vector from atom 2 to atom 1, r2 is the vector from
	atom 2 to atom 3, and Ka is the bond angle energy constant.

**/
double		md_angular_forces(Bmolgroup* molgroup, double Kangle, int wrap)
{
	if ( Kangle <= 0 ) return 0;
	
	double			d1len2, d2len2, dot, fac, dcosa, c11, c12, c22, energy(0);
	Vector3<double>	d1, d2, Force;
	Bangle*			angle;
	
	for ( angle = molgroup->angle; angle; angle = angle->next ) {
		if ( wrap ) {
			d1 = vector3_difference_PBC(angle->atom2->coord, angle->atom1->coord, molgroup->box);
			d2 = vector3_difference_PBC(angle->atom2->coord, angle->atom3->coord, molgroup->box);
		} else {
			d1 = angle->atom2->coord - angle->atom1->coord;
			d2 = angle->atom2->coord - angle->atom3->coord;
		}
		d1len2 = d1.length2();
		d2len2 = d2.length2();
		dot = d1.scalar(d2);
		fac = 1/sqrt(d1len2*d2len2);
		dcosa = cos(angle->a) - dot*fac;
		energy += Kangle*dcosa*dcosa;
		c12 = 2*Kangle*fac*dcosa;
		c11 = c12*dot/d1len2;
		c22 = c12*dot/d2len2;
		Force = (d1 * c11) - (d2 * c12);
		angle->atom1->F += Force;
		Force = (d2 * c22) - (d1 * c12);
		angle->atom3->F += Force;
	}

	return energy;
}

/**
@brief 	Calculates the non-bonded forces and energy.
@param 	*molgroup		molecular structure.
@param 	*md				molecular dynamics structure.
@return double			total non-bonded energy.

	The energy is defined as a Lennard-Jones term for the Van der Waals 
	interactions based on a reference length, ro, and a Coulomb term for 
	electrostatic interactions based on the atomic charges, q1 and q2: 
		Enb = Kvdw*((1/12)*(|ro|/|r|)^12 - (1/6)*(|ro|/|r|)^6) + Kelec*q1*q2/|r|
	The force is the derivative of the energy:
		Fnb = -Kvdw*((|ro|/|r|)^12 - (|ro|/|r|)^6)*r/|r| + Kelec*q1*q2*r/|r|^3
	where r is the distance vector, Kvdw is the Van der Waals energy constant
	and Kelec is the electrostatic energy constant.

**/
double		md_nonbonded_forces(Bmolgroup* molgroup, Bmd* md)
{
	if ( md->Kvdw < 0 ) md->Kvdw = 0;
	if ( md->Kelec < 0 ) md->Kelec = 0;
	if ( md->Kvdw <= 0 && md->Kelec <= 0 ) return 0;
	
	int				dononbond;
	long			i, ii, x, y, z, xx, yy, zz, ix, iy, iz;
	Vector3<double>	box = molgroup->box;
	Vector3<double>	sampling(md->cutoff, md->cutoff, md->cutoff);
	Vector3<int>	size((int) (box[0]/sampling[0] + 0.001), 
		(int) (box[1]/sampling[1] + 0.001), (int) (box[2]/sampling[2] + 0.001));
	size = size.max(1);
	for ( i=0; i<3; i++ ) sampling[i] = box[i]/size[i] + 0.001;
	long			boxsize = (long) size.volume();
	Latom			*latom, *latom2;
	
	Latom**			alist = molgroup_atom_mesh_lists(molgroup, size, sampling);
	
	md->Evdw = md->Eelec = 0;
	
	for ( z=0; z<size[2]; z++ ) {
		for ( y=0; y<size[1]; y++ ) {
			for ( x=0; x<size[0]; x++ ) {
				i = (z*size[1] + y)*size[0] + x;
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
								for ( latom2 = alist[ii]; latom2; latom2 = latom2->next ) {
									dononbond = 1;
									if ( latom2->atom == latom->atom ) dononbond = 0;	// Same atom
									else if ( latom2->atom->next == latom->atom ) dononbond = 0;	// Bonded atoms
									else if ( latom2->atom == latom->atom->next ) dononbond = 0;	// Bonded atoms
									else if ( latom2->atom->next ) {	// Atom 2 bonds away
										if ( latom2->atom->next->next == latom->atom ) dononbond = 0;	
									} else if ( latom->atom->next ) {	// Atom 2 bonds away
										if ( latom2->atom == latom->atom->next->next ) dononbond = 0;
									}
									if ( dononbond ) {
										atom_nonbonded_forces(latom->atom, latom2->atom, md, box);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for ( i=0; i<boxsize; i++ ) kill_list((char *) alist[i], sizeof(Latom));
	delete[] alist;
	
	return md->Eelec+md->Evdw;
}

/**
@brief 	Calculates the non-bonded forces and energy between two atoms.
@param 	*atom			central atom.
@param 	*atom2			second atom.
@param 	*md				molecular dynamics structure.
@param 	box				dynamics box.
@return double			total non-bonded energy.

	The energy is defined as a Lennard-Jones term for the Van der Waals 
	interactions based on a reference length, ro, and a Coulomb term for 
	electrostatic interactions based on the atomic charges, q1 and q2: 
		Enb = Kvdw*((1/12)*(|ro|/|r|)^12 - (1/6)*(|ro|/|r|)^6) + Kelec*q1*q2/|r|
	The force is the derivative of the energy:
		Fnb = -Kvdw*((|ro|/|r|)^12 - (|ro|/|r|)^6)*r/|r| + Kelec*q1*q2*r/|r|^3
	where r is the distance vector, Kvdw is the Van der Waals energy constant
	and Kelec is the electrostatic energy constant.

**/
int			atom_nonbonded_forces(Batom* atom, Batom* atom2, Bmd* md, Vector3<double> box)
{
	if ( md->Kvdw < 0 ) md->Kvdw = 0;
	if ( md->Kelec < 0 ) md->Kelec = 0;
	if ( md->Kvdw <= 0 && md->Kelec <= 0 ) return 0;
	
	Vector3<double>	d, Force;
	double			rd2, rd6, rd12, dist2, invdist2, fac;
	Bbondtype*		bt = NULL;
	
	if ( md->wrap )
		d = vector3_difference_PBC(atom2->coord, atom->coord, box);
	else
		d = atom2->coord - atom->coord;
	dist2 = d.length2();
	if ( dist2 < md->cutoff*md->cutoff ) {
		invdist2 = 1/dist2;
		if ( md->Kvdw ) {
			bt = md_find_bond_type(atom, atom2, md->bond);
			if ( bt && bt->vdwdist > 0 ) {
				rd2 = bt->vdwdist*bt->vdwdist*invdist2;
				rd6 = rd2*rd2*rd2;
				rd12 = rd6*rd6;
				md->Evdw += md->Kvdw*(md->VdWcoeff1*rd12 - md->VdWcoeff2*rd6);
				fac = -md->Kvdw*(rd12 - rd6)*invdist2;
				Force = d * fac;
				atom->F += Force;
			}
		}
		if ( md->Kelec ) {
			fac = md->Kelec*atom->chrg*atom2->chrg*sqrt(invdist2);
			md->Eelec += fac;
			fac *= invdist2;
			Force = d * fac;
			atom->F += Force;
		}
	}
	
	return 0;
}

/**
@brief 	Calculates the atomic forces and energy resulting from a single point force.
@param 	*molgroup		molecular structure.
@param 	point			center of point force.
@param 	Kpoint			point force constant.
@param 	decay			energy decay with distance.
@return double			point force energy.

	The energy is defined as an exponential decay over distance from the 
	center of the point force:
		Ep = Kp * exp(-decay*dist)
	The force is the derivative of the energy:
		Fp = Kp * decay * dir * exp(-decay*dist)
	where Kp is the point force constant, dist is the distance of the atom 
	from the center of the point force, decay is the energy decay with distance
	from the point force center, and dir is the normalized direction vector
	pointing from the point force center to the atom, indicating the direction
	of force.
**/
double		md_point_force(Bmolgroup* molgroup, Vector3<double> point, double Kpoint, double decay)
{
	if ( Kpoint <= 0 ) return 0;
	
	double			dist, en1, energy = 0;
	Vector3<double>	dir;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				dir = atom->coord - point;
				dist = dir.length();
				en1 = Kpoint*exp(-decay*dist);
				energy += en1;
				atom->F += dir * (decay*en1/dist);
			}
		}
	}
	
	return energy;
}
