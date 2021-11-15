/**
@file	model_mechanics.cpp
@brief	Functions to do molecular mechanics
@author Bernard Heymann
@date	Created: 20010828
@date	Modified: 20210318
**/

#include "model_mechanics.h"
#include "model_poly.h"
#include "model_shell.h"
#include "model_transform.h"
#include "model_util.h"
#include "math_util.h"
#include "random_numbers.h"
#include "symmetry.h"
#include "Vector3.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#define ETERMS	16

/**
@brief 	Minimizes the energy of a model.
@param 	*model			model structure.
@param 	&md				model parameters.
@param 	mm_type			type of mechanics: 0=minimization, 1=dynamics
@param 	max_iter		number of minimization iterations.
@param 	max_shift		maximum shift per iteration.
@param 	velocitylimit	limit on velocity per time step.
@return double			final energy.

	The inclusion of energy terms is based on positive K-constants in the 
	model parameter structure.
	Only the first model in the linked list is used.
**/
double		model_mechanics(Bmodel* model, Bmodparam& md, int mm_type, int max_iter,
				double max_shift, double velocitylimit)
{
	if ( !model->select() ) {
		cerr << "Error: No model selected!" << endl << endl;
		return 0;
	}

	if ( !model->poly ) model_poly_generate(model);

	if ( max_shift < 1 ) max_shift = 1e30;
	
	Vector3<int>	ksize;
	Bimage*			map = NULL;
	
	if ( md.Kmap && model->mapfile().length() > 1 ) {
		map = read_img(model->mapfile(), 1, model->image_number());
		map->rescale_to_avg_std(0, 1);
		if ( md.sigma ) {
			if ( md.sigma < map->sampling(0)[0] ) md.sigma = map->sampling(0)[0];
			ksize = Vector3<int>((int) (6*md.sigma/map->sampling(0)[0]), (int) (6*md.sigma/map->sampling(0)[1]), (int) (6*md.sigma/map->sampling(0)[2]));
		}
	}

	if ( verbose ) {
		if ( mm_type ) {
			cout << "Dynamics for model " << model->identifier() << ":" << endl;
			cout << "Time step:                      " << md.timestep << endl;
			cout << "Kfriction:                      " << md.Kfriction << endl;
		} else {
			cout << "Minimizing model " << model->identifier() << ":" << endl;
		}
		cout << "Kdistance:                      " << md.Kdistance << endl;
		cout << "Klink:                          " << md.Klink << endl;
		cout << "Kangle:                         " << md.Kangle << endl;
		cout << "Kpolyangle:                     " << md.Kpolyangle << endl;
		cout << "Kpolygon:                       " << md.Kpolygon << endl;
		cout << "Kpolyplane:                     " << md.Kpolyplane << endl;
		cout << "Kpoint and decay:               " << md.Kpoint << " (" << md.pointdecay << ")" << endl;
		cout << "Kradial and radius:             " << md.Kradial << " (" << md.radius << ")" << endl;
		cout << "Kplane:                         " << md.Kplane << endl;
		if ( md.guide ) {
			cout << "Guide model:                    " << md.guide->identifier() << endl;
			cout << "Kguide:                         " << md.Kguide << endl;
		}
		if ( map ) {
			cout << "Map:                            " << model->mapfile() << endl;
			cout << "Kmap:                           " << md.Kmap << endl;
			cout << "Sigma:                          " << md.sigma << " A" << endl;
			cout << "Kernel size:                    " << ksize << endl;	
		}
		cout << "Maximum shift per iteration:    " << max_shift << " A" << endl;
		cout << "Model mass:                     " << model_mass(model) << " Da" << endl << endl;
	}
	

	long			i, iter(0);
	int				h[ETERMS];
	double			E[ETERMS], E0[ETERMS], Epot0(0);
	for ( i=0; i<ETERMS; i++ ) {
		E[i] = E0[i] = 0;
		h[i] = 0;
	}
	if ( md.Klink > 0 ) h[0] = 1;
	if ( md.Kangle > 0 ) h[1] = 1;
	if ( md.Kpolyangle > 0 ) h[2] = 1;
	if ( md.Kpolygon > 0 ) h[3] = 1;
	if ( md.Kpolyplane > 0 ) h[4] = 1;
	if ( md.Kdistance > 0 ) h[5] = 1;
	if ( md.Kelec > 0 ) h[6] = 1;
	if ( md.Kpoint > 0 ) h[7] = 1;
	if ( md.Kradial > 0 ) h[8] = 1;
	if ( md.Kplane > 0 ) h[9] = 1;
	if ( md.Kguide > 0 ) h[10] = 1;
	if ( md.Kmap ) h[11] = 1;
	
	Bstring*		head_list = NULL;
	Bstring*		head;
	string_add(&head_list, "Elink     ");
	string_add(&head_list, "Eangle    ");
	string_add(&head_list, "Epolyangle");
	string_add(&head_list, "Epolygon  ");
	string_add(&head_list, "Epolyplane");
	string_add(&head_list, "Edistance ");
	string_add(&head_list, "Eelectrostatic ");
	string_add(&head_list, "Epoint    ");
	string_add(&head_list, "Eradial   ");
	string_add(&head_list, "Eplane    ");
	string_add(&head_list, "Eguide    ");
	string_add(&head_list, "Emap      ");
	

	model_calculate_deviations(model, md);
	
//	Vector3<double>	refcom = model_center_of_mass(model);
	
	double			gyrad0 = model_gyration_radius(model);
	
//	if ( verbose & VERB_PROCESS ) {
	if ( verbose ) {
		cout << "Iter";
		for ( i=0, head = head_list; i<ETERMS && head; i++, head = head->next ) if ( h[i] )
			cout << tab << *head;
		cout << "\tEpotential";
		if ( mm_type ) cout << "\tEkinetic";
		cout << endl;
	}
	for ( iter=1; iter<=max_iter; iter++ ) {
		model_zero_forces(model);
		E[0] = md.Elink = model_link_energy(model, md);
		E[1] = md.Eangle = model_angle_energy(model, md);
		E[2] = md.Epolyangle = model_polygon_angle_energy(model, md.Kpolyangle);
		E[3] = md.Epolygon = model_polygon_energy(model, md.Kpolygon);
		E[4] = md.Epolyplane = model_polygon_plane_energy(model, md.Kpolyplane);
		E[5] = md.Edistance = model_distance_energy(model, md);
//		E[5] = md.Edistance = model_grid_distance_energy(model, md);
//		E[5] = md.Edistance = model_neighbor_distance_energy(model, md);
		E[6] = md.Eelec = model_electrostatic_energy(model, md);
		E[7] = md.Epoint = model_point_force(model, md.point, md.Kpoint, md.pointdecay);
		E[8] = md.Eradial = model_radial_energy(model, md.point, md.radius, md.Kradial);
		E[9] = md.Eplane = model_neighbor_plane_energy(model, md.Kplane);
//		E[10] = md.Eguide = model_polyhedron_guide_energy(model, md.guide, md.Kguide);
		E[10] = md.Eguide = model_guide_energy(model, md);
		E[11] = md.Emap = model_map_energy(model, map, md.Kmap, md.sigma);
		for ( i=0, md.Epot=0; i<ETERMS; i++ ) md.Epot += E[i];
		if ( !isfinite(md.Epot) ) {
			cerr << endl << "***** Warning: The system exploded!!! *****" << endl << endl;
			return md.Epot;
		}
		if ( iter == 1 ) {
			for ( i=0; i<ETERMS; i++ ) E0[i] = E[i];
			Epot0 = md.Epot;
		}
		if ( mm_type < 1 ) model_minimize(model, max_shift);
		else md.Ekin = model_verlet(model, md.timestep, md.Kfriction, velocitylimit);
//		model_shift(model, refcom - model_center_of_mass(model));
//		if ( verbose & VERB_PROCESS ) {
		if ( verbose ) {
			cout << iter;
			for ( i=0; i<ETERMS; i++ ) if ( h[i] )
				cout << tab << E[i];
			cout << tab << md.Epot;
			if ( mm_type ) cout << tab << md.Ekin;
			if ( verbose < VERB_PROCESS ) {
				cout << "\r" << flush;
			} else {
				cout << endl;
			}
		}
	}

	double			gyrad = model_gyration_radius(model);
	
	if ( verbose ) {
		cout << endl << endl;
		cout << "EnergyTerm\tStart    \tEnd      \tChange" << endl;
		for ( i=0, head = head_list; i<ETERMS && head; i++, head = head->next ) if ( h[i] ) {
			cout << *head << tab << E0[i] << tab << E[i] << tab << E[i] - E0[i] << endl;
		}
		cout << "Epotential\t" << Epot0 << tab << md.Epot << tab << md.Epot - Epot0 << endl;
		if ( mm_type ) cout << "Ekinetic\t\t" << md.Ekin << endl;
		cout << "Gyration radius\t" << gyrad0 << tab << gyrad << tab << gyrad - gyrad0 << endl;
		cout << endl;
	}

	delete map;
	
	model_calculate_deviations(model, md);
		
	return md.Epot+md.Ekin;
}

/**
@brief 	Move components random distances down the energy gradient.
@param 	*model		model structure.
@param 	max_shift	maximum shift for each component.
@return int			number of components.

	The distance of movement is limited to the maximum shift.
	Only the first model in the linked list is used.

**/
int			model_minimize(Bmodel* model, double max_shift)
{
	int					n;
//	Bmodel*				mp;
	Bcomponent*			comp;
	
	double				irm = 1.0/get_rand_max();
	Vector3<double>		shift;
	
	for ( n=0, comp = model->comp; comp; comp = comp->next, n++ ) if ( comp->select() ) {
		if ( comp->type()->mass() > 0 ) comp->force(comp->force() / comp->type()->mass());
		shift = comp->force() * (random()*irm);
		shift = shift.min(max_shift);
		shift = shift.max(-max_shift);
		comp->shift(shift);
		comp->FOM(comp->force().length());
	}
	
	return n;
}

/**
@brief 	Model dynamics using the velocity verlet integrator.
@param 	*model			model structure.
@param 	timestep		dynamics time step.
@param 	Kfriction		friction coefficient.
@param 	velocitylimit	limit on velocity per time step.
@return double			kinetic energy.

	Leapfrog integration for any coordinate x, velocity vx and force Fx:
		x(t+1) = x(t) + vx(t+1) * dt
		vx(t+1) = (Fx(t) * dt/m + vx(t)) * kf
		where
			kf:	friction constant (1=no friction)
			dt:	time step
			m: atomic mass
	The velocity is limited each time step to damp chaotic oscillations.
	Only the first model in the linked list is used.
**/
double		model_verlet(Bmodel* model, double timestep, double Kfriction, double velocitylimit)
{
	double				vel, Ekin(0);
	double				minvel = -velocitylimit, maxvel = velocitylimit;	// Maximum velocity per time step
	Bcomponent*			comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		comp->velocity(comp->velocity() * ( 1 + (timestep/comp->type()->mass())));
		comp->velocity(comp->velocity() * Kfriction);
		comp->velocity(vector3_scalar_range(comp->velocity(), minvel, maxvel));
		comp->shift(comp->velocity() * timestep);
//		if ( wrap ) comp->location(vector3_set_PBC(comp->location(), box);
		vel = comp->velocity().length();
		Ekin += 0.5*comp->type()->mass()*vel*vel;
		comp->FOM(vel);
	}
	
	return Ekin;
}

/*
@brief 	Calculates the harmonic potential between components.
@param 	*link		link.
@param 	d0			reference distance between components.
@param 	Kd			distance force constant.
@return double		link energy.

	The harmonic energy function:
		E = Kd*(d - d0)^2
	The force is zero when d = reference link length (d0).
**/
double		component_harmonic_potential(Bcomponent* comp1, Bcomponent* comp2, double d0, double Kd)
{
	if ( d0 <= 0 || Kd <= 0 ) return 0;
	
	double			d, dev, fac, E(0);
	Vector3<double>	v, F;

	v = comp1->location() - comp2->location();
	d = v.length();
	dev = d - d0;
	E = Kd*dev*dev;
	if ( d ) {
		fac = -2*Kd*dev/d;
		F = v * fac;
		comp1->force(comp1->force() + F);
		comp2->force(comp2->force() - F);
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG component_harmonic_potential: c1=" << comp1->identifier() << 
			" c2=" << comp2->identifier() << " Kd=" << Kd << " d0=" << d0 << " E=" << E << endl;

	return E;
}

/*
@brief	Calculates the soft sphere potential between components.
@param	comp1		first component.
@param	comp2		second component.
@param 	d0			reference distance between components.
@param 	Kd			distance force constant.
@return	double		distance energy.

	The soft sphere potential is given by:
		E = Kd*(d0/d)^12
	The potential is set to zero for d > 3*d0;
**/
double		component_soft_potential(Bcomponent* comp1, Bcomponent* comp2, double d0, double Kd)
{
	if ( d0 <= 0 || Kd <= 0 ) return 0;
	
	double			d2, d02(d0*d0), fac, E(0);
	double			rd2, rd6, rd12;
	Vector3<double>	v, F;

	v = comp2->location() - comp1->location();	// Points from 1 to 2
	d2 = v.length2();
	if ( d2 < 9*d02 ) {
		rd2 = d02/d2;
		rd6 = rd2*rd2*rd2;
		rd12 = rd6*rd6;
		E = Kd*rd12;
		fac = -12*Kd*rd12/d2;
		F = v * fac;				// Force on 1 in direction of 2
		comp1->force(comp1->force() + F);
		comp2->force(comp2->force() - F);
	}

	return E;
}

/*
@brief	Calculates the Lennard-Jones potential between components.
@param	comp1		first component.
@param	comp2		second component.
@param 	d0			reference distance between components.
@param 	Kd			distance force constant.
@return	double		distance energy.

	The Lennard-Jones potential is given by:
		E = Kd*((d0/d)^12 - 2*(d0/d)^6)
	The force is zero when d = d0.
	The potential is set to zero for d > 3*d0;
**/
double		component_lennard_jones_potential(Bcomponent* comp1, Bcomponent* comp2, double d0, double Kd)
{
	if ( d0 <= 0 || Kd <= 0 ) return 0;
	
	double			d2, d02(d0*d0), fac, E(0);
	double			rd2, rd6, rd12;
	Vector3<double>	v, F;

	v = comp2->location() - comp1->location();	// Points from 1 to 2
	d2 = v.length2();
	if ( d2 < 9*d02 ) {
		rd2 = d02/d2;
		rd6 = rd2*rd2*rd2;
		rd12 = rd6*rd6;
		E = Kd*(rd12 - 2*rd6);
		fac = -12*Kd*(rd12 - rd6)/d2;
		F = v * fac;				// Force on 1 in direction of 2
		comp1->force(comp1->force() + F);
		comp2->force(comp2->force() - F);
	}
	
	return E;
}

/*
@brief	Calculates the Morse potential between components.
@param	comp1		first component.
@param	comp2		second component.
@param 	d0			reference distance between components.
@param 	Kd			distance force constant.
@return	double		distance energy.

	The Morse potential is given by:
		E = Kd*((1-exp((a/d0)*(d0-d)))^2 - 1)
	The width of the energy well is given by a, set here to 6.
	The force is zero when d = d0.
	The potential is set to zero for d > 3*d0;
**/
double		component_morse_potential(Bcomponent* comp1, Bcomponent* comp2, double d0, double Kd)
{
	if ( d0 <= 0 || Kd <= 0 ) return 0;
	
	double			d, fac, E(0);
	double			ef, ef1, a = 6/d0;
	Vector3<double>	v, F;

	v = comp2->location() - comp1->location();	// Points from 1 to 2
	d = v.length();
	if ( d < 3*d0 ) {
		ef = exp(a*(d0 - d));
		ef1 = 1 - ef;
		E = Kd*(ef1*ef1 - 1);
		fac = 2*a*Kd*ef1*ef/d;
		F = v * fac;	// Force on 1 in direction of 2
		comp1->force(comp1->force() + F);
		comp2->force(comp2->force() - F);
	}

	return E;
}

/*
@brief	Calculates the potential between two links to a component.
@param	comp		component.
@param	comp1		first component linked to.
@param	comp2		second component linked to.
@param 	a0			reference angle between linked components.
@param 	Ka			angular force constant.
@return	double		angular energy.

	The angular potential is given by:
		E = Ka*(a - a0)^2
	The force is zero when a = a0.
**/
double		component_angular_potential(Bcomponent* comp, Bcomponent* comp1, Bcomponent* comp2, double a0, double Ka)
{
	if ( Ka <= 0 ) return 0;
	
	double			v1len, v2len, da, fac, E(0);
	Vector3<double>	v1, v2, n;
	
	v1 = comp1->location() - comp->location();
	v2 = comp2->location() - comp->location();
	n = v1.cross(v2);
	n.normalize();
	v1len = v1.normalize();
	v2len = v2.normalize();
	da = acos(v1.scalar(v2)) - a0;
	fac = -2*Ka*da;
	E = Ka*da*da*v1len*v2len;
	comp1->force(comp1->force() - n.cross(v1) * fac * v1len);
	comp2->force(comp2->force() + n.cross(v2) * fac * v2len);

	return E;
}

/*
@brief	Calculates the potential between a component and a point.
@param	comp		component.
@param 	point		reference point = radial center.
@param 	d0			reference distance between component and point.
@param 	Kradial		radial force constant.
@return	double		radial energy.

	The radial potential is given by:
		E = Kd*(d - d0)^2
	The force is zero when d = d0.
**/
double		component_radial_potential(Bcomponent* comp, Vector3<double> point, double d0, double Kradial)
{
	if ( Kradial <= 0 ) return 0;
	
	double			d, dev, fac, E(0);
	Vector3<double>	v, F;

	v = comp->location() - point;	// Points from center to component
	d = v.length();
	dev = d - d0;
	E = Kradial*dev*dev;
	if ( d ) {
		fac = -2*Kradial*dev/d;
		F = v * fac;
		comp->force(comp->force() + F);
	}

	return E;
}

double		component_distance_potential(Bcomponent* comp1, Bcomponent* comp2, double Kd, int type)
{
	double			rd(1), E(0);

	rd = comp1->radius() + comp2->radius();

	switch ( type ) {
		case 1:
			E += component_harmonic_potential(comp1, comp2, rd, Kd);
			break;
		case 2:
			E += component_soft_potential(comp1, comp2, rd, Kd);
			break;
		case 3:
			E += component_lennard_jones_potential(comp1, comp2, rd, Kd);
			break;
		case 4:
			E += component_morse_potential(comp1, comp2, rd, Kd);
			break;
//		default:
	}
	
	return E;
}

double		component_electrostatic_potential(Bcomponent* comp1, Bcomponent* comp2,
				double q1, double q2, double Ke, double cutoff)
{
	if ( Ke <= 0 ) return 0;

	double			fac(Ke*q1*q2);
	
	if ( fac == 0 ) return 0;

	double			d, d2, E(0);
	Vector3<double>	v, F;

	v = comp1->location() - comp2->location();
	d2 = v.length2();
	d = sqrt(d2);
	if ( d < cutoff ) {
		fac /= d;
		E = fac;
		fac /= d2;
		F = v * fac;
		comp1->force(comp1->force() + F);
		comp2->force(comp2->force() - F);
	}

	return E;
}

/**
@brief 	Calculates the electrostatic potentials between components.
@param 	*model		model structure.
@param 	&md			model parameters with distance interactions specifications.
@return double		electrostatic energy.

	Electrostaic potential:
		E = Ke*q1*q2/d
	Only the first model in the linked list is used.
**/
double		model_electrostatic_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kelec <= 0 ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_electrostatic_energy: Ke=" << md.Kelec << endl;
	
	double			E(0), q1, q2;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp->next; comp = comp->next ) if ( comp->select() && comp->type()->index() >= 0 ) {
		for ( comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() && comp2->type()->index() >= 0 ) {
			q1 = comp->type()->charge();
			q2 = comp2->type()->charge();
			E += component_electrostatic_potential(comp, comp2, q1, q2, md.Kelec, md.cutoff);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_electrostatic_energy: c1=" << comp->identifier() <<
					" c2=" << comp2->identifier() << endl;
		}
	}

	return E;
}

/**
@brief 	Calculates the distance-related potentials between components.
@param 	*model		model structure.
@param 	&md			model parameters with distance interactions specifications.
@return double		distance energy.

	Distance potential types:
		0	none
		1	harmonic - only for explicit links
		2	soft
		3	Lennard-Jones
		4	Morse
	Only the first model in the linked list is used.
**/
double		model_distance_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kdistance <= 0 ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_distance_energy: Kd=" << md.Kdistance << endl;
	
	long			i, j;
	double			rd, Kd, E(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp->next; comp = comp->next ) if ( comp->select() && comp->type()->index() >= 0 ) {
		for ( comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() && comp2->type()->index() >= 0 ) {
//			i = comp->type()->index()*md.ntype + comp2->type()->index();
//			Kd = md.Kdistance * md.Kd[i];
			i = comp->type()->index();
			j = comp2->type()->index();
			Blinktype&		lt = md.linktype[i][j];
			rd = lt.distance();
			Kd = md.Kdistance * lt.Kdistance();
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_distance_energy: c1=" << comp->identifier() << 
					" c2=" << comp2->identifier() << " type=" << lt.select() <<
					" Kd=" << Kd << " d0=" << rd << endl;
			switch ( lt.select() ) {
				case 1:
					E += component_harmonic_potential(comp, comp2, rd, Kd);
					break;
				case 2:
					E += component_soft_potential(comp, comp2, rd, Kd);
					break;
				case 3:
					E += component_lennard_jones_potential(comp, comp2, rd, Kd);
					break;
				case 4:
					E += component_morse_potential(comp, comp2, rd, Kd);
					break;
//				default:
			}
		}
	}

	return E;
}

/**
@brief 	Calculates the non-bonded forces and energy.
@param 	*model		model structure.
@param 	&md			model parameters with distance interactions specifications.
@return double		distance energy.


	Distance potential types:
		0	none
		1	harmonic - only for explicit links
		2	soft
		3	Lennard-Jones
		4	Morse

**/
double		model_grid_distance_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kdistance <= 0 ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_distance_energy: Kd=" << md.Kdistance << endl;

	Vector3<long> 	gsize;
	Vector3<double>	gori, gsam(md.cutoff, md.cutoff, md.cutoff);
	vector<vector<Bcomponent*>>	grid = model_component_grid(model, gsize, gori, gsam);

	int				dodist;
	long			i, ii, x, y, z, xx, yy, zz, ix, iy, iz;
	
	md.Edistance = md.Eelec = 0;
	
	for ( i=z=0; z<gsize[2]; ++z ) {
		for ( y=0; y<gsize[1]; ++y ) {
			for ( x=0; x<gsize[0]; ++x, ++i ) {
				for ( Bcomponent* comp1: grid[i] ) {
					for ( zz=z-1; zz<=z+1; ++zz ) {
						iz = zz;
						if ( iz < 0 ) iz += gsize[2];
						if ( iz >= gsize[2]) iz -= gsize[2];
						for ( yy=y-1; yy<=y+1; ++yy ) {
							iy = yy;
							if ( iy < 0 ) iy += gsize[1];
							if ( iy >= gsize[1] ) iy -= gsize[1];
							for ( xx=x-1; xx<=x+1; ++xx ) {
								ix = xx;
								if ( ix < 0 ) ix += gsize[0];
								if ( ix >= gsize[0] ) ix -= gsize[0];
								ii = (iz*gsize[1] + iy)*gsize[0] + ix;
								for ( Bcomponent* comp2: grid[ii] ) {
									dodist = 1;
									if ( comp1 == comp2 ) dodist = 0;	// Same atom
//									else if ( latom2->atom->next == latom->atom ) dononbond = 0;	// Bonded atoms
//									else if ( latom2->atom == latom->atom->next ) dononbond = 0;	// Bonded atoms
//									else if ( latom2->atom->next ) {	// Atom 2 bonds away
//										if ( latom2->atom->next->next == latom->atom ) dononbond = 0;
//									} else if ( latom->atom->next ) {	// Atom 2 bonds away
//										if ( latom2->atom == latom->atom->next->next ) dononbond = 0;
//									}
									if ( dodist ) {
										md.Edistance += component_distance_potential(comp1, comp2, md.Kdistance, 4);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	return md.Eelec+md.Edistance;
}

/**
@brief 	Calculates the distance-related potentials between components and their neighbors.
@param 	*model		model structure.
@param 	&md			model parameters with distance interactions specifications.
@return double		distance energy.

	Distance potential types:
		0	none
		1	harmonic - only for explicit links
		2	soft
		3	Lennard-Jones
		4	Morse
	Only the first model in the linked list is used.
**/
double		model_neighbor_distance_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kdistance <= 0 ) return 0;
	
	if ( md.linktype.size() < 1 ) {
		cerr << "The link type distance matrix is not defined!" << endl;
		bexit(-1);
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG model_neighbor_distance_energy: Kd=" << md.Kdistance << endl;
		cout << "DEBUG model_neighbor_distance_energy: lt.size=" << md.linktype.size() << endl;
	}
	
	long			i, j, n;
	double			rd, Kd, E(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		i = comp->type()->index();
//		if ( i1 >= 0 ) {
			for ( n=0; n<comp->neighbor.size() && comp->neighbor[n]; n++ ) {
				comp2 = comp->neighbor[n];
				j = comp2->type()->index();
				if ( ( comp->identifier() < comp2->identifier() ) && comp2->select() ) {
//					i = i1*md.ntype + i2;
					Blinktype&		lt = md.linktype[i][j];
					rd = lt.distance();
					Kd = md.Kdistance * lt.Kdistance();
					switch ( lt.select() ) {
						case 1:
							E += component_harmonic_potential(comp, comp2, rd, Kd);
							break;
						case 2:
							E += component_soft_potential(comp, comp2, rd, Kd);
							break;
						case 3:
							E += component_lennard_jones_potential(comp, comp2, rd, Kd);
							break;
						case 4:
							E += component_morse_potential(comp, comp2, rd, Kd);
							break;
//						default:
					}
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG model_neighbor_distance_energy: c1=" << comp->identifier() <<
							" c2=" << comp2->identifier() << " type=" << lt.select() <<
							" Kd=" << Kd << " d0=" << rd << " E=" << E << endl;
				}
			}
//		}
	}

	return E;
}

/**
@brief 	Calculates the soft sphere potential between components.
@param 	*model		model structure.
@param 	Kd			distance force constant.
@param 	d0			reference distance between components.
@return double		distance energy.

	The soft sphere potential is given by:
		E = Kd*(d0/d)^12
	The potential is set to zero for d > 3*d0;
	Only the first model in the linked list is used.
**/
double		model_soft_sphere_energy(Bmodel* model, double Kd, double d0)
{
	if ( Kd <= 0 ) return 0;
	
	double			E(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp->next; comp = comp->next )
		for ( comp2 = comp->next; comp2; comp2 = comp2->next )
			E += component_soft_potential(comp, comp2, d0, Kd);

	return E;
}

/**
@brief 	Calculates the Lennard-Jones potential between components.
@param 	*model		model structure.
@param 	Kd			distance force constant (Kd).
@param 	d0			reference distance between components (d0).
@return double		distance energy.

	The Lennard-Jones potential is given by:
		E = Kd*((d0/d)^12 - 2*(d0/d)^6)
	The force is zero when d = d0.
	The potential is set to zero for d > 3*d0;
	Only the first model in the linked list is used.
**/
double		model_lennard_jones_energy(Bmodel* model, double Kd, double d0)
{
	if ( Kd <= 0 ) return 0;
	
	double			E(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp->next; comp = comp->next )
		for ( comp2 = comp->next; comp2; comp2 = comp2->next )
			E += component_lennard_jones_potential(comp, comp2, d0, Kd);

	return E;
}

/**
@brief 	Calculates the Morse potential between components.
@param 	*model		model structure.
@param 	Kd			distance force constant (Kd).
@param 	d0			reference distance between components (d0).
@return double		distance energy.

	The Morse potential is given by:
		E = Kd*((1-exp((a/d0)*(d0-d)))^2 - 1)
	The width of the energy well is given by a, set here to 6.
	The force is zero when d = d0.
	The potential is set to zero for d > 3*d0;
	Only the first model in the linked list is used.
**/
double		model_morse_energy(Bmodel* model, double Kd, double d0)
{
	if ( Kd <= 0 ) return 0;
	
	double			E(0);
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp->next; comp = comp->next )
		for ( comp2 = comp->next; comp2; comp2 = comp2->next )
			E += component_morse_potential(comp, comp2, d0, Kd);

	return E;
}

/**
@brief 	Calculates the model link energy.
@param 	*model		model structure.
@param 	&md			model parameters with distance interactions specifications.
@return double		distance energy.

	The link potential is a harmonic function:
		E = Kl*(d - l0)^2
	The force is zero when d = reference link length (l0).
	Only the first model in the linked list is used.
**/
double		model_link_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Klink <= 0 ) return 0;
	
	int				i, j;
	double			E(0);
	Blink*			link;
	
	for ( link = model->link; link; link = link->next ) {
//		i = link->comp[0]->type()->index()*md.ntype + link->comp[1]->type()->index();
//		E += component_harmonic_potential(link->comp[0], link->comp[1],
//				md.llen[i], md.Klink * md.Kl[i]);
		i = link->comp[0]->type()->index();
		j = link->comp[1]->type()->index();
		Blinktype&	lt = md.linktype[i][j];
		E += component_harmonic_potential(link->comp[0], link->comp[1],
				lt.length(), md.Klink * lt.Klength());
	}

	return E;
}

/**
@brief 	Calculates the model link energy.
@param 	*model		model structure.
@param 	Klink		link force constant (Kl).
@return double		link energy.

	The link potential is a harmonic function:
		E = Kl*(d - l0)^2
	The force is zero when d = reference link length (l0).
	Only the first model in the linked list is used.
**/
double		model_link_energy(Bmodel* model, double Klink)
{
	if ( Klink <= 0 ) return 0;
	
	double			E(0);
	Blink*			link;
	
	for ( link = model->link; link; link = link->next )
		E += component_harmonic_potential(link->comp[0], link->comp[1], link->length(), Klink);

	return E;
}

/**
@brief 	Calculates the model angular energy.
@param 	*model		model structure.
@param 	&md			model parameters with reference angle specifications.
@return double		angular energy.

	The angular potential is a harmonic function:
		E = Ka*(a - a0)^2
	The force is zero when a = reference angle (a0).
	Only the first model in the linked list is used.
**/
double		model_angle_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kangle <= 0 ) return 0;
	
	int				i, j, k;
	double			E(0);

	Bcomponent*		comp;
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		i = comp->type()->index();
		for ( j=1; j<comp->link.size() && comp->link[j]; j++ ) {
			comp1 = comp->link[j];
			j = comp1->type()->index();
			for ( k=0; k<j && comp->link[k]; k++ ) {
				comp2 = comp->link[k];
				k = comp2->type()->index();
				Bangletype&		at = md.angletype[i][j][k];
//				i = (comp->type()->index()*md.ntype + comp1->type()->index())*md.ntype + comp2->type()->index();
//				E += component_angular_potential(comp, comp1, comp2, md.angle[i], md.Kangle * md.Ka[i]);
				E += component_angular_potential(comp, comp1, comp2, at.angle(), md.Kangle * at.Kangle());
			}
		}
	}
	
	return E;
}


/**
@brief 	Calculates the link angular forces and energy.
@param 	*model			model structure.
@param 	Kpolyangle		angle force constant (Ka).
@return double			total link angle energy.

	The energy is defined as a harmonic function around the reference 
	link angle, a0:
		Ea = Ka*(cos(a0)-r1*r2/(|r1|*|r2|))^2
		Ea = Ka*(a0 - a)^2
	The force is the derivative of the energy on the first and last atoms:
		Fa1 = 2*Ka*(cos(a0)-r1*r2/(|r1|*|r2|))/(|r1|*|r2|) * ((r1*r2/|r1|)*r1-r2)
		Fa3 = 2*Ka*(cos(a0)-r1*r2/(|r1|*|r2|))/(|r1|*|r2|) * ((r1*r2/|r2|)*r2-r1)
	where r1 is the vector from atom 2 to atom 1, r2 is the vector from
	atom 2 to atom 3, and Ka is the link angle force constant.
	Only the first model in the linked list is used.
**/
double		model_polygon_angle_energy(Bmodel* model, double Kpolyangle)
{
	if ( Kpolyangle <= 0 ) return 0;
	
	int				i;
	double			a0, E(0);
	Bpolygon*		poly;
	Bcomponent		*comp, *comp1, *comp2;

	if ( !model->poly ) model_poly_generate(model);
	
	for ( poly = model->poly; poly; poly = poly->next ) if ( poly->closed() ) {
		a0 = M_PI*(1 - 2.0L/poly->size());
		for ( i=0, comp1 = poly->comp[poly->size()-1]; i<poly->comp.size() && poly->comp[i]; comp1 = comp, i++ ) {
			comp = poly->comp[i];
			comp2 = poly->comp[i+1];
			if ( !comp2 ) comp2 = poly->comp[0];
			E += component_angular_potential(comp, comp1, comp2, a0, Kpolyangle);
		}
	}

	return E;
}

/**
@brief 	Calculates the polygon regularity.
@param 	*model		model structure.
@param 	linklength	reference link length (l0).
@param 	Kpoly		polygon regularity constant (Kp).
@return double		polygon energy.

	Given the distances of all the vertices from the polygon center, the
	polygon regularity is defined as the deviation of these distances
	from the reference distance:
		E = Kp*(d - d0)
	where the reference distance is calculated from the link length:
		d0 = l0/sqrt(2*(1-cos(2PI/n))
	where n is the polygon order.
	Only the first model in the linked list is used.
**/
double		model_polygon_energy(Bmodel* model, double linklength, double Kpoly)
{
	if ( Kpoly <= 0 ) return 0;
	
	int				i;
	double			dr, dist, dev, fac, E(0);
	Vector3<double>	d, polycen, dcen;
	Bpolygon*		poly;
	
	if ( !model->poly ) model_poly_generate(model);

	for ( poly = model->poly; poly; poly = poly->next ) if ( poly->closed() ) {
		for ( i=0, polycen = 0; i<poly->comp.size() && poly->comp[i]; i++ )
			polycen += poly->comp[i]->location();	// Center of polygon
		polycen /= i;
		if ( linklength > 0 ) {
			dr = linklength/sqrt(2*(1 - cos(TWOPI/i)));	// Reference distance to polygon center
		} else {
			for ( i=0, dcen = 0; i<poly->comp.size() && poly->comp[i]; i++ )
				dcen += poly->comp[i]->location().distance(polycen);	// Distance to center
			dcen /= i;
			dr = dcen.length();	// Reference distance to polygon center
		}
		for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
			d = poly->comp[i]->location() - polycen;
			dist = d.length();
			dev = dist - dr;
			E += Kpoly*dev*dev;
			fac = -2*Kpoly*dev/dist;
			poly->comp[i]->force(poly->comp[i]->force() + d*fac);
		}
	}
	
	return E;
}

/**
@brief 	Calculates the polygon regularity energy and forces.
@param 	*model		model structure.
@param 	Kpoly		polygon regularity constant (Kp).
@return double		polygon energy.

	Given the distances of all the vertices from the polygon center, the
	polygon regularity is defined as the deviation of these distances
	from the average distance:
		E = Kp*(d - d0)
	Only the first model in the linked list is used.
**/
double		model_polygon_energy(Bmodel* model, double Kpoly)
{
	return model_polygon_energy(model, 0, Kpoly);
}

/**
@brief 	Calculates the polygon planarity energy and forces.
@param 	*model		model structure.
@param 	Kpolyplane	polygon planarity constant.
@return double		polygon planar energy.

	A plane is fit through the polygon vertices and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	The polygon plane energy is calculated as a harmonic deviation from 
	the fitted plane.
	Only the first model in the linked list is used.
**/
double		model_polygon_plane_energy(Bmodel* model, double Kpolyplane)
{
	if ( Kpolyplane <= 0 ) return 0;
	
    long     	i;
	double				dev, fac, E(0), offset;
	Vector3<double>		normal;
	Bpolygon*			poly;

	if ( !model->poly ) model_poly_generate(model);

	for ( poly = model->poly; poly; poly = poly->next ) if ( poly->closed() ) {
		normal = component_plane(poly->comp, offset);
		for ( i=0; i<poly->comp.size() && poly->comp[i]; i++ ) {
			dev = poly->comp[i]->location().scalar(normal) - offset;
			E += Kpolyplane*dev*dev;
			fac = -2*Kpolyplane*dev;
			poly->comp[i]->force(poly->comp[i]->force() + normal*fac);
		}
	}
	
	return E;
}

/**
@brief 	Calculates the planarity energy and forces with respect to neighbors.
@param 	*model		model structure.
@param 	Kplane		neighbor planarity constant.
@return double		neighbor planar energy.

	A plane is fit through the neigbor locations and the normal calculated from:
		n•p = d
	where n is the normal vector, p is a point in the plane, and d is the offset.
	The deviation of a component location from the neighbor plane is calculated
	and converted to a harmonic energy and force.
	Only the first model in the linked list is used.
**/
double		model_neighbor_plane_energy(Bmodel* model, double Kplane)
{
	if ( Kplane <= 0 ) return 0;
	
	double				dev, fac, E(0), offset;
	Vector3<double>		v, normal, center;
	Bcomponent*			comp;

	for ( comp = model->comp; comp; comp = comp->next ) {
		normal = component_plane(comp->neighbor, offset);
		dev = comp->location().scalar(normal) - offset;
		E += Kplane*dev*dev;
		fac = -2*Kplane*dev;
		comp->force(comp->force() + normal*fac);
/*		for ( i=0; i<MAXLINK && comp->neighbor[i]; i++ ) {
			v = comp->neighbor[i]->location() - comp->location();
			v -= normal * normal.scalar(v);
			dev = v.normalize() - 50;
			E += Kplane*dev*dev;
			fac = -2*Kplane*dev;
			comp->force(comp->force() + v*fac);
		}*/
	}
	
	return E;
}

/**
@brief 	Calculates the forces and energy resulting from a single point force.
@param 	*model		model structure.
@param 	point		center of point force.
@param 	Kpoint		point force constant.
@param 	decay		energy decay with distance.
@return double		point force energy.

	The energy is defined as an exponential decay over distance from the 
	center of the point force:
		Ep = Kp * exp(-decay*dist)
	The force is the derivative of the energy:
		Fp = Kp * decay * dir * exp(-decay*dist)
	where Kp is the point force constant, dist is the distance of the component 
	from the center of the point force, decay is the energy decay with distance
	from the point force center, and dir is the normalized direction vector
	pointing from the point force center to the component, indicating the direction
	of force.
	Only the first model in the linked list is used.
**/
double		model_point_force(Bmodel* model, Vector3<double> point, double Kpoint, double decay)
{
	if ( Kpoint <= 0 ) return 0;
	
	double			dist, en1, energy(0);
	Vector3<double>	dir;
	Bmodel*			mp;
	Bcomponent*		comp;

	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			dir = comp->location() - point;
			dist = dir.length();
			en1 = Kpoint*exp(-decay*dist);
			energy += en1;
			comp->force(comp->force() + dir * (decay*en1/dist));
		}
	}
	
	return energy;
}

/**
@brief 	Calculates the potential between a components and a point.
@param 	*model		model structure.
@param 	point		reference point = radial center.
@param 	d0			reference distance between component and point.
@param 	Kradial		radial force constant.
@return double		energy.

	The radial potential is given by:
		E = Kd*(d - d0)^2
	The force is zero when d = d0.
	Only the first model in the linked list is used.
**/
double		model_radial_energy(Bmodel* model, Vector3<double> point, double d0, double Kradial)
{
	if ( Kradial <= 0 ) return 0;
	
	double			E(0);
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next )
		E += component_radial_potential(comp, point, d0, Kradial);
	
	return E;
}

/**
@brief 	Calculates the potential with respect to a guide model.
@param 	*model			model structure.
@param 	&md				model aparameters.
@return double			energy.

	Only the first model in the linked list is used.
**/
double		model_guide_energy(Bmodel* model, Bmodparam& md)
{
	if ( md.Kguide <= 0 ) return 0;
	if ( !md.guide ) return 0;
		
	double			d, fac(2*md.Kguide), E(0);
	Vector3<double>	v, F;
	Bcomponent*		comp;
	Bcomponent*		gcomp;

	for ( comp = model->comp; comp; comp = comp->next ) {
		gcomp = md.guide->closest_component(comp->location());
		if ( gcomp ) {
			v = gcomp->location() - comp->location();
			d = v.length();
			if ( d > md.cutoff ) {
				d -= md.cutoff;
				E += md.Kguide*d*d;
				F = v * fac;
				comp->force(comp->force() + F);
			}
		}
	}

	return E;
}

/**
@brief 	Calculates the potential with respect to a guide polyhedron.
@param 	*model			model structure.
@param 	*gmod			guide polyhedron model.
@param 	Kguide			guide polyhedron force constant.
@return double			energy.

	Only the first model in the linked list is used.
**/
double		model_polyhedron_guide_energy(Bmodel* model, Bmodel* gmod, double Kguide)
{
	if ( Kguide <= 0 ) return 0;
	if ( !gmod ) return 0;
		
	int				n;
	double			d, len, fac, E(0);
	Vector3<double>	center, F;
	Bcomponent*		comp;
	Bcomponent*		gcomp;

	// Calculate the center of the vertices in the guide model
	for ( n=0, gcomp = gmod->comp; gcomp; gcomp = gcomp->next, n++ )
		center += gcomp->location();
	
	if ( n ) center /= n;
	
	// Calculate the component vectors in the guide model
	for ( gcomp = gmod->comp; gcomp; gcomp = gcomp->next )
		gcomp->force(gcomp->location());
//		gcomp->force(gcomp->location() - center);
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		d = model_inside_outside(comp->location(), gmod, 0, 1);
		if ( !isfinite(d) ) {
			cerr << "d not finite: " << d << endl;
			bexit(-1);
		}
		E += Kguide*d*d;
		len = comp->location().length();
		if ( len ) {
			fac = 2*Kguide*d/len;
			F = comp->location() * fac;
			comp->force(comp->force() + F);
		}
	}

	return E;
}

/**
@brief 	Energy calculation of a model into a map.  
@param 	*model		model structure.
@param 	*map		map.
@param 	Kmap		map force constant.
@return double		energy.

	The map must be possitive density.
	Only the first model in the linked list is used.
**/
double		model_map_energy(Bmodel* model, Bimage* map, double Kmap)
{
	if ( !model ) return 0;
	if ( !map ) return 0;
	if ( Kmap == 0 ) return 0;
	
	map->change_type(Float);
	
	long				i, n(0);
	long				nx(map->sizeX()), nxy(nx*map->sizeY());
	double				E(0), scale(1/map->standard_deviation());
	Vector3<long>		coor;
	Vector3<double>		invunit(1/map->sampling(0));
	Vector3<double>		g;
	Bcomponent*			comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		coor = map->image->image_coordinates(comp->location());
		i = map->index(coor, 0);
		if ( coor[0] <= 0 ) g[0] = -coor[0];
		else if ( coor[0] >= map->sizeX()-1 ) g[0] = map->sizeX() - coor[0];
		else g[0] = (*map)[i+1] - (*map)[i-1];
		if ( coor[1] <= 0 ) g[1] = -coor[1];
		else if ( coor[1] >= map->sizeY()-1 ) g[1] = map->sizeY() - coor[1];
		else g[1] = (*map)[i+nx] - (*map)[i-nx];
		if ( coor[2] <= 0 ) g[2] = -coor[2];
		else if ( coor[2] >= map->sizeZ()-1 ) g[2] = map->sizeZ() - coor[2];
		else g[2] = (*map)[i+nxy] - (*map)[i-nxy];
		n++;
		g *= invunit;
		E += Kmap*g.length2();
		g *= 2*Kmap;
		comp->force(comp->force() + g);
	}
	
	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Model outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap*((E/n + map->average())*scale + 10);
	
	return E;
}

double		model_map_gradient_energy(Bmodel* model, Bimage* map, double Kmap)
{
	if ( !model ) return 0;
	if ( !map ) return 0;
	if ( Kmap == 0 ) return 0;
	
	if ( map->channels() != 3 ) {
		cerr << "Error: A gradient map is required!" << endl;
		bexit(-1);
	}
	
	map->change_type(Float);
	
	long				i, n(0);
	double				E(0);
	Vector3<long>		loc;
	Vector3<double>		invunit(0.5*Kmap/map->sampling(0)[0], 0.5*Kmap/map->sampling(0)[1], 0.5*Kmap/map->sampling(0)[2]);
	Vector3<double>		F;
	Bcomponent*			comp;
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		loc = map->image->image_coordinates(comp->location());
		if ( map->within_boundaries(loc) ) {
			i = map->index(loc, 0);
			F = map->vector3(i);
			E += F.length2();
			comp->force(comp->force() + F * invunit);
			n++;
		}
	}
	
	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Model outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap*((E/n + map->average())/map->standard_deviation() + 10);
	
	return E;
}

double		model_map_energy(Bmodel* model, Bimage* map, double Kmap, double sigma)
{
	if ( !model ) return 0;
	if ( !map ) return 0;
	if ( Kmap == 0 ) return 0;
	
	if ( map->channels() == 3 ) return model_map_gradient_energy(model, map, Kmap);
	
	if ( sigma < map->sampling(0)[0] ) return model_map_energy(model, map, Kmap);
	
	map->change_type(Float);
	
	Vector3<double>		pixsigma(sigma/map->sampling(0)[0], sigma/map->sampling(0)[1], sigma/map->sampling(0)[2]);
//	Vector3<double>		pixsigma2 = pixsigma*pixsigma;
	Vector3<int>		ksize = Vector3<int>((int) (6*pixsigma[0]), (int) (6*pixsigma[1]), (int) (6*pixsigma[2]));
	ksize = ksize.max(3);
	Vector3<int>		hksize = ksize/2;
	double*				kernel = new double[(long)ksize.volume()];
	
	long				i, k, n(0);
	long				x, y, z, ix, iy, iz, kx, ky, kz;
	double				dx, dy, dz, v, w, E1, E(0);
	double				sigma2 = sigma*sigma;
	double				ffac = Kmap/sigma2;
	Vector3<double>		F;
	Bcomponent*			comp;
	
	for ( i=0, z=0; z<ksize[2]; z++ ) {
		dz = (z - hksize[2])/pixsigma[2];
		for ( y=0; y<ksize[1]; y++ ) {
			dy = (y - hksize[1])/pixsigma[1];
			for ( x=0; x<ksize[0]; x++, i++ ) {
				dx = (x - hksize[0])/pixsigma[0];
				kernel[i] = exp(-0.5*(dx*dx + dy*dy + dz*dz));
			}
		}
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		x = (long) (comp->location()[0]/map->sampling(0)[0] + map->image->origin()[0] + 0.5);
		y = (long) (comp->location()[1]/map->sampling(0)[1] + map->image->origin()[1] + 0.5);
		z = (long) (comp->location()[2]/map->sampling(0)[2] + map->image->origin()[2] + 0.5);
		E1 = w =0;
		F = 0;
		for ( iz = z - hksize[2], kz=0; kz<ksize[2]; iz++, kz++ ) if ( iz>=0 && iz<map->sizeZ() ) {
			dz = (kz - hksize[2])*map->sampling(0)[2];
			for ( iy = y - hksize[1], ky=0; ky<ksize[1]; iy++, ky++ ) if ( iy>=0 && iy<map->sizeY() ) {
				dy = (ky - hksize[1])*map->sampling(0)[1];
				for ( ix = x - hksize[0], kx=0; kx<ksize[0]; ix++, kx++ ) if ( ix>=0 && ix<map->sizeX() ) {
					dx = (kx - hksize[0])*map->sampling(0)[0];
					i = map->index(0, x, y, z, 0);
					k = (kz*ksize[1] + ky)*ksize[0] + kx;
					v = (*map)[i]*kernel[k];
					w += kernel[k];
					E1 += v;
					F += Vector3<double>(v*dx, v*dy, v*dz);
				}
			}
		}
		if ( w ) {
			E1 /= w;		// Average density, correction for kernel volume
			F /= w;			// Average force
			comp->density(E1);
			if ( comp->type()->mass() ) {
				E1 *= comp->type()->mass();
				F *= comp->type()->mass();
			}
			E -= E1;
			comp->force(comp->force() + F * ffac);
			n++;
		}
	}
	
	if ( n < 1 ) {
		if ( verbose & VERB_PROCESS )
			cout << "Warning: Model outside map boundaries!" << endl;
		E = 1e100;
	} else
		E = Kmap*E/n;
	
	delete[] kernel;
	
	return E;
}

/**
@brief 	Zeroes the component forces.
@param 	*model		model.
@return int			number of components.

	Only the first model in the linked list is used.
**/
int			model_zero_forces(Bmodel* model)
{
	int				n(0);
	Bcomponent*		comp;

	for ( comp = model->comp; comp; comp = comp->next, n++ )
		comp->force(Vector3<float>(0,0,0));
	
	return n;
}

/**
@brief 	Calculates the deviations from regularity in a polyhedral model.
@param 	*model		model structure.
@return int			0.

	Link deviation:
		d = sqrt(sum((l - l0)^2)/n)
	where l is the length and l0 is the reference length.
	Angle deviation:
		d = sqrt(sum((a - a0)^2)/n)
	where a is the angle between two links and a0 is the reference angle:
		a0 = PI*(1 - 2/v)
	where v is the number of vertices in the associated polygon.
	Only the first model in the linked list is used.
**/
int			model_calculate_deviations(Bmodel* model)
{
	int					i, n;
	double				ll, d, dd, a, aa, ar, ad;
	Vector3<double>		d1, d2;
	Bcomponent			*comp1, *comp2, *comp3;
	Blink*				link;
	Bpolygon*			poly;
	
	for ( n=0, ll=dd=0, link = model->link; link; link = link->next, n++ ) {
		d = link->comp[0]->location().distance(link->comp[1]->location());
		ll += d;
		d -= link->length();
		dd += d*d;
	}
	
	if ( n ) {
		ll /= n;
		dd = sqrt(dd/n);
	}
	
	if ( verbose )
		cout << "Link average and deviation:     " << ll << " " << dd << " A" << endl;
	
	for ( n=0, aa=ad=0, poly = model->poly; poly; poly = poly->next ) if ( poly->closed() ) {
		ar = M_PI*(1 - 2.0L/poly->size());
		for ( i=0, comp1 = poly->comp[poly->size()-1]; i<poly->comp.size() && poly->comp[i]; comp1 = comp2, i++ ) {
			comp2 = poly->comp[i];
			comp3 = poly->comp[i+1];
			if ( !comp3 ) comp3 = poly->comp[0];
			d1 = comp2->location() - comp1->location();
			d2 = comp2->location() - comp3->location();
			a = d1.angle(d2);
			aa += a;
			a -= ar;
			ad += a*a;
			n++;
		}
	}
	
	if ( n ) {
		aa /= n;
		ad = sqrt(ad/n);
	}
	
	if ( verbose )
		cout << "Angle average and deviation:    " << aa*180.0/M_PI << " " << ad*180.0/M_PI << " degrees" << endl << endl;

	return 0;
}

int			model_calculate_deviations(Bmodel* model, Bmodparam& md)
{
	int					i, j, k, n;
	int					ii, jj, kk;
	double				ll, d, dd, a, aa, ar, ad;
	Vector3<double>		d1, d2;
	Bcomponent			*comp, *comp1, *comp2;
	Blink*				link;
	Bpolygon*			poly;
	
	if ( verbose )
		cout << "        \tCount\tAverage\tDeviation" << endl;
	
	for ( n=0, ll=dd=0, link = model->link; link; link = link->next, n++ ) {
		i = link->comp[0]->type()->index();
		j = link->comp[1]->type()->index();
		d = link->comp[0]->location().distance(link->comp[1]->location());
		Blinktype&		lt = md.linktype[i][j];
		ll += d;
		d -= lt.length();
		dd += d*d;
	}
	
	if ( n ) {
		ll /= n;
		dd = sqrt(dd/n);
	}
	
	if ( verbose )
		cout << "Links    \t" << n << tab << ll << tab << dd << endl;

	for ( n=0, aa=ad=0, comp = model->comp; comp; comp = comp->next ) {
		ii = comp->type()->index();
		for ( j=1; j<comp->link.size(); j++ ) {
			comp1 = comp->link[j];
			jj = comp1->type()->index();
			for ( k=0; k<j; k++ ) {
				comp2 = comp->link[k];
				kk = comp2->type()->index();
				Bangletype&		at = md.angletype[ii][jj][kk];
//				i = (comp->type()->index()*md.ntype + comp1->type()->index())*md.ntype + comp2->type()->index();
				if ( at.angle() ) {
					d1 = comp1->location() - comp->location();
					d2 = comp2->location() - comp->location();
					a = d1.angle(d2);
					aa += a;
					a -= at.angle();
					ad += a*a;
					n++;
				}
			}
		}
	}
	
	if ( n ) {
		aa /= n;
		ad = sqrt(ad/n);
	}
	
	if ( verbose )
		cout << "Angles   \t" << n << tab << aa*180.0/M_PI << tab << ad*180.0/M_PI << endl;

//	cout << "Polygons: " << model->polygon_count() << endl;
	
	for ( n=0, aa=ad=0, poly = model->poly; poly; poly = poly->next ) if ( poly->closed() ) {
		ar = M_PI*(1 - 2.0L/poly->size());
		comp1 = poly->comp[0];
		for ( i=0; i<poly->comp.size(); i++ ) {
			comp = poly->comp[i];
			if ( i+1 < poly->comp.size() ) comp2 = poly->comp[i+1];
			else comp2 = poly->comp[0];
			d1 = comp1->location() - comp->location();
			d2 = comp2->location() - comp->location();
			a = d1.angle(d2);
			aa += a;
			a -= ar;
			ad += a*a;
			n++;
			comp1 = comp;
		}
	}

	if ( n ) {
		aa /= n;
		ad = sqrt(ad/n);
	}
	
	if ( verbose )
		cout << "Polygon angles\t" << n << tab << aa*180.0/M_PI << tab << ad*180.0/M_PI << endl << endl;

	return 0;
}


double			model_average_linklength(Bmodel* model)
{
	int					n(0);
	double				linklength(0);
	Blink*				link;
	
	for ( link = model->link; link; link = link->next, n++ )
		linklength += link->length();
	
	if ( n ) linklength /= n;
	
	return linklength;
}

/**
@brief 	Regularizes a model.
@param 	*model		model structure.
@param 	max_iter	maximum number of iterations.
@param 	distance	reference distance.
@param 	Kdistance	distance strength constant.
@param 	Klink		link strength constant.
@param 	Kpolyangle	angle strength constant.
@param 	Kpolygon	polygon regularity constant.
@param 	Kpolyplane	polygon planarity constant.
@param 	Kpoint		force away from the center-of-mass.
@param 	decay		point force decay constant.
@return int			0.

	Only the first model in the linked list is used.
**/
int			model_regularize(Bmodel* model, int max_iter, double distance, 
				double Kdistance, double Klink, double Kpolyangle, double Kpolygon, 
				double Kpolyplane, double Kpoint, double decay)
{
	if ( !model->select() ) {
		cout << "No model selected for regularization!" << endl << endl;
		return 0;
	}
	
	double			linklength = model_average_linklength(model);
	
	double			max_shift = 0.1*linklength;
	
	if ( max_shift < 1 ) max_shift = 0.1*distance;
	
	if ( max_shift < 1 ) max_shift = 10;
	
	if ( verbose ) {
		cout << "Regularizing a model:" << endl;
		cout << "Reference distance:             " << distance << " A" << endl;
		cout << "Kdistance:                      " << Kdistance << endl;
		cout << "Link length:                    " << linklength << " A" << endl;
		cout << "Klink:                          " << Klink << endl;
		cout << "Kpolyangle:                     " << Kpolyangle << endl;
		cout << "Kpolygon:                       " << Kpolygon << endl;
		cout << "Kpolyplane:                     " << Kpolyplane << endl;
		cout << "Kpoint:                         " << Kpoint << " (" << decay << ")" << endl;
		cout << "Maximum shift per iteration:    " << max_shift << " A" << endl << endl;
	}

	long			iter(0);
	double			Edist0(0), Elink0(0), Eangle0(0), Epoly0(0), Eplane0(0), Epoint0(0), E0(0);
	double			Edist(0), Elink(0), Eangle(0), Epoly(0), Eplane(0), Epoint(0), E(0);

	model_calculate_deviations(model);
	
	Vector3<double>	refcom = model_center_of_mass(model);
	
	if ( verbose & VERB_PROCESS )
		cout << "Iter\tEdist\tElink\tEangle\tEpoly\tEplane\tEpoint\tE" << endl;
	for ( iter=1; iter<=max_iter; iter++ ) {
		model_zero_forces(model);
		Edist = model_lennard_jones_energy(model, Kdistance, distance);
		Elink = model_link_energy(model, Klink);
		Eangle = model_polygon_angle_energy(model, Kpolyangle);
		Epoly = model_polygon_energy(model, Kpolygon);
		Eplane = model_polygon_plane_energy(model, Kpolyplane);
		Epoint = model_point_force(model, refcom, Kpoint, decay);
		E = Edist + Elink + Eangle + Epoly + Eplane + Epoint;
		if ( iter == 1 ) {
			Edist0 = Edist;
			Elink0 = Elink;
			Eangle0 = Eangle;
			Epoly0 = Epoly;
			Eplane0 = Eplane;
			Epoint0 = Epoint;
			E0 = E;
		}
		model_minimize(model, max_shift);
		model_shift(model, refcom - model_center_of_mass(model));
		if ( verbose & VERB_PROCESS )
			cout << iter << tab << Edist << tab << Elink << tab << Eangle << tab 
				<< Epoly << tab << Eplane << tab << Epoint << tab << E << endl;
	}
	
	if ( verbose ) {
		cout << endl << "\tEdistance\tElinklength\tEangledev\tEpolyreg\tEpolyplane\tEpointdist\tE" << endl;
		cout << "Start:\t" << Edist0 << tab << Elink0 << tab << Eangle0 << tab << 
			Epoly0 << tab << Eplane0 << tab << Epoint0 << tab << E0 << endl;
		cout << "End:\t" << Edist << tab << Elink << tab << Eangle << tab << 
			Epoly << tab << Eplane << tab << Epoint << tab << E << endl;
		cout << "Change:\t" << Edist - Edist0 << tab << Elink - Elink0 << tab << 
			Eangle - Eangle0 << tab << Epoly - Epoly0 << tab << Eplane - Eplane0 << tab << 
			Epoint - Epoint0 << tab << E - E0 << endl << endl;
	}
	
	model_calculate_deviations(model);
	
	return 0;
}

