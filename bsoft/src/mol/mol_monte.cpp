/**
@file	mol_monte.cpp
@brief	Functions for a monte carlo metroplis algorithm to energy minimize molecular positions.
@author Bernard Heymann
@date	Created: 20041230
@date	Modified: 20180226
**/

#include "rwmolecule.h"
#include "rwimg.h"
#include "rwmd.h"
#include "mol_md.h"
#include "mol_monte.h"
#include "mol_bonds.h"
#include "mol_map_energy.h"
#include "mol_compare.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "Matrix.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

int			display_energy_headers(Bmd* md)
{
	cout << "#";
	if ( md->Kbond ) cout << "\tEbond";
	if ( md->Kangle ) cout << "\tEangle";
	if ( md->Kvdw ) cout << "\tEvdw";
	if ( md->Kelec ) cout << "\tEelec";
	if ( md->Ksep ) cout << "\tEsep";
	if ( md->Kpoint ) cout << "\tEpoint";
	if ( md->Kmap ) cout << "\tEmap";
	cout << "\tEpot\tdEpot" << endl;
	
	return 0;
}

int			display_energy(long iter, Bmd* md, double dE)
{
	cout << iter;
	if ( md->Kbond ) cout << tab << md->Ebond;
	if ( md->Kangle ) cout << tab << md->Eangle;
	if ( md->Kvdw ) cout << tab << md->Evdw;
	if ( md->Kelec ) cout << tab << md->Eelec;
	if ( md->Ksep ) cout << tab << md->Esep;
	if ( md->Kpoint ) cout << tab << md->Epoint;
	if ( md->Kmap ) cout << tab << md->Emap;
	cout << tab << md->Epot << tab << dE << endl;
	
	return 0;
}

/**
@brief 	Monte Carlo fit of a molecule to a map.  
@param 	*molgroup			molecule group.
@param 	*md					molecular dynamics structure.
@param 	*map				map.
@param 	beta				equivalent of 1/kT.
@param 	max_angle			maximum allowed angular step size.
@param 	max_shift			maximum allowed shift.
@param 	max_iter			maximum number of iterations.
@param 	rigid				type of atom grouping.
@fn 	(Efunc)(Bmolgroup*, Bimage*, Bmd*)	energy function.
@fn 	(Tfunc)(Bmolgroup*, double, double)		transformation function.
@return	Bmolgroup*			transformed coordinates.
**/
Bmolgroup*	monte_carlo_metropolis(Bmolgroup* molgroup, Bmd* md, Bimage* map, 
				double beta, double max_angle, double max_shift, 
				long max_iter, int rigid, 
				double (Efunc)(Bmolgroup*, Bimage*, Bmd*),
				int (Tfunc)(Bmolgroup*, double, double))
{
	int				done(0), accept, select_min(0);
	long			i, isel(0), cycle(0), naccept(0);
	double			r, E, pE, bestE, Z(0), Ze, Zi, Emax;
	double			irm = 1.0/get_rand_max();
	Vector3<double>	gloc = molgroup_center_of_mass(molgroup);
	Bmolgroup*		mgc = molgroup_list_copy(molgroup);
	Bmolgroup*		mg = NULL;
	Bmolgroup*		mgold = NULL;
	Bmolecule*		mol;
	
	double*			Elist = new double[max_iter+1];
	
	long			nmg = count_list((char *) molgroup);
	long			nmol = count_list((char *) molgroup->mol);
	long			nbond = count_list((char *) molgroup->bond);
	Vector3<double>*	iloc = new Vector3<double>[nmol];

	if ( verbose & VERB_RESULT ) {
		cout << "Monte carlo metropolis energy minimization:" << endl;
		cout << "Molecule file:                  " << molgroup->filename << endl;
		cout << "Number of molecule groups:      " << nmg << endl;
		cout << "Number of molecules:            " << nmol << endl;
		if ( map ) cout << "Map file:                       " << map->file_name() << endl;
		cout << "Coordinate minima and maxima:   " << molgroup->min << " " << molgroup->max << endl;
		cout << "Bounding box:                   " << molgroup->box << endl;
		cout << "Kbond:                          " << md->Kbond << endl;
		cout << "Kangle:                         " << md->Kangle << endl;
		cout << "Kelectrostatic:                 " << md->Kelec << endl;
		cout << "Kvdw:                           " << md->Kvdw << endl;
		cout << "Kmap:                           " << md->Kmap << endl;
		cout << "Kseparation:                    " << md->Ksep << endl;
		cout << "Kpoint:                         " << md->Kpoint << " (" << md->pointdecay << ")" << endl;
		cout << "location:                       " << md->point << endl;
		cout << "beta constant:                  " << beta << endl;
		cout << "Cutoff distance:                " << md->cutoff << endl;
		cout << "Separation distance:            " << md->sepdist << endl;
		cout << "Maximum shift:                  " << max_shift << endl;
		cout << "Maximum angular increment:      " << max_angle*180.0/M_PI << endl;
		cout << "Steps along bond:               " << md->bondsteps << endl;
		cout << "Number of bonds:                " << nbond << endl;
		if ( rigid < 1 ) cout << "Treating the whole ensemble as a rigid body" << endl;
		else if ( rigid == 1 ) cout << "Treating each file with a molecule group as a rigid body" << endl;
		else if ( rigid == 2 ) cout << "Treating each molecule as a rigid body" << endl;
		else cout << "Modifying each atom's coordinates" << endl;
		cout << endl;
		md_calculate_deviations(molgroup, md->wrap);
	}
	
	if ( rigid < 1 ) {
		for ( mg = mgc; mg; mg = mg->next ) mg->select = "ALL";
	} else if ( rigid == 1 ) {
		for ( mg = mgc; mg; mg = mg->next ) {
			mg->sel = 0;
			mg->select = "FILE";
		}
	} else if ( rigid == 2 ) {
		for ( i=0, mol = mgc->mol; mol; mol = mol->next, ++i ) {
			mol->sel = 0;
			iloc[i] = mol_center_of_mass(mol);
		}
		mgc->select = "MOL";
	} else mgc->select = "ATOM";

	pE = (Efunc)(mgc, map, md);

	bestE = E = pE;
	Elist[0] = pE;
	
	if ( verbose & VERB_RESULT ) {
		display_energy_headers(md);
		display_energy(0, md, E-pE);
	}
	while ( !done && ( cycle < max_iter ) ) {
		cycle++;
		mgold = molgroup_list_copy(mgc);
		if ( rigid == 1 ) {
			if ( select_min ) {
				for ( i=0, Emax=-1e30, mg=mgc; mg; mg=mg->next, i++ ) {
					if ( Emax < mg->fom ) {
						Emax = mg->fom;
						isel = i;
					}
				}
			} else {
				isel = (long) ((nmg-0.001)*irm*random());
			}
			for ( i=0, mg=mgc; mg; mg=mg->next, i++ ) {
				if ( i == isel ) mg->sel = 1;
				else mg->sel = 0;
			}
		} else if ( rigid == 2 ) {
			if ( select_min ) {
				for ( i=0, Emax=-1e30, mol=mgc->mol; mol; mol=mol->next, i++ ) {
					if ( Emax < mol->fom ) {
						Emax = mol->fom;
						isel = i;
					}
				}
			} else {
				isel = (long) ((nmol-0.001)*irm*random());
			}
			for ( i=0, mol=mgc->mol; mol; mol=mol->next, i++ ) {
				if ( i == isel ) mol->sel = 1;
				else mol->sel = 0;
			}
		}
		(Tfunc)(mgc, max_angle, max_shift);
		E = (Efunc)(mgc, map, md);
		Elist[cycle] = E;
//		display_energy(cycle, md, E-pE);
		if ( bestE > E ) {
			bestE = E;
			molgroup_list_kill(molgroup);
			molgroup = molgroup_list_copy(mgc);
			if ( verbose & VERB_RESULT )
				display_energy(cycle, md, E-pE);
		}
		accept = 1;
		if ( E > 1e30 ) accept = 0;		// High energy indicates the molecule might be outside the box
		if ( E > pE ) {
			r = random()*irm;
			if ( exp(beta*(pE-E)) < r ) accept = 0;
		}
		if ( accept ) {
			pE = E;
			naccept++;
			molgroup_list_kill(mgold);
		} else {
			molgroup_list_kill(mgc);
			mgc = mgold;
		}
//		molgroup_test_if_within_box(mgc, mgc->min, mgc->max);
	}
	
	molgroup_list_kill(mgc);
	
	for ( i=0, Z=Ze=0; i<=cycle; i++ ) {
		Zi = exp(-beta*Elist[i]);
		Z += Zi;
		if ( Elist[i] - bestE < 0.001 ) Ze += Zi;
	}
	if ( Z <= 0 ) Z = 1;
	
	delete[] Elist;
	
	md->Epot = bestE;
	
	gloc = molgroup_center_of_mass(molgroup) - gloc;
	
	if ( verbose & VERB_RESULT ) {
		cout << "Total shift:                        " << gloc << endl;
		for ( i=0, mol = molgroup->mol; mol; mol = mol->next, ++i ) {
			iloc[i] = mol_center_of_mass(mol) - iloc[i];
			E = mol_map_energy(mol, map, md->Kmap);
			cout << tab << mol->id << tab << iloc[i] << tab << E << endl;
		}
		cout << "Fraction of perturbations accepted: " << naccept*1.0/cycle << endl;
		cout << "Final energy:                       " << bestE << endl;
		cout << "Partition function:                 " << Z << endl;
		cout << "Probability of final orientation:   " << Ze/Z << endl << endl;
		md_calculate_deviations(molgroup, md->wrap);
	}
	
	delete[] iloc;
	
	return molgroup;
}

/**
@brief 	Generates multiple versions of a molecule at different locations.  
@param 	*molgroup		molecule group.
@param 	*pmask			mask to get limit grid positions.
@param 	grid_sampling	sampling for grid search.
@param 	filename		output base file name.
@return Bmolgroup*		linked list of molecule groups.

	The locations of the molecules are defined on the intersection of a 
	grid and a user-defined mask.
	The numbered output filename is also set.
	Note: the molecule group minima and maxima should be set to the 
	volume to be searched before calling this function.

**/
Bmolgroup*	molgroup_generate_masked_grid_list(Bmolgroup* molgroup, Bimage* pmask,
				Vector3<double> grid_sampling, Bstring filename)
{
	Bstring			basename("t.pdb");
	long			i, x, y, z, n(0);
	Vector3<double>	loc;
	Vector3<double>	bbmin = molgroup->min + grid_sampling;
	Vector3<double>	bbmax = molgroup->max - grid_sampling;
	Vector3<double>	com = molgroup_center_of_mass(molgroup);
	Bmolgroup*		mglist = NULL;
	Bmolgroup*		mgc = NULL;
	
	if ( verbose ) {
		cout << "Generating molecules on a masked grid:" << endl;
		cout << "Bounding box:                   " << bbmin << " " << bbmax << endl;
		cout << "Grid sampling:                  " << grid_sampling << endl;
		if ( pmask ) cout << "Mask file:                      " << pmask->file_name() << endl;
	}
	
	int				accept_flag;
	
	for ( loc[2]=bbmin[2]; loc[2]<=bbmax[2]; loc[2]+=grid_sampling[2] ) {
		for ( loc[1]=bbmin[1]; loc[1]<=bbmax[1]; loc[1]+=grid_sampling[1] ) {
			for ( loc[0]=bbmin[0]; loc[0]<=bbmax[0]; loc[0]+=grid_sampling[0] ) {
				accept_flag = 1;
				if ( pmask ) {
					x = (long) (loc[0]/pmask->sampling(0)[0] + pmask->image->origin()[0] + 0.5);
					y = (long) (loc[1]/pmask->sampling(0)[1] + pmask->image->origin()[1] + 0.5);
					z = (long) (loc[2]/pmask->sampling(0)[2] + pmask->image->origin()[2] + 0.5);
					i = pmask->index(0, x, y, z, 0);
					if ( i < 0 ) accept_flag = 0;
					else if ( (*pmask)[i] < 0.5 ) accept_flag = 0;
				}
				if ( accept_flag ) {
					n++;
					if ( mgc ) {
						mgc->next = molgroup_copy(molgroup);
						mgc = mgc->next;
					} else {
						mgc = molgroup_copy(molgroup);
						mglist = mgc;
					}
					molgroup_coor_shift(mgc, (loc - com));
					molgroup_set_box_to_map_boundaries(mgc, pmask);
					if ( filename.length() )
						mgc->filename = filename.pre_rev('.') + Bstring(n, "_%06d.") + filename.post_rev('.');
					else
						mgc->filename = basename.pre_rev('.') + Bstring(n, "_%06d.") + basename.post_rev('.');
				}
			}
		}
	}

	if ( verbose )
		cout << "Molecules generated:            " << n << endl << endl;
		
	return mglist;
}

Bmodel*		molgroup_generate_masked_grid_list(Bmolgroup* molgroup, 
				Vector3<double> grid_sampling, Bimage* pmask)
{
	long			i, x, y, z, j(0);
	Vector3<double>	loc;
	Vector3<double>	bbmin = molgroup->min + grid_sampling;
	Vector3<double>	bbmax = molgroup->max - grid_sampling;
	string			id("Grid");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	string			comptype("1");
	Bcomptype*		ct = model->add_type(comptype);
	ct->file_name(molgroup->filename.str());
	
	if ( verbose ) {
		cout << "Generating models on a masked grid:" << endl;
		cout << "Molecule file name:             " << molgroup->filename << endl;
		cout << "Bounding box:                   " << bbmin << " " << bbmax << endl;
		cout << "Grid sampling:                  " << grid_sampling << endl;
		if ( pmask ) cout << "Mask file:                      " << pmask->file_name() << endl;
	}
	
	int				accept_flag;
	
	for ( loc[2]=bbmin[2]; loc[2]<=bbmax[2]; loc[2]+=grid_sampling[2] ) {
		for ( loc[1]=bbmin[1]; loc[1]<=bbmax[1]; loc[1]+=grid_sampling[1] ) {
			for ( loc[0]=bbmin[0]; loc[0]<=bbmax[0]; loc[0]+=grid_sampling[0] ) {
				accept_flag = 1;
				if ( pmask ) {
					x = (long) (loc[0]/pmask->sampling(0)[0] + pmask->image->origin()[0] + 0.5);
					y = (long) (loc[1]/pmask->sampling(0)[1] + pmask->image->origin()[1] + 0.5);
					z = (long) (loc[2]/pmask->sampling(0)[2] + pmask->image->origin()[2] + 0.5);
					i = pmask->index(0, x, y, z, 0);
					if ( i < 0 ) accept_flag = 0;
					else if ( (*pmask)[i] < 0.5 ) accept_flag = 0;
				}
				if ( accept_flag ) {
//					comp = (Bcomponent *) add_item((char **) &comp, sizeof(Bcomponent));
//					if ( !model->comp ) model->comp = comp;
//					comp->identifier(Bstring(++id, "%d"));
					if ( comp ) comp = comp->add(++j);
					else comp = model->add_component(++j);
					comp->type(ct);
					comp->location(loc);
					comp->view()[2] = 1;
					comp->select(1);
				}
			}
		}
	}

	if ( verbose )
		cout << "Components generated:           " << id << endl << endl;
	
	return model;
}

/**
@brief 	Generates multiple versions of a molecule at different locations.  
@param 	*molgroup		molecule group.
@param 	angle_step		angular step size in radians.
@param 	filename		output base file name.
@param 	whole			treat the whole ensemble as a rigid body.
@return Bmolgroup*		linked list of molecule groups.

	The molecule group is rotated in place to give all orientations
	with a given angle step size between the views.
	The numbered output filename is also set.

**/
Bmolgroup*	molgroup_generate_orientation_list(Bmolgroup* molgroup,
				double angle_step, Bstring filename, int whole)
{
	Bstring			basename("t.pdb");
	int				n(0), nphi;
	double			theta, phi, phi_step;
	Bmolgroup*		mglist = NULL;
	Bmolgroup*		mgc = NULL;
	Bmolecule*		mol;
	Transform		t;
	t.origin = molgroup_center_of_mass(molgroup);

	if ( verbose & VERB_FULL )
		cout << "Theta\tAxis.x\tAxis.y" << endl;
	for ( theta=0; theta<M_PI+angle_step/2; theta+=angle_step ) {
		nphi = (int) (2*M_PI*sin(theta)/angle_step + 0.5);	// Number of views at this theta
		phi_step = 3*M_PI;
		if ( nphi ) phi_step = M_PI*2.0/nphi;
		for ( phi=-M_PI; phi<M_PI-phi_step/2; phi+=phi_step ) {
			n++;
			t.axis = Vector3<double>(cos(phi), sin(phi), 0);
			t.angle = theta;
			if ( verbose & VERB_FULL )
				cout << theta*180.0/M_PI << tab << t.axis[0] << tab << t.axis[1] << endl;
			if ( mglist ) {
				mgc->next = molgroup_copy(molgroup);
				mgc = mgc->next;
			} else mglist = mgc = molgroup_copy(molgroup);
			if ( whole ) {
				molgroup_coor_rotate(mgc, t);
			} else {
				for ( mol = mgc->mol; mol; mol = mol->next ) {
					t.origin = mol_center_of_mass(mol);
					mol_coor_rotate(mol, t);
				}
			}
			if ( filename.length() )
				mgc->filename = filename.pre_rev('.') + Bstring(n, "_%06d.") + filename.post_rev('.');
			else
				mgc->filename = basename.pre_rev('.') + Bstring(n, "_%06d.") + basename.post_rev('.');
		}
	}
		
	if ( verbose )
		cout << "Molecules generated:            " << n << endl << endl;
		
	return mglist;
}

Bmodel*		molgroup_generate_orientation_list(Bmolgroup* molgroup, double angle_step)
{
	int				i(0), nphi;
	double			theta, phi, phi_step;
	Vector3<double>	com = molgroup_center_of_mass(molgroup);
	string			id("Orientations");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	string			comptype("1");
	Bcomptype*		ct = model->add_type(comptype);
	ct->file_name(molgroup->filename.str());

	for ( theta=0; theta<M_PI+angle_step/2; theta+=angle_step ) {
		nphi = (int) (2*M_PI*sin(theta)/angle_step + 0.5);	// Number of views at this theta
		phi_step = 3*M_PI;
		if ( nphi ) phi_step = M_PI*2.0/nphi;
		for ( phi=-M_PI; phi<M_PI-phi_step/2; phi+=phi_step ) {
//			comp = (Bcomponent *) add_item((char **) &comp, sizeof(Bcomponent));
//			if ( !model->comp ) model->comp = comp;
//			comp->identifier(Bstring(++id, "%d"));
			if ( comp ) comp = comp->add(++i);
			else model->add_component(++i);
			comp->type(ct);
			comp->location(com);
			comp->view(View2<float>(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), 0));
			comp->select(1);
		}
	}
		
	if ( verbose )
		cout << "Components generated:           " << id << endl << endl;
		
	return model;
}

/**
@brief 	Monte Carlo fit of a molecule to a list of locations and orientations.  
@param 	*molgroup		molecule group.
@param 	*model			model parameters (modified).
@param 	*md				molecular dynamics structure.
@param 	*map			map.
@param 	beta			equivalent of 1/kT.
@param 	max_angle		maximum allowed angular step size.
@param 	max_shift		maximum allowed shift.
@param 	max_iter		maximum number of iterations.
@param 	rigid			selection of mode.
@return Bmolgroup*		new molecule group list.

	For each model, the molecule group is first transformed to the model 
	location and orientation. This is then refined with the Monte Carlo
	Metropolis algorithm and the model parameters updated.

**/
int			mcm_molecule_list(Bmolgroup* molgroup, Bmodel* model, Bmd* md, Bimage* map, double beta, 
				double max_angle, double max_shift, long max_iter, int rigid)
{
	int				i;
	Bmolgroup*		mc = NULL;
	Bcomponent*		comp = NULL;
	Quaternion		q;
	Transform		t;
	t.origin = molgroup_center_of_mass(molgroup);
	
	model->mapfile(map->file_name());
	
	for ( i=1, comp = model->comp; comp; comp = comp->next, i++ ) {
//		comp->type = model_add_type_by_filename(model, molgroup->filename, 0);
		comp->type(model->add_type(molgroup->filename.str(), 0));
		mc = molgroup_copy(molgroup);
//		q = quaternion_from_view(comp->view);
		q = comp->view().quaternion();
//		t = transform_from_quaternion(q);
		t = Transform(q);
		t.origin = molgroup_center_of_mass(mc);
		t.trans = comp->location() - t.origin;
		molgroup_coor_rotate(mc, t);
		if ( map ) molgroup_set_box_to_map_boundaries(mc, map);
		mc = monte_carlo_metropolis(mc, md, map, beta, max_angle, max_shift, 
			max_iter, rigid, monte_rigid_body_fit_energy, molgroup_rigid_body_transform);
		comp->location(molgroup_center_of_mass(mc));
		t = molgroup_find_transformation(mc, molgroup);
//		comp->view = view_from_angle_and_axis3(t.angle, t.axis);
		comp->view(View2<float>(t.angle, t.axis));
		comp->FOM(exp(-0.14*md->Epot));
		molgroup_kill(mc);
		if ( verbose )
			cout << i << ":\t" << comp->location() << tab << comp->view() << tab << md->Epot << endl;
	}
	
	return i;
}

/**
@brief 	Monte Carlo fit of a set of molecule groups.  
@param 	*molgroup		molecule group.
@param 	*model			model parameters (modified).
@param 	*md				molecular dynamics structure.
@param 	*map			map.
@param 	beta			equivalent of 1/kT.
@param 	max_angle		maximum allowed angular step size.
@param 	max_shift		maximum allowed shift.
@param 	max_iter		maximum number of iterations.
@param 	rigid			selection of mode.
@return	Bmolgroup*		new molecule group list.

	For each model, the molecule group is first transformed to the model 
	location and orientation. The whole ensemble is then refined with the 
	Monte Carlo Metropolis algorithm and the model parameters updated.

**/
int			mcm_molecule_groups(Bmolgroup* molgroup, Bmodel* model, Bmd* md, Bimage* map, double beta, 
				double max_angle, double max_shift, long max_iter, int rigid)
{
	int				i = 1;
	Bmolgroup*		mglist = NULL;
	Bmolgroup*		mc = NULL;
	Bcomponent*		comp = NULL;
	Quaternion		q;
	Transform		t;
	t.origin = molgroup_center_of_mass(molgroup);

	model->mapfile(map->file_name());
	
//	Bstring			filename;
	
	for ( comp = model->comp; comp; comp = comp->next ) {
//		comp->type = model_add_type_by_filename(model, molgroup->filename, 0);
		comp->type(model->add_type(molgroup->filename.str(), 0));
		if ( !mglist ) {
			mglist = mc = molgroup_copy(molgroup);
		} else {
			mc->next = molgroup_copy(molgroup);
			mc = mc->next;
		}
//		q = quaternion_from_view(comp->view);
		q = comp->view().quaternion();
//		t = transform_from_quaternion(q);
		t = Transform(q);
		t.origin = molgroup_center_of_mass(mc);
		t.trans = comp->location() - t.origin;
		molgroup_coor_rotate(mc, t);
		if ( map ) molgroup_set_box_to_map_boundaries(mc, map);
//		filename = Bstring(i++, "test_%03d.pdb");
//		write_molecule(filename, mc);
	}
	
	mglist = monte_carlo_metropolis(mglist, md, map, beta, max_angle, max_shift, 
			max_iter, rigid, monte_rigid_body_fit_energy, molgroup_rigid_body_transform);
	
	for ( i=1, comp = model->comp, mc = mglist; comp && mc; comp = comp->next, mc = mc->next, i++ ) {
		comp->location(molgroup_center_of_mass(mc));
		t = molgroup_find_transformation(mc, molgroup);
//		comp->view = view_from_angle_and_axis3(t.angle, t.axis);
		comp->view(View2<float>(t.angle, t.axis));
		comp->FOM(exp(-0.14*md->Epot));
		if ( verbose )
			cout << i << ":\t" << comp->location() << tab << comp->view() << tab << md->Epot << endl;
	}

	molgroup_list_kill(mglist);

	return i;
}

/**
@brief 	Sets the box in a molgroup to that defined by the map boundaries.  
@param 	*molgroup		molecule group.
@param 	*map			map to get boundaries from.
@return  int			0.
**/
int			molgroup_set_box_to_map_boundaries(Bmolgroup* molgroup, Bimage* map)
{
	if ( !map ) return 0;
	
	molgroup->min[0] = -map->image->origin()[0]*map->sampling(0)[0];
	molgroup->min[1] = -map->image->origin()[1]*map->sampling(0)[1];
	molgroup->min[2] = -map->image->origin()[2]*map->sampling(0)[2];
	molgroup->max[0] = molgroup->min[0] + map->sizeX()*map->sampling(0)[0];
	molgroup->max[1] = molgroup->min[1] + map->sizeY()*map->sampling(0)[1];
	molgroup->max[2] = molgroup->min[2] + map->sizeZ()*map->sampling(0)[2];
	molgroup->box = molgroup->max - molgroup->min;
	
	return 0;
}

/**
@brief 	Tests if a molecule overlaps with a defined box.  
@param 	*molgroup		molecule group.
@param 	min				start of box.
@param 	max				end of box.
@return long			0.
**/
long		molgroup_test_if_within_box(Bmolgroup* molgroup, Vector3<double> min, Vector3<double> max)
{
	Bmolecule*		mol;	
	Bresidue*		res;
	Batom*			atom;
	long			i, natom(0), natom_in(0);
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				i = ( atom->coord[0] >= min[0] && atom->coord[0] <= max[0] );
				i += ( atom->coord[1] >= min[1] && atom->coord[1] <= max[1] );
				i += ( atom->coord[2] >= min[2] && atom->coord[2] <= max[2] );
				if ( i > 2 ) natom_in++;
				natom++;
			}
		}
	}
	
	if ( verbose )
		cout << "Fraction within bounding box:   " << natom_in*1.0/natom << endl;
	
	return natom_in;
}

/**
@brief 	Calculates the potential energy for rigid body fitting.
@param 	*molgroup		molecular structure.
@param 	*map			density map.
@param 	*md				molecular dynamics structure.
@return double			potential energy.

	The energy is the sum of the overlap, map, and point force energies.

**/
double		monte_rigid_body_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md)
{
	md->Epot = molgroup_atom_overlap(molgroup, md);
	
	if ( map && md->Kmap ) {
		md->Emap = molgroup_map_energy(molgroup, map, md->Kmap);
		md->Epot += md->Emap;
	}
	
	if ( md->Kpoint ) {
		md->Epoint = md->Kpoint*md->point.distance(molgroup_center_of_mass(molgroup));
		md->Epot += md->Epoint;
	}
	
	if ( md->Kbond ) {
		md->Ebond = md_bond_forces(molgroup, md->Kbond, md->wrap);
		md->Epot += md->Ebond;
	}
	
	return md->Epot;
}

/**
@brief 	Calculates the potential energy for fitting atoms to a map.
@param 	*molgroup		molecular structure.
@param 	*map			density map.
@param 	*md				molecular dynamics structure.
@return double			potential energy.

	The energy is the sum of the bond, angle, and map energies.

**/
double		monte_atom_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md)
{
	md_zero_forces(molgroup);
	
	md->Ebond = md_bond_forces(molgroup, md->Kbond, md->wrap);
	
	md->Eangle = md_angular_forces(molgroup, md->Kangle, md->wrap);
	
	md_nonbonded_forces(molgroup, md);
	
	md->Epoint = 0;
	if ( md->Kpoint )
		md->Epoint = md->Kpoint*md->point.distance(molgroup_center_of_mass(molgroup));
	
	md->Emap = 0;
	if ( map ) md->Emap = molgroup_map_energy(molgroup, map, md->Kmap);
	
	md->Epot = md->Ebond + md->Eangle + md->Evdw + md->Eelec + md->Epoint + md->Emap;
	
	return md->Epot;
}

/**
@brief 	Calculates the potential energy for fitting bonds to a map.
@param 	*molgroup		molecular structure.
@param 	*map			density map.
@param 	*md				molecular dynamics structure.
@return double			potential energy.

	The energy is the sum of the bond, angle, and map energies.

**/
double		monte_bond_fit_energy(Bmolgroup* molgroup, Bimage* map, Bmd* md)
{
	md_zero_forces(molgroup);
	
	md->Ebond = md_bond_forces(molgroup, md->Kbond, md->wrap);
	
	md->Eangle = md_angular_forces(molgroup, md->Kangle, md->wrap);
	
	md_nonbonded_forces(molgroup, md);
	
	md->Epoint = 0;
	if ( md->Kpoint )
		md->Epoint = md->Kpoint*md->point.distance(molgroup_center_of_mass(molgroup));
	
	md->Emap = 0;
	if ( map ) md->Emap = molgroup_bond_fit_map_energy(molgroup, map, md->Kmap, md->bondsteps);
	
	md->Epot = md->Ebond + md->Eangle + md->Evdw + md->Eelec + md->Epoint + md->Emap;
	
	return md->Epot;
}


/**
@brief 	Calculates an energy term based on atom overlap.
@param 	*molgroup		molecular structure.
@param 	*md				molecular dynamics structure.
@return double			total overlap energy.

	The energy is defined as linear decay to the reference separation distance
	and zero beyond:
		Esep = Ksep * (1 - d/dsep)  for  d < dsep, zero otherwise

**/
double		molgroup_atom_overlap(Bmolgroup* molgroup, Bmd* md)
{
	if ( md->Ksep <= 0 ) return 0;
	
	long			n(0), i, x, y, z;
	Vector3<double>	box = molgroup->box;
	Vector3<double>	sampling(md->sepdist, md->sepdist, md->sepdist);
	Vector3<int>	size;
	size[0] = (int) (box[0]/sampling[0] + 0.001);
	size[1] = (int) (box[1]/sampling[1] + 0.001);
	size[2] = (int) (box[2]/sampling[2] + 0.001);
	size = size.max(1);
	sampling[0] = box[0]/size[0] + 0.001;
	sampling[1] = box[1]/size[1] + 0.001;
	sampling[2] = box[2]/size[2] + 0.001;
	long			boxsize = (long) size.volume();

	Bmolgroup*		mg;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;

	Bimage*			p = new Bimage(Float, TSimple, size, 1);
	float*			data = (float *) p->data_pointer();
	p->next = new Bimage(Float, TSimple, p->size(), p->images());
	float*			dataone = (float *) p->next->data_pointer();
	
	for ( mg = molgroup; mg; mg = mg->next ) {
		for ( mol = mg->mol; mol; mol = mol->next ) {
			for ( res = mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					x = (long) ((atom->coord[0] - mg->min[0])/sampling[0]);
					y = (long) ((atom->coord[1] - mg->min[1])/sampling[1]);
					z = (long) ((atom->coord[2] - mg->min[2])/sampling[2]);
					if ( x>=0 && x<size[0] && y>=0 && y<size[1] && z>=0 && z<size[2] ) {
						i = p->index(x, y, z);
						dataone[i] = 1;
					}
				}
			}
			if ( !mg->select.contains("FILE") ) {
				for ( i=0; i<boxsize; i++ ) {
					data[i] += dataone[i];
					dataone[i] = 0;
				}
			}
		}
		if ( mg->select.contains("FILE") ) {
			for ( i=0; i<boxsize; i++ ) {
				data[i] += dataone[i];
				dataone[i] = 0;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		Bstring			filename("overlap.map");
		write_img(filename, p, 0);
	}
	
	md->Esep = 0;
	
	for ( i=n=0; i<boxsize; i++ ) {
		if ( data[i] ) {
			n++;
			if ( data[i] > 1 ) md->Esep += data[i] - 1;
		}
	}
	
	md->Esep *= md->Ksep * md->sepdist/n;
	
	delete p;
	
	return md->Esep;
}

/**
@brief 	Randomly transforms a molecule or molecule group.
@param 	*molgroup		molecular structure.
@param 	max_angle		maximum rotation angle.
@param 	shift_std		gaussian length for shift vector.
@return int				0.

	The transformation is calculted as a random angular rotation and a
	random shift. The shift is sampled from a random vector with a
	gaussian length distribution.

**/
int			molgroup_rigid_body_transform(Bmolgroup* molgroup, double max_angle, double shift_std)
{
	double			irm = 1.0/get_rand_max();
	Bmolecule*		mol;

	Transform		t;
	
	t.axis = vector3_random(-1, 1);
	t.axis.normalize();
	
	t.angle = max_angle*(random()*irm*2 - 1);
	
	t.trans = vector3_random_gaussian(0, shift_std);
	
	if ( molgroup->select.contains("ALL") ) {
		t.origin = molgroup_center_of_mass(molgroup);
		molgroup_coor_rotate(molgroup, t);
	} else if ( molgroup->select.contains("FILE") ) {
		for ( ; molgroup; molgroup=molgroup->next ) if ( molgroup->sel ) {
			t.origin = molgroup_center_of_mass(molgroup);
			molgroup_coor_rotate(molgroup, t);
			molgroup->sel = 0;
		}
	} else {
		for ( mol=molgroup->mol; mol; mol=mol->next ) if ( mol->sel ) {
			t.origin = mol_center_of_mass(mol);
			mol_coor_rotate(mol, t);
			mol->sel = 0;
		}
	}
	
	return 0;
}

/**
@brief 	Move atoms random distances down the energy gradient.
@param 	*molgroup		molecular structure.
@param 	max_angle		(not used).
@param 	max_shift		maximum shift for each atom.
@return double			0.

	The distance of movement is limited to the maximum shift.

**/
int			molgroup_move_atoms_down_energy(Bmolgroup* molgroup, double max_angle, double max_shift)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	double			irm = 1.0/get_rand_max();
	Vector3<double>	shift;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				shift = atom->F * (random()*irm);
				shift = shift.min(max_shift);
				atom->coord += shift;
			}
		}
	}
	
	return 0;
}

