/**
@file	mol_alpha.cpp
@brief	Functions to make and analyze alpha helices.
@author Bernard Heymann
@date	Created: 20050315
@date	Modified: 20111030
**/

#include "rwmolecule.h"
#include "mol_alpha.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "Bimage.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "Matrix3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
Records for the template residue in an alpha helix:
ATOM     39  N   ALA S   7       1.378  -0.646  -0.926  1.00  0.00          N  0.000
ATOM     40  CA  ALA S   7       2.303   0.000  -0.000  1.00  0.00          C  0.000
ATOM     41  C   ALA S   7       1.599   0.860   1.028  1.00  0.00          C  0.000
ATOM     42  O   ALA S   7       2.020   1.002   2.171  1.00  0.00          O  0.000
ATOM     43  CB  ALA S   7       3.260   0.864  -0.824  1.00  0.00          C  0.000
ATOM     44  H   ALA S   7       1.279  -0.316  -1.865  1.00  0.00          H  0.000
*/

#define NTEMPATOM	6		// Number of atoms defined in template
#define	DAXIS		2.303	// Distance between Ca and helical axis
#define HRISE		1.5		// Helical rise per residue
#define HANGLE		100		// Helical rotation per residue in degrees

struct AlphaTemplate {
	char atomtype[8];
	double coord[3];
} alpha_template[] = {
	{"N",  {1.378, -0.646, -0.926}},
	{"CA", {2.303,  0.000,  0.000}},
	{"C",  {1.599,  0.860,  1.028}},
	{"O",  {2.020,  1.002,  2.171}},
	{"CB", {3.260,  0.864, -0.824}},
	{"H",  {1.279, -0.316, -1.865}}
} ;

/**
@brief 	Generates an alpha helix of the desired length.  
@param 	length				number of alanines to generate.
@return  	int						0.
**/
Bmolgroup*	molgroup_generate_alpha_helix(int length)
{
	Bmolgroup*	molgroup = molgroup_init();
	
	molgroup->mol = mol_generate_alpha_helix(length);
	
	return molgroup;
}

/**
@brief 	Generates an alpha helix of the desired length.  
@param 	length				number of alanines to generate.
@return  	int						0.
**/
Bmolecule*	mol_generate_alpha_helix(int length)
{
	int				i, j, n;
	char			molname[8] = "A";
	char			restype[8] = "ALA";
	Bmolecule*		mol = NULL;
	mol = molecule_add(&mol, molname);
	
	Bsecondary*		sec = (Bsecondary *) add_item((char **) &mol->sec, sizeof(Bsecondary));
	Bresidue*		res = NULL;
	Batom*			atom = NULL;
	Matrix3			mat(1);
	Vector3<double>	shift, axis(0,0,1);
	
	if ( verbose & VERB_PROCESS )
		cout << "Generating an aplha helix of length " << length << endl << endl;
	
	for ( i=n=0; i<length; i++ ) {
		res = residue_add(&mol->res, restype);
		res->num = i + 1;
		mat = Matrix3(axis, HANGLE*i*M_PI/180.0);
		shift = Vector3<double>(0, 0, HRISE*i);
		for ( j=0; j<NTEMPATOM; j++ ) {
			atom = atom_add(&res->atom, alpha_template[j].atomtype);
			atom->coord[0] = alpha_template[j].coord[0];
			atom->coord[1] = alpha_template[j].coord[1];
			atom->coord[2] = alpha_template[j].coord[2];
			atom->coord = mat * atom->coord;
			atom->coord += shift;
			atom->num = ++n;
			atom->sel = 1;
		}
	}
	
	sec->num = 1;
	sec->id[0] = 'A';
	sec->first = mol->res;
	sec->last = res;
	
	return mol;
}

/**
@brief 	Sets a residue range to an alpha helix in all molecules.  
@param 	*molgroup		molecule group.
@param 	helix_start			first residue in helix.
@param 	helix_end			last residue in helix.
@return  	int						0.
**/
int			molgroup_set_alpha_helix(Bmolgroup* molgroup, int helix_start, int helix_end)
{
	Bmolecule*	mol;
	Bsecondary*	sec;
	Bresidue*	res1, *res2;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( res1 = mol->res; res1 && res1->num < helix_start; res1 = res1->next ) ;
		if ( res1 ) {
			for ( res2 = res1; res2->next && res2->next->num < helix_end; res2 = res2->next ) ;
			if ( res2->next ) {
				if ( res2->next->num == helix_end ) res2 = res2->next;
				if ( res2 ) {
					sec = (Bsecondary *) add_item((char **) &mol->sec, sizeof(Bsecondary));
					sec->first = res1;
					sec->last = res2;
				}
			}
		}
	}
	
	return 0;
}


/**
@brief 	Calculates the centers and orientations of all alpha helices.  
@param 	*molgroup		molecule group.
@return  	int						0.
**/
int			molgroup_find_helical_axes(Bmolgroup* molgroup)
{
	Bmolecule*		mol;
	Bsecondary*		sec;
	int				i;
	Vector3<double>	c, a;
	
	if ( verbose & VERB_PROCESS )
		cout << "Mol\tHelix\tFirst\tLast\tLength\tCenter\t\t\tAxis" << endl;
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for ( i=1, sec = mol->sec; sec; sec = sec->next, i++ ) {
			c = alpha_find_center(sec->first, sec->last);
			a = alpha_find_orientation(sec->first, sec->last);
			if ( verbose & VERB_PROCESS )
				cout << mol->id << tab << i << tab << sec->first->num << tab << sec->last->num << tab << 
					sec->last->num - sec->first->num + 1 << tab << 
					c[0] << tab << c[1] << tab << c[2] << tab << 
					a[0] << tab << a[1] << tab << a[2] << endl;
		}
	}
	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	return 0;
}

/**
@brief 	Calculates the center of an alpha helix.  
@param 	*resfirst		first residue in helix.
@param 	*reslast		last residue in helix.
@return  	Vector3<double>			center.
**/
Vector3<double>		alpha_find_center(Bresidue* resfirst, Bresidue* reslast)
{
	int				i, n;
	Vector3<double>	ca[3];
	Bresidue*		res;
	Batom*			atom;
	
	for ( n=1, res = resfirst; res && res != reslast; res = res->next ) n++;
	for ( i=0, res = resfirst; res && i<n/2-1; res = res->next ) i++;
	
	for ( i=0; res && i<3; res = res->next, i++ ) {
		for ( atom = res->atom; atom && strncmp(atom->type, "CA", 2); atom = atom->next ) ;
		ca[i] = atom->coord;
	}
	
	Vector3<double>		center = point_on_helix_axis(ca[0], ca[1], ca[2]);
	
	if ( verbose & VERB_FULL )
		cout << "Helical center:                 " << center << " (" << n/2 << ")" << endl;
	
	return center;
}

/**
@brief 	Calculates the orientation of an alpha helix.  
@param 	*resfirst		first residue in helix.
@param 	*reslast		last residue in helix.
@return  	Vector3<double>			axis.
**/
Vector3<double>	alpha_find_orientation(Bresidue* resfirst, Bresidue* reslast)
{
	int				i, nca;
	Bresidue*		res;
	Batom*			atom;
	Vector3<double>	ca[1000], v1, v2, va, vas;
	
	for ( nca=0, res = resfirst; res != reslast->next; res = res->next )
		for ( atom = res->atom; atom; atom = atom->next )
			if ( strncmp(atom->type, "CA", 2) == 0 )
				ca[nca++] = atom->coord;
	
	for ( i=0; i<nca-3; i++ ) {
		v1 = (ca[i] + ca[i+2]) - (ca[i+1] * 2);
		v2 = (ca[i+1] + ca[i+3]) - (ca[i+2] * 2);
		v1.normalize();
		v2.normalize();
		va = v1.cross(v2);
		va.normalize();
		vas += va;
	}
	
	vas.normalize();
	
	i = nca/2;
	va = point_on_helix_axis(ca[i-1], ca[i], ca[i+1]);
	
	if ( verbose & VERB_FULL ) {
		cout << "Helical center:                 " << va << endl;
		cout << "Helical axis:                   " << vas << endl << endl;
	}
	
	return vas;
}

/**
@brief 	Calculates the orientation of an alpha helix.  
@param 	*mol			molecule.
@param 	set_std				rotate and shift to a standard orientation.
@return  	Vector3<double>			axis.
**/
Vector3<double>		mol_find_alpha_orientation(Bmolecule* mol, int set_std)
{
	if ( verbose & VERB_PROCESS ) {
		cout << "Finding the helical axis:" << endl;
		cout << "Molecule:                       " << mol->id << endl;
	}
	
	int				i, nca;
	Bresidue*		res;
	Batom*			atom;
	double			dca;
	Transform		t;
	Vector3<double>	ca[1000], v1, v2, va, vas, d, ds, axis(0,0,1);
	
	for ( nca=0, res = mol->res; res; res = res->next )
		for ( atom = res->atom; atom; atom = atom->next )
			if ( strncmp(atom->type, "CA", 2) == 0 )
				ca[nca++] = atom->coord;
	
	for ( i=0; i<nca-3; i++ ) {
		v1 = (ca[i] + ca[i+2]) - (ca[i+1] * 2);
		v2 = (ca[i+1] + ca[i+3]) - (ca[i+2] * 2);
		v1.normalize();
		v2.normalize();
		va = v1.cross(v2);
		va.normalize();
		vas += va;
		d = ((va * 1.5) + (ca[i+1] - ca[i+2])) / (v2 - v1);
		ds += d;
		if ( verbose & VERB_FULL )
			cout << va << " " << d << endl;
	}
	
	vas.normalize();
	d = ds / (nca-3.0);
	dca = 0;
	i = 0;
	if ( isfinite(d[0]) && fabs(d[0]) > 0.5 ) {
		dca += d[0];
		i++;
	}
	if ( isfinite(d[1]) && fabs(d[1]) > 0.5 ) {
		dca += d[1];
		i++;
	}
	if ( isfinite(d[2]) && fabs(d[2]) > 0.5 ) {
		dca += d[2];
		i++;
	}
	dca /= i;
	
	i = nca/2;
	v1 = point_on_helix_axis(ca[i-1], ca[i], ca[i+1]);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Helical center:                 " << v1 << endl;
		cout << "Helical axis:                   " << vas << endl;
		cout << "Ca-axis distance:               " << dca << " A" << endl << endl;
	}
	
	if ( set_std ) {
//		i = nca/2;
//		v1 = vector3_subtract(vector3_add(ca[i-1], ca[i+1]), vector3_scale(ca[i], 2));
//		v1 = vector3_normalize(v1);
	
		t.axis = vas.cross(axis);
		t.axis.normalize();
		t.origin = v1;
		t.angle = vas.angle(axis);
	
		t.trans = -t.origin;
		if ( t.angle > 1e-20 )
			mol_coor_rotate(mol, t);
		else
			mol_coor_shift(mol, t.trans);

		for ( i=0, res = mol->res; res && i<nca/2; res = res->next ) i++;
		for ( atom = res->atom; atom && strncmp(atom->type, "CA", 2); atom = atom->next ) ;
	
		t.angle = -atan2(atom->coord[1], atom->coord[0]);
		t.origin = 0;
		t.axis = Vector3<double>(0,0,1);
		t.trans = 0;
	
		mol_coor_rotate(mol, t);
	}
	
	return vas;
}

/**
@brief 	Calculates a point on an alpha helix corresponding to a Ca atom.  

	The point on the alpha helix corresponding to the second Ca atom is returned.

@param 	ca1				first Ca atom.
@param 	ca2				second Ca atom.
@param 	ca3				third Ca atom.
@return  	Vector3<double>					point on alpha helix.
**/
Vector3<double>		point_on_helix_axis(Vector3<double> ca1, Vector3<double> ca2, Vector3<double> ca3)
{
	Vector3<double>		v = (ca1 + ca3) - (ca2 * 2.0);
	
	v.normalize();
	v = (v * DAXIS) + ca2;
	
	return v;
}

long		find_closest_mean(long k, Vector3<double>* mc, Vector3<double>* ma, Vector3<double> c, Vector3<double> a, double afac)
{
	long		i, iclose = -1;
	double		d, dmin = 1e30, angle;
	
	for ( i=0; i<k; i++ ) {
		d = c.distance2(mc[i]);
		angle = afac*a.angle(ma[i]);
		d = sqrt(d + angle*angle);
		if ( dmin > d ) {
			dmin = d;
			iclose = i;
		}
	}
	
	return iclose;
}

/**
@brief 	Cluster a set of alpha helices to consolidate helices.  
@param 	*molgroup		molecule group with alpha helices (deallocated).
@return  	Bmolgroup*				new molecule group.
**/
Bmolgroup*	molgroup_consolidate_alpha(Bmolgroup* molgroup)
{
	long		i, j, x, y, z, nsec;
	double		sampling = 3;
	Bmolecule*	mol, *newmol;
	Bsecondary*	sec;
	
	for ( nsec=0, mol = molgroup->mol; mol; mol = mol->next )
		for ( sec = mol->sec; sec; sec = sec->next ) nsec++;
		
	if ( nsec < 1 ) {
		cerr << "Error: No helices found!" << endl;
		return NULL;
	}
	
	// Convert molecules to position and orientation vectors
	Vector3<double>*	c = new Vector3<double>[nsec];
	Vector3<double>*	a = new Vector3<double>[nsec];
	
	for ( i=0, mol = molgroup->mol; mol; mol = mol->next ) {
		for ( sec = mol->sec; sec; sec = sec->next, i++ ) {
			c[i] = alpha_find_center(sec->first, sec->last);
			a[i] = alpha_find_orientation(sec->first, sec->last);
		}
	}
	
	molgroup->box = molgroup->max - molgroup->min;
	
	// Find the number of density peaks as the number of means and their starting points
	Bimage*		p = new Bimage(Integer, TSimple, 
					(long) (molgroup->box[0]/sampling+1),
					(long) (molgroup->box[1]/sampling+1), 
					(long) (molgroup->box[2]/sampling+1), 1);
	p->sampling(sampling, sampling, sampling);
	
	for ( i=0; i<nsec; i++ ) {
		x = (long) ((c[i][0] - molgroup->min[0])/p->sampling(0)[0]);
		y = (long) ((c[i][1] - molgroup->min[1])/p->sampling(0)[1]);
		z = (long) ((c[i][2] - molgroup->min[2])/p->sampling(0)[2]);
		j = p->index(0, x, y, z, 0);
		p->set(j, (*p)[j] + 1); 
	}

	Bimage* 	pseg = p->regions(1, 0);

	delete p;
	
//	write_img("t.tif", pseg, 0);

	long		k = (long) pseg->maximum();
	int*		n = new int[k];
	Vector3<double>*	mc = new Vector3<double>[k];
	Vector3<double>*	ma = new Vector3<double>[k];
	Vector3<double>*	mcn = new Vector3<double>[k];
	Vector3<double>*	man = new Vector3<double>[k];
	Vector3<double>*	min = new Vector3<double>[k];
	Vector3<double>*	max = new Vector3<double>[k];
	
	for ( i=0; i<k; i++ ) n[k] = 0;
	
	for ( i=z=0; z<pseg->sizeZ(); z++ ) {
		for ( y=0; y<pseg->sizeY(); y++ ) {
			for ( x=0; x<pseg->sizeX(); x++, i++ ) {
				if ( (*pseg)[i] > 0 ) {
					j = (long)(*pseg)[i] - 1;
					n[j]++;
					mc[j] += Vector3<double>(x*pseg->sampling(0)[0] + molgroup->min[0],
						y*pseg->sampling(0)[1] + molgroup->min[1], z*pseg->sampling(0)[2] + molgroup->min[2]);
				}
			}
		}
	}
	
	delete pseg;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Consolidating a set of helices:" << endl;
		cout << "Number of helix groups:         " << k << endl;
	}
	
	for ( j=0; j<k; j++ )
		if ( n[j] ) mc[j] /= n[j];
		
	// Do k-means iterations
	long			it = 0;
	double			afac = 0;
	double			d = 1e30;
	while ( d > 1e-10 ) {
		it++;
		for ( j=0; j<k; j++ ) {
			mcn[j] = man[j] = 0;
			min[j] = Vector3<double>(1e30, 1e30, 1e30);
			max[j] = Vector3<double>(-1e30, -1e30, -1e30);
			n[j] = 0;
		}
		for ( i=0; i<nsec; i++ ) {
			j = find_closest_mean(k, mc, ma, c[i], a[i], afac);
			mcn[j] += c[i];
			man[j] += a[i];
			min[j] = min[j].min(c[i]);
			max[j] = max[j].max(c[i]);
			n[j]++;
		}
		d = 0;
		if ( verbose & VERB_FULL )
			cout << "Iteration " << it << ":" << endl;
		for ( j=0; j<k; j++ ) {
			if ( n[j] ) {
				mcn[j] /= n[j];
				man[j] /= n[j];
			} else {
				mcn[j] = vector3_random(molgroup->max, molgroup->min);
				man[j] = vector3_random(-1.0, 1.0);
				man[j].normalize();
			}
			if ( verbose & VERB_FULL )
				cout << j+1 << tab << n[j] << tab << 
					mcn[j][0] << tab << mcn[j][1] << tab << mcn[j][2] << tab << 
					man[j][0] << tab << man[j][1] << tab << man[j][2] << endl;
			d += mcn[j].distance2(mc[j]);
			d += man[j].distance2(ma[j]);
			mc[j] = mcn[j];
			ma[j] = man[j];
		}
		d = sqrt(d/k);
		afac = 10;
		if ( verbose )
			cout << it << tab << d << endl;
	}
	
	molgroup_kill(molgroup);
	
	// Generate a new set of helices
	molgroup = molgroup_init();
	Bmolecule*		molalpha = NULL;
	View			view;
	Vector3<double>	origin;
	char			molname = 'A';
	
	for ( j=0; j<k; j++ ) if ( n[j] ) {
		d = max[j].distance(min[j]);
		molalpha = mol_generate_alpha_helix(14 + (int)(d/1.5));
		mol_shift_to_center_of_mass(molalpha);
//		view = View(ma[j][0], ma[j][1], ma[j][2], 0);
		view = View(-ma[j][0], -ma[j][1], ma[j][2], 0);
		newmol = mol_rotate_to_view(molalpha, view, origin, mc[j]);
		newmol->id[0] = molname;
		if ( !molgroup->mol ) molgroup->mol = newmol;
		else mol->next = newmol;
		mol = newmol;
		molname++;
		if ( molname > 'Z' ) molname = 'A';
		molecule_kill(molalpha);
	}
	
	delete[] c;
	delete[] a;
	delete[] n;
	delete[] mc;
	delete[] ma;
	delete[] mcn;
	delete[] man;
	delete[] min;
	delete[] max;
	
	return molgroup;
}

