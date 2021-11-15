/**
@file	mol_symmetry.cpp
@brief	Library routines used for symmetry operations on atomic coordinates
@author Bernard Heymann
@date	Created: 20021020
@date	Modified: 20180226
**/

#include "rwmolecule.h"
#include "rwatomprop.h"
#include "rwresprop.h"
#include "mol_compare.h"
#include "mol_symmetry.h"
#include "mol_transform.h"
#include "mol_util.h"
#include "seq_align.h"
#include "Matrix3.h"
#include "matrix_linear.h"
#include "symmetry.h"
#include "random_numbers.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Generates all symmetry-related coordinates.
@param 	*molgroup			molecule group structure.
@param 	*sym				symmetry structure.
@param 	ref_view			reference view.
@return int					0.
**/
int 		molgroup_apply_point_group(Bmolgroup* molgroup, Bsymmetry& sym, View ref_view)
{
	int 			i, j;
	int 			nunits(sym.order()), nmol;
	
	if ( ref_view.vector_size() < 1e-10 ) ref_view = View(0, 0, 1, 0);
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Applying symmetry " << sym.label() << ":" << endl;
		cout << "Number of asymmetric units:     " << nunits << endl;
		cout << "Reference symmetry axis:        {" << 
				ref_view[0] << "," << ref_view[1] << "," << ref_view[2] << "}" << endl;
		cout << "Reference rotation angle:       " <<  
				ref_view.angle()*180/M_PI << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Applying symmetry " << sym.label() << endl << endl;
		
//	if ( verbose & VERB_FULL )
//		sym_show_matrices(sym);
	
	Bmolgroup*		new_molgroup = molgroup;
	Bmolecule		*mol;
	Bresidue		*res;
	Batom			*atom;
	Bbond*			bond;
	Bangle*			angle;
	
	// Initial user-specified rotation to get the map in a standard orientation
//	if ( fabs(ref_view.a) > 1e-30 )
//		molgroup_coor_rotate(molgroup, ref_view, shift);
	
	Matrix3			ref_mat = ref_view.matrix();
	Matrix3			mat(1);
	Vector3<double>	new_axis;

	nmol = 0;
	for ( mol = molgroup->mol; mol; mol = mol->next ) nmol++;

	for ( i=0; i<sym.operations(); i++ ) {
		new_axis = ref_mat * sym[i].axis();
		for ( j=1; j<sym[i].order(); j++ ) {
			mat = Matrix3(new_axis, j*TWOPI*1.0L/sym[i].order());
			if ( verbose & VERB_FULL )
				matrix3_show_hp(mat);
			new_molgroup->next = molgroup_copy(molgroup);
			new_molgroup = new_molgroup->next;
			for ( mol = new_molgroup->mol; mol; mol = mol->next ) {
				for( res = mol->res; res; res = res->next ) {
					for ( atom = res->atom; atom; atom = atom->next ) {
						atom->coord = mat * atom->coord;
					}
				}
			}
		}
		mol = molgroup->mol;
		bond = molgroup->bond;
		angle = molgroup->angle;
		for ( new_molgroup = molgroup->next; new_molgroup; new_molgroup = new_molgroup->next ) {
			if ( mol ) {
				for ( ; mol->next; mol = mol->next ) ;
				mol->next = new_molgroup->mol;
			} else {
				mol = molgroup->mol = new_molgroup->mol;
			}
			new_molgroup->mol = NULL;
			if ( bond ) {
				for ( ; bond->next; bond = bond->next ) ;
				bond->next = new_molgroup->bond;
			} else {
				bond = molgroup->bond = new_molgroup->bond;
			}
			new_molgroup->bond = NULL;
			if ( angle ) {
				for ( ; angle->next; angle = angle->next ) ;
				angle->next = new_molgroup->angle;
			} else {
				angle = molgroup->angle = new_molgroup->angle;
			}
			new_molgroup->angle = NULL;
		}
		molgroup_list_kill(molgroup->next);
		molgroup->next = NULL;
		new_molgroup = molgroup;
		molgroup_atom_renumber(molgroup, 1);
		nmol *= sym[i].order();
	}
	
	return 0;
}

/**
@brief	Generates a helix from the given parameters.
@param 	*molgroup			molecule group structure.
@param 	ref_view			reference view.
@param 	helix_rise			helical rise.
@param 	helix_angle			helical rotation angle.
@param 	gen_down			number of asymmetric units generated upwards.
@param 	gen_up				number of asymmetric units generated downwards.
@return int					0.
**/
int			molgroup_generate_helix(Bmolgroup* molgroup, View ref_view, 
				double helix_rise, double helix_angle, int gen_down, int gen_up)
{
	int 			i;
	
	if ( ref_view.vector_size() < 1e-10 ) ref_view = View(0, 0, 1, 0);
	
	if ( gen_down < 0 ) gen_down = 0;
	if ( gen_up < 0 ) gen_up = 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Applying helical symmetry:" << endl;
		cout << "Number of asymmetric units:     " << gen_down + 1 + gen_up << endl;
		cout << "Helical rise:                   " << helix_rise << endl;
		cout << "Helical rotation angle:         " << helix_angle*180/M_PI << endl;
		cout << "Reference symmetry axis:        {" << 
				ref_view[0] << "," << ref_view[1] << "," << ref_view[2] << "}" << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Applying helical symmetry" << endl << endl;
		
	Bmolgroup*		mg = NULL;
	Bmolecule		*mol;
	Bresidue		*res;
	Batom			*atom;
	
	// Initial user-specified rotation to get the map in a standard orientation
//	if ( fabs(ref_view.a) > 1e-30 )
//		molgroup_coor_rotate(molgroup, ref_view, shift);
	
//	Matrix3			ref_mat = ref_view.matrix();
	Matrix3			mat(1);
	Vector3<double>	axis(ref_view[0], ref_view[1], ref_view[2]);
	Vector3<double>	shift;

	for ( mg = molgroup, i=-gen_down; i<gen_up; i++, mg = mg->next )
		mg->next = molgroup_copy(molgroup);

	for ( mg = molgroup, i=-gen_down; i<=gen_up; i++, mg = mg->next ) {
		shift[2] = i*helix_rise;
		mat = Matrix3(axis, i*helix_angle);
		if ( verbose & VERB_FULL )
			matrix3_show_hp(mat);
		for ( mol = mg->mol; mol; mol = mol->next ) {
			for( res = mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					atom->coord = mat * atom->coord;
					atom->coord += shift;
				}
			}
		}
	}

	return 0;
}


/**
@brief	Generates all symmetry-related coordinates using PDB file matrices.
@param 	*molgroup 		molecule group structure.
@param 	&filename		PDB file with matrices.
@return int 				0.

	Uses the SMTRY records.

**/
int			molgroup_apply_symmetry_from_pdb(Bmolgroup* molgroup, Bstring& filename)
{
	if ( !filename.contains(".pdb") && !filename.contains(".PDB") &&
			!filename.contains(".ent") && !filename.contains(".ENT") )
		return -1;

    ifstream		fpdb(filename.c_str());
    if ( fpdb.fail() ) return -1;
    
	Matrix3			mat[1000];
	Vector3<double>	trans[1000];
    char			aline[100], tag1[12], tag2[12], tag3[12], tag4[12];
	int				n(0), i(0);
	
	while ( fpdb.getline(aline, 100) && fpdb.good() ) {
		if ( strncmp(aline, "REMARK 290   SMTRY", 18) == 0 ) {
			if ( i==0 ) sscanf(aline, "%s %s %s %s %lf %lf %lf %lf",
				tag1, tag2, tag3, tag4, &mat[n][0][0], &mat[n][0][1], &mat[n][0][2], &trans[n][0]);
			else if ( i==1 ) sscanf(aline, "%s %s %s %s %lf %lf %lf %lf",
				tag1, tag2, tag3, tag4, &mat[n][1][0], &mat[n][1][1], &mat[n][1][2], &trans[n][1]);
			else sscanf(aline, "%s %s %s %s %lf %lf %lf %lf",
				tag1, tag2, tag3, tag4, &mat[n][2][0], &mat[n][2][1], &mat[n][2][2], &trans[n][2]);
			i++; 
			if ( strstr(aline, "SMTRY3") ) {
				n++;
				i = 0;
			}
		}
	}
	
    fpdb.close();
	
	Bmolecule		*mol, *new_mol;
	Bresidue		*res;
	Batom			*atom;
	int				nmol(0), m;
	for ( mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( verbose & VERB_PROCESS )
		cout << "Extracting " << n << " matrices from " << filename << " and applying:" << endl;

	for ( i=1; i<n; i++ ) {
		if ( verbose & VERB_PROCESS ) {
			cout << "Matrix " << i+1 << ":" << endl;
			cout << mat[i] << endl;
			cout << "Shift:                          " << trans[i] << endl;
		}
		for ( m=0, mol = molgroup->mol; m<nmol; mol = mol->next, m++ ) {
			new_mol = mol_copy_and_add_to_molgroup(molgroup, mol);
			for( res = new_mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					atom->coord = mat[i] * atom->coord;
					atom->coord += trans[i];
				}
			}
		}
	}

	return 0;
}

/**
@brief	Generates all symmetry-related coordinates using PDB file matrices.
@param 	*molgroup 			molecule group structure.
@param 	&filename			PDB file with matrices.
@return int 					0.

	Uses the BIOMT records.

**/
int			molgroup_apply_matrices_from_pdb(Bmolgroup* molgroup, Bstring& filename)
{
	if ( !filename.contains(".pdb") && !filename.contains(".PDB") &&
			!filename.contains(".ent") && !filename.contains(".ENT") )
		return -1;

    ifstream		fpdb(filename.c_str());
    if ( fpdb.fail() ) return -1;
    
	Matrix3			mat[1000];
	Vector3<double>	trans[1000];
    char			aline[100], tag1[12], tag2[12], tag3[12], tag4[12];
	int				n(0), i(0);
	
	while ( fpdb.getline(aline, 100) && fpdb.good() ) {
		if ( strncmp(aline, "REMARK 350   BIOMT", 18) == 0 ) {
			if ( i==0 ) sscanf(aline, "%s %s %s %s %lf %lf %lf %lf",
				tag1, tag2, tag3, tag4, &mat[n][0][0], &mat[n][0][1], &mat[n][0][2], &trans[n][0]);
			else if ( i==1 ) sscanf(aline, "%s %s %s %s %lf %lf %lf %lf",
				tag1, tag2, tag3, tag4, &mat[n][1][0], &mat[n][1][1], &mat[n][1][2], &trans[n][1]);
			else sscanf(aline, "%s %s %s %s %lf %lf %lf %lf", 
				tag1, tag2, tag3, tag4, &mat[n][2][0], &mat[n][2][1], &mat[n][2][2], &trans[n][2]);
			i++;
			if ( strstr(aline, "BIOMT3") ) {
				n++;
				i = 0;
			}
		}
	}
	
    fpdb.close();

	Bmolecule		*mol, *new_mol;
	Bresidue		*res;
	Batom			*atom;
	int				nmol(0), m;
	for ( mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( verbose & VERB_PROCESS )
		cout << "Extracting " << n << " matrices from " << filename << " and applying:" << endl;

	for ( i=1; i<n; i++ ) {
		if ( verbose & VERB_PROCESS ) {
			cout << "Matrix " << i+1 << ":" << endl;
			cout << mat[i] << endl;
			cout << "Shift:                          " << trans[i] << endl;
		}
		for ( m=0, mol = molgroup->mol; m<nmol; mol = mol->next, m++ ) {
			new_mol = mol_copy_and_add_to_molgroup(molgroup, mol);
			for( res = new_mol->res; res; res = res->next ) {
				for ( atom = res->atom; atom; atom = atom->next ) {
					atom->coord = mat[i] * atom->coord;
					atom->coord += trans[i];
				}
			}
		}
	}

	return 0;
}

struct MassCOM {
	double	mass;
	Vector3<double>	com;
} ;

static int  QsortMassCOM(const void *x, const void *y)
{
	MassCOM*		c1 = (MassCOM *) x;
	MassCOM*		c2 = (MassCOM *) y;
	
	if ( fabs((c1->mass - c2->mass)/(c1->mass + c2->mass)) > 0.1 ) {
		if ( c1->mass < c2->mass ) return 1;
		else return -1;
	} else {
		if ( c1->com.length() < c2->com.length() ) return 1;
		else return -1;		// Sort from high to low
	}
}

/**
@brief	Searches for the standard view based on point group symmetry.
@param 	*molgroup		molecule group structure.
@param 	*sym			point group symmetry.
@param 	ref_view		reference view.
@return int				0.

	The molecule group is first analyzed to identify the different
	chains and calculate their centers-of-mass and weights.
	The overall center-of-mass defines a point on at least the
	major symmetry axis (cyclic symmetries), or the likely intersection
	of symmetry axes.
	Note: This function does a reasonable job of orienting the molecule
	group, but it may be off by up to an angstrom!!!

**/
int 		molgroup_find_standard_view(Bmolgroup* molgroup, Bsymmetry& sym, View ref_view)
{
	molgroup_shift_to_center_of_mass(molgroup);
	
	int				g, i, j, k, nmol, nunits(sym.order()), nmassif;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	for ( nmol = 0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( nmol < nunits ) {
		cerr << "Error: Not enough molecules (" << nmol << ") to fit into " << nunits << " asymmetric units!" << endl << endl;
		return -1;
	}

	// Get all chains' centers of mass and weight
	double			themass, themass2, maxmass(0);
	MassCOM*		c = new MassCOM[nmol];
	for ( i=0; i<nmol; i++ ) c[i].com = c[i].mass = 0;
	
	nmassif = 0;
    for ( i=0, mol = molgroup->mol; mol; mol = mol->next, i++ ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				c[i].com += atom->coord * atom->mass;
				c[i].mass += atom->mass;
			}
	    }
		if ( c[i].mass ) {
			c[i].com /= c[i].mass;
			nmassif++;
		}
		if ( maxmass < c[i].mass ) maxmass = c[i].mass;
	}
	
	// Sort the chains based on weight and distance from COM
	qsort((void *) c, nmol, sizeof(MassCOM), (int (*)(const void *, const void *)) QsortMassCOM);
	
	if ( verbose & VERB_FULL ) {
		cout << endl << "Ordered list of centers-of-mass:" << endl;
		cout << "    x\t    y\t    z\tMass (Da)" << endl;
		for ( i=0; i<nmol; i++ )
			cout << c[i].com[0] << tab << c[i].com[1] << tab << c[i].com[2] << tab << c[i].mass << endl;
		cout << endl << endl;
	}
	
	// Determine the number of groups of chains, each group likely symmetry-related
	// A group is defined as a set of chains (number = symmetry order), 
	// that differ little in weight (stdev < 5%)
	int				ngroup(0);
	Vector3<double>	com, translate;
	
	if ( verbose &  VERB_FULL )
		cout << endl << "Group\tMass\tStDev" << endl;
	for ( i=0; i<nmol-nunits+1; i+=nunits ) {   // Loop over all groups of size defined by symmetry
		themass = themass2 = 0;
		com = 0;
		for ( j=i; j<i+nunits; j++ ) {
			themass += c[j].mass;
			themass2 += c[j].mass*c[j].mass;
			com += c[j].com;	// Center of group on symmetry axis
		}
		themass /= nunits;
		themass2 = themass2/nunits - themass*themass;
		if ( themass2 > 0 ) themass2 = sqrt(themass2);
		else themass2 = 0;
		if ( themass > 100 && themass2/themass < 0.05 ) {
			ngroup++;
			com *= -1.0/nunits;	// Translation to symmetry axis
			translate += com;
		}
		if ( verbose &  VERB_FULL )
			cout << ngroup << tab << themass << tab << themass2 << endl;
	}
	translate /= ngroup;
	
	// Subtract the global center-of-mass from the chain COM's
	com = 0;
	for ( i=0; i<nmassif; i++ ) com += c[i].com;
	com /= nmassif;
	for ( i=0; i<nmassif; i++ ) c[i].com -= com;
	
	// Generate pairs of vectors and calculate normals
	int				n(0), nvec = (ngroup*nunits*(nunits-1)*(nunits-2))/6;
	double			a;
	Vector3<double>	axis, v1, v2, d;
	Vector3<double>*	v = new Vector3<double>[nvec];
	
	if ( verbose &  VERB_FULL )
		cout << endl << "Group\tv1\tv2\tv3\tx\ty\tz" << endl;
	for ( g=0; g<ngroup; g++ ) {
		// Calculate normal vectors
		for ( i=2+nunits*g; i<nunits*(g+1); i++ ) {
			for ( j=1+nunits*g; j<i; j++ ) {
				v1 = c[i].com - c[j].com;
				for ( k=nunits*g; k<j; k++ ) {
					v2 = c[i].com - c[k].com;
//					v[n] = v1.cross(v2) * c[i].mass;
					v[n] = v1.cross(v2);
					v[n].normalize();
					if ( verbose &  VERB_FULL )
						cout << g << tab << k << tab << j << tab << i << tab 
							<< v[n][0] << tab << v[n][1] << tab << v[n][2] << endl;
					n++;
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_find_standard_view: nvec=" << nvec << " n=" << n << endl;
		
	// Find a common normal when comparing pairs of vectors and assign it as the major axis
	int				vbest = -1, m = 100;
	int				nx = 2*m, ny = 2*m, nz = 2;
	int				sv = nx * ny * nz;
	int*			nv = new int[sv];
	for ( i=0; i<sv; i++ ) nv[i] = 0;
	
	for ( i=1; i<n; i++ ) {
		for ( j=0; j<i; j++ ) {
			v1 = v[i];
			v1.normalize();
			a = fabs(v[i].angle(v[j]));
			if ( a > M_PI_2 ) {
				a = M_PI - a;
				v1 = -v1;
			}
			v1 *= m;
			if ( a < 0.05 ) {
				k = (int) (nx*(v1[1]+m)) + (int) (v1[0]+m);
				if ( k > nx*ny ) k = 0;
				if ( v1[2] > 0 ) k += nx*ny;
				nv[k]++;
				if ( vbest < nv[k] ) {
					vbest = nv[k];
					axis = v1;
				}
			}
		}
	}
	delete[] nv;
	
	for ( i=0; i<n; i++ ) {
		a = fabs(axis.angle(v[i]));
		if ( a > M_PI_2 ) {
			a = M_PI - a;
			v[i] = -v[i];
		}
		if ( a < 0.05 ) axis += v[i];
	}
	
	axis.normalize();
	
	// Find the 2-fold on the x-axis as the sum vector of two COM's
	// that has a difference vector that is most parallel with the major axis
	double		ang(0), ang_min = 1000;
	if ( sym.point() > 200 ) {
		for ( i=1; i<nunits; i++ ) {
			for ( j=0; j<i; j++ ) {
				d = c[i].com - c[j].com;
				ang = axis.angle(d);
				if ( fabs(ang_min) > fabs(ang) ) {
					ang_min = ang;
					v2 = c[i].com + c[j].com;
					v2.normalize();
				}
			}
		}
	}
	
	delete[] c;
	delete[] v;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Angle with major axis: " << ang_min*180/M_PI << endl;
		cout << "Dihedral axis:         " << v2 << endl;
		cout << "Angle between axes:    " << (axis.angle(v2))*180/M_PI << endl;
	}
	
	if ( verbose & VERB_LABEL ) {
		cout << "Translate:                      " << translate << endl;
		cout << "Major symmetry axis:            " << axis << endl;
		if ( sym.point() > 200 )
			cout << "Minor symmetry axis:            " << v2 << endl;
		cout << endl;
	}
	
	// The first matrix rotates it to the major symmetry axis
	// The second does the in-plane rotation
	View			view(-axis[0], -axis[1], axis[2], 0);
	Matrix3			mat = view.matrix();
	if ( sym.point() > 200 ) {
		v1 = 0;
		v1[0] = 1;	
		v2 = mat * v2;
		ang = fmod(v1.angle(v2), M_PI*2.0/sym[0].order());
		view= View(0, 0, 1, ang);
		mat = view.matrix() * mat;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Dihedral axis:                      " << v2 << endl;
		cout << "Angle: " << ang*180.0/M_PI << endl;
	}
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				atom->coord += translate;
				atom->coord = mat * atom->coord;
			}
		}
	}
	
	return 0;
}

/**
@brief	Searches for the standard view based on point group symmetry.
@param 	*molgroup 	molecule group structure.
@param 	*sym		point group symmetry.
@param 	ref_view	reference view (default should be 0,0,1,0).
@param 	*simat		residue similarity matrix.
@return int 		0.

	Each pair of chains in the molecule groupis tested for sequence
	identity to find symmetry-related molecules. For each pair of matched
	molecules, the transformation to superimpose the one onto the other
	is determined and the symmetry axis and translation calculated.
	The collection of symmetry axes are clustered with a k-means 
	algorithm and the predominant class assigned to the major 
	symmetry axis. For dihedral point groups, a minor axis is also
	assigned (randomly at this time). The molecule group is then
	transformed to orient it with the major axis on {0,0,1} and
	the minor axis on {1,0,0}, and the symmetry center at {0,0,0}.
	Note: This function has not been extensively tested with all
	point groups!!!

**/
int 		molgroup_orient_to_standard_view(Bmolgroup* molgroup, Bsymmetry& sym,
				View ref_view, Bresidue_matrix* simat)
{
	random_seed();
	
	if ( sym.point() < 102 ) return 0;
	
	long			i, j, k, m, nmol, nres_cut, cut_max(30);
	long			nunit(sym.order()), ngroup, nid, n2=0, nn=0, no;
	Vector3<double>	a2[1000], an[1000];
	Bmolecule*		mol, *mol2, *symmol;
	Bresidue*		res;
	Batom*  		atom;
	
	for ( i=m=0; i<sym.operations(); i++ )
		if ( m < sym[i].order() ) m = sym[i].order();
	
	double			a, refang = TWOPI/m;
	
	for ( nmol = 0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	ngroup = nmol/nunit;
	
	if ( ngroup < 1 ) {
		cerr << "Error: Not enough molecules (" << nmol << ") to fit into " << nunit << " asymmetric units!" << endl << endl;
		return -1;
	}

	molgroup_shift_to_center_of_mass(molgroup);

	if ( verbose & VERB_PROCESS ) {
		cout << "Orienting to the standard view for symmetry " << sym.label() << endl;
		cout << "Number of molecules:            " << nmol << endl;
		cout << "Number of groups:               " << ngroup << endl;
	}
	
	int*			g = new int[nmol];	// Group indices
	for ( i=0; i<nmol; i++ ) g[i] = 0;
	
	Transform		t;
	int				offset, notfound;
	
	// Calculate a matrix of rotations between pairs of molecules 
	for ( k=0, i=1; k<nmol; k++, i++, i=(i>ngroup)?1:i ) g[k] = i;
	for ( i=0, mol = molgroup->mol; mol->next; mol = mol->next, i++ ) {
		nres_cut = mol->nres/2;
		if ( nres_cut > cut_max ) for ( j=i+1, mol2 = mol->next; mol2; mol2 = mol2->next, j++ ) {
			// if two sequences can be aligned, find the transformation to superimpose them
			offset = seq_find_best_offset(mol, mol2, nid, simat);
			if ( nid > cut_max || nid >= nres_cut ) {
				g[j] = g[i];
				t = mol_find_transformation(mol, mol2, offset);
				symmol = molecule_copy(mol);
				mol_coor_rotate(symmol, t);
				t.fom = mol_calculate_rmsd(mol2, symmol);
				molecule_kill(symmol);
				angle_set_negPI_to_PI(t.angle);
				a = t.angle/refang;
				no = (long) (a + 0.5);
				if ( fabs(t.angle - M_PI) < 0.1 ) {
					for ( k=0, notfound=1; k<n2 && notfound; k++ ) {
						if ( a2[k].angle(t.axis) < 0.1 ) {
							a2[k] += t.axis;
							notfound = 0;
						}
					}
					if ( notfound ) {
						a2[n2] = t.axis;
						n2++;
					}
				}
				if ( fabs(a - no) < 0.03 ) {
					for ( k=0, notfound=1; k<nn && notfound; k++ ) {
						if ( an[k].angle(t.axis) < 0.1 ) {
							an[k] += t.axis;
							notfound = 0;
						}
					}
					if ( notfound ) {
						an[nn] = t.axis;
						nn++;
					}
				}
				m++;
			}
		}
	}

	if ( verbose & VERB_PROCESS ) {
		cout << "Num\tMol\tGroup\tLength" << endl;
		for ( i=0, mol = molgroup->mol; mol; mol = mol->next, i++ )
			cout << i+1 << tab << mol->id << tab << g[i] << tab << mol->nres << endl;
		cout << endl;
	}
	
	delete[] g;
			
	double			d, dm;
	Vector3<double>	axis1, axis2;
	
	if ( sym.point() == 202 || sym.point() == 320 || sym.point() == 532 ) {
		// Pick the longest 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			d = a2[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = a2[k];
			}
		}
		// Pick the longest perpendicular 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			a = axis1.angle(a2[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = a2[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = a2[k];
				}
			}
		}
	} else if ( sym.point() == 432 ) {
		// Pick the longest 4-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			d = an[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = an[k];
			}
		}
		// Pick the longest perpendicular 4-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			a = axis1.angle(an[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = an[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = an[k];
				}
			}
		}
	} else {
		// Pick the longest n-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			d = an[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = an[k];
			}
		}
		// Pick the longest perpendicular 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			a = axis1.angle(a2[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = a2[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = a2[k];
				}
			}
		}
	}
	
	if ( axis1.length() < 0.9 ) {
		cerr << "Error: Major axis not found!" << endl;
	}
	
	if ( sym.point() > 200 && axis2.length() < 0.9 ) {
		cerr << "Error: Minor axis not found!" << endl;
	}
	
	axis1.normalize();
	axis2.normalize();
	
	
	// Translation required is the negative of the COM
//	Vector3<double>	translate = -(molgroup_center_of_mass(molgroup));
	
	if ( verbose & VERB_LABEL ) {
		if ( sym.point() < 200 ) {
			cout << "Symmetry axis:                  " << axis1 << endl;
		} else {
			cout << "Major symmetry axis:            " << axis1 << endl;
			cout << "Minor symmetry axis:            " << axis2 << endl;
		}
	}
	
	// The first matrix rotates it to the major symmetry axis
	// The second does the in-plane rotation
	Vector3<double>	v1(0,0,1), v2(1,0,0);
	Matrix3         mat = Matrix3(axis1, v1);
	if ( sym.point() > 200 ) {
		v1 = mat * axis2;
		mat = Matrix3(v1, v2) * mat;
	}
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
//				atom->coord += translate;
				atom->coord = mat * atom->coord;
			}
		}
	}
	
	return 0;
}

/**
@brief	Calculates the RMSD for symmetry in the standard orientation.
@param 	*molgroup		molecule group structure.
@param 	*sym			point group symmetry.
@return double			RMSD.

	The centers-of-mass for all the molecules are calculated as a reduced
	representation of the molecule group. All symmetry operations are 
	imposed on the centers, with the RMSD defined as minimum distance
	between an original cneter and a transformed center.

**/
double		molgroup_symmetry_RMSD(Bmolgroup* molgroup, Bsymmetry& sym)
{
	int				i, j, k, l, n, nmol;
	double			d, mind, R(0);
	Vector3<double>	coms;
	Matrix3			mat(1);
	
	Bmolecule*	mol;
	
	for ( nmol=0, mol=molgroup->mol; mol; mol=mol->next ) nmol++;
	
	Vector3<double>*	com = new Vector3<double>[nmol];
	
	for ( i=0, mol=molgroup->mol; mol; mol=mol->next, i++ )
		com[i] = mol_center_of_mass(mol);
	
	for ( i=n=0; i<sym.operations(); i++ ) {
		for ( j=1; j<sym[i].order(); j++ ) {
			mat = Matrix3(sym[i].axis(), j*TWOPI/sym[i].order());
			for ( k=0; k<nmol; k++ ) {
				coms = mat * com[k];
				for ( mind=1e30, l=0; l<nmol; l++ ) {
					d = coms.distance(com[l]);
					if ( mind > d ) mind = d;
				}
				R += mind*mind;
				n++;
			}
		}
	}
	
	R = sqrt(R/n);
	
	if ( verbose ) 
		cout << "Symmetry RMSD:                  " << R << " A (" << n << ")" << endl << endl;
	
	delete[] com;

	return R;
}

/**
@brief	Calculates the B factors from symmetry-related molecules.
@param 	*molgroup		molecule group structure.
@param 	*sym			point group symmetry.
@return double			RMSD.

	The centers-of-mass for all the molecules are calculated as a reduced
	representation of the molecule group. All symmetry operations are 
	imposed on the centers, with the RMSD defined as minimum distance
	between an original cneter and a transformed center.

**/
double		molgroup_symmetry_B(Bmolgroup* molgroup, Bsymmetry& sym)
{
	long			i, j, nmol(0), nsym(0);
	double			d, mind, R(0);
	Vector3<double>	coms, coor;
	Matrix3			m(1);
	
	Bmolecule*		mol, *molr;
	Bresidue*		res, *resr;
	Batom*  		atom, *atomr;
	
	for ( nmol=0, mol=molgroup->mol; mol; mol=mol->next ) nmol++;
	
	if ( verbose )
		cout << "Calculating B-factors for symmetry " << sym.label() << endl;
	
	Vector3<double>*	com = new Vector3<double>[nmol];
//	Matrix3*			mat = symmetry_get_all_matrices(sym, nsym);
	vector<Matrix3>		mat = sym.matrices();
	nsym = mat.size();
	
	for ( i=0, mol=molgroup->mol; mol; mol=mol->next, i++ ) {
		com[i] = mol_center_of_mass(mol);
		for ( mind=1e30, j=0; j<nsym; j++ ) {
			coms = mat[j] * com[i];
			d = com[0].distance(coms);
			if ( mind > d ) {
				mind = d;
				mol->sel = j;
			}
		}
		cout << mol->id << endl << mat[mol->sel] << endl;
	}
	
	// Initialize the vectors for the first molecule
	molr = molgroup->mol;
	for( res = molr->res; res; res = res->next )
		for ( atom = res->atom; atom; atom = atom->next ) {
			atom->b = 0;
			atom->vel = Vector3<double>(0,0,0);
		}
	
	// Accumulate the sums and square sums in the first molecule atoms
	for ( i=0, mol=molgroup->mol; mol; mol=mol->next, i++ ) {
		for( res = mol->res, resr = molr->res; res; res = res->next, resr = resr->next ) {
			for ( atom = res->atom, atomr = resr->atom; atom; atom = atom->next, atomr = atomr->next ) {
//				m = mat[i].transpose();
				coor = mat[mol->sel] * atom->coord;
				atomr->vel += coor;
				atomr->b += coor.length2();
			}
		}
	}
	
	// Calculate the B factors in the first molecule
	for( res = molr->res; res; res = res->next )
		for ( atom = res->atom; atom; atom = atom->next )
			atom->b = (atom->b - atom->vel.length2()/nmol)/nmol;
	
	// Distribute the B factors to all molecules
	for ( mol=molgroup->mol->next; mol; mol=mol->next )
		for( res = mol->res, resr = molr->res; res; res = res->next, resr = resr->next )
			for ( atom = res->atom, atomr = resr->atom; atom; atom = atom->next, atomr = atomr->next )
				atom->b = atomr->b;
	
	delete[] com;
//	delete[] mat;

	return R;
}

/**
@brief 	Generates unit cells from a set of coordinates.
@param 	*molgroup	molecule group.
@param 	uc			unit cell dimensions.
@param 	number		number of unit cells in each lattice direction.
@return int 		0, <0 if error.

	The input molecule group is duplicated to generate the requested number
	of copies in each lattice direction.

**/
int 		molgroup_generate_crystal(Bmolgroup* molgroup, UnitCell uc, Vector3<int> number)
{
	if ( number.volume() < 2 ) return 0;
	
	if ( !uc.check() ) {
		cerr << "Error: Please specify the unit cell!" << endl;
		return -1;
	}
	
	int				i, x, y, z, nmol(0);
	Vector3<double>	d;
	Matrix3			mat(uc.a(), uc.b()*cos(uc.gamma()), uc.c()*cos(uc.beta()),
						0, uc.b()*sin(uc.gamma()), uc.c()*(cos(uc.alpha()) - cos(uc.beta())*cos(uc.gamma()))/sin(uc.gamma()),
						0, 0, uc.volume()/(uc.a()*uc.b()*sin(uc.gamma())));
	Bmolecule		*mol, *newmol;
	Bresidue		*res, *newres;
	Batom			*atom, *newatom;
	
	for ( nmol=0, mol = molgroup->mol; mol; mol = mol->next ) nmol++;
	
	if ( verbose ) {
		cout << "Generating new unit cells for " << nmol << " molecules:" << endl;
		cout << "Number:                         " << number << " = " << (int)number.volume() << endl;
	}
	
	for ( z=0; z<number[2]; z++ ) {
		for ( y=0; y<number[1]; y++ ) {
			for ( x=0; x<number[0]; x++ ) {
				if ( verbose )
					cout << "Generating unit cell:           " << x << " " << y << " " << z << endl;
				d = Vector3<double>(x,y,z);
				d = mat * d;
				if ( x+y+z > 0 ) for ( i=0, mol = molgroup->mol; i<nmol; mol = mol->next, i++ ) { 
					newmol = mol_copy_and_add_to_molgroup(molgroup, mol);
					for ( res = mol->res, newres = newmol->res; res; res = res->next, newres = newres->next ) {
						for ( atom = res->atom, newatom = newres->atom; atom; atom = atom->next, newatom = newatom->next ) {
							newatom->coord = atom->coord + d;
						}
					}
				}
			}
		}
	}
	
	return 0;
}

