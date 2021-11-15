/**
@file	molecule_to_map.cpp
@brief	Functions to calculate a 3D map from atomic coordinates 
@author Bernard Heymann
@date	Created: 19970914
@date	Modified: 20210216
**/
 
#include "molecule_to_map.h"
#include "scatter.h"
#include "Complex.h"
#include "UnitCell.h"
#include "linked_list.h"
#include "utilities.h"
#include <fstream>
 
// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Constants
#define MAXSID  	1024		// Maximum grid sidelength
#define MAXSCAT 	100			// Maximum number of atomic scattering data points
#define MAXRAD 		100			// Maximum radial points at 0.1 A/point for calculating atomic potential curves

// Internal function prototypes
int 		mol_to_image(Bmolgroup* molgroup, Bimage* p, Batomtype* atompar, 
				double resolution, double Bfactor, int wrap, int gextype);
int			mol_to_structure_factors(Bmolgroup* molgroup, Bimage* p, 
				Batomtype* atompar, double resolution, int wrap, double Bfactor);

/**
@brief 	Calculates a 3D density map from a set of atomic coordinates.
@param 	*molgroup	set of molecules with atomic coordinates.
@param 	origin		3-valued origin vector (angstrom).
@param 	size		3-valued size vector (voxels).
@param 	sampling	sampling/voxel size (angstrom/voxel).
@param 	resolution	resolution (angstrom).
@param 	Bfactor		global B-factor to use, if 0, use individual atom B-factors
@param 	wrap		wrapping flag.
@param 	gextype		type of gaussian used: 0 = single, 1 = atomic potential
@param 	spacegroup	crystal space group.
@param 	unit_cell	6-valued vector of unit cell parameters.
@return Bimage*		the new map.

	A 3D map is calculated from atomic coordinates by placing a gaussian
	sphere at each set of atomic coordinates. The resolution is set as
	twice the sigma coefficient of the gaussian function. The amplitude
	of the gaussian function is set so that the total density calculated
	equals the atomic mass. The resultant map therefore has the density
	units of Dalton/voxel.
	The statistics of the new image is calculated.

**/
Bimage* 	img_from_molecule(Bmolgroup* molgroup, Vector3<double> origin, 
				Vector3<long> size, Vector3<double> sampling, double resolution, double Bfactor,
				int wrap, int gextype, int spacegroup, UnitCell unit_cell)
{
	Batomtype*	atompar = NULL;
	Bstring		paramfile;
	if ( gextype) atompar = get_atom_properties(paramfile);

	if ( sampling.volume() < 1e-6 )
		sampling[0] = sampling[1] = sampling[2] = 1;
	else
		sampling = sampling.max(0.01);
	
	if ( resolution < sampling[0] ) resolution = sampling[0];
	if ( resolution < 0.01 ) resolution = 2*sampling[0];
	
	if ( size[0] < 1 ) size[0] = (long) ((molgroup->max[0] + 3*resolution)/sampling[0] + origin[0]);
	if ( size[1] < 1 ) size[1] = (long) ((molgroup->max[1] + 3*resolution)/sampling[1] + origin[1]);
	if ( size[2] < 1 ) size[2] = (long) ((molgroup->max[2] + 3*resolution)/sampling[2] + origin[2]);
	size = size.min(MAXSID);
	size = size.max(1);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_sf_from_molecule: size = " << size << endl;
		
	Bimage* 	p = new Bimage(Float, TSimple, size, 1);
	p->unit_cell(unit_cell);
	p->origin(origin);
	p->sampling(sampling);
	p->space_group(spacegroup);
	
	int 		side = (long) floor((3.0/2.0)*resolution/p->sampling(0)[0]+1);
    if ( side > MAXSID ) side = MAXSID;
	
	if ( verbose & VERB_PROCESS ) {
		if ( gextype ) cout << "Calculating atomic potential:" << endl;
		else cout << "Calculating density by simple gaussian expansion:" << endl;
	    cout << "Dimensions:                     " << p->size() << " voxels" << endl;
    	cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
    	if ( !gextype )
			cout << "Resolution:                     " << resolution << " A" << endl;
		cout << "Global B-factor:                " << Bfactor << " A2" << endl;
	    cout << "Wrap:                           " << wrap << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Calculating density" << endl << endl;
	
	mol_to_image(molgroup, p, atompar, resolution, Bfactor, wrap, gextype);
	
	if ( gextype) kill_list((char *) atompar, sizeof(Batomtype));
	
	p->statistics();
	
	return p;
}

/*
@brief 	Calculates a 3D density map from a set of atomic coordinates.
@param 	*molgroup 	set of molecules with atomic coordinates.
@param 	*p			image structure for the new map.
@param	*atompar	atom type parameters.
@param 	resolution	resolution (angstrom).
@param 	Bfactor		global B-factor to use, if 0, use individual atom B-factors
@param 	wrap		wrapping flag.
@param 	gextype		type of gaussian used: 0 = single, 1 = atomic potential
@return int			0.

	A 3D map is calculated from atomic coordinates by placing a gaussian
	sphere at each set of atomic coordinates. The resolution is set as
	twice the sigma coefficient of the gaussian function. The amplitude
	of the gaussian function is set so that the total density calculated
	equals the atomic mass. The resultant map therefore has the density
	units of Dalton/voxel.
	The content and statistics of the new image is not checked.

**/
int 		mol_to_image(Bmolgroup* molgroup, Bimage* p, Batomtype* atompar, 
				double resolution, double Bfactor, int wrap, int gextype)
{
	long   			h, i, j, irad, nsig, t;
	long			x, y, z, ix, iy, iz;
	Vector3<int> 	lo, hi;
	Vector3<double>	cgrid;
	double			sigma, sigma2, avgsig, dx, dy, dz, dist2;
	double			sigma2_global = resolution*resolution/4 + Bfactor/(8*M_PI*M_PI);
	double			interval = p->sizeX()*p->sampling(0)[0]*0.5/MAXRAD;		// Sampling for radial potential curve calculations
	if ( interval > 0.1 ) interval = 0.1;
	double			sam2 = p->sampling(0)[0]*p->sampling(0)[0];
	double			voxel_size = p->sampling(0).volume();
	double			thisamp, thisdens, thisatom, summass(0), totmass(0);
	long   			side = (long) floor((3.0/2.0)*resolution/p->sampling(0)[0]+1);
    if ( side > MAXSID ) side = MAXSID;
	if ( gextype ) {
		side = (long) (MAXRAD*interval/p->sampling(0)[0] - 1);
		if ( side*p->sampling(0)[0] > 5 ) side = (long) (5/p->sampling(0)[0] + 1);
	}
	if ( side > p->sizeX() ) side = p->sizeX();
	double			radius2 = side*side*sam2;
	
	if ( verbose & VERB_PROCESS ) {
    	cout << "Sigma:                          " << sqrt(sigma2_global) << " A" << endl;
    	cout << "Cube edge length:               " << 2*side+1 << " voxels" << endl;
	}
	
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	double*			apot = NULL;
	if ( atompar ) apot = get_potential_curves(atompar, interval);
	
	nsig = 0;
	avgsig = 0;
    for ( h=0, mol = molgroup->mol; mol; mol = mol->next ) if ( mol->sel ) {
		h++;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				cgrid[0] = atom->coord[0]/p->sampling(0)[0] + p->image->origin()[0];
				cgrid[1] = atom->coord[1]/p->sampling(0)[1] + p->image->origin()[1];
				cgrid[2] = atom->coord[2]/p->sampling(0)[2] + p->image->origin()[2];
				lo[0] = (long) (cgrid[0] + 0.5 - side);
				hi[0] = (long) (cgrid[0] + 0.5 + side);
				lo[1] = (long) (cgrid[1] + 0.5 - side);
				hi[1] = (long) (cgrid[1] + 0.5 + side);
				lo[2] = (long) (cgrid[2] + 0.5 - side);
				hi[2] = (long) (cgrid[2] + 0.5 + side);
				if ( !wrap ) {
					if ( lo[0] < 0 ) lo[0] = 0;
					if ( lo[1] < 0 ) lo[1] = 0;
					if ( lo[2] < 0 ) lo[2] = 0;
					if ( hi[0] > (long)p->sizeX() - 1 ) hi[0] = p->sizeX() - 1;
					if ( hi[1] > (long)p->sizeY() - 1 ) hi[1] = p->sizeY() - 1;
					if ( hi[2] > (long)p->sizeZ() - 1 ) hi[2] = p->sizeZ() - 1;
				}
				if ( Bfactor > 0 ) sigma2 = sigma2_global;
				else sigma2 = sigma2_global + atom->b/(8*M_PI*M_PI);	// B-factors
//				if ( sigma2 < 0.01 ) sigma2 = 0.01;
				sigma = sqrt(sigma2);
				nsig++;
				avgsig += sigma;
				if ( verbose & VERB_DEBUG )
					cout << atom->coord << " " << atom->type << " " << atom->mass << endl;
				thisatom = 0;
				thisamp = atom->mass*voxel_size/(sqrt(TWOPI*sigma2)*TWOPI*sigma2);
				sigma2 = -1.0/(2.0*sigma2);
				ix = (long) (cgrid[0] + 0.5);
				iy = (long) (cgrid[1] + 0.5);
				iz = (long) (cgrid[2] + 0.5);
				i = p->index(0, ix, iy, iz, 0);
				for ( z=lo[2]; z<=hi[2]; z++ ) {
					dz = (cgrid[2] - z)*(cgrid[2] - z)*sam2;
					for ( y=lo[1]; y<=hi[1]; y++ ) {
						dy = (cgrid[1] - y)*(cgrid[1] - y)*sam2;
						for ( x=lo[0]; x<=hi[0]; x++ ) {
							dx = (cgrid[0] - x)*(cgrid[0] - x)*sam2;
							dist2 = dx + dy + dz;
							if ( dist2 <= radius2 ) {
								ix = x;
								iy = y;
								iz = z;
								if ( wrap ) {
									while ( ix < 0 ) ix += p->sizeX();    /* This wraps the density */
									while ( iy < 0 ) iy += p->sizeY();
									while ( iz < 0 ) iz += p->sizeZ();
									while ( ix >= p->sizeX() ) ix -= p->sizeX();
									while ( iy >= p->sizeY() ) iy -= p->sizeY();
									while ( iz >= p->sizeZ() ) iz -= p->sizeZ();
								}
								j = p->index(ix, iy, iz);
								thisdens = 0;
								if ( gextype ) {
									t = atom->tnum;
									if ( j == i ) {
										thisdens = apot[t*MAXRAD];
									} else {
										irad = (long) (sqrt(dist2)/interval + 0.5);
										if ( irad < MAXRAD ) thisdens = apot[t*MAXRAD+irad];
									}
								} else {
									thisdens = thisamp*exp(dist2*sigma2);
								}
								p->add(j, thisdens);
								thisatom += thisdens;
							}
						}
					}
				}
				totmass += thisatom;
				summass += atom->mass;
			}
	    }
		if ( verbose & VERB_TIME )
			cout << "Molecules done:                 " << h+1 << "\r" << flush;
	}
	
	double		scale = POTPREFAC;
	if ( p->sizeX() > 1 ) scale *= p->sampling(0)[0];
	if ( p->sizeY() > 1 ) scale *= p->sampling(0)[1];
	if ( p->sizeZ() > 1 ) scale *= p->sampling(0)[2];
	
	if ( gextype ) {
		p->rescale(scale, 0);
		totmass *= scale;
	}

	if ( verbose & VERB_PROCESS ) {
	    cout << "Average sigma:                  " << avgsig/nsig << " A" << endl;
	    cout << "Total mass (sum):               " << summass << " Da" << endl;
    	cout << "Total mass (density):           " << totmass << " Da" << endl;
    	cout << "               Ratio:           " << totmass/summass << endl;
    	cout << "Total volume:                   " << totmass/RHO << " A3" << endl << endl;
	}
	
	if  ( apot ) delete[] apot;
	
	return 0;
}

/**
@brief 	Compares reference and calculated maps and calculates an occupancy
	for every atom in the molecule set.
@param 	*molgroup 	set of molecules with atomic coordinates.
@param 	*pcalc		map calculated from the set of molecules.
@param 	*p			reference map.
@return int			0.
**/
int 		compare_mol_map(Bmolgroup* molgroup, Bimage* pcalc, Bimage* p)
{
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	
	long   			j, natom(0);
	long			x, y, z;
	double			min = 1e30, max = -1e30, Qavg(0), Qsd(0);
	double			scale = (pcalc->maximum() - pcalc->minimum())/(p->maximum() - p->minimum());
	
    for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) {
				x = (long) floor(atom->coord[0]/pcalc->sampling(0)[0] + pcalc->image->origin()[0] + 0.5);
				y = (long) floor(atom->coord[1]/pcalc->sampling(0)[1] + pcalc->image->origin()[1] + 0.5);
				z = (long) floor(atom->coord[2]/pcalc->sampling(0)[2] + pcalc->image->origin()[2] + 0.5);
				while ( x < 0 ) x += pcalc->sizeX();
				while ( x >= pcalc->sizeX() ) x -= pcalc->sizeX();
				while ( y < 0 ) y += pcalc->sizeY();
				while ( y >= pcalc->sizeY() ) y -= pcalc->sizeY();
				while ( z < 0 ) z += pcalc->sizeZ();
				while ( z >= pcalc->sizeZ() ) z -= pcalc->sizeZ();
				j = pcalc->index(x, y, z);
				if ( (*pcalc)[j] ) atom->q = scale*(*p)[j]/(*pcalc)[j];
				if ( atom->q > 1 ) atom->q = 1;
				if ( min > atom->q ) min = atom->q;
				if ( max < atom->q ) max = atom->q;
				Qavg += atom->q;
				Qsd += atom->q*atom->q;
				atom->b = atom->q;
				natom++;
			}
		}
	}
	Qavg /= natom;
	Qsd = sqrt(Qsd/natom - Qavg*Qavg);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Comparison of molecule " << mol->id << " with map " << p->file_name() << endl;
    	cout << "Origin (original):              " << p->image->origin() << endl;
    	cout << "Origin (calculated):            " << pcalc->image->origin() << endl;
		cout << "Scale:                          " << scale << endl;
		cout << "Min & max:                      " << min << " " << max << endl;
		cout << "Average and standard deviation: " << Qavg << " " << Qsd << endl << endl;
	} else {
		cout << "Comparison of molecule " << mol->id << " with map " << p->file_name() << endl;
		cout << "Average and standard deviation: " << Qavg << " " << Qsd << endl << endl;
	}
		
	return 0;
}



/**
@brief 	Calculates a 3D set of structure factors from a set of atomic coordinates.
@param 	*molgroup	set of molecules with atomic coordinates.
@param 	origin		3-valued origin vector (voxels).
@param 	size		3-valued size vector (voxels).
@param 	sampling	sampling/voxel size (angstrom/voxel).
@param 	resolution	resolution (angstrom).
@param 	spacegroup	crystal space group.
@param 	unit_cell	6-valued vector of unit cell parameters.
@param 	wrap		0=cut atoms outside box, 1=wrap coordinates within unit cell.
@param 	Bfactor		constant for decay curve.
@param 	&paramfile	parameter file with scattering coefficients.
@return Bimage*		the new structure factors.

	All structure factors within a given resolution are calculated from
	all the selected atomic coordinates. The coordinates are fractionalized
	to fit into the given size box. If the size of the box is not given,
	it defaults to:
		x_size = (max(x_coor) - min(x_coor)) / x_sampling
	The atomic scattering profiles are read from the STAR database as the
	amplitudes and B-factors of reciprocal space gaussians. For each profile,
	a lookup table is calculated to speed up further calculations.
	The statistics of the new image is calculated.

**/
Bimage*		img_sf_from_molecule(Bmolgroup* molgroup,
				Vector3<double> origin, Vector3<long> size, Vector3<double> sampling, double resolution, 
				int spacegroup, UnitCell unit_cell, int wrap, double Bfactor, Bstring& paramfile)
{
	Batomtype*	atompar = get_atom_properties(paramfile);

	if ( sampling.volume() < 1e-6 )
		sampling[0] = sampling[1] = sampling[2] = 1;
	else
		sampling = sampling.max(0.01);
	
	if ( resolution < sampling[0] ) resolution = sampling[0];
	if ( resolution < 0.01 ) resolution = 2*sampling[0];
	
	if ( size[0] < 1 ) size[0] = (long) ((molgroup->max[0] - molgroup->min[0])/sampling[0]);
	if ( size[1] < 1 ) size[1] = (long) ((molgroup->max[1] - molgroup->min[1])/sampling[1]);
	if ( size[2] < 1 ) size[2] = (long) ((molgroup->max[2] - molgroup->min[2])/sampling[2]);
	size = size.min(MAXSID);
	size = size.max(1);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG img_sf_from_molecule: size = " << size << endl;
		
	Bimage* 	p = new Bimage(Float, TComplex, size, 1);
	p->unit_cell(unit_cell);
	p->origin(origin);
	p->sampling(sampling);
	p->fourier_type(Standard);
	p->space_group(spacegroup);
	
	double		smax = 2.0/p->real_size().length();
	if ( resolution > 0.01 ) smax = 1.0/resolution;
	
	if ( ( verbose & VERB_LABEL ) || ( verbose & VERB_PROCESS ) )
		cout << "Calculating structure factors from atomic coordinates" << endl;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Size:                           " << p->size() << " voxels" << endl;
		cout << "Origin:                         " << p->image->origin() << endl;
		cout << "Sampling:                       " << p->sampling(0) << " A/voxel" << endl;
		cout << "Resolution:                     " << resolution << " A" << endl;
		cout << "Reciprocal resolution:          " << smax << " /A" << endl;
	}
	
	mol_to_structure_factors(molgroup, p, atompar, resolution, wrap, Bfactor);
	
	atom_type_kill(atompar);
	
	p->statistics();
	
	return p;
}

double		one_sf(Bmolgroup* molgroup, Bimage* p, long i, double s, double scale, double* scat)
{
	if ( s <= 0 ) return 0;
	
	long			x, y, z;
	long			h, k, l, t, si;
	double			fraction, fraction1, f, phi;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Complex<float>	temp;
	
	p->coordinates(i, x, y, z);
	if ( x < (p->sizeX() + 1)/2 ) h = x; else h = (long)x - p->sizeX();
	if ( y < (p->sizeY() + 1)/2 ) k = y; else k = (long)y - p->sizeY();
	if ( z < (p->sizeZ() + 1)/2 ) l = z; else l = (long)z - p->sizeZ();
	
	si = (long) s;
	fraction = s - si;
	fraction1 = 1 - fraction;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next ) 
				if ( atom->sel ) {
					t = atom->tnum*MAXSCAT + si;
					f = atom->q*(fraction1*scat[t] + fraction*scat[t+1]);
					phi = MIN2PI*(h*atom->vel[0] + 
							k*atom->vel[1] + l*atom->vel[2]);
					temp = Complex<float>(cos(phi), sin(phi));
					p->set(i, p->complex(i) + f*temp);
				}
		}
	}
	
	p->set(i, p->complex(i) * scale);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mol_to_structure_factors: " << x << " " << y << " " << z << " " <<
			(p->complex(i)).real() << " " << (p->complex(i)).imag() << endl;

	return 0;
}

/*
@brief 	Calculates a 3D set of structure factors from a set of atomic coordinates.
@param 	*molgroup 	set of molecules with atomic coordinates.
@param 	*p			image structure for the new structure factors.
@param	*atompar	atom type parameters.
@param 	resolution	resolution (angstrom).
@param 	wrap		0=cut atoms outside box, 1=wrap coordinates within unit cell.
@param 	Bfactor		constant for decay curve.
@return int			0.

	All structure factors within a given resolution are calculated from
	all the selected atomic coordinates. The coordinates are fractionalized
	to fit into the given size box. If the size of the box is not given,
	it defaults to:
		x_size = (max(x_coor) - min(x_coor)) / x_sampling
	The atomic scattering profiles are read from the STAR database as the
	amplitudes and B-factors of reciprocal space gaussians. For each profile,
	a lookup table is calculated to speed up further calculations.
	The content and statistics of the new image is not checked.

**/
int			mol_to_structure_factors(Bmolgroup* molgroup, Bimage* p,
				Batomtype* atompar, double resolution, int wrap, double Bfactor)
{
	long 			i, j, h, k, l, x, y, z, yi, zi, n, t, nmax;
	long 			nscat, nmol(0), natom(0), nsel(0);
    double			s, smax, sx, sy, sz;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*  		atom;
	Batomtype*  	at;

	// Set the resolution limit
	smax = 1;
	if ( resolution ) smax = 1.0/resolution;
	resolution = 1.0/smax;
    
	double			recip_interval = 1.1*smax/MAXSCAT;
	double			inv_recip_int = 1/recip_interval;
	Vector3<double> ori = p->sampling(0)*p->image->origin();
	Vector3<double> fcoor;
	
	// Convert to fractional coordinates
	UnitCell		unit_cell = p->unit_cell();
	unit_cell.a(p->real_size()[0]);
	unit_cell.b(p->real_size()[1]);
	unit_cell.c(p->real_size()[2]);
	p->unit_cell(unit_cell);
	Matrix3			frac_mat = unit_cell.skew_matrix();

	double*			scat = get_scattering_curves(atompar, Bfactor, recip_interval, nscat);
    
    for ( nmol = 0, mol = molgroup->mol; mol; mol = mol->next ) if ( mol->sel ) {
		nmol++;
		for( res = mol->res; res; res = res->next ) {
			for ( atom = res->atom; atom; atom = atom->next, natom++ ) {
				atom->tnum = -1;
				for ( t=0, at = atompar; at && atom->tnum < 0; at = at->next, t++ )
					if ( strncmp(atom->el, at->el, 2) == 0 )
						atom->tnum = t;
				if ( atom->tnum < 0 ) atom->sel = 0;
				if ( atom->q <= 0 ) atom->sel = 0;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG mol_to_structure_factors: Atom = " << atom->el << endl;
				fcoor = atom->coord + ori;
				atom->vel = frac_mat * fcoor;
				if ( wrap == 1 ) {
					while ( atom->vel[0] < 0 ) atom->vel[0] += 1;
					while ( atom->vel[1] < 0 ) atom->vel[1] += 1;
					while ( atom->vel[2] < 0 ) atom->vel[2] += 1;
					while ( atom->vel[0] > 1 ) atom->vel[0] -= 1;
					while ( atom->vel[1] > 1 ) atom->vel[1] -= 1;
					while ( atom->vel[2] > 1 ) atom->vel[2] -= 1;
				} else if ( wrap > 1 ) {
					if ( atom->vel[0] < 0 ) atom->sel = 0;
					if ( atom->vel[1] < 0 ) atom->sel = 0;
					if ( atom->vel[0] > 1 ) atom->sel = 0;
					if ( atom->vel[1] > 1 ) atom->sel = 0;
					while ( atom->vel[2] < 0 ) atom->vel[2] += 1;
					while ( atom->vel[2] > 1 ) atom->vel[2] -= 1;
				} else {
					if ( atom->vel[0] < 0 ) atom->sel = 0;
					if ( atom->vel[1] < 0 ) atom->sel = 0;
					if ( atom->vel[2] < 0 ) atom->sel = 0;
					if ( atom->vel[0] > 1 ) atom->sel = 0;
					if ( atom->vel[1] > 1 ) atom->sel = 0;
					if ( atom->vel[2] > 1 ) atom->sel = 0;
				}
				if ( atom->sel ) nsel++;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG mol_to_structure_factors: " 
						<< atom->el << " " << atom->vel << endl;
			}
		}
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Molecules and atoms:            " << nmol << " " << natom << " (" << nsel << ")" << endl;
		cout << "Unit cell:                      " << 
			unit_cell.a() << " " << unit_cell.b() << " " << unit_cell.c() << " " << 
			unit_cell.alpha()*180/M_PI << " " << unit_cell.beta()*180/M_PI << " " << unit_cell.gamma()*180/M_PI << endl;
		cout << "Fractionalization matrix:       " << endl << frac_mat << endl;
	}
	
	// Calculate the structure factors
	// Default space group = 1
	// Estimate the number of structure factors from the spherical volume or circular area
	n = 0;
	nmax = (long) ((4.0*M_PI/3.0)*smax*smax*smax/(frac_mat[0][0]*frac_mat[1][1]*frac_mat[2][2]));
	if ( p->sizeZ() < 2 )
		nmax = (long) (M_PI*smax*smax/(frac_mat[0][0]*frac_mat[1][1]));
	// Factor required to calculate the estimated completion from the number of
	// structure factors already calculated
//	pfac = 100.0/(nmax/2 + M_PI*smax*smax/(2*frac_mat[1][1]*frac_mat[2][2]));
//	if ( p->sizeZ() < 2 )
//		pfac = 100.0/(nmax/2 + smax/frac_mat[1][1]);
	if ( verbose & VERB_PROCESS )
		cout << "Reflections estimated:          " << nmax << endl;
		
	// Scaled by the potential prefactor K divided by the volume
	double		scale = POTPREFAC/p->size().volume();
	if ( p->sizeX() > 1 ) scale /= p->sampling(0)[0];
	if ( p->sizeY() > 1 ) scale /= p->sampling(0)[1];
	if ( p->sizeZ() > 1 ) scale /= p->sampling(0)[2];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG mol_to_structure_factors: scale = " << scale << endl;
	
	long			imgsize = (long) p->size().volume();
	long			hx = p->sizeX()/2, hy = p->sizeY()/2, hz = p->sizeZ()/2;
	long			hx1 = (p->sizeX()+1)/2, hy1 = (p->sizeY()+1)/2, hz1 = (p->sizeZ()+1)/2;
	Complex<float>	temp;
	
	p->next = new Bimage(Float, TSimple, p->size(), p->images());
	float*			fom = (float *) p->next->data_pointer();
	
	// Set up FOM block with s values
	for ( l=-hz; l<hz1; l++ ) {
		z = ( l < 0 ) ? l + p->sizeZ() : l;
		sz = l*frac_mat[2][2];
		sz *= sz;
		for ( k=-hy; k<hy1; k++ ) {
			y = ( k < 0 ) ? k + p->sizeY() : k;
			sy = k*frac_mat[1][1];
			sy *= sy;
			for ( h=-hx; h<=0; h++ ) {
				x = ( h < 0 ) ? h + p->sizeX() : h;
				sx = h*frac_mat[0][0];
				sx *= sx;
				s = sqrt(sx + sy + sz);
				if ( s < smax ) {
					i = p->index(x, y, z);
					s *= inv_recip_int;
					fom[i] = s;
					n++;
				}
			}
		}
	}
	
	fom[0] = 1e-10;		// Make F0 non-zero
	
#ifdef HAVE_GCD
	dispatch_apply(imgsize, dispatch_get_global_queue(0, 0), ^(size_t i){
		if ( fom[i] ) one_sf(molgroup, p, i, fom[i], scale, scat);
	});
#else
#pragma omp parallel for
	for ( i=0; i<imgsize; i++ ) if ( fom[i] )
		one_sf(molgroup, p, i, fom[i], scale, scat);
#endif

//	cout << "first half done" << endl;
	
	// Fill in the second half from the first
	for ( l=-hz; l<hz1; l++ ) {
		z = ( l < 0 ) ? l + p->sizeZ() : l;
		zi = ( z == 0 ) ? 0 : p->sizeZ() - z;
		sz = l*frac_mat[2][2];
		sz *= sz;
		for ( k=-hy; k<hy1; k++ ) {
			y = ( k < 0 ) ? k + p->sizeY() : k;
			yi = ( y == 0 ) ? 0 : p->sizeY() - y;
			sy = k*frac_mat[1][1];
			sy *= sy;
			for ( h=1; h<hx1; h++ ) {
				sx = h*frac_mat[0][0];
				sx *= sx;
				s = sqrt(sx + sy + sz);
				if ( s < smax ) {
					i = p->index(p->sizeX() - h, yi, zi);
					j = p->index(h, y, z);
					p->set(j, (p->complex(i)).conj());
					n++;
				}
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Reflections calculated:         " << n << endl << endl;

	delete[] scat;

	return 0;
}

double*		get_potential_curves(Batomtype* atompar, double interval)
{
	long			i, j, n, t;
	double			interval2 = -4*M_PI*M_PI*interval*interval, prefac, expfac;
	Batomtype*  	at;
	
    // Get the atomic properties
	for ( n=0, at = atompar; at; at = at->next, n++ ) ;
	
    // Calculate the atomic scattering profiles
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_potential_curves: n=" << n << " interval=" << interval << endl;
		
	double*		apot = new double[MAXRAD*n];
	for ( i=0; i<MAXRAD*n; i++ ) apot[i] = 0;
	
	for ( t=0, at = atompar; at; at = at->next, t++ ) {
		apot[t*MAXRAD] = at->sfc;
		for ( j=0; j<5; j++ ) {
			prefac = M_PI*sqrt(M_PI)*at->sfa[j]/(at->sfb[j]*sqrt(at->sfb[j]));
			expfac = interval2/at->sfb[j];
			for ( i=0; i<MAXRAD; i++ )
				apot[t*MAXRAD+i] += prefac*exp(expfac*i*i);
		}
    }
	
	if ( verbose & VERB_FULL ) {
		cout << "Scattering coefficients:" << endl;
		for ( at = atompar; at; at = at->next )
			cout << tab << at->el << "(" << at->z << ")";
		cout << endl;
		for ( j=0; j<5; j++ ) {
			cout << "a" << j+1;
			for ( at = atompar; at; at = at->next )
				cout << tab << at->sfa[j];
			cout << endl;
		}
		for ( j=0; j<5; j++ ) {
			cout << "b" << j+1;
			for ( at = atompar; at; at = at->next )
				cout << tab << at->sfb[j];
			cout << endl;
		}
		cout << "c";
		for ( at = atompar; at; at = at->next )
			cout << tab << at->sfc;
		cout << endl << "Atomic radial potential curves:" << endl << "s";
		for ( at = atompar; at; at = at->next )
			cout << tab << at->el << "(" << at->z << ")";
		cout << endl;
	    for ( i=0; i<MAXRAD; i++ ) {
			cout << i*interval;
			for ( t=0, at = atompar; at; at = at->next, t++ )
				cout << tab << apot[t*MAXRAD+i];
			cout << endl;
		}
		cout << endl;
	}
	
	return apot;
}

double*		get_scattering_curves(Batomtype* atompar, double Bfactor,
				double recip_interval, long& nscat)
{
	long			i, j, n, t;
	double			s2, fact = -Bfactor/4;
	Batomtype*  	at;
	
	if ( recip_interval < 1e-10 ) {
		cerr << "Error: Reciprocal space interval not defined!" << endl;
		return NULL;
	}
	
    double			recip_int2 = recip_interval*recip_interval/4;
	
    // Get the atomic properties
	for ( n=0, at = atompar; at; at = at->next, n++ ) ;
	
    // Calculate the atomic scattering profiles
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_scattering_curves: nscat = " << n << endl;
		
	double*			scat = new double[MAXSCAT*n];
	for ( i=0; i<MAXSCAT*n; i++ ) scat[i] = 0;
	
    for ( i=0; i<MAXSCAT; i++ ) {
    	s2 = i*i*recip_int2;
		for ( t=0, at = atompar; at; at = at->next, t++ ) {
    		scat[t*MAXSCAT+i] = at->sfc;
			for ( j=0; j<5; j++ )
				scat[t*MAXSCAT+i] += at->sfa[j]*exp(-at->sfb[j]*s2);
    		scat[t*MAXSCAT+i] *= exp(fact*s2);
	    }
    }
	
	if ( verbose & VERB_FULL ) {
		cout << "Scattering coefficients:" << endl;
		for ( at = atompar; at; at = at->next )
			cout << tab << at->el << "(" << at->z << ")";
		cout << endl;
		for ( j=0; j<5; j++ ) {
			cout << "a" << j+1;
			for ( at = atompar; at; at = at->next )
				cout << tab << at->sfa[j];
			cout << endl;
		}
		for ( j=0; j<5; j++ ) {
			cout << "b" << j+1;
			for ( at = atompar; at; at = at->next )
				cout << tab << at->sfb[j];
			cout << endl;
		}
		cout << "c";
		for ( at = atompar; at; at = at->next )
			cout << tab << at->sfc;
		cout << endl << "Scattering curves:" << endl << "s";
		for ( at = atompar; at; at = at->next )
			cout << tab << at->el << "(" << at->z << ")";
		cout << endl;
	    for ( i=0; i<MAXSCAT; i++ ) {
			s2 = i*i*recip_int2;
			cout << sqrt(s2);
			for ( t=0, at = atompar; at; at = at->next, t++ )
				cout << tab << scat[t*MAXSCAT+i];
			cout << endl;
		}
		cout << endl;
	}
	
	nscat = n;
	
	return scat;
}

map<string, vector<double>>	get_scattering_curves(map<string,Bcomptype>& types,
			double Bfactor, double recip_interval)
{
	if ( recip_interval < 1e-10 ) recip_interval = 0.01;
	
	long			i;
	double			s2, fact = -Bfactor/4;
	double			recip_int2 = recip_interval*recip_interval/4;
	
    // Calculate the atomic scattering profiles
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_scattering_curves: nscat = " << types.size() << endl;

	map<string, vector<double>>	scat = calculate_scattering_curves(
				types, recip_interval, 1.0);

	// Impose the B-factor
	for ( auto c: scat ) {
		i = 0;
		for ( auto v: c.second ) {
			s2 = i*i*recip_int2;
			v *= exp(fact*s2);
			i++;
		}
	}
	
	return scat;
}

