/**
@file	molecule_to_map.h 
@brief	Header file for functions to calculate a 3D map from atomic coordinates 
@author Bernard Heymann 
@date	Created: 19970914
@date	Modified: 20150128
**/
 
#include "rwmolecule.h"
#include "rwimg.h"
#include "rwatomprop.h"
#include "symmetry.h"

#define	POTPREFAC	47.87801	// Atomic potential prefactor K = h^2 / (2*PI*e*mo)

// Function prototypes
Bimage* 	img_from_molecule(Bmolgroup* molgroup, Vector3<double> origin, 
				Vector3<long> size, Vector3<double> sam, 
				double resolution, double Bfactor,
				int wrap, int gextype, int spacegroup, UnitCell unit_cell);
int 		compare_mol_map(Bmolgroup* molgroup, Bimage* pcalc, Bimage* pimg);
Bimage*		img_sf_from_molecule(Bmolgroup* molgroup,
				Vector3<double> origin, Vector3<long> size, 
				Vector3<double> sam, double resolution, 
				int spacegroup, UnitCell unit_cell, 
				int wrap, double Bfactor, Bstring& paramfile);
double*		get_potential_curves(Batomtype* atompar, double interval);
double*		get_scattering_curves(Batomtype* atompar, double Bfactor,
				double recip_interval, long& nscat);

