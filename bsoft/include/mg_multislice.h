/**
@file	mg_multislice.h
@brief	Generates and manipulates projects for multislice calculations 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20221024
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "rwmolecule.h"

// Function prototypes
Bimage*		img_calc_wave_propagator(Vector3<long> size, Vector3<double> sam,
				double thickness, double volt);
Bimage*		img_calc_multi_slice(Bimage* pgrate, double thickness, double volt);
Bimage*		img_calc_multi_slice_correct(Bimage* pgrate, double thickness, double volt);
Bimage*		img_calc_potential(Bmolgroup* molgroup, Vector3<long> size, Vector3<double> origin,
				Vector3<double> sam, double thickness, double resolution, double Bfactor, 
				Bstring& paramfile, int type);
int			img_calc_phase_grating(Bimage* p, double volt);
int			img_apply_complex_CTF(Bimage* p, CTFparam& cp);
Bproject*	project_generate(int nfield, int nmg, int npart,
				Vector3<double>  pixel_size, double img_origin,
				CTFparam& cp, double def_min, double def_max, double dose,
				double tsigma, Bstring& fieldbase, Bstring& mgbase, Bstring& partbase,
				int fieldnumber, int mgnumber, int partnumber);
Bproject*	project_generate_asu(Bstring& symmetry_string,
				Vector3<double>  pixel_size, double img_origin,
				double theta_step, double phi_step,
				CTFparam& cp, double defocus, double dose, Bstring& mgbase, Bstring& partbase);
int			project_generate_potential(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution, 
				double Bfactor, int pottype, Bstring& paramfile);
int			project_generate_image(Bproject* project, double thickness, double resolution);
int			project_apply_distortions(Bproject* project, int poisson, double gauss, double kmtf);
