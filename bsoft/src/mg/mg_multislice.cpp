/**
@file	mg_multislice.cpp
@brief	Generates and manipulates projects for multislice calculations 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20221024
**/

#include "mg_processing.h"
#include "mg_ctf.h"
#include "molecule_to_map.h"
#include "mol_transform.h"
#include "mol_edit.h"
#include "Complex.h"
#include "Vector3.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the wave propagation function between slices.
@param 	size		size of projection image (z = 1).
@param 	sam			pixel size in x and y, slice thickness in z.
@param 	thickness	slice thickness (in angstrom).
@param 	volt		acceleration voltage (volt).
@return Bimage*	 	wave propagation function image.
**/
Bimage*		img_calc_wave_propagator(Vector3<long> size, Vector3<double> sam,
				double thickness, double volt)
{
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->sampling(sam);
	
	long			i, h, k, x, y;
	double			wavelength = electron_wavelength(volt);
	double			arg, sx2, sy2;
	double			afactor = M_PI*wavelength*thickness;
	double			nfactor = 1.0/(p->sizeX()*p->sizeY());	// Normalization for passing through FFT's
	double			xscale2 = 1.0/(p->sizeX()*p->sizeX()*sam[0]*sam[0]);
	double			yscale2 = 1.0/(p->sizeY()*p->sizeY()*sam[1]*sam[1]);
	
	for ( i=y=0; y<p->sizeY(); ++y ) {
		k = y;
		if ( k > (p->sizeY()-1)/2 ) k -= p->sizeY();
		sy2 = k*k*yscale2;
		for ( x=0; x<p->sizeX(); ++x, ++i ) {
			h = x;
			if ( h > (p->sizeX()-1)/2 ) h -= p->sizeX();
			sx2 = h*h*xscale2;
			arg = afactor*(sx2 + sy2);
			p->set(i, Complex<double>(cos(arg), sin(arg)) * nfactor);
		}
	}
	
//	write_img("propagator.map", p);
	
	return p;
}

/**
@brief 	Simulates the electron imaging process using a multi-slice approach.
@param 	*pgrate			phase grating multi-image.
@param 	thickness		slice thickness (in angstrom).
@param 	volt	 		acceleration voltage (volt).
@return Bimage*	 		simulated projection image transform.

	The passage of the electron beam is simulated as the interaction of
	a planar wave with successive planar phase gratings spaced at regular 
	intervals, with the wave propagated between the 2D gratings. The
	phase gratings are derived from slabs of the atomic potential 
	calculated from the atomic structure using scattering profiles.
	Note: The final product is a 2D transform of the exit wave.
	Reference: Cowley, J. M. (1995) Diffraction Physics. 3rd Rev. Ed. 
		Elsevier Science, Amsterdam.

**/
Bimage*		img_calc_multi_slice(Bimage* pgrate, double thickness, double volt)
{
	Vector3<long>	size(pgrate->size());
	
	// wavelength & size - loop one slice -> propagator (one 2D image)
	Bimage*			prop = img_calc_wave_propagator(size, pgrate->sampling(0), thickness, volt);
	
	// For all slices: Multiply previous exitwave with phasegrating - FFT - 
	// multiply with propagator - FFT-1 - next slice -> exitwave (one 2D image)
	long   			n;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->label(pgrate->label());
	p->fourier_type(NoTransform);
	p->origin(pgrate->image->origin());
	p->sampling(pgrate->sampling(0));
	
	long			slice_size(size[0]*size[1]);
	Complex<double>	cv(1);
	for ( n=0; n<slice_size; ++n ) p->set(n, cv);
	
	Bimage* 		ptemp;
	for ( n=0; n<pgrate->images(); ++n ) {
//		cout << n << endl;
		ptemp = pgrate->extract(n);					// Get phase grating slice
		p->complex_product(ptemp);				// Multiply with current wave
		delete ptemp;
		p->fft(FFTW_FORWARD, 0);		// Normalization in propagator
		if ( n < pgrate->images()-1 ) {
			p->complex_product(prop);			// Multiply with propagator
			p->fft(FFTW_BACKWARD, 0);	// Normalization in propagator
		}
	}
	
	delete prop;
	
	return p;
}

Bimage*		img_calc_multi_slice_correct(Bimage* pgrate, double thickness, double volt)
{
	Vector3<long>	size(pgrate->size());
	
	// wavelength & size - loop one slice -> propagator (one 2D image)
	Bimage*			prop = img_calc_wave_propagator(size, pgrate->sampling(0), thickness, volt);
	
	// For all slices: Multiply previous exitwave with phasegrating - FFT -
	// multiply with propagator - FFT-1 - next slice -> exitwave (one 2D image)
	long   i;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->label(pgrate->label());
	p->fourier_type(Standard);
	p->origin(pgrate->image->origin());
	p->sampling(pgrate->sampling(0));
	
//	long	slice_size = size[0]*size[1];
//	Complex<double>	cv(1);
//	for ( i=0; i<slice_size; i++ ) p->set(i, cv);
	
	Bimage* 		ptemp;
	for ( i=0; i<pgrate->images(); i++ ) {
		ptemp = pgrate->extract(i);					// Get phase grating slice
		ptemp->fft(FFTW_FORWARD, 0);		// Normalization in propagator
		p->add(ptemp);					// Add to current wave
		delete ptemp;
		if ( i < pgrate->images()-1 )
			p->complex_product(prop);			// Multiply with propagator
	}
	
	delete prop;

	p->fft(FFTW_BACKWARD, 0, Real);	// Normalization in propagator

	return p;
}

/**
@brief 	Calculates the atomic potential.
@param 	*molgroup	set of molecules.
@param 	size		size of projection image (z = 1).
@param 	origin		origin in x and y.
@param 	sam			voxel size.
@param 	thickness	slice thickness (in angstrom).
@param 	resolution	resolution limit (angstrom).
@param 	Bfactor		overall temperature factor.
@param 	&paramfile	parameter file for atomic scattering coefficients.
@param 	type		type of potential calculation: 0=reciprocal space, 1=real space, 2=gaussian
@return Bimage*	 	complex potential image.
**/
Bimage*		img_calc_potential(Bmolgroup* molgroup, Vector3<long> size, Vector3<double> origin,
				Vector3<double> sam, double thickness, double resolution, double Bfactor,
				Bstring& paramfile, int type)
{
	// Atoms & scat cross-sections - Loop slices - calc 2D reciprocal space structure factors 
	long			i, j, k;
	int				nslices(0);
	Bmolgroup**		slice_molgroup = molgroup_split_into_slices(molgroup, thickness, nslices);
	
	if ( verbose )
		cout << "Calculating the atomic potential" << endl;
		
	int				spacegroup(0);
	UnitCell		unit_cell(sam[0]*size[0],sam[1]*size[1],sam[2]*size[2],M_PI_2,M_PI_2,M_PI_2);
	Bimage*			ptemp = NULL;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, nslices);
	long   imgsize = p->sizeX()*p->sizeY()*p->channels();
	p->label(molgroup->comment.str());
	p->unit_cell(unit_cell);
	p->sampling(sam[0], sam[1], thickness);
	p->space_group(spacegroup);
	
	for ( i=0; i<nslices; i++ ) {
		if ( type ) {
			origin[2] = -slice_molgroup[i]->min[2]/sam[2] - 0.5;
			size[2] = (int) (thickness/sam[2] + 1);
			ptemp = img_from_molecule(slice_molgroup[i], origin, size, sam, 
					resolution, 0.001, 0, 2 - type, spacegroup, unit_cell);
			ptemp->project('z', 1);
			ptemp->simple_to_complex();
		} else {
//			origin[2] = -slice_molgroup[i]->min[2]/sam[2] - 0.5;
//			origin[2] = -i*thickness;
//			origin[2] = 0;
			size[2] = 1;
			sam[2] = thickness;
			ptemp = img_sf_from_molecule(slice_molgroup[i], origin, size, sam, 
					resolution, spacegroup, unit_cell, 2, Bfactor, paramfile);
			ptemp->fft(FFTW_BACKWARD, 0, Real);	// No normalization - get actual potential at the beam
		}
		for ( j=0, k=i*imgsize; j<imgsize; j++, k++ ) p->set(k, (*ptemp)[j]);
		delete ptemp;
		p->origin(i, origin);
		if ( verbose )
			cout << "Slice " << i+1 << " done" << endl << endl;
	}

	p->fourier_type(NoTransform);
	
	for ( i=0; i<nslices; i++ )
		molgroup_kill(slice_molgroup[i]);
	
	delete[] slice_molgroup;
	
	return p;
}

/**
@brief 	Calculates the phase grating approximation from the atomic potential.
@param 	*p				atomic potential image (modified).
@param 	volt			acceleration voltage (volt).
@return int 				0.

	All calculations are complex.

**/
int			img_calc_phase_grating(Bimage* p, double volt)
{
	p->simple_to_complex();
	
	long   			i, j, n;
	long   			imgsize = p->sizeX()*p->sizeY();
//	double			wavelength = electron_wavelength(volt);
	double			arg;
//	double			sigma = 0.0208886 * wavelength *
//						sqrt(1 + 5.8886579e-4/(wavelength*wavelength));	// According to Roar Kilaas
	double			sigma = 1e-10*TWOPI*ECHARGE/(PLANCK*LIGHTSPEED*beta(volt));
	
//	if ( verbose & VERB_FULL )
	if ( verbose )
		cout << "Sigma (interaction factor):     " << sigma << endl;

	for ( i=n=0; n<p->images(); n++ ) {
		for ( j=0; j<imgsize; i++, j++ ) {
			arg = -sigma*(p->complex(i)).real();
			p->set(i, Complex<double>(cos(arg), sin(arg)));
		}
	}
	
	return 0;
}

/**
@brief 	Applies a complex CTF function to a Fourier transform.
@param 	*p				complex Fourier transform (modified).
@param 	def_avg			defocus minimum (angstrom).
@param 	def_dev			defocus maximum (angstrom).
@param 	ast_angle 		astigmatism angle (radians).
@param 	volts 			acceleration voltage (volts).
@param 	Cs				spherical aberration (angstrom).
@param 	Cc				chromatic aberration (angstrom).
@param 	amp_shift		amplitude contrast phase shift (radian).
@param 	alpha			beam source size (radians).
@param 	energy_spread	effective energy spread (fraction: typically 10^-5)
@return int 			0.


	The CTF is applied as a multiplication with a complex number:
		new_datum.re = datum.re*amp_fac*cos(dphi) - datum.im*phi_fac*sin(dphi)
		new_datum.re = datum.re*phi_fac*sin(dphi) + datum.im*amp_fac*cos(dphi)
	Both input and output are complex transforms.


**/
int			img_apply_complex_CTF(Bimage* p, CTFparam& cp)
{
	
	long			i, n, x, y, z;
	double			xx, yy, zz, sx2, sy2, sz2, s2, env;
	Complex<double>	ctf;
	Vector3<double>	freq_scale(1.0/p->real_size());
	
	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) )
		cout << "Applying a CTF to a complex Fourier transform" << endl;

	if ( verbose & VERB_PROCESS ) {
		cp.show();
	}
	
	for ( i=n=0; n<p->images(); n++ ) {
		for ( z=0; z<p->sizeZ(); z++ ) {
			zz = z;
			if ( z > (p->sizeZ() - 1)/2 ) zz -= p->sizeZ();
			sz2 = freq_scale[2]*zz;
			sz2 *= sz2;
			for ( y=0; y<p->sizeY(); y++ ) {
				yy = y;
				if ( y > (p->sizeY() - 1)/2 ) yy -= p->sizeY();
				sy2 = freq_scale[1]*yy;
				sy2 *= sy2;
				for ( x=0; x<p->sizeX(); x++, i++ ) {
					xx = x;
					if ( x > (p->sizeX() - 1)/2 ) xx -= p->sizeX();
					sx2 = freq_scale[0]*xx;
					sx2 *= sx2;
					s2 = sx2 + sy2 + sz2;
					ctf = cp.calculate_complex(s2, atan2(yy,xx));
					env = cp.partial_coherence_and_energy_spread(s2);
					p->set(i, ((*p)[i] * ctf) * env);
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates a project for multislice simulations.
@param 	nfield			number of fields-of-view.
@param 	nmg				number of micrographs per field-of-view.
@param 	npart			number of particles per micrograph.
@param 	pixel_size		micrograph pixel size.
@param 	img_origin		image origin within the simulation box.
@param 	cp				CTF parameters).
@param 	def_min			defocus minimum (angstrom).
@param 	def_max			defocus maximum (angstrom).
@param 	dose			electron dose (e/angstrom^2).
@param 	tsigma			translation standard deviation (pixels).
@param 	&fieldbase		field base name.
@param 	&mgbase			micrograph base name.
@param 	&partbase		particle image base name.
@param 	fieldnumber		field-of-view number.
@param 	mgnumber		micrograph number.
@param 	partnumber		particle number.
@return Bproject*			project structure.
**/
Bproject*	project_generate(int nfield, int nmg, int npart,
				Vector3<double> pixel_size, double img_origin,
				CTFparam& cp, double def_min, double def_max, double dose,
				double tsigma, Bstring& fieldbase, Bstring& mgbase, Bstring& partbase,
				int fieldnumber, int mgnumber, int partnumber)
{
	random_seed();
	
	double			def_range = def_max - def_min;
	if ( def_range < 0 ) def_range = 0;
	
	int				i, j, k, block = 0;
	double			irm = 1.0/get_rand_max(), irm2 = 2*irm;
	Vector3<float>	t;
	View*			v = new View[npart];
	Bproject*		project = new Bproject;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bstring			field_id, mg_id;
	
	if ( verbose ) {
		cout << "Generating a new project:" << endl;
		cout << "Number of fields:               " << nfield << endl;
		cout << "Number of micrographs:          " << nmg << endl;
		cout << "Number of particles:            " << npart << endl;
	}

	for ( i=0; i<nfield; i++ ) {
		field_id = fieldbase + Bstring(fieldnumber, "_%06d");
		if ( project->field ) field = field_add(&field, field_id);
		else field = field_add(&project->field, field_id);
		for ( k=0; k<npart; k++ ) {
			v[k] = View(random()*irm2 - 1, random()*irm2 - 1, random()*irm2 - 1, M_PI*(random()*irm2 - 1));
			v[k].normalize();
		}
		for ( j=0; j<nmg; j++ ) {
			mg_id = mgbase + Bstring(fieldnumber, "_%03d") + Bstring(mgnumber, "_%03d");
			if ( field->mg ) mg = micrograph_add(&mg, mg_id);
			else mg = micrograph_add(&field->mg, mg_id);
			mg->block = block;
			block++;
			mg->fpart = partbase + Bstring(fieldnumber, "_%03d") + Bstring(mgnumber, "_%03d.spi");
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			mg->ctf->update(cp);
			mg->ctf->defocus_average(def_range*random()*irm + def_min);
			mg->pixel_size = pixel_size;
			mg->dose = dose;
			mg->box_size[0] = mg->box_size[1] = (int) (2*img_origin);
			mg->box_size[2] = 1;
			mg->ctf->zero(1);
			mg->ctf->baseline(0, 1);
			mg->ctf->envelope(0, 3);
			mg->ctf->envelope(1, -M_PI*M_PI*mg->ctf->alpha()*mg->ctf->alpha()*mg->ctf->defocus_average()*mg->ctf->defocus_average());
			for ( k=0; k<npart; k++ ) {
				if ( mg->part ) part = particle_add(&part, k+partnumber);
				else part = particle_add(&mg->part, k+partnumber);
				if ( tsigma > 0 ) t = vector3_xy_random_gaussian(0.0, (double)tsigma);
				part->view = v[k];
				part->loc[0] = part->loc[1] = img_origin;
				part->ori = t + img_origin;
				part->ori[2] = 0;
			}
			mgnumber++;
		}
		fieldnumber++;
	}	
	
	delete[] v;
	
	return project;
}

/**
@brief 	Generates a project for multislice calculations of an asymmetric unit.
@param 	&symmetry_string   symmetry designation.
@param 	pixel_size		micrograph pixel size.
@param 	img_origin		image origin within the simulation box.
@param 	theta_step		step size in theta (radians).
@param 	phi_step		step size in phi (radians).
@param 	cp				CTF parameters).
@param 	defocus			defocus minimum (angstrom).
@param 	dose			electron dose (e/angstrom^2).
@param 	&mgbase			micrograph base name.
@param 	&partbase		particle image base name.
@return Bproject*			project structure.
**/
Bproject*	project_generate_asu(Bstring& symmetry_string,
				Vector3<double> pixel_size, double img_origin,
				double theta_step, double phi_step,
				CTFparam& cp, double defocus, double dose, Bstring& mgbase, Bstring& partbase)
{
	int				k, npart;
	Bproject*		project = new Bproject;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;

	Bsymmetry 	sym(symmetry_string);
	View*		view = asymmetric_unit_views(sym, theta_step, phi_step, 0);
	View*		v;
	for ( v=view, npart=0; v; v=v->next ) npart++;

	field = field_add(&project->field, mgbase);
	mg = micrograph_add(&field->mg, mgbase);
	mg->block = 0;
	mg->fpart = partbase + ".spi";
	if ( !mg->ctf ) mg->ctf = new CTFparam;
	mg->ctf->update(cp);
	mg->ctf->defocus_average(defocus);
	mg->pixel_size = pixel_size;
	mg->dose = dose;
	mg->box_size[0] = mg->box_size[1] = (int) (2*img_origin);
	mg->box_size[2] = 1;
	for ( k=0, v=view; k<npart && v; k++, v=v->next ) {
		if ( mg->part ) part = particle_add(&part, k+1);
		else part = particle_add(&mg->part, k+1);
		part->view = *v;
		part->loc = part->ori = Vector3<float>(img_origin, img_origin, 0);
	}
	
	kill_list((char *) view, sizeof(View));
	
	return project;
}

/**
@brief 	Generates potential images using a multislice calculation.
@param 	*molgroup		molecule group structure.
@param 	*water			block of water as solvent.
@param 	*project		project structure with parameters.
@param 	&fieldname		selected field name (if "" do all).
@param 	&mgname			selected micrograph (if "" do all).
@param 	partselect		selected particle ( if <1 do all).
@param 	size			size of simulation block (angstrom).
@param 	thickness		thickness of slices for the multislice calculation (angstrom).
@param 	resolution		resolution for the multislice calculation (angstrom).
@param 	Bfactor			B-factor to apply to the multislice calculation (angstrom^2).
@param 	pottype			type of potential to calculate (???).
@param 	&paramfile		parameter file for scattering curves (???).
@return int				0.
**/
int			project_generate_potential(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int pottype, Bstring& paramfile)
{
	random_seed();
	
	Vector3<double>			sam, box;
	sam = project->field->mg->pixel_size;
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];
	
	if ( size.volume() < 1 ) {
		size[0] = (long) (molgroup->box[0]/sam[0]);
		size[1] = (long) (molgroup->box[1]/sam[1]);
		size[2] = (long) (molgroup->box[2]/sam[2]);
	}
	if ( water ) {
		if ( molgroup->box.length2() < water->box.length2() ) {
			size[0] = (long) (water->box[0]/sam[0]);
			size[1] = (long) (water->box[1]/sam[1]);
			size[2] = (long) (water->box[2]/sam[2]);
		}
		box = sam * size;
		water->box = Vector3<float>(box[0], box[1], box[2]);
	}
	box = sam * size;
	molgroup->box = Vector3<float>(box[0], box[1], box[2]);
	
	Vector3<double>		sampling3(sam[0], sam[1], sam[2]);
	
	if ( verbose )
		cout << "Simulation box:                 " << molgroup->box << " = " << molgroup->box.volume() << " A^3" << endl;
	
	int				j;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bmolgroup*		mol_rot = NULL;
	Bmolgroup*		mol_sim = NULL;
	Bimage			*ppot;
	Bstring			potname;
	double			ranshift = molgroup->box[0];
	Vector3<float>	t, rot_ori;
	Vector3<double>	ori;
		
	for ( field = project->field; field; field = field->next ) 
		if ( fieldname.length() < 1 || field->id == fieldname ) {
		if ( verbose )
			cout << "Calculating potentials for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next )
			if ( mgname.length() < 1 || mg->id == mgname ) {
			if ( verbose )
				cout << "Calculating potentials for micrograph " << mg->id << endl;
			for ( part = mg->part; part; part = part->next )
				if ( partselect < 1 || part->id == partselect ) {
//				potname = number_filename(mg->fpart, part->id, 3);
				potname = mg->fpart;
				potname = potname.pre_rev('.') + Bstring(part->id, "_%03d.") + potname.post_rev('.');
				if ( access(potname.c_str(), F_OK) != 0 ) { // Skip if the file already exists
					if ( verbose )
						cout << "Calculating potentials for particle " << part->id << endl;
					t = (part->ori - part->loc) * sampling3;
//					cout << "origin=" << origin << " shift=" << t << endl;
					mol_rot = molgroup_rotate_from_view(molgroup, part->view, rot_ori, t);
					if ( water ) {
						mol_sim = molgroup_copy(water);
						molgroup_coor_shift_PBC(mol_sim, vector3_random(-ranshift, ranshift));
						molgroup_insert(mol_sim, mol_rot, 2);
						molgroup_kill(mol_rot);
					} else {
						mol_sim = mol_rot;
					}
					ori = {part->ori[0], part->ori[1], part->ori[2]};
					ppot = img_calc_potential(mol_sim, size, ori, sam, thickness,
						resolution, Bfactor, paramfile, pottype);
					ppot->complex_to_real();
					for ( j=0; j<ppot->images(); j++ ) {
						ppot->image[j].view(part->view);
						ppot->image[j].origin(part->ori);
					}
					part->fpart = potname;
					write_img(potname, ppot, 0);
					molgroup_kill(mol_sim);
					delete ppot;
				}
			}
		}
	}
	
	return 0;
}

/**
@brief 	Generates final images from a multislice calculation.
@param 	*project		project structure with parameters.
@param 	thickness			thickness of slices for the multislice calculation (angstrom).
@param 	resolution		resolution for the multislice calculation (angstrom).
@return int						0.
**/
int			project_generate_image(Bproject* project, double thickness, double resolution)
{
	Bstring			potname = project->field->mg->fpart;
	potname = potname.pre_rev('.') + Bstring(project->field->mg->part->id, "_%03d.") + potname.post_rev('.');
	
	Bimage*			ppot = read_img(potname, 0, 0);
	Vector3<long>	size(ppot->sizeX(), ppot->sizeY(), 1);
	delete ppot;
	
	Vector3<double>	realsize(size);
	realsize *= project->field->mg->pixel_size;
	
	if ( resolution < 2*project->field->mg->pixel_size[0] )
		resolution = 2*project->field->mg->pixel_size[0];
	
	long			i, j, k, npart;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;
	Bparticle*		part = NULL;
	Bimage			*p, *pone;
	int				spacegroup = 1;
		
	long	imgsize = 2*size[0]*size[1];
	UnitCell		unit_cell(realsize[0], realsize[1], realsize[2], M_PI_2, M_PI_2, M_PI_2);
	
	for ( field = project->field; field; field = field->next ) {
		if ( verbose )
			cout << "Calculating images for field " << field->id << endl;
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose )
				cout << "Calculating images for micrograph " << mg->id << endl;
			for ( npart=0, part = mg->part; part; part = part->next, npart++ ) ;
			p = new Bimage(Float, TComplex, size[0], size[1], 1, npart);
			p->unit_cell(unit_cell);
			p->space_group(spacegroup);
			p->sampling(mg->pixel_size);
			for ( i=0, part = mg->part; part; part = part->next, i++ ) {
				if ( verbose )
					cout << "Calculating image for particle " << part->id << endl;
				potname = mg->fpart;
				potname = potname.pre_rev('.') + Bstring(part->id, "_%03d.") + potname.post_rev('.');
				ppot = read_img(potname, 1, -1);
				ppot->simple_to_complex();
				img_calc_phase_grating(ppot, mg->ctf->volt());
				pone = img_calc_multi_slice(ppot, thickness, mg->ctf->volt());
				for ( j=i*imgsize, k=0; k<imgsize; j++, k++ ) p->set(j, (*pone)[k]);
				p->image[i] = ppot->image[0];
				delete ppot;
				delete pone;
			}
			if ( mg->ctf->defocus_average() )
				img_apply_complex_CTF(p, *mg->ctf);
			p->fft(FFTW_BACKWARD, 2, Real);
			p->complex_to_intensities();
			p->statistics();
			write_img(mg->fpart, p, 0);
			delete p;
		}
	}
	
	return 0;
}

/**
@brief 	Applies imaging distortions to the final images from a multislice calculation.
@param 	*project		project structure with parameters.
@param 	poisson			flag to add Poisson noise.
@param 	gauss			width of gaussian noise to add (0=no noise).
@param 	kmtf			mass transfer decay constant (0=no decay).
@return int				0.
**/
int			project_apply_distortions(Bproject* project, int poisson, double gauss, double kmtf)
{
	double			dose_per_pixel;
	Bfield* 		field = NULL;
	Bmicrograph*	mg = NULL;

	Bimage*			p;
	Bstring			insert("_tf.");				// Transfer functions(s)
	if ( poisson || gauss ) insert = "_tfn.";	// Transfer functions(s) and noise
	
	for ( field = project->field; field; field = field->next ) {
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( verbose )
				cout << "Applying distortions to micrograph " << mg->id << endl;
			p = read_img(mg->fpart, 1, -1);
			if ( p->compound_type() == TComplex ) {
				p->fft(FFTW_BACKWARD, 2, Real);
				p->complex_to_intensities();
			} else {
				p->change_type(Float);
			}
			dose_per_pixel = mg->dose*mg->pixel_size[0]*mg->pixel_size[1];
			if ( dose_per_pixel )
				p->rescale(dose_per_pixel/p->average(), 0);
			if ( poisson )
//				img_add_poisson_noise(p);
				p->noise_poisson(p->average());
			if ( gauss < 1000 )
//				img_add_gaussian_noise(p, gauss);
				p->noise_gaussian(0, p->standard_deviation()/sqrt(gauss));
			if ( kmtf > 0 )
				p->fspace_weigh_B_factor(4*kmtf);
			mg->fpart = mg->fpart.pre_rev('.') + insert + mg->fpart.post_rev('.');
			write_img(mg->fpart, p, 0);
			delete p;
		}
	}
	
	return 0;
}


