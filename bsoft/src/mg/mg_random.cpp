/**
@file	mg_random.cpp
@brief	Functions for applying some randomization to parameters
@author Bernard Heymann
@date	Created: 20010206
@date	Modified: 20140225
**/

#include "mg_random.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Introduces random errors into micrograph defocus averages.

	The random deviations are obtained from a Gaussian distribution
	generator.

@param 	*project	project parameter structure.
@param 	std			standard deviation of random error (angstroms).
@return int 				0.
**/
int			project_random_defocus(Bproject* project, double std)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	long			n(0);
	double			rnum, rdev(0);

	if ( verbose & VERB_PROCESS ) {
		cout << "Introducing random errors into defocus average:" << endl;
		cout << "Target standard deviation:      " << std << " A" << endl;
	}
	
	random_seed();
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next ) {
			if ( !mg->ctf ) mg->ctf = new CTFparam;
			do {
				rnum = random_gaussian(0, std);
			} while ( 2*fabs(rnum) > mg->ctf->defocus_average() );
			mg->ctf->defocus_average(mg->ctf->defocus_average() + rnum);
			n++;
			rdev += rnum*rnum;
		}
	
	rdev = sqrt(rdev/n);
	
	if ( verbose & VERB_PROCESS )
		cout << "Standard deviation:             " << rdev << " A" << endl;

	return 0;
}

/**
@brief 	Introduces random errors into particle origins.

	The random deviations are obtained from a Gaussian distribution
	generator.

@param 	*project	project parameter structure.
@param 	std			standard deviation of random error (pixels).
@return int 				0.
**/
int			project_random_origins(Bproject* project, double std)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	long				n(0);
	double				ori_dev(0);
	Vector3<float>		rvec;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Introducing random errors into origins:" << endl;
		cout << "Target standard deviation:      " << std << " pixels" << endl;
	}

	random_seed();
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				rvec = vector3_xy_random_gaussian(0.0, std);
				part->ori += rvec;
				n++;
				ori_dev += rvec.length2();
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					rvec = vector3_xy_random_gaussian(0.0, std);
					part->ori += rvec;
					n++;
					ori_dev += rvec.length2();
				}
			}
		}
	}
	
	ori_dev = sqrt(ori_dev/n);
	
	if ( verbose & VERB_PROCESS )
		cout << "Standard deviation:             " << ori_dev << " pixels" << endl << endl;

	return 0;
}

/**
@brief 	Replace the original particle orientations by random views.

	The random views are generated in the ranges {[0,1], [0,1], [0,1], [-PI,PI]}.

@param 	*project	project parameter structure.
@return int					0.
**/
int			project_random_views(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;

	random_seed();	 // Seed random generator = process ID
	
	if ( verbose & VERB_PROCESS )
		cout << "Replacing particle views with random views" << endl << endl;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( verbose & VERB_FULL )
					cout << "Original view: " << part->view << endl;
				part->view = random_view();
				if ( verbose & VERB_FULL )
					cout << "Random view: " << part->view << endl;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( verbose & VERB_FULL )
						cout << "Original view: " << part->view << endl;
					part->view = random_view();
					if ( verbose & VERB_FULL )
						cout << "Random view: " << part->view << endl;
				}
			}
		}
	}

	return 0;
}

/**
@brief 	Replace the original particle orientations by random side views for a helix.

	The random views are generated in the ranges {[0,1], [0,1], [0,1], [-PI,PI]}.

@param 	*project	project parameter structure.
@return int					0.
**/
int			project_random_helical_views(Bproject* project)
{
	if ( !project ) return 0;
	
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	Euler				euler;
	double				irm = TWOPI/get_rand_max();

	random_seed();	 // Seed random generator = process ID
	
	if ( verbose & VERB_PROCESS )
		cout << "Replacing particle views with random helical views" << endl << endl;
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( verbose & VERB_FULL )
					cout << "Original view: " << part->view << endl;
				euler = Euler(part->view);
				euler[2] = random()*irm;
				part->view = euler.view();
				if ( verbose & VERB_FULL )
					cout << "Random view: " << part->view << endl;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( verbose & VERB_FULL )
						cout << "Original view: " << part->view << endl;
					euler = Euler(part->view);
					euler[1] = M_PI_2;
					euler[2] = random()*irm;
					part->view = euler.view();
					if ( verbose & VERB_FULL )
						cout << "Random view: " << part->view << endl;
				}
			}
		}
	}

	return 0;
}

double		particle_random_view(Bparticle* part, double std)
{
	Vector3<double>		vec = part->view.vector3();
	Vector3<double>		axis = vector3_random(-1.0, 1.0);
	
	axis.normalize();
	axis = axis.cross(vec);
	axis.normalize();
	
	double				angle = random_gaussian(0.0, std);
	Matrix3				mat = Matrix3(axis, angle);
	
	axis = mat * vec;
	part->view = View(axis[0], axis[1], axis[2], part->view.angle());
	
	return axis.angle(vec);
}

/**
@brief 	Introduces random errors into particle views.

	The random deviations are obtained from a Gaussian distribution
	generator.

@param 	*project	project parameter structure.
@param 	std			standard deviation of random error (radians).
@return int 				0.
**/
int			project_random_views(Bproject* project, double std)
{
	if ( !project ) return 0;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	long				n(0);
	double				angle, ang_dev(0);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Introducing random errors into views:" << endl;
		cout << "Target standard deviation:      " << std*180.0/M_PI << " degrees" << endl;
	}

	random_seed();
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
//				angle = particle_random_view(part, std);
				angle = random_view_error(part->view, std);
				n++;
				ang_dev += angle*angle;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
//					angle = particle_random_view(part, std);
					angle = random_view_error(part->view, std);
					n++;
					ang_dev += angle*angle;
				}
			}
		}
	}
	
	ang_dev = sqrt(ang_dev/n);
	
	if ( verbose & VERB_PROCESS )
		cout << "Standard deviation:             " << ang_dev*180.0/M_PI << " degrees" << endl << endl;

	return 0;
}

/**
@brief 	Introduces random errors into particle magnifications.

	The random deviations are obtained from a Gaussian distribution
	generator.

@param 	*project	project parameter structure.
@param 	std			standard deviation of random error (fraction).
@return int 				0.
**/
int			project_random_magnification(Bproject* project, double std)
{
	if ( !project ) return 0;
	
	Bfield* 		field;
	Bmicrograph*	mg;
	Bparticle*		part;
	long			n(0);
	double			rnum, rdev(0);
	
	if ( std > 0.5 ) std = 0.5;

	if ( verbose & VERB_PROCESS ) {
		cout << "Introducing random errors into defocus average:" << endl;
		cout << "Target standard deviation:      " << std << " A" << endl;
	}
	
	random_seed();
	
	for ( field = project->field; field; field = field->next )
		for ( mg = field->mg; mg; mg = mg->next )
			for ( part = mg->part; part; part = part->next ) {
				do {
					rnum = random_gaussian(0, std);
				} while ( fabs(rnum) > 0.5 );
				part->mag += rnum;
				n++;
				rdev += rnum*rnum;
			}
	
	rdev = sqrt(rdev/n);
	
	if ( verbose & VERB_PROCESS )
		cout << "Standard deviation:             " << rdev << " A" << endl << endl;

	return 0;
}

/**
@author Eduardo Sanz-Garcia and David Belnap
@brief 	Replaces the original orientations by symmetry-related views.

	For each particle all the symmetry-related views are generated and
	one randomly selected.

@param 	*project	project parameter structure.
@param 	*sym		symmetry structure.
@return int					0.
**/
int			project_random_symmetry_related_views(Bproject* project, Bsymmetry& sym)
{
	if ( !project ) return 0;
	
	Bfield* 			field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
	Bparticle*			part;
	View*				v;
	double				r;			// Stores a random value in range [0,1) or [0,n)
	int					y;			// Stores a random integer in range [0,n-1]
	int					i;
	double				irm = 1.0L/(get_rand_max() + 1.0L);
	int			 		nviews(sym.order());	// Number of symmetry-related views

	random_seed();	// Seed random generator = process ID
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Substituting original orientation by randomly assigning one of the symmetry-related views:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
	}
	
	if ( project->select ) {
		for ( rec = project->rec; rec; rec = rec->next ) {
			for ( part = rec->part; part; part = part->next ) {
				if ( verbose & VERB_FULL )
					cout << "Original view vector: " << part->view << endl;

				v = symmetry_get_all_views(sym, part->view); // Get symmetry-related views
				
				r = nviews*irm*random();
				y = (int) r;	// y is a random integer in the range [0,n-1]

				for ( i=0; i<y ; i++, v=v->next) ;
				part->view = *v;
				part->view.next = NULL;
				
				kill_list((char *) v, sizeof(View));
				
				if ( verbose & VERB_FULL )
					cout << "View chosen #" << i << ": " << part->view << endl;
			}
		}
	} else {
		for ( field = project->field; field; field = field->next ) {
			for ( mg = field->mg; mg; mg = mg->next ) {
				for ( part = mg->part; part; part = part->next ) {
					if ( verbose & VERB_FULL )
						cout << "Original view vector: " << part->view << endl;

					v = symmetry_get_all_views(sym, part->view); // Get symmetry-related views
				
					r = nviews*irm*random();
					y = (int) r;	// y is a random integer in the range [0,n-1]

					for ( i=0; i<y ; i++, v=v->next) ;
					part->view = *v;
					part->view.next = NULL;
					
					kill_list((char *) v, sizeof(View));
					
					if ( verbose & VERB_FULL )
						cout << "View chosen #" << i << ": " << part->view << endl;
				}
			}
		}
	}

	return 0;
}

