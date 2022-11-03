/**
@file	symmetry.cpp
@brief	General symmetry functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20210116
**/

#include "symmetry.h"
#include "random_numbers.h"
#include "Euler.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Symmetry constructor from a symmetry identifier.
@param 	sym			string containing point group identifier.

	The point group symmetries are identified by the following strings:
		C<n>		cyclic point group of order n.
		D<n>			dihedral point group of order n.
		T			tetrahedral point group.
		O			octahedral/cubic point group.
		I			icosahedral/dodecahedral point group.
		H<r>,<a>,<d>	helical symmetry with rise r, rise angle a and dyad d (1/2).
	For the higher symmetries the following adjustments are available:
		T-3			no three-fold operator.
		O-2			no two-fold operator.
		O-3			no three-fold operator.
		O-4			no four-fold operator.
		I-2			no two-fold operator.
		I-3			no three-fold operator.
		I-5			no five-fold operator.
		I90			90 degrees rotated around z-axis.
		I90-3			90 degrees rotated around z-axis and no three-fold operator.
	If the point group string is empty, the default is C1 (asymmetric).

**/
Bsymmetry::Bsymmetry(Bstring sym)
{
	int 			theorder, helix_dyad(1);
	double			helix_rise = 1, helix_angle = M_PI;
	Matrix3 		mat(0,-1,0,1,0,0,0,0,1); // 90 degree rotation for I90
	
	clean_symstring(sym);
	
 	if ( verbose & VERB_FULL )
		cout << endl << "Getting symmetry operators for " << lbl << endl;
	
	if ( lbl[0] == 'C' ) {									// Cyclic point groups
		if ( lbl == "Cs" ) theorder = 1;					// Reflection
		else theorder = lbl.substr(1,10).integer();
		pnt = 100 + theorder;
		op.push_back(Bsymop(0,0,1,theorder,M_PI*2.0L/theorder));
	} else if ( lbl[0] == 'D' ) {							// Dihedral point groups
		theorder = lbl.substr(1,10).integer();
		pnt = 200 + theorder;
		op.push_back(Bsymop(0,0,1,theorder,M_PI*2.0L/theorder));
		op.push_back(Bsymop(1,0,0,2,M_PI));					// Rotate around 2-fold x-axis
	} else if ( lbl[0] == 'T' ) {							// Tetrahedral point group
		pnt = 320;
		if ( !lbl.contains("-3") )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		op.push_back(Bsymop(0,0,1,2,M_PI));					// Rotate around 2-fold z-axis
		op.push_back(Bsymop(1,0,0,2,M_PI));					// Rotate around 2-fold x-axis
	} else if ( lbl[0] == 'O' ) {							// Octahedral point group
		pnt = 432;
		if ( !lbl.contains("-3") )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		if ( lbl.contains("-4") ) {
			op.push_back(Bsymop(1,-1,0,2,M_PI));			// Rotate around 2-fold axis
		} else if ( lbl.contains("-2") ) {
			op.push_back(Bsymop(0,0,1,4,M_PI_2));			// Rotate around 4-fold z-axis
		} else {
			op.push_back(Bsymop(0,0,1,4,M_PI_2));			// Rotate around 4-fold z-axis
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Rotate around 2-fold x-axis
		}
	} else if ( lbl[0] == 'I' ) {							// Icosahedral point group
		pnt = 532;
		if ( !lbl.contains("-3") )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		if ( lbl.contains("-5") ) {
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Rotate around 2-fold x-axis
		} else if ( lbl.contains("-2") ) {
			op.push_back(Bsymop(1,-1,0,2,M_PI));			// Rotate around 2-fold axis
		} else {
			op.push_back(Bsymop(1,1.0L/GOLDEN,GOLDEN,2,M_PI));	// Rotate around 2-fold
		}
		if ( !lbl.contains("-5") ) {
			op.push_back(Bsymop(1.0L/GOLDEN,1,0,5,M_PI*2.0L/5.0));	// Rotate around 5-fold
		}
		if ( !lbl.contains("-2") ) {
			op.push_back(Bsymop(0,0,1,2,M_PI));				// Rotate around 2-fold z-axis
		}
		if ( lbl.contains("90") )
			transform(mat);
	} else if ( lbl[0] == 'H' ) {							// Helical symmetry
		sscanf(lbl.c_str(), "H%lf,%lf,%d", &helix_rise, &helix_angle, &helix_dyad);
		pnt = 600;
		op.push_back(Bsymop(0,0,1,1,helix_angle*M_PI/180.0L));	// Helical axis
		op[0].shift(helix_rise);							// Rise in angstrom
		if ( helix_dyad == 2 )
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Dyad axis
	} else {
		error_show("Error in Bsymmetry::Bsymmetry", __FILE__, __LINE__);
		cerr << "This is not a valid symmetry designation: " << lbl << endl;
		return;
	}
	
	for ( auto it = op.begin(); it != op.end(); ++it )
		it->normalize();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bsymmetry::Bsymmetry: lbl=" << lbl << endl;
}

/**
@brief 	Returns the label for helical symmetry.
@param 	helix_rise		helical rise (angstroms).
@param	helix_angle		helical rotation angle (radians).
@param	dyad_axis		presence of dyad axis (1/2).
@param	cyclic			cyclic symmetry.
@param	seam_shift		fractional shift for seamed helices.
@return string			label.

	Thge symmetry order is defined as the product of all the individual
	orders of the symmetry operations, or alternatively, the number of views.

**/
string		symmetry_helical_label(double helix_rise, double helix_angle,
								   int dyad_axis, int cyclic, double seam_shift)
{
	string			label = "H" + to_string(helix_rise);
	label += "," + to_string(helix_angle*180.0/M_PI);
	label += "," + to_string(dyad_axis);
	label += "," + to_string(cyclic);
	label += "," + to_string(seam_shift);

	return label;
}

/**
@brief 	Returns an asymmetric unit reference point.
@param 	&sym		symmetry structure.
@return View			reference view.
**/
View		view_symmetry_reference(Bsymmetry& sym)
{
	View			ref;
 
	if ( sym.point() < 200 ) {
		ref[0] = 1;
		ref[2] = 0;
	} else if ( sym.point() < 300 ) {
		ref[0] = sin(M_PI/4);
		ref[2] = cos(M_PI/4);
	} else if ( sym.point() == 320 ) {
		ref[0] = sin(M_PI/4);
		ref[2] = cos(M_PI/4);
	} else if ( sym.point() == 432 ) {
		ref[0] = sin(M_PI/8);
		ref[2] = cos(M_PI/8);
	} else if ( sym.point() == 532 ) {
		ref[0] = sin(M_PI/18);
		ref[2] = cos(M_PI/18);
	} else if ( sym.point() == 600 ) {	// What is it for helical symmetry?
	}

	return ref;
}

/**
@brief 	Rotation matrix to orient a symmetry axis on the z axis.
@param 	&sym		symmetry structure.
@param 	axis		desired symmetry axis order.
@param 	axis_flag	view modifier.
@return Matrix3		new rotation matrix.
**/
Matrix3		symmetry_rotate_to_axis(Bsymmetry& sym, long axis, long axis_flag)
{
	double			sqrt2(sqrt(2));
	Matrix3			mat(1), mat2(1);

	if ( axis_flag ) mat2 = Matrix3(Vector3<double>(0,0,1), M_PI_2);
	
	if ( sym.point() < 200 ) {							// Cyclic
	} else if ( sym.point() < 300 ) {					// Dihedral
		if ( axis_flag ) mat2 = Matrix3(Vector3<double>(0,0,1), M_PI/sym[0].order());
		if ( axis > 2 ) mat = Matrix3(Vector3<double>(0,1,0), M_PI_2);
	} else if ( sym.point() < 400 ) {					// Tetrahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(-1/sqrt2,1/sqrt2,0), atan(sqrt2));
	} else if ( sym.point() < 500 ) {					// Octahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(-1/sqrt2,1/sqrt2,0), atan(sqrt2));
		if ( axis == 2 ) mat = Matrix3(Vector3<double>(0,1,0), M_PI/4);
	} else {											// Icosahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(0,1,0), atan(1/(GOLDEN*sqrt(3))));
		if ( axis == 5 ) mat = Matrix3(Vector3<double>(1,0,0), atan(1/GOLDEN));
	}	
		
	mat = mat2*mat;

	return mat;
}

/**
@brief 	Get all symmetry axes.
@param 	&sym				symmetry structure.
@return vector<Vector3<double>>	array of axes.
**/
vector<Vector3<double>>	symmetry_get_axes(Bsymmetry& sym)
{
	long			i, j, k, m, ns;
	Vector3<double>	axis1;
	Matrix3			mat;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG symmetry_get_axes: label=" << sym.label() << endl;
	
	if ( sym.point() < 300 ) {					// Cyclic and dihedral
		axis1[2] = 1;
	} else if ( sym.point() < 400 ) {					// Tetrahedral
		axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
	} else if ( sym.point() < 500 ) {					// Octahedral
		if ( sym.label().contains("-4") ) {
			axis1[2] = 1;
		} else if ( sym.label().contains("-3") ) {
			axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
		} else {
			cerr << "Error: The symmetry designation must be either O-3 or O-4!" << endl;
			bexit(-1);
		}
	} else if ( sym.point() < 600 ) {						// Icosahedral
		if ( sym.label().contains("-5") ) {
			axis1[1] = 1/sqrt(2+GOLDEN);
			axis1[2] = 1/sqrt(3-GOLDEN);
		} else if ( sym.label().contains("-3") ) {
			axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
		} else {
			cerr << "Error: The symmetry designation must be either I-3 or I-5!" << endl;
			bexit(-1);
		}
	}
	
	vector<Vector3<double>>	axis;
	
	axis.push_back(axis1);

	for ( i=0, k=1, ns=1; i<sym.operations(); i++ ) {
		for ( j=1; j<sym[i].order(); j++ ) {
			mat = Matrix3(sym[i].axis(), j*TWOPI/sym[i].order());
			for ( m=0; m<ns; m++, k++ ) axis.push_back(mat * axis[m]);
		}
		ns = k;
	}
	
	if ( verbose & VERB_FULL )
		for ( i=0; i<axis.size(); i++ )
			cout << "Axis[" << i+1 << "]: " << axis[i] << endl;
	
	return axis;
}

/**
@brief 	Get all symmetry-related views of one given view.
@param 	&sym		symmetry structure.
@param 	asu_view	asymmetric unit vector and rotation angle.
@return View*			linked list of views.

	The number of views generated for a point group symmetry is
	calculated as the product of the order fields in the symmetry
	structure.

**/
View*		symmetry_get_all_views(Bsymmetry& sym, View asu_view)
{
	int				i, j, k;
	double			angle;
	Quaternion		q, qv; 
	
	View*			view = NULL;
	View*			v;
	View*			vn = (View *) add_item((char **) &view, sizeof(View));
	*vn = asu_view;
	vn->normalize();
	vn->next = NULL;
	
	if ( verbose & VERB_FULL ) {
		cout << "Getting all the symmetric views:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "View:                           " << asu_view << endl;
	}
	
	int				nview(1);
	for ( i=0; i<sym.operations(); i++ ) {
		for ( j=1; j<sym[i].order(); j++ ) {
			angle = j*TWOPI*1.0L/sym[i].order();
			q = Quaternion(sym[i].axis(), angle);
			for ( k=0, v=view; k<nview; k++, v=v->next ) {
				qv = v->quaternion();
				qv = q * qv;
//				vn = (View *) add_item((char **) &view, sizeof(View));
				vn = (View *) add_item((char **) &vn, sizeof(View));
				*vn = View(qv);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG symmetry_get_all_views: " << vn << endl;
			}
		}
		nview *= sym[i].order();
	}
	
	if ( verbose & VERB_FULL )
		for ( i=1, v=view; v; v=v->next, i++ )
			cout << "View " << i << ": " << tab << *v << endl;
	
	return view;
}

View*		symmetry_get_all_views(Bsymmetry& sym, View* views)
{
	View*		v;
	View*		vt;
	View*		vc = NULL;
	View*		vn = NULL;
	
	for ( v=views; v; v=v->next ) {
		vt = symmetry_get_all_views(sym, *v);
		if ( vc ) vc->next = vt;
		else vn = vc = vt;
		while ( vc->next ) vc = vc->next;
	}
	
	return vn;
}

/**
@brief 	Get all symmetry-related views of one given view.
@param 	&sym			symmetry structure.
@return vector<Matrix3>		array of matrices.

	The number of views generated for a point group symmetry is
	calculated as the product of the order fields in the symmetry
	structure.

**/
vector<Matrix3>	symmetry_get_all_matrices(Bsymmetry& sym)
{
	vector<Matrix3>	mat = sym.matrices();
	
	if ( verbose & VERB_FULL ) {
		cout << "Getting all the symmetric matrices:" << endl;
		cout << "Symmetry:                       " << sym.label() << endl;
		cout << "Number of matrices:             " << mat.size() << endl;
	}
	
	if ( verbose & VERB_FULL ) {
		for ( int i=1; i<=mat.size(); ++i ) {
			cout << "Matrix " << i << ":" << endl;
			matrix3_show_hp(mat[i-1]);
		}
	}
	
	return mat;
}

/**
@brief 	Initializes a well-distributed set of views in an asymmetric unit.
@param 	&sym		symmetry structure.
@param 	theta_step	angular step size from primary symmetry axis (radians).
@param 	phi_step	angular step size around primary symmetry axis (radians).
@param 	flag		flag: 0=half, 1=full, 2=no in-plane.
@return View* 		a linked list of views.

	A set of views is calculated with tesselation within each asymmetric
	unit such that the views are well-distributed.
 	Flag bits:
		1: both halves of the asymmetric unit are covered.
 		2: no in-plane rotations are applied.

**/
View* 		asymmetric_unit_views(Bsymmetry& sym, double theta_step, double phi_step, int flag)
{
	if ( theta_step < M_PI/1800 ) theta_step = M_PI/180;
	if ( phi_step < M_PI/1800 ) phi_step = M_PI/180;

	char			full = flag & 1;
	char			inplanerot = !(flag & 2);
	int 			i, j, n = 0, ntheta, nphi, nrphi, istart;
	double			theta, phi, max_theta, theta_start, phi_start, phi_end;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG asymmetric_unit_views: flag=" << flag << " full=" << full << endl;

	max_theta = M_PI_2 + 1e-6;			// Most symmetries limited to upper half
	ntheta = (int) (max_theta/theta_step + 0.5);
	nphi = (int) (M_PI*2.0/(sym[0].order()*phi_step) + 0.5);
	
	if ( sym.point() < 200 ) {
		if ( full ) {
			max_theta = M_PI + 1e-6;
			ntheta += ntheta;
		}
	} else if ( sym.point() == 432 ) {
		max_theta = M_PI_4 + 1e-6;
		ntheta = (int) (max_theta/theta_step + 0.5);
	} else if ( sym.point() == 532 ) {
		max_theta = atan(1/(GOLDEN+1)) + 1e-6;
		ntheta = (int) (max_theta/theta_step + 2.5);
		nphi = (int) (atan(1/GOLDEN)/phi_step + 1.5);
	} else if ( sym.point() == 600 ) {
		ntheta = 1;
		nphi = (int) (M_PI*2.0/phi_step + 1);
	}
	
	if ( ntheta < 1 ) ntheta = 1;
	if ( nphi < 1 ) nphi = 1;
		
	if ( verbose & (VERB_PROCESS | VERB_LABEL) ) {
		if ( full )
			cout << "Getting all the asymmetric unit views for symmetry " << sym.label() << endl;
		else
			cout << "Getting half the asymmetric unit views for symmetry " << sym.label() << endl;
		cout << "Theta and phi step sizes:       " << 
				theta_step*180/M_PI << " " << phi_step*180/M_PI << " degrees" << endl;
		cout << "Theta and phi steps:            " << ntheta << " " << nphi << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG asymmetric_unit_views: Views allocated = " << 2*ntheta*nphi << endl;
	
	View*			view = NULL;	
	View*			v = NULL;
	
	if ( sym.point() < 500 ) {						// Top view for all but icosahedral and helical symmetry
		v = (View *) add_item((char **) &view, sizeof(View));
		(*v)[2] = 1;
		n = 1;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG asymmetric_unit_views: Adding the top view" << endl;
	}
	
	if ( sym.point() > 100 && sym.point() < 200 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = istart = 0;
			if ( nrphi ) {
				phi = M_PI*1.0/(sym[0].order()*nrphi);
				istart = -nrphi;
				if ( sym.point() < 102 ) nrphi -= 1;
			}
			for ( j=istart; j<=nrphi; j++ ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				*v = View(sin(theta)*cos(j*phi), sin(theta)*sin(j*phi), cos(theta), 0);
				if ( inplanerot ) (*v)[3] = j*phi;
				n++;
			}
		}
		if ( full && v->z() > -0.999999 ) {
			v = (View *) add_item((char **) &view, sizeof(View));
			(*v)[2] = -1;
			n++;
		}
	} else if ( sym.point() > 200 && sym.point() < 300 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = 0;
			if ( nrphi ) phi = M_PI*1.0/(sym[0].order()*nrphi);
			istart = 0;
			if ( full ) istart = -nrphi;
			for ( j=istart; j<=nrphi; j++ ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				*v = View(sin(theta)*cos(j*phi), sin(theta)*sin(j*phi), cos(theta), 0);
				if ( inplanerot ) (*v)[3] = j*phi;
				n++;
			}
		}
	} else if ( sym.point() == 320 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			phi_end = theta;
			if ( theta > M_PI_4 ) phi_end = M_PI_2 - theta;
			if ( full ) phi_start = -theta;
			if ( full && theta > M_PI_4 ) phi_start = theta - M_PI_2;
			for ( phi=phi_start; phi<=phi_end; phi+=phi_step ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				*v = View(sin(theta), sin(phi), cos(theta), 0);
				n++;
			}
		}
	} else if ( sym.point() == 432 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			if ( full ) phi_start = -theta;
			for ( phi=phi_start; phi<=theta; phi+=phi_step ) {
				v = (View *) add_item((char **) &view, sizeof(View));
				*v = View(sin(theta), sin(phi), cos(theta), 0);
				n++;
			}
		}
	} else if ( sym.point() == 532 ) {
		if ( sym.label().contains("I90") ) {
			theta_start = 0;
			if ( full ) theta_start = -ntheta*theta_step;
			for ( theta=theta_start; theta<=max_theta; theta+=theta_step ) {
				for ( phi=0; phi<=GOLDEN*(max_theta - fabs(theta)); phi+=phi_step ) {
					v = (View *) add_item((char **) &view, sizeof(View));
					*v = View(tan(phi), tan(theta), cos(theta), 0);
					n++;
				}
			}
		} else {
			for ( theta=0; theta<=max_theta; theta+=theta_step ) {
				nrphi = (int) (GOLDEN*(max_theta - theta)/phi_step);
				phi_start = 0;
				if ( full ) phi_start = -nrphi*phi_step;
				for ( phi=phi_start; phi<=GOLDEN*(max_theta - theta); phi+=phi_step ) {
					v = (View *) add_item((char **) &view, sizeof(View));
					*v = View(tan(theta), tan(phi), cos(theta), 0);
					n++;
				}
			}
		}
	} else if ( sym.point() == 600 ) {					// Helical symmetry
		for ( phi = 0; phi < TWOPI-0.001; phi += phi_step ) {
			v = (View *) add_item((char **) &view, sizeof(View));
			*v = View(cos(phi), sin(phi), 0, 0);
			if ( inplanerot ) (*v)[3] = phi;
			n++;
		}
	} else {
		cerr << "Warning: Symmetry type " << sym.point() << " not supported!" << endl;
	}

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << "Views generated:                " << n << endl << endl;

	for ( v=view; v; v = v->next ) v->normalize();
	
	if ( verbose & VERB_FULL ) {
		cout << "View\tx\ty\tz\ta" << endl;
		for ( v=view, i=1; v; v = v->next, i++ )
			cout << i << tab << *v << endl;
		cout << endl;
	}
	
	return view;
}

View* 		asymmetric_unit_views(Bsymmetry& sym, double theta_step, double phi_step, double alpha_step, int flag)
{
	View*		view = asymmetric_unit_views(sym, theta_step, phi_step, flag);
	
	View*		view2 = view_list_expand_angles(view, -M_PI, M_PI - alpha_step/2, alpha_step);
	
	kill_list((char *) view, sizeof(View));

	return view2;
}

/**
@brief 	Initializes a set of views around the z-axis for helical projection.
@param 	&sym		symmetry structure.
@param 	side_ang	starting angle (radians).
@param 	theta_step	angular step size perpendicular to equator (radians).
@param 	phi_step	angular step size around equator (radians).
@return View* 		a set of 4-value views.

	A set of views is calculated corresponding to views around the z-axis
	including some tilting to account for oblique views.

**/
View*		side_views(Bsymmetry& sym, double side_ang, double theta_step, double phi_step)
{
	if ( side_ang < 0 ) side_ang = -side_ang;
	if ( theta_step <= 0 ) {
		side_ang = 0;
		theta_step = M_PI/180;
	} else {
		side_ang = theta_step*floor(side_ang/theta_step);
	}
	if ( phi_step <= 0 ) phi_step = M_PI_2;	
	
	int				order = sym[0].order();
	if ( order < 1 ) order = 1;
	
	double			phi, theta;
	int				n;
	Euler			euler;
	View*			view = NULL;
	View*			v = NULL;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Getting all the asymmetric unit side views for symmetry " << sym.label() << endl;
		cout << "Step size around equator:       " << phi_step*180.0/M_PI << " degrees" << endl;
		if ( side_ang ) {
			cout << "Deviation from equator:         " << side_ang*180.0/M_PI << " degrees" << endl;
			cout << "Step size from equator:         " << theta_step*180.0/M_PI << " degrees" << endl;
		}
	}
	
	for ( n = 0, phi = -M_PI*1.0L/order; phi < M_PI*1.0L/order; phi += phi_step ) {
		for ( theta = M_PI_2 - side_ang; theta <= M_PI_2 + side_ang; theta += theta_step ) {
			v = (View *) add_item((char **) &view, sizeof(View));
			euler = Euler(0, theta, phi);
			*v = euler.view();
			n++;
//			cout << "theta=" << theta*180.0/M_PI << " phi=" << phi*180.0/M_PI << endl;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of views:                " << n << endl << endl;
	
	return view;
}

View*		side_views(Bsymmetry& sym, double side_ang, double theta_step, double phi_step, double alpha_step)
{
	View*		view = side_views(sym, side_ang, theta_step, phi_step);
	
	View*		view2 = view_list_expand_angles(view, -M_PI, M_PI - alpha_step/2, alpha_step);
	
	kill_list((char *) view, sizeof(View));

	return view2;
}

/**
@brief 	Change the views to those in the asymmetric unit.
@param 	&sym		symmetry structure.
@param 	*view		list of views (replaced).
@return int			0.

	The view is replaced with the one in the standard asymmetric unit.

**/
int 		change_views_to_asymmetric_unit(Bsymmetry& sym, View* view)
{
	if ( sym.point() < 102 ) return 0;
	
	View*			v;

	for ( v=view; v; v=v->next )
		*v = find_asymmetric_unit_view(sym, *v);
	
	return 0;
}

/**
@brief 	Finds the corresponding view in the asymmetric unit.
@param 	&sym		symmetry structure.
@param 	theview		view.
@return View			the asymmetric unit view.

	The asymmetric unit view is found and the the new view with the 
	link from the old view is returned.

**/
View 		find_asymmetric_unit_view(Bsymmetry& sym, View theview)
{
	if ( sym.point() < 102 ) return theview;
	if ( sym.point() >= 600 ) return theview;
	
	View*			v, bv(theview), tv;
	View*			view = symmetry_get_all_views(sym, theview);
	
	int				found(0);
	double			tol(1e-10);
	double			lim = tan(M_PI*1.0L/sym[0].order());
	
	for ( v=view; v && !found; v=v->next )
			if ( ( v->x() - tol >= 0 ) || ( v->x() + tol >= 0 && v->y() >= 0 ) ) {
		tv.x(0); tv.y(0); tv.z(0);
		if ( sym.point() < 200 ) {
			if ( sym[0].order() == 2 ) {
				tv = *v;
			} else {
				if ( fabs(v->y()) <= v->x()*lim + tol )
					tv = *v;
			}
		} else if ( sym.point() > 200 && sym.point() < 300 ) {
			if ( v->z() + tol >= 0 ) {
				if ( sym[0].order() == 2 ) {
					tv = *v;
				} else {
					if ( fabs(v->y()) <= v->x()*lim + tol ) tv = *v;
				}
			}
		} else if ( sym.point() == 320 ) {
			if ( fabs(v->y()) <= v->x() + tol && fabs(v->y()) <= v->z() + tol )
				tv = *v;
		} else if ( sym.point() == 432 ) {
			if ( fabs(v->y()) <= v->x() + tol && v->x() <= v->z() + tol )
				tv = *v;
		} else if ( sym.point() == 532 ) {
			if ( sym.label().contains("I90") ) {
				if ( v->z() + tol >= 0 && fabs(v->y()) <= (v->z()/GOLDEN - v->x())/GOLDEN + tol )
					tv = *v;
			} else {
				if ( v->z() + tol >= 0 && fabs(v->y()) <= v->z()/GOLDEN - v->x()*GOLDEN + tol )
					tv = *v;
			}
		} else if ( sym.point() >= 600 ) {	// What is it for helical symmetry?
			tv = theview;
		}
		if ( tv.vector_size() > 0.9 ) {
			bv = tv;
			found++;
		}
	}
	
	if ( !found ) {
		cerr << "ASU view not found: " << theview << endl;
//		for ( v=view; v; v=v->next )
//			cerr << *v << endl;
	} else {
		if ( verbose & VERB_STATS ) {
			cout << "Old view:\t" << theview << endl;
			cout << "New view:\t" << bv << endl;
		}
	}

	kill_list((char *) view, sizeof(View));
	
	bv.next = theview.next;
	
	return bv;
}

/**
@brief 	Finds the closest symmetric match between two views.
@param 	&sym		symmetry structure.
@param 	view_ref	reference view.
@param 	view		test view.
@return View			matched symmetric version of test view.

	A list of symmetry-related views of the test view is searched
	for the closest to the reference view.
	The matched view is returned.

**/
View 		find_closest_symmetric_view(Bsymmetry& sym, View view_ref, View view)
{
	if ( sym.point() < 102 ) return view;
	if ( sym.point() >= 600 ) return view;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG find_closest_symmetric_view: " << sym.label() << endl;
	
	double			r, rb(1e30);
	View*			v, bv(view);
	View*			viewsym = symmetry_get_all_views(sym, view);
	
	for ( v=viewsym; v; v=v->next ) {
//		r = v->residual(view_ref);
		r = v->angle(view_ref);
		if ( rb > r ) {
			rb = r;
			bv = *v;
		}
	}
	
	kill_list((char *) viewsym, sizeof(View));
	
	bv.next = view.next;
	
	return bv;
}

/**
@brief 	Returns a reference view for each asymmetric unit.
@param 	&sym		symmetry structure.
@return View*			reference views.
**/
View* 		reference_asymmetric_unit_views(Bsymmetry& sym)
{
	View			ref;
	
	if ( sym.point() < 200 ) {
		ref[0] = 1;
		ref[2] = 0;
	} else if ( sym.point() < 300 ) {
		ref[0] = ref[2] = 1;
	} else if ( sym.point() == 320 ) {
		ref[0] = ref[1] = cos(M_PI/8.0);
		ref[2] = sin(M_PI/8.0);
	} else if ( sym.point() == 432 ) {
		ref[0] = cos(M_PI/3.0);
		ref[2] = sin(M_PI/3.0);
	} else if ( sym.point() == 532 ) {
		if ( sym.label().contains("I90") ) {
			ref[1] = 1;
			ref[2] = GOLDEN + 1;
		} else {
			ref[0] = 1;
			ref[2] = GOLDEN + 1;
		}
	} else if ( sym.point() == 600 ) {	// What is it for helical symmetry?
	}

	ref.normalize();

	View*			rv = symmetry_get_all_views(sym, ref);

	return rv;
}

/**
@brief 	Returns a randomly picked symmetry view.
@param	&asu_view	asymmetric unit view.
@param 	&sym		symmetry structure.
@return View			picked view.
**/
View 		random_symmetric_view(View& asu_view, Bsymmetry& sym)
{
	View*			v;
	View*			view = symmetry_get_all_views(sym, asu_view);
	
	long			n(sym.order());
	
	long			j, i(n*random()*0.99999L/get_rand_max());

	for ( j=0, v=view; v && j<i; v=v->next, j++ ) ;
	
	View			vr(*v);
	
	kill_list((char *) view, sizeof(View));
	
	return vr;
}

/**
@brief 	Returns an asymmetric unit index number.
@param 	theview		view to test.
@param 	&sym		symmetry structure.
@return int			view number.
**/
int 		test_asymmetric_unit_view(View theview, Bsymmetry& sym)
{
	int				i, n = 0;
	double			a, amin = TWOPI;

	View*			rv = reference_asymmetric_unit_views(sym);
	View*			v;
	
	for ( i=0, v=rv; v; v=v->next, i++ ) {
		a = theview.angle(*v);
		if ( amin > a ) {
			amin = a;
			n = i;
		}
	}
	
	kill_list((char *) rv, sizeof(View));
		
	return n;
}

/**
@brief 	Show symmetry matrices.
@param 	&sym		symmetry structure.
@return int 			number of symmetry matrices.
**/
int			sym_show_matrices(Bsymmetry& sym)
{
	vector<Matrix3>	mat = sym.matrices();
	int				n(mat.size());

	cout << "\nSymmetry matrices (" << n << "):" << endl;
	for ( int i=0; i<n; ++i ) {
		cout << "Matrix " << i+1 << ":" << endl;
		matrix3_show_hp(mat[i]);
	}
	cout << endl;
	
	return n;
}

/**
@brief 	Show symmetry matrices associated with each symmetry operator.
@param 	&sym		symmetry structure.
@return int 			number of symmetry matrices.
**/
int			sym_show_operational_matrices(Bsymmetry& sym)
{
	int			i;

	for ( i=0; i<sym.operations(); ++i ) {
		cout << "Operation " << i+1 << ":" << endl;
		cout << "Axis and angle: " << sym[i].axis() << " " << 360.0/sym[i].order() << endl;
		cout << sym[i].matrix() << endl;
	}
	cout << endl;
	
	return sym.operations();
}

/**
@brief 	Show PDB format symmetry matrices associated with each symmetry operator.
@param 	&sym		symmetry structure.
@return int 			number of symmetry matrices.
**/
int			sym_show_pdb_matrices(Bsymmetry& sym)
{
	int			i, j, k;
	double		t(0);

	vector<Matrix3>	mat = sym.matrices();
	long			nsym = mat.size();

	for ( i=0; i<nsym; i++ ) {
		for ( j=0; j<3; j++ ) {
			cout << "BIOMT" << j+1 << setw(4) << i+1 << setprecision(6);
			for ( k=0; k<3; k++ )
					cout << setw(10) << mat[i][j][k];
			cout << setprecision(6) << setw(15) << t << endl;
		}
	}
	
	return nsym;
}


