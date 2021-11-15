/**
@file	View.cpp
@brief	View functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20160708
**/

#include "View.h"
#include "View2.h"
#include "linked_list.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, View v) {
	output.setf(ios::fixed, ios::floatfield);
	output << v.vector3() << tab << setw(10) << v.angle()*180.0/M_PI;
	return output;
}

//template <typename Type>
//ostream& operator<<(ostream& output, View2<Type> v) {
ostream& operator<<(ostream& output, View2<float> v) {
	output.setf(ios::fixed, ios::floatfield);
	output << v.vector3() << tab << setw(10) << v.angle()*180.0/M_PI;
	return output;
}

ostream& operator<<(ostream& output, View2<double> v) {
	output.setf(ios::fixed, ios::floatfield);
	output << v.vector3() << tab << setw(10) << v.angle()*180.0/M_PI;
	return output;
}

/**
@brief 	Displays a linked list of views.
@param 	*v				the linked list of views.
@return int				number of views.
**/
int			show_views(View* v)
{
	int			n(0);
	
	for ( ; v; v = v->next, n++ )
		cout << *v << endl;
	
	return n;
}

/**
@brief 	Copies views from a list into an array.
@param 	*v				the linked list of views.
@param 	n				number of views.
@return View*			view array.
**/
View*		view_array(View* v, long& n)
{
	long		i;
	
	n = count_list((char *)v);

	View*		varr = new View[n];

	for ( i=0; v; v=v->next, i++ ) varr[i] = *v;
	
	return varr;
}

/**
@brief 	Initializes a set of views tilted around the y axis.
@param 	ang_min			starting angle (radians).
@param 	ang_max			ending angle (radians).
@param 	ang_step		angular step size (radians).
@param 	axis			tilt axis angle (radians).
@return View* 			a set of 4-value views.

	A set of views is calculated corresponding to tilted views imaged
	during tomography. The tilt axis angle is taken as a counter-clockwise
	rotation from the x-axis.

**/
View*		tilt_views(double ang_min, double ang_max, double ang_step, double axis)
{
	if ( ang_max < ang_min ) swap(ang_max, ang_min);
	if ( ang_step < 0 ) ang_step = -ang_step;
	
	double			ang;
	int				n(0);
	View*			view = NULL;
	View*			v = NULL;
	
	for ( n=0, ang=ang_min; ang<ang_max+0.5*ang_step; ang+=ang_step, n++ ) {
		v = (View *) add_item((char **) &view, sizeof(View));
//		*v = view_from_tilt_and_axis(ang, axis);
		*v = View(ang, axis);
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Getting " << n << " views for tilting from " << ang_min*180.0/M_PI << 
			" to " << ang_max*180.0/M_PI << " by step " << ang_step*180.0/M_PI <<
			" around axis " << axis*180.0/M_PI << endl;
	
	return view;
}

/**
@brief 	Calculates a random view.
@return View			random view.

	A random seed should already have been generated.

**/
View		random_view()
{
	double			irm = 2.0/get_rand_max();
	View			v;
	
	v[0] = random()*irm - 1;
	v[1] = random()*irm - 1;
	v[2] = random()*irm - 1;
	v[3] = M_PI*(random()*irm - 1);
	v.normalize();
	
	return v;
}

/**
@brief 	Generates a random reslicing 3x3 rotation matrix.
@return View 			new view.

	The view represents any one or more 90 degree rotations,
	randomly chosen.

**/
View		view_random_reslice()
{
	double		irm = 1.0/get_rand_max();
	View		v;
	
	int			rv = (int) (5.999*irm*random());
	int			ra = (int) (3.999*irm*random());
	switch ( rv ) {
		case 0: v = View(1,0,0,ra*M_PI_2); break;
		case 1: v = View(0,1,0,ra*M_PI_2); break;
		case 2: v = View(0,0,1,ra*M_PI_2); break;
		case 3: v = View(-1,0,0,ra*M_PI_2); break;
		case 4: v = View(0,-1,0,ra*M_PI_2); break;
		case 5: v = View(0,0,-1,ra*M_PI_2); break;
	}
	
	return v;
}

/**
@brief 	Calculates a set of random views.
@param 	nviews			number of views.
@return View*			list of random views.
**/
View*		random_views(int nviews)
{
	if ( nviews < 1 ) {
		error_show("Error in random_views: No views defined!", __FILE__, __LINE__);
		return NULL;
	}
	
	random_seed();
	
	int				i;
	View*			view = NULL;
	View*			v = NULL;
	
	for ( i=0; i<nviews; i++ ) {
		v = (View *) add_item((char **) &v, sizeof(View));
		if ( !view ) view = v;
		*v = random_view();
	}
	
	return view;
}

/**
@brief 	Calculates a new view with a random error from the given view.
@param 	&v				given view (modified).
@param 	std				standard deviation of gaussian distribution.
@return double			random angle.

	The rotation between the given and new views is defined by the 
	gaussian distributed random angle and a random axis.
	A random seed should already have been generated.

**/
double		random_view_error(View& v, double std)
{
	double			angle = random_gaussian(0.0, std);

	Vector3<double>	axis = vector3_random(-1.0, 1.0);
	axis.normalize();

	Quaternion		q(v.quaternion()*Quaternion(axis, angle));
	
	v = View(q);
	
	return angle;
}

/**
@brief 	Generates a list of views within an angular distance from the input view.
@param 	theview			the input view.
@param 	theta_step		theta step size (radians).
@param 	phi_step		phi step size (radians).
@param 	alpha_step		alpha step size (radians).
@param 	view_angle_limit	angular distance limit from view vector (radians).
@param 	alpha_angle_limit	angular distance limit from view rotation angle (radians).
@return View*	 		a list of views.

	The list of views forms a 3D search grid in orientation space.

**/
View*		views_within_limits(View theview, double theta_step, double phi_step, 
					double alpha_step, double view_angle_limit, double alpha_angle_limit)
{
	double			theta, phi, alpha;
	double			theta_start = (theta_step)? -theta_step*floor(view_angle_limit/theta_step): 0;
	double			phi_start = (phi_step)? -phi_step*floor(view_angle_limit/phi_step): 0;
	double			alpha_start = theview.angle();
	double			alpha_end = alpha_start;
	
	if ( alpha_step > 0 && alpha_angle_limit > 0 ) {
		alpha = alpha_step*floor(alpha_angle_limit/alpha_step);
		alpha_start -= alpha;
		alpha_end += alpha;
	}
	
	if ( alpha_step <= 0 ) alpha_step = 1;
	
	Vector3<double>	ref_axis(0,0,1);
		
	theview.normalize();
	
	Vector3<double>	vv(theview[0], theview[1], theview[2]);
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	View*			view = NULL;	
	View*			v = NULL;
	Quaternion    	q1, q;
	Quaternion    	qv(0, theview[0], theview[1], theview[2]);
	
	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_within_limits: The View: " << theview << endl; 
		cout << "DEBUG views_within_limits: theta_step: "<< theta_step*180/M_PI 
			<< ", phi_step: " << phi_step*180/M_PI
			<< ", alpha_step: " << alpha_step*180/M_PI << endl;
		cout << "DEBUG views_within_limits: Axis1: " << axis1 << endl;
		cout << "DEBUG views_within_limits: Axis2: " << axis2 << endl;
	}
	
	if ( theta_step < 1e-6 ) theta_step = 1e6;
	if ( phi_step < 1e-6 ) phi_step = 1e6;
	if ( alpha_step < 1e-6 ) alpha_step = 1e6;
	for ( theta = theta_start; theta <= view_angle_limit; theta += theta_step ) {
//		q1 = quaternion_from_angle_and_axis3(theta, axis1);
		q1 = Quaternion(axis1, theta);
		for ( phi = phi_start; phi <= view_angle_limit; phi += phi_step ) {
//			q = quaternion_from_angle_and_axis3(phi, axis2);
			q = Quaternion(axis2, phi);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha = alpha_start; alpha <= alpha_end; alpha += alpha_step ) {
				v = (View *) add_item((char **) &v, sizeof(View));
				if ( !view ) view = v;
				*v = View(q[1], q[2], q[3], alpha);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_within_limits: " << *v << " - " << theview.angle(*v)*180/M_PI << endl;
			}
		}
	}
	
	return view;
}

list<View2<float>>	views_within_limits2(View2<float> theview, double theta_step, double phi_step,
					double alpha_step, double view_angle_limit, double alpha_angle_limit)
{
	double			theta, phi, alpha;
	double			theta_start = (theta_step)? -theta_step*floor(view_angle_limit/theta_step): 0;
	double			phi_start = (phi_step)? -phi_step*floor(view_angle_limit/phi_step): 0;
	double			alpha_start = theview.angle();
	double			alpha_end = alpha_start;
	
	if ( alpha_step > 0 && alpha_angle_limit > 0 ) {
		alpha = alpha_step*floor(alpha_angle_limit/alpha_step);
		alpha_start -= alpha;
		alpha_end += alpha;
	}
	
	if ( alpha_step <= 0 ) alpha_step = 1;
	
	Vector3<double>	ref_axis(0,0,1);
		
	theview.normalize();
	
	Vector3<double>	vv(theview[0], theview[1], theview[2]);
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	list<View2<float>>	v;
	Quaternion    		q1, q;
	Quaternion    		qv(0, theview[0], theview[1], theview[2]);
	
	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_within_limits: The View: " << theview << endl;
		cout << "DEBUG views_within_limits: theta_step: "<< theta_step*180/M_PI
			<< ", phi_step: " << phi_step*180/M_PI
			<< ", alpha_step: " << alpha_step*180/M_PI << endl;
		cout << "DEBUG views_within_limits: Axis1: " << axis1 << endl;
		cout << "DEBUG views_within_limits: Axis2: " << axis2 << endl;
	}
	
	if ( theta_step < 1e-6 ) theta_step = 1e6;
	if ( phi_step < 1e-6 ) phi_step = 1e6;
	if ( alpha_step < 1e-6 ) alpha_step = 1e6;
	for ( theta = theta_start; theta <= view_angle_limit; theta += theta_step ) {
		q1 = Quaternion(axis1, theta);
		for ( phi = phi_start; phi <= view_angle_limit; phi += phi_step ) {
			q = Quaternion(axis2, phi);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha = alpha_start; alpha <= alpha_end; alpha += alpha_step ) {
//				v = (View *) add_item((char **) &v, sizeof(View));
//				if ( !view ) view = v;
//				*v = View(q[1], q[2], q[3], alpha);
				v.push_back(View2<float>(q[1], q[2], q[3], alpha));
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_within_limits: " << v.back() << " - " << theview.angle(v.back())*180/M_PI << endl;
			}
		}
	}
	
	return v;
}


/**
@brief 	Generates a 3x3 grid of views around an input view.
@param 	theview			the input view.
@param 	alpha_step		angular step size around view vector (radians).
@return View*	 			a list of views.
**/
View*		views_for_refinement(View theview, double alpha_step)
{
	View*			view = NULL;	

	if ( alpha_step <= 0 ) {
		view = (View *) add_item((char **) &view, sizeof(View));
		*view = theview;
		return view;
	}
	
	double			alpha1, alpha2, alpha3;
	Vector3<double>	ref_axis(0,0,1);
		
	theview.normalize();
	
	Vector3<double>	vv(theview[0], theview[1], theview[2]);
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	View*			v = NULL;
	Quaternion    	q1, q;
	Quaternion    	qv(0, theview[0], theview[1], theview[2]);

	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_for_refinement: The View: " << theview << endl; 
		cout << "DEBUG views_for_refinement: alpha_step: "<< alpha_step*180/M_PI << endl;
		cout << "DEBUG views_for_refinement: Axis1: " << axis1 << endl;
		cout << "DEBUG views_for_refinement: Axis2: " << axis2 << endl;
	}
	
	for ( alpha1 = -alpha_step; alpha1 <= alpha_step; alpha1 += alpha_step ) {
		q1 = Quaternion(axis1, alpha1);
		for ( alpha2 = -alpha_step; alpha2 <= alpha_step; alpha2 += alpha_step ) {
			q = Quaternion(axis2, alpha2);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha3 = -alpha_step; alpha3 <= alpha_step; alpha3 += alpha_step ) {
				v = (View *) add_item((char **) &v, sizeof(View));
				if ( !view ) view = v;
				*v = View(q[1], q[2], q[3], theview.angle() + alpha3);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_for_refinement: " << *v << " - " << theview.angle(*v)*180/M_PI << endl;
			}
		}
	}
	
	return view;
}

/**
@brief 	Generates a 3x3 grid of views around an input view.
@param 	theview			the input view.
@param 	alpha_step1		angular step size for view vector in one direction (radians).
@param 	alpha_step2		angular step size for view vector in second direction (radians).
@param 	alpha_step3		angular step size around view vector (radians).
@param 	max_alpha3		maximum for step 3.
@return View*	 			a list of views.
**/
View*		views_for_refinement(View theview, double alpha_step1, double alpha_step2,
				double alpha_step3, double max_alpha3)
{
	View*			view = NULL;	

	if ( alpha_step1 <= 0 && alpha_step2 <= 0 && alpha_step3 <= 0 ) {
		view = (View *) add_item((char **) &view, sizeof(View));
		*view = theview;
		return view;
	}

	if ( alpha_step1 < 0 ) alpha_step1 = 0;
	if ( alpha_step2 < 0 ) alpha_step2 = 0;
	if ( alpha_step3 <= 0 ) alpha_step3 = max_alpha3 = 0;
	
	double			tol(M_PI/1800.0);	// Tolerance for small angles: 0.1 degrees
	double			alpha1, alpha2, alpha3;
	Vector3<double>	ref_axis(0,0,1);
		
	theview.normalize();
	
	Vector3<double>	vv(theview[0], theview[1], theview[2]);
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	View*			v = NULL;
	Quaternion    	q1, q;
	Quaternion    	qv(0, theview[0], theview[1], theview[2]);

	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_for_refinement: The View: " << theview << endl; 
		cout << "DEBUG views_for_refinement: alpha_steps: " << alpha_step1*180/M_PI
			<< " " << alpha_step2*180/M_PI << " " << alpha_step3*180/M_PI << endl;
		cout << "DEBUG views_for_refinement: Axis1: " << axis1 << endl;
		cout << "DEBUG views_for_refinement: Axis2: " << axis2 << endl;
	}
	
	for ( alpha1 = -alpha_step1; alpha1 <= alpha_step1; alpha1 += alpha_step1 ) {
//		q1 = quaternion_from_angle_and_axis3(alpha1, axis1);
		q1 = Quaternion(axis1, alpha1);
		for ( alpha2 = -alpha_step2; alpha2 <= alpha_step2; alpha2 += alpha_step2 ) {
//			q = quaternion_from_angle_and_axis3(alpha2, axis2);
			q = Quaternion(axis2, alpha2);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha3 = -max_alpha3; alpha3 <= max_alpha3; alpha3 += alpha_step3 ) {
				v = (View *) add_item((char **) &v, sizeof(View));
				if ( !view ) view = v;
				*v = View(q[1], q[2], q[3], theview.angle() + alpha3);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_for_refinement: " << *v << " - " << theview.angle(*v)*180/M_PI << endl;
				if ( alpha_step3 < tol ) alpha3 += 1;
			}
			if ( alpha_step2 < tol ) alpha2 += 1;
		}
		if ( alpha_step1 < tol ) alpha1 += 1;
	}
	
	return view;
}

/**
@brief 	Expands each view to several views with different rotation angles.
@param 	*views			view list.
@param 	amin			minimum angle.
@param 	amax			maximum angle.
@param 	astep			angular step.
@return View*			new view list.

	The new angles are added to the existing angles of the view.

**/
View*		view_list_expand_angles(View* views, double amin, double amax, double astep)
{
	double		a;
	View*		newviews = NULL;
	View*		v = NULL;
	View*		vn = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_list_expand_angles: amin=" << amin << " amax=" << amax << " astep=" << astep << endl;
	
	for ( v = views; v; v = v->next ) {
		if ( verbose & VERB_FULL )
			cout << "Old view: " << v << ": " << *v << endl;
		for ( a = amin; a <= amax; a += astep ) {
			vn = (View *) add_item((char **) &vn, sizeof(View));
			if ( !newviews ) newviews = vn;
			*vn = *v;
			(*vn)[3] += a;
			if ( verbose & VERB_FULL )
				cout << vn << *vn << endl;
		}
	}
	
	if ( verbose & VERB_FULL )
		show_views(newviews);
	
	return newviews;
}

/**
@brief 	Replaces a view list with a subset of it.
@param 	**view_list		view list.
@param 	start			offset of first view of subset.
@param 	size			number of views in subset.
@return int				number of views selected, <0 on error.

	The old view list is destroyed.

**/
int			view_list_subset(View** view_list, int start, int size)
{
	if ( size < 1 ) {
		cerr << "Warning: No subset selected - view list remains the same!" << endl;
		return 0;
	}
	
	int			i;
	View*		view;
	View*		view_next = NULL;
	
	if ( start > 0 ) {
		for ( i=0, view = *view_list; view && i<start-1; view = view->next, i++ ) ;
		if ( view ) {
			view_next = view->next;
			view->next = NULL;
		}
	
		kill_list((char *) *view_list, sizeof(View));
	
		*view_list = view_next;
	}

	if ( !*view_list ) {
		cerr << "Error: No views selected!" << endl;
		return -1;
	}
		
	for ( i=0, view = *view_list; view && i<size-1; view = view->next, i++ ) ;
	
	if ( view ) {
		i++;
		view_next = view->next;
		view->next = NULL;
		kill_list((char *) view_next, sizeof(View));
	}
	
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of views selected:       " << i << endl;
	
	return i;
}


