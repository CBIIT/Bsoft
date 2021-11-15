/**
@file	marker.cpp
@brief	Functions for dealing with markers and linkers
@author Bernard Heymann
@date	Created: 20020619
@date	Modified: 20210112
**/

#include "rwmodel.h"
#include "marker.h"
#include "matrix_linear.h"
#include "Matrix3.h"
#include "Matrix.h"
#include "ps_marker.h"
#include "qsort_functions.h"
#include "linked_list.h"
#include "utilities.h"
#include "timer.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Local function prototypes
int 		markers_findOneToOne(Bmarker* refset, Bmarker* applyset);
double*		markers_RoughMatch(Bmarker* setA, Bmarker* setB, double dist_cutoff);
int 		markers_sort_rel_distance(Bmarker* set, Bmarker* ref);

/**
@brief 	Calculating marker statistics.
@param 	*markers	list of markers.
@return long		number of markers (<0 means failure).
**/
long		marker_stats(Bmarker* markers)
{
	long			n(0);
	Bmarker*		mark;
	Vector3<double>	min(1e30, 1e30, 1e30), max(-1e30, -1e30, -1e30);
	Vector3<double>	avg, std;

	for ( mark = markers; mark; mark = mark->next ) if ( mark->sel ) {
		min = min.min(mark->loc);
		max = max.max(mark->loc);
		avg += mark->loc;
		std += mark->loc * mark->loc;
		n++;
	}
	
	avg /= n;
	std = ((std/n) - (avg*avg)).square_root();
	
	cout << "Number of selected markers:     " << n << endl;
	cout << "Minimum:                        " << min << endl;
	cout << "Maximum:                        " << max << endl;
	cout << "Average:                        " << avg << endl;
	cout << "Standard deviation:             " << std << endl;
	
	return n;
}

/**
@brief 	Calculating the minimum marker location.
@param 	*markers		list of markers.
@return Vector3<double>	location minimum.
**/
Vector3<double>	marker_minimum(Bmarker* markers)
{
	Bmarker*		mark;
	Vector3<double>	min(1e30, 1e30, 1e30);

	for ( mark = markers; mark; mark = mark->next )
		if ( mark->sel )
			min = min.min(mark->loc);
	
	return min;
}

/**
@brief 	Calculating the maximum marker location.
@param 	*markers		list of markers.
@return Vector3<double>	location maximum.
**/
Vector3<double>	marker_maximum(Bmarker* markers)
{
	Bmarker*		mark;
	Vector3<double>	max(-1e30, -1e30, -1e30);

	for ( mark = markers; mark; mark = mark->next )
		if ( mark->sel )
			max = max.max(mark->loc);
	
	return max;
}

/**
@brief 	Calculating the marker location range.
@param 	*markers		list of markers.
@return Vector3<double>	location range.
**/
Vector3<double>	marker_range(Bmarker* markers)
{
	return marker_maximum(markers) - marker_minimum(markers);
}

/**
@brief 	Copies a set of markers.
@param 	*mark			marker list.
@param 	sel_flag		flag for selection.
@return Bmarker*			new marker list.
**/
Bmarker*	markers_copy(Bmarker* mark, int sel_flag)
{
	Bmarker*		newmark = NULL;
	Bmarker*		m = NULL;
	Bmarker*		m2 = NULL;
	long	nmark(0);
	
	for ( m = mark; m; m = m->next ) if ( !sel_flag || m->sel ) {
		m2 = (Bmarker *) add_item((char **) &m2, sizeof(Bmarker));
		if ( !newmark ) newmark = m2;
		m2->id = m->id;
		m2->img_num = m->img_num;
		m2->loc = m->loc;
		m2->err = m->err;
		m2->res = m->res;
		m2->fom = m->fom;
		m2->sel = m->sel;
		nmark++;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Marker model initialized with " << nmark << " markers" << endl << endl;
	
	return newmark;
}

/**
@brief 	Copies a set of markers.
@param 	*mark			marker list.
@return Bmarker*			new marker list.
**/
Bmarker*	markers_copy(Bmarker* mark)
{
	return markers_copy(mark, 0);
}

/**
@brief 	Copies a set of markers based on selection.
@param 	*mark			marker list.
@return Bmarker*				new marker list.
**/
Bmarker*	markers_copy_selected(Bmarker* mark)
{
	return markers_copy(mark, 1);
}

/**
@brief 	Adds a set of markers to another.
@param 	**mark			marker list.
@param 	*mark2			marker list to add.
@param 	mindist			minimum distance between markers.
@param 	sel					selection flag to set for second list.
@return Bmarker*				new marker list.

	The second list of markers are checked to eliminate those close to
	markers in the first set. The remaining markers are then added and 
	their selection flags set if the variable sel is greater than zero.

**/
long		markers_add(Bmarker** mark, Bmarker* mark2, double mindist, int sel)
{
	int				id(0);
	Bmarker*		m;
	Bmarker*		m2;

	for ( m = *mark; m; m = m->next ) if ( id < m->id ) id = m->id;
	for ( m2 = mark2; m2; m2 = m2->next ) m2->sel = 1;
	
	if ( *mark ) {
		for ( m = *mark; m; m = m->next )
			for ( m2 = mark2; m2; m2 = m2->next )
				if ( m->loc.distance(m2->loc) < mindist )
					m2->sel = 0;
		m = *mark;
		for ( m2 = mark2; m2; m2 = m2->next ) if ( m2->sel ) {
			m = (Bmarker *) add_item((char **) &m, sizeof(Bmarker));
			m->id = ++id;
			m->img_num = m2->img_num;
			m->loc = m2->loc;
			m->err = m2->err;
			m->res = m2->res;
			m->fom = m2->fom;
			m->sel = m2->sel;
			if ( sel > 0 ) m->sel = sel;
		}
		kill_list((char *) mark2, sizeof(Bmarker));
	} else *mark = mark2;
	
	long			nmark = count_list((char *) *mark);
	
	return nmark;
}

/**
@brief 	Deletes markers not selected.
@param 	**mark			pointer to marker list.
@return long			number of remaining markers.
**/
long	markers_delete_non_selected(Bmarker** mark)
{
	Bmarker*		m = markers_copy_selected(*mark);
	kill_list((char *) *mark, sizeof(Bmarker));
	*mark = m;
	
	return count_list((char *) m);
}
/*
long	markers_delete_non_selected(Bmarker** mark)
{
	Bmarker*		m = NULL;
	long			nmark(0);
	
	for ( m = *mark; m;  ) {
		if ( m->sel < 1 ) {
			m = (Bmarker *) remove_item((char **)mark, (char *)m, sizeof(Bmarker));
		} else {
			m = m->next;
			nmark++;
		}
	}
	
	return nmark;
}
*/

/**
@brief 	Shifting markers.
@param 	*markers	linked list of markers.
@param 	shift		shift vector.
@return int			error code (<0 means failure).
**/
int			marker_shift(Bmarker* markers, Vector3<double> shift)
{
	Bmarker*		mark;
	
	for ( mark = markers; mark; mark = mark->next )
		mark->loc += shift;
	
	return 0;
}

/**
@brief 	Scaling markers.
@param 	*markers	linked list of markers.
@param 	scale		scale vector.
@return int			error code (<0 means failure).
**/
int			marker_scale(Bmarker* markers, Vector3<double> scale)
{
	Bmarker*		mark;
	
	for ( mark = markers; mark; mark = mark->next ) {
		mark->loc *= scale;
		mark->err *= scale;
	}
	
	return 0;
}

/**
@brief 	Transforming markers.
@param 	*markers	linked list of markers.
@param 	t			transform structure.
@return int			number of markers.
**/
int			marker_transform(Bmarker* markers, Transform t)
{
	int				nmark(0);
	Bmarker*		mark;
	
	Matrix3			mat = Matrix3(t.axis, t.angle);
	mat = t.scale * mat;

	for ( mark = markers; mark; mark = mark->next, nmark++ ) {
		mark->loc = t * mark->loc;
		mark->err = mat * mark->err;
	}
	
	return nmark;
}

/**
@brief 	Compare marker sets and calculates the root-mean-square deviation.
@param 	*mark1		linked list of markers.
@param 	*mark2		linked list of markers.
@return double		RMSD.
**/
double		marker_compare(Bmarker* mark1, Bmarker* mark2)
{
	int				nmark(0);
	double			R(0);
	Bmarker*		m1;
	Bmarker*		m2;
	
	for ( m1 = mark1, m2 = mark2; m1 && m2; m1 = m1->next, m2 = m2->next, nmark++ ) {
		m2->err = m1->loc - m2->loc;
		m2->res = m2->err.length();
		R += m2->res * m2->res;
	}
	
	if ( nmark ) R = sqrt(R/nmark);
	
	return R;
}

/**
@brief 	Determines the rotation matrix between two sets of markers.
@param 	*markers1		linked list of markers.
@param 	*markers2		linked list of markers.
@param 	origin1	origin of first set.
@param 	origin2	origin of second set.
@return Matrix3					rotation matrix.

	The function calculates the rotation matrix required to superimposed 
	the first set onto the second set.
	Requirement: The two marker sets must share a significant number of markers.

**/
Matrix3		marker_find_matrix(Bmarker* markers1, Bmarker* markers2, 
				Vector3<double> origin1, Vector3<double> origin2)
{
	int				i, j, n(0);
	Bmarker*		mark1;
	Bmarker*		mark2;
	Vector3<double>	loc1, loc2;
	Matrix3			mat(1);
	vector<double>	x(4), y(4), bx(4), by(4), bz(4);
	Matrix			a(4,4);
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	for ( mark1 = markers1; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = markers2; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			loc1 = mark1->loc - origin1;
			loc2 = mark2->loc - origin2;
			for ( i=0; i<3; i++ ) {
				x[i] = loc1[i];
				y[i] = loc2[i];
			}
			x[3] = y[3] = 1;
			for ( i=0; i<4; i++ ) {
				for ( j=0; j<4; j++ )
					a[i][j] += x[i]*x[j];
				bx[i] += x[i]*y[0];
				by[i] += x[i]*y[1];
				bz[i] += x[i]*y[2];
			}
			n++;
		}
	}
	
	if ( n < 3 ) {
		cerr << "Error in marker_find_matrix: Too few shared markers! (" << n << ")" << endl << endl;
		return mat;
	}
	
//	cout << a;
	a.LU_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	a.multiply_in_place(bz);
//	cout << a;
	
	for ( i=0; i<3; i++ ) {
		mat[0][i] = bx[i];
		mat[1][i] = by[i];
		mat[2][i] = bz[i];
	}
	
	return mat;
}

/**
@brief 	Determines the rotation matrix between two sets of markers.
@param 	*markers1		linked list of markers.
@param 	*markers2		linked list of markers.
@param 	origin			origin of rotation.
@return Transform		transform structure.

	The function calculates the rotation matrix required to superimposed 
	the first set onto the second set.
	Requirement: The two marker sets must share a significant number of markers.

**/
Transform	marker_find_transform(Bmarker* markers1, Bmarker* markers2, Vector3<double> origin)
{
	int				i, j, n(0);
	Bmarker*		mark1;
	Bmarker*		mark2;
	Vector3<double>	loc1, loc2;
	Matrix3			mat(1);
	Matrix			a(4,4);
//	double			x[4], y[4], bx[4], by[4], bz[4];
	vector<double>	x(4), y(4), bx(4), by(4), bz(4);
	Transform		t;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG marker_find_transform: origin=" << origin << endl;
	
	for ( i=0; i<4; i++ ) bx[i] = by[i] = bz[i] = 0;
	
	for ( mark1 = markers1; mark1; mark1 = mark1->next ) if ( mark1->sel ) {
		for ( mark2 = markers2; mark2 && mark2->id != mark1->id; mark2 = mark2->next ) ;
		if ( mark2 && mark2->sel ) {
			loc1 = mark1->loc - origin;
			loc2 = mark2->loc - origin;
			for ( i=0; i<3; i++ ) {
				x[i] = loc1[i];
				y[i] = loc2[i];
			}
			x[3] = y[3] = 1;
			for ( i=0; i<4; i++ ) {
				for ( j=0; j<4; j++ )
					a[i][j] += x[i]*x[j];
				bx[i] += x[i]*y[0];
				by[i] += x[i]*y[1];
				bz[i] += x[i]*y[2];
			}
			n++;
		}
	}
	
	if ( n < 3 ) {
		cerr << "Error in marker_find_transform: Too few shared markers! (" << n << ")" << endl << endl;
		return t;
	}

	if ( verbose & VERB_DEBUG )
		cout << a << endl;
		
	t = transform_matrix_solve(a, bx, by, bz, 0);

	t.origin = origin;

	return t;
}

/**
@brief 	Finds the rotation and shift between two sets of vectors.
@param 	*set1		first coordinate set.
@param 	*set2		second coordinate set.
@param 	tolerance 	acceptable deviation for fit.
@return Transform	structure with shift, scale, rotation angle, and R factor.

	A large number of angles is tested and the shift calculated at each angle.
	The best angle is selected based on the sum of the standard deviations of
	the x, y and z shifts.
	Outliers are iteratively removed until the shift standard deviations drop
	below a tolerance value.
	The angle returned is the rotation applied to the first set to obtain
	an estimate of the second set.
	Requirement: The two sets must have the same number of points.

**/
Transform	markers_find_rottrans(Bmarker* set1, Bmarker* set2, double tolerance)
{
	Transform       t;
	long 	n = count_list((char *) set1);
	
	if ( n != count_list((char *) set2) ) {
		cerr << "Error: The two sets must have the same number of points!" << endl;
		return t;
	}
	
	// Set pointers to point arrays
	Bmarker*		m1;
	Bmarker*		m2;
	
	int				i, j(0), goodpoints = n, lastgoodpoints = n + 1, imaxdev; 
	double			xshift, yshift, xscale, yscale, phi, maxdev;
	double			xstdev, ystdev, bestxstdev = tolerance*100, bestystdev = tolerance*100;
	double			phi_min = -M_PI, phi_max = M_PI, phi_step = M_PI/10000;
	double			cos_phi, sin_phi;
	double			bestphi(0), bestxshift(0), bestyshift(0);
	double			bestxscale = 1, bestyscale = 1;
	double			R, bestR = 1e37;
	double			dx, dy, dd;
	
	double*			s1 = new double[n];
	double*			s2 = new double[n];
	
	// Data points to ignore
	int* 			ignore = new int[n];

	for ( i=0; i<n; i++ ) s1[i] = s2[i] = ignore[i] = 0;
	
//	if ( verbose & VERB_DEBUG )
	if ( verbose & VERB_FULL ) {
		cout << "Coordinate sets:" << endl;
		for ( m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next )
			cout << m1->loc[0] << "," << m1->loc[1] << " - " << m2->loc[0] << "," << m2->loc[1] << endl;
		cout << endl << "Fitting the coordinate sets:" << endl;
	}
	
	while ( ( bestxstdev > tolerance || bestystdev > tolerance ) 
			&& goodpoints > n/2 && goodpoints < lastgoodpoints ) {
		lastgoodpoints = goodpoints;
		bestR = 1e30; 
		for ( phi=phi_min; phi<=phi_max; phi+=phi_step ) {
			cos_phi = cos(phi);
			sin_phi = sin(phi);
			for ( i=j=0, m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next, i++ ) {
				if ( !ignore[i] ) {
					s1[j] = m1->loc[0]*cos_phi - m1->loc[1]*sin_phi;
					s2[j] = m2->loc[0];
					j++;
				}
			}
			linear_least_squares(0, j-1, s1, s2, xshift, xscale);
			for ( i=j=0, xstdev = 0, m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next, i++ ) {
				if ( !ignore[i] ) { 
					dx = m2->loc[0] - (xscale*s1[j] + xshift); 
					xstdev += dx*dx;
					s1[j] = m1->loc[0]*sin_phi + m1->loc[1]*cos_phi;
					s2[j] = m2->loc[1];
					j++;
				}
			}
			linear_least_squares(0, j-1, s1, s2, yshift, yscale);
			for ( i=j=0, ystdev = 0, m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next, i++ ) {
				if ( !ignore[i] ) { 
					dy = m2->loc[1] - (yscale*s1[j] + yshift); 
					ystdev += dy*dy;
					j++;
				}
			}
			R = 1e37;
			if ( j > 0 ) {
				R = sqrt((xstdev + ystdev)/j); 
				xstdev = sqrt(xstdev/j); 
				ystdev = sqrt(ystdev/j); 
			}
			if ( bestR > R ) { 
				bestxstdev = xstdev;
				bestystdev = ystdev;
				bestR = R; 
				bestphi = phi; 
				bestxshift = xshift; 
				bestyshift = yshift; 
				bestxscale = xscale; 
				bestyscale = yscale; 
			}
		}
	 
		if ( bestxscale < 0 && bestyscale < 0 && fabs(bestphi) > M_PI_2 ) {
			bestxscale = -bestxscale; 
			bestyscale = -bestyscale; 
			bestphi = bestphi + M_PI;
			if ( bestphi > M_PI ) bestphi -= TWOPI; 
		}
		
		// Calculate the coordinates for the second set using the first set 
		imaxdev = 0; 
		maxdev = 0; 
		cos_phi = cos(bestphi);
		sin_phi = sin(bestphi);
		for ( i=j=0, m1=set1, m2=set2; m1; m1=m1->next, m2=m2->next, i++ ) {
			if ( !ignore[i] ) { 
				dx = m2->loc[0] - (bestxscale*(m1->loc[0]*cos_phi - 
					m1->loc[1]*sin_phi) + bestxshift); 
				dy = m2->loc[1] - (bestyscale*(m1->loc[0]*sin_phi + 
					m1->loc[1]*cos_phi) + bestyshift); 
				dd = dx*dx + dy*dy;
				if ( maxdev < dd ) {
					maxdev = dd;
					imaxdev = i;
				}
				j++;
				if ( verbose & VERB_DEBUG ) {
					cout << "DEBUG markers_find_rottrans: " << dx << "," << dy << endl;
				}
			}
		}
		ignore[imaxdev] = 1;
		goodpoints--;
		phi_min = bestphi - 0.1;
		phi_max = bestphi + 0.1;
		if ( goodpoints < n/4 )
			cerr << "Error: Number of points used (" << goodpoints << ") too small!" << endl;
		if ( verbose >= VERB_FULL ) 
			cout << " " <<  
					j << " " << bestxshift << " " << bestyshift << " " << bestxscale << " " << bestyscale << " " <<  
					bestphi*180.0/M_PI << " (" << bestxstdev << " " << bestystdev << " " << bestR << ")" << endl; 
	}
	if ( verbose >= VERB_FULL ) 
		cout << "Fit: " <<  
			j << " " << bestxshift << " " << bestyshift << " " << bestxscale << " " << bestyscale << " " <<  
					bestphi*180.0/M_PI << " (" << bestxstdev << " " << bestystdev << " " << bestR << ")" << endl; 
	if ( verbose >= VERB_PROCESS ) cout << endl;
	
	delete[] s1;
	delete[] s2;
	delete[] ignore;

	t.trans[0] = bestxshift; 
	t.trans[1] = bestyshift; 
	t.scale[0] = bestxscale; 
	t.scale[1] = bestyscale; 
	t.angle = bestphi; 
	t.fom = bestR; 

	return t; 
}

/**
@author  Samuel Payne
@brief 	This function finds the transform parameters between sets of 
	coordinates.  These sets are expected to have a large intersection, 
	but do not have to be 1 to 1.
@param 	**set		set of coordinates from images to be aligned.
@param 	nseries			number of images in the set.
@param 	refset			reference set number.
@param 	size	size of images/frames of sets.
@return Transform*			list with shift, scale, rotation angle, and R factor.

	Gets a one to one mapping of the two coordinate sets using
	markers_find_rottrans to compare pairs of sets. 
	A postscript file is made to graphically show the deviations of each 
	point from its calculated position. 
	Assumptions: This assumes each set has been sorted according to
	distance from the center of the set.  	

**/
Transform* 	markers_map_and_find_transform(Bmarker** set, int nseries, int refset, Vector3<long> size)
{
	int				j, matched_pairs;
	Transform*      tlist = NULL;
	Transform*      t = NULL;
	Bstring			filename;
	
	if ( verbose & VERB_LABEL )
		cout << "X-shift  Y-shift  Z-shift  X-scale  Y-scale  Z-scale    Angle" << endl;
	 
	for ( j=0; j<nseries; j++ ) {
		t = (Transform *) add_item((char **) &t, sizeof(Transform));
		if ( !tlist ) tlist = t;
		*t = Transform();
		if ( j != refset ) {
			matched_pairs = markers_findOneToOne(set[refset],set[j]);
			if ( matched_pairs != 0 ) {
				*t = markers_find_rottrans(set[refset], set[j], 2);
				if ( t->fom ) {
					filename = Bstring(j, "final_errors%d.ps");
					ps_marker_errors(filename, set[refset], set[j], *t, size, 20);
					if ( verbose & VERB_LABEL )
						cout << t->trans << " " << t->scale << " " << t->angle*180/M_PI << endl;
				}
				if ( verbose & VERB_DEBUG ) {
					filename = Bstring(j, "final_match%d.ps");
					ps_marker_match(set[refset], set[j], filename);
				}
			} else {
				error_show("Error in markers_map_and_find_transform", __FILE__, __LINE__);
				cerr <<"Not able to match micrographs " << refset+1 << " and " << j+1 << endl;
				return tlist;
			}
		}
	}
	
	return tlist;
}

/*
@author  Samuel Payne
@brief 	Finds a one to one mapping of two coordinate sets.  The sets should 
	have a large intersection to begin with.
@param 	*refset		coordinate set serving as reference
@param 	*applyset	coordinate set to be matched to refset
@return int matched_pairs			the number of points in the 1 to 1 mapping	

	To get a rough matching, we attempt to match every point in B with 
	a point in A. (Looping through a few points in A as well.)	The 
	transformation parameters are recorded for the rough match with 
	the lowest error. Finally it goes through an exaustive brute force 
	matching with the estimates of the transform parameters and 
	produces the one to one mapping.  It rewrites the point* sets so 
	that the matching points have the same index.	
	(set1[0] matches set2[0], set1[1] matches set2[1], ...) 

**/
int 	markers_findOneToOne(Bmarker* refset, Bmarker* applyset)
{
	long 			num_matches=0;
	double 			error;
	Vector3<double>	calc;
	double			besterror = 100, phi(0), xshift(0), yshift(0), zshift(0);
	double*			match_params = NULL;
	Bmarker			*m, *mr, *mA, *mB, *setA, *setB, *bestB;
	
	long	nA = count_list((char *) refset);
	long	nB = count_list((char *) applyset);
	Vector3<double>	min = refset->loc;
	Vector3<double>	max = refset->loc;
	
//	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG markers_findOneToOne: nA=" << nA << " nB=" << nB << endl;
		
	for ( m=refset; m; m=m->next ) {
		min = min.min(m->loc);
		max = max.max(m->loc);
	}
	
	double			refsize = max.distance(min);
	Vector3<double>	refcenter = (max + min) * 0.5;
	
	double 			dist_cutoff = refsize *.03;
	double 			dist_err_cutoff = 0.5 * refsize;
	
	min = applyset->loc;
	max = applyset->loc;
	for ( m=applyset; m; m=m->next ) {
		min = min.min(m->loc);
		max = max.max(m->loc);
	}
	
//	double			applysize = max.distance(min);
	Vector3<double>	applycenter = (max + min) * 0.5;
	
	for ( mA=refset; mA; mA=mA->next ) {
		if ( mA->loc.distance(refcenter) < refsize/4 ) {
			setA = (Bmarker *) copy_list((char *) refset, sizeof(Bmarker));
			markers_sort_rel_distance(setA, mA);
			for ( mB=applyset; mB; mB=mB->next ) {
				if ( mB->loc.distance(applycenter) < refsize/4 ) {
					setB = (Bmarker *) copy_list((char *) applyset, sizeof(Bmarker));
					markers_sort_rel_distance(setB, mB);
					match_params = markers_RoughMatch(setA, setB, dist_cutoff);
					if ( (match_params[0] < besterror) && (match_params[4] > nA/3) ) {
						besterror = match_params[0];
						phi = match_params[1];
						xshift = match_params[2];
						yshift = match_params[3];
					}//end if besterr
					delete[] match_params;
					kill_list((char *) setB, sizeof(Bmarker));
				}
			}
			kill_list((char *) setA, sizeof(Bmarker));
		}
	}
	
//	if ( besterror > 99 ) return 0; //no matches with an error less than 100 are acceptable

//	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG markers_findOneToOne: Final match found. Parameter estimates:" << endl;
		cout << "\tangle " << phi*180/M_PI << " xshift " << xshift << " yshift " << yshift << endl;
//	}
	
	num_matches = 0;
	setA = (Bmarker *) copy_list((char *) refset, sizeof(Bmarker));
	setB = (Bmarker *) copy_list((char *) applyset, sizeof(Bmarker));
	m = applyset;
	mr = refset;
	for ( mA=setA; mA; mA=mA->next ) {
		calc[0] = mA->loc[0]* cos(phi) + mA->loc[1] * sin(phi) + xshift;
		calc[1] = -1*mA->loc[0]* sin(phi) + mA->loc[1] * cos(phi) + yshift;
		calc[2] = zshift;
		bestB = NULL;
		besterror = dist_err_cutoff; 
		for ( mB=setB; mB; mB=mB->next ) {
			error = mB->loc.distance(calc);
			if ( error < besterror ) {
				besterror = error;
				bestB = mB;
			}
		}
		if ( bestB ) {
			copy_item((char *)m, (char *)bestB, sizeof(Bmarker));
			copy_item((char *)mr, (char *)mA, sizeof(Bmarker));
			m = m->next;
			mr = mr->next;
			num_matches++;
		}
	}
			
	kill_list((char *) setA, sizeof(Bmarker));
	kill_list((char *) setB, sizeof(Bmarker));

	kill_list((char *) m, sizeof(Bmarker));
	kill_list((char *) mr, sizeof(Bmarker));
	
	return num_matches;
}

/*
@author  Samuel Payne
@brief	Produces a match of two sets and estimate the transform parameters.
@param	Bmarker* setA	the first coordinate set.
@param	Bmarker* setB	the second coordinate set.
@param	double dist_cutoff		max difference allowed in matching d_rel variables
@return	double*					the transform parameters, the set error, and size.

	This finds the rough matching of two sets of points. These sets are 
	assumed to be sorted relative to their distance to a relative point.  
	It then computes a histogram of all the angular deviations.  Taking 
	the largest bin as the true angle (not averaging all values as that 
	is subject to outliers). It removes all matches whose angular 
	deviations are further from the true angle than the angle tolerance 
	allows.  Finally it calculates the x and y shift parameters, error and 
	set size.
**/
double*		markers_RoughMatch(Bmarker* setA, Bmarker* setB, double dist_cutoff)
{
	int 			a, i, j, bins = 360, max;
	double 			best_phi=0, phi, ang_tol = M_PI/10, xshift, yshift, error = 1e37;
    double			scale = bins*1.0/(TWOPI);
	double			dx, dy, sumdx_sq, sumdy_sq, sumdx, sumdy;
	
	long	nA = count_list((char *) setA);
//	long	nB = count_list((char *) setB);
	
	Bmarker			*mA, *mB;
	double			d;

	Bmarker**		match = new Bmarker*[nA];
	int*			histo = new int[bins];
	double*			params = new double[5];
	double*			phi_devs = new double[nA];
	
	/* The loop starts at a=1, b=1 because the first points have a zero 
		distances. They always match each other (a trivial solution).  they 
		also have a zero angle, which is not characteristic of the set in 
		general. so it should be excluded. */
	for ( a=1, mA=setA->next, mB=setB->next; mA && mB; mA=mA->next, mB=mB->next, a++ ) {
//	for( a=1, b=1; a<nA && b<nB; a++, b++ ) {
		/* These while loops find an 'a' and 'b' that match each other according to the
			relative distance variable.  while the next point of b matches better than the
			current than advance the b counter. same for a.  This strategy works because the
			sets are sorted according to the d_rel variable (relative distance).  If the difference 
			between the distance vars is less than a cutoff the match is accepted. */
		d = mA->err.length2();
		while ( mB->next &&
				( fabs(d - mB->err.length2()) > fabs(d - mB->next->err.length2()) ) ) {
			mB = mB->next;
		}
		d = mB->err.length2();
		while ( mA->next &&
				( fabs(mA->err.length2() - d) > fabs(mA->next->err.length2() - d)) ) {
			mA = mA->next;
			a++;
		}
		if ( fabs(mA->err.length2() - mB->err.length2()) < dist_cutoff ) {
			/* Here a histogram of the angular deviations is set up.  There are 360 
				bins which are meant to correspond to the 360 degrees in a unit circle.
				The angles are shifted between 0 and 2PIs so they can be scaled into a 
				bin easily.  The additional 0.5 shift before casting to int ensures that 
				values slightly under a bin are put in bin they are closest to.  In the 
				bottom else j==bin, therefore j = 360, which corresponds to 360 degrees, 
				the same as 0 degrees.  therefore it is put into histo[0]. */
			match[a] = mB;
			phi = angle_set_negPI_to_PI(mA->fom - mB->fom);
			phi_devs[a] = phi;
	   		j = (int) (scale*(phi)+ 0.5);
			if ( j < bins ) histo[j]++;
			else histo[0]++; 
		}
	}
	
	//find the Max bin value.  that is an initial estimate of the angle
	for ( max=0, i=0; i<bins; i++ ) {
		if ( histo[i] > max ) {
			best_phi = i; //now holds the degree value of the angle
			max = histo[i];
		}
	}
	best_phi = best_phi / scale;  //now holds the RADIAN value
	xshift = yshift = j = 0;
	for ( a=1, mA=setA; mA && a<nA; a++, mA=mA->next ) {
		if ( match[a] ) {
			phi = fabs(phi_devs[a] - best_phi);
			if ( phi > M_PI ) phi -= TWOPI;
			if ( fabs(phi) > ang_tol ) {
				match[a] = NULL;
			} else {
				dx = match[a]->loc[0] - mA->loc[0]*cos(best_phi) - mA->loc[1]*sin(best_phi);
				dy = match[a]->loc[1] + mA->loc[0]*sin(best_phi) - mA->loc[1]*cos(best_phi);
				xshift += dx;
				yshift += dy;
				j++; //here j is acting as a count variable
			}
		}
	}
	xshift = xshift/j;
	yshift = yshift/j;
	sumdx = sumdy = sumdx_sq = sumdy_sq = j = 0;
	
	for ( a=1, mA=setA; mA && a<nA; a++, mA=mA->next ) {
		if ( match[a] ) {
			dx = match[a]->loc[0] - mA->loc[0]*cos(best_phi) - mA->loc[1]*sin(best_phi) - xshift;
			dy = match[a]->loc[1] + mA->loc[0]*sin(best_phi) - mA->loc[1]*cos(best_phi) - yshift;
			sumdx += dx;
			sumdy += dy;
			sumdx_sq += dx*dx;
			sumdy_sq += dy*dy;				
			j++; //here j is acting as a count var
		}
	}
	
	if ( j ) {
		sumdx = sumdx/j;
		sumdy = sumdy/j;
		error = sqrt((sumdx_sq + sumdy_sq)/j - sumdx*sumdx - sumdy*sumdy);
		if ( !isfinite(error) ) error = 1e37;
	}
	
	params[0] = error;
	params[1] = best_phi;
	params[2] = xshift;
	params[3] = yshift;
	params[4] = j; //num_matches
	
	if ( verbose & VERB_DEBUG ) {
		/*if ( error < 100 && j > nA/5 )*/
		if ( max > 15 ){
			cout << endl << "DEBUG markers_RoughMatch:" << endl;
			cout << "error " << error << " angle " << best_phi << " xshift " <<
					xshift << " yshift " << yshift << " count " << j << endl;
			cout << "DEBUG markers_RoughMatch: angular deviation histogram bins" << endl;
			for ( a=0; a<bins; a++ ) {
				if ( histo[a] != 0 )
					cout << "bin " << a << " value " << histo[a] << endl;
			}
		}
	}
	
	delete[] histo;
	delete[] match;
	delete[] phi_devs;
	
	return params;
}

/**
@author  Samuel Payne
@brief 	This method constricts a marker set to target_num points.
@param 	target_num		number of points desired in set.
@param 	sign			sign to direct sorting.
@param 	**set		pointer to set containing all the points.
@return long 				number of points selected.

	This finds the target_num darkest points in the set.  The extreme 
	grey scale value associated with each point is held in value field
	of the Point structure.  This is copied into an array and sorted. The 
	cutoff value is the target_num index in the array.  A new array of 
	points is created containing points below the cutoff and then returned.

**/
long 		markers_limit(int target_num, int sign, Bmarker** set)
{
	long 	n = count_list((char *) set);
	
	if ( target_num < 1 || n < 1 ) return -1;
	if ( target_num >= n ) return 0;

	long			i;
	Bmarker*		m, *mp;
	double*			toSort = new double[n];
	
	for ( i=0, m=*set; i<n && m; i++, m=m->next ) toSort[i] = m->fom;
	
	if ( sign < 0 )
		qsort((void *) toSort, n, sizeof(double), QsortSmallToLargeFloat);
	else
		qsort((void *) toSort, n, sizeof(double), QsortLargeToSmallFloat);
				
//	if ( verbose & VERB_DEBUG ) 
		cout << "DEBUG markers_limit: n=" << n << " target_num=" << target_num << endl;
	
	long			count(0);
	double			cutoff = toSort[target_num];
	
	delete[] toSort;
	
	for ( m=mp=*set; m; ) { cout << "Marker " << m->id << " Previous " << mp->id << endl;
		if ( sign*m->fom < sign*cutoff ) {
			cout << "Removing " << m->id << " from " << (*set)->id << endl;
			if ( *set == m ) *set = m->next;
			m = (Bmarker *) remove_item((char **) &mp, (char *) m, sizeof(Bmarker));
			if ( m ) cout << "New marker " << m->id << " Previous " << mp->id << " Set " << (*set)->id << endl;
		} else {
			count++;
			mp = m;
			m = m->next;
		}
	}
	
	if ( verbose & VERB_LABEL ) 
		cout << "Limiting the size of the coordinate set from " << n << " to " << 
			target_num << " (threshold = " << cutoff << ")" << endl;

	return count;
}

/*
@author  Samuel Payne
@brief 	Sorts a set of markers based on the distance to a reference marker.
@param 	*set		coordinate set.
@param 	*ref		reference point.
@return int 				0.

	The FOM is set to the angle of rotation of the vector connecting
	the reference marker to the marker.	

**/
int 		markers_sort_rel_distance(Bmarker* set, Bmarker* ref)
{
	Bmarker*	m, *mp;
	
	for ( m=set; m; m=m->next ) {
		m->err = m->loc - ref->loc;
		m->fom = m->err.length2();
	}
	
	long		count = 1;
	
	while ( count ) {
		count = 0;
		for ( m=mp=set; m; mp=m, m=m->next ) {
			if ( m->fom < mp->fom ) {
				mp->next = m->next;
				m->next = mp;
				m = mp;
				count++;
			}
		}
	}
	
	for ( m=set; m; m=m->next )
		m->fom = atan2(m->err[1], m->err[0]);
	
	return 0;
}

/**
@brief 	Convert markers to a model.
@param 	*markers	linked list of markers.
@param 	origin		marker origin.
@param 	sam			scaling marker coordinates.
@return Bmodel*		new model.
**/
Bmodel*		model_from_markers(Bmarker* markers, Vector3<double> origin, Vector3<double> sam)
{
	if ( !markers ) return NULL;
	
	Bmarker*		mark;
	string			id("Markers");
	Bmodel*			model = new Bmodel(id);
	Bcomponent*		comp = NULL;
	string			comptype("MRK");
//	Bstring			id;

	if ( verbose ) 
		cout << "Generating a model from markers" << endl << endl;

//	model->identifier("Markers");
	
	for ( mark = markers; mark; mark = mark->next ) {
//		id = Bstring(mark->id, "%d");
//		comp = component_add(&comp, id);
//		if ( !model->comp ) model->comp = comp;
		if ( comp ) comp = comp->add(mark->id);
		else model->comp = comp = new Bcomponent(mark->id);
//		comp->type = model_add_type_by_id(model, comptype);
		comp->type(model->add_type(comptype));
		comp->location(sam * (mark->loc - origin));
		comp->radius(sam[0]);
	}
	
	return model;
}

/**
@brief 	Fits a plane to a set of markers.
@param 	*markers		linked list of markers.
@param 	origin			marker origin.
@return Vector3<double>	normal to the plane.

	The plane is given by:
		a*x + b*y + c*z = d
	The fit is assessed by calculating the distance of each marker to
	the plane:
		R = sqrt(sum(|a*x + b*y + c*z - d|^2)/n)

**/
Vector3<double>	marker_plane(Bmarker* markers, Vector3<double> origin)
{	
    long     		i, j, n;
	double			offset, d, R;
	vector<double>	b(3);
	Matrix			a(3,3);
	Bmarker*		mark;
	Vector3<double>	loc, normal, sum;
	
	for ( i=0; i<3; i++ ) b[i] = 0;
	sum = 0;
	for ( n=0, mark = markers; mark; mark = mark->next, n++ ) {
		loc = mark->loc - origin;
		sum += loc;
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<=i; j++ )
				a[j][i] += loc[i]*loc[j];
			b[i] += loc[i];
		}
	}
	sum /= n;
	for ( i=1; i<3; i++ )
		for ( j=0; j<i; j++ )
			a[i][j] = a[j][i];
	offset = 0;
	if ( a[0][0] == 0 ) {
		normal = Vector3<double>(1, 0, 0);
	} else if ( a[1][1] == 0 ) {
		normal = Vector3<double>(0, 1, 0);
	} else if ( a[2][2] == 0 ) {
		normal = Vector3<double>(0, 0, 1);
	} else {
		a.LU_decomposition(b);
		normal = Vector3<double>(b[0], b[1], b[2]);
		normal.normalize();
		offset = sum.scalar(normal);
		cout << b[0] << tab << b[1] << tab << b[2] << tab << offset << endl;
	}
	n = 0;
	R = 0;
	cout << "Marker\tDistance" << endl;
	for ( mark = markers; mark; mark = mark->next ) {
		d = (mark->loc - origin).scalar(normal) - offset;
		cout << mark->id << tab << d << endl;
		R += d*d;
		n++;
	}
	R = sqrt(R/n);
	cout << "Solution: " << normal[0] << " x + " << normal[1] << " y + " << 
		normal[2] << " z = " << offset << "   R = " << R << endl;
	
	return normal;
}

/**
@brief 	Sorts a set of markers by ID number.
@param 	**mark		pointer to marker list.
@return Bmarker*			new pointer to marker list.
**/
Bmarker*	markers_sort_by_id(Bmarker** mark)
{
	if ( ! *mark ) return NULL;
	
	long	i, max_id(0);
	Bmarker*		m;
	
	for ( m = *mark; m; m = m->next ) if ( max_id <= m->id ) max_id = m->id + 1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG markers_sort_by_id: max_id = " << max_id << endl;
	
	Bmarker**		mlist = new Bmarker*[max_id];
	for ( i=0; i<max_id; i++ ) mlist[i] = NULL;
	
	for ( m = *mark; m; m = m->next ) mlist[m->id] = m;
	
	m = *mark = NULL;
	
	for ( i=0; i<max_id; i++ ) {
		if ( mlist[i] ) {
			if ( m ) {
				m->next = mlist[i];
				m = m->next;
			} else {
				m = mlist[i];
				*mark = m;
			}
		}
	}
	
	m->next = NULL;
	
	delete[] mlist;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG markers_sort_by_id: *mark = " << *mark << endl;
	
	return *mark;
}

/**
@brief 	Renumbers a set of markers.
@param 	**mark		pointer to marker list.
@return long		number of markers.
**/
long	markers_renumber(Bmarker* mark)
{
	if ( ! mark ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG markers_renumber" << endl;
	
	long	nmark(0);
	Bmarker*		m;
	
	for ( m = mark; m; m = m->next )
		m->id = ++nmark;
		
	return nmark;
}

/**
@brief 	Centers a set of markers.
@param 	**mark			pointer to marker list.
@param 	center	center to shift to.
@return long			number of markers.

	The marker center is calculated and the all markers shifted by the
	difference with the target center.

**/
long	markers_center(Bmarker* mark, Vector3<double> center)
{
	if ( ! mark ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG markers_center" << endl;
	
	long	nmark;
	Bmarker*		m;
	Vector3<double>	mark_center, shift;
	
	for ( nmark=0, m = mark; m; m = m->next, nmark++ )
		mark_center += m->loc;
	
	if ( nmark ) {
		mark_center /= nmark;
		shift = center - mark_center;
		for ( m = mark; m; m = m->next ) m->loc += shift;
	}
	
	return nmark;
}

/**
@brief 	Generates a marker set at a given radius with z=0.
@param 	rad				radius.
@param 	ainc				angular increment.
@return int						0, <0 on error.

	Markers are generated at the given radius and angular increments.

**/
Bmarker*	marker_set_at_radius(double rad, double ainc)
{
	long		n;
	double		a;
	Bmarker*	marklist = NULL;
	Bmarker*	mark = NULL;
	
	for ( n=0, a=-M_PI; a<M_PI; a+=ainc ) {
		mark = (Bmarker *) add_item((char **) &mark, sizeof(Bmarker));
		if ( !marklist ) marklist = mark;
		n++;
		mark->id = n;
		mark->loc = Vector3<double>(rad*cos(a), rad*sin(a), 0);
	}
	
	return marklist;
}

/**
@brief 	Generates a bild file with the two marker sets.
@param 	&filename		bild format file name.
@param 	*mark1			marker list 1.
@param 	*mark2			marker list 2.
@return int						0, <0 on error.

	A sphere is drawn for every marker.
	The first set is red and the second set is blue.

**/
int			marker_sets_to_bild(Bstring& filename, Bmarker* mark1, Bmarker* mark2)
{
	Bmarker*		mark;	
	
	if ( verbose )
		cout << "Generating a bild file:         " << filename << endl << endl;
	
	ofstream		fbld(filename.c_str());
	if ( fbld.fail() ) return -1;
	
	fbld << endl;
	fbld << ".color red" << endl;
	for ( mark = mark1; mark; mark = mark->next ) if ( mark->sel ) {
		fbld << ".sphere " << mark->loc[0] << " " << mark->loc[1] << " " << mark->loc[2] << " " << 10.0 << endl;
	}
	
	fbld << endl;
	fbld << ".color blue" << endl;
	for ( mark = mark2; mark; mark = mark->next ) if ( mark->sel ) {
		fbld << ".sphere " << mark->loc[0] << " " << mark->loc[1] << " " << mark->loc[2] << " " << 10.0 << endl;
	}
	
	fbld << endl;
	
	fbld.close();
		
	return 0;
}


