/**
@file	mg_subtomo.h
@brief	Header file for functions finding particles (subtomograms) in a tomogram
@author	Juha Huiskonen
@author	Bernard Heymann
@date	Created:  20071010
@date	Modified: 20120124 (BH)
@date	Modified: 20150108 (BH) - incorporated into Bsoft
@date	Modified: 20150720 (BH)
**/

#include "rwimg.h"
#include "View.h"
#include "symmetry.h"
#include "mg_processing.h"

// Sphere structure
#ifndef _Sphere_
#define _Sphere_
struct Sphere {
        double A,B,C;	// Sphere coordinates
        double R;	// Sphere radius
} ;
#endif

// Function prototypes

long	reconstruction_search_subtomo(Breconstruction* rec,
			Bimage* p, Bimage* ptemp, Bimage* pmask, Bimage* pmask2,
			double alpha_step, double theta_step, double phi_step,
			double alpha_limit, double thetaphi_limit,
			double hires, double lores, double shiftlimit, double mindist,
			double threshold, int maxhits, Vector3<long> bin, Bsymmetry& sym,
			int refinepeaks, Bstring ccmax_file);

long	reconstruction_refine_subtomo(Breconstruction* rec,
			Bimage* p, Bimage* ptemp, Bimage* pmask, Bimage* pmask2,
			double alpha_step, double theta_step, double phi_step,
			double alpha_limit, double thetaphi_limit,
			double hires, double lores, double shiftlimit,
			double shiftlimitz, double shiftlimitxy, double mindist,
			Vector3<long> bin, Bsymmetry& sym, int iters, int refinepeaks, Bstring ccmax_file);

Bparticle*	img_find_refine_peaks(Bimage* pcc, View view, double shift_limit,
				double shift_along, double shift_orthogonal, double mindist,
				double threshold, int maxhits, int refinepeaks);

double		img_find_peak_subtomo(Bimage* p, View view, double shift,
				double shift_along, double shift_orthogonal);




Vector3<double>  closest_point_line(Vector3<double> point, Vector3<double> endpoint1, Vector3<double> endpoint2);
double  		closest_point_line_distance2(Vector3<double> point, Vector3<double> endpoint1, Vector3<double> endpoint2);
//double  		closest_point_plane_distance2(Vector3<double> point, Vector3<double> point2, View view);
double  		closest_point_disc_distance2(Vector3<double> point, Vector3<double> point2, View view, double radius);

//Vector3<double>	closest_point_plane(Vector3<double> point, Vector3<double> point2, View vieW);

Sphere		locations_fit_sphere(Bparticle* part, int N, double Nstop);

