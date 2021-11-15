/**
@file	marker.h
@brief	Header file for image parameter functions
@author Bernard Heymann
@date	Created: 20051226
@date	Modified: 20210112
**/

#include "rwmodel.h"
#include "Transform.h"

#ifndef _MarkerLinker_
#define _MarkerLinker_
/************************************************************************
@Object: struct Bmarker
@Description:
	Structure for a feature marker.
@Features:
	3D coordinates for the center of the marker.
*************************************************************************/
struct Bmarker {
	Bmarker*		next;			// Next marker in list
	int				id;				// Marker identifier
	int				img_num;		// Marker image number
	Vector3<float>	loc;			// Coordinates in the micrograph or tomogram
	Vector3<float>	err;			// Error in coordinates
	float			res;			// Residual
	float			fom;			// FOM for marker
	int				sel;			// Selection flag
} ;

#endif

// Function prototypes
long		marker_stats(Bmarker* markers);
Vector3<double>	marker_minimum(Bmarker* markers);
Vector3<double>	marker_maximum(Bmarker* markers);
Vector3<double>	marker_range(Bmarker* markers);
Bmarker*	markers_copy(Bmarker* mark);
Bmarker*	markers_copy_selected(Bmarker* mark);
long		markers_add(Bmarker** mark, Bmarker* mark2, double mindist, int sel);
long		markers_delete_non_selected(Bmarker** mark);
int			marker_shift(Bmarker* markers, Vector3<double> shift);
int			marker_scale(Bmarker* markers, Vector3<double> scale);
int			marker_transform(Bmarker* markers, Transform t);
double		marker_compare(Bmarker* mark1, Bmarker* mark2);
Matrix3		marker_find_matrix(Bmarker* markers1, Bmarker* markers2,
				Vector3<double> origin1, Vector3<double> origin2);
Transform	marker_find_transform(Bmarker* markers1, Bmarker* markers2, Vector3<double> origin);
Transform	markers_find_rottrans(Bmarker* set1, Bmarker* set2, double tolerance);
Transform* 	markers_map_and_find_transform(Bmarker** set, int nseries, int refset, Vector3<long> size);
long 		markers_limit(int target_num, int sign, Bmarker** set);
Bmodel*		model_from_markers(Bmarker* markers, Vector3<double> origin, Vector3<double> sam);
Vector3<double>	marker_plane(Bmarker* markers, Vector3<double> origin);
Bmarker*	markers_sort_by_id(Bmarker** mark);
long		markers_renumber(Bmarker* mark);
long		markers_center(Bmarker* mark, Vector3<double> center);
Bmarker*	marker_set_at_radius(double rad, double ainc);
int			marker_sets_to_bild(Bstring& filename, Bmarker* mark1, Bmarker* mark2);

