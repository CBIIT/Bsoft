/**
@file	mg_marker.h
@brief	Header file for micrograph marker functions
@author Bernard Heymann
@date	Created: 20150212
@date	Modified: 20150212
**/

#include "marker.h"
#include "mg_processing.h"

// Function prototypes
long		project_marker_lists(Bproject* project);
long		project_marker_in_particle(Bproject* project);
long		project_marker_in_particle_image(Bproject* project, double marker_radius);

