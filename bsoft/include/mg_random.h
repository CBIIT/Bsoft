/**
@file	mg_random.h
@brief	Header file for applying some randomization to parameters
@author Bernard Heymann
@date	Created: 20000426
@date	Modified: 20101116
**/

#include "mg_processing.h"
#include "symmetry.h"


// Function prototypes
int			project_random_defocus(Bproject* project, double std);
int			project_random_origins(Bproject* project, double std);
int			project_random_views(Bproject* project);
int			project_random_helical_views(Bproject* project);
int			project_random_views(Bproject* project, double std);
int			project_random_magnification(Bproject* project, double std);
int			project_random_symmetry_related_views(Bproject* project, Bsymmetry& sym);



