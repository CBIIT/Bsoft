/**
@file	rwmodel_bild.h
@brief	Header file for reading and writing Chimera BILD files
@author Bernard Heymann
@date	Created: 20140706
@date	Modified: 20210305
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_bild(Bstring* file_list);
int			write_model_bild(Bstring& filename, Bmodel* model);
int			model_to_bild_orientations(Bstring& filename, Bmodel* model, int vec_type, int color_type);
int			model_to_bild_view_sphere(Bstring& filename, Bmodel* model, int color_type);
int			model_to_bild_force_vectors(Bstring& filename, Bmodel* model, int color_type);
int			model_to_bild_view_polygons(Bstring& filename, Bmodel* model, int order, int color_type);
int			model_to_bild_polygons(Bstring& filename, Bmodel* model, int color_type);
int			model_to_bild_neighbor_planes(Bstring& filename, Bmodel* model, int color_type);


