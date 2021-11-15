/**
@file	mg_helix.h
@brief	Header file for functions to process helical data.
@author Bernard Heymann
@date	Created: 20061110
@date	Modified: 20151227
**/

#include "mg_processing.h"
#include "rwimg.h"

// Function prototypes
double*		filament_profile(Bfilnode* fnode, Bimage* p, long img_num, int id, double width, long& n);
double		filaments_center(Bfilament* fillist, Bimage* p, long img_num, int filament_width);
double		project_center_filaments(Bproject* project, int filament_width);
int			project_filament_powerspectrum(Bproject* project, int pad, int rotated, Bstring& path);
Bimage*		project_filament_density(Bproject* project, int filament_width);
int			mg_generate_layer_lines(Bmicrograph* mg, int rad_lim);
int			img_mask_layer_lines(Bimage* p, Blayerline* layer_line, float helix_axis, float width);
double*		img_extract_layer_line(Bimage* p, Blayerline* line, float helix_axis, int length);
int			mg_extract_show_layer_lines(Bmicrograph* mg, int length, int show);

