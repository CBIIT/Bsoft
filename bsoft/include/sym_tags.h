/**
@file	sym_tags.h
@brief	All STAR file format tags for symmetry parameters
@author Bernard Heymann
@date	Created: 20000419
@date	Modified: 20080924
**/

// Do not change the constant names because they are referenced in code
// Changes in the tag strings will not affect program execution
//  but old data files will be uninterpretable

#define	SYMMETRY				"symmetry"

// Symmetry information
#define SYMMETRY_TYPE			"symmetry.cell_setting"
#define SYMMETRY_NUMBER 		"symmetry.Int_Tables_number"
#define SYMMETRY_NAME			"symmetry.space_group_name_H-M"
#define SYMMETRY_EQUIVID		"symmetry_equiv.id"
#define SYMMETRY_EQUIVXYZ		"symmetry_equiv.pos_as_xyz"
#define SYMMETRY_POINT_GROUP 	"symmetry.point_group"
#define SYMMETRY_PG_NUMBER	 	"symmetry.point_group_number"
#define SYMMETRY_AXIS_ORDER 	"symmetry_axis.order"
#define SYMMETRY_AXIS_X 		"symmetry_axis.x"
#define SYMMETRY_AXIS_Y 		"symmetry_axis.y"
#define SYMMETRY_AXIS_Z 		"symmetry_axis.z"

// Crystallographic unit cell parameters
#define UNIT_CELL_A				"cell.length_a"
#define UNIT_CELL_B				"cell.length_b"
#define UNIT_CELL_C				"cell.length_c"
#define UNIT_CELL_ALPHA 		"cell.angle_alpha"
#define UNIT_CELL_BETA			"cell.angle_beta"
#define UNIT_CELL_GAMMA			"cell.angle_gamma"


