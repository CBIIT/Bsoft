/**
@file	model_tags.h
@brief	All STAR and XML file format tags for models
@author Bernard Heymann
@date	Created: 20000419
@date	Modified: 20210219
**/

// Do not change the constant names because they are referenced in code
// Changes in the tag strings will not affect program execution
//  but old data files will be uninterpretable

#define COMMENT					"comment"
#define	ID						"id"

// Multi-molecular models
#define MODEL					"model"
#define MODEL_ID				"model.id"
#define MODEL_TYPE				"model.type_id"
#define MODEL_HAND				"model.hand"
#define MODEL_SYM				"model.point_group"
#define MODEL_MAP_FILENAME		"model.map_file_name"
#define MODEL_MAP_NUMBER		"model.map_number"
#define MODEL_FOM				"model.fom"
#define MODEL_SELECT			"model.select"

// Multi-molecular model components
#define COMPONENT				"component"
#define COMPONENT_ID			"component.id"
#define COMPONENT_TYPE_ID		"component.type_id"
#define COMPONENT_X				"component.x"
#define COMPONENT_Y				"component.y"
#define COMPONENT_Z				"component.z"
#define COMPONENT_VIEW_X		"component.view_x"
#define COMPONENT_VIEW_Y		"component.view_y"
#define COMPONENT_VIEW_Z		"component.view_z"
#define COMPONENT_VIEW_ANGLE	"component.view_angle"
#define COMPONENT_RADIUS		"component.radius"
#define COMPONENT_RED			"component.red"
#define COMPONENT_GREEN			"component.green"
#define COMPONENT_BLUE			"component.blue"
#define COMPONENT_ALPHA			"component.alpha"
#define COMPONENT_DENSITY		"component.density"
#define COMPONENT_FOM			"component.fom"
#define COMPONENT_SELECT		"component.select"

// Multi-molecular model component types
#define COMPTYPE				"component_type"
#define COMPTYPE_ID				"component_type.id"
#define COMPTYPE_FILENAME		"component_type.file_name"
#define COMPTYPE_MAP_NUMBER		"component_type.map_number"	// Deprecated
#define COMPTYPE_NUMBER			"component_type.number"
#define COMPTYPE_MASS			"component_type.mass"
#define COMPTYPE_FOM			"component_type.fom"
#define COMPTYPE_SELECT			"component_type.select"

// Multi-molecular model component connectors
#define COMPLINK				"component_link"
#define COMPLINK_ID				"component_link.id"
#define COMPLINK_1				"component_link.component_1"
#define COMPLINK_2				"component_link.component_2"
#define COMPLINK_ANGLE			"component_link.rotation_angle"
#define COMPLINK_RADIUS			"component_link.radius"
#define COMPLINK_LENGTH			"component_link.length"
#define COMPLINK_RED			"component_link.red"
#define COMPLINK_GREEN			"component_link.green"
#define COMPLINK_BLUE			"component_link.blue"
#define COMPLINK_ALPHA			"component_link.alpha"
#define COMPLINK_FOM			"component_link.fom"
#define COMPLINK_SELECT			"component_link.select"

// Multi-molecular model component connector types
#define LINKTYPE				"link_type"
#define LINKTYPE_ID1			"link_type.component_1"
#define LINKTYPE_ID2			"link_type.component_2"
#define LINKTYPE_LENGTH			"link_type.length"
#define LINKTYPE_DISTANCE		"link_type.distance"
#define LINKTYPE_KLENGTH		"link_type.Klength"
#define LINKTYPE_KDISTANCE		"link_type.Kdistance"
#define LINKTYPE_FOM			"link_type.fom"
#define LINKTYPE_SELECT			"link_type.select"

// Multi-molecular model component connector types
#define ANGLETYPE				"angle_type"
#define ANGLETYPE_ID1			"angle_type.component_1"
#define ANGLETYPE_ID2			"angle_type.component_2"
#define ANGLETYPE_ID3			"angle_type.component_3"
#define ANGLETYPE_ANGLE			"angle_type.angle"
#define ANGLETYPE_KANGLE		"angle_type.Kangle"
#define ANGLETYPE_FOM			"angle_type.fom"
#define ANGLETYPE_SELECT		"angle_type.select"


