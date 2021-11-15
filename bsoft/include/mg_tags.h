/**
@file	mg_tags.h
@brief	All STAR and XML file format tags for micrographs and reconstructions
@author Bernard Heymann
@date	Created: 20000419
@date	Modified: 20210728
**/

// Do not change the constant names because they are referenced in code
// Changes in the tag strings will not affect program execution
//  but old data files will be uninterpretable

// Data block tags or identifiers: Program block tags are in upper case
//#define IMAGE					"image"
//#define STRUCFAC				"structure_factors"
//#define MOLECULE				"molecule"
//#define SEQUENCE				"sequence"
//#define POLARFOURIER 			"PFT"
//#define THREEDIMRECON 		"3DR"

#define COMMENT					"comment"
#define	ID						"id"

#define	PROJECT					"project"
#define	FIELD					"field"
#define	FIELD_ID				"field.id"
#define	FIELD_SELECT			"field.select"
#define	FIELD_FOM				"field.fom"

// Map tags
#define MAP						"map.3D_reconstruction"
#define MAP_ID					"map.3D_reconstruction.id"
#define MAP_REFERENCE			"map.reference.file_name"
#define MAP_RECONSTRUCTION		"map.3D_reconstruction.file_name"
#define	MAP_TRANSFORM_FILE		"map.3D_reconstruction_fourier_transform.file_name"
#define	MAP_POWERSPEC_FILE		"map.3D_reconstruction_powerspectrum.file_name"
#define MAP_SIZE_X				"map.3D_reconstruction.size_x"
#define MAP_SIZE_Y				"map.3D_reconstruction.size_y"
#define MAP_SIZE_Z				"map.3D_reconstruction.size_z"
#define	MAP_ORIGIN_X			"map.3D_reconstruction.origin_x"
#define	MAP_ORIGIN_Y			"map.3D_reconstruction.origin_y"
#define	MAP_ORIGIN_Z			"map.3D_reconstruction.origin_z"
#define MAP_SCALE_X				"map.3D_reconstruction.scale_x"
#define MAP_SCALE_Y				"map.3D_reconstruction.scale_y"
#define MAP_SCALE_Z				"map.3D_reconstruction.scale_z"
#define MAP_VOXEL_SIZE			"map.3D_reconstruction.voxel_size" // Deprecated
#define MAP_VOXEL_SIZE_X		"map.3D_reconstruction.voxel_size_x"
#define MAP_VOXEL_SIZE_Y		"map.3D_reconstruction.voxel_size_y"
#define MAP_VOXEL_SIZE_Z		"map.3D_reconstruction.voxel_size_z"
#define MAP_SELECT				"map.3D_reconstruction.select"
#define MAP_FOM					"map.3D_reconstruction.fom"
#define MAP_MAGNIFICATION		"map.magnification"
#define	MAP_VIEW_X				"map.view_x"
#define	MAP_VIEW_Y				"map.view_y"
#define	MAP_VIEW_Z				"map.view_z"
#define	MAP_VIEW_ANGLE			"map.view_angle"
#define MAP_BACK_RWEIGHT		"map.back_projection.rweight"
#define MAP_MODEL				"map.model"
#define MAP_SYMMETRY			"map.symmetry"

// Parameters for micrographs
#define MICROGRAPH					"micrograph"
#define	MICROGRAPH_FILE				"micrograph.file_name"
#define	MICROGRAPH_FRAMES_FILE		"micrograph_frames.file_name"
#define	MICROGRAPH_FRAME			"micrograph_frame.number"
#define	MICROGRAPH_FRAME_PIXEL_X	"micrograph_frame.pixel_size_x"
#define	MICROGRAPH_FRAME_PIXEL_Y	"micrograph_frame.pixel_size_y"
#define	MICROGRAPH_FRAME_PIXEL_Z	"micrograph_frame.pixel_size_z"
#define	MICROGRAPH_FRAME_SHIFT_X	"micrograph_frame.shift_x"
#define	MICROGRAPH_FRAME_SHIFT_Y	"micrograph_frame.shift_y"
#define	MICROGRAPH_FRAME_SELECT		"micrograph_frame.select"
#define	MICROGRAPH_FRAME_FOM		"micrograph_frame.fom"
#define	MICROGRAPH_PARTICLE_FILE	"micrograph_particle.file_name"	// deprecated
#define	MICROGRAPH_FILAMENT_FILE	"micrograph_filament.file_name"	// deprecated
#define	MICROGRAPH_TRANSFORM_FILE	"micrograph_fourier_transform.file_name"
#define	MICROGRAPH_POWERSPEC_FILE	"micrograph_powerspectrum.file_name"
#define	MICROGRAPH_ID				"micrograph.id"
#define	MICROGRAPH_FIELD_ID			"micrograph.field_id"
#define	MICROGRAPH_NUMBER			"micrograph.number"
#define	MICROGRAPH_SELECT			"micrograph.select"
#define	MICROGRAPH_FOM				"micrograph.fom"
#define	MICROGRAPH_MAGNIFICATION	"micrograph.magnification"
#define	MICROGRAPH_SAMPLING			"micrograph.sampling"
#define	MICROGRAPH_PIXEL			"micrograph.pixel_size"	// Deprecated
#define	MICROGRAPH_PIXEL_X			"micrograph.pixel_size_x"
#define	MICROGRAPH_PIXEL_Y			"micrograph.pixel_size_y"
#define	MICROGRAPH_PIXEL_Z			"micrograph.pixel_size_z" // ???
#define	MICROGRAPH_UNITS			"micrograph.units" 		// Supercede by general tag
#define	MICROGRAPH_DOSE				"micrograph.electron_dose"
#define	MICROGRAPH_EXPOSURE			"micrograph.exposure"
#define	MICROGRAPH_INTENSITY		"micrograph.intensity"
#define	MICROGRAPH_WATER_RING		"micrograph.water_ring_index"
#define	MICROGRAPH_ORIGIN_X 		"micrograph.origin_x"
#define	MICROGRAPH_ORIGIN_Y 		"micrograph.origin_y"
#define	MICROGRAPH_ORIGIN_Z 		"micrograph.origin_z"
#define MICROGRAPH_SCALE_X			"micrograph.scale_x"
#define MICROGRAPH_SCALE_Y			"micrograph.scale_y"
#define MICROGRAPH_SCALE_Z			"micrograph.scale_z"
#define	MICROGRAPH_TILT_AXIS		"micrograph.tilt_axis"
#define	MICROGRAPH_TILT_ANGLE		"micrograph.tilt_angle"
#define	MICROGRAPH_LEVEL_ANGLE		"micrograph.level_angle"
#define MICROGRAPH_ROT_ANGLE		"micrograph.rotation_angle"	// in-plane rotation
#define	MICROGRAPH_VIEW_X			"micrograph.view_x"
#define	MICROGRAPH_VIEW_Y			"micrograph.view_y"
#define	MICROGRAPH_VIEW_Z			"micrograph.view_z"
#define	MICROGRAPH_VIEW_ANGLE		"micrograph.view_angle"
#define MICROGRAPH_MATRIX_1_1		"micrograph.matrix_1_1"
#define MICROGRAPH_MATRIX_1_2		"micrograph.matrix_1_2"
#define MICROGRAPH_MATRIX_1_3		"micrograph.matrix_1_3"
#define MICROGRAPH_MATRIX_2_1		"micrograph.matrix_2_1"
#define MICROGRAPH_MATRIX_2_2		"micrograph.matrix_2_2"
#define MICROGRAPH_MATRIX_2_3		"micrograph.matrix_2_3"
#define MICROGRAPH_MATRIX_3_1		"micrograph.matrix_3_1"
#define MICROGRAPH_MATRIX_3_2		"micrograph.matrix_3_2"
#define MICROGRAPH_MATRIX_3_3		"micrograph.matrix_3_3"
#define MICROGRAPH_HVEC_X			"micrograph.h_x"
#define MICROGRAPH_HVEC_Y			"micrograph.h_y"
#define MICROGRAPH_HVEC_Z			"micrograph.h_z"
#define MICROGRAPH_KVEC_X			"micrograph.k_x"
#define MICROGRAPH_KVEC_Y			"micrograph.k_y"
#define MICROGRAPH_KVEC_Z			"micrograph.k_z"
#define MICROGRAPH_LVEC_X			"micrograph.l_x"
#define MICROGRAPH_LVEC_Y			"micrograph.l_y"
#define MICROGRAPH_LVEC_Z			"micrograph.l_z"
#define MICROGRAPH_HELIX_AXIS		"micrograph.helix_axis_angle"
#define MICROGRAPH_HELIX_RISE		"micrograph.helix_subunit_rise"
#define MICROGRAPH_HELIX_ANGLE		"micrograph.helix_subunit_angle"
#define MICROGRAPH_HELIX_RADIUS		"micrograph.helix_radius"
#define	MICROGRAPH_VOLTAGE			"micrograph.voltage"	// deprecated
#define	MICROGRAPH_CTF_CS			"micrograph.ctf.Cs"	// deprecated
#define	MICROGRAPH_CTF_CC			"micrograph.ctf.Cc"	// deprecated
#define	MICROGRAPH_CTF_ALPHA		"micrograph.ctf.alpha"	// deprecated
#define	MICROGRAPH_CTF_DE			"micrograph.ctf.energy_spread"	// deprecated
#define	MICROGRAPH_CTF_AMP_CONT		"micrograph.ctf.amp_contrast"	// deprecated
#define	MICROGRAPH_CTF_ZERO			"micrograph.ctf.first_zero"	// deprecated
#define	MICROGRAPH_CTF_DEF_AVG		"micrograph.ctf.defocus_average"	// deprecated
#define	MICROGRAPH_CTF_DEF_DEV		"micrograph.ctf.defocus_deviation"	// deprecated
#define MICROGRAPH_CTF_DEF_MIN		"micrograph.ctf.defocus_min"	// deprecated
#define MICROGRAPH_CTF_DEF_MAX		"micrograph.ctf.defocus_max"	// deprecated
#define	MICROGRAPH_CTF_AST_ANG		"micrograph.ctf.astigmatism_angle"	// deprecated
#define	MICROGRAPH_CTF_BASELINE		"micrograph.ctf.baseline"	// deprecated
#define	MICROGRAPH_CTF_ENVELOPE 	"micrograph.ctf.envelope"	// deprecated
#define	MICROGRAPH_BOX_RADIUS		"micrograph.box_radius"	// deprecated
#define	MICROGRAPH_BOX_RADIUS_X		"micrograph.box_radius_x"	// deprecated
#define	MICROGRAPH_BOX_RADIUS_Y		"micrograph.box_radius_y"	// deprecated
#define	MICROGRAPH_BOX_RADIUS_Z		"micrograph.box_radius_z"	// deprecated
#define	MICROGRAPH_BAD				"micrograph.bad"	// deprecated
#define	MICROGRAPH_BAD_RADIUS		"micrograph.bad_radius"	// deprecated
#define	MICROGRAPH_BAD_X			"micrograph.bad_x"	// deprecated
#define	MICROGRAPH_BAD_Y			"micrograph.bad_y"	// deprecated
#define	MICROGRAPH_BAD_Z			"micrograph.bad_z"	// deprecated
#define MICROGRAPH_MARKER_RADIUS	"micrograph.marker_radius"	// deprecated
#define MICROGRAPH_MARKER_ID		"micrograph.marker_id"	// deprecated
#define MICROGRAPH_MARKER_X			"micrograph.marker_x"	// deprecated
#define MICROGRAPH_MARKER_Y			"micrograph.marker_y"	// deprecated
#define MICROGRAPH_MARKER_Z			"micrograph.marker_z"	// deprecated
#define MICROGRAPH_MARKER_ERROR_X	"micrograph.marker_error_x"	// deprecated
#define MICROGRAPH_MARKER_ERROR_Y	"micrograph.marker_error_y"	// deprecated
#define MICROGRAPH_MARKER_ERROR_Z	"micrograph.marker_error_z"	// deprecated
#define MICROGRAPH_MARKER_FOM		"micrograph.marker_fom"	// deprecated
#define MICROGRAPH_FILAMENT_WIDTH	"micrograph.filament_width"	// deprecated
#define MICROGRAPH_FILNODE_RADIUS	"micrograph.filament_node_radius"	// deprecated

// Parameters for the CTF
#define	CTF							"ctf"
#define	CTF_VOLTAGE					"ctf.voltage"
#define	CTF_FOCAL					"ctf.focal_length"
#define	CTF_APERTURE				"ctf.objective_aperture"
#define	CTF_CS						"ctf.Cs"
#define	CTF_CC						"ctf.Cc"
#define	CTF_ALPHA					"ctf.alpha"
#define	CTF_DE						"ctf.energy_spread"
#define	CTF_AMP						"ctf.amp_contrast"
#define	CTF_AMP_SHIFT				"ctf.amplitude_phase_shift"
#define	CTF_ZERO					"ctf.first_zero"
#define	CTF_DEF_AVG					"ctf.defocus_average"
#define	CTF_DEF_DEV					"ctf.defocus_deviation"
#define CTF_DEF_MIN					"ctf.defocus_min"
#define CTF_DEF_MAX					"ctf.defocus_max"
#define	CTF_AST_ANG					"ctf.astigmatism_angle"
#define	CTF_BASELINE				"ctf.baseline"
#define	CTF_ENVELOPE				"ctf.envelope"

// Parameters for each individual particle image
#define PARTICLE				"particle"
#define	PARTICLE_FILE			"particle.file_name"
#define	PARTICLE_NUMBER			"particle.number"			// Deprecated
#define	PARTICLE_ID				"particle.id"
#define	PARTICLE_GROUP			"particle.group_id"
#define	PARTICLE_MG_ID			"particle.micrograph_id"
#define	PARTICLE_MG_X			"particle.micrograph_x"	// Deprecated
#define	PARTICLE_MG_Y			"particle.micrograph_y"	// Deprecated
#define	PARTICLE_MG_Z			"particle.micrograph_z"	// Deprecated
#define	PARTICLE_X				"particle.x"
#define	PARTICLE_Y				"particle.y"
#define	PARTICLE_Z				"particle.z"
#define	PARTICLE_X_ORIGIN		"particle.x_origin"		// Deprecated
#define	PARTICLE_Y_ORIGIN		"particle.y_origin"		// Deprecated
#define	PARTICLE_Z_ORIGIN		"particle.z_origin"		// Deprecated
#define	PARTICLE_ORIGIN_X		"particle.origin_x"
#define	PARTICLE_ORIGIN_Y		"particle.origin_y"
#define	PARTICLE_ORIGIN_Z		"particle.origin_z"
#define	PARTICLE_PSI			"particle.psi"
#define	PARTICLE_THETA			"particle.theta"
#define	PARTICLE_PHI			"particle.phi"
#define	PARTICLE_OMEGA			"particle.omega"
#define	PARTICLE_VIEW_X			"particle.view_x"
#define	PARTICLE_VIEW_Y			"particle.view_y"
#define	PARTICLE_VIEW_Z			"particle.view_z"
#define	PARTICLE_VIEW_ANGLE		"particle.view_angle"
#define	PARTICLE_MAGNIF			"particle.magnification" // Deprecated
#define	PARTICLE_PIXEL			"particle.pixel_size"	// Deprecated
#define	PARTICLE_PIXEL_X		"particle.pixel_size_x"
#define	PARTICLE_PIXEL_Y		"particle.pixel_size_y"
#define	PARTICLE_PIXEL_Z		"particle.pixel_size_z"
#define	PARTICLE_DEFOCUS		"particle.defocus"
#define	PARTICLE_DEF_DEV		"particle.defocus_deviation"
#define	PARTICLE_AST_ANG		"particle.astigmatism_angle"
#define PARTICLE_SELECT 		"particle.select"
#define	PARTICLE_FOM			"particle.fom"
#define	PARTICLE_FOM_CV			"particle.fom_crossvalidation"
#define	PARTICLE_FOM_SNR		"particle.signal_to_noise"
#define	PARTICLE_FOM_AVG		"particle.fom_average"
#define	PARTICLE_FOM_STD		"particle.fom_stdev"
#define	PARTICLE_FOM1			"particle.fom1"
#define	PARTICLE_FOM2			"particle.fom2"
#define	PARTICLE_FOM3			"particle.fom3"
#define	PARTICLE_FOM4			"particle.fom4"
#define	PARTICLE_FOM5			"particle.fom5"
#define	PARTICLE_FOM6			"particle.fom6"
#define	PARTICLE_FOM7			"particle.fom7"
#define	PARTICLE_FOM8			"particle.fom8"
#define	PARTICLE_FOM9			"particle.fom9"
#define PARTICLE_HANDA_FOM		"particle.handa_fom"
#define PARTICLE_HANDB_FOM		"particle.handb_fom"
#define	PARTICLE_CC				"particle.cc"
#define	PARTICLE_CC_AVG			"particle.cc_avg"
#define	PARTICLE_PFT_CC			"particle.pft_cc"
#define	PARTICLE_PRJ_CC			"particle.prj_cc"
#define	PARTICLE_CMP_CC			"particle.cmp_cc"
#define	PARTICLE_RFACTORAB		"particle.rfactorab"
#define	PARTICLE_COVERAGE		"particle.coverage"
#define	PARTICLE_DENSITY		"particle.density"
#define	PARTICLE_BOX_SIZE		"particle.box_size"
#define	PARTICLE_BOX_SIZE_X		"particle.box_size_x"
#define	PARTICLE_BOX_SIZE_Y		"particle.box_size_y"
#define	PARTICLE_BOX_SIZE_Z		"particle.box_size_z"
#define	PARTICLE_BOX_RADIUS		"particle.box_radius"		// Deprecated
#define	PARTICLE_BOX_RADIUS_X	"particle.box_radius_x"		// Deprecated
#define	PARTICLE_BOX_RADIUS_Y	"particle.box_radius_y"		// Deprecated
#define	PARTICLE_BOX_RADIUS_Z	"particle.box_radius_z"		// Deprecated
#define	PARTICLE_BAD			"particle.bad"
#define	PARTICLE_BAD_RADIUS		"particle.bad_radius"
#define	PARTICLE_BAD_X			"particle.bad_x"
#define	PARTICLE_BAD_Y			"particle.bad_y"
#define	PARTICLE_BAD_Z			"particle.bad_z"

// Tag for class average section
#define	CLASS_AVERAGE			"class_average"

// Parameters for filaments
#define	FILAMENT				"filament"
#define	FILAMENT_FILE			"filament.file_name"
#define	FILAMENT_ID				"filament.id"
#define	FILAMENT_NODE			"filament.node"
#define	FILAMENT_NODE_ID		"filament.node_id"
#define	FILAMENT_NODE_X			"filament.x"
#define	FILAMENT_NODE_Y			"filament.y"
#define	FILAMENT_NODE_Z			"filament.z"
#define FILAMENT_WIDTH			"filament.width"
#define FILNODE_RADIUS			"filament.node_radius"

// Orientation parameter sets
#define ORIENT_ID				"orient.id"
#define ORIENT_ORIGIN_X			"orient.origin_x"
#define ORIENT_ORIGIN_Y			"orient.origin_y"
#define ORIENT_ORIGIN_Z			"orient.origin_z"
#define ORIENT_VIEW_X			"orient.view_x"
#define ORIENT_VIEW_Y			"orient.view_y"
#define ORIENT_VIEW_Z			"orient.view_z"
#define ORIENT_VIEW_ANGLE		"orient.view_angle"
#define ORIENT_FOM				"orient.fom"
#define ORIENT_SELECT			"orient.select"

// Marker parameters
#define MARKER					"marker"
#define MARKER_RADIUS			"marker.radius"
#define MARKER_ID				"marker.id"
#define MARKER_X				"marker.x"
#define MARKER_Y				"marker.y"
#define MARKER_Z				"marker.z"
#define MARKER_ERROR_X			"marker.error_x"
#define MARKER_ERROR_Y			"marker.error_y"
#define MARKER_ERROR_Z			"marker.error_z"
#define MARKER_IMAGE			"marker.image_number"
#define MARKER_RESIDUAL			"marker.residual"
#define MARKER_FOM				"marker.fom"
#define MARKER_SELECT			"marker.select"

// Structure factors
#define REFLEX					"refln"
#define REFLEX_RADIUS			"refln.radius"
#define REFLEX_X 				"refln.location_x"
#define REFLEX_Y 				"refln.location_y"
#define REFLEX_Z 				"refln.location_z"
#define REFLEX_H 				"refln.index_h"
#define REFLEX_K 				"refln.index_k"
#define REFLEX_L 				"refln.index_l"
#define REFLEX_AMP  			"refln.F_meas_au"
//#define REFLEX_SIGAMP 		"refln.F_sigma"
#define REFLEX_SIGAMP 			"refln.F_meas_sigma_au"
#define REFLEX_PHI  			"refln.phase_meas"
#define REFLEX_SIGPHI 			"refln.phase_sigma"
#define REFLEX_FOM  			"refln.weight"
#define REFLEX_STATUS 			"refln.status"

// Layer lines
#define LAYERLINE				"layer_line"
#define LAYERLINE_NUMBER		"layer_line.number"
#define LAYERLINE_ORDER			"layer_line.bessel_order"
#define LAYERLINE_DISTANCE		"layer_line.distance"
#define LAYERLINE_FREQ			"layer_line.frequency"
#define LAYERLINE_AMP			"layer_line.amplitude"
#define LAYERLINE_FOM			"layer_line.fom"
#define LAYERLINE_SELECT		"layer_line.select"
