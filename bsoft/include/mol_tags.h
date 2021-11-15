/**
@file	mol_tags.h
@brief	All STAR file format tags for molecules and molecular parameters
@author Bernard Heymann
@date	Created: 20000419
@date	Modified: 20190603
**/

// Do not change the constant names because they are referenced in code
// Changes in the tag strings will not affect program execution
//  but old data files will be uninterpretable

// Crystallographic unit cell
#define	CELL_A					"cell.length_a"
#define	CELL_B					"cell.length_b"
#define	CELL_C					"cell.length_c"
#define	CELL_ALPHA				"cell.angle_alpha"
#define	CELL_BETA				"cell.angle_beta"
#define	CELL_GAMMA				"cell.angle_gamma"

// Atom type properties
#define ATOM_TYPE				"atom_type"
#define ATOM_TYPE_SYMBOL		"atom_type.symbol"
#define ATOM_TYPE_NUMBER	 	"atom_type.atom_number"
#define ATOM_TYPE_OXIDATION 	"atom_type.oxidation_number"
#define ATOM_TYPE_MASS			"atom_type.mass"
#define ATOM_TYPE_RADIUS_BOND 	"atom_type.radius_bond"
#define ATOM_TYPE_RADIUS_VDW 	"atom_type.radius_contact"
#define ATOM_TYPE_SCAT_A1 		"atom_type.scat_Cromer_Mann_a1"
#define ATOM_TYPE_SCAT_A2 		"atom_type.scat_Cromer_Mann_a2"
#define ATOM_TYPE_SCAT_A3 		"atom_type.scat_Cromer_Mann_a3"
#define ATOM_TYPE_SCAT_A4 		"atom_type.scat_Cromer_Mann_a4"
#define ATOM_TYPE_SCAT_A5 		"atom_type.scat_Cromer_Mann_a5"
#define ATOM_TYPE_SCAT_B1 		"atom_type.scat_Cromer_Mann_b1"
#define ATOM_TYPE_SCAT_B2 		"atom_type.scat_Cromer_Mann_b2"
#define ATOM_TYPE_SCAT_B3 		"atom_type.scat_Cromer_Mann_b3"
#define ATOM_TYPE_SCAT_B4		"atom_type.scat_Cromer_Mann_b4"
#define ATOM_TYPE_SCAT_B5		"atom_type.scat_Cromer_Mann_b5"
#define ATOM_TYPE_SCAT_C 		"atom_type.scat_Cromer_Mann_c"
#define ATOM_TYPE_COUNT	 		"atom_type.count"

// Bond rules
#define BOND_RULE_SYMBOL		"bond_rule.symbol"
#define BOND_RULE_BONDS			"bond_rule.bonds"

// Bond type properties
#define BOND_TYPE_SYMBOL1		"bond_type.symbol1"
#define BOND_TYPE_SYMBOL2		"bond_type.symbol2"
#define BOND_TYPE_LENGTH		"bond_type.length"
#define BOND_TYPE_VDWDIST		"bond_type.VdW_distance"

// Angle type properties
#define ANGLE_TYPE_SYMBOL1		"angle_type.symbol1"
#define ANGLE_TYPE_SYMBOL2		"angle_type.symbol2"
#define ANGLE_TYPE_SYMBOL3		"angle_type.symbol3"
#define ANGLE_TYPE_ANGLE		"angle_type.angle"

// Specific atom parameters
#define ATOM_PDB				"atom_site.group_PDB"
#define ATOM_RESNUMBER  		"atom_site.label_seq_id"
#define ATOM_RESINSERT  		"atom_site.pdbx_PDB_ins_code"
#define ATOM_NUMBER 			"atom_site.id"
#define ATOM_SEGMENT 			"atom_site.auth_asym_id"
#define ATOM_LABEL  			"atom_site.group_PDB"
#define ATOM_SYMBOL 			"atom_site.type_symbol"
#define ATOM_TOPTYPE			"atom_site.label_atom_id"
#define ATOM_RESIDUE 			"atom_site.label_comp_id" 
#define ATOM_CHAIN  			"atom_site.label_asym_id"
#define ATOM_ENTITY  			"atom_site.label_entity_id" 
#define ATOM_X  				"atom_site.Cartn_x" 
#define ATOM_Y  				"atom_site.Cartn_y" 
#define ATOM_Z  				"atom_site.Cartn_z" 
#define ATOM_FRACT_X  			"atom_site.fract_x" 
#define ATOM_FRACT_Y  			"atom_site.fract_y" 
#define ATOM_FRACT_Z  			"atom_site.fract_z" 
#define ATOM_OCCUPANCY  		"atom_site.occupancy" 
#define ATOM_BFACTOR 			"atom_site.B_iso_or_equiv" 
#define ATOM_CHARGE 			"atom_site.charge" 
#define ATOM_FOOTNOTE_ID 		"atom_site.footnote_id"

// Protein residue properties
#define RESPROP_CODON			"residue_property.codon"
#define RESPROP_CODE1 			"residue_property.code1"	
#define RESPROP_CODE3 			"residue_property.code3" 	
#define RESPROP_MASS 			"residue_property.mass"
#define RESPROP_VOLUME  		"residue_property.volume"
#define RESPROP_EXTENSION 		"residue_property.extension"
#define RESPROP_EXTSTDEV 		"residue_property.extstdev"
#define RESPROP_HYDROPHOBICITY  "residue_property.hydrophobicity"
#define RESPROP_CHARGE  		"residue_property.charge"
#define RESPROP_H  				"residue_property.H"
#define RESPROP_C  				"residue_property.C"
#define RESPROP_N  				"residue_property.N"
#define RESPROP_O  				"residue_property.O"
#define RESPROP_S  				"residue_property.S"
#define RESMATRIX1_CODE1 		"residue1_matrix.code1"	
#define RESMATRIX2_CODE1 		"residue2_matrix.code1"	
#define RESMATRIX_SUBSTITUTION  "residue_matrix.substitution"	
#define RESCOIL_CODE1			"residue_coil.code1"
#define RESCOIL_POS_A			"residue_coil.position_a"
#define RESCOIL_POS_B			"residue_coil.position_b"
#define RESCOIL_POS_C			"residue_coil.position_c"
#define RESCOIL_POS_D			"residue_coil.position_d"
#define RESCOIL_POS_E			"residue_coil.position_e"
#define RESCOIL_POS_F			"residue_coil.position_f"
#define RESCOIL_POS_G			"residue_coil.position_g"

// Nucleic acid and protein sequences
#define MOLECULE_NAME 			"molecule.name"
#define MOLECULE_LENGTH 		"molecule.length"
#define MOLECULE_SEQUENCE 		"molecule.sequence"

// Positions and orientations for multi-molecular solutions to fitting
#define MOLECULE_FILENAME		"molecule.file_name"
#define MOLECULE_MAP_FILE		"molecule.map_file_name"
#define MOLECULE_MAP_NUMBER		"molecule.map_number"
#define MOLECULE_X				"molecule.x"
#define MOLECULE_Y				"molecule.y"
#define MOLECULE_Z				"molecule.z"
#define MOLECULE_VIEW_X			"molecule.view_x"
#define MOLECULE_VIEW_Y			"molecule.view_y"
#define MOLECULE_VIEW_Z			"molecule.view_z"
#define MOLECULE_VIEW_ANGLE		"molecule.view_angle"
#define MOLECULE_FOM			"molecule.fom"
#define MOLECULE_SELECT			"molecule.select"


