/**
@file	mol_tags.h
@brief	All STAR file format tags for molecules and molecular parameters
@author 	Bernard Heymann
@date	Created: 20000419
@date	Modified: 20220929
**/

// Do not change the constant names because they are referenced in code
// Changes in the tag strings will not affect program execution
//  but old data files will be uninterpretable

// Entry identifier
#define ENTRY_ID				"pdbx_database_status.entry_id"

// Molecules/chains
#define MOL_ID					"struct_asym.id"

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

// Helix range parameters
#define HELIX_ID				"struct_conf.pdbx_PDB_helix_id"
#define HELIX_CHAIN				"struct_conf.beg_label_asym_id"
#define HELIX_RESIDUE1			"struct_conf.beg_label_comp_id"
#define HELIX_RESIDUE2			"struct_conf.end_label_comp_id"
#define HELIX_RESNUM1			"struct_conf.beg_label_seq_id"
#define HELIX_RESNUM2			"struct_conf.end_label_seq_id"

// Sheet range parameters
#define SHEET_ID				"struct_sheet_range.sheet_id"
#define SHEET_STRAND_ID			"struct_sheet_range.id"
#define SHEET_CHAIN				"struct_sheet_range.beg_label_asym_id"
#define SHEET_RESIDUE1			"struct_sheet_range.beg_label_comp_id"
#define SHEET_RESIDUE2			"struct_sheet_range.end_label_comp_id"
#define SHEET_RESNUM1			"struct_sheet_range.beg_label_seq_id"
#define SHEET_RESNUM2			"struct_sheet_range.end_label_seq_id"
#define SHEET_ORDER_ID			"struct_sheet_order.sheet_id"
#define SHEET_ORDER_STRAND_ID	"struct_sheet_order.range_id_2"
#define SHEET_ORDER_SENSE		"struct_sheet_order.sense"

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

// Chemical models
#define CHEMICAL_ID				"chem_comp.id"
#define CHEMICAL_NAME			"chem_comp.name"
#define CHEMICAL_PDB			"chem_comp.pdbx_type"
#define CHEMICAL_FOMRULA		"chem_comp.formula"
#define CHEMICAL_WEIGHT			"chem_comp.formula_weight"
#define CHEMICAL_ONE			"chem_comp.one_letter_code"
#define CHEMICAL_THREE			"chem_comp.three_letter_code"

// Chemical atoms
#define CHEMICAL_ATOM_RES		"chem_comp_atom.comp_id"
#define CHEMICAL_ATOM_ID		"chem_comp_atom.atom_id"
#define CHEMICAL_ATOM_ID_ALT	"chem_comp_atom.alt_atom_id"
#define CHEMICAL_ATOM_SYMBOL	"chem_comp_atom.type_symbol"
#define CHEMICAL_ATOM_CHARGE	"chem_comp_atom.charge"
#define CHEMICAL_ATOM_X			"chem_comp_atom.model_Cartn_x"
#define CHEMICAL_ATOM_Y			"chem_comp_atom.model_Cartn_y"
#define CHEMICAL_ATOM_Z			"chem_comp_atom.model_Cartn_z"
#define CHEMICAL_ATOM_XI		"chem_comp_atom.pdbx_model_Cartn_x_ideal"
#define CHEMICAL_ATOM_YI		"chem_comp_atom.pdbx_model_Cartn_y_ideal"
#define CHEMICAL_ATOM_ZI		"chem_comp_atom.pdbx_model_Cartn_z_ideal"
#define CHEMICAL_ATOM_IDP		"chem_comp_atom.pdbx_component_atom_id"
#define CHEMICAL_ATOM_RESP		"chem_comp_atom.pdbx_component_comp_id"
#define CHEMICAL_ATOM_NUMBER	"chem_comp_atom.pdbx_ordinal"

