/**
@file	rwmgEMX.cpp
@brief	Reads and writes micrograph exchange files
@author Bernard Heymann
@date	Created: 20130123
@date	Modified: 20200128
**/

#ifdef HAVE_XML
#include "rwxml.h"
#endif

#include "mg_processing.h"
#include "mg_img_proc.h"
#include "mg_tomography.h"
#include "rwmg.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"

#define	CONV_SWITCH	0	// Bad convention if 1

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
#ifdef HAVE_XML
int			emx_to_project(xmlDocPtr doc, Bproject* project, Bstring& filename);
int 		project_to_emx(Bproject* project, xmlDocPtr doc, int mg_select, int rec_select);

Vector3<int>	emx_get_integer_vector2(xmlNodePtr parent, const char* tag)
{
	Vector3<int>		vec;
	
	xmlNodePtr			node = xml_find_node(parent, tag);

	if ( !node ) return vec;
	
	vec[0] = xml_get_integer(node, "X");
	vec[1] = xml_get_integer(node, "Y");
	
	return vec;
}

Vector3<double>	emx_get_real_vector2(xmlNodePtr parent, const char* tag)
{
	Vector3<double>		vec;
	
	xmlNodePtr			node = xml_find_node(parent, tag);

	if ( !node ) return vec;
	
	vec[0] = xml_get_real(node, "X");
	vec[1] = xml_get_real(node, "Y");
	
	return vec;
}

int			emx_get_transformation(xmlNodePtr parent, Matrix3& mat, Vector3<double>& origin)
{
	xmlNodePtr			node = xml_find_node(parent, "transformationMatrix");
	
	if ( !node ) return -1;
	
	mat[0][0] = xml_get_real(node, "t11");
	mat[0][1] = xml_get_real(node, "t12");
	mat[0][2] = xml_get_real(node, "t13");
	mat[1][0] = xml_get_real(node, "t21");
	mat[1][1] = xml_get_real(node, "t22");
	mat[1][2] = xml_get_real(node, "t23");
	mat[2][0] = xml_get_real(node, "t31");
	mat[2][1] = xml_get_real(node, "t32");
	mat[2][2] = xml_get_real(node, "t33");

	origin[0] = xml_get_real(node, "t14");
	origin[1] = xml_get_real(node, "t24");
	origin[2] = xml_get_real(node, "t34");

	Matrix3				mattr = mat.transpose();

#if CONV_SWITCH
// 20150218 - convention T
	origin = mattr * origin;
	origin = -origin;
#else
// 20150323 - convention T-1
	mat = mattr;
#endif

	return 0;
}

xmlNodePtr	emx_set_real_vector2(xmlNodePtr parent, const char* tag, double v1, double v2, const char* unit)
{
	xmlNodePtr			node = xmlNewChild(parent, NULL, BAD_CAST tag, NULL);
	xmlNodePtr			node2 = xml_set_real(node, "X", v1, "%g");
	xmlNewProp(node2, BAD_CAST "unit", BAD_CAST unit);
	node2 = xml_set_real(node, "Y", v2, "%g");
	xmlNewProp(node2, BAD_CAST "unit", BAD_CAST unit);

	return node;
}

xmlNodePtr	emx_set_real_unit(xmlNodePtr parent, const char* tag, double v, const char* unit)
{
	xmlNodePtr			node = xml_set_real(parent, tag, v, "%g");
	
	if ( unit ) xmlNewProp(node, BAD_CAST "unit", BAD_CAST unit);

	return node;
}

xmlNodePtr	emx_set_transformation(xmlNodePtr parent, Matrix3 mat, Vector3<double> origin)
{
#if CONV_SWITCH
// 20150218 - convention T
	origin = mat * origin;
	origin = -origin;
#else
// 20150323 - convention T-1
	mat = mat.transpose();
#endif
	
	xmlNodePtr			node = xmlNewChild(parent, NULL, BAD_CAST "transformationMatrix", NULL);
	
	xml_set_real(node, "t11", mat[0][0], "%g");
	xml_set_real(node, "t12", mat[0][1], "%g");
	xml_set_real(node, "t13", mat[0][2], "%g");
	xml_set_real(node, "t21", mat[1][0], "%g");
	xml_set_real(node, "t22", mat[1][1], "%g");
	xml_set_real(node, "t23", mat[1][2], "%g");
	xml_set_real(node, "t31", mat[2][0], "%g");
	xml_set_real(node, "t32", mat[2][1], "%g");
	xml_set_real(node, "t33", mat[2][2], "%g");

	emx_set_real_unit(node, "t14", origin[0], "px");
	emx_set_real_unit(node, "t24", origin[1], "px");
	emx_set_real_unit(node, "t34", origin[2], "px");

	return node;
}

xmlNodePtr	emx_set_transformation(xmlNodePtr parent, View& view, Vector3<double> origin)
{
	Matrix3				mat = view.matrix();

	return emx_set_transformation(parent, mat, origin);
}

#endif

/**
@brief 	Reading micrograph parameters from XML files.
@param 	&filename	file name.
@param 	*project	initialized project structure.
@param 	&xsdfile	XML schema to validate against (can be empty).
@return int			error code (<0 means failure).
**/
#if defined(HAVE_XML)
int			read_project_emx(Bstring& filename, Bproject* project, Bstring& xsdfile)
{
	if ( filename.length() < 1 ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return -1;
	}
	
	int				err(0);
    xmlDocPtr		doc;
    xmlNodePtr		emx_node;
	
	detect_and_fix_carriage_return(filename.c_str());

	doc = xmlParseFile(filename.c_str());
	if (doc == NULL) return -1;

	if ( verbose & VERB_PROCESS )
		cout << "# Reading EMX file:               " << filename << endl;

	if ( xsdfile.length() ) xml_validate(doc, xsdfile);

	emx_node = xmlDocGetRootElement(doc);
	if ( !emx_node ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		xmlFreeDoc(doc);
		return -1;
	}
	
	if ( xmlStrcmp(emx_node->name, (const xmlChar *) "EMX") ) {
		cerr << "Error: The document " << filename << " is not an EMX file!" << endl;
		xmlFreeDoc(doc);
		return -1;
	}

	xmlChar*		content = xmlGetProp(emx_node, BAD_CAST "version");
	if ( verbose ) {
		if ( content ) cout << "EMX version " << content << endl;
		else cout << "EMX version ?" << endl;
	}
		
//		xmlDocDump(stderr, doc);
//		xmlElemDump(stderr, doc, node);
				
	err = emx_to_project(doc, project, filename);

	xmlFreeDoc(doc);

	xmlCleanupParser();
	
	project_calculate_angles(project);
	
	return err;
}

int			read_project_emx(Bstring* file_list, Bproject* project, Bstring& xsdfile)
{
	if ( !file_list ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return -1;
	}
	
	// Get the list of filenames
	int				err(0);
	Bstring*		thisfile = NULL;
	
	if ( verbose & VERB_DEBUG_STAR ) {
		cout << "DEBUG read_project_emx: EMX filenames:" << endl;
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		err = read_project_emx(*thisfile, project, xsdfile);
	}
	
	project_calculate_angles(project);
	
	return err;
}
#else
int			read_project_emx(Bstring& filename, Bproject* project, Bstring& xsdfile)
{
	cerr << "Error: EMX files are not supported!" << endl << endl;
	
	return -1;
}
#endif

/**
@brief 	Writing micrograph parameters to an EMX file.
@param 	&filename		file name.
@param 	*project		project structure.
@param 	mg_select		flag to only write selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return int				error code (<0 means failure).
**/
#if defined(HAVE_XML)
int			write_project_emx(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
	int					i, err(0);
	Bstring				name;

    xmlDocPtr			doc;
	
	Bfield*				field;
	Bmicrograph*		mg;
//	Breconstruction*	rec;
	char				format[32];
	
	snprintf(format, 32, "_%%0%dd.", project->split);
	
	if ( project->split ) {
		for ( i=0, field = project->field; field; field=field->next ) {
			for ( mg=field->mg; mg; mg=mg->next, i++ ) {
				if ( project->split == 9 && mg->id.length() > 0 ) {
					name = mg->id + ".xml";
					name =	name.replace(' ', '_');
				} else {
					name = filename.pre_rev('.') + Bstring(i, format) + filename.post_rev('.');
				}
				if ( verbose & VERB_PROCESS )
					cout << "# Writing file:                 " << name << endl;
				doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
				if ( doc == NULL ) {
					cerr << "Error: The XML document tree was not created!" << endl;
					return -1;
				}
				err = project_to_emx(project, doc, i, 0);
				if ( err == 0 ) {
					if ( verbose & VERB_PROCESS )
						cout << "# Writing EMX file:               " << filename << endl;
					err = xmlSaveFormatFile(name.c_str(), doc, 1);
				}
				xmlFreeDoc(doc);
			}
		}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "# Writing file:                 " << filename << endl;
		doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
		if ( !doc ) {
			cerr << "Error: The XML document tree was not created!" << endl;
			return -1;
		}
		err = project_to_emx(project, doc, 0, 0);
		if ( err == 0 ) {
			if ( verbose & VERB_PROCESS )
				cout << "# Writing EMX file:               " << filename << endl;
			err = xmlSaveFormatFile(filename.c_str(), doc, 1);
		}
		xmlFreeDoc(doc);
	}
	
	return err;
}
#else
int			write_project_emx(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
	cerr << "Error: EMX files are not supported!" << endl << endl;
	
	return -1;
}
#endif

#ifdef HAVE_XML


/*
@brief 	Converts micrograph data from an EMX document to a project structure.
@param 	xmlDocPtr doc	XML document pointer.
@param 	*project		project structure.
@param 	&filename		file name to resolve paths.
@return int				error code (<0 means failure).

	The function sets up the project hierarchy from a XML file.
	All angles given in the document are assumed to be in degrees, 
	and are converted here only to radians.

**/
int			emx_to_project(xmlDocPtr doc, Bproject* project, Bstring& filename)
{
	int					err(0);
	
	xmlNodePtr			emx_node = xmlDocGetRootElement(doc);
	xmlNodePtr			mg_node = NULL;
	xmlNodePtr			part_node = NULL;
	xmlNodePtr			part_mg_node = NULL;
	xmlNodePtr			node = NULL;

	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
//	Breconstruction*	rec = NULL;
	Bparticle*			part = NULL;
	
	Bstring				id, mg_file, part_file, prev_mg_file;
	long				i, pid, img_num;
	double				defU, defV;
	Vector3<double>		origin, center;
	Matrix3				mat;

	for ( i=0, node = emx_node->xmlChildrenNode; node; node = node->next ) {
		mg_file = 0;
		part_file = 0;
		if ( !xmlStrcmp(node->name, BAD_CAST "micrograph") ) {
			mg_node = node;
			mg_file = (char *) xmlGetProp(mg_node, BAD_CAST "fileName");
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG emx_to_project: mg=" << mg_file << endl;
			id = Bstring(++i, "%d");
			if ( mg_file.length() ) id = mg_file.pre_rev('.');
			if ( id.length() < 1 ) id = Bstring(i, "%d");;
//			if ( !field ) field = field_add(&project->field, id);
			if ( prev_mg_file != mg_file ) field = field_add(&project->field, id);
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG emx_to_project: id=" << id << endl;
			mg = micrograph_add(&field->mg, id);
			mg->block = i;
			mg->fmg = mg_file;
			mg->img_num = xml_get_integer_attribute(mg_node, "index") - 1;
			if ( mg->img_num < 0 ) {
				if ( verbose & VERB_FULL )
					cerr << "Warning: " << mg->id << " index number is less than 1" << endl;
				mg->img_num = 0;
			}
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG emx_to_project: mg: mg=" << mg_file << " img_num=" << mg->img_num << endl;

			center = micrograph_get_nominal_origin(mg);

			emx_get_transformation(mg_node, mg->matrix, origin);
			mg->origin = center + origin;
				
			mg->select = xml_get_integer(mg_node, "activeFlag");
			mg->fom = xml_get_real(mg_node, "fom");

			if ( xml_find_node(mg_node, "pixelSpacing") )
				mg->pixel_size = emx_get_real_vector2(mg_node, "pixelSpacing");
			
			if ( xml_check_for_node(mg_node, "defocusU") || xml_check_for_node(mg_node, "acceleratingVoltage") ) {
				mg->ctf = new CTFparam;
				mg->ctf->volt(1e3 * xml_get_real(mg_node, "acceleratingVoltage"));
				mg->ctf->amp_shift(asin(xml_get_real(mg_node, "amplitudeContrast")));
				mg->ctf->Cs(1e7 * xml_get_real(mg_node, "cs"));
				defU = 10*xml_get_real(mg_node, "defocusU");
				defV = 10*xml_get_real(mg_node, "defocusV");
				if ( defU < 100 ) defU *= 1e3;
				if ( defV < 100 ) defV *= 1e3;
				mg->ctf->defocus_average((defU + defV)/2);
				mg->ctf->defocus_deviation(fabs(defU - defV)/2);
				mg->ctf->astigmatism_angle(xml_get_real(node, "defocusUAngle") * M_PI/180.0);
			}
			prev_mg_file = mg_file;
		} else if ( !xmlStrcmp(node->name, BAD_CAST "particle") ) {
			part_node = node;
			part_file = (char *) xmlGetProp(part_node, BAD_CAST "fileName");
			pid = xml_get_integer_attribute(part_node, "index");
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG emx_to_project: part=" << part_file << " pid=" << pid << endl;
			
			part_mg_node = xml_find_node(part_node, "micrograph");
			if ( part_mg_node) {
				mg = NULL;
				mg_file = (char *) xmlGetProp(part_mg_node, BAD_CAST "fileName");
//				cout << "***" << mg_file << "***" << endl;
				if ( mg_file.length() ) {
					img_num = xml_get_integer_attribute(part_mg_node, "index") - 1;
					if ( img_num < 0 ) img_num = 0;
//					cout << "***" << mg_file << "***" << img_num << "***" << endl;
//					cout << "***" << field->mg->fmg << "***" << field->mg->img_num << "***" << endl;
					for ( field = project->field; field; field = field->next ) {
						for ( mg = field->mg; mg; mg = mg->next ) {
							if ( mg->fmg == mg_file && mg->img_num == img_num ) break;
						}
						if ( mg ) break;
					}
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG emx_to_project: part: mg=" << mg->fmg << " img_num=" << mg->img_num << endl;
				} else {
					cerr << "Warning: No micrograph file name specified!" << endl;
				}
			}
			
			if ( !mg ) {
				if ( verbose & VERB_FULL )
					cerr << "Warning: No micrograph for particle " << pid << " found!" << endl;
				id = Bstring(++i, "%d");
				if ( !field ) {
					if ( verbose )
						cout << "Generating a new field-of-view record with id: " << id << endl;
					field = field_add(&project->field, id);
				}
				if ( verbose )
					cout << "Generating a new micrograph record with id: " << id << endl;
				mg = micrograph_add(&field->mg, id);
			}
				
//			mg->fpart = (char *) xmlGetProp(part_node, BAD_CAST "fileName");

			if ( xml_find_node(part_node, "boxSize") )
				mg->box_size = emx_get_integer_vector2(part_node, "boxSize");
			
			if ( verbose & VERB_DEBUG )
				cout << "Particle " << pid << endl;
			part = particle_add(&mg->part, pid);

			part->fpart = part_file;

			if ( mg->box_size.volume() < 1 )
				mg->box_size = particle_get_box_size(part);

			if ( xml_find_node(part_node, "pixelSpacing") )
				part->pixel_size = emx_get_real_vector2(part_node, "pixelSpacing");
			
			if ( part->pixel_size[0] <= 0 )
				part->pixel_size = Vector3<double>(1,1,1);
			
			part->loc = emx_get_real_vector2(part_node, "centerCoord");

			emx_get_transformation(part_node, mat, origin);
			part->ori = origin + mg->box_size/2;
			part->view = View(mat);
			
			defU = 10*xml_get_real(part_node, "defocusU");
			defV = 10*xml_get_real(part_node, "defocusV");
			if ( defU < 100 ) defU *= 1e3;
			if ( defV < 100 ) defV *= 1e3;
			part->def = (defU + defV)/2;
			part->dev = fabs(defU - defV)/2;
			part->ast = xml_get_real(part_node, "defocusUAngle") * M_PI/180.0;
			
			part->sel = xml_get_integer(part_node, "activeFlag");
			part->fom[0] = xml_get_real(part_node, "fom");
		}
	}
	
	return err;
}


/*
@brief	Converts micrograph data from an project structure to a XML document.
@param	Bproject*		project parameter structure.
@param	xmlDocPtr doc	XML document.
@param 	mg_select		flag to only write selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return	int				error code (<0 means failure).

	The function packes new micrograph data back into an existing XML
	document from which the old project parameters were obtained.
	All angles used internally are in radians, and are converted here
	to degrees for output. The angles within the micrograph structure
	are left in radians.
	The stylesheet written into the document is "Bsoft_micrograph.xsl".
**/
int 		project_to_emx(Bproject* project, xmlDocPtr doc, int mg_select, int rec_select)
{
	int					err(0);
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Bparticle*			part = NULL;

	xmlNodePtr			root_node = xmlNewDocPI(doc, BAD_CAST "xml-stylesheet",
							BAD_CAST "href=\"emx.xsl\" type=\"text/xsl\"");
    xmlDocSetRootElement(doc, root_node);
	
	xmlNodePtr			emx_node = xmlNewNode(NULL, BAD_CAST "EMX");
	
	xmlAddSibling(root_node, emx_node);

	xmlNewProp(emx_node, BAD_CAST "version", BAD_CAST "1.0");
	
	Bstring				cmnt(" Written by Bsoft ");
	cmnt += BVERSION;
	xmlNodePtr			cmnt_node = xmlNewComment(BAD_CAST cmnt.c_str());
	xmlAddChild(emx_node, cmnt_node);
	
	xmlNodePtr			mg_node = NULL;
	xmlNodePtr			part_node = NULL;
	xmlNodePtr			node = NULL;
	
	double				angle;
	Vector3<double>		center;
	Bstring				filename;
	
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			if ( mg->fmg.length() ) filename = mg->fmg;
			else filename = "dummy.mrc";
			mg_node = xmlNewChild(emx_node, NULL, BAD_CAST "micrograph", NULL);
			xmlNewProp(mg_node, BAD_CAST "fileName", BAD_CAST filename.c_str());
			xml_set_integer_attribute(mg_node, "index", mg->img_num+1, "%d");
			
			center = micrograph_get_nominal_origin(mg);
			
			emx_set_real_vector2(mg_node, "pixelSpacing", mg->pixel_size[0], mg->pixel_size[1], "A/px");

			emx_set_transformation(mg_node, mg->matrix, mg->origin - center);

			if ( mg->ctf ) {
				emx_set_real_unit(mg_node, "acceleratingVoltage", mg->ctf->volt()/1e3, "kV");
				emx_set_real_unit(mg_node, "amplitudeContrast", sin(mg->ctf->amp_shift()), NULL);
				emx_set_real_unit(mg_node, "cs", mg->ctf->Cs()/1e7, "mm");
				emx_set_real_unit(mg_node, "defocusU", (mg->ctf->defocus_average() + mg->ctf->defocus_deviation())/10, "nm");
				emx_set_real_unit(mg_node, "defocusV", (mg->ctf->defocus_average() - mg->ctf->defocus_deviation())/10, "nm");
				angle = mg->ctf->astigmatism_angle()*180.0/M_PI;
				if ( angle < 0 ) angle += 180;
				emx_set_real_unit(mg_node, "defocusUAngle", angle, "deg");
			}
			
			xml_set_integer(mg_node, "activeFlag", mg->select, "%d");
			xml_set_real(mg_node, "fom", mg->fom, "%g");
		}
	}
	
	for ( field=project->field; field; field=field->next ) {
		for ( mg=field->mg; mg; mg=mg->next ) {
			if ( mg->fmg.length() ) filename = mg->fmg;
			else filename = "dummy.mrc";
			center = mg->box_size/2;
			for ( part=mg->part; part; part=part->next ) {
				part_node = xmlNewChild(emx_node, NULL, BAD_CAST "particle", NULL);
				if ( part->fpart.length() )
					xmlNewProp(part_node, BAD_CAST "fileName", BAD_CAST part->fpart.c_str());
				else
					xmlNewProp(part_node, BAD_CAST "fileName", BAD_CAST mg->fpart.c_str());
				xml_set_integer_attribute(part_node, "index", part->id, "%d");

				node = xmlNewChild(part_node, NULL, BAD_CAST "micrograph", NULL);
				xmlNewProp(node, BAD_CAST "fileName", BAD_CAST filename.c_str());
				xml_set_integer_attribute(node, "index", mg->img_num+1, "%d");

				emx_set_real_vector2(part_node, "pixelSpacing", part->pixel_size[0], part->pixel_size[1], "A/px");

				emx_set_real_vector2(part_node, "boxSize", mg->box_size[0], mg->box_size[1], "px");
				emx_set_real_vector2(part_node, "centerCoord", part->loc[0], part->loc[1], "px");

				emx_set_transformation(part_node, part->view, part->ori - center);

				if ( part->def ) {
					emx_set_real_unit(part_node, "defocusU", (part->def + part->dev)/10, "nm");
					emx_set_real_unit(part_node, "defocusV", (part->def - part->dev)/10, "nm");
					angle = part->ast*180.0/M_PI;
					if ( angle < 0 ) angle += 180;
					emx_set_real_unit(part_node, "defocusUAngle", angle, "deg");
				}
				
				xml_set_integer(part_node, "activeFlag", part->sel, "%d");
				xml_set_real(part_node, "fom", part->fom[0], "%g");
			}
		}
	}
	
	return err;
}

#endif

