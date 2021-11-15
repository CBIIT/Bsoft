/**
@file	rwmgXML.cpp
@brief	Reads and writes micrograph XML files
@author Bernard Heymann
@date	Created: 20050920
@date	Modified: 20200202
**/

#ifdef HAVE_XML
#include "rwxml.h"
#endif

#include "mg_tags.h"
#include "mg_processing.h"
#include "rwmg.h"
#include "file_util.h"
#include "linked_list.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
#ifdef HAVE_XML
int			xml_to_project(xmlDocPtr doc, Bproject* project, Bstring& filename);
int 		project_to_xml(Bproject* project, xmlDocPtr doc, int mg_select, int rec_select);
Bparticle*	particle_from_xml(xmlNodePtr node, Bparticle** partlist, FOMType fom_tag[NFOM]);
xmlNodePtr	particle_to_xml(Bparticle* part, xmlNodePtr parent, int euler_flag, int omega_flag, FOMType fom_tag[NFOM]);
#endif

/**
@brief 	Reading micrograph parameters from XML files.
@param 	&filename	file name.
@param 	*project	initialized project structure.
@return int					error code (<0 means failure).
**/
#if defined(HAVE_XML)
int			read_project_xml(Bstring& filename, Bproject* project)
{
	if ( filename.length() < 1 ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return -1;
	}
	
	int				err(0);
    xmlDocPtr		doc;
    xmlNodePtr		node;
	
	detect_and_fix_carriage_return(filename.c_str());

	doc = xmlParseFile(filename.c_str());
	if (doc == NULL) return -1;

	if ( verbose & VERB_PROCESS )
		cout << "# Reading XML file:               " << filename << endl;

	node = xmlDocGetRootElement(doc);
	if ( node == NULL ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		xmlFreeDoc(doc);
		return -1;
	}
	
	if ( xmlStrcmp(node->name, (const xmlChar *) PROJECT) ) {
		cerr << "Error: The document " << filename << " is not a micrograph parameter file!" << endl;
		xmlFreeDoc(doc);
		return -1;
	}
		
//		xmlDocDump(stderr, doc);
//		xmlElemDump(stderr, doc, node);
				
	err = xml_to_project(doc, project, filename);

	xmlFreeDoc(doc);

	xmlCleanupParser();
	
	return err;
}

int			read_project_xml(Bstring* file_list, Bproject* project)
{
	if ( !file_list ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return -1;
	}
	
	// Get the list of filenames
	int				err(0);
	Bstring*		thisfile = NULL;
	
	if ( verbose & VERB_DEBUG_STAR ) {
		cout << "DEBUG read_project_xml: XML filenames:" << endl;
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		err = read_project_xml(*thisfile, project);
	}
	
	return err;
}

Bparticle*	read_particle_xml(Bstring& filename, Bparticle** partlist, FOMType fom_tag[NFOM])
{
	if ( filename.length() < 1 ) {
		error_show("No file name found!", __FILE__, __LINE__);
		return NULL;
	}
	
    xmlDocPtr		doc;
    xmlNodePtr		node;
	
	doc = xmlParseFile(filename.c_str());
	if (doc == NULL) return NULL;

	if ( verbose & VERB_FULL )
		cout << "# Reading XML file:               " << filename << endl;

	node = xmlDocGetRootElement(doc);
	if ( node == NULL ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		xmlFreeDoc(doc);
		return NULL;
	}
	
	if ( xmlStrcmp(node->name, (const xmlChar *) PARTICLE) ) {
		cerr << "Error: The document " << filename << " is not a particle parameter file!" << endl;
		xmlFreeDoc(doc);
		return NULL;
	}
		
	Bparticle*		part = particle_from_xml(node, partlist, fom_tag);

	xmlFreeDoc(doc);

	return part;
}

#else
int			read_project_xml(Bstring& filename, Bproject* project)
{
	cerr << "Error: XML files are not supported!" << endl << endl;
	
	return -1;
}
#endif

/**
@brief 	Writing micrograph parameters to a XML file.
@param 	&filename		file name.
@param 	*project		project structure.
@param 	mg_select		flag to only write selected micrographs.
@param 	rec_select		flag to only convert selected reconstructions.
@return int				error code (<0 means failure).
**/
#if defined(HAVE_XML)
int			write_project_xml(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
	int					i, err(0);
	Bstring				name;

    xmlDocPtr			doc;
	
	Bfield*				field;
	Bmicrograph*		mg;
	Breconstruction*	rec;
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
				err = project_to_xml(project, doc, i, 0);
				if ( err == 0 )
					err = xmlSaveFormatFile(name.c_str(), doc, 1);
				xmlFreeDoc(doc);
			}
		}
		for ( i=0, rec = project->rec; rec; rec = rec->next ) {
			if ( project->split == 9 && rec->id.length() > 0 ) {
				name = rec->id + ".xml";
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
			err = project_to_xml(project, doc, 0, i);
			if ( err == 0 )
				err = xmlSaveFormatFile(name.c_str(), doc, 1);
			xmlFreeDoc(doc);
		}
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "# Writing file:                 " << filename << endl;
		doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
		if ( doc == NULL ) {
			cerr << "Error: The XML document tree was not created!" << endl;
			return -1;
		}
		err = project_to_xml(project, doc, 0, 0);
		if ( err == 0 )
			err = xmlSaveFormatFile(filename.c_str(), doc, 1);
		xmlFreeDoc(doc);
	}
	
	return err;
}

int			write_particle_xml(Bstring& filename, Bparticle* part, int euler_flag, int omega_flag, FOMType fom_tag[NFOM])
{
    xmlDocPtr		doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
	if (doc == NULL) return -1;

	if ( verbose & VERB_PROCESS )
		cout << "# Writing XML file:               " << filename << endl;

	xmlNodePtr		node = xmlNewDocPI(doc, BAD_CAST "xml-stylesheet",
						BAD_CAST "href=\"Bsoft_particle.xsl\" type=\"text/xsl\"");
    xmlDocSetRootElement(doc, node);
	
	xmlNodePtr		part_node = particle_to_xml(part, NULL, euler_flag, omega_flag, fom_tag);
//	cout << "part_node=" << part_node << endl;

	xmlAddSibling(node, part_node);
	
	if ( xmlSaveFormatFile(filename.c_str(), doc, 1) < 0 )
		cerr << "Error: XML file " << filename << " not written!" << endl;
	
	xmlFreeDoc(doc);

	return 0;
}
#else
int			write_project_xml(Bstring& filename, Bproject* project, int mg_select, int rec_select)
{
	cerr << "Error: XML files are not supported!" << endl << endl;
	
	return -1;
}
#endif

#ifdef HAVE_XML

CTFparam*		ctf_from_xml(xmlNodePtr node)
{
	if ( !xml_check_for_node(node, CTF_DEF_AVG) ) return NULL;
	
	Bstring				s;
	CTFparam*			ctf = new CTFparam;
	
	ctf->volt(xml_get_real(node, CTF_VOLTAGE));
	if ( xml_find_node(node, CTF_AMP_SHIFT) )
		ctf->amp_shift(xml_get_real(node, CTF_AMP_SHIFT));
	else
		ctf->amp_shift(asin(xml_get_real(node, CTF_AMP)));
	ctf->defocus_average(xml_get_real(node, CTF_DEF_AVG));
	ctf->defocus_deviation(xml_get_real(node, CTF_DEF_DEV));
	ctf->astigmatism_angle(xml_get_real(node, CTF_AST_ANG) * M_PI/180.0);
	ctf->Cs(xml_get_real(node, CTF_CS));
	ctf->Cc(xml_get_real(node, CTF_CC));
	ctf->alpha(xml_get_real(node, CTF_ALPHA));
	ctf->dE(xml_get_real(node, CTF_DE));
	ctf->zero(1);
	s = xml_copy_string(node, CTF_BASELINE);
	ctf->parse_baseline_equation(s);
	s = xml_copy_string(node, CTF_ENVELOPE);
	ctf->parse_envelope_equation(s);

	return ctf;
}

Bframe*		frame_from_xml(xmlNodePtr node, Bframe** framelist)
{
	xmlChar*		frame_id;
	Bframe*			frame = NULL;
	
	frame_id = xmlGetProp(node, BAD_CAST ID);
	frame = frame_add(framelist, atoi((char *) frame_id));
	xmlFree(frame_id);
	if ( !frame ) return frame;
	
	frame->shift[0] = xml_get_real(node, MICROGRAPH_FRAME_SHIFT_X);
	frame->shift[1] = xml_get_real(node, MICROGRAPH_FRAME_SHIFT_Y);
	frame->fom = xml_get_real(node, MICROGRAPH_FRAME_FOM);
	frame->sel = xml_get_integer(node, MICROGRAPH_FRAME_SELECT);

	return frame;
}

Bparticle*	particle_from_xml(xmlNodePtr node, Bparticle** partlist, FOMType fom_tag[NFOM])
{
	int 			f, euler_flag, omega_flag;
	Bstring			tag;
	xmlChar*		part_id;
	Bparticle*		part = NULL;
	Euler			euler;
	
	part_id = xmlGetProp(node, BAD_CAST ID);
	part = particle_add(partlist, atoi((char *) part_id));
	xmlFree(part_id);
	if ( !part ) return part;
	
	part->fpart = (char *) xmlGetProp(node, BAD_CAST PARTICLE_FILE);
	part->group = xml_get_integer(node, PARTICLE_GROUP);
	part->def = xml_get_real(node, PARTICLE_DEFOCUS);
	part->dev = xml_get_real(node, PARTICLE_DEF_DEV);
	part->ast = xml_get_real(node, PARTICLE_AST_ANG);
	part->mag = xml_get_real(node, PARTICLE_MAGNIF);
	part->loc[0] = xml_get_real(node, PARTICLE_X);
	part->loc[1] = xml_get_real(node, PARTICLE_Y);
	part->loc[2] = xml_get_real(node, PARTICLE_Z);
	part->pixel_size[0] = xml_get_real(node, PARTICLE_PIXEL_X);
	part->pixel_size[1] = xml_get_real(node, PARTICLE_PIXEL_Y);
	part->pixel_size[2] = xml_get_real(node, PARTICLE_PIXEL_Z);
	part->ori[0] = xml_get_real(node, PARTICLE_ORIGIN_X);
	part->ori[1] = xml_get_real(node, PARTICLE_ORIGIN_Y);
	part->ori[2] = xml_get_real(node, PARTICLE_ORIGIN_Z);
	euler_flag = xml_check_for_node(node, PARTICLE_THETA);
	if ( euler_flag == 0 ) {
		part->view[0] = xml_get_real(node, PARTICLE_VIEW_X);
		part->view[1] = xml_get_real(node, PARTICLE_VIEW_Y);
		part->view[2] = xml_get_real(node, PARTICLE_VIEW_Z);
		part->view[3] = xml_get_real(node, PARTICLE_VIEW_ANGLE)*M_PI/180.0;
	} else {
		omega_flag = xml_check_for_node(node, PARTICLE_OMEGA);
		if ( omega_flag )
			euler[0] = -xml_get_real(node, PARTICLE_OMEGA)*M_PI/180.0;
		else
			euler[0] = xml_get_real(node, PARTICLE_PSI)*M_PI/180.0;
		euler[1] = xml_get_real(node, PARTICLE_THETA)*M_PI/180.0;
		euler[2] = xml_get_real(node, PARTICLE_PHI)*M_PI/180.0;
		part->view = euler.view();
	}
	for ( f=0; f<NFOM; f++ ) if ( fom_tag[f] ) {
		tag = get_fom_tag(fom_tag[f]);
		part->fom[f] = xml_get_real(node, tag.c_str());
	}
	part->sel = xml_get_integer(node, PARTICLE_SELECT);

	return part;
}

Bfilament*	filament_from_xml(xmlNodePtr node, Bfilament** fillist)
{
	xmlChar*		fil_id;
	xmlNodePtr		fnode_node = NULL;
	Bfilament*		fil = NULL;
	Bfilnode*		fnode = NULL;
	
	fil_id = xmlGetProp(node, BAD_CAST ID);
	fil = filament_add(fillist, atoi((char *) fil_id));
	xmlFree(fil_id);
	for ( fnode_node = node->xmlChildrenNode; fnode_node; fnode_node = fnode_node->next ) {
		if ( !xmlStrcmp(fnode_node->name, BAD_CAST FILAMENT_NODE) ) {
			fil_id = xmlGetProp(node, BAD_CAST ID);
			fnode = filament_node_add(&fil->node, atoi((char *) fil_id));
			xmlFree(fil_id);
			fnode->loc[0] = xml_get_real(node, FILAMENT_NODE_X);
			fnode->loc[1] = xml_get_real(node, FILAMENT_NODE_Y);
			fnode->loc[2] = xml_get_real(node, FILAMENT_NODE_Z);
		}
	}

	return fil;
}

Bbadarea*	badarea_from_xml(xmlNodePtr node, Bbadarea** badlist)
{
	Bbadarea*		bad = (Bbadarea *) add_item((char **) badlist, sizeof(Bbadarea));
	
	bad->loc[0] = xml_get_real(node, PARTICLE_BAD_X);
	bad->loc[1] = xml_get_real(node, PARTICLE_BAD_Y);
	bad->loc[2] = xml_get_real(node, PARTICLE_BAD_Z);

	return bad;
}

Bmarker*	marker_from_xml(xmlNodePtr node, Bmarker** marklist)
{
	xmlChar*		mark_id;
	Bmarker*		mark = (Bmarker *) add_item((char **) marklist, sizeof(Bmarker));

	mark_id = xmlGetProp(node, BAD_CAST ID);
	mark->id = atoi((char *) mark_id);
	xmlFree(mark_id);
	mark->fom = 1;
	mark->loc[0] = xml_get_real(node, MARKER_X);
	mark->loc[1] = xml_get_real(node, MARKER_Y);
	mark->loc[2] = xml_get_real(node, MARKER_Z);
	mark->err[0] = xml_get_real(node, MARKER_ERROR_X);
	mark->err[1] = xml_get_real(node, MARKER_ERROR_Y);
	mark->err[2] = xml_get_real(node, MARKER_ERROR_Z);
	mark->res = xml_get_real(node, MARKER_RESIDUAL);
	mark->fom = xml_get_real(node, MARKER_FOM);
	mark->sel = xml_get_integer(node, MARKER_SELECT);

	return mark;
}

Bstrucfac*	strucfac_from_xml(xmlNodePtr node, Bstrucfac** sflist)
{
	Bstrucfac*		sf = (Bstrucfac *) add_item((char **) sflist, sizeof(Bstrucfac));
	
	sf->fom = sf->sel = 1;
	sf->index[0] = xml_get_integer(node, REFLEX_H);
	sf->index[1] = xml_get_integer(node, REFLEX_K);
	sf->index[2] = xml_get_integer(node, REFLEX_L);
	sf->loc[0] = xml_get_real(node, REFLEX_X);
	sf->loc[1] = xml_get_real(node, REFLEX_Y);
	sf->loc[2] = xml_get_real(node, REFLEX_Z);
	sf->amp = xml_get_real(node, REFLEX_AMP);
	sf->sigamp = xml_get_real(node, REFLEX_SIGAMP);
	sf->phi = xml_get_real(node, REFLEX_PHI);
	sf->sigphi = xml_get_real(node, REFLEX_SIGPHI);
	sf->fom = xml_get_real(node, REFLEX_FOM);
	sf->sel = xml_get_integer(node, REFLEX_STATUS);

	return sf;
}

Blayerline*	layerline_from_xml(xmlNodePtr node, Blayerline** linelist)
{
	Blayerline*		line = (Blayerline *) add_item((char **) linelist, sizeof(Blayerline));
	
	line->fom = line->sel = 1;
	line->number = xml_get_integer(node, LAYERLINE_NUMBER);
	line->order = xml_get_integer(node, LAYERLINE_ORDER);
	line->distance = xml_get_real(node, LAYERLINE_DISTANCE);
	line->freq = xml_get_real(node, LAYERLINE_FREQ);
	line->amp = xml_get_real(node, LAYERLINE_AMP);
	line->fom = xml_get_real(node, LAYERLINE_FOM);
	line->sel = xml_get_integer(node, LAYERLINE_SELECT);

	return line;
}

/*
@brief 	Converts micrograph data from a XML document to a project structure.
@param 	xmlDocPtr doc	XML document pointer.
@param 	*project		project structure.
@param 	&filename		file name to resolve paths.
@return int				error code (<0 means failure).

	The function sets up the project hierarchy from a XML file.
	All angles given in the XML document are assumed to be in degrees, 
	and are converted here only to radians.

**/
int			xml_to_project(xmlDocPtr doc, Bproject* project, Bstring& filename)
{
	int					err(0);
	
	Bfield*				field = project->field;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;

	int					i, f;
	xmlChar				*xstring, *id;
	Bstring				bid;
	Bstring				s;

    xmlNodePtr			project_node = xmlDocGetRootElement(doc);
	xmlNodePtr			child_node = NULL;
	xmlNodePtr			mg_node = NULL;
	xmlNodePtr			node = NULL;
	
	// Make sure FOM tags are properly set up
	if ( xml_check_for_node_in_tree(project_node, PARTICLE_FOM) ) project->fom_tag[0] = FOM;
	if ( xml_check_for_node_in_tree(project_node, PARTICLE_FOM_CV) ) project->fom_tag[1] = FOM_CV;
	if ( xml_check_for_node_in_tree(project_node, PARTICLE_HANDA_FOM) ) project->fom_tag[1] = FOM_HAND_A;
	if ( xml_check_for_node_in_tree(project_node, PARTICLE_HANDB_FOM) ) project->fom_tag[2] = FOM_HAND_B;
	
	for ( f=i=0; f<NFOM; f++ ) i += project->fom_tag[f];
	if ( i < 2 ) project->fom_tag[0] = FOM;
	
	for ( i=0, child_node = project_node->xmlChildrenNode; child_node; child_node = child_node->next ) {
		if ( project->comment.empty() ) {
			if ( !xmlStrcmp(child_node->name, BAD_CAST COMMENT) )
				project->comment = (char *) xmlNodeGetContent(child_node);
		}
		if ( !project->reference ) {
			if ( !xmlStrcmp(child_node->name, BAD_CAST MAP_REFERENCE) ) {
				xstring = xmlNodeGetContent(child_node);
				string_add(&project->reference, (char *) xstring);
				xmlFree(xstring);
			}
		}
        if ( !xmlStrcmp(child_node->name, BAD_CAST FIELD) ) {
			id = xmlGetProp(child_node, BAD_CAST ID);
			bid = (char *) id;
			for ( field=project->field; field; field=field->next )
				if ( field->id == bid ) break;
			if ( !field ) field = field_add(&project->field, bid);
			xmlFree(id);
			for ( mg_node = child_node->xmlChildrenNode; mg_node; mg_node = mg_node->next ) {
				if ( !xmlStrcmp(mg_node->name, BAD_CAST MICROGRAPH) ) {
					id = xmlGetProp(mg_node, BAD_CAST ID);
					bid = (char *) id;
					mg = micrograph_add(&field->mg, bid);
					xmlFree(id);
					mg->block = i++;
					mg->fmg = xml_get_string(mg_node, MICROGRAPH_FILE);
					mg->fframe = xml_get_string(mg_node, MICROGRAPH_FRAMES_FILE);
					mg->fpart = xml_get_string(mg_node, PARTICLE_FILE);
					mg->ffil = xml_get_string(mg_node, FILAMENT_FILE);
					mg->fft = xml_get_string(mg_node, MICROGRAPH_TRANSFORM_FILE);
					mg->fps = xml_get_string(mg_node, MICROGRAPH_POWERSPEC_FILE);
					mg->img_num = xml_get_integer(mg_node, MICROGRAPH_NUMBER);
					mg->select = xml_get_integer(mg_node, MICROGRAPH_SELECT);
					mg->fom = xml_get_real(mg_node, MICROGRAPH_FOM);
					mg->magnification = xml_get_real(mg_node, MICROGRAPH_MAGNIFICATION);
					mg->sampling = xml_get_real(mg_node, MICROGRAPH_SAMPLING);
					mg->pixel_size[0] = xml_get_real(mg_node, MICROGRAPH_PIXEL_X);
					mg->pixel_size[1] = xml_get_real(mg_node, MICROGRAPH_PIXEL_Y);
					mg->frame_pixel_size[0] = xml_get_real(mg_node, MICROGRAPH_FRAME_PIXEL_X);
					mg->frame_pixel_size[1] = xml_get_real(mg_node, MICROGRAPH_FRAME_PIXEL_Y);
					mg->dose = xml_get_real(mg_node, MICROGRAPH_DOSE);
					mg->exposure = xml_get_real(mg_node, MICROGRAPH_EXPOSURE);
					mg->intensity = xml_get_real(mg_node, MICROGRAPH_INTENSITY);
					mg->wri = xml_get_real(mg_node, MICROGRAPH_WATER_RING);
					mg->tilt_axis = xml_get_real(mg_node, MICROGRAPH_TILT_AXIS);
					mg->tilt_angle = xml_get_real(mg_node, MICROGRAPH_TILT_ANGLE);
					mg->rot_angle = xml_get_real(mg_node, MICROGRAPH_ROT_ANGLE);
					mg->origin[0] = xml_get_real(mg_node, MICROGRAPH_ORIGIN_X);
					mg->origin[1] = xml_get_real(mg_node, MICROGRAPH_ORIGIN_Y);
					mg->origin[2] = xml_get_real(mg_node, MICROGRAPH_ORIGIN_Z);
					mg->scale[0] = xml_get_real(mg_node, MICROGRAPH_SCALE_X);
					mg->scale[1] = xml_get_real(mg_node, MICROGRAPH_SCALE_Y);
					mg->scale[2] = xml_get_real(mg_node, MICROGRAPH_SCALE_Z);
					mg->matrix[0][0] = xml_get_real(mg_node, MICROGRAPH_MATRIX_1_1);
					mg->matrix[0][1] = xml_get_real(mg_node, MICROGRAPH_MATRIX_1_2);
					mg->matrix[0][2] = xml_get_real(mg_node, MICROGRAPH_MATRIX_1_3);
					mg->matrix[1][0] = xml_get_real(mg_node, MICROGRAPH_MATRIX_2_1);
					mg->matrix[1][1] = xml_get_real(mg_node, MICROGRAPH_MATRIX_2_2);
					mg->matrix[1][2] = xml_get_real(mg_node, MICROGRAPH_MATRIX_2_3);
					mg->matrix[2][0] = xml_get_real(mg_node, MICROGRAPH_MATRIX_3_1);
					mg->matrix[2][1] = xml_get_real(mg_node, MICROGRAPH_MATRIX_3_2);
					mg->matrix[2][2] = xml_get_real(mg_node, MICROGRAPH_MATRIX_3_3);
					mg->hvec[0] = xml_get_real(mg_node, MICROGRAPH_HVEC_X);
					mg->hvec[1] = xml_get_real(mg_node, MICROGRAPH_HVEC_Y);
					mg->hvec[2] = xml_get_real(mg_node, MICROGRAPH_HVEC_Z);
					mg->kvec[0] = xml_get_real(mg_node, MICROGRAPH_KVEC_X);
					mg->kvec[1] = xml_get_real(mg_node, MICROGRAPH_KVEC_Y);
					mg->kvec[2] = xml_get_real(mg_node, MICROGRAPH_KVEC_Z);
					mg->lvec[0] = xml_get_real(mg_node, MICROGRAPH_LVEC_X);
					mg->lvec[1] = xml_get_real(mg_node, MICROGRAPH_LVEC_Y);
					mg->lvec[2] = xml_get_real(mg_node, MICROGRAPH_LVEC_Z);
					mg->helix_axis = xml_get_real(mg_node, MICROGRAPH_HELIX_AXIS);
					mg->helix_rise = xml_get_real(mg_node, MICROGRAPH_HELIX_RISE);
					mg->helix_angle = xml_get_real(mg_node, MICROGRAPH_HELIX_ANGLE);
					mg->helix_radius = xml_get_real(mg_node, MICROGRAPH_HELIX_RADIUS);
					if ( mg->sampling < 1e3 ) mg->sampling *= 1e4;
					if ( mg->magnification < 1e3 ) mg->magnification *= 1e3;
					mg->ctf = ctf_from_xml(mg_node);
					mg->tilt_axis *= M_PI/180.0;
					mg->tilt_angle *= M_PI/180.0;
					mg->rot_angle  *= M_PI/180.0;
					if ( xml_check_for_node(mg_node, PARTICLE_BOX_RADIUS) )
						mg->box_size[0] = mg->box_size[1] = mg->box_size[2] = 
							2 * xml_get_integer(mg_node, PARTICLE_BOX_RADIUS);
					if ( xml_check_for_node(mg_node, PARTICLE_BOX_RADIUS_X) ) {
						mg->box_size[0] = 2 * xml_get_integer(mg_node, PARTICLE_BOX_RADIUS_X);
						mg->box_size[1] = 2 * xml_get_integer(mg_node, PARTICLE_BOX_RADIUS_Y);
						mg->box_size[2] = 2 * xml_get_integer(mg_node, PARTICLE_BOX_RADIUS_Z);
					}
					if ( xml_check_for_node(mg_node, PARTICLE_BOX_SIZE) )
						mg->box_size[0] = mg->box_size[1] = mg->box_size[2] =
							xml_get_integer(mg_node, PARTICLE_BOX_SIZE);
					if ( xml_check_for_node(mg_node, PARTICLE_BOX_SIZE_X) ) {
						mg->box_size[0] = xml_get_integer(mg_node, PARTICLE_BOX_SIZE_X);
						mg->box_size[1] = xml_get_integer(mg_node, PARTICLE_BOX_SIZE_Y);
						mg->box_size[2] = xml_get_integer(mg_node, PARTICLE_BOX_SIZE_Z);
					}
					mg->filament_width = xml_get_real(mg_node, FILAMENT_WIDTH);
					mg->fil_node_radius = xml_get_real(mg_node, FILNODE_RADIUS);
					mg->bad_radius = xml_get_real(mg_node, PARTICLE_BAD_RADIUS);
					mg->sf_radius = xml_get_real(mg_node, REFLEX_RADIUS);
					mg->mark_radius = xml_get_real(mg_node, MARKER_RADIUS);
					for ( node = mg_node->xmlChildrenNode; node; node = node->next ) {
						if ( !xmlStrcmp(node->name, BAD_CAST MICROGRAPH_FRAME) )
							frame_from_xml(node, &mg->frame);
						if ( !xmlStrcmp(node->name, BAD_CAST PARTICLE) )
							particle_from_xml(node, &mg->part, project->fom_tag);
						if ( !xmlStrcmp(node->name, BAD_CAST FILAMENT) )
							filament_from_xml(node, &mg->fil);
						if ( !xmlStrcmp(node->name, BAD_CAST PARTICLE_BAD) )
							badarea_from_xml(node, &mg->bad);
						if ( !xmlStrcmp(node->name, BAD_CAST MARKER) )
							marker_from_xml(node, &mg->mark);
						if ( !xmlStrcmp(node->name, BAD_CAST REFLEX) )
							strucfac_from_xml(node, &mg->sf);
						if ( !xmlStrcmp(node->name, BAD_CAST LAYERLINE) )
							layerline_from_xml(node, &mg->layer);
					}
				}
			}
		} else if ( !xmlStrcmp(child_node->name, BAD_CAST MAP) ) {
			id = xmlGetProp(child_node, BAD_CAST ID);
			bid = (char *) id;
			rec = reconstruction_add(&rec, bid);
			if ( !project->rec ) project->rec = rec;
			xmlFree(id);
			rec->block = i++;
			rec->select = 1;
			rec->frec = xml_get_string(child_node, MAP_RECONSTRUCTION);
			rec->fpart = xml_get_string(child_node, PARTICLE_FILE);
			rec->ffil = xml_get_string(child_node, FILAMENT_FILE);
			rec->fft = xml_get_string(child_node, MAP_TRANSFORM_FILE);
			rec->fps = xml_get_string(child_node, MAP_POWERSPEC_FILE);
			s = xml_copy_string(child_node, MAP_MODEL);
			if ( s.length() ) string_add(&rec->model, s);
			rec->select = xml_get_integer(child_node, MAP_SELECT);
			rec->fom = xml_get_real(child_node, MAP_FOM);
			rec->origin[0] = xml_get_real(child_node, MAP_ORIGIN_X);
			rec->origin[1] = xml_get_real(child_node, MAP_ORIGIN_Y);
			rec->origin[2] = xml_get_real(child_node, MAP_ORIGIN_Z);
			rec->scale[0] = xml_get_real(child_node, MAP_SCALE_X);
			rec->scale[1] = xml_get_real(child_node, MAP_SCALE_Y);
			rec->scale[2] = xml_get_real(child_node, MAP_SCALE_Z);
			rec->voxel_size[0] = xml_get_real(child_node, MAP_VOXEL_SIZE_X);
			rec->voxel_size[1] = xml_get_real(child_node, MAP_VOXEL_SIZE_Y);
			rec->voxel_size[2] = xml_get_real(child_node, MAP_VOXEL_SIZE_Z);
			rec->symmetry = xml_get_string(child_node, MAP_SYMMETRY);
			if ( xml_check_for_node(child_node, CTF_DEF_AVG) )
				rec->ctf = ctf_from_xml(child_node);
			if ( xml_check_for_node(child_node, PARTICLE_BOX_RADIUS) )
				rec->box_size[0] = rec->box_size[1] = rec->box_size[2] =
					2 * xml_get_integer(child_node, PARTICLE_BOX_RADIUS);
			if ( xml_check_for_node(child_node, PARTICLE_BOX_RADIUS_X) ) {
				rec->box_size[0] = 2 * xml_get_integer(child_node, PARTICLE_BOX_RADIUS_X);
				rec->box_size[1] = 2 * xml_get_integer(child_node, PARTICLE_BOX_RADIUS_Y);
				rec->box_size[2] = 2 * xml_get_integer(child_node, PARTICLE_BOX_RADIUS_Z);
			}
			if ( xml_check_for_node(child_node, PARTICLE_BOX_SIZE) )
				rec->box_size[0] = rec->box_size[1] = rec->box_size[2] =
					xml_get_integer(child_node, PARTICLE_BOX_SIZE);
			if ( xml_check_for_node(child_node, PARTICLE_BOX_SIZE_X) ) {
				rec->box_size[0] = xml_get_integer(child_node, PARTICLE_BOX_SIZE_X);
				rec->box_size[1] = xml_get_integer(child_node, PARTICLE_BOX_SIZE_Y);
				rec->box_size[2] = xml_get_integer(child_node, PARTICLE_BOX_SIZE_Z);
			}
			rec->filament_width = xml_get_real(child_node, FILAMENT_WIDTH);
			rec->fil_node_radius = xml_get_real(child_node, FILNODE_RADIUS);
			rec->bad_radius = xml_get_real(child_node, PARTICLE_BAD_RADIUS);
			rec->sf_radius = xml_get_real(child_node, REFLEX_RADIUS);
			rec->mark_radius = xml_get_real(child_node, MARKER_RADIUS);
			rec->view[0] = xml_get_real(child_node, MAP_VIEW_X);
			rec->view[1] = xml_get_real(child_node, MAP_VIEW_Y);
			rec->view[2] = xml_get_real(child_node, MAP_VIEW_Z);
			rec->view[3] = xml_get_real(child_node, MAP_VIEW_ANGLE)*M_PI/180.0;
			for ( node = child_node->xmlChildrenNode; node; node = node->next ) {
				if ( !xmlStrcmp(node->name, BAD_CAST PARTICLE) )
					particle_from_xml(node, &rec->part, project->fom_tag);
				if ( !xmlStrcmp(node->name, BAD_CAST FILAMENT) )
					filament_from_xml(node, &rec->fil);
				if ( !xmlStrcmp(node->name, BAD_CAST PARTICLE_BAD) )
					badarea_from_xml(node, &rec->bad);
				if ( !xmlStrcmp(node->name, BAD_CAST MARKER) )
					marker_from_xml(node, &rec->mark);
				if ( !xmlStrcmp(node->name, BAD_CAST REFLEX) )
					strucfac_from_xml(node, &rec->sf);
			}
		} else if ( !xmlStrcmp(child_node->name, BAD_CAST CLASS_AVERAGE) ) {
			for ( node = child_node->xmlChildrenNode; node; node = node->next )
				if ( !xmlStrcmp(node->name, BAD_CAST PARTICLE) )
					particle_from_xml(node, &project->class_avg, project->fom_tag);
		}
	}
	
	return err;
}

int			ctf_to_xml(CTFparam* ctf, xmlNodePtr parent)
{
	if ( !ctf ) return -1;
	
	Bstring		s;
	
	xml_set_real(parent, CTF_VOLTAGE, ctf->volt(), "%g");
	xml_set_real(parent, CTF_CS, ctf->Cs(), "%g");
	xml_set_real(parent, CTF_CC, ctf->Cc(), "%g");
	xml_set_real(parent, CTF_ALPHA, ctf->alpha(), "%g");
	xml_set_real(parent, CTF_DE, ctf->dE(), "%g");
	xml_set_real(parent, CTF_AMP_SHIFT, ctf->amp_shift(), "%g");
	xml_set_real(parent, CTF_DEF_AVG, ctf->defocus_average(), "%g");
	xml_set_real(parent, CTF_DEF_DEV, ctf->defocus_deviation(), "%g");
	xml_set_real(parent, CTF_AST_ANG, ctf->astigmatism_angle()*180.0/M_PI, "%g");
	xml_set_real(parent, CTF_ZERO, ctf->zero(1), "%g");
	if ( ctf->baseline_type() ) {
		s = ctf->baseline_equation();
		xmlNewChild(parent, NULL, BAD_CAST CTF_BASELINE, BAD_CAST s.c_str());
	}
	if ( ctf->envelope(0) ) {
		s = ctf->envelope_equation();
		xmlNewChild(parent, NULL, BAD_CAST CTF_ENVELOPE, BAD_CAST s.c_str());
	}
	
	return 1;
}

int			frame_to_xml(Bframe* framelist, xmlNodePtr parent)
{
	int				nframe(0);
	Bframe*			frame = NULL;
	xmlNodePtr		frame_node = NULL;

	for ( frame=framelist; frame; frame=frame->next, nframe++ ) {
		frame_node = xmlNewChild(parent, NULL, BAD_CAST MICROGRAPH_FRAME, NULL);
		xml_set_integer_attribute(frame_node, ID, frame->id, "%d");
		xml_set_real(frame_node, MICROGRAPH_FRAME_SHIFT_X, frame->shift[0], "%7.3f");
		xml_set_real(frame_node, MICROGRAPH_FRAME_SHIFT_Y, frame->shift[1], "%7.3f");
		xml_set_real(frame_node, MICROGRAPH_FRAME_FOM, frame->fom, "%7.4lf");
		xml_set_integer(frame_node, MICROGRAPH_FRAME_SELECT, frame->sel, "%4d");
	}
	
	return nframe;
}

xmlNodePtr	particle_to_xml(Bparticle* part, xmlNodePtr parent, int euler_flag, int omega_flag, FOMType fom_tag[NFOM])
{
	int				f;
	Bstring			tag;
	Euler			euler;

	xmlNodePtr		part_node = NULL;

	if ( parent ) part_node = xmlNewChild(parent, NULL, BAD_CAST PARTICLE, NULL);
	else part_node = xmlNewNode(NULL, BAD_CAST PARTICLE);
//	cout << "part_node=" << part_node << endl;

	if ( part->fpart.length() )
		xmlNewProp(part_node, BAD_CAST PARTICLE_FILE, BAD_CAST part->fpart.c_str());
	xml_set_integer_attribute(part_node, ID, part->id, "%d");
	xml_set_integer(part_node, PARTICLE_GROUP, part->group, "%d");
	
	if ( part->def ) {
		xml_set_real(part_node, PARTICLE_DEFOCUS, part->def, "%g");
		xml_set_real(part_node, PARTICLE_DEF_DEV, part->dev, "%g");
		xml_set_real(part_node, PARTICLE_AST_ANG, part->ast, "%g");
	}
	
	xml_set_real(part_node, PARTICLE_MAGNIF, part->mag, "%g");
	xml_set_real(part_node, PARTICLE_X, part->loc[0], "%7.2f");
	xml_set_real(part_node, PARTICLE_Y, part->loc[1], "%7.2f");
	xml_set_real(part_node, PARTICLE_Z, part->loc[2], "%7.2f");
	xml_set_real(part_node, PARTICLE_PIXEL_X, part->pixel_size[0], "%7.4f");
	xml_set_real(part_node, PARTICLE_PIXEL_Y, part->pixel_size[1], "%7.4f");
	xml_set_real(part_node, PARTICLE_PIXEL_Z, part->pixel_size[2], "%7.4f");
	xml_set_real(part_node, PARTICLE_ORIGIN_X, part->ori[0], "%7.3f");
	xml_set_real(part_node, PARTICLE_ORIGIN_Y, part->ori[1], "%7.3f");
	xml_set_real(part_node, PARTICLE_ORIGIN_Z, part->ori[2], "%7.3f");
	
	if ( euler_flag < 1 ) {
		xml_set_real(part_node, PARTICLE_VIEW_X, part->view[0], "%7.4lf");
		xml_set_real(part_node, PARTICLE_VIEW_Y, part->view[1], "%7.4lf");
		xml_set_real(part_node, PARTICLE_VIEW_Z, part->view[2], "%7.4lf");
		xml_set_real(part_node, PARTICLE_VIEW_ANGLE, part->view.angle()*180.0/M_PI, "%7.2lf");
	} else {
		euler = Euler(part->view);
		if ( omega_flag ) euler[0] = -euler.psi();
		xml_set_real(part_node, PARTICLE_PHI, euler.phi()*180.0/M_PI, "%7.2f");
		xml_set_real(part_node, PARTICLE_THETA, euler.theta()*180.0/M_PI, "%7.2f");
		xml_set_real(part_node, PARTICLE_PSI, euler.psi()*180.0/M_PI, "%7.2f");
	}
	
	for ( f=0; f<NFOM; f++ ) if ( fom_tag[f] ) {
		tag = get_fom_tag(fom_tag[f]);
		if ( tag.length() > 1 ) {
			xml_set_real(part_node, tag.c_str(), part->fom[f], "%7.4f");
		}
	}
	
	xml_set_integer(part_node, PARTICLE_SELECT, part->sel, "%4d");

	return part_node;
}

int			filament_to_xml(Bfilament* fillist, xmlNodePtr parent)
{
	int				nfil(0);
	Bfilament*		fil;
	Bfilnode*		fnode = NULL;
	xmlNodePtr		fil_node = NULL;
	xmlNodePtr		fnode_node = NULL;
	
	for ( fil=fillist; fil; fil=fil->next, nfil++ ) {
		fil_node = xmlNewChild(parent, NULL, BAD_CAST FILAMENT, NULL);
		xml_set_integer_attribute(fil_node, ID, fil->id, "%5d");
		for ( fnode=fil->node; fnode; fnode=fnode->next ) {
			fnode_node = xmlNewChild(fil_node, NULL, BAD_CAST FILAMENT_NODE, NULL);
			xml_set_integer_attribute(fnode_node, ID, fnode->id, "%5d");
			xml_set_real(fnode_node, FILAMENT_NODE_X, fnode->loc[0], "%7.2f");
			xml_set_real(fnode_node, FILAMENT_NODE_Y, fnode->loc[1], "%7.2f");
			xml_set_real(fnode_node, FILAMENT_NODE_Z, fnode->loc[2], "%7.2f");
		}
	}

	return nfil;
}

int			badarea_to_xml(Bbadarea* badlist, xmlNodePtr parent)
{
	int				nbad(0);
	Bbadarea*		bad;
	xmlNodePtr		bad_node = NULL;
	
	for ( bad=badlist; bad; bad=bad->next, nbad++ ) {
		bad_node = xmlNewChild(parent, NULL, BAD_CAST PARTICLE_BAD, NULL);
		xml_set_real(bad_node, PARTICLE_BAD_X, bad->loc[0], "%7.2f");
		xml_set_real(bad_node, PARTICLE_BAD_Y, bad->loc[1], "%7.2f");
		xml_set_real(bad_node, PARTICLE_BAD_Z, bad->loc[2], "%7.2f");
	}

	return nbad;
}

int			marker_to_xml(Bmarker* marklist, xmlNodePtr parent)
{
	int				nmark(0);
	Bmarker*		mark;
	xmlNodePtr		mark_node = NULL;
	
	for ( mark=marklist; mark; mark=mark->next, nmark++ ) {
		mark_node = xmlNewChild(parent, NULL, BAD_CAST MARKER, NULL);
		xml_set_integer_attribute(mark_node, ID, mark->id, "%5d");
		xml_set_real(mark_node, MARKER_X, mark->loc[0], "%7.2f");
		xml_set_real(mark_node, MARKER_Y, mark->loc[1], "%7.2f");
		xml_set_real(mark_node, MARKER_Z, mark->loc[2], "%7.2f");
		xml_set_real(mark_node, MARKER_ERROR_X, mark->err[0], "%7.2f");
		xml_set_real(mark_node, MARKER_ERROR_Y, mark->err[1], "%7.2f");
		xml_set_real(mark_node, MARKER_ERROR_Z, mark->err[2], "%7.2f");
		xml_set_real(mark_node, MARKER_RESIDUAL, mark->res, "%7.4f");
		xml_set_real(mark_node, MARKER_FOM, mark->fom, "%7.4f");
		xml_set_integer(mark_node, MARKER_SELECT, mark->sel, "%d");
	}

	return nmark;
}

int			strucfac_to_xml(Bstrucfac* sflist, xmlNodePtr parent)
{
	int				nsf(0);
	Bstrucfac*		sf = NULL;
	xmlNodePtr		sf_node = NULL;
	
	for ( sf=sflist; sf; sf=sf->next, nsf++ ) {
		sf_node = xmlNewChild(parent, NULL, BAD_CAST REFLEX, NULL);
		xml_set_integer(sf_node, REFLEX_H, sf->index[0], "%4d");
		xml_set_integer(sf_node, REFLEX_K, sf->index[1], "%4d");
		xml_set_integer(sf_node, REFLEX_L, sf->index[2], "%4d");
		xml_set_real(sf_node, REFLEX_X, sf->loc[0], "%7.2f");
		xml_set_real(sf_node, REFLEX_Y, sf->loc[1], "%7.2f");
		xml_set_real(sf_node, REFLEX_Z, sf->loc[2], "%7.2f");
		xml_set_real(sf_node, REFLEX_AMP, sf->amp, "%7.2f");
		xml_set_real(sf_node, REFLEX_SIGAMP, sf->sigamp, "%7.2f");
		xml_set_real(sf_node, REFLEX_PHI, sf->phi, "%7.2f");
		xml_set_real(sf_node, REFLEX_SIGPHI, sf->sigphi, "%7.2f");
		xml_set_real(sf_node, REFLEX_FOM, sf->fom, "%7.4f");
		xml_set_integer(sf_node, REFLEX_STATUS, sf->sel, "%4d");
	}

	return nsf;
}

int			layerline_to_xml(Blayerline* linelist, xmlNodePtr parent)
{
	int				nll(0);
	Blayerline*		line = NULL;
	xmlNodePtr		ll_node = NULL;
	
	for ( line=linelist; line; line=line->next, nll++ ) {
		ll_node = xmlNewChild(parent, NULL, BAD_CAST LAYERLINE, NULL);
		xml_set_integer(ll_node, LAYERLINE_NUMBER, line->number, "%4d");
		xml_set_integer(ll_node, LAYERLINE_ORDER, line->order, "%4d");
		xml_set_real(ll_node, LAYERLINE_DISTANCE, line->distance, "%7.2f");
		xml_set_real(ll_node, LAYERLINE_FREQ, line->freq, "%7.4f");
		xml_set_real(ll_node, LAYERLINE_AMP, line->amp, "%7.2f");
		xml_set_real(ll_node, LAYERLINE_FOM, line->fom, "%7.2f");
		xml_set_integer(ll_node, LAYERLINE_SELECT, line->sel, "%4d");
	}
	
	return nll;
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
int 		project_to_xml(Bproject* project, xmlDocPtr doc, int mg_select, int rec_select)
{
	int					err(0);
	
	Bfield*				field = NULL;
	Bmicrograph*		mg = NULL;
	Breconstruction*	rec = NULL;
	Bparticle*			part;
	Bstring*			data = NULL;

	int					i, j, sel, npart(0);

	xmlNodePtr			node = xmlNewDocPI(doc, BAD_CAST "xml-stylesheet", 
							BAD_CAST "href=\"Bsoft_micrograph.xsl\" type=\"text/xsl\"");
    xmlDocSetRootElement(doc, node);
	
	xmlNodePtr			project_node = xmlNewNode(NULL, BAD_CAST PROJECT);
//	xmlDocSetRootElement(doc, project_node);
	
	xmlAddSibling(node, project_node);
	
	xmlNodePtr			field_node = NULL;
	xmlNodePtr			mg_node = NULL;
	xmlNodePtr			rec_node = NULL;
	
	xmlNewChild(project_node, NULL, BAD_CAST COMMENT, BAD_CAST project->comment.c_str());
	
	for ( data = project->reference; data; data=data->next )
		xmlNewChild(project_node, NULL, BAD_CAST MAP_REFERENCE, BAD_CAST data->c_str());

	for ( i=0, field=project->field; field; field=field->next ) {
		sel = 0;
		if ( project->split ) {
			if ( mg_select < 0 ) sel = 1;
			else for ( j=i, mg=field->mg; mg; mg=mg->next, j++ ) if ( j==mg_select ) sel = 1;
		} else {
			if ( mg_select < 1 ) sel = 1;
			else for ( mg=field->mg; mg; mg=mg->next ) if ( mg->select ) sel = 1;
		}
		if ( sel ) {
			field_node = xmlNewChild(project_node, NULL, BAD_CAST FIELD, NULL);
			xmlNewProp(field_node, BAD_CAST ID, BAD_CAST field->id.c_str());
		}
		for ( mg=field->mg; mg; mg=mg->next, i++ ) {
			sel = 0;
			if ( project->split ) {
				if ( i==mg_select || mg_select < 0 ) sel = 1;
			} else {
				if ( mg_select < 1 || mg->select > 0 ) sel = 1;
			}
			if ( sel ) {
				mg_node = xmlNewChild(field_node, NULL, BAD_CAST MICROGRAPH, NULL);
				xmlNewProp(mg_node, BAD_CAST ID, BAD_CAST mg->id.c_str());
				node = xml_set_integer(mg_node, MICROGRAPH_NUMBER, mg->img_num, "%d");
				if ( mg->fmg.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST MICROGRAPH_FILE, BAD_CAST mg->fmg.c_str());
				if ( mg->fframe.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST MICROGRAPH_FRAMES_FILE, BAD_CAST mg->fframe.c_str());
				if ( mg->fpart.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST PARTICLE_FILE, BAD_CAST mg->fpart.c_str());
				if ( mg->ffil.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST FILAMENT_FILE, BAD_CAST mg->ffil.c_str());
				if ( mg->fft.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST MICROGRAPH_TRANSFORM_FILE, BAD_CAST mg->fft.c_str());
				if ( mg->fps.length() > 0 )
					xmlNewChild(mg_node, NULL, BAD_CAST MICROGRAPH_POWERSPEC_FILE, BAD_CAST mg->fps.c_str());
				xml_set_integer(mg_node, MICROGRAPH_SELECT, mg->select, "%d");
				xml_set_real(mg_node, MICROGRAPH_FOM, mg->fom, "%g");
				xml_set_real(mg_node, MICROGRAPH_MAGNIFICATION, mg->magnification, "%g");
				xml_set_real(mg_node, MICROGRAPH_SAMPLING, mg->sampling, "%g");
				xml_set_real(mg_node, MICROGRAPH_PIXEL_X, mg->pixel_size[0], "%g");
				xml_set_real(mg_node, MICROGRAPH_PIXEL_Y, mg->pixel_size[1], "%g");
				xml_set_real(mg_node, MICROGRAPH_FRAME_PIXEL_X, mg->frame_pixel_size[0], "%g");
				xml_set_real(mg_node, MICROGRAPH_FRAME_PIXEL_Y, mg->frame_pixel_size[1], "%g");
				xml_set_real(mg_node, MICROGRAPH_DOSE, mg->dose, "%g");
				xml_set_real(mg_node, MICROGRAPH_EXPOSURE, mg->exposure, "%g");
				xml_set_real(mg_node, MICROGRAPH_INTENSITY, mg->intensity, "%g");
				xml_set_real(mg_node, MICROGRAPH_WATER_RING, mg->wri, "%g");
				xml_set_real(mg_node, MICROGRAPH_TILT_AXIS, mg->tilt_axis*180.0/M_PI, "%g");
				xml_set_real(mg_node, MICROGRAPH_TILT_ANGLE, mg->tilt_angle*180.0/M_PI, "%g");
				xml_set_real(mg_node, MICROGRAPH_ROT_ANGLE, mg->rot_angle*180.0/M_PI, "%g");
				xml_set_real(mg_node, MICROGRAPH_ORIGIN_X, mg->origin[0], "%g");
				xml_set_real(mg_node, MICROGRAPH_ORIGIN_Y, mg->origin[1], "%g");
				xml_set_real(mg_node, MICROGRAPH_ORIGIN_Z, mg->origin[2], "%g");
				xml_set_real(mg_node, MICROGRAPH_SCALE_X, mg->scale[0], "%g");
				xml_set_real(mg_node, MICROGRAPH_SCALE_Y, mg->scale[1], "%g");
				xml_set_real(mg_node, MICROGRAPH_SCALE_Z, mg->scale[2], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_1_1, mg->matrix[0][0], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_1_2, mg->matrix[0][1], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_1_3, mg->matrix[0][2], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_2_1, mg->matrix[1][0], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_2_2, mg->matrix[1][1], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_2_3, mg->matrix[1][2], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_3_1, mg->matrix[2][0], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_3_2, mg->matrix[2][1], "%g");
				xml_set_real(mg_node, MICROGRAPH_MATRIX_3_3, mg->matrix[2][2], "%g");
				xml_set_real(mg_node, MICROGRAPH_HVEC_X, mg->hvec[0], "%f");
				xml_set_real(mg_node, MICROGRAPH_HVEC_Y, mg->hvec[1], "%f");
				xml_set_real(mg_node, MICROGRAPH_HVEC_Z, mg->hvec[2], "%f");
				xml_set_real(mg_node, MICROGRAPH_KVEC_X, mg->kvec[0], "%f");
				xml_set_real(mg_node, MICROGRAPH_KVEC_Y, mg->kvec[1], "%f");
				xml_set_real(mg_node, MICROGRAPH_KVEC_Z, mg->kvec[2], "%f");
				xml_set_real(mg_node, MICROGRAPH_LVEC_X, mg->lvec[0], "%f");
				xml_set_real(mg_node, MICROGRAPH_LVEC_Y, mg->lvec[1], "%f");
				xml_set_real(mg_node, MICROGRAPH_LVEC_Z, mg->lvec[2], "%f");
				xml_set_real(mg_node, MICROGRAPH_HELIX_AXIS, mg->helix_axis*180.0/M_PI, "%f");
				xml_set_real(mg_node, MICROGRAPH_HELIX_RISE, mg->helix_rise, "%f");
				xml_set_real(mg_node, MICROGRAPH_HELIX_ANGLE, mg->helix_angle*180.0/M_PI, "%f");
				xml_set_real(mg_node, MICROGRAPH_HELIX_RADIUS, mg->helix_radius, "%f");
				xml_set_real(mg_node, PARTICLE_BOX_SIZE_X, mg->box_size[0], "%g");
				xml_set_real(mg_node, PARTICLE_BOX_SIZE_Y, mg->box_size[1], "%g");
				xml_set_real(mg_node, PARTICLE_BOX_SIZE_Z, mg->box_size[2], "%g");
				xml_set_real(mg_node, FILAMENT_WIDTH, mg->filament_width, "%g");
				xml_set_real(mg_node, FILNODE_RADIUS, mg->fil_node_radius, "%g");
				xml_set_real(mg_node, PARTICLE_BAD_RADIUS, mg->bad_radius, "%g");
				xml_set_real(mg_node, REFLEX_RADIUS, mg->sf_radius, "%g");
				xml_set_real(mg_node, MARKER_RADIUS, mg->mark_radius, "%g");
				frame_to_xml(mg->frame, mg_node);
				ctf_to_xml(mg->ctf, mg_node);
				for ( part=mg->part; part; part=part->next, npart++ )
					particle_to_xml(part, mg_node, project->euler_flag, project->omega_flag, project->fom_tag);
				filament_to_xml(mg->fil, mg_node);
				badarea_to_xml(mg->bad, mg_node);
				marker_to_xml(mg->mark, mg_node);
				strucfac_to_xml(mg->sf, mg_node);
				layerline_to_xml(mg->layer, mg_node);
				if ( verbose & VERB_FULL )
					cout << "Writing field \"" << field->id << "\", micrograph \"" 
						<< mg->id << "\" with " << npart << " particles" << endl;
			}
		}
	}
	
	for ( rec = project->rec; rec; rec = rec->next ) {
		sel = 0;
		if ( project->split ) {
			if ( i==rec_select || rec_select < 0 ) sel = 1;
		} else {
			if ( rec_select < 1 || rec->select > 0 ) sel = 1;
		}
		if ( sel ) {
			rec_node = xmlNewChild(project_node, NULL, BAD_CAST MAP, NULL);
			xmlNewProp(rec_node, BAD_CAST ID, BAD_CAST rec->id.c_str());
			if ( rec->frec.length() > 0 )
				xmlNewChild(rec_node, NULL, BAD_CAST MAP_RECONSTRUCTION, BAD_CAST rec->frec.c_str());
			if ( rec->fft.length() > 0 )
				xmlNewChild(rec_node, NULL, BAD_CAST MAP_TRANSFORM_FILE, BAD_CAST rec->fft.c_str());
			if ( rec->fps.length() > 0 )
				xmlNewChild(rec_node, NULL, BAD_CAST MAP_POWERSPEC_FILE, BAD_CAST rec->fps.c_str());
			if ( rec->fpart.length() > 0 && !rec->part->fpart.length() )
				xmlNewChild(rec_node, NULL, BAD_CAST PARTICLE_FILE, BAD_CAST rec->fpart.c_str());
			if ( rec->ffil.length() > 0 )
				xmlNewChild(rec_node, NULL, BAD_CAST FILAMENT_FILE, BAD_CAST rec->ffil.c_str());
			// Note: Should be a list
			if ( rec->model )
				xmlNewChild(rec_node, NULL, BAD_CAST MAP_MODEL, BAD_CAST rec->model->c_str());
			xml_set_integer(rec_node, MAP_SELECT, rec->select, "%d");
			xml_set_real(rec_node, MAP_FOM, rec->fom, "%g");
			xml_set_real(rec_node, MAP_ORIGIN_X, rec->origin[0], "%g");
			xml_set_real(rec_node, MAP_ORIGIN_Y, rec->origin[1], "%g");
			xml_set_real(rec_node, MAP_ORIGIN_Z, rec->origin[2], "%g");
			xml_set_real(rec_node, MAP_SCALE_X, rec->scale[0], "%g");
			xml_set_real(rec_node, MAP_SCALE_Y, rec->scale[1], "%g");
			xml_set_real(rec_node, MAP_SCALE_Z, rec->scale[2], "%g");
			xml_set_real(rec_node, MAP_VOXEL_SIZE_X, rec->voxel_size[0], "%g");
			xml_set_real(rec_node, MAP_VOXEL_SIZE_Y, rec->voxel_size[1], "%g");
			xml_set_real(rec_node, MAP_VOXEL_SIZE_Z, rec->voxel_size[2], "%g");
			if ( rec->symmetry.length() > 0 )
				xmlNewChild(rec_node, NULL, BAD_CAST MAP_SYMMETRY, BAD_CAST rec->symmetry.c_str());
			ctf_to_xml(rec->ctf, rec_node);
			xml_set_real(rec_node, PARTICLE_BOX_SIZE_X, rec->box_size[0], "%g");
			xml_set_real(rec_node, PARTICLE_BOX_SIZE_Y, rec->box_size[1], "%g");
			xml_set_real(rec_node, PARTICLE_BOX_SIZE_Z, rec->box_size[2], "%g");
			xml_set_real(rec_node, PARTICLE_BAD_RADIUS, rec->bad_radius, "%g");
			xml_set_real(rec_node, FILAMENT_WIDTH, rec->filament_width, "%g");
			xml_set_real(rec_node, FILNODE_RADIUS, rec->fil_node_radius, "%g");
			xml_set_real(rec_node, REFLEX_RADIUS, rec->sf_radius, "%g");
			xml_set_real(rec_node, MARKER_RADIUS, rec->mark_radius, "%g");
			xml_set_real(rec_node, MAP_VIEW_X, rec->view[0], "%g");
			xml_set_real(rec_node, MAP_VIEW_Y, rec->view[1], "%g");
			xml_set_real(rec_node, MAP_VIEW_Z, rec->view[2], "%g");
			xml_set_real(rec_node, MAP_VIEW_ANGLE, rec->view.angle()*180.0/M_PI, "%g");
			xml_set_real(rec_node, FILNODE_RADIUS, rec->fil_node_radius, "%g");
			xml_set_real(rec_node, FILNODE_RADIUS, rec->fil_node_radius, "%g");
			for ( part=rec->part; part; part=part->next, npart++ )
				particle_to_xml(part, rec_node, project->euler_flag, project->omega_flag, project->fom_tag);
			filament_to_xml(rec->fil, rec_node);
			badarea_to_xml(rec->bad, rec_node);
			marker_to_xml(rec->mark, rec_node);
			strucfac_to_xml(rec->sf, rec_node);
		}
	}
	
	if ( project->class_avg ) {
		xmlNewChild(project_node, NULL, BAD_CAST CLASS_AVERAGE, NULL);
		for ( part=project->class_avg; part; part=part->next, npart++ )
			particle_to_xml(part, node, project->euler_flag, project->omega_flag, project->fom_tag);
	}
	
	return err;
}

#endif

