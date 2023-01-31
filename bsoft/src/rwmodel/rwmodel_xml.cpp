/**
@file	rwmodel_xml.cpp
@brief	Library routines to read and write XML model parameters
@author 	Bernard Heymann
@date	Created: 20081029
@date	Modified: 20221123
**/

#ifdef HAVE_XML
#include "rwxml.h"
#endif

#include "rwmodel.h"
#include "model_links.h"
#include "model_util.h"
#include "model_tags.h"
#include "file_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads XML model parameters.
@param 	*file_list	list of model parameter file names.
@return Bmodel*		model parameters.
**/
#if defined(HAVE_XML)
Bmodel*		read_model_xml(Bstring* file_list)	
{
	if ( !file_list ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return  NULL;
	}
	
	// Get the list of filenames
	long			i;
	Bstring*		thisfile = NULL;
	
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomptype*		comptype = NULL;
	Bcomponent*		comp = NULL;
	Bcomponent*		comp2 = NULL;
	Blink*			link = NULL;
	Bstring			id, path;
	
    xmlDocPtr		doc;
    xmlNodePtr		root_node;
	xmlNodePtr		comment_node;
    xmlNodePtr		model_node;
    xmlNodePtr		node;
	
	if ( verbose & VERB_DEBUG_XML ) {
		cout << "DEBUG read_xml: XML filenames: " << endl;
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		if ( verbose )
			cout << "Reading file:                   " << *thisfile << endl;
		detect_and_fix_carriage_return(thisfile->c_str());
		doc = xmlParseFile(thisfile->c_str());
		if ( doc == NULL ) return  NULL;

		path = thisfile->pre_rev('/');

		if ( verbose & VERB_PROCESS )
			cout << "# Reading XML file:               " << *thisfile << endl;

		root_node = xmlDocGetRootElement(doc);
		if ( root_node == NULL ) {
			error_show(thisfile->c_str(), __FILE__, __LINE__);
			xmlFreeDoc(doc);
			return  NULL;
		}
		
//		cout << "XML root name: " << root_node->name << endl;
	
		if ( xmlStrcmp(root_node->name, BAD_CAST "model_file") ) {
			cerr <<"Error: The document " << *thisfile << " is not a model file!" << endl;
			xmlFreeDoc(doc);
			return  NULL;
		}
		
//		xmlDocDump(stderr, doc);
//		xmlElemDump(stderr, doc, node);
		
		for ( comment_node = root_node->xmlChildrenNode; comment_node && 
				xmlStrcmp(comment_node->name, BAD_CAST COMMENT); comment_node = comment_node->next ) ;
		
		for ( i=0, model_node = root_node->xmlChildrenNode; model_node; model_node = model_node->next )
				if ( !xmlStrcmp(model_node->name, BAD_CAST MODEL) ) {
			id = 0;
			id = (char *) xmlGetProp(model_node, BAD_CAST ID);
			if ( id.length() < 1 ) id = Bstring(i, "%d");
//			mp = model_add(&mp, id);
//			if ( !model ) model = mp;
			if ( mp ) mp = mp->add(id.str());
			else model = mp = new Bmodel(id.str());
			if ( !xmlStrcmp(model_node->name, BAD_CAST COMMENT) )
				mp->comment((char *) xmlNodeGetContent(model_node));
			else if ( comment_node )
				mp->comment((char *) xmlNodeGetContent(comment_node));
			if ( mp->identifier().size() < 1 ) mp->identifier(to_string(i));
			mp->model_type(xml_get_string(model_node, MODEL_TYPE));
			mp->symmetry(xml_get_string(model_node, MODEL_SYM));
			mp->mapfile(find_file(xml_get_string(model_node, MODEL_MAP_FILENAME), path).str());
			mp->handedness(xml_get_integer(model_node, MODEL_HAND));
			mp->image_number(xml_get_integer(model_node, MODEL_MAP_NUMBER));
			mp->FOM(xml_get_real(model_node, MODEL_FOM));
			mp->select(xml_get_integer(model_node, MODEL_SELECT));
			comptype = NULL;
			comp = NULL;
			for ( node = model_node->xmlChildrenNode; node; node = node->next ) {
				if ( !xmlStrcmp(node->name, BAD_CAST COMPTYPE) ) {
					id = (char *) xmlGetProp(node, BAD_CAST ID);
					comptype = mp->add_type(id);
					if ( !mp->type ) mp->type = comptype;
					comptype->file_name(xml_get_string(node, COMPTYPE_FILENAME));
					comptype->image_number(xml_get_integer(node, COMPTYPE_NUMBER));
					comptype->mass(xml_get_real(node, COMPTYPE_MASS));
					comptype->FOM(xml_get_real(node, COMPTYPE_FOM));
					comptype->select(xml_get_integer(node, COMPTYPE_SELECT));
				}
				if ( !xmlStrcmp(node->name, BAD_CAST COMPONENT) ) {
					id = (char *) xmlGetProp(node, BAD_CAST ID);
//					comp = component_add(&comp, id);
//					if ( !mp->comp ) mp->comp = comp;
//					comp = mp->add_component(id);
					if ( comp ) comp = comp->add(id.str());
					else model->comp = comp = new Bcomponent(id.str());
					id = xml_get_string(node, COMPONENT_TYPE_ID);
//					comp->type = model_add_type_by_id(mp, id);
					comp->type(mp->add_type(id.str()));
					comp->location()[0] = xml_get_real(node, COMPONENT_X);
					comp->location()[1] = xml_get_real(node, COMPONENT_Y);
					comp->location()[2] = xml_get_real(node, COMPONENT_Z);
					comp->view()[0] = xml_get_real(node, COMPONENT_VIEW_X);
					comp->view()[1] = xml_get_real(node, COMPONENT_VIEW_Y);
					comp->view()[2] = xml_get_real(node, COMPONENT_VIEW_Z);
					comp->view()[3] = xml_get_real(node, COMPONENT_VIEW_ANGLE) * M_PI/180.0;
					comp->radius(xml_get_real(node, COMPONENT_RADIUS));
					comp->color()[0] = xml_get_real(node, COMPONENT_RED);
					comp->color()[1] = xml_get_real(node, COMPONENT_GREEN);
					comp->color()[2] = xml_get_real(node, COMPONENT_BLUE);
					comp->color()[3] = xml_get_real(node, COMPONENT_ALPHA);
					comp->density(xml_get_real(node, COMPONENT_DENSITY));
					comp->FOM(xml_get_real(node, COMPONENT_FOM));
					comp->select(xml_get_integer(node, COMPONENT_SELECT));
					comp->description(xml_get_string(node, COMPONENT_DESCRIPTION));
				}
			}
			link = NULL;
			for ( node = model_node->xmlChildrenNode; node; node = node->next ) {
				if ( !xmlStrcmp(node->name, BAD_CAST COMPLINK) ) {
					id = xml_get_string(node, COMPLINK_1);
					for ( comp = mp->comp; comp && comp->identifier() != id.str(); comp = comp->next ) ;
					if ( !comp )
						cerr << "Error: Component " << id << " not found! (1)" << endl;
					id = xml_get_string(node, COMPLINK_2);
					for ( comp2 = mp->comp; comp2 && comp2->identifier() != id.str(); comp2 = comp2->next ) ;
					if ( !comp2 )
						cerr << "Error: Component " << id << " not found! (2)" << endl;
					link = link_add(&link, comp, comp2, 1, 1);
					if ( !mp->link ) mp->link = link;
					link->angle(xml_get_real(node, COMPLINK_ANGLE) * M_PI/180.0);
					link->radius(xml_get_real(node, COMPLINK_RADIUS));
					link->length(xml_get_real(node, COMPLINK_LENGTH));
					link->color()[0] = xml_get_real(node, COMPLINK_RED);
					link->color()[1] = xml_get_real(node, COMPLINK_GREEN);
					link->color()[2] = xml_get_real(node, COMPLINK_BLUE);
					link->color()[3] = xml_get_real(node, COMPLINK_ALPHA);
					link->FOM(xml_get_real(node, COMPLINK_FOM));
					link->select(xml_get_integer(node, COMPLINK_SELECT));
				}
			}
			for ( link=mp->link; link; link=link->next )
				link->length(link->comp[0]->location().distance(link->comp[1]->location()));
		}

		xmlFreeDoc(doc);

		xmlCleanupParser();
	}
	
//	model_list_setup_links(model);
	models_process(model, model_setup_links);
	
	return model;
}
#else
Bmodel*		read_model_xml(Bstring* file_list)	
{
	cerr << "Error: XML files are not supported!" << endl << endl;
	
	return  NULL;
}
#endif
	
/**
@brief 	Writes XML model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			0, <0 on error.
**/
#if defined(HAVE_XML)
int			write_model_xml(Bstring& filename, Bmodel* model)	
{
	int				i, err = 0;

	xmlNodePtr		comment_node;
    xmlNodePtr		model_node;
    xmlNodePtr		type_node;
    xmlNodePtr		comp_node;
    xmlNodePtr		link_node;
	
	Bmodel*			mp = NULL;
	Bcomptype*		comptype = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
	
	if ( verbose & VERB_PROCESS )
		cout << "# Writing file:                 " << filename << endl;
		
	xmlDocPtr		doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
	if ( doc == NULL ) {
		cerr << "Error: The XML document tree was not created!" << endl;
		return  -1;
	}

	xmlNodePtr		root_node = xmlNewDocPI(doc, BAD_CAST "xml-stylesheet", 
						BAD_CAST "href=\"Bsoft_model.xsl\" type=\"text/xsl\"");

	root_node = xmlNewNode(NULL, BAD_CAST "model_file");
    xmlDocSetRootElement(doc, root_node);

	comment_node = xmlNewChild(root_node, NULL, BAD_CAST COMMENT, NULL);
	xmlNodeSetContent(comment_node, BAD_CAST model->comment().c_str());
	
	for ( i=1, mp = model; mp; mp = mp->next, i++ ) {
		model_node = xmlNewChild(root_node, NULL, BAD_CAST MODEL, NULL);
		if ( mp->identifier().length() < 1 ) mp->identifier(Bstring(i, "%d"));
		xmlNewProp(model_node, BAD_CAST ID, BAD_CAST mp->identifier().c_str());
		if ( mp->model_type().length() > 0 )
			xmlNewChild(model_node, NULL, BAD_CAST MODEL_TYPE, BAD_CAST mp->model_type().c_str());
		if ( mp->symmetry().length() > 0 )
			xmlNewChild(model_node, NULL, BAD_CAST MODEL_SYM, BAD_CAST mp->symmetry().c_str());
		xml_set_integer(model_node, MODEL_HAND, mp->handedness(), "%d");
		if ( mp->mapfile().length() > 0 )
			xmlNewChild(model_node, NULL, BAD_CAST MODEL_MAP_FILENAME, BAD_CAST mp->mapfile().c_str());
		xml_set_integer(model_node, MODEL_MAP_NUMBER, mp->image_number(), "%d");
		xml_set_real(model_node, MODEL_FOM, mp->FOM(), "%g");
		xml_set_integer(model_node, MODEL_SELECT, mp->select(), "%d");
		for ( comptype = mp->type; comptype; comptype = comptype->next ) {
			type_node = xmlNewChild(model_node, NULL, BAD_CAST COMPTYPE, NULL);
			xmlNewProp(type_node, BAD_CAST ID, BAD_CAST comptype->identifier().c_str());
			xmlNewChild(type_node, NULL, BAD_CAST COMPTYPE_FILENAME, BAD_CAST comptype->file_name().c_str());
			xml_set_integer(type_node, COMPTYPE_NUMBER, comptype->image_number(), "%4d");
			xml_set_real(type_node, COMPTYPE_MASS, comptype->mass(), "%12.2f");
			xml_set_real(type_node, COMPTYPE_FOM, comptype->FOM(), "%7.4f");
			xml_set_integer(type_node, COMPTYPE_SELECT, comptype->select(), "%4d");
		}
		for ( comp = mp->comp; comp; comp = comp->next ) {
			comp_node = xmlNewChild(model_node, NULL, BAD_CAST COMPONENT, NULL);
			xmlNewProp(comp_node, BAD_CAST ID, BAD_CAST comp->identifier().c_str());
			xmlNewChild(comp_node, NULL, BAD_CAST COMPONENT_TYPE_ID, BAD_CAST comp->type()->identifier().c_str());
			xml_set_real(comp_node, COMPONENT_X, comp->location()[0], "%8.3f");
			xml_set_real(comp_node, COMPONENT_Y, comp->location()[1], "%8.3f");
			xml_set_real(comp_node, COMPONENT_Z, comp->location()[2], "%8.3f");
			xml_set_real(comp_node, COMPONENT_VIEW_X, comp->view()[0], "%7.4lf");
			xml_set_real(comp_node, COMPONENT_VIEW_Y, comp->view()[1], "%7.4lf");
			xml_set_real(comp_node, COMPONENT_VIEW_Z, comp->view()[2], "%7.4lf");
			xml_set_real(comp_node, COMPONENT_VIEW_ANGLE, comp->view()[3] * 180.0/M_PI, "%7.2lf");
			xml_set_real(comp_node, COMPONENT_RADIUS, comp->radius(), "%8.3f");
			xml_set_real(comp_node, COMPONENT_RED, comp->color()[0], "%5.3f");
			xml_set_real(comp_node, COMPONENT_GREEN, comp->color()[1], "%5.3f");
			xml_set_real(comp_node, COMPONENT_BLUE, comp->color()[2], "%5.3f");
			xml_set_real(comp_node, COMPONENT_ALPHA, comp->color()[3], "%5.3f");
			xml_set_real(comp_node, COMPONENT_DENSITY, comp->density(), "%8.3f");
			xml_set_real(comp_node, COMPONENT_FOM, comp->FOM(), "%8.3f");
			xml_set_integer(comp_node, COMPONENT_SELECT, comp->select(), "%4d");
			xmlNewChild(type_node, NULL, BAD_CAST COMPONENT_DESCRIPTION, BAD_CAST comp->description().c_str());
		}
		for ( link = mp->link; link; link = link->next ) {
			link_node = xmlNewChild(model_node, NULL, BAD_CAST COMPLINK, NULL);
			xmlNewChild(link_node, NULL, BAD_CAST COMPLINK_1, BAD_CAST link->comp[0]->identifier().c_str());
			xmlNewChild(link_node, NULL, BAD_CAST COMPLINK_2, BAD_CAST link->comp[1]->identifier().c_str());
			xml_set_real(link_node, COMPLINK_ANGLE, link->angle() * 180.0/M_PI, "%7.2f");
			xml_set_real(link_node, COMPLINK_RADIUS, link->radius(), "%8.3f");
			xml_set_real(link_node, COMPLINK_LENGTH, link->length(), "%8.3f");
			xml_set_real(link_node, COMPLINK_RED, link->color()[0], "%5.3f");
			xml_set_real(link_node, COMPLINK_GREEN, link->color()[1], "%5.3f");
			xml_set_real(link_node, COMPLINK_BLUE, link->color()[2], "%5.3f");
			xml_set_real(link_node, COMPLINK_ALPHA, link->color()[3], "%5.3f");
			xml_set_real(link_node, COMPLINK_FOM, link->FOM(), "%8.3f");
			xml_set_integer(link_node, COMPLINK_SELECT, link->select(), "%4d");
		}
	}
	
	if ( err == 0 ) err = xmlSaveFormatFile(filename.c_str(), doc, 1);
	
	xmlFreeDoc(doc);

	return err;
}
#else
int			write_model_xml(Bstring& filename, Bmodel* model)	
{
	cerr << "Error: XML files are not supported!" << endl << endl;
	
	return 0;
}
#endif

