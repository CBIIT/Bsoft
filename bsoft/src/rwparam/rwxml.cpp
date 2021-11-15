/**
@file	rwxml.cpp
@brief	Reads to and writes from XML structures
@author Bernard Heymann
@date	Created: 20050920
@date	Modified: 20210615
**/

#ifdef HAVE_XML

#include "rwxml.h"
#include "string_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Checks for an XML node.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return int 		1 if found, 0 if not.

	The children of the parent XML parent node is traversed to find the 
	node with the given XML tag.

**/
int			xml_check_for_node(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node;
	
	for ( node = parent->xmlChildrenNode; node; node=node->next )
        if ( !xmlStrcmp(node->name, BAD_CAST tag) )
			return 1;
	
	return 0;
}

/**
@brief 	Checks recursively for an XML node in a tree.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return int 		1 if found, 0 if not.

	The tree below the parent XML node is traversed to find a 
	node with the given XML tag.

**/
int			xml_check_for_node_in_tree(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node;
	
	int				result(0);
	
	for ( node = parent->xmlChildrenNode; node && result==0; node=node->next ) {
        if ( !xmlStrcmp(node->name, BAD_CAST tag) )
			result = 1;
		else
			result = xml_check_for_node_in_tree(node, tag);
	}
	
	return result;
}

/**
@brief 	Finds a child node with a given tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return xmlNodePtr 	node pointer.

	The children of the parent XML parent node is traversed to find the 
	node with the given XML tag.

**/
xmlNodePtr	xml_find_node(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node = NULL;
	
	for ( node = parent->xmlChildrenNode; node; node = node->next )
        if ( !xmlStrcmp(node->name, BAD_CAST tag) )
			return node;
	
	if ( verbose & VERB_FULL )
		cerr << "Node " << tag << " of parent " << parent->name << " not found!" << endl;
	
	return node;
}

/**
@brief 	Copies a string value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return Bstring 		content string, empty on error.

	The children of the parent XML node is traversed to obtain the 
	string associated with the given XML tag.

**/
Bstring		xml_copy_string(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node = xml_find_node(parent, tag);
	
	Bstring 		s;
	
	if ( !node ) return s;
	
	xmlChar*		content = xmlNodeGetContent(node);

	if ( !content ) return s;
	
	s = (char *)content;
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_copy_string: tag=" << tag << " string=" << s << " (" << s.length() << ")" << endl;
	
	return s;
}

/**
@brief 	Gets a string attribute value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML attribute tag.
@return string 		the string value, empty if the tag doesn't exist.

	The children of the parent XML node is traversed to obtain the 
	integer value associated with the given XML tag.

**/
string		xml_get_string_attribute(xmlNodePtr parent, const char* tag)
{
	xmlChar*		content = xmlGetProp(parent, BAD_CAST tag);
	
	string			value("");
	
	if ( !content ) return value;
	
	value = (char *)content;
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_string_attribute: tag=" << tag << " value=" << value << endl;
	
	return value;
}

/**
@brief 	Gets a integer attribute value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML attribute tag.
@return long 		the integer value, 0 if the tag doesn't exist.

	The children of the parent XML node is traversed to obtain the 
	integer value associated with the given XML tag.

**/
long		xml_get_integer_attribute(xmlNodePtr parent, const char* tag)
{
	xmlChar*		content = xmlGetProp(parent, BAD_CAST tag);
	
	if ( !content ) return 0;
	
	long			value = atol((char *)content);
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_integer_attribute: tag=" << tag << " value=" << value << endl;
	
	return value;
}

/**
@brief 	Gets a real attribute value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML attribute tag.
@return double 		the real value, 0 if the tag doesn't exist.

	The children of the parent XML node is traversed to obtain the 
	integer value associated with the given XML tag.

**/
double		xml_get_real_attribute(xmlNodePtr parent, const char* tag)
{
	xmlChar*		content = xmlGetProp(parent, BAD_CAST tag);
	
	if ( !content ) return 0;
	
	double			value = atof((char *)content);
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_real_attribute: tag=" << tag << " value=" << value << endl;
	
	return value;
}

/**
@brief 	Returns a string value associated with an XML node.
@param 	node		XML node.
@return string 		the string.

**/
string		xml_get_string(xmlNodePtr node)
{
	string			s;

	if ( !node ) return s;

	xmlChar*		content = xmlNodeGetContent(node);

	if ( !content ) return s;
	
	s = (char *) content;
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_string: s=" << s << " (" << s.length() << ")" << endl;
	
	return s;
}

/**
@brief 	Returns a string value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return string 		the string.

	The children of the parent XML node is traversed to obtain the 
	string associated with the given XML tag.

**/
string		xml_get_string(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node = xml_find_node(parent, tag);
	
	return			xml_get_string(node);
}

/**
@brief 	Gets a integer value associated with an XML node.
@param 	node		XML node.
@return long 			the integer value, 0 if the tag doesn't exist.

**/
long		xml_get_integer(xmlNodePtr node)
{
	long			value(0);
	
	if ( !node ) return value;
	
	xmlChar*		content = xmlNodeGetContent(node);

	if ( !content ) return value;
	
	value = atol((char *)content);
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_integer: value=" << value << endl;
	
	return value;
}

/**
@brief 	Gets a integer value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@return long 			the integer value, 0 if the tag doesn't exist.

	The children of the parent XML node is traversed to obtain the
	integer value associated with the given XML tag.

**/
long		xml_get_integer(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node = xml_find_node(parent, tag);
	
	return xml_get_integer(node);
}

/**
@brief 	Gets a floating point value associated with an XML node.
@param 	node		XML node.
@return double 		the floating point value, 0 if the tag doesn't exist.

**/
double		xml_get_real(xmlNodePtr node)
{
	double			value(0);
	
	if ( !node ) return value;
	
	xmlChar*		content = xmlNodeGetContent(node);

	if ( !content ) return value;
	
	value = atof((char *)content);
	
	xmlFree(content);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_get_real: value=" << value << endl;
	
	return value;
}

/**
@brief 	Gets a floating point value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag string.
@return double 		the floating point value, 0 if the tag doesn't exist.

	The children of the parent XML parent node is traversed to obtain the
	floating point value associated with the given XML tag.

**/
double		xml_get_real(xmlNodePtr parent, const char* tag)
{
	xmlNodePtr		node = xml_find_node(parent, tag);
	
	return xml_get_real(node);
}

/**
@brief 	Sets an integer value associated with an XML attribute tag.
@param 	parent		parent XML node.
@param	*tag		an XML attribute tag.
@param 	value		tag value.
@param	*format		string content format.
@return xmlAttrPtr 	new attribute, ? if the tag doesn't exist.

	A new child node is added with the given XML tag and its content 
	set as a string with the given format. 

**/
xmlAttrPtr	xml_set_integer_attribute(xmlNodePtr parent, const char* tag, long value, const char* format)
{
	xmlChar			s[256];
	
	xmlStrPrintf(s, 256, (xmlFormat) format, value);
	
	xmlAttrPtr		attr = xmlNewProp(parent, BAD_CAST tag, BAD_CAST s);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_set_integer_attribute: tag=" << tag << " string=" << s << " (" << xmlStrlen(s) << ")" << endl;
	
	return attr;
}

/**
@brief 	Sets an real value associated with an XML attribute tag.
@param 	parent		parent XML node.
@param	*tag		an XML attribute tag.
@param 	value		tag value.
@param	*format		string content format.
@return xmlAttrPtr 	new attribute, ? if the tag doesn't exist.

	A new child node is added with the given XML tag and its content 
	set as a string with the given format. 

**/
xmlAttrPtr	xml_set_real_attribute(xmlNodePtr parent, const char* tag, double value, const char* format)
{
	xmlChar			s[256];
	
	xmlStrPrintf(s, 256, (xmlFormat) format, value);
	
	xmlAttrPtr		attr = xmlNewProp(parent, BAD_CAST tag, BAD_CAST s);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_set_real_attribute: tag=" << tag << " string=" << s << " (" << xmlStrlen(s) << ")" << endl;
	
	return attr;
}

/**
@brief 	Sets an integer value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@param 	value		tag value.
@param	*format		string content format.
@return xmlNodePtr 	new child node, ? if the tag doesn't exist.

	A new child node is added with the given XML tag and its content 
	set as a string with the given format. 

**/
xmlNodePtr	xml_set_integer(xmlNodePtr parent, const char* tag, long value, const char* format)
{
	xmlChar			s[256];
	
	xmlStrPrintf(s, 256, (xmlFormat) format, value);
	
	xmlNodePtr		node = xmlNewChild(parent, NULL, BAD_CAST tag, BAD_CAST s);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_set_integer: tag=" << tag << " string=" << s << " (" << xmlStrlen(s) << ")" << endl;
	
	return node;
}

/**
@brief 	Sets a floating point value associated with an XML tag.
@param 	parent		parent XML node.
@param	*tag		an XML tag.
@param 	value		tag value.
@param	*format		string content format.
@return xmlNodePtr 	new child node, ? if the tag doesn't exist.

	A new child node is added with the given XML tag and its content 
	set as a string with the given format. 

**/
xmlNodePtr	xml_set_real(xmlNodePtr parent, const char* tag, double value, const char* format)
{
	xmlChar			s[256];
	
	xmlStrPrintf(s, 256, (xmlFormat) format, value);
	
	xmlNodePtr		node = xmlNewChild(parent, NULL, BAD_CAST tag, BAD_CAST s);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_set_real: tag=" << tag << " string=" << s << " (" << xmlStrlen(s) << ")" << endl;
	
	return node;
}

/**
@brief 	Validates an XML document based on a schema.
@param 	doc			document pointer.
@param 	&xsdfile	schema file name.
@return int 		1 if valid, 0 if not, <0 on error.

	A new child node is added with the given XML tag and its content 
	set as a string with the given format. 

**/
int 		xml_validate(const xmlDocPtr doc, Bstring& xsdfile)
{
	if ( xsdfile.length() < 1 ) return 0;
	
	xmlDocPtr schema_doc = xmlParseFile(xsdfile.c_str());
	
    if ( !schema_doc ) {
        cerr << "Error: the schema cannot be loaded or is not well-formed" << endl;
        return -1;
    }
	
    xmlSchemaParserCtxtPtr parser_ctxt = xmlSchemaNewDocParserCtxt(schema_doc);
    if ( !parser_ctxt ) {
        cerr << "Error: unable to create a parser context for the schema" << endl;
        xmlFreeDoc(schema_doc);
        return -2;
    }
	
    xmlSchemaPtr schema = xmlSchemaParse(parser_ctxt);
    if ( !schema ) {
        cerr << "Error: the schema is not valid" << endl;
        xmlSchemaFreeParserCtxt(parser_ctxt);
        xmlFreeDoc(schema_doc);
        return -3;
    }
	
    xmlSchemaValidCtxtPtr valid_ctxt = xmlSchemaNewValidCtxt(schema);
    if ( !valid_ctxt ) {
        cerr << "Error: unable to create a validation context for the schema" << endl;
        xmlSchemaFree(schema);
        xmlSchemaFreeParserCtxt(parser_ctxt);
        xmlFreeDoc(schema_doc);
        return -4; 
    }
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG xml_validate: file=" << xsdfile << endl;
	
    int 	is_valid = (xmlSchemaValidateDoc(valid_ctxt, doc) == 0);
	
	if ( verbose ) {
		cout << "File is ";
		if ( is_valid ) cout << "valid";
		else cout << "not valid";
		cout << " (" << xsdfile << ")" << endl;
	}
	
	xmlSchemaFreeValidCtxt(valid_ctxt);
	xmlSchemaFree(schema);
	xmlSchemaFreeParserCtxt(parser_ctxt);
	xmlFreeDoc(schema_doc);
	
    /* force the return value to be non-negative on success */
    return is_valid ? 1 : 0;
}

JSvalue		json_from_xml(const char* data)
{
	JSvalue			js(JSobject);
	string			tag, val;

	xmlChar*		chr = xmlStrdup((xmlChar*)data);
	xmlDocPtr		doc = xmlParseDoc(chr);
	xmlNodePtr		node = xmlDocGetRootElement(doc);
	
	for ( node = node->xmlChildrenNode; node; node = node->next ) {
		if ( verbose & VERB_FULL )
			cout << xmlGetProp(node, BAD_CAST "name") << tab << xmlNodeGetContent(node) << endl;
//		tag = (char *)xmlGetProp(node, BAD_CAST "name");
//		js[tag] = atof((char *)xmlNodeGetContent(node));
		tag = xml_get_string_attribute(node, "name");
		val = xml_get_string(node);
		if ( check_for_number(val) )
			js[tag] = xml_get_real(node);
		else
			js[tag] = xml_get_integer(node);
	}

	return js;
}

#endif

