/**
@file	rwxml.h
@brief	Reads and writes micrograph XML files
@author Bernard Heymann
@date	Created: 20050920
@date	Modified: 20210615
**/

#include <libxml/parser.h>
#include <libxml/xmlschemas.h>
#include "Bstring.h"
#include "json.h"

#if defined(OLD_XML)
typedef const xmlChar* xmlFormat;
#else
typedef const char* xmlFormat;
#endif

// Function prototypes
int			xml_check_for_node(xmlNodePtr parent, const char* tag);
int			xml_check_for_node_in_tree(xmlNodePtr parent, const char* tag);
xmlNodePtr	xml_find_node(xmlNodePtr parent, const char* tag);
Bstring		xml_copy_string(xmlNodePtr parent, const char* tag);
string		xml_get_string_attribute(xmlNodePtr parent, const char* tag);
long		xml_get_integer_attribute(xmlNodePtr parent, const char* tag);
double		xml_get_real_attribute(xmlNodePtr parent, const char* tag);
string		xml_get_string(xmlNodePtr node);
string		xml_get_string(xmlNodePtr parent, const char* tag);
long		xml_get_integer(xmlNodePtr node);
long		xml_get_integer(xmlNodePtr parent, const char* tag);
double		xml_get_real(xmlNodePtr node);
double		xml_get_real(xmlNodePtr parent, const char* tag);
xmlAttrPtr	xml_set_integer_attribute(xmlNodePtr parent, const char* tag, long value, const char* format);
xmlAttrPtr	xml_set_real_attribute(xmlNodePtr parent, const char* tag, double value, const char* format);
xmlNodePtr	xml_set_integer(xmlNodePtr parent, const char* tag, long value, const char* format);
xmlNodePtr	xml_set_real(xmlNodePtr parent, const char* tag, double value, const char* format);
int 		xml_validate(const xmlDocPtr doc, Bstring& xsdfile);
JSvalue		json_from_xml(const char* data);

