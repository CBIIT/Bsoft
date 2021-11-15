/**
@file	rwFSC_XML.cpp
@brief	Library routines to read and write FSC curves
@author Bernard Heymann
@date	Created: 20121216
@date	Modified: 20140423
**/

#ifdef HAVE_XML
#include "rwxml.h"
#endif

#include "ps_plot.h" 
#include "utilities.h"

#define FSCX	"Resolution (A-1)"
#define	FSCY	"Correlation Coefficient"

// Declaration of the global variables 
extern int 	verbose;		// Output verbosity

#ifdef HAVE_XML
/**
@brief 	Reads an XML file with an FSC curve.
@param 	&filename	input file name.
@param	&xsdfile	XML schema to validate against (can be empty).
@return Bplot* 		new plot.

	A new plot with two columns is generated.

**/
Bplot*		xml_read_fsc(Bstring &filename, Bstring &xsdfile)
{
	if ( verbose & VERB_PROCESS )
		cout << "# Reading XML file:               " << filename << endl;

    xmlDocPtr			doc = xmlParseFile(filename.c_str());

	if ( xsdfile.length() ) xml_validate(doc, xsdfile);

    xmlNodePtr			fsc_node = xmlDocGetRootElement(doc);
	xmlNodePtr			coor_node;

	if ( fsc_node == NULL ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		xmlFreeDoc(doc);
		return NULL;
	}

	if ( xmlStrcmp(fsc_node->name, (const xmlChar *) "fsc") ) {
		cerr << "Error: The document " << filename << " is not an FSC curve!" << endl;
		xmlFreeDoc(doc);
		return NULL;
	}

	int					i, ncol(2), nrow(0);
	
	for ( nrow=0, coor_node = fsc_node->xmlChildrenNode; coor_node; coor_node = coor_node->next )
		if ( !xmlStrcmp(coor_node->name, BAD_CAST "coordinate") ) nrow++;
	
	Bstring				title = (char *) xmlGetProp(fsc_node, BAD_CAST "title");
	Bplot*				plot = new Bplot(1, nrow, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label((char *) xmlGetProp(fsc_node, BAD_CAST "xaxis"));
	plot->page(0).column(1).label((char *) xmlGetProp(fsc_node, BAD_CAST "yaxis"));
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(1/hi_res);
//	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(1).label("Resolution(A)");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).inc(0.1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);

	for ( i=0, coor_node = fsc_node->xmlChildrenNode; coor_node; coor_node = coor_node->next ) {
		if ( !xmlStrcmp(coor_node->name, BAD_CAST "coordinate") ) {
			(*plot)[i] = xml_get_real(coor_node, "x");
			(*plot)[nrow+i] = xml_get_real(coor_node, "y");
			i++;
		}
	}

	xmlFreeDoc(doc);

	xmlCleanupParser();

	return plot;
}

/**
@brief 	Reads an XML file with an FSC curve.
@param 	&filename	input file name.
@return Bplot*		new plot.

	A new plot with two columns is generated.

**/
Bplot*		xml_read_fsc(Bstring &filename)
{
	Bstring		xsdfile;
	return xml_read_fsc(filename, xsdfile);
}

/**
@brief 	Writes one or more XML files from a FSC plot.
@param 	&filename	output file name.
@param 	*plot		resolution plot information.
@return void

	The plot must have 2 columns for each page as follows:
	1.	Spatial frequency or Resolution (A-1)
	2.	FSC or Correlation Coefficient
	A plot with multiple pages gnerates multiple, numbered output files.

**/
void		xml_write_fsc(Bstring& filename, Bplot* plot)
{
	long		n, i, j, k, cx, cy, nrows = plot->rows();
	Bstring		filename1(filename);

    xmlDocPtr			doc;
	xmlNodePtr			node;
	xmlNodePtr			fsc_node;
	xmlNodePtr			coor_node;

	for ( n=0; n<plot->pages(); n++ ) {
		if ( plot->pages() > 1 ) filename1 = filename.pre_rev('.') + Bstring(n+1, "_%03ld.xml");
		if ( verbose & VERB_PROCESS )
			cout << "# Writing an XML file:            " << filename1 << endl << endl;
		cx = plot->page(n).column(0).number();
		cy = plot->page(n).column(1).number();
		doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
		node = xmlNewDocPI(doc, BAD_CAST "xml-stylesheet", 
							BAD_CAST "href=\"Bsoft_micrograph.xsl\" type=\"text/xsl\"");
		xmlDocSetRootElement(doc, node);
		fsc_node = xmlNewNode(NULL, BAD_CAST "fsc");
		xmlAddSibling(node, fsc_node);
		xmlNewProp(fsc_node, BAD_CAST "title", BAD_CAST "FSC Plot");
		xmlNewProp(fsc_node, BAD_CAST "xaxis", BAD_CAST FSCX);
		xmlNewProp(fsc_node, BAD_CAST "yaxis", BAD_CAST FSCY);
		for ( i=cx*nrows, j=cy*nrows, k=0; k<nrows; i++, j++, k++ ) {
			coor_node = xmlNewChild(fsc_node, NULL, BAD_CAST "coordinate", NULL);
			node = xml_set_real(coor_node, "x", (*plot)[i], "%g");
			node = xml_set_real(coor_node, "y", (*plot)[j], "%g");
		}
		xmlSaveFormatFile(filename1.c_str(), doc, 1);
		xmlFreeDoc(doc);
	}
}

#else

/**
@brief 	Reads an XML file with an FSC curve.
@param 	&filename	input file name.
@return Bplot*		new plot.

	A new plot with two columns is generated.

**/
Bplot*		xml_read_fsc(Bstring& filename)
{
	if ( verbose & VERB_PROCESS )
		cout << "# Reading XML file:               " << filename << endl;

    ifstream			fx;
    fx.open(filename.c_str());
    if ( fx.fail() ) return NULL;
	
	int					i, j, ncol(2), nrow(0);
	float				x, y;
    char				aline[MAXLINELEN];

 	while ( !fx.eof() ) {
		fx.getline(aline, MAXLINELEN);
		if ( strstr(aline, "<coordinate>") ) nrow++;
	}
	
	Bstring				title = "FSC Plot";
	Bplot*				plot = new Bplot(1, nrow, ncol);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(ncol);
	for ( i=0; i<ncol; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label(FSCX);
	plot->page(0).column(1).label(FSCY);
//	plot->page(0).axis(1).min(0);
//	plot->page(0).axis(1).max(1/hi_res);
//	plot->page(0).axis(1).inc(0.1/hi_res);
	plot->page(0).axis(1).flags(1);
	plot->page(0).axis(1).label("Resolution(A)");
	plot->page(0).axis(3).min(0);
	plot->page(0).axis(3).max(1);
	plot->page(0).axis(3).inc(0.1);
	plot->page(0).column(1).type(2);
	plot->page(0).column(0).axis(1);
	plot->page(0).column(1).axis(3);

	fx.seekg(0, ios::beg);

 	for ( i=0, j=nrow; !fx.eof();  ) {
		fx.getline(aline, MAXLINELEN);
		if ( sscanf(aline, "<x>%f</x>", &x) ) (*plot)[i] = x;
		if ( sscanf(aline, "<y>%f</y>", &y) ) { (*plot)[j] = y; i++; j++; }
	}
	
	fx.close();

	return plot;
}

void		xml_write_fsc(Bstring& filename, Bplot* plot)
{
	long		n, i, j, k, cx, cy, nrows = plot->rows();
	Bstring		filename1(filename);
	ofstream	fx;

	for ( n=0; n<plot->pages(); n++ ) {
		if ( plot->pages() > 1 ) filename1 = filename.pre_rev('.') + Bstring(n+1, "_%03ld.xml");
		cx = plot->page(n).column(0).number();
		cy = plot->page(n).column(1).number();
		fx.open(filename1.c_str());
		fx << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
		fx << "<fsc title=\"FSC Plot\" xaxis=\"" << FSCX << "\" yaxis=\"" << FSCY << "\">" << endl;
		for ( i=cx*nrows, j=cy*nrows, k=0; k<nrows; i++, j++, k++ ) {
			fx << "<coordinate>" << endl;
			fx << "\t<x>" << (*plot)[i] << "</x>" << endl;
			fx << "\t<y>" << (*plot)[j] << "</y>" << endl;
			fx << "</coordinate>" << endl;
		}
		fx << "</fsc>" << endl;
		fx.close();
	}
}

#endif


