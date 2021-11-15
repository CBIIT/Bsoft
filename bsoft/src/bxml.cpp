/**
@file	bxml.cpp
@brief	Reads, writes and validates XML files
@author Bernard Heymann
@date	Created: 20140423
@date	Modified: 20140424
**/

#include "rwxml.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistence
const char* use[] = {
" ",
"Usage: bxml [options] input.xml [input.xml]",
"-------------------------------------------",
"Reads, writes and validates XML files.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Input:",
"-validate mg.xsd         Validate XML input file(s) with this schema.",
" ",
"Output:",
"-output file.xml         Output XML format file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	Bstring			xsdfile;					// XML schema file
	Bstring			outfile;					// Output XML format file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "validate" )
			xsdfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No XML files specified!" << endl;
		bexit(-1);
	}

	Bstring*		thisfile;
	xmlDocPtr		doc;
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		doc = xmlParseFile(thisfile->c_str());
		if ( !doc ) {
			cerr << "Error: " << *thisfile << " not parsed!" << endl;
			bexit(-1);
		} else {
			if ( verbose & VERB_PROCESS )
				cout << "# Reading file:                   " << *thisfile << endl;
			if ( xsdfile.length() ) xml_validate(doc, xsdfile);
		}
		if ( outfile.length() ) xmlSaveFormatFile(outfile.c_str(), doc, 1);
		xmlFreeDoc(doc);
	}

	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

