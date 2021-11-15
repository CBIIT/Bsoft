/**
@file	json.cpp
@author	Bernard Heymann
@date	20160409 - 20180328

**/

#include "json.h"

ostream&	operator<<(ostream& output, JStype t) {
	switch ( t ) {
		case JSstring: cout << "JSstring"; break;
		case JSbool: cout << "JSbool"; break;
		case JSinteger: cout << "JSinteger"; break;
		case JSreal: cout << "JSreal"; break;
		case JSobject: cout << "JSobject"; break;
		case JSarray: cout << "JSarray"; break;
		case JSnull: cout << "JSnull"; break;
		default: cout << "undefined";
	}
	return output;
}

int		indent(0);

void	indent_line(ostream& output) {
	for ( int i=0; i<indent; ++i ) output << "\t";
}

ostream&	operator<<(ostream& output, JSvalue& jsv) {
//	cout << "***** " << jsv.type() << " *****" << endl;
	if ( jsv.type() == JSstring ) {
		output << "\"" << jsv.value() << "\"";
	} else if ( jsv.type() == JSbool ) {
		output << ((jsv.boolean())? "true": "false");
	} else if ( jsv.type() == JSinteger ) {
		output << jsv.integer();
	} else if ( jsv.type() == JSreal ) {
		output << jsv.real();
	} else if ( jsv.type() == JSobject ) {
		output << "{";
		if ( jsv.size() ) output << endl;
		indent++;
		for ( auto it=jsv.object_begin(); it!=jsv.object_end(); ++it ) {
			if ( it != jsv.object_begin() ) output << ", " << endl;
			indent_line(output);
			output << "\"" << it->first << "\": " << it->second;
		}
		indent--;
		if ( jsv.size() ) {
			output << endl;
			indent_line(output);
		}
		output << "}";
	} else if ( jsv.type() == JSarray ) {
		output << "[";
		if ( jsv.size() > 1 ) {
			output << endl;
			indent++;
			indent_line(output);
		}
		for ( auto it=jsv.begin(); it!=jsv.end(); ++it ) {
			if ( it != jsv.begin() ) output << ", ";
			output << *it;
		}
		if ( jsv.size() > 1 ) {
			output << endl;
			indent--;
			indent_line(output);
		}
		output << "]";
	} else if ( jsv.type() == JSnull ) {
		output << "null";
	}
	return output;
}





