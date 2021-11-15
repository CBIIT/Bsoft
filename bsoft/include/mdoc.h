/**
@file	mdoc.h
@author	Bernard Heymann
@date	20190109
**/

#include <iomanip>
#include <sstream>

#include "json.h"
#include "utilities.h"

using namespace std;

#ifndef _MDOC_


class MDOCparser {
private:
	string		s, f;
	void		get_clean_line(ifstream& fmdoc) {
		if ( fmdoc.eof() ) return;
		getline(fmdoc, s);
		s.erase(remove(s.end()-1, s.end(), '\r'), s.end());
//		cout << s << endl;
	}
	void		read_object(JSvalue& obj) {
		if ( s.length() < 3 ) return;
		size_t			i(s.find_first_of("="));
		string			tag(s.substr(0,i-1));
		i = s.find_first_not_of(" \t", i+1);
		obj[tag] = parse_string(s.substr(i));
	}
	JSvalue		read_section(ifstream& fmdoc) {
		JSvalue			obj(JSobject);
		s = s.substr(1,s.length()-2);
//		cout << "Reading section: " << s << endl;
		while ( !fmdoc.eof() && s[0] != '[' ) {
			read_object(obj);
			get_clean_line(fmdoc);
		}
		return obj;
	}
	JSvalue		parse_string(string s) {
		// Remove any quotes
		s.erase( remove( s.begin(), s.end(), '\"' ), s.end() );
    	istringstream	ss(s);
    	vector<string>	sv((istream_iterator<string>(ss)), istream_iterator<string>());
		int				ft(0);
		for ( auto it = s.begin(); it!=s.end(); ++it ) {
			if ( isalpha(*it) ) return JSvalue(s);
			else if ( *it == '.' ) ft = 1;
		}
		if ( sv.size() > 1 ) {
//			cout << "#" << s << "#" << tab << sv.size() << endl;
			if ( ft ) {
				vector<double>	d(sv.size());
				try { for ( size_t i=0; i<sv.size(); ++i ) d[i] = stod(sv[i]); }
				catch (...) { fail("No conversion to a real number."); }
				return JSvalue(d);
			} else {
				vector<long>	l(sv.size());
				try { for ( size_t i=0; i<sv.size(); ++i ) l[i] = stod(sv[i]); }
				catch (...) { fail("No conversion to an integer."); }
				return JSvalue(l);
			}
		} else {
			if ( ft ) {
				double			d;
				try { d = stod(s); }
				catch (...) { fail("No conversion to a real number."); }
				return JSvalue(d);
			} else {
				long			l;
				try { l = stol(s); }
				catch (...) { fail("No conversion to an integer."); }
				return JSvalue(l);
			}
		}
	}
	void		fail(string msg) {
		cerr << "MDOC parser error: " << msg << endl;
		exit(-1);
	}
public:
	MDOCparser() {};
	MDOCparser(string filename) { f = filename; }
	JSvalue		parse() {
		if ( f.length() < 1 ) fail("No filename specified!");
		return parse(f);
	}
	JSvalue		parse(string filename) {
		f = filename;
//		cout << "Reading " << f << endl;
		ifstream	fmdoc(f);
		if ( fmdoc.fail() ) fail("File " + f + " not opened");
	
		JSvalue			root(JSobject);
		JSvalue			sections(JSarray);
		
		get_clean_line(fmdoc);

		while ( !fmdoc.eof() ) {
			if ( s[0] == '[' ) sections.push_back(read_section(fmdoc));
			else {
				read_object(root);
				get_clean_line(fmdoc);
			}
//			cout << "***" << s << endl;
		}
		
		root["Sections"] = sections;

//		cout << "Number of sections: " << sections.size() << endl;
		
		fmdoc.close();
	
//		cout << "done " << f << endl;
	
		return root;
	}
};

#define _MDOC_
#endif
