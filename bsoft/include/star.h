/**
@file	star.h
@author	Bernard Heymann
@date	20151106 - 20220103
**/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <regex>
#include <algorithm>

#include "string_util.h"
#include "utilities.h"

using namespace std;

#ifndef _STAR_

/**
@class 	BstarLoop
@brief	Structure for a table of entries in a STAR database block.

	Each item consists of a string or a zero-delimited sequence of strings.
	The tag map values indicate the order of the row elements and must be set up properly.
***/
class BstarLoop {
private:
	map<string, int>			tg;	// Vector of item identifiers
	vector< vector<string> >	d;	// Table of data strings
public:
	BstarLoop() { }
	BstarLoop(ifstream& fstar) {
		string			s;
		int				ntag(0);
	
		while ( !fstar.eof() ) {
			s = get_clean_line(fstar);
//			cout << "=" << s << "=" << endl;
			if ( !s.length() ) break;
			if ( s[0] == '#' ) break;
			if ( s[0] == '_' ) {		// Strip off the underscore
				s = s.substr(1, s.find_first_of(' ')-1);
				tg[s] = ntag;
				ntag++;
			} else {
				vector<string>	v = splitn(s, ntag);
				if ( v[ntag-1].size() < 1 ) break;
				d.push_back(v);
			}
		}
	}
	int				write(ofstream& fstar) {
		int				err(0);
		
		vector<string>	vs(tg.size());

		for ( auto s: tg )
			vs[s.second] = s.first;		// Get the correct order of the tags
/*
		fstar << endl << "loop_" << endl;
	
		for ( auto s: vs )
			fstar << "_" << s << endl;	// Lead with an underscore

		for ( auto r: d ) {
			for ( auto s: r )
				fstar << setw(8) << s << " ";
			fstar << endl;
		}

		fstar << endl;
*/
		fstar << "\nloop_\n";
	
		for ( auto s: vs )
			fstar << "_" << s << "\n";	// Lead with an underscore

		for ( auto r: d ) {
			for ( auto s: r )
				fstar << setw(8) << s << " ";
			fstar << "\n";
		}

		fstar << "\n";

		return err;
	}
	map<string, int>&	tags() { return tg; }
	void			tags(map<string, int>& t) { tg = t; }
	vector< vector<string> >&	data() { return d; }
	long			tag(const string& t) { return tg.at(t); }
	long			find(const string& t) {
		auto itr = tg.find(t);
		if ( itr == tg.end() ) return -1;
		return itr->second;
	}
	long			columns() { return tg.size(); }
	long			rows() { return d.size(); }
	void			add_row(vector<string> vs) { d.push_back(vs); }
	vector<string>& add_row(long i) { vector<string> vs(i); d.push_back(vs); return d.back(); }
	vector<string>&	operator[](long i) { return d[i]; }
	string			value(long i, const string& t) {
		if ( tg.find(t) == tg.end() ) return "?";
		return d[i][tg[t]];
	}
	long			integer(long i, const string& t) {
		if ( tg.find(t) == tg.end() ) return 0;
		return stol(d[i][tg[t]]);
	}
	long			integer(long i, long j) {
		try { return stol(d[i][j]); }
		catch (...) { return 0; }
	}
	double			real(long i, const string& t) {
		if ( tg.find(t) == tg.end() ) return 0;
		return stod(d[i][tg[t]]);
	}
	double			real(long i, long j) {
		try { return stod(d[i][j]); }
		catch (...) { return 0; }
	}
	void			erase(const string& t) {
		long			j(tg[t]);
		tg.erase(t);
		for ( auto r: d ) r.erase(r.begin()+j);
	}
	void		replace_tag(string oldtag, string newtag) {
		if ( tg.find(oldtag) != tg.end() ) {
			tg[newtag] = tg[oldtag];
			tg.erase(oldtag);
		}
	}
	void			show_tags() {
		for ( auto i = tg.begin(); i != tg.end(); ++i )
			cout << "-" << i->first << "-" << tab << i->second << endl;
	}
} ;


/**
@class Bstar_block
@brief	Structure for a data block with multiple items in a STAR database.

	Each data block contains a set set of unique items defined by tags 
	and with the associate data as single or multiple values.
	The order of the items in the data block are important and preserved.
***/
class BstarBlock {
private:
	string				tg;		// Block identifier
	string				fn;		// File name for this block
	string				com;	// Comment following tag-value pair
	map<string,string>	it;		// List of data items
	vector<BstarLoop>	lp;		// List of loops
public:
	BstarBlock() {}
	BstarBlock(const string& s) {
		if ( regex_search(s, regex("^data_")) ) tg = s.substr(5);
		else tg = s;
		//	cout << "block " << tag << endl;
	}
	string			read(ifstream& fstar) {
		string			s, t;
		while ( !fstar.eof() && !regex_search(s, regex("^data_")) ) {
			s = get_clean_line(fstar);
//			cout << "=" << s << "=" << endl;
			if ( s.length() && s[0] != '#' ) {
				if ( s[0] == '_' ) {		// single item
					size_t		i(s.find_first_of(' '));
					if ( i != string::npos ) {
						t = s.substr(1, i-1);	// Strip off the underscore
						it[t] = quote_or_not(s.substr(i));
					} else {
						t = s.substr(1);		// Strip off the underscore
					}
				} else if ( s[0] == ';' ) {	// multiline
					string			v;
					if ( s.length() > 1 ) v = s.substr(1);
					while ( !fstar.eof() ) {
						s = get_clean_line(fstar);
						if ( s[0] == ';' ) break;
						v += s;
					};
					it[t] = v;
				} else if ( regex_search(s, regex("^loop_")) ) {
					lp.push_back(BstarLoop(fstar));
				}
			}
		}
		return s;
	}
	int				write(ofstream& fstar) {
		int						err(0), line_length(80);
/*
		fstar << endl << "data_" << tg << endl << endl;

		for ( auto ii: it ) {
			fstar << left << "_" << setw(39) << ii.first;	// Lead with an underscore
			string		s = ii.second;
			if ( s.length() < line_length ) {
				if ( s[0] != '"' && s.find_first_of(' ') != string::npos )
					fstar << " \"" << s << "\"" << endl;
				else fstar << " " << s << endl;
			} else {
				fstar << endl << ";" << endl;
				for ( size_t i = 0; i<s.length(); i += line_length )
					fstar << s.substr(i, line_length) << endl;
				fstar << ";" << endl;
			}
		}
*/
		fstar << "\n" << "data_" << tg << "\n\n";

		for ( auto ii: it ) {
			fstar << left << "_" << setw(39) << ii.first;	// Lead with an underscore
			string		s = ii.second;
			if ( s.length() < line_length ) {
				if ( s[0] != '"' && s.find_first_of(' ') != string::npos )
					fstar << " \"" << s << "\"\n";
				else fstar << " " << s << "\n";
			} else {
				fstar << "\n;\n";
				for ( size_t i = 0; i<s.length(); i += line_length )
					fstar << s.substr(i, line_length) << "\n";
				fstar << ";\n";
			}
		}

		for ( auto il: lp )
			err += il.write(fstar);
	
		return err;
	}
	int				write(string filename) {
		int				err(0);
		cout << "writing " << filename << endl;
		string			comment = "# Written by Bsoft\n";
		
		fn = filename;
	
		ofstream		fstar(filename.c_str());
		if ( fstar.fail() ) {
			cerr << "Error: Not able to write " << filename << endl;
			return -1;
		}

		fstar << comment << "\n";

		write(fstar);
	
		fstar.close();
	
		return err;
	}
	string			tag() { return tg; }
	void			file_name(string s) { fn = s; }
	string			file_name() { return fn; }
	string&			operator[](const string& t) { return it[t]; }
	bool			exists(const string& t) {
		if ( it.find(t) == it.end() ) return 0;
		return 1;
	}
	bool			exists_loop(const string& t) {
		for ( auto il: lp )
			if ( il.find(t) >= 0 ) return 1;
		return 0;
	}
	string			at(const string& t) {
//		cout << "-" << t << "-" << endl;
		if ( it.find(t) == it.end() ) return "?";
		return it.at(t);
	}
	long			integer(const string& t) {
//		cout << "-" << t << "-" << endl;
		if ( it.find(t) == it.end() ) return 0;
		return stol(it.at(t));
	}
	double			real(const string& t) {
//		cout << "-" << t << "-" << endl;
		if ( it.find(t) == it.end() ) return 0;
//		return stod(it.at(t));
		try { return stod(it[t]); }
		catch (...) { return 0; }
	}
	map<string,string>&	items() { return it; }
	void			erase(const string& t) {
		if ( it.find(t) != it.end() ) it.erase(t);
	}
	vector<BstarLoop>&	loops() { return lp; }
	BstarLoop&		loop(long i) {
		long			j(0);
		for ( auto il = lp.begin(); il != lp.end(); ++il )
			if ( j==i ) return *il;
		return lp.front();
	}
	BstarLoop&		add_loop() {
		lp.push_back(BstarLoop());
		return lp.back();
	}
	void		replace_tag(string oldtag, string newtag) {
		if ( it.find(oldtag) != it.end() ) {
			it[newtag] = it[oldtag];
			it.erase(oldtag);
		} else {
			for ( auto il: lp ) il.replace_tag(oldtag, newtag);
		}
	}
	void			show_tags() {
		for ( auto i = it.begin(); i != it.end(); ++i )
			cout << "-" << i->first << "-" << tab << i->second << endl;
	}
} ;

/**
@class Bstar
@brief	Overall STAR class to hold all data blocks and options for I/O..

	The split flag allows the user to output data blocks in separate files
	in stead of one big file.
	The line length field allows the user to output long lines without
	wrapping it around.
	The comments are ignored but output to the a new file - this can be used
	to document the history of the file.
	The STAR database is a hierarchy consisting of blocks, each with a set
	of items.
***/
class Bstar {
private:
	size_t 				lin;	// Length of lines in output file
	string				com;	// List of comments before 1st data block
	vector<BstarBlock>	b;		// List of data blocks
public:
	Bstar() { lin = 80; }
	Bstar(string filename) { read(filename); }
	/**
	@brief 	Reads paramaters and data into a STAR data base from a file.
	@param	filename	a file name.

		Every data block is read separately and comments are preserved as far 
		as possible.

	**/
	int			read(string filename) {
//		cout << "reading " << filename << endl;
   		ifstream		fstar(filename.c_str());
		if ( fstar.fail() ) {
//			error_show(thisfile->c_str(), __FILE__, __LINE__);
			return -1;
		}
		string			s;
		int				comment_add(b.size()==0);

		s = get_clean_line(fstar);

		while ( !fstar.eof() ) {
//			cout << "-" << s << "-" << endl;
			if ( regex_search(s, regex("^data_")) ) {
//				cout << "parsing block " << s << endl;
				b.push_back(BstarBlock(s));
				b.back().file_name(filename);
				s = b.back().read(fstar);
				comment_add = 0;
			} else {
				if ( comment_add )
					com = com + s + "\n";
				s = get_clean_line(fstar);
			}
		}

		fstar.close();
//		cout << "done " << filename << endl;
		
		return 0;
	}
	/**
	@brief 	Writes a STAR data base to a STAR format file.
	@param	filename	file name.
	@return int			error code (<0 means failure).

	**/
	int			write(string filename) {
		int				err(0);
	
//		cout << "writing " << filename << endl;

//		if ( com.length() < 1 ) com = "# Written by Bsoft";
		if ( com.length() < 1 ) com = command_line_time2();
		
		ofstream		fstar(filename.c_str());
		if ( fstar.fail() ) {
			cerr << "Error: Not able to write " << filename << endl;
			return -1;
		}

		fstar << com << "\n";

		for ( auto ib: b )
			err += ib.write(fstar);
		
		fstar.close();
	
		return err;
	}
	int			write(string filename, int split) {
		if ( split == 0 ) return write(filename);

		int				err(0), i(0);
		string			blockname;

		for ( auto ib: b ) {
			if ( split == 9 ) {
				if ( ib.tag().length() > 0 ) {
					blockname = ib.tag();
//					remove_spaces(blockname);
					remove_if(blockname.begin(), blockname.end(), ::isspace);
					blockname += ".star";
				} else {
					blockname = ".star";
					blockname = insert(blockname, ++i, 4);
				}
			} else {
				blockname = insert(filename, ++i, 4);
		    }
			err += ib.write(blockname);
		}
		
		return err;
	}
	void		line_length(long i) { lin = i; }
	long		line_length() { return lin; }
	void		comment(const string& s) { com = s; }
	string&		comment() { return com; }
	long		list_comments(long len) {
		if ( len < 10 ) len = 1000000;
		stringstream 	ss(com);
		string 			to;
		long			n(0);
		while( getline(ss,to,'\n') ) {
			if ( to[0] == '#' && to.length() > 1 ) {
				cout << to.substr(2,len) <<endl;
				n++;
			}
		}
		return n;
	}
	vector<BstarBlock>&	blocks() { return b; }
	bool	exists(string t) {
		for ( auto ib: b )
			if ( ib.exists(t) ) return 1;
		return 0;
	}
	BstarBlock&	block(long i) {
		long			j(0);
		for ( auto ib = b.begin(); ib != b.end(); ++ib, ++j )
			if ( j==i ) return *ib;
		return b.front();
	}
	BstarBlock&	find(string t) {
		for ( auto &ib: b )
			if ( ib.exists(t) ) return ib;
		return b.front();
	}
	BstarBlock&	add_block(const string& s) {
		b.push_back(BstarBlock(s));
		return b.back();
	}
	long		number_of_blocks(string t) {
		long			n(0);
		for ( auto ib: b )
			if ( ib.exists(t) ) n++;
		return n;
	}
	void		replace_tag(string oldtag, string newtag) {
		for ( auto ib: b ) ib.replace_tag(oldtag, newtag);
	}
	void		erase(const string& s) {
		for ( auto ib = b.begin(); ib != b.end(); ++ib ) {
			if ( ib->tag() == s ) {
				b.erase(ib);
				break;
			}
		}
	}
} ;
#define _STAR_
#endif


