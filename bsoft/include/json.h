/**
@file	json.h
@author	Bernard Heymann
@date	20160409 - 20191219

**/

#include <cstddef>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <type_traits>

using namespace std;

#ifndef _JSON_
#define _JSON_

enum JStype {
	JSstring,
	JSbool,
	JSinteger,
	JSreal,
	JSobject,
	JSarray,
	JSnull
};

// Prototypes
class 		JSvalue;
ostream&	operator<<(ostream& output, JStype t);
ostream&	operator<<(ostream& output, JSvalue& jsv);

/*
	The string s contains eiher the tag if the type is a pair
	or the value if the type is a number or string.
	The list v contains elements that could be:
		a single value if type is a pair
		an array of values if type is an array
		an array of pairs if type is an object.
	The last three types are the barewords: true, false and null
	Type		Str/Val	Vec/Map
	JSstring	value	-
	JSbool		t/f		-
	JSinteger	long	-
	JSreal		double	-
	JSobject	-		map
	JSarray		-		vec
	JSnull		-		-
*/
class JSvalue {					// 104
private:
	JStype					t;	// 4 - aligned on 8
	string					s;	// 24
	bool					b;	// 1 - aligned on 8
	long					i;	// 8
	double					r;	// 8
	vector<JSvalue>			v;	// 24
	map<string, JSvalue>	m;	// 24
	void			clean_string(string val) {
		if ( val.front() == '"' ) {
			 if ( val.back() == '"' ) s = val.substr(1,val.length() - 2);
			 else s = val.substr(1);
		} else s = val;
	}
	void			fail(string msg) {
		cerr << "JSON error! " << msg << " (" << t << ")" << endl;
		show_value(*this);
		cerr << "array:" << endl;
		for ( auto it = v.begin(); it != v.end(); ++it )
			show_value(*it);
		cerr << "object:" << endl;
		for ( auto it = m.begin(); it != m.end(); ++it ) {
			cerr << "tag: " << it->first << endl;
			show_value(it->second);
		}
//		exit(-1);
	}
	void			show_value(JSvalue& val) {
		cerr << "string: \"" << s << "\"" << endl;
		cerr << "boolean: " << b << endl;
		cerr << "integer: " << i << endl;
		cerr << "real:    " << r << endl;
	}
public:
	JSvalue() : t(JSnull), b(0), i(0), r(0) { }
	JSvalue(JStype tp) : t(tp) { }
	JSvalue(bool tf) : t(JSbool), b(tf) { }
	JSvalue(const char* val) : t(JSstring) { clean_string(val); }
	JSvalue(const string& val) : t(JSstring) { clean_string(val); }
	JSvalue(const int& val) : t(JSinteger) { i = val; }
	JSvalue(const unsigned int& val) : t(JSinteger) { i = val; }
	JSvalue(const long& val) : t(JSinteger) { i = val; }
	JSvalue(const unsigned long& val) : t(JSinteger) { i = val; }
	JSvalue(const float& val) : t(JSreal) { r = val; }
	JSvalue(const double& val) : t(JSreal) { r = val; }
	template <typename T>
	JSvalue(const vector<T>& arr) : t(JSarray) {
		v.clear();
		for ( auto it = arr.begin(); it != arr.end(); ++it )
			v.push_back(*it);
	}
//	JSvalue(const map<string, JSvalue>& mp) : t(JSobject) { m = mp; }
	template <typename T>
	JSvalue(const map<string, T>& mp) : t(JSobject) {
		m.clear();
		for ( auto it = mp.begin(); it != mp.end(); ++it )
			m[it->first] = JSvalue(it->second);
	}
	const JStype&			type() const { return t; }
	// Size
	size_t					size() {
		if ( t == JSarray ) return v.size();
		else if ( t == JSobject ) return m.size();
		else if ( t == JSstring ) return s.size();
		else if ( t == JSnull ) return 0;
		return 1;
	}
	void					clear() {
		if ( t == JSarray ) return v.clear();
		else if ( t == JSobject ) return m.clear();
		else if ( t == JSstring ) return s.clear();
	}
	long					memory() {
		long			n(sizeof(JSvalue));
		if ( t == JSarray ) for ( auto it: v ) n += it.memory();
		else if ( t == JSobject ) for ( auto it: m ) n += it.first.size() + it.second.memory();
		else if ( t == JSstring ) n += s.size();
		return n;
	}
	// Access
	string&					value() {
		if ( t != JSstring ) fail("value: Expected a string.");
		return s;
	}
	bool&					boolean() {
		if ( t != JSbool ) fail("boolean: Expected a boolean.");
		return b;
	}
	long&					integer() {
		if ( t != JSinteger ) {
			if ( t != JSreal ) fail("integer: Expected an integer.");
			else i = r;
		}
		return i;
	}
	double&					real() {
		if ( t != JSreal ) {
			if ( t != JSinteger ) fail("real: Expected a real number.");
			else r = i;
		}
		return r;
	}
	// Array access
	vector<JSvalue>&		array() {
		if ( t != JSarray ) fail("array: Expected an array.");
		return v;
	}
	JSvalue&				operator[](long j) {
		if ( t != JSarray ) fail("Expected an array.");
		return v.at(j);
	}
	vector<long>			array_integer() {
		if ( t != JSarray ) fail("array_integer: Expected an array.");
		vector<long> 		vec;
		for ( auto it = v.begin(); it != v.end(); ++it )
			vec.push_back(it->integer());
		return vec;
	}
	vector<double>			array_real() {
		if ( t != JSarray ) fail("array_real: Expected an array.");
		vector<double> 		vec;
		for ( auto it = v.begin(); it != v.end(); ++it )
			vec.push_back(it->real());
		return vec;
	}
	JSvalue&					front() {
		if ( t != JSarray ) fail("front: Expected an array.");
		return v.front();
	}
	JSvalue&					back() {
		if ( t != JSarray ) fail("back: Expected an array.");
		return v.back();
	}
	vector<JSvalue>::iterator	begin() { return v.begin(); }
	vector<JSvalue>::iterator	end() { return v.end(); }
	vector<JSvalue>::iterator	erase(long j) {
		if ( t != JSarray ) fail("erase: Expected an array.");
		return v.erase(v.begin()+j);
	}
	vector<JSvalue>::iterator	erase(vector<JSvalue>::iterator first, vector<JSvalue>::iterator last) {
		if ( t != JSarray ) fail("erase: Expected an array.");
		return v.erase(v.begin(), v.end());
	}
	void					push_back(const JSvalue& val) {
		if ( t != JSarray ) fail("push_back: Expected an array.");
		v.push_back(val);
	}
	// Object access
	map<string, JSvalue>&	object() {
		if ( t != JSobject ) fail("Expected an object.");
		return m;
	}
	JSvalue&				operator[](string tag) {
		if ( t != JSobject ) fail("Expected an object with tag " + tag + ".");
		if ( m.find(tag) == m.end() ) m[tag] = JSvalue();
		return m[tag];
	}
	JSvalue&				object(string tag) {
		if ( t != JSobject ) fail("object: Expected an object with tag " + tag + ".");
		if ( m.find(tag) == m.end() ) fail("Tag " + tag + " not found.");
		return m[tag];
	}
	JSvalue&				at(string tag) {
		if ( t != JSobject ) fail("at: Expected an object with tag " + tag + ".");
		if ( m.find(tag) == m.end() ) fail("Tag " + tag + " not found.");
		return m[tag];
	}
	bool					exists(string tag) {
		if ( t != JSobject ) fail("exists: Expected an object with tag " + tag + ".");
		if ( m.find(tag) == m.end() ) return false;
		return true;
	}
	long					erase(string tag) {
		if ( t != JSobject ) fail("erase: Expected an object with tag " + tag + ".");
		return m.erase(tag);
	}
	map<string, JSvalue>::iterator	object_begin() { return m.begin(); }
	map<string, JSvalue>::iterator	object_end() { return m.end(); }
	map<string, JSvalue>::iterator	erase(map<string, JSvalue>::iterator first, map<string, JSvalue>::iterator last) {
		if ( t != JSobject ) fail("erase: Expected an object.");
		return m.erase(first, last);
	}
	void				append(JSvalue& obj) {
		if ( t != JSobject ) fail("append: Expected an object.");
		if ( obj.t != JSobject ) fail("append: Expected an object.");
		m.insert(obj.object_begin(), obj.object_end());
	}
	// Simple type access
	bool&				boolean(string& tag) {
		return object(tag).boolean();
	}
	long&				integer(string& tag) {
		return object(tag).integer();
	}
	double&				real(string& tag) {
		return object(tag).real();
	}
	vector<JSvalue>&	array(string& tag) {
		return object(tag).array();
	}
	// Relational operators
	bool				operator==(string val) {
		return s == val;
	}
	bool				operator!=(string val) {
		return s != val;
	}
	// Computational operators
	void				operator+=(double d) {
		if ( t == JSinteger ) i += d;
		else if ( t == JSreal ) r += d;
		else if ( t == JSarray )
			for ( auto it = this->begin(); it != this->end(); ++it )
				*it += d;
	}
	void				operator*=(double d) {
		if ( t == JSinteger ) i *= d;
		else if ( t == JSreal ) r *= d;
		else if ( t == JSarray )
			for ( auto it = this->begin(); it != this->end(); ++it )
				*it *= d;
	}
	void				operator/=(double d) {
		if ( d == 0 ) {
			cerr << "Error: Divisor is zero!" << endl;
			return;
		}
		if ( t == JSinteger ) i /= d;
		else if ( t == JSreal ) r /= d;
		else if ( t == JSarray )
			for ( auto it = this->begin(); it != this->end(); ++it )
				*it /= d;
	}
	// Queries
	vector<JSvalue*>	operator()(string& jsonpath) {
//		cout << "Parsing path: " << jsonpath << endl;
		vector<JSvalue*>		qarr;
		if ( jsonpath.find_first_of(".") < 2 ) qdots(jsonpath, qarr);
		else qbrackets(jsonpath, qarr);
//		cout << "Query size: " << qarr.size() << endl;
		return qarr;
	}
	void				qdots(string& jsonpath, vector<JSvalue*>& qarr) {
		//		cout << "jsonpath=" << jsonpath << endl;
		size_t			b1(jsonpath.find_first_of("."));
		size_t			b2(jsonpath.find_first_of("["));
		//		cout << b1 << "\t" << b2 << endl;
		if ( b2 < b1 ) {
			b1 = b2;
			b2 = jsonpath.find_first_of("]", b1);
		} else {
			size_t			b3(jsonpath.find_first_of(".", b1+1));
			if ( b3 < b2 ) b2 = b3;
		}
		if ( b1 >= string::npos ) {
			qarr.push_back(this);
			return;
		}
		//		cout << b1 << "\t" << b2 << endl;
		if ( b2 >= string::npos ) b2 = jsonpath.length();
		string			tag(jsonpath.substr(b1+1,b2-b1-1));
		string			nupath(jsonpath.substr(b2));	// new path starts with '.' or '['
		//		cout << "tag=" << tag << endl;
		//		cout << "nupath=" << nupath << endl;
		if ( t == JSarray ) {
			if ( tag == "*" ) {
				for ( auto it = begin(); it != end(); ++it )
					it->qdots(nupath, qarr);
			} else {
				try { i = stol(tag); }
				catch (...) { fail("No conversion to an integer."); }
				v[i].qdots(nupath, qarr);
			}
		} else if ( t == JSobject ) {
			if ( tag == "*" ) {
				for ( auto it = object_begin(); it != object_end(); ++it )
					it->second.qdots(nupath, qarr);
			} else {
				if ( tag[0] == '\'' )
					tag = tag.substr(1, tag.find_last_of("'") - 1);
				if ( tag[0] == '"' )
					tag = tag.substr(1, tag.find_last_of("\"") - 1);
				m[tag].qdots(nupath, qarr);
			}
		} else fail("Tag " + tag + " not found.");
	}
	void				qbrackets(string& jsonpath, vector<JSvalue*>& qarr) {
//		cout << "jsonpath=" << jsonpath << endl;
		size_t			b1(jsonpath.find_first_of("["));
		if ( b1 == string::npos ) {
			qarr.push_back(this);
			return;
		}
		size_t			b2(jsonpath.find_first_of("]"));
		string			tag(jsonpath.substr(b1+1,b2-b1-1));
		string			nupath(jsonpath.substr(b2+1));
//		cout << "tag=" << tag << endl;
		if ( t == JSarray ) {
			if ( tag == "*" ) {
				for ( auto it = begin(); it != end(); ++it )
					it->qbrackets(nupath, qarr);
			} else {
				try { i = stol(tag); }
				catch (...) { fail("No conversion to an integer."); }
				v[i].qbrackets(nupath, qarr);
			}
		} else if ( t == JSobject ) {
			if ( tag == "*" ) {
				for ( auto it = object_begin(); it != object_end(); ++it )
					it->second.qbrackets(nupath, qarr);
			} else {
				if ( tag[0] == '\'' )
					tag = tag.substr(1, tag.find_last_of("'") - 1);
				if ( tag[0] == '"' )
					tag = tag.substr(1, tag.find_last_of("\"") - 1);
				m[tag].qbrackets(nupath, qarr);
			}
		} else fail("Tag " + tag + " not found.");
	}
	void					write(string filename) {
		ofstream		fjs(filename);
		fjs << *this << endl;
		fjs.close();
	}
};

/*
	The following incorrect syntax is ignored:
	Any content after the top-level closing bracket
	Any fields before or after commas that only contain white space
	Any white space within numeric fields or leading zeroes
*/
class JSparser {
private:
	string				content;
	string::iterator	it;
	void		read(string filename) {
//		cout << "Reading " << filename << endl;
		ifstream	f(filename);
		if ( f.fail() ) {
			cerr << "Error: File " << filename << " not opened" << endl;
			return;
		}
		f.seekg(0, ios::end);
		content.resize(f.tellg());
//		cout << "Length = " << content.size() << endl;
		f.seekg(0, ios::beg);
		f.read(&content[0], content.size());
		f.close();
	}
	void		fail(string msg) {
		long			i;
		cerr << "JSON parser error! " << msg << " (";
		for ( i=0; it!=content.end() && i<10 ; ++it, ++i ) cerr << *it;
		if ( i >= 10 ) cerr << "...";
		cerr << ")" << endl;
		exit(-1);
	}
public:
	JSparser() {};
	JSparser(string filename) { read(filename); }
	JSparser(const char* buf, size_t n) { content = string(buf, n); }
	JSvalue		parse() {
		JSvalue		val;
		for ( it=content.begin(); it!=content.end(); ++it ) {
			if ( *it == '{' ) return parseJSobject();
			else if ( *it == '[' ) return parseJSarray();
			else if ( !isspace(*it) ) break;
		}
		fail("No { or [ found at the start.");
		return val;
	}
	JSvalue		parse(string filename) {
		read(filename);
		return parse();
	}
private:
	JSvalue		parseJSvalue() {
		JSvalue		val;
		for ( ; it!=content.end() && *it != '}' && *it != ']' && *it != ','; ++it ) {
			if ( *it == '{' ) {
				return parseJSobject();
			} else if ( *it == '[' ) {
				return parseJSarray();
			} else if ( *it == '"' ) {
				return parseJSstring();
			} else if ( isdigit(*it) || *it == '.'  || *it == '+'  || *it == '-' ) {
				return parseJSnumber();
			} else if ( *it == 't' ||  *it == 'f' ||  *it == 'n' ) {
				return parseJSbare();
			} else if ( !isspace(*it) ) {
				fail("Found an unquoted string.");
			}
		}
		return val;
	}
	string		parseString() {
		string		s;
		for ( ; it!=content.end() && *it != '"'; ++it )
			if ( !isspace(*it) ) fail("Unquoted string.");
		if ( *it != '"' ) fail("No opening quote on a string.");
		for ( ++it; it!=content.end(); ++it ) {
			if ( *it == '\\' ) { s += *it; ++it; s += *it; }
			else if ( *it == '"' ) break;
			else s += *it;
//			cout << *it;
		}
//		cout << s << " \"" << *it << "\"" << endl;
		if ( *it != '"' ) fail("No closing quote on a string.");
		++it;
		return s;
	}
	JSvalue		parseJSstring() {
		string			s(parseString());
		return JSvalue(s);
	}
	JSvalue		parseJSnumber() {
		string			s;
		int				ft(0);
		for ( ; it!=content.end() && *it != '}' && *it != ']'
				&& *it != ',' && !isspace(*it); ++it ) {
			if ( *it == '.' || *it == 'e' ) ft = 1;
			if ( *it == 'x' ) fail("Numbers cannot be hex");
			s += *it;
		}
		long			l;
		double			d;
		if ( ft ) {
			try { d = stod(s); }
			catch (...) { fail("No conversion to a real number."); }
			return JSvalue(d);
		} else {
			try { l = stol(s); }
			catch (...) { fail("No conversion to an integer."); }
			return JSvalue(l);
		}
	}
	JSvalue		parseJSbare() {
		string			s;
		for ( ; it!=content.end() && !isspace(*it) && *it != ',' && *it != '}'; ++it )
			s += *it;
		if ( s.compare("true") == 0 ) {
			return JSvalue(true);
		} else if ( s.compare("false") == 0 ) {
			return JSvalue(false);
		} else if ( s.compare("null") == 0 ) {
			return JSvalue();
		} else fail("Bare word " + s + " not allowed!");
		return JSvalue();
	}
	JSvalue		parseJSarray() {
		it++;
		JSvalue			arr(JSarray);
		for ( ; it!=content.end() && *it != ']'; ++it ) {
//			cout << *it << endl;
			if ( !isspace(*it) && *it != '}' && *it != ',' )
				arr.push_back(parseJSvalue());
//			cout << "parseJSarray:'" << *it << "'" << endl;
			if ( *it == ':' ) fail("Unexpected colon.");
			if ( it==content.end() || *it == ']' ) break;
		}
		if ( *it != ']' ) fail("No closing bracket on an array.");
		it++;
//		cout << "Parsing array" << *it << endl;
		return arr;
	}
	JSvalue		parseJSobject() {
		it++;
		string					tag;
		JSvalue					obj(JSobject);
		for ( ; it!=content.end() && *it != '}'; ++it ) {
			if ( !isspace(*it) && *it != ']' && *it != ',' ) {
				tag = parseString();
				for ( ; it!=content.end() && isspace(*it); ++it ) ;
				if ( *it != ':' ) fail("Expected a colon.");
				it++;
				obj[tag] = parseJSvalue();
			}
			if ( it==content.end() || *it == '}' ) break;
		}
		if ( *it != '}' ) fail("No closing bracket on an object.");
		it++;
		return obj;
	}
};

#endif	// _JSON_




