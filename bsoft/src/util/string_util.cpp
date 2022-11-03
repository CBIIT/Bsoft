/**
@file	string_util.cpp
@author	Bernard Heymann
@date	20160911 - 20220101

**/


#include "string_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


string		command_line(int argc, char** argv)
{
	string		cl(argv[0]);
	
	for ( long i=1; i<argc; ++i ) {
		cl += " ";
		cl += argv[i];
	}

	return cl;
}

string		base(string& filename) {
	return filename.substr(0, filename.rfind("."));
}

string		extension(string filename)
{
	string		ext(filename.substr(filename.rfind(".")+1));
        
	transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        
	return ext;
}

string		replace_extension(string& filename, string ext)
{
	string		newname = filename.substr(0, filename.rfind(".")+1);
	
	return newname + ext;
}

string		insert(string& filename, string ins)
{
	string		newname = filename.substr(0, filename.rfind("."));
	string		ext = filename.substr(filename.rfind(".")+1);
	
	newname += ins + ext;
	
	return newname;
}

string		insert(string& filename, int i, int n)
{
	string		ins = to_string(i);
	if ( ins.size() < n ) ins = string(n-ins.size(),'0') + ins;
	
	string		newname = filename.substr(0, filename.rfind("."));
	string		ext = filename.substr(filename.rfind(".")+1);
	
	newname += ins + ext;
	
	return newname;
}

void		remove_spaces(string& s)
{
	remove_if(s.begin(), s.end(), ::isspace);
}
	
string		remove_spaces2(string& s)
{
	long		i, j(0);
	string		ns(s);
	
	for ( i=0; i<s.size(); ++i )
		if ( !isspace(s[i]) )
			ns[j++] = s[i];
	
	ns.resize(j);
	
	return ns;
}

/**
@brief 	Finds the parameter file path.
@param 	&filename	the parameter filename.
@return Bstring				the full path and file name.

	The parameter file path should primarily be defined in the environmental
	variable "BPARAM".
	Otherwise a default is used which may or may not be valid.
	The string returned has a length containing the path plus one byte.

**/
string		parameter_file_path(string& filename)
{
	string		path;
	
	if ( getenv("BPARAM") ) {
		path = getenv("BPARAM");
 		if ( path[-1] != ']'  && path[-1] != '/')
 			path += "/";
	} else
		path = "/usr/local/bsoft/parameters/";
	path += filename;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG parameter_file_path: " << path << endl;
	
	return path;
}

string 		quote_or_not(const string &s)
{
	if ( s.length() < 2 ) return s;
	long				i, j, k, m(s.size()), f(0);
	const char*			cs = s.c_str();
	for ( i=j=k=0; i<m; ++i ) {
		if ( cs[i] == '"' ) f = 1 - f;
		if ( !f && isspace(cs[i]) ) {
			if ( k ) break;
		} else {
			if ( !k ) j = i;
			k++;
		}
	}
    return s.substr(j, k);
}


/**
@brief  Splits a string on whitespace.
@param  &s			string to be split.
@return vector<string>  vector of strings.
**/
vector<string> split(const string &s)
{
    istringstream		ss(s);
    vector<string>		sv((istream_iterator<string>(ss)),
                                 istream_iterator<string>());
    return sv;
}

/**
@brief  Splits a string on whitespace.
@param  &s			string to be split.
@param  n 			number of strings to expect.
@return vector<string>  	vector of strings.

	Keeps whitespace within quoted strings intact.
	Excess beyond n strings is discarded.
	
**/
vector<string> splitn(const string &s, long n)
{
	vector<string>		sv(n);
//	bool				q(0);
	long				i, j, k, l, m(s.size());
	const char*			cs = s.c_str();
	for ( i=j=k=l=0; i<m && l<n; ++i ) {
		if ( cs[i] == '"' ) {
			for ( j = ++i; cs[i] != '"' && i<m; ++i );
			k = i-j;
			i++;
		}
		if ( isspace(cs[i]) ) {	// End of string
			if ( k ) {
				sv[l++] = s.substr(j, k);
				k = 0;
			}
		} else {
			if ( !k ) j = i;	// Start of string
			k++;
		}
	}
	if ( l < n && k ) sv[l] = s.substr(j, k);	// Final string
    return sv;
}
/*
vector<string> splitn(const string &s, long n)
{
	vector<string>		sv(n);
	bool				q(0);
	long				i, j, k, l, m(s.size());
	const char*			cs = s.c_str();
	for ( i=j=k=l=0; i<m && l<n; ++i ) {
		if ( cs[i] == '"' ) q = 1 - q;	// Toggle on quotes
		if ( !q && isspace(cs[i]) ) {	// End of string if not within quotes
			if ( k ) {
				sv[l++] = s.substr(j, k);
				k = 0;
			}
		} else {
			if ( !k ) j = i+q;	// Start of string
			k++;
		}
	}
	if ( l < n && k ) sv[l] = s.substr(j, k);	// Final string
    return sv;
}
*/
/**
@brief 	Splits a string using a given delimiter.
@param	&s				string to be split.
@param	delim			delimiter.
@return	vector<string>		vector of strings.
**/
vector<string> split(const string &s, char delim)
{
    stringstream 	ss(s);
    string 			item;
    vector<string> 	tokens;
	
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
	
    return tokens;
}

/**
@brief 	Concatenates a vector of strings.
@param	&vs				vector of strings.
@return	string			concatenated strings.
**/
string	 concatenate(const vector<string>& vs)
{
	long			i(0), len(0);
    
    for ( auto s: vs ) len += s.size() + 1;
    
     string 		cs(len,' ');
 
     for ( auto s: vs ) {
		cs.replace(i, s.size(), s);
		i += s.size() + 1;
	}

    return cs;
}

/**
@brief 	Concatenates a vector of integers into a comma-separated array.
@param	&v				vector of integers.
@return	string			concatenated strings.
**/
string	 concatenate(const vector<long>& v)
{
	string			s;
	
	if ( v.size() < 1 ) return s;
	
	s = to_string(v[0]);
	
	for ( long i=1; i<v.size(); ++i ) s += "," + to_string(v[i]);
	
	return s;
}

/**
@brief 	Concatenates a vector of real values into a comma-separated array.
@param	&v				vector of real values.
@return	string			concatenated strings.
**/
string	 concatenate(const vector<double>& v)
{
	string			s;
	
	if ( v.size() < 1 ) return s;
	
	s = to_string(v[0]);
	
	for ( long i=1; i<v.size(); ++i ) s += "," + to_string(v[i]);
	
	return s;
}


/**
@brief 	Parses a comma-separated array of integers into a vector.
@param	vecstr			string.
@return	vector<long>		vector of integers.
**/
vector<long>	parse_integer_vector(string vecstr)
{
	vector<string>	sv(split(vecstr, ','));
	
	vector<long>	vec;
	
	for ( auto it = sv.begin(); it != sv.end(); ++it )
		vec.push_back(stol(*it));
	
	return vec;
}

/**
@brief 	Parses a comma-separated array of real values into a vector.
@param	vecstr			string.
@return	vector<double>		vector of real values.
**/
vector<double>	parse_real_vector(string vecstr)
{
	vector<string>	sv(split(vecstr, ','));
	
	vector<double>	vec;
	
	for ( auto it = sv.begin(); it != sv.end(); ++it )
		vec.push_back(stod(*it));
	
	return vec;
}

string		local_time()
{
	chrono::system_clock::time_point t = chrono::system_clock::now();

	time_t		tt = chrono::system_clock::to_time_t(t);
	
	string		st = ctime(&tt);
	
	st.erase(st.length()-1);
	
//  cout << "Time: " << st << "." << endl;
	
  return st;
}

string	get_clean_line(ifstream& fstar)
{
	string			s;
	
	getline(fstar, s);
	
	int				n(s.size());

	if ( !n ) return s;

	int				i(0), j(0), k(0);

	// Find first non-whitespace
	while ( isspace(s[j]) ) j++;

	for ( i=0; j<n; ++i, ++j ) {
		if ( isspace(s[j]) ) {
			s[i] = ' ';		// All control characters become spaces
		} else {
			k = i;			// Set to the last non-space character
			s[i] = s[j];
		}
	}
	k++;
	
	s.resize(k);
	
//	cout << "-" << s << "-" << endl;

	return s;
}

long		to_integer(string s)
{
	try { return stol(s); }
	catch (...) { }
	return 0;
}

double		to_real(string s)
{
	try { return stod(s); }
	catch (...) { }
//	catch (...) { cout << "Exception: " << s << endl; }
	return 0;
}

// Check if string or number
bool		check_for_number(string& s)
{
	long		k, digit(0);
	for ( k=0; k<s.size(); k++ ) {
		if ( isalpha(s[k]) ) return 0;
		if ( isdigit(s[k]) ) digit++;
	}
	if ( digit ) return 1;
	return 0;
}

// Check if string, integer or real
int			check_for_type(string& s)
{
	long		k, digit(0), period(0);
	for ( k=0; k<s.size(); k++ ) {
		if ( isalpha(s[k]) ) return 0;
		if ( isdigit(s[k]) ) digit++;
		if ( s[k] == '.' ) period++;
	}
	if ( digit ) {
		if ( period ) return 2;
		return 1;
	}
	return 0;
}

