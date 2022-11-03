/**
@file	string_util.h
@author	Bernard Heymann
@date	20160911 - 20220101

**/


#ifndef _STRINGUTIL_

#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <complex>

using namespace std;

// Function prototypes
string		command_line(int argc, char** argv);
string		base(string& filename);
string		extension(string filename);
string		replace_extension(string& filename, string ext);
string		insert(string& filename, string ins);
string		insert(string& filename, int i, int n);
void		remove_spaces(string& s);
string		remove_spaces2(string& s);
string		parameter_file_path(string& filename);

string 		quote_or_not(const string &s);

/**
@brief  Splits a string on whitespace.
@param  &s                      string to be split.
@return vector<string>  vector of strings.
**/
vector<string> split(const string &s);
vector<string> splitn(const string &s, long n);

/**
@brief 	Splits a string using a given delimiter.
@param	&s				string to be split.
@param	delim			delimiter.
@return	vector<string>	vector of strings.
**/
vector<string> split(const string &s, char delim);

string	 concatenate(const vector<string>& vs);
string	 concatenate(const vector<long>& v);
string	 concatenate(const vector<double>& v);

vector<long>	parse_integer_vector(string vecstr);
vector<double>	parse_real_vector(string vecstr);

string		local_time();

string		get_clean_line(ifstream& fstar);

long		to_integer(string s);

double		to_real(string s);

bool		check_for_number(string& s);

int			check_for_type(string& s);

#define _STRINGUTIL_
#endif
