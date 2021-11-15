/**
@file	Bstring.h 
@brief	Header file for the string class 
@author Bernard Heymann 
@date	Created: 20051020
@date	Modified: 20190208
**/

#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

#ifndef _Bstring_
#define _Bstring_

#define MAX_FORMAT_SIZE	64

/************************************************************************
@Object: class Bstring
@Description:
	Bsoft linkable string class.
@Features:
	The variables are a public pointer to the next string and a private 
	pointer to a character string.
	The string data are allocated internally.
*************************************************************************/
class Bstring {
public:
	Bstring*	next;
private:
	char*		data;
public:
	Bstring() { data = NULL; next = NULL; }
	Bstring(const Bstring& s);
	Bstring(const string& s);
	Bstring(const Bstring& s, long start, long len);
	Bstring(const char* c);
	Bstring(const char* c, long start, long len);
	Bstring(const char c, long size);
	string		str() {
		string s;
		if ( data ) s = data;
		return s;
	}
	template <typename T>
	Bstring(T value, const char* format) {
		data = NULL;
		char*	temp = new char[MAX_FORMAT_SIZE];
		::snprintf(temp, MAX_FORMAT_SIZE, format, value);
		*this = temp;
		delete[] temp;
		next = NULL;
	}
	~Bstring() { delete[] data; data = NULL; }
	Bstring& operator=(const Bstring& s);
	Bstring& operator=(const string& s);
	Bstring& operator=(const char* c);
	Bstring& operator+=(const Bstring& s);
	Bstring operator+(const Bstring& s);
	Bstring operator+(const Bstring& s) const;
	Bstring& operator+=(const char c);
	Bstring operator+(const char c);
	bool operator==(const Bstring& s);
	bool operator!=(const Bstring& s);
	bool operator>(const Bstring& s);
	bool operator<(const Bstring& s);
	bool operator>=(const Bstring& s);
	bool operator<=(const Bstring& s);
	Bstring operator<<(int n);
	Bstring operator>>(int n);
	char& operator[](long i);
	const char* c_str() { return (const char*)data; }
	bool	empty() { return data == NULL; }
	long	length() const;
	long	count(const char c);
	long	compare(const Bstring& s);
	long	compare_value(const Bstring& s);
	Bstring	common(const Bstring& s);
	void	fill(const char c);
	void	place(long i, const Bstring& s);
	Bstring lower();
	Bstring upper();
	Bstring substr(long start, long len);
	Bstring no_space();
	Bstring no_lead_space();
	Bstring alnum();
	Bstring left(long len);
	Bstring right(long len);
	Bstring pre(char c);
	Bstring post(char c);
	Bstring pre_rev(char c);
	Bstring post_rev(char c);
	Bstring within(char c1, char c2);
	Bstring replace(char cold, char cnew);
	Bstring remove(const char c);
	Bstring swap(const long i, const long j);
	Bstring erase(const long i);
	Bstring truncate(const long n);
	Bstring base();
	Bstring extension();
	Bstring canonical(int n);
	long	integer();
	double	real();
	long	index(const char c) { return index(c, 0); }
	long	index(const char c, long start);
	long	find(const Bstring& s) { return find(s, 0); }
	long	find(const Bstring& s, long start);
	long	rfind(const Bstring& s) { return rfind(s, 0); }
	long	rfind(const Bstring& s, long start);
	bool	contains(const Bstring& s);
	Bstring insert(long pos, const Bstring& s);
	Bstring* split();
	Bstring* split(const Bstring& delim);
	vector<long>	split_into_integers(const Bstring& delim);
	vector<double>	split_into_doubles(const Bstring& delim);
	friend	Bstring operator+(const char* c, const Bstring& s);
};

//template Bstring::Bstring(long, const char*);

#endif

istream& operator>>(istream& input, Bstring& s);
ostream& operator<<(ostream& output, Bstring& s);

// Function prototypes 
Bstring*	string_add(Bstring** list, const char* string);
Bstring*	string_add(Bstring** list, Bstring& string);
Bstring		string_catenate(Bstring* list);
int			string_kill(Bstring* list);
Bstring		parameter_file_path(Bstring filename);
int			string_sort(Bstring* slist, int descending, int value);

