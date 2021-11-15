/**
@file	Bstring.cpp
@brief	Functions in the string class 
@author Bernard Heymann 
@date	Created: 20051020
@date	Modified: 20210115
**/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "Bstring.h"
#include "utilities.h"

// Definition of the global variables 
extern int 	verbose;		// Level of output to the screen

/**
	Constructor
*/

Bstring::Bstring(const Bstring& s)
{
	data = NULL;
	long	len = s.length();
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, s.data, len);
		data[len] = 0;
	}
	next = NULL;
}

Bstring::Bstring(const string& s)
{
	data = NULL;
	long	len = s.length();
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, s.data(), len);
		data[len] = 0;
	}
	next = NULL;
}

Bstring::Bstring(const Bstring& s, long start, long len)
{
	data = NULL;
	if ( start < 0 ) start = 0;
	long	end = start + len - 1;
	long	inlen = s.length();
	if ( end >= inlen ) end = inlen - 1;
	if ( start > end ) start = end;
	len = end - start + 1;
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, &s.data[start], len);
		data[len] = 0;
	}
	next = NULL;
}

Bstring::Bstring(const char* c)
{
	data = NULL;
	if ( c ) {
		long	len = ::strlen(c);
		if ( len ) {
			data = new char[len + 1];
			::strncpy(data, c, len);
			data[len] = 0;
		}
	}
	next = NULL;
}

Bstring::Bstring(const char* c, long start, long len)
{
	data = NULL;
	if ( c ) {
		if ( start < 0 ) start = 0;
		long	end = start + len - 1;
		long	inlen = ::strlen(c);
		if ( end >= inlen ) end = inlen - 1;
		if ( start > end ) start = end;
		len = end - start + 1;
		if ( len ) {
			data = new char[len + 1];
			::strncpy(data, c+start, len);
			data[len] = 0;
		}
	}
	next = NULL;
}

Bstring::Bstring(const char c, long size)
{
	data = NULL;
	if ( size ) {
		data = new char[size + 1];
		for ( long i=0; i<size; i++ ) data[i] = c;
		data[size] = 0;
	}
	next = NULL;
}

Bstring& Bstring::operator=(const Bstring& s)
{
	if ( this == &s ) return *this;
	if ( data ) delete[] data;	// Delete existing string
	data = NULL;
	long	len = s.length();
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, s.data, len);			// Copy new string
		data[len] = 0;
	}
	return *this;
}

Bstring& Bstring::operator=(const string& s)
{
	if ( data ) delete[] data;	// Delete existing string
	data = NULL;
	long	len = s.length();
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, s.c_str(), len);			// Copy new string
		data[len] = 0;
	}
	return *this;
}

Bstring& Bstring::operator=(const char* c)
{
	if ( data ) delete[] data;	// Delete existing string
	data = NULL;
	long	len(0);
	if ( c ) len = ::strlen(c);
	if ( len ) {
		data = new char[len + 1];
		::strncpy(data, c, len);
		data[len] = 0;
	}
	return *this;
}

Bstring& Bstring::operator+=(const Bstring& s)
{
	*this = *this + s;
	return *this;
}

Bstring Bstring::operator+(const Bstring& s)
{
	long	l1 = length();
	long	l2 = s.length();
	long	len = l1 + l2;
	Bstring	ns('\0', len);
	char*	ptr = ns.data;
	if ( data ) ptr = ::stpncpy(ns.data, data, l1);
	if ( s.data ) ::strncpy(ptr, s.data, l2);
	ns.data[len] = 0;
	return ns;
}

Bstring Bstring::operator+(const Bstring& s) const
{
	long	l1 = length();
	long	l2 = s.length();
	long	len = l1 + l2;
	Bstring	ns('\0', len);
	char*	ptr = ns.data;
	if ( data ) ptr = ::stpncpy(ns.data, data, l1);
	if ( s.data ) ::strncpy(ptr, s.data, l2);
	ns.data[len] = 0;
	return ns;
}

Bstring& Bstring::operator+=(const char c)
{
	*this = *this + c;
	return *this;
}

Bstring Bstring::operator+(const char c)
{
	long	len = length();
	long	nlen = len + 1;
	Bstring	ns('\0', nlen);
	if ( data ) ::strcpy(ns.data, data);
	ns.data[len] = c;
	return ns;
}

bool Bstring::operator==(const Bstring& s)
{
	return ( compare(s) == 0 );
}

bool Bstring::operator!=(const Bstring& s)
{
	return ( compare(s) != 0 );
}

bool Bstring::operator>(const Bstring& s)
{
	return ( compare(s) > 0 );
}

bool Bstring::operator<(const Bstring& s)
{
	return ( compare(s) < 0 );
}

bool Bstring::operator>=(const Bstring& s)
{
	return ( compare(s) >= 0 );
}

bool Bstring::operator<=(const Bstring& s)
{
	return ( compare(s) <= 0 );
}

Bstring Bstring::operator<<(int n)
{
	Bstring		ns(*this, n, length() - n);
	ns += substr(0, n);
	return ns;
}

Bstring Bstring::operator>>(int n)
{
	Bstring		ns(*this, length() - n, n);
	ns += substr(0, length() - n);
	return ns;
}

char& Bstring::operator[](long i)
{
	if ( !data ) {
		cerr << "Error in Bstring::operator[]: String length = 0!" << endl;
		bexit(-1);
	}
	long		len = length();
	if ( i < 0 ) {
		if ( len ) while ( i < 0 ) i += len;
		else i = 0;
	} else {
		while ( i >= len ) i -= len;
	}
	return data[i];
}

long Bstring::length() const
{
	if ( !data ) return 0;
	return ::strlen(data);
}

long Bstring::count(const char c)
{
	long		n(0);
	if ( data )
		for ( char* s = data; *s; s++ ) if ( *s == c ) n++;
	return n;
}

long Bstring::compare(const Bstring& s)
{
	if ( !data && !s.data ) return 0;
	if ( !data ) return -1;
	if ( !s.data ) return 1;
	return ::strcmp(data, s.data);
}

long Bstring::compare_value(const Bstring& s)
{
	if ( !data || !s.data ) {
		cerr << "Error in Bstring::compare_value! (data = " << &data << " s.data = " << &(s.data) << endl;
		bexit(-1);
	}
	long		ld = length() - s.length();
	if ( ld ) return ld;
	return ::strcmp(data, s.data);
}

Bstring	Bstring::common(const Bstring& s)
{
	long		n(length());
	if ( data && s.data )
		for ( n=0; data[n] && s.data[n] && data[n] == s.data[n]; ++n ) ;
	return substr(0, n);
}

void Bstring::fill(const char c)
{
	if ( data ) 
		for ( char* s = data; *s; s++ ) *s = c;
}

void Bstring::place(long i, const Bstring& s)
{
	long		len(s.length());
	if ( data ) strncpy(&data[i], s.data, len);
}

Bstring Bstring::lower()
{
	Bstring		ns(*this);
	char		*s, *sn;
	for ( s = data, sn = ns.data; s && *s; s++, sn++ ) *sn = tolower(*s);
	return ns;
}

Bstring Bstring::upper()
{
	Bstring		ns(*this);
	char		*s, *sn;
	for ( s = data, sn = ns.data; *s; s++, sn++ ) *sn = toupper(*s);
	return ns;
}

Bstring Bstring::substr(long start, long len)
{
	long			copylen = length() - start;
	if ( copylen > len ) copylen = len;
	Bstring	ns(*this, start, copylen);
	return ns;
}

Bstring Bstring::no_space()
{
	long	i, j;
	for ( i=j=0; data[i]; i++ ) if ( !isspace(data[i]) ) j++;
	Bstring			ns('\0', j);
	for ( i=j=0; data[i]; i++ ) if ( !isspace(data[i]) ) ns.data[j++] = data[i];
	return ns;
}

Bstring Bstring::no_lead_space()
{
	long	i, j;
	for ( i=j=0; isspace(data[i]); i++ ) j++;
	Bstring	ns(*this, j, length());
	return ns;
}

Bstring Bstring::alnum()
{
	long	i, j;
	for ( i=j=0; data[i]; i++ ) if ( isalnum(data[i]) ) j++;
	Bstring			ns('\0', j);
	for ( i=j=0; data[i]; i++ ) if ( isalnum(data[i]) ) ns.data[j++] = data[i];
	return ns;
}

Bstring Bstring::left(long len)
{
	long			copylen = length();
	if ( copylen > len ) copylen = len;
	Bstring			ns(' ', copylen);
	strncpy(ns.data, data, copylen);
	return ns;
}

Bstring Bstring::right(long len)
{
	long			copylen = length();
	if ( copylen > len ) copylen = len;
	long			start = len - copylen;
	Bstring			ns(' ', len);
	strncpy(ns.data + start, data, copylen);
	return ns;
}

/*	Return the string before the first character.
	If the character is not in the string, return the whole string.
	If the character is the first in the string, return an empty string.
*/
Bstring Bstring::pre(char c)
{
	Bstring		ns;
	long		i, len = length();
	for ( i=0; i<len && data[i] != c; i++ ) ;
	if ( i > 0 ) ns = Bstring(*this, 0, i);
	return ns;
}

/*	Return the string after the first character.
	If the character is not in the string, return an empty string.
	If the character is the last in the string, return an empty string.
*/
Bstring Bstring::post(char c)
{
	Bstring		ns;
	long		i, len = length(), j = len;
	for ( i=0; i<len && data[i] != c; i++ ) ;
	if ( i < len ) {
		i++;
		j = len - i;
		ns = Bstring(*this, i, j);
	}
//	} else i = 0;
//	if ( j > 0 ) ns = Bstring(*this, i, j);
	return ns;
}

/*	Return the string before the last character.
	If the character is not in the string, return the whole string.
	If the character is the first in the string, return an empty string.
*/
Bstring Bstring::pre_rev(char c)
{
	Bstring		ns;
	long		i, j = -1;
	char*		s;
	for ( i=0, s = data; *s; i++, s++ ) if ( *s == c ) j = i;
	if ( j < 0 ) j = i;
	if ( j > 0 ) ns = Bstring(*this, 0, j);
	return ns;
}

/*	Return the string after the last character.
	If the character is not in the string, return an empty string.
	If the character is the last in the string, return an empty string.
*/
Bstring Bstring::post_rev(char c)
{
	Bstring		ns;
	long		i, j = -1;
	char*		s;
	for ( i=0, s = data; *s; i++, s++ ) if ( *s == c ) j = i;
	if ( j > -1 ) {
		j++;
		i -= j;
		if ( i > 0 ) ns = Bstring(*this, j, i);
	}
	return ns;
}

Bstring Bstring::within(char c1, char c2)
{
	Bstring		ns;
	long		i, j;
	char*		s;
	for ( i=0, s = data; *s && *s != c1; i++, s++ ) ;
	i++;
	if ( *s ) s++;
	if ( *s ) {
		for ( j=0; *s && *s != c2; j++, s++ ) ;
		if ( *s && j > 0 ) ns = Bstring(*this, i, j);
	}
	return ns;
}

Bstring Bstring::replace(char cold, char cnew)
{
	Bstring		ns(*this);
	for ( char* s = ns.data; *s; s++ )
		if ( *s == cold ) *s = cnew;
	return ns;
}

Bstring Bstring::remove(const char c)
{
	long		nc = count(c);
	if ( nc < 1 ) return *this;
	long		n = length() - nc;
	Bstring		ns(' ', n);
	char		*s, *sn;
	for ( s = data, sn = ns.data; *s; s++ )
		if ( *s != c ) *(sn++) = *s;
	return ns;
}

Bstring Bstring::swap(const long i, const long j)
{
	Bstring		ns(*this);
	ns.data[i] = data[j];
	ns.data[j] = data[i];
	return ns;
}

Bstring Bstring::erase(const long i)
{
	long		j, n = length() - 1;
	Bstring		ns(' ', n);
	char		*s, *sn;
	for ( j=0, s = data, sn = ns.data; *s; s++, j++ )
		if ( j != i ) *(sn++) = *s;
	return ns;
}

Bstring Bstring::truncate(const long n)
{
	long		nlen = length() - n;
	Bstring		ns(' ', nlen);
	char		*s, *sn;
	long		i;
	for ( i=0, s = data, sn = ns.data; i<nlen; i++ )
		*(sn++) = *(s++);
	return ns;
}

/*	Removes both the leading path and extension from a filename string.
*/
Bstring Bstring::base()
{
	Bstring		ns(*this);
	if ( ns.contains("/") ) ns = ns.post_rev('/');
	return ns.pre_rev('.');	
}

Bstring Bstring::extension()
{
	Bstring		ns = pre(',');
	if ( ns.find("#", 0) > -1 )
		ns = "raw";
	else if ( ns.find(":", 0) > -1 )
		ns = ns.post(':');
	else
		ns = ns.post_rev('.');
	return ns.lower();	
}

Bstring Bstring::canonical(int n)
{
	int			i;
	Bstring		s(*this);
	Bstring		cs(*this);
	for ( i=0; i<s.length(); i+=n ) {
		s = s << n;
		if ( cs > s ) cs = s;
	}	
	return cs;
}

long	Bstring::integer()
{
	if ( !data ) return 0;
	return strtol(data, NULL, 10);
}

double	Bstring::real()
{
	if ( !data ) return 0;
	return strtod(data, NULL);
}

/* The offset of a search character is returned if found, otherwise a negative value is returned */
long	Bstring::index(const char c, long start)
{
	long		offset;
	char*		aptr = &data[start];
	for ( offset = start; *aptr && *aptr != c; aptr++, offset++ ) ;
	if ( *aptr != c ) offset = -1;
	return offset;
}

/* The offset of a search string is returned if found, otherwise a negative value is returned */
long	Bstring::find(const Bstring& s, long start)
{
	long		len = s.length();
	long		offset;
	char*		aptr = &data[start];
	for ( offset = start; *aptr; aptr++, offset++ )
		if ( ::strncmp(s.data, aptr, len) == 0 ) break;
	if ( *aptr == 0 ) offset = -1;
	return offset;
}

/* The offset of a search string is returned if found, otherwise a negative value is returned */
long	Bstring::rfind(const Bstring& s, long start)
{
	long		i, len = s.length();
	long		offset = -1;
	char*		aptr = &data[start];
	for ( i = start; *aptr; aptr++, i++ )
		if ( ::strncmp(s.data, aptr, len) == 0 ) offset = i;
	return offset;
}

bool	Bstring::contains(const Bstring& s)
{
	if ( !data ) return 0;
	return strstr(data, s.data) != NULL;
}

Bstring Bstring::insert(long pos, const Bstring& s)
{
	Bstring		ns = substr(0, pos);
	Bstring		ne = substr(pos, length() - pos);
	return ns + s + ne;
}

Bstring* Bstring::split()
{
	if ( !data ) return NULL;
	long		len = length();
	if ( len < 1 ) return NULL;
	long		i;
	Bstring		tstr;
	Bstring*	list = NULL;
	Bstring*	astring = NULL;
	char*		sptr = data;
	char*		aptr = data;
	for ( i=0; i<=len; i++, aptr++ ) {
		if ( isspace(*sptr) ) sptr++;
		else if ( *aptr==0 || isspace(*aptr) ) {
			tstr = substr(sptr - data, aptr - sptr);
			if ( tstr.length() ) {
				astring = string_add(&astring, tstr);
				if ( !list ) list = astring;
				sptr = aptr;
			}
		}
	}
	return list;
}

Bstring* Bstring::split(const Bstring& delim)
{
	long		dlen = delim.length();
	if ( dlen < 1 ) return split();
	long		len = length();
	if ( len < 1 ) return NULL;
	long		i;
	Bstring		tstr;
	Bstring*	list = NULL;
	Bstring*	astring = NULL;
	char*		sptr = data;
	char*		aptr = data;
	for ( i=0; i<=len; i++, aptr++ ) {
		if ( *aptr==0 || ::strncmp(aptr, delim.data, dlen) == 0 ) {
			tstr = substr(sptr - data, aptr - sptr);
			if ( tstr.length() ) {
				astring = string_add(&astring, tstr);
				if ( !list ) list = astring;
				sptr = aptr + dlen;
				aptr += dlen - 1;
				i += dlen - 1;
			}
		}
	}
	return list;
}

vector<long>	Bstring::split_into_integers(const Bstring& delim)
{
	vector<long>	array;
	
	Bstring*		list = split(delim);
	if ( !list ) return array;
	
	Bstring*		astring = NULL;
	for ( astring=list; astring; astring=astring->next )
		array.push_back(astring->integer());
	
	string_kill(list);
	
	return array;
}

vector<double>	Bstring::split_into_doubles(const Bstring& delim)
{
	vector<double>	array;
	
	Bstring*		list = split(delim);
	if ( !list ) return array;
	
	Bstring*		astring = NULL;
	for ( astring=list; astring; astring=astring->next )
		array.push_back(astring->real());

	string_kill(list);

	return array;
}

Bstring operator+(const char* c, const Bstring& s)
{
	Bstring ns(c);
	return ns += s;
}

istream& operator>>(istream& input, Bstring& s) {
	string		str;
	input >> str;
	s = str;
	return input;
}

ostream& operator<<(ostream& output, Bstring& s) {
	if ( s.length() ) output << left << s.c_str();
	return output;
}

/**
@brief 	Adds a string to a linked list.
@param 	**list 		the string linked list.
@param 	&string		string - unchanged.
@return Bstring*			new string structure.
**/
Bstring*	string_add(Bstring** list, const char* string)
{
	Bstring 		this_string(string);
	return string_add(list, this_string);
}

Bstring*	string_add(Bstring** list, Bstring& string)
{
	Bstring* 		this_string = *list;
	
	if ( string.empty() ) return this_string;
	
	Bstring* 		new_string = new Bstring;
	
	*new_string = string;
	
	if ( !this_string )
		*list = new_string;
	else {
		while ( this_string->next ) this_string = this_string->next;
		this_string->next = new_string;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG string_add: " << *new_string << " (" << new_string->length() << ")" << endl;
	
	return new_string;
}

/**
@brief 	Catenates a string linked list into one string.
@param 	*list 		the string linked list.
@return Bstring		new string.
**/
Bstring		string_catenate(Bstring* list)
{
	long		i, len(0);
	Bstring*	str;
	
	for ( str = list; str; str = str->next )
		len += str->length();
	
	Bstring		cat = Bstring(' ', len);

	for ( i=0, str = list; str && i<len; i+=str->length(), str = str->next )
		cat.place(i, *str);
	
	return cat;
}

/**
@brief 	Kills a string linked list.
@param 	*list 		the string linked list.
@return int			0.
**/
int			string_kill(Bstring* list)
{
	if ( !list ) return 0;
	
	Bstring		*s, *s2;
	
	for ( s = list; s; ) {
		s2 = s->next;
//		*s = 0;
		delete s;
		s = s2;
	}
	
	return 0;
}

/**
@brief 	Finds the parameter file path.
@param 	&filename	the parameter filename.
@return Bstring		the full path and file name.

	The parameter file path should primarily be defined in the environmental
	variable "BPARAM".
	Otherwise a default is used which may or may not be valid.
	The string returned has a length containing the path plus one byte.

**/
Bstring		parameter_file_path(Bstring filename)
{
	Bstring		path;

	if ( access(filename.c_str(), R_OK) == 0 )
		return filename;
	
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

/**
@brief 	Utility function for sorting string lexically in qsort.
@param 	*x		first string pointer.
@param 	*y		second string pointer.
@return int		-1 if x < y and 1 otherwise.
**/
int			QsortSmallToLargeString(const void *x, const void *y)
{
	Bstring* 	s1 = (Bstring*) x;
	Bstring* 	s2 = (Bstring*) y;
	if( *s1 < *s2 ) return -1;
	else return 1;
}

/**
@brief	Utility function for sorting string lexically in qsort.
@param	*x		first string pointer.
@param	*y		second string pointer.
@return				-1 if x < y and 1 otherwise.
**/
int			QsortLargeToSmallString(const void *x, const void *y)
{
	Bstring* 	s1 = (Bstring*) x;
	Bstring* 	s2 = (Bstring*) y;
	if( *s1 > *s2 ) return -1;
	else return 1;
}

/**
@brief	Utility function for sorting string values in qsort.
@param	*x		first string pointer.
@param	*y		second string pointer.
@return				-1 if x < y and 1 otherwise.
**/
int			QsortSmallToLargeStringValue(const void *x, const void *y)
{
	Bstring* 	s1 = (Bstring*) x;
	Bstring* 	s2 = (Bstring*) y;
	if( s1->compare_value(*s2) < 0 ) return -1;
	else return 1;
}

/**
@brief	Utility function for sorting string values in qsort.
@param	*x		first string pointer.
@param	*y		second string pointer.
@return				-1 if x < y and 1 otherwise.
**/
int			QsortLargeToSmallStringValue(const void *x, const void *y)
{
	Bstring* 	s1 = (Bstring*) x;
	Bstring* 	s2 = (Bstring*) y;
	if( s1->compare_value(*s2) > 0 ) return -1;
	else return 1;
}

/**
@brief 	Sorts a list of strings in ascending (descending) lexical or value order.
@param 	*slist		string list.
@param 	descending		flag to do a descending sort.
@param 	value			sort in value order.
@return int					0.

	The string list is first convert to an array, quicksorted, and
	written back into the string list.

**/
int			string_sort(Bstring* slist, int descending, int value)
{
	int			i, n(0);
	Bstring*	s;
	
	for ( s = slist; s; s = s->next ) n++;
	
	Bstring*	ssort = new Bstring[n];

	for ( i=0, s = slist; s; s = s->next, i++ ) ssort[i] = *s;
	
	if ( value ) {
		if ( descending ) qsort((void *) ssort, n, sizeof(Bstring), QsortLargeToSmallStringValue);
		else qsort((void *) ssort, n, sizeof(Bstring), QsortSmallToLargeStringValue);
	} else {
		if ( descending ) qsort((void *) ssort, n, sizeof(Bstring), QsortLargeToSmallString);
		else qsort((void *) ssort, n, sizeof(Bstring), QsortSmallToLargeString);
	}
	
	for ( i=0, s = slist; s; s = s->next, i++ ) *s = ssort[i];
	
	delete[] ssort;

	return 0;
}
