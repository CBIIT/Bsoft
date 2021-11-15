/**
@file	bfile.cpp
@brief	Program to read file contents directly and poke single values
@author Bernard Heymann
@date	Created: 19980129
@date	Modified: 20210430
**/

#include "options.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
int 		show_block(size_t start, size_t size, unsigned char* buf);
int 		show_block_hex(size_t start, size_t size, unsigned char* buf);
int 		show_char(size_t start, size_t size, unsigned char* buf);
int 		swap_block(size_t size, size_t nbytes, unsigned char* buf);
int 		show_table(size_t start, size_t size, unsigned char* buf);
int 		show_structure(size_t start, size_t size, int swap, char* strdef, unsigned char* buf);
int 		show_eer(size_t start, size_t size, unsigned char* buf);
int 		poke_value(fstream* fin, size_t loc, char typechar, double value);
int			find_string(const char* string, size_t start, size_t size, unsigned char* buf);
int 		replace(fstream* fin, char* rstring1, char* rstring2);
char* 		string2byte(char* string);
size_t	string_bytes(char* string);

const char* use[] = {
" ",
"Usage: bfile [options] infile [outfile]",
"---------------------------------------",
"Examines and modifies binary files.",
" ",
"Actions:",
"-blocks            Character value blocks (16 per line).",
"-hexadecimal       Character hexadecimal value blocks (16 per line).",
"-characters        Character blocks (64 per line).",
"-format bbsif      Output formatted according to a sequence of data types.",
"-find astring      Find the given string.",
"-poke 3124,20.4,f  Poking at a location a value of a specific data type.",
"-replace \"\\n,\\r\"   Replace all occurrences of the first string by the second.",
"-eer               Show EER decompression details.",
" ",
"Parameters:",
"-verbose 1         Verbosity of output.",
"-bytes 128         Number of bytes to read (default 64).",
"-start 646         Starting position to read from (default 0).",
"-reverse 2         Reverse byte order for (2) bytes.",
"-VAX               Convert VAX floating point numbers.",
" ",
"Data types:",
"    b: byte or unsigned char (1 byte)",
"    c: signed char (1 byte)",
"    u: word or unsigned short (2 bytes)",
"    s: short (2 bytes)",
"    j: unsigned integer (4 bytes)",
"    i: integer (4 bytes)",
"    f: float (4 bytes)",
"    d: double (8 bytes)",
" ",
NULL
};

int			main(int argc, char **argv)
{
	// Initialize variables
    size_t			buflen(64);   		/* Default buffer length	*/
    size_t			start(0);	    	/* Default start of buffer	*/
	size_t 			nbytes(0); 			/* Number of bytes to swap */
    int     	    setvax(0);
    int     	    setblock(0);
    int     	    sethex(0);
    int     	    setchar(0);
    int				seteer(0);		// Show EER details
    char*     	    strdef = new char[128];
	Bstring			find;				/* String to find */
	int 			poke(0);			/* Flag for poking values */
	size_t			loc(0);				/* Location for poking values */
	char 			typechar = 'u';		/* Default uchar for poking */
	double			pokevalue(0);		/* Value for poking */
	char			rstring1[200] = ""; // String to replace
	char			rstring2[200] = ""; // New string
	size_t			i;
      
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "bytes" ) buflen = curropt->value.integer();
		if ( curropt->tag == "start" ) start = curropt->value.integer();
		if ( curropt->tag == "reverse" ) nbytes = curropt->value.integer();
		if ( curropt->tag == "VAX" ) setvax = 1;
		if ( curropt->tag == "blocks" ) setblock = 1;
		if ( curropt->tag == "hexadecimal" ) sethex = 1;
		if ( curropt->tag == "characters" ) setchar = 1;
		if ( curropt->tag == "format" ) sscanf(curropt->value.c_str(), "%s", strdef);
		if ( curropt->tag == "find" ) find = curropt->value;
		if ( curropt->tag == "eer" ) seteer = 1;
        if ( curropt->tag == "poke" ) sscanf(curropt->value.c_str(), "%ld,%lf,%c", &loc, &pokevalue, &typechar);
		if ( curropt->tag == "replace" ) {
			memcpy(rstring1, &curropt->value[0], curropt->value.length()+1);
			if ( strlen(rstring1) < 1 )
				cerr << "-replace: Two strings must be specified!" << endl;
			else {
				i = strcspn(rstring1, ",");
				strcpy(rstring2, rstring1+i+1);
				rstring1[i] = 0;
			}
		}
	}
	option_kill(option);
        
	fstream*		fin = new fstream(argv[optind], fstream::in | fstream::out);
	if ( fin->fail() ) {
		cerr << "Error: File " << argv[optind] << " not read!" << endl << endl;
		bexit(-1);
	} else optind++;
    
	if ( strlen(rstring1) > 0 ) {
		replace(fin, rstring1, rstring2);
		bexit(0);
	}
	
	if ( poke ) {
		poke_value(fin, loc, typechar, pokevalue);
		start = loc;
		buflen = 8;
	}
	
	fin->seekg (0, ios::end);
	size_t	fsize = (size_t)fin->tellg() - start;
	if ( buflen > fsize ) buflen = fsize;
	
    size_t	size = buflen*sizeof(unsigned char);
    unsigned char*	buf = new unsigned char[size];

	fin->seekg (start, ios::beg);
	fin->read((char *)buf, size);
	fin->close();
	delete fin;
    
	if ( setblock ) {
		show_block(start, size, buf);
		bexit(0);
	}
	
	if ( sethex ) {
		show_block_hex(start, size, buf);
		bexit(0);
	}
	
	if ( setchar ) {
		show_char(start, size, buf);
		bexit(0);
	}
	
	if ( strlen(strdef) ) {
		show_structure(start, size, nbytes, strdef, buf);
		bexit(0);
	}
	
	if ( find.length() ) {
		find_string(find.c_str(), start, size, buf);
//		delete[] find;
		find = 0;
		bexit(0);
	}

	if ( seteer ) {
		show_eer(start, size, buf);
		bexit(0);
	}

    if ( setvax ) {
	    cout << "Getting VAX floating point" << endl;
		for ( i=0; i<size; i+=4 ) vax2ieee((unsigned char *)(buf+i), 0);
	}
	
    if ( nbytes ) swap_block(size, nbytes, buf);
    
	if ( verbose ) show_table(start, size, buf);
	
	if ( optind < argc ) {
		ofstream	fout(argv[optind]);
		if ( fout.fail() ) {
			cerr << "Error: Unable to open " << argv[optind] << endl;
			bexit(-1);
		}
		fout.write((char *)buf, size);
		fout.close();
	}
	
	delete[] buf;
    
	bexit(0);
}

int 		show_block(size_t start, size_t size, unsigned char* buf)
{
	size_t 			i, j;
	
	for ( i=0; i<size; i+=16 ) {
		cout << setw(10) << i+start << ": ";
		for ( j=i; j<i+16; j++ ) cout << setw(4) << (int) buf[j];
		cout << "  |";
		for ( j=i; j<i+16; j++ ) {
			if ( isprint(buf[j]) ) cout << buf[j];
			else cout << " ";
		}
		cout << "|" << endl;
	}
	
	return 0;
}

int 		show_block_hex(size_t start, size_t size, unsigned char* buf)
{
	size_t 			i, j;
	
	for ( i=0; i<size; i+=16 ) {
		cout << setw(10) << i+start << ": ";
		for ( j=i; j<i+16; j++ ) cout << hex << setw(3) << (int) buf[j];
		cout << "  |";
		for ( j=i; j<i+16; j++ ) {
			if ( isprint(buf[j]) ) cout << buf[j];
			else cout << " ";
		}
		cout << "|" << endl;
	}
	
	return 0;
}

int 		show_char(size_t start, size_t size, unsigned char* buf)
{
	size_t 			i, j;
	
	for ( i=0; i<size; i+=64 ) {
		cout << setw(10) << i+start << ": ";
		for ( j=i; j<i+64; j++ ) {
			if ( isprint(buf[j]) ) cout << buf[j];
			else cout << " ";
		}
		cout << endl;
	}
	
	return 0;
}

int 		swap_block(size_t size, size_t nbytes, unsigned char* buf)
{
	size_t 					i, j;
	unsigned char			t[8];
	
    cout << "Swapping byte order for " << nbytes << "-byte types" << endl;
	for ( i=0; i<size; i+=nbytes ) {
		memcpy(t, buf+i, nbytes);
		for( j=0; j<nbytes; j++ ) buf[i+j] = t[nbytes-1-j];
	}
	
	return 0;
}

int 		show_table(size_t start, size_t size, unsigned char* buf)
{
    unsigned char   bits[9], thechar;
	size_t 	i, j;
    short* 			inshort = (short *) buf;
    int*			inint = (int *) buf;
    long*			inlong = (long *) buf;
    float*			infloat = (float *) buf;
    double*			indouble = (double *) buf;

    cout << "  loc      char   bits    char  short     int       long       float           double" << endl;
    for ( i=0; i<size; i++ ) {
    	for ( j=0; j<8; j++) {
    	    bits[j] = '0';
    	    if ( ( buf[i] << j ) & 0x80 )
    	    	bits[j] = '1';
    	}
		bits[8] = '\0';
    	thechar = (unsigned char)buf[i];
    	if ( !isprint(thechar) ) thechar = ' ';
    	cout << setw(10) << i+start << setw(4) << thechar << setw(10) << bits << setw(5) << (int)buf[i];
    	if ( i%2 == 0 ) cout << setw(8) << inshort[i/2];
    	if ( i%4 == 0 ) cout << setw(10) << inint[i/4];
    	if ( i%4 == 0 ) cout << setw(24) << infloat[i/4];
    	if ( i%8 == 1 ) cout << setw(36) << inlong[i/8];
    	if ( i%8 == 1 ) cout << setw(24) << indouble[i/8];
    	cout << endl;
    }
	
	return 0;
}

int 		show_structure(size_t start, size_t size, int swap, char* strdef, unsigned char* buf)
{
	size_t 			i, j;
	unsigned char*			uchar;
	signed char*			schar;
	unsigned short*			ushort;
	signed short*			sshort;
	int*					integer;
	float*					floating;
	double*					doublefloat;
	
    cout << "loc\ttype\tvalue" << endl;
	i = j = 0;
	while ( i<size ) {
		switch ( strdef[j] ) {
			case 'b':
				uchar = (unsigned char *) (buf + i);
				cout << i+start << tab << strdef[j] << tab << *uchar << endl;
				i += sizeof(char);
				break;
			case 'c':
				schar = (signed char *) (buf + i);
				cout << i+start << tab << strdef[j] << tab << *schar << endl;
				i += sizeof(char);
				break;
			case 'u':
				ushort = (unsigned short *) (buf + i);
				if ( swap ) swapbytes((unsigned char *) ushort, sizeof(short));
				cout << i+start << tab << strdef[j] << tab << *ushort << endl;
				i += sizeof(short);
				break;
			case 's':
				sshort = (short *) (buf + i);
				if ( swap ) swapbytes((unsigned char *) sshort, sizeof(short));
				cout << i+start << tab << strdef[j] << tab << *sshort << endl;
				i += sizeof(short);
				break;
			case 'i':
				integer = (int *) (buf + i);
				if ( swap ) swapbytes((unsigned char *) integer, sizeof(int));
				cout << i+start << tab << strdef[j] << tab << *integer << endl;
				i += sizeof(int);
				break;
			case 'f':
				floating = (float *) (buf + i);
				if ( swap ) swapbytes((unsigned char *) floating, sizeof(float));
				cout << i+start << tab << strdef[j] << tab << *floating << endl;
				i += sizeof(float);
				break;
			case 'd':
				doublefloat = (double *) (buf + i);
				if ( swap ) swapbytes((unsigned char *) doublefloat, sizeof(double));
				cout << i+start << tab << strdef[j] << tab << *doublefloat << endl;
				i += sizeof(double);
				break;
			default: break;
		}
		j++;
		if ( j > strlen(strdef) ) j = 0;
	}
	
	return 0;
}

int			show_bits(unsigned long val)
{
	bitset<64> b(val);
	cout << val << tab << b << endl;
	
	return 0;
}

int			show_code(unsigned long code)
{
	bitset<11> b(code);
	cout << code << tab << (code & 127) << tab << b << endl;
	
	return 0;
}

int			shift_right(unsigned long& bits, long n)
{
	if ( bits < ULLONG_MAX ) {
		bits >>= n;
		return 0;
	}
	
	bits /= 2;
	bits >>= (n-1);

	return 0;
}

int 		show_eer(size_t start, size_t size, unsigned char* buf)
{
	long		pxbits(7);					// Value size for a pixel
	long		subpx(2);					// Pixel division depth
	long		codebits(pxbits+2*subpx);	// Full size of a code element
	long		maxcode((1<<codebits)-1);	// Maximum code element value covering all bits
	long		maxval((1<<pxbits)-1);		// Maximum number of pixels in a code element
	long		n(0), xbits, rbits(64), npx(0), ne(0);
	unsigned long	code(0);
	unsigned long	rawcc(size);
	unsigned long*	rawcp = (unsigned long *)buf;
	unsigned long	bits = *rawcp;	// Read the first word
	rawcp++;	// Adjust the input pointer and counter
	rawcc -= sizeof(long);
//	show_bits(bits);
//	cout << maxcode << tab << maxval << endl;
	
	while ( rawcc > 0 ) {
		code = bits & maxcode;
		show_bits(bits);
		show_code(code);
//		bits >>= codebits;
		shift_right(bits, codebits);
		rbits -= codebits;
		if ( rbits <= 0 ) {
			bits = *rawcp;	// Read the next word
			rawcp++;		// Adjust the input pointer and counter
			rawcc -= sizeof(long);
			cout << "rbits: " << rbits << endl;
			show_bits(bits);
			xbits = -rbits;
			if ( xbits > 0 && xbits <= codebits ) {	// Add the missing part
				code |= ((bits & ((1<<xbits)-1)) << (codebits-xbits));
//				bits >>= xbits;
				shift_right(bits, xbits);
			}
			show_bits(bits);
			show_code(code);
			rbits += 64;
		}
//		show_code(code);
		npx = code & maxval;
		n += npx;
		buf += npx;
		if ( npx < maxval ) {		// Count an electron
			*buf = 1;
			buf++;
			n++;
			ne++;
		}
//		cout << npx << tab << n << tab << rbits << tab << ne << endl;
	}

	return 0;
}

int			find_string(const char* string, size_t start, size_t size, unsigned char* buf)
{
	size_t 			i, j, n(0), ssize = strlen(string);
	
	cout << "Searching for: " << string << endl;
	for ( i=0; i<size - ssize + 1; i++, buf++ ) {
		if ( strncasecmp(string, (char *)buf, ssize) == 0 ) {
			cout << setw(10) << i+start << ":";
			for ( j=0; j<64; j++ ) cout << buf[j];
			cout << endl;
			n++;
		}
	}
	cout << "Found:                          " << n << endl;
	
	return 0;
}

int 		poke_value(fstream* fin, size_t loc, char typechar, double value)
{
	unsigned char	bval = (unsigned char)value;
	signed char		cval = (signed char)value;
	unsigned short	uval = (unsigned short)value;
	short			sval = (short)value;
	int 			ival = (int)value;
	double 			dval = (double)value;
	
	fin->seekg(loc, ios::beg);
	switch ( typechar ) {
		case 'b':
			fin->write((char *)&bval, sizeof(char));
			break;
		case 'c':
			fin->write((char *)&cval, sizeof(char));
			break;
		case 'u':
			fin->write((char *)&uval, sizeof(short));
			break;
		case 's':
			fin->write((char *)&sval, sizeof(short));
			break;
		case 'i':
			fin->write((char *)&ival, sizeof(int));
			break;
		case 'f':
			fin->write((char *)&value, sizeof(float));
			break;
		case 'd':
			fin->write((char *)&dval, sizeof(double));
			break;
		default: break;
	}
	
	return 0;
}

int 		replace(fstream* fin, char* rstring1, char* rstring2)
{
	size_t 		i, j, m, n, nr(0);
	
	char*				fromstr = string2byte(rstring1);
	char*				tostr = string2byte(rstring2);
	
	if ( verbose )
		cout << "Replacing \"" << rstring1 << "\" with \"" << rstring2 << "\"" << endl;
	
//	if ( strlen(tostr) < 1 ) strcpy(tostr, " ");
				
	if ( strlen(fromstr) < 1 ) {
		cerr << "Error: The search string has zero length!" << endl;
		return -1;
	}
	
	fin->seekg(0, ios::end);
	size_t 		filesize = fin->tellg();
	fin->seekg(0, ios::beg);
	
	if ( verbose )
		cout << "File size:                      " << filesize << endl;
	
	char*		buf = new char[filesize];
	char*		newbuf = new char[filesize];
	
	fin->read(buf, filesize);
	
	i = j = nr = 0;
	m = string_bytes(fromstr);
	n = string_bytes(tostr);
	while ( i < filesize ) {
		if ( memcmp(buf+i, fromstr, m) == 0 ) {
			memcpy(newbuf+j, tostr, n);
			i += m;
			j += n;
			nr++;
		} else {
			newbuf[j] = buf[i];
			i++;
			j++;
		}
	}
	
	fin->seekg(0, ios::beg);
	
	fin->write(newbuf, filesize);
	
	fin->close();
	
	if ( verbose )
		cout << "Number of replacements:         " << nr << endl << endl;
	
	return 0;
}

char* 		string2byte(char* string)
{
	size_t 		i, j, setcntl(0);
	char*				bytes = new char[strlen(string)+1];
	
	j = 0;
	for ( i=0; i<strlen(string); i++ ) {
		bytes[j] = 0;
		if ( setcntl ) {
			if ( string[i] == 'a' ) bytes[j] = 7;
			if ( string[i] == 'b' ) bytes[j] = 8;
			if ( string[i] == 't' ) bytes[j] = 9;
			if ( string[i] == 'n' ) bytes[j] = 10;
			if ( string[i] == 'v' ) bytes[j] = 11;
			if ( string[i] == 'f' ) bytes[j] = 12;
			if ( string[i] == 'r' ) bytes[j] = 13;
			if ( string[i] == '\\' ) bytes[j] = 28;
			setcntl = 0;
			j++;
		} else if ( string[i] == '\\' ) {
			setcntl = 1;
		} else {
			bytes[j] = string[i];
			j++;
		}
	}
	
	delete[] bytes;
	
	return bytes;
}

size_t	string_bytes(char* string)
{
	size_t 		i(0), n(0);
	
	for ( i=0; i<strlen(string); i++ ) if ( string[i] != '\\' ) n++;
	
	return n;
}

