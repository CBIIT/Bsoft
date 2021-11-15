/**
@file	bfix.cpp
@brief	Program to fix end-of-line characters in ascii files
@author Bernard Heymann
@date	Created: 19990206
@date	Modified: 20120306
**/

#include "utilities.h"
#include "options.h"
#include "timer.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototype
long		replace_a_character(unsigned char* buf, long size, int chin, int chout);
long		delete_a_character(unsigned char* buf, long size, int c);

const char* use[] = {
" ",
"Usage: bfix [options] file.in file.out",
"--------------------------------------",
"Replaces single characters in a file.",
"If only one file name is given, the output is sent to stdout.",
"If a \".\" character is used for output, the input file name is used",
"(without the path).",
"Note: replace and delete are mutually exclusive - the last option wins",
" ",
"Actions:",
"-n                 Replace <cr> with <nl> (To UNIX style)",
"-c                 Replace <nl> with <cr> (To MAC style)",
"-replace 9,32      Replace all character values 9 with character values 32",
"-delete 13         Delete all character values 13",
" ",
"Parameters:",
"-verbose 7         Verbosity of output",
" ",
NULL
};

int			main(int argc, char **argv)
{
	int				char_out(0), char_in(0);
	Bstring			filename;
    
	int 			replace_flag(0);		// Flag for character replacement
      
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "verbose" ) verbose = get_option_verbose(curropt->value);
		if ( curropt->tag == "c" ) { replace_flag = 1; char_in = 10; char_out = 13; }
		if ( curropt->tag == "n" ) { replace_flag = 2; char_in = 13; char_out = 10; }
		if ( curropt->tag == "replace" ) {
			if ( curropt->values(char_in, char_out) < 1 ) {
				cerr << "-replace: At least one character must be specified!" << endl;
				bexit(-1);
			} else
				replace_flag = 3;
        }
		if ( curropt->tag == "delete" ) {
			if ( ( char_in = curropt->value.integer() ) < 1 ) {
				cerr << "-delete: A character value must be specified!" << endl;
				bexit(-1);
			} else
				replace_flag = 4;
        }
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	filename = argv[optind++];
	
	ifstream		fin(filename.c_str());
	if ( fin.fail() ) {
		cerr << "Error: File " << filename << " not read!" << endl;
		bexit(-1);
	}
    
	// Get the file size to determine the space to allocate
	fin.seekg (0, ios::end);
	long			filesize = fin.tellg();
	fin.seekg (0, ios::beg);
		
	if ( verbose )
		cout << "Input file size:                " << filesize << endl;
	
    if ( verbose & VERB_LABEL ) {
		switch ( replace_flag ) {
    		case 1: cout << "Converting all NL in " << filename << " to CR" << endl; break;
			case 2: cout << "Converting all CR in " << filename << " to NL" << endl; break;
    		case 3: cout << "Replacing all " << char_in << " in " << filename << " to " << char_out << endl; break;
    		case 4: cout << "Deleting all " << char_in << " in " << filename << endl; break;
		}
	}
	
	// Add an additional 20 % to the file size to allow for expansion
	long			bufsize = (long) (filesize*1.2);
    unsigned char*	buf = new unsigned char[bufsize];
	
	fin.read((char *)buf, bufsize);
//	if ( fin.fail() )
//		cerr << "Error: File " << filename << " not read!" << endl;
	
	fin.close();
	
	if ( optind < argc ) {
		filename = argv[optind];
		if ( filename.contains(".") ) {
			if ( filename.contains("/") )
				filename = filename.post_rev('/');
		}
	} else {
		filename = 0;
		if ( verbose & VERB_LABEL ) cout << endl;
	}
    
	if ( replace_flag ) {
		if ( replace_flag < 4 )
			filesize = replace_a_character(buf, filesize, char_in, char_out);
		else
			filesize = delete_a_character(buf, filesize, char_in);
	}
	
	ofstream		fout;
    if ( filename.length() ) {
		fout.open(filename.c_str());
		if ( fout.fail() ) {
			cerr << "Error: File " << filename << " not opened!" << endl;
			bexit(-1);
		}
		if ( verbose ) {
			cout << "Writing:                        " << filename << endl;
			cout << "Output file size:               " << filesize << endl;
		}
    	fout.write((char *)buf, filesize);
    	fout.close();
	} else
		cout << buf << endl;
	
    delete[] buf;
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

long		replace_a_character(unsigned char* buf, long size, int chin, int chout)
{
	long		i, n;
	
	for ( i=n=0; i<size; i++ ) {
		if ( buf[i] == chin ) {
			buf[i] = chout;
			n++;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Characters replaced:            " << n << endl << endl;

	return i;
}

long		delete_a_character(unsigned char* buf, long size, int c)
{
	long		i, j, n;
	
	for ( i=j=n=0; i<size; i++ ) {
		if ( buf[i] != c )
			buf[j++] = buf[i];
		else
			n++;
	}
	
	buf[j] = 0;
	
	if ( verbose & VERB_PROCESS )
		cout << "Characters deleted:             " << n << endl << endl;

	return j;
}

