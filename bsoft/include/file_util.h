/**
@file	file_util.h 
@brief	Header file for file checking functions 
@author Bernard Heymann 
@date	Created: 20010601
@date	Modified: 20210413
**/

#include "Bstring.h"
#include <fstream>

#ifndef _filetype_ 
/**
@enum	FileType
@brief	File type enumeration.

The types are mostly based on Bsoft objects.
Additional types enumerated above 10 are for other packages.
**/
enum FileType { 
	Unknown_FileType = 0,	// Indeterminate
	Image = 1,				// Image file 
	Micrograph = 2,			// Micrograph parameter file 
	Molecule = 3,			// Molecular coordinate file 
	Model = 4,				// Model parameter file
	MgRelion = 11			// Relion micrograph
} ;
#define _filetype_ 
#endif 

// Function prototypes 
int			fread_large(unsigned char* aptr, size_t pagesize, size_t offset, ifstream* fimg);
int			fread_large(unsigned char* aptr, size_t pagesize, size_t offset, ifstream& fimg);
Bstring		find_file(Bstring filename, Bstring path, int flag=0);
string		find_file(string filename, string path, int flag=0);
vector<string>	file_list(string path);
vector<string>	file_list(string path, string ext);
FileType	file_type(const char* filename);
FileType	file_type(Bstring& filename);
int			detect_and_fix_carriage_return(const char* filename);

