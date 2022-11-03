/**
@file	file_util.cpp
@brief	Library functions for file checking 
@author Bernard Heymann 
@date	Created: 20070101
@date	Modified: 20210413
**/

#include "file_util.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>
#include <dirent.h>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bstring	test_access(Bstring filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG test_access: " << filename << endl;
	
	if ( access(filename.c_str(), F_OK) == 0 ) return filename;

	Bstring		empty;
	
	return empty;
}

string	test_access(string filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG test_access: " << filename << endl;
	
	if ( access(filename.c_str(), F_OK) == 0 ) return filename;

	string		empty;
	
	return empty;
}

/**
@brief 	Searches for a file within or without the given path.
@param 	filename	file name to search for.
@param 	path		path to search in.
@param	flag		if not found: bit 4 = warn; bit 5 = delete file name.
@return Bstring		found file, empty if not found and delete_flag set.

	The input filename is first tested for access.
	If not found, the filename without its original path is tested.
	If not found, the filename with the given path is tested.
	If not found, an error is reported and the original filename returned 
	unless the delete flag is set.

**/
Bstring		find_file(Bstring filename, Bstring path, int flag)
{
	if ( filename.length() < 1 ) return filename;

	if ( filename == "?" ) return filename;

	if ( filename.empty() ) {
		filename = "?";
		return filename;
	}

	if ( !path.empty() )
		if ( path[-1] != '/' ) path += "/";

	long		i;
	Bstring		foundfile = test_access(filename);
	Bstring		testfile, testpath;

	if ( foundfile.empty() && path.length() )
		foundfile = test_access(path + filename);

	if ( foundfile.length() ) {
		return foundfile;
	} else if ( filename.contains("://") ) {
		i = 3 + filename.find("://", 0);
		i = filename.find("/", i);
		foundfile = test_access(filename.substr(i, filename.length()));
	} else {
		if ( filename.contains("/") ) {
			testfile = filename.post_rev('/');
			testpath = filename.pre_rev('/') + "/";
			foundfile = test_access(testfile);
		} else {
			testfile = filename;
		}
		
		if ( foundfile.empty() && path.length() )
			foundfile = test_access(path + testfile);

		if ( foundfile.empty() && testpath.length() )
			foundfile = test_access(testpath + testfile);
		
		if ( foundfile.empty() && path.length() && testpath.length() ) {
			foundfile = test_access(path + testpath + testfile);
		}
	}
	
	if ( foundfile.empty() ) {
		if ( verbose && ( flag & 8 ) )
			cerr << "Warning: File " << filename << " not found!" << endl;
		if ( flag & 16 ) {
			if ( verbose && ( flag & 8 ) )
				cerr << tab << " deleting file name!" << endl;
		} else {
			foundfile = filename;
		}
	}
	
	return foundfile;
}

string		find_file(string filename, string path, int flag)
{
	if ( filename.length() < 1 ) return filename;

	if ( filename == "?" ) return filename;

	if ( filename.length() < 1 ) {
		filename = "?";
		return filename;
	}

	if ( path.length() )
		if ( path[-1] != '/' ) path += "/";

	long		i;
	string		foundfile = test_access(filename);
	string		testfile, testpath;
	
	if ( foundfile.length() ) {
		return foundfile;
	} else if ( filename.find("://") != string::npos ) {
		i = 3 + filename.find("://", 0);
		i = filename.find("/", i);
		foundfile = test_access(filename.substr(i, filename.length()));
	} else {
		if ( filename.find("/") != string::npos ) {
			testfile = filename.substr(filename.rfind('/'));
			testpath = filename.substr(0, filename.rfind('/')) + "/";
			foundfile = test_access(testfile);
		} else {
			testfile = filename;
		}
		
		if ( foundfile.length() < 1 && path.length() )
			foundfile = test_access(path + testfile);

		if ( foundfile.length() < 1 && testpath.length() )
			foundfile = test_access(testpath + testfile);
		
		if ( foundfile.length() < 1 && path.length() && testpath.length() ) {
			foundfile = test_access(path + testpath + testfile);
		}
	}
	
	if ( foundfile.length() < 1 ) {
		if ( verbose && ( flag & 8 ) )
			cerr << "Warning: File " << filename << " not found!" << endl;
		if ( flag & 16 ) {
			foundfile = filename;
			if ( verbose && ( flag & 8 ) )
				cerr << tab << " deleting file name!" << endl;
		}
	}
	
	return foundfile;
}

/**
@brief 	Returns a list of files in the requested directory.
@param 	&path			directory path.
@return vector<string>		list of file names.

**/
vector<string>	file_list(string path)
{
	vector<string>	files;
	if ( path.length() < 1 ) return files;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG file_list: " << path << endl;

	DIR 			*dir;
	struct dirent 	*entry;
	string			fn;
	
	if ( ( dir = opendir(path.c_str()) ) ) {
		while ( ( entry = readdir(dir) ) ) {
			fn = entry->d_name;
			if( fn[0] != '.' ) {
				fn = path + "/" + fn;
				files.push_back(fn);
			}
		}
		closedir(dir);
	}

	sort(files.begin(), files.end());

	return files;
}

/**
@brief 	Returns a list of files in the requested directory.
@param 	&path			directory path.
@param 	&ext			file name extension.
@return vector<string>		list of file names.

**/
vector<string>	file_list(string path, string ext)
{
	if ( ext.length() < 1 ) return file_list(path);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG file_list: " << path << tab << ext << endl;
	
	DIR 			*dir;
	struct dirent 	*entry;
	string			fn;
	vector<string>	files;
	
	if ( ( dir = opendir(path.c_str()) ) ) {
		while ( ( entry = readdir(dir) ) ) {
			fn = entry->d_name;
			if( fn[0] != '.' && extension(fn) == ext ) {
				fn = path + "/" + fn;
				files.push_back(fn);
			}
		}
		closedir(dir);
	}
	
	sort(files.begin(), files.end());
	
	return files;
}


/**
@brief 	Checks the file type using the extension and contents.
@param 	*filename		file name.
@return FileType			enumerated file type.

	The file extension is the main determinant of the file type.
	File formats with multiple types (such as the STAR and PDB formats)
	are distinguished based on content.

**/
FileType	file_type(const char* filename)
{
	Bstring		name(filename);
	return file_type(name);
}

FileType	file_type(Bstring& filename)
{
	FileType	type = Unknown_FileType;

	if ( filename.length() < 1 ) return type;
	
	ifstream	f;
	char		aline[MAXLINELEN];
	char*		aptr;
	string		sline;
	int			i, check_flag(0);

	Bstring		cleanname = filename;
	Bstring		ext;
	
	if ( filename.contains(":") ) {
		cleanname = filename.pre(':');
		ext = filename.post(':');
	} else {
		ext = cleanname.extension();
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG file_type: ext = " << ext << endl;
	
	if ( ext == "star" ) {
		f.open(cleanname.c_str());
		if ( f.fail() ) {
			error_show(cleanname.c_str(), __FILE__, __LINE__);
			return type;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG file_type: file " << cleanname << " opened" << endl;
		while ( !f.eof() && type == Unknown_FileType ) {
			getline(f, sline);
//			cout << "-" << sline << "-" << endl;
			for ( aptr = (char *) sline.c_str(); isspace(aptr[0]); aptr++ ) ;
			if ( strncmp(aptr, "data_", 5) == 0 ) check_flag = 1;
//			cout << check_flag << ": " << aptr << endl;
			if ( check_flag && aptr[0] != '#' ) {
				if ( strstr(aptr, "micrograph.id") ) type = Micrograph;
				if ( strstr(aptr, "map.3D_reconstruction.id") ) type = Micrograph;
				if ( strstr(aptr, "map.reference.file_name") ) type = Micrograph;
				if ( strstr(aptr, "_rln") ) type = MgRelion;	// Relion STAR
				if ( strstr(aptr, "model.id") ) type = Model;
				if ( strstr(aptr, "molecule.name") ) type = Molecule;
				if ( strstr(aptr, "image.file_name") ) type = Image;
			}
		}
		f.close();
	} else if ( ext == "xml" ) {
		f.open(cleanname.c_str());
		if ( f.fail() ) {
			error_show(cleanname.c_str(), __FILE__, __LINE__);
			return type;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG file_type: file " << cleanname << " opened" << endl;
		while ( !f.eof() && type == Unknown_FileType ) {
			getline(f, sline);
			for ( aptr = (char *) sline.c_str(); isspace(aptr[0]); aptr++ ) ;
			if ( strncmp(aptr, "<project", 8) == 0 ) type = Micrograph;
			if ( strncmp(aptr, "<model", 6) == 0 ) type = Model;
		}
		f.close();
	} else if ( ext == "emx" ) {
		type = Micrograph;
	} else if ( ext == "imod" ) {
		type = Micrograph;
	} else if ( ext == "mdoc" ) {
		type = Micrograph;
	} else if ( ext.contains("pdb") || ext.contains("ent") ) {
		type = Molecule;
		f.open(cleanname.c_str());
		if ( f.fail() ) {
			error_show(cleanname.c_str(), __FILE__, __LINE__);
			return type;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG file_type: file " << cleanname << " opened" << endl;
		while ( f.getline(aline, MAXLINELEN) && type == Unknown_FileType && strncmp(aline, "ATOM", 4) ) {
			if ( strncmp(aline, "REMARK", 6) == 0 && strstr(aline, "Model") ) type = Model;
		}
		f.close();
	} else if ( ext.contains("cmm") ) type = Model;
	else if ( ext.contains("v3d") ) type = Model;
    else if ( ext.contains("cif") ) type = Molecule;
	else if ( ext.contains("txt") ) type = Molecule;
	else if ( ext.contains("aln") ) type = Molecule;
	else if ( ext.contains("embl") ) type = Molecule;
	else if ( ext.contains("fasta") ) type = Molecule;
	else if ( ext.contains("gb") || ext.contains("gen") ) type = Molecule;
    else if ( ext.contains("gro") ) type = Molecule;
	else if ( ext.contains("phylip") ) type = Molecule;
	else if ( ext.contains("pir") ) type = Molecule;
	else if ( ext.contains("wh") || ext.contains("wah") ) type = Molecule;
	else if ( ext.contains("raw") ) type = Image;
	else if ( ext.contains("asc") || ext.contains( "txt") ) type = Image;
	else if ( ext.contains("bcr") ) type = Image;
	else if ( ext.contains("pic") ) type = Image;
	else if ( ext.contains("brx" ) || ext.contains("brix" ) ) type = Image;
	else if ( ext.contains("dat" ) ) type = Image;
	else if ( ext.contains("ccp") || ext.contains("map") ) type = Image;
	else if ( ext.contains("di") ) type = Image;
	else if ( ext.contains("dm") ) type = Image;
	else if ( ext.contains("dx") ) type = Image;
	else if ( ext.contains("omap" ) || ext.contains("dsn6") || ext.contains("dn6") ) type = Image;
	else if ( ext.contains("eer") ) type = Image;
	else if ( ext.contains("em") ) type = Image;
	else if ( ext.contains("pot") ) type = Image;
	else if ( ext.contains("grd") ) type = Image;
	else if ( ext.contains("hkl") ) type = Image;
	else if ( ext.contains("img") || ext.contains("hed") ) type = Image;
	else if ( ext.contains("ip") ) type = Image;
	else if ( ext.contains("jpg") || ext.contains("jpeg") ) type = Image;
	else if ( ext.contains("krn") ) type = Image;
	else if ( ext.contains("mif") ) type = Image;
	else if ( ext.contains("mff") ) type = Image;
	else if ( ext.contains("mrc") || ext == "stk" ) type = Image;
	else if ( ext == "st" || ext.contains("ali") || ext == ("rec") ) type = Image;	// IMOD extensions
	else if ( ext.contains("nd2") ) type = Image;
	else if ( ext.contains("pif") || ext.contains("sf") ) type = Image;
	else if ( ext.contains("bp") || ext.contains("bq") ) type = Image;
	else if ( ext.contains("png") ) type = Image;
	else if ( ext.contains("pbm") || ext.contains("pgm") || ext.contains("ppm") ) type = Image;
	else if ( ext.contains("ser") ) type = Image;
	else if ( ext.contains("sit") ) type = Image;
	else if ( ext.contains("spe") ) type = Image;
	else if ( ext.contains("spi") ) {		// Spider image and parameter files
		f.open(cleanname.c_str());
		if ( f.fail() ) {
			error_show(cleanname.c_str(), __FILE__, __LINE__);
			return type;
		}
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG file_type: file " << cleanname << " opened" << endl;
		for ( i=0; i<5 && type == Unknown_FileType; i++ ) {
			f.getline(aline, MAXLINELEN);
			aptr = aline;
			while ( aptr[0] && isspace(aptr[0]) ) aptr++;	// Remove leading spaces
			if ( aptr[0] == ';' ) type = Micrograph;
		}
		f.close();
		if ( type == Unknown_FileType ) type = Image;
	} else if ( ext.contains("spm") || ext.contains("sup") || ext == "f" ) type = Image;
	else if ( ext.contains("tif") ) type = Image;
	else if ( ext.contains("tga") ) type = Image;
	else if ( ext.contains("xpl") || ext.contains("rfl") ) type = Image;
	else cerr << "Error: File type with extension \"" << ext << "\" not supported!" << endl;

	return type;
}

/**
@brief 	Reads blocks of memeory no larger than 1Gb.
@param 	*aptr		pointer to pre-allocated memory.
@param 	pagesize	size of pre-allocated memory.
@param 	offset		offset in file.
@param 	*fimg		file pointer.
@return int			0, <0 if error.

	Each block is packed in sequence into the pre-allocated memory provided.

**/
int			fread_large(unsigned char* aptr, size_t pagesize, size_t offset, ifstream* fimg)
{
	size_t			j, readsize, pagemax = 1073741824;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG fread_large: pagesize=" << pagesize << " offset=" << offset << endl;
	
	fimg->seekg(offset, ios::beg);
	for ( j=0; j<pagesize; j+=pagemax ) {
		readsize = pagesize - j;
		if ( readsize > pagemax ) readsize = pagemax;
		fimg->read((char *)aptr, readsize);
		if ( fimg->fail() ) return -1;
		aptr += readsize;
	}
	
	return 0;
}

int			fread_large(unsigned char* aptr, size_t pagesize, size_t offset, ifstream& fimg)
{
	size_t			j, readsize, pagemax = 1073741824;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG fread_large: pagesize=" << pagesize << " offset=" << offset << endl;
	
	fimg.seekg(offset, ios::beg);
	for ( j=0; j<pagesize; j+=pagemax ) {
		readsize = pagesize - j;
		if ( readsize > pagemax ) readsize = pagemax;
		fimg.read((char *)aptr, readsize);
		if ( fimg.fail() ) return -1;
		aptr += readsize;
	}
	
	return 0;
}



/**
@brief 	Detects carriage returns in text files and converts them to new-lines.
@param 	*filename	file name.
@return int			0, <0 if error.

	The first line is read and if any carriage returns are found, the whole
	file is scanned and carriage returns converted to new-lines.

**/
int			detect_and_fix_carriage_return(const char* filename)
{
	long		max_line_len = 10240;
	char		aline[max_line_len];
	
	ifstream		fin(filename);
	if ( fin.fail() ) {
		error_show(filename, __FILE__, __LINE__);
		return -1;
	}
	
	fin.getline(aline, max_line_len);

	fin.close();

	if ( !strchr(aline, '\r') ) return 0;

	if ( access(filename, W_OK) ) {
		error_show("Error in detect_and_fix_carriage_return", __FILE__, __LINE__);
		cerr << "Warning: Carriage returns in " << filename << "! Fix with bfix." << endl;
		return -1;
	}
	
	cerr << "Warning: Fixing carriage returns in " << filename << endl;

	fin.open(filename);
	
	ofstream		fout("t.t");
	char			*pin, *pout;
	
	while ( !fin.eof() ) {
		fin.getline(aline, max_line_len);
		for ( pin = pout = aline; *pin != 0; pin++ ) {
			if ( *pin != '\r' ) {
				*pout = *pin;
				pout++;
			}
		}
		*pout = 0;
		fout << aline << endl;
	}
	
	fin.close();
	fout.close();
	
	rename("t.t", filename);
	
	return 0;
}

