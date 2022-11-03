/**
@file	rwstar.cpp
@brief	Parse and interpret STAR format files
@author Bernard Heymann
@date	Created: 19990605
@date	Modified: 20190819
**/

#include "rwstar.h"
#include "linked_list.h"
#include "file_util.h"
#include "utilities.h"
#include <iostream>
#include <stdarg.h>

using namespace std;

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes

char*		clean_line(char* string)
{
	unsigned int	i, j, start(0), last = strlen(string);
	
	if ( last < 1 ) return string;
	
	while ( isspace(string[start]) && start < strlen(string) ) start++;	// Find first non-whitespace
		
	char*		aptr = string;
	
	j = 0;
	for ( i=start; i<strlen(string); i++ ) {
		if ( string[i] != '\n' && string[i] != '\r' ) {	// Ignore newlines
			if ( string[i] < '\r' ) aptr[j] = ' ';		// All control characters become spaces
			else aptr[j] = string[i];
			if ( aptr[j] != ' ' ) last = j;
			j++;
		}
	}
	
	last++;
	aptr[last] = 0;
		
	return aptr;
}

string		clean_line(string s)
{
	size_t		i, j(0), start(0), last(s.length());
	
	if ( last < 1 ) return s;
	
//	while ( isspace(s[start]) && start < s.length() ) start++;	// Find first non-whitespace
	s.erase(0, s.find_first_not_of(" \t"));
		
	for ( i=start; i<s.length(); i++ ) {
		if ( s[i] != '\n' && s[i] != '\r' ) {	// Ignore newlines
			if ( s[i] < '\r' ) s[j] = ' ';		// All control characters become spaces
			else s[j] = s[i];
			if ( s[j] != ' ' ) last = j;
			j++;
		}
	}
	
	last++;
	s[last] = 0;
		
	return s;
}

int			no_square_brackets(Bstar_old* star)
{
	Bstar_block*	block;
	Bstar_item* 	item;
	
	if ( verbose & VERB_FULL )
		cout << "Removing square brackets from tags" << endl;
		
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			item->tag = item->tag.replace('[', '_');
			item->tag = item->tag.remove(']');
		}
	}
	
	return 0;
}

bool		tag_compare(Bstring item_tag, Bstring tag)
{
	if ( item_tag == tag ) return true;
	
	tag.replace('.', '_');
	
	if ( item_tag == tag ) return true;
	
	return false;
}

/**
@brief 	Creates a STAR data base.
@return Bstar_old*			the new STAR data base.

	A STAR structure is allocated.
	This function should be called before reading a STAR file, or before
	composing a STAR database for writing.

**/
Bstar_old*		init_star()
{
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG init_star: Initializing" << endl;
	
	// The main database structure and pointers to the data blocks
	Bstar_old*		star = new Bstar_old;
	star->split = 0;
	star->line_length = LINELENGTH;		// Length of line for output
	star->block = NULL;
	
	return star;
}

/**
@brief 	Destroys a STAR data base.
@param 	*star			the STAR data base.
@return int				0.

	A STAR data base structure and all of the data blocks and items are freed.

**/
int			kill_star(Bstar_old* star)
{
	if ( star == NULL ) return 0;
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG kill_star: Deleting STAR data base" << endl;
	
	Bstar_block*	block, *block2;
	
	for ( block=star->block; block; ) {
		block2 = block->next;
		kill_block(block);
		block = block2;
	}
	
	star->comment = 0;
	
	delete star;
	
	return 0;
}

/**
@brief 	Destroys a STAR data block.
@param	*block	pointer to the STAR data block.
@return int					0.

	A data block in a STAR data base structure and all of the items 
	associated with that data block are freed.

**/
int 		kill_block(Bstar_block* block)
{
	if ( block == NULL ) return 0;
	
	Bstar_item* 	item, *item2;
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG kill_block: Deleting block " << block->tag << endl;
	
	for ( item = block->item; item; ) {
		item2 = item->next;
		kill_item(item);
		item = item2;
	}
	
	block->tag = 0;
	
	block->comment = 0;
	
	block->filename = 0;
	
	delete[] (char*) block;
	
	return 0;
}

/**
@brief 	Destroys an item in a data block in a STAR data base.
@param	*item	the item.
@return int					0.

	An item and all of the items referenced are freed. 

**/
int 		kill_item(Bstar_item* item)
{		
	if ( item == NULL ) return 0;
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG kill_item: Deleting item " << item->tag << endl;
	
	string_kill(item->data);
	
	item->tag = 0;
	item->comment = 0;
	
	delete[] (char*) item;
	
	return 0;
}

/**
@brief 	Reads paramaters and data into a STAR data base from a list of files.
@param	*filename	a list of file names separated by commas.
@param	*star				an existing STAR data base.
@return int						error code (<0 means failure).

	Every data block is read separately and comments are preserved as far 
	as possible.

**/
int			read_star(const char* filename, Bstar_old* star)
{
	int				err(0);
	
	// Get the list of filenames
	Bstring			files(filename);
	Bstring*		file_list = files.split(",");
	
	err = read_star(file_list, star);
	
	string_kill(file_list);
	
	return err;
}

int			read_star(Bstring& filename, Bstar_old* star)
{
	Bstring			onefile(filename);
	return read_star(&onefile, star);
}

int			read_star(Bstring* file_list, Bstar_old* star)
{
	int				err(0);
	
	if ( !file_list ) {
		error_show("No file names found!", __FILE__, __LINE__);
		return -1;
	}
	
	// Get the list of filenames
	Bstring*		thisfile = NULL;
	
	if ( verbose & VERB_DEBUG_STAR ) {
		cout << "DEBUG read_star: STAR filenames: " << endl;
		for ( thisfile = file_list; thisfile; thisfile = thisfile->next )
			cout << " " << *thisfile;
		cout << endl;
	}
	
	// Parse every file
    ifstream*		fstar = new ifstream;
    char*       	aline = new char[MAXLINELEN+4];
	string			sline;
	memset(aline, 0, MAXLINELEN+4);
	char			*aptr;
	string			comment;
	int 			comment_done(0), comment_add(1);
	Bstar_block*	block;
	
	for ( thisfile = file_list; thisfile; thisfile = thisfile->next ) {
		if ( detect_and_fix_carriage_return(thisfile->c_str()) ) {
			cerr << "Warning: Carriage returns not fixed in " << *thisfile << endl;
//			error_show(thisfile->c_str(), __FILE__, __LINE__);
//			return -1;
		}

		if ( verbose & VERB_DEBUG )
			cout << "Opening " << *thisfile << endl;
	    fstar->open(thisfile->c_str());
	    if ( fstar->fail() ) {
			error_show(thisfile->c_str(), __FILE__, __LINE__);
			return -1;
		}
		
		if ( verbose & VERB_PROCESS )
			cout << "# Reading STAR file:              " << *thisfile << endl;
		
		comment_done = 0;
		while ( !comment_done ) {
			getline(*fstar, sline);
//			cout << "-" << sline << "-" << endl;
			if ( fstar->eof() ) {
				comment_done = 1;
			} else {
				sline.erase(0, sline.find_first_not_of(" \t"));
				if ( sline.size() > 4 && sline.compare(0, 5, "data_", 5) == 0 )
					comment_done = 1;
				else if ( comment_add ) comment = comment + sline + "\n";
			}
		}
		
		comment_add = 0;
		
//		cout << "comment length = " << comment.length() << endl;
		
		strncpy(aline, sline.c_str(), MAXLINELEN);
		aptr = clean_line(aline);	// Points to either a data block or nothing
		
		if ( strncmp(aptr, "data_", 5) == 0 ) {
			if ( !star->block ) star->block = read_block(fstar, aptr, *thisfile);
			else {
				for ( block = star->block; block->next; block = block->next ) ;
				block->next = read_block(fstar, aptr, *thisfile);
			}
		} else
			cerr << "Error: No data blocks in file " << *thisfile << endl;
		
		fstar->close();
	}
	
	delete fstar;
	delete[] aline;

	for ( block = star->block; block->next; block = block->next )
		block->next->number = block->number + 1;

	if ( comment.length() ) star->comment = comment;
	else star->comment = "\n";
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG read_star: Transferring comments:" << star->comment << endl << endl;
	
	no_square_brackets(star);
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG read_star: Done!" << endl;
	
	return err;
}

/**
@brief 	Reads paramaters and data into a STAR data block from an open file.
@param 	*fstar 			an open STAR format file.
@param 	*aptr			a pointer to the current line in the file.
@param 	&filename		file containing this block.
@return Bstar_block*	new block read.

	A block defines a unit of parameters or a unit of data. Every data 
	block is read separately and comments are preserved as far as possible.

**/
Bstar_block* read_block(ifstream* fstar, char* aptr, Bstring& filename)
{
	Bstar_block*	block = NULL;
	add_item(reinterpret_cast<char **>(&block), sizeof(Bstar_block));
	
	block->filename = filename;
	
	Bstar_item*		item = NULL;
	Bstar_item*		new_item = NULL;
	
	int				loop_counter(0), loop = -1;
	
	// Copy the block tag
	aptr += 5;
	if ( strlen(aptr) > 0 )
		block->tag = aptr;
	else
		block->tag = "";
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "# Reading block: " << block->tag << endl;
				
	char*       	aline = new char[MAXLINELEN+4];
	memset(aline, 0, MAXLINELEN+4);
	
	while ( fstar->getline(aline, MAXLINELEN) ) {
		aptr = clean_line(aline);
		if ( strncmp(aptr, "data_", 5) == 0 ) {
			block->next = read_block(fstar, aptr, filename);
		} else {
			if ( aptr[0] == '_' || aptr[0] == '#' || aptr[0] == ';' ) {	// Single item read
				new_item = read_single_item(fstar, aptr);
				loop = -1;
			} else if ( strncmp(aptr, "loop_", 5) == 0 ) { 	// A loop block
				new_item = read_loop_items(fstar);
				loop = ++loop_counter;
			}
			if ( new_item ) {
				if ( item ) item->next = new_item;
				else {
					block->item = item = new_item;
					item->loop = loop;
				}
				while ( item->next ) {
					item = item->next;
					item->loop = loop;
				}
				new_item = NULL;
			}
		}
	}
	
	delete[] aline;
	
	return block;
}

/**
@brief 	Reads a single-valued item or comment line into a STAR data block 
	from an open file.
@param 	*fstar 			an open STAR format file.
@param 	*aline			a pointer to the current line in the file.
@return Bstar_item* 	new STAR item.

	All tags with single values and outside loops are interpreted here.
	A single item is defined where the first non-space character on a
	line is an underscore. A comment is defined by a '#' or ';' as the
	first character on the line.
	Note: The loop flag field of the STAR item is equal to -1.

**/
Bstar_item* read_single_item(ifstream* fstar, char* aline)
{
	Bstar_item*		item = NULL;
	add_item(reinterpret_cast<char **>(&item), sizeof(Bstar_item));
	item->type = StringItem;
	
	int				k, digit(0), alpha(0);
	Bstring			svalue;
	char*			aptr = aline;
	char			temp[100];

	if ( aptr[0] == '#' && strlen(aptr) < 2 ) return 0;
	if ( aptr[0] == '_' && strlen(aptr) < 4 ) return 0;
		
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG single item: Starting" << endl;

	item->loop = -1;
	
	if ( aptr[0] == '#' || aptr[0] == ';' ) {	// A comment
		item->tag = "comment";
	} else {
		sscanf(aptr, "%s%n", temp, &k);			// Read tag
		item->tag = &temp[1];					// Omit the initial underscore
		aptr = clean_line(aline + k);			// Start of data item
		if ( strlen(aptr) < 1 ) {				// Next line
			fstar->getline(aline, MAXLINELEN);
			aptr = clean_line(aline);
		}
	}

	if ( aptr[0] == '#' ) {						// Comment
		aptr++;
		svalue = aptr;
	} else if ( aptr[0] == '\'' ) {				// String delimiters
		aptr++;
		k = strcspn(aptr, "\'");
		aptr[k] = 0;
		svalue = aptr;
	} else if ( aptr[0] == '\"' ) {
		while ( aptr[0] == '\"' ) aptr++;		// Take care of multiple quotes
		k = strcspn(aptr, "\"");
		aptr[k] = 0;
		svalue = aptr;
	} else if ( aptr[0] == ';' ) {
		aptr++;
		k = strcspn(aptr, "\n");
		aptr[k] = 0;
		if ( k > 0 ) svalue = aptr;
		while ( fstar->getline(aline, MAXLINELEN) && aline[0] != ';' ) {
			k = strcspn(aline, "\n");
			aline[k] = 0;
			if ( k > 0 ) svalue += aline;
		}
	} else {
		k = strcspn(aptr, "#");
		aptr[k] = 0;
		svalue = aptr;
	}
	
	if ( svalue.length() < 1 ) {
		cerr << "Warning: Tag " << item->tag << " value not read!" << endl;
		kill_item(item);
		item = NULL;
	}
	
	svalue = svalue.remove('\"');
	
	string_add(&item->data, svalue);
	item->maxlen = item->data->length();
	for ( k=0; k<svalue.length(); k++ ) {	// Check if string or number
		if ( isalpha(svalue[k]) ) alpha++;
		else if ( isdigit(svalue[k]) ) digit++;
	}
	if ( alpha < 1 && digit > 0 ) item->type = NumberItem;
	
	if ( verbose & VERB_DEBUG_STAR ) {
		cout << "DEBUG single item: Tag: " << item->tag << ", maxlen=" << item->maxlen << endl;
		cout << "DEBUG single item: Data: " << *(item->data) << endl;
	}
	
	return item;
}

/**
@brief 	Reads a loop structure with multiple items into a STAR data block 
	from an open file.
@param 	*fstar 			an open STAR format file.
@return Bstar_item* 	first new item in list.

	The loop is read line by line, checking to get every column value in a row.
	A row may extend over multiple lines, as long as it contains the number values
	specified by the number of tags at the beginning of the loop.
	A multiple line string value must be enclosed in ";" as the first character
	in the lines before and after the string.
	The loop ends with an empty line or when too few values occur in a row.
	An empty line at the end of a loop is required.
	Note: The loop flag field of the STAR item is equal to the item index
	of the first item in the loop.

**/
Bstar_item* read_loop_items(ifstream* fstar)
{
	Bstar_item* 	item_list = NULL;
	Bstar_item* 	item = NULL;
	
    char*       	aline = new char[MAXLINELEN+4];
	memset(aline, 0, MAXLINELEN+4);
    char*       	aval = new char[MAXLINELEN+4];
	memset(aval, 0, MAXLINELEN+4);
	char*			aptr = aline;
	int 			i, j, k, loopdone(0), number(0);
	int 			digit, alpha;
	int 			nitems(0);
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG loop: Starting" << endl;

	while ( fstar->getline(aline, MAXLINELEN) &&
				( aptr = clean_line(aline) ) &&
				aptr[0] == '_' ) {			// Find all the tags
		item = reinterpret_cast<Bstar_item *> (add_item(reinterpret_cast<char **>(&item_list), sizeof(Bstar_item)));
		sscanf(aptr, "%s%n", aval, &k);		// Read tag
		item->tag = &aval[1];				// Omit the initial underscore
		item->type = NumberItem;			// Default is NumberItem --> check later if it should be StringItem
		nitems++;
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG loop: number of tags=" << nitems << endl;
	
	if ( nitems < 1 ) {
		cerr << "Error: No tags in the loop!" << endl;
		return NULL;
	}

	// Strings to hold the values temporarily
	Bstring*	sval = new Bstring[nitems];
	Bstring**	data = new Bstring*[nitems];
	Bstring**	datalist = new Bstring*[nitems];
	memset(sval, 0, nitems*sizeof(Bstring));
	memset(data, 0, nitems*sizeof(Bstring*));
	memset(datalist, 0, nitems*sizeof(Bstring*));

	loopdone = number = 0;
	if ( strlen(aptr) < 1 || aptr[0] == '\n' || aptr[0] == '#' )	// Check for blank line
			loopdone = 1;
	while ( !loopdone ) {				// Get the number of rows in the loop
		if ( verbose & VERB_DEBUG_STAR )
			cout << "DEBUG loop: starting " << number << " " << aline << endl;
// I cannot take out commas - it's not generally only used as delimiter
//		for ( j=0; j<strlen(aptr); j++ ) if ( aptr[j] == ',' ) aptr[j] = ' ';
		for ( i=0; i<nitems && !loopdone; i++ ) {
			aptr = clean_line(aptr);
			if ( strlen(aptr) < 1 ) {							// Loop columns extend over multiple lines
				if ( fstar->getline(aline, MAXLINELEN) )			// Check for end-of-file
					aptr = clean_line(aline);
				else  loopdone = 1;
			}
			if ( aline[0] == ';' ) {							// Multi-line string items
				aptr = clean_line(aline + 1);
				j = strcspn(aptr, "\n");
				aptr[j] = 0;
				sval[i] = aptr;
				while ( fstar->getline(aline, MAXLINELEN) && aline[0] != ';' ) {
					j = strcspn(aline, "\n");
					aline[j] = 0;
					sval[i] = sval[i] + aline;
				}
			} else {
				if ( sscanf(aptr, "%s", aval) < 1 ) {			// Short items in loop
					if ( fstar->getline(aline, MAXLINELEN) )
						aptr = clean_line(aline);
					if ( sscanf(aptr, "%s", aval) < 1 ) 		// Check for end-of-loop
						loopdone = 1;
					else sval[i] = aval;
				} else sval[i] = aval;
			}
			if ( ( aptr = strchr(aptr, ' ') ) == NULL )
				aptr = aline + strlen(aline);
		}
		if ( verbose & VERB_DEBUG_STAR ) {
			cout << "DEBUG loop: finishing " << number << endl;
			for ( i=0; i<nitems; i++ ) cout << tab << sval[i] << endl;
			cout << endl;
		}
		if ( i == nitems && !loopdone ) {
			number++;
			for ( j=0, item=item_list; j<nitems; j++, item=item->next ) {
				sval[j] = sval[j].remove('\"');
				if ( item->maxlen < (int)sval[j].length() ) item->maxlen = sval[j].length();
				data[j] = string_add(&data[j], sval[j]);
				if ( !datalist[j] ) datalist[j] = data[j];
				digit = alpha = 0;
				for ( k=0; k<(int)sval[j].length(); k++ ) {		// Check if it is a string
					if ( isalpha(sval[j][k]) ) alpha++;
					else if ( isdigit(sval[j][k]) ) digit++;
				}
				if ( alpha > 0 ) item->type = StringItem;
			}
			if ( fstar->getline(aline, MAXLINELEN) )			// Check for end-of-file
				aptr = clean_line(aline);
			else loopdone = 1;
		} else loopdone = 1;								// Check for end-of-loop
		if ( strlen(aptr) < 1 || aptr[0] == '\n' || aptr[0] == '#' )	// Check for blank line
			loopdone = 1;
	}
	
	for ( i=0, item=item_list; i<nitems && item; i++, item=item->next )		// Assign the items
		item->data = datalist[i];

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG loop: Columns=" << nitems << " Rows=" << number << endl;

	delete[] aline;
	delete[] aval;
	delete[] sval;	
	delete[] data;
	delete[] datalist;
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG loop: memory freed" << endl;
	
	return item_list;
}

/**
@brief 	Writes a STAR data base to one or more STAR format files.
@param	*filename	the base file name (can be NULL if star->split == 9).
@param	*star		the STAR data base.
@return int			error code (<0 means failure).

	The STAR data base structure contains a flag (star->split) to indicate
	whether one or multiple files should be written.  In the case of multiple
	files, the base name is taken from the input file name, with an underscore
	and a number appended. The length of the number is determined by the
	star->split variable.

**/
int			write_star(const char* filename, Bstar_old* star)
{
	Bstring		thefile(filename);
	return write_star(thefile, star);
}

int			write_star(Bstring& filename, Bstar_old* star)
{
	int				err(0);
	
	if ( star->split < 9 )
		if ( filename.length() < 1 ) return -1;
	
	ofstream*		fstar = new ofstream;
	Bstring			default_comment = "# Written by Bsoft";
	
	if ( star->comment.length() < 1 ) star->comment = default_comment;
	
	if ( star->split == 0 ) {
	    fstar->open(filename.c_str());
		if ( fstar->fail() ) {
			cerr << "Error: Not able to write " << filename << endl;
			return -1;
		}
	
		if ( verbose & VERB_PROCESS )
			cout << "# Writing file:                 " << filename << endl;
	
		*fstar << star->comment << endl;
	} else {
		if ( verbose & VERB_PROCESS )
			cout << "Splitting STAR file (" << star->split << "):" << endl;
	}
	
	int				i;
	Bstar_block*	block;
	Bstring			blockname;
	char			format[32];
	
	snprintf(format, 32, "_%%0%dd.", star->split);
	
	for ( i=0, block = star->block; block; block = block->next, i++ ) if ( block->item ) {
		if ( star->split ) {
			if ( star->split == 9 ) {
				if ( block->tag.length() > 0 ) {
					blockname = block->tag;
					blockname = blockname.replace(' ', '_');
					blockname += ".star";
				} else {
					blockname = Bstring(i+1, "%09d.star");
				}
			} else {
				blockname = filename.pre_rev('.') + Bstring(i+1, format) + filename.post_rev('.');
		    }
			fstar->open(blockname.c_str());
			if ( fstar->fail() ) {
				cerr << "Error: Not able to write " << blockname << endl;
				return -1;
			}
			if ( verbose & VERB_PROCESS )
				cout << "# Writing file:                 " << blockname << endl;
			*fstar << star->comment << endl;
		}
		err = write_block(fstar, block, star->line_length);
		if ( star->split ) fstar->close();
	}

	if ( star->split == 0 ) fstar->close();
	
	delete fstar;
		
	return err;
}

/**
@brief 	Writes paramaters and data from a STAR data block into an open file.
@param 	*fstar 		an open STAR format file.
@param 	*block		a data block.
@param 	linelength	output maximum line length.
@return int			error code (<0 on error).

	A block defines a unit of parameters or a unit of data.

**/
int			write_block(ofstream* fstar, Bstar_block* block, int linelength)
{
	int				err(0);
	long			j, m, size, nitems, pos;
	Bstring			tag, tstr;
	Bstar_item*		item, *item2;
	Bstring**		dptr;
	
	if ( block->tag.empty() ) block->tag = "";
	*fstar << endl << "data_" << block->tag << endl << endl;
	
//	cout << block->item->tag << endl;
	
	for ( item = block->item; item; item = item->next ) {
		if ( item->tag.contains("comment") ) {
			tag = "\n# ";
		} else if ( item->tag.length() > 0 ) {
			if ( item->tag[0] == '_' ) tag = item->tag;
			else  tag = "_" + item->tag;
		}
		if ( item->loop < 0 ) {		// Write single value items
			if ( item->data->length() > linelength - item->tag.length() - 8 ) {	// Multi-line items
				if ( tag.length() > 0 )
					*fstar << tag << endl;
				*fstar << ";" << endl;
				for ( j=0; j<item->data->length(); j += linelength ) {
					size = item->data->length() - j;
					if ( size > linelength ) size = linelength;
					if ( size ) {
						tstr = item->data->substr(j, size);
						*fstar << tstr << endl;
					}
				}
				*fstar << ";" << endl;
			} else {	// Single items fitting on one or two lines
				if ( tag.length() > 0 ) {
					if ( item->type == StringItem ) {
						if ( item->data->length() > linelength - tag.length() - 8 ) {
							*fstar << tag << endl << "\"" << item->data << "\"" << endl;
						} else if ( strchr(item->data->c_str(), ' ') ) {
							*fstar << setw(40) << tag << " \"" << *(item->data) << "\"" << endl;
						} else
							*fstar << setw(40) << tag << " " << *(item->data) << endl;
					} else {
						*fstar << setw(40) << tag << " " << *(item->data) << endl;
					}
				} else {
					*fstar << *(item->data) << endl;	// No-tag data
				}
			}
		} else if ( item->data ) {			// Write a loop block
			*fstar << endl << "loop_" << endl;
			for ( nitems=0, item2=item; item2 && item2->loop == item->loop; item2=item2->next, nitems++ ) {
				if ( item2->tag[0] == '_' ) tag = item2->tag;
				else tag = "_" + item2->tag;
				*fstar << tag << endl;
			}
			dptr = new Bstring*[nitems];
			for ( j=0, item2=item; j<nitems; item2=item2->next, j++ )
				dptr[j] = item2->data;
			while ( dptr[0] ) {
				pos = 0;
				tstr = 0;
				for ( j=0, item2=item; j<nitems; item2=item2->next, j++ ) {
					if ( item2->maxlen > linelength ) {	// Multi-line items
						tstr += ";\n";
/*						for ( m=0; m<item2->maxlen; m += linelength ) {
							size = item2->maxlen - m;
							if ( size > linelength )
								size = linelength;
							if ( size > dptr[j]->length() - m )
								size = dptr[j]->length() - m;
							if ( size ) {
								tstr += dptr[j]->substr(m, size) + "\n";
							}
						}*/
						for ( m=0; m<dptr[j]->length(); m += linelength ) {
							size = dptr[j]->length() - m;
							if ( size > linelength )
								size = linelength;
							if ( size ) {
								tstr += dptr[j]->substr(m, size) + "\n";
							}
						}
						if ( j == nitems - 1 )
							tstr += ";";
						else
							tstr += ";\n";
						pos = 0;
					} else {
						if ( item2->type == NumberItem ) {	// Right justify
							tstr += " " + dptr[j]->right(item2->maxlen);
						} else {	// Write strings as is
//							tstr += *dptr[j] + " ";
							if ( dptr[j]->contains("'") ) *dptr[j] = "\"" + *dptr[j] + "\"";
							if ( item2->maxlen < dptr[j]->length() ) item2->maxlen = dptr[j]->length();
							tstr += " " + dptr[j]->left(item2->maxlen);
//							pos += dptr[j]->length() + 1;
						}
						pos += item2->maxlen + 1;
					}
					if ( pos + item2->maxlen > linelength && j < nitems - 1 ) {
						tstr += "\n";
						pos = 0;
					}
					dptr[j] = dptr[j]->next;
				}
				*fstar << tstr << endl;
			}
			for ( j=1; item->next && j<nitems; j++ ) item = item->next;
			delete[] dptr;
			*fstar << endl;
		}
	}
	*fstar << endl;
	
	return err;
}

/**
@brief 	Puts a set of strings and time in the main comment of the STAR data base.
@param	*star 		the STAR data base.
@param 	n			the number of strings.
@param	**strings	an array of strings.
@return int 		string length of the new comment.

	This is designed to pack the command line into a string followed by
	a second string for the time.

**/
int 		star_update_comment(Bstar_old* star, int n, char** strings)
{
	time_t		ti = time(NULL);
	
	int			i; 
	Bstring		new_comment("# ");
	
	for ( i=0; i<n; i++ )
		new_comment = new_comment + strings[i] + " ";
	
	new_comment = new_comment + "\n# " +  asctime(localtime(&ti)) + "\n";
	
	if ( verbose & VERB_DEBUG )
		cout << endl << new_comment << endl;
	 
	star->comment += new_comment;
	
	return star->comment.length();
}

/**
@brief 	Lists the command lines in the STAR header.
@param	*star 		the STAR data base.
@param	len 		maximum line length, infinite if zero.
@return long 		the number of command lines.

**/
long		star_list_comments(Bstar_old* star, long len)
{
	if ( len < 1 ) len = 10000000000;
	
	long		n(0);
	Bstring*	s;
	Bstring		s1;
	
	Bstring*	strarr = star->comment.split("\n");
	
	for ( s = strarr; s; s = s->next ) {
		if ( (*s)[0] == '#' ) {
			s1 = s->substr(2,len);
			cout << s1 << endl;
			n++;
		}
	}

	string_kill(strarr);

	return n;
}

/**
@brief 	Sets the maximum string lengths for each item in each block.
@param	*star 		the STAR data base.
@return int 		0.

	This is designed to clean up after creating a STAR database.

**/
int			star_set_string_lengths(Bstar_old* star)
{
	Bstar_block*	block;
	Bstar_item* 	item;
	Bstring*		data;
	
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			for ( data = item->data; data; data = data->next )
				if ( item->maxlen < data->length() )
					item->maxlen = data->length();
		}
	}
	
	return 0;
}

/**
@brief 	Replaces one tag with another in the STAR data base.
@param	*star 		the STAR database.
@param	*tag		old tag.
@param	*newtag	new	tag.
@return int 		0.

	The item with a given tag has the tag replaced with a new one.
	If an item with the new tag exists, it is deleted first.

**/
int 		item_change_tag(Bstar_old* star, const char* tag, const char* newtag)
{
	Bstar_block*	block;
	Bstar_item* 	item;
	
	if ( verbose & VERB_FULL )
		cout << "Changing all tags \"" << tag << "\" to \"" << newtag << "\"" << endl;
		
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag ) {
				item_delete_from_block(block, newtag);	// First make sure any other newtag is deleted
				item->tag = newtag;
			}
		}
	}
	
	return 0;
}

/**
@brief 	Prints the list of tags in the STAR data base.
@param	*star 		the STAR database
@return int 		number of tags.
**/
int 		show_tags(Bstar_old* star)
{
	int 			ntags(0);
	Bstring			tag;
	Bstar_block*	block;
	Bstar_item* 	item;
	
	for ( block = star->block; block; block = block->next ) {
		cout << endl << "data_" << block->tag << endl;
		cout << endl << "loop_" << endl;
		cout << "_star_tag" << endl;
		cout << "_star_type" << endl;
		for ( item = block->item; item; item = item->next ) {
			tag = "\"_" + item->tag + "\"";
			cout << setw(40) << item->tag;
			if ( item->type == StringItem )
				cout << "String" << endl;
			else
				cout << "Number" << endl;
			ntags++;
		}
		cout << endl;
	}
	cout << endl;
	
	return ntags;
}

/**
@brief 	Gets the STAR item associated with a tag in a STAR data block.
@param	*block		a STAR data block.
@param	*tag		a STAR tag string.
@return int			the STAR item index, -1 if the tag doesn't exist.

	The items in the STAR data block are traversed to find the item
	associated with a STAR tag defined in a header file.

**/
int			item_index(Bstar_block* block, const char* tag)
{
	int				i;
	Bstar_item* 	item = NULL;

	for ( i=0, item = block->item; item && !tag_compare(item->tag, tag); item = item->next ) i++;

	if ( !item ) i = -1;
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_index: " << tag << " " << i << endl;
	
	return i;
}

/**
@brief 	Gets the first block associated with a tag in a STAR data base.
@param	*star 		the STAR data base.
@param	*tag		a STAR tag string.
@return Bstar_block* the block, NULL if not found.
**/
Bstar_block*	block_find_with_tag(Bstar_old* star, const char* tag)
{
	int				found(0);
	Bstar_block*	block = NULL;
	Bstar_item* 	item = NULL;
	
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag ) {
				found = 1;
				break;
			}
		}
		if ( found ) break;
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_first_block: " << tag << " " << &block << endl;
	
	return block;
}


/**
@brief 	Gets the STAR item associated with a tag in a STAR data block.
@param	*block	a STAR data block.
@param	*tag		a STAR tag string.
@return Bstar_item* the STAR item, NULL if the tag doesn't exist.

	The items in the STAR data block are traversed to find the item
	associated with a STAR tag defined in a header file.

**/
/*Bstar_item* item_find(Bstar_block* block, const char* tag)
{
	Bstar_item* 	item = NULL;

	for ( item = block->item; item && !tag_compare(item->tag, tag); item = item->next ) ;

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_find: " << tag << " " << &item << endl;
	
	return item;
}
*/
Bstar_item* item_find(Bstar_block* block, Bstring tag)
{
	Bstar_item* 	item = NULL;

	for ( item = block->item; item && !tag_compare(item->tag, tag); item = item->next ) ;

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_find: " << tag << " " << &item << endl;
	
	return item;
}

/**
@brief 	Gets the STAR item associated with a tag in a STAR data block.
@param	*block		a STAR data block.
@param	*tag		a STAR tag string.
@return Bstar_item* the STAR item, NULL if the tag doesn't exist.

	The items in the STAR data block are traversed to find the item
	associated with a STAR tag defined in a header file.

**/
Bstar_item* item_find_or_make(Bstar_block* block, const char* tag)
{
	Bstar_item* 	item = item_find(block, tag);

//	if ( item ) {
//		string_kill(item->data);
//		item->data = NULL;
//	} else {
	if ( !item ) {
		item = reinterpret_cast<Bstar_item *> (add_item(reinterpret_cast<char **>(&block->item), sizeof(Bstar_item)));
		item->tag = tag;
	}

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_find_or_make: " << tag << " " << &item << endl;
	
	return item;
}

/**
@brief 	Gets the number of values associated with a tag in a STAR data base.
@param	*star 		the STAR data base.
@param	*tag		a STAR tag string.
@return long	the number of values.

	All STAR data blocks are traversed to count the number of
	values associated with a STAR tag defined in a header file.

**/
long	item_get_number(Bstar_old* star, const char* tag)
{
	long	number(0);
	Bstar_block*	block;
	Bstar_item* 	item;
	
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag )
				number += count_list((char *)item->data);
		}
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_number: " << tag << " " << &number << endl;
	
	return number;
}

/**
@brief 	Gets the number of values associated with a tag in a STAR data block.
@param	*block	block in the STAR database.
@param	*tag		a STAR tag string.
@return long	the number of values.

	The STAR data block is traversed to count the number of
	values associated with a STAR tag defined in a header file.

**/
long	item_get_number_for_block(Bstar_block* block, const char* tag)
{
	long	number(0);
	Bstar_item* 	item;
	
	for ( item = block->item; item; item = item->next ) {
		if ( item->tag == tag )
			number += count_list((char *)item->data);
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_number_for_block: " << tag << " " << &number << endl;
	
	return number;
}

/**
@brief 	Gets a string value associated with a tag in a STAR data block.
@param	*block	block in the STAR database.
@param	*tag		a STAR tag string.
@return char* 		the string value, NULL if the tag doesn't exist.

	The STAR data block is traversed to obtain the string
	value associated with a STAR tag defined in a header file.

**/
char* 		item_get_string(Bstar_block* block, const char* tag)
{
	int				i;
	Bstar_item* 	item;
	char* 			string = NULL;

	for ( item = block->item; item && !string; item = item->next ) {
		if ( item->tag == tag ) {
			string = new char[item->data->length()+1];
			for ( i=0; i<item->data->length(); i++ ) string[i] = (*item->data)[i];
			string[item->data->length()] = 0;
		}
	}

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_string: " << tag << " " << string << endl;
	
	return string;
}

/**
@brief 	Copies a string value associated with a tag in a STAR data block.
@param 	*string		destination string - must exist on stack or allocated.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@return int 		0.

	The STAR data blocks is traversed to obtain the string
	value associated with a STAR tag defined in a header file.

**/
int 		item_copy_string(char* string, Bstar_block* block, const char* tag)
{
	Bstar_item* 	item;

	for ( item = block->item; item; item = item->next ) {
		if ( item->tag == tag )
			strcpy(string, item->data->c_str());
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_copy_string: " << tag << " " << string << endl;
	
	return 0;
}

/**
@brief 	Gets an integer value associated with a tag in a STAR data block.
@param	*block	block in the STAR database.
@param	*tag		a STAR tag string.
@return int 		the integer value, 0 if the tag doesn't exist.

	The STAR data block is traversed to obtain the first integer
	value associated with a STAR tag defined in a header file.

**/
int 		item_get_integer(Bstar_block* block, const char* tag)
{
	Bstar_item* 	item;
	int 			value(0);

	for ( item = block->item; item; item = item->next ) {
		if ( item->tag == tag ) {
			if ( *(item->data) == "?" || *(item->data) == "." ) value = 0;
			else value = item->data->integer();
		}
	}

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_integer: " << tag << " " << value << endl;
	
	return value;
}

/**
@brief 	Gets a floating point value associated with a tag in a STAR data block.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@return float 		the floating point value, 0 if the tag doesn't exist.

	The STAR data block is traversed to obtain the floating point
	value associated with a STAR tag defined in a header file.

**/
float 		item_get_float(Bstar_block* block, const char* tag)
{
	Bstar_item* 	item;
	float 			value(0);

	for ( item = block->item; item; item = item->next ) {
		if ( item->tag == tag ) {
			if ( *(item->data) == "?" || *(item->data) == "." ) 
				value = strtod ("NAN()",NULL);
			else value = item->data->real();
		}
	}

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_get_float: " << tag << " " << value << endl;
	
	return value;
}

/**
@brief 	Writes a string into a data item associated with a specific data block
	and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	*string		string value.
@return int 		error code.
**/
int 		item_put_string(Bstar_block* block, const char* tag, char* string)
{
	Bstring 		this_string(string);
	return item_put_string(block, tag, this_string);
}

int 		item_put_string(Bstar_block* block, const char* tag, Bstring& string)
{
	if ( !block ) return -1;

	long	len = string.length();

	if ( len < 1 ) return 0;
	
	Bstar_item* 	item = item_find_or_make(block, tag);

	item->type = StringItem;
	item->loop = -1;
	
	string_add(&item->data, string);
	
	if ( item->maxlen < len ) item->maxlen = len;

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_string: " << tag << " " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Writes a list of strings into a data item associated with 
	a specific data block and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	*list		list of strings.
@return int			error code.
**/
int 		item_put_string_list(Bstar_block* block, const char* tag, Bstring* list)
{
	if ( !block ) return -1;

	Bstar_item* 	item = item_find_or_make(block, tag);
	Bstring*		data = NULL;
	Bstring*		string = NULL;

	item->type = StringItem;
	
	for ( string = list; string; string = string->next ) {
		data = string_add(&data, *string);
		if ( !item->data ) item->data = data;
		if ( item->maxlen < string->length() ) item->maxlen = string->length();
	}

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_string_list: " << tag << " " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Writes a integer into a data item associated with a specific data block
	and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	value		integer value.
@param	*format		string format.
@return int			error code.
**/
int 		item_put_integer(Bstar_block* block, const char* tag, int value, const char* format)
{
	if ( !block ) return -1;
	
	char			string[32] = "";
	snprintf(string, 32, format, value);

	long	len = strlen(string);

	Bstar_item* 	item = item_find_or_make(block, tag);

	item->type = NumberItem;
	item->loop = -1;
	
	string_add(&item->data, string);

	if ( item->maxlen < len ) item->maxlen = len;

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_integer: " << tag << " " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Writes a floating point value into a data item associated with a 
	specific data block and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	value		integer value.
@param	*format		string format.
@return int			error code.
**/
int 		item_put_float(Bstar_block* block, const char* tag, float value, const char* format)
{
	if ( !block ) return -1;
	
	char			string[32] = "";
	snprintf(string, 32, format, value);
	if ( isnan(value) ) snprintf(string, 32, ".");

	long	len = strlen(string);

	Bstar_item* 	item = item_find_or_make(block, tag);

	item->type = NumberItem;
	item->loop = -1;
	
	string_add(&item->data, string);

	if ( item->maxlen < len ) item->maxlen = len;

	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_float: " << tag << " " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Writes a list of string or numeric values into a data item associated with 
	a specific data block and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	*list		linked list.
@param 	offset		offset of structure element.
@param	*format		string format.
@return int			error code.
**/
int 		item_put_list(Bstar_block* block, const char* tag, char* list, 
				size_t offset, const char* format)
{
	if ( !block ) return -1;

	Bstar_item* 	item = item_find_or_make(block, tag);
	Bstring*		data = item->data;
	if ( data ) while ( data->next ) data = data->next;

	item->type = NumberItem;
	if ( strstr(format, "s") ) item->type = StringItem;

	char			string[32] = "";
	char**			it;

	for ( it = (char **) list; it; it = (char **) *it ) {
		if ( strstr(format, "s") )
			snprintf(string, 32, format, (char *)it + offset);
		else if ( strstr(format, "d") )
			snprintf(string, 32, format, *((int*)((char *)it + offset)));
		else if ( strstr(format, "lf") )
			snprintf(string, 32, format, *((double*)((char *)it + offset)));
		else
			snprintf(string, 32, format, *((float*)((char *)it + offset)));
		if ( strstr(string, "nan") ) snprintf(string, 32, ".");
		data = string_add(&data, string);
		if ( !item->data ) item->data = data;
		if ( item->maxlen < strlen(string) ) item->maxlen = strlen(string);
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_list: " << tag << " (" << offset << ") " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Writes a list of floating point numbers into a data item associated 
	with a specific data block and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	number		the number of values.
@param 	*value		a list of floating point values.
@param	*format		string format.
@return int			0.

	NaN values are taken as missing numbers and indicated by a '.' in the
	STAR file.

**/
int 		item_put_float_array(Bstar_block* block, const char* tag, int number, 
				float* value, const char* format)
{
	if ( !block || number < 1 || !value ) return -1;

	Bstar_item* 	item = item_find_or_make(block, tag);
	Bstring*		data = item->data;
	if ( data ) while ( data->next ) data = data->next;

	item->type = NumberItem;

	char			string[32] = "";
	int				i;

	for ( i=0; i<number; i++ ) {
		snprintf(string, 32, format, value[i]);
		if ( strstr(string, "nan") ) snprintf(string, 32, ".");
		data = string_add(&data, string);
		if ( !item->data ) item->data = data;
		if ( item->maxlen < strlen(string) ) item->maxlen = strlen(string);
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_float_array: " << tag << " " << *(item->data) << endl;
		
	return 0;
}

/**
@brief 	Writes a list of angular values in degrees into a data item associated with 
	a specific data block and tag in a STAR data base.
@param	*block		block in the STAR database.
@param	*tag		a STAR tag string.
@param 	*list		linked list.
@param 	offset		offset of structure element.
@param	*format		string format.
@return int			error code.

	Each angle in radians is first converted to degrees before writing it
	into a string.

**/
int 		item_put_angle_list(Bstar_block* block, const char* tag, char* list, 
				size_t offset, const char* format)
{
	if ( !block ) return -1;

	Bstar_item* 	item = item_find_or_make(block, tag);
	Bstring*		data = item->data;
	if ( data ) while ( data->next ) data = data->next;

	item->type = NumberItem;

	double			angle;
	char			string[32] = "";
	char**			it;

	for ( it = (char **) list; it; it = (char **) *it ) {
		if ( strstr(format, "lf") )
			angle = *((double*)((char *)it + offset));
		else
			angle = *((float*)((char *)it + offset));
		snprintf(string, 32, format, angle*180.0/M_PI);
		if ( strstr(string, "nan") ) snprintf(string, 32, ".");
		data = string_add(&data, string);
		if ( !item->data ) item->data = data;
		if ( item->maxlen < strlen(string) ) item->maxlen = strlen(string);
	}
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG item_put_float_list: " << tag << " " << *(item->data) << endl;
	
	return 0;
}

/**
@brief 	Lists all items associated with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param 	&tag		tag for items to be listed.
@return int 				0.

	The item with a given tag is listed to standard output as an end-of-line
	delimited array.

**/
int 		item_list(Bstar_old* star, Bstring& tag)
{
	Bstar_block*	block;
	Bstar_item* 	item;
	Bstring*		data;
	
	if ( verbose & VERB_LABEL )
		cout << "Listing all items associated with tag \"" << tag << "\"" << endl;

	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag ) {
				for ( data = item->data; data; data = data->next )
					cout << " " << *data;
			}
		}
	}
	cout << endl;
	
	return 0;
}

/**
@brief 	Deletes all items associated with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param	*tag		tag for items to be deleted.
@return int 				0.

	The item with a given tag is deleted in all blocks it is found
	and the item pointers are rearranged to fill in the gap.

**/
int 		item_delete_all(Bstar_old* star, const char* tag)
{
	Bstar_block*	block;
	
	if ( verbose & VERB_FULL )
		cout << "Deleting all items associated with tag \"" << tag << "\"" << endl;
	
	for ( block = star->block; block; block = block->next )
		item_delete_from_block(block, tag);
	
	return 0;
}

/**
@brief 	Deletes an item associated with a given tag from a block in the STAR data base.
@param	*block	block in the STAR database.
@param	*tag		tag for items to be deleted.
@return int 		0.

	The item with a given tag is deleted in the specified block it is found
	and the item pointers are rearranged to fill in the gap.

**/
int 		item_delete_from_block(Bstar_block* block, const char* tag)
{
	Bstar_item* 	item, *item2;
	
	if ( verbose & VERB_FULL )
		cout << "Deleting tag \"" << tag << "\" from block" << endl;
	for ( item = item2 = block->item; item;  ) {
		if ( item->tag == tag ) {
			if ( item == block->item ) {
				block->item = item2 = item->next;
				kill_item(item);
				item = item2;
			} else {
				item2->next = item->next;
				kill_item(item);
				item = item2->next;
			}
		} else {
			item2 = item;
			item = item->next;
		}
	}
	
	return 0;
}

/**
@brief 	Deletes all blocks with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param 	&tag		tag for blocks to be deleted.
@return int 		number of blocks deleted.

	The blocks containing a given tag are deleted and the
	block pointers are rearranged to fill in the gap.

**/
int 		block_delete(Bstar_old* star, Bstring& tag)
{
	int				del, ndel(0);
	Bstar_block*	block, *block2;
	Bstar_item* 	item;
	
	if ( verbose & VERB_LABEL )
		cout << "Deleting all blocks associated with tag \"" << tag << "\"" << endl;

	for ( block = block2 = star->block; block;  ) {
		for ( del = 0, item = block->item; item && del == 0; item = item->next )
			if ( item->tag == tag ) del = 1;
		if ( del ) {
			if ( block == star->block ) {
				star->block = block2 = block->next;
				kill_block(block);
				block = block2;
			} else {
				block2->next = block->next;
				kill_block(block);
				block = block2->next;
			}
			ndel++;
		} else {
			block2 = block;
			block = block->next;
		}
	}

	if ( verbose & VERB_LABEL )
		cout << "Blocks deleted:                 " << ndel << endl;
	
	return ndel;
}

/**
@brief 	Scales and shifts all items associated with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param 	&tag		tag for items to be modified.
@param 	iscale 		multiplier.
@param 	ishift 		value added.
@return int 		total number of values changed.

	The item must be integer and is modified as:
		new_value = old_value*scale + shift.

**/
int 		item_integer_scale_shift(Bstar_old* star, Bstring& tag, int iscale, int ishift)
{
	int 			idata, total(0);
	char			format[32];
	Bstar_block*	block;
	Bstar_item* 	item;
	Bstring*		data;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Modifying all integer items associated with tag \"" << tag << "\"" << endl;
		cout << "Scale and shift:                " << iscale << " " << ishift << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Modifying all integer items associated with tag \"" << tag << "\"" << endl;
		
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag ) {
				snprintf(format, 32, "%%%dd", item->maxlen);
				for ( data = item->data; data; data = data->next ) {
					idata = data->integer()*iscale + ishift;
					data[0] = Bstring(idata, format);
					total++;
				}
			}
		}
	}
	
	return total;
}

/**
@brief 	Scales and shifts all items associated with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param 	&tag		tag for items to be modified.
@param 	scale 		multiplier.
@param 	shift 		value added.
@return int 		total number of values changed.

	The item must be numeric and is modified as:
		new_value = old_value*scale + shift.

**/
int 		item_float_scale_shift(Bstar_old* star, Bstring& tag, float scale, float shift)
{
	int 			total(0), type, len;
	float			fdata, maxpre;
	char			format[32];
	Bstar_block*	block;
	Bstar_item* 	item;
	Bstring*		data;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Modifying all floating point items associated with tag \"" << tag << "\"" << endl;
		cout << "Scale and shift:                " << scale << " " << shift << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Modifying all floating point items associated with tag \"" << tag << "\"" << endl;
	
	
	for ( block = star->block; block; block = block->next ) {
		for ( item = block->item; item; item = item->next ) {
			if ( item->tag == tag ) {
				type = item_get_format(item, format);
				if ( type < 2 ) snprintf(format, 32, "%%%d.0f", item->maxlen);
				maxpre = pow(10.0, 1.0*format[1])*scale + shift;
				if ( maxpre < 1 ) maxpre = 1;
				format[1] = (char) (log(maxpre)/log(10.0) + 1) + '0';
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG item_float_scale_shift: format = " << format << endl;
				for ( data = item->data; data; data = data->next ) {
					fdata = data->real()*scale + shift;
					data[0] = Bstring(fdata, format);
					len = (*data).length();
					if ( item->maxlen < len ) item->maxlen = len;
					total++;
				}
			}
		}
	}
	
	return total;
}

/**
@brief 	Assigns a loop identification number to items.
@param	*block 		block in the STAR database.
@param 	loop		loop identifier to use.
@param 	n			number of patterns to test for.
@param 	... (tag_pattern)	tag pattern to match to set the loop identifier.
@return int 		number assigned.

	The items in the data block are rearranged so that the item assigned to 
	the loop follows the other loop items.

**/
int 		loop_set_identifier(Bstar_block* block, int loop, int n, ...)
{
	va_list			ap;
	char*			tag_pattern;
	int				i, nloop(0);
	Bstar_item* 	item;
	
	va_start(ap, n);
	for ( i=0; i<n; i++ ) {
		tag_pattern = va_arg(ap, char *);
		for ( item = block->item; item; item = item->next ) {
			if ( strstr(item->tag.c_str(), tag_pattern) == item->tag.c_str() ) {
				item->loop = loop;
				nloop++;
			}
		}
	}
	va_end(ap);
	
	if ( verbose & VERB_DEBUG_STAR )
		cout << "DEBUG loop_set_identifier: loop=" << loop << " number=" << nloop << endl;

	return i;
}

/**
@brief 	Gets the format from the item.
@param	*item 		STAR item.
@param 	*format		pointer to pre-allocated format string (modified).
@return int 		data type: 0=string, 1=integer, 2=float.

	Returns the format in the given format string.

**/
int			item_get_format(Bstar_item* item, char* format)
{
	int				type(0), pre, post, maxpre(0), maxpost(0);
	Bstring*		data;
	
	if ( item->type == StringItem )
		strcpy(format, "%s");
	else {
		type = 1;
		for ( data = item->data; data; data = data->next ) {
			if ( data->find(".") > -1 && data->length() > 1 ) {
				type = 2;
				pre = (data->pre('.')).length();
				post = (data->post('.')).length();
				if ( maxpre < pre ) maxpre = pre;
				if ( maxpost < post ) maxpost = post;
			}
			if ( item->maxlen < data->length() ) item->maxlen = data->length();
		}
		if ( type == 2 )
			snprintf(format, 32, "%%%d.%df", maxpre + maxpost + 1, maxpost);
		else
			snprintf(format, 32, "%%%dd", item->maxlen);
	}
	
	return type;
}

