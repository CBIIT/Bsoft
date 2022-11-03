/**
@file	rwstar.h
@brief	A header file for STAR format files
@author Bernard Heymann
@date	Created: 19990605
@date	Modified: 20190819
**/

#include "Bstring.h"

// Constants
#define LINELENGTH			80		// Length of line for output

#ifndef _STAR_OLD_
/************************************************************************
@Object: enum StarType
@Description:
	Types of elements in a STAR database.
@Features:
	Comments are ignored but passed to output files.
	Numbers are either integer or floating point.
*************************************************************************/
enum StarType {
	CommentItem = 0,
	StringItem = 1,
	NumberItem = 2
} ;

/************************************************************************
@Object: struct Bstar_item
@Description:
	Structure for a single entry or a single column in a STAR database.
@Features:
	Each item consists of a string or a zero-delimited sequence of strings.
*************************************************************************/
struct Bstar_item {
	Bstar_item*		next;			// Pointer to next item
	Bstring			tag;			// Item identifier
	StarType		type;			// Item type
	int 			loop; 			// Loop identifier
	int 			maxlen; 		// Maximum length of a string
	Bstring			comment;		// Comment following tag-value pair
	Bstring*		data;			// Pointer to numbers or strings
} ;

/************************************************************************
@Object: struct Bstar_block
@Description:
	Structure for a data block with multiple items in a STAR database.
@Features:
	Each data block contains a set set of unique items defined by tags 
	and with the associate data as single or multiple values.
	The order of the items in the data block are important and preserved.
*************************************************************************/
struct Bstar_block {
	Bstar_block*	next;			// Pointer to next block
	int				number;			// Block number (classic)
	Bstring			tag;			// Block identifier
	Bstring			comment;		// Comment following tag-value pair
	Bstring			filename;		// File containing this block
	Bstar_item* 	item;			// Pointer to items
} ;

/************************************************************************
@Object: struct Bstar_old
@Description:
	Overall STAR structure to hold all data blocks and options for I/O..
@Features:
	The split flag allows the user to output data blocks in separte files
	in stead of one big file.
	The line length field allows the user to output long lines without
	wrapping it around.
	The comments are ignored but output to the a new file - this can be used
	to document the history of the file.
	The STAR database is a hierarchy consisting of blocks, each with a set
	of items.
*************************************************************************/
struct Bstar_old {
	int 			split;			// Flag to indicate writing each data block in its own file
	int 			line_length;	// Length of lines in output file
	Bstring			comment;		// List of comments before 1st data block
	Bstar_block*	block;			// Pointer to first block
} ;
#define _STAR_OLD_
#endif

// Function prototypes
Bstar_old*		init_star();
int			kill_star(Bstar_old* star);
int 		kill_block(Bstar_block* block);
int 		kill_item(Bstar_item* item);
int			read_star(const char* filename, Bstar_old* star);
int			read_star(Bstring& filename, Bstar_old* star);
int			read_star(Bstring* file_list, Bstar_old* star);
Bstar_block* read_block(ifstream* fstar, char* aptr, Bstring& filename);
Bstar_item* read_single_item(ifstream* fstar, char* aline);
Bstar_item* read_loop_items(ifstream* fstar);
int			write_star(const char* filename, Bstar_old* star);
int			write_star(Bstring& filename, Bstar_old* star);
int			write_block(ofstream* fstar, Bstar_block* block, int linelength);
int 		star_update_comment(Bstar_old* star, int n, char** strings);
long		star_list_comments(Bstar_old* star, long len);
int			star_set_string_lengths(Bstar_old* star);
int 		item_change_tag(Bstar_old* star, const char* tag, const char* newtag);
int 		show_tags(Bstar_old* star);
int			item_index(Bstar_block* block, const char* tag);
Bstar_block*	block_find_with_tag(Bstar_old* star, const char* tag);
//Bstar_item* item_find(Bstar_block* block, const char* tag);
Bstar_item* item_find(Bstar_block* block, Bstring tag);
Bstar_item* item_find_or_make(Bstar_block* block, const char* tag);
long		item_get_number(Bstar_old* star, const char* tag);
long		item_get_number_for_block(Bstar_block* block, const char* tag);
char* 		item_get_string(Bstar_block* block, const char* tag);
int 		item_copy_string(char* string, Bstar_block* block, const char* tag);
int 		item_get_integer(Bstar_block* block, const char* tag);
float 		item_get_float(Bstar_block* block, const char* tag);
int 		item_put_string(Bstar_block* block, const char* tag, char* string);
int 		item_put_string(Bstar_block* block, const char* tag, Bstring& string);
int 		item_put_string_list(Bstar_block* block, const char* tag, Bstring* list);
int 		item_put_integer(Bstar_block* block, const char* tag, int value, const char* format);
int 		item_put_float(Bstar_block* block, const char* tag, float value, const char* format);
int 		item_put_list(Bstar_block* block, const char* tag, char* list, 
				size_t offset, const char* format);
int 		item_put_float_array(Bstar_block* block, const char* tag, int number, 
				float* value, const char* format);
int 		item_put_angle_list(Bstar_block* block, const char* tag, char* list, 
				size_t offset, const char* format);
int 		item_list(Bstar_old* star, Bstring& tag);
int 		item_delete_all(Bstar_old* star, const char* tag);
int 		item_delete_from_block(Bstar_block* block, const char* tag);
int 		block_delete(Bstar_old* star, Bstring& tag);
int 		item_integer_scale_shift(Bstar_old* star, Bstring& tag, int iscale, int ishift);
int 		item_float_scale_shift(Bstar_old* star, Bstring& tag, float scale, float shift);
int 		loop_set_identifier(Bstar_block* block, int loop, int n, ...);
int			item_get_format(Bstar_item* item, char* format);


