/**
@file	bstar.cpp
@brief	Examines STAR format files
@author Bernard Heymann
@date	Created: 20010409
@date	Modified: 20210214
**/

#include "star.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


int 		show_tags(Bstar2& star);
int 		item_delete_all(Bstar2& star, string& tag);
int 		item_scale_shift(Bstar2& star, string& tag, double scale, double shift, int flag);
int 		item_list(Bstar2& star, string& tag);


// Usage assistance
const char* use[] = {
" ",
"Usage: bstar [options] input.star [input.star]",
"----------------------------------------------",
"Examines and modifies STAR format files.",
" ",
"Actions:",
"-header 85               List all the header lines limited to the indicated length.",
"-split 3                 Split the data blocks into individual files, inserting digits before the extension.",
"-list time_unit          List all items with this STAR tag.",
"-delete resolution       Delete all items with this STAR tag.",
"-blockdelete symmetry    Delete all data blocks containing this STAR tag.",
"-integerscale select,1,3 Scale and shift all integer items with this STAR tag.",
"-floatscale psi,1.5,-5.1 Scale and shift all floating point items with this STAR tag.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-linelength 120          Output line length (default 100).",
" ",
"Output:",
"-output file.star        Output STAR format file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	long			header(0);				// Show header lines using this length
	int				split(0);				// Flag to split into blocks
	int				showflag(1);			// Default show tags
	string			tag_to_list;			// STAR item which should be listed
	string			tag_to_be_deleted;		// STAR item which should be deleted
	string			block_to_be_deleted;	// STAR item in block which should be deleted
	string			tag_integer;			// Integer STAR item which should be scaled and shifted
	int				iscale(1);
	int				ishift(0);
	string			tag_float;				// Floating point STAR item which should be scaled and shifted
	double			scale(1);
	double			shift(0);
	string			outstar;				// Output STAR format file
	
	// Initialize the STAR database
	Bstar2			star;
	star.line_length(120);				// Set the output line length
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	Bstring*		strarr;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "header" )
			if ( ( header = curropt->value.integer() ) < 1 )
				cerr << "-header: A line length must be specified!" << endl;
		strarr = curropt->value.split(",");
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("dat") ) split = 9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 4 ) split = 4;
		}
		if ( curropt->tag == "linelength" ) {
			if ( curropt->value.integer() < 1 )
				cerr << "-linelength: A line length must be specified!" << endl;
			else
				star.line_length(curropt->value.integer());
		}
		if ( curropt->tag == "list" ) {
			showflag = 0;
			tag_to_list = curropt->value.c_str();
			if ( tag_to_list.size() < 1 )
				cerr << "-list: A STAR tag must be specified!" << endl;
		}
		if ( curropt->tag == "delete" ) {
			showflag = 0;
			tag_to_be_deleted = curropt->value.c_str();
			if ( tag_to_be_deleted.length() < 1 )
				cerr << "-delete: A STAR tag must be specified!" << endl;
		}
		if ( curropt->tag == "blockdelete" ) {
			showflag = 0;
			block_to_be_deleted = curropt->value.c_str();
			if ( block_to_be_deleted.length() < 1 )
				cerr << "-blockdelete: A STAR tag must be specified!" << endl;
		}
		if ( curropt->tag == "integerscale" ) {
			showflag = 0;
			tag_integer = strarr->c_str();
			if ( tag_integer.length() < 1 )
				cerr << "-integerscale: A STAR tag must be specified!" << endl;
			else {
				if ( strarr->next ) iscale = strarr->next->integer();
				if ( strarr->next->next ) ishift = strarr->next->next->integer();
			}
		}
		if ( curropt->tag == "floatscale" ) {
			showflag = 0;
			tag_float = strarr->c_str();
			if ( tag_float.length() < 1 )
				cerr << "-floatscale: A STAR tag must be specified!" << endl;
			else {
				if ( strarr->next ) scale = strarr->next->real();
				if ( strarr->next->next ) shift = strarr->next->next->real();
			}
		}
		if ( curropt->tag == "output" )
			outstar = curropt->filename().c_str();
		string_kill(strarr);
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	// Read all the STAR format files
	while ( optind < argc )
		star.read(argv[optind++]);
	
	if ( header )
		star.list_comments(header);
	
	if ( tag_to_list.length() > 0 )
		item_list(star, tag_to_list);
	
	if ( tag_to_be_deleted.length() > 0 )
		item_delete_all(star, tag_to_be_deleted);
	
	if ( block_to_be_deleted.length() > 0 )
		star.erase(block_to_be_deleted);
	
	if ( tag_integer.length() > 0 )
		item_scale_shift(star, tag_integer, iscale, ishift, 0);
	
	if ( tag_float.length() > 0 )
		item_scale_shift(star, tag_float, scale, shift, 1);

	if ( showflag && verbose ) show_tags(star);

	// write an output STAR format file if a name is given
    if ( outstar.length() ) {
		star.comment(star.comment() + command_line_time().c_str());
		star.write(outstar, split);
	}
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

/**
@brief 	Prints the list of tags in the STAR data base.
@param	&star 		the STAR database
@return int 			number of tags.
**/
int 		show_tags(Bstar2& star)
{
	int 			ntags(0);
	string			tag;
	
	for ( auto ib: star.blocks() ) {
		cout << endl << "data_" << ib.tag() << endl;
		cout << endl << "loop_" << endl;
		cout << "_star_tag" << endl;
		cout << "_star_type" << endl;
		for ( auto it: ib.items() ) {
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG show_tags: tag = " << it.first << tab << it.second << endl;
			tag = "\"_" + it.first + "\"";
			cout << left << setw(40) << tag;
			if ( check_for_number(it.second) )
				cout << "Number" << endl;
			else
				cout << "String" << endl;
			ntags++;
		}
		for ( auto il: ib.loops() ) {
			map<string, int>&	t = il.tags();
			for ( auto it: t ) {
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG show_tags: loop tag = " << it.first << tab << it.second << endl;
				tag = "\"_" + it.first + "\"";
				cout << left << setw(40) << tag;
				if ( check_for_number(il[0][it.second]) )
					cout << "Number" << endl;
				else
					cout << "String" << endl;
				ntags++;
			}
		}
		cout << endl;
	}
	cout << endl;
	
	return ntags;
}

/**
@brief 	Deletes all items associated with a given tag from the STAR data base.
@param	*star 		the STAR database.
@param	*tag		tag for items to be deleted.
@return int 				0.

	The item with a given tag is deleted in all blocks it is found
	and the item pointers are rearranged to fill in the gap.

**/
int 		item_delete_all(Bstar2& star, string& tag)
{
	if ( verbose & VERB_FULL )
		cout << "Deleting all items associated with tag \"" << tag << "\"" << endl;
	
	for ( auto ib: star.blocks() ) {
		if ( ib.exists(tag) ) ib.erase(tag);
		else for ( auto il: ib.loops() ) {
			if ( il.find(tag) >= 0 ) il.erase(tag);
		}
	}
	
	return 0;
}

/**
@brief 	Gets the format from the item.
@param	&s			string.
@param 	*format		pointer to pre-allocated format string (modified).
@return int 			data type: 0=string, 1=integer, 2=float.

	Returns the format in the given format string.

**/
int			string_get_format(string& s, char* format)
{
	int				type(0), pre, post;
	
	if ( check_for_number(s) ) {
		auto it = s.find('.');
		if ( it == string::npos ) {
			snprintf(format, 32, "%%%dd", 10);
		} else {
			pre = it;
			post = s.size() - it - 1;
			snprintf(format, 32, "%%%d.%df", pre + post + 1, post);
		}
	} else {
		strcpy(format, "%s");
	}
		
	return type;
}

/**
@brief 	Scales and shifts all items associated with a given tag from the STAR data base.
@param	&star 		the STAR database.
@param 	&tag		tag for items to be modified.
@param 	scale 		multiplier.
@param 	shift 		value added.
@param 	flag 		0=integer, 1=real.
@return int 			total number of values changed.

	The item must be integer and is modified as:
		new_value = old_value*scale + shift.

**/
int 		item_scale_shift(Bstar2& star, string& tag, double scale, double shift, int flag)
{
	long 			t(flag+1), j, idata, rdata, total(0);
//	char			format[32];

	if ( verbose & VERB_PROCESS ) {
		cout << "Modifying all integer items associated with tag \"" << tag << "\"" << endl;
		cout << "Scale and shift:                " << scale << " " << shift << endl;
	} else if ( verbose & VERB_LABEL )
		cout << "Modifying all integer items associated with tag \"" << tag << "\"" << endl;

//	if ( flag ) snprintf(format, 32, "%g");
//	else snprintf(format, 32, "%d");

	for ( auto ib: star.blocks() ) {
		for ( auto it: ib.items() ) {
			if ( it.first == tag ) {
//				t = check_for_type(it.second);
				if ( t == 1 ) {
					idata = to_integer(it.second)*scale + shift;
					it.second = to_string(idata);
				} else if ( t == 2 ) {
					rdata = to_real(it.second)*scale + shift;
					it.second = to_string(rdata);
				}
				if ( t ) total++;
			}
		}
		for ( auto il: ib.loops() ) {
			if ( ( j = il.find(tag) ) >= 0 ) {
//				t = check_for_type(il.data()[j][0]);
				for ( auto ir: il.data() ) {
					if ( t == 1 ) {
						idata = to_integer(ir[j])*scale + shift;
						ir[j] = to_string(idata);
					} else if ( t == 2 ) {
						rdata = to_real(ir[j])*scale + shift;
						ir[j] = to_string(rdata);
					}
					if ( t ) total++;
				}
			}
		}
	}

	return total;
}


/**
@brief 	Lists all items associated with a given tag from the STAR data base.
@param	&star 		the STAR database.
@param 	&tag		tag for items to be listed.
@return int 			0.

	The item with a given tag is listed to standard output as an end-of-line
	delimited array.

**/
int 		item_list(Bstar2& star, string& tag)
{
	long			j;
	
	if ( verbose & VERB_LABEL )
		cout << "Listing all items associated with tag \"" << tag << "\"" << endl;

	for ( auto ib: star.blocks() ) {
		for ( auto it: ib.items() ) {
			if ( it.first == tag )
				cout << " " << it.second;
		}
		for ( auto il: ib.loops() ) {
			if ( ( j = il.find(tag) ) >= 0 ) {
				for ( auto ir: il.data() ) {
					cout << " " << ir[j];
				}
			}
		}
	}

	cout << endl;
	
	return 0;
}
