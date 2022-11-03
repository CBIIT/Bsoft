/**
@file	options.cpp
@brief	Functions to handle general options 
@author Bernard Heymann 
@date	Created: 20010613
@date	Modified: 20220722
**/
 
#include "options.h" 
#include "utilities.h"

#define MAX_TAG_LEN	20

// Declaration of global variables
extern int 		verbose;		// Level of output to the shell
extern string	command;		// Command line
extern int		thread_limit;	// Thread limit

/**
@brief 	Parses command line arguments based on a template.
@param 	*use[]			usage list of strings, the template.
@param 	argc			number of command line arguments.
@param 	*argv[]			array of command line argument strings.
@param 	&optind			first argument after option list.
@return Boption* 		linked list of tag-value pairs.

	The usage list is parsed to find the desired option tag.
	Each option in the usage list must have the following format:
		-<tag> <value,value,...>
	The first character on the line must be '-', the tag may not
	extend beyond the 17'th character and the value must start
	before the 19'th character.
	A partial input tag is tolerated as long as it is unambiguous.
	An ambiguous or unknown tag or a tag requiring a value but
	without one causes program abortion.
	Special tags:
		-verbose 3      sets the verbosity level for all programs.
		-help           returns the usage information and quits.

**/
Boption*	get_option_list(const char* use[], int argc, char* argv[], int& optind)
{
	if ( argc < 2 ) {
		usage(use, 0);
		return NULL;
	}
	
	cmd_line(argc, argv);
	
	int			i, v(1);
	Bstring		usetag, tag, value;
	Boption*	option = NULL;
	Boption*	curropt = NULL;
	Boption*	nuopt = NULL;
	
	for ( i=1; i<argc; i++ ) {
		if ( strncmp(argv[i], "-verbose", strlen(argv[i])) == 0 && i<argc-1 ) {
			if ( argv[i+1][0] != '-' )
				verbose = get_option_verbose(argv[i+1]);
		}
		if ( strncmp(argv[i], "-threads", strlen(argv[i])) == 0 ) {
			thread_limit = atoi(argv[++i]);
			if ( thread_limit < 1 )
				cerr << "Error: At least one thread must be specified!" << endl;
		}
		if ( strncmp(argv[i], "-help", 5) == 0 ) {
			usage(use, 1);
			return NULL;
		}
	}

	if ( verbose ) {
		Bstring ct = command_line_time();
		cout << ct << endl;
	}

	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG get_option_list:";
		for ( i=0; i<argc; i++ ) cout << " " << argv[i];
		cout << endl;
	}
	
	optind = 1;
	for ( i=1; i<argc && v; i+=v ) {
//		cout << i << tab << argv[i] << endl;
		if ( argv[i][0] == '-' ) {
			nuopt = new Boption(argv[i], argv[i+1], use);
			if ( !option ) option = curropt = nuopt;
			else {
				curropt->next = nuopt;
				curropt = curropt->next;
			}
			if ( curropt->value.length() ) v = 2;
			else v = 1;
			optind = i + v;
			if ( verbose & VERB_DEBUG ) {
				if ( v == 1 ) cout << "DEBUG get_option_list: tag=" << curropt->tag << endl;
				else cout << "DEBUG get_option_list: tag=" << curropt->tag << " value=" << curropt->value << endl;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG get_option_list: optind=" << optind;
		if ( optind < argc ) cout << " arg=" << argv[optind] << endl;
		cout << endl;
	}
	
	if ( option->errors() ) {
 		cerr << "Option errors: " << option->errors() << endl;
		bexit(-1);
	}
	
	return option;
}

/**
@brief 	Deallocates a linked list of option structures.
@param 	*option		linked list of tag-value pairs.
@return int			0.
**/
int			option_kill(Boption* option)
{
	Boption		*curropt, *curropt2;
	
	for ( curropt = option; curropt; ) {
		curropt2 = curropt->next;
		curropt->tag = 0;
		curropt->value = 0;
		delete curropt;
		curropt = curropt2;
	}
	
	return 0;
}

/**
@brief 	Sets the verbose option.
@param 	&optarg		verbosity level.
@return int			verbosity number, <0 on error.

	The verbosity level is defined by the following constants setting 
	particular bits:
		VERB_NONE		0		No output
		VERB_RESULT		1		Program results
		VERB_LABEL		2		Function information
		VERB_PROCESS	4		Selected processing information
		VERB_STATS		8		Statistical information on objects
		VERB_FULL		16		All processing information
		VERB_TIME		32		Timing information
		VERB_MEMORY 	64		Memory allocation and freeing
		VERB_DEBUG		256 		Debugging information
		VERB_DEBUG_STAR	512 		STAR code debugging information
		VERB_DEBUG_DM	1024 	Digital Micrograph format debugging information
		VERB_DEBUG_ND2	2048 	Nikon ND2 format debugging information

**/
int 		get_option_verbose(char* optarg)
{
	Bstring		arg(optarg);
	return get_option_verbose(arg);
}

int 		get_option_verbose(Bstring& optarg)
{
	verbose = optarg.integer();
	
	Bstring			arg = optarg.lower();
	
	if ( arg.contains("non") ) verbose |= VERB_NONE;
	if ( arg.contains("res") ) verbose |= VERB_RESULT;
	if ( arg.contains("lab") ) verbose |= VERB_LABEL;
	if ( arg.contains("pro") ) verbose |= VERB_PROCESS;
	if ( arg.contains("sta") ) verbose |= VERB_STATS;
	if ( arg.contains("ful") ) verbose |= VERB_FULL;
	if ( arg.contains("tim") ) verbose |= VERB_TIME;
	if ( arg.contains("mem") ) verbose |= VERB_MEMORY;
	if ( arg.contains("deb") ) verbose |= VERB_DEBUG;
	if ( arg.contains("star") ) verbose |= VERB_DEBUG_STAR;
	if ( arg.contains("dm") ) verbose |= VERB_DEBUG_DM;
	if ( arg.contains("nd2") ) verbose |= VERB_DEBUG_ND2;
	if ( arg.contains("eer") ) verbose |= VERB_DEBUG_EER;
	
	return verbose; 
}

/**
@brief 	Returns the adjusted mass based on the added character.
@param 	&optarg		mass with added character.
@return double		adjusted mass.
**/
double		get_option_mass(Bstring& optarg)
{
	double			mass = optarg.real();
	
	if ( optarg.contains("k") || optarg.contains("K") )
		mass *= 1e3;
	else if ( optarg.contains("m") || optarg.contains("M") )
		mass *= 1e6;
	else if ( optarg.contains("g") || optarg.contains("G") )
		mass *= 1e9;
	
	return mass;
}


