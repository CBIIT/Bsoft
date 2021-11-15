/**
@file	linked_list.cpp
@brief	Generalized linked list functions
@author Bernard Heymann
@date	Created: 20031203
@date	Modified: 20150206
**/

#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Adds an item to a linked list.
@param 	**list			pointer to first item in the list.
@param 	size			size of item.
@return char* 			new item.

	If the list is not defined, the new item becomes the first in the list.
	Otherwise, the list is traversed to find the end and the new item appended.
	Any structure with a pointer to itself as a first element can be used in a
	linked list. However, a linked list can only consist of one type of structure.

**/
char*		add_item(char** list, unsigned long size)
{
	char**		curr_item = (char **) *list;
	char*		new_item = new char[size];
	char*		next_item = NULL;
	
	memset(new_item, 0, size);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG add_item: " << &list << " " << &curr_item << " " << &new_item << endl;
	
	if ( !curr_item )
		*list = new_item;
	else {
		while ( ( next_item = *curr_item ) ) curr_item = (char **) next_item;
		*curr_item = new_item;
	}
	
	return new_item;
}

/**
@brief 	Appends an item to a linked list.
@param 	**list			pointer to first item in the list.
@param 	*item			item to append.
@return char* 			item.

	If the list is not defined, the new item becomes the first in the list.
	Otherwise, the list is traversed to find the end and the new item appended.
	Any structure with a pointer to itself as a first element can be used in a
	linked list. However, a linked list can only consist of one type of structure.

**/
char*		append_item(char** list, char* item)
{
	char**		curr_item = (char **) *list;
	char*		next_item = NULL;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG append_item: " << &list << " " << &item << endl;
	
	if ( !curr_item )
		*list = item;
	else {
		while ( ( next_item = *curr_item ) ) curr_item = (char **) next_item;
		*curr_item = item;
	}
	
	return item;
}

/**
@brief 	Copies an item to an existing item without changing the link.
@param 	*toitem			item to copy to.
@param 	*fromitem		item to copy from.
@param 	size			size of item.
@return char* 			the resultant item.

	The items must be the same type.

**/
char*		copy_item(char* toitem, char* fromitem, unsigned long size)
{
	char*		link = *((char **) toitem);
	
	memcpy(toitem, fromitem, size);
	
	memcpy(toitem, &link, sizeof(char*));
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG copy_item: %" << &toitem << " " << &link << " " << &(*((char **) toitem)) << endl;
	
	return toitem;
}

/**
@brief 	Finds the given item and deletes it from the linked list.
@param 	**list			pointer to first item in the linked list.
@param 	*item			item to be deleted.
@param 	size			size of item.
@return char* 			item after the one removed.

	If the item is the first in the list, the list pointer is set to point
	to the next item.
	Otherwise, the list is traversed to find the item, the previous item's
	pointer is set to the next item, and the current item deallocated.

**/
char*		remove_item(char** list, char* item, unsigned long size)
{
	if ( !(*list) || !item ) return NULL;
	
	char**		curr_item = (char **) *list;	// First item pointer
	char**		next_item = (char **) item;		// Item to be removed pointer
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG remove_item: " << &list << " " << &curr_item << " " << &item << endl;
		cout << "DEBUG remove_item: " << list << " " << curr_item << " " << item << endl;
	}
	
	if ( *list == item ) {
		*list = *curr_item;
	} else {
		for ( ; curr_item && *curr_item != item; curr_item = (char **) *curr_item ) ;
		*curr_item = *next_item;
	}
	
	delete[] item;
	
	return *curr_item;
}

/**
@brief 	Finds the given item and replaces it with the new item.
@param 	**list			pointer to first item in the linked list.
@param 	*item			item to be replaced.
@param 	*new_item		new item.
@return char* 			new item.

	If the item is the first in the list, the list pointer is set to point
	to the new item.
	Otherwise, the list is traversed to find the item, the previous item's
	pointer is set to the new item, the new item's pointer is set to that
	of the old item, and the old item is deallocated.

**/
char*		replace_item(char** list, char* item, char* new_item)
{
	if ( !(*list) || !item ) return NULL;
	
	char**		curr_item = (char **) *list;	// First item pointer
//	char**		next_item = (char **) item;		// Item to be replaced pointer
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG replace_item: " << &list << " " << &curr_item << " " << &item << " " << &new_item << endl;
		cout << "DEBUG replace_item: " << list << " " << curr_item << " " << item << " " << new_item << endl;
	}
	
	if ( *list == item ) {
		*list = new_item;
	} else {
		for ( ; curr_item && *curr_item != item; curr_item = (char **) *curr_item ) ;
		*curr_item = new_item;
	}
	
	*new_item = *item;
	
	delete[] item;
	item = NULL;
	
	return new_item;
}

/**
@brief 	Generates a complete copy of a linked list.
@param 	*list		linked list.
@param 	size		size of list item.
@return char*		new list.
**/
char*		copy_list(char* list, unsigned long size)
{
	if ( !list ) return 0;
	
	char**		item;
	char*		newlist = NULL;
	char*		newitem = NULL;
	
	for ( item = (char **) list; item; item = (char **) *item ) {
		newitem = add_item(&newitem, size);
		if ( !newlist ) newlist = newitem;
		memcpy(newitem, item, size);
		memset(newitem, 0, sizeof(char *));
	}
			
	return newlist;
}
	

/**
@brief 	Frees all the items in a linked list.
@param 	*list		first item in the linked list.
@param 	size		size of item.
@return long 		number of items deallocated.

	The list is traversed, setting a pointer to the next item before
	deallocating the current item.

**/
long 		kill_list(char* list, unsigned long size)
{
	if ( !list ) return 0;
	
	long		n(0);
	char**		item;
	
	for ( item = (char **) list; item; n++ ) {
		list = *item;
		delete[] item;
		item = (char **) list;
	}
			
	return n;
}

/**
@brief 	Counts the number of items in a linked list.
@param 	*list		first item in the linked list.
@return long 		number of items in the list.
**/
long 		count_list(char* list)
{
	if ( !list ) return 0;
	
	long		n(0);
	char**		item;
	
	for ( item = (char **) list; item; item = (char **) *item ) n++;
			
	return n;
}

/**
@brief 	Reverse the order of items in a linked list.
@param 	*list		first item in the linked list.
@return long 		number of items in the list.
**/
long 		reverse_list(char** list)
{
	if ( !(*list) ) return 0;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reverse_list: " << &list << endl;

	long		n(0);			
	char**		item;
	char**		itemp = NULL;	
	char*		link = NULL;

	for ( n=0, item = (char **) *list; item; item = (char **) link, n++ ) {
		link = *((char **) item);
		if ( *list == (char *) item ) *item = NULL;
		else *item = (char *) itemp;
		itemp = item;
	}
	
	*list = (char *) itemp;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG reverse_list: " << n << " items reversed" << endl;
	
	return n;
}


