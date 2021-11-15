/**
@file	linked_list.h
@brief	Header file for generalized linked list functions
@author Bernard Heymann
@date	Created: 20031203
@date	Modified: 20150206
**/

#include <stdlib.h>
#include <string.h>

//Function prototypes
char*		add_item(char** list, unsigned long size);
char*		append_item(char** list, char* item);
char*		copy_item(char* toitem, char* fromitem, unsigned long size);
char*		remove_item(char** list, char* item, unsigned long size);
char*		replace_item(char** list, char* item, char* new_item);
char*		copy_list(char* list, unsigned long size);
long 		kill_list(char* list, unsigned long size);
long 		count_list(char* list);
long 		reverse_list(char** list);

