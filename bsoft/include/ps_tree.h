/**
@file	ps_tree.h 
@brief	Header file for writing a tree into a postcsript file.
@author Bernard Heymann 
@date	Created: 20010722
@date	Modified: 20120416
**/
 
#include "rwtree.h"
#include "ps_plot.h"
#include "Bstring.h"

// Function prototypes
int 		ps_draw_tree(Bstring& filename, Bnode* root, float scale);
int 		ps_draw_node(ofstream* fps, Bnode* node, float scale, float x, float y);
