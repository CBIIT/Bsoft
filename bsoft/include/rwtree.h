/**
@file	rwtree.h 
@brief	Header file for reading (and writing) trees.
@author Bernard Heymann 
@date	Created: 20010722
@date	Modified: 20010805
**/
 

#ifndef _Bnode_
/************************************************************************
@Object: struct Bnode
@Description:
	Unrooted tree node structure.
@Features:
	An unrooted tree consists of a set of trifurcated internal nodes and
	a set of end or tip nodes. The choice of the root node is arbitrary,
	although it is always defined as the entry point into the tree.
	The branch length is the distance from the node closer to the root
	of the tree (thus inherent directionality which must be reset if
	the tree is rerooted). The plot angle is just a convenient field
	for tree display.
*************************************************************************/
struct Bnode {
	Bnode *one, *two, *three; 	// Node pointers, two are 0 for tip nodes
	int number;					// Number of tip nodes in branch, 1 for a tip node
	char label[100];			// Node label, '.' denotes an internal node
	float length; 				// Branch length
	float angle;				// Plot angle
};
#define _Bnode_
#endif

// Function prototypes
Bnode* 		read_tree(char* filename);
Bnode*		tree_read_node(char** handle, Bnode* parent);
int 		tree_calculate_angle(Bnode* node, int* itip, float tip_angle);
int 		tree_delete_node(Bnode* node);
int 		tree_rotate(Bnode* node, float angle);
int 		tree_show(Bnode* node);
