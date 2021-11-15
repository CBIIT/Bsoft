/**
@file	rwtree.cpp
@brief	Functions for reading (and writing) trees.
@author Bernard Heymann
@date	Created: 20010722
@date	Modified: 20110810
**/
 
#include "rwtree.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a phylip format tree file.

	Reads the tree file using a recursive approach.

@param 	*filename		tree file name.
@return Bnode* 				a pointer to the root node of the tree.
**/
Bnode* 		read_tree(char* filename)
{
	FILE        *ftree;
	if ( ( ftree = fopen(filename, "r") ) == NULL ) return NULL;
	
	if ( verbose & VERB_LABEL )
		cout << "Reading file:                   " << filename << endl;
	
	fseek(ftree, 0, SEEK_END);
	long		filesize = ftell(ftree);
	fseek(ftree, 0, SEEK_SET);
	
	char*		buffer = (char *) malloc(filesize);
	if ( fread( buffer, filesize, 1, ftree ) < 1 ) {
		fclose(ftree);
		return NULL;
	}
	
	fclose(ftree);
		
	char*			aptr = buffer;
	while ( *aptr != '(' ) aptr++;
	unsigned int	i, j(0);
	for ( i=0; i<strlen(aptr); i++ ) if ( !isspace(aptr[i]) ) {
		aptr[j] = aptr[i];
		j++;
	}
	aptr[j] = 0;

	if ( verbose & VERB_DEBUG )
		cout << buffer << endl << endl;
	
	// The root node is trifurcated because the tree is "unrooted"
	Bnode*			root = tree_read_node(&aptr, NULL);
	
	int				itip(0);
	tree_calculate_angle(root, &itip, M_PI*2.0/root->number);
	
	if ( verbose & VERB_PROCESS ) {
		tree_show(root);
		cout << "Number of tree tips:            " << root->number << endl;
	}
	
	return root;
}

/**
@brief 	Reads a node into a tree structure.
@param	**handle		pointer to the buffer holding the phylip format tree.
@param 	*one			one node.
@return Bnode* 				a pointer to the new node of the tree.
**/
Bnode*		tree_read_node(char** handle, Bnode* one)
{
	Bnode*		node = new Bnode;
	memset(node, 0, sizeof(Bnode));

	node->one = one;
	node->label[0] = '.';
	int 		i;
	char		astring[128];
	
	// Every internal node starts with a '(' and is delimited by a ','
	if ( **handle == '(' ) {
		*handle += 1;
		node->two = tree_read_node(handle, node);
		*handle += 1;
		node->three = tree_read_node(handle, node);
		node->number += node->two->number + node->three->number;
		*handle += 1;
		if ( !one ) {
			node->one = tree_read_node(handle, node);
			node->number += node->one->number;
		}
	} else {
		node->number = 1;							// Count the tip node
		if ( **handle == ',' ) *handle += 1;
		for ( i=0; **handle != ':' && **handle != 0; i++, *handle += 1 )
				node->label[i] = **handle;			// Get the node label
	}
	
	*handle += 1;
	
	// Get the branch length
	if ( strlen(*handle) > 0 ) {
		for ( i=0; **handle != ',' && **handle != ')'; i++, *handle += 1 )
				astring[i] = **handle;
		astring[i] = 0;
		sscanf(astring, "%f", &node->length);
		if ( node->length == 0 && strlen(astring) > 0 ) node->length = 0.000001;
	}
	
	if ( **handle == ')' ) *handle += 1;
	
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG tree_read_node: " << &node << " " << &node->one << " " 
			<< &node->two << " " << &node->three << " " << node->number 
			<< " " << node->label << " " << node->length << endl;
	
	return node;
}

/**
@brief 	Calculates an angle for each node for plotting.
@param 	*node			node.
@param 	*itip			the number of tip nodes passed.
@param 	tip_angle 	the angular increment for each tip node.
@return int 				0.
**/
int 		tree_calculate_angle(Bnode* node, int* itip, float tip_angle)
{
	if ( node->two ) {
		tree_calculate_angle(node->two, itip, tip_angle);
		tree_calculate_angle(node->three, itip, tip_angle);
		node->angle = (node->two->angle + node->three->angle)/2;
	} else {
		node->angle = (*itip)*tip_angle;
		*itip += 1;
	}
	
	if ( node->length == 0 )		// Third branch of the root node
		tree_calculate_angle(node->one, itip, tip_angle);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG tree_calculate_angle: " << *itip << " " << &node << " " 
			<< &node->one << " " << &node->two << " " << &node->three << " " 
			<< node->number << " " << node->label << " " << node->length << " " << node->angle << endl;
	
	return 0;
}
	
/**
@brief 	Deletes a node with all daughter nodes.
@param 	*node			node.
@return int 				0.
**/
int 		tree_delete_node(Bnode* node)
{
	if ( node->two ) {
		tree_delete_node(node->two);
		tree_delete_node(node->three);
	}
	
	if ( node->length == 0 )		// Third branch of the root node
		tree_delete_node(node->one);
	
	delete node;
	
	return 0;
}

/**
@brief 	Adds an angle to all nodes.

	The angle must be specified in radians.

@param 	*node			node.
@param 	angle 		rotation angle (radians).
@return int 				0.
**/
int 		tree_rotate(Bnode* node, float angle)
{
	if ( node->length == 0 )		// Third branch of the root node
		tree_rotate(node->one, angle);
	
	if ( node->two ) {
		tree_rotate(node->two, angle);
		tree_rotate(node->three, angle);
	}
	
	node->angle += angle;
	
	return 0;
}

/**
@brief 	Shows a node with all daughter nodes.
@param 	*node			node.
@return int 				0.
**/
int 		tree_show(Bnode* node)
{
	cout << " - " << node->label << " (" << node->length << ")";
	
	if ( node->length == 0 )		// Third branch of the root node
		tree_show(node->one);
	
	if ( node->two ) {
		tree_show(node->two);
		tree_show(node->three);
	}
	
	if ( node->number < 2 ) cout << endl;
	
	return 0;
}
