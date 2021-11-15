/**
@file	ps_tree.cpp
@brief	Postscript tree output functions.
@author Bernard Heymann
@date	Created: 20010722
@date	Modified: 20120416
**/
 
#include "ps_plot.h"
#include "ps_tree.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Draws a tree into a postscript file.
@param 	&filename	postscript file name.
@param 	*root		root node.
@param 	scale 	 	scale for branch lengths.
@return int 			0.

Requirements:
	The tree must be well-formed, i.e., it must have branch lengths 
	and angles defined.

**/
int 		ps_draw_tree(Bstring& filename, Bnode* root, float scale)
{
	ofstream*	fps = ps_open_and_init(filename, filename, 1, 612, 792);
	
	if ( verbose & VERB_LABEL )
		cout << "Writing file:                   " << filename << endl;
	
	*fps << "%%%%Page: 1 1" << endl;

	*fps << "/Helvetica findfont 10 scalefont setfont" << endl;
	*fps << "1 setlinecap\n1 setlinejoin\n2 setlinewidth" << endl;
	
	ps_draw_node(fps, root, scale, 300, 400);

	*fps << "showpage" << endl;
	
	ps_close(fps);
	
	return 0;
}

/**
@brief 	Draws a tree node into a postscript file.
@param 	*fps		file pointer.
@param 	*node 		node.
@param 	scale 		scale for branch lengths
@param 	x 			x-coordinate for the node.
@param 	y 			y-coordinate for the node.
@return int 			0.

Requirements:
	The tree must be well-formed, i.e., it must have branch lengths 
	and angles defined.

**/
int 		ps_draw_node(ofstream* fps, Bnode* node, float scale, float x, float y)
{
	float		dx = scale*node->length*cos(node->angle);
	float		dy = scale*node->length*sin(node->angle);
	
	*fps << x << " " << y << " moveto ";
	x += dx;
	y += dy;
	*fps << x << " " << y << " lineto stroke" << endl;
	
	if ( node->two ) {
		ps_draw_node(fps, node->two, scale, x, y);
		ps_draw_node(fps, node->three, scale, x, y);
	} else {
		if ( node->angle < M_PI_2 || node->angle > 3*M_PI_2 )
			*fps << x + 5 << " " << y << " moveto ";
		else
			*fps << x - 5 << " (" << node->label << ") stringwidth pop sub " << y << " moveto ";
		*fps << "(" << node->label << ") show" << endl;
	}
	
	if ( node->length == 0 )		// Third branch of the root node
		ps_draw_node(fps, node->one, scale, x, y);
	
	return 0;
}



