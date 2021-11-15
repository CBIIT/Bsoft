/**
@file	rwPostScript.cpp
@brief	Functions for (reading and) writing postscript image files
@author Bernard Heymann
@date	Created: 20010614
@date 	Modified: 20171124
**/

#include "ps_plot.h" 
#include "rwPostScript.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Writing a postscript image format.
@param	*p			the image structure.
@return	int 		0.

Every sub-image is drawn on its own page
	One-dimensional images are written as graphs.
	Three-dimensional images are montaged with all the slices of each
	sub-image arrayed on a page. 

**/
int 		writePostScriptImage(Bimage* p)
{
	if ( p->sizeZ()*p->sizeY() < 2 ) {
		Bplot* 		plot = p->plot();
		ps_plot(p->file_name(), plot);
		delete plot;
		return 0;		
	}
	
	p->change_type(UCharacter);
	
	long			i = 0, j = 0, n, x, y, z, c, nc, nr;
	long			line_length = 32;
	long			width = 600, height = 800;
	float			scale = 1;
	unsigned char*	data = p->data_pointer();
	
	long	cols = (int) (sqrt((p->sizeZ() - 1)*width*1.0/height) + 1);
	long	rows = (int) ((p->sizeZ() - 1)*1.0/cols + 1);
	
	if ( cols*p->sizeX() > width ) scale = width*1.0/(cols*p->sizeX());
	if ( scale*rows*p->sizeY() > height ) scale = height*1.0/(rows*p->sizeY());
	
	if ( p->sizeZ()*p->sizeY() < 2 ) {
		width = 600;
		height = 800;
	}
	
	ofstream*		fimg = ps_open_and_init(p->file_name(), p->file_name(), p->images(), width, height);
	
	*fimg << "%Command line: " << p->label() << endl;
	*fimg << "/Helvetica findfont 10 scalefont setfont" << endl;
	for ( n=0; n<p->images(); n++ ) {
		*fimg << "%%Page: Image " << n+1 << endl;
		*fimg << "gsave" << endl;
			if ( p->compound_type() == TRGB )
				*fimg << "/DeviceRGB setcolorspace" << endl;
			else
				*fimg << "/DeviceGray setcolorspace" << endl;
			*fimg << scale << " " << scale << " scale" << endl;
			for ( z=0; z<p->sizeZ(); z++ ) {
				*fimg << "gsave" << endl;
				nr = z/cols;
				nc = z - cols*nr;
				*fimg << nc*p->sizeX() << " " << nr*p->sizeY() << " translate" << endl;
				*fimg << "<<" << endl;
				*fimg << "	/ImageType 1" << endl;
				*fimg << "	/Width " << p->sizeX() << " /Height " << p->sizeY() << endl;
				*fimg << "	/BitsPerComponent 8" << endl;
				if ( p->compound_type() == TRGB )
					*fimg << "	/Decode [0 1 0 1 0 1]" << endl;
				else
					*fimg << "	/Decode [0 1]" << endl;
				*fimg << "	/ImageMatrix [1 0 0 1 0 0]" << endl;
				*fimg << "	/DataSource currentfile /ASCIIHexDecode filter" << endl;
				*fimg << ">>" << endl;
				*fimg << "%%BeginData: " << p->channels()*p->sizeX()*p->sizeY() << " hex bytes" << endl;
				*fimg << "image" << endl;
				*fimg << hex << right << setfill('0');
				j = 0;
				for ( y=0; y<p->sizeY(); y++ ) {
					for ( x=0; x<p->sizeX(); x++ ) {
						for ( c=0; c<p->channels(); c++ ) {
//							i = (((n*p->z + z)*p->sizeY() + y)*p->sizeX() + x)*p->channels() + c;
							i = p->index(c,x,y,z,n);
							if ( j >= line_length ) {
								*fimg << endl;
								j = 0;
							}
							*fimg << setw(2) << (int)data[i];
							j++;
						}
					}
				}
				*fimg << ">";
				if ( j <= line_length ) *fimg << endl;
				*fimg << dec << "%%EndData" << endl;
				if ( p->sizeZ() > 1 ) {
					*fimg << "1 1 0 setrgbcolor" << endl;
					*fimg << "0 0.9 moveto" << endl << "(" << z << ") show" << endl;
				}
				*fimg << "grestore" << endl;
			}
			*fimg << "showpage" << endl;
		*fimg << "grestore" << endl;
	}

	ps_close(fimg);
	
	return 0;
}

