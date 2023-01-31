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

int 		writePostScriptSlice(ofstream* fimg, Bimage* p, long n, long z, double tx, double ty, long m);

/**
@brief	Writing a postscript image.
@param	*p			the image structure.
@return	int 			0.

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
	
	long			n, z, nc, nr;
	long			width(600), height(800);
	double			scale(1), tx, ty;
	
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
				nr = z/cols;
				nc = z - cols*nr;
				tx = nc*p->sizeX();
				ty = nr*p->sizeY();
		 		writePostScriptSlice(fimg, p, n, z, tx, ty, z);
			}
			*fimg << "showpage" << endl;
		*fimg << "grestore" << endl;
	}

	ps_close(fimg);
	
	return 0;
}

/**
@brief	Writing the slices arranged in a postscript image.
@param	*p			the image structure.
@param	&t			array of slice locations.
@return	int 			0.

Every sub-image is drawn on its own page
	One-dimensional images are written as graphs.
	Three-dimensional images are montaged with all the slices of each
	sub-image arrayed on a page.

**/
int 		writePostScriptImage(Bimage* p, vector<Vector3<double>>& t)
{
	p->change_type(UCharacter);
	
	long			i, n, z;
	long			width(600), height(800);
	double			xmin(1e30), xmax(-1e30);
	double			scale(0.8/sqrt(1.0*t.size()));
	
	for ( auto x: t ) {
		if ( xmin > x[0] ) xmin = x[0];
		if ( xmax < x[0] ) xmax = x[0];
	}

	for ( auto& x: t ) x = (x/(xmax - xmin) * p->sizeX() + width/2) * sqrt(t.size());

	ofstream*		fimg = ps_open_and_init(p->file_name(), p->file_name(), p->images(), width, height);
	
	*fimg << "%Command line: " << p->label() << endl;
	*fimg << "/Helvetica findfont 10 scalefont setfont" << endl;
	for ( i=n=0; n<p->images(); ++n ) {
		*fimg << "%%Page: Image " << n+1 << endl;
		*fimg << "gsave" << endl;
			if ( p->compound_type() == TRGB )
				*fimg << "/DeviceRGB setcolorspace" << endl;
			else
				*fimg << "/DeviceGray setcolorspace" << endl;
			*fimg << scale << " " << scale << " scale" << endl;
			*fimg << "/Helvetica findfont " << 10/scale << " scalefont setfont" << endl;
			for ( z=0; z<p->sizeZ(); ++z, ++i ) {
		 		writePostScriptSlice(fimg, p, n, z, t[i][0], t[i][1], n+1);
			}
		*fimg << "grestore" << endl;
	}

	*fimg << "showpage" << endl;

	ps_close(fimg);
	
	return 0;
}

/**
@brief	Writing an image slice to a postscript file.
@param	*fimg		output stream.
@param	*p			the image structure.
@param	n			sub-image.
@param	z			slice within the sub-image.
@param	tx			translation in x.
@param	ty			translation in y.
@param	m			slice number.
@return	int 			0.

Every sub-image is drawn on its own page
	One-dimensional images are written as graphs.
	Three-dimensional images are montaged with all the slices of each
	sub-image arrayed on a page.

**/
int 		writePostScriptSlice(ofstream* fimg, Bimage* p, long n, long z, double tx, double ty, long m)
{
	long			i, j(0), x, y, c;
	long			line_length(32);
	unsigned char*	data = p->data_pointer();
	
	cout << tx << tab << ty << endl;
	
	*fimg << "gsave" << endl;
	*fimg << tx << " " << ty << " translate" << endl;
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
	for ( y=0; y<p->sizeY(); y++ ) {
		for ( x=0; x<p->sizeX(); x++ ) {
			for ( c=0; c<p->channels(); c++ ) {
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
	*fimg << "1 1 0 setrgbcolor" << endl;
	*fimg << "10 10 moveto" << endl << "(" << m << ") show" << endl;
	*fimg << "grestore" << endl;

	return 0;
}
