/**
@file	rwDX.cpp
@brief	Functions for reading and writing OpenDX files
@author Bernard Heymann
@date	Created: 20150319
@date	Modified: 20150319
**/

#include "rwDX.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Reading an OpenDX image format.
@param	*p			the image structure.
@param 	readdata	flag to activate reading of image data.
@return	int			error code (<0 means failure).
# Comments
          object 1 class gridpositions counts nx ny nz
          origin xmin ymin zmin
          delta hx 0.0 0.0
          delta 0.0 hy 0.0 
          delta 0.0 0.0 hz
          object 2 class gridconnections counts nx ny nz
          object 3 class array type double rank 0 times n
          u(0,0,0) u(0,0,1) u(0,0,2)
          ...
          u(0,0,nz-3) u(0,0,nz-2) u(0,0,nz-1)
          u(0,1,0) u(0,1,1) u(0,1,2)
          ...
          u(0,1,nz-3) u(0,1,nz-2) u(0,1,nz-1)
          ...
          u(0,ny-3,nz-3) u(0,ny-2,nz-2) u(0,ny-1,nz-1)
          u(1,0,0) u(1,0,1) u(1,0,2)
          ...
          attribute "dep" string "positions"
          object "regular positions regular connections" class field
          component "positions" value 1
          component "connections" value 2
          component "data" value 3
**/
int 	readDX(Bimage* p, int readdata)
{
    ifstream		fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readDX: Start reading text image" << endl;
	
	Bstring			str, tag;
	long			i(0), j(0), n(0);
	Vector3<long>	size;
	Vector3<double>	ori, sam(1,1,1);
	
	// Set up default parameters
	p->size(1, 1, 1);
	p->images(1);
	p->channels(1);
	
	p->data_type(Float);

	while ( !fimg.eof() ) {
		fimg >> str;
		if ( tag == "counts" && n < 1 ) {
			if ( i < 3 ) size[i++] = str.integer();
			else {
				n = size[0]*size[1]*size[2];
				p->size(size);
				if ( readdata ) p->data_alloc();
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readDX: size=" << p->size() << endl;
			}
		}
		if ( tag == "origin" && i < 3 )
			ori[i++] = str.integer();
		if ( tag == "delta" && i < 3 ) {
			if ( i == j-1 ) sam[i] = str.real();
//			cout << tag << tab << i << tab << j << endl;
			i++;
		}
		if ( readdata && tag == "times" && i < n )
			p->set(i++, str.real());
		if ( str == "counts" || str == "origin" || str == "delta" || str == "times" ) {
			tag = str;
			i = 0;
			if ( str == "delta" ) j++;
			if ( str == "times" ) fimg >> str;
		}
	}
	
	p->origin(-ori);
	p->sampling(sam);
	if ( readdata ) {
		if ( p->sizeZ() > 1 ) str = "zyx";
		else if ( p->sizeY() > 1 ) str = "yxz";
		p->reslice(str);
	}
	
	fimg.close();
	
	return 0;
}

/**
@brief	Writing a DX image format.
@param	*p			the image structure.
@return	int					error code (<0 means failure).
**/
int 	writeDX(Bimage* p)
{
    ofstream        fimg;
    fimg.open(p->file_name());
    if ( fimg.fail() ) return -1;
	
	long			i, j, x, y, z, n = p->size()[0]*p->size()[1]*p->size()[2];
	
	fimg << "# Written by Bsoft" << endl;
	fimg << "# " << p->label() << endl;
	
	fimg << "object 1 class gridpositions counts " << p->sizeX()
		<< " " << p->sizeY() << " " << p->sizeZ() << endl;
	fimg << "origin " << -p->image->origin()[0]
		<< " " << -p->image->origin()[1] << " " << -p->image->origin()[2] << endl;
	fimg << "delta " << p->sampling(0)[0] << " 0.0 0.0" << endl;
    fimg << "delta 0.0 " << p->sampling(0)[1] << " 0.0" << endl;
    fimg << "delta 0.0 0.0 " << p->sampling(0)[2] << endl;

	fimg << "object 2 class gridconnections counts " << p->sizeX()
		<< " " << p->sizeY() << " " << p->sizeZ() << endl;
	
	fimg << "object 3 class array type double rank 0 times " << n;	
	for ( j=x=0; x<p->sizeX(); x++ ) {
		for ( y=0; y<p->sizeY(); y++ ) {
			for ( z=0; z<p->sizeZ(); z++, j++ ) {
				if ( j%3 ) fimg << " ";
				else fimg << endl;
				i = p->index(x,y,z);
				fimg << (*p)[i];
			}
		}
	}
	fimg << endl;

	fimg << "attribute \"dep\" string \"positions\"" << endl;
	fimg << "object \"regular positions regular connections\" class field" << endl;
	fimg << "component \"positions\" value 1" << endl;
	fimg << "component \"connections\" value 2" << endl;
	fimg << "component \"data\" value 3" << endl;
	
	fimg << endl;
	
	fimg.close();
	
	return 0;
}

